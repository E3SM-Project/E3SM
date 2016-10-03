"""Module containing tools for dealing with Fortran namelists.

The public interface consists of the following functions:
- `character_literal_to_string`
- `compress_literal_list`
- `expand_literal_list`
- `fortran_namelist_base_value`
- `is_valid_fortran_name`
- `is_valid_fortran_namelist_literal`
- `literal_to_python_value`
- `merge_literal_lists`
- `parse`
- `string_to_character_literal`

In addition, the `Namelist` class represents a namelist held in memory.

For the moment, only a subset of namelist syntax is supported; specifically, we
assume that only variables of intrinsic type are used, and indexing/co-indexing
of arrays to set a portion of a variable is not supported. (However, null values
and repeated values may be used to set or fill a variable as indexing would.)

We also always assume that a period (".") is the decimal separator, not a comma
(","). We also assume that the file encoding is UTF-8 or some compatible format
(e.g. ASCII).

Otherwise, most Fortran syntax rules implemented here are compatible with
Fortran 2008 (which is largely the same as previous standards, and will be
similar to Fortran 2015). The only exceptions should be cases where (a) we
deliberately prohibit "troublesome" behavior that would be allowed by the
standard, or (b) we rely on conventions shared by all major compilers.

One important convention is that newline characters can be used to denote the
end of a record. This makes them equivalent to spaces at most locations in a
Fortran namelist, except that newlines also end comments, and they are ignored
entirely within strings.

While the treatment of comments in this module is standard, it may be somewhat
surprising. Namelist comments are only allowed in two situations:

(1) As the only thing on a line (aside from optional indentation with spaces).
(2) Immediately after a "value separator" (the space, newline, comma, or slash
after a value).

This implies that all lines except for the last are syntax errors, in this
example:

```
&group_name! This is not a valid comment because it's after the group name.
foo ! Neither is this, because it's between a name and an equals sign.
= 2 ! Nor this, because it comes between the value and the following comma.
, bar = ! Nor this, because it's between an equals sign and a value.
2! Nor this, because it should be separated from the value by a comma or space.
bazz = 3 ! Nor this, because it comes between the value and the following slash.
/! This is fine, but technically it is outside the namelist, not a comment.
```

However, the above would actually be valid if all the "comments" were removed.
The Fortran standard is not clear about whether whitespace is allowed after
inline comments and before subsequent non-whitespace text (!), but this module
allows such whitespace, to preserve the sanity of both implementors and users.

The Fortran standard only applies to the interior of namelist groups, and not to
text between one namelist group and the next. This module assumes that namelist
groups are separated by (optional) whitespace and comments, and nothing else.
"""

###############################################################################
#
# Lexer/parser design notes
#
# The bulk of the complexity of this module is in the `_NamelistParser` object.
# Lexing, parsing, and translation of namelist data is all performed in a single
# pass (though it would be possible to use separate stages if needed). The style
# is that of a recursive descent parser, i.e. the functions correspond roughly
# to concepts in the Fortran namelist grammar, and top-down parsing is used.
# Parsing is done left-to-right with no backtracking.
#
# The most important attributes of a `_NamelistParser` are the input text
# itself (`_text`), and the current position in the text (`_pos`). The position
# is only changed via the `_advance` method, which also maintains line and
# column numbers for error-reporting purposes. The `_settings` attribute
# holds the final output, i.e. the variable name-value pairs.
#
# Parsing errors are signaled by one of two exceptions. The first is
# `_NamelistParseError`, which always signals an unrecoverable error. This is
# caught and translated to a user-visible error in `parse`. The second is
# `_NamelistEOF`, which may or may not represent a true error. During parsing of
# a standard namelist, it is treated in the same manner as
# `_NamelistParseError`, unless it occurs outside of any namelist group, in
# which case the `parse_namelist` method will catch it and return normally.
#
# The non-standard "groupless" format complicates things significantly by
# allowing an end-of-file at any location where a '/' would normally be. This is
# the reason for most of the `allow_eof` flags and related logic, since any
# `_NamelistEOF` exceptions raised must be caught and dealt with.
#
###############################################################################

# Disable these because of doctest, and because we don't typically follow the
# (rather specific) pylint naming conventions.
# pylint: disable=line-too-long,too-many-lines,invalid-name

import re

# Disable these because this is our standard setup
# pylint: disable=wildcard-import,unused-wildcard-import

from CIME.XML.standard_module_setup import *
from CIME.utils import expect

logger = logging.getLogger(__name__)

# Fortran syntax regular expressions.
# Variable names.
FORTRAN_NAME_REGEX = re.compile(r"^[a-z][a-z0-9_]{0,62}$", re.IGNORECASE)
FORTRAN_LITERAL_REGEXES = {}
# Integer literals.
_int_re_string = r"(\+|-)?[0-9]+"
FORTRAN_LITERAL_REGEXES['integer'] = re.compile("^" + _int_re_string + "$")
# Real/complex literals.
_ieee_exceptional_re_string = r"inf(inity)?|nan(\([^)]+\))?"
_float_re_string = r"((\+|-)?([0-9]+(\.[0-9]*)?|\.[0-9]+)([ed]?%s)?|%s)" % \
                   (_int_re_string, _ieee_exceptional_re_string)
FORTRAN_LITERAL_REGEXES['real'] = re.compile("^" + _float_re_string + "$",
                                             re.IGNORECASE)
FORTRAN_LITERAL_REGEXES['complex'] = re.compile(r"^\([ \n]*" +
                                                _float_re_string +
                                                r"[ \n]*,[ \n]*" +
                                                _float_re_string +
                                                r"[ \n]*\)$", re.IGNORECASE)
# Character literals.
_char_single_re_string = r"'[^']*(''[^']*)*'"
_char_double_re_string = r'"[^"]*(""[^"]*)*"'
FORTRAN_LITERAL_REGEXES['character'] = re.compile("^(" +
                                                  _char_single_re_string + "|" +
                                                  _char_double_re_string +
                                                  ")$")
# Logical literals.
FORTRAN_LITERAL_REGEXES['logical'] = re.compile(r"^\.?[tf][^=/ \n]*$",
                                                re.IGNORECASE)
# Repeated value prefix.
FORTRAN_REPEAT_PREFIX_REGEX = re.compile(r"^[0-9]*[1-9]+[0-9]*\*")


def is_valid_fortran_name(string):
    """Check that a variable name is allowed in Fortran.

    The rules are:
    1. The name must start with a letter.
    2. All characters in a name must be alphanumeric (or underscores).
    3. The maximum name length is 63 characters.

    >>> is_valid_fortran_name("")
    False
    >>> is_valid_fortran_name("a")
    True
    >>> is_valid_fortran_name("A")
    True
    >>> is_valid_fortran_name("2")
    False
    >>> is_valid_fortran_name("_")
    False
    >>> is_valid_fortran_name("abc#123")
    False
    >>> is_valid_fortran_name("aLiBi_123")
    True
    >>> is_valid_fortran_name("A" * 64)
    False
    >>> is_valid_fortran_name("A" * 63)
    True
    """
    return FORTRAN_NAME_REGEX.search(string) is not None


def fortran_namelist_base_value(string):
    r"""Strip off whitespace and repetition syntax from a namelist value.

    >>> fortran_namelist_base_value("")
    ''
    >>> fortran_namelist_base_value("f")
    'f'
    >>> fortran_namelist_base_value("6*")
    ''
    >>> fortran_namelist_base_value("6*f")
    'f'
    >>> fortran_namelist_base_value(" \n6* \n")
    ''
    >>> fortran_namelist_base_value("\n 6*f\n ")
    'f'
    """
    # Strip leading/trailing whitespace.
    string = string.strip(" \n")
    # Strip off repeated value prefix.
    if FORTRAN_REPEAT_PREFIX_REGEX.search(string) is not None:
        string = string[string.find('*') + 1:]
    return string


def character_literal_to_string(literal):
    """Convert a Fortran character literal to a Python string.

    This function assumes (without checking) that `literal` is a valid literal.

    >>> character_literal_to_string("'blah'")
    'blah'
    >>> character_literal_to_string('"blah"')
    'blah'
    >>> character_literal_to_string("'don''t'")
    "don't"
    >>> character_literal_to_string('"' + '""Hello!""' + '"')
    '"Hello!"'
    """
    # Figure out whether a quote or apostrophe is the delimiter.
    delimiter = None
    for char in literal:
        if char in ("'", '"'):
            delimiter = char
    # Find left and right edges of the string, extract middle.
    left_pos = literal.find(delimiter)
    right_pos = literal.rfind(delimiter)
    new_literal = literal[left_pos+1:right_pos]
    # Replace escaped quote and apostrophe characters.
    return new_literal.replace(delimiter * 2, delimiter)


def string_to_character_literal(string):
    r"""Convert a Python string to a Fortran character literal.

    This function always uses double quotes (") as the delimiter.

    >>> string_to_character_literal('blah')
    '"blah"'
    >>> string_to_character_literal("'blah'")
    '"\'blah\'"'
    >>> string_to_character_literal('She said "Hi!".')
    '"She said ""Hi!""."'
    """
    string = string.replace('"', '""')
    return '"' + string + '"'

def is_valid_fortran_namelist_literal(type_, string):
    r"""Determine whether a literal is valid in a Fortran namelist.

    Note that kind parameters are *not* allowed in namelists, which simplifies
    this check a bit. Internal whitespace is allowed for complex and character
    literals only. BOZ literals and compiler extensions (e.g. backslash escapes)
    are not allowed.

    Null values, however, are allowed for all types. This means that passing in
    a string containing nothing but spaces and newlines will always cause
    `True` to be returned. Repetition (e.g. `5*'a'`) is also allowed, including
    repetition of null values.

    Detailed rules and examples follow.

    Integers: Must be a sequence of one or more digits, with an optional sign.

    >>> is_valid_fortran_namelist_literal("integer", "")
    True
    >>> is_valid_fortran_namelist_literal("integer", " ")
    True
    >>> is_valid_fortran_namelist_literal("integer", "\n")
    True
    >>> is_valid_fortran_namelist_literal("integer", "5*")
    True
    >>> is_valid_fortran_namelist_literal("integer", "1")
    True
    >>> is_valid_fortran_namelist_literal("integer", "5*1")
    True
    >>> is_valid_fortran_namelist_literal("integer", " 5*1")
    True
    >>> is_valid_fortran_namelist_literal("integer", "5* 1")
    False
    >>> is_valid_fortran_namelist_literal("integer", "5 *1")
    False
    >>> is_valid_fortran_namelist_literal("integer", "a")
    False
    >>> is_valid_fortran_namelist_literal("integer", " 1")
    True
    >>> is_valid_fortran_namelist_literal("integer", "1 ")
    True
    >>> is_valid_fortran_namelist_literal("integer", "1 2")
    False
    >>> is_valid_fortran_namelist_literal("integer", "0123456789")
    True
    >>> is_valid_fortran_namelist_literal("integer", "+22")
    True
    >>> is_valid_fortran_namelist_literal("integer", "-26")
    True
    >>> is_valid_fortran_namelist_literal("integer", "2A")
    False
    >>> is_valid_fortran_namelist_literal("integer", "1_8")
    False
    >>> is_valid_fortran_namelist_literal("integer", "2.1")
    False
    >>> is_valid_fortran_namelist_literal("integer", "2e6")
    False

    Reals:
    - For fixed-point format, there is an optional sign, followed by an integer
    part, or a decimal point followed by a fractional part, or both.
    - Scientific notation is allowed, with an optional, case-insensitive "e" or
    "d" followed by an optionally-signed integer exponent. (Either the "e"/"d"
    or a sign must be present to separate the number from the exponent.)
    - The (case-insensitive) strings "inf", "infinity", and "nan" are allowed.
    NaN values can also contain additional information in parentheses, e.g.
    "NaN(x1234ABCD)".

    >>> is_valid_fortran_namelist_literal("real", "")
    True
    >>> is_valid_fortran_namelist_literal("real", "a")
    False
    >>> is_valid_fortran_namelist_literal("real", "1")
    True
    >>> is_valid_fortran_namelist_literal("real", " 1")
    True
    >>> is_valid_fortran_namelist_literal("real", "1 ")
    True
    >>> is_valid_fortran_namelist_literal("real", "1 2")
    False
    >>> is_valid_fortran_namelist_literal("real", "+1")
    True
    >>> is_valid_fortran_namelist_literal("real", "-1")
    True
    >>> is_valid_fortran_namelist_literal("real", "1.")
    True
    >>> is_valid_fortran_namelist_literal("real", "1.5")
    True
    >>> is_valid_fortran_namelist_literal("real", ".5")
    True
    >>> is_valid_fortran_namelist_literal("real", "+.5")
    True
    >>> is_valid_fortran_namelist_literal("real", ".")
    False
    >>> is_valid_fortran_namelist_literal("real", "+")
    False
    >>> is_valid_fortran_namelist_literal("real", "1e6")
    True
    >>> is_valid_fortran_namelist_literal("real", "1e-6")
    True
    >>> is_valid_fortran_namelist_literal("real", "1e+6")
    True
    >>> is_valid_fortran_namelist_literal("real", ".5e6")
    True
    >>> is_valid_fortran_namelist_literal("real", "1e")
    False
    >>> is_valid_fortran_namelist_literal("real", "1D6")
    True
    >>> is_valid_fortran_namelist_literal("real", "1q6")
    False
    >>> is_valid_fortran_namelist_literal("real", "1+6")
    True
    >>> is_valid_fortran_namelist_literal("real", "1.6.5")
    False
    >>> is_valid_fortran_namelist_literal("real", "1._8")
    False
    >>> is_valid_fortran_namelist_literal("real", "1,5")
    False
    >>> is_valid_fortran_namelist_literal("real", "inf")
    True
    >>> is_valid_fortran_namelist_literal("real", "INFINITY")
    True
    >>> is_valid_fortran_namelist_literal("real", "NaN")
    True
    >>> is_valid_fortran_namelist_literal("real", "nan(x56)")
    True
    >>> is_valid_fortran_namelist_literal("real", "nan())")
    False

    Complex numbers:
    - A pair of real numbers enclosed by parentheses, and separated by a comma.
    - Any number of spaces or newlines may be placed before or after each real.

    >>> is_valid_fortran_namelist_literal("complex", "")
    True
    >>> is_valid_fortran_namelist_literal("complex", "()")
    False
    >>> is_valid_fortran_namelist_literal("complex", "(,)")
    False
    >>> is_valid_fortran_namelist_literal("complex", "( ,\n)")
    False
    >>> is_valid_fortran_namelist_literal("complex", "(a,2.)")
    False
    >>> is_valid_fortran_namelist_literal("complex", "(1.,b)")
    False
    >>> is_valid_fortran_namelist_literal("complex", "(1,2)")
    True
    >>> is_valid_fortran_namelist_literal("complex", "(-1.e+06,+2.d-5)")
    True
    >>> is_valid_fortran_namelist_literal("complex", "(inf,NaN)")
    True
    >>> is_valid_fortran_namelist_literal("complex", "(  1. ,  2. )")
    True
    >>> is_valid_fortran_namelist_literal("complex", "( \n \n 1. \n,\n 2.\n)")
    True
    >>> is_valid_fortran_namelist_literal("complex", " (1.,2.)")
    True
    >>> is_valid_fortran_namelist_literal("complex", "(1.,2.) ")
    True

    Character sequences (strings):
    - Must begin and end with the same delimiter character, either a single
    quote (apostrophe), or a double quote (quotation mark).
    - Whichever character is used as a delimiter must not appear in the
    string itself, unless it appears in doubled pairs (e.g. '''' or "'" are the
    two ways of representing a string containing a single apostrophe).
    - Note that newlines cannot be represented in a namelist character literal
    since they are interpreted as an "end of record", but they are allowed as
    long as they don't come between one of the aforementioned double pairs of
    characters.

    >>> is_valid_fortran_namelist_literal("character", "")
    True
    >>> is_valid_fortran_namelist_literal("character", "''")
    True
    >>> is_valid_fortran_namelist_literal("character", " ''")
    True
    >>> is_valid_fortran_namelist_literal("character", "'\n'")
    True
    >>> is_valid_fortran_namelist_literal("character", "''\n''")
    False
    >>> is_valid_fortran_namelist_literal("character", "'''")
    False
    >>> is_valid_fortran_namelist_literal("character", "''''")
    True
    >>> is_valid_fortran_namelist_literal("character", "'''Cookie'''")
    True
    >>> is_valid_fortran_namelist_literal("character", "'''Cookie''")
    False
    >>> is_valid_fortran_namelist_literal("character", "'\"'")
    True
    >>> is_valid_fortran_namelist_literal("character", "'\"\"'")
    True
    >>> is_valid_fortran_namelist_literal("character", '""')
    True
    >>> is_valid_fortran_namelist_literal("character", '"" ')
    True
    >>> is_valid_fortran_namelist_literal("character", '"\n"')
    True
    >>> is_valid_fortran_namelist_literal("character", '""\n""')
    False
    >>> is_valid_fortran_namelist_literal("character", '""' + '"')
    False
    >>> is_valid_fortran_namelist_literal("character", '""' + '""')
    True
    >>> is_valid_fortran_namelist_literal("character", '"' + '""Cookie""' + '"')
    True
    >>> is_valid_fortran_namelist_literal("character", '""Cookie""' + '"')
    False
    >>> is_valid_fortran_namelist_literal("character", '"\'"')
    True
    >>> is_valid_fortran_namelist_literal("character", '"\'\'"')
    True

    Logicals:
    - Must contain a (case-insensitive) "t" or "f".
    - This must be either the first nonblank character, or the second following
    a period.
    - The rest of the string is ignored, but cannot contain a comma, newline,
    equals sign, slash, or space (except that trailing spaces are allowed and
    ignored).

    >>> is_valid_fortran_namelist_literal("logical", "")
    True
    >>> is_valid_fortran_namelist_literal("logical", "t")
    True
    >>> is_valid_fortran_namelist_literal("logical", "F")
    True
    >>> is_valid_fortran_namelist_literal("logical", ".T")
    True
    >>> is_valid_fortran_namelist_literal("logical", ".f")
    True
    >>> is_valid_fortran_namelist_literal("logical", " f")
    True
    >>> is_valid_fortran_namelist_literal("logical", " .t")
    True
    >>> is_valid_fortran_namelist_literal("logical", "at")
    False
    >>> is_valid_fortran_namelist_literal("logical", ".TRUE.")
    True
    >>> is_valid_fortran_namelist_literal("logical", ".false.")
    True
    >>> is_valid_fortran_namelist_literal("logical", ".TEXAS$")
    True
    >>> is_valid_fortran_namelist_literal("logical", ".f=")
    False
    >>> is_valid_fortran_namelist_literal("logical", ".f/1")
    False
    >>> is_valid_fortran_namelist_literal("logical", ".t\nted")
    False
    >>> is_valid_fortran_namelist_literal("logical", ".Fant astic")
    False
    >>> is_valid_fortran_namelist_literal("logical", ".t2 ")
    True
    """
    expect(type_ in FORTRAN_LITERAL_REGEXES,
           "Invalid Fortran type for a namelist: %r" % str(type_))
    # Strip off whitespace and repetition.
    string = fortran_namelist_base_value(string)
    # Null values are always allowed.
    if string == '':
        return True
    return FORTRAN_LITERAL_REGEXES[type_].search(string) is not None


def literal_to_python_value(literal, type_=None):
    r"""Convert a Fortran literal string to a Python value.

    This function assumes that the input contains a single value, i.e.
    repetition syntax is not used. The type can be specified by passing a string
    as the `type_` argument, or if this option is not provided, this function
    will attempt to autodetect the variable type.

    Note that it is not possible to be certain whether a literal like "123" is
    intended to represent an integer or a floating-point value, however, nor can
    we be certain of the precision that will be used to hold this value in
    actual Fortran code. We also cannot use the optional information in a NaN
    float, so this will cause the function to throw an error if that information
    is present (e.g. a string like "NAN(1234)" will cause an error).

    The Python type of the return value is as follows for different `type_`
    arguments:
    "character" - `str`
    "complex" - `complex`
    "integer" - `int`
    "logical" - `bool`
    "real" - `float`

    If a null value is input (i.e. the empty string), `None` will be returned.

    >>> literal_to_python_value("'She''s a winner!'")
    "She's a winner!"
    >>> literal_to_python_value("1")
    1
    >>> literal_to_python_value("1.")
    1.0
    >>> literal_to_python_value(" (\n 1. , 2. )\n ")
    (1+2j)
    >>> literal_to_python_value(".true.")
    True
    >>> literal_to_python_value("Fortune")
    False
    >>> literal_to_python_value("bacon")
    Traceback (most recent call last):
    ...
    SystemExit: ERROR: 'bacon' is not a valid literal for any Fortran type.
    >>> literal_to_python_value("1", type_="real")
    1.0
    >>> literal_to_python_value("bacon", type_="logical")
    Traceback (most recent call last):
    ...
    SystemExit: ERROR: 'bacon' is not a valid literal of type 'logical'.
    >>> literal_to_python_value("1", type_="booga")
    Traceback (most recent call last):
    ...
    SystemExit: ERROR: Invalid Fortran type for a namelist: 'booga'
    >>> literal_to_python_value("2*1")
    Traceback (most recent call last):
    ...
    SystemExit: ERROR: Cannot use repetition syntax in literal_to_python_value
    >>> literal_to_python_value("")
    >>> literal_to_python_value("-1.D+10")
    -10000000000.0
    >>> literal_to_python_value("nan(1234)")
    Traceback (most recent call last):
    ...
    ValueError: invalid literal for float(): nan(1234)
    """
    expect(FORTRAN_REPEAT_PREFIX_REGEX.search(literal) is None,
           "Cannot use repetition syntax in literal_to_python_value")
    # Handle null value.
    if fortran_namelist_base_value(literal) == '':
        return None
    if type_ is None:
        # Autodetect type.
        for test_type in ('character', 'complex', 'integer', 'logical', 'real'):
            if is_valid_fortran_namelist_literal(test_type, literal):
                type_ = test_type
                break
        expect(type_ is not None,
               "%r is not a valid literal for any Fortran type." % str(literal))
    else:
        # Check that type is valid.
        expect(is_valid_fortran_namelist_literal(type_, literal),
               "%r is not a valid literal of type %r." %
               (str(literal), str(type_)))
    # Conversion for each type.
    if type_ == 'character':
        return character_literal_to_string(literal)
    elif type_ == 'complex':
        literal = literal.lstrip(' \n(').rstrip(' \n)')
        real_part, _, imag_part = literal.partition(',')
        return complex(float(real_part), float(imag_part))
    elif type_ == 'integer':
        return int(literal)
    elif type_ == 'logical':
        literal = literal.lstrip(' \n.')
        return literal[0] in 'tT'
    elif type_ == 'real':
        literal = literal.lower().replace('d', 'e')
        return float(literal)


def expand_literal_list(literals):
    """Expands a list of literal values to get rid of repetition syntax.

    >>> expand_literal_list([])
    []
    >>> expand_literal_list(['true'])
    ['true']
    >>> expand_literal_list(['1', '2', 'f*', '3*3', '5'])
    ['1', '2', 'f*', '3', '3', '3', '5']
    >>> expand_literal_list([u'2*f*'])
    [u'f*', u'f*']
    """
    expanded = []
    for literal in literals:
        if FORTRAN_REPEAT_PREFIX_REGEX.search(literal) is not None:
            num, _, value = literal.partition('*')
            expanded += int(num) * [value]
        else:
            expanded.append(literal)
    return expanded


def compress_literal_list(literals):
    """Uses repetition syntax to shorten a literal list.

    >>> compress_literal_list([])
    []
    >>> compress_literal_list(['true'])
    ['true']
    >>> compress_literal_list(['1', '2', 'f*', '3', '3', '3', '5'])
    ['1', '2', 'f*', '3*3', '5']
    >>> compress_literal_list([u'f*', u'f*'])
    [u'2*f*']
    """
    compressed = []
    if len(literals) == 0:
        return compressed
    # Start with the first literal.
    old_literal = literals[0]
    num_reps = 1
    for literal in literals[1:]:
        if literal == old_literal:
            # For each new literal, if it matches the old one, it increases the
            # number of repetitions by one.
            num_reps += 1
        else:
            # Otherwise, write out the previous literal and start tracking the
            # new one.
            rep_str = str(num_reps) + '*' if num_reps > 1 else ''
            if isinstance(old_literal, basestring):
                compressed.append(rep_str + old_literal)
            else:
                compressed.append(rep_str + str(old_literal))
            old_literal = literal
            num_reps = 1
    rep_str = str(num_reps) + '*' if num_reps > 1 else ''
    if isinstance(old_literal, basestring):
        compressed.append(rep_str + old_literal)
    else:
        compressed.append(rep_str + str(old_literal))
    return compressed 

def merge_literal_lists(default, overwrite):
    """Merge two lists of literal value strings.

    The `overwrite` values have higher precedence, so will overwrite the
    `default` values. However, if `overwrite` contains null values, or is
    shorter than `default` (and thus implicitly ends in null values), the
    elements of `default` will be used where `overwrite` is null.

    >>> merge_literal_lists([], [])
    []
    >>> merge_literal_lists(['true'], ['false'])
    ['false']
    >>> merge_literal_lists([], ['false'])
    ['false']
    >>> merge_literal_lists(['true'], [''])
    ['true']
    >>> merge_literal_lists([], [''])
    ['']
    >>> merge_literal_lists(['true'], [])
    ['true']
    >>> merge_literal_lists(['true'], [])
    ['true']
    >>> merge_literal_lists(['3*false', '3*true'], ['true', '4*', 'false'])
    ['true', '2*false', '2*true', 'false']
    """
    merged = []
    default = expand_literal_list(default)
    overwrite = expand_literal_list(overwrite)
    for default_elem, elem in zip(default, overwrite):
        if elem == '':
            merged.append(default_elem)
        else:
            merged.append(elem)
    def_len = len(default)
    ovw_len = len(overwrite)
    if ovw_len < def_len:
        merged[ovw_len:def_len] = default[ovw_len:def_len]
    else:
        merged[def_len:ovw_len] = overwrite[def_len:ovw_len]
    return compress_literal_list(merged)


def parse(in_file=None, text=None, groupless=False, convert_tab_to_space=True):
    """Parse a Fortran namelist.

    The `in_file` argument must be either a `str` or `unicode` object containing
    a file name, or a text I/O object with a `read` method that returns the text
    of the namelist.

    Alternatively, the `text` argument can be provided, in which case it must be
    the text of the namelist itself.

    The `groupless` argument changes namelist parsing in two ways:

    1. `parse` allows an alternate file format where no group names or slashes
       are present. In effect, the file is parsed as if an invisible, arbitrary
       group name was prepended, and an invisible slash was appended. However,
       if any group names actually are present, the file is parsed normally.
    2. The return value of this function is not a `Namelist` object. Instead a
       single, flattened dictionary of name-value pairs is returned.

    The `convert_tab_to_space` option can be used to force all tabs in the file
    to be converted to spaces, and is on by default. Note that this will usually
    allow files that use tabs as whitespace to be parsed. However, the
    implementation of this option is crude; it converts *all* tabs in the file,
    including those in character literals. (Note that there are many characters
    that cannot be passed in via namelist in any standard way, including '\n',
    so it is already a bad idea to assume that the namelist will preserve
    whitespace in strings, aside from simple spaces.)

    The return value, if `groupless=False`, is a `Namelist` object.

    All names and values returned are ultimately unicode strings. E.g. a value
    of "6*2" is returned as that string; it is not converted to 6 copies of the
    Python integer `2`. Null values are returned as the empty string ("").
    """
    expect(in_file is not None or text is not None,
           "Must specify an input file or text to the namelist parser.")
    expect(in_file is None or text is None,
           "Cannot specify both input file and text to the namelist parser.")
    if isinstance(in_file, str) or isinstance(in_file, unicode):
        logger.debug("Reading namelist at: %s", in_file)
        with open(in_file) as in_file_obj:
            text = in_file_obj.read()
    elif in_file is not None:
        logger.debug("Reading namelist from file object")
        text = in_file.read()
    if convert_tab_to_space:
        text = text.replace('\t', ' ')
    try:
        namelist_dict = _NamelistParser(text, groupless).parse_namelist()
    except (_NamelistEOF, _NamelistParseError) as error:
        # Deal with unexpected EOF or other parsing errors.
        expect(False, str(error))
    if groupless:
        return namelist_dict
    else:
        return Namelist(namelist_dict)


class Namelist(object):

    """Class representing a Fortran namelist.

    Public methods:
    __init__
    delete_variable
    get_group_names
    get_value
    get_variable_names
    get_variable_value
    merge_nl
    set_variable_value
    write
    """

    def __init__(self, groups=None):
        """Construct a new `Namelist` object.

        The `groups` argument is a dictionary associating group names to
        dictionaries of name/value pairs. If omitted, an empty namelist object
        is created.

        Unless you are deliberately creating an empty `Namelist`, it is easier/
        safer to use `parse` than to directly call this constructor.
        """
        self._groups = {}
        if groups is not None:
            for group_name in groups:
                group_lc = group_name.lower()
                self._groups[group_lc] = {}
                for variable_name in groups[group_name]:
                    variable_lc = variable_name.lower()
                    self._groups[group_lc][variable_lc] = \
                                        groups[group_name][variable_name]

    def get_group_names(self):
        """Return a list of all groups in the namelist.

        >>> Namelist().get_group_names()
        []
        >>> sorted(parse(text='&foo / &bar /').get_group_names())
        [u'bar', u'foo']
        """
        return self._groups.keys()

    def get_variable_names(self, group_name):
        """Return a list of all variables in the given namelist group.

        If the specified group is not in the namelist, returns an empty list.

        >>> Namelist().get_variable_names('foo')
        []
        >>> x = parse(text='&foo bar=,bazz=true,bang=6*""/')
        >>> sorted(x.get_variable_names('fOo'))
        [u'bang', u'bar', u'bazz']
        """
        group_name = group_name.lower()
        if group_name not in self._groups:
            return []
        return self._groups[group_name].keys()

    def get_variable_value(self, group_name, variable_name):
        """Return the value of the specified variable.

        This function always returns a non-empty list containing strings. If the
        specified `group_name` or `variable_name` is not present, `[u'']` is
        returned.

        >>> Namelist().get_variable_value('foo', 'bar')
        [u'']
        >>> parse(text='&foo bar=1,2 /').get_variable_value('foo', 'bazz')
        [u'']
        >>> parse(text='&foo bar=1,2 /').get_variable_value('foO', 'Bar')
        [u'1', u'2']
        """
        group_name = group_name.lower()
        variable_name = variable_name.lower()
        if group_name not in self._groups or \
           variable_name not in self._groups[group_name]:
            return [u'']
        return self._groups[group_name][variable_name]

    def get_value(self, variable_name):
        """Return the value of a uniquely-named variable.

        This function is similar to `get_variable_value`, except that it does
        not require a `group_name`, and it requires that the `variable_name` be
        unique across all groups.

        >>> parse(text='&foo bar=1 / &bazz bar=1 /').get_value('bar')
        Traceback (most recent call last):
        ...
        SystemExit: ERROR: Namelist.get_value: Variable %s is present in multiple groups: [u'bazz', u'foo']
        >>> parse(text='&foo bar=1 / &bazz /').get_value('Bar')
        [u'1']
        >>> parse(text='&foo / &bazz /').get_value('bar')
        [u'']
        """
        variable_name = variable_name.lower()
        possible_groups = [group_name for group_name in self._groups
                           if variable_name in self._groups[group_name]]
        expect(len(possible_groups) <= 1,
               "Namelist.get_value: Variable %s is present in multiple groups: "
               + str(possible_groups))
        if possible_groups:
            return self._groups[possible_groups[0]][variable_name]
        else:
            return [u'']

    def set_variable_value(self, group_name, variable_name, value):
        """Set the value of the specified variable.

        >>> x = parse(text='&foo bar=1 /')
        >>> x.set_variable_value('foo', 'bar', [u'2'])
        >>> x.set_variable_value('foo', 'bazz', [u'3'])
        >>> x.set_variable_value('Brack', 'baR', [u'4'])
        >>> x.get_variable_value('foo', 'bar')
        [u'2']
        >>> x.get_variable_value('foo', 'bazz')
        [u'3']
        >>> x.get_variable_value('brack', 'bar')
        [u'4']
        """
        group_name = group_name.lower()
        variable_name = variable_name.lower()
        if group_name not in self._groups:
            self._groups[group_name] = {}
        self._groups[group_name][variable_name] = value

    def delete_variable(self, group_name, variable_name):
        """Delete a variable from a specified group.

        If the specified group or variable does not exist, this is a no-op.

        >>> x = parse(text='&foo bar=1 /')
        >>> x.delete_variable('FOO', 'BAR')
        >>> x.delete_variable('foo', 'bazz')
        >>> x.delete_variable('brack', 'bazz')
        >>> x.get_variable_names('foo')
        []
        >>> x.get_variable_names('brack')
        []
        """
        group_name = group_name.lower()
        variable_name = variable_name.lower()
        if group_name in self._groups and \
           variable_name in self._groups[group_name]:
            del self._groups[group_name][variable_name]

    def merge_nl(self, other, overwrite=False):
        """Merge this namelist object with another.

        Values in the invoking (`self`) `Namelist` will take precedence over
        values in the `other` `Namelist`, unless `overwrite=True` is passed in,
        in which case `other` values take precedence.

        >>> x = parse(text='&foo bar=1 bazz=,2 brat=3/')
        >>> y = parse(text='&foo bar=2 bazz=3*1 baker=4 / &foo2 barter=5 /')
        >>> x.merge_nl(y)
        >>> sorted(x.get_group_names())
        [u'foo', u'foo2']
        >>> sorted(x.get_variable_names('foo'))
        [u'baker', u'bar', u'bazz', u'brat']
        >>> sorted(x.get_variable_names('foo2'))
        [u'barter']
        >>> x.get_value('bar')
        [u'1']
        >>> x.get_value('bazz')
        [u'1', u'2', u'1']
        >>> x.get_value('brat')
        [u'3']
        >>> x.get_value('baker')
        [u'4']
        >>> x.get_value('barter')
        [u'5']

        >>> x = parse(text='&foo bar=1 bazz=,2 brat=3/')
        >>> y = parse(text='&foo bar=2 bazz=3*1 baker=4 / &foo2 barter=5 /')
        >>> x.merge_nl(y, overwrite=True)
        >>> sorted(x.get_group_names())
        [u'foo', u'foo2']
        >>> sorted(x.get_variable_names('foo'))
        [u'baker', u'bar', u'bazz', u'brat']
        >>> sorted(x.get_variable_names('foo2'))
        [u'barter']
        >>> x.get_value('bar')
        [u'2']
        >>> x.get_value('bazz')
        [u'3*1']
        >>> x.get_value('brat')
        [u'3']
        >>> x.get_value('baker')
        [u'4']
        >>> x.get_value('barter')
        [u'5']
        """
        # Pretty simple strategy: go through the entire other namelist, and
        # merge all values with this one's.
        for group_name in other.get_group_names():
            for variable_name in other.get_variable_names(group_name):
                self_val = self.get_variable_value(group_name, variable_name)
                other_val = other.get_variable_value(group_name, variable_name)
                if overwrite:
                    merged_val = merge_literal_lists(self_val, other_val)
                else:
                    merged_val = merge_literal_lists(other_val, self_val)
                self.set_variable_value(group_name, variable_name, merged_val)

    def write(self, out_file, groups=None, append=False, format_='nml'):
        """Write a Fortran namelist to a file.

        As with `parse`, the `out_file` argument can be either a file name, or a
        file object with a `write` method that accepts unicode. If specified,
        the `groups` argument specifies a subset of all groups to write out.

        If `out_file` is a file name, and `append=True` is passed in, the
        namelist will be appended to the named file instead of overwriting it.
        The `append` option has no effect if a file object is passed in.

        The `format_` option can be either 'nml' (namelist) or 'rc', and
        specifies the file format. Formats other than 'nml' may not support all
        possible output values.
        """
        expect(format_ in ('nml', 'rc'),
               "Namelist.write: unexpected output format %r" % str(format_))
        if isinstance(out_file, str) or isinstance(out_file, unicode):
            logger.debug("Writing namelist to: %s", out_file)
            flag = 'a' if append else 'w'
            with open(out_file, flag) as file_obj:
                self._write(file_obj, groups, format_)
        else:
            logger.debug("Writing namelist to file object")
            self._write(out_file, groups, format_)

    def _write(self, out_file, groups, format_):
        """Unwrapped version of `write` assuming that a file object is input."""
        if groups is None:
            groups = self._groups.keys()
        if format_ == 'nml':
            equals = ' ='
        elif format_ == 'rc':
            equals = ':'
        for group_name in sorted(group.lower() for group in groups):
            if format_ == 'nml':
                out_file.write("&%s\n" % group_name)
            group = self._groups[group_name]
            for name in sorted(group.keys()):
                values = group[name]
                # To prettify things for long lists of values, build strings
                # line-by-line.
                if values[0] == "True" or values[0] == "False":
                    values[0] = values[0].replace("True",".true.").replace("False",".false.")
                lines = ["  %s%s %s" % (name, equals, values[0])]
                for value in values[1:]:
                    if value == "True" or value == "False":
                        value = value.replace("True",".true.").replace("False",".false.")
                    if len(lines[-1]) + len(value) <= 77:
                        lines[-1] += ", " + value
                    else:
                        lines[-1] += ",\n"
                        lines.append("      " + value)
                lines[-1] += "\n"
                for line in lines:
                    out_file.write(line)
            if format_ == 'nml':
                out_file.write("/\n")


class _NamelistEOF(Exception):

    """Exception thrown for an unexpected end-of-file in a namelist.

    This is an internal helper class, and should never be raised in a context
    where it would be visible to a user. (Typically it should be caught and
    converted to some other error, or ignored.)
    """

    def __init__(self, message=None):
        """Create a `_NamelistEOF`, optionally using an error message."""
        super(_NamelistEOF, self).__init__()
        self._message = message

    def __str__(self):
        """Get an error message suitable for display."""
        string = "Unexpected end of file encountered in namelist."
        if self._message is not None:
            string += " (%s)" % self._message
        return string


class _NamelistParseError(Exception):

    """Exception thrown when namelist input has a syntax error.

    This is an internal helper class, and should never be raised in a context
    where it would be visible to a user. (Typically it should be caught and
    converted to some other error, or ignored.)
    """

    def __init__(self, message=None):
        """Create a `_NamelistParseError`, optionally using an error message."""
        super(_NamelistParseError, self).__init__()
        self._message = message

    def __str__(self):
        """Get an error message suitable for display."""
        string = "Error in parsing namelist"
        if self._message is not None:
            string += ": %s" % self._message
        return string


class _NamelistParser(object): # pylint:disable=too-few-public-methods

    """Class to validate and read from Fortran namelist input.

    This is intended to be an internal helper class and should not be used
    directly. Use the `parse` function in this module instead.
    """

    def __init__(self, text, groupless=False):
        """Create a `_NamelistParser` given text to parse in a string."""
        # Current location within the file.
        self._pos = 0
        self._line = 1
        self._col = 0
        # Text and its size.
        self._text = unicode(text)
        self._len = len(self._text)
        # Dictionary with group names as keys, and dictionaries of variable
        # name-value pairs as values. (Or a single flat dictionary if
        # `groupless=True`.)
        self._settings = {}
        self._groupless = groupless

    def _line_col_string(self):
        r"""Return a string specifying the current line and column number.

        >>> x = _NamelistParser('abc\nd\nef')
        >>> x._advance(5)
        >>> x._line_col_string()
        'line 2, column 1'
        """
        return "line %s, column %s" % (self._line, self._col)

    def _curr(self):
        """Return the character at the current position."""
        return self._text[self._pos]

    def _next(self):
        """Return the character at the next position.

        >>> _NamelistParser(' ')._next()
        Traceback (most recent call last):
            ...
        _NamelistEOF: Unexpected end of file encountered in namelist. (At line 1, column 1)
        """
        # If at the end of the file, we should raise _NamelistEOF. The easiest
        # way to do this is to just advance.
        if self._pos == self._len - 1:
            self._advance()
        return self._text[self._pos+1]

    def _advance(self, nchars=1, check_eof=False):
        r"""Advance the parser's current position by `nchars` characters.

        The `nchars` argument must be non-negative. If the end of file is
        reached, an exception is thrown, unless `check_eof=True` is passed. If
        `check_eof=True` is passed, the position is advanced past the end of the
        file (`self._pos == `self._len`), and a boolean is returned to signal
        whether or not the end of the file was reached.

        >>> _NamelistParser('abcd')._advance(-1)
        Traceback (most recent call last):
            ...
        AssertionError: _NamelistParser attempted to 'advance' backwards
        >>> x = _NamelistParser('abc\nd\nef')
        >>> (x._pos, x._line, x._col)
        (0, 1, 0)
        >>> x._advance(0)
        >>> (x._pos, x._line, x._col)
        (0, 1, 0)
        >>> x._advance(2)
        >>> (x._pos, x._line, x._col)
        (2, 1, 2)
        >>> x._advance(1)
        >>> (x._pos, x._line, x._col)
        (3, 1, 3)
        >>> x._advance(1)
        >>> (x._pos, x._line, x._col)
        (4, 2, 0)
        >>> x._advance(3)
        >>> (x._pos, x._line, x._col)
        (7, 3, 1)
        >>> x._advance(1)
        Traceback (most recent call last):
            ...
        _NamelistEOF: Unexpected end of file encountered in namelist. (At line 3, column 2)
        >>> _NamelistParser('abc\n')._advance(4)
        Traceback (most recent call last):
            ...
        _NamelistEOF: Unexpected end of file encountered in namelist. (At line 2, column 0)
        >>> x = _NamelistParser('ab')
        >>> x._advance(check_eof=True)
        False
        >>> x._curr()
        u'b'
        >>> x._advance(check_eof=True)
        True
        """
        assert nchars >= 0, \
            "_NamelistParser attempted to 'advance' backwards"
        new_pos = min(self._pos + nchars, self._len)
        consumed_text = self._text[self._pos:new_pos]
        self._pos = new_pos
        lines = consumed_text.count('\n')
        self._line += lines
        # If we started a new line, set self._col to be relative to the start of
        # the current line.
        if lines > 0:
            self._col = -(consumed_text.rfind('\n') + 1)
        self._col += len(consumed_text)
        end_of_file = new_pos == self._len
        if check_eof:
            return end_of_file
        elif end_of_file:
            raise _NamelistEOF(message="At "+self._line_col_string())

    def _eat_whitespace(self, allow_initial_comment=False):
        r"""Advance until the next non-whitespace character.

        Returns a boolean representing whether anything was eaten. Note that
        this function also skips past new lines containing comments. Comments in
        the current line will be skipped if `allow_initial_comment=True` is
        passed in.

        >>> x = _NamelistParser(' \n a ')
        >>> x._eat_whitespace()
        True
        >>> x._curr()
        u'a'
        >>> x._eat_whitespace()
        False
        >>> x._advance()
        >>> x._eat_whitespace()
        Traceback (most recent call last):
            ...
        _NamelistEOF: Unexpected end of file encountered in namelist. (At line 2, column 3)
        >>> x = _NamelistParser(' \n! blah\n ! blah\n a')
        >>> x._eat_whitespace()
        True
        >>> x._curr()
        u'a'
        >>> x = _NamelistParser('! blah\n a')
        >>> x._eat_whitespace()
        False
        >>> x._curr()
        u'!'
        >>> x = _NamelistParser(' ! blah\n a')
        >>> x._eat_whitespace()
        True
        >>> x._curr()
        u'!'
        >>> x = _NamelistParser(' ! blah\n a')
        >>> x._eat_whitespace(allow_initial_comment=True)
        True
        >>> x._curr()
        u'a'
        """
        eaten = False
        comment_allowed = allow_initial_comment
        while True:
            while self._curr() in (' ', '\n'):
                comment_allowed |= self._curr() == '\n'
                eaten = True
                self._advance()
            # Note the reliance on short-circuit `and` here.
            if not (comment_allowed and self._eat_comment()):
                break
        return eaten

    def _eat_comment(self):
        r"""If currently positioned at a '!', advance past the comment's end.

        Only works properly if not currently inside a comment or string. Returns
        a boolean representing whether anything was eaten.

        >>> x = _NamelistParser('! foo\n ! bar\na ! bazz')
        >>> x._eat_comment()
        True
        >>> x._curr()
        u' '
        >>> x._eat_comment()
        False
        >>> x._eat_whitespace()
        True
        >>> x._eat_comment()
        True
        >>> x._curr()
        u'a'
        >>> x._advance(2)
        >>> x._eat_comment()
        Traceback (most recent call last):
            ...
        _NamelistEOF: Unexpected end of file encountered in namelist. (At line 3, column 8)
        >>> x = _NamelistParser('! foo\n')
        >>> x._eat_comment()
        Traceback (most recent call last):
            ...
        _NamelistEOF: Unexpected end of file encountered in namelist. (At line 2, column 0)
        """
        if self._curr() != '!':
            return False
        newline_pos = self._text[self._pos:].find('\n')
        if newline_pos == -1:
            # This is the last line.
            self._advance(self._len - self._pos)
        else:
            # Advance to the next line.
            self._advance(newline_pos)
            # Advance to the first character of the next line.
            self._advance()
        return True

    def _expect_char(self, chars):
        """Raise an error if the wrong character is present.

        Does not return anything, but raises a `_NamelistParseError` if `chars`
        does not contain the character at the current position.

        >>> x = _NamelistParser('ab')
        >>> x._expect_char('a')
        >>> x._advance()
        >>> x._expect_char('a')
        Traceback (most recent call last):
            ...
        _NamelistParseError: Error in parsing namelist: expected 'a' but found 'b' at line 1, column 1
        >>> x._expect_char('ab')
        """
        if self._curr() not in chars:
            if len(chars) == 1:
                char_description = repr(str(chars))
            else:
                char_description = "one of the characters in %r" % str(chars)
            raise _NamelistParseError("expected %s but found %r at %s" %
                                      (char_description, str(self._curr()),
                                       self._line_col_string()))

    def _parse_namelist_group_name(self):
        r"""Parses and returns a namelist group name at the current position.

        >>> _NamelistParser('abc')._parse_namelist_group_name()
        Traceback (most recent call last):
            ...
        _NamelistParseError: Error in parsing namelist: expected '&' but found 'a' at line 1, column 0
        >>> _NamelistParser('&abc')._parse_namelist_group_name()
        Traceback (most recent call last):
            ...
        _NamelistEOF: Unexpected end of file encountered in namelist. (At line 1, column 4)
        >>> _NamelistParser('&abc ')._parse_namelist_group_name()
        u'abc'
        >>> _NamelistParser('&abc\n')._parse_namelist_group_name()
        u'abc'
        >>> _NamelistParser('&abc/ ')._parse_namelist_group_name()
        Traceback (most recent call last):
            ...
        _NamelistParseError: Error in parsing namelist: 'abc/' is not a valid variable name at line 1, column 5
        >>> _NamelistParser('&abc= ')._parse_namelist_group_name()
        Traceback (most recent call last):
            ...
        _NamelistParseError: Error in parsing namelist: 'abc=' is not a valid variable name at line 1, column 5
        >>> _NamelistParser('& ')._parse_namelist_group_name()
        Traceback (most recent call last):
            ...
        _NamelistParseError: Error in parsing namelist: '' is not a valid variable name at line 1, column 1
        """
        self._expect_char("&")
        self._advance()
        return self._parse_variable_name(allow_equals=False)

    def _parse_variable_name(self, allow_equals=True):
        r"""Parses and returns a variable name at the current position.

        The `allow_equals` flag controls whether '=' can denote the end of the
        variable name; if it is `False`, only white space can be used for this
        purpose.

        >>> _NamelistParser('abc')._parse_variable_name()
        Traceback (most recent call last):
            ...
        _NamelistEOF: Unexpected end of file encountered in namelist. (At line 1, column 3)
        >>> _NamelistParser('abc ')._parse_variable_name()
        u'abc'
        >>> _NamelistParser('ABC ')._parse_variable_name()
        u'abc'
        >>> _NamelistParser('abc\n')._parse_variable_name()
        u'abc'
        >>> _NamelistParser('abc=')._parse_variable_name()
        u'abc'
        >>> _NamelistParser('abc, ')._parse_variable_name()
        Traceback (most recent call last):
            ...
        _NamelistParseError: Error in parsing namelist: 'abc,' is not a valid variable name at line 1, column 4
        >>> _NamelistParser(' ')._parse_variable_name()
        Traceback (most recent call last):
            ...
        _NamelistParseError: Error in parsing namelist: '' is not a valid variable name at line 1, column 0
        """
        old_pos = self._pos
        separators = (' ', '\n', '=') if allow_equals else (' ', '\n')
        while self._curr() not in separators:
            self._advance()
        text = self._text[old_pos:self._pos]
        if not is_valid_fortran_name(text):
            raise _NamelistParseError("%r is not a valid variable name at %s" %
                                      (str(text), self._line_col_string()))
        return text.lower()

    def _parse_character_literal(self):
        """Parse and return a character literal (a string).

        Position on return is the last character of the string; we avoid
        advancing past that in order to avoid potential EOF errors.

        >>> _NamelistParser('"abc')._parse_character_literal()
        Traceback (most recent call last):
            ...
        _NamelistEOF: Unexpected end of file encountered in namelist. (At line 1, column 4)
        >>> _NamelistParser('"abc" ')._parse_character_literal()
        u'"abc"'
        >>> _NamelistParser("'abc' ")._parse_character_literal()
        u"'abc'"
        >>> _NamelistParser("*abc* ")._parse_character_literal()
        Traceback (most recent call last):
            ...
        _NamelistParseError: Error in parsing namelist: *abc* is not a valid character literal at line 1, column 4
        >>> _NamelistParser("'abc''def' ")._parse_character_literal()
        u"'abc''def'"
        >>> _NamelistParser("'abc''' ")._parse_character_literal()
        u"'abc'''"
        >>> _NamelistParser("'''abc' ")._parse_character_literal()
        u"'''abc'"
        """
        delimiter = self._curr()
        old_pos = self._pos
        self._advance()
        while True:
            while self._curr() != delimiter:
                self._advance()
            # Avoid end-of-file condition.
            if self._pos == self._len - 1:
                break
            # Doubled delimiters are escaped.
            if self._next() == delimiter:
                self._advance(2)
            else:
                break
        text = self._text[old_pos:self._pos+1]
        if not is_valid_fortran_namelist_literal("character", text):
            raise _NamelistParseError("%s is not a valid character literal at %s" %
                                      (text, self._line_col_string()))
        return text

    def _parse_complex_literal(self):
        """Parse and return a complex literal.

        Position on return is the last character of the string; we avoid
        advancing past that in order to avoid potential EOF errors.

        >>> _NamelistParser('(1.,2.')._parse_complex_literal()
        Traceback (most recent call last):
            ...
        _NamelistEOF: Unexpected end of file encountered in namelist. (At line 1, column 6)
        >>> _NamelistParser('(1.,2.) ')._parse_complex_literal()
        u'(1.,2.)'
        >>> _NamelistParser("(A,B) ")._parse_complex_literal()
        Traceback (most recent call last):
            ...
        _NamelistParseError: Error in parsing namelist: '(A,B)' is not a valid complex literal at line 1, column 4
        """
        old_pos = self._pos
        while self._curr() != ')':
            self._advance()
        text = self._text[old_pos:self._pos+1]
        if not is_valid_fortran_namelist_literal("complex", text):
            raise _NamelistParseError("%r is not a valid complex literal at %s"
                                      % (str(text), self._line_col_string()))
        return text

    def _look_ahead_for_equals(self, pos):
        r"""Look ahead to see if the next whitespace character is '='.

        The `pos` argument is the position in the text to start from while
        looking. This function returns a boolean.

        >>> _NamelistParser('=')._look_ahead_for_equals(0)
        True
        >>> _NamelistParser('a \n=')._look_ahead_for_equals(1)
        True
        >>> _NamelistParser('')._look_ahead_for_equals(0)
        False
        >>> _NamelistParser('a=')._look_ahead_for_equals(0)
        False
        """
        for test_pos in range(pos, self._len):
            if self._text[test_pos] not in (' ', '\n'):
                if self._text[test_pos] == '=':
                    return True
                else:
                    break
        return False

    def _parse_literal(self, allow_name=False, allow_eof_end=False):
        r"""Parse and return a variable value at the current position.

        The basic strategy is this:
        - If a value starts with an apostrophe/quotation mark, parse it as a
        character value (string).
        - If a value starts with a left parenthesis, parse it as a complex
        number.
        - Otherwise, read until the next value separator (comma, space, newline,
        or slash).

        If the argument `allow_name=True` is passed in, we allow the possibility
        that the current position is at the start of the variable name in a new
        name-value pair. In this case, `None` is returned, and the current
        position remains unchanged.

        If the argument `allow_eof_end=True` is passed in, we allow end-of-file
        to mark the end of a literal.

        >>> _NamelistParser('"abc" ')._parse_literal()
        u'"abc"'
        >>> _NamelistParser("'abc' ")._parse_literal()
        u"'abc'"
        >>> _NamelistParser('"abc"')._parse_literal()
        Traceback (most recent call last):
            ...
        _NamelistEOF: Unexpected end of file encountered in namelist. (At line 1, column 5)
        >>> _NamelistParser('"abc"')._parse_literal(allow_eof_end=True)
        u'"abc"'
        >>> _NamelistParser('(1.,2.) ')._parse_literal()
        u'(1.,2.)'
        >>> _NamelistParser('(1.,2.)')._parse_literal()
        Traceback (most recent call last):
            ...
        _NamelistEOF: Unexpected end of file encountered in namelist. (At line 1, column 7)
        >>> _NamelistParser('(1.,2.)')._parse_literal(allow_eof_end=True)
        u'(1.,2.)'
        >>> _NamelistParser('5 ')._parse_literal()
        u'5'
        >>> _NamelistParser('6.9 ')._parse_literal()
        u'6.9'
        >>> _NamelistParser('inf ')._parse_literal()
        u'inf'
        >>> _NamelistParser('nan(booga) ')._parse_literal()
        u'nan(booga)'
        >>> _NamelistParser('.FLORIDA$ ')._parse_literal()
        u'.FLORIDA$'
        >>> _NamelistParser('hamburger ')._parse_literal()
        Traceback (most recent call last):
            ...
        _NamelistParseError: Error in parsing namelist: expected literal value, but got 'hamburger' at line 1, column 9
        >>> _NamelistParser('5,')._parse_literal()
        u'5'
        >>> _NamelistParser('5\n')._parse_literal()
        u'5'
        >>> _NamelistParser('5/')._parse_literal()
        u'5'
        >>> _NamelistParser(',')._parse_literal()
        u''
        >>> _NamelistParser('6*5 ')._parse_literal()
        u'6*5'
        >>> _NamelistParser('6*(1., 2.) ')._parse_literal()
        u'6*(1., 2.)'
        >>> _NamelistParser('6*"a" ')._parse_literal()
        u'6*"a"'
        >>> _NamelistParser('6*')._parse_literal()
        Traceback (most recent call last):
            ...
        _NamelistEOF: Unexpected end of file encountered in namelist. (At line 1, column 2)
        >>> _NamelistParser('6*')._parse_literal(allow_eof_end=True)
        u'6*'
        >>> _NamelistParser('foo= ')._parse_literal()
        Traceback (most recent call last):
            ...
        _NamelistParseError: Error in parsing namelist: expected literal value, but got 'foo=' at line 1, column 4
        >>> _NamelistParser('5,')._parse_literal(allow_name=True)
        u'5'
        >>> x = _NamelistParser('foo= ')
        >>> x._parse_literal(allow_name=True)
        >>> x._curr()
        u'f'
        >>> _NamelistParser('6*foo= ')._parse_literal(allow_name=True)
        Traceback (most recent call last):
            ...
        _NamelistParseError: Error in parsing namelist: expected literal value, but got '6*foo=' at line 1, column 6
        >>> x = _NamelistParser('foo = ')
        >>> x._parse_literal(allow_name=True)
        >>> x._curr()
        u'f'
        >>> x = _NamelistParser('foo\n= ')
        >>> x._parse_literal(allow_name=True)
        >>> x._curr()
        u'f'
        >>> _NamelistParser('')._parse_literal(allow_eof_end=True)
        u''
        """
        # Deal with empty input string.
        if allow_eof_end and self._pos == self._len:
            return u''
        # Deal with a repeated value prefix.
        old_pos = self._pos
        if FORTRAN_REPEAT_PREFIX_REGEX.search(self._text[self._pos:]):
            allow_name = False
            while self._curr() != '*':
                self._advance()
            if self._advance(check_eof=allow_eof_end):
                # In case the file ends with the 'r*' form of null value.
                return self._text[old_pos:]
        prefix = self._text[old_pos:self._pos]
        # Deal with delimited literals.
        if self._curr() in ('"', "'"):
            literal = self._parse_character_literal()
            self._advance(check_eof=allow_eof_end)
            return prefix + literal
        if self._curr() == '(':
            literal = self._parse_complex_literal()
            self._advance(check_eof=allow_eof_end)
            return prefix + literal
        # Deal with non-delimited literals.
        new_pos = self._pos
        separators = [' ', '\n', ',', '/']
        if allow_name:
            separators.append('=')
        while new_pos != self._len and self._text[new_pos] not in separators:
            new_pos += 1
        if not allow_eof_end and new_pos == self._len:
            # At the end of the file, give up by throwing an EOF.
            self._advance(self._len)
        # If `allow_name` is set, we need to check and see if the next non-blank
        # character is '=', and return `None` if so.
        if allow_name and self._look_ahead_for_equals(new_pos):
            return
        self._advance(new_pos - self._pos, check_eof=allow_eof_end)
        text = self._text[old_pos:self._pos]
        if not any(is_valid_fortran_namelist_literal(type_, text)
                   for type_ in ("integer", "logical", "real")):
            raise _NamelistParseError("expected literal value, but got %r at %s"
                                      % (str(text), self._line_col_string()))
        return text

    def _expect_separator(self, allow_eof=False):
        r"""Advance past the current value separator.

        This function raises an error if we are not positioned at a valid value
        separator. It returns `False` if the end-of-namelist ('/') was
        encountered, in which case this function will leave the current position
        at the '/'. This function returns `True` otherwise, and skips to the
        location of the next non-whitespace character.

        If `allow_eof=True` is passed to this function, the meanings of '/' and
        the end-of-file are reversed. That is, an exception will be raised if a
        '/' is encountered, but the end-of-file will cause `False` to be
        returned rather than `True`. (An end-of-file after a ',' will be taken
        to be part of the next separator, and will not cause `False` to be
        returned.)

        >>> x = _NamelistParser("\na")
        >>> x._expect_separator()
        True
        >>> x._curr()
        u'a'
        >>> x = _NamelistParser(" a")
        >>> x._expect_separator()
        True
        >>> x._curr()
        u'a'
        >>> x = _NamelistParser(",a")
        >>> x._expect_separator()
        True
        >>> x._curr()
        u'a'
        >>> x = _NamelistParser("/a")
        >>> x._expect_separator()
        False
        >>> x._curr()
        u'/'
        >>> x = _NamelistParser("a")
        >>> x._expect_separator()
        Traceback (most recent call last):
            ...
        _NamelistParseError: Error in parsing namelist: expected one of the characters in ' \n,/' but found 'a' at line 1, column 0
        >>> x = _NamelistParser(" , a")
        >>> x._expect_separator()
        True
        >>> x._curr()
        u'a'
        >>> x = _NamelistParser(" / a")
        >>> x._expect_separator()
        False
        >>> x._curr()
        u'/'
        >>> x = _NamelistParser(" , ! Some stuff\n a")
        >>> x._expect_separator()
        True
        >>> x._curr()
        u'a'
        >>> x = _NamelistParser(" , ! Some stuff\n ! Other stuff\n a")
        >>> x._expect_separator()
        True
        >>> x._curr()
        u'a'
        >>> _NamelistParser("")._expect_separator(allow_eof=True)
        False
        >>> x = _NamelistParser(" ")
        >>> x._expect_separator(allow_eof=True)
        False
        >>> x = _NamelistParser(" ,")
        >>> x._expect_separator(allow_eof=True)
        True
        >>> x = _NamelistParser(" / ")
        >>> x._expect_separator(allow_eof=True)
        Traceback (most recent call last):
            ...
        _NamelistParseError: Error in parsing namelist: found group-terminating '/' in file without group names at line 1, column 1
        """
        errstring = "found group-terminating '/' in file without group names at "
        # Deal with the possibility that we are already at EOF.
        if allow_eof and self._pos == self._len:
            return False
        # Must actually be at a value separator.
        self._expect_char(' \n,/')
        try:
            self._eat_whitespace()
            if self._curr() == '/':
                if allow_eof:
                    raise _NamelistParseError(errstring +
                                              self._line_col_string())
                else:
                    return False
        except _NamelistEOF:
            if allow_eof:
                return False
            else:
                raise
        try:
            if self._curr() == ',':
                self._advance()
                self._eat_whitespace(allow_initial_comment=True)
        except _NamelistEOF:
            if not allow_eof:
                raise
        return True

    def _parse_name_and_values(self, allow_eof_end=False):
        r"""Parse and return a variable name and values assigned to that name.

        The return value of this function is a tuple containing (a) the name of
        the variable in a string, and (b) a list of the variable's values. Null
        values are represented by the empty string.

        If `allow_eof_end=True`, the end of the sequence of values might come
        from an empty string rather than a slash. (This is used for the
        alternate file format in "groupless" mode.)

        >>> _NamelistParser("foo='bar' /")._parse_name_and_values()
        (u'foo', [u"'bar'"])
        >>> _NamelistParser("foo ='bar' /")._parse_name_and_values()
        (u'foo', [u"'bar'"])
        >>> _NamelistParser("foo=\n'bar' /")._parse_name_and_values()
        (u'foo', [u"'bar'"])
        >>> _NamelistParser("foo 'bar' /")._parse_name_and_values()
        Traceback (most recent call last):
            ...
        _NamelistParseError: Error in parsing namelist: expected '=' but found "'" at line 1, column 4
        >>> _NamelistParser("foo='bar','bazz' /")._parse_name_and_values()
        (u'foo', [u"'bar'", u"'bazz'"])
        >>> _NamelistParser("foo=,,'bazz',6*/")._parse_name_and_values()
        (u'foo', [u'', u'', u"'bazz'", u'6*'])
        >>> _NamelistParser("foo='bar' 'bazz' foo2='ban'")._parse_name_and_values()
        (u'foo', [u"'bar'", u"'bazz'"])
        >>> _NamelistParser("foo= foo2='ban' ")._parse_name_and_values()
        Traceback (most recent call last):
            ...
        _NamelistParseError: Error in parsing namelist: expected literal value, but got "foo2='ban'" at line 1, column 15
        >>> _NamelistParser("foo=,,'bazz',6* ")._parse_name_and_values(allow_eof_end=True)
        (u'foo', [u'', u'', u"'bazz'", u'6*'])
        >>> _NamelistParser("foo=")._parse_name_and_values()
        Traceback (most recent call last):
            ...
        _NamelistEOF: Unexpected end of file encountered in namelist. (At line 1, column 4)
        >>> _NamelistParser("foo=")._parse_name_and_values(allow_eof_end=True)
        (u'foo', [u''])
        >>> _NamelistParser("foo= ")._parse_name_and_values(allow_eof_end=True)
        (u'foo', [u''])
        >>> _NamelistParser("foo=2")._parse_name_and_values(allow_eof_end=True)
        (u'foo', [u'2'])
        >>> _NamelistParser("foo=1,2")._parse_name_and_values(allow_eof_end=True)
        (u'foo', [u'1', u'2'])
        >>> _NamelistParser("foo=1,")._parse_name_and_values(allow_eof_end=True)
        (u'foo', [u'1', u''])
        """
        name = self._parse_variable_name()
        self._eat_whitespace()
        self._expect_char("=")
        try:
            self._advance()
            self._eat_whitespace()
        except _NamelistEOF:
            # If we hit the end of file, return a name assigned to a null value.
            if allow_eof_end:
                return name, [u'']
            else:
                raise
        # Expect at least one literal, even if it's a null value.
        values = [self._parse_literal(allow_eof_end=allow_eof_end)]
        # While we haven't reached the end of the namelist group...
        while self._expect_separator(allow_eof=allow_eof_end):
            # see if we can parse a literal (we might get a variable name)...
            literal = self._parse_literal(allow_name=True,
                                          allow_eof_end=allow_eof_end)
            if literal is None:
                break
            # and if it really is a literal, add it.
            values.append(literal)
        return name, values

    def _parse_namelist_group(self):
        r"""Parse an entire namelist group, adding info to `self._settings`.

        This function assumes that we start at the beginning of the group name
        (e.g. '&'), and will return at the end of the namelist group ('/').

        >>> x = _NamelistParser("&group /")
        >>> x._parse_namelist_group()
        >>> x._settings
        {u'group': {}}
        >>> x._curr()
        u'/'
        >>> x = _NamelistParser("&group\n foo='bar','bazz'\n,, foo2=2*5\n /")
        >>> x._parse_namelist_group()
        >>> x._settings
        {u'group': {u'foo': [u"'bar'", u"'bazz'", u''], u'foo2': [u'2*5']}}
        >>> x = _NamelistParser("&group\n foo='bar','bazz'\n,, foo2=2*5\n /", groupless=True)
        >>> x._parse_namelist_group()
        >>> x._settings
        {u'foo': [u"'bar'", u"'bazz'", u''], u'foo2': [u'2*5']}
        >>> x._curr()
        u'/'
        >>> x = _NamelistParser("&group /&group /")
        >>> x._parse_namelist_group()
        >>> x._advance()
        >>> x._parse_namelist_group()
        Traceback (most recent call last):
            ...
        _NamelistParseError: Error in parsing namelist: Namelist group 'group' encountered twice.
        >>> x = _NamelistParser("&group foo='bar', foo='bazz' /")
        >>> x._parse_namelist_group()
        >>> x._settings
        {u'group': {u'foo': [u"'bazz'"]}}
        >>> x = _NamelistParser("&group foo='bar', foo= /")
        >>> x._parse_namelist_group()
        >>> x._settings
        {u'group': {u'foo': [u"'bar'"]}}
        >>> x = _NamelistParser("&group foo='bar', foo= /", groupless=True)
        >>> x._parse_namelist_group()
        >>> x._settings
        {u'foo': [u"'bar'"]}
        """
        group_name = self._parse_namelist_group_name()
        if not self._groupless:
            # Make sure that this is the first time we've seen this group.
            if group_name in self._settings:
                raise _NamelistParseError("Namelist group %r encountered twice."
                                          % str(group_name))
            self._settings[group_name] = {}
        self._eat_whitespace()
        while self._curr() != '/':
            name, values = self._parse_name_and_values()
            if self._groupless:
                if name in self._settings:
                    values = merge_literal_lists(self._settings[name], values)
                self._settings[name] = values
            else:
                group = self._settings[group_name]
                if name in group:
                    values = merge_literal_lists(group[name], values)
                group[name] = values

    def parse_namelist(self):
        r"""Parse the contents of an entire namelist file.

        Returned information is a dictionary of dictionaries, mapping variables
        first by namelist group name, then by variable name.

        >>> _NamelistParser("").parse_namelist()
        {}
        >>> _NamelistParser(" \n!Comment").parse_namelist()
        {}
        >>> _NamelistParser(" &group /").parse_namelist()
        {u'group': {}}
        >>> _NamelistParser("! Comment \n &group /! Comment\n ").parse_namelist()
        {u'group': {}}
        >>> _NamelistParser("! Comment \n &group /! Comment ").parse_namelist()
        {u'group': {}}
        >>> _NamelistParser("&group1\n foo='bar','bazz'\n,, foo2=2*5\n / &group2 /").parse_namelist()
        {u'group1': {u'foo': [u"'bar'", u"'bazz'", u''], u'foo2': [u'2*5']}, u'group2': {}}
        >>> _NamelistParser("!blah \n foo='bar','bazz'\n,, foo2=2*5\n ", groupless=True).parse_namelist()
        {u'foo': [u"'bar'", u"'bazz'", u''], u'foo2': [u'2*5']}
        >>> _NamelistParser("!blah \n foo='bar','bazz'\n,, foo2=2*5,6\n ", groupless=True).parse_namelist()
        {u'foo': [u"'bar'", u"'bazz'", u''], u'foo2': [u'2*5', u'6']}
        >>> _NamelistParser("!blah \n foo='bar'", groupless=True).parse_namelist()
        {u'foo': [u"'bar'"]}
        >>> _NamelistParser("foo='bar', foo='bazz'", groupless=True).parse_namelist()
        {u'foo': [u"'bazz'"]}
        >>> _NamelistParser("foo='bar', foo=", groupless=True).parse_namelist()
        {u'foo': [u"'bar'"]}
        """
        # Return empty dictionary for empty files.
        if self._len == 0:
            return self._settings
        # Remove initial whitespace and comments, and return empty dictionary if
        # that's all we have.
        try:
            self._eat_whitespace(allow_initial_comment=True)
        except _NamelistEOF:
            return self._settings
        # Handle case with no namelist groups.
        if self._groupless and self._curr() != '&':
            while self._pos < self._len:
                name, values = self._parse_name_and_values(allow_eof_end=True)
                if name in self._settings:
                    values = merge_literal_lists(self._settings[name], values)
                self._settings[name] = values
            return self._settings
        # Loop over namelist groups in the file.
        while True:
            self._parse_namelist_group()
            # After each group, try to move forward to the next one. If we run
            # out of text, return what we've found.
            try:
                self._advance()
                self._eat_whitespace(allow_initial_comment=True)
            except _NamelistEOF:
                return self._settings
