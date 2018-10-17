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
import collections

# Disable these because this is our standard setup
# pylint: disable=wildcard-import,unused-wildcard-import

from CIME.XML.standard_module_setup import *
from CIME.utils import expect, string_in_list
import six

logger = logging.getLogger(__name__)

# Fortran syntax regular expressions.
# Variable names.
#FORTRAN_NAME_REGEX = re.compile(r"(^[a-z][a-z0-9_]{0,62})(\([+-]?\d*:?[+-]?\d*:?[+-]?\d*\))?$", re.IGNORECASE)
FORTRAN_NAME_REGEX = re.compile(r"""(^[a-z][a-z0-9_@]{0,62})                            #  The variable name
                                  (\(                                                   # begin optional index expression
                                  (([+-]?\d+)                                           # Single valued index
                                  |                                                     # or
                                  (([+-]?\d+)?:([+-]?\d+)?:?([+-]?\d+)?))               # colon seperated triplet
                                  \))?\s*$"""                                           # end optional index expression
                                , re.IGNORECASE | re.VERBOSE)

FORTRAN_LITERAL_REGEXES = {}
# Integer literals.
_int_re_string = r"(\+|-)?[0-9]+"
FORTRAN_LITERAL_REGEXES['integer'] = re.compile("^" + _int_re_string + "$")
# Real/complex literals.
_ieee_exceptional_re_string = r"inf(inity)?|nan(\([^)]+\))?"
_float_re_string = r"((\+|-)?([0-9]+(\.[0-9]*)?|\.[0-9]+)([ed]?{})?|{})".format(_int_re_string, _ieee_exceptional_re_string)
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
    4. We only handle a single dimension !!!

    >>> is_valid_fortran_name("")
    False
    >>> is_valid_fortran_name("a")
    True
    >>> is_valid_fortran_name("A")
    True
    >>> is_valid_fortran_name("A(4)")
    True
    >>> is_valid_fortran_name("A(::)")
    True
    >>> is_valid_fortran_name("A(1:2:3)")
    True
    >>> is_valid_fortran_name("A(1::)")
    True
    >>> is_valid_fortran_name("A(:-2:)")
    True
    >>> is_valid_fortran_name("A(1::+3)")
    True
    >>> is_valid_fortran_name("A(1,3)")
    False
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

def get_fortran_name_only(full_var):
    """ remove array section if any and return only the variable name
    >>> get_fortran_name_only('foo')
    'foo'
    >>> get_fortran_name_only('foo(3)')
    'foo'
    >>> get_fortran_name_only('foo(::)')
    'foo'
    >>> get_fortran_name_only('foo(1::)')
    'foo'
    >>> get_fortran_name_only('foo(:+2:)')
    'foo'
    >>> get_fortran_name_only('foo(::-3)')
    'foo'
    >>> get_fortran_name_only('foo(::)')
    'foo'
    """
    m = FORTRAN_NAME_REGEX.search(full_var)
    return m.group(1)

def get_fortran_variable_indices(varname, varlen=1, allow_any_len=False):
    """ get indices from a fortran namelist variable as a triplet of minindex, maxindex and step

    >>> get_fortran_variable_indices('foo(3)')
    (3, 3, 1)
    >>> get_fortran_variable_indices('foo(1:2:3)')
    (1, 2, 3)
    >>> get_fortran_variable_indices('foo(::)', varlen=4)
    (1, 4, 1)
    >>> get_fortran_variable_indices('foo(::2)', varlen=4)
    (1, 4, 2)
    >>> get_fortran_variable_indices('foo(::)', allow_any_len=True)
    (1, -1, 1)
    """
    m = FORTRAN_NAME_REGEX.search(varname)
    (minindex, maxindex, step) = (1, varlen, 1)

    if m.group(4) is not None:
        minindex = int(m.group(4))
        maxindex = minindex
        step = 1

    elif m.group(5) is not None:
        if m.group(6) is not None:
            minindex = int(m.group(6))
        if m.group(7) is not None:
            maxindex = int(m.group(7))
        if m.group(8) is not None:
            step = int(m.group(8))

    if allow_any_len and maxindex == minindex:
        maxindex = -1

    expect(step != 0,"Step size 0 not allowed")

    return (minindex, maxindex, step)

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
           "Invalid Fortran type for a namelist: {!r}".format(str(type_)))
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
    >>> shouldRaise(ValueError, literal_to_python_value, "nan(1234)")
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
               "{!r} is not a valid literal for any Fortran type.".format(str(literal)))
    else:
        # Check that type is valid.
        expect(is_valid_fortran_namelist_literal(type_, literal),
               "{!r} is not a valid literal of type {!r}.".format(str(literal), str(type_)))
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
    >>> expand_literal_list(['2*f*'])
    ['f*', 'f*']
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
    ['1', '2', 'f*', '3', '3', '3', '5']
    >>> compress_literal_list(['f*', 'f*'])
    ['f*', 'f*']
    """
    compressed = []
    if len(literals) == 0:
        return compressed
    # for right now do not compress
    do_compression = False
    if do_compression:
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
                if isinstance(old_literal, six.string_types):
                    compressed.append(rep_str + old_literal)
                else:
                    compressed.append(rep_str + str(old_literal))
                old_literal = literal
                num_reps = 1
        rep_str = str(num_reps) + '*' if num_reps > 1 else ''
        if isinstance(old_literal, six.string_types):
            compressed.append(rep_str + old_literal)
        else:
            compressed.append(rep_str + str(old_literal))
        return compressed
    else:
        for literal in literals:
            if isinstance(literal, six.string_types):
                compressed.append(literal)
            else:
                compressed.append(str(literal))
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
    ['true', 'false', 'false', 'true', 'true', 'false']
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
    if isinstance(in_file, six.string_types):
        logger.debug("Reading namelist at: {}".format(in_file))
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


def shouldRaise(eclass, method, *args, **kw):
    """
    A helper function to make doctests py3 compatible
    http://python3porting.com/problems.html#running-doctests
    """
    try:
        method(*args, **kw)
    except:
        e = sys.exc_info()[1]
        if not isinstance(e, eclass):
            raise
        return
    raise Exception("Expected exception %s not raised" %
                    str(eclass))



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
                expect(group_name is not None, " Got None in groups {}".format(groups))
                self._groups[group_name] = collections.OrderedDict()
                for variable_name in groups[group_name]:
                    self._groups[group_name][variable_name] = groups[group_name][variable_name]

    def clean_groups(self):
        self._groups = collections.OrderedDict()

    def get_group_names(self):
        """Return a list of all groups in the namelist.

        >>> Namelist().get_group_names()
        []
        >>> sorted(parse(text='&foo / &bar /').get_group_names())
        ['bar', 'foo']
        """
        return list(self._groups.keys())

    def get_variable_names(self, group_name):
        """Return a list of all variables in the given namelist group.

        If the specified group is not in the namelist, returns an empty list.

        >>> Namelist().get_variable_names('foo')
        []
        >>> x = parse(text='&foo bar=,bazz=true,bazz(2)=fred,bang=6*""/')
        >>> sorted(x.get_variable_names('fOo'))
        ['bang', 'bar', 'bazz', 'bazz(2)']
        >>> x = parse(text='&foo bar=,bazz=true,bang=6*""/')
        >>> sorted(x.get_variable_names('fOo'))
        ['bang', 'bar', 'bazz']
        >>> x = parse(text='&foo bar(::)=,bazz=false,bazz(2)=true,bazz(:2:)=6*""/')
        >>> sorted(x.get_variable_names('fOo'))
        ['bar(::)', 'bazz', 'bazz(2)', 'bazz(:2:)']
        """
        gn = string_in_list(group_name,self._groups)
        if not gn:
            return []
        return list(self._groups[gn].keys())

    def get_variable_value(self, group_name, variable_name):
        """Return the value of the specified variable.

        This function always returns a non-empty list containing strings. If the
        specified `group_name` or `variable_name` is not present, `['']` is
        returned.

        >>> Namelist().get_variable_value('foo', 'bar')
        ['']
        >>> parse(text='&foo bar=1,2 /').get_variable_value('foo', 'bazz')
        ['']
        >>> parse(text='&foo bar=1,2 /').get_variable_value('foO', 'Bar')
        ['1', '2']
        """
        gn = string_in_list(group_name,self._groups)
        if gn:
            vn = string_in_list(variable_name,self._groups[gn])
            if vn:
                return self._groups[gn][vn]
        return ['']

    def get_value(self, variable_name):
        """Return the value of a uniquely-named variable.

        This function is similar to `get_variable_value`, except that it does
        not require a `group_name`, and it requires that the `variable_name` be
        unique across all groups.

        >>> parse(text='&foo bar=1 / &bazz bar=1 /').get_value('bar')  # doctest: +ELLIPSIS
        Traceback (most recent call last):
        ...
        SystemExit: ERROR: Namelist.get_value: Variable {} is present in multiple groups: ...
        >>> parse(text='&foo bar=1 / &bazz /').get_value('Bar')
        ['1']
        >>> parse(text='&foo bar(2)=1 / &bazz /').get_value('Bar(2)')
        ['1']
        >>> parse(text='&foo / &bazz /').get_value('bar')
        ['']
        """
        possible_groups = []
        vn = None
        for group_name in self._groups:
            vnt = string_in_list(variable_name, self._groups[group_name])
            if vnt:
                vn = vnt
                possible_groups.append(group_name)
        expect(len(possible_groups) <= 1,
               "Namelist.get_value: Variable {} is present in multiple groups: "
               + str(possible_groups))
        if possible_groups:
            return self._groups[possible_groups[0]][vn]
        else:
            return ['']

    def set_variable_value(self, group_name, variable_name, value, var_size=1):
        """Set the value of the specified variable.

        >>> x = parse(text='&foo bar=1 /')
        >>> x.get_variable_value('foo', 'bar')
        ['1']
        >>> x.set_variable_value('foo', 'bar(2)', ['3'], var_size=4)
        >>> x.get_variable_value('foo', 'bar')
        ['1', '3']
        >>> x.set_variable_value('foo', 'bar(1)', ['2'])
        >>> x.get_variable_value('foo', 'bar')
        ['2', '3']
        >>> x.set_variable_value('foo', 'bar', ['1'])
        >>> x.get_variable_value('foo', 'bar')
        ['1', '3']
        >>> x.set_variable_value('foo', 'bazz', ['3'])
        >>> x.set_variable_value('Brack', 'baR', ['4'])
        >>> x.get_variable_value('foo', 'bazz')
        ['3']
        >>> x.get_variable_value('brack', 'bar')
        ['4']
        >>> x.set_variable_value('foo', 'red(2:6:2)', ['2', '4', '6'], var_size=12)
        >>> x.get_variable_value('foo', 'red')
        ['', '2', '', '4', '', '6']
        """
        minindex, maxindex, step = get_fortran_variable_indices(variable_name, var_size)
        variable_name = get_fortran_name_only(variable_name)

        expect(minindex > 0, "Indices < 1 not supported in CIME interface to fortran namelists... lower bound={}".format(minindex))
        gn = string_in_list(group_name,self._groups)
        if not gn:
            gn = group_name
            self._groups[gn] = {}

        tlen = 1
        vn = string_in_list(variable_name, self._groups[gn])
        if vn:
            tlen = len(self._groups[gn][vn])
        else:
            vn = variable_name
            tlen = 1
            self._groups[gn][vn] = ['']

        if minindex > tlen:
            self._groups[gn][vn].extend(['']*(minindex-tlen-1))

        for i in range(minindex, maxindex+2*step, step):
            while len(self._groups[gn][vn]) < i:
                self._groups[gn][vn].append('')
            self._groups[gn][vn][i-1] = value.pop(0)
            if len(value) == 0:
                break

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
        gn = string_in_list(group_name,self._groups)
        if gn:
            vn=string_in_list(variable_name,self._groups[gn])
            if vn:
                del self._groups[gn][vn]

    def merge_nl(self, other, overwrite=False):
        """Merge this namelist object with another.

        Values in the invoking (`self`) `Namelist` will take precedence over
        values in the `other` `Namelist`, unless `overwrite=True` is passed in,
        in which case `other` values take precedence.

        >>> x = parse(text='&foo bar=1 bazz=,2 brat=3/')
        >>> y = parse(text='&foo bar=2 bazz=3*1 baker=4 / &foo2 barter=5 /')
        >>> y.get_value('bazz')
        ['1', '1', '1']
        >>> x.merge_nl(y)
        >>> sorted(x.get_group_names())
        ['foo', 'foo2']
        >>> sorted(x.get_variable_names('foo'))
        ['baker', 'bar', 'bazz', 'brat']
        >>> sorted(x.get_variable_names('foo2'))
        ['barter']
        >>> x.get_value('bar')
        ['1']
        >>> x.get_value('bazz')
        ['1', '2', '1']
        >>> x.get_value('brat')
        ['3']
        >>> x.get_value('baker')
        ['4']
        >>> x.get_value('barter')
        ['5']
        >>> x = parse(text='&foo bar=1 bazz=,2 brat=3/')
        >>> y = parse(text='&foo bar=2 bazz=3*1 baker=4 / &foo2 barter=5 /')
        >>> x.merge_nl(y, overwrite=True)
        >>> sorted(x.get_group_names())
        ['foo', 'foo2']
        >>> sorted(x.get_variable_names('foo'))
        ['baker', 'bar', 'bazz', 'brat']
        >>> sorted(x.get_variable_names('foo2'))
        ['barter']
        >>> x.get_value('bar')
        ['2']
        >>> x.get_value('bazz')
        ['1', '1', '1']
        >>> x.get_value('brat')
        ['3']
        >>> x.get_value('baker')
        ['4']
        >>> x.get_value('barter')
        ['5']
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
                self.set_variable_value(group_name, variable_name, merged_val,
                                        var_size=len(merged_val))

    def get_group_variables(self, group_name):
        group_variables = {}
        group = self._groups[group_name]
        for name in sorted(group.keys()):
            value = group[name][0]
            group_variables[name] = value
        return group_variables

    def write(self, out_file, groups=None, append=False, format_='nml', sorted_groups=True,
              skip_comps=None, atm_cpl_dt=None, ocn_cpl_dt=None):
        """Write a the output data (normally fortran namelist) to the  out_file

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
        expect(format_ in ('nml', 'rc', 'nmlcontents', 'nuopc'),
               "Namelist.write: unexpected output format {!r}".format(str(format_)))
        if isinstance(out_file, six.string_types):
            logger.debug("Writing namelist to: {}".format(out_file))
            flag = 'a' if append else 'w'
            with open(out_file, flag) as file_obj:
                if format_ == 'nuopc':
                    self._write_nuopc(file_obj, groups, sorted_groups=sorted_groups,
                                      skip_comps=skip_comps, atm_cpl_dt=atm_cpl_dt, ocn_cpl_dt=ocn_cpl_dt)
                else:
                    self._write(file_obj, groups, format_, sorted_groups=sorted_groups)
        else:
            logger.debug("Writing namelist to file object")
            if format_ == 'nuopc':
                self._write_nuopc(out_file, groups, sorted_groups=sorted_groups,
                                  skip_comps=skip_comps, atm_cpl_dt=atm_cpl_dt, ocn_cpl_dt=ocn_cpl_dt)
            else:
                self._write(out_file, groups, format_, sorted_groups=sorted_groups)

    def _write(self, out_file, groups, format_, sorted_groups):
        """Unwrapped version of `write` assuming that a file object is input."""
        if groups is None:
            groups = list(self._groups.keys())
        if format_ == 'nml' or format_ == 'nmlcontents':
            equals = ' ='
        elif format_ == 'rc':
            equals = ':'
        if (sorted_groups):
            group_names = sorted(group for group in groups)
        else:
            group_names = groups
        for group_name in group_names:
            if format_ == 'nml':
                out_file.write("&{}\n".format(group_name))
            group = self._groups[group_name]
            for name in sorted(group.keys()):
                values = group[name]

                # @ is used in a namelist to put the same namelist variable in multiple groups
                # in the write phase, all characters in the namelist variable name after
                # the @ and including the @ should be removed
                if "@" in name:
                    name = re.sub('@.+$', "", name)

                # To prettify things for long lists of values, build strings
                # line-by-line.
                if values[0] == "True" or values[0] == "False":
                    values[0] = values[0].replace("True",".true.").replace("False",".false.")
                lines = ["  {}{} {}".format(name, equals, values[0])]
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
            if format_ == 'nmlcontents':
                out_file.write("\n")


    def _write_nuopc(self, out_file, groups, sorted_groups, skip_comps, atm_cpl_dt, ocn_cpl_dt):
        """Unwrapped version of `write` assuming that a file object is input."""

        if groups is None:
            groups = self._groups.keys()

        if (sorted_groups):
            group_names = sorted(group for group in groups)
        else:
            group_names = groups

        for group_name in group_names:
            if "_attributes" not in group_name and "nuopc_" not in group_name:
                continue

            if "_attributes" in group_name:
                out_file.write("{}::\n".format(group_name))

            group = self._groups[group_name]
            for name in sorted(group.keys()):
                values = group[name]
                if "component_list" in name:
                    for skip_comp in skip_comps:
                        if skip_comp in values[0]:
                            values[0] = values[0].replace(skip_comp,"")

                # @ is used in a namelist to put the same namelist variable in multiple groups
                # in the write phase, all characters in the namelist variable name after
                # the @ and including the @ should be removed
                if "@" in name:
                    name = re.sub('@.+$', "", name)

                equals = " ="
                if group_name == 'nuopc_runseq':
                    equals = '::\n       '
                elif "_var" in group_name:
                    equals = ':'

                # To prettify things for long lists of values, build strings
                # line-by-line.
                if values[0] == "True" or values[0] == "False":
                    values[0] = values[0].replace("True",".true.").replace("False",".false.")

                if "_attribute" in group_name:
                    lines = ["     {}{} {}".format(name, equals, values[0])]
                else:
                    lines = ["{}{} {}".format(name, equals, values[0])]

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
                    line = line.replace('"','')
                    # remove un-needed entries from the nuopc_runseq based
                    # on the prognostic_comps and skip_comps lists
                    if group_name == 'nuopc_runseq':
                        run_entries = line.splitlines()
                        newline = ""
                        for run_entry in run_entries:
                            print_entry = True
                            for skip_comp in skip_comps:
                                if "@" not in run_entry:
                                    if skip_comp in run_entry:
                                        print_entry = False
                                        logger.info("Writing nuopc_runseq, skipping {}".format(run_entry))
                                    if skip_comp.lower().strip() in run_entry:
                                        print_entry = False
                                        logger.info("Writing nuopc_runseq, skipping {}".format(run_entry))
                            if print_entry:
                                if "@atm_cpl_dt" in run_entry:
                                    run_entry = run_entry.replace("atm_cpl_dt",atm_cpl_dt)
                                if "@ocn_cpl_dt" in run_entry:
                                    run_entry = run_entry.replace("ocn_cpl_dt",ocn_cpl_dt)
                                newline += run_entry + "\n"
                        out_file.write(newline)
                    else:
                        out_file.write(line)

            if "_attribute" in group_name or "runseq" in group_name:
                out_file.write("::\n\n")

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
            string += " ({})".format(self._message)
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
            string += ": {}".format(self._message)
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
        self._text = str(text)
        self._len = len(self._text)
        # Dictionary with group names as keys, and dictionaries of variable
        # name-value pairs as values. (Or a single flat dictionary if
        # `groupless=True`.)
        self._settings = collections.OrderedDict()
        # Fortran allows setting a particular index of an array
        # such as foo(2)='k'
        # this dict is set to that value if used.
        self._groupless = groupless

    def _line_col_string(self):
        r"""Return a string specifying the current line and column number.

        >>> x = _NamelistParser('abc\nd\nef')
        >>> x._advance(5)
        >>> x._line_col_string()
        'line 2, column 1'
        """
        return "line {}, column {}".format(self._line, self._col)

    def _curr(self):
        """Return the character at the current position."""
        return self._text[self._pos]

    def _next(self):
        """Return the character at the next position.

        >>> shouldRaise(_NamelistEOF, _NamelistParser(' ')._next)

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
        >>> shouldRaise(_NamelistEOF, x._advance, 1)

        >>> shouldRaise(_NamelistEOF, _NamelistParser('abc\n')._advance, 4)

        >>> x = _NamelistParser('ab')
        >>> x._advance(check_eof=True)
        False
        >>> x._curr()
        'b'
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
            raise _NamelistEOF(message=None)

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
        'a'
        >>> x._eat_whitespace()
        False
        >>> x._advance()
        >>> shouldRaise(_NamelistEOF, x._eat_whitespace)

        >>> x = _NamelistParser(' \n! blah\n ! blah\n a')
        >>> x._eat_whitespace()
        True
        >>> x._curr()
        'a'
        >>> x = _NamelistParser('! blah\n a')
        >>> x._eat_whitespace()
        False
        >>> x._curr()
        '!'
        >>> x = _NamelistParser(' ! blah\n a')
        >>> x._eat_whitespace()
        True
        >>> x._curr()
        '!'
        >>> x = _NamelistParser(' ! blah\n a')
        >>> x._eat_whitespace(allow_initial_comment=True)
        True
        >>> x._curr()
        'a'
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
        ' '
        >>> x._eat_comment()
        False
        >>> x._eat_whitespace()
        True
        >>> x._eat_comment()
        True
        >>> x._curr()
        'a'
        >>> x._advance(2)
        >>> shouldRaise(_NamelistEOF, x._eat_comment)

        >>> x = _NamelistParser('! foo\n')
        >>> shouldRaise(_NamelistEOF, x._eat_comment)

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
        >>> shouldRaise(_NamelistParseError, x._expect_char, 'a')

        >>> x._expect_char('ab')
        """
        if self._curr() not in chars:
            if len(chars) == 1:
                char_description = repr(str(chars))
            else:
                char_description = "one of the characters in {!r}".format(str(chars))
            raise _NamelistParseError("expected {} but found {!r}".format(char_description, str(self._curr())))

    def _parse_namelist_group_name(self):
        r"""Parses and returns a namelist group name at the current position.

        >>> shouldRaise(_NamelistParseError, _NamelistParser('abc')._parse_namelist_group_name)

        >>> shouldRaise(_NamelistEOF, _NamelistParser('&abc')._parse_namelist_group_name)

        >>> _NamelistParser('&abc ')._parse_namelist_group_name()
        'abc'
        >>> _NamelistParser('&abc\n')._parse_namelist_group_name()
        'abc'
        >>> shouldRaise(_NamelistParseError, _NamelistParser('&abc/ ')._parse_namelist_group_name)

        >>> shouldRaise(_NamelistParseError, _NamelistParser('&abc= ')._parse_namelist_group_name)

        >>> shouldRaise(_NamelistParseError, _NamelistParser('& ')._parse_namelist_group_name)

        """
        self._expect_char("&")
        self._advance()
        return self._parse_variable_name(allow_equals=False)

    def _parse_variable_name(self, allow_equals=True):
        r"""Parses and returns a variable name at the current position.

        The `allow_equals` flag controls whether '=' can denote the end of the
        variable name; if it is `False`, only white space can be used for this
        purpose.

        >>> shouldRaise(_NamelistEOF, _NamelistParser('abc')._parse_variable_name)

        >>> _NamelistParser('foo(2)= ')._parse_variable_name()
        'foo(2)'
        >>> _NamelistParser('abc ')._parse_variable_name()
        'abc'
        >>> _NamelistParser('ABC ')._parse_variable_name()
        'abc'
        >>> _NamelistParser('abc\n')._parse_variable_name()
        'abc'
        >>> _NamelistParser('abc%fred\n')._parse_variable_name()
        'abc%fred'
        >>> _NamelistParser('abc(2)@fred\n')._parse_variable_name()
        'abc(2)@fred'
        >>> _NamelistParser('abc(1:2:3)\n')._parse_variable_name()
        'abc(1:2:3)'
        >>> _NamelistParser('abc=')._parse_variable_name()
        'abc'
        >>> try:
        ...     _NamelistParser('abc(1,2) ')._parse_variable_name()
        ...     raise AssertionError("_NamelistParseError not raised")
        ... except _NamelistParseError:
        ...    pass
        >>> try:
        ...     _NamelistParser('abc, ')._parse_variable_name()
        ...     raise AssertionError("_NamelistParseError not raised")
        ... except _NamelistParseError:
        ...    pass
        >>> try:
        ...    _NamelistParser(' ')._parse_variable_name()
        ...    raise AssertionError("_NamelistParseError not raised")
        ... except _NamelistParseError:
        ...    pass
        >>> _NamelistParser('foo+= ')._parse_variable_name()
        'foo'
        """
        old_pos = self._pos
        separators = (' ', '\n', '=', '+') if allow_equals else (' ', '\n')
        while self._curr() not in separators:
            self._advance()
        text = self._text[old_pos:self._pos]
        if '(' in text:
            expect(')' in text,"Parsing error ")
        elif ')' in text:
            expect(False,"Parsing error ")

        # @ is used in a namelist to put the same namelist variable in multiple groups
        # in the write phase, all characters in the namelist variable name after
        # the @ and including the @ should be removed
        if "%" in text:
            text_check = re.sub('%.+$', "", text)
        elif "@" in text:
            text_check = re.sub('@.+$', "", text)
        else:
            text_check = text

        if not is_valid_fortran_name(text_check):
            if re.search(r".*\(.*\,.*\)", text_check):
                err_str = "Multiple dimensions not supported in CIME namelist variables {!r}".format(str(text))
            else:
                err_str = "{!r} is not a valid variable name".format(str(text))
            raise _NamelistParseError(err_str)
        name = text.lower()

        return name

    def _parse_character_literal(self):
        """Parse and return a character literal (a string).

        Position on return is the last character of the string; we avoid
        advancing past that in order to avoid potential EOF errors.

        >>> shouldRaise(_NamelistEOF, _NamelistParser('"abc')._parse_character_literal)

        >>> _NamelistParser('"abc" ')._parse_character_literal()
        '"abc"'
        >>> _NamelistParser("'abc' ")._parse_character_literal()
        "'abc'"
        >>> shouldRaise(_NamelistParseError, _NamelistParser("*abc* ")._parse_character_literal)

        >>> _NamelistParser("'abc''def' ")._parse_character_literal()
        "'abc''def'"
        >>> _NamelistParser("'abc''' ")._parse_character_literal()
        "'abc'''"
        >>> _NamelistParser("'''abc' ")._parse_character_literal()
        "'''abc'"
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
            raise _NamelistParseError("{} is not a valid character literal".format(text))
        return text

    def _parse_complex_literal(self):
        """Parse and return a complex literal.

        Position on return is the last character of the string; we avoid
        advancing past that in order to avoid potential EOF errors.

        >>> shouldRaise(_NamelistEOF, _NamelistParser('(1.,2.')._parse_complex_literal)

        >>> _NamelistParser('(1.,2.) ')._parse_complex_literal()
        '(1.,2.)'
        >>> shouldRaise(_NamelistParseError, _NamelistParser("(A,B) ")._parse_complex_literal)

        """
        old_pos = self._pos
        while self._curr() != ')':
            self._advance()
        text = self._text[old_pos:self._pos+1]
        if not is_valid_fortran_namelist_literal("complex", text):
            raise _NamelistParseError("{!r} is not a valid complex literal".format(str(text)))
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

    def _look_ahead_for_plusequals(self, pos):
        r"""Look ahead to see if the next two non-whitespace character are '+='.

        The `pos` argument is the position in the text to start from while
        looking. This function returns a boolean.

        >>> _NamelistParser('+=')._look_ahead_for_plusequals(0)
        True
        >>> _NamelistParser('a \n+=')._look_ahead_for_plusequals(1)
        True
        >>> _NamelistParser('')._look_ahead_for_plusequals(0)
        False
        >>> _NamelistParser('a+=')._look_ahead_for_plusequals(0)
        False
        """
        for test_pos in range(pos, self._len):
            if self._text[test_pos] not in (' ', '\n'):
                if self._text[test_pos] == '+':
                    return self._look_ahead_for_equals(test_pos + 1)
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
        '"abc"'
        >>> _NamelistParser("'abc' ")._parse_literal()
        "'abc'"
        >>> shouldRaise(_NamelistEOF, _NamelistParser('"abc"')._parse_literal)

        >>> _NamelistParser('"abc"')._parse_literal(allow_eof_end=True)
        '"abc"'
        >>> _NamelistParser('(1.,2.) ')._parse_literal()
        '(1.,2.)'
        >>> shouldRaise(_NamelistEOF, _NamelistParser('(1.,2.)')._parse_literal)

        >>> _NamelistParser('(1.,2.)')._parse_literal(allow_eof_end=True)
        '(1.,2.)'
        >>> _NamelistParser('5 ')._parse_literal()
        '5'
        >>> _NamelistParser('6.9 ')._parse_literal()
        '6.9'
        >>> _NamelistParser('inf ')._parse_literal()
        'inf'
        >>> _NamelistParser('nan(booga) ')._parse_literal()
        'nan(booga)'
        >>> _NamelistParser('.FLORIDA$ ')._parse_literal()
        '.FLORIDA$'
        >>> shouldRaise(_NamelistParseError, _NamelistParser('hamburger ')._parse_literal)

        >>> _NamelistParser('5,')._parse_literal()
        '5'
        >>> _NamelistParser('5\n')._parse_literal()
        '5'
        >>> _NamelistParser('5/')._parse_literal()
        '5'
        >>> _NamelistParser(',')._parse_literal()
        ''
        >>> _NamelistParser('6*5 ')._parse_literal()
        '6*5'
        >>> _NamelistParser('6*(1., 2.) ')._parse_literal()
        '6*(1., 2.)'
        >>> _NamelistParser('6*"a" ')._parse_literal()
        '6*"a"'
        >>> shouldRaise(_NamelistEOF, _NamelistParser('6*')._parse_literal)

        >>> _NamelistParser('6*')._parse_literal(allow_eof_end=True)
        '6*'
        >>> shouldRaise(_NamelistParseError, _NamelistParser('foo= ')._parse_literal)

        >>> shouldRaise(_NamelistParseError, _NamelistParser('foo+= ')._parse_literal)

        >>> _NamelistParser('5,')._parse_literal(allow_name=True)
        '5'
        >>> x = _NamelistParser('foo= ')
        >>> x._parse_literal(allow_name=True)
        >>> x._curr()
        'f'
        >>> x = _NamelistParser('foo+= ')
        >>> x._parse_literal(allow_name=True)
        >>> x._curr()
        'f'
        >>> shouldRaise(_NamelistParseError, _NamelistParser('6*foo= ')._parse_literal, allow_name=True)

        >>> shouldRaise(_NamelistParseError, _NamelistParser('6*foo+= ')._parse_literal, allow_name=True)

        >>> x = _NamelistParser('foo = ')
        >>> x._parse_literal(allow_name=True)
        >>> x._curr()
        'f'
        >>> x = _NamelistParser('foo\n= ')
        >>> x._parse_literal(allow_name=True)
        >>> x._curr()
        'f'
        >>> _NamelistParser('')._parse_literal(allow_eof_end=True)
        ''
        """
        # Deal with empty input string.
        if allow_eof_end and self._pos == self._len:
            return ''
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
            separators.append('+')
        while new_pos != self._len and self._text[new_pos] not in separators:
            # allow commas if they are inside ()
            if self._text[new_pos] == '(':
                separators.remove(',')
            elif self._text[new_pos] == ')':
                separators.append(',')
            new_pos += 1

        if not allow_eof_end and new_pos == self._len:
            # At the end of the file, give up by throwing an EOF.
            self._advance(self._len)
        # If `allow_name` is set, we need to check and see if the next non-blank
        # character is '=' or the next two are '+=', and return `None` if so.
        if allow_name and self._look_ahead_for_equals(new_pos):
            return
        elif allow_name and self._look_ahead_for_plusequals(new_pos):
            return

        self._advance(new_pos - self._pos, check_eof=allow_eof_end)
        text = self._text[old_pos:self._pos]
        if not any(is_valid_fortran_namelist_literal(type_, text)
                   for type_ in ("integer", "logical", "real")):
            raise _NamelistParseError("expected literal value, but got {!r}".format(str(text)))
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
        'a'
        >>> x = _NamelistParser(" a")
        >>> x._expect_separator()
        True
        >>> x._curr()
        'a'
        >>> x = _NamelistParser(",a")
        >>> x._expect_separator()
        True
        >>> x._curr()
        'a'
        >>> x = _NamelistParser("/a")
        >>> x._expect_separator()
        False
        >>> x._curr()
        '/'
        >>> x = _NamelistParser("a")
        >>> shouldRaise(_NamelistParseError, x._expect_separator)

        >>> x = _NamelistParser(" , a")
        >>> x._expect_separator()
        True
        >>> x._curr()
        'a'
        >>> x = _NamelistParser(" / a")
        >>> x._expect_separator()
        False
        >>> x._curr()
        '/'
        >>> x = _NamelistParser(" , ! Some stuff\n a")
        >>> x._expect_separator()
        True
        >>> x._curr()
        'a'
        >>> x = _NamelistParser(" , ! Some stuff\n ! Other stuff\n a")
        >>> x._expect_separator()
        True
        >>> x._curr()
        'a'
        >>> _NamelistParser("")._expect_separator(allow_eof=True)
        False
        >>> x = _NamelistParser(" ")
        >>> x._expect_separator(allow_eof=True)
        False
        >>> x = _NamelistParser(" ,")
        >>> x._expect_separator(allow_eof=True)
        True
        >>> x = _NamelistParser(" / ")
        >>> shouldRaise(_NamelistParseError, x._expect_separator, allow_eof=True)

        """
        errstring = "found group-terminating '/' in file without group names"
        # Deal with the possibility that we are already at EOF.
        if allow_eof and self._pos == self._len:
            return False
        # Must actually be at a value separator.
        self._expect_char(' \n,/')
        try:
            self._eat_whitespace()
            if self._curr() == '/':
                if allow_eof:
                    raise _NamelistParseError(errstring)
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
        the variable in a string, (b) a list of the variable's values, and
        (c) whether or not to add the found value to existing variable. Null
        values are represented by the empty string.

        If `allow_eof_end=True`, the end of the sequence of values might come
        from an empty string rather than a slash. (This is used for the
        alternate file format in "groupless" mode.)

        >>> _NamelistParser("foo='bar' /")._parse_name_and_values()
        ('foo', ["'bar'"], False)
        >>> _NamelistParser("foo(3)='bar' /")._parse_name_and_values()
        ('foo(3)', ["'bar'"], False)
        >>> _NamelistParser("foo ='bar' /")._parse_name_and_values()
        ('foo', ["'bar'"], False)
        >>> _NamelistParser("foo=\n'bar' /")._parse_name_and_values()
        ('foo', ["'bar'"], False)
        >>> shouldRaise(_NamelistParseError, _NamelistParser("foo 'bar' /")._parse_name_and_values)

        >>> _NamelistParser("foo='bar','bazz' /")._parse_name_and_values()
        ('foo', ["'bar'", "'bazz'"], False)
        >>> _NamelistParser("foo=,,'bazz',6*/")._parse_name_and_values()
        ('foo', ['', '', "'bazz'", '6*'], False)
        >>> _NamelistParser("foo='bar' 'bazz' foo2='ban'")._parse_name_and_values()
        ('foo', ["'bar'", "'bazz'"], False)
        >>> _NamelistParser("foo='bar' 'bazz' foo2(2)='ban'")._parse_name_and_values()
        ('foo', ["'bar'", "'bazz'"], False)
        >>> shouldRaise(_NamelistParseError, _NamelistParser("foo= foo2='ban' ")._parse_name_and_values)

        >>> _NamelistParser("foo=,,'bazz',6* ")._parse_name_and_values(allow_eof_end=True)
        ('foo', ['', '', "'bazz'", '6*'], False)
        >>> _NamelistParser("foo(3)='bazz'")._parse_name_and_values(allow_eof_end=True)
        ('foo(3)', ["'bazz'"], False)
        >>> shouldRaise(_NamelistEOF, _NamelistParser("foo=")._parse_name_and_values)

        >>> _NamelistParser("foo=")._parse_name_and_values(allow_eof_end=True)
        ('foo', [''], False)
        >>> _NamelistParser("foo= ")._parse_name_and_values(allow_eof_end=True)
        ('foo', [''], False)
        >>> _NamelistParser("foo=2")._parse_name_and_values(allow_eof_end=True)
        ('foo', ['2'], False)
        >>> _NamelistParser("foo=1,2")._parse_name_and_values(allow_eof_end=True)
        ('foo', ['1', '2'], False)
        >>> _NamelistParser("foo(1:2)=1,2,3 ")._parse_name_and_values(allow_eof_end=True)
        Traceback (most recent call last):
        ...
        SystemExit: ERROR: Too many values for array foo(1:2)
        >>> _NamelistParser("foo=1,")._parse_name_and_values(allow_eof_end=True)
        ('foo', ['1', ''], False)
        >>> _NamelistParser("foo+=1")._parse_name_and_values(allow_eof_end=True)
        ('foo', ['1'], True)
        """
        name = self._parse_variable_name()
        addto = False  # This keeps track of whether += existed

        self._eat_whitespace()
        # check to see if we have a "+="
        if self._curr() == '+':
            self._advance()
            addto=True  # tell parser that we want to add to dictionary values
        self._expect_char("=")
        try:
            self._advance()
            self._eat_whitespace()
        except _NamelistEOF:
            # If we hit the end of file, return a name assigned to a null value.
            if allow_eof_end:
                return name, [''], addto
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
        (minindex, maxindex, step) = get_fortran_variable_indices(name,allow_any_len=True)
        if (minindex > 1 or maxindex > minindex or step > 1) and maxindex > 0:
            arraylen =max(0,1 + ((maxindex - minindex)/step))
            expect(len(values) <= arraylen, "Too many values for array {}".format(name))

        return name, values, addto

    def _parse_namelist_group(self):
        r"""Parse an entire namelist group, adding info to `self._settings`.

        This function assumes that we start at the beginning of the group name
        (e.g. '&'), and will return at the end of the namelist group ('/').

        >>> x = _NamelistParser("&group /")
        >>> x._parse_namelist_group()
        >>> x._settings
        OrderedDict([('group', {})])
        >>> x._curr()
        '/'
        >>> x = _NamelistParser("&group\n foo='bar','bazz'\n,, foo2=2*5\n /")
        >>> x._parse_namelist_group()
        >>> x._settings
        OrderedDict([('group', {'foo': ["'bar'", "'bazz'", ''], 'foo2': ['5', '5']})])
        >>> x = _NamelistParser("&group\n foo='bar','bazz'\n,, foo2=2*5\n /", groupless=True)
        >>> x._parse_namelist_group()
        >>> x._settings
        OrderedDict([('foo', ["'bar'", "'bazz'", '']), ('foo2', ['5', '5'])])
        >>> x._curr()
        '/'
        >>> x = _NamelistParser("&group /&group /")
        >>> x._parse_namelist_group()
        >>> x._advance()
        >>> shouldRaise(_NamelistParseError, x._parse_namelist_group)

        >>> x = _NamelistParser("&group foo='bar', foo='bazz' /")
        >>> x._parse_namelist_group()
        >>> x._settings
        OrderedDict([('group', {'foo': ["'bazz'"]})])
        >>> x = _NamelistParser("&group foo='bar', foo= /")
        >>> x._parse_namelist_group()
        >>> x._settings
        OrderedDict([('group', {'foo': ["'bar'"]})])
        >>> x = _NamelistParser("&group foo='bar', foo= /", groupless=True)
        >>> x._parse_namelist_group()
        >>> x._settings
        OrderedDict([('foo', ["'bar'"])])
        >>> x = _NamelistParser("&group foo='bar', foo+='baz' /", groupless=True)
        >>> x._parse_namelist_group()
        >>> x._settings
        OrderedDict([('foo', ["'bar'", "'baz'"])])
        >>> x = _NamelistParser("&group foo+='bar' /", groupless=True)
        >>> x._parse_namelist_group()
        >>> x._settings
        OrderedDict([('foo', ["'bar'"])])
        >>> x = _NamelistParser("&group foo='bar', foo+='baz' /")
        >>> x._parse_namelist_group()
        >>> x._settings
        OrderedDict([('group', {'foo': ["'bar'", "'baz'"]})])
        >>> x = _NamelistParser("&group foo+='bar' /")
        >>> x._parse_namelist_group()
        >>> x._settings
        OrderedDict([('group', {'foo': ["'bar'"]})])
        """
        group_name = self._parse_namelist_group_name()
        if not self._groupless:
            # Make sure that this is the first time we've seen this group.
            if group_name in self._settings:
                raise _NamelistParseError("Namelist group {!r} encountered twice.".format(str(group_name)))
            self._settings[group_name] = {}
        self._eat_whitespace()
        while self._curr() != '/':
            name, values, addto = self._parse_name_and_values()
            dsettings = []
            if self._groupless:
                if name in self._settings:
                    dsettings = self._settings[name]
                    if addto:
                        values = self._settings[name] + values
                if not addto:
                    values = merge_literal_lists(dsettings, values)
                self._settings[name] = values
            else:
                group = self._settings[group_name]
                if name in group:
                    dsettings = group[name]
                    if addto:
                        values = group[name] + values
                if not addto:
                    values = merge_literal_lists(dsettings, values)
                group[name] = values

    def parse_namelist(self):
        r"""Parse the contents of an entire namelist file.

        Returned information is a dictionary of dictionaries, mapping variables
        first by namelist group name, then by variable name.

        >>> _NamelistParser("").parse_namelist()
        OrderedDict()
        >>> _NamelistParser(" \n!Comment").parse_namelist()
        OrderedDict()
        >>> _NamelistParser(" &group /").parse_namelist()
        OrderedDict([('group', {})])
        >>> _NamelistParser("! Comment \n &group /! Comment\n ").parse_namelist()
        OrderedDict([('group', {})])
        >>> _NamelistParser("! Comment \n &group /! Comment ").parse_namelist()
        OrderedDict([('group', {})])
        >>> _NamelistParser("&group1\n foo='bar','bazz'\n,, foo2=2*5\n / &group2 /").parse_namelist()
        OrderedDict([('group1', {'foo': ["'bar'", "'bazz'", ''], 'foo2': ['5', '5']}), ('group2', {})])
        >>> _NamelistParser("!blah \n foo='bar','bazz'\n,, foo2=2*5\n ", groupless=True).parse_namelist()
        OrderedDict([('foo', ["'bar'", "'bazz'", '']), ('foo2', ['2*5'])])
        >>> _NamelistParser("!blah \n foo='bar','bazz'\n,, foo2=2*5,6\n ", groupless=True).parse_namelist()
        OrderedDict([('foo', ["'bar'", "'bazz'", '']), ('foo2', ['2*5', '6'])])
        >>> _NamelistParser("!blah \n foo='bar'", groupless=True).parse_namelist()
        OrderedDict([('foo', ["'bar'"])])
        >>> _NamelistParser("foo='bar', foo(3)='bazz'", groupless=True).parse_namelist()
        OrderedDict([('foo', ["'bar'"]), ('foo(3)', ["'bazz'"])])
        >>> _NamelistParser("foo(2)='bar'", groupless=True).parse_namelist()
        OrderedDict([('foo(2)', ["'bar'"])])
        >>> _NamelistParser("foo(2)='bar', foo(3)='bazz'", groupless=True).parse_namelist()
        OrderedDict([('foo(2)', ["'bar'"]), ('foo(3)', ["'bazz'"])])
        >>> _NamelistParser("foo='bar', foo='bazz'", groupless=True).parse_namelist()
        OrderedDict([('foo', ["'bazz'"])])
        >>> _NamelistParser("foo='bar'\n foo+='bazz'", groupless=True).parse_namelist()
        OrderedDict([('foo', ["'bar'", "'bazz'"])])
        >>> _NamelistParser("foo='bar', foo='bazz'", groupless=True).parse_namelist()
        OrderedDict([('foo', ["'bazz'"])])
        >>> _NamelistParser("foo='bar', foo=", groupless=True).parse_namelist()
        OrderedDict([('foo', ["'bar'"])])
        >>> _NamelistParser("foo='bar', 'bazz'\n foo+='ban'", groupless=True).parse_namelist()
        OrderedDict([('foo', ["'bar'", "'bazz'", "'ban'"])])
        >>> _NamelistParser("foo+='bar'", groupless=True).parse_namelist()
        OrderedDict([('foo', ["'bar'"])])
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
                name, values, addto = self._parse_name_and_values(allow_eof_end=True)
                if name in self._settings:
                    if addto:
                        values = self._settings[name] + values
                    else:
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
