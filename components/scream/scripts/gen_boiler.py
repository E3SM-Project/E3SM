from utils import expect, get_git_toplevel_dir

from collections import OrderedDict
import pathlib, re, os

#
# Global hardcoded data
#

# Templates: maps piece name to generic file text
FILE_TEMPLATES = {
    "cxx_bfb_unit_impl": lambda phys, sub, gen_code:
"""#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/{physics}/{physics}_functions.hpp"
#include "physics/{physics}/{physics}_functions_f90.hpp"

#include "{physics}_unit_tests_common.hpp"

namespace scream {{
namespace {physics} {{
namespace unit_test {{

template <typename D>
struct UnitWrap::UnitTest<D>::{test_data_struct} {{

{gen_code}

}};

}} // namespace unit_test
}} // namespace {physics}
}} // namespace scream

namespace {{

TEST_CASE("{sub}_bfb", "[{physics}]")
{{
  using TestStruct = scream::{physics}::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::{test_data_struct};

  TestStruct::run_bfb();
}}

}} // empty namespace
""".format(physics=phys, sub=sub, test_data_struct=get_data_test_struct_name(sub), gen_code=gen_code),

###############################################################################

    "cxx_func_impl": lambda phys, sub, gen_code:
"""#ifndef {phys_upper}_{sub_upper}_IMPL_HPP
#define {phys_upper}_{sub_upper}_IMPL_HPP

#include "{physics}_functions.hpp" // for ETI only but harmless for GPU

namespace scream {{
namespace {physics} {{

/*
 * Implementation of {physics} {sub}. Clients should NOT
 * #include this file, but include {physics}_functions.hpp instead.
 */

template<typename S, typename D>
{gen_code}

}} // namespace {physics}
}} // namespace scream

#endif
""".format(physics=phys, sub=sub, gen_code=gen_code, phys_upper=phys.upper(), sub_upper=sub.upper()),

###############################################################################

    "cxx_eti": lambda phys, sub, gen_code: gen_code
}

# piece map. maps the name of a piece of boilerplate that needs to be genereate to:
#   (filepath (relative to cxx_root))
FILEPATH, FILECREATE, INSERT_REGEX, ID_SELF_BEGIN_REGEX, ID_SELF_END_REGEX, DESC = range(6)
PIECES = OrderedDict([
    ("f90_c2f_bind", (
        lambda phys, sub, gb: "{}_iso_c.f90".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "f90_c2f_bind"),
        lambda phys, sub, gb: re.compile(r"^\s*end\s+module\s{}_iso_c".format(phys)), # put at end of module
        lambda phys, sub, gb: get_subroutine_begin_regex(sub + "_c"), # sub_c begin
        lambda phys, sub, gb: get_subroutine_end_regex(sub + "_c"),    # sub_c end
        lambda *x           : "The c to f90 fortran subroutine(<name>_c)"
    )),

    ("f90_f2c_bind"  , (
        lambda phys, sub, gb: "{}_iso_f.f90".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "f90_f2c_bind"),
        lambda phys, sub, gb: re.compile(r"^\s*end\s+interface"), # put at end of interface
        lambda phys, sub, gb: get_subroutine_begin_regex(sub + "_f"), # sub_f begin
        lambda phys, sub, gb: get_subroutine_end_regex(sub + "_f"),   # sub_f begin
        lambda *x           : "The f90 to c fortran subroutine(<name>_f)"
    )),

    ("cxx_c2f_bind_decl"  , (
        lambda phys, sub, gb: "{}_functions_f90.cpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_bind_decl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(comment='extern "C" : end _c decls'), # reqs special comment
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub + "_c"), # cxx_c decl
        lambda phys, sub, gb: re.compile(r".*;\s*$"),                   # ;
        lambda *x           : "The c to f90 c function declaration(<name>_c)"
    )),

    ("cxx_c2f_glue_decl"  , (
        lambda phys, sub, gb: "{}_functions_f90.hpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_glue_decl"),
        lambda phys, sub, gb: re.compile(r'^\s*extern\s+"C"'), # put before _f decls
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub), # cxx(data) decl
        lambda phys, sub, gb: re.compile(r".*;\s*"),             # ;
        lambda *x           : "The cxx to c function declaration(<name>(Data))"
    )),

    ("cxx_c2f_glue_impl" , (
        lambda phys, sub, gb: "{}_functions_f90.cpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_glue_impl"),
        lambda phys, sub, gb: re.compile(r"^\s*// end _c impls"), # reqs special comment
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub), # cxx(data)
        lambda phys, sub, gb: get_cxx_close_block_regex(at_line_start=True), # terminating }
        lambda *x           : "The cxx to c function implementation(<name>(Data))"
    )),

    ("cxx_c2f_data"  , (
        lambda phys, sub, gb: "{}_functions_f90.hpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_data"),
        lambda phys, sub, gb: re.compile(r"^\s*// Glue functions to call fortran"),  # reqs special comment
        lambda phys, sub, gb: get_cxx_struct_begin_regex(get_data_struct_name(sub)), # struct Sub
        lambda phys, sub, gb: get_cxx_close_block_regex(semicolon=True),             # terminating };
        lambda *x           : "The cxx data struct definition(struct Data)"
    )),

    ("cxx_f2c_bind_decl"  , (
        lambda phys, sub, gb: "{}_functions_f90.hpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_f2c_bind_decl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(comment="end _f function decls"), # reqs special comment
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub + "_f"), # cxx_f decl
        lambda phys, sub, gb: re.compile(r".*;\s*$"),                   # ;
        lambda *x           : "The f90 to cxx function declaration(<name>_f)"
    )),

    ("cxx_f2c_bind_impl"  , (
        lambda phys, sub, gb: "{}_functions_f90.cpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_f2c_bind_impl"),
        lambda phys, sub, gb: get_namespace_close_regex(phys),          # insert at end of namespace
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub + "_f"),      # cxx_f
        lambda phys, sub, gb: get_cxx_close_block_regex(at_line_start=True), # terminating }
        lambda *x           : "The f90 to cxx function implementation(<name>_f)"
    )),

    ("cxx_func_decl", (
        lambda phys, sub, gb: "{}_functions.hpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_func_decl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(semicolon=True, comment="struct Functions"), # end of struct, reqs special comment
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub, static=True), # cxx decl
        lambda phys, sub, gb: re.compile(r".*;\s*$"),            # ;
        lambda *x           : "The cxx kokkos function declaration(<name>)"
    )),

    ("cxx_incl_impl", (
        lambda phys, sub, gb: "{}_functions.hpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_incl_impl"),
        lambda phys, sub, gb: re.compile(r"^\s*#\s*endif\s+//\s*KOKKOS_ENABLE_CUDA"), # insert at end of impl includes, reqs special comment
        lambda phys, sub, gb: re.compile(r'^\s*#\s*include\s+"{}"'.format(get_piece_data(phys, sub, "cxx_func_impl", FILEPATH, gb))),
        lambda phys, sub, gb: re.compile(r".*"),
        lambda *x           : "The include of *impl.hpp file at bottom of main hpp"
    )),

    ("cxx_func_impl", (
        lambda phys, sub, gb: "{}_{}_impl.hpp".format(phys, sub),
        lambda phys, sub, gb: create_template(phys, sub, gb, "cxx_func_impl"),
        lambda phys, sub, gb: get_namespace_close_regex(phys),   # insert at end of namespace
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub, template="Functions<S,D>"), # cxx begin
        lambda phys, sub, gb: get_cxx_close_block_regex(at_line_start=True), # terminating }
        lambda *x           : "The cxx kokkos function stub implementation(<name>)"
    )),

    ("cxx_bfb_unit_decl", (
        lambda phys, sub, gb: "tests/{}_unit_tests_common.hpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_bfb_unit_decl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(semicolon=True), # insert at end of test struct
        lambda phys, sub, gb: get_cxx_struct_begin_regex(get_data_test_struct_name(sub)), # struct decl
        lambda phys, sub, gb: re.compile(r".*;\s*$"), # end of struct decl
        lambda *x           : "The cxx unit test struct declaration"
    )),

    ("cxx_bfb_unit_impl", (
        lambda phys, sub, gb: "tests/{}_{}_tests.cpp".format(phys, sub),
        lambda phys, sub, gb: create_template(phys, sub, gb, "cxx_bfb_unit_impl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(semicolon=True, at_line_start=True), # insert of end of struct
        lambda phys, sub, gb: get_cxx_function_begin_regex("run_bfb", static=True),  # run_bfb
        lambda phys, sub, gb: get_cxx_close_block_regex(comment="run_bfb"), # } // run_bfb # reqs special comment
        lambda *x           : "The cxx bfb unit test implementation"
    )),

    ("cxx_eti", (
        lambda phys, sub, gb: "{}_{}.cpp".format(phys, sub),
        lambda phys, sub, gb: create_template(phys, sub, gb, "cxx_eti"),
        lambda phys, sub, gb: re.compile(".*"), # insert at top of file
        lambda phys, sub, gb: re.compile(".*"), # start at top of file
        lambda phys, sub, gb: get_namespace_close_regex("scream"), #end of file
        lambda *x           : "The cxx explicit template instatiation file"
    )),

    ("cmake_impl_eti", (
        lambda phys, sub, gb: "CMakeLists.txt",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cmake_impl_eti"),
        lambda phys, sub, gb: re.compile(r".*[)]\s*#\s*{} ETI SRCS".format(phys.upper())), # insert at end of ETI src list, reqs special comment
        lambda phys, sub, gb: re.compile(r".*{}".format(get_piece_data(phys, sub, "cxx_eti", FILEPATH, gb))),
        lambda phys, sub, gb: re.compile(".*"),
        lambda *x           : "Make cmake aware of the ETI file if not cuda build"
    )),

    ("cmake_unit_test", (
        lambda phys, sub, gb: "tests/CMakeLists.txt",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cmake_unit_test"),
        lambda phys, sub, gb: re.compile(r".*[)]\s*#\s*{}_TESTS_SRCS".format(phys.upper())), # insert at end of test src list, reqs special comment
        lambda phys, sub, gb: re.compile(r".*{}".format(os.path.basename(get_piece_data(phys, sub, "cxx_bfb_unit_impl", FILEPATH, gb)))),
        lambda phys, sub, gb: re.compile(".*"),
        lambda *x           : "Make cmake aware of the unit test"
    )),

])

# physics map. maps the name of a physics packages containing the original fortran subroutines to:
#   (path-to-origin, path-to-cxx-src)
ORIGIN_FILE, CXX_ROOT, INIT_CODE = range(3)
PHYSICS = {
    "p3"   : (
        "components/cam/src/physics/cam/micro_p3.F90",
        "components/scream/src/physics/p3",
        "p3_init();"
    ),
    "shoc" : (
        "components/cam/src/physics/cam/shoc.F90",
        "components/scream/src/physics/shoc",
        "shoc_init(REPLACE_ME, true);"
    ),
}

# A good set of arg data for unit testing
UT_ARG_DATA = [
    ("foo1", "real", "in", ("shcol",)),
    ("foo2", "real", "in", ("shcol",)),
    ("bar1", "real", "in", ("shcol","nlev")),
    ("bar2", "real", "in", ("shcol","nlev")),
    ("bak1", "real", "in", ("shcol","nlevi")),
    ("bak2", "real", "in", ("shcol","nlevi")),
    ("gag", "real", "in", None),
    ("baz", "real", "inout", ("shcol",)),
    ("bag", "integer", "in", ("shcol",)),
    ("bab1", "integer", "out", None),
    ("bab2", "integer", "out", None),
    ("val", "logical", "in", None),
    ("shcol", "integer", "in", None),
    ("nlev", "integer", "in", None),
    ("nlevi", "integer", "in", None),
    ("ball1", "integer", "out", ("shcol",)),
    ("ball2", "integer", "out", ("shcol",)),
]

UT_ARG_DATA2 = [
    ("foo1", "real", "in", ("shcol",)),
    ("foo2", "real", "in", ("shcol",)),
    ("bar1", "real", "in", ("shcol","nlev")),
    ("bar2", "real", "in", ("shcol","nlev")),
    ("bak1", "real", "in", ("shcol","nlevi")),
    ("bak2", "real", "in", ("shcol","nlevi")),
    ("tracerd1", "real", "in", ("shcol","nlev", "ntracers")),
    ("tracerd2", "real", "in", ("shcol","nlev", "ntracers")),
    ("gag", "real", "in", None),
    ("baz", "real", "inout", ("shcol",)),
    ("bag", "integer", "in", ("shcol",)),
    ("bab1", "integer", "out", None),
    ("bab2", "integer", "out", None),
    ("val", "logical", "in", None),
    ("shcol", "integer", "in", None),
    ("nlev", "integer", "in", None),
    ("nlevi", "integer", "in", None),
    ("ntracers", "integer", "in", None),
    ("ball1", "integer", "out", ("shcol",)),
    ("ball2", "integer", "out", ("shcol",)),
]

#
# Free functions
#

###############################################################################
def get_piece_description(piece):
###############################################################################
    return PIECES[piece][DESC]()

###############################################################################
def get_subroutine_begin_regex(name):
###############################################################################
    """
    >>> bool(get_subroutine_begin_regex("fake_sub").match("subroutine fake_sub("))
    True
    >>> bool(get_subroutine_begin_regex("fake_sub").match("subroutine fake_sub_2("))
    False
    >>> bool(get_subroutine_begin_regex("fake_sub").match("  subroutine  fake_sub ("))
    True
    >>> bool(get_subroutine_begin_regex("fake_sub").match("  subroutine  fake_sub ( one, two )"))
    True
    >>> bool(get_subroutine_begin_regex("fake_sub").match("! subroutine fake_sub("))
    False
    >>> bool(get_subroutine_begin_regex("fake_sub").match("subroutine fake_sub"))
    False
    """
    subroutine_begin_regex_str = r"^\s*subroutine\s+{}\s*[(]".format(name)
    return re.compile(subroutine_begin_regex_str)

###############################################################################
def get_subroutine_end_regex(name):
###############################################################################
    """
    >>> bool(get_subroutine_end_regex("fake_sub").match("end subroutine fake_sub"))
    True
    >>> bool(get_subroutine_end_regex("fake_sub").match("end subroutine fake_sub_2"))
    False
    >>> bool(get_subroutine_end_regex("fake_sub").match("  end  subroutine  fake_sub "))
    True
    >>> bool(get_subroutine_end_regex("fake_sub").match("!end  subroutine  fake_sub "))
    False
    """
    subroutine_end_regex_str = r"^\s*end\s+subroutine\s+{}\s*$".format(name)
    return re.compile(subroutine_end_regex_str)

###############################################################################
def get_cxx_function_begin_regex(name, static=False, template=None):
###############################################################################
    """
    >>> bool(get_cxx_function_begin_regex("fake_sub").match("void fake_sub("))
    True
    >>> bool(get_cxx_function_begin_regex("fake_sub").match("void fake_sub_2("))
    False
    >>> bool(get_cxx_function_begin_regex("fake_sub").match("  void  fake_sub ("))
    True
    >>> bool(get_cxx_function_begin_regex("fake_sub").match("  void  fake_sub ( one, two )"))
    True
    >>> bool(get_cxx_function_begin_regex("fake_sub").match("// void fake_sub("))
    False
    >>> bool(get_cxx_function_begin_regex("fake_sub").match("void fake_sub"))
    False
    >>> bool(get_cxx_function_begin_regex("fake_sub", static=True).match("static void fake_sub("))
    True
    >>> bool(get_cxx_function_begin_regex("fake_sub", template="Foo<T>").match("void Foo<T>::fake_sub("))
    True
    """
    static_regex_str = r"static\s+" if static else ""
    template_regex_str = r"{}::".format(template) if template else ""
    function_begin_regex_str = r"^\s*{}void\s+{}{}\s*[(]".format(static_regex_str, template_regex_str, name)
    return re.compile(function_begin_regex_str)

###############################################################################
def get_cxx_close_block_regex(semicolon=False, comment=None, at_line_start=False):
###############################################################################
    """
    >>> bool(get_cxx_close_block_regex().match("}"))
    True
    >>> bool(get_cxx_close_block_regex(at_line_start=True).match("}"))
    True
    >>> bool(get_cxx_close_block_regex().match(" } "))
    True
    >>> bool(get_cxx_close_block_regex(at_line_start=True).match(" }; "))
    False
    >>> bool(get_cxx_close_block_regex(at_line_start=True).match(" } "))
    False
    >>> bool(get_cxx_close_block_regex(comment="hi").match(" } "))
    False
    >>> bool(get_cxx_close_block_regex(comment="hi").match(" } // hi"))
    True
    >>> bool(get_cxx_close_block_regex(comment="hi").match("} // hi  "))
    True
    >>> bool(get_cxx_close_block_regex(semicolon=True).match(" } "))
    False
    >>> bool(get_cxx_close_block_regex(semicolon=True).match(" } ; "))
    True
    >>> bool(get_cxx_close_block_regex(semicolon=True).match("};"))
    True
    >>> bool(get_cxx_close_block_regex(semicolon=True, comment="hi", at_line_start=True).match("};  // hi"))
    True
    >>> bool(get_cxx_close_block_regex(semicolon=True, comment="hi", at_line_start=True).match("};  // hi there"))
    False
    """
    semicolon_regex_str = r"\s*;" if semicolon else ""
    line_start_regex_str = "" if at_line_start else r"\s*"
    comment_regex_str   = r"\s*//\s*{}".format(comment) if comment else ""
    close_block_regex_str = re.compile(r"^{}}}{}{}\s*$".format(line_start_regex_str, semicolon_regex_str, comment_regex_str))
    return re.compile(close_block_regex_str)

###############################################################################
def get_namespace_close_regex(namespace):
###############################################################################
    """
    >>> bool(get_namespace_close_regex("foo").match(" } // namespace foo"))
    True
    >>> bool(get_namespace_close_regex("foo").match(" } // namespace foo_bar"))
    False
    """
    return get_cxx_close_block_regex(comment=r"namespace\s+{}".format(namespace))

###############################################################################
def get_cxx_struct_begin_regex(struct):
###############################################################################
    """
    >>> bool(get_cxx_struct_begin_regex("Foo").match("struct Foo {"))
    True
    >>> bool(get_cxx_struct_begin_regex("Foo").match("struct Foo"))
    True
    >>> bool(get_cxx_struct_begin_regex("Foo").match("struct FooBar"))
    False
    """
    struct_regex_str = r"^\s*struct\s+{}([\W]|$)".format(struct)
    return re.compile(struct_regex_str)

###############################################################################
def get_data_struct_name(sub):
###############################################################################
    """
    >>> get_data_struct_name("my_sub_name")
    'MySubNameData'
    >>> get_data_struct_name("sub")
    'SubData'
    """
    return "".join([item.capitalize() for item in sub.split("_")]) + "Data"

###############################################################################
def get_data_test_struct_name(sub):
###############################################################################
    """
    >>> get_data_test_struct_name("my_sub_name")
    'TestMySubName'
    >>> get_data_test_struct_name("update_prognostics_implicit")
    'TestUpdatePrognosticsImplicit'
    """
    return "Test{}".format(get_data_struct_name(sub)[:-4])

###############################################################################
def get_supported_pieces():
###############################################################################
    return PIECES.keys()

###############################################################################
def get_supported_physics():
###############################################################################
    return PHYSICS.keys()

###############################################################################
def get_piece_data(physics, sub, piece_name, piece_data, gb):
###############################################################################
    return PIECES[piece_name][piece_data](physics, sub, gb)

###############################################################################
def get_physics_data(physics_name, physics_data):
###############################################################################
    return PHYSICS[physics_name][physics_data]

###############################################################################
def expect_exists(physics, sub, gb, piece):
###############################################################################
    filepath = gb.get_path_for_piece_file(physics, sub, piece)
    expect(filepath.exists(), "For generating {}'s {} for phyiscs {}, expected file {} to already exist".\
           format(sub, piece, physics, filepath))
    return False # File was not created

###############################################################################
def create_template(physics, sub, gb, piece, force=False, force_arg_data=None):
###############################################################################
    """
    Create a file based on a template if it doesn't exist. Return True if a file was created.

    >>> gb = GenBoiler(["linear_interp"], ["cxx_func_impl"], dry_run=True)
    >>> create_template("shoc", "linear_interp", gb, "cxx_func_impl", force=True, force_arg_data=UT_ARG_DATA) #doctest: +ELLIPSIS
    Would create file .../components/scream/src/physics/shoc/shoc_linear_interp_impl.hpp with contents:
    #ifndef SHOC_LINEAR_INTERP_IMPL_HPP
    #define SHOC_LINEAR_INTERP_IMPL_HPP
    <BLANKLINE>
    #include "shoc_functions.hpp" // for ETI only but harmless for GPU
    <BLANKLINE>
    namespace scream {
    namespace shoc {
    <BLANKLINE>
    /*
     * Implementation of shoc linear_interp. Clients should NOT
     * #include this file, but include shoc_functions.hpp instead.
     */
    <BLANKLINE>
    template<typename S, typename D>
    KOKKOS_FUNCTION
    void Functions<S,D>::linear_interp(const uview_1d<const Spack>& foo1, const uview_1d<const Spack>& foo2, const uview_1d<const Spack>& bar1, const uview_1d<const Spack>& bar2, const uview_1d<const Spack>& bak1, const uview_1d<const Spack>& bak2, const Spack& gag, const uview_1d<Spack>& baz, const uview_1d<const Int>& bag, Int& bab1, Int& bab2, const bool& val, const Int& shcol, const Int& nlev, const Int& nlevi, const uview_1d<Int>& ball1, const uview_1d<Int>& ball2)
    {
      // TODO
      // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
    }
    <BLANKLINE>
    } // namespace shoc
    } // namespace scream
    <BLANKLINE>
    #endif
    <BLANKLINE>
    True
    """
    filepath = gb.get_path_for_piece_file(physics, sub, piece)
    if not filepath.exists() or force:
        expect(piece in FILE_TEMPLATES,
               "{} does not exist and there is no template for generating files for piece {}".format(filepath, piece))

        gen_code = getattr(gb, "gen_{}".format(piece))(physics, sub, force_arg_data=force_arg_data)
        contents = FILE_TEMPLATES[piece](physics, sub, gen_code)
        if gb.dry_run():
            print("Would create file {} with contents:\n{}".format(filepath, contents))
        else:
            with filepath.open("w") as fd:
                fd.write(contents)

            print("SUCCESS (Generated from scratch)")

        return True
    else:
        return False

###############################################################################
def remove_comments_and_ws(contents):
###############################################################################
    """
    >>> teststr = '''
    ... module mymod
    ...   subroutine foo(a, b, &{0}
    ...                c, d, e,&
    ... !bad{0}
    ... &f)
    ...
    ...     real, intent(in) :: a, b, & !go
    ...                 c, d, e, f{0}
    ...
    ...   ! hi
    ...   !hi ! there{0}
    ... !hi ! there
    ...   end subroutine foo{0}
    ... end module mymod
    ... '''.format(" ")
    >>> print(remove_comments_and_ws(teststr))
    module mymod
    subroutine foo(a, b, &
    c, d, e,&
    &f)
    real, intent(in) :: a, b, &
    c, d, e, f
    end subroutine foo
    end module mymod
    """
    new_lines = []
    comment_regex = re.compile(r"^([^!]*)")
    for line in contents.splitlines():
        m = comment_regex.match(line)
        if m is not None:
            line = m.groups()[0].strip()
        else:
            line = line.strip()

        if line != "":
            new_lines.append(line)

    return "\n".join(new_lines)

###############################################################################
def resolve_line_continuations(contents):
###############################################################################
    """
    >>> teststr = '''
    ... module mymod
    ...   subroutine foo(a, b, &
    ...                c, d, e,&{0}
    ... !bad
    ... &f){0}
    ...
    ...     real, intent(in) :: a, b, & !go{0}
    ...                 c, d, e, f
    ...
    ...   ! hi
    ...   !hi ! there
    ... !hi ! there{0}
    ...   end subroutine foo
    ... end module mymod{0}
    ... '''.format("  ")
    >>> print(resolve_line_continuations(teststr))
    module mymod
    subroutine foo(a, b, c, d, e,f)
    real, intent(in) :: a, b, c, d, e, f
    end subroutine foo
    end module mymod
    """
    # Must remove comments and whitespace for the alg below to work
    contents = remove_comments_and_ws(contents)

    new_lines = []
    contination = False
    for line in contents.splitlines():
        line = line.lstrip("&") # always remove leading &, they are optional
        if line != "":
            if contination:
                new_lines[-1] += line.rstrip("&")
            else:
                new_lines.append(line.rstrip("&"))

            contination = line.endswith("&")

    return "\n".join(new_lines)

###############################################################################
def normalize_f90(contents):
###############################################################################
    # We do not currently attempt to preprocess the contents
    return resolve_line_continuations(contents).lower()

###############################################################################
def split_top_commas(line):
###############################################################################
    """
    Split a line along commas except for commas within parentheses
    """
    raw_splits = [item.strip() for item in line.strip().split(",")]
    top_splits = []
    balanced = True
    for raw_split in raw_splits:
        if balanced:
            top_splits.append(raw_split)
        else:
            top_splits[-1] += ",{}".format(raw_split)

        balanced = top_splits[-1].count("(") == top_splits[-1].count(")")

    return top_splits

###############################################################################
def get_arg_order(line):
###############################################################################
    """
    Given a line of fortran declaring a subroutine, return the arg names in order

    >>> get_arg_order("subroutine p3_set_tables( mu_r_user, revap_user,vn_user, vm_user )")
    ['mu_r_user', 'revap_user', 'vn_user', 'vm_user']
    """
    args_raw = line.rstrip(")").split("(", maxsplit=1)[-1]
    return [item.strip() for item in args_raw.split(",") if item]

ARG_NAME, ARG_TYPE, ARG_INTENT, ARG_DIMS = range(4)
###############################################################################
def parse_f90_args(line):
###############################################################################
    """
    Given a line of fortran code declaring an argument[s], return [(argname, argtype, intent, dims)]
    Anywhere you see "arg_data" in this program, it's referring to a list of data produced by
    this function. Anywhere you see "arg_datum", it's referring to a single item in this list.

    >>> parse_f90_args('integer, intent(in) :: kts, kte, kbot')
    [('kts', 'integer', 'in', None), ('kte', 'integer', 'in', None), ('kbot', 'integer', 'in', None)]
    >>> parse_f90_args('real(rtype),intent(inout ), dimension(kts:kte) :: pres,dpres,  dz ')
    [('pres', 'real', 'inout', ('kts:kte',)), ('dpres', 'real', 'inout', ('kts:kte',)), ('dz', 'real', 'inout', ('kts:kte',))]
    >>> parse_f90_args('logical (btype), intent( in) ::do_predict_nc')
    [('do_predict_nc', 'logical', 'in', None)]
    >>> parse_f90_args('real(rtype),intent(inout), dimension( kts:kte, its: ite) :: dz')
    [('dz', 'real', 'inout', ('kts:kte', 'its:ite'))]
    >>> parse_f90_args('real(rtype),intent(inout), dimension(3) :: dz')
    [('dz', 'real', 'inout', ('3',))]
    >>> parse_f90_args('real(rtype),intent(inout), dimension(3,4) :: dz')
    [('dz', 'real', 'inout', ('3', '4'))]
    >>> parse_f90_args('real(rtype), dimension(3,4),intent(inout) :: dz')
    [('dz', 'real', 'inout', ('3', '4'))]
    >>> parse_f90_args('real(rtype), intent(in) :: x1(ncol,km1), y1(ncol , km1 )')
    [('x1', 'real', 'in', ('ncol', 'km1')), ('y1', 'real', 'in', ('ncol', 'km1'))]
    >>> parse_f90_args('real(rtype), intent(in) :: x1(ncol,km1,ntracers)')
    [('x1', 'real', 'in', ('ncol', 'km1', 'ntracers'))]
    """
    expect(line.count("::") == 1, "Expected line format 'type-info :: names' for: {}".format(line))
    metadata_str, names_str = line.split("::")
    names_dims = split_top_commas(names_str)
    metadata   = split_top_commas(metadata_str)

    argtype = metadata[0].split("(")[0].strip()
    intent, dims = None, None
    for metadatum in metadata:
        if metadatum.startswith("intent"):
            expect(intent is None, "Multiple intents in line: {}".format(line))
            intent = metadatum.split("(")[-1].rstrip(")").strip()
        elif metadatum.startswith("dimension"):
            expect(dims is None, "Multiple dimensions in line: {}".format(line))
            dims_raw = metadatum.split("(")[-1].rstrip(")").strip()
            dims = tuple(item.replace(" ", "") for item in dims_raw.split(","))

    names = []
    for name_dim in names_dims:
        if "(" in name_dim:
            name, dims_raw = name_dim.split("(")
            dims_raw = dims_raw.rstrip(")").strip()
            dims_check = tuple(item.replace(" ", "") for item in dims_raw.split(","))
            expect(dims is None or dims_check == dims, "Inconsistent dimensions in line: {}".format(line))
            dims = dims_check
            names.append(name.strip())
        else:
            names.append(name_dim.strip())

    return [(name, argtype, intent, dims) for name in names]

###############################################################################
def parse_origin(contents, subs):
###############################################################################
    r"""
    Returns a map of subname->[(argname, argtype, intent, dims)]

    >>> teststr = '''
    ...
    ...   SUBROUTINE p3_get_tables(mu_r_user, revap_user, &
    ...           tracerd, vn_user, vm_user)
    ...     ! This can be called after p3_init_b.
    ...     implicit none
    ...     real(rtype), dimension(150), intent(out) :: mu_r_user
    ...     real(rtype), dimension(300,10), intent(out) :: vn_user, vm_user, revap_user
    ...     real(rtype), dimension(300,10,42), intent(out) :: tracerd
    ...     mu_r_user(:) = mu_r_table(:)
    ...     revap_user(:,:) = revap_table(:,:)
    ...     vn_user(:,:) = vn_table(:,:)
    ...     vm_user(:,:) = vm_table(:,:)
    ...
    ...    return
    ...
    ...   end SUBROUTINE p3_get_tables
    ...
    ...   subroutine p3_set_tables( mu_r_user, revap_user,vn_user, vm_user )
    ...     ! This can be called instead of p3_init_b.
    ...     implicit none
    ...     real(rtype), dimension(300,10), intent(in) :: vn_user, vm_user, revap_user
    ...     real(rtype), dimension(150), intent(in) :: mu_r_user
    ...     mu_r_table(:) = mu_r_user(:)
    ...     revap_table(:,:) = revap_user(:,:)
    ...     vn_table(:,:) = vn_user(:,:)
    ...     vm_table(:,:) = vm_user(:,:)
    ...
    ...    return
    ...
    ...   END SUBROUTINE p3_set_tables
    ...
    ...   SUBROUTINE p3_init_b()
    ...     implicit none
    ...     integer                      :: i,ii,jj,kk
    ...     real(rtype)                         :: lamr,mu_r,dm,dum1,dum2,dum3,dum4,dum5,  &
    ...          dd,amg,vt,dia
    ...
    ...     ! AaronDonahue: Switching to table ver 4 means switching to a constand mu_r,
    ...     ! so this section is commented out.
    ...     do i = 1,150
    ...   END SUBROUTINE p3_init_b
    ... '''
    >>> print("\n".join([str(item) for item in sorted(parse_origin(teststr, ["p3_get_tables", "p3_init_b"]).items())]))
    ('p3_get_tables', [('mu_r_user', 'real', 'out', ('150',)), ('revap_user', 'real', 'out', ('300', '10')), ('tracerd', 'real', 'out', ('300', '10', '42')), ('vn_user', 'real', 'out', ('300', '10')), ('vm_user', 'real', 'out', ('300', '10'))])
    ('p3_init_b', [])
    """
    begin_regexes = [get_subroutine_begin_regex(sub) for sub in subs]
    arg_decl_regex = re.compile(r"^.+intent\s*[(]\s*(in|out|inout)\s*[)]")

    contents = normalize_f90(contents)

    db = {}
    active_sub = None
    arg_order = []
    arg_decls = []
    for line in contents.splitlines():
        begin_match = None
        for sub, begin_regex in zip(subs, begin_regexes):
            begin_match = begin_regex.match(line)
            if begin_match is not None:
                expect(active_sub is None, "subroutine {} was still active when {} began".format(active_sub, sub))
                active_sub = sub
                arg_order = get_arg_order(line)

        if active_sub:
            decl_match = arg_decl_regex.match(line)
            if decl_match is not None:
                arg_decls.extend(parse_f90_args(line))

            end_regex = get_subroutine_end_regex(active_sub)
            end_match = end_regex.match(line)
            if end_match is not None:
                expect(active_sub not in db, "Found multiple matches for {}".format(active_sub))
                expect(len(arg_order) == len(arg_decls),
                       "Number of decls:\n{}\nDid not match arg list: {}".format(arg_decls, arg_order))

                # we need our decls to be ordered based on arg list order
                ordered_decls = []
                for arg in arg_order:
                    found = False
                    for arg_decl in arg_decls:
                        if arg_decl[ARG_NAME] == arg:
                            ordered_decls.append(arg_decl)
                            found = True
                            break

                    expect(found, "Could not find decl for arg {} in\n{}".format(arg, arg_decls))

                db[active_sub] = ordered_decls
                active_sub = None
                arg_decls = []

    return db

C_TYPE_MAP = {"real" : "c_real", "integer" : "c_int", "logical" : "c_bool"}
###############################################################################
def gen_arg_f90_decl(argtype, intent, dims, names):
###############################################################################
    """
    Generate f90 argument declaration based on the input data

    >>> gen_arg_f90_decl("real", "in", ("10", "150"), ["foo", "bar"])
    'real(kind=c_real) , intent(in), dimension(10, 150) :: foo, bar'
    >>> gen_arg_f90_decl("real", "out", ("10", "150"), ["foo", "bar"])
    'real(kind=c_real) , intent(out), dimension(10, 150) :: foo, bar'
    >>> gen_arg_f90_decl("real", "out", ("10", "150", "42"), ["foo", "bar"])
    'real(kind=c_real) , intent(out), dimension(10, 150, 42) :: foo, bar'
    >>> gen_arg_f90_decl("logical", "in", None, ["biz", "baz"])
    'logical(kind=c_bool) , value, intent(in) :: biz, baz'
    >>> gen_arg_f90_decl("integer", "inout", None, ["barg"])
    'integer(kind=c_int) , intent(inout) :: barg'
    >>> gen_arg_f90_decl("integer", "out", None, ["barg"])
    'integer(kind=c_int) , intent(out) :: barg'
    """
    expect(argtype in C_TYPE_MAP, "Unrecognized argtype for C_TYPE_MAP: {}".format(argtype))
    c_type = C_TYPE_MAP[argtype]
    value  = ", value" if dims is None and intent == "in" else ""
    intent_s = ", intent({})".format(intent)
    dimension_s = ", dimension({})".format(", ".join(dims)) if dims is not None else ""
    names_s = ", ".join(names)
    return "{argtype}(kind={c_type}) {value}{intent}{dimension} :: {names}".\
        format(argtype=argtype, c_type=c_type, value=value, intent=intent_s, dimension=dimension_s, names=names_s)

CXX_TYPE_MAP = {"real" : "Real", "integer" : "Int", "logical" : "bool"}
###############################################################################
def get_cxx_type(arg_datum):
###############################################################################
    """
    Based on arg datum, give c++ bridge type

    >>> get_cxx_type(("foo", "real", "in", ("100",)))
    'Real*'
    >>> get_cxx_type(("foo", "real", "in", None))
    'Real'
    >>> get_cxx_type(("foo", "real", "inout", ("100",)))
    'Real*'
    >>> get_cxx_type(("foo", "real", "inout", None))
    'Real*'
    >>> get_cxx_type(("foo", "real", "out", ("100",)))
    'Real*'
    >>> get_cxx_type(("foo", "real", "out", None))
    'Real*'
    >>> get_cxx_type(("foo", "integer", "inout", None))
    'Int*'
    """
    is_ptr = arg_datum[ARG_DIMS] is not None or arg_datum[ARG_INTENT] != "in"
    arg_type = arg_datum[ARG_TYPE]
    expect(arg_type in CXX_TYPE_MAP, "Unrecognized argtype for CXX_TYPE_MAP: {}".format(arg_type))
    arg_cxx_type = CXX_TYPE_MAP[arg_type]
    return "{}{}".format(arg_cxx_type, "*" if is_ptr else "")

KOKKOS_TYPE_MAP = {"real" : "Spack", "integer" : "Int", "logical" : "bool"}
###############################################################################
def get_kokkos_type(arg_datum):
###############################################################################
    """
    Based on arg datum, give c++ kokkos type

    Note: We can only guess at the correct types, especially whether an argument
    should be packed data or not!

    >>> get_kokkos_type(("foo", "real", "in", ("100",)))
    'const uview_1d<const Spack>&'
    >>> get_kokkos_type(("foo", "real", "in", None))
    'const Spack&'
    >>> get_kokkos_type(("foo", "real", "inout", ("100",)))
    'const uview_1d<Spack>&'
    >>> get_kokkos_type(("foo", "real", "inout", None))
    'Spack&'
    >>> get_kokkos_type(("foo", "real", "out", ("100",)))
    'const uview_1d<Spack>&'
    >>> get_kokkos_type(("foo", "real", "out", None))
    'Spack&'
    >>> get_kokkos_type(("foo", "integer", "inout", None))
    'Int&'
    """
    is_const  = arg_datum[ARG_INTENT] == "in"
    is_view   = arg_datum[ARG_DIMS] is not None
    base_type = "{}{}".format("const " if is_const else "", KOKKOS_TYPE_MAP[arg_datum[ARG_TYPE]])

    # We assume 1d even if the f90 array is 2d since we assume c++ will spawn a kernel
    # over one of the dimensions
    return "const uview_1d<{}>&".format(base_type) if is_view else "{}&".format(base_type)

###############################################################################
def gen_arg_cxx_decls(arg_data, kokkos=False):
###############################################################################
    """
    Get all arg decls for a set of arg data

    >>> gen_arg_cxx_decls([("foo", "real", "in", ("100",)), ("bar", "real", "in", None)])
    ['Real* foo', 'Real bar']
    >>> gen_arg_cxx_decls([("foo", "real", "in", ("100",)), ("bar", "real", "in", None)], kokkos=True)
    ['const uview_1d<const Spack>& foo', 'const Spack& bar']
    """
    arg_names    = [item[ARG_NAME] for item in arg_data]
    get_type     = get_kokkos_type if kokkos else get_cxx_type
    arg_types    = [get_type(item) for item in arg_data]
    arg_sig_list = ["{} {}".format(arg_type, arg_name) for arg_name, arg_type in zip(arg_names, arg_types)]
    return arg_sig_list

###############################################################################
def needs_transpose(arg_data):
###############################################################################
    """
    Based on data, does this sub need to have its data transposed when going between languages?

    >>> needs_transpose([("foo", "real", "in", ("100",)), ("bar", "real", "in", None)])
    False
    >>> needs_transpose([("foo", "real", "in", ("100","50")), ("bar", "real", "in", None)])
    True
    >>> needs_transpose([("foo", "real", "in", ("100","50","42")), ("bar", "real", "in", None)])
    True
    >>> needs_transpose([("foo", "real", "in", None), ("bar", "real", "in", None)])
    False
    """
    for arg_datum in arg_data:
        arg_dims = arg_datum[ARG_DIMS]
        if arg_dims is not None and len(arg_dims) > 1:
            return True

    return False

###############################################################################
def gen_cxx_data_args(physics, arg_data):
###############################################################################
    """
    Based on data, generate unpacking of Data struct args

    >>> gen_cxx_data_args("shoc", UT_ARG_DATA)
    ['d.foo1', 'd.foo2', 'd.bar1', 'd.bar2', 'd.bak1', 'd.bak2', 'd.gag', 'd.baz', 'd.bag', '&d.bab1', '&d.bab2', 'd.val', 'd.shcol()', 'd.nlev()', 'd.nlevi()', 'd.ball1', 'd.ball2']

    >>> gen_cxx_data_args("shoc", UT_ARG_DATA2)
    ['d.foo1', 'd.foo2', 'd.bar1', 'd.bar2', 'd.bak1', 'd.bak2', 'd.tracerd1', 'd.tracerd2', 'd.gag', 'd.baz', 'd.bag', '&d.bab1', '&d.bab2', 'd.val', 'd.shcol', 'd.nlev', 'd.nlevi', 'd.ntracers', 'd.ball1', 'd.ball2']
    """
    all_dims = group_data(arg_data)[3]
    func_call = "()" if is_sugar_compatible(physics, arg_data) else ""
    args_needs_ptr = [item[ARG_DIMS] is None and item[ARG_INTENT] != "in" for item in arg_data]
    arg_names      = [item[ARG_NAME] for item in arg_data]
    arg_dim_call   = [item[ARG_NAME] in all_dims for item in arg_data]
    args = ["{}d.{}{}".format("&" if need_ptr else "", arg_name, func_call if dim_call else "")
            for arg_name, need_ptr, dim_call in zip(arg_names, args_needs_ptr, arg_dim_call)]
    return args

###############################################################################
def gen_arg_f90_decls(arg_data):
###############################################################################
    r"""
    Generate f90 argument declarations, will attempt to group these together if possible.

    >>> print("\n".join(gen_arg_f90_decls(UT_ARG_DATA)))
    real(kind=c_real) , intent(in), dimension(shcol) :: foo1, foo2
    real(kind=c_real) , intent(in), dimension(shcol, nlev) :: bar1, bar2
    real(kind=c_real) , intent(in), dimension(shcol, nlevi) :: bak1, bak2
    real(kind=c_real) , value, intent(in) :: gag
    real(kind=c_real) , intent(inout), dimension(shcol) :: baz
    integer(kind=c_int) , intent(in), dimension(shcol) :: bag
    integer(kind=c_int) , intent(out) :: bab1, bab2
    logical(kind=c_bool) , value, intent(in) :: val
    integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi
    integer(kind=c_int) , intent(out), dimension(shcol) :: ball1, ball2
    """
    metadata = OrderedDict()
    for name, argtype, intent, dims in arg_data:
        metatuple = (argtype, intent, dims)
        metadata.setdefault(metatuple, []).append(name)

    result = []
    for metatuple, names in metadata.items():
        result.append(gen_arg_f90_decl(*metatuple, names))

    return result

###############################################################################
def has_arrays(arg_data):
###############################################################################
    for _, _, _, dims in arg_data:
        if dims is not None:
            return True

    return False

###############################################################################
def gen_struct_members(arg_data):
###############################################################################
    r"""
    Gen cxx code for data struct members

    >>> print("\n".join(gen_struct_members(UT_ARG_DATA)))
    // Inputs
    Real *foo1, *foo2, *bar1, *bar2, *bak1, *bak2;
    Real gag;
    Int *bag;
    bool val;
    Int shcol, nlev, nlevi;
    <BLANKLINE>
    // Inputs/Outputs
    Real *baz;
    <BLANKLINE>
    // Outputs
    Int bab1, bab2;
    Int *ball1, *ball2;
    <BLANKLINE>
    """
    metadata = {} # intent -> (type, is_ptr) -> names
    for name, argtype, intent, dims in arg_data:
        metadata.setdefault(intent, OrderedDict()).setdefault((argtype, dims is not None), []).append(name)

    intent_order = ( ("in", "Inputs"), ("inout", "Inputs/Outputs"), ("out", "Outputs") )
    result = []
    for intent, comment in intent_order:
        if intent in metadata:
            result.append("// {}".format(comment))
            type_map = metadata[intent]
            for type_info, names in type_map.items():
                type_name, is_ptr = type_info
                decl_str = CXX_TYPE_MAP[type_name]
                decl_str += " {};".format(", ".join(["{}{}".format("*" if is_ptr else "", name) for name in names]))
                result.append(decl_str)

            result.append("")

    return result

###############################################################################
def group_data(arg_data, filter_out_intent=None):
###############################################################################
    r"""
    Given data, return ([fst_dims], [snd_dims], [trd_dims], [all-dims], [scalars], {dims->[real_data]}, {dims->[int_data]})

    >>> print("\n".join([str(item) for item in group_data(UT_ARG_DATA2)]))
    ['shcol']
    ['nlev', 'nlevi']
    ['ntracers']
    ['shcol', 'nlev', 'nlevi', 'ntracers']
    [('gag', 'Real'), ('bab1', 'Int'), ('bab2', 'Int'), ('val', 'bool')]
    OrderedDict([(('shcol',), ['foo1', 'foo2', 'baz']), (('shcol', 'nlev'), ['bar1', 'bar2']), (('shcol', 'nlevi'), ['bak1', 'bak2']), (('shcol', 'nlev', 'ntracers'), ['tracerd1', 'tracerd2'])])
    OrderedDict([(('shcol',), ['bag', 'ball1', 'ball2'])])
    """
    scalars  = []

    fst_dims = []
    snd_dims = []
    trd_dims = []

    for name, argtype, _, dims in arg_data:
        if dims is not None:
            expect(len(dims) >= 1 and len(dims) <= 3,
                   "Only 1d-3d data is supported, {} has too many dims: {}".format(name, len(dims)))

            if dims[0] not in fst_dims:
                fst_dims.append(dims[0])
            if len(dims) > 1 and dims[1] not in snd_dims:
                snd_dims.append(dims[1])
            if len(dims) > 2 and dims[2] not in trd_dims:
                trd_dims.append(dims[2])

    all_dims = list(OrderedDict([(item, None) for item in (fst_dims + snd_dims + trd_dims)]))
    real_data = OrderedDict()
    int_data = OrderedDict()

    for name, argtype, intent, dims in arg_data:
        if filter_out_intent is None or intent != filter_out_intent:
            if dims is None:
                if name not in all_dims:
                    scalars.append( (name, CXX_TYPE_MAP[argtype]))
                else:
                    expect(argtype == "integer", "Expected dimension {} to be of type integer".format(name))

            elif argtype == "integer":
                int_data.setdefault(dims, []).append(name)

            else:
                real_data.setdefault(dims, []).append(name)

    return fst_dims, snd_dims, trd_dims, all_dims, scalars, real_data, int_data

###############################################################################
def is_sugar_compatible(physics, arg_data):
###############################################################################
    """
    Can we use PhysicsTestData syntax sugar

    >>> is_sugar_compatible("shoc", UT_ARG_DATA)
    True
    >>> is_sugar_compatible("shoc", UT_ARG_DATA2)
    False
    >>> is_sugar_compatible("p3", UT_ARG_DATA)
    False
    """
    fst_dims, snd_dims, trd_dims, all_dims, _, _, int_data = group_data(arg_data)
    sugar_compatible = physics == "shoc"
    if len(all_dims) > 3 or len(fst_dims) > 1 or len(snd_dims) > 2 or len(trd_dims) > 1:
        sugar_compatible = False
    if int_data and (len(int_data) > 1 or (fst_dims[0],) not in int_data):
        sugar_compatible = False

    return sugar_compatible

###############################################################################
def gen_struct_api_sugar(physics, struct_name, arg_data):
###############################################################################
    """
    Return struct contents using syntactic sugar
    """
    _, _, _, all_dims, scalars, real_data, int_data = group_data(arg_data)
    ik_reals, ij_reals, i_reals, td_reals, i_ints = [], [], [], [], []
    for dims, reals in real_data.items():
        if len(dims) == 1:
            expect(not i_reals, "Multiple sets of i_reals?")
            i_reals = reals

        if len(dims) == 2:
            if not ik_reals:
                ik_reals = reals
            elif not ij_reals:
                ij_reals = reals
            else:
                expect(False, "Multiple sets of 2d reals?")

        if len(dims) == 3:
            expect(not td_reals, "Multiple sets of 3d reals?")
            td_reals = reals

    i_ints = list(int_data.values())[0]

    result = []
    dim_args = [(item, "Int") for item in all_dims if item is not None]
    cons_args = dim_args + scalars
    result.append("{struct_name}({cons_args}) :".\
                  format(struct_name=struct_name,
                         cons_args=", ".join(["{} {}_".format(argtype, name) for name, argtype in cons_args])))
    parent_call = "  PhysicsTestData({}".format(", ".join(["{}_".format(name) for name, _ in dim_args]))
    for item in (ik_reals, ij_reals, i_reals, i_ints, td_reals):
        if len(item) > 0:
            parent_call += ", {{{}}}".format(", ".join(["&{}".format(name) for name in item]))

    parent_call += ")"
    if scalars:
        parent_call += ", {}".format(", ".join(["{0}({0}_)".format(name) for name, _ in scalars]))

    parent_call += " {}"
    result.append(parent_call)
    result.append("")

    if physics != "p3":
        if len(scalars) == 0:
            result.append("SHOC_NO_SCALAR({}, {})".format(struct_name, len(dim_args)))
        else:
            result.append("SHOC_SCALARS({}, {}, {}, {})".format(struct_name, len(dim_args), len(scalars),
                                                                ", ".join([name for name, _ in scalars])))
    else:
        expect(False, "p3 sugar is not supported for now") # TODO

    return result

###############################################################################
def gen_struct_api_generic(physics, struct_name, arg_data):
###############################################################################
    _, _, _, all_dims, scalars, real_data, int_data = group_data(arg_data)

    result = []
    dim_args = [(item, "Int") for item in all_dims if item is not None]
    cons_args = dim_args + scalars
    result.append("{struct_name}({cons_args}) :".\
                  format(struct_name=struct_name,
                         cons_args=", ".join(["{} {}_".format(argtype, name) for name, argtype in cons_args])))

    dim_cxx_vec = []
    real_vec = []
    int_vec = []
    for dims, reals in real_data.items():
        dim_cxx_vec.append("{{ {} }}".format(", ".join(["{}_".format(item) for item in dims])))
        real_vec.append("{{ {} }}".format(", ".join(["&{}".format(item) for item in reals])))

    for dims, ints in int_data.items():
        dim_cxx_vec.append("{{ {} }}".format(", ".join(["{}_".format(item) for item in dims])))
        int_vec.append("{{ {} }}".format(", ".join(["&{}".format(item) for item in ints])))

    parent_call = "  PhysicsTestDataGeneric({{{}}}, {{{}}}, {{{}}})".format(", ".join(dim_cxx_vec), ", ".join(real_vec), ", ".join(int_vec))

    parent_call += ", {}".format(", ".join(["{0}({0}_)".format(name) for name, _ in cons_args]))

    parent_call += " {}"
    result.append(parent_call)
    result.append("")

    result.append("PTDG_STD_DEF({}, {}, {});".\
                  format(struct_name, len(cons_args), ", ".join([name for name, _ in cons_args])))

    return result

###############################################################################
def gen_struct_api(physics, struct_name, arg_data):
###############################################################################
    r"""
    Given data, generate code for data struct api

    >>> print("\n".join(gen_struct_api("shoc", "DataSubName", UT_ARG_DATA)))
    DataSubName(Int shcol_, Int nlev_, Int nlevi_, Real gag_, Int bab1_, Int bab2_, bool val_) :
      PhysicsTestData(shcol_, nlev_, nlevi_, {&bar1, &bar2}, {&bak1, &bak2}, {&foo1, &foo2, &baz}, {&bag, &ball1, &ball2}), gag(gag_), bab1(bab1_), bab2(bab2_), val(val_) {}
    <BLANKLINE>
    SHOC_SCALARS(DataSubName, 3, 4, gag, bab1, bab2, val)

    >>> print("\n".join(gen_struct_api("shoc", "DataSubName", UT_ARG_DATA2)))
    DataSubName(Int shcol_, Int nlev_, Int nlevi_, Int ntracers_, Real gag_, Int bab1_, Int bab2_, bool val_) :
      PhysicsTestDataGeneric({{ shcol_ }, { shcol_, nlev_ }, { shcol_, nlevi_ }, { shcol_, nlev_, ntracers_ }, { shcol_ }}, {{ &foo1, &foo2, &baz }, { &bar1, &bar2 }, { &bak1, &bak2 }, { &tracerd1, &tracerd2 }}, {{ &bag, &ball1, &ball2 }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_), ntracers(ntracers_), gag(gag_), bab1(bab1_), bab2(bab2_), val(val_) {}
    <BLANKLINE>
    PTDG_STD_DEF(DataSubName, 8, shcol, nlev, nlevi, ntracers, gag, bab1, bab2, val);
    """
    sugar_compatible = is_sugar_compatible(physics, arg_data)

    if sugar_compatible:
        return gen_struct_api_sugar(physics, struct_name, arg_data)
    else:
        return gen_struct_api_generic(physics, struct_name, arg_data)

###############################################################################
def find_insertion(lines, insert_regex):
###############################################################################
    """
    Find the index that matches insert_regex. If not found, return None

    >>> lines = ["foo", "bar", "baz", "bag"]
    >>> find_insertion(lines, re.compile("baz"))
    2
    >>> find_insertion(lines, re.compile("ball"))
    >>>
    """
    for idx, line in enumerate(lines):
        has_match = bool(insert_regex.match(line))
        if has_match:
            return idx

    return None

###############################################################################
def check_existing_piece(lines, begin_regex, end_regex):
###############################################################################
    """
    Check to see if an existing block/piece of code existing in a file, given its
    starting and ending regex. Returns None if not found

    >>> lines = ["foo", "bar", "baz", "bag"]
    >>> check_existing_piece(lines, re.compile("foo"), re.compile("bag"))
    (0, 4)
    >>> check_existing_piece(lines, re.compile("foo"), re.compile("foo"))
    (0, 1)
    >>> check_existing_piece(lines, re.compile("zxzc"), re.compile("foo"))
    """
    begin_idx = None
    end_idx   = None

    for idx, line in enumerate(lines):
        begin_match = bool(begin_regex.match(line))
        end_match   = bool(end_regex.match(line))

        if begin_match:
            expect(begin_idx is None,
                   "Found multiple begin matches for pattern '{}' before end pattern '{}' was found".\
                   format(begin_regex.pattern, end_regex.pattern))

            begin_idx = idx

        if end_match and begin_idx is not None:
            end_idx = idx
            break

    if begin_idx is not None:
        expect(end_idx is not None,
               "Found no ending match for begin pattern '{}' starting on line {} and searching end pattern '{}'".\
               format(begin_regex.pattern, begin_idx, end_regex.pattern))

    return None if begin_idx is None else (begin_idx, end_idx+1)

#
# Main classes
#

###############################################################################
class GenBoiler(object):
###############################################################################

    ###########################################################################
    def __init__(self,
                 subs        = None,
                 pieces      = get_supported_pieces(),
                 physics     = get_supported_physics(),
                 overwrite   = False,
                 kernel      = False,
                 source_repo = get_git_toplevel_dir(),
                 target_repo = get_git_toplevel_dir(),
                 dry_run     = False,
                 verbose     = False):
    ###########################################################################
        expect(source_repo is not None, "Must either run from a valid repo or provide a --source-repo")
        expect(target_repo is not None, "Must either run from a valid repo or provide a --target-repo")

        normalized_source_repo = get_git_toplevel_dir(repo=source_repo)
        expect(normalized_source_repo is not None, "source repo {} is not a valid repo".format(source_repo))
        source_repo = normalized_source_repo

        normalized_target_repo = get_git_toplevel_dir(repo=target_repo)
        expect(normalized_target_repo is not None, "target repo {} is not a valid repo".format(target_repo))
        target_repo = normalized_target_repo

        # configuration
        self._subs        = subs
        self._pieces      = pieces
        self._physics     = physics
        self._overwrite   = overwrite
        self._kernel      = kernel
        self._source_repo = pathlib.Path(source_repo).resolve()
        self._target_repo = pathlib.Path(target_repo).resolve()
        self._dry_run     = dry_run
        self._verbose     = verbose

        # internals
        self._db          = {}

        # TODO: support subroutine rename?
        # TODO: support smart line wrapping?
        # TODO: support generation in main physics file

        if not self._pieces:
            self._pieces = get_supported_pieces()

    ###########################################################################
    def _get_db(self, phys):
    ###########################################################################
        if phys in self._db:
            return self._db[phys]
        else:
            origin_file = self._source_repo / get_physics_data(phys, ORIGIN_FILE)
            expect(origin_file.exists(), "Missing origin file for physics {}: {}".format(phys, origin_file))
            db = parse_origin(origin_file.open().read(), self._subs)
            self._db[phys] = db
            if self._verbose:
                print("For physics {}, found:")
                for sub in self._subs:
                    if sub in db:
                        print("  For subroutine {}, found args:")
                        for name, argtype, intent, dims in db[sub]:
                            print("    name:{} type:{} intent:{} dims:({})".\
                                  format(name, argtype, intent, ",".join(dims) if dims else "scalar"))
            return db

    ###########################################################################
    def _get_arg_data(self, phys, sub):
    ###########################################################################
        phys_db = self._get_db(phys)
        expect(sub in phys_db, "No data for subroutine {} in physics {}".format(sub, phys))
        return phys_db[sub]

    ###########################################################################
    def dry_run(self):
    ###########################################################################
        return self._dry_run

    ###############################################################################
    def get_path_for_piece_file(self, physics, sub, piece):
    ###############################################################################
        root_dir = pathlib.Path(get_physics_data(physics, CXX_ROOT))
        filepath = self._target_repo / root_dir / get_piece_data(physics, sub, piece, FILEPATH, self)
        return filepath

    #
    # Individual piece generator methods
    #

    ###########################################################################
    def gen_f90_c2f_bind(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_f90_c2f_bind("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
          subroutine fake_sub_c(foo1, foo2, bar1, bar2, bak1, bak2, gag, baz, bag, bab1, bab2, val, shcol, nlev, nlevi, ball1, ball2) bind(C)
            use shoc, only : fake_sub
        <BLANKLINE>
            real(kind=c_real) , intent(in), dimension(shcol) :: foo1, foo2
            real(kind=c_real) , intent(in), dimension(shcol, nlev) :: bar1, bar2
            real(kind=c_real) , intent(in), dimension(shcol, nlevi) :: bak1, bak2
            real(kind=c_real) , value, intent(in) :: gag
            real(kind=c_real) , intent(inout), dimension(shcol) :: baz
            integer(kind=c_int) , intent(in), dimension(shcol) :: bag
            integer(kind=c_int) , intent(out) :: bab1, bab2
            logical(kind=c_bool) , value, intent(in) :: val
            integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi
            integer(kind=c_int) , intent(out), dimension(shcol) :: ball1, ball2
        <BLANKLINE>
            call fake_sub(foo1, foo2, bar1, bar2, bak1, bak2, gag, baz, bag, bab1, bab2, val, shcol, nlev, nlevi, ball1, ball2)
          end subroutine fake_sub_c
        """
        arg_data = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        arg_names = ", ".join([item[ARG_NAME] for item in arg_data])
        arg_decls = gen_arg_f90_decls(arg_data)
        phys_mod = "micro_p3" if phys == "p3" else phys
        result = \
"""  subroutine {sub}_c({arg_names}) bind(C)
    use {phys_mod}, only : {sub}

    {arg_decls}

    call {sub}({arg_names})
  end subroutine {sub}_c""".format(sub=sub, arg_names=arg_names, phys_mod=phys_mod, arg_decls="\n    ".join(arg_decls))

        return result

    ###########################################################################
    def gen_f90_f2c_bind(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_f90_f2c_bind("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
          subroutine fake_sub_f(foo1, foo2, bar1, bar2, bak1, bak2, gag, baz, bag, bab1, bab2, val, shcol, nlev, nlevi, ball1, ball2) bind(C)
            use iso_c_binding
        <BLANKLINE>
            real(kind=c_real) , intent(in), dimension(shcol) :: foo1, foo2
            real(kind=c_real) , intent(in), dimension(shcol, nlev) :: bar1, bar2
            real(kind=c_real) , intent(in), dimension(shcol, nlevi) :: bak1, bak2
            real(kind=c_real) , value, intent(in) :: gag
            real(kind=c_real) , intent(inout), dimension(shcol) :: baz
            integer(kind=c_int) , intent(in), dimension(shcol) :: bag
            integer(kind=c_int) , intent(out) :: bab1, bab2
            logical(kind=c_bool) , value, intent(in) :: val
            integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi
            integer(kind=c_int) , intent(out), dimension(shcol) :: ball1, ball2
          end subroutine fake_sub_f
        """
        arg_data = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        arg_names = ", ".join([item[ARG_NAME] for item in arg_data])
        arg_decls = gen_arg_f90_decls(arg_data)
        result = \
"""  subroutine {sub}_f({arg_names}) bind(C)
    use iso_c_binding

    {arg_decls}
  end subroutine {sub}_f""".format(sub=sub, arg_names=arg_names, arg_decls="\n    ".join(arg_decls))

        return result

    ###########################################################################
    def gen_cxx_c2f_bind_decl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_c2f_bind_decl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
        void fake_sub_c(Real* foo1, Real* foo2, Real* bar1, Real* bar2, Real* bak1, Real* bak2, Real gag, Real* baz, Int* bag, Int* bab1, Int* bab2, bool val, Int shcol, Int nlev, Int nlevi, Int* ball1, Int* ball2);
        <BLANKLINE>
        """
        arg_data = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        arg_decls = gen_arg_cxx_decls(arg_data)
        result = "void {sub}_c({arg_sig});\n".format(sub=sub, arg_sig=", ".join(arg_decls))
        return result

    ###########################################################################
    def gen_cxx_c2f_glue_decl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_c2f_glue_decl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
        void fake_sub(FakeSubData& d);
        """
        struct_name      = get_data_struct_name(sub)
        result = "void {sub}({struct_name}& d);".format(sub=sub, struct_name=struct_name)
        return result

    ###########################################################################
    def gen_cxx_c2f_glue_impl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_c2f_glue_impl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
        void fake_sub(FakeSubData& d)
        {
          shoc_init(d.nlev(), true);
          d.transpose<ekat::TransposeDirection::c2f>();
          fake_sub_c(d.foo1, d.foo2, d.bar1, d.bar2, d.bak1, d.bak2, d.gag, d.baz, d.bag, &d.bab1, &d.bab2, d.val, d.shcol(), d.nlev(), d.nlevi(), d.ball1, d.ball2);
          d.transpose<ekat::TransposeDirection::f2c>();
        }

        >>> print(gb.gen_cxx_c2f_glue_impl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA2))
        void fake_sub(FakeSubData& d)
        {
          shoc_init(d.nlev, true);
          d.transpose<ekat::TransposeDirection::c2f>();
          fake_sub_c(d.foo1, d.foo2, d.bar1, d.bar2, d.bak1, d.bak2, d.tracerd1, d.tracerd2, d.gag, d.baz, d.bag, &d.bab1, &d.bab2, d.val, d.shcol, d.nlev, d.nlevi, d.ntracers, d.ball1, d.ball2);
          d.transpose<ekat::TransposeDirection::f2c>();
        }
        """
        arg_data         = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        arg_data_args    = ", ".join(gen_cxx_data_args(phys, arg_data))
        need_transpose   = needs_transpose(arg_data)
        transpose_code_1 = "\n  d.transpose<ekat::TransposeDirection::c2f>();" if need_transpose else ""
        transpose_code_2 = "\n  d.transpose<ekat::TransposeDirection::f2c>();" if need_transpose else ""
        data_struct      = get_data_struct_name(sub)
        init_code        = get_physics_data(phys, INIT_CODE)
        func_call        = "()" if is_sugar_compatible(phys, arg_data) else ""
        init_code        = init_code.replace("REPLACE_ME", "d.nlev{}".format(func_call))

        result = \
"""void {sub}({data_struct}& d)
{{
  {init_code}{transpose_code_1}
  {sub}_c({arg_data_args});{transpose_code_2}
}}""".format(sub=sub, data_struct=data_struct, init_code=init_code, transpose_code_1=transpose_code_1, transpose_code_2=transpose_code_2, arg_data_args=arg_data_args)
        return result

    ###########################################################################
    def gen_cxx_c2f_data(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_c2f_data("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
        struct FakeSubData : public PhysicsTestData {
          // Inputs
          Real *foo1, *foo2, *bar1, *bar2, *bak1, *bak2;
          Real gag;
          Int *bag;
          bool val;
          Int shcol, nlev, nlevi;
        <BLANKLINE>
          // Inputs/Outputs
          Real *baz;
        <BLANKLINE>
          // Outputs
          Int bab1, bab2;
          Int *ball1, *ball2;
        <BLANKLINE>
          FakeSubData(Int shcol_, Int nlev_, Int nlevi_, Real gag_, Int bab1_, Int bab2_, bool val_) :
            PhysicsTestData(shcol_, nlev_, nlevi_, {&bar1, &bar2}, {&bak1, &bak2}, {&foo1, &foo2, &baz}, {&bag, &ball1, &ball2}), gag(gag_), bab1(bab1_), bab2(bab2_), val(val_) {}
        <BLANKLINE>
          SHOC_SCALARS(FakeSubData, 3, 4, gag, bab1, bab2, val)
        };
        """
        arg_data         = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        struct_members   = "\n  ".join(gen_struct_members(arg_data))
        any_arrays       = has_arrays(arg_data)
        struct_name      = get_data_struct_name(sub)
        base_class       = "PhysicsTestData" if is_sugar_compatible(phys, arg_data) else "PhysicsTestDataGeneric"
        inheritance      = " : public {}".format(base_class) if any_arrays else ""
        api              = "\n  " + "\n  ".join(gen_struct_api(phys, struct_name, arg_data) if any_arrays else "")

        result = \
"""struct {struct_name}{inheritance} {{
  {struct_members}{api}
}};""".format(struct_name=struct_name, inheritance=inheritance, struct_members=struct_members, api=api)
        return result

    ###########################################################################
    def gen_cxx_f2c_bind_decl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_f2c_bind_decl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
        void fake_sub_f(Real* foo1, Real* foo2, Real* bar1, Real* bar2, Real* bak1, Real* bak2, Real gag, Real* baz, Int* bag, Int* bab1, Int* bab2, bool val, Int shcol, Int nlev, Int nlevi, Int* ball1, Int* ball2);
        """
        arg_data  = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        arg_decls = gen_arg_cxx_decls(arg_data)

        return "void {sub}_f({arg_sig});".format(sub=sub, arg_sig=", ".join(arg_decls))

    ###########################################################################
    def gen_cxx_f2c_bind_impl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_f2c_bind_impl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
        void fake_sub_f(Real* foo1, Real* foo2, Real* bar1, Real* bar2, Real* bak1, Real* bak2, Real gag, Real* baz, Int* bag, Int* bab1, Int* bab2, bool val, Int shcol, Int nlev, Int nlevi, Int* ball1, Int* ball2)
        {
          // TODO
        }
        <BLANKLINE>
        """
        decl = self.gen_cxx_f2c_bind_decl(phys, sub, force_arg_data=force_arg_data).rstrip(";")

        # TODO - is it possible to fill-out some or all of the implementation?
        result = \
"""{decl}
{{
  // TODO
}}
""".format(decl=decl)
        return result

    ###########################################################################
    def gen_cxx_func_decl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_func_decl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
          KOKKOS_FUNCTION
          static void fake_sub(const uview_1d<const Spack>& foo1, const uview_1d<const Spack>& foo2, const uview_1d<const Spack>& bar1, const uview_1d<const Spack>& bar2, const uview_1d<const Spack>& bak1, const uview_1d<const Spack>& bak2, const Spack& gag, const uview_1d<Spack>& baz, const uview_1d<const Int>& bag, Int& bab1, Int& bab2, const bool& val, const Int& shcol, const Int& nlev, const Int& nlevi, const uview_1d<Int>& ball1, const uview_1d<Int>& ball2);
        """
        arg_data = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        arg_decls = gen_arg_cxx_decls(arg_data, kokkos=True)

        return "  KOKKOS_FUNCTION\n  static void {sub}({arg_sig});".format(sub=sub, arg_sig=", ".join(arg_decls))

    ###########################################################################
    def gen_cxx_incl_impl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_incl_impl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
        # include "shoc_fake_sub_impl.hpp"
        """
        impl_path = get_piece_data(phys, sub, "cxx_func_impl", FILEPATH, self)
        return '# include "{}"'.format(impl_path)

    ###########################################################################
    def gen_cxx_func_impl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_func_impl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
        KOKKOS_FUNCTION
        void Functions<S,D>::fake_sub(const uview_1d<const Spack>& foo1, const uview_1d<const Spack>& foo2, const uview_1d<const Spack>& bar1, const uview_1d<const Spack>& bar2, const uview_1d<const Spack>& bak1, const uview_1d<const Spack>& bak2, const Spack& gag, const uview_1d<Spack>& baz, const uview_1d<const Int>& bag, Int& bab1, Int& bab2, const bool& val, const Int& shcol, const Int& nlev, const Int& nlevi, const uview_1d<Int>& ball1, const uview_1d<Int>& ball2)
        {
          // TODO
          // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
        }
        """
        decl = self.gen_cxx_func_decl(phys, sub, force_arg_data=force_arg_data).rstrip(";").replace("void ", "void Functions<S,D>::").replace("static ", "")
        decl = "\n".join(line.strip() for line in decl.splitlines())

        # I don't think any intelligent guess at an impl is possible here
        result = \
"""{decl}
{{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}}""".format(decl=decl)
        return result

    ###########################################################################
    def gen_cxx_bfb_unit_decl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_bfb_unit_decl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
            struct TestFakeSub;
        """
        test_struct = get_data_test_struct_name(sub)
        return "    struct {};".format(test_struct)

    ###########################################################################
    def gen_cxx_bfb_unit_impl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_bfb_unit_impl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
          static void run_bfb()
          {
            FakeSubData f90_data[] = {
              // TODO
            };
        <BLANKLINE>
            static constexpr Int num_runs = sizeof(f90_data) / sizeof(FakeSubData);
        <BLANKLINE>
            // Generate random input data
            for (auto& d : f90_data) {
              d.randomize();
            }
        <BLANKLINE>
            // Create copies of data for use by cxx. Needs to happen before fortran calls so that
            // inout data is in original state
            FakeSubData cxx_data[] = {
              // TODO
            };
        <BLANKLINE>
            // Assume all data is in C layout
        <BLANKLINE>
            // Get data from fortran
            for (auto& d : f90_data) {
              // expects data in C layout
              fake_sub(d);
            }
        <BLANKLINE>
            // Get data from cxx
            for (auto& d : cxx_data) {
              d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
              fake_sub_f(d.foo1, d.foo2, d.bar1, d.bar2, d.bak1, d.bak2, d.gag, d.baz, d.bag, &d.bab1, &d.bab2, d.val, d.shcol(), d.nlev(), d.nlevi(), d.ball1, d.ball2);
              d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
            }
        <BLANKLINE>
            // Verify BFB results, all data should be in C layout
            for (Int i = 0; i < num_runs; ++i) {
              FakeSubData& d_f90 = f90_data[i];
              FakeSubData& d_cxx = cxx_data[i];
              REQUIRE(d_f90.bab1 == d_cxx.bab1);
              REQUIRE(d_f90.bab2 == d_cxx.bab2);
              for (Int k = 0; k < d_f90.total(d_f90.baz); ++k) {
                REQUIRE(d_f90.total(d_f90.baz) == d_cxx.total(d_cxx.baz));
                REQUIRE(d_f90.baz[k] == d_cxx.baz[k]);
                REQUIRE(d_f90.total(d_f90.baz) == d_cxx.total(d_cxx.ball1));
                REQUIRE(d_f90.ball1[k] == d_cxx.ball1[k]);
                REQUIRE(d_f90.total(d_f90.baz) == d_cxx.total(d_cxx.ball2));
                REQUIRE(d_f90.ball2[k] == d_cxx.ball2[k]);
              }
        <BLANKLINE>
            }
          } // run_bfb
        """
        arg_data = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        data_struct = get_data_struct_name(sub)
        has_array = has_arrays(arg_data)
        need_transpose = needs_transpose(arg_data)
        arg_data_args    = ", ".join(gen_cxx_data_args(phys, arg_data))

        gen_random = "" if not has_array else \
"""

    // Generate random input data
    for (auto& d : f90_data) {
      d.randomize();
    }"""

        c2f_transpose_code = "" if not need_transpose else \
"""
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout"""
        f2c_transpose_code = "" if not need_transpose else \
"""
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout"""

        _, _, _, _, scalars, real_data, int_data = group_data(arg_data, filter_out_intent="in")
        check_scalars, check_arrays = "", ""
        for scalar in scalars:
            check_scalars += "      REQUIRE(d_f90.{name} == d_cxx.{name});\n".format(name=scalar[0])

        all_data = OrderedDict(real_data)
        for k, v in int_data.items():
            if k in all_data:
                all_data[k].extend(v)
            else:
                all_data[k] = v

        for _, data in all_data.items():
            check_arrays += "      for (Int k = 0; k < d_f90.total(d_f90.{}); ++k) {{\n".format(data[0])
            for datum in data:
                check_arrays += "        REQUIRE(d_f90.total(d_f90.{orig}) == d_cxx.total(d_cxx.{name}));\n".format(orig=data[0], name=datum)
                check_arrays += "        REQUIRE(d_f90.{name}[k] == d_cxx.{name}[k]);\n".format(name=datum)

            check_arrays += "      }\n"

        result = \
"""  static void run_bfb()
  {{
    {data_struct} f90_data[] = {{
      // TODO
    }};

    static constexpr Int num_runs = sizeof(f90_data) / sizeof({data_struct});{gen_random}

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    {data_struct} cxx_data[] = {{
      // TODO
    }};

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {{
      // expects data in C layout
      {sub}(d);
    }}

    // Get data from cxx
    for (auto& d : cxx_data) {{{c2f_transpose_code}
      {sub}_f({arg_data_args});{f2c_transpose_code}
    }}

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {{
      {data_struct}& d_f90 = f90_data[i];
      {data_struct}& d_cxx = cxx_data[i];
{check_scalars}{check_arrays}
    }}
  }} // run_bfb""".format(data_struct=data_struct,
                          sub=sub,
                          gen_random=gen_random,
                          c2f_transpose_code=c2f_transpose_code,
                          f2c_transpose_code=f2c_transpose_code,
                          arg_data_args=arg_data_args,
                          check_scalars=check_scalars,
                          check_arrays=check_arrays)

        return result

    ###########################################################################
    def gen_cxx_eti(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_eti("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
        #include "shoc_fake_sub_impl.hpp"
        <BLANKLINE>
        namespace scream {
        namespace shoc {
        <BLANKLINE>
        /*
         * Explicit instantiation for doing fake_sub on Reals using the
         * default device.
         */
        <BLANKLINE>
        template struct Functions<Real,DefaultDevice>;
        <BLANKLINE>
        } // namespace shoc
        } // namespace scream
        <BLANKLINE>
        """
        include_file = get_piece_data(phys, sub, "cxx_func_impl", FILEPATH, self)

        result = \
"""#include "{include_file}"

namespace scream {{
namespace {phys} {{

/*
 * Explicit instantiation for doing {sub} on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

}} // namespace {phys}
}} // namespace scream
""".format(sub=sub, include_file=include_file, phys=phys)

        return result

    ###########################################################################
    def gen_cmake_impl_eti(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cmake_impl_eti("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
            shoc_fake_sub.cpp
        """
        eti_src = get_piece_data(phys, sub, "cxx_eti", FILEPATH, self)
        return "    {}".format(eti_src)

    ###########################################################################
    def gen_cmake_unit_test(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cmake_unit_test("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
            shoc_fake_sub_tests.cpp
        """
        test_src = os.path.basename(get_piece_data(phys, sub, "cxx_bfb_unit_impl", FILEPATH, self))
        return "    {}".format(test_src)

    #
    # Main methods
    #

    ###########################################################################
    def gen_piece(self, phys, sub, piece, force_arg_data=None, force_file_lines=None):
    ###########################################################################
        r"""
        Generate code for a specific piece for a physics+subroutine

        >>> gb = GenBoiler(dry_run=True)
        >>> force_file_lines = [
        ... "fake_line_before_1",
        ... "fake_line_before_2",
        ... "void fake_sub(FakeSubData& d)",
        ... "{",
        ... "  // bad line",
        ... "}",
        ... "fake_line_after_1",
        ... "fake_line_after_2",
        ... ]
        >>> gb.gen_piece("shoc", "fake_sub", "cxx_c2f_glue_impl", force_arg_data=UT_ARG_DATA, force_file_lines=force_file_lines)
        In file shoc_functions_f90.cpp, would replace:
        void fake_sub(FakeSubData& d)
        {
          // bad line
        }
        <BLANKLINE>
        WITH:
        void fake_sub(FakeSubData& d)
        {
          shoc_init(d.nlev(), true);
          d.transpose<ekat::TransposeDirection::c2f>();
          fake_sub_c(d.foo1, d.foo2, d.bar1, d.bar2, d.bak1, d.bak2, d.gag, d.baz, d.bag, &d.bab1, &d.bab2, d.val, d.shcol(), d.nlev(), d.nlevi(), d.ball1, d.ball2);
          d.transpose<ekat::TransposeDirection::f2c>();
        }

        >>> force_file_lines = [
        ... "fake_line_before_1",
        ... "fake_line_before_2",
        ... "// end _c impls",
        ... "fake_line_after_1",
        ... "fake_line_after_2",
        ... ]
        >>> gb.gen_piece("shoc", "fake_sub", "cxx_c2f_glue_impl", force_arg_data=UT_ARG_DATA, force_file_lines=force_file_lines)
        In file shoc_functions_f90.cpp, at line 2, would insert:
        void fake_sub(FakeSubData& d)
        {
          shoc_init(d.nlev(), true);
          d.transpose<ekat::TransposeDirection::c2f>();
          fake_sub_c(d.foo1, d.foo2, d.bar1, d.bar2, d.bak1, d.bak2, d.gag, d.baz, d.bag, &d.bab1, &d.bab2, d.val, d.shcol(), d.nlev(), d.nlevi(), d.ball1, d.ball2);
          d.transpose<ekat::TransposeDirection::f2c>();
        }

        >>> force_file_lines[2:2] = ["void fake_sub(FakeSubData& d)", "{", "}"]
        >>> print("\n".join(force_file_lines))
        fake_line_before_1
        fake_line_before_2
        void fake_sub(FakeSubData& d)
        {
        }
        // end _c impls
        fake_line_after_1
        fake_line_after_2

        >>> force_file_lines = [
        ... "fake_line_before_1",
        ... "fake_line_before_2",
        ... "void fake_sub_c();",
        ... "fake_line_after_1",
        ... "fake_line_after_2",
        ... ]
        >>> gb.gen_piece("shoc", "fake_sub", "cxx_c2f_bind_decl", force_arg_data=UT_ARG_DATA, force_file_lines=force_file_lines)
        In file shoc_functions_f90.cpp, would replace:
        void fake_sub_c();
        <BLANKLINE>
        WITH:
        void fake_sub_c(Real* foo1, Real* foo2, Real* bar1, Real* bar2, Real* bak1, Real* bak2, Real gag, Real* baz, Int* bag, Int* bab1, Int* bab2, bool val, Int shcol, Int nlev, Int nlevi, Int* ball1, Int* ball2);
        """
        if force_arg_data is None: # don't want unit tests printing this
            print("===============================================================================")
            print("Trying to generate piece {} for subroutine {} for physics {}\n".format(piece, sub, phys))

        base_filepath, was_filegen, insert_regex, self_begin_regex, self_end_regex, _ \
            = [item(phys, sub, self) for item in PIECES[piece]]

        filepath = self.get_path_for_piece_file(phys, sub, piece)

        if was_filegen:
            # If freshly generated file, we're done
            pass
        else:
            orig_lines = force_file_lines if force_file_lines else filepath.open().read().splitlines()
            needs_rewrite = False
            gen_lines  = getattr(self, "gen_{}".format(piece))(phys, sub, force_arg_data=force_arg_data).splitlines()

            # Check to see if piece already exists
            try:
                existing_piece_line_range = check_existing_piece(orig_lines, self_begin_regex, self_end_regex)
            except SystemExit as e:
                expect(False, "Problem parsing file {} for existing piece {}: {}".format(filepath, piece, e))

            if existing_piece_line_range is not None:
                # Replace existing
                if self._dry_run:
                    print("In file {}, would replace:\n{}\n\nWITH:\n{}".\
                          format(base_filepath, "\n".join(orig_lines[slice(*existing_piece_line_range)]), "\n".join(gen_lines)))
                elif not self._overwrite:
                    print("Already detected piece {} for subroutine {} in file {}, code:\n{}".\
                          format(piece, sub, base_filepath, "\n".join(orig_lines[slice(*existing_piece_line_range)])))
                else:
                    orig_lines[slice(*existing_piece_line_range)] = gen_lines
                    needs_rewrite = True
            else:
                # Look for place to insert and insert
                insert_line = find_insertion(orig_lines, insert_regex)
                expect(insert_line is not None,
                       "Could not find place to insert generated code for {} in file {} based on regex {}".\
                       format(piece, filepath, insert_regex))

                if self._dry_run:
                    print("In file {}, at line {}, would insert:\n{}".format(base_filepath, insert_line, "\n".join(gen_lines)))
                else:
                    orig_lines[insert_line:insert_line] = gen_lines
                    needs_rewrite = True

            if needs_rewrite:
                with filepath.open("w") as fd:
                    fd.write("\n".join(orig_lines))

                print("SUCCESS")

    ###########################################################################
    def gen_boiler(self):
    ###########################################################################
        all_success = True
        for sub in self._subs:
            for phys in self._physics:
                for piece in self._pieces:
                    try:
                        self.gen_piece(phys, sub, piece)
                    except SystemExit as e:
                        print("Warning: failed to generate subroutine {} piece {} for physics {}, error: {}".\
                              format(sub, piece, phys, e))
                        all_success = False

        return all_success
