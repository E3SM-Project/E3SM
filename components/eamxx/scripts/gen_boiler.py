from utils import expect
from git_utils import get_git_toplevel_dir

from collections import OrderedDict
import re
from pathlib import Path

#
# Global hardcoded data
#

# Templates: maps piece name to generic file text
FILE_TEMPLATES = {
    "cxx_bfb_unit_impl": lambda phys, sub, gen_code:
f"""#include "catch2/catch.hpp"

#include "share/eamxx_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/{phys}/{phys}_functions.hpp"
#include "physics/{phys}/{phys}_functions_f90.hpp"

#include "{phys}_unit_tests_common.hpp"

namespace scream {{
namespace {phys} {{
namespace unit_test {{

template <typename D>
struct UnitWrap::UnitTest<D>::{get_data_test_struct_name(sub)} {{

{gen_code}

}};

}} // namespace unit_test
}} // namespace {phys}
}} // namespace scream

namespace {{

TEST_CASE("{sub}_bfb", "[{phys}]")
{{
  using TestStruct = scream::{phys}::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::{get_data_test_struct_name(sub)};

  TestStruct::run_bfb();
}}

}} // empty namespace
""",

###############################################################################

    "cxx_func_impl": lambda phys, sub, gen_code:
f"""#ifndef {phys.upper()}_{sub.upper()}_IMPL_HPP
#define {phys.upper()}_{sub.upper()}_IMPL_HPP

#include "{phys}_functions.hpp" // for ETI only but harmless for GPU

namespace scream {{
namespace {phys} {{

/*
 * Implementation of {phys} {sub}. Clients should NOT
 * #include this file, but include {phys}_functions.hpp instead.
 */

template<typename S, typename D>
{gen_code}

}} // namespace {phys}
}} // namespace scream

#endif
""",

###############################################################################

    "cxx_eti": lambda phys, sub, gen_code: gen_code
}

# piece map. maps the name of a piece of boilerplate that needs to be genereate to:
#   (filepath (relative to cxx_root))
FILEPATH, FILECREATE, INSERT_REGEX, ID_SELF_BEGIN_REGEX, ID_SELF_END_REGEX, DESC = range(6)
PIECES = OrderedDict([
    ("f90_c2f_bind", (
        lambda phys, sub, gb: f"{phys}_iso_c.f90",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "f90_c2f_bind"),
        lambda phys, sub, gb: re.compile(fr"^\s*end\s+module\s{phys}_iso_c"), # put at end of module
        lambda phys, sub, gb: get_subroutine_begin_regex(sub + "_c"), # sub_c begin
        lambda phys, sub, gb: get_subroutine_end_regex(sub + "_c"),    # sub_c end
        lambda *x           : "The c to f90 fortran subroutine(<name>_c)"
    )),

    ("f90_f2c_bind"  , (
        lambda phys, sub, gb: f"{phys}_iso_f.f90",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "f90_f2c_bind"),
        lambda phys, sub, gb: re.compile(r"^\s*end\s+interface"), # put at end of interface
        lambda phys, sub, gb: get_subroutine_begin_regex(sub + "_f"), # sub_f begin
        lambda phys, sub, gb: get_subroutine_end_regex(sub + "_f"),   # sub_f begin
        lambda *x           : "The f90 to c fortran subroutine(<name>_f)"
    )),

    ("cxx_c2f_bind_decl"  , (
        lambda phys, sub, gb: f"{phys}_functions_f90.cpp",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_bind_decl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(comment='extern "C" : end _c decls'), # reqs special comment
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub + "_c"), # cxx_c decl
        lambda phys, sub, gb: re.compile(r".*;\s*$"),                   # ;
        lambda *x           : "The c to f90 c function declaration(<name>_c)"
    )),

    ("cxx_c2f_glue_decl"  , (
        lambda phys, sub, gb: f"{phys}_functions_f90.hpp",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_glue_decl"),
        lambda phys, sub, gb: re.compile(r'^\s*extern\s+"C"'), # put before _f decls
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub), # cxx(data) decl
        lambda phys, sub, gb: re.compile(r".*;\s*"),             # ;
        lambda *x           : "The cxx to c function declaration(<name>(Data))"
    )),

    ("cxx_c2f_glue_impl" , (
        lambda phys, sub, gb: f"{phys}_functions_f90.cpp",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_glue_impl"),
        lambda phys, sub, gb: re.compile(r"^\s*// end _c impls"), # reqs special comment
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub), # cxx(data)
        lambda phys, sub, gb: get_cxx_close_block_regex(at_line_start=True), # terminating }
        lambda *x           : "The cxx to c function implementation(<name>(Data))"
    )),

    ("cxx_c2f_data"  , (
        lambda phys, sub, gb: f"{phys}_functions_f90.hpp",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_data"),
        lambda phys, sub, gb: re.compile(r"^\s*// Glue functions to call fortran"),  # reqs special comment
        lambda phys, sub, gb: get_cxx_struct_begin_regex(get_data_struct_name(sub)), # struct Sub
        lambda phys, sub, gb: get_cxx_close_block_regex(semicolon=True),             # terminating };
        lambda *x           : "The cxx data struct definition(struct Data)"
    )),

    ("cxx_f2c_bind_decl"  , (
        lambda phys, sub, gb: f"tests/infra/{phys}_test_data.hpp",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_f2c_bind_decl"),
        lambda phys, sub, gb: get_plain_comment_regex(comment="end _host function decls"), # reqs special comment
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub + "_host"), # cxx_host decl
        lambda phys, sub, gb: re.compile(r".*;\s*$"),                   # ;
        lambda *x           : "The f90 to cxx function declaration(<name>_host)"
    )),

    ("cxx_f2c_bind_impl"  , (
        lambda phys, sub, gb: f"tests/infra/{phys}_test_data.cpp",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_f2c_bind_impl"),
        lambda phys, sub, gb: get_namespace_close_regex(phys),          # insert at end of namespace
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub + "_host"),      # cxx_f
        lambda phys, sub, gb: get_cxx_close_block_regex(at_line_start=True), # terminating }
        lambda *x           : "The f90 to cxx function implementation(<name>_host)"
    )),

    ("cxx_func_decl", (
        lambda phys, sub, gb: f"{phys}_functions.hpp",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_func_decl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(semicolon=True, comment="struct Functions"), # end of struct, reqs special comment
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub, static=True), # cxx decl
        lambda phys, sub, gb: re.compile(r".*;\s*$"),            # ;
        lambda *x           : "The cxx kokkos function declaration(<name>)"
    )),

    ("cxx_incl_impl", (
        lambda phys, sub, gb: f"{phys}_functions.hpp",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_incl_impl"),
        lambda phys, sub, gb: re.compile(r"^\s*#\s*endif\s+//\s*GPU"), # insert at end of impl includes, reqs special comment
        lambda phys, sub, gb: re.compile(r'^\s*#\s*include\s+"{}"'.format(get_piece_data(phys, sub, "cxx_func_impl", FILEPATH, gb))),
        lambda phys, sub, gb: re.compile(r".*"),
        lambda *x           : "The include of *impl.hpp file at bottom of main hpp"
    )),

    ("cxx_func_impl", (
        lambda phys, sub, gb: f"impl/{phys}_{sub}_impl.hpp",
        lambda phys, sub, gb: create_template(phys, sub, gb, "cxx_func_impl"),
        lambda phys, sub, gb: get_namespace_close_regex(phys),   # insert at end of namespace
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub, template="Functions<S,D>"), # cxx begin
        lambda phys, sub, gb: get_cxx_close_block_regex(at_line_start=True), # terminating }
        lambda *x           : "The cxx kokkos function stub implementation(<name>)"
    )),

    ("cxx_bfb_unit_decl", (
        lambda phys, sub, gb: f"tests/{phys}_unit_tests_common.hpp",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_bfb_unit_decl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(semicolon=True), # insert at end of test struct
        lambda phys, sub, gb: get_cxx_struct_begin_regex(get_data_test_struct_name(sub)), # struct decl
        lambda phys, sub, gb: re.compile(r".*;\s*$"), # end of struct decl
        lambda *x           : "The cxx unit test struct declaration"
    )),

    ("cxx_bfb_unit_impl", (
        lambda phys, sub, gb: f"tests/{phys}_{sub}_tests.cpp",
        lambda phys, sub, gb: create_template(phys, sub, gb, "cxx_bfb_unit_impl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(semicolon=True, at_line_start=True), # insert of end of struct
        lambda phys, sub, gb: get_cxx_function_begin_regex("run_bfb", static=True),  # run_bfb
        lambda phys, sub, gb: get_cxx_close_block_regex(comment="run_bfb"), # } // run_bfb # reqs special comment
        lambda *x           : "The cxx bfb unit test implementation"
    )),

    ("cxx_eti", (
        lambda phys, sub, gb: f"eti/{phys}_{sub}.cpp",
        lambda phys, sub, gb: create_template(phys, sub, gb, "cxx_eti"),
        lambda phys, sub, gb: re.compile(".*"), # insert at top of file
        lambda phys, sub, gb: re.compile(".*"), # start at top of file
        lambda phys, sub, gb: get_namespace_close_regex("scream"), #end of file
        lambda *x           : "The cxx explicit template instatiation file"
    )),

    ("cmake_impl_eti", (
        lambda phys, sub, gb: "CMakeLists.txt",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cmake_impl_eti"),
        lambda phys, sub, gb: re.compile(fr".*[)]\s*#\s*{phys.upper()} ETI SRCS"), # insert at end of ETI src list, reqs special comment
        lambda phys, sub, gb: re.compile(fr".*{get_piece_data(phys, sub, 'cxx_eti', FILEPATH, gb)}"),
        lambda phys, sub, gb: re.compile(".*"),
        lambda *x           : "Make cmake aware of the ETI file if not cuda build"
    )),

    ("cmake_unit_test", (
        lambda phys, sub, gb: "tests/CMakeLists.txt",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cmake_unit_test"),
        lambda phys, sub, gb: re.compile(fr".*[)]\s*#\s*{phys.upper()}_TESTS_SRCS"), # insert at end of test src list, reqs special comment
        lambda phys, sub, gb: re.compile(fr".*{Path(get_piece_data(phys, sub, 'cxx_bfb_unit_impl', FILEPATH, gb)).name}"),
        lambda phys, sub, gb: re.compile(".*"),
        lambda *x           : "Make cmake aware of the unit test"
    )),

])

# physics map. maps the name of a physics packages containing the original fortran subroutines to:
#   (path-to-origin, path-to-cxx-src)
ORIGIN_FILES, CXX_ROOT, INIT_CODE = range(3)
PHYSICS = {
    "p3"   : (
        ("components/eam/src/physics/cam/micro_p3.F90",),
        "components/eamxx/src/physics/p3",
        "p3_init();"
    ),
    "shoc" : (
        ("components/eam/src/physics/cam/shoc.F90",),
        "components/eamxx/src/physics/shoc",
        "shoc_init(d.nlev, true);"
    ),
    "dp" : (
        ("components/eam/src/control/apply_iop_forcing.F90", "components/eam/src/dynamics/se/se_iop_intr_mod.F90", "components/eam/src/control/iop_data_mod.F90", "components/eam/src/control/history_iop.F90"),
        "components/eamxx/src/physics/dp",
        "dp_init(d.plev, true);"
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
    ("tracerd1", "real", "in", ("shcol","nlev", "ntracers")),
    ("tracerd2", "real", "in", ("shcol","nlev", "ntracers")),
    ("gag", "real", "in", None),
    ("baz", "real", "inout", ("shcol",)),
    ("bag", "integer", "in", ("shcol",)),
    ("bab1", "integer", "out", None),
    ("bab2", "integer", "out", None),
    ("val", "logical", "in", None),
    ("vals", "logical", "in", ("shcol",)),
    ("shcol", "integer", "in", None),
    ("nlev", "integer", "in", None),
    ("nlevi", "integer", "in", None),
    ("ntracers", "integer", "in", None),
    ("ball1", "integer", "out", ("shcol",)),
    ("ball2", "integer", "out", ("shcol",)),
]

UT_ARG_DATA_ALL_SCALAR = [
    ("foo1", "real", "in", None),
    ("foo2", "real", "in", None),
    ("bar1", "real", "inout", None),
    ("bar2", "real", "inout", None),
    ("baz1", "real", "out", None),
    ("baz2", "real", "out", None),
    ("gag1", "integer", "in", None),
    ("gag2", "integer", "in", None),
    ("gal1", "integer", "inout", None),
    ("gal2", "integer", "inout", None),
    ("bal1", "integer", "out", None),
    ("bal2", "integer", "out", None),
    ("bit1", "logical", "in", None),
    ("bit2", "logical", "in", None),
    ("gut1", "logical", "inout", None),
    ("gut2", "logical", "inout", None),
    ("gat1", "logical", "out", None),
    ("gat2", "logical", "out", None),
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
    subroutine_begin_regex_str = fr"^\s*subroutine\s+{name}\s*[(]"
    return re.compile(subroutine_begin_regex_str)

###############################################################################
def get_function_begin_regex(name):
###############################################################################
    """
    >>> get_function_begin_regex("fake_sub").match("function fake_sub(foo, bar) result(baz)").groups()[-1]
    'baz'
    >>> get_function_begin_regex("fake_sub").match("pure function fake_sub(foo, bar) result(baz)").groups()[-1]
    'baz'
    >>> get_function_begin_regex("fake_sub").match("pure function fake_sub() result(baz)").groups()[-1]
    'baz'
    >>> get_function_begin_regex("fake_sub").match("  pure  function  fake_sub ( foo, bar )   result (  baz)").groups()[-1]
    'baz'
    >>> bool(get_function_begin_regex("fake_sub").match("function fake_sub2(foo, bar) result(baz)"))
    False
    >>> bool(get_function_begin_regex("fake_sub").match("! function fake_sub(foo, bar) result(baz)"))
    False
    >>> bool(get_function_begin_regex("fake_sub").match("end function fake_sub"))
    False
    """
    function_begin_regex_str = fr"^\s*((pure\s+)?function)\s+{name}\s*[(].*result\s*[(]\s*([^) ]+)"
    return re.compile(function_begin_regex_str)

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
    >>> bool(get_subroutine_end_regex("fake_sub").match("end function fake_sub"))
    True
    >>> bool(get_subroutine_end_regex("fake_sub").match("end function fake_sub_2"))
    False
    """
    subroutine_end_regex_str = fr"^\s*end\s+(subroutine|function)\s+{name}\s*$"
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
    template_regex_str = fr"{template}::" if template else ""
    function_begin_regex_str = fr"^\s*{static_regex_str}void\s+{template_regex_str}{name}\s*[(]"
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
    comment_regex_str   = fr"\s*//\s*{comment}" if comment else ""
    close_block_regex_str = re.compile(fr"^{line_start_regex_str}}}{semicolon_regex_str}{comment_regex_str}\s*$")
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
    return get_cxx_close_block_regex(comment=fr"namespace\s+{namespace}")

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
    struct_regex_str = fr"^\s*struct\s+{struct}([\W]|$)"
    return re.compile(struct_regex_str)

###############################################################################
def get_plain_comment_regex(comment):
###############################################################################
    comment_regex_str   = fr"^\s*//\s*{comment}"
    return re.compile(comment_regex_str)

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
    return f"Test{get_data_struct_name(sub)[:-4]}"

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
    expect(filepath.exists(),
           f"For generating {sub}'s {piece} for phyiscs {physics}, expected file {filepath} to already exist")
    return False # File was not created

###############################################################################
def create_template(physics, sub, gb, piece, force=False, force_arg_data=None):
###############################################################################
    """
    Create a file based on a template if it doesn't exist. Return True if a file was created.

    >>> gb = GenBoiler(["linear_interp"], ["cxx_func_impl"], dry_run=True)
    >>> create_template("shoc", "linear_interp", gb, "cxx_func_impl", force=True, force_arg_data=UT_ARG_DATA) #doctest: +ELLIPSIS
    Would create file .../components/eamxx/src/physics/shoc/impl/shoc_linear_interp_impl.hpp with contents:
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
    void Functions<S,D>::linear_interp(const uview_1d<const Spack>& foo1, const uview_1d<const Spack>& foo2, const uview_1d<const Spack>& bar1, const uview_1d<const Spack>& bar2, const uview_1d<const Spack>& bak1, const uview_1d<const Spack>& bak2, const uview_1d<const Spack>& tracerd1, const uview_1d<const Spack>& tracerd2, const Spack& gag, const uview_1d<Spack>& baz, const uview_1d<const Int>& bag, Int& bab1, Int& bab2, const bool& val, const uview_1d<const bool>& vals, const Int& shcol, const Int& nlev, const Int& nlevi, const Int& ntracers, const uview_1d<Int>& ball1, const uview_1d<Int>& ball2)
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
               f"{filepath} does not exist and there is no template for generating files for piece {piece}")

        gen_code = getattr(gb, f"gen_{piece}")(physics, sub, force_arg_data=force_arg_data)
        contents = FILE_TEMPLATES[piece](physics, sub, gen_code)
        if gb.dry_run():
            print(f"Would create file {filepath} with contents:\n{contents}")
        else:
            with filepath.open("w", encoding="utf-8") as fd:
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
            top_splits[-1] += f",{raw_split}"

        balanced = top_splits[-1].count("(") == top_splits[-1].count(")")

    return top_splits

###############################################################################
def get_arg_order(line):
###############################################################################
    """
    Given a line of fortran declaring a subroutine, return the arg names in order

    >>> get_arg_order("subroutine p3_set_tables( mu_r_user, revap_user,vn_user, vm_user )")
    ['mu_r_user', 'revap_user', 'vn_user', 'vm_user']
    >>> get_arg_order("function p3_set_tables( mu_r_user, revap_user,vn_user, vm_user ) result(bar)")
    ['mu_r_user', 'revap_user', 'vn_user', 'vm_user', 'bar']
    >>> get_arg_order("pure function p3_set_tables( mu_r_user, revap_user,vn_user, vm_user ) result( bar)")
    ['mu_r_user', 'revap_user', 'vn_user', 'vm_user', 'bar']
    >>> get_arg_order("pure function p3_set_tables(mu_r_user,revap_user,vn_user,vm_user) result(bar)")
    ['mu_r_user', 'revap_user', 'vn_user', 'vm_user', 'bar']
    """
    tokens = line.split()
    if tokens[0] == "function" or tokens[1] == "function":
        first_paren = False
        first_paren_contents = ""
        for c in line:
            if c == "(":
                expect(not first_paren, f"Bad line, multiple opening parens: {line}")
                first_paren = True
            elif c == ")":
                break
            elif first_paren:
                first_paren_contents += c

        basic_args = [item.strip() for item in first_paren_contents.split(",") if item.strip()]

        basic_args.append(line.rstrip(")").split(")")[-1].split("(")[-1].strip())
        return basic_args

    else:
        # subroutine logic
        args_raw = line.rstrip(")").split("(", maxsplit=1)[-1]
        return [item.strip() for item in args_raw.split(",") if item.strip()]

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
    >>> parse_f90_args('type(element_t), intent(inout) :: elem(:)')
    [('elem', 'type::element_t', 'inout', (':',))]
    >>> parse_f90_args('character*(max_path_len), intent(out), optional ::  iopfile_out')
    [('iopfile_out', 'type::string', 'out', None)]
    """
    expect(line.count("::") == 1, f"Expected line format 'type-info :: names' for: {line}")
    metadata_str, names_str = line.split("::")
    names_dims = split_top_commas(names_str)
    metadata   = split_top_commas(metadata_str)

    argtoken = metadata[0]
    argtype = argtoken.split("(")[0].strip()
    if argtype == "type":
        expect("(" in argtoken, f"Undefined type for {argtoken}")
        argtype += ("::" + argtoken.split("(")[1].strip().rstrip(")"))
    elif argtype == "character*":
        argtype = "type::string"

    intent, dims = None, None
    for metadatum in metadata:
        if metadatum.startswith("intent"):
            expect(intent is None, f"Multiple intents in line: {line}")
            intent = metadatum.split("(")[-1].rstrip(")").strip()
        elif metadatum.startswith("dimension"):
            expect(dims is None, f"Multiple dimensions in line: {line}")
            dims_raw = metadatum.split("(")[-1].rstrip(")").strip()
            dims = tuple(item.replace(" ", "") for item in dims_raw.split(","))

    names = []
    for name_dim in names_dims:
        if "(" in name_dim:
            name, dims_raw = name_dim.split("(")
            dims_raw = dims_raw.rstrip(")").strip()
            dims_check = tuple(item.replace(" ", "") for item in dims_raw.split(","))
            expect(dims is None or dims_check == dims, f"Inconsistent dimensions in line: {line}")
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
    ...
    ... function impli_srf_stress_term(shcol, rho_zi_sfc, uw_sfc, vw_sfc, u_wind_sfc, v_wind_sfc) result (ksrf)
    ...   !intent-ins
    ...   integer,     intent(in) :: shcol
    ...
    ...   !air density at interfaces [kg/m3]
    ...   real(rtype), intent(in) :: rho_zi_sfc(shcol)
    ...   !vertical zonal momentum flux at surface [m3/s3]
    ...   real(rtype), intent(in) :: uw_sfc(shcol)
    ...   !vertical meridional momentum flux at surface [m3/s3]
    ...   real(rtype), intent(in) :: vw_sfc(shcol)
    ...   !zonal wind [m/s]
    ...   real(rtype), intent(in) :: u_wind_sfc(shcol)
    ...   !meridional wind [m/s]
    ...   real(rtype), intent(in) :: v_wind_sfc(shcol)
    ...
    ...   !function return value
    ...   real(rtype) :: ksrf(shcol)
    ...
    ...   return foo
    ...  end function impli_srf_stress_term
    ...
    ...  subroutine advance_iop_forcing(scm_dt, ps_in, &             ! In
    ...                    u_in, v_in, t_in, q_in, t_phys_frc,&    ! In
    ...                    u_update, v_update, t_update, q_update) ! Out
    ...
    ...    ! Input arguments
    ...    real(r8), intent(in) :: ps_in             ! surface pressure [Pa]
    ...    real(r8), intent(in) :: u_in(plev)        ! zonal wind [m/s]
    ...    real(r8), intent(in) :: v_in(plev)        ! meridional wind [m/s]
    ...    real(r8), intent(in) :: t_in(plev)        ! temperature [K]
    ...    real(r8), intent(in) :: q_in(plev,pcnst)  ! q tracer array [units vary]
    ...    real(r8), intent(in) :: t_phys_frc(plev)  ! temperature forcing from physics [K/s]
    ...    real(r8), intent(in) :: scm_dt            ! model time step [s]
    ...
    ...    ! Output arguments
    ...    real(r8), intent(out) :: t_update(plev)      ! updated temperature [K]
    ...    real(r8), intent(out) :: q_update(plev,pcnst)! updated q tracer array [units vary]
    ...    real(r8), intent(out) :: u_update(plev)      ! updated zonal wind [m/s]
    ...    real(r8), intent(out) :: v_update(plev)      ! updated meridional wind [m/s]
    ...
    ...  end subroutine advance_iop_forcing
    ...
    ...  subroutine iop_setinitial(elem)
    ...    type(element_t), intent(inout) :: elem(:)
    ...  end subroutine iop_setinitial
    ... '''
    >>> print("\n".join([str(item) for item in sorted(parse_origin(teststr, ["p3_get_tables", "p3_init_b"]).items())]))
    ('p3_get_tables', [('mu_r_user', 'real', 'out', ('150',)), ('revap_user', 'real', 'out', ('300', '10')), ('tracerd', 'real', 'out', ('300', '10', '42')), ('vn_user', 'real', 'out', ('300', '10')), ('vm_user', 'real', 'out', ('300', '10'))])
    ('p3_init_b', [])
    >>> print("\n".join([str(item) for item in parse_origin(teststr, ["impli_srf_stress_term"]).items()]))
    ('impli_srf_stress_term', [('shcol', 'integer', 'in', None), ('rho_zi_sfc', 'real', 'in', ('shcol',)), ('uw_sfc', 'real', 'in', ('shcol',)), ('vw_sfc', 'real', 'in', ('shcol',)), ('u_wind_sfc', 'real', 'in', ('shcol',)), ('v_wind_sfc', 'real', 'in', ('shcol',)), ('ksrf', 'real', 'out', ('shcol',))])
    >>> print("\n".join([str(item) for item in parse_origin(teststr, ["advance_iop_forcing"]).items()]))
    ('advance_iop_forcing', [('plev', 'integer', 'in', None), ('pcnst', 'integer', 'in', None), ('scm_dt', 'real', 'in', None), ('ps_in', 'real', 'in', None), ('u_in', 'real', 'in', ('plev',)), ('v_in', 'real', 'in', ('plev',)), ('t_in', 'real', 'in', ('plev',)), ('q_in', 'real', 'in', ('plev', 'pcnst')), ('t_phys_frc', 'real', 'in', ('plev',)), ('u_update', 'real', 'out', ('plev',)), ('v_update', 'real', 'out', ('plev',)), ('t_update', 'real', 'out', ('plev',)), ('q_update', 'real', 'out', ('plev', 'pcnst'))])
    >>> print("\n".join([str(item) for item in parse_origin(teststr, ["iop_setinitial"]).items()]))
    ('iop_setinitial', [('elem', 'type::element_t', 'inout', (':',))])
    """
    begin_sub_regexes  = [get_subroutine_begin_regex(sub) for sub in subs]
    begin_func_regexes = [get_function_begin_regex(sub)   for sub in subs]
    arg_decl_regex = re.compile(r"^.+intent\s*[(]\s*(in|out|inout)\s*[)]")

    contents = normalize_f90(contents)

    db = {}
    active_sub = None
    result_name = None
    arg_order = []
    arg_decls = []
    for line in contents.splitlines():
        for sub, begin_sub_regex, begin_func_regex in zip(subs, begin_sub_regexes, begin_func_regexes):
            begin_sub_match = begin_sub_regex.match(line)
            begin_func_match = begin_func_regex.match(line)
            if begin_sub_match is not None:
                expect(active_sub is None, f"subroutine {active_sub} was still active when {sub} began")
                active_sub = sub
                arg_order = get_arg_order(line)
            elif begin_func_match is not None:
                expect(active_sub is None, f"subroutine {active_sub} was still active when {sub} began")
                active_sub = sub
                arg_order = get_arg_order(line)
                result_name = begin_func_match.groups()[-1]

        if active_sub:
            decl_match = arg_decl_regex.match(line)
            if decl_match is not None:
                arg_decls.extend(parse_f90_args(line))
            elif result_name:
                result_decl_regex = re.compile(fr".+::\s*{result_name}([^\w]|$)")
                result_decl_match = result_decl_regex.match(line)
                if result_decl_match is not None:
                    line = line.replace("::", " , intent(out) ::")
                    arg_decls.extend(parse_f90_args(line))

            end_regex = get_subroutine_end_regex(active_sub)
            end_match = end_regex.match(line)
            if end_match is not None:
                expect(active_sub not in db, f"Found multiple matches for {active_sub}")
                expect(len(arg_order) == len(arg_decls),
                       f"Number of decls:\n{arg_decls}\nDid not match arg list: {arg_order}")

                # we need our decls to be ordered based on arg list order
                ordered_decls = []
                for arg in arg_order:
                    found = False
                    for arg_decl in arg_decls:
                        if arg_decl[ARG_NAME] == arg:
                            ordered_decls.append(arg_decl)
                            found = True
                            break

                    expect(found, f"Could not find decl for arg {arg} in\n{arg_decls}")

                # Dim resolution. Arrays with global dims must have the
                # dim as an input in the converted code.
                global_ints_to_insert = []
                arg_names = set()
                for arg_datum in ordered_decls:
                    arg_name = arg_datum[ARG_NAME]
                    arg_names.add(arg_name)

                for arg_datum in ordered_decls:
                    arg_dims = arg_datum[ARG_DIMS]
                    if arg_dims is not None:
                        for arg_dim in arg_dims:
                            if not arg_dim.isdigit() and arg_dim not in arg_names and arg_dim != ":":
                                global_ints_to_insert.append((arg_dim, "integer", "in", None))
                                arg_names.add(arg_dim)

                db[active_sub] = global_ints_to_insert + ordered_decls
                active_sub = None
                result_name = None
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
    >>> gen_arg_f90_decl('type::element_t', 'inout', (':',), ["foo"])
    'type(c_ptr) , intent(inout), dimension(:) :: foo'
    """
    value  = ", value" if dims is None and intent == "in" else ""
    intent_s = f", intent({intent})"
    dimension_s = f", dimension({', '.join(dims)})" if dims is not None else ""
    names_s = ", ".join(names)

    if is_custom_type(argtype):
        return f"type(c_ptr) {intent_s}{dimension_s} :: {names_s}"
    else:
        expect(argtype in C_TYPE_MAP, f"Unrecognized argtype for C_TYPE_MAP: {argtype}")
        c_type = C_TYPE_MAP[argtype]
        return f"{argtype}(kind={c_type}) {value}{intent_s}{dimension_s} :: {names_s}"

###############################################################################
def is_custom_type(arg_type):
###############################################################################
    return arg_type.startswith("type::")

CXX_TYPE_MAP = {"real" : "Real", "integer" : "Int", "logical" : "bool"}
###############################################################################
def get_cxx_scalar_type(arg_type):
###############################################################################
    if is_custom_type(arg_type):
        arg_cxx_type = arg_type.split("::")[-1]
    else:
        expect(arg_type in CXX_TYPE_MAP, f"Unrecognized argtype for CXX_TYPE_MAP: {arg_type}")
        arg_cxx_type = CXX_TYPE_MAP[arg_type]

    return arg_cxx_type

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
    >>> get_cxx_type(('elem', 'type::element_t', 'inout', (':',)))
    'element_t*'
    """
    is_ptr = arg_datum[ARG_DIMS] is not None or arg_datum[ARG_INTENT] != "in"
    arg_type = arg_datum[ARG_TYPE]
    arg_cxx_type = get_cxx_scalar_type(arg_type)
    return f"{arg_cxx_type}{'*' if is_ptr else ''}"

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
    >>> get_kokkos_type(('elem', 'type::element_t', 'inout', (':',)))
    'const uview_1d<element_t>&'
    """
    is_const  = arg_datum[ARG_INTENT] == "in"
    is_view   = arg_datum[ARG_DIMS] is not None
    arg_type  = arg_datum[ARG_TYPE]
    if is_custom_type(arg_type):
        kokkos_type = arg_type.split("::")[-1]
    else:
        kokkos_type = KOKKOS_TYPE_MAP[arg_type]

    base_type = f"{'const ' if is_const else ''}{kokkos_type}"

    # We assume 1d even if the f90 array is 2d since we assume c++ will spawn a kernel
    # over one of the dimensions
    return f"const uview_1d<{base_type}>&" if is_view else f"{base_type}&"

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
    arg_sig_list = [f"{arg_type} {arg_name}" for arg_name, arg_type in zip(arg_names, arg_types)]
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
def split_by_intent(arg_data):
###############################################################################
    """
    Take arg data and split into three lists of names based on intent: [inputs], [intouts], [outputs]

    >>> split_by_intent(UT_ARG_DATA)
    (['foo1', 'foo2', 'bar1', 'bar2', 'bak1', 'bak2', 'tracerd1', 'tracerd2', 'gag', 'bag', 'val', 'vals', 'shcol', 'nlev', 'nlevi', 'ntracers'], ['baz'], ['bab1', 'bab2', 'ball1', 'ball2'])
    """
    inputs, inouts, outputs = [], [], []
    for name, _, intent, _ in arg_data:
        if intent == "in":
            inputs.append(name)
        elif intent == "inout":
            inouts.append(name)
        elif intent == "out":
            outputs.append(name)
        else:
            expect(False, f"Unhandled intent: {intent}")

    return inputs, inouts, outputs

###############################################################################
def split_by_type(arg_data):
###############################################################################
    """
    Take arg data and split into three lists of names based on type: [reals], [ints], [logicals]

    >>> split_by_type(UT_ARG_DATA)
    (['foo1', 'foo2', 'bar1', 'bar2', 'bak1', 'bak2', 'tracerd1', 'tracerd2', 'gag', 'baz'], ['bag', 'bab1', 'bab2', 'shcol', 'nlev', 'nlevi', 'ntracers', 'ball1', 'ball2'], ['val', 'vals'])
    """
    reals, ints, logicals = [], [], []
    for name, argtype, _, _ in arg_data:
        if argtype == "real":
            reals.append(name)
        elif argtype == "integer":
            ints.append(name)
        elif argtype == "logical":
            logicals.append(name)
        elif is_custom_type(argtype):
            pass
        else:
            expect(False, f"Unhandled argtype: {argtype}")

    return reals, ints, logicals

###############################################################################
def split_by_scalar_vs_view(arg_data):
###############################################################################
    """
    Take arg data and split into two lists of names based on scalar/not-scalar: [scalars] [non-scalars]
    """
    scalars, non_scalars = [], []
    for name, _, _, dims in arg_data:
        if dims is not None:
            non_scalars.append(name)
        else:
            scalars.append(name)

    return scalars, non_scalars

###############################################################################
def gen_cxx_data_args(physics, arg_data):
###############################################################################
    """
    Based on data, generate unpacking of Data struct args

    >>> gen_cxx_data_args("shoc", UT_ARG_DATA)
    ['d.foo1', 'd.foo2', 'd.bar1', 'd.bar2', 'd.bak1', 'd.bak2', 'd.tracerd1', 'd.tracerd2', 'd.gag', 'd.baz', 'd.bag', '&d.bab1', '&d.bab2', 'd.val', 'd.vals', 'd.shcol', 'd.nlev', 'd.nlevi', 'd.ntracers', 'd.ball1', 'd.ball2']
    """
    all_dims = group_data(arg_data)[3]
    args_needs_ptr = [item[ARG_DIMS] is None and item[ARG_INTENT] != "in" for item in arg_data]
    arg_names      = [item[ARG_NAME] for item in arg_data]
    arg_dim_call   = [item[ARG_NAME] in all_dims for item in arg_data]
    args = [f"{'&' if need_ptr else ''}d.{arg_name}"
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
    real(kind=c_real) , intent(in), dimension(shcol, nlev, ntracers) :: tracerd1, tracerd2
    real(kind=c_real) , value, intent(in) :: gag
    real(kind=c_real) , intent(inout), dimension(shcol) :: baz
    integer(kind=c_int) , intent(in), dimension(shcol) :: bag
    integer(kind=c_int) , intent(out) :: bab1, bab2
    logical(kind=c_bool) , value, intent(in) :: val
    logical(kind=c_bool) , intent(in), dimension(shcol) :: vals
    integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi, ntracers
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
    Real *foo1, *foo2, *bar1, *bar2, *bak1, *bak2, *tracerd1, *tracerd2;
    Real gag;
    Int *bag;
    bool val;
    bool *vals;
    Int shcol, nlev, nlevi, ntracers;
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
            result.append(f"// {comment}")
            type_map = metadata[intent]
            for type_info, names in type_map.items():
                type_name, is_ptr = type_info
                decl_str = get_cxx_scalar_type(type_name)
                decl_str += f" {', '.join(['{}{}'.format('*' if is_ptr else '', name) for name in names])};"
                result.append(decl_str)

            result.append("")

    return result

###############################################################################
def group_data(arg_data, filter_out_intent=None, filter_scalar_custom_types=False):
###############################################################################
    r"""
    Given data, return ([fst_dims], [snd_dims], [trd_dims], [all-dims], [scalars], {dims->[real_data]}, {dims->[int_data]}, {dims->[bool_data]})

    >>> print("\n".join([str(item) for item in group_data(UT_ARG_DATA)]))
    ['shcol']
    ['nlev', 'nlevi']
    ['ntracers']
    ['shcol', 'nlev', 'nlevi', 'ntracers']
    [('gag', 'Real'), ('bab1', 'Int'), ('bab2', 'Int'), ('val', 'bool')]
    OrderedDict([(('shcol',), ['foo1', 'foo2', 'baz']), (('shcol', 'nlev'), ['bar1', 'bar2']), (('shcol', 'nlevi'), ['bak1', 'bak2']), (('shcol', 'nlev', 'ntracers'), ['tracerd1', 'tracerd2'])])
    OrderedDict([(('shcol',), ['bag', 'ball1', 'ball2'])])
    OrderedDict([(('shcol',), ['vals'])])
    >>> print("\n".join([str(item) for item in group_data(UT_ARG_DATA_ALL_SCALAR)]))
    []
    []
    []
    []
    [('foo1', 'Real'), ('foo2', 'Real'), ('bar1', 'Real'), ('bar2', 'Real'), ('baz1', 'Real'), ('baz2', 'Real'), ('gag1', 'Int'), ('gag2', 'Int'), ('gal1', 'Int'), ('gal2', 'Int'), ('bal1', 'Int'), ('bal2', 'Int'), ('bit1', 'bool'), ('bit2', 'bool'), ('gut1', 'bool'), ('gut2', 'bool'), ('gat1', 'bool'), ('gat2', 'bool')]
    OrderedDict()
    OrderedDict()
    OrderedDict()
    """
    scalars  = []

    fst_dims = []
    snd_dims = []
    trd_dims = []

    for name, argtype, _, dims in arg_data:
        if dims is not None:
            expect(len(dims) >= 1 and len(dims) <= 3,
                   f"Only 1d-3d data is supported, {name} has too many dims: {len(dims)}")

            if dims[0] not in fst_dims:
                fst_dims.append(dims[0])
            if len(dims) > 1 and dims[1] not in snd_dims:
                snd_dims.append(dims[1])
            if len(dims) > 2 and dims[2] not in trd_dims:
                trd_dims.append(dims[2])

    all_dims = list(OrderedDict([(item, None) for item in (fst_dims + snd_dims + trd_dims)]))
    real_data = OrderedDict()
    int_data = OrderedDict()
    bool_data = OrderedDict()

    for name, argtype, intent, dims in arg_data:
        if filter_out_intent is None or intent != filter_out_intent:
            if dims is None:
                if not (is_custom_type(argtype) and filter_scalar_custom_types):
                    if name not in all_dims:
                        scalars.append( (name, get_cxx_scalar_type(argtype)))
                    else:
                        expect(argtype == "integer", f"Expected dimension {name} to be of type integer")

            elif argtype == "integer":
                int_data.setdefault(dims, []).append(name)

            elif argtype == "real":
                real_data.setdefault(dims, []).append(name)

            elif argtype == "logical":
                bool_data.setdefault(dims, []).append(name)

    return fst_dims, snd_dims, trd_dims, all_dims, scalars, real_data, int_data, bool_data

###############################################################################
def gen_struct_api(physics, struct_name, arg_data):
###############################################################################
    r"""
    Given data, generate code for data struct api

    >>> print("\n".join(gen_struct_api("shoc", "DataSubName", UT_ARG_DATA)))
    DataSubName(Int shcol_, Int nlev_, Int nlevi_, Int ntracers_, Real gag_, Int bab1_, Int bab2_, bool val_) :
      PhysicsTestData({{ shcol_ }, { shcol_, nlev_ }, { shcol_, nlevi_ }, { shcol_, nlev_, ntracers_ }, { shcol_ }, { shcol_ }}, {{ &foo1, &foo2, &baz }, { &bar1, &bar2 }, { &bak1, &bak2 }, { &tracerd1, &tracerd2 }}, {{ &bag, &ball1, &ball2 }}, {{ &vals }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_), ntracers(ntracers_), gag(gag_), bab1(bab1_), bab2(bab2_), val(val_) {}
    <BLANKLINE>
    PTD_STD_DEF(DataSubName, 8, shcol, nlev, nlevi, ntracers, gag, bab1, bab2, val);
    """
    _, _, _, all_dims, scalars, real_data, int_data, bool_data = group_data(arg_data, filter_scalar_custom_types=True)

    result = []
    dim_args = [(item, "Int") for item in all_dims if item is not None]
    cons_args = dim_args + scalars
    result.append("{struct_name}({cons_args}) :".\
                  format(struct_name=struct_name,
                         cons_args=", ".join(["{} {}_".format(argtype, name) for name, argtype in cons_args])))

    dim_cxx_vec = []
    real_vec = []
    int_vec = []
    bool_vec = []
    for data, data_vec in zip([real_data, int_data, bool_data], [real_vec, int_vec, bool_vec]):
        for dims, items in data.items():
            dim_cxx_vec.append(f"{{ {', '.join(['{}_'.format(item) for item in dims])} }}")
            data_vec.append(f"{{ {', '.join(['&{}'.format(item) for item in items])} }}")

    parent_call = f"  PhysicsTestData({{{', '.join(dim_cxx_vec)}}}, {{{', '.join(real_vec)}}}"
    if int_vec or bool_vec:
        parent_call += f", {{{', '.join(int_vec)}}}"
    if bool_vec:
        parent_call += f", {{{', '.join(bool_vec)}}}"
    parent_call += ")"

    parent_call += f", {', '.join(['{0}({0}_)'.format(name) for name, _ in cons_args])}"

    parent_call += " {}"
    result.append(parent_call)
    result.append("")

    result.append("PTD_STD_DEF({}, {}, {});".\
                  format(struct_name, len(cons_args), ", ".join([name for name, _ in cons_args])))

    return result

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
                   f"Found multiple begin matches for pattern '{begin_regex.pattern}' before end pattern '{end_regex.pattern}' was found")

            begin_idx = idx

        if end_match and begin_idx is not None:
            end_idx = idx
            break

    if begin_idx is not None:
        expect(end_idx is not None,
               "Found no ending match for begin pattern '{begin_regex.pattern}' starting on line {begin_idx} and searching end pattern '{end_regex.pattern}'")

    return None if begin_idx is None else (begin_idx, end_idx+1)

###############################################################################
def get_data_by_name(arg_data, arg_name, data_idx):
###############################################################################
    for name, a, b, c in arg_data:
        if name == arg_name:
            return [name, a, b, c][data_idx]

    expect(False, f"Name {arg_name} not found")

###############################################################################
def get_rank_map(arg_data, arg_names):
###############################################################################
    # Create map of rank -> [args]
    rank_map = {}
    for arg in arg_names:
        dims = get_data_by_name(arg_data, arg, ARG_DIMS)
        rank = len(dims)
        if rank in rank_map:
            rank_map[rank].append(arg)
        else:
            rank_map[rank] = [arg]

    return rank_map

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
        expect(normalized_source_repo is not None, f"source repo {source_repo} is not a valid repo")
        source_repo = normalized_source_repo

        normalized_target_repo = get_git_toplevel_dir(repo=target_repo)
        expect(normalized_target_repo is not None, f"target repo {target_repo} is not a valid repo")
        target_repo = normalized_target_repo

        # configuration
        self._subs        = subs
        self._pieces      = pieces
        self._physics     = physics
        self._overwrite   = overwrite
        self._kernel      = kernel
        self._source_repo = Path(source_repo).resolve()
        self._target_repo = Path(target_repo).resolve()
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
        if phys not in self._db:
            origin_files = get_physics_data(phys, ORIGIN_FILES)
            self._db[phys] = {}
            for origin_file in origin_files:
                origin_file = self._source_repo / origin_file
                expect(origin_file.exists(), f"Missing origin file for physics {phys}: {origin_file}")
                db = parse_origin(origin_file.open(encoding="utf-8").read(), self._subs)
                self._db[phys].update(db)
                if self._verbose:
                    print(f"For physics {phys}, found:")
                    for sub in self._subs:
                        if sub in db:
                            print(f"  For subroutine {sub}, found args:")
                            for name, argtype, intent, dims in db[sub]:
                                print("    name:{} type:{} intent:{} dims:({})".\
                                      format(name, argtype, intent, ",".join(dims) if dims else "scalar"))

        return self._db[phys]

    ###########################################################################
    def _get_arg_data(self, phys, sub):
    ###########################################################################
        phys_db = self._get_db(phys)
        expect(sub in phys_db, f"No data for subroutine {sub} in physics {phys}")
        return phys_db[sub]

    ###########################################################################
    def dry_run(self):
    ###########################################################################
        return self._dry_run

    ###############################################################################
    def get_path_for_piece_file(self, physics, sub, piece):
    ###############################################################################
        root_dir = Path(get_physics_data(physics, CXX_ROOT))
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
          subroutine fake_sub_c(foo1, foo2, bar1, bar2, bak1, bak2, tracerd1, tracerd2, gag, baz, bag, bab1, bab2, val, vals, shcol, nlev, nlevi, ntracers, ball1, ball2) bind(C)
            use shoc, only : fake_sub
        <BLANKLINE>
            real(kind=c_real) , intent(in), dimension(shcol) :: foo1, foo2
            real(kind=c_real) , intent(in), dimension(shcol, nlev) :: bar1, bar2
            real(kind=c_real) , intent(in), dimension(shcol, nlevi) :: bak1, bak2
            real(kind=c_real) , intent(in), dimension(shcol, nlev, ntracers) :: tracerd1, tracerd2
            real(kind=c_real) , value, intent(in) :: gag
            real(kind=c_real) , intent(inout), dimension(shcol) :: baz
            integer(kind=c_int) , intent(in), dimension(shcol) :: bag
            integer(kind=c_int) , intent(out) :: bab1, bab2
            logical(kind=c_bool) , value, intent(in) :: val
            logical(kind=c_bool) , intent(in), dimension(shcol) :: vals
            integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi, ntracers
            integer(kind=c_int) , intent(out), dimension(shcol) :: ball1, ball2
        <BLANKLINE>
            call fake_sub(foo1, foo2, bar1, bar2, bak1, bak2, tracerd1, tracerd2, gag, baz, bag, bab1, bab2, val, vals, shcol, nlev, nlevi, ntracers, ball1, ball2)
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
          subroutine fake_sub_f(foo1, foo2, bar1, bar2, bak1, bak2, tracerd1, tracerd2, gag, baz, bag, bab1, bab2, val, vals, shcol, nlev, nlevi, ntracers, ball1, ball2) bind(C)
            use iso_c_binding
        <BLANKLINE>
            real(kind=c_real) , intent(in), dimension(shcol) :: foo1, foo2
            real(kind=c_real) , intent(in), dimension(shcol, nlev) :: bar1, bar2
            real(kind=c_real) , intent(in), dimension(shcol, nlevi) :: bak1, bak2
            real(kind=c_real) , intent(in), dimension(shcol, nlev, ntracers) :: tracerd1, tracerd2
            real(kind=c_real) , value, intent(in) :: gag
            real(kind=c_real) , intent(inout), dimension(shcol) :: baz
            integer(kind=c_int) , intent(in), dimension(shcol) :: bag
            integer(kind=c_int) , intent(out) :: bab1, bab2
            logical(kind=c_bool) , value, intent(in) :: val
            logical(kind=c_bool) , intent(in), dimension(shcol) :: vals
            integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi, ntracers
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
        void fake_sub_c(Real* foo1, Real* foo2, Real* bar1, Real* bar2, Real* bak1, Real* bak2, Real* tracerd1, Real* tracerd2, Real gag, Real* baz, Int* bag, Int* bab1, Int* bab2, bool val, bool* vals, Int shcol, Int nlev, Int nlevi, Int ntracers, Int* ball1, Int* ball2);
        <BLANKLINE>
        """
        arg_data = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        arg_decls = gen_arg_cxx_decls(arg_data)
        result = f"void {sub}_c({', '.join(arg_decls)});\n"
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
        result = f"void {sub}({struct_name}& d);"
        return result

    ###########################################################################
    def gen_cxx_c2f_glue_impl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_c2f_glue_impl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
        void fake_sub(FakeSubData& d)
        {
          shoc_init(d.nlev, true);
          d.transpose<ekat::TransposeDirection::c2f>();
          fake_sub_c(d.foo1, d.foo2, d.bar1, d.bar2, d.bak1, d.bak2, d.tracerd1, d.tracerd2, d.gag, d.baz, d.bag, &d.bab1, &d.bab2, d.val, d.vals, d.shcol, d.nlev, d.nlevi, d.ntracers, d.ball1, d.ball2);
          d.transpose<ekat::TransposeDirection::f2c>();
        }
        <BLANKLINE>
        <BLANKLINE>
        """
        arg_data         = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        arg_data_args    = ", ".join(gen_cxx_data_args(phys, arg_data))
        need_transpose   = needs_transpose(arg_data)
        transpose_code_1 = "\n  d.transpose<ekat::TransposeDirection::c2f>();" if need_transpose else ""
        transpose_code_2 = "\n  d.transpose<ekat::TransposeDirection::f2c>();" if need_transpose else ""
        data_struct      = get_data_struct_name(sub)
        init_code        = get_physics_data(phys, INIT_CODE)

        result = \
f"""void {sub}({data_struct}& d)
{{
  {init_code}{transpose_code_1}
  {sub}_c({arg_data_args});{transpose_code_2}
}}

"""
        return result

    ###########################################################################
    def gen_cxx_c2f_data(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_c2f_data("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
        struct FakeSubData : public PhysicsTestData {
          // Inputs
          Real *foo1, *foo2, *bar1, *bar2, *bak1, *bak2, *tracerd1, *tracerd2;
          Real gag;
          Int *bag;
          bool val;
          bool *vals;
          Int shcol, nlev, nlevi, ntracers;
        <BLANKLINE>
          // Inputs/Outputs
          Real *baz;
        <BLANKLINE>
          // Outputs
          Int bab1, bab2;
          Int *ball1, *ball2;
        <BLANKLINE>
          FakeSubData(Int shcol_, Int nlev_, Int nlevi_, Int ntracers_, Real gag_, Int bab1_, Int bab2_, bool val_) :
            PhysicsTestData({{ shcol_ }, { shcol_, nlev_ }, { shcol_, nlevi_ }, { shcol_, nlev_, ntracers_ }, { shcol_ }, { shcol_ }}, {{ &foo1, &foo2, &baz }, { &bar1, &bar2 }, { &bak1, &bak2 }, { &tracerd1, &tracerd2 }}, {{ &bag, &ball1, &ball2 }}, {{ &vals }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_), ntracers(ntracers_), gag(gag_), bab1(bab1_), bab2(bab2_), val(val_) {}
        <BLANKLINE>
          PTD_STD_DEF(FakeSubData, 8, shcol, nlev, nlevi, ntracers, gag, bab1, bab2, val);
        };
        <BLANKLINE>
        <BLANKLINE>
        """
        arg_data         = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        struct_members   = "\n  ".join(gen_struct_members(arg_data))
        any_arrays       = has_arrays(arg_data)
        struct_name      = get_data_struct_name(sub)
        inheritance      = " : public PhysicsTestData" if any_arrays else ""
        api              = "\n  " + "\n  ".join(gen_struct_api(phys, struct_name, arg_data) if any_arrays else "")

        result = \
f"""struct {struct_name}{inheritance} {{
  {struct_members}{api}
}};

"""
        return result

    ###########################################################################
    def gen_cxx_f2c_bind_decl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_f2c_bind_decl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
        void fake_sub_f(Real* foo1, Real* foo2, Real* bar1, Real* bar2, Real* bak1, Real* bak2, Real* tracerd1, Real* tracerd2, Real gag, Real* baz, Int* bag, Int* bab1, Int* bab2, bool val, bool* vals, Int shcol, Int nlev, Int nlevi, Int ntracers, Int* ball1, Int* ball2);
        """
        arg_data  = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        arg_decls = gen_arg_cxx_decls(arg_data)

        return f"void {sub}_host({', '.join(arg_decls)});"

    ###########################################################################
    def gen_cxx_f2c_bind_impl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_f2c_bind_impl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
        void fake_sub_f(Real* foo1, Real* foo2, Real* bar1, Real* bar2, Real* bak1, Real* bak2, Real* tracerd1, Real* tracerd2, Real gag, Real* baz, Int* bag, Int* bab1, Int* bab2, bool val, bool* vals, Int shcol, Int nlev, Int nlevi, Int ntracers, Int* ball1, Int* ball2)
        {
          // TODO
        }
        <BLANKLINE>
        >>> print(gb.gen_cxx_f2c_bind_impl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA_ALL_SCALAR))
        void fake_sub_f(Real foo1, Real foo2, Real* bar1, Real* bar2, Real* baz1, Real* baz2, Int gag1, Int gag2, Int* gal1, Int* gal2, Int* bal1, Int* bal2, bool bit1, bool bit2, bool* gut1, bool* gut2, bool* gat1, bool* gat2)
        {
        #if 0
          using PF = Functions<Real, DefaultDevice>;
        <BLANKLINE>
          using Spack   = typename PF::Spack;
          using view_1d = typename PF::view_1d<Real>;
          using iview_1d = typename PF::view_1d<Int>;
          using bview_1d = typename PF::view_1d<bool>;
        <BLANKLINE>
          view_1d t_d("t_d", 4);
          const auto t_h = Kokkos::create_mirror_view(t_d);
        <BLANKLINE>
          iview_1d it_d("it_d", 4);
          const auto it_h = Kokkos::create_mirror_view(it_d);
        <BLANKLINE>
          bview_1d bt_d("bt_d", 4);
          const auto bt_h = Kokkos::create_mirror_view(bt_d);
        <BLANKLINE>
          Real local_bar1(*bar1), local_bar2(*bar2);
          Int local_gal1(*gal1), local_gal2(*gal2);
          bool local_gut1(*gut1), local_gut2(*gut2);
          Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
            Spack bar1_(local_bar1), bar2_(local_bar2), foo1_(foo1), foo2_(foo2), baz1_(), baz2_();
            Int bal1_(), bal2_(), gal1_(local_gal1), gal2_(local_gal2);
            bool gat1_(), gat2_(), gut1_(local_gut1), gut2_(local_gut2);
            PF::fake_sub(foo1_, foo2_, bar1_, bar2_, baz1_, baz2_, gag1, gag2, gal1_, gal2_, bal1_, bal2_, bit1, bit2, gut1_, gut2_, gat1_, gat2_);
            t_d(0) = bar1_[0];
            t_d(1) = bar2_[0];
            t_d(2) = baz1_[0];
            t_d(3) = baz2_[0];
            it_d(0) = bal1_;
            it_d(1) = bal2_;
            it_d(2) = gal1_;
            it_d(3) = gal2_;
            bt_d(0) = gat1_;
            bt_d(1) = gat2_;
            bt_d(2) = gut1_;
            bt_d(3) = gut2_;
          });
          Kokkos::deep_copy(t_h, t_d);
          Kokkos::deep_copy(it_h, it_d);
          Kokkos::deep_copy(bt_h, bt_d);
          *bar1 = t_h(0);
          *bar2 = t_h(1);
          *baz1 = t_h(2);
          *baz2 = t_h(3);
          *bal1 = it_h(0);
          *bal2 = it_h(1);
          *gal1 = it_h(2);
          *gal2 = it_h(3);
          *gat1 = bt_h(0);
          *gat2 = bt_h(1);
          *gut1 = bt_h(2);
          *gut2 = bt_h(3);
        #endif
        <BLANKLINE>
        }
        <BLANKLINE>
        """
        arg_data  = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        arg_names = [item[ARG_NAME] for item in arg_data]
        decl      = self.gen_cxx_f2c_bind_decl(phys, sub, force_arg_data=force_arg_data).rstrip(";")

        impl = ""
        if has_arrays(arg_data):
            #
            # Steps:
            # 1) Set up typedefs
            # 2) Sync to device
            # 3) Unpack view array
            # 4) Get nk_pack and policy
            # 5) Get subviews
            # 6) Call fn
            # 7) Sync back to host
            #
            inputs, inouts, outputs = split_by_intent(arg_data)
            reals, ints, bools   = split_by_type(arg_data)
            scalars, views = split_by_scalar_vs_view(arg_data)
            all_inputs = inputs + inouts
            all_outputs = inouts + outputs

            vreals = list(sorted(set(reals) & set(views)))
            vints  = list(sorted(set(ints)  & set(views)))
            vbools = list(sorted(set(bools) & set(views)))

            sreals = list(sorted(set(reals) & set(scalars)))
            sints  = list(sorted(set(ints)  & set(scalars)))
            sbools = list(sorted(set(bools) & set(scalars)))

            ivreals = list(sorted(set(vreals) & set(all_inputs)))
            ivints  = list(sorted(set(vints)  & set(all_inputs)))
            ivbools = list(sorted(set(vbools) & set(all_inputs)))

            ovreals = list(sorted(set(vreals) & set(all_outputs)))
            ovints  = list(sorted(set(vints)  & set(all_outputs)))
            ovbools = list(sorted(set(vbools) & set(all_outputs)))

            isreals = list(sorted(set(sreals) & set(all_inputs)))
            isints  = list(sorted(set(sints)  & set(all_inputs)))
            isbools = list(sorted(set(sbools) & set(all_inputs)))

            osreals = list(sorted(set(sreals) & set(all_outputs)))
            osints  = list(sorted(set(sints)  & set(all_outputs)))
            osbools = list(sorted(set(sbools) & set(all_outputs)))

            #
            # 1) Set up typedefs
            #

            # set up basics
            impl += "#if 0\n" # There's no way to guarantee this code compiles
            impl += "  using SHF        = Functions<Real, DefaultDevice>;\n"
            impl += "  using Scalar     = typename SHF::Scalar;\n"
            impl += "  using Spack      = typename SHF::Spack;\n"
            impl += "  using KT         = typename SHF::KT;\n"
            impl += "  using ExeSpace   = typename KT::ExeSpace;\n"
            impl += "  using MemberType = typename SHF::MemberType;\n\n"

            prefix_list  = ["", "i", "b"]
            type_list    = ["Real", "Int", "bool"]
            ktype_list   = ["Spack", "Int", "bool"]

            # make necessary view types. Anything that's an array needs a view type
            for view_group, prefix_char, typename in zip([vreals, vints, vbools], prefix_list, type_list):
                if view_group:
                    rank_map = get_rank_map(arg_data, view_group)
                    for rank in rank_map:
                        if typename == "Real" and rank > 1:
                            # Probably this should be packed data
                            impl += f"  using {prefix_char}view_{rank}d = typename SHF::view_{rank}d<Spack>;\n"
                        else:
                            impl += f"  using {prefix_char}view_{rank}d = typename SHF::view_{rank}d<{typename}>;\n"

            impl += "\n"

            #
            # 2) Sync to device. Do ALL views, not just inputs
            #

            for input_group, prefix_char, typename in zip([vreals, vints, vbools], prefix_list, type_list):
                if input_group:
                    rank_map = get_rank_map(arg_data, input_group)

                    for rank, arg_list in rank_map.items():
                        impl += f"  static constexpr Int {prefix_char}num_arrays_{rank} = {len(arg_list)};\n"
                        impl += f"  std::vector<{prefix_char}view_{rank}d> {prefix_char}temp_d_{rank}({prefix_char}num_arrays_{rank});\n"
                        for rank_itr in range(rank):
                            dims = [get_data_by_name(arg_data, arg_name, ARG_DIMS)[rank_itr] for arg_name in arg_list]
                            impl += f"  std::vector<int> {prefix_char}dim_{rank}_{rank_itr}_sizes = {{{', '.join(dims)}}};\n"
                        dim_vectors = [f"{prefix_char}dim_{rank}_{rank_itr}_sizes" for rank_itr in range(rank)]
                        funcname = "ekat::host_to_device" if (typename == "Real" and rank > 1) else "ScreamDeepCopy::copy_to_device"
                        impl += f"  {funcname}({{{', '.join(arg_list)}}}, {', '.join(dim_vectors)}, {prefix_char}temp_d_{rank});\n\n"

            #
            # 3) Unpack view array
            #

            for input_group, prefix_char, typename in zip([vreals, vints, vbools], prefix_list, type_list):
                if input_group:
                    rank_map = get_rank_map(arg_data, input_group)

                    for rank, arg_list in rank_map.items():
                        impl += f"  {prefix_char}view_{rank}d\n"
                        for idx, input_item in enumerate(arg_list):
                            impl += f"    {input_item}_d({prefix_char}temp_d_{rank}[{idx}]){';' if idx == len(arg_list) - 1 else ','}\n"
                        impl += "\n"


            #
            # 4) Get nk_pack and policy, launch kernel
            #
            impl += "  const Int nk_pack = ekat::npack<Spack>(nlev);\n"
            impl += "  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);\n"
            impl += "  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {\n"
            impl += "    const Int i = team.league_rank();\n\n"

            #
            # 5) Get subviews
            #
            for view_group, prefix_char, typename in zip([vreals, vints, vbools], prefix_list, type_list):
                if view_group:
                    for view_arg in view_group:
                        dims = get_data_by_name(arg_data, view_arg, ARG_DIMS)
                        if "shcol" in dims:
                            if len(dims) == 1:
                                impl += f"    const Scalar {view_arg}_s = {view_arg}_d(i);\n"
                            else:
                                impl += f"    const auto {view_arg}_s = ekat::subview({view_arg}_d, i);\n"

                    impl += "\n"

            #
            # 6) Call fn
            #
            kernel_arg_names = []
            for arg_name in arg_names:
                if arg_name in views:
                    if "shcol" in dims:
                        kernel_arg_names.append(f"{arg_name}_s")
                    else:
                        kernel_arg_names.append(f"{arg_name}_d")
                else:
                    kernel_arg_names.append(arg_name)

            impl += f"    SHF::{sub}({', '.join(kernel_arg_names)});\n"
            impl +=  "  });\n"

            #
            # 7) Sync back to host
            #
            for output_group, prefix_char, typename in zip([ovreals, ovints, ovbools], prefix_list, type_list):
                if output_group:
                    rank_map = get_rank_map(arg_data, output_group)

                    for rank, arg_list in rank_map.items():
                        impl += f"  std::vector<{prefix_char}view_{rank}d> {prefix_char}tempout_d_{rank}({prefix_char}num_arrays_{rank});\n"
                        for rank_itr in range(rank):
                            dims = [get_data_by_name(arg_data, arg_name, ARG_DIMS)[rank_itr] for arg_name in arg_list]
                            impl += f"  std::vector<int> {prefix_char}dim_{rank}_{rank_itr}_out_sizes = {{{', '.join(dims)}}};\n"
                        dim_vectors = [f"{prefix_char}dim_{rank}_{rank_itr}_out_sizes" for rank_itr in range(rank)]
                        funcname = "ekat::device_to_host" if (typename == "Real" and rank > 1) else "ScreamDeepCopy::copy_to_host"
                        impl += f"  {funcname}({{{', '.join(arg_list)}}}, {', '.join(dim_vectors)}, {prefix_char}tempout_d_{rank});\n\n"

            impl += "#endif\n"

        else:
            inputs, inouts, outputs = split_by_intent(arg_data)
            reals, ints, logicals   = split_by_type(arg_data)
            all_inputs = inputs + inouts
            all_outputs = inouts + outputs

            oreals = list(sorted(set(all_outputs) & set(reals)))
            oints  = list(sorted(set(all_outputs) & set(ints)))
            obools = list(sorted(set(all_outputs) & set(logicals)))

            ireals = list(sorted(set(all_inputs) & set(reals)))

            ioreals = list(sorted(set(inouts) & set(reals)))
            ioints  = list(sorted(set(inouts) & set(ints)))
            iobools = list(sorted(set(inouts) & set(logicals)))

            ooreals = list(sorted(set(outputs) & set(reals)))

            # set up basics
            impl += "#if 0\n" # There's no way to guarantee this code compiles
            impl += "  using PF = Functions<Real, DefaultDevice>;\n"
            impl += "\n"
            impl += "  using Spack   = typename PF::Spack;\n"

            prefix_list  = ["", "i", "b"]
            type_list    = ["Real", "Int", "bool"]
            ktype_list   = ["Spack", "Int", "bool"]

            # make necessary view types
            for output_group, prefix_char, typename in zip([oreals, oints, obools], prefix_list, type_list):
                if output_group:
                    impl += f"  using {prefix_char}view_1d = typename PF::view_1d<{typename}>;\n"

            impl += "\n"

            # make output views for host and device
            for output_group, prefix_char in zip([oreals, oints, obools], prefix_list):
                if output_group:
                    impl += f'  {prefix_char}view_1d {prefix_char}t_d("{prefix_char}t_d", {len(output_group)});\n'
                    impl += f"  const auto {prefix_char}t_h = Kokkos::create_mirror_view({prefix_char}t_d);\n"
                    impl += "\n"

            # inout data must be derefenced before the kernel
            for io_group, typename in zip([ioreals, ioints, iobools], type_list):
                if io_group:
                    impl += f"  {typename} {', '.join(['local_{0}(*{0})'.format(item) for item in io_group])};\n"

            # start a kernel
            impl += "  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {\n"

            # Declare output temporaries for non-reals. We need temporaries for all reals, including
            # inputs to convert them to Spacks. This code will need to be manually corrected for reals that should
            # not be packed (like dt)
            for output_group, typename in zip([list(ireals) + list(ooreals), oints, obools], ktype_list):
                if output_group:
                    impl += f"    {typename} "
                    temp_cons = []
                    for item in output_group:
                        if item in inouts:
                            temp_cons.append(f"{item}_(local_{item})")
                        elif item in outputs:
                            temp_cons.append(f"{item}_()")
                        else:
                            temp_cons.append(f"{item}_({item})")

                    impl += f"{', '.join(temp_cons)};\n"

            # Make cxx call
            kernel_arg_names = []
            for arg_name in arg_names:
                if arg_name in inputs and arg_name not in reals:
                    kernel_arg_names.append(arg_name)
                else:
                    kernel_arg_names.append(arg_name + "_")

            impl += f"    PF::{sub}({', '.join(kernel_arg_names)});\n"

            # Load output data into views
            for output_group, prefix_char in zip([oreals, oints, obools], prefix_list):
                for idx, item in enumerate(output_group):
                    if output_group == oreals:
                        impl += f"    {prefix_char}t_d({idx}) = {item}_[0];\n"
                    else:
                        impl += f"    {prefix_char}t_d({idx}) = {item}_;\n"

            # finish kernel
            impl += "  });\n"

            # copy outputs back to host
            for output_group, prefix_char in zip([oreals, oints, obools], prefix_list):
                if output_group:
                    impl += f"  Kokkos::deep_copy({prefix_char}t_h, {prefix_char}t_d);\n"

            # copy from views into pointer args
            for output_group, prefix_char in zip([oreals, oints, obools], prefix_list):
                for idx, item in enumerate(output_group):
                    impl += f"  *{item} = {prefix_char}t_h({idx});\n"

            impl += "#endif\n"

        result = \
f"""{decl}
{{
{impl}
}}
"""
        return result

    ###########################################################################
    def gen_cxx_func_decl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_func_decl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
          KOKKOS_FUNCTION
          static void fake_sub(const uview_1d<const Spack>& foo1, const uview_1d<const Spack>& foo2, const uview_1d<const Spack>& bar1, const uview_1d<const Spack>& bar2, const uview_1d<const Spack>& bak1, const uview_1d<const Spack>& bak2, const uview_1d<const Spack>& tracerd1, const uview_1d<const Spack>& tracerd2, const Spack& gag, const uview_1d<Spack>& baz, const uview_1d<const Int>& bag, Int& bab1, Int& bab2, const bool& val, const uview_1d<const bool>& vals, const Int& shcol, const Int& nlev, const Int& nlevi, const Int& ntracers, const uview_1d<Int>& ball1, const uview_1d<Int>& ball2);
        """
        arg_data = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        arg_decls = gen_arg_cxx_decls(arg_data, kokkos=True)

        return f"  KOKKOS_FUNCTION\n  static void {sub}({', '.join(arg_decls)});"

    ###########################################################################
    def gen_cxx_incl_impl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_incl_impl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
        # include "impl/shoc_fake_sub_impl.hpp"
        """
        impl_path = get_piece_data(phys, sub, "cxx_func_impl", FILEPATH, self)
        return f'# include "{impl_path}"'

    ###########################################################################
    def gen_cxx_func_impl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_func_impl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
        KOKKOS_FUNCTION
        void Functions<S,D>::fake_sub(const uview_1d<const Spack>& foo1, const uview_1d<const Spack>& foo2, const uview_1d<const Spack>& bar1, const uview_1d<const Spack>& bar2, const uview_1d<const Spack>& bak1, const uview_1d<const Spack>& bak2, const uview_1d<const Spack>& tracerd1, const uview_1d<const Spack>& tracerd2, const Spack& gag, const uview_1d<Spack>& baz, const uview_1d<const Int>& bag, Int& bab1, Int& bab2, const bool& val, const uview_1d<const bool>& vals, const Int& shcol, const Int& nlev, const Int& nlevi, const Int& ntracers, const uview_1d<Int>& ball1, const uview_1d<Int>& ball2)
        {
          // TODO
          // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
        }
        """
        decl = self.gen_cxx_func_decl(phys, sub, force_arg_data=force_arg_data).rstrip(";").replace("void ", "void Functions<S,D>::").replace("static ", "")
        decl = "\n".join(line.strip() for line in decl.splitlines())

        # I don't think any intelligent guess at an impl is possible here
        result = \
f"""{decl}
{{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}}"""
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
        return f"    struct {test_struct};"

    ###########################################################################
    def gen_cxx_bfb_unit_impl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_bfb_unit_impl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
          static void run_bfb()
          {
            auto engine = setup_random_test();
        <BLANKLINE>
            FakeSubData f90_data[] = {
              // TODO
            };
        <BLANKLINE>
            static constexpr Int num_runs = sizeof(f90_data) / sizeof(FakeSubData);
        <BLANKLINE>
            // Generate random input data
            // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
            for (auto& d : f90_data) {
              d.randomize(engine);
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
              fake_sub_f(d.foo1, d.foo2, d.bar1, d.bar2, d.bak1, d.bak2, d.tracerd1, d.tracerd2, d.gag, d.baz, d.bag, &d.bab1, &d.bab2, d.val, d.vals, d.shcol, d.nlev, d.nlevi, d.ntracers, d.ball1, d.ball2);
              d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
            }
        <BLANKLINE>
            // Verify BFB results, all data should be in C layout
            if (SCREAM_BFB_TESTING) {
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
            }
          } // run_bfb
        >>> print(gb.gen_cxx_bfb_unit_impl("shoc", "fake_sub", force_arg_data=UT_ARG_DATA_ALL_SCALAR))
          static void run_bfb()
          {
            auto engine = setup_random_test();
        <BLANKLINE>
            FakeSubData f90_data[max_pack_size] = {
              // TODO
            };
        <BLANKLINE>
            static constexpr Int num_runs = sizeof(f90_data) / sizeof(FakeSubData);
        <BLANKLINE>
            // Generate random input data
            // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
            for (auto& d : f90_data) {
              d.randomize(engine);
            }
        <BLANKLINE>
            // Create copies of data for use by cxx and sync it to device. Needs to happen before fortran calls so that
            // inout data is in original state
            view_1d<FakeSubData> cxx_device("cxx_device", max_pack_size);
            const auto cxx_host = Kokkos::create_mirror_view(cxx_device);
            std::copy(&f90_data[0], &f90_data[0] + max_pack_size, cxx_host.data());
            Kokkos::deep_copy(cxx_device, cxx_host);
        <BLANKLINE>
            // Get data from fortran
            for (auto& d : f90_data) {
              fake_sub(d);
            }
        <BLANKLINE>
            // Get data from cxx. Run fake_sub from a kernel and copy results back to host
            Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
              const Int offset = i * Spack::n;
        <BLANKLINE>
              // Init pack inputs
              Spack bar1, bar2, foo1, foo2;
              for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
                bar1[s] = cxx_device(vs).bar1;
                bar2[s] = cxx_device(vs).bar2;
                foo1[s] = cxx_device(vs).foo1;
                foo2[s] = cxx_device(vs).foo2;
              }
        <BLANKLINE>
              // Init outputs
              Spack baz1(0), baz2(0);
        <BLANKLINE>
        <BLANKLINE>
              Functions::fake_sub(foo1, foo2, bar1, bar2, baz1, baz2, cxx_device(0).gag1, cxx_device(0).gag2, cxx_device(0).gal1, cxx_device(0).gal2, cxx_device(0).bal1, cxx_device(0).bal2, cxx_device(0).bit1, cxx_device(0).bit2, cxx_device(0).gut1, cxx_device(0).gut2, cxx_device(0).gat1, cxx_device(0).gat2);
        <BLANKLINE>
              // Copy spacks back into cxx_device view
              for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
                cxx_device(vs).bar1 = bar1[s];
                cxx_device(vs).bar2 = bar2[s];
                cxx_device(vs).baz1 = baz1[s];
                cxx_device(vs).baz2 = baz2[s];
              }
        <BLANKLINE>
            });
        <BLANKLINE>
            Kokkos::deep_copy(cxx_host, cxx_device);
        <BLANKLINE>
            // Verify BFB results
            if (SCREAM_BFB_TESTING) {
              for (Int i = 0; i < num_runs; ++i) {
                FakeSubData& d_f90 = f90_data[i];
                FakeSubData& d_cxx = cxx_host[i];
                REQUIRE(d_f90.bar1 == d_cxx.bar1);
                REQUIRE(d_f90.bar2 == d_cxx.bar2);
                REQUIRE(d_f90.baz1 == d_cxx.baz1);
                REQUIRE(d_f90.baz2 == d_cxx.baz2);
                REQUIRE(d_f90.gal1 == d_cxx.gal1);
                REQUIRE(d_f90.gal2 == d_cxx.gal2);
                REQUIRE(d_f90.bal1 == d_cxx.bal1);
                REQUIRE(d_f90.bal2 == d_cxx.bal2);
                REQUIRE(d_f90.gut1 == d_cxx.gut1);
                REQUIRE(d_f90.gut2 == d_cxx.gut2);
                REQUIRE(d_f90.gat1 == d_cxx.gat1);
                REQUIRE(d_f90.gat2 == d_cxx.gat2);
              }
            }
          } // run_bfb
        """
        arg_data = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        data_struct = get_data_struct_name(sub)
        has_array = has_arrays(arg_data)
        need_transpose = needs_transpose(arg_data)
        arg_data_args    = ", ".join(gen_cxx_data_args(phys, arg_data))

        gen_random = \
"""

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize(engine);
    }"""

        _, _, _, _, scalars, real_data, int_data, bool_data = group_data(arg_data, filter_out_intent="in")
        check_scalars, check_arrays = "", ""
        for scalar in scalars:
            check_scalars += f"        REQUIRE(d_f90.{scalar[0]} == d_cxx.{scalar[0]});\n"

        if has_array:
            c2f_transpose_code = "" if not need_transpose else \
"""
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout"""
            f2c_transpose_code = "" if not need_transpose else \
"""
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout"""

            all_data = OrderedDict(real_data)
            for type_data in [int_data, bool_data]:
                for k, v in type_data.items():
                    if k in all_data:
                        all_data[k].extend(v)
                    else:
                        all_data[k] = v

            for _, data in all_data.items():
                check_arrays += f"        for (Int k = 0; k < d_f90.total(d_f90.{data[0]}); ++k) {{\n"
                for datum in data:
                    check_arrays += f"          REQUIRE(d_f90.total(d_f90.{data[0]}) == d_cxx.total(d_cxx.{datum}));\n"
                    check_arrays += f"          REQUIRE(d_f90.{datum}[k] == d_cxx.{datum}[k]);\n"

                check_arrays += "        }\n"

        if has_array:
            result = \
"""  static void run_bfb()
  {{
    auto engine = setup_random_test();

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
    if (SCREAM_BFB_TESTING) {{
      for (Int i = 0; i < num_runs; ++i) {{
        {data_struct}& d_f90 = f90_data[i];
        {data_struct}& d_cxx = cxx_data[i];
{check_scalars}{check_arrays}
      }}
    }}
  }} // run_bfb""".format(data_struct=data_struct,
                          sub=sub,
                          gen_random=gen_random,
                          c2f_transpose_code=c2f_transpose_code,
                          f2c_transpose_code=f2c_transpose_code,
                          arg_data_args=arg_data_args,
                          check_scalars=check_scalars,
                          check_arrays=check_arrays)
        else:
            inputs, inouts, outputs = split_by_intent(arg_data)
            reals                   = split_by_type(arg_data)[0]
            all_inputs              = inputs + inouts
            all_outputs             = inouts + outputs

            ireals  = list(sorted(set(all_inputs) & set(reals)))
            oreals  = list(sorted(set(all_outputs) & set(reals)))
            ooreals = list(sorted(set(outputs) & set(reals)))

            spack_init = ""
            if ireals:
                spack_init = \
"""// Init pack inputs
      Spack {ireals};
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {{
        {ireal_assigns}
      }}
""".format(ireals=", ".join(ireals), ireal_assigns="\n        ".join(["{0}[s] = cxx_device(vs).{0};".format(ireal) for ireal in ireals]))

            spack_output_init = ""
            if ooreals:
                spack_output_init = \
f"""// Init outputs
      Spack {', '.join(['{}(0)'.format(ooreal) for ooreal in ooreals])};
"""

            scalars = group_data(arg_data)[4]
            func_call = f"Functions::{sub}({', '.join([(scalar if scalar in reals else 'cxx_device(0).{}'.format(scalar)) for scalar, _ in scalars])});"

            spack_output_to_dview = ""
            if oreals:
                spack_output_to_dview = \
"""// Copy spacks back into cxx_device view
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {{
        {}
      }}
""".format("\n        ".join(["cxx_device(vs).{0} = {0}[s];".format(oreal) for oreal in oreals]))

            result = \
"""  static void run_bfb()
  {{
    auto engine = setup_random_test();

    {data_struct} f90_data[max_pack_size] = {{
      // TODO
    }};

    static constexpr Int num_runs = sizeof(f90_data) / sizeof({data_struct});{gen_random}

    // Create copies of data for use by cxx and sync it to device. Needs to happen before fortran calls so that
    // inout data is in original state
    view_1d<{data_struct}> cxx_device("cxx_device", max_pack_size);
    const auto cxx_host = Kokkos::create_mirror_view(cxx_device);
    std::copy(&f90_data[0], &f90_data[0] + max_pack_size, cxx_host.data());
    Kokkos::deep_copy(cxx_device, cxx_host);

    // Get data from fortran
    for (auto& d : f90_data) {{
      {sub}(d);
    }}

    // Get data from cxx. Run {sub} from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {{
      const Int offset = i * Spack::n;

      {spack_init}
      {spack_output_init}

      {func_call}

      {spack_output_to_dview}
    }});

    Kokkos::deep_copy(cxx_host, cxx_device);

    // Verify BFB results
    if (SCREAM_BFB_TESTING) {{
      for (Int i = 0; i < num_runs; ++i) {{
        {data_struct}& d_f90 = f90_data[i];
        {data_struct}& d_cxx = cxx_host[i];
{check_scalars}      }}
    }}
  }} // run_bfb""".format(data_struct=data_struct,
                          sub=sub,
                          gen_random=gen_random,
                          check_scalars=check_scalars,
                          spack_init=spack_init,
                          spack_output_init=spack_output_init,
                          func_call=func_call,
                          spack_output_to_dview=spack_output_to_dview)

        return result

    ###########################################################################
    def gen_cxx_eti(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cxx_eti("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
        #include "impl/shoc_fake_sub_impl.hpp"
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
f"""#include "{include_file}"

namespace scream {{
namespace {phys} {{

/*
 * Explicit instantiation for doing {sub} on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

}} // namespace {phys}
}} // namespace scream
"""

        return result

    ###########################################################################
    def gen_cmake_impl_eti(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cmake_impl_eti("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
            eti/shoc_fake_sub.cpp
        """
        eti_src = get_piece_data(phys, sub, "cxx_eti", FILEPATH, self)
        return f"    {eti_src}"

    ###########################################################################
    def gen_cmake_unit_test(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        >>> gb = GenBoiler([])
        >>> print(gb.gen_cmake_unit_test("shoc", "fake_sub", force_arg_data=UT_ARG_DATA))
            shoc_fake_sub_tests.cpp
        """
        test_src = Path(get_piece_data(phys, sub, "cxx_bfb_unit_impl", FILEPATH, self)).name
        return f"    {test_src}"

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
          shoc_init(d.nlev, true);
          d.transpose<ekat::TransposeDirection::c2f>();
          fake_sub_c(d.foo1, d.foo2, d.bar1, d.bar2, d.bak1, d.bak2, d.tracerd1, d.tracerd2, d.gag, d.baz, d.bag, &d.bab1, &d.bab2, d.val, d.vals, d.shcol, d.nlev, d.nlevi, d.ntracers, d.ball1, d.ball2);
          d.transpose<ekat::TransposeDirection::f2c>();
        }
        <BLANKLINE>

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
          shoc_init(d.nlev, true);
          d.transpose<ekat::TransposeDirection::c2f>();
          fake_sub_c(d.foo1, d.foo2, d.bar1, d.bar2, d.bak1, d.bak2, d.tracerd1, d.tracerd2, d.gag, d.baz, d.bag, &d.bab1, &d.bab2, d.val, d.vals, d.shcol, d.nlev, d.nlevi, d.ntracers, d.ball1, d.ball2);
          d.transpose<ekat::TransposeDirection::f2c>();
        }
        <BLANKLINE>

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
        void fake_sub_c(Real* foo1, Real* foo2, Real* bar1, Real* bar2, Real* bak1, Real* bak2, Real* tracerd1, Real* tracerd2, Real gag, Real* baz, Int* bag, Int* bab1, Int* bab2, bool val, bool* vals, Int shcol, Int nlev, Int nlevi, Int ntracers, Int* ball1, Int* ball2);
        """
        if force_arg_data is None: # don't want unit tests printing this
            print("===============================================================================")
            print(f"Trying to generate piece {piece} for subroutine {sub} for physics {phys}\n")

        base_filepath, was_filegen, insert_regex, self_begin_regex, self_end_regex, _ \
            = [item(phys, sub, self) for item in PIECES[piece]]

        filepath = self.get_path_for_piece_file(phys, sub, piece)

        if was_filegen:
            # If freshly generated file, we're done
            pass
        else:
            orig_lines = force_file_lines if force_file_lines else filepath.open(encoding="utf-8").read().splitlines()
            needs_rewrite = False
            gen_lines  = getattr(self, f"gen_{piece}")(phys, sub, force_arg_data=force_arg_data).splitlines()

            # Check to see if piece already exists
            try:
                existing_piece_line_range = check_existing_piece(orig_lines, self_begin_regex, self_end_regex)
            except SystemExit as e:
                expect(False, f"Problem parsing file {filepath} for existing piece {piece}: {e}")

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
                with filepath.open("w", encoding="utf-8") as fd:
                    fd.write("\n".join(orig_lines) + "\n")

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
                        print(f"Warning: failed to generate subroutine {sub} piece {piece} for physics {phys}, error: {e}")
                        all_success = False

        return all_success
