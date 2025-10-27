from utils import expect
from git_utils import get_git_toplevel_dir

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
#include "physics/{phys}/{phys}_functions.hpp"
#include "physics/{phys}/tests/infra/{phys}_test_data.hpp"

#include "{phys}_unit_tests_common.hpp"

#include <ekat_pack.hpp>
#include <ekat_team_policy_utils.hpp>

namespace scream {{
namespace {phys} {{
namespace unit_test {{

template <typename D>
struct UnitWrap::UnitTest<D>::{get_data_test_struct_name(sub)} : public UnitWrap::UnitTest<D>::Base {{

{gen_code}

}};

}} // namespace unit_test
}} // namespace {phys}
}} // namespace scream

namespace {{

TEST_CASE("{sub}_bfb", "[{phys}]")
{{
  using TestStruct = scream::{phys}::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::{get_data_test_struct_name(sub)};

  TestStruct t;
  t.run_bfb();
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

{gen_code}

}} // namespace {phys}
}} // namespace scream

#endif
""",

###############################################################################

    "cxx_eti": lambda phys, sub, gen_code: gen_code
}

# piece map. maps the name of a piece of boilerplate that needs to be genereate to:
#   filepath (relative to cxx_root))
#   file creation / file exists check
#   insert regex : A regex that, when matched, indicates the generated code should be inserted here
#   id self begin: A regex that, when matched, indicates the generated for the piece starts here
#   id self end: A regex that, when matched, indicates the generated for the piece ends here
#   desc: A description of this piece
FILEPATH, FILECREATE, INSERT_REGEX, ID_SELF_BEGIN_REGEX, ID_SELF_END_REGEX, DESC = range(6)
PIECES = dict([
    ("f90_c2f_bind", (
        lambda phys, sub, gb: f"tests/infra/{phys}_c2f_bridge.f90",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "f90_c2f_bind"),
        lambda phys, sub, gb: re.compile(fr"^\s*end\s+module\s{phys}_c2f_bridge"), # put at end of module
        lambda phys, sub, gb: get_subroutine_begin_regex(sub + "_bridge_f"), # sub_f begin
        lambda phys, sub, gb: get_subroutine_end_regex(sub + "_bridge_f"),    # sub_f end
        lambda *x           : "The c to f90 fortran subroutine implementation(<name>_bridge_f)"
    )),

    ("f90_f2c_bind"  , (
        lambda phys, sub, gb: f"{phys}_f2c_bridge.f90",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "f90_f2c_bind"),
        lambda phys, sub, gb: re.compile(r"^\s*end\s+interface"), # put at end of interface
        lambda phys, sub, gb: get_subroutine_begin_regex(sub + "_bridge_c"), # sub_c begin
        lambda phys, sub, gb: get_subroutine_end_regex(sub + "_bridge_c"),   # sub_c begin
        lambda *x           : "The f90 to c fortran subroutine declaration(<name>_bridge_c)"
    )),

    ("cxx_f2c_bind_impl"  , (
        lambda phys, sub, gb: f"tests/infra/{phys}_test_data.cpp",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_f2c_bind_impl"),
        lambda phys, sub, gb: get_plain_comment_regex(comment='extern "C" : end _c decls'), # reqs special comment
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub + "_bridge_c"),      # sub_bridge_c begin
        lambda phys, sub, gb: get_cxx_close_block_regex(at_line_start=True), # terminating }
        lambda *x           : "The f90 to cxx function declaration and implementation(<name>_bridge_c)"
    )),

    ("cxx_c2f_bind_decl"  , (
        lambda phys, sub, gb: f"tests/infra/{phys}_test_data.cpp",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_bind_decl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(comment='extern "C" : end _f decls'), # reqs special comment
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub + "_bridge_f"), # c decl
        lambda phys, sub, gb: re.compile(r".*;\s*$"),                   # ;
        lambda *x           : "The c to f90 c function declaration(<name>_bridge_f)"
    )),

    ("cxx_c2f_glue_decl"  , (
        lambda phys, sub, gb: f"tests/infra/{phys}_test_data.hpp",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_glue_decl"),
        lambda phys, sub, gb: re.compile(r"^\s*// End glue function decls"),  # reqs special comment
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub + "_f"), # cxx(data) decl
        lambda phys, sub, gb: re.compile(r".*;\s*"),             # ;
        lambda *x           : "The cxx to f90 function declaration(<name>_f(Data))"
    )),

    ("cxx_c2f_glue_impl" , (
        lambda phys, sub, gb: f"tests/infra/{phys}_test_data.cpp",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_glue_impl"),
        lambda phys, sub, gb: re.compile(r"^\s*// end glue impls"), # reqs special comment
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub + "_f"), # cxx(data)
        lambda phys, sub, gb: get_cxx_close_block_regex(at_line_start=True), # terminating }
        lambda *x           : "The cxx to f90 function implementation(<name>_f(Data))"
    )),

    ("cxx_t2cxx_glue_decl" , (
        lambda phys, sub, gb: f"tests/infra/{phys}_test_data.hpp",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_t2cxx_glue_decl"),
        lambda phys, sub, gb: re.compile(r"^\s*// End glue function decls"),  # reqs special comment
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub), # cxx(data)
        lambda phys, sub, gb: re.compile(r".*;\s*"),             # ;
        lambda *x           : "The CXX test data to CXX function declaration(<name>(Data))"
    )),

    ("cxx_t2cxx_glue_impl" , (
        lambda phys, sub, gb: f"tests/infra/{phys}_test_data.cpp",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_t2cxx_glue_impl"),
        lambda phys, sub, gb: re.compile(r"^\s*// end glue impls"), # reqs special comment
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub), # cxx(data)
        lambda phys, sub, gb: get_cxx_close_block_regex(at_line_start=True), # terminating }
        lambda *x           : "The CXX test data to CXX function implementation(<name>(Data))"
    )),

    ("cxx_c2f_data"  , (
        lambda phys, sub, gb: f"tests/infra/{phys}_test_data.hpp",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_data"),
        lambda phys, sub, gb: re.compile(r"^\s*// Glue functions for host"),  # reqs special comment
        lambda phys, sub, gb: get_cxx_struct_begin_regex(get_data_struct_name(sub)), # struct Sub
        lambda phys, sub, gb: get_cxx_close_block_regex(semicolon=True),             # terminating };
        lambda *x           : "The cxx data struct definition(struct Data)"
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
        lambda phys, sub, gb: f"tests/infra/{phys}_unit_tests_common.hpp",
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_bfb_unit_decl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(semicolon=True, comment="UnitWrap"), # Insert at end of UnitWrap struc
        lambda phys, sub, gb: get_cxx_struct_begin_regex(get_data_test_struct_name(sub)), # struct decl
        lambda phys, sub, gb: re.compile(r".*;\s*$"), # end of struct decl
        lambda *x           : "The cxx unit test struct declaration"
    )),

    ("cxx_bfb_unit_impl", (
        lambda phys, sub, gb: f"tests/{phys}_{sub}_tests.cpp",
        lambda phys, sub, gb: create_template(phys, sub, gb, "cxx_bfb_unit_impl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(semicolon=True, at_line_start=True), # insert of end of struct
        lambda phys, sub, gb: get_cxx_function_begin_regex("run_bfb"),  # run_bfb
        lambda phys, sub, gb: get_cxx_close_block_regex(comment="run_bfb"), # } // run_bfb # reqs special comment
        lambda *x           : "The cxx bfb unit test implementation"
    )),

    ("cxx_eti", (
        lambda phys, sub, gb: f"eti/{phys}_{sub}.cpp",
        lambda phys, sub, gb: create_template(phys, sub, gb, "cxx_eti", force=True),
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
#   (path-to-origin, path-to-cxx-src, init-code)
ORIGIN_FILES, CXX_ROOT, INIT_CODE, FINALIZE_CODE, COLS_DIMNAME, UNPACKED = range(6)
PHYSICS = {
    "p3"   : (
        ("components/eam/src/physics/cam/micro_p3.F90",),
        "components/eamxx/src/physics/p3",
        "p3_init();",
        "",
        "its:ite",
        False
    ),
    "shoc" : (
        ("components/eam/src/physics/cam/shoc.F90",),
        "components/eamxx/src/physics/shoc",
        "shoc_init(d.nlev, true);",
        "",
        "shcol",
        False
    ),
    "dp" : (
        (
            "components/eam/src/control/apply_iop_forcing.F90",
            "components/eam/src/dynamics/se/se_iop_intr_mod.F90",
            "components/eam/src/control/iop_data_mod.F90",
            "components/eam/src/control/history_iop.F90"
        ),
        "components/eamxx/src/physics/dp",
        "dp_init(d.plev, true);",
        ""
        "",
        False
    ),
    "gw" : (
        (
            "components/eam/src/physics/cam/gw/gw_common.F90",
            "components/eam/src/physics/cam/gw/gw_convect.F90",
            "components/eam/src/physics/cam/gw/gw_diffusion.F90",
            "components/eam/src/physics/cam/gw/gw_oro.F90",
            "components/eam/src/physics/cam/gw/gw_utils.F90",
            "components/eam/src/physics/cam/gw/gw_front.F90"
        ),
        "components/eamxx/src/physics/gw",
        "gw_common_init(); // Might need more specific init",
        "gw_finalize_cxx();",
        "ncol",
        True
    ),
}

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
    Get a regex object to match a fortran subroutine with a certain name
    """
    subroutine_begin_regex_str = fr"^\s*subroutine\s+{name}\s*[(]"
    return re.compile(subroutine_begin_regex_str)

###############################################################################
def get_function_begin_regex(name):
###############################################################################
    """
    Get a regex to match a fortran function with a certain name
    """
    function_begin_regex_str = fr"^\s*((pure\s+)?function)\s+{name}\s*[(].*result\s*[(]\s*([^) ]+)"
    return re.compile(function_begin_regex_str)

###############################################################################
def get_subroutine_end_regex(name):
###############################################################################
    """
    Get a regex to match the end of a fortran function with a certain name
    """
    subroutine_end_regex_str = fr"^\s*end\s+(subroutine|function)\s+{name}\s*$"
    return re.compile(subroutine_end_regex_str)

###############################################################################
def get_cxx_function_begin_regex(name, static=False, template=None):
###############################################################################
    """
    Get a regex to match a cxx function with a certain name and properties
    """
    static_regex_str = r"static\s+" if static else ""
    template_regex_str = fr"{template}::" if template else ""
    function_begin_regex_str = fr"^\s*{static_regex_str}void\s+{template_regex_str}{name}\s*[(]"
    return re.compile(function_begin_regex_str)

###############################################################################
def get_cxx_close_block_regex(semicolon=False, comment=None, at_line_start=False):
###############################################################################
    """
    Get a regex to match the closing of a block with certain properties
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
    Get a regex to match the closing of a namespace
    """
    return get_cxx_close_block_regex(comment=fr"namespace\s+{namespace}")

###############################################################################
def get_cxx_struct_begin_regex(struct):
###############################################################################
    """
    Get a regex to match the beginning of a struct with a certain name
    """
    struct_regex_str = fr"^\s*struct\s+{struct}([\W]|$)"
    return re.compile(struct_regex_str)

###############################################################################
def get_plain_comment_regex(comment):
###############################################################################
    """
    Get a regex to match a C comment with a certain text
    """
    comment_regex_str   = fr"^\s*//\s*{comment}"
    return re.compile(comment_regex_str)

###############################################################################
def get_data_struct_name(sub):
###############################################################################
    """
    For a function name, return the corresponding struct name
    """
    return "".join([item.capitalize() for item in sub.split("_")]) + "Data"

###############################################################################
def get_data_test_struct_name(sub):
###############################################################################
    """
    For a function name, return the corresponding test struct name
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
    """
    filepath = gb.get_path_for_piece_file(physics, sub, piece)
    if not filepath.exists() or force:
        expect(piece in FILE_TEMPLATES,
               f"{filepath} does not exist and there is no template for generating files for piece {piece}")

        gen_code = getattr(gb, f"gen_{piece}")(physics, sub, force_arg_data=force_arg_data)
        contents = FILE_TEMPLATES[piece](physics, sub, gen_code)
        if gb.dry_run():
            print(f"Would create file {filepath.relative_to(Path.cwd().parent.parent)} with contents:\n{contents}")
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
    Remove comments and whitespaces from fortran code
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
def normalize_f90(contents):
###############################################################################
    """
    Normalize fortran by removing comments, whitespaces, resolving line continuations,
    and lowercasing everything.
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

    return ("\n".join(new_lines)).lower()

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
    all_dims = []
    for name_dim in names_dims:
        if "(" in name_dim:
            name, dims_raw = name_dim.split("(")
            dims_raw = dims_raw.rstrip(")").strip()
            dims_check = tuple(item.replace(" ", "") for item in dims_raw.split(","))
            all_dims.append(dims_check)
            names.append(name.strip())
        else:
            all_dims.append(dims)
            names.append(name_dim.strip())

    return [(name, argtype, intent, dims) for name, dims in zip(names, all_dims)]

###############################################################################
def parse_origin(contents, subs):
###############################################################################
    """
    Returns a map of subname->[(argname, argtype, intent, dims)]
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
                            dim_scalars = extract_dim_scalars(arg_dim)
                            for arg_dim in dim_scalars:
                                if arg_dim not in arg_names:
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
    """
    is_ptr = arg_datum[ARG_DIMS] is not None or arg_datum[ARG_INTENT] != "in"
    arg_type = arg_datum[ARG_TYPE]
    arg_cxx_type = get_cxx_scalar_type(arg_type)
    return f"{arg_cxx_type}{'*' if is_ptr else ''}"

###############################################################################
def get_kokkos_type(arg_datum, col_dim, unpacked=False):
###############################################################################
    """
    Based on arg datum, give c++ kokkos type

    Note: We can only guess at the correct types, especially whether an argument
    should be packed data or not!
    """
    expect(col_dim is not None, "Kokkos must know col_dim name")

    is_const  = arg_datum[ARG_INTENT] == "in"
    dims      = arg_datum[ARG_DIMS]
    # Remove the cols dim, we parallelize over cols
    if dims is not None and col_dim in dims:
        dims = list(dims)
        dims.remove(col_dim)

    is_view   = bool(dims)
    arg_type  = arg_datum[ARG_TYPE]
    if is_custom_type(arg_type):
        kokkos_type = arg_type.split("::")[-1]
    else:
        kokkos_type = CXX_TYPE_MAP[arg_type]
        if kokkos_type == "Real" and not unpacked:
            kokkos_type = "Spack"

    base_type = f"{'const ' if is_const else ''}{kokkos_type}"

    return f"const uview_{len(dims)}d<{base_type}>&" if is_view else f"{base_type}&"

###############################################################################
def gen_arg_cxx_decls(arg_data, kokkos=False, unpacked=False, col_dim=None):
###############################################################################
    """
    Get all arg decls for a set of arg data. kokkos flag tells us to use C++/Kokkos
    types instead of C types.
    """
    arg_names    = [item[ARG_NAME] for item in arg_data]
    if kokkos:
        arg_types = [get_kokkos_type(item, col_dim, unpacked=unpacked) for item in arg_data]
    else:
        arg_types = [get_cxx_type(item) for item in arg_data]

    arg_sig_list = [(f"{arg_type} {arg_name}", arg_datum[ARG_INTENT])
                    for arg_name, arg_type, arg_datum in zip(arg_names, arg_types, arg_data)]

    # For kokkos functions, we will almost always want the team and we don't want
    # the col_dim
    if kokkos:
        arg_sig_list.insert(0, ("const MemberType& team", "in"))
        for arg_sig, arg_intent in arg_sig_list:
            if arg_sig.split()[-1] == col_dim:
                expect(arg_intent == "in", f"col_dim {col_dim} wasn't an input, {arg_intent}?")
                arg_sig_list.remove((arg_sig, arg_intent))
                break

    result = []

    # For permanent sigs, we want them to look nice. We may want to order these
    # by intent and scalar vs array, but for now, just mimic the fortran order.
    if kokkos:
        intent_map = {"in" : "Inputs", "inout" : "Inputs/Outputs", "out" : "Outputs"}
        curr = None
        for arg_sig, arg_intent in arg_sig_list:
            if arg_intent != curr:
                fullname = intent_map[arg_intent]
                result.append(f"// {fullname}")
                curr = arg_intent

            result.append(arg_sig)

    else:
        result = [arg_sig for arg_sig, _ in arg_sig_list]

    return result

###############################################################################
def split_by_intent(arg_data):
###############################################################################
    """
    Take arg data and split into three lists of names based on intent: [inputs], [intouts], [outputs]
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
def gen_cxx_data_args(arg_data):
###############################################################################
    """
    Based on data, generate unpacking of Data struct args
    """
    all_dims = group_data(arg_data)[0]
    args_needs_ptr = [item[ARG_DIMS] is None and item[ARG_INTENT] != "in" for item in arg_data]
    arg_names      = [item[ARG_NAME] for item in arg_data]
    arg_dim_call   = [item[ARG_NAME] in all_dims for item in arg_data]
    args = [f"{'&' if need_ptr else ''}d.{arg_name}"
            for arg_name, need_ptr, dim_call in zip(arg_names, args_needs_ptr, arg_dim_call)]
    return args

###############################################################################
def gen_arg_f90_decls(arg_data):
###############################################################################
    """
    Generate f90 argument declarations, will attempt to group these together if possible.
    """
    metadata = {}
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
    """
    Return if arg_data contains any array data
    """
    for _, _, _, dims in arg_data:
        if dims is not None:
            return True

    return False

###############################################################################
def gen_struct_members(arg_data):
###############################################################################
    """
    Gen cxx code for data struct members
    """
    metadata = {} # intent -> (type, is_ptr) -> names
    for name, argtype, intent, dims in arg_data:
        metadata.setdefault(intent, {}).setdefault((argtype, dims is not None), []).append(name)

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
def extract_dim_tokens(dim):
###############################################################################
    """
    Given a dimensions spec, extract the tokens
    """
    return [item.replace(" ", "") for item in dim.split(":") if item.strip() != ""]

###############################################################################
def extract_dim_scalars(dim):
###############################################################################
    """
    Given a dimensions spec, extract the scalar variables involved
    """
    tokens = extract_dim_tokens(dim)
    result = []
    for item in tokens:
        # Strip off the negative if it's there
        token = item.lstrip("-")
        # Skip hardcoded ints
        if not token.isdigit() and token not in result:
            result.append(token)

    return result

###############################################################################
def group_data(arg_data, filter_out_intent=None, filter_scalar_custom_types=False):
###############################################################################
    """
    Given data, return ([all-dims], [scalars], {dims->[real_data]}, {dims->[int_data]}, {dims->[bool_data]})
    """
    scalars  = []

    all_dims = []

    for name, argtype, _, dims in arg_data:
        if dims is not None:
            for dim in dims:
                dscalars = extract_dim_scalars(dim)
                for dscalar in dscalars:
                    if dscalar not in all_dims:
                        all_dims.append(dscalar)

    real_data = {}
    int_data = {}
    bool_data = {}

    for name, argtype, intent, dims in arg_data:
        if filter_out_intent is None or intent != filter_out_intent:
            if dims is None:
                if not (is_custom_type(argtype) and filter_scalar_custom_types):
                    if name not in all_dims:
                        scalars.append( (name, get_cxx_scalar_type(argtype)))
                    else:
                        expect(argtype == "integer", f"Expected dimension {name} to be of type integer")
                        expect(intent == "in", f"Expected dimension {name} to be intent in")

            elif argtype == "integer":
                int_data.setdefault(dims, []).append(name)

            elif argtype == "real":
                real_data.setdefault(dims, []).append(name)

            elif argtype == "logical":
                bool_data.setdefault(dims, []).append(name)

    return all_dims, scalars, real_data, int_data, bool_data

###############################################################################
def get_list_of_lists(items, indent):
###############################################################################
    result = "{\n"
    for item in items:
        result += f"{indent}{{{item}}},\n"
    result = result.rstrip(",\n")
    result += f"\n{indent[0:-2]}}}"

    return result

###############################################################################
def convert_to_cxx_dim(dim, add_underscore=False, from_d=False):
###############################################################################
    """
    Convert a fortran dim, potentially containing a range, to an item count
    """
    tokens = extract_dim_tokens(dim)
    uns = "_" if add_underscore else ""
    obj = "d." if from_d else ""

    # null case, could not determine anything
    if not tokens:
        return ""

    # case 1, single token
    elif len(tokens) == 1:
        expect(not tokens[0].startswith("-"), f"Received weird negative fortran dim: '{dim}'")
        return obj + tokens[0] + uns

    # case 2, multiple tokens
    elif len(tokens) == 2:
        first_token = tokens[0]
        second_token = tokens[1]
        first_int = None
        second_int = None
        try:
            first_int = int(first_token)
        except ValueError:
            pass

        try:
            second_int = int(second_token)
        except ValueError:
            pass

        # case 2.1, both tokens are ints
        if first_int is not None and second_int is not None:
            return str(second_int - (first_int - 1))

        # case 2.2, first token is an int
        elif first_int is not None:
            if first_int <= 0:
                return f"{obj}{second_token}{uns} + {1 + abs(first_int)}"
            elif first_int == 1:
                return second_token
            else:
                return f"{obj}{second_token}{uns} - {first_int - 1}"

        # case 2.3, second token is an int
        elif second_int is not None:
            expect(False, f"Received weird fortran range with the 2nd token as int: '{dim}'")

        # case 2.4, first token is negative
        elif first_token.startswith("-"):
            if first_token.strip("-") == second_token:
                return f"{obj}{second_token}{uns}*2 + 1"
            else:
                return f"{obj}{first_token.strip('-')}{uns} + {obj}{second_token}{uns} + 1"

        else:
            return f"{obj}{second_token}{uns} - {obj}{first_token}{uns}"

    else:
        expect(False, f"Received weird fortran range with more than 2 tokens: '{tokens}'")

###############################################################################
def gen_struct_api(struct_name, arg_data):
###############################################################################
    """
    Given data, generate code for data struct api
    """
    all_dims, scalars, real_data, int_data, bool_data = group_data(arg_data, filter_scalar_custom_types=True)

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
            dim_cxx_vec.append(f"{', '.join([convert_to_cxx_dim(item, True) for item in dims])}")
            data_vec.append(f"{', '.join(['&{}'.format(item) for item in items])}")

    parent_call = "  PhysicsTestData("
    parent_call += get_list_of_lists(dim_cxx_vec, "      ")
    parent_call += ",\n    "
    parent_call += get_list_of_lists(real_vec, "      ")

    if int_vec:
        parent_call += ",\n    "
        parent_call += get_list_of_lists(int_vec, "      ")
    if bool_vec:
        parent_call += ",\n    "
        parent_call += get_list_of_lists(bool_vec, "      ")

    parent_call += "),\n"

    parent_call += f"    {', '.join(['{0}({0}_)'.format(name) for name, _ in cons_args])}"

    result.append(parent_call)
    result.append("{}")
    result.append("")

    result.append("PTD_STD_DEF({}, {}, {});".\
                  format(struct_name, len(cons_args), ", ".join([name for name, _ in cons_args])))

    return result

###############################################################################
def find_insertion(lines, insert_regex):
###############################################################################
    """
    Find the index that matches insert_regex. If not found, return None
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

    rank_map_sorted = dict(sorted(rank_map.items()))

    return rank_map_sorted

PREFIX_MAP = {"Real" : "r", "Int" : "i", "bool" : "b"}
###############################################################################
def get_view_type(typename, rank):
###############################################################################
    """
    Return a device view type name based on base type and rank
    """
    expect(typename in PREFIX_MAP, f"Unknown typename {typename}")
    return f"view{rank}d{PREFIX_MAP[typename]}_d"

###############################################################################
def has_uniform_sizes(arg_data, rank, arg_list):
###############################################################################

    for rank_itr in range(rank):
        dims = [convert_to_cxx_dim(get_data_by_name(arg_data, arg_name, ARG_DIMS)[rank_itr]) for arg_name in arg_list]
        if len(set(dims)) > 1:
            return False

    return True

###############################################################################
def get_htd_dth_call(arg_data, rank, arg_list, typename, is_output=False, f2c=False):
###############################################################################
    result = ""

    view_type = get_view_type(typename, rank)
    prefix_char = PREFIX_MAP[typename]
    arg_list_d = arg_list if f2c else [f"d.{item}" for item in arg_list]
    arg_list_v = [f"{item}_d" for item in arg_list]
    vec_name = f"vec{rank}d{prefix_char}_{'out' if is_output else 'in'}"
    if is_output:
        result += f"  std::vector<{view_type}> {vec_name} = {{{', '.join(arg_list_v)}}};\n"
    else:
        result += f"  std::vector<{view_type}> {vec_name}({len(arg_list)});\n"
    funcname = "ekat::device_to_host" if is_output else "ekat::host_to_device"

    if has_uniform_sizes(arg_data, rank, arg_list):
        dims = [convert_to_cxx_dim(get_data_by_name(arg_data, arg_list[0], ARG_DIMS)[rank_itr], from_d=not f2c) for rank_itr in range(rank)]
        result += f"  {funcname}({{{', '.join(arg_list_d)}}}, {', '.join(dims)}, {vec_name});\n\n"

    else:
        for rank_itr in range(rank):
            dims = [convert_to_cxx_dim(get_data_by_name(arg_data, arg_name, ARG_DIMS)[rank_itr], from_d=not f2c) for arg_name in arg_list]
            result += f"  std::vector<int> {vec_name}_{rank_itr}_sizes = {{{', '.join(dims)}}};\n"

        dim_vectors = [f"{vec_name}_{rank_itr}_sizes" for rank_itr in range(rank)]
        result += f"  {funcname}({{{', '.join(arg_list_d)}}}, {', '.join(dim_vectors)}, {vec_name});\n\n"

    return result

###############################################################################
def gen_glue_impl(phys, sub, arg_data, arg_names, col_dim, f2c=False, unpacked=False):
###############################################################################
    """
    Generate code that takes a TestData struct and unpacks it to call the CXX
    version of the physics function.
    """
    init_code = get_physics_data(phys, INIT_CODE)
    init_code = f"  {init_code}\n\n" if init_code else ""
    final_code = get_physics_data(phys, FINALIZE_CODE)
    final_code = f"  {final_code}\n" if final_code else ""
    obj = "" if f2c else "d."

    impl = "#if 0\n" # There's no way to guarantee this code compiles

    if not f2c:
        impl += init_code

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
        all_inputs  = inputs + inouts
        all_outputs = inouts + outputs

        iscalars = list(sorted(set(all_inputs) & set(scalars)))
        oscalars = list(sorted(set(all_outputs) & set(scalars)))

        oviews = list(sorted(set(all_outputs) & set(views)))

        vreals = list(sorted(set(reals) & set(views)))
        vints  = list(sorted(set(ints)  & set(views)))
        vbools = list(sorted(set(bools) & set(views)))

        sreals = list(sorted(set(reals) & set(scalars)))
        sints  = list(sorted(set(ints)  & set(scalars)))
        sbools = list(sorted(set(bools) & set(scalars)))

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
        # 1) Set up typedefs (or just have these at the top of file so they can be shared?)
        #

        # set up basics

        type_list    = ["Real", "Int", "bool"]
        impl += "  // create device views and copy\n"

        #
        # 2) Sync to device. Do ALL views, not just inputs
        #

        for input_group, typename in zip([vreals, vints, vbools], type_list):
            if input_group:
                rank_map = get_rank_map(arg_data, input_group)

                for rank, arg_list in rank_map.items():
                    impl += get_htd_dth_call(arg_data, rank, arg_list, typename, f2c=f2c)

        #
        # 3) Unpack view array
        #

        for input_group, typename in zip([vreals, vints, vbools], type_list):
            prefix_char = PREFIX_MAP[typename]
            if input_group:
                rank_map = get_rank_map(arg_data, input_group)

                for rank, arg_list in rank_map.items():
                    view_type = get_view_type(typename, rank)
                    impl += f"  {view_type}\n"
                    for idx, input_item in enumerate(arg_list):
                        impl += f"    {input_item}_d(vec{rank}d{prefix_char}_in[{idx}]){';' if idx == len(arg_list) - 1 else ','}\n"
                    impl += "\n"


        #
        # 4) Get nk_pack and policy, unpack scalars, and launch kernel
        #
        if unpacked:
            impl += f"  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy({obj}{col_dim}, {obj}nlev);\n\n"
        else:
            impl += f"  const Int nk_pack = ekat::npack<Spack>({obj}nlev);\n"
            impl += f"  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy({obj}{col_dim}, nk_pack);\n\n"

        if scalars:
            if not f2c:
                impl += "  // unpack data scalars because we do not want the lambda to capture d\n"
                for input_group, typename in zip([isreals, isints, isbools], type_list):
                    if input_group:
                        for arg in input_group:
                            if arg not in oscalars and arg != col_dim:
                                impl += f"  const {typename} {arg} = {obj}{arg};\n"

            # We use 0-rank views to handle output scalars
            for output_group, typename in zip([osreals, osints, osbools], type_list):
                if output_group:
                    view_type = get_view_type(typename, 0)
                    hview_type = view_type[0:-2] + "_h"
                    for arg in output_group:
                        impl += f'  {hview_type} {arg}_h("{arg}_h");\n'
                        if arg in iscalars:
                            impl += f'  {arg}_h() = {obj}{arg};\n'
                        impl += f'  {view_type} {arg}_d = Kokkos::create_mirror_view_and_copy(DefaultDevice(), {arg}_h);\n'

            impl += "\n"

        impl += "  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {\n"
        impl += "    const Int i = team.league_rank();\n\n"

        #
        # 5) Get subviews
        #
        impl += "    // Get single-column subviews of all inputs, shouldn't need any i-indexing\n"
        impl += "    // after this.\n"

        for view_arg in views:
            dims = get_data_by_name(arg_data, view_arg, ARG_DIMS)
            if col_dim in dims:
                if len(dims) == 1:
                    pass
                else:
                    impl += f"    const auto {view_arg}_c = ekat::subview({view_arg}_d, i);\n"

        impl += "\n"

        #
        # 6) Call fn
        #
        kernel_arg_names = ["team"]
        for arg_name in arg_names:
            if arg_name in views:
                dims = get_data_by_name(arg_data, arg_name, ARG_DIMS)
                if len(dims) == 1 and col_dim in dims:
                    kernel_arg_names.append(f"{arg_name}_d(i)")
                else:
                    kernel_arg_names.append(f"{arg_name}_c")
            elif arg_name in oscalars:
                kernel_arg_names.append(f"{arg_name}_d()")
            else:
                if arg_name != col_dim:
                    kernel_arg_names.append(arg_name)

        joinstr = ',\n      '
        impl += f"    SHF::{sub}(\n      {joinstr.join(kernel_arg_names)});\n"
        impl +=  "  });\n\n"

        #
        # 7) Sync back to host
        #
        if oscalars:
            impl += "  // Get outputs back, start with scalars\n"
            for arg in oscalars:
                impl += f"  Kokkos::deep_copy({arg}_h, {arg}_d);\n"
                impl += f"  {obj}{arg} = {arg}_h();\n"

            impl += "\n"

        if oviews:
            impl += "  // Now get arrays\n"
            for output_group, typename in zip([ovreals, ovints, ovbools], type_list):
                if output_group:
                    rank_map = get_rank_map(arg_data, output_group)

                    for rank, arg_list in rank_map.items():
                        impl += get_htd_dth_call(arg_data, rank, arg_list, typename, is_output=True, f2c=f2c)

    else:
        expect(False, "Not yet supported")

    if not f2c:
        impl += final_code

    impl += "#endif"

    return impl

#
# Main classes
#

###############################################################################
class GenBoiler(object):
###############################################################################

    ###########################################################################
    def __init__(self,
                 physics,
                 subs        = None,
                 pieces      = get_supported_pieces(),
                 overwrite   = False,
                 kernel      = False,
                 source_repo = get_git_toplevel_dir(),
                 target_repo = get_git_toplevel_dir(),
                 f2c         = False,
                 unpacked    = False,
                 col_dim     = None,
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
        self._subs        = [] if subs is None else subs
        self._pieces      = pieces
        self._physics     = physics
        self._overwrite   = overwrite
        self._kernel      = kernel
        self._source_repo = Path(source_repo).resolve()
        self._target_repo = Path(target_repo).resolve()
        self._unpacked    = unpacked if unpacked else get_physics_data(physics, UNPACKED)
        self._col_dim     = get_physics_data(physics, COLS_DIMNAME) if col_dim is None else col_dim
        self._dry_run     = dry_run
        self._verbose     = verbose

        # internals
        self._db          = {}

        if not self._pieces:
            self._pieces = get_supported_pieces()

        # If user did not request f2c bridging, don't generate them
        if self._pieces == get_supported_pieces() and not f2c:
            self._pieces = [piece for piece in self._pieces if "f2c" not in piece]

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
    def gen_cxx_c2f_bind_decl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        In C, generate the C to F90 brige declaration. The definition will be in fortran
        """
        arg_data = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        arg_decls = gen_arg_cxx_decls(arg_data, unpacked=self._unpacked)
        result = f"void {sub}_bridge_f({', '.join(arg_decls)});\n"
        return result

    ###########################################################################
    def gen_f90_c2f_bind(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        In F90, generate the C to F90 bridge function.
        """
        arg_data = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        arg_names = ", ".join([item[ARG_NAME] for item in arg_data])
        arg_decls = gen_arg_f90_decls(arg_data)
        phys_mod = "micro_p3" if phys == "p3" else phys
        result = \
"""  subroutine {sub}_bridge_f({arg_names}) bind(C)
    use {phys_mod}, only : {sub}

    {arg_decls}

    call {sub}({arg_names})
  end subroutine {sub}_bridge_f""".format(sub=sub, arg_names=arg_names, phys_mod=phys_mod, arg_decls="\n    ".join(arg_decls))

        return result

    ###########################################################################
    def gen_f90_f2c_bind(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        In F90, generate the F90 to C bridge declaration. The definition will be in C
        """
        arg_data = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        arg_names = ", ".join([item[ARG_NAME] for item in arg_data])
        arg_decls = gen_arg_f90_decls(arg_data)
        result = \
"""  subroutine {sub}_bridge_c({arg_names}) bind(C)
    use iso_c_binding

    {arg_decls}
  end subroutine {sub}_bridge_c""".format(sub=sub, arg_names=arg_names, arg_decls="\n    ".join(arg_decls))

        return result

    ###########################################################################
    def gen_cxx_c2f_glue_decl(self, phys, sub, force_arg_data=None): # pylint: disable=W0613
    ###########################################################################
        """
        In C, generate the C to F90 bridge declaration. This version takes the test
        data struct.
        """
        struct_name = get_data_struct_name(sub)
        result = f"void {sub}_f({struct_name}& d);"
        return result

    ###########################################################################
    def gen_cxx_c2f_glue_impl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        In C, generate the C to F90 bridge function. This version takes the test
        data struct.
        """
        arg_data         = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        arg_data_args    = ", ".join(gen_cxx_data_args(arg_data))
        transition_code_1 = "d.transition<ekat::TransposeDirection::c2f>();"
        transition_code_2 = "d.transition<ekat::TransposeDirection::f2c>();"
        data_struct      = get_data_struct_name(sub)
        init_code        = get_physics_data(phys, INIT_CODE)

        result = \
f"""void {sub}_f({data_struct}& d)
{{
  {transition_code_1}
  {init_code}
  {sub}_c({arg_data_args});
  {transition_code_2}
}}

"""
        return result

    ###########################################################################
    def gen_cxx_t2cxx_glue_decl(self, phys, sub, force_arg_data=None): # pylint: disable=W0613
    ###########################################################################
        """
        In C, generate the C to CXX bridge declaration. This version takes the test
        data struct.
        """
        struct_name = get_data_struct_name(sub)
        result = f"void {sub}({struct_name}& d);"
        return result

    ###########################################################################
    def gen_cxx_t2cxx_glue_impl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        In C, generate the C to CXX bridge implementation. This version takes the test
        data struct.
        """
        arg_data  = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        arg_names = [item[ARG_NAME] for item in arg_data]
        decl      = self.gen_cxx_t2cxx_glue_decl(phys, sub, force_arg_data=force_arg_data).rstrip(";")

        impl = gen_glue_impl(phys, sub, arg_data, arg_names, self._col_dim, unpacked=self._unpacked)

        result = \
f"""{decl}
{{
{impl}
}}
"""
        return result


    ###########################################################################
    def gen_cxx_c2f_data(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        In CXX, generate the the test data struct.
        """
        arg_data         = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        struct_members   = "\n  ".join(gen_struct_members(arg_data))
        any_arrays       = has_arrays(arg_data)
        struct_name      = get_data_struct_name(sub)
        inheritance      = " : public PhysicsTestData" if any_arrays else ""
        api              = "\n  " + "\n  ".join(gen_struct_api(struct_name, arg_data) if any_arrays else "")

        result = \
f"""struct {struct_name}{inheritance} {{
  {struct_members}{api}
}};

"""

        result = "\n".join([item.rstrip() for item in result.splitlines()])

        return result

    ###########################################################################
    def gen_cxx_f2c_bind_impl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        In C, generate the F90 to C bridge declatation and implementation.
        """
        arg_data  = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        arg_names = [item[ARG_NAME] for item in arg_data]
        arg_decls = gen_arg_cxx_decls(arg_data, unpacked=self._unpacked)
        decl      = f"void {sub}_bridge_c({', '.join(arg_decls)})"

        impl = gen_glue_impl(phys, sub, arg_data, arg_names, self._col_dim, f2c=True)

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
        In CXX, generate the function declaration in the main physics header
        """
        arg_data = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        arg_decls = gen_arg_cxx_decls(arg_data, kokkos=True, unpacked=self._unpacked, col_dim=self._col_dim)

        arg_decls_str = ("\n    ".join([item if item.startswith("//") else f"{item}," for item in arg_decls])).rstrip(",")

        return f"  KOKKOS_FUNCTION\n  static void {sub}(\n    {arg_decls_str});"

    ###########################################################################
    def gen_cxx_incl_impl(self, phys, sub, force_arg_data=None): # pylint: disable=W0613
    ###########################################################################
        """
        In CXX, generate the include of the implementation header
        """
        impl_path = get_piece_data(phys, sub, "cxx_func_impl", FILEPATH, self)
        return f'# include "{impl_path}"'

    ###########################################################################
    def gen_cxx_func_impl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        In CXX, generate a stub implemenation of the function
        """
        decl = self.gen_cxx_func_decl(phys, sub, force_arg_data=force_arg_data).rstrip(";").replace("void ", "void Functions<S,D>::").replace("static ", "")
        decls = [line.strip() for line in decl.splitlines()]
        for idx, item in enumerate(decls):
            if idx > 1:
                decls[idx] = f"  {item.strip()}"
            else:
                decls[idx] = item.strip()

        decl = "\n".join(decls)
        # I don't think any intelligent guess at an impl is possible here
        result = \
f"""template<typename S, typename D>
{decl}
{{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}}"""
        return result

    ###########################################################################
    def gen_cxx_bfb_unit_decl(self, phys, sub, force_arg_data=None): # pylint: disable=W0613
    ###########################################################################
        """
        In CXX, within the common test struct, declare the test struct for this function.
        """
        test_struct = get_data_test_struct_name(sub)
        return f"    struct {test_struct};"

    ###########################################################################
    def gen_cxx_bfb_unit_impl(self, phys, sub, force_arg_data=None):
    ###########################################################################
        """
        In CXX, generate the BFB unit test implementation
        """
        arg_data = force_arg_data if force_arg_data else self._get_arg_data(phys, sub)
        data_struct = get_data_struct_name(sub)
        has_array = has_arrays(arg_data)

        gen_random = \
"""

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }"""

        _, scalars, real_data, int_data, bool_data = group_data(arg_data, filter_out_intent="in")
        check_scalars, check_arrays, scalar_comments = "", "", ""
        for scalar in scalars:
            check_scalars += f"        REQUIRE(d_baseline.{scalar[0]} == d_test.{scalar[0]});\n"

        all_dims, input_scalars, _, _, _ = group_data(arg_data, filter_out_intent="out")
        all_scalar_inputs = all_dims + [scalar_name for scalar_name, _ in input_scalars]
        scalar_comments = "// " + ", ".join(all_scalar_inputs)

        if has_array:
            all_data = dict(real_data)
            for type_data in [int_data, bool_data]:
                for k, v in type_data.items():
                    if k in all_data:
                        all_data[k].extend(v)
                    else:
                        all_data[k] = v

            for _, data in all_data.items():
                for datum in data:
                    check_arrays += f"        REQUIRE(d_baseline.total(d_baseline.{data[0]}) == d_test.total(d_test.{datum}));\n"

                check_arrays += f"        for (Int k = 0; k < d_baseline.total(d_baseline.{data[0]}); ++k) {{\n"
                for datum in data:
                    check_arrays += f"          REQUIRE(d_baseline.{datum}[k] == d_test.{datum}[k]);\n"

                check_arrays += "        }\n"

        if has_array:
            result = \
"""  void run_bfb()
  {{
    auto engine = Base::get_engine();

    // Set up inputs
    {data_struct} baseline_data[] = {{
      // TODO
      {scalar_comments}
      {data_struct}(),
    }};

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof({data_struct});{gen_random}

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    {data_struct} test_data[] = {{
      // TODO
      {data_struct}(baseline_data[0]),
    }};

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {{
      for (auto& d : baseline_data) {{
        d.read(Base::m_ifile);
      }}
    }}

    // Get data from test
    for (auto& d : test_data) {{
      if (this->m_baseline_action == GENERATE) {{
        {sub}_f(d);
      }}
      else {{
        {sub}(d);
      }}
    }}

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {{
      for (Int i = 0; i < num_runs; ++i) {{
        {data_struct}& d_baseline = baseline_data[i];
        {data_struct}& d_test = test_data[i];
{check_scalars}{check_arrays}      }}
    }}
    else if (this->m_baseline_action == GENERATE) {{
      for (Int i = 0; i < num_runs; ++i) {{
        test_data[i].write(Base::m_ofile);
      }}
    }}
  }} // run_bfb""".format(data_struct=data_struct,
                          scalar_comments=scalar_comments,
                          sub=sub,
                          gen_random=gen_random,
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
""".format(ireals=", ".join(ireals), ireal_assigns="\n        ".join(["{0}[s] = test_device(vs).{0};".format(ireal) for ireal in ireals]))

            spack_output_init = ""
            if ooreals:
                spack_output_init = \
f"""// Init outputs
      Spack {', '.join(['{}(0)'.format(ooreal) for ooreal in ooreals])};
"""

            scalars = group_data(arg_data)[1]
            func_call = f"Functions::{sub}({', '.join([(scalar if scalar in reals else 'test_device(0).{}'.format(scalar)) for scalar, _ in scalars])});"

            spack_output_to_dview = ""
            if oreals:
                spack_output_to_dview = \
"""// Copy spacks back into test_device view
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {{
        {}
      }}
""".format("\n        ".join(["test_device(vs).{0} = {0}[s];".format(oreal) for oreal in oreals]))

            result = \
"""  void run_bfb()
  {{
    auto engine = Base::get_engine();

    {data_struct} baseline_data[max_pack_size] = {{
      // TODO
    }};

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof({data_struct});{gen_random}

    // Create copies of data for use by test and sync it to device. Needs to happen before read calls so that
    // inout data is in original state
    view_1d<{data_struct}> test_device("test_device", max_pack_size);
    const auto test_host = Kokkos::create_mirror_view(test_device);
    std::copy(&baseline_data[0], &baseline_data[0] + max_pack_size, test_host.data());
    Kokkos::deep_copy(test_device, test_host);

    // Get data from fortran
    for (auto& d : baseline_data) {{
      {sub}(d);
    }}

    // Get data from test. Run {sub} from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {{
      const Int offset = i * Spack::n;

      {spack_init}
      {spack_output_init}

      {func_call}

      {spack_output_to_dview}
    }});

    Kokkos::deep_copy(test_host, test_device);

    // Verify BFB results
    if (SCREAM_BFB_TESTING) {{
      for (Int i = 0; i < num_runs; ++i) {{
        {data_struct}& d_baseline = baseline_data[i];
        {data_struct}& d_test = test_host[i];
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
    def gen_cxx_eti(self, phys, sub, force_arg_data=None): # pylint: disable=W0613
    ###########################################################################
        """
        In CXX, generate the ETI file
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
    def gen_cmake_impl_eti(self, phys, sub, force_arg_data=None): # pylint: disable=W0613
    ###########################################################################
        """
        In CMake, add the ETI file to the list of sources for the phys library
        """
        eti_src = get_piece_data(phys, sub, "cxx_eti", FILEPATH, self)
        return f"    {eti_src}"

    ###########################################################################
    def gen_cmake_unit_test(self, phys, sub, force_arg_data=None): # pylint: disable=W0613
    ###########################################################################
        """
        In CMake, add the bfb unit test to the list of sources for the test
        """
        test_src = Path(get_piece_data(phys, sub, "cxx_bfb_unit_impl", FILEPATH, self)).name
        return f"    {test_src}"

    #
    # Main methods
    #

    ###########################################################################
    def gen_piece(self, phys, sub, piece, force_arg_data=None, force_file_lines=None):
    ###########################################################################
        """
        Generate code for a specific piece for a physics+subroutine
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
            for piece in self._pieces:
                try:
                    self.gen_piece(self._physics, sub, piece)
                except SystemExit as e:
                    print(f"Warning: failed to generate subroutine {sub} piece {piece} for physics {self._physics}, error:\n{e}")
                    all_success = False

        print("ALL_SUCCESS" if all_success else "THERE WERE FAILURES")

        return all_success

if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(parse_f90_args, globals())
