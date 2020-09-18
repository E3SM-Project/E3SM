from utils import expect

from collections import OrderedDict
import pathlib, re

#
# Global hardcoded data
#

# Templates: maps piece name to generic file text
FILE_TEMPLATES = {
    "cxx_func_impl": lambda phys, sub, gen_code:
"""
#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/{physics}/{physics}_functions.hpp"
#include "physics/{physics}/{phyiscs}_functions_f90.hpp"

#include "{physics}_unit_tests_common.hpp"

namespace scream {
namespace {physics} {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::{test_data_struct} {

{gen_code}

}; // struct {test_data_struct}

} // namespace unit_test
} // namespace {physics}
} // namespace scream

namespace {

TEST_CASE({sub}_bfb, "[{physics}_functions]")
{
  using TestStruct = scream::{physics}::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::{test_data_struct};

  TestStruct::run_bfb();
}

} // empty namespace
""".format(physics=phys, test_data_struct=get_data_test_struct_name(sub), gen_code=gen_code),

    "cxx_bfb_unit_impl": lambda phys, sub, gen_code:
"""
#ifndef {phys_upper}_{sub_upper}_IMPL_HPP
#define {phys_upper}_{sub_upper}_IMPL_HPP

#include "{physics}_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace {physics} {

/*
 * Implementation of {physics} {sub}. Clients should NOT
 * #include this file, but include {physics}_functions.hpp instead.
 */

{gen_code}

} // namespace p3
} // namespace scream

#endif
""".format(physics=phys, test_data_struct=get_data_test_struct_name(sub), gen_code=gen_code, phys_upper=phys.upper(), sub_upper=sub.upper()),
}

# piece map. maps the name of a piece of boilerplate that needs to be genereate to:
#   (filepath (relative to cxx_root))
FILEPATH, FILECREATE, INSERT_REGEX, ID_SELF_BEGIN_REGEX, ID_SELF_END_REGEX = range(5)
PIECES = {
    "f90_c2f_bind"  : (
        lambda phys, sub, gb: "{}_iso_c.f90".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "f90_c2f_bind"),
        lambda phys, sub, gb: re.compile(r"^\s*end\s+module\s{}_iso_c".format(phys)),
        lambda phys, sub, gb: get_subroutine_begin_regex(sub + "_c"),
        lambda phys, sub, gb: get_subroutine_end_regex(sub + "_c"),
    ),
    "f90_f2c_bind"  : (
        lambda phys, sub, gb: "{}_iso_f.f90".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "f90_f2c_bind"),
        lambda phys, sub, gb: re.compile(r"^\s*end\s+module\s{}_iso_f".format(phys)),
        lambda phys, sub, gb: get_subroutine_begin_regex(sub + "_f"),
        lambda phys, sub, gb: get_subroutine_end_regex(sub + "_f"),
    ),
    "cxx_c2f_bind_decl"  : (
        lambda phys, sub, gb: "{}_functions_f90.cpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_bind_decl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(comment='extern "C" : end _c decls'),
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub + "_c"),
        lambda phys, sub, gb: re.compile(r".*;\s*$"),
    ),
    "cxx_c2f_bind_impl"  : (
        lambda phys, sub, gb: "{}_functions_f90.cpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_bind_impl"),
        lambda phys, sub, gb: re.compile(r"^\s*// end _c impls"),
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub),
        lambda phys, sub, gb: re.compile(r".*[}]\s*$"),
    ),
    "cxx_c2f_data"  : (
        lambda phys, sub, gb: "{}_functions_f90.hpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_data"),
        lambda phys, sub, gb: get_namespace_close_regex(phys),
        lambda phys, sub, gb: get_cxx_struct_begin_regex(get_data_struct_name(sub)),
        lambda phys, sub, gb: re.compile(r".*(namespace|struct)"),
    ),
    "cxx_f2c_bind_decl"  : (
        lambda phys, sub, gb: "{}_functions_f90.hpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_f2c_bind_decl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(comment="end _f function decls"),
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub + "_f"),
        lambda phys, sub, gb: re.compile(r".*;\s*$"),
    ),
    "cxx_f2c_bind_impl"  : (
        lambda phys, sub, gb: "{}_functions_f90.cpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_f2c_bind_impl"),
        lambda phys, sub, gb: get_namespace_close_regex(phys),
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub + "_f"),
        lambda phys, sub, gb: get_cxx_close_block_regex(comment="end {}".format(sub + "_f")),
    ),
    "cxx_func_hpp"  : (
        lambda phys, sub, gb: "{}_functions.hpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_func_hpp"),
        lambda phys, sub, gb: get_cxx_close_block_regex(semicolon=True, comment="struct Functions"),
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub),
        lambda phys, sub, gb: re.compile(r".*;\s*$"),
    ),
    "cxx_func_impl" : (
        lambda phys, sub, gb: "{}_{}_impl.hpp".format(phys, sub),
        lambda phys, sub, gb: create_template(phys, sub, gb, "cxx_func_impl"),
        lambda phys, sub, gb: get_namespace_close_regex(phys),
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub),
        lambda phys, sub, gb: re.compile(r".*{\s*$"),
    ),
    "cxx_bfb_unit_impl"  : (
        lambda phys, sub, gb: "tests/{}_{}_tests.cpp".format(phys, sub),
        lambda phys, sub, gb: create_template(phys, sub, gb, "cxx_bfb_unit_impl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(semicolon=True, comment="struct {}".format(get_data_test_struct_name(sub))),
        lambda phys, sub, gb: get_cxx_function_begin_regex("run_bfb"),
        lambda phys, sub, gb: get_cxx_close_block_regex(comment="run_bfb"),
    ),
    "cxx_bfb_unit_decl" : (
        lambda phys, sub, gb: "tests/{}_unit_tests_common.hpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_bfb_unit_decl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(semicolon=True),
        lambda phys, sub, gb: get_cxx_struct_begin_regex(get_data_test_struct_name(sub)),
        lambda phys, sub, gb: re.compile(r".*;\s*$"),
    ),
}

# physics map. maps the name of a physics packages containing the original fortran subroutines to:
#   (path-to-origin, path-to-cxx-src)
ORIGIN_FILE, CXX_ROOT = range(2)
PHYSICS = {
    "p3"   : (
        "components/cam/src/physics/cam/micro_p3.F90",
        "components/scream/src/physics/p3",
    ),
    "shoc" : (
        "components/cam/src/physics/cam/shoc.F90",
        "components/scream/src/physics/shoc",
    ),
}

#
# Free functions
#

###############################################################################
def get_subroutine_begin_regex(name):
###############################################################################
    subroutine_begin_regex_str = r"^\s*subroutine\s+{}\s*[(]".format(name)
    return re.compile(subroutine_begin_regex_str)

###############################################################################
def get_subroutine_end_regex(name):
###############################################################################
    subroutine_end_regex_str = r"^\s*end\s+subroutine\s+{}\s*$".format(name)
    return re.compile(subroutine_end_regex_str)

###############################################################################
def get_cxx_function_begin_regex(name):
###############################################################################
    function_begin_regex_str = r"^\s*void\s+{}\s*[(]".format(name)
    return re.compile(function_begin_regex_str)

###############################################################################
def get_cxx_close_block_regex(semicolon=False, comment=None):
###############################################################################
    semicolon_regex_str = r"\s*;" if semicolon else ""
    comment_regex_str   = r"\s*//\s*{}".format(comment) if comment else ""
    close_block_regex_str = re.compile(r"^\s*}{}{}\s*$".format(semicolon_regex_str, comment_regex_str))
    return re.compile(close_block_regex_str)

###############################################################################
def get_namespace_close_regex(namespace):
###############################################################################
    return get_cxx_close_block_regex(comment=r"namespace\s+{}".format(namespace))

###############################################################################
def get_cxx_struct_begin_regex(struct):
###############################################################################
    struct_regex_str = r"^\s*struct\s+{}([\W]|$)".format(struct)
    return re.compile(struct_regex_str)

###############################################################################
def get_data_struct_name(sub):
###############################################################################
    return "".join([item.capitalize() for item in sub.split("_")])

###############################################################################
def get_data_test_struct_name(sub):
###############################################################################
    return "Test{}".format(get_data_struct_name(sub))

###############################################################################
def get_supported_pieces():
###############################################################################
    return PIECES.keys()

###############################################################################
def get_supported_physics():
###############################################################################
    return PHYSICS.keys()

###############################################################################
def get_piece_data(physics, sub, piece_name, piece_data):
###############################################################################
    return PIECES[piece_name][piece_data](physics, sub)

###############################################################################
def get_physics_data(physics_name, physics_data):
###############################################################################
    return PHYSICS[physics_name][physics_data]

###############################################################################
def expect_exists(physics, sub, piece, gb):
###############################################################################
    filepath = gb.get_path_for_piece_file(physics, sub, piece, gb)
    expect(filepath.exists(), "For generating {}'s {} for phyiscs {}, expected file {} to already exist".\
           format(sub, piece, physics, filepath))

###############################################################################
def create_template(physics, sub, piece, gb):
###############################################################################
    filepath = gb.get_path_for_piece_file(physics, sub, piece)
    if not filepath.exists():
        expect(piece in FILE_TEMPLATES,
               "{} does not exist and there is no template for generating files for piece {}".format(filepath, piece))

        gen_code = getattr(gb, "gen_{}".format(piece))(physics, sub)
        contents = FILE_TEMPLATES[piece](physics, sub, gen_code)
        with filepath.open("w") as fd:
            fd.write(contents)

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

ARG_NAME, ARG_TYPE, ARG_INTENT, ARG_DIMS = range(4)
###############################################################################
def parse_f90_args(line):
###############################################################################
    """
    Given a line of fortran code declaring an argument[s], return [(argname, argtype, intent, dims)]

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
    """
    expect(line.count("::") == 1, "Expected line format 'type-info :: names' for: {}".format(line))
    metadata, names_str = line.split("::")
    names = [name.strip() for name in names_str.split(",")]
    metadata_raw = [item.strip() for item in metadata.strip().split(",")]
    metadata = []
    balanced = True
    for metadatum_raw in metadata_raw:
        if balanced:
            metadata.append(metadatum_raw)
        else:
            metadata[-1] += ",{}".format(metadatum_raw)

        balanced = metadata[-1].count("(") == metadata[-1].count(")")

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

    return [(name, argtype, intent, dims) for name in names]

###############################################################################
def parse_origin(contents, subs):
###############################################################################
    """
    Returns a map of subname->[(argname, argtype, intent, dims)]

    >>> teststr = '''
    ...
    ...   SUBROUTINE p3_get_tables(mu_r_user, revap_user, vn_user, vm_user)
    ...     ! This can be called after p3_init_b.
    ...     implicit none
    ...     real(rtype), dimension(150), intent(out) :: mu_r_user
    ...     real(rtype), dimension(300,10), intent(out) :: vn_user, vm_user, revap_user
    ...     mu_r_user(:) = mu_r_table(:)
    ...     revap_user(:,:) = revap_table(:,:)
    ...     vn_user(:,:) = vn_table(:,:)
    ...     vm_user(:,:) = vm_table(:,:)
    ...
    ...    return
    ...
    ...   end SUBROUTINE p3_get_tables
    ...
    ...   subroutine p3_set_tables(mu_r_user, revap_user, vn_user, vm_user)
    ...     ! This can be called instead of p3_init_b.
    ...     implicit none
    ...     real(rtype), dimension(150), intent(in) :: mu_r_user
    ...     real(rtype), dimension(300,10), intent(in) :: vn_user, vm_user, revap_user
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
    >>> sorted(parse_origin(teststr, ["p3_get_tables", "p3_init_b"]).items())
    [('p3_get_tables', [('mu_r_user', 'real', 'out', ('150',)), ('vn_user', 'real', 'out', ('300', '10')), ('vm_user', 'real', 'out', ('300', '10')), ('revap_user', 'real', 'out', ('300', '10'))]), ('p3_init_b', [])]
    """
    begin_regexes = [get_subroutine_begin_regex(sub) for sub in subs]
    arg_decl_regex = re.compile(r"^.+intent\s*[(]\s*(in|out|inout)\s*[)]")

    contents = normalize_f90(contents)

    db = {}
    active_sub = None
    arg_decls = []
    for line in contents.splitlines():
        begin_match = None
        for sub, begin_regex in zip(subs, begin_regexes):
            begin_match = begin_regex.match(line)
            if begin_match is not None:
                expect(active_sub is None, "subroutine {} was still active when {} began".format(active_sub, sub))
                active_sub = sub

        if active_sub:
            decl_match = arg_decl_regex.match(line)
            if decl_match is not None:
                arg_decls.extend(parse_f90_args(line))

            end_regex = get_subroutine_end_regex(active_sub)
            end_match = end_regex.match(line)
            if end_match is not None:
                expect(active_sub not in db, "Found multiple matches for {}".format(active_sub))
                db[active_sub] = arg_decls
                active_sub = None
                arg_decls = []

    return db

C_TYPE_MAP = {"real" : "c_real", "integer" : "c_int", "logical" : "c_bool"}
###############################################################################
def gen_arg_f90_decl(argtype, intent, dims, names):
###############################################################################
    expect(argtype in C_TYPE_MAP, "Unrecognized argtype: {}".format(argtype))
    c_type = C_TYPE_MAP[argtype]
    value  = ", value" if dims is None and intent == "in" else ""
    intent_s = ", intent({})".format(intent)
    dimension_s = ", dimension({})".format(", ".join(dims)) if dims is not None else ""
    names_s = ", ".join(names)
    return "{argtype}(kind={c_type}) {value}{intent}{dimension} :: {names}".\
        format(argtype=argtype, c_type=c_type, value=value, intent=intent_s, dimension=dimension_s, names=names_s)

###############################################################################
def gen_arg_f90_decls(arg_data):
###############################################################################
    """
    Generate f90 argument declarations, will attempt to group these together if possible.
    """
    metadata = OrderedDict()
    for name, argtype, intent, dims in arg_data:
        metatuple = (argtype, intent, dims)
        if metatuple in metadata:
            metadata[metatuple].append(name)
        else:
            metadata[metatuple] = [name]

    result = []
    for metatuple, names in metadata.items():
        result.append(gen_arg_f90_decl(*metatuple, names))

    return result

#
# Main classes
#

###############################################################################
class GenBoiler(object):
###############################################################################

    ###########################################################################
    def __init__(self, subs, pieces, physics, overwrite, kernel, source_repo, target_repo, dry_run):
    ###########################################################################
        self._subs        = subs
        self._pieces      = pieces
        self._physics     = physics
        self._overwrite   = overwrite
        self._kernel      = kernel
        self._source_repo = pathlib.Path(source_repo).resolve()
        self._target_repo = pathlib.Path(target_repo).resolve()
        self._dry_run     = dry_run
        self._db          = {}

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
            return db

    ###########################################################################
    def _get_arg_data(self, phys, sub):
    ###########################################################################
        phys_db = self._get_db(phys)
        expect(sub in phys_db, "No data for subroutine {} in physics {}".format(sub, phys))
        return phys_db[sub]

    ###############################################################################
    def get_path_for_piece_file(self, physics, sub, piece):
    ###############################################################################
        root_dir = pathlib.Path(get_physics_data(physics, CXX_ROOT))
        filepath = self._target_repo / root_dir / get_piece_data(physics, sub, piece, FILEPATH)
        return filepath

    ###########################################################################
    def gen_piece(self, phys, sub, piece):
    ###########################################################################
        filepath, filegen, insert_regex, self_begin_regex, self_end_regex \
            = [item(phys, sub, self) for item in PIECES[piece]]

        # TODO:
        # Gen file if not there
        # If freshly generated file, we're done
        # If file already exists...
        #   Check to see if piece already exists
        #   If yes:
        #     find the lines it occupies and overwrite them with fresh generation
        #   If no:
        #     fine insertion point and insert fresh generation
        #   write the new lines

    ###########################################################################
    def gen_f90_c2f_bind(self, phys, sub):
    ###########################################################################
        arg_data = self._get_arg_data(phys, sub)
        arg_names = ", ".join([item[ARG_NAME] for item in arg_data])
        arg_decls = gen_arg_f90_decls(arg_data)
        phys_mod = "micro_p3" if phys == "p3" else phys
        result = \
"""
  subroutine {sub}_c({arg_names}) bind(C)
    use {phys_mod}, only : {sub}

    {arg_decls}

    call {sub}({arg_names})
  end subroutine {sub}_c
""".format(sub=sub, arg_names=arg_names, phys_mod=phys_mod, arg_decls="\n    ".join(arg_decls))

        return result

    ###########################################################################
    def gen_boiler(self):
    ###########################################################################
        all_success = True
        for sub in self._subs:
            for phys in self._physics:
                for piece in self._pieces:
                    try:
                        gen_piece(phys, sub, piece)
                    except Exception as e:
                        print("Warning: failed to generate subroutine {} piece {} for physics {}, error: {}".\
                              format(sub, piece, phys, e))
                        all_success = False

        return all_success
