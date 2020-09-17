from utils import expect

import pathlib, re

#
# Global hardcoded data
#

# piece map. maps the name of a piece of boilerplate that needs to be genereate to:
#   (filepath (relative to cxx_root))
FILEPATH, INSERT_REGEX, ID_SELF_BEGIN_REGEX, ID_SELF_END_REGEX  = range(4)
PIECES = {
    "f90_c2f_bind"  : (
        "micro_$physics_iso_c.f90",
        r"^\s*end\s+module\smicro_$physics_iso_c",
    ),
    "f90_f2c_bind"  : (
        "micro_$physics_iso_f.f90",
        r"^\s*end\s+module\smicro_$physics_iso_c",
    ),
    "cxx_c2f_bind"  : (),
    "cxx_f2c_bind"  : (),
    "cxx_func_hpp"  : (),
    "cxx_func_impl" : (),
    "cxx_data"      : (),
    "cxx_bfb_unit"  : (),
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
def get_supported_pieces():
###############################################################################
    return PIECES.keys()

###############################################################################
def get_supported_physics():
###############################################################################
    return PHYSICS.keys()

###############################################################################
def get_piece_data(piece_name, piece_data):
###############################################################################
    return PIECES[piece_name][piece_data]

###############################################################################
def get_physics_data(physics_name, physics_data):
###############################################################################
    return PHYSICS[physics_name][physics_data]

###############################################################################
def remove_comments_and_ws(contents):
###############################################################################
    """
    >>> teststr = '''
    ... module mymod
    ...   subroutine foo(a, b, &
    ...                c, d, e,&
    ... !bad
    ... &f)
    ...
    ...     real, intent(in) :: a, b, & !go
    ...                 c, d, e, f
    ...
    ...   ! hi
    ...   !hi ! there
    ... !hi ! there
    ...   end subroutine foo
    ... end module mymod
    ... '''
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
    ...                c, d, e,&
    ... !bad
    ... &f)
    ...
    ...     real, intent(in) :: a, b, & !go
    ...                 c, d, e, f
    ...
    ...   ! hi
    ...   !hi ! there
    ... !hi ! there
    ...   end subroutine foo
    ... end module mymod
    ... '''
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
def parse_f90_args(line):
###############################################################################
    """
    Given a line of fortran code declaring an argument[s], return [(argname, argtype, intent, is_array)]

    >>> parse_f90_args('integer, intent(in) :: kts, kte, kbot')
    [('kts', 'integer', 'in', False), ('kte', 'integer', 'in', False), ('kbot', 'integer', 'in', False)]
    >>> parse_f90_args('real(rtype),intent(inout ), dimension(kts:kte) :: pres,dpres,  dz ')
    [('pres', 'real', 'inout', True), ('dpres', 'real', 'inout', True), ('dz', 'real', 'inout', True)]
    >>> parse_f90_args('logical (btype), intent( in) ::do_predict_nc')
    [('do_predict_nc', 'logical', 'in', False)]
    """
    expect(line.count("::") == 1, "Expected line format 'type-info :: names' for: {}".format(line))
    metadata, names_str = line.split("::")
    names = [name.strip() for name in names_str.split(",")]
    metadata = [item.strip() for item in metadata.strip().split(",")]
    argtype = metadata[0].split("(")[0].strip()
    intent = None
    for metadatum in metadata:
        if metadatum.startswith("intent"):
            expect(intent is None, "Multiple intents in line: {}".format(line))
            intent = metadatum.split("(")[-1].rstrip(")").strip()

    is_array = "dimension" in " ".join(metadata)

    return [(name, argtype, intent, is_array) for name in names]

###############################################################################
def parse_origin(contents, subs):
###############################################################################
    """
    Returns a map of subname->[(argname, argtype, intent, is_array)]

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
    [('p3_get_tables', [('mu_r_user', 'real', 'out', True), ('vn_user', 'real', 'out', True), ('vm_user', 'real', 'out', True), ('revap_user', 'real', 'out', True)]), ('p3_init_b', [])]
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
    def gen_piece(self, sub, phys, piece):
    ###########################################################################
        db = self._get_db(phys)
        expect(sub in db, "No subroutine {} in physics {}".format(sub, phys))

    ###########################################################################
    def gen_boiler(self):
    ###########################################################################
        all_success = True
        for sub in self._subs:
            for phys in self._physics:
                for piece in self._pieces:
                    try:
                        gen_piece(sub, phys, piece)
                    except Exception as e:
                        print("Warning: failed to generate subroutine {} piece {} for physics {}, error: {}".\
                              format(sub, piece, phys, e))
                        all_success = False

        return all_success
