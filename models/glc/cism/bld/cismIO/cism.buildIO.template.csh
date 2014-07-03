# -----------------------------------------------------------------------------
# NOTE: If you are viewing this script within the bld subdirectory of the cism
# code directory, please note that this is not a complete script. Instead, it
# is embedded in a script that is created by cism.cpl7.template (in the parent
# directory). That is where some variables are defined (cism_confIOdir,
# sourcemod_dir). This is done because cism.cpl7.template has access to the
# CASEROOT and CASEBUILD environment variables, whereas this script (which is
# meant to be run as a standalone script -- NOT part of the cesm build) does
# not necessarily know the values of these variables.
#
# If you are viewing this script from within your CASE directory, then the
# above note does not apply.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# generate new glimmer _io.F90 files
# -----------------------------------------------------------------------------
cd $cism_confIOdir

foreach file (glide glint glint_mbal)
  set file_varsdef = ${file}_vars.def
  set file_ioF90 = ${file}_io.F90
  if (-f ${file_varsdef}) then
    # ---------------------------------------------------------------------------
    #  create new _io.F90 file using the glimmer python script
    # ---------------------------------------------------------------------------
    $PYTHON generate_ncvars.py $file_varsdef ncdf_template.F90.in

    if (-f ${file_ioF90}) then
      # ---------------------------------------------------------------------------
      #  compare new _io.F90 file with current version in the objdir (if it exists)
      #  if different, copy the new one to the objdir
      # ---------------------------------------------------------------------------
      cp ${file_ioF90} ${sourcemod_dir}/${file_ioF90}
    else
      # ---------------------------------------------------------------------------
      #  if new _io.F90 file not created for some reason, exit
      # ---------------------------------------------------------------------------
      echo ERROR: glimmer python script failed to produce new file: ${file_ioF90}
      exit 2
    endif

  else
    echo ERROR: missing glimmer variable definition file: ${file_varsdef}
    exit 2
  endif
end


