#! /bin/csh -f

#==============================================================================
# Purpose: Regenerate Buildnml_Prestage/pop2.input_data_list after
#          changing OCN_TRACER_MODULES in the \$case directory
#
# Usage:   Execute regen_list.csh in the \$case directory
#==============================================================================

set    srcdir       = $CODEROOT/ocn/pop2
set    my_path      = $CASEROOT/SourceMods/src.pop2

setenv OCN_PRESTAGE    TRUE
setenv INPUT_TEMPLATES $srcdir/input_templates
setenv INPUT           $EXEROOT/ocn/input


#------------------------------------------------------------
# source the \$case environment scripts
#------------------------------------------------------------
foreach file ( `\ls env_*` )
source $file
end

#------------------------------------------------------------
# regenerate the Buildnml_Prestage/pop2.input_data_list file
#------------------------------------------------------------

cat $INPUT_TEMPLATES/${OCN_GRID}_inputdata >&! $CASEROOT/Buildnml_Prestage/pop2.input_data_list
foreach module ( `echo $OCN_TRACER_MODULES` )
$srcdir/input_templates/ocn.${module}.setup.csh ccsm_prestage $CASEROOT/Buildnml_Prestage/pop2.input_data_list
end
