#! /bin/csh -f

#==============================================================================
# Purpose: Regenerate Buildnml_Prestage/pop.input_data_list after
#          changing OCN_TRACER_MODULES in the $case directory
#
# Usage:   Execute regen_list.csh in the $case directory
#==============================================================================

set CASEROOT           = `./xmlquery CASEROOT            -value`
set SRCROOT            = `./xmlquery SRCROOT             -value`
set EXEROOT            = `./xmlquery EXEROOT             -value`
set OCNGRID            = `./xmlquery OCNGRID             -value`
set OCN_TRACER_MODULES = `./xmlquery OCN_TRACER_MODULES  -value`

set srcdir          = $SRCROOT/components/pop
set my_path         = $CASEROOT/SourceMods/src.pop

set INPUT_TEMPLATES = $srcdir/input_templates
set INPUT           = $EXEROOT/ocn/input

#------------------------------------------------------------
# regenerate the $CASEROOT/Buildconf/pop.input_data_list file
#------------------------------------------------------------

cat $INPUT_TEMPLATES/${OCN_GRID}_inputdata >&! $CASEROOT/Buildnml_Prestage/pop.input_data_list
foreach module ( `echo $OCN_TRACER_MODULES` )
    $srcdir/input_templates/ocn.${module}.setup.csh ccsm_prestage $CASEROOT/Buildnml_Prestage/pop.input_data_list
end
