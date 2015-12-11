#! /bin/csh -f

#==============================================================================
# Purpose: Regenerate Buildnml_Prestage/pop.input_data_list after
#          changing OCN_TRACER_MODULES in the $case directory
#
# Usage:   Execute regen_list.csh in the $case directory
#==============================================================================

set xmlquery_data=`./xmlquery -s JGFSEP CASEROOT SRCROOT EXEROOT OCNGRID OCN_TRACER_MODULES`
set CASEROOT=`echo $xmlquery_data | awk -F'JGFSEP' '{print $1}'`
set SRCROOT=`echo $xmlquery_data | awk -F'JGFSEP' '{print $2}'`
set EXEROOT=`echo $xmlquery_data | awk -F'JGFSEP' '{print $3}'`
set OCNGRID=`echo $xmlquery_data | awk -F'JGFSEP' '{print $4}'`
set OCN_TRACER_MODULES=`echo $xmlquery_data | awk -F'JGFSEP' '{print $5}'`

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
