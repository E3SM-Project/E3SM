================================================================================
SVN $Id: ChangeLog 68792 2015-03-10 00:57:28Z mvertens $
SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/scripts/trunk_tags/scripts4_150309/ChangeLog $
================================================================================
This file describes what tags were created and why
================================================================================

================================================================================
Originator: mvertens
Date: 9 Mar 2015
Model: scripts
Version: scripts4_150309
One-line: Major scripts refactoring
	
The following new directory structure will be implemented as part of the scripts
	refactoring and introduction of the new cime/ directory

$CCSMROOT/	
  cime/  <=================
      /components
               /data_comps
               /stub_comps
               /xcpl_comps (previously dead)
      /driver_cpl
      /externals
               /gptl   (brought in as subtree in git)
               /mct    (brought in as subtree in git)
               /pio    (brought in as subtree in git)
               /CMake  (brought in as subtree in git)
               /genf90 (brought in as subtree in git) 
      /machines
      /scripts
      /share
               /csm_share
               /esmf_wrf_timemgr
               /timing
       /tools
               /cprnc
               /load_balancing_tool
               /mapping
               /unit_testing
               /validation_testing
       /utils
	      /per5lib
	      /pythonlib
	
  components/ <==============
      /aquap
      /cam
      /cice
      /cism
      /clm
      /pop2
      /rtm
      /ww3

	
	
=> SVN externals - there will be no SVN_EXTERNAL_DIRECTORIES under cime/
D       SVN_EXTERNAL_DIRECTORIES

=> The following files were changed to account for the directory refactoring 
   as well as the new testing changes that accompanied this tag	
M       create_clone
M       create_newcase
M       create_test
	
	
=> The following tests were modified, removed or added
   - See the following README file for a full documentation of test functionality
A       ccsm_utils/Testcases/README
D       ccsm_utils/Testcases/ERB_script
D       ccsm_utils/Testcases/ERH_script
D       ccsm_utils/Testcases/OEM_build.csh
D       ccsm_utils/Testcases/OEM_script
D       ccsm_utils/Testcases/P4A_script
D       ccsm_utils/Testcases/PRS_build.csh
D       ccsm_utils/Testcases/PRS_script
A       ccsm_utils/Testcases/APT_build.csh
M       ccsm_utils/Testcases/APT_script
M       ccsm_utils/Testcases/CME_build.csh
M       ccsm_utils/Testcases/CME_script
M       ccsm_utils/Testcases/ERI_script
A       ccsm_utils/Testcases/ERP_build.csh
A       ccsm_utils/Testcases/ERP_script
M       ccsm_utils/Testcases/ERS_script
R       ccsm_utils/Testcases/ERT_script
M       ccsm_utils/Testcases/ICP_build.csh
M       ccsm_utils/Testcases/ICP_script
M       ccsm_utils/Testcases/LAR_script
M       ccsm_utils/Testcases/NCK_build.csh
M       ccsm_utils/Testcases/NCK_script
M       ccsm_utils/Testcases/NCR_build.csh
M       ccsm_utils/Testcases/NCR_script
M       ccsm_utils/Testcases/NOC_build.csh
M       ccsm_utils/Testcases/NOC_script
M       ccsm_utils/Testcases/OCP_build.csh
M       ccsm_utils/Testcases/OCP_script
M       ccsm_utils/Testcases/PEA_build.csh
M       ccsm_utils/Testcases/PEA_script
M       ccsm_utils/Testcases/PEM_build.csh
M       ccsm_utils/Testcases/PEM_script
A       ccsm_utils/Testcases/PET_build.csh
M       ccsm_utils/Testcases/PET_script
M       ccsm_utils/Testcases/PFS_script
A       ccsm_utils/Testcases/PMT_build.csh
A       ccsm_utils/Testcases/PMT_script
M       ccsm_utils/Testcases/SBN_script
M       ccsm_utils/Testcases/SEQ_script
M       ccsm_utils/Testcases/SMS_script
M       ccsm_utils/Testcases/SSP_script
M       ccsm_utils/Testcases/STA_script
M       ccsm_utils/Testcases/config_tests.xml
M       ccsm_utils/Testcases/tests_build.csh

=> The testlist manager was renamed and moved to the top level scripts directory
   - Note that manage_testlists now has a new category argument to deal with the split of test
     lists as seen below - to see the new arguments issue the command
     > manage_testslists -help	 
A       manage_testlists
D       ccsm_utils/Testlistxml/manage_xml_entries
	  
=> Testlist and test mods were moved to component directories -
   - The original testlist.xml file has been split into the following component specific testlists
     exist in the component directories
        /cime/scripts/ccsm_utils/Testlistxml/testlist_allactive.xml
        /cime/driver_cpl/cesmtest/testlist_drv.xml
        /components/cam/cesmtest/testlist_cam.xml
        /components/clm/cesmtest/testlist_clm.xml
        /components/cice/cesmtest/testlist_cice.xml
        /components/pop2/cesmtest/testlist_pop.xml
        /components/cism/cesmtest/testlist_cism.xml
	
=> As a result the following has occurred in the scripts/ directory	
D       ccsm_utils/Testlistxml/testlist.xml
A       ccsm_utils/Testlistxml/testlist_allactive.xml
A       ccsm_utils/Testlistxml/testmods_dirs/allactive
	
=> The following directories were moved to $CCSMROOT/components/cam/cesmtest/testmods_dirs  
D       ccsm_utils/Testlistxml/testmods_dirs/cam
D       ccsm_utils/Testlistxml/testmods_dirs/cam/cam4_port
D       ccsm_utils/Testlistxml/testmods_dirs/cam/cam4_port/user_nl_cam
D       ccsm_utils/Testlistxml/testmods_dirs/cam/cam5_port
D       ccsm_utils/Testlistxml/testmods_dirs/cam/cam5_port/user_nl_cam
D       ccsm_utils/Testlistxml/testmods_dirs/cam/cosp
D       ccsm_utils/Testlistxml/testmods_dirs/cam/cosp/xmlchange_cmnds
	
=> The following directories were moved to $CCSMROOT/components/cism/cesmtest/testmods_dirs  
D       ccsm_utils/Testlistxml/testmods_dirs/cism
D       ccsm_utils/Testlistxml/testmods_dirs/cism/apply_to_multiinstance
D       ccsm_utils/Testlistxml/testmods_dirs/cism/apply_to_multiinstance/README
D       ccsm_utils/Testlistxml/testmods_dirs/cism/apply_to_multiinstance/shell_commands
D       ccsm_utils/Testlistxml/testmods_dirs/cism/oneway
D       ccsm_utils/Testlistxml/testmods_dirs/cism/oneway/README
D       ccsm_utils/Testlistxml/testmods_dirs/cism/oneway/xmlchange_cmnds
D       ccsm_utils/Testlistxml/testmods_dirs/cism/override_glc_frac
D       ccsm_utils/Testlistxml/testmods_dirs/cism/override_glc_frac/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/cism/override_glc_frac/user_nl_cism
D       ccsm_utils/Testlistxml/testmods_dirs/cism/test_coupling
D       ccsm_utils/Testlistxml/testmods_dirs/cism/test_coupling/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/cism/test_coupling/user_nl_cism
D       ccsm_utils/Testlistxml/testmods_dirs/cism/trilinos
D       ccsm_utils/Testlistxml/testmods_dirs/cism/trilinos/README
D       ccsm_utils/Testlistxml/testmods_dirs/cism/trilinos/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/cism/trilinos/shell_commands
D       ccsm_utils/Testlistxml/testmods_dirs/cism/trilinos/user_nl_cism
	
=> The following directories were moved to $CCSMROOT/components/clm/cesmtest/testmods_dirs  
D       ccsm_utils/Testlistxml/testmods_dirs/clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/NoVSNoNI
D       ccsm_utils/Testlistxml/testmods_dirs/clm/NoVSNoNI/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/NoVSNoNI/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/SNICARFRC
D       ccsm_utils/Testlistxml/testmods_dirs/clm/SNICARFRC/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/SNICARFRC/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/USUMB
D       ccsm_utils/Testlistxml/testmods_dirs/clm/USUMB/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/USUMB/xmlchange_cmnds
D       ccsm_utils/Testlistxml/testmods_dirs/clm/af_bias_v5
D       ccsm_utils/Testlistxml/testmods_dirs/clm/af_bias_v5/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/af_bias_v5/user_nl_datm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/af_bias_v5/xmlchange_cmnds
D       ccsm_utils/Testlistxml/testmods_dirs/clm/allActive
D       ccsm_utils/Testlistxml/testmods_dirs/clm/allActive/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ciso
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ciso/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ciso/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/cn_conly
D       ccsm_utils/Testlistxml/testmods_dirs/clm/cn_conly/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/cn_conly/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/crop
D       ccsm_utils/Testlistxml/testmods_dirs/clm/crop/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/crop/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/cropMonthOutput
D       ccsm_utils/Testlistxml/testmods_dirs/clm/cropMonthOutput/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/cropMonthOutput/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/crop_trans_f10
D       ccsm_utils/Testlistxml/testmods_dirs/clm/crop_trans_f10/README
D       ccsm_utils/Testlistxml/testmods_dirs/clm/crop_trans_f10/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/crop_trans_f10/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/crop_trans_sville
D       ccsm_utils/Testlistxml/testmods_dirs/clm/crop_trans_sville/README
D       ccsm_utils/Testlistxml/testmods_dirs/clm/crop_trans_sville/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/crop_trans_sville/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/decStart
D       ccsm_utils/Testlistxml/testmods_dirs/clm/decStart/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/decStart/xmlchange_cmnds
D       ccsm_utils/Testlistxml/testmods_dirs/clm/default
D       ccsm_utils/Testlistxml/testmods_dirs/clm/default/shell_commands
D       ccsm_utils/Testlistxml/testmods_dirs/clm/default/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/drydepnomegan
D       ccsm_utils/Testlistxml/testmods_dirs/clm/drydepnomegan/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/drydepnomegan/xmlchange_cmnds
D       ccsm_utils/Testlistxml/testmods_dirs/clm/edTest
D       ccsm_utils/Testlistxml/testmods_dirs/clm/edTest/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC/user_nl_cpl
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_changeFlags
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_changeFlags/README
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_changeFlags/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_changeFlags/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_changeFlags/xmlchange_cmnds
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_decrease
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_decrease/README
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_decrease/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_decrease/user_nl_cism
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_increase
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_increase/README
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_increase/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_increase/user_nl_cism
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_long
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_long/README
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_long/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_long/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn_reduceOutput
D       ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn_reduceOutput/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn_reduceOutput/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn_reduceOutput/xmlchange_cmnds
D       ccsm_utils/Testlistxml/testmods_dirs/clm/irrig_o3_reduceOutput
D       ccsm_utils/Testlistxml/testmods_dirs/clm/irrig_o3_reduceOutput/README
D       ccsm_utils/Testlistxml/testmods_dirs/clm/irrig_o3_reduceOutput/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/irrig_o3_reduceOutput/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/monthly
D       ccsm_utils/Testlistxml/testmods_dirs/clm/monthly/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/monthly/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/monthly/user_nl_cpl
D       ccsm_utils/Testlistxml/testmods_dirs/clm/monthly_noinitial
D       ccsm_utils/Testlistxml/testmods_dirs/clm/monthly_noinitial/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/monthly_noinitial/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/o3
D       ccsm_utils/Testlistxml/testmods_dirs/clm/o3/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/o3/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/oldhyd
D       ccsm_utils/Testlistxml/testmods_dirs/clm/oldhyd/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/pts
D       ccsm_utils/Testlistxml/testmods_dirs/clm/pts/README
D       ccsm_utils/Testlistxml/testmods_dirs/clm/pts/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/pts/xmlchange_cmnds
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLA
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLA/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLA/xmlchange_cmnds
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLB
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLB/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLB/xmlchange_cmnds
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROA
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROA/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROA/xmlchange_cmnds
D       ccsm_utils/Testlistxml/testmods_dirs/clm/reduceOutput
D       ccsm_utils/Testlistxml/testmods_dirs/clm/reduceOutput/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/rootlit
D       ccsm_utils/Testlistxml/testmods_dirs/clm/rootlit/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/rootlit/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subset
D       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subset/README
D       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subset/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetEarly
D       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetEarly/README
D       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetEarly/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetEarly/xmlchange_cmnds
D       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetLate
D       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetLate/README
D       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetLate/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetLate/xmlchange_cmnds
D       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetMid
D       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetMid/README
D       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetMid/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetMid/xmlchange_cmnds
D       ccsm_utils/Testlistxml/testmods_dirs/clm/vrtlay
D       ccsm_utils/Testlistxml/testmods_dirs/clm/vrtlay/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/vrtlay/user_nl_clm
	
=> The following directories were moved to $CCSMROOT/components/pop/cesmtest/testmods_dirs  
D       ccsm_utils/Testlistxml/testmods_dirs/pop
D       ccsm_utils/Testlistxml/testmods_dirs/pop/ecosys
D       ccsm_utils/Testlistxml/testmods_dirs/pop/ecosys/shell_commands
D       ccsm_utils/Testlistxml/testmods_dirs/pop/ecosys/user_nl_pop2
D       ccsm_utils/Testlistxml/testmods_dirs/pop/ecosys_restore_gx3v7
D       ccsm_utils/Testlistxml/testmods_dirs/pop/ecosys_restore_gx3v7/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/pop/ecosys_restore_gx3v7/user_nl_pop2
	
- The following directories were moved to $CCSMROOT/components/rtm/cesmtest/testmods_dirs  
D       ccsm_utils/Testlistxml/testmods_dirs/rtm
D       ccsm_utils/Testlistxml/testmods_dirs/rtm/rtmOff
D       ccsm_utils/Testlistxml/testmods_dirs/rtm/rtmOff/user_nl_rtm
D       ccsm_utils/Testlistxml/testmods_dirs/rtm/rtmOnFloodOnEffvelOff
D       ccsm_utils/Testlistxml/testmods_dirs/rtm/rtmOnFloodOnEffvelOff/user_nl_rtm
D       ccsm_utils/Testlistxml/testmods_dirs/rtm/rtmOnFloodOnEffvelOn
D       ccsm_utils/Testlistxml/testmods_dirs/rtm/rtmOnFloodOnEffvelOn/user_nl_rtm
D       ccsm_utils/Testlistxml/testmods_dirs/rtm/rtmOnIceOff
D       ccsm_utils/Testlistxml/testmods_dirs/rtm/rtmOnIceOff/user_nl_rtm
D       ccsm_utils/Testlistxml/testmods_dirs/rtm/rtmOnIceOn
D       ccsm_utils/Testlistxml/testmods_dirs/rtm/rtmOnIceOn/user_nl_rtm

=> The following empty directories were removed 
D       ccsm_utils/Tools/atm
D       ccsm_utils/Tools/glc
D       ccsm_utils/Tools/ice
D       ccsm_utils/Tools/lnd
D       ccsm_utils/Tools/lnd/clm
D       ccsm_utils/Tools/ocn

=> A big change was removing the translation from xml to environment variables by
   introducing new extensive use of the xmlquery command in the scripts. As a result
 - ccsm_getenv is no longer needed and has been moved to $CCSMROOT/cime/machines	
 - xml2env is no longer needed and was removed	
D       ccsm_utils/Tools/ccsm_getenv
D       ccsm_utils/Tools/xml2env

=> getTiming.csh was translated to a new perl routine getTiming which is now in Tools/
  - as part of this, the functionality in getTiming2.pl was migrated as a function
    in  getTimin
  - perf_summary.pl is not used and was removed
  - the directory timing/ was no longer needed and was removed	
	
A       ccsm_utils/Tools/getTiming
D       ccsm_utils/Tools/timing/getTiming.csh
D       ccsm_utils/Tools/timing/getTiming2.pl
D       ccsm_utils/Tools/timing/perf_summary.pl
D       ccsm_utils/Tools/timing

=> ccsm_check_locked files was translated from csh to perl and
   and was renamed to check_lockedfiles
D       ccsm_utils/Tools/ccsm_check_lockedfiles
A       ccsm_utils/Tools/check_lockedfiles

=> the following files were added/removed as part of the testing refactor to
   enable the comparison of all component history files	both for test 
   functionality and testing refactor
M       ccsm_utils/Tools/component_compare.sh
A       ccsm_utils/Tools/component_compare_move.sh
A       ccsm_utils/Tools/component_compare_test.sh
A       ccsm_utils/Tools/component_compgen_baseline.sh
M       ccsm_utils/Tools/component_generate.sh
A       ccsm_utils/Tools/component_write_comparefail.pl
A       ccsm_utils/Tools/testcase_setup
D       ccsm_utils/Tools/testcase_begin
D       ccsm_utils/Tools/testcase_end
D       ccsm_utils/Tools/testcase_env.csh
D       ccsm_utils/Tools/testcase_setup.csh

=> the following file was removed as part of the testing refactor
D       ccsm_utils/Tools/component_gen_comp
D       ccsm_utils/Tools/hist_compare.csh

=> SetupTools.pm was changed significantly
 - the function expand_env_var is gone - all expansion is now done ONLY
   on xml variables
 - a new function expand_xml_var has been introduced
M       ccsm_utils/Tools/SetupTools.pm

=> the following files were removed - their functionality was added as a
   new function in SetupTools, "create_namelist_infile"
 - note that all component xxx.buildnml scripts NO LONGER CALL either
   user_nl_add or user_nlcreate but instead are now perl scripts and call
   SetupTools::create_namelist_infile	
D       ccsm_utils/Tools/user_nl_add
D       ccsm_utils/Tools/user_nlcreate

=> config_compsets.xml was modified to bring in new waccm cam5 compsets and to 
   replace any apostrophes in the xml element attributes with &apos;	
M       ccsm_utils/Case.template/config_compsets.xml

=> config_grid.xml had a bug fix for ne30np4 in trying to reference a non-existent
   1.9x1.25 grid	
M       ccsm_utils/Case.template/config_grid.xml

=> the following variables were removed from config_definition.xml
 - MODEL_GEN_COMP	
   the following variables were added in config_definition.xml
 - POP_CPPDEFS, CICE_CPPDEFS, RUN_REFDIR
   in particular, RUN_REFDIR was introduce to be able to now have a new location
   for spun-up initial conditions different that ccsm4_init	 
M       ccsm_utils/Case.template/config_definition.xml

=> removed functino _write_env from ConfigCase.pm since are no longer translating to
   environment variables	
M       ccsm_utils/Case.template/ConfigCase.pm

=> a new main build script, cesm_build.pl, is now rewritten in perl 	
   this is in perl - but the environment for now is still needed for the makefile
   hence the wrapper script, cesm_build.csh is needed	
A       ccsm_utils/Tools/cesm_build.pl
A       ccsm_utils/Tools/cesm_build.csh
D       ccsm_utils/Tools/cesm_buildexe
D       ccsm_utils/Tools/cesm_buildstart

=> removed reference to ccsm_getenv and now call xmlquery - BUT STILL CSH
M       ccsm_utils/Tools/cesm_clean_build
M       ccsm_utils/Tools/cesm_postrun_setup
M       ccsm_utils/Tools/cesm_prerun_setup
M       ccsm_utils/Tools/cesm_prestage
M       ccsm_utils/Tools/cesm_submit

=> NOW REWRITTEN IN PERL FROM CSH
M       ccsm_utils/Tools/cesm_setup
M       ccsm_utils/Tools/check_case
M       ccsm_utils/Tools/check_input_data
M       ccsm_utils/Tools/preview_namelists

# removed reference to environment for CASEBASEID and made it an input variable	
M       ccsm_utils/Tools/compare_namelists.pl

# removed reference to XML/Lite and use xmlquery rather than xmlvars hash	
M       ccsm_utils/Tools/create_production_test
M       ccsm_utils/Tools/cs.status
	
# removed reference to absolute_path and now using perl abs_path
M       ccsm_utils/Tools/xmlchange
	
# now being called directly in ccsm_utils/Tools rather than $CASEROOT/Tools	
M       ccsm_utils/Tools/st_archive
	
# refactored to be utilized in a majority of the scripts
M       ccsm_utils/Tools/xmlquery

M       ccsm_utils/Tools/cs.submit
# The following were deleted and will be generated in postprocessing
  using upcoming tools to read the batch info.
M       ccsm_utils/Tools/tseries_generate.run
M       ccsm_utils/Tools/tseries_generate.submit

================================================================================
Originator: jshollen
Date: 27 Feb 2015
Model: scripts
Version: scripts4_150227
One-line: Updated testreporter to handle expected fails. 

M       ccsm_utils/Tools/testreporter.pl

================================================================================
Originator: santos
Date: 04 Feb 2015
Model: scripts
Version: scripts4_150204
One-line: Test f09_g16.B1850C5CN on yellowstone, mira, bluewaters

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: santos
Date: 02 Feb 2015
Model: scripts
Version: scripts4_150202
One-line: Update WACCM5 compsets with CLM4.5, CLUBB

M       ccsm_utils/Case.template/config_compsets.xml
        - All WACCM5 compsets now use CLM4.5
        - Added WACCM5+CLUBB compset (this will not actually work correctly
          until a future CAM tag fixes issues in vertical diffusion).

================================================================================
Originator: erik
Date: 30 Jan 2015
Model: scripts
Version: scripts4_150130
One-line:  Correct name of compset to FGC5L45BGC in prebeta testlist

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: 29 Jan 2015
Model: scripts
Version: scripts4_150129a
One-line: Added stampede prealpha test list 

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: 29 Jan 2015
Model: scripts
Version: scripts4_150129
One-line: Mira test list updates

M       ChangeLog

================================================================================
Originator: fvitt
Date: 27 Jan 2015
Model: scripts
Version: scripts4_150127
One-line: added chemistry compsets with CLM4.5, waccm-x with ionosphere and fixed waccm-sc compsets

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Testlistxml/testlist.xml
Add compsets:
  FSTRATMAM3L45 - CAM-Chem, CAM5 phys, CLM4.5
  BWTC5L45CCMIR2 - WACCM, CAM5 phys, CLM4.5
  FWXI - WACCMX with enhanced ionosphere
Corrected WACCM-SC compsets to use "waccm_sc" chemistry

M       ccsm_utils/Testcases/ICP_build.csh
Reverted the change introduced in scripts4_150112

================================================================================
Originator: erik
Date: 22 Jan 2015
Model: scripts
Version: scripts4_150122
One-line:  Add more CAM5-CLM4.5 BGC compsets for B and F and convert more prealpha/prebeta
           tests to use CAM5-CLM4.5 compsets

M       ccsm_utils/Case.template/config_compsets.xml - Add some new CAM5-CLM45 B and F coupled compsets:
                 BG1850C5L45BGC, FGC5L45BGC, FG1850C5L45BGC, change FGHISTCN to FGHISTC5L45BGC
M       ccsm_utils/Testlistxml/manage_xml_entries ---- Add comment that gives mv command
M       ccsm_utils/Testlistxml/testlist.xml ---------- Update prebeta and prealpha test lists
              to use more CAM5-CLM4.5-BGC compsets

================================================================================
Originator: jedwards
Date: 20 Jan 2015
Model: scripts
Version: scripts4_150120
One-line: fix an error in determining if the run completed

M ccsm_utils/Tools/cesm_postrun_setup

================================================================================
Originator: aliceb
Date: 14 Jan 2015
Model: scripts
Version: scripts4_150114
One-line: updated create_newcase help text

M       create_newcase

================================================================================
Originator: jshollen
Date: 12 Jan 2015
Model: scripts
Version: scripts4_150112c
One-line: added 3-test alpha test list for all beta test machines

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: 12 Jan 2015
Model: scripts
Version: scripts4_150112b
One-line: added T62_g16.G,GIAF to yellowstone for CICE 

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================

Originator: jshollen
Date: 12 Jan 2015
Model: scripts
Version: scripts4_150112a
One-line: putting T31_T31 fix into trunk

M       create_newcase

================================================================================
Originator: dbailey
Date: 12 Jan 2015
Model: scripts
Version: scripts4_150112
One-line: CICE5 decomps renamed

blkrobin -> sectrobin
blkcart -> sectcart

ccsm_utils/Testcases/ICP_build.csh

================================================================================
Originator: sacks
Date: 10 Jan 2015
Model: scripts
Version: scripts4_150110
One-line: Change cism testmods to work with multi-instance

user_nl_cism files need to be copied to _0001 and _0002 for multi-instance cases

This change allows NCK tests with the test_coupling testmod to pass.

A       ccsm_utils/Testlistxml/testmods_dirs/cism/apply_to_multiinstance
A       ccsm_utils/Testlistxml/testmods_dirs/cism/apply_to_multiinstance/README
A       ccsm_utils/Testlistxml/testmods_dirs/cism/apply_to_multiinstance/shell_commands
A       ccsm_utils/Testlistxml/testmods_dirs/cism/override_glc_frac/include_user_mods
A       ccsm_utils/Testlistxml/testmods_dirs/cism/test_coupling/include_user_mods
A       ccsm_utils/Testlistxml/testmods_dirs/cism/trilinos/include_user_mods

================================================================================
Originator: erik
Date: 08 Jan 2015
Model: scripts
Version: scripts4_150108
One-line: Put clm4_5/clm5_0 library into shared build area

Send bldroot and compspec into $comp.buildexe.csh script and for clm specify 
if shared library should be used (for clm4_5/clm5_0) or clm4_0 non-shared
library should be used.

M       ccsm_utils/Tools/cesm_buildexe

Only works with Machines_150108, and backup_clm40_f09_fglcmask_n02_clm4_5_1_r101
(or clm4_5_1_r104) and after.

================================================================================
Originator: sacks
Date: 30 Dec 2014
Model: scripts
Version: scripts4_141230a
One-line: tweak aux_clm tests

M       ccsm_utils/Testlistxml/testmods_dirs/clm/crop_trans_f10/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/crop_trans_sville/user_nl_clm

================================================================================
Originator: sacks
Date: 30 Dec 2014
Model: scripts
Version: scripts4_141230
One-line: add transient crop compset, add transient crop tests in aux_clm45 test list

M       ccsm_utils/Case.template/config_compsets.xml
A       ccsm_utils/Testlistxml/testmods_dirs/clm/crop_trans_f10/user_nl_clm
A       ccsm_utils/Testlistxml/testmods_dirs/clm/crop_trans_f10/include_user_mods
A       ccsm_utils/Testlistxml/testmods_dirs/clm/crop_trans_f10/README
A       ccsm_utils/Testlistxml/testmods_dirs/clm/crop_trans_f10
A       ccsm_utils/Testlistxml/testmods_dirs/clm/crop_trans_sville/user_nl_clm
A       ccsm_utils/Testlistxml/testmods_dirs/clm/crop_trans_sville/include_user_mods
A       ccsm_utils/Testlistxml/testmods_dirs/clm/crop_trans_sville/README
A       ccsm_utils/Testlistxml/testmods_dirs/clm/crop_trans_sville
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: 23 Dec 2014
Model: scripts
Version: scripts4_141223a
One-line: Remove ne120 tests from all prebeta lists save yellowstone intel, bluewaters,
remove ne240 tests from prebeta

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: sacks
Date: 23 Dec 2014
Model: scripts
Version: scripts4_141223
One-line: point clm test to new file

M       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subset/user_nl_clm

================================================================================
Originator: aliceb
Date: 19 Dec 2014
Model: scripts
Version: scripts4_141219
One-line: test list changes: remove clm4.5 tests from Janus until intel15
works, shorten titan test list, remove ESMF tests from systems without working
ESMF. 

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: aliceb
Date: 12 Dec 2014
Model: scripts
Version: scripts4_141212
One-line: bugfix for lt_archive.sh as described in bugzilla #2083

M       ccsm_utils/Tools/lt_archive.sh

================================================================================
Originator: andre
Date: 09 Dec 2014
Model: scripts
Version: scripts4_141209a
One-line: bugfix for scripts4_141208: fix error checking for userdefined config_machines.xml files.

M       ccsm_utils/Case.template/ConfigCase.pm

================================================================================
Originator: jshollen
Date: 09 Dec 2014
Model: scripts
Version: scripts4_141209
One-line: Added DTEST, ETEST to goldbach nag for cice testing. 

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: andre
Date: 08 Dec 2014
Model: scripts
Version: scripts4_141208
One-line: Configure and run cases and tests on user defined machines without changing sandboxes.

Add the ability to configure and run cases and tests on user defined machines
without manually changing every sandbox. Look in ~/.cesm for config_compilers.xml,
config_machines.xml, config_pes.xml and mkbatch.xxxx. If the user defined machine
is found in those files, print an UNSUPPORTED machine warning and use the configuration,
otherwise look for the machine configuration in the standard location.

M       ccsm_utils/Tools/SetupTools.pm
M       ccsm_utils/Tools/testcase_setup.csh
M       ccsm_utils/Case.template/ConfigCase.pm
M       create_newcase

================================================================================
Originator: mlevy
Date: 05 Dec 2014
Model: scripts
Version: scripts4_141205
One-line: POP2 changed namelist options for ecosys_restore test

M       ccsm_utils/Testlistxml/testmods_dirs/pop/ecosys_restore_gx3v7/user_nl_pop2

================================================================================
Originator: sacks
Date: 02 Dec 2014
Model: scripts
Version: scripts4_141202
One-line: Update perl5lib, rework clm testmods to use new += syntax for hist_fincl

M       SVN_EXTERNAL_DIRECTORIES
M       ccsm_utils/Testlistxml/testmods_dirs/clm/crop/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn_reduceOutput/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/irrig_o3_reduceOutput/include_user_mods
M       ccsm_utils/Testlistxml/testmods_dirs/clm/irrig_o3_reduceOutput/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/monthly/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/o3/user_nl_clm

================================================================================
Originator: jedwards
Date: 01 Dec 2014
Model: scripts
Version: scripts4_141201b
One-line: Fix issues with timing file and test completion

M  ccsm_utils/Tools/timing/getTiming.csh
M  ccsm_utils/Tools/cesm_postrun_setup



================================================================================
Originator: jedwards
Date: 01 Dec 2014
Model: scripts
Version: scripts4_141201a
One-line: Add resubmit feature for max_cpltime_step namelist option in drvseq5_1_04
               Change the cs.status script to give one test per line with sub tests indented.


M	     ccsm_utils/Tools/cs.status
M          ccsm_utils/Tools/cesm_postrun_setup


================================================================================
Originator: sacks
Date: 01 Dec 2014
Model: scripts
Version: scripts4_141201
One-line: add test to aux_clm45 test list

Add an intel version of the o3 test while the pgi test is failing, due to a pgi
compiler bug

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: sacks
Date: 30 Nov 2014
Model: scripts
Version: scripts4_141130
One-line: remove test from aux_clm45 test list

removed the temporary SMS test that was added in the last tag

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: sacks
Date: 29 Nov 2014
Model: scripts
Version: scripts4_141129
One-line: tweak aux_clm45 test list

- Change an ERH test to ERI (the ERH test is failing, due to what seems to be a
  scripts problem)

- Add a temporary SMS test to cover the functionality of the currently-failing
  ERS_Ly5 test

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: sacks
Date: 25 Nov 2014
Model: scripts
Version: scripts4_141125
One-line: add some aux_clm45 tests that include ozone

Add one test, and change the irrigOn_reduceOutput tests to also include ozone.

Note that these new tests will only work with a new version of CLM (not yet
committed, but will be committed soon).

M       ccsm_utils/Testlistxml/testlist.xml
A       ccsm_utils/Testlistxml/testmods_dirs/clm/irrig_o3_reduceOutput
A       ccsm_utils/Testlistxml/testmods_dirs/clm/irrig_o3_reduceOutput/README
A       ccsm_utils/Testlistxml/testmods_dirs/clm/irrig_o3_reduceOutput/include_user_mods
A       ccsm_utils/Testlistxml/testmods_dirs/clm/irrig_o3_reduceOutput/user_nl_clm
A       ccsm_utils/Testlistxml/testmods_dirs/clm/o3
A       ccsm_utils/Testlistxml/testmods_dirs/clm/o3/include_user_mods
A       ccsm_utils/Testlistxml/testmods_dirs/clm/o3/user_nl_clm

================================================================================
Originator: sacks
Date: 24 Nov 2014
Model: scripts
Version: scripts4_141124a
One-line: move golbach_intel aux_clm40 tests to yellowstone

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: mvertens, sacks
Date: 24 Nov 2014
Model: scripts
Version: scripts4_141124
One-line: compset and test changes for ED refactor; move golbach_intel aux_clm45 to yellowstone

These changes came from Mariana, from branch refactor_koven

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Testlistxml/testlist.xml
M       ccsm_utils/Testlistxml/testmods_dirs/clm/edTest/user_nl_clm

================================================================================
Originator: sacks
Date: 18 Nov 2014
Model: scripts
Version: scripts4_141118
One-line: misc. changes to support CISM2, rework testing of cism

*** Change available coupled compsets for CISM2: add BG1850C5L45BGCIS2, remove
    some others that likely would not be used; make CISM_USE_TRILINOS default to
    FALSE always, since this is no longer needed for running CISM2, and doesn't
    work on many machines; change default run length for TG compsets to 5 years
    instead of 10
M       ccsm_utils/Case.template/config_compsets.xml

*** Change documentation for CISM_USE_TRILINOS
M       ccsm_utils/Case.template/config_definition.xml

*** Add some grid combos that use the CISM 4 km and 20 km grids
M       ccsm_utils/Case.template/config_grid.xml

*** Major overhaul of aux_glc test list; tweak cism tests in prealpha and
    prebeta test suites; add cism-test_coupling testmod to a number of tests:
    this turns on frequent dynamics in CISM, so that we can exercise most of the
    CISM code in a few-day test run
M       ccsm_utils/Testlistxml/testlist.xml
D       ccsm_utils/Testlistxml/testmods_dirs/cism/gradient_margin_2
D       ccsm_utils/Testlistxml/testmods_dirs/cism/gradient_margin_2/user_nl_cism
A  +    ccsm_utils/Testlistxml/testmods_dirs/cism/test_coupling

*** Change test mods in accordance with the change in CISM_USE_TRILINOS default
M       ccsm_utils/Testlistxml/testmods_dirs/cism/trilinos/README
A  +    ccsm_utils/Testlistxml/testmods_dirs/cism/trilinos/shell_commands
D       ccsm_utils/Testlistxml/testmods_dirs/cism/no_trilinos
D       ccsm_utils/Testlistxml/testmods_dirs/cism/no_trilinos/README
D       ccsm_utils/Testlistxml/testmods_dirs/cism/no_trilinos/xmlchange_cmnds

================================================================================
Originator: sacks
Date: 14 Nov 2014
Model: scripts
Version: scripts4_141114
One-line: Tweak some aux_clm45 tests

I came across some notes from Sam Levis saying that, to properly test crop
restarts, the restart file should be written after at least 2 years have
elapsed. This led me to tweak the crop restart tests, some of which had been
ERS_Lm25. I lengthened one of them to satisfy Sam's point; for the other two, I
have chosen restart intervals to test restarts in the middle of the first year
and in the middle of the second year, figuring that that might trigger some
different logic.

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: sacks
Date: 12 Nov 2014
Model: scripts
Version: scripts4_141112a
One-line: Clean up the testlist using the new -cleanxml option

This removes duplicate tests, sorts entries and consolidates duplicate nodes in
the xml file

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: sacks
Date: 12 Nov 2014
Model: scripts
Version: scripts4_141112
One-line: Add -cleanxml option to manage_xml_entries

From the help message:
    This option sorts the xml file, consolidates duplicate nodes, and removes duplicate entries.

    It should be run after making any manual edits to the xml file, so that the next users
    of this file are starting from a clean state. This will minimize merge conflicts.

    Note that sorting is automatically done after other manipulations performed by this
    program, so this mode of operation only needs to be performed after making manual edits.

*** Add -cleanxml option, and do some refactoring
M       ccsm_utils/Testlistxml/manage_xml_entries

*** Add some tests for the new functionality
M       ccsm_utils/Testlistxml/test_manage_xml_entries/test_manage_xml_entries.t
A  +    ccsm_utils/Testlistxml/test_manage_xml_entries/testlist.aux_glc_subset.xml
A  +    ccsm_utils/Testlistxml/test_manage_xml_entries/testlist.aux_glc_subset_mangled.xml
A  +    ccsm_utils/Testlistxml/test_manage_xml_entries/testlist.original.15Oct2013.cleaned.xml

================================================================================
Originator: jshollen
Date: 11 Nov 2014
Model: scripts
Version: scripts4_141111
One-line: Switch ERH,ERB Ld3 to Ld7

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: andre
Date: 07 Nov 2014
Model: scripts
Version: scripts4_141107
One-line: Add aux_pop_obgc_se test suite and pop testmods
Depends: cesm_pop_2_1_20141031

A       ccsm_utils/Testlistxml/testmods_dirs/pop
A       ccsm_utils/Testlistxml/testmods_dirs/pop/ecosys
A       ccsm_utils/Testlistxml/testmods_dirs/pop/ecosys/user_nl_pop2
A       ccsm_utils/Testlistxml/testmods_dirs/pop/ecosys/shell_commands
A       ccsm_utils/Testlistxml/testmods_dirs/pop/ecosys_restore_gx3v7
A       ccsm_utils/Testlistxml/testmods_dirs/pop/ecosys_restore_gx3v7/include_user_mods
A       ccsm_utils/Testlistxml/testmods_dirs/pop/ecosys_restore_gx3v7/user_nl_pop2
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: 05 Nov 2014
Model: scripts
Version: scripts4_141105a
One-line:  Add GECO test to yellowstone pgi

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jedwards
Date: 05 Nov 2014
Model: scripts
Version: scripts4_141105
One-line:  resolves bugzilla 2084

M   ccsm_utils/Tools/archive_metadata.sh


================================================================================
Originator: sacks
Date: 31 Oct 2014
Model: scripts
Version: scripts4_141031
One-line: Allow recursive testmods and user_mods

Main change is to allow recursive testmods and user_mods (i.e., a mods directory
that includes one or more other mods directories). This is done by having an
include_user_mods text file in the mods directory, which contains paths
(relative or absolute) to other mods directories that should be applied.

I have applied this new functionality to the clm testmods.

Other misc. changes in this tag:

- use shell_commands instead of xmlchange_commands in mods directories (but
  still allow xmlchange_cmnds for backwards compatibility)

- for clm-default: apply multi-instance testmods in a more robust way that
  allows for recursive testmods

- remove a test from the aux_clm45 test list, for which there was already a
  similar test that had a testmods directory

*** Main changes
A  +    ccsm_utils/Tools/UserModsTools.pm
M       create_newcase

*** Apply changes to clm testmods
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/NoVSNoNI/include_user_mods
M       ccsm_utils/Testlistxml/testmods_dirs/clm/NoVSNoNI/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/SNICARFRC/include_user_mods
M       ccsm_utils/Testlistxml/testmods_dirs/clm/SNICARFRC/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/ciso/include_user_mods
M       ccsm_utils/Testlistxml/testmods_dirs/clm/ciso/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/cn_conly/include_user_mods
M       ccsm_utils/Testlistxml/testmods_dirs/clm/cn_conly/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/crop/include_user_mods
M       ccsm_utils/Testlistxml/testmods_dirs/clm/crop/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/cropMonthOutput/include_user_mods
M       ccsm_utils/Testlistxml/testmods_dirs/clm/cropMonthOutput/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/decStart/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/decStart/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/default/shell_commands
M       ccsm_utils/Testlistxml/testmods_dirs/clm/default/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/default/user_nl_clm_0001
D       ccsm_utils/Testlistxml/testmods_dirs/clm/default/user_nl_clm_0002
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/drydepnomegan/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/drydepnomegan/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC/include_user_mods
M       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_changeFlags/include_user_mods
M       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_changeFlags/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_changeFlags/user_nl_cpl
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_decrease/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_decrease/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_decrease/user_nl_cpl
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_increase/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_increase/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_increase/user_nl_cpl
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_long/include_user_mods
M       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_long/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_long/user_nl_cpl
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn_reduceOutput/include_user_mods
M       ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn_reduceOutput/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/monthly/include_user_mods
M       ccsm_utils/Testlistxml/testmods_dirs/clm/monthly/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/monthly_noinitial/include_user_mods
M       ccsm_utils/Testlistxml/testmods_dirs/clm/monthly_noinitial/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/monthly_noinitial/user_nl_cpl
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/pts
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLA/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLA/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLA/xmlchange_cmnds
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLB/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLB/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLB/xmlchange_cmnds
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROA/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROA/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROA/xmlchange_cmnds
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/rootlit/include_user_mods
M       ccsm_utils/Testlistxml/testmods_dirs/clm/rootlit/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subset
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetEarly/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetEarly/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetLate/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetLate/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetMid/include_user_mods
D       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetMid/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/vrtlay/include_user_mods
M       ccsm_utils/Testlistxml/testmods_dirs/clm/vrtlay/user_nl_clm

*** Unrelated change (remove test from aux_clm: see above)
M       ccsm_utils/Testlistxml/testlist.xml


================================================================================

Originator: mickelso
Date: 30 Oct 2014
Model: scripts
Version: scripts4_141030
One-line: Added a couple of prebeta tests for babbage

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================

Originator: sacks
Date: 27 Oct 2014
Model: scripts
Version: scripts4_141027a
One-line: Tweak aux_glc test list

M       ccsm_utils/Testlistxml/testlist.xml
A       ccsm_utils/Testlistxml/testmods_dirs/cism/gradient_margin_2
A       ccsm_utils/Testlistxml/testmods_dirs/cism/gradient_margin_2/user_nl_cism

================================================================================
Originator: erik
Date: 27 Oct 2014
Model: scripts
Version: scripts4_141027
One-line: Backout L45 coupled tests to CLM40 for non-existant grids

A few new coupled L45 tests were at resolutions without datasets: ne120, f45-1850
and f02. Move these tests back to using CLM40 until new datasets are provided.

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jedwards
Date: 24 Oct 2014
Model: scripts
Version: scripts4_141024c
One-line: correct cs.submit when using testroot option

M    ccsm_utils/Tools/cs.submit

================================================================================
Originator: jedwards
Date: 24 Oct 2014
Model: scripts
Version: scripts4_141024b
One-line:  Needed to fix test to match name change from CCSM to CESM and
	get timing tool to match name of timing summary file.


M            64739   ccsm_utils/Tools/timing/getTiming.csh
M            64739   ccsm_utils/Testcases/ERT_script
M            64739   ccsm_utils/Testcases/PRS_script
M            64739   ccsm_utils/Testcases/STA_script
M            64739   ccsm_utils/Testcases/ERS_script
M            64739   ccsm_utils/Testcases/SSP_script

================================================================================
Originator: jedwards
Date: 24 Oct 2014
Model: scripts
Version: scripts4_141024a
One-line: refactor cs.status - improve output format

M   ccsm_utils/Tools/cs.status

================================================================================
Originator: jedwards
Date: 24 Oct 2014
Model: scripts
Version: scripts4_141024
One-line: refactor cs.submit - now only resubmits tests whose status is not PASS, PEND or RUN

M   ccsm_utils/Tools/cs.submit


================================================================================
Originator: sacks
Date: 23 Oct 2014
Model: scripts
Version: scripts4_141023
One-line: Re-remove some aux_glc tests

Re-remove the FG and IG tests that I had removed in scripts4_141017a... that
change was lost in scripts4_141022

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: sacks
Date: 22 Oct 2014
Model: scripts
Version: scripts4_141022g
One-line: add a 4km test to the aux_glc test suite

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: sacks
Date: 22 Oct 2014
Model: scripts
Version: scripts4_141022f
One-line: allow easy creation of cases using 4-km Greenland grid

M       ccsm_utils/Case.template/config_grid.xml

================================================================================
Originator: sacks
Date: 22 Oct 2014
Model: scripts
Version: scripts4_141022e
One-line: allow 4-km Greenland grid

M       ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: sacks
Date: 22 Oct 2014
Model: scripts
Version: scripts4_141022d
One-line: change CISM2 PEA test in aux_glc testlist so that it builds

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: sacks
Date: 22 Oct 2014
Model: scripts
Version: scripts4_141022c
One-line: change description for the CISM_USE_TRILINOS xml variable

M       ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: sacks
Date: 22 Oct 2014
Model: scripts
Version: scripts4_141022b
One-line: fix DLND_CPLHIST for TG cases that point to 20TR output

Some 20TR cases had been accidentally renamed to HIST rather than 20TR when
these compsets were renamed; this led to pointers to non-existent files

M       ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: sacks
Date: 22 Oct 2014
Model: scripts
Version: scripts4_141022a
One-line: remove CISM2P/CISM2S distinction

Now that trilinos is not required in order to run CISM2 in parallel, there is no
longer much need for a CISM2S compset. So I am removing this, and removing the P
vs. S in the CISM2 compset names.

In order to get reasonable PE layouts for CISM2 compsets, this requires
Machines_141022.

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: muszala
Date: 22 Oct 2014
Model: scripts
Version: scripts4_141022
One-line: aux_clm45 testlist reorganization

M       ccsm_utils/Case.template/config_compsets.xml
-- xml doesn't like '--' in a comment line
  -<!---P compsets -- PORT (CAM's Parallel Offline Radiation Tool) -->
  +<!---P compsets PORT (CAM's Parallel Offline Radiation Tool) -->
M       ccsm_utils/Testlistxml/testlist.xml
-- yearly reorg of aux_clm45 tests.  add gnu+yellowstone, remove
-- pgi+goldbach, balance tests between machine/compiler
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/cropMonthOutput
-- adds capability to just dump monthly output to crop tests
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROB
-- not required anymore
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROB/user_nl_clm
-- not required anymore
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROB/xmlchange_cmnds
-- not required anymore
D       ccsm_utils/Testlistxml/testmods_dirs/clm/voc
-- not required anymore
D       ccsm_utils/Testlistxml/testmods_dirs/clm/voc/user_nl_clm
-- not required anymore
D       ccsm_utils/Testlistxml/testmods_dirs/clm/voc/user_nl_cpl
-- not required anymore

================================================================================
Originator: sacks
Date: 21 Oct 2014
Model: scripts
Version: scripts4_141021
One-line: allow 'xmlchange_commands' in addition to 'xmlchange_cmds' in testmods

Main purpose is to allow xmlchange_commands in addition to xmlchange_cmds in
testmods. Eventually, I'd hope to do away with xmlchange_cmnds entirely, but for
now that's supported for backwards compatibility. Motivation: cmnds is a
non-standard abbreviation, and would be easy to get wrong (also the
documentation was wrong, saying to use xmlchange_cmds).

In doing this, I have also:

(1) Removed handling of testmods from create_test, instead leveraging the
existing -user_mods_dir functionality in create_newcase, in order to remove
duplicated code.

(2) Removed the -user_testmods_dir option of create_test, which didn't appear to
be working correctly anyway (it seems like it was supposed to allow directories
outside of the main testmods area, but it didn't).

M       create_test
M       create_newcase


================================================================================
Originator: sacks
Date: 18 Oct 2014
Model: scripts
Version: scripts4_141018
One-line: add some no_trilinos tests for aux_glc test list

M       ccsm_utils/Testlistxml/testlist.xml
A       ccsm_utils/Testlistxml/testmods_dirs/cism/no_trilinos
A       ccsm_utils/Testlistxml/testmods_dirs/cism/no_trilinos/README
A       ccsm_utils/Testlistxml/testmods_dirs/cism/no_trilinos/xmlchange_cmnds

================================================================================
Originator: jshollen
Date: 17 Oct 2014
Model: scripts
Version: scripts4_141017c
One-line: fix cesm_submit for SBN, regular tests, plain cases

M       ccsm_utils/Tools/cesm_submit

================================================================================
Originator: aliceb
Date: 17 Oct 2014
Model: scripts
Version: scripts4_141017b
One-line: updates for the tseries generation support in CESM scripts

(Note: you must manually edit the BSUB -P entry in the tseries_generate.run 
script to accurately reflect your project number. This will change with 
the integration into the new batch system that Jay is working on.)

M       ccsm_utils/Tools/cesm_tseries_generator.py
M       ccsm_utils/Tools/tseries_generate.submit
M       ccsm_utils/Tools/tseries_generate.run
M       ccsm_utils/Tools/pythonlib/cesmEnvLib.py
M       ccsm_utils/Tools/cesm_postrun_setup
M       ccsm_utils/Case.template/config_archive.xml
M       ccsm_utils/Case.template/config_definition.xml
M       create_newcase

================================================================================

Originator: sacks
Date: 17 Oct 2014
Model: scripts
Version: scripts4_141017a
One-line: Remove some aux_glc tests

Removed FG tests (redundant with BG tests), and a couple of redundant IG tests

(Note: running this through manage_xml_entries also resulted in a rearrangement
of other tests)

M       ccsm_utils/Testlistxml/testlist.xml


================================================================================
Originator: jshollen
Date: 17 Oct 2014
Model: scripts
Version: scripts4_141017
One-line: Fix cesm_submit bug for regular cases

M       ccsm_utils/Tools/cesm_submit

================================================================================
Originator: jedwards
Date: 16 Oct 2014
Model: scripts
Version: scripts4_141016
One-line: rework of getTiming2.pl for compatibility with new gptl lib

M ccsm_utils/Tools/timing/getTiming2.pl


================================================================================
Originator: jedwards
Date: 15 Oct 2014
Model: scripts
Version: scripts4_141015
One-line: Update GPTL timer library, add profile_papi_enable option
                requires driver, timer and Machines tag updates.  
	
M ccsm_utils/Tools/timing/getTiming.csh
M ccsm_utils/Tools/timing/getTiming2.pl
M ccsm_utils/Case.template/config_definition.xml


================================================================================
Originator: erik
Date: 09 Oct 2014
Model: scripts
Version: scripts4_141009
One-line:  Fix several namelist-compare issues, change 20TR=>HIST, work on PDAY
           work on test-lists

Fix bugs: 2024, 2035, 2037

Move scripts4_140819_nlcomparefix to trunk.

Change compsets with 20TR to HIST. Change fixed PDAY to 2000, but change 1850PDAY 
to PINDPDAY. Fix some issues with SBN tests. Move some aux_clm45 tests to use
new test-mods to turn on drydep namelist without the megan namelist. And move
some prebeata and prealpha tests from using CLM40 to using CLM45.

Fix a bunch of issues with create_test for namelist-compare. And have namelist-compare
cases retain config-test options so the compare will be more valid.

----------- Add tests for running with drydep namelist and without megan namelist
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/drydepnomegan
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/drydepnomegan/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/drydepnomegan/xmlchange_cmnds

M       ccsm_utils/Case.template/config_compsets.xml - Rename 20TR to HIST, work on PDAY
           compsets
M       ccsm_utils/Testcases/SBN_script ----- Make sure some env variables are set
M       ccsm_utils/Testlistxml/testlist.xml - Change testlists, add more CLM45 coupled
           tests, rename 20TR compsets to HIST
M       ccsm_utils/Tools/cesm_submit -------- a fix for SBN
M       ccsm_utils/Tools/testcase_end ------- set base_name for both BASEGEN and BASECMP (bug 2035)
M       ccsm_utils/Tools/testcase_env.csh --- unset BASECMP_NAME
M       create_test --- Fix some errors, that use diagnostics found, set undefined variables.
                        Change namelist tests so they retain test options in comparison.
                        Change warn to print which seems to help a panic that happens in
                        debug mode (bug 2024). Fix namelist compare to compare to another case that
                        isn't identical (bug 2037).

================================================================================
Originator: sacks
Date: 07 Oct 2014
Model: scripts
Version: scripts4_141007
One-line: Modify IS2 tests in aux_glc testlist

(1) Remove hopper tests. I may want to add these back in later (especially to
    test gnu), but for now I'm not maintaining these.

(2) Add SMS_D_Ly1.f09_g16_gl10.TGIS2 tests that use trilinos, both for
    yellowstone-intel and titan-pgi. I may want to lengthen these tests, but for
    now these 1-year tests should be sufficient to make sure trilinos basically
    works.

(3) Change one yellowstone-intel IS2 test to yellowstone-pgi: yellowstone-pgi
    should work now that we're no longer using trilinos by default (we don't
    have trilinos for pgi on yellowstone).

M       ccsm_utils/Testlistxml/testlist.xml
A       ccsm_utils/Testlistxml/testmods_dirs/cism/trilinos
A       ccsm_utils/Testlistxml/testmods_dirs/cism/trilinos/README
A       ccsm_utils/Testlistxml/testmods_dirs/cism/trilinos/user_nl_cism

================================================================================
Originator: fvitt
Date: 03 Oct 2014
Model: scripts
Version: scripts4_141003
One-line: Added PORT and CCMI chemistry compsets and renamed FSSOA to FSTRATSOA

The new compsets:

PORT (Parallel Offline Radiation Tool) for CAM4 and CAM5 configurations
use stub surface models:
  PC4  shortname: P_2000_CAM4  longname: 2000_CAM4%PORT_SLND_SICE_SOCN_SROF_SGLC_SWAV
  PC5  shortname: P_2000_CAM5  longname: 2000_CAM5%PORT_SLND_SICE_SOCN_SROF_SGLC_SWAV

CCMI chemistry (waccm and cam-chem) REFC1 compsets:
  alias: FWMC4L40CCMIR1    shortname: F_1950-2010_CCMI_REFC1_WACCM_MA           longname: FRC1_CAM4%WCMA_CLM40%SP_CICE%PRES_DOCN%DOM_RTM_SGLC_SWAV
  alias: FWTC4L40CCMIR1    shortname: F_1950-2010_CCMI_REFC1_WACCM_TSMLT        longname: FRC1_CAM4%WTSM_CLM40%SP_CICE%PRES_DOCN%DOM_RTM_SGLC_SWAV
  alias: FTSC4L40CCMIR1    shortname: F_1950-2010_CCMI_REFC1_TROP_STRAT         longname: FRC1_CAM4%SSOA_CLM40%SP_CICE%PRES_DOCN%DOM_RTM_SGLC_SWAV
  alias: FWMC4L40CCMIR1SD  shortname: F_1975-2010_CCMI_REFC1SD_WACCM_MA         longname: SDC1_CAM4%WCMA_CLM40%SP_CICE%PRES_DOCN%DOM_RTM_SGLC_SWAV
  alias: FWTC4L40CCMIR1SD  shortname: F_1975-2010_CCMI_REFC1SD_WACCM_TSMLT      longname: SDC1_CAM4%WTSM_CLM40%SP_CICE%PRES_DOCN%DOM_RTM_SGLC_SWAV
  alias: FTSC4L40CCMIR1SD  shortname: F_1975-2010_CCMI_REFC1SD_TROP_STRAT       longname: SDC1_CAM4%SSOA_CLM40%SP_CICE%PRES_DOCN%DOM_RTM_SGLC_SWAV

CCMI chemistry (waccm and cam-chem) REFC2 compsets:
  alias: BWMC4L40CCMIR2    shortname: B_1950-2100_CCMI_REFC2_RCP6.0_WACCM_MA    longname: C2R6_CAM4%WCMA_CLM40%SP_CICE_POP2_RTM_SGLC_SWAV
  alias: BWMC4L40CCMIS2R85 shortname: B_1950-2100_CCMI_SENC2_RCP8.5_WACCM_MA    longname: C2R8_CAM4%WCMA_CLM40%SP_CICE_POP2_RTM_SGLC_SWAV
  alias: BWTC4L40CCMIR2    shortname: B_1950-2100_CCMI_REFC2_RCP6.0_WACCM_TSMLT longname: C2R6_CAM4%WTSM_CLM40%SP_CICE_POP2_RTM_SGLC_SWAV
  alias: BTSC4L40CCMIR2    shortname: B_1950-2100_CCMI_REFC2_RCP6.0_TROP_STRAT  longname: C2R6_CAM4%SSOA_CLM40%SP_CICE_POP2_RTM_SGLC_SWAV
  alias: BTSC4L40CCMIS2R45 shortname: B_2004-2100_CCMI_SENC2_RCP4.5_TROP_STRAT  longname: C2R4_CAM4%SSOA_CLM40%SP_CICE_POP2_RTM_SGLC_SWAV

Renamed FSSOA to FSTRATSOA:
  FSTRATSOA   shortname: F_2000_STRATSOA

M      ccsm_utils/Case.template/config_compsets.xml
 - added PORT and CCMI chemistry compsets and renamed FSSOA to FSTRATSOA

M      ccsm_utils/Testlistxml/testlist.xml
 - added tests for compsets listed above

A      ccsm_utils/Testlistxml/testmods_dirs/cam/cam4_port/user_nl_cam
A      ccsm_utils/Testlistxml/testmods_dirs/cam/cam4_port
A      ccsm_utils/Testlistxml/testmods_dirs/cam/cam5_port/user_nl_cam
A      ccsm_utils/Testlistxml/testmods_dirs/cam/cam5_port
 - inputs for PORT tests are specified in user_nl_cam

================================================================================
Originator: tcraig
Date: 30 Sep 2014
Model: scripts
Version: scripts4_140930
One-line: add preview_namelist call after prestaging for cice5,
 also make sure different parts of scripts are in the correct directory.

M       ccsm_utils/Tools/cesm_setup
M       ccsm_utils/Tools/cesm_buildnml
M       ccsm_utils/Tools/cesm_prestage

================================================================================
Originator: mlevy
Date: 25 Sep 2014
Model: scripts
Version: scripts4_140925
One-line: add C1DECO and G1D compsets

M       ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: sacks
Date: 24 Sep 2014
Model: scripts
Version: scripts4_140924a
One-line: another tweak to CLM's irrigation testmods

change secondary hist file to annual average

Note that this and the previous tag will change answers for CLM hist file
comparisons for any tests that use this testmods directory (currently there are
two such tests in the aux_clm45 test suite).

M       ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn_reduceOutput/user_nl_clm

================================================================================
Originator: sacks
Date: 24 Sep 2014
Model: scripts
Version: scripts4_140924
One-line: tweak CLM's irrigation testmods

add QIRRIG output to primary and secondary (1-d) history files, and use
double-precision output

M       ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn_reduceOutput/user_nl_clm

================================================================================
Originator: dbailey
Date: 18 Sep 2014
Model: scripts
Version: scripts4_140918
One-line: Remove cam5=.true. and add cam4=.true. for CICE.

Requires cice4_0_20140918

ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: jshollen
Date: 17 Sept 2014
Model: scripts
Version: scripts4_140917
One-line: create_test will set the machine automatically from text test lists.

M       create_test

================================================================================
Originator: sacks
Date: 16 Sept 2014
Model: scripts
Version: scripts4_140916c
One-line: Rename CLM_UPDATE_GLC_AREAS to GLC_TWO_WAY_COUPLING, and tweak aux_glc test list

GLC_TWO_WAY_COUPLING now has a more general purpose than the old
CLM_UPDATE_GLC_AREAS, and is needed by both cism and clm. 

Also, changed an aux_glc test to switch the value of GLC_TWO_WAY_COUPLING, and
removed an aux_glc test that has been CFAILing for a long time (a pgi esmf test
on yellowstone, which is not supported).

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Testlistxml/testlist.xml
A  +    ccsm_utils/Testlistxml/testmods_dirs/cism/oneway
M       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_changeFlags/xmlchange_cmnds


================================================================================
Originator: muszala
Date: 16 Sept 2014
Model: scripts
Version: scripts4_140916b
One-line: rename some clm tests

When moving tests from goldbach to yellowstone I misnamed some tests
as clm_aux45.  They are now correctly named aux_clm45.

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: 16 Sept 2014
Model: scripts
Version: scripts4_140916
One-line: Tests added to ensure we're testing grids/compsets in the validation guide
B2013WSCCN on Janus shortened to 3 days

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: 15 Sept 2014
Model: scripts
Version: scripts4_140915
One-line: Add STA short-term archiving test.

M       ccsm_utils/Tools/testcase_end
M       ccsm_utils/Testcases/config_tests.xml
A  +    ccsm_utils/Testcases/STA_script
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: sacks
Date: 12 Sept 2014
Model: scripts
Version: scripts4_140912
One-line: fix checks of COMPARE_MEMORY and COMPARE_THROUGHPUT in testcase_end

This fixes two problems that were introduced in scripts4_140515:

(1) The memcomp and tputcomp tests weren't ever being run

(2) For test types that did not define COMPARE_MEMORY and/or COMPARE_THROUGHPUT,
the test script was exiting abnormally before it completed. This meant that
model_gen_comp was not being run, and the TestStatus.out file was not being
copied to the baseline directory.

M       ccsm_utils/Tools/testcase_end

================================================================================
Originator: aliceb
Date: 11 Sept 2014
Model: scripts
Version: scripts4_140911a
One-line: st_archive minor fix for multi-instance restarts

M       st_archive

================================================================================


Originator: jshollen
Date: 11 Sept 2014
Model: scripts
Version: scripts4_140911
One-line: cloned hopper prebeta test list for edison

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: muszala
Date: 10 Sept 2014
Model: scripts
Version: scripts4_140910
One-line: move goldbach pgi tests to yellowstone

M       testlist.xml

================================================================================
Originator: sacks
Date: 08 Sept 2014
Model: scripts
Version: scripts4_140908
One-line: restore ChangeLog entry for scripts4_140905

M       ChangeLog

================================================================================

Originator: aliceb
Date: 05 Sept 2014
Model: scripts
Version: scripts4_140905c
One-line: create_newcase typo fix

M       create_newcase

================================================================================

Originator: aliceb
Date: 05 Sept 2014
Model: scripts
Version: scripts4_140905b
One-line: create_newcase minor changes for $CASEROOT/Tools/pythonlib

M       create_newcase

================================================================================

Originator: aliceb
Date: 05 Sept 2014
Model: scripts
Version: scripts4_140905a
One-line: first iteration of integrating the pyReshaper code into CESM scripts

 M      ccsm_utils/Tools
X       ccsm_utils/Tools/perl5lib
A  +    ccsm_utils/Tools/README.post_process
A  +    ccsm_utils/Tools/cesm_tseries_generator.py
A  +    ccsm_utils/Tools/tseries_generate.submit
A  +    ccsm_utils/Tools/tseries_generate.run
A  +    ccsm_utils/Tools/pythonlib
A  +    ccsm_utils/Tools/pythonlib/cesmEnvLib.py
M       ccsm_utils/Tools/st_archive
M       ccsm_utils/Tools/cesm_postrun_setup
M       ccsm_utils/Tools/xmlchange
M       ccsm_utils/Case.template/config_archive.xsd
M       ccsm_utils/Case.template/config_archive.xml
M       ccsm_utils/Case.template/ConfigCase.pm
MM      ccsm_utils/Case.template/config_definition.xml
 M      ccsm_utils/Testcases
M       ChangeLog
M       create_newcase

================================================================================
Originator: sacks
Date: 05 Sept 2014
Model: scripts
Version: scripts4_140905
One-line: add some comments in ProjectTools.pm

M       ccsm_utils/Tools/ProjectTools.pm

================================================================================

Originator: sacks
Date: 03 Sept 2014
Model: scripts
Version: scripts4_140903b
One-line: fix clean_build for glc

For glc, there are a lot of build things that are not in the obj
subdirectory. In particular, we need to at least remove CMakeCache.txt
and CMakeFiles in some cases (e.g., when changing the compiler - this
is needed to get a PEA test to pass). To be safe, simply blow away the
whole bld/glc directory.

M       ccsm_utils/Tools/cesm_clean_build

================================================================================
Originator: jshollen
Date: 03 Sept 2014
Model: scripts
Version: scripts4_140903a
One-line: add NCR test to yellowstone prebeta testlist for each compiler.

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: tcraig
Date: 03 Sept 2014
Model: scripts
Version: scripts4_140903
One-line: add NCR test, add "D" pecount (_P) feature

- add NCR test to compare multi-instance fully concurrent setup against
  single instance case run sequentially.  should use _P4x1D option to
  create test (ie. NCR_P4x1D.T31_g37.B1850CN.yellowstone_intel)
- minor mod to NCK to compare all three instance, not just two for passing
- add "D" feature to create newcase for pecount (_P) option.  D will
  distribute all models on unique tasks (ie. fully concurrent).

M       ccsm_utils/Testcases/NCK_script
M       ccsm_utils/Testcases/config_tests.xml
A       ccsm_utils/Testcases/NCR_script
M       ccsm_utils/Testcases/NCK_build.csh
A       ccsm_utils/Testcases/NCR_build.csh
M       create_newcase

================================================================================
Originator: sacks
Date: 02 Sept 2014
Model: scripts
Version: scripts4_140902
One-line: fix component_gen_comp to skip sharedlibroot directories

M       ccsm_utils/Tools/component_gen_comp

================================================================================
Originator: sacks
Date: 29 Aug 2014
Model: scripts
Version: scripts4_140829
One-line: Add $PROJECT xml variable for more robust specification of project /
          account numbers

These changes set a PROJECT xml variable that specifies the project / account
number for your case, for use on machines that need one. ProjectTools.pm takes
over some functionality that used to be replicated (with variations) in the
mkbatch files in Machines, for getting a project number from the environment
($PROJECT, $ACCOUNT, ~/.cesm_proj or ~/.ccsm_proj). The original implementation
of these changes was done by Pat Worley.

In addition, in create_test, I now escape any appearances of '$' in
sharedlibroot variable passed to create_newcase (this is needed to prevent an
unwanted expansion of $PROJECT). This change appears to have fixed the problem
where CESMSCRATCHROOT could not have unresolved variables.

Note that I have added PROJECT and PROJECT_REQUIRED to env_case.xml. They
fundamentally need to be locked by cesm_setup time, and currently need to be
locked at create_newcase time because they are used to create the long-term
archive batch script, which is created by create_newcase. But even putting aside
the issue of the long-term archive batch script, there did not appear to be an
appropriate place to put these variables so that they would be locked after
create_newcase but before cesm_setup (env_mach_pes.xml felt inappropriate) - and
putting them in env_case.xml will be more robust in case directories that are
currently created during cesm_setup are changed in the future to be created
during create_newcase (because $PROJECT is used for some directory paths on some
machines).

Finally, note that some machines (maybe just titan) used to have logic to allow
getting a default project if none was set. That logic has been removed, so you
now need to explicitly set a project on titan - either through the -project
argument, or by setting a $PROJECT or $ACCOUNT environment variable, or through
a ~/.cesm_proj or ~/.ccsm_proj file.

A       ccsm_utils/Tools/ProjectTools.pm
M       ccsm_utils/Case.template/config_definition.xml
M       create_clone
M       create_newcase
M       create_test

================================================================================
Originator: sacks
Date: 28 Aug 2014
Model: scripts
Version: scripts4_140828
One-line: Use $DIN_LOC_ROOT in some CLM testmods, rather than
          hard-coding yellowstone directories

M       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetEarly/README
M       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetEarly/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetLate/README
M       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetLate/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetMid/README
M       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetMid/user_nl_clm

================================================================================
Originator: mlevy
Date: 25 Aug 2014
Model: scripts
Version: scripts4_140825
One-line: Add B_ONEDIM (B1D) compset; same as C1D, but fully coupled.

M       ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: erik
Date: 14 Aug 2014
Model: scripts
Version: scripts4_140814a
One-line: Add support_level to grids and compsets, delete or add tests for configurations
          fix ERH test, add B and F compsets with CLM45BGC, change compsets for klindsay
          to use CLM45BGC

Fix several bugs:

2024 nlcompareonly option to create_test not working correctly for reporting
2019 ERH tests don't save the base env_run.xml, so have trouble when resubmitted...
2018 Failed tests in cesm1_3_beta11 needed for CLM
2005 Remove untested named compsets and grids
1999 T85_g16 has inconsistent land domain and surface datasets

M       create_test ---- Improve documentation, add notes about changes that
           cause trouble (bug 2024). Allow regular expression matching for
             -xml_compiler and -xml_category. Run testBuildSubmit for nlcompareonly
           option.
M       create_newcase - Inform user of support level and add it to -list
           option, also abort if bad aption given to -list
M       ccsm_utils/Case.template/config_compsets.xml -- Add support_level to
           compsets that can't be tested. Add B and F compsets with CLM45BGC
           change BPRP compsets for klindsay to use CLM45BGC
           compsets added: BC5L45BGC, B1850C5L45BGC, BPIPDC5L45BGC, 
                           BRCP85C5L45BGC, F1850C5L45BGC, FC5L45BGC
           compsets rm: F20TR, F20C5TRCN, F1850CNMAM3, I20TRCRUCLM45CN
           compsets mv: B1850BPRPCLM45=>B1850BPRPL45BGC, B1850BPRPC5CLM45=>B1850BPRPC5L45BGC
                        F20C5TR=>F1850PDC5L45BGC
M       ccsm_utils/Case.template/config_grid.xml --- Remove some unused
           grids: T05, T21, T85_f09_g16, f45_g37, f05_f05, ne120_f09_g16,
                  ne120_f02_t12, ne240_g16, ne240_f02_t12, f09_g16_rx1,
                  ne30_f09_g16_rx1
           Add support_level to grids
M       ccsm_utils/Testcases/ERH_script -- Save base env_run so rerunable
M       ccsm_utils/Testcases/SBN_script -- Set COMPARE_BASELINE=FALSE so
          won't do history or log file comparison, only namelist comparison
M       ccsm_utils/Testlistxml/testlist.xml - Add CAM5/CLM45BGC B and F tests
          add tests with CLM45BGC to: mira, edison, janus
          add some yellowstone_gnu tests to aux_clm45
M       ccsm_utils/Tools/cs.status ---- Add logic for nobatch, and abort if
          filename NOT correctly renamed
M       ccsm_utils/Tools/cs.submit ---- Abort if filename NOT correctly renamed
M       ccsm_utils/Tools/testcase_end - Do NOT delete the TestStatus.nlcomp file, 
          so that the test can be "rerun"

================================================================================
Originator: jedwards
Date: 14 Aug 2014
Model: scripts
Version: scripts4_140814
One-line: resolve env and xml vars used in user_nl_*

M   ccsm_utils/Tools/user_nl_add


================================================================================
Originator: jshollen
Date: 13 Aug 2014
Model: scripts
Version: scripts4_140813
One-line: create_production_test needed to write an env_test.xml file.

M       ccsm_utils/Tools/create_production_test

================================================================================
Originator: jshollen
Date: 5 Aug 2014
Model: scripts
Version: scripts4_140805
One-line: added CSL timing test list.

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: muszala
Date: 4 Aug 2014
Model: scripts
Version: scripts4_140804
One-line: add support for cruncep v5 files and add anomaly forcing tests

MM      ccsm_utils/Case.template/config_definition.xml
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/af_bias_v5
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: sacks
Date: 1 Aug 2014
Model: scripts
Version: scripts4_140801
One-line: fix component_gen_comp to allow regex metacharacters in the testid

M       ccsm_utils/Tools/component_gen_comp

================================================================================
Originator: jshollen
Date: 31 Jul 2014
Model: scripts
Version: scripts4_140731
One-line: added B1850BPRPCLM45, B1850BPRPC5CLM45 compsets for Keith Lindsay's CSL timings

M       ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: mlevy
Date: 28 Jul 2014
Model: scripts
Version: scripts4_140728
One-line: Rename C_ONED compset to C_ONEDIM and OCN_ONED env_run variable to
          OCN_ONEDIM

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: jshollen
Date: 21 Jul 2014
Model: scripts
Version: scripts4_140722
One-line: create_newcase will now print out where compsets/grids are tested, requires
perl5lib_140722, create_test will now sanely refuse to create duplicate tests. 

 M      .
M       create_test
M       SVN_EXTERNAL_DIRECTORIES
M       create_newcase

================================================================================
Originator: jshollen
Date: 21 Jul 2014
Model: scripts
Version: scripts4_140721
One-line: cs.status enhancements, will do throughput performance comparisons
as well as calculate total cost.  See script help for details.

M       ccsm_utils/Tools/cs.status

================================================================================
Originator: mlevy
Date: 14 Jul 2014
Model: scripts
Version: scripts4_140714
One-line: New compset for POP (single column); easier double precision output

M       ccsm_utils/Case.template/config_compsets.xml
        - Add C1D (C_ONED) compset
M       ccsm_utils/Case.template/config_definition.xml
        - Add OCN_ONED variable to env_run.xml, POP_TAVG_R8 to env_build.xml

================================================================================
Originator: sacks
Date: 8 Jul 2014
Model: scripts
Version: scripts4_140708
One-line: Remove untested cism-related compsets and grids, or add tests for them

*** Remove untested EGCN and EG1850CN, and unneeded EG1850CNTEST (now there are now EGCN compsets)
M       ccsm_utils/Case.template/config_compsets.xml

*** Remove untested T31_T31_gl10, f09_f09_gl10, f19_g16_gl10
M       ccsm_utils/Case.template/config_grid.xml

*** Remove tests of EG1850CNTEST from prerelease list, change some B20TRCN to
    BG20TRCN tests in prebeta, and add a BG20TRCN test on edison
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: santos
Date: 2 Jul 2014
Model: scripts
Version: scripts4_140702
One-line: Fix check_input_data verbose option.

M       ccsm_utils/Tools/check_input_data
        - Fix hack to stop printing data list files by default.

================================================================================
Originator: santos
Date: 1 Jul 2014
Model: scripts
Version: scripts4_140701
One-line: Add messing tests to prebeta.

M       ccsm_utils/Testlistxml/testlist.xml
         - Add tests for the following compsets:
           - F_1850-PDAY_CAM5
           - F_1996_WACCMX
           - F_2000_WACCM5
           - F_ADIAB_PHYS
           - F_IDEAL_PHYS

================================================================================
Originator: erik
Date: 29 Jun 2014
Model: scripts
Version: scripts4_140629
One-line: CLM script changes from clm4_6_00 to scripts trunk

Move clm scripts branch tag to trunk. Adds in ED and CLM50 compsets and tests.
Fix the issue with IGCN for cesm1_3_alpha11a (CLM_UPDATE_GLC_AREAS should be FALSE for clm40).

Fixes issue with T85 and T341 grids (bug 1999). Fix whitespace indentation in create_newcase.

Add some more tests in for CLM. Add a new testlist aux_scripts for scripts testing.

D       ccsm_utils/Testlistxml/testmods_dirs/clm/edTestGb -- Remove special test directory
              for for ED on goldbach, no longer needed.

M       ccsm_utils/Case.template/config_compsets.xml -- Change names for ED compsets a bit,
              set CLM_UPDATE_GLC_AREAS to FALSE for clm40, TRUE for clm45 only.
M       ccsm_utils/Case.template/config_grid.xml ------ Fix T85 and T341 grids (bug 1999)
M       ccsm_utils/Testlistxml/testlist.xml ----------- Add aux_script test list, add some CLM
              tests to prebeta, prealpha test lists

M       ccsm_utils/Testlistxml/testmods_dirs/clm/default/user_nl_clm_0001 - Remove BUILDHEAT/Qanth
M       ccsm_utils/Testlistxml/testmods_dirs/clm/default/user_nl_clm_0002 - Remove BUILDHEAT/Qanth
M       ccsm_utils/Testlistxml/testmods_dirs/clm/edTest/user_nl_clm ------- Remove special ED params
             file now in CLM build-namelist namelist_defaults.
M       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLA/xmlchange_cmnds --- Remove setting of ROOTPE
M       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLB/xmlchange_cmnds --- Remove setting of ROOTPE
M       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROA/xmlchange_cmnds --- Remove setting of ROOTPE
M       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROB/xmlchange_cmnds --- Remove setting of ROOTPE

M       ccsm_utils/Tools/testcase_env.csh - Use baseline_name_gen/cmp in Uppercase! (bug 1982)
M       ccsm_utils/Tools/testcase_end ----- Remove un-needed "set extension_list"

M       create_newcase ---- Change formatting (indentation)

================================================================================
Originator: tcraig
Date: 26 Jun 2014
Model: scripts
Version: scripts4_140626
One-line: Update for new cpl_seq_options, remove ocn_tight_coupling

M       ccsm_utils/Tools/check_case
M       ccsm_utils/Tools/timing/getTiming2.pl
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: jshollen
Date: 24 Jun 2014
Model: scripts
Version: scripts4_140624
One-line: New test for Mike, 2 runs, 2nd run gets half the tasks and double the 
threads, BFB flag TRUE for env_run and usr_nl_pop

M       ccsm_utils/Testlistxml/testlist.xml
M       ccsm_utils/Testcases/config_tests.xml
A       ccsm_utils/Testcases/OEM_script
A       ccsm_utils/Testcases/OEM_build.csh

================================================================================
Originator: jedwards
Date: 24 Jun 2014
Model: scripts
Version: scripts4_140624
One-line: correct test issue when running compare without Generate

M   ccsm_utils/Tools/testcase_end


================================================================================
Originator: jshollen
Date: 23 Jun 2014
Model: scripts
Version: scripts4_140623
One-line: Have xmlchange make backup copies of all files it changes,
in order to prevent file corruption due to filesystem problems. 

M       ccsm_utils/Tools/xmlchange

================================================================================
Originator: santos
Date: 28 May 2014
Model: scripts
Version: scripts4_140528a
One-line: Shorten FSD and WACCM-5 tests.

M       ccsm_utils/Testlistxml/testlist.xml
        - Add "_Ld3" to FSDSMAM, FSDSSOA, FSDWSF, and B1850W5CN tests.

================================================================================
Originator: mlevy
Date: 28 May 2014
Model: scripts
Version: scripts4_140528
One-line: Fix two bugs when using Darwin tracer module in POP2

1) create_newcase was setting OCN_TRACER_MODULE_OPT to " darwin" rather than
   "darwin" and the leading space was causing problems in other scripts

2) CDARWIN, GDARWIN, and G1850DARWIN should all be using drof rather than rtm
   for the runoff component

M       ccsm_utils/Case.template/config_compsets.xml
M       create_newcase
================================================================================
Originator: sacks
Date: 20 May 2014
Model: scripts
Version: scripts4_140520a
One-line: Rename CLM's fpftdyn to flanduse_timeseries in tests; update perl5lib

  Depends on an upcoming CLM tag (currently slated for clm4_5_75), but this
  dependency only matters if you are running the aux_clm tests.

M       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetEarly/README
M       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetEarly/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetLate/README
M       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetLate/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetMid/README
M       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetMid/user_nl_clm

*** Update perl5lib to a version that does similar renames for its unit tests
M       SVN_EXTERNAL_DIRECTORIES

================================================================================
Originator: jedwards
Date: 20 May 2014
Model: scripts
Version: scripts4_140520
One-line: Get cpl logs from case directory instead of DOUT_S_ROOT


M       create_clone
M       ccsm_utils/Tools/testcase_end
M       ccsm_utils/Testcases/ERB_script
M       ccsm_utils/Testcases/ERH_script
M       ccsm_utils/Testcases/SSP_script
M       ccsm_utils/Testcases/LAR_script
M       ccsm_utils/Testcases/ERI_script


================================================================================
Originator: jshollen
Date: 16 May 2014
Model: scripts
Version: scripts4_140516
One-line: Test script fixes.

M       create_test
M       ccsm_utils/Tools/cesm_postrun_setup
M       ccsm_utils/Tools/testcase_end

================================================================================
Originator: santos
Date: 15 May 2014
Model: scripts
Version: scripts4_140515d
One-line: Don't have check_input_data always print file lists.

M       ccsm_utils/Tools/check_input_data
        - Add verbose option to control printing of input_data_list files
          found. This output is rarely useful and frequently confuses
          users, so cesm_prestage shouldn't be producing it.

================================================================================
Originator: jedwards
Date: 15 May 2014
Model: scripts
Version: scripts4_140515c
One-line: Deleted some needed code in the previous (15a) commit

M        ccsm_utils/Tools/cesm_setup


================================================================================
Originator: santos
Date: 15 May 2014
Model: scripts
Version: scripts4_140515b
One-line: Clarify descriptions of the SSTICE variables.

M       ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: jedwards
Date: 15 May 2014
Model: scripts
Version: scripts4_140515a
One-line: Need to setenv for testcase_env.csh

M            ccsm_utils/Tools/testcase_env.csh
M            ccsm_utils/Tools/cesm_setup

================================================================================
Originator: jedwards
Date: 15 May 2014
Model: scripts
Version: scripts4_140515
One-line: The env_test.xml needs to be read for tests, but not production

M            60267   ccsm_utils/Tools/xml2env
M            60267   ccsm_utils/Tools/ccsm_getenv
M            60267   ccsm_utils/Tools/testcase_begin
M            60267   ccsm_utils/Tools/testcase_env.csh
M            60267   ccsm_utils/Tools/testcase_end
M            60267   ccsm_utils/Tools/testcase_setup.csh


================================================================================
Originator: jedwards
Date: 12 May 2014
Model: scripts
Version: scripts4_140512a
One-line: Fix error in create_test (sharedlibroot option didn't work)

M    create_test

================================================================================
Originator: jedwards
Date: 12 May 2014
Model: scripts
Version: scripts4_140512
One-line: remove bad aquap definition, fix scam test  

M           ccsm_utils/Case.template/config_compsets.xml
M           ccsm_utils/Testlistxml/testlist.xml


================================================================================
Originator: mlevy
Date: 9 May 2014
Model: scripts
Version: scripts4_140509c
One-line: Fix for bug 1882 (typo in declaration of tri-grids)

M            60101   ccsm_utils/Case.template/config_grid.xml

================================================================================
Originator: aliceb
Date: 9 May 2014
Model: scripts
Version: scripts4_140509b
One-line: updates to st_archive, config_archive.xml, create_newcase and added config_archive.xsd

M       ccsm_utils/Tools/st_archive
A       ccsm_utils/Case.template/config_archive.xsd
M       ccsm_utils/Case.template/config_archive.xml
M       create_newcase

================================================================================
Originator: mlevy
Date: 9 May 2014
Model: scripts
Version: scripts4_140509
One-line: $CASE.test would return non-zero exit status due to bad csh script

M            60088   ccsm_utils/Tools/testcase_end

================================================================================
Originator: jedwards
Date: 7 May 2014
Model: scripts
Version: scripts4_140507a
One-line: Fixes IOP tests 

M            60022   ccsm_utils/Tools/testcase_begin
M            60022   ccsm_utils/Tools/testcase_end

================================================================================
Originator: jedwards
Date: 7 May 2014
Model: scripts
Version: scripts4_140507
One-line: allow each component to build with or without omp independently

M   ccsm_utils/Tools/cesm_buildexe


================================================================================
Originator: jedwards
Date: 6 May 2014
Model: scripts
Version: scripts4_140506
One-line: clean up and rearrange cesm_postrun_setup

M     ccsm_utils/Tools/cesm_postrun_setup


================================================================================
Originator: jshollen
Date: 2 May 2014
Model: scripts
Version: scripts4_140502a
One-line: prevent variable expansion when finding CESMSCRATCHROOT in
create_test

M       create_test

================================================================================
Originator: sacks
Date: 2 May 2014
Model: scripts
Version: scripts4_140502
One-line: add and tweak tests in aux_clm45 and aux_glc test lists

M       ccsm_utils/Testlistxml/testlist.xml
M       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_changeFlags/README
M       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_changeFlags/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_changeFlags/xmlchange_cmnds
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_decrease
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_increase
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_long
M       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetEarly/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetLate/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetMid/user_nl_clm

================================================================================
Originator: mlevy
Date: 30 Apr 2014
Model: scripts
Version: scripts4_140430
One-line: Add prebeta tests for 3 new compsets

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: mlevy
Date: 29 Apr 2014
Model: scripts
Version: scripts4_140429b
One-line: Add 3 new compsets (B + RCP + ecosys)

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: jshollen
Date: 29 Apr 2014
Model: scripts
Version: scripts4_140429a
One-line: Fix SSP, ERH tests

M       ccsm_utils/Testcases/ERH_script
M       ccsm_utils/Testcases/SSP_script

================================================================================
Originator: santos
Date: 29 Apr 2014
Model: scripts
Version: scripts4_140429
One-line: Added new SC-WACCM B compsets.

M       ccsm_utils/Case.template/config_compsets.xml
        - Add new SC-WACCM compsets (1850, 55TR, RCP26, RCP45, and RCP85).
        - Change nuclear winter compset to a generic black carbon compset.
        - Correct 1850 WACCM5 use_case.
        - Correct descriptions of GEOS compsets.

ccsm_utils/Testlistxml/testlist.xml
        - Add tests for new compsets.

================================================================================
Originator: fvitt
Date: 24 Apr 2014
Model: scripts
Version: scripts4_140424
One-line: Added new compsets FSDSSOA and FSDSMAM

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: 23 Apr 2014
Model: scripts
Version: scripts4_140423a
One-line: Forgot Changelog for previous trunk tag

================================================================================
Originator: jshollen
Date: 23 Apr 2014
Model: scripts
Version: scripts4_140423
One-line: Fixed endless recursion bug in ConfigCase.pm, fixed ERI test script to work with new st_archive script. 

M       ccsm_utils/Case.template/ConfigCase.pm
M       ccsm_utils/Testcases/ERI_script

================================================================================
Originator: jshollen
Date: 18 Apr 2014
Model: scripts
Version: scripts4_140418
One-line: cfg_ref set_machine call needs to be set to yellowstone for namelist comparison tests
to work 

M       create_test

================================================================================
Originator: jshollen
Date: 16 Apr 2014
Model: scripts
Version: scripts4_140416
One-line: With Jim's help, changes to build,clean build, and test scripts to fix the 
shared builds

M       ccsm_utils/Tools/preview_namelists
D       ccsm_utils/Tools/clean_build
M       ccsm_utils/Tools/cesm_buildexe
M       ccsm_utils/Tools/cesm_clean_build
M       ccsm_utils/Testcases/CME_build.csh
M       ccsm_utils/Testcases/PEA_build.csh
M       ccsm_utils/Testcases/NCK_build.csh

================================================================================
Originator: aliceb
Date: 14 Apr 2014
Model: scripts
Version: scripts4_140414b
One-line: updates to st_archive, config_archive.xml, and config_definition.xml in support of time series generation.

M       ccsm_utils/Tools/st_archive
M       ccsm_utils/Case.template/config_archive.xml
M       ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: jedwards
Date: 14 Apr 2014
Model: scripts
Version: scripts4_140414
One-line: rename SCRATCHROOT to CESMSCRATCHROOT

M    ccsm_utils/Tools/ccsm_getenv
M    ccsm_utils/Case.template/config_definition.xml
M    create_test

================================================================================
Originator: aliceb
Date: 09 Apr 2014
Model: scripts
Version: scripts4_140409
One-line:  Updates to st_archive script usage and group definitions in config_definition.xml

M       ccsm_utils/Tools/st_archive
M       ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: aliceb
Date: 03 Apr 2014
Model: scripts
Version: scripts4_140403
One-line:  Updates to support new st_archive perl script

M       create_clone
A       ccsm_utils/Tools/st_archive
D       ccsm_utils/Tools/st_archive.sh
M       ccsm_utils/Tools/cesm_postrun_setup
M       ccsm_utils/Tools/xmlchange
A       ccsm_utils/Case.template/config_archive.xml
M       ccsm_utils/Case.template/ConfigCase.pm
M       ccsm_utils/Case.template/config_definition.xml
M       create_newcase

================================================================================
Originator: erik
Date: 02 Apr 2014
Model: scripts
Version: scripts4_140402b
One-line:  Always send MODEL_GEN_COMP

M       ccsm_utils/Tools/cesm_setup

================================================================================
Originator: jshollen
Date: 02 Apr 2014
Model: scripts
Version: scripts4_140402a
One-line:  Fixed merge issue, code cleanup

M       create_test

================================================================================
Originator: jshollen
Date: 02 Apr 2014
Model: scripts
Version: scripts4_140402
One-line:  Shared builds use new $SCRATCHROOT, multi-instance builds should now be fixed

M       create_test
M       ccsm_utils/Tools/cesm_buildexe
MM      ccsm_utils/Case.template/config_definition.xml
M       create_newcase

================================================================================
Originator: erik
Date: 30 Mar 2014
Model: scripts
Version: scripts4_140330
One-line:  Fix bug when model_gen_comp NOT used (fixes bug 1961)

Set MODEL_GEN_COMP to UNSET if not used, and always pass it, and test it
against UNSET to see if it should be used.

M       ccsm_utils/Tools/testcase_env.csh
M       ccsm_utils/Tools/cesm_setup
M       ccsm_utils/Tools/testcase_end
M       ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: jshollen
Date: 28 Mar 2014
Model: scripts
Version: scripts4_140328
One-line:  Fix shared build bug in csm_share, added proper shared build cleaning 
to .clean_build scripts. 

M       ccsm_utils/Tools/cesm_buildexe
M       ccsm_utils/Tools/cesm_clean_build
M       ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: erik
Date: 20 Mar 2014
Model: scripts
Version: scripts4_140320b
One-line:  Fix some reporting issues with -model_gen_comp option

M       ccsm_utils/Tools/testcase_end
M       ccsm_utils/Tools/component_compare.sh
M       ccsm_utils/Tools/component_gen_comp
M       ccsm_utils/Tools/component_generate.sh

================================================================================
Originator: erik
Date: 20 Mar 2014
Model: scripts
Version: scripts4_140320
One-line:  Seperate out clm4_0 tests from aux_clm to aux_clm40 and aux_clm45
           Also fix path for USUMB test and add model_gen_comp option to create_test

M       ccsm_utils/Testlistxml/testlist.xml --- split and remove aux_clm into aux_clm40 and aux_clm45
M       ccsm_utils/Testlistxml/testmods_dirs/clm/USUMB/xmlchange_cmnds

   Add -model_gen_comp option to create_test (bug 1922)
M       ccsm_utils/Tools/testcase_end
M       ccsm_utils/Tools/cesm_setup
M       ccsm_utils/Tools/testcase_setup.csh
M       create_test
M       ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: jshollen
Date: 17 Mar 2014
Model: scripts
Version: scripts4_140317
One-line:  Bug fixes for shared library builds, fixes threaded/non-threaded
and mpi-serial builds. 

M       create_test
M       ccsm_utils/Tools/cesm_buildexe

Performing status on external item at 'ccsm_utils/Tools/perl5lib':
================================================================================
Originator: jshollen
Date: 13 Mar 2014
Model: scripts
Version: scripts4_140313a
One-line:  Fix typos in resolution strings, add missing ICE & OCN DOMAiN_FILE 
entries for same per Brian, add mira prebeta test list. 

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: erik
Date: 13 Mar 2014
Model: scripts
Version: scripts4_140313
One-line:  Fix my screwup in the testlist

M   ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: erik
Date: 12 Mar 2014
Model: scripts
Version: scripts4_140312
One-line: Start adding in CLM50 compsets and tests, fix USUMB test
          so uses DIN_LOC_ROOT_CLMFORC, add more metadata to -list grids, remove WRF grid.

M       ccsm_utils/Case.template/config_grid.xml -- Remove WRF wus12 grid
M       create_newcase -- Add more metadata to -list grids (bug 1932)
M       ccsm_utils/Testlistxml/testmods_dirs/clm/USUMB/xmlchange_cmnds - set
            DIN_LOC_ROOT_CLMFORC (bug 1936)
M       ccsm_utils/Case.template/config_compsets.xml - Add CLM50 compsets
M       ccsm_utils/Testlistxml/testlist.xml ---------- Add CLM50 tests

================================================================================
Originator: jshollen
Date: 11 March 2014
Model: scripts
Version: scripts4_140311
One-line:  CESM library builds (mct, gptl, pio, csm_share) can now be shared,
needs Machines4_140311

M       create_test
M      ccsm_utils/Tools
M       ccsm_utils/Tools/ccsm_getenv
M       ccsm_utils/Tools/cesm_buildexe
MM      ccsm_utils/Case.template/config_definition.xml
M       create_newcase

================================================================================
Originator: sacks
Date: 5 March 2014
Model: scripts
Version: scripts4_140305
One-line: Modify an aux_glc test

M       ccsm_utils/Testlistxml/testmods_dirs/cism/override_glc_frac/user_nl_cism

================================================================================
Originator: sacks
Date: 4 March 2014
Model: scripts
Version: scripts4_140304
One-line: Add an aux_glc test

A       ccsm_utils/Testlistxml/testmods_dirs/cism
A       ccsm_utils/Testlistxml/testmods_dirs/cism/override_glc_frac
A       ccsm_utils/Testlistxml/testmods_dirs/cism/override_glc_frac/user_nl_cism
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: 24 Feb 2014
Model: scripts
Version: scripts4_140224a
One-line: Check for the existence of XML::LibXML in the machine's perl installation, 
warn and exit if not found. 

M       create_test
M       ccsm_utils/Case.template/ConfigCase.pm
M       create_newcase

================================================================================
Originator: erik
Date: 24 Feb 2014
Model: scripts
Version: scripts4_140224
One-line: Turn on transient CO2 streams for transient I compsets

Turn DATM_CO2_TSERIES on for datm as well as other variables needed to 
turn on transient CO2 for both historical and RCP scenarios.

M       ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: sacks
Date: 21 Feb 2014
Model: scripts
Version: scripts4_140221a
One-line: rename GLC_UPDATE_CLM_AREAS to CLM_UPDATE_GLC_AREAS and move to run_component_clm

Erik points out that, since this option is specific to CLM, this is where it belongs.

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: sacks
Date: 21 Feb 2014
Model: scripts
Version: scripts4_140221
One-line: Add GLC_UPDATE_CLM_AREAS env_run.xml variable

This will control whether CLM responds to dynamic areas coming from CISM

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: jshollen
Date: 20 Feb 2014
Model: scripts
Version: scripts4_140220
One-line: added ERS_D A,X compset tests to yellowstone and goldbach. 

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: muszala
Date: 14 Feb 2014
Model: scripts
Version: scripts4_140214b
One-line: remove _Mmpich for nag on goldbach

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: muszala
Date: 14 Feb 2014
Model: scripts
Version: scripts4_140214a
One-line: forgot to add ./xmlchange ROOTPE_ATM=0 for pts. mode runs

M       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLA/xmlchange_cmnds
M       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLB/xmlchange_cmnds
M       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROA/xmlchange_cmnds
M       ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROB/xmlchange_cmnds

================================================================================
Originator: muszala
Date: 14 Feb 2014
Model: scripts
Version: scripts4_140214
One-line: move pts. mode runs to testmods.  add ED functionality

pts. mode is being deprecated in CLM for all science use.  retain
tests for SCAM.  Need to change NTASKS* in xmlchange to keep
clm4_0 pts. mode running

! ED compsets and add CLM45CN back in for ED testing
M       ccsm_utils/Case.template/config_compsets.xml
! redo pts. mode tests to use testmods
M       ccsm_utils/Testlistxml/testlist.xml
! clean up, retab and remove all pts. mode functionality
M       create_newcase
! testmods additions for ed and pts tests
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/edTest
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/edTest/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/edTestGb
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/edTestGb/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLA
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLA/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLA/xmlchange_cmnds
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLB
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLB/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/ptsRLB/xmlchange_cmnds
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROA
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROA/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROA/xmlchange_cmnds
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROB
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROB/user_nl_clm
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/ptsROB/xmlchange_cmnds

================================================================================
Originator: jshollen
Date: 10 Feb 2014
Model: scripts
Version: scripts4_140210
One-line: cesm_setup Perl syntax fix, triggered by perl 5.18.0

M       ccsm_utils/Tools/cesm_setup

================================================================================
Originator: erik
Date: 9 Feb 2014
Model: scripts
Version: scripts4_140209
One-line: Fix argument to CLM build-namelist for CLM irrig test

M       ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn_reduceOutput/xmlchange_cmnds

================================================================================
Originator: erik
Date: 8 Feb 2014
Model: scripts
Version: scripts4_140208
One-line: Change clm tests so vsoilc is in user_nl_clm

M       ccsm_utils/Testlistxml/testmods_dirs/clm/NoVSNoNI/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/NoVSNoNI/xmlchange_cmnds
M       ccsm_utils/Testlistxml/testmods_dirs/clm/rootlit/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/rootlit/xmlchange_cmnds

================================================================================
Originator: erik
Date: 7 Feb 2014
Model: scripts
Version: scripts4_140207
One-line: Move snicarfrc setting to namelist from commandline and irrig from
          namelist to command line, switch frankfurt to goldbach

M       ccsm_utils/Testlistxml/testlist.xml -- change frankfurt to goldbach

M       ccsm_utils/Testlistxml/testmods_dirs/clm/SNICARFRC/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/SNICARFRC/xmlchange_cmnds
M       ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn_reduceOutput/user_nl_clm
A       ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn_reduceOutput/xmlchange_cmnds

================================================================================
Originator: sacks
Date: 4 Feb 2014
Model: scripts
Version: scripts4_140204
One-line: Add some history fields to some aux_clm tests

M       ccsm_utils/Testlistxml/testmods_dirs/clm/monthly/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/monthly_noinitial/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_changeFlags/user_nl_clm

================================================================================
Originator: sacks
Date: 29 Jan 2014
Model: scripts
Version: scripts4_140129
One-line: Add an aux_glc test

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: santos
Date: 24 Jan 2014
Model: scripts
Version: scripts4_140124
One-line: Add new 60 layer AMIP compset.

M       ccsm_utils/Case.template/config_compsets.xml
        - Add new compset, FAMIPC5L60.

M       ccsm_utils/Testlistxml/manage_xml_entries
        - Correct minor errors in the usage message.

M       ccsm_utils/Testlistxml/testlist.xml
        - Add test for new FAMIPC5L60 compset to Mira prebeta tests.

================================================================================
Originator: jshollen
Date: 21 Jan 2014
Model: scripts
Version: scripts4_140121
One-line: ERS_IOP4{c,p}.ne30_g16.F.yellowstone_intel changed res to ne30_ne30

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: 15 Jan 2014
Model: scripts
Version: scripts4_140115b
One-line: New method for ConfigCase to resolve all env vars, needed recursive
function. 

M       ccsm_utils/Case.template/ConfigCase.pm

================================================================================
Originator: tcraig
Date: 15 Jan 2014
Model: scripts
Version: scripts4_140115a
One-line: Add AVGHIST env variables for coupler, update check_exactrestart

check_exactrestart.pl mods modify the comparison to skip lines in
both file1 and file2 that are not comm_diag lines.  this is needed
for the new self-documenting coupler merge routines.

M       ccsm_utils/Tools/check_exactrestart.pl
M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Testcases/config_tests.xml
================================================================================
Originator: jshollen
Date: 15 Jan 2014
Model: scripts
Version: scripts4_140115
One-line: New method for ConfigCase.pm, resolves all variables in all xml
files. 

M       ccsm_utils/Case.template/ConfigCase.pm

================================================================================

Originator: erik
Date: 14 Jan 2014
Model: scripts
Version: scripts4_140114b
One-line: Remove PTCLM and link_dirtree, new USUMB test directory, add DATM_CO2_TSERIES
          remove BUILDHEAT, Qanth from CLM tests

Move transCO2_n03_scripts4_140114 branch tag to trunk.

Fix bugs: 1901, 1900, 1437, 979

M       ccsm_utils/Case.template/config_compsets.xml --- Set DATM_CO2_TSERIES (bug 979)
           when needed (comment out lines that turns it on by default for transient compsets)
           Fix IRCP85CLM45BGC so uses rcp8.5 rather than rcp4.5 (bug 1901)
M       ccsm_utils/Case.template/config_definition.xml - Add DATM_CO2_TSERIES
M       ccsm_utils/Case.template/config_grid.xml ------- Add g16_g16 grid

********** Update USUMB test data (produced by PTCLM)
M       ccsm_utils/Testlistxml/testmods_dirs/clm/USUMB/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/USUMB/xmlchange_cmnds

********** Remove PTCLM external and link_dirtree (bug 1437) from scripts
M       SVN_EXTERNAL_DIRECTORIES
D       link_dirtree

The PTCLM external is now part of CLM.

================================================================================
Originator: sacks
Date: 14 Jan 2014
Model: scripts
Version: scripts4_140114
One-line: Add tests to aux_clm test suite

Add some tests to test various edge cases for transient runs

A       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetMid
A       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetMid/user_nl_clm
A       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetMid/xmlchange_cmnds
A       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetMid/README
A       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetEarly
A       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetEarly/user_nl_clm
A       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetEarly/xmlchange_cmnds
A       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetEarly/README
A       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetLate
A       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetLate/user_nl_clm
A       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetLate/xmlchange_cmnds
A       ccsm_utils/Testlistxml/testmods_dirs/clm/tropicAtl_subsetLate/README
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: 13 Jan 2014
Model: scripts
Version: scripts4_140113
One-line: *** WRF compset, grid, and test removal ***

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: 02 Jan 2014
Model: scripts
Version: scripts4_140102b
One-line: Fix eos testlist typo

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: 02 Jan 2014
Model: scripts
Version: scripts4_140102a
One-line: Fix PEA mpi-serial build. 

M       ccsm_utils/Testcases/PEA_build.csh

================================================================================
Originator: sacks
Date: 02 Jan 2014
Model: scripts
Version: scripts4_140102
One-line: Fix improperly defined compset, I20TRCLM45BGC (bug 1869)

M       ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: jshollen
Date: 31 Dec 2013
Model: scripts
Version: scripts4_131231
One-line: removed SEQ_PFC tests from yellowstone prebeta lists, the
SEQ_IOP_PFC tests already create these

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: 30 Dec 2013
Model: scripts
Version: scripts4_131230
One-line: added eos prebeta testlist. 

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jedwards
Date: 03 Dec 2013
Model: scripts
Version: scripts4_131203
One-line: added support for mpilib in config_compilers.xml

M  ccsm_utils/Tools/cesm_setup
M  ccsm_utils/Tools/SetupTools.pm

================================================================================
Originator: jshollen
Date: 27 Nov 2013
Model: scripts
Version: scripts4_131127
One-line: added ability to allow/disallow pe layouts for specific tests in 
config_pes.xml

M       ccsm_utils/Case.template/ConfigCase.pm
M       create_newcase

================================================================================
Originator: sacks
Date: 26 Nov 2013
Model: scripts
Version: scripts4_131126a
One-line: Fix typo in CLM aux test list

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: sacks
Date: 26 Nov 2013
Model: scripts
Version: scripts4_131126
One-line: Tweak aux clm test list

Add a testmods directory for CLM, and change the testmods directories used by
two aux_clm tests (both IG tests).

A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_changeFlags
M  +    ccsm_utils/Testlistxml/testmods_dirs/clm/glcMEC_changeFlags/user_nl_clm
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: fvitt
Date: 22 Nov 2013
Model: scripts
Version: scripts4_131122
One-line: Update to lt_archive for pleiades

M  ccsm_utils/Tools/lt_archive.sh
   -  added "copy_dirs_local" mode for pleiades -- this is more efficient
      then the scp method.  This mode can be used if archiving is done via
      the NAS data analysis nodes which have both the /nobackup and lou disks
      mounted locally.

================================================================================
Originator: erik
Date: 19 Nov 2013
Model: scripts
Version: scripts4_131119
One-line: Also copy SourceMods for -user_mods_dir option

M       create_newcase -- also copy files under SourceMods if they exist and
           are in the new case. And print what doing under the verbose option.

================================================================================
Originator: erik
Date: 18 Nov 2013
Model: scripts
Version: scripts4_131118
One-line: add -user_mods_dir option to create_newcase

See bug 1868 for a full discussion this. But, adding an option to look for
user_nl_* files and xmlchange_cmnds file to execute to create_newcase similar
to the user_testmod_dir option in create_test, but allowing you to do this 
as a general user setting up a case. Eventually, we'll change create_test to
make use of the option in create_newcase rather than replicating it.

M       create_newcase -- add -user_mods_dir option

================================================================================
Originator: jedwards
Date: 14 Nov 2013
Model: scripts
Version: scripts4_131114
One-line: add netcdf4p and 4c to IOP tests 

M            ccsm_utils/Tools/testcase_begin
M            ccsm_utils/Testlistxml/testlist.xml


================================================================================
Originator: jshollen
Date: 13 Nov 2013
Model: scripts
Version: scripts4_131113a
One-line: Removed doc directory from scripts, now an external. 

D       doc

================================================================================
Originator: jshollen
Date: 13 Nov 2013
Model: scripts
Version: scripts4_131113
One-line: manage_xml_entries file path bug fix, script should now be usable
from any directory

M       ccsm_utils/Testlistxml/manage_xml_entries

================================================================================
Originator: jshollen
Date: 07 Nov 2013
Model: scripts
Version: scripts4_131107b
One-line: manage_xml_entries test add bug fix, doc fixes. Found stray nldir in
testlist.xml

manage_xml_entries was not adding tests due to incomplete XPath query. 

M       ccsm_utils/Testlistxml/manage_xml_entries
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: sacks
Date: 07 Nov 2013
Model: scripts
Version: scripts4_131107a
One-line: change default testid to include yymmdd; tweak CLM aux test list

*** Change default testid from hhmmss to yymmdd-hhmmss
M       create_test

*** Change CME_Ly5.f10_f10.I1850CLM45BGC.yellowstone_intel.clm-default to
    CME_Ly4.f10_f10.I1850CLM45BGC.yellowstone_intel.clm-monthly to avoid running
    out of time; add
    ERS_D_Mmpich.f10_f10.ICLM45BGCCROP.frankfurt_nag.clm-allActive
A       ccsm_utils/Testlistxml/testmods_dirs/clm/allActive
A       ccsm_utils/Testlistxml/testmods_dirs/clm/allActive/user_nl_clm
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: 07 Nov 2013
Model: scripts
Version: scripts4_131107
One-line: replace query_tests, add/make_xml_tests with manage_xml_entries, add
NOC test to yellowstone intel prebeta

A  + ccsm_utils/Testlistxml/test_manage_xml_entries
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/testlist.eosadded.xml
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/frankfurtnagprealpha.orig.txt
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/testlist.5by5amazonremoved.15Oct2013.xml
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/aux_clm_short.original.txt
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/testlist.testwithwithouttestmods.xml
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/testlist.xcompsetremoved.15Oct2013.xml
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/testlist.nagremoved.xml
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/testmods_dupes_test
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/testlist.testmods.clmdecStartremoved.xml
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/test_manage_xml_entries.t
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/frankfurt_nag_prealpha
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/testlist.fakecompsetXA.added.xml
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/testlist.olympusremoved.xml
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/testlist.auxclmshortremoved.15Oct2013.xml
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/aux_clm_short.07Oct2013
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/testlist.testmodstesting.xml
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/testlist.SEQ_PFCremoved.15Oct2013.xml
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/testlist.original.15Oct2013.xml
A  + ccsm_utils/Testlistxml/test_manage_xml_entries/testlist.adds.existing.xml
D    ccsm_utils/Testlistxml/add_xml_entries
A  + ccsm_utils/Testlistxml/manage_xml_entries
D    ccsm_utils/Testlistxml/make_xml_entries
M    ccsm_utils/Testlistxml/testlist.xml
A  + ccsm_utils/Testlistxml/README
M    ccsm_utils/Testcases
D    query_tests
================================================================================
Originator: tcraig
Date: 06 Nov 2013
Model: scripts
Version: scripts4_131106
One-line: fix NCK test and add NOC test

M       ccsm_utils/Testcases/NCK_script
A       ccsm_utils/Testcases/NOC_script
M       ccsm_utils/Testcases/config_tests.xml
A       ccsm_utils/Testcases/NOC_build.csh

================================================================================
Originator: mvertens
Date: 30 Oct 2013
Model: scripts
Version: scripts4_131030
One-line: merged in changes from scripts/branches/clm_controlMod_cpp (Ben Andre)

This removes all CPP flags from CLM and replaces them with run-time namelist settings	

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Testlistxml/testmods_dirs/clm/rootlit/xmlchange_cmnds
M       ccsm_utils/Testlistxml/testmods_dirs/clm/NoVSNoNI/xmlchange_cmnds
M       ccsm_utils/Testlistxml/testmods_dirs/clm/SNICARFRC/xmlchange_cmnds
M       ccsm_utils/Testcases/config_tests.xml
M       ccsm_utils/Testcases/SSP_script

================================================================================
Originator: jshollen 
Date: 28 Oct 2013
Model: scripts
Version: scripts4_131028
One-line: single-test create_test bug fix, namelistCompare should be using the
casebaseid

M       create_test

================================================================================
Originator: jedwards (via Cecile Hannay)
Date: 24 Oct 2013
Model: scripts
Version: scripts4_131024
One-line: 1. Change in the initial condition for FAMIP_CAM5_CN at ne30 (to use a better spunup condition)
2. Change in the initial condition for F_1850-PDAY_CAM5 at ne30 and f09 (to use a better spunup condition)
3. Add time-varying PFT from RCP4.5 in F_1850-PDAY_CAM5  ( to be more consistent with what is made in the atmosphere).

	M ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: sacks
Date: 18 Oct 2013
Model: scripts
Version: scripts4_131018
One-line: add a test to aux_clm test list

Add

ERS_Ly5.f10_f10.ICLM45BGCCROP.yellowstone_intel.clm-irrigOn_reduceOutput

which catches a restart problem that has gone undetected until now.

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: sacks
Date: 3 Oct 2013
Model: scripts
Version: scripts4_131003
One-line: tweak aux_clm test list

Turn irrigation on for a long test instead of a short one (because for
the short test, it's likely that irrigation never actually happens);
in particular, use irrigation for this test:
PET_P15x2_Ly3.f10_f10.ICLM45BGCCROP.yellowstone_intel.clm-irrigOn_reduceOutput
instead of:
ERI_N2.f19_g16.ICRUCLM45BGCCROP.yellowstone_intel

A  +    scripts/ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn_reduceOutput
M  +    scripts/ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn_reduceOutput/user_nl_clm
D       scripts/ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn
D       scripts/ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn/user_nl_clm
D       scripts/ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn/user_nl_clm_0001
D       scripts/ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn/user_nl_clm_0002
M       scripts/ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: sacks
Date: 2 Oct 2013
Model: scripts
Version: scripts4_131002
One-line: add user_nl_clm files for use in multi-instance tests

Note that this changes answers for any multi-instance test that uses
the clm-irrigOn testmods directory; it changes behavior (but not
answers) for multi-instance tests using the clm-default testmods directory.

A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/default/user_nl_clm_0001
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/default/user_nl_clm_0002
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn/user_nl_clm_0001
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn/user_nl_clm_0002

================================================================================
Originator: jedwards
Date: 1 Oct 2013
Model: scripts
Version: scripts4_131001a
One-line: fix namelist compare in create_test (sometimes tested wrong case)

M     create_test


================================================================================
Originator: sacks
Date: 1 Oct 2013
Model: scripts
Version: scripts4_131001
One-line: Remove some duplicate aux_clm edison tests

Six tests were in both the aux_clm_ys_intel and aux_clm_ys_pgi test lists. In
all cases, I removed the test from the aux_clm_ys_intel list.

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: sacks
Date: 30 Sep 2013
Model: scripts
Version: scripts4_130930a
One-line:  add aux_clm test lists for edison

Add temporary aux_clm test lists for edison, for while yellowstone is
down. (These were added in a somewhat kludgey way, since they are
intended to be temporary.) Also remove some duplicated sections in
testlist.xml.

M       ccsm_utils/Testlistxml/testlist.xml


================================================================================
Originator: sacks
Date: 30 Sep 2013
Model: scripts
Version: scripts4_130930
One-line:  add an estimated run cost associated with STOP_OPTION

M       ccsm_utils/Tools/cesm_setup

================================================================================
Originator: erik
Date: 29 Sep 2013
Model: scripts
Version: scripts4_130929
One-line:  Update PTCLM, add more I1PT compsets, fix PRS_build add reduceOutput to most frankfurt I tests

Update PTCLM

M       ccsm_utils/Case.template/config_compsets.xml - Add BGC, and CN I1PT compsets
M       ccsm_utils/Testlistxml/testlist.xml ---------- Add reduceOutput to most frankfurt I tests
M       ccsm_utils/Testcases/PRS_build.csh ----------- Reset build to first PE layout after done

================================================================================
Originator: mai
Date: 25 Sep 2013
Model: scripts
Version: scripts4_130925
One-line: add option to keep data on disk when running long-term archiver

M            51592   ccsm_utils/Tools/lt_archive.sh
M            51592   ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: sacks
Date: Sep 24 2013
Model: scripts
Version: scripts4_130924
One-line: Fix PEM test

Fix PEM test so that it passes the first time around, rather than requiring
resubmission. Previously, the job script was incorrectly set up to use the
smaller number of processors, rather than the larger number.

M       ccsm_utils/Testcases/PEM_build.csh

================================================================================
Originator: jshollen
Date: Sep 20 2013
Model: scripts
Version: scripts4_130920
One-line: Move frankfurt nag ERI to ERS, move ne30_g16.F1850 to ne30_ne30 
for yellowstone. 

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jedwards
Date: Sep 19 2013
Model: scripts
Version: scripts4_130919
One-line: Move taskmaker.pl to Machines

D      ccsm_utils/Tools/taskmaker.pl
M      create_newcase

================================================================================
Originator: santos
Date: Sep 18 2013
Model: scripts
Version: scripts4_130918
One-line: Add new SCAM compset, ARM97 with CAM5.

M       ccsm_utils/Case.template/config_compsets.xml
         - Add new compset, F_ARM97_SCAM5.

M       ccsm_utils/Testlistxml/testlist.xml
         - Change NAG SCAM test to be a prealpha test of the new compset,
           and also test it with PGI on yellowstone.

M       create_test
         - Fix some typos in the help message.

================================================================================
Originator: jshollen
Date: Sep 17 2013
Model: scripts
Version: scripts4_130917
One-line: add bluewaters test list, minor fix in ConfigCase.pm

M       ccsm_utils/Case.template/ConfigCase.pm
M       ccsm_utils/Testlistxml/add_xml_entries
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: sacks
Date: Sep 16 2013
Model: scripts
Version: scripts4_130916
One-line: Rework PET and PEM tests in aux_clm test list

Get better test coverage of PET tests, and to a lesser extent PEM tests, in the
aux_clm test list. Explicitly specify PE layout of PET tests, both to get more
reasonable PE layouts and to force multiple threads even if config_pes says not
to. Include some longer PET tests.

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: erik
Date: Sep 15 2013
Model: scripts
Version: scripts4_130915
One-line: Add PRS test that does restart testing from startup with changed PE count

PRS test halves the number of tasks and threads and does a restart test. Thus
making sure that writing and read of restart files is NOT sensitive to PE count.

M       ccsm_utils/Testcases/config_tests.xml
A       ccsm_utils/Testcases/PRS_script
A       ccsm_utils/Testcases/PRS_build.csh

================================================================================
Originator: jshollen
Date: Sep 13 2013
Model: scripts
Version: scripts4_130913
One-line: add_xml_entry: now handles testmods correctly, now sorts testlist.xml

D       ccsm_utils/Testlistxml/sort_xml_entries
M       ccsm_utils/Testlistxml/add_xml_entries

================================================================================
Originator: erik
Date: Sep 12 2013
Model: scripts
Version: scripts4_130912
One-line: Fix parsing bug, remove rootlitfrac from NoVSNoNI test
          correct and add more settings for USUMB test case

M       create_test
MM      ccsm_utils/Testlistxml/testmods_dirs/clm/USUMB/xmlchange_cmnds
M       ccsm_utils/Testlistxml/testmods_dirs/clm/NoVSNoNI/user_nl_clm
 M      ccsm_utils/Testlistxml/testmods_dirs/clm/decStart/xmlchange_cmnds

================================================================================
Originator: erik
Date: Sep 11 2013
Model: scripts
Version: scripts4_130911b
One-line: Fix NoVSNoNI and rootlit tests, make sure all numa and smallville tests 
          use compsets with CROP, update PTCLM

M       ccsm_utils/Testlistxml/testmods_dirs/clm/rootlit/xmlchange_cmnds
M       ccsm_utils/Testlistxml/testmods_dirs/clm/NoVSNoNI/xmlchange_cmnds
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: Sep 11 2013
Model: scripts
Version: scripts4_130911
One-line: fixed testlist, had deleted all testmods inadvertently.  (use of add_xml_entries had deleted them)

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: Sep 10 2013
Model: scripts
Version: scripts4_130910
One-line: add frankfurt nag prebeta testlist. 

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: Sep 9 2013
Model: scripts
Version: scripts4_130909
One-line: add frankfurt nag prealpha testlist. 

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: tcraig
Date: Sep 5 2013
Model: scripts
Version: scripts4_130905
One-line: merge cplupa branch

- remove sno component
- add support for new glc-ocn/ice and rof-ice coupling
- update timing tool for driver modifications
- echo command line to README.case for create_test
svn merge $SVNREPO/scripts/trunk_tags/scripts4_130830a $SVNREPO/scripts/branch_tags/cplupa_tags/cplupa_n04_scripts4_130830a

M       create_test
M       ccsm_utils/Tools/timing/getTiming2.pl
M       ccsm_utils/Tools/cesm_buildexe
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: sacks
Date: Aug 30 2013
Model: scripts
Version: scripts4_130830a
One-line: tweak aux_clm test list

ERI.f19_g16.IG1850CN.yellowstone_pgi -> ERI_D.f19_g16.IG1850CN.yellowstone_pgi
- so we have a _D test for IG CLM4, given that the frankfurt one is CFAILing

PET_D_P1x30.ne30_g16.ICN.yellowstone_intel -> PET_D_P4x30.ne30_g16.ICN.yellowstone_intel
- 1x30 was running out of wallclock time

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: Aug 30 2013
Model: scripts
Version: scripts4_130830
One-line: cesm_prestage quoting fix, create_test now properly handles no xml
tests found. 

M       create_test
M       ccsm_utils/Tools/cesm_prestage

================================================================================
Originator: muszala
Date: Aug 29 2013
Model: scripts
Version: scripts4_130829
One-line: add _Mmpich to frankfurt nag tests and reduceOutput option for clm

A       ccsm_utils/Testlistxml/testmods_dirs/clm/reduceOutput
A       ccsm_utils/Testlistxml/testmods_dirs/clm/reduceOutput/user_nl_clm
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: santos
Date: Aug 28 2013
Model: scripts
Version: scripts4_130828a
One-line: Print short names from "create_newcase -list compsets".

M       create_newcase
         - Print compset short names in list_compsets.

================================================================================
Originator: santos
Date: Aug 28 2013
Model: scripts
Version: scripts4_130828
One-line: Clarify message from cesm_prestage.

M       ccsm_utils/Tools/cesm_prestage
         - Clarify message so that it doesn't seem to say that input_data_list
           files are missing.

================================================================================
Originator: muszala
Date: Aug 16 2013
Model: scripts
Version: scripts4_130816
One-line: remove and modify tests so they work with new ch4 params files in clm

D       ccsm_utils/Testlistxml/testmods_dirs/clm/ch4_set2_ciso
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ch4_set2_ciso/user_nl_clm
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ch4_set3_pftroot
D       ccsm_utils/Testlistxml/testmods_dirs/clm/ch4_set3_pftroot/user_nl_clm
M       ccsm_utils/Testlistxml/testmods_dirs/clm/rootlit/user_nl_clm
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: erik
Date: Aug 09 2013
Model: scripts
Version: scripts4_130809b
One-line: Do xmlchange_cmds BEFORE cesm_setup, remove some compsets move into tests

Remove SNCR from compsets (never worked anyway). And remove NoVS and NoVSNoNi
from compsets and move into testmods_dirs. For USUMB/decStart tests, make sure
./xmlchange has "./" in front of it. Work on formatting (lining up columns)
in testlist.xml

M       create_test
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Testlistxml/testlist.xml
M       ccsm_utils/Testlistxml/testmods_dirs/clm/USUMB/xmlchange_cmnds
M       ccsm_utils/Testlistxml/testmods_dirs/clm/decStart/xmlchange_cmnds

A       ccsm_utils/Testlistxml/testmods_dirs/clm/rootlit/xmlchange_cmnds
A  +    ccsm_utils/Testlistxml/testmods_dirs/clm/NoVSNoNI
A       ccsm_utils/Testlistxml/testmods_dirs/clm/NoVSNoNI/xmlchange_cmnds
A       ccsm_utils/Testlistxml/testmods_dirs/clm/SNICARFRC
A       ccsm_utils/Testlistxml/testmods_dirs/clm/SNICARFRC/user_nl_clm
A       ccsm_utils/Testlistxml/testmods_dirs/clm/SNICARFRC/xmlchange_cmnds


================================================================================
Originator: jedwards
Date: Aug 09 2013
Model: scripts
Version: scripts4_130809
One-line: adjust aprun output from taskmaker with numa flags for intel compiler (on edison)

M ccsm_utils/Tools/taskmaker.pl


================================================================================
Originator: sacks
Date: Aug 02 2013
Model: scripts
Version: scripts4_130802a
One-line: tweak clm aux test lists, fix make_xml_entries

Change / add the following clm aux tests (< is old, > is new):

*** add / tweak some threading tests
< PET.f10_f10.I1850CLM45BGC.yellowstone_pgi.clm-default
< PET.f10_f10.ICLM45BGC.yellowstone_intel
> PET_PT_D.f10_f10.I1850CLM45BGC.yellowstone_pgi.clm-ciso
> PET_PT.f10_f10.ICLM45BGC.yellowstone_intel
> PET_PT_D.f19_g16.ICLM45BGCDVCROP.yellowstone_pgi.clm-crop

*** c13c14only was the same as ciso, so I deleted it
< ERS.f10_f10.I1850CLM45BGC.frankfurt_nag.clm-c13c14only
> ERS.f10_f10.I1850CLM45BGC.frankfurt_nag.clm-ciso

*** CISM hasn't been ported to NAG, so move IG tests elsewhere; note
    that the frankfurt pgi test is currently failing for the same
    reason as a bunch of other frankfurt pgi tests (a netcdf problem)
< ERS_D.f19_g16.IGRCP26CN.frankfurt_nag
< ERS_Lm3.f19_g16.IGRCP60CN.frankfurt_nag
> ERS_D.f19_g16.IGRCP26CN.frankfurt_pgi
> ERS_Lm3.f19_g16.IGRCP60CN.yellowstone_intel




M       ccsm_utils/Testlistxml/testlist.xml
D       ccsm_utils/Testlistxml/testmods_dirs/clm/c13c14only
D       ccsm_utils/Testlistxml/testmods_dirs/clm/c13c14only/user_nl_clm

*** work correctly with testmod dir
M       ccsm_utils/Testlistxml/make_xml_entries

================================================================================
Originator: erik
Date: Aug 02 2013
Model: scripts
Version: scripts4_130802
One-line: Sync up testlist and compsets for "I" compsets

Make sure every "I" compset has a test for it (except I20TRCRUCLM45CN).
And remove unneeded "I" compsets. Add in a 1x1_US-UMB test and a SNICAR_FRC test for aux_clm.
Also add a 1x1_US-UMB test for prebeta, and prebeta test on yellowstone_intel for ICLM45CRUBGC. 
Make a IG1850, f19_g16 test on yellowstone_pgi be for IG1850CN instead.

Consistently use alias names for "I" compset tests. A few were using snames.
(move I_2000_CRUFRC_CN ICRUCN [ERS_D hcru], move I_2000_CRUFRC_CLM45_BGC to ICRUCLM45BGC [ERS_D hcru], 
I_1948-2004 to I4804[ERB ne30])

Update PTCLM version.

M       ccsm_utils/Case.template/config_compsets.xml
A       ccsm_utils/Testlistxml/testmods_dirs/clm/USUMB
A       ccsm_utils/Testlistxml/testmods_dirs/clm/USUMB/user_nl_clm
A       ccsm_utils/Testlistxml/testmods_dirs/clm/USUMB/xmlchange_cmnds
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: Jul 31 2013
Model: scripts
Version: scripts4_130731
One-line: testreporter bug fix in HTTP POST method. 

M       ccsm_utils/Tools/testreporter.pl

================================================================================
Originator: muszala
Date: Jul 30 2013
Model: scripts
Version: scripts4_130730
One-line: Reorg. of aux_clm tests.  Add some clm test parameters to testmods.

M       ccsm_utils/Case.template/config_compsets.xml
A       ccsm_utils/Testlistxml/testmods_dirs/clm/c13c14only
A       ccsm_utils/Testlistxml/testmods_dirs/clm/c13c14only/user_nl_clm
A       ccsm_utils/Testlistxml/testmods_dirs/clm/decStart
A       ccsm_utils/Testlistxml/testmods_dirs/clm/decStart/user_nl_clm
A       ccsm_utils/Testlistxml/testmods_dirs/clm/decStart/xmlchange_cmnds
A       ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn
A       ccsm_utils/Testlistxml/testmods_dirs/clm/irrigOn/user_nl_clm
M       ccsm_utils/Testlistxml/testlist.xml
	
================================================================================
Originator: jedwards
Date: Jul 18 2013
Model: scripts
Version: scripts4_130718
One-line: Update perl5 lib external, multiinstance and thread tests 

M    SVN_EXTERNAL_DIRECTORIES
M    ccsm_utils/Tools/cesm_setup
M    ccsm_utils/Testcases/PET_script
	
================================================================================
Originator: mlevy
Date: Jul 15 2013
Model: scripts
Version: scripts4_130715
One-line:  add support for C_MPAS compset using 120km ocean grid

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/config_definition.xml
================================================================================
Originator: jedwards
Date: Jul 09 2013
Model: scripts
Version: scripts4_130709a
One-line: create_test now appends testmods options to test name, query_tests 
now appends testmods to testlists when -outputlist option invoked.  

M       create_test
M       query_tests

================================================================================
Originator: jedwards
Date: Jul 09 2013
Model: scripts
Version: scripts4_130709
One-line: fix variable expansion in SetupTools.pm (wasn't expanding env vars
in perl)

M       ccsm_utils/Tools/SetupTools.pm

================================================================================
Originator: jshollen
Date: Jul 03 2013
Model: scripts
Version: scripts4_130703
One-line:  add SMS.T42_T42.FARM95C4 test to yellowstone intel prebeta

M      ccsm_utils/Testlistxml/testlist.xml
================================================================================
Originator: jshollen
Date: Jul 01 2013
Model: scripts
Version: scripts4_130701
One-line:  Fix create_test typo, add ERS_Ly5.f10_f10.I20TRCRUCLM45BGC to
aux_clm testlist. 

M       create_test
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: mvertens
Date: Jun 27 2013
Model: scripts
Version: scripts4_130627
One-line: Put in capability to have both xmlchange_cmds file and user_nl file in
	testmods_dir (which was nl_dirs)
	
M       create_newcase
M       create_test
M       query_tests
X       ccsm_utils/Tools/perl5lib
X       ccsm_utils/Tools/lnd/clm/PTCLM
M       ccsm_utils/Case.template/config_definition.xml
	
A  +    ccsm_utils/Testlistxml/testmods_dirs
A       ccsm_utils/Testlistxml/testmods_dirs/cam
A       ccsm_utils/Testlistxml/testmods_dirs/cam/cosp
A       ccsm_utils/Testlistxml/testmods_dirs/cam/cosp/xmlchange_cmnds
	
D       ccsm_utils/Testlistxml/nl_dirs
D       ccsm_utils/Testlistxml/nl_dirs/rtm
D       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOnFloodOnEffvelOff
D       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOnFloodOnEffvelOff/user_nl_rtm
D       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOnFloodOnEffvelOn
D       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOnFloodOnEffvelOn/user_nl_rtm
D       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOnIceOff
D       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOnIceOff/user_nl_rtm
D       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOff
D       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOff/user_nl_rtm
D       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOnIceOn
D       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOnIceOn/user_nl_rtm
D       ccsm_utils/Testlistxml/nl_dirs/clm
D       ccsm_utils/Testlistxml/nl_dirs/clm/monthly_noinitial
D       ccsm_utils/Testlistxml/nl_dirs/clm/monthly_noinitial/user_nl_clm
D       ccsm_utils/Testlistxml/nl_dirs/clm/monthly_noinitial/user_nl_cpl
D       ccsm_utils/Testlistxml/nl_dirs/clm/vrtlay
D       ccsm_utils/Testlistxml/nl_dirs/clm/vrtlay/user_nl_clm
D       ccsm_utils/Testlistxml/nl_dirs/clm/crop
D       ccsm_utils/Testlistxml/nl_dirs/clm/crop/user_nl_clm
D       ccsm_utils/Testlistxml/nl_dirs/clm/oldhyd
D       ccsm_utils/Testlistxml/nl_dirs/clm/oldhyd/user_nl_clm
D       ccsm_utils/Testlistxml/nl_dirs/clm/cn_conly
D       ccsm_utils/Testlistxml/nl_dirs/clm/cn_conly/user_nl_clm
D       ccsm_utils/Testlistxml/nl_dirs/clm/default
D       ccsm_utils/Testlistxml/nl_dirs/clm/default/user_nl_clm
D       ccsm_utils/Testlistxml/nl_dirs/clm/voc
D       ccsm_utils/Testlistxml/nl_dirs/clm/voc/user_nl_clm
D       ccsm_utils/Testlistxml/nl_dirs/clm/voc/user_nl_cpl
D       ccsm_utils/Testlistxml/nl_dirs/clm/ch4_set2_ciso
D       ccsm_utils/Testlistxml/nl_dirs/clm/ch4_set2_ciso/user_nl_clm
D       ccsm_utils/Testlistxml/nl_dirs/clm/ch4_set3_pftroot
D       ccsm_utils/Testlistxml/nl_dirs/clm/ch4_set3_pftroot/user_nl_clm
D       ccsm_utils/Testlistxml/nl_dirs/clm/glcMEC
D       ccsm_utils/Testlistxml/nl_dirs/clm/glcMEC/user_nl_clm
D       ccsm_utils/Testlistxml/nl_dirs/clm/glcMEC/user_nl_cpl
D       ccsm_utils/Testlistxml/nl_dirs/clm/monthly
D       ccsm_utils/Testlistxml/nl_dirs/clm/monthly/user_nl_clm
D       ccsm_utils/Testlistxml/nl_dirs/clm/monthly/user_nl_cpl
D       ccsm_utils/Testlistxml/nl_dirs/clm/rootlit
D       ccsm_utils/Testlistxml/nl_dirs/clm/rootlit/user_nl_clm
D       ccsm_utils/Testlistxml/nl_dirs/clm/ciso
D       ccsm_utils/Testlistxml/nl_dirs/clm/ciso/user_nl_clm
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: santos
Date: Jun 19 2013
Model: scripts
Version: scripts4_130619
One-line: Fix getTiming2.pl grep command that needed an extra escape (\).

M       ccsm_utils/Tools/timing/getTiming2.pl
         - Fix bug introduced in scripts4_130612.

================================================================================
Originator: jedwards
Date: Jun 13 2013
Model: scripts
Version: scripts4_130613
One-line: turn off runoff grid and add ocn/ice domain files for ne240.  
	Fix bug in create_clone (ignoring -mach_dir arg) have 
	check_exactrestart.pl pass S compset cases with no lines to compare.

M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Tools/check_exactrestart.pl
M       create_clone
	
================================================================================
Originator: santos
Date: Jun 12 2013
Model: scripts
Version: scripts4_130612
One-line: Explicitly use extended regex in getTiming2.pl grep command

M       ccsm_utils/Tools/timing/getTiming2.pl

================================================================================
Originator: jshollen
Date: Jun 06 2013
Model: scripts
Version: scripts4_130610
One-line: Re-Removed unneeded section from cism namelist documentation

M      doc/modelnl/nldef2html_cism

================================================================================
Originator: jshollen
Date: Jun 06 2013
Model: scripts
Version: scripts4_130606
One-line: add aux_waccm tests

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: Jun 03 2013
Model: scripts
Version: scripts4_130603
One-line:  creat__test was not setting baseline_root properly. 

M       create_test
	   
================================================================================
Originator: aliceb
Date: May 30 2013
Model: scripts
Version: scripts4_130530
One-line:  updating archive_metadata.sh to include user namelists (user_nl_xxx) (Alice)
	   Correct SST dataset for 1850-pday compsets (Jim)
	
M       ccsm_utils/Case.template/config_compsets.xml
M	ccsm_utils/Tools/archive_metadata.sh

================================================================================
Originator: erik
Date: May 29 2013
Model: scripts
Version: scripts4_130529
One-line:  Update PTCLM, remove VIC temp files, set RTM to either 1850/2000 startup
           fix 1x1_urbanc_alpha 

Update PTCLM

M       SVN_EXTERNAL_DIRECTORIES -- Update PTCLM
M       ccsm_utils/Case.template/config_compsets.xml --- Fix
            some settings for 1x1_urbanc_alpha. Add settings for
            1850 or 2000 startup year for RTM
M       ccsm_utils/Case.template/config_definition.xml -- add note
            about where restart files for DOUT_S_SAVE_INT_REST_FILES
M       ccsm_utils/Testlistxml/testlist.xml -- Remove user_nl_clm
            from VIC tests

------------ Delete the VIC temporary datasets
D       ccsm_utils/Testlistxml/nl_dirs/clm/vic
D       ccsm_utils/Testlistxml/nl_dirs/clm/vic/vic_f09
D       ccsm_utils/Testlistxml/nl_dirs/clm/vic/vic_f09/user_nl_clm
D       ccsm_utils/Testlistxml/nl_dirs/clm/vic/vic_f19
D       ccsm_utils/Testlistxml/nl_dirs/clm/vic/vic_f19/user_nl_clm
D       ccsm_utils/Testlistxml/nl_dirs/clm/vic/vic_vrtlay
D       ccsm_utils/Testlistxml/nl_dirs/clm/vic/vic_vrtlay/user_nl_clm

================================================================================
Originator: santos
Date: May 28 2013
Model: scripts
Version: scripts4_130528a
One-line: Short term archiver no longer deletes data if "mv" errors.

M       ccsm_utils/Tools/st_archive.sh
         - If data couldn't be moved out of a temporary directory, this
           script would end up deleting it all. Now it doesn't do that.

================================================================================
Originator: jshollen
Date: May 28 2013
Model: scripts
Version: scripts4_130528
One-line: query-tests now outputs xml test lists. 

M       query_tests

================================================================================
Originator: sacks
Date: May 24 2013
Model: scripts
Version: scripts4_130524
One-line: update documentation: note cmake requirement for cism

M       doc/usersguide/introduction.xml

================================================================================
Originator: jshollen
Date: May 23 2013
Model: scripts
Version: scripts4_130523
One-line: create_test documentation changes requested by Bill and Tony.

M      create_test

================================================================================
Originator: sacks
Date: May 22 2013
Model: scripts
Version: scripts4_130522a
One-line: Fix some I compsets

(1) Fix I_1850_CLM45_CN_GLC_CISM1 (IG1850CLM45CN) so that it is truly
    using CLM4.5 (changes answers for this compset)

(2) Fix sname of IG4804: change I_1948-2004_CLM45_GLC_CISM1 to
    I_1948-2004_GLC_CISM1

(3) Fix sname: change I_RCP2.6_CLML45_CN to I_RCP2.6_CLM45_CN (note
    extra "L" before the 45 in the old name)

M       ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: jshollen
Date: May 22 2013
Model: scripts
Version: scripts4_130522
One-line: Testreporter update for namelist test summaries. 

M       ccsm_utils/Tools/testreporter.pl

================================================================================
Originator: erik
Date: May 20 2013
Model: scripts
Version: scripts4_130520
One-line: I CRU compset years change

Change 1850 CRUNCEP years to 1901-1920
       2000 CRUNCEP years to 1991-2010
       20TR CRUNCEP years to 1901-1920 with align year of 1901.
       RCP  CRUNCEP years to 1991-2010 with align year of 2005.

Change REFDATE for I1850CLM40CRUCN_f09_g16_clm4500_c130514 to 1122

M       ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: erik
Date: May 17 2013
Model: scripts
Version: scripts4_130517
One-line:  Work on SSP test and add I1850CLM40CRUCN_f09_g16_clm4500_c130514 ref-case for I1850CLM40CRUCN

M       ccsm_utils/Tools/timing/getTiming2.pl -------- fix divide by zero issue
M       ccsm_utils/Case.template/config_compsets.xml - add new compsets for BGCDV and BGCCROP
               add I1850CLM40CRUCN_f09_g16_clm4500_c130514 ref-case for I1850CLM40CRUCN
M       ccsm_utils/Testlistxml/testlist.xml ---------- add tests for new compsets
M       ccsm_utils/Testcases/config_tests.xml -------- Work on SSP test
M       ccsm_utils/Testcases/SSP_script -------------- Get SSP test working

================================================================================
Originator: jedwards
Date: May 15 2013
Model: scripts
Version: scripts4_130515
One-line:  Improve non-standard namelist handeling in create_test, remove a test with no grdi defined.

	M create_test
	M ccsm_utils/Testlistxml/testlist.xml

	
================================================================================
Originator: jshollen
Date: May 13 2013
Model: scripts
Version: scripts4_130513b
One-line:  change f09_g16 TG to ERS, add shortlist.glc.auxtests as
aux_glc_short. 

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jshollen
Date: May 13 2013
Model: scripts
Version: scripts4_130513a
One-line:  compare_namelist fix. 

M       ccsm_utils/Tools/compare_namelists.pl


================================================================================
Originator: jshollen
Date: May 13 2013
Model: scripts
Version: scripts4_130513
One-line:  testreporter xml bug fix, compare_namelist new output allows for
eaiser summarization of output. 

M       ccsm_utils/Tools/compare_namelists.pl
M       ccsm_utils/Tools/testreporter.pl

================================================================================
Originator: jshollen
Date: May 11 2013
Model: scripts
Version: scripts4_130511
One-line:  Create_test sets the baseline tag in testspec.xml, testreporter has
more robust xml reporting. 

M       create_test
M       ccsm_utils/Tools/testreporter.pl

================================================================================
Originator: erik
Date: May 10 2013
Model: scripts
Version: scripts4_130510
One-line:  Fix more create_test issues, add more tests for CLM45/DV/CROP to test lists

M       create_test --- Fix cesmroot option and problem that would abort on setting
           STOP_N when it was undef.
M       ccsm_utils/Case.template/config_compsets.xml -- add ICRUCLM45BGCTEST for
           testing CLM45 with CRUNCEP forcing on machines with only 2002-2003 data
MM      ccsm_utils/Testlistxml/sort_xml_entries -- fix typo so would work
M       ccsm_utils/Testlistxml/testlist.xml --- Add CLM45/DV/CROP tests
M       ccsm_utils/Testcases/SSP_script -- Change final mode to ref case hybrid
M       ccsm_utils/Testcases/SBN_script -- Fix bug 1687

================================================================================
Originator: jshollen
Date: May 9 2013
Model: scripts
Version: scripts4_130509a
One-line:  create_test now uses config_machines baselineroot if none is set, better
error checking for namelist comparisons.  

M       create_test

================================================================================
Originator: erik
Date: May 9 2013
Model: scripts
Version: scripts4_130509
One-line:  Fix some bugs with create_test, fix SBN test, add SSP, fix S1850 compsets
           remove asphaltjungle grid

   Fix -help, and -mach machine_compiler in create_test, make sure baselineroot set for compare/generate
   Fix the S1850 compsets
   Fix SBN test add a SSP test for CLM45 spinup
   Add and improve documentation of I compsets.
   Work with testlist for I compsets a bit.
   Remove asphaltjungle grid

M       create_test
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Testlistxml/testlist.xml
M       ccsm_utils/Testcases/config_tests.xml

A       ccsm_utils/Case.template/config_compsets.xsl - So can view the XML compset file
            to check for errors.
A       ccsm_utils/Testcases/SSP_script ---- CLM45 spinup test

================================================================================
Originator: mvertens
Date: May 7 2013
Model: scripts
Version: scripts4_130508

One-line: more updates to config_compsets.xml and testlist.xml that fix SFAILS in
	  namelists tests

M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Testlistxml/testlist.xml
M       README
	
================================================================================
Originator: mvertens
Date: May 7 2013
Model: scripts
Version: scripts4_130507b
One-line: more updates to config_compsets.xml and testlist.xml that fix SFAILS in
	  namelists tests
	
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Testlistxml/testlist.xml
	
================================================================================
Originator: mvertens
Date: May 7 2013
Model: scripts
Version: scripts4_130507a
One-line: updates to config_compsets.xml and testlist.xml that fix SFAILS in
	  namelists tests

M       ccsm_utils/Case.template/config_compsets.xml
D       ccsm_utils/Testlistxml/sort_file.pl
A  +    ccsm_utils/Testlistxml/sort_xml_entries
M       ccsm_utils/Testlistxml/add_xml_entries
M       ccsm_utils/Testlistxml/testlist.xml
	
================================================================================
Originator: jshollen
Date: May 7 2013
Model: scripts
Version: scripts4_130507
One-line: testlist update, bug fixes for compare_namelists and create_test,
          query_tests docs updates. 

M       create_test
M       ccsm_utils/Tools/compare_namelists.pl
M       ccsm_utils/Testlistxml/testlist.xml
M       query_tests

================================================================================
Originator: mvertens
Date: May 6 2013
Model: update for
Version: scripts4_130506
One-line: Updates to create_tables for new compsets and grids

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
A       doc/modelnl/nl_clm45.html
M       doc/modelnl/nl_drv.html
M       doc/modelnl/nl_cism.html
D       doc/modelnl/nldef2html_clm
M       doc/modelnl/machines.html
M       doc/modelnl/xmldef2html_grid
M       doc/modelnl/compsets.html
M       doc/modelnl/nldef2html_cism
A       doc/modelnl/nl_clm40.html
M       doc/modelnl/env_run.html
M       doc/modelnl/env_case.html
M       doc/modelnl/grid.html
M       doc/modelnl/nl_cam.html
M       doc/modelnl/create_tables
M       doc/modelnl/nl_rtm.html
A       doc/modelnl/nldef2html_clm40
M       doc/modelnl/env_pesetup.html
M       doc/modelnl/xmldef2html_compsets
A       doc/modelnl/nldef2html_clm45
M       doc/modelnl/nl_clm.html
M       doc/modelnl/env_build.html
M       doc/modelnl/nl_pop2.html
M       create_newcase
	
================================================================================
Originator: tcraig
Date: May 3 2013
Model: scripts
Version: scripts4_130503a
One-line: modify rx1 grid/compset match
M       ccsm_utils/Case.template/config_grid.xml
================================================================================

Originator: mlevy
Date: May 3 2013
Model: scripts
Version: scripts4_130503
One-line: set B[1850,20TR][BPRP,BDRD] compsets to startup rather than hybid
          (new POP ecosystem means no spun-up cases to start from)

M       ccsm_utils/Case.template/config_compsets.xml

================================================================================

Originator: jshollen
Date: May 2 2013
Model: scripts
Version: scripts4_130502
One-line: create_test: fix bugs found in 130501e

M       create_test

================================================================================

Originator: jshollen
Date: May 1 2013
Model: scripts
Version: scripts4_130501e
One-line: create_test: update documentation, fix bugs

M       create_test

================================================================================

Originator: jedwards
Date: May 1 2013
Model: scripts
Version: scripts4_130501d
One-line: fix preview_namelist script again

M ccsm_utils/Tools/preview_namelists

===================================================================================

Originator: jedwards
Date: May 1 2013
Model: scripts
Version: scripts4_130501c
One-line: fix preview_namelist script for multiinstance cases

M ccsm_utils/Tools/preview_namelists


================================================================================

Originator: erik
Date: May 1 2013
Model: scripts
Version: scripts4_130501b
One-line: Move aux_clm_int tests to aux_clm for NAG, remove many tests that fail

Remove the anoxia_wt tests and remove setting of shape_fluxprof as it was removed.
Remove all but one PET_PT test for CLM45 as it's known NOT to work with threading.
Change aux_clm_int tests to aux_clm for NAG compiler.
Change ERS_Ln48 tests to either SMS or ERS_Ld3.

M       ccsm_utils/Testlistxml/nl_dirs/clm/ch4_set2_ciso/user_nl_clm
M       ccsm_utils/Testlistxml/nl_dirs/clm/ch4_set3_pftroot/user_nl_clm
D       ccsm_utils/Testlistxml/nl_dirs/clm/anoxia_wtsat
D       ccsm_utils/Testlistxml/nl_dirs/clm/anoxia_wtsat/user_nl_clm
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================

Originator: jedwards
Date: May 1 2013
Model: scripts
Version: scripts4_130501
One-line: add performance tests to testlist.xml, add outputlist option to query_tests, update add_xml_entries to sort and format output

M      ccsm_utils/Testlistxml/add_xml_entries
M      ccsm_utils/Testlistxml/testlist.xml
M      query_tests

================================================================================

Originator: jedwards
Date: Apr 30 2013
Model: scripts
Version: scripts4_130430b
One-line: add build and tests for cam aquaplanet cases

M           ccsm_utils/Case.template/config_compsets.xml
M           ccsm_utils/Case.template/config_grid.xml
M           ccsm_utils/Case.template/config_definition.xml
M           ccsm_utils/Testlistxml/testlist.xml



================================================================================

Originator: mvertens
Date: Apr 30 2013
Model: scripts
Version: scripts4_130430a
One-line: updated query for testlists

D       query_testlist
A  +    query_tests

================================================================================
	
Originator: mvertens
Date: Apr 30 2013
Model: scripts
Version: scripts4_130430
One-line: changes for cam ideal physics and adiabatic

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Testlistxml/testlist.xml
M       create_newcase
	
================================================================================
	
Originator: jshollen
Date: Apr 26 2013
Model: scripts
Version: scripts4_130426b
One-line: create_test compare namelist compare bug fix.  

M       create_test

================================================================================

Originator: santos
Date: Apr 26 2013
Model: scripts
Version: scripts4_130426a
One-line: Change error message mentioning "cpl.log" to one mentioning "cesm.log"

M       ccsm_utils/Tools/cesm_postrun_setup
         - Check cpl.log to see if a case ran, but mention cesm.log in error
           messages.

================================================================================
Originator: jedwards
Date: Apr 26 2013
Model: scripts
Version: scripts4_130426
One-line: remove IOP.T31 tests, fix some compset parsing 

M              ccsm_utils/Case.template/config_compsets.xml
M              ccsm_utils/Case.template/ConfigCase.pm
M              ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: jedwards
Date: Apr 25 2013
Model: scripts
Version: scripts4_130425
One-line: fix scam use_case in config_compsets.xml

M      config_compsets.xml


================================================================================
Originator: jshollen
Date: Apr 24 2013
Model: scripts
Version: scripts4_130424b
One-line: Jim's compset fixes, update SCAM iop dataset. 

M       SVN_EXTERNAL_DIRECTORIES
M       ccsm_utils/Case.template/config_compsets.xml
M       create_newcase

================================================================================
Originator: sacks
Date: Apr 24 2013
Model: scripts
Version: scripts4_130424a
One-line: Removed unneeded section from cism namelist documentation

M       doc/modelnl/nldef2html_cism
	
================================================================================
Originator: jedwards
Date: Apr 24 2013
Model: scripts
Version: scripts4_130424
One-line: Fixed SCAM (FARM95) test inputfile, check compset longnames against def
	     remove some unsupported tests, fix a bad compset name
	
M               ccsm_utils/Tools/cesm_postrun_setup
M               ccsm_utils/Case.template/config_compsets.xml
M               ccsm_utils/Testlistxml/testlist.xml
M               create_newcase


================================================================================
Originator: mvertens
Date: Apr 23 2013
Model: scripts
Version: scripts4_130423a
One-line: Fixed SCAM (FARM95) test and added it to pre-alpha test for yellowstone

M       ccsm_utils/Tools/cesm_setup
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/ConfigCase.pm
M       ccsm_utils/Case.template/config_definition.xml
	  - above are all fixes for SCAM to work
	
M       ccsm_utils/Testlistxml/testlist.xml
  	  - new prelpha yellowstone test for FARM95 
	
M       create_newcase
	  - help comments just slightly changed
	
================================================================================
Originator: jshollen
Date: Apr 23 2013
Model: scripts
Version: scripts4_130423
One-line: Update test list, namelist compare log file for test suite runs, 
create_test bug fix. 

M       create_test
M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: mvertens
Date: Apr 23 2013
Model: scripts
Version: scripts4_130423
One-line: Fixed SCAM (FARM95) test and added it to pre-alpha test for yellowstone

M       ccsm_utils/Tools/cesm_setup
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/ConfigCase.pm
M       ccsm_utils/Case.template/config_definition.xml
	  - above are all fixes for SCAM to work
	
M       ccsm_utils/Testlistxml/testlist.xml
  	  - new prelpha yellowstone test for FARM95 
	
M       create_newcase
	  - help comments just slightly changed
	
================================================================================
Originator: mvertens
Date: Apr 22 2013
Model: scripts
Version: scripts4_130422b
One-line: Update test list 

M       query_testlist
M       ccsm_utils/Testlistxml/add_xml_entries
M       ccsm_utils/Testlistxml/make_xml_entries
	   - added more functionality to add_xml_entries, make_xml_entries and 
   	     query_testlist  - call each of these with the help flag for more 
 	     documentation
	
M       ccsm_utils/Case.template/config_compsets.xml
	   - added new case for BGCDV 
	
M       ccsm_utils/Case.template/config_grid.xml
	   - fixed bug for T21_T21
	
M       ccsm_utils/Testlistxml/nl_dirs/clm/monthly/user_nl_cpl
A       ccsm_utils/Testlistxml/nl_dirs/clm/monthly_noinitial
A       ccsm_utils/Testlistxml/nl_dirs/clm/monthly_noinitial/user_nl_clm
A       ccsm_utils/Testlistxml/nl_dirs/clm/monthly_noinitial/user_nl_cpl
M       ccsm_utils/Testlistxml/testlist.xml
           -  added the following tests
        ERS_Ly5.f10_f10.I20TRCLM45CN4Me.yellowstone_intel.clm/monthly_noinitial	
        SMS_D.f09_g16.TG.yellowstone_intel
        SMS_D_Ly1.f09_g16_gl10.TGIS2.yellowstone_intel
        ERI.f19_f19.FC5PM.yellowstone_intel
        ERS.f19_g16.BNUKE_C4WBC_L40CN.yellowstone_intel
        ERS.f19_f19.FGEOS_C4WSF_L40CN.yellowstone_intel
        ERI.f19_f19.F1850C5PM.yellowstone_intel
        ERI.f19_f19.FSSOA.yellowstone_intel
        ERI.f19_f19.FC5CLUBB.yellowstone_intel


	
================================================================================
Originator: santos
Date: Apr 22 2013
Model: scripts
Version: scripts4_130422a
One-line: Fix SCAM compset, add references to cesm.log file

M       ccsm_utils/Tools/st_archive.sh
         - Change ccsm.log to cesm.log

M       ccsm_utils/Tools/cesm_postrun_setup
         - Users are referred to cesm.log rather than cpl.log in case of error.

M       ccsm_utils/Case.template/config_compsets.xml
         - Change F_ARM95_CAM5 to CAM4, since we only have a CAM4 use_case.

================================================================================
Originator: jedwards
Date: Apr 22 2013
Model: scripts
Version: scripts4_130422
One-line: correct update of  testlist.xml including nldir field

M     ccsm_utils/Testlistxml/testlist.xml


================================================================================
Originator: fischer
Date: Apr 19 2013
Model: scripts
Version: scripts4_130419a
One-line: rename homme to se, fix T42_T42 masks

M       SVN_EXTERNAL_DIRECTORIES
  perl5lib needed updating to replace homme with se

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_definition.xml
  rename homme to se

M       ccsm_utils/Case.template/config_grid.xml
  replace gx1v6 mask with usgs

================================================================================
Originator: jshollen
Date: Apr 19 2013
Model: scripts
Version: scripts4_130419
One-line: Testlist updates: Jim's clm updates, remove _IOP for T31_g37
yellowstone, add frankfurt prebeta

M       ccsm_utils/Testlistxml/testlist.xml

================================================================================
Originator: santos
Date: Apr 18 2013
Model: scripts
Version: scripts4_130418
One-line: Fixed problem with GEOS+WCSF compset.

M       ccsm_utils/Case.template/config_compsets.xml
         - Change some "%WCCM" matches back to "%WC". This fixes
           a problem where the compset with GEOS and WCSF had the
           wrong value for "-nlev" in CAM_CONFIG_OPTS.

================================================================================
Originator: mvertens
Date: Apr 17 2013
Model: scripts
Version: scripts4_130417
One-line: Fixed problem with docn changes for AMIP compsets
	This bug was introduced in scripts4_130414b
	
Verified that running create_test -nlcompareonly versus cesm1_2_alpha06e gave
	identical namelists other than for the F compsets with transient docn
	forcing

M       query_testlist
M       ccsm_utils/Tools/compare_namelists.pl
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
	
================================================================================
Originator: fvitt
Date: Apr 16 2013
Model: scripts
Version: scripts4_130416b
One-line: Restored "-chem waccm_mozart" CAM_CONFIG_OPTS for compsets %WCBC and %WCMX

M       ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: jshollen
Date: Apr 16 2013
Model: scripts
Version: scripts4_130416a
One-line: Stop testcase.test scripts from overwriting TestStatus,
TestStatus.out, needed for namelist comparison testing. 

M       ccsm_utils/Tools/testcase_begin
M       ccsm_utils/Tools/testcase_end

================================================================================
Originator: mvertens
Date: Apr 16 2013
Model: scripts
Version: scripts4_130416
One-line: updates to aux_clm tests

M       ccsm_utils/Tools/testcase_end
M       ccsm_utils/Testlistxml/add_xml_entries
M       ccsm_utils/Testlistxml/testlist.xml
	
================================================================================
Originator: jshollen
Date: Apr 15 2013
Model: scripts
Version: scripts4_130415
One-line: merge branches/a06e_scripts4_130415, remove text-based test lists,
test system bug fixes.  

M       ccsm_utils/Tools
M       ccsm_utils/Tools/testreporter.pl
M       ccsm_utils/Tools/taskmaker.pl
M       ccsm_utils/Tools/ccsm_getenv
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Testlistxml/testlist.xml
M       ccsm_utils/Testcases
D       ccsm_utils/Testlists
D       ccsm_utils/Testlists/titan.glc.auxtest
D       ccsm_utils/Testlists/lynx_pgi.prerelease
D       ccsm_utils/Testlists/bluefire.prerelease
D       ccsm_utils/Testlists/frankfurt_pgi.prerelease
D       ccsm_utils/Testlists/bluefire.science.auxtest
D       ccsm_utils/Testlists/hadley.auxtest
D       ccsm_utils/Testlists/bluefire.pop2.auxtest
D       ccsm_utils/Testlists/titan.cam.auxtest
D       ccsm_utils/Testlists/frankfurt_lahey.prerelease
D       ccsm_utils/Testlists/eastwind.rasm.auxtest
D       ccsm_utils/Testlists/garnet.rasm.auxtest
D       ccsm_utils/Testlists/yellowstone_intel.glc.auxtest
D       ccsm_utils/Testlists/lynx_intel.prealpha
D       ccsm_utils/Testlists/B01
D       ccsm_utils/Testlists/yellowstone_intel.clm.auxtest
D       ccsm_utils/Testlists/B41
D       ccsm_utils/Testlists/B42
D       ccsm_utils/Testlists/frankfurt_intel.prealpha
D       ccsm_utils/Testlists/B43
D       ccsm_utils/Testlists/B44
D       ccsm_utils/Testlists/janus_intel.prebeta
D       ccsm_utils/Testlists/hopper_pgi.glc.auxtest
D       ccsm_utils/Testlists/newmachine.port.auxtest
D       ccsm_utils/Testlists/lynx_gnu.prealpha
D       ccsm_utils/Testlists/titan.prebeta
D       ccsm_utils/Testlists/frankfurt_intel.clm.auxtest
D       ccsm_utils/Testlists/intrepid.prerelease
D       ccsm_utils/Testlists/olympus.rasm.auxtest
D       ccsm_utils/Testlists/janus_intel.prealpha
D       ccsm_utils/Testlists/lynx_pathscale.prealpha
D       ccsm_utils/Testlists/brutus_pm.auxtest
D       ccsm_utils/Testlists/brutus_po.auxtest
D       ccsm_utils/Testlists/hopper.prebeta
D       ccsm_utils/Testlists/evergreen.rasm.auxtest
D       ccsm_utils/Testlists/bluefire.esmf.auxtest
D       ccsm_utils/Testlists/yellowstone_pgi.glc.auxtest
D       ccsm_utils/Testlists/hopper.prerelease
D       ccsm_utils/Testlists/lynx_pgi.prealpha
D       ccsm_utils/Testlists/bluefire.prealpha
D       ccsm_utils/Testlists/yellowstone_pgi.clm.auxtest
D       ccsm_utils/Testlists/frankfurt_pgi.prealpha
D       ccsm_utils/Testlists/intrepid.prebeta
D       ccsm_utils/Testlists/yellowstone_intel.prebeta
D       ccsm_utils/Testlists/brutus_im.auxtest
D       ccsm_utils/Testlists/lynx_intel.prebeta
D       ccsm_utils/Testlists/hopper_gnu.glc.auxtest
D       ccsm_utils/Testlists/brutus_io.auxtest
D       ccsm_utils/Testlists/frankfurt_pgi.clm.auxtest
D       ccsm_utils/Testlists/yellowstone_intel.prealpha
D       ccsm_utils/Testlists/bluefire.cice1.auxtest
D       ccsm_utils/Testlists/bluefire.cice2.auxtest
D       ccsm_utils/Testlists/allIcompsetsRes.clm.auxtest
D       ccsm_utils/Testlists/A01
D       ccsm_utils/Testlists/A02
D       ccsm_utils/Testlists/bluefire.drv.auxtest
D       ccsm_utils/Testlists/A03
D       ccsm_utils/Testlists/C01
D       ccsm_utils/Testlists/frankfurt_intel.prerelease
D       ccsm_utils/Testlists/A04
D       ccsm_utils/Testlists/C02
D       ccsm_utils/Testlists/E01
D       ccsm_utils/Testlists/A05
D       ccsm_utils/Testlists/C03
D       ccsm_utils/Testlists/C04
D       ccsm_utils/Testlists/C41
D       ccsm_utils/Testlists/C05
D       ccsm_utils/Testlists/C42
D       ccsm_utils/Testlists/bluefire.cam.auxtest
D       ccsm_utils/Testlists/C43
D       ccsm_utils/Testlists/bluefire.rasm.auxtest
D       ccsm_utils/Testlists/C44
D       ccsm_utils/Testlists/C45
D       ccsm_utils/Testlists/shortlist.glc.auxtest
D       ccsm_utils/Testlists/C46
D       ccsm_utils/Testlists/C47
D       ccsm_utils/Testlists/chugach.rasm.auxtest
D       ccsm_utils/Testlists/namelists.frankfurt
D       ccsm_utils/Testlists/shortlist.clm.auxtest
D       ccsm_utils/Testlists/S01
D       ccsm_utils/Testlists/raptor.rasm.auxtest
D       ccsm_utils/Testlists/yellowstone_pgi.prebeta
D       ccsm_utils/Testlists/frankfurt_lahey.prebeta
D       ccsm_utils/Testlists/titan.prerelease
D       ccsm_utils/Testlists/W01
D       ccsm_utils/Testlists/lynx_pgi.prebeta
D       ccsm_utils/Testlists/bluefire.prebeta
D       ccsm_utils/Testlists/yellowstone_pgi.prealpha
M       create_newcase

================================================================================
Originator: mvertens
Date: Apr 14 2013
Model: scripts
Version: scripts4_130414b
One-line: updated config_compsets.xml to have sname correspond to older longname
	for backwards compatibility, added query_testlist functionality and
	added capability to have testlists point to nl_dirs for the clm_aux tests
        migrated aux_clm2 and aux_rtm2 into aux_clm tests
	 
A       query_testlist
M       create_test
M       ccsm_utils/Case.template/config_grid.xml
A  +    ccsm_utils/Testlistxml/add_xml_entries
D       ccsm_utils/Testlistxml/add_entry.pl
D       ccsm_utils/Testlistxml/parse_file.pl
A       ccsm_utils/Testlistxml/make_xml_entries
M       ccsm_utils/Testlistxml/testlist.xml
A       ccsm_utils/Testlistxml/nl_dirs
A       ccsm_utils/Testlistxml/nl_dirs/rtm
A       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOnFloodOnEffvelOff
A       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOnFloodOnEffvelOff/user_nl_rtm
A       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOnFloodOnEffvelOn
A       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOnFloodOnEffvelOn/user_nl_rtm
A       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOnIceOff
A       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOnIceOff/user_nl_rtm
A       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOff
A       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOff/user_nl_rtm
A       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOnIceOn
A       ccsm_utils/Testlistxml/nl_dirs/rtm/rtmOnIceOn/user_nl_rtm
A       ccsm_utils/Testlistxml/nl_dirs/clm
A       ccsm_utils/Testlistxml/nl_dirs/clm/anoxia_wtsat
A       ccsm_utils/Testlistxml/nl_dirs/clm/anoxia_wtsat/user_nl_clm
A       ccsm_utils/Testlistxml/nl_dirs/clm/vic
A       ccsm_utils/Testlistxml/nl_dirs/clm/vic/vic_f09
A       ccsm_utils/Testlistxml/nl_dirs/clm/vic/vic_f09/user_nl_clm
A       ccsm_utils/Testlistxml/nl_dirs/clm/vic/vic_f19
A       ccsm_utils/Testlistxml/nl_dirs/clm/vic/vic_f19/user_nl_clm
A       ccsm_utils/Testlistxml/nl_dirs/clm/vic/vic_vrtlay
A       ccsm_utils/Testlistxml/nl_dirs/clm/vic/vic_vrtlay/user_nl_clm
A       ccsm_utils/Testlistxml/nl_dirs/clm/vrtlay
A       ccsm_utils/Testlistxml/nl_dirs/clm/vrtlay/user_nl_clm
A       ccsm_utils/Testlistxml/nl_dirs/clm/crop
A       ccsm_utils/Testlistxml/nl_dirs/clm/crop/user_nl_clm
A       ccsm_utils/Testlistxml/nl_dirs/clm/oldhyd
A       ccsm_utils/Testlistxml/nl_dirs/clm/oldhyd/user_nl_clm
A       ccsm_utils/Testlistxml/nl_dirs/clm/cn_conly
A       ccsm_utils/Testlistxml/nl_dirs/clm/cn_conly/user_nl_clm
A       ccsm_utils/Testlistxml/nl_dirs/clm/default
A       ccsm_utils/Testlistxml/nl_dirs/clm/default/user_nl_clm
A       ccsm_utils/Testlistxml/nl_dirs/clm/voc
A       ccsm_utils/Testlistxml/nl_dirs/clm/voc/user_nl_clm
A       ccsm_utils/Testlistxml/nl_dirs/clm/voc/user_nl_cpl
A       ccsm_utils/Testlistxml/nl_dirs/clm/ch4_set2_ciso
A       ccsm_utils/Testlistxml/nl_dirs/clm/ch4_set2_ciso/user_nl_clm
A       ccsm_utils/Testlistxml/nl_dirs/clm/ch4_set3_pftroot
A       ccsm_utils/Testlistxml/nl_dirs/clm/ch4_set3_pftroot/user_nl_clm
A       ccsm_utils/Testlistxml/nl_dirs/clm/glcMEC
A       ccsm_utils/Testlistxml/nl_dirs/clm/glcMEC/user_nl_clm
A       ccsm_utils/Testlistxml/nl_dirs/clm/glcMEC/user_nl_cpl
A       ccsm_utils/Testlistxml/nl_dirs/clm/monthly
A       ccsm_utils/Testlistxml/nl_dirs/clm/monthly/user_nl_clm
A       ccsm_utils/Testlistxml/nl_dirs/clm/monthly/user_nl_cpl
A       ccsm_utils/Testlistxml/nl_dirs/clm/rootlit
A       ccsm_utils/Testlistxml/nl_dirs/clm/rootlit/user_nl_clm
A       ccsm_utils/Testlistxml/nl_dirs/clm/ciso
A       ccsm_utils/Testlistxml/nl_dirs/clm/ciso/user_nl_clm

	
================================================================================
Originator: jshollen
Date: Apr 12 2013
Model: scripts
Version: scripts4_130412a
One-line: merge nlcomparefixes_scripts4_130410c_tags/nlcomparefixes00_scripts4_130410c
getting nlcompareonly and create_test namelist comparisons working. 

M       create_test
M       ccsm_utils/Tools
M       ccsm_utils/Tools/testcase_begin
M       ccsm_utils/Tools/cs.status
M       ccsm_utils/Tools/taskmaker.pl
M       ccsm_utils/Tools/ccsm_getenv
M       ccsm_utils/Tools/testcase_end
M       ccsm_utils/Tools/cs.submit
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/config_definition.xml
A  +    ccsm_utils/Testlistxml/sort_file.pl
M       ccsm_utils/Testlistxml/testlist.xml
M       ccsm_utils/Testcases
M       ccsm_utils/Testlists/A01

================================================================================
Originator: sacks
Date: Apr 12 2013
Model: scripts
Version: scripts4_130412
One-line: Fix CISM_USE_TRILINOS for CISM2S compsets

M       ccsm_utils/Case.template/config_compsets.xml
	
================================================================================
Originator: jedwards
Date: Apr 10 2013
Model: scripts
Version: scripts4_130410c
One-line: Update testlist.xml, fix parsing tools

M            ccsm_utils/Testlistxml/add_entry.pl
M            ccsm_utils/Testlistxml/parse_file.pl
M            ccsm_utils/Testlistxml/testlist.xml
	


================================================================================
Originator: fischer
Date: Apr 10 2013
Model: scripts
Version: scripts4_130410b
One-line: Fix T42_T42 grid

M       ccsm_utils/Case.template/config_grid.xml
   Typo for T42_T42 grid USGS should be usgs

================================================================================
Originator: jedwards, mvertens, jshollen
Date: Apr 10 2013
Model: scripts
Version: scripts4_130410a
One-line: Add support for namelist consistancy checking in create_test

M            ccsm_utils/Tools/testcase_end
M            ccsm_utils/Tools/cs.submit
M            ccsm_utils/Case.template/config_grid.xml
M            ccsm_utils/Testlistxml/testlist.xml
M            create_test
M            create_newcase


================================================================================
Originator: tcraig
Date: Apr 10 2013
Model: scripts
Version: scripts4_130410
One-line: Add OCP test, remove CSMDATA from check_input_data
	
M       ccsm_utils/Tools/check_input_data
M       ccsm_utils/Testcases/OCP_script
M       ccsm_utils/Testcases/OCP_build.csh
M       ccsm_utils/Testcases/config_tests.xml
M       ccsm_utils/Testcases/ICP_build.csh
================================================================================
Originator: mlevy
Date: Apr 9 2013
Model: scripts
Version: scripts4_130409
One-line: Rolling back T31_g37 mapping files to pre-Apr 4 tag. Also, updating
          config_definition to note that "cart3d" is now default for vect_map.

M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/config_definition.xml

================================================================================

Originator: mvertens
Date: Apr 8 2013
Model: scripts
Version: scripts4_130408b
One-line: updated scripts trunk to branch_tags/newcompsets2_tags/newcompsets2_01_scripts4_130405a
          this is the branch tag that will be used to create clm4_0_71 trunk tag
	
M       create_newcase
M       ccsm_utils/Case.template/config_compsets.xml

	Ran namelists.yellowstone tests to verify that all namelist changes were expected
	
================================================================================
Originator: jshollen
Date: Apr 8 2013
Model: scripts
Version: scripts4_130408a
One-line: merge in branch_tags/xmltests_tags/xmltests03_scripts4_130408, new
create_test, xml test lists, added/deleted tests for Chris.

M       create_test
M      ccsm_utils/Tools
M       ccsm_utils/Tools/testcase_begin
A  +    ccsm_utils/Tools/cs.status
M      ccsm_utils/Tools/taskmaker.pl
M      ccsm_utils/Tools/ccsm_getenv
A  +    ccsm_utils/Tools/cs.submit
M      ccsm_utils/Case.template/config_definition.xml
A  +    ccsm_utils/Testlistxml
A  +    ccsm_utils/Testlistxml/add_entry.pl
A  +    ccsm_utils/Testlistxml/parse_file.pl
A  +    ccsm_utils/Testlistxml/testlist.xml
M      ccsm_utils/Testcases
M       ccsm_utils/Testlists/B43
M       ccsm_utils/Testlists/C04
M       ccsm_utils/Testlists/C05
M       ChangeLog
M       create_newcase
D       create_test_suite


================================================================================
Originator: santos
Date: Apr 8 2013
Model: scripts
Version: scripts4_130408
One-line: Add CARMA compsets, change "WACM" to "WCCM" in compset names.

M       ccsm_utils/Case.template/config_compsets.xml
         - Change "WACM" to "WCCM". This means that all WACCM compsets
           now match "CAM[45]%WC", and this can be used in matching rules.
         - Add compset "BNUKE_C4WBC_L40CN", a hypothetical nuclear winter
           scenario that uses WACCM and the stratospheric black carbon CARMA
           model.
         - Add compset "FGEOS_C4WSF_L40CN", a specified-dynamics case with
           WACCM sulfur chemistry and the sulfate CARMA model.
         - Used matching rules and "additive" nature of CAM_CONFIG_OPTS to
           simplify and generalize some entries.
         - Fix CCSM_CO2_PPMV setting for WACCM-X cases.
         - Remove duplicate setting of RUN_STARTDATE for RCP._CAM compsets.

M       ccsm_utils/Case.template/config_grid.xml
         - Add grids needed by CAM tests that were removed in a previous tag.
            - These include horizontal resolutions T05, T21, and T42.

M       ccsm_utils/Testlists/B43
         - Add tests for new compsets.

================================================================================
Originator: fischer
Date: 5 Apr 2013
Model: scripts
Version: scripts4_130405a
One-line: add prescribed mam and clubb compsets, CCSM_CO2_PPMV settings

M       ccsm_utils/Case.template/config_compsets.xml
 - add new 1850 and 2000 prescribed mam compsets
 - add new CLUBB compset
 - add new F_2000_STRATSOA (FSTRATSOA) compset
 - remove F_TROP_STRAT_CHEM_CN
 - change reference to "waccm_mozart_v1" to "waccm_mozart"
 - set CCSM_CO2_PPMV to 0.000001 for all *_WACCM, *_SD, *_MOZART, _MOZSOA, *_STRAT* compsets
 - turn off "CN" for compset FMOZSOA
 - configure FMOZSOA compset with cam4 phys
 - add age-of-air tracers to *_MOZMAM, *_MOZSOA, *_STRATMAM* compsets
 - set CICE_NAMELIST_OPTS="cam5=.true." for  *_STRATMAM*, *_MOZMAM compsets
 - the hybrid startup type set up for B_1850_WACCM5_CN needs a 70-level CAM IC file
 - new build-namelist use_case files for compsets *_MOZMAM, *_MOZSOA, *_STRATMAM*

M       ccsm_utils/Testlists/B43
 - Add ERI.f19_f19.FC5PM for prescribed MAM test

================================================================================
Originator: sacks
Date: 05 April 2013
Model: scripts
Version: scripts4_130405
One-line: Use new, working initial conditions for BG1850CN @ f09

Created a new refcase directory, bg40.1850.track1.1deg.006b, which is
similar to b40.1850.track1.1deg.006, but with spun-up CISM initial
conditions and CLM initial conditions that agree with the newer CLM surface
datasets for IG runs.

*** Use new refcase
M       ccsm_utils/Case.template/config_compsets.xml

*** Add a BG1850CN f09 test
M       ccsm_utils/Testlists/E01

*** Remove documentation that BG1850CN is broken at f09
M       doc/modelnl/compsets.html


================================================================================
Originator: mnlevy
Date: 04 April 2013
Model: scripts
Version: scripts4_130404
One-line: Changed mapping files for f09_g16, f19_g16, T62_g16, T62_g37, and
          T31_g37. All use aave for ocn2atm[F/S] and atm2ocnF, bilin for
					atm2ocnS, and patch for atm2ocnV. Also set default for VECT_MAP
					to cart3d (though I didn't update any documentation about this...)

M       ccsm_utils/Case.template/config_grid.xml

================================================================================
Originator: sacks
Date: 03 April 2013
Model: scripts
Version: scripts4_130403
One-line: add CISM_OBSERVED_IC in env_run.xml
	
M       ccsm_utils/Case.template/config_definition.xml
	
================================================================================
Originator: mvertens
Date: 30 Mar 2013
Model: scripts
Version: scripts4_130330
One-line: turned on RTM by default for F compsets (RTM_MODE is set to ACTIVE) and
	added user_compset to create_newcase options
	

M       ccsm_utils/Case.template/config_compsets.xml
M       create_newcase
	
	Put in check in create_newcase so that if there is no LND2ROF_FMAPNAME
	file (i.e. LND2ROF_FMAPNAME = idmap) then RTM_MODE will be set to NULL 
	
        Also - added a new argument user_compset so that the user can specify
	a new compset longname on the command line  
	
	
================================================================================
Originator: tcraig
Date: 25 Mar 2013
Model: scripts
Version: scripts4_130326
One-line: add wave model, CISM_GRID env variable, other mods

- add wave models (swav,xwav,ww3) in scripts, compsets, timer output, tests, env variables
- merge $SVNREPO/scripts/trunk_tags/scripts4_130222 $SVNREPO/scripts/branch_tags/ww3bq_tags/ww3bq01_scripts4_130222
- add B1850CLM45CN, B1850CLM45CNF, and I1850CLM45CNF compsets
- add some clm45 and flood tests
- add back GLC_GRID, GLC_NX, GLC_NY.  add new env variable CISM_GRID
- cleanup testcase_end indentation and fix missing end
- add AWAV compset and WW3 datm mode, ww3a grid
- add DICE%COPY and DOCN%COPY and some data model NULL modes
- add RTM%FLOOD option
- add XROF_FLOOD_MODE env variable
- initial migration of MAPTYPE back to scripts, but commented out, still in driver

M       create_test
M       ccsm_utils/Tools/testcase_begin
M       ccsm_utils/Tools/lt_archive.sh
M       ccsm_utils/Tools/st_archive.sh
M       ccsm_utils/Tools/cesm_setup
M       ccsm_utils/Tools/timing/getTiming2.pl
M       ccsm_utils/Tools/taskmaker.pl
M       ccsm_utils/Tools/ccsm_getenv
M       ccsm_utils/Tools/testcase_end
M       ccsm_utils/Tools/cesm_buildexe
M       ccsm_utils/Tools/cesm_clean_build
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/ConfigCase.pm
M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Testcases/PEM_auto_pes_file
M       ccsm_utils/Testcases/NCK_build.csh
M       ccsm_utils/Testcases/PEM_build.csh
M       ccsm_utils/Testcases/PET_script
M       ccsm_utils/Testcases/SEQ_script
M       ccsm_utils/Testlists/C01
M       sample_pes_file.xml
M       create_newcase

================================================================================
Originator: mvertens
Date: 21 Mar 2013
Model: scripts
Version: scripts4_130321
One-line: fixed for minor bug that appeared in scripts4_130320

M       ccsm_utils/Case.template/config_grid.xml
	
================================================================================
Originator: mvertens
Date: 20 Mar 2013
Model: scripts
Version: scripts4_130320
One-line: fixes found for bugs that appeared in scripts4_130318
	
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Testlists/B41
M       ccsm_utils/Testlists/B42
M       ccsm_utils/Testlists/bluefire.cice2.auxtest
M       ccsm_utils/Testlists/C02
M       ccsm_utils/Testlists/namelists.frankfurt
M       create_newcase
	
================================================================================
Originator: mvertens
Date: 18 Mar 2013
Model: scripts
Version: scripts4_130318
One-line: Scripts refactoring for compsets and grids
	
D       sample_grid_file.xml
D       sample_compset_file.xml
D       ccsm_utils/Case.template/config_grid.xsl
D       ccsm_utils/Case.template/config_definition.xsl
D       ccsm_utils/Case.template/config_compsets.xsl
	
M       create_newcase
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/ConfigCase.pm
M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Tools/cesm_buildnml
	
A       ccsm_utils/Testlists/namelists.frankfurt
	
================================================================================
Originator: erik
Date: 15 Mar 2013
Model: scripts
Version: scripts4_130315c
One-line: Seperate out CLM intel and pgi test lists

D       ccsm_utils/Testlists/yellowstone.interactiveonly.clm.auxtest -- just delete
D       ccsm_utils/Testlists/yellowstone.clm.auxtest -- split into intel/pgi
D       ccsm_utils/Testlists/frankfurt.clm.auxtest ---- split into intel/pgi

New test names that are either wholely intel lists or pgi lists.

A  +    ccsm_utils/Testlists/frankfurt_intel.clm.auxtest
A  +    ccsm_utils/Testlists/yellowstone_intel.clm.auxtest
A  +    ccsm_utils/Testlists/yellowstone_pgi.clm.auxtest
A  +    ccsm_utils/Testlists/frankfurt_pgi.clm.auxtest

Add some comments about how to run to top of file.

M       ccsm_utils/Testlists/allIcompsetsRes.clm.auxtest

================================================================================
Originator: erik
Date: 15 Mar 2013
Model: scripts
Version: scripts4_130315b
One-line: Fix error in clm auxtest list, and add urban single point 
          RES_COMPSET_MATCH for CLM45 compsets

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Testlists/yellowstone.clm.auxtest

================================================================================
Originator: erik
Date: 15 Mar 2013
Model: scripts
Version: scripts4_130315
One-line: Change RES_COMPSET_MATCH for I_2000_1PTFRC, and I_RCP*_CN compsets

Fix bug 1643. RES_COMPSET_MATCH for some I compsets was not being used
because create_newcase was changed to require an exact match.

M       ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: jshollen
Date: 13 Mar 2013
Model: scripts
Version: scripts4_130313
One-line: Remove Windows carriage returns from st_archive.sh

M       ccsm_utils/Tools/st_archive.sh

================================================================================
Originator: jshollen
Date: 12 Mar 2013
Model: scripts
Version: scripts4_130312
One-line: move frankfurt prebeta test lists to prealpha

D       ccsm_utils/Testlists/frankfurt_pgi.prebeta
A  +    ccsm_utils/Testlists/frankfurt_intel.prealpha
A  +    ccsm_utils/Testlists/frankfurt_pgi.prealpha
D       ccsm_utils/Testlists/frankfurt_intel.prebeta


================================================================================
Originator: jshollen
Date: 08 Mar 2013
Model: scripts
Version: scripts4_130308
One-line: updated create_test_suite to also compare *.rc in offline namelist compare

M       create_test_suite

================================================================================
Originator: jedwards
Date: 06 Mar 2013
Model: scripts
Version: scripts4_130306
One-line: updated st_archive script for multi-instance cases (From Tim Hoar)

M   ccsm_utils/Tools/st_archive.sh


===========================
Originator: jshollen
Date: 04 Mar 2013
Model: scripts
Version: scripts4_130304
One-line: added CECO T62_g16 restart test

M       ccsm_utils/Testlists/C01

================================================================================
Originator: mlevy
Date: 01 Mar 2013
Model: scripts
Version: scripts4_130301
One-line: Update preview_namelist to copy all multi-instance namelists to
          CaseDocs (without reintroducing the "copy *_in* bug)

M       ccsm_utils/Tools/preview_namelists

================================================================================
Originator: erik
Date: 27 Feb 2013
Model: scripts
Version: scripts4_130227b
One-line: Combine CLM40 and CLM45 test lists, update PTCLM, CLM initial files for CLM40 compsets only

Update PTCLM. Now basic usages works when NOT creating input files.

++++++++++++++ CLM Test lists that combine clm4_0 and clm4_5
A       ccsm_utils/Testlists/yellowstone.interactiveonly.clm.auxtest
A       ccsm_utils/Testlists/yellowstone.clm.auxtest
A       ccsm_utils/Testlists/frankfurt.clm.auxtest
A       ccsm_utils/Testlists/shortlist.clm.auxtest

++++++++++++++ Delete CLM test lists that are specific for clm4_0 or clm4_5
++++++++++++++ Also delete unused test lists for titan, intrepid, lynx, bluefire.
D       ccsm_utils/Testlists/titan.clm.auxtest
D       ccsm_utils/Testlists/intrepid.clm.auxtest
D       ccsm_utils/Testlists/bluefire.clm.trans.auxtest
D       ccsm_utils/Testlists/lynx.clm.auxtest
D       ccsm_utils/Testlists/bluefire.clm40.auxtest
D       ccsm_utils/Testlists/bluefire.clm45.auxtest
D       ccsm_utils/Testlists/lynx.clm.interactive.auxtest
D       ccsm_utils/Testlists/yellowstone.clm40.auxtest
D       ccsm_utils/Testlists/yellowstone.clm45.auxtest
D       ccsm_utils/Testlists/shortlist.clm40.auxtest
D       ccsm_utils/Testlists/shortlist.clm45.auxtest
D       ccsm_utils/Testlists/yellowstone.clm40.interactive.auxtest
D       ccsm_utils/Testlists/frankfurt.clm40.auxtest
D       ccsm_utils/Testlists/bluefire.clm.interactive.auxtest
D       ccsm_utils/Testlists/yellowstone.clm45.interactive.auxtest
D       ccsm_utils/Testlists/frankfurt.clm45.auxtest

++++++++++++++ Use CLM initial files for CLM4_0 CN compsets only
++++++++++++++ Add CLM45 compsets to all I compsets test list
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Testlists/allIcompsetsRes.clm.auxtest

================================================================================
Originator: jedwards
Date: 27 Feb 2013
Model: scripts
Version: scripts4_130227
One-line: Improve output of compare_namelists.pl

M    ccsm_utils/Tools/compare_namelists.pl

===========================
Originator: jedwards
Date: 26 Feb 2013
Model: scripts
Version: scripts4_130226
One-line: Improve output of compare_namelists.pl

M    ccsm_utils/Tools/compare_namelists.pl

================================================================================
Originator: sacks
Date: 22 Feb 2013
Model: scripts
Version: scripts4_130222
One-line: Rework component_gen_comp to handle multiple history files
          per component

Now handles, e.g., both .h0 and .h1 file for CLM

M       ccsm_utils/Tools/component_compare.sh
M       ccsm_utils/Tools/component_gen_comp
M       ccsm_utils/Tools/component_generate.sh

================================================================================
Originator: fvitt
Date: 21 Feb 2013
Model: scripts
Version: scripts4_130221a
One-line: Corrected hsi check in copy_dirs_hsi section of lt_archive

M       ccsm_utils/Tools/lt_archive.sh

================================================================================
Originator: jshollen
Date: 21 Feb 2013
Model: scripts
Version: scripts4_130221
One-line: namelist compare fixed for the test scripts. 

M       ccsm_utils/Tools/testcase_end

================================================================================
Originator: jshollen
Date: 20 Feb 2013
Model: scripts
Version: scripts4_130220
One-line: add offline namelist comparison to create_test_suite

M       create_test_suite

================================================================================
Originator: jshollen
Date: 07 Feb 2013
Model: scripts
Version: scripts4_130207
One-line: updated testreporter, reports on namelist compares. 

M       ccsm_utils/Tools/testreporter.pl

================================================================================
Originator: sacks
Date: 06 Feb 2013
Model: scripts
Version: scripts4_130206
One-line: changes to support cism2

- add cism2 compsets, and rename cism1 compsets
- add trilinos support in build
- add additional glc grid support
- add SourceMods/src.cism/glimmer-cism subdirectory, needed now that
  glimmer-cism code is built as a separate library from the
  cesm-specific code
- rework test lists to test cism2 compsets

M       create_newcase
M       ccsm_utils/Tools/cesm_buildexe
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/config_definition.xml
A       ccsm_utils/Testlists/titan.glc.auxtest
A       ccsm_utils/Testlists/yellowstone_intel.glc.auxtest
A       ccsm_utils/Testlists/hopper_pgi.glc.auxtest
A       ccsm_utils/Testlists/yellowstone_pgi.glc.auxtest
D       ccsm_utils/Testlists/bluefire.glc.auxtest
A       ccsm_utils/Testlists/hopper_gnu.glc.auxtest
M       ccsm_utils/Testlists/A01
M       ccsm_utils/Testlists/shortlist.glc.auxtest
M       doc/modelnl/nldef2html_cism

================================================================================
Originator: jedwards
Date: 05 Feb 2013
Model: scripts
Version: scripts4_130205a
One-line: Remove redundant test output

M       ccsm_utils/Tools/testcase_end

================================================================================
riginator: jedwards
Date: 05 Feb 2013
Model: scripts
Version: scripts4_130205
One-line:  add compare_namelists.pl to compare namelists against baseline

A       ccsm_utils/Tools/compare_namelists.pl
M      ccsm_utils/Tools/testcase_end
M      create_newcase
	

===========================
Originator: jshollen
Date: 04 Feb 2013
Model: scripts
Version: scripts4_130204
One-line:  New version of testreporter.pl, uploads test results to the testdb.

M       ccsm_utils/Tools/testreporter.pl

================================================================================
Originator: erik
Date: 30 Jan 2013
Model: scripts
Version: scripts4_130130
One-line:  Add user_nl_dir to testnames, add RTM_BLDNML_OPTS/RTM_NAMELIST_OPTS

M       create_test --- add user_nl_dir to testname
M       SVN_EXTERNAL_DIRECTORIES -- Update PTCLM
M       ccsm_utils/Case.template/config_compsets.xml -- add placeholder for U compset
              add RTM_NAMELIST_OPTS for CLM45 compsets, fix I20TRCLM45 compset
M       ccsm_utils/Case.template/config_definition.xml - add RTM_BLDNML_OPTS/RTM_NAMELIST_OPTS

Needs rtm1_0_17 or later

================================================================================
Originator: jedwards
Date: 28 Jan 2013
Model: scripts
Version: scripts4_130128
One-line: add namelist consistantcy check to tests

M      ccsm_utils/Tools/testcase_begin
M      ccsm_utils/Tools/testcase_end

================================================================================
Originator: sacks
Date: 27 Jan 2013
Model: scripts
Version: scripts4_130127
One-line: Change check_exactrestart.pl to fail when no comm_diag lines are found

M       ccsm_utils/Tools/check_exactrestart.pl

================================================================================
Originator: sacks
Date: 25 Jan 2013
Model: scripts
Version: scripts4_130125
One-line: update component_gen_comp for yellowstone defaults;
          use cprnc in user's path if possible

M       ccsm_utils/Tools/component_compare.sh
M       ccsm_utils/Tools/component_gen_comp

================================================================================
Originator: jshollen
Date: 09 Jan 2013
Model: scripts
Version: scripts4_130109
One-line: cesm_setup now tracks BASELINE_ROOT for tests, fixes NCK test 

M       ccsm_utils/Tools/cesm_setup
M       ccsm_utils/Case.template/config_definition.xml
M       create_newcase

================================================================================
Originator: sacks
Date: 08 Jan 2013
Model: scripts
Version: scripts4_130108
One-line: Update perl5lib to allow '%glc' in data model field names

 M      .
M       SVN_EXTERNAL_DIRECTORIES

================================================================================
Originator: sacks
Date: 07 Jan 2013
Model: scripts
Version: scripts4_130107
One-line: Fix documentation

M       doc/usersguide/introduction.xml

================================================================================
Originator: mmvertens
Date: 03 Jan 2013
Model: scripts
Version: scripts4_130103b
One-line: New compsets and clm aux tests for clm4_5 (coming into clm4_0_60)
	  Results should be bit-for-bit for standard clm4_0 
	
M       ccsm_utils/Case.template/config_compsets.xml
D       ccsm_utils/Testlists/yellowstone.clm.auxtest
D       ccsm_utils/Testlists/yellowstone.clm.interactive.auxtest
A       ccsm_utils/Testlists/bluefire.clm40.auxtest
D       ccsm_utils/Testlists/frankfurt.clm.auxtest
A       ccsm_utils/Testlists/bluefire.clm45.auxtest
A       ccsm_utils/Testlists/yellowstone.clm40.auxtest
D       ccsm_utils/Testlists/bluefire.clm.auxtest
A       ccsm_utils/Testlists/yellowstone.clm45.auxtest
A       ccsm_utils/Testlists/shortlist.clm40.auxtest
A       ccsm_utils/Testlists/shortlist.clm45.auxtest
A       ccsm_utils/Testlists/yellowstone.clm40.interactive.auxtest
A       ccsm_utils/Testlists/frankfurt.clm40.auxtest
A       ccsm_utils/Testlists/yellowstone.clm45.interactive.auxtest
A       ccsm_utils/Testlists/frankfurt.clm45.auxtest
D       ccsm_utils/Testlists/shortlist.clm.auxtest

================================================================================
Originator: mlevy
Date: 03 Jan 2013
Model: scripts
Version: scripts4_130103
One-line: Bug-fix for f25_f25 resolution (OCN_GRID was incorrect)

M       ccsm_utils/Case.template/config_grid.xml

================================================================================
Originator: jedwards
Date: 01 Jan 2013
Model: scripts
Version: scripts4_130101
One-line: Refactor check_exactrestart.pl to give more relivent information on failure

M ccsm_utils/Tools/check_exactrestart.pl


================================================================================
Originator: santos
Date: 27 Dec 2012
Model: scripts
Version: scripts4_121227
One-line: Simple fix for fix WACCM5 cases

M       ccsm_utils/Case.template/config_compsets.xml
        - Add "-phys cam5" back to WACCM5 match. This
          would be unnecessary due to the default, except
          that earlier "WACCM" matches currently set "-phys cam4",
          so that needs to be overridden.

================================================================================
Originator: fvitt
Date: 19 Dec 2012
Model: scripts
Version: scripts4_121219
One-line: Updates for pleiades

M       ccsm_utils/Tools/lt_archive.sh

================================================================================
Originator: mvertens
Date: 18 Dec 2012
Model: scripts
Version: scripts4_121218b
One-line: Updates for clm auxilliary tests

M       ccsm_utils/Testlists/yellowstone.clm.auxtest

================================================================================
Originator: tcraig
Date: 18 Dec 2012
Model: scripts
Version: scripts4_121218
One-line: Update machine name parsing for timing file

M       ccsm_utils/Tools/timing/getTiming2.pl	
================================================================================
Originator: tcraig
Date: 17 Dec 2012
Model: scripts
Version: scripts4_121217b
One-line: Update cost calculation for timing output

M       ccsm_utils/Tools/timing/getTiming.csh
M       ccsm_utils/Tools/timing/getTiming2.pl
M       ccsm_utils/Case.template/config_definition.xml
================================================================================
Originator: jshollen
Date: 17 Dec 2012
Model: scripts
Version: scripts4_121217a
One-line: Merge changes from release branch, config_compsets & release test changes

A       ccsm_utils/Testlists/B44
A       ccsm_utils/Testlists/C47
M       ChangeLog
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xsl
M       ccsm_utils/Testlists/B01
M       ccsm_utils/Testlists/B42
M       ccsm_utils/Testlists/B43
M       ccsm_utils/Testlists/C01
M       ccsm_utils/Testlists/C02
M       ccsm_utils/Testlists/C03
M       ccsm_utils/Testlists/C05
M       ccsm_utils/Testlists/C41
M       ccsm_utils/Testlists/E01
M       ccsm_utils/Testlists/bluefire.cam.auxtest
M       ccsm_utils/Testlists/bluefire.clm.auxtest
M       ccsm_utils/Testlists/bluefire.prebeta
M       ccsm_utils/Testlists/lynx.clm.auxtest
M       ccsm_utils/Testlists/lynx_intel.prebeta
M       ccsm_utils/Testlists/lynx_pgi.prebeta
M       ccsm_utils/Testlists/titan.clm.auxtest
M       ccsm_utils/Testlists/yellowstone.clm.auxtest

================================================================================
Originator: santos
Date: 17 Dec 2012
Model: scripts
Version: scripts4_121217
One-line: Fix bug where data sets get pulled into CaseDocs in Intel tests.

M       ccsm_utils/Tools/preview_namelists
        - Change "*_in*" to just "*_in", preventing datasets with "*_intel"
          from being pulled into CaseDocs.

================================================================================
Originator: erik
Date: 13 Dec 2012
Model: scripts
Version: scripts4_121213
One-line: Bring changes to xmlquery from cesm1_1_0 release branch onto trunk

Add CROP compset and add clm test lists for yellowstone

M       ccsm_utils/Tools/xmlquery
M       ccsm_utils/Tools/SetupTools.pm
M       ccsm_utils/Case.template/config_compsets.xml
A       ccsm_utils/Testlists/yellowstone.clm.auxtest
A       ccsm_utils/Testlists/yellowstone.clm.interactive.auxtest

================================================================================
Originator: jedwards
Date: 11 Dec 2012
Model: scripts
Version: scripts4_121211
One-line: remove invalid test from lists (waccm not in cam-se)

M   ccsm_utils/Testlists/B42

================================================================================
Originator: santos
Date: 7 Dec 2012
Model: scripts
Version: scripts4_121207b
One-line: Add WACCM5 compsets.

M       ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: jshollen
Date: 7 Dec 2012
Model: scripts
Version: scripts4_121207
One-line: cesm_setup can now regenerate .test script

M       ccsm_utils/Case.template/ConfigCase.pm
M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Tools/cesm_setup
M       create_newcase
M       create_test

================================================================================
Originator: jshollen
Date: 3 Dec 2012
Model: scripts
Version: scripts4_121203
One-line: clone yellowstone intel test list for pgi. 

A  +    ccsm_utils/Testlists/yellowstone_intel.prealpha
A  +    ccsm_utils/Testlists/yellowstone_intel.prebeta
A       ccsm_utils/Testlists/yellowstone_pgi.prealpha
A       ccsm_utils/Testlists/yellowstone_pgi.prebeta
D       ccsm_utils/Testlists/yellowstone.prealpha
D       ccsm_utils/Testlists/yellowstone.prebeta

================================================================================
Originator: fvitt

Date: 29 Nov 2012
Model: scripts
Version: scripts4_121129
One-line: Changes for pleiades

M       ccsm_utils/Tools/taskmaker.pl

================================================================================
Originator: mvertens
	
Date: 27 Nov 2012
Model: scripts
Version: scripts4_121127
One-line: changes for CRU NCEP forcing and compsets

M       ccsm_utils/Case.template/config_grid.xsl
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Testlists/titan.clm.auxtest
M       ccsm_utils/Testlists/lynx.clm.auxtest
M       ccsm_utils/Testlists/bluefire.clm.auxtest
	
================================================================================
Originator: tcraig
Date: 25 Nov 2012
Model: scripts
Version: scripts4_121125
One-line: fix some test issues

 - add RTM_MODE and RTM_FLOOD_MODE
 - fix problems with xmlchange and testing
 - fix IOP test issues with multiple builds
 - modify perl5lib for drv map names all lower case

M       SVN_EXTERNAL_DIRECTORIES
M       ccsm_utils/Tools/testcase_begin
M       ccsm_utils/Tools/xmlchange
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/config_definition.xml
================================================================================
Originator: sacks
Date: 15 Nov 2012
Model: scripts
Version: scripts4_121115
One-line: Fix parsing of BEG_COMPSET_MATCH

Previously, BEG_COMPSET_MATCHes with vertical bars were not being
treated properly. I have fixed this. It appears that this only affects
the parsing of two compsets:

B_2000_MOZSOA_CN 
F_2000_MOZSOA_CN

And it only changes settings that are relevant for F_2000_MOZSOA_CN,
not B_2000_MOZSOA_CN.


M       create_newcase

================================================================================
Originator: mlevy
Date: 14 Nov 2012
Model: scripts
Version: scripts4_121114
One-line: New domain.ocn.ne* files

M       ccsm_utils/Case.template/config_grid.xml
================================================================================
Originator: tcraig
Date: 05 Nov 2012
Model: scripts
Version: scripts4_121105
One-line: Fix NCK test

M       ccsm_utils/Testcases/NCK_script
================================================================================
Originator: tcraig
Date: 01 Nov 2012
Model: scripts
Version: scripts4_121101b
One-line: Fix create_production_test

M       ccsm_utils/Tools/testcase_env.csh
		
================================================================================
Originator: jedwards
Date: 01 Nov 2012
Model: scripts
Version: scripts4_121101a
One-line: Fix issues found in last commit
M               create_clone
M               create_production_test
	
================================================================================
Originator: jedwards
Date: 01 Nov 2012
Model: scripts
Version: scripts4_121101
One-line: Create SetupTools.pm to avoid duplication of several resused perl functions
M               create_clone
M               ccsm_utils/Tools/cesm_setup
M               ccsm_utils/Tools/xmlquery
A               ccsm_utils/Tools/SetupTools.pm
M               create_newcase


================================================================================
Originator: tcraig
Date: 30 Oct 2012
Model: scripts
Version: scripts4_121030
One-line: ICP fix
	
M       ccsm_utils/Testcases/ICP_build.csh
	
================================================================================
Originator: tcraig
Date: 29 Oct 2012
Model: scripts
Version: scripts4_121029
One-line: add test_build, rename to cesm_setup and env_mach_pes.xml, update userdefined, change F compset to ahve rof grid null

svn merge $SVNREPO/scripts/trunk_tags/scripts4_121025 $SVNREPO/scripts/branch_tags/tbld_tags/tbld05_scripts4_121025
 - add test_build script to generate cesm executables up front for tests
 - refactor tests for new test build, update/remove some tests
 - rename setup to cesm_setup, rename env_pesetup.xml to env_mach_pes.xml
 - change F compset rof grid default to null
 - modify cs.submit for test build and other future extensions to it
 - update userdefined setup
 - modify how mpilib and mpi-serial defaults are set, add _M confopts
	
M       create_test
M       create_clone
M       ccsm_utils/Tools/testcase_begin
M       ccsm_utils/Tools/ccsm_check_lockedfiles
M       ccsm_utils/Tools/archive_metadata.sh
A  +    ccsm_utils/Tools/cesm_setup
M       ccsm_utils/Tools/cesm_prestage
M       ccsm_utils/Tools/cesm_buildstart
M       ccsm_utils/Tools/xml2env
MM      ccsm_utils/Tools/taskmaker.pl
MM      ccsm_utils/Tools/ccsm_getenv
D       ccsm_utils/Tools/setup
M       ccsm_utils/Tools/testcase_end
M       ccsm_utils/Tools/xmlchange
M       ccsm_utils/Tools/testcase_setup.csh
X       ccsm_utils/Tools/lnd/clm/PTCLM
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/ConfigCase.pm
MM      ccsm_utils/Case.template/config_definition.xml
D       ccsm_utils/Testcases/PMT_script
D       ccsm_utils/Testcases/SEQ_auto_pes_file
D       ccsm_utils/Testcases/PET_auto_pes_file
M       ccsm_utils/Testcases/PEA_script
D       ccsm_utils/Testcases/PST_script
M       ccsm_utils/Testcases/ERB_script
A  +    ccsm_utils/Testcases/tests_build.csh
D       ccsm_utils/Testcases/PMT_auto_pes_file
M       ccsm_utils/Testcases/ICP_script
M       ccsm_utils/Testcases/NCK_script
M       ccsm_utils/Testcases/ERH_script
D       ccsm_utils/Testcases/PST_auto_pes_file
M       ccsm_utils/Testcases/PEM_script
M       ccsm_utils/Testcases/APT_script
M       ccsm_utils/Testcases/config_tests.xml
M       ccsm_utils/Testcases/CME_script
A  +    ccsm_utils/Testcases/CME_build.csh
A  +    ccsm_utils/Testcases/PEA_build.csh
M       ccsm_utils/Testcases/ERI_script
A  +    ccsm_utils/Testcases/NCK_build.csh
A  +    ccsm_utils/Testcases/ICP_build.csh
A  +    ccsm_utils/Testcases/PEM_build.csh
M       ccsm_utils/Testcases/PET_script
M       ccsm_utils/Testcases/SEQ_script
M       ccsm_utils/Testlists/titan.clm.auxtest
M       ccsm_utils/Testlists/eastwind.rasm.auxtest
M       ccsm_utils/Testlists/B01
M       ccsm_utils/Testlists/B41
M       ccsm_utils/Testlists/B42
M       ccsm_utils/Testlists/B43
M       ccsm_utils/Testlists/lynx.clm.auxtest
M       ccsm_utils/Testlists/olympus.rasm.auxtest
M       ccsm_utils/Testlists/evergreen.rasm.auxtest
M       ccsm_utils/Testlists/bluefire.glc.auxtest
M       ccsm_utils/Testlists/bluefire.clm.auxtest
M       ccsm_utils/Testlists/bluefire.cice2.auxtest
M       ccsm_utils/Testlists/A01
M       ccsm_utils/Testlists/C01
M       ccsm_utils/Testlists/shortlist.clm.auxtest
M       ccsm_utils/Testlists/S01
M       create_newcase
M       create_test_suite

================================================================================
Originator: jedwards
Date: 25 Oct 2012
Model: scripts
Version: scripts4_121025
One-line: make ne120 rof r05, change intrepid tests from f19 to ne30

M             ccsm_utils/Case.template/config_grid.xml
M             ccsm_utils/Testlists/C05

================================================================================
Originator: tcraig
Date: 23 Oct 2012
Model: scripts
Version: scripts4_121023
One-line: Fix ICP test
	
M       ccsm_utils/Tools/cesm_prerun_setup
M       ccsm_utils/Testcases/ICP_script
================================================================================
Originator: sacks
Date: 22 Oct 2012
Model: scripts
Version: scripts4_121022
One-line: add "expand all help" and "collapse all help" buttons on
          namelist documentation pages

*** Add needed functions
M       doc/modelnl/showinfo.js

*** Add buttons to all pages
M       doc/modelnl/nldef2html_cam
M       doc/modelnl/nldef2html_cice
M       doc/modelnl/nldef2html_cism
M       doc/modelnl/nldef2html_clm
M       doc/modelnl/nldef2html_datm
M       doc/modelnl/nldef2html_dice
M       doc/modelnl/nldef2html_dlnd
M       doc/modelnl/nldef2html_docn
M       doc/modelnl/nldef2html_drof
M       doc/modelnl/nldef2html_drv
M       doc/modelnl/nldef2html_pop2
M       doc/modelnl/nldef2html_rtm
M       doc/modelnl/xmldef2html_compsets
M       doc/modelnl/xmldef2html_env_build
M       doc/modelnl/xmldef2html_env_case
M       doc/modelnl/xmldef2html_env_pesetup
M       doc/modelnl/xmldef2html_env_run
M       doc/modelnl/xmldef2html_grid
M       doc/modelnl/xmldef2html_machines

================================================================================
Originator: jshollen
Date: 15 Oct 2012
Model: scripts
Version: scripts4_121015a
One-line: add ERS.ne30_g16.B1850C5CN to ys prebeta list, make yellowstone prealpha/prebeta lists. 

A       ccsm_utils/Testlists/yellowstone.prealpha
A       ccsm_utils/Testlists/B43
A       ccsm_utils/Testlists/yellowstone.prebeta

================================================================================
Originator: jedwards
Date: 15 Oct 2012
Model: scripts
Version: scripts4_121015
One-line: add a hybrid startup for ne30_g16.B1850C5CN

M   ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: sacks
Date: 12 Oct 2012
Model: scripts
Version: scripts4_121012
One-line: Break cism namelist documentation into more groups

(Also, many doc changes from others, not documented here)

M       doc/modelnl/nldef2html_cism

================================================================================
Originator: sacks
Date: 09 Oct 2012
Model: scripts
Version: scripts4_121009b
One-line: Remove test from bluefire.glc.auxtest list

I had added a test temporarily; I can remove it now that the main ERS
tests are passing again.

M       ccsm_utils/Testlists/bluefire.glc.auxtest

================================================================================
Originator: jedwards
Date: 09 Oct 2012
Model: scripts
Version: scripts4_121009a
One-line: Add back capability to specify mpilib in create_test and
create_test_suite (needed for yellowstone)

M     create_test
M     create_test_suite


================================================================================
Originator: tcraig
Date: 09 Oct 2012
Model: scripts
Version: scripts4_121009
One-line: Fix use of EXEROOT aprun option in taskmaker.pl

M       ccsm_utils/Tools/taskmaker.pl
	
================================================================================
Originator: erik
Date: 8 Oct 2012
Model: scripts
Version: scripts4_121008b
One-line: Fix call to expandXMLVars for check_input_data when var needs to be expanded

M       ccsm_utils/Tools/check_input_data

================================================================================
Originator: sacks
Date: 8 Oct 2012
Model: scripts
Version: scripts4_121008
One-line: update perl5lib

 M      .
M       SVN_EXTERNAL_DIRECTORIES

================================================================================
Originator: jshollen
Date: 4 Oct 2012
Model: scripts
Version: scripts4_121004
One-line: Merge & update hopper's prebeta test list with jaguarpf 

A       ccsm_utils/Testlists/B42
M       ccsm_utils/Testlists/hopper.prebeta

================================================================================
Originator: sacks
Date: 3 Oct 2012
Model: scripts
Version: scripts4_121003
One-line: add note that STOP_OPTION & REST_OPTION must be nyears for _GLC cases
          
M       ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: jim
Date: 2 Oct 2012
Model: scripts
Version: scripts4_121002a
One-line: Allow for variables of the form $ENV{Variable_Name} in config_machines.xml
M  ccsm_utils/Case.template/ConfigCase.pm

================================================================================
Originator: erik
Date: 2 Oct 2012
Model: scripts
Version: scripts4_121002
One-line: Work with clm testlists, update B1850 case for T31_g37, add expandXMLVar
          call for check_input_data if needed, update perl5lib

M       ccsm_utils/Tools/check_input_data
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Testlists/bluefire.clm.trans.auxtest
M       ccsm_utils/Testlists/allIcompsetsRes.clm.auxtest

================================================================================
Originator: tcraig
Date: 1 Oct 2012
Model: scripts
Version: scripts4_121001a
One-line: fix testid issue with IOP test

M       create_test
M       ccsm_utils/Tools/testcase_begin
M       ccsm_utils/Tools/testcase_env.csh
================================================================================
Originator: sacks
Date: 1 Oct 2012
Model: scripts
Version: scripts4_121001
One-line: point to new forcing data locations for TG compsets

new locations have files with naming convention changed to match
what's expected by new dlnd

M       ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: sacks
Date: 30 Sep 2012
Model: scripts
Version: scripts4_120930
One-line: add test in glc test list

M       ccsm_utils/Testlists/bluefire.glc.auxtest

================================================================================
Originator: santos
Date: 28 Sep 2012
Model: scripts
Version: scripts4_120928f
One-line: WACCM-X 1996 compset, fixed setting of *_OPTS in env_build.xml

NOTE: This was intended to be scripts4_120928d, but, by mistake, no change
      was committed in that tag.

M       ccsm_utils/Case.template/config_compsets.xml
 - Add WACCM-X 1996 solar minimum compset.
M       create_newcase
 - Change handling of GEN_COMPSET_MATCH to better handle cases where multiple
   options are specified multiple times.
 - New WACCM-X specification relies on this behavior, but no change to existing
   compsets.

================================================================================
Originator: sacks
Date: 28 Sep 2012
Model: scripts
Version: scripts4_120928e
One-line: add note that BG1850CN is broken at f09

M       ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: sacks
Date: 28 Sep 2012
Model: scripts
Version: scripts4_120928c
One-line: change default GLC grid to gland5UM

This will change answers for all _GLC compsets!

M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/config_definition.xml

M       ccsm_utils/Testlists/bluefire.glc.auxtest
M       ccsm_utils/Testlists/A01
M       ccsm_utils/Testlists/E01
M       ccsm_utils/Testlists/shortlist.glc.auxtest
        - remove explicit _gl5UM designation from grids

================================================================================
Originator: sacks
Date: 28 Sep 2012
Model: scripts
Version: scripts4_120928b
One-line: change grid name for gl5UM, fix glc test lists

M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Testlists/bluefire.glc.auxtest
M       ccsm_utils/Testlists/shortlist.glc.auxtest

================================================================================
Originator: tcraig
Date: 28 Sep 2012
Model: scripts
Version: scripts4_120928
One-line: lots of little updates

       - modify some testcases for better coverage
       - gx3v7 grid switched to domain.ocn.gx3v7.120323.nc (db)
       - get rid of cesm_buildnml
       - fix duplication of code in preview_namelists
       - update preview_namelists, add -verbose and modify LID setting
       - update create_production_test (je)
       - fix compiler validation check in create_newcase (mv)
       - remove rpointer.ocn.init from cesm_prestage logic (ml)
       - update IOP test
       - fix ERT, remove cpl log file check, hist file only
       - fix check_exactrestart.pl to not fail automatically for no comm_diag 
         output (needed for S tests)
       - update grids and compsets for glc and rof grid support
       - add a few more ETEST cases

M       create_test
M       create_clone
M       ccsm_utils/Tools/testcase_begin
M       ccsm_utils/Tools/testcase_env.csh
M       ccsm_utils/Tools/preview_namelists
M       ccsm_utils/Tools/cesm_prerun_setup
M       ccsm_utils/Tools/cesm_buildnml
M       ccsm_utils/Tools/check_exactrestart.pl
M       ccsm_utils/Tools/cesm_prestage
M       ccsm_utils/Tools/cesm_buildstart
M       ccsm_utils/Tools/testcase_end
M       ccsm_utils/Tools/create_production_test
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Testcases/PMT_script
M       ccsm_utils/Testcases/SMS_script
M       ccsm_utils/Testcases/CME_script
M       ccsm_utils/Testcases/PEA_script
M       ccsm_utils/Testcases/PST_script
M       ccsm_utils/Testcases/ERB_script
M       ccsm_utils/Testcases/NCK_script
M       ccsm_utils/Testcases/ICP_script
M       ccsm_utils/Testcases/LAR_script
M       ccsm_utils/Testcases/ERH_script
M       ccsm_utils/Testcases/ERI_script
M       ccsm_utils/Testcases/PEM_script
M       ccsm_utils/Testcases/SBN_script
M       ccsm_utils/Testcases/APT_script
M       ccsm_utils/Testcases/P4A_script
M       ccsm_utils/Testcases/PFS_script
M       ccsm_utils/Testcases/PET_script
M       ccsm_utils/Testcases/SEQ_script
M       ccsm_utils/Testcases/ERS_script
M       ccsm_utils/Testcases/ERT_script
M       ccsm_utils/Testlists/B01
D       ccsm_utils/Testlists/bluefire.posttag.old
D       ccsm_utils/Testlists/bluefire.pretag.old
M       ccsm_utils/Testlists/bluefire.glc.auxtest
M       ccsm_utils/Testlists/A01
M       ccsm_utils/Testlists/C01
M       ccsm_utils/Testlists/C02
M       ccsm_utils/Testlists/E01
M       ccsm_utils/Testlists/shortlist.glc.auxtest
M       create_newcase
M       create_test_suite
================================================================================
Originator: mvertens
Date: 24 Sep 2012
Model: scripts
Version: scripts4_120924a
One-line: setup fixes so that Macros and user_nl_xxx can be reset without
	  invoking setup -clean
	
M       ccsm_utils/Tools/setup
	
================================================================================
Originator: tcraig
Date: 24 Sep 2012
Model: scripts
Version: scripts4_120924
One-line: change D and E tests to DTEST and ETEST

M       ccsm_utils/Testlists/B01
M       ccsm_utils/Testlists/B41
M       ccsm_utils/Testlists/newmachine.port.auxtest
M       ccsm_utils/Testlists/C01
M       ccsm_utils/Testlists/C02
M       ccsm_utils/Testlists/E01
M       ccsm_utils/Testlists/C03
M       ccsm_utils/Testlists/C04
M       ccsm_utils/Testlists/C05
================================================================================
Originator: sacks
Date: 23 Sep 2012
Model: scripts
Version: scripts4_120923a
One-line: update shortlist.glc.auxtest

M       ccsm_utils/Testlists/shortlist.glc.auxtest
        - update SMS_D IG test to be consistent with change in bluefire.glc.auxtest

================================================================================
Originator: sacks
Date: 23 Sep 2012
Model: scripts
Version: scripts4_120923
One-line: update test lists, minor bug fixes

M       ccsm_utils/Tools/component_gen_comp
        - bug fix for multi-instance

M       ccsm_utils/Case.template/config_compsets.xml
        - fix BG1850CN to use gland10 (this setting was accidentally
          removed at some point)

M       ccsm_utils/Testlists/bluefire.glc.auxtest
        - add TG CME and ERS_E test
        - replace SMS IGLONG test with CME_Ly5 IG test
        - remove unnecessary ERS_Lm3 test, change SMS_D test to use
          IG20TR rather than IG (so that we still have an IG20TR test)
        - fix syntax in ERS20y -> ERS_Ly20

M       ccsm_utils/Testlists/E01
        - add TG CME test, change some B tests to BG

================================================================================
Originator: tcraig
Date: 21 Sep 2012
Model: scripts
Version: scripts4_120921a
One-line: update tests to support variable length confopts

- this will change some test results
- delete several tests in favor of _L option
- update and clean up some test lists
- update clean_build
- add DTEST, ETEST
- add DOCN_SOME_FILENAME
- add DROF env variables

M       ccsm_utils/Tools/testcase_env.csh
M       ccsm_utils/Tools/cesm_clean_build
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_definition.xml
D       ccsm_utils/Testcases/ERU_script
D       ccsm_utils/Testcases/ERS3y_script
M       ccsm_utils/Testcases/ERB_script
M       ccsm_utils/Testcases/ICP_script
M       ccsm_utils/Testcases/ERH_script
D       ccsm_utils/Testcases/ERS211d_script
D       ccsm_utils/Testcases/SMS2_script
D       ccsm_utils/Testcases/ERP_script
M       ccsm_utils/Testcases/PFS_script
D       ccsm_utils/Testcases/SMS6_script
M       ccsm_utils/Testcases/ERT_script
D       ccsm_utils/Testcases/CME10y_script
M       ccsm_utils/Testcases/config_tests.xml
M       ccsm_utils/Testcases/CME_script
D       ccsm_utils/Testcases/ERS2_script
M       ccsm_utils/Testcases/LAR_script
M       ccsm_utils/Testcases/ERI_script
D       ccsm_utils/Testcases/ERS6_script
D       ccsm_utils/Testcases/ERI44y_script
M       ccsm_utils/Testcases/P4A_script
D       ccsm_utils/Testcases/ERS20y_script
D       ccsm_utils/Testcases/CME4_script
D       ccsm_utils/Testcases/ERS48s_script
M       ccsm_utils/Testcases/SEQ_script
M       ccsm_utils/Testcases/ERS_script
A  +    ccsm_utils/Testlists/frankfurt_pgi.prebeta
M       ccsm_utils/Testlists/titan.clm.auxtest
D       ccsm_utils/Testlists/franklin.prebeta
A  +    ccsm_utils/Testlists/frankfurt_pgi.prerelease
D       ccsm_utils/Testlists/edinburgh_pgi.prerelease
M       ccsm_utils/Testlists/hadley.auxtest
A  +    ccsm_utils/Testlists/frankfurt_lahey.prerelease
M       ccsm_utils/Testlists/eastwind.rasm.auxtest
M       ccsm_utils/Testlists/B41
D       ccsm_utils/Testlists/edinburgh_pgi.prebeta
D       ccsm_utils/Testlists/edinburgh_intel.prerelease
D       ccsm_utils/Testlists/franklin.prerelease
M       ccsm_utils/Testlists/newmachine.port.auxtest
D       ccsm_utils/Testlists/atlas.auxtest
M       ccsm_utils/Testlists/lynx.clm.auxtest
M       ccsm_utils/Testlists/olympus.rasm.auxtest
A  +    ccsm_utils/Testlists/frankfurt.clm.auxtest
M       ccsm_utils/Testlists/evergreen.rasm.auxtest
D       ccsm_utils/Testlists/edinburgh_intel.prebeta
M       ccsm_utils/Testlists/bluefire.posttag.old
M       ccsm_utils/Testlists/bluefire.pretag.old
M       ccsm_utils/Testlists/bluefire.glc.auxtest
M       ccsm_utils/Testlists/bluefire.clm.auxtest
D       ccsm_utils/Testlists/edinburgh_lahey.prebeta
M       ccsm_utils/Testlists/bluefire.cice1.auxtest
M       ccsm_utils/Testlists/bluefire.cice2.auxtest
D       ccsm_utils/Testlists/edinburgh.clm.auxtest
A  +    ccsm_utils/Testlists/frankfurt_intel.prebeta
M       ccsm_utils/Testlists/bluefire.clm.interactive.auxtest
D       ccsm_utils/Testlists/edinburgh_lahey.prerelease
M       ccsm_utils/Testlists/A01
M       ccsm_utils/Testlists/A02
M       ccsm_utils/Testlists/bluefire.drv.auxtest
M       ccsm_utils/Testlists/A03
M       ccsm_utils/Testlists/C01
A  +    ccsm_utils/Testlists/frankfurt_intel.prerelease
M       ccsm_utils/Testlists/C02
M       ccsm_utils/Testlists/A04
M       ccsm_utils/Testlists/E01
M       ccsm_utils/Testlists/A05
M       ccsm_utils/Testlists/C03
M       ccsm_utils/Testlists/C04
M       ccsm_utils/Testlists/C41
M       ccsm_utils/Testlists/C05
M       ccsm_utils/Testlists/C42
D       ccsm_utils/Testlists/chester.auxtest
M       ccsm_utils/Testlists/C43
M       ccsm_utils/Testlists/bluefire.rasm.auxtest
M       ccsm_utils/Testlists/C44
M       ccsm_utils/Testlists/shortlist.glc.auxtest
M       ccsm_utils/Testlists/C46
A  +    ccsm_utils/Testlists/frankfurt_lahey.prebeta
================================================================================
Originator: jshollen
Date: 21 Sep 2012
Model: scripts
Version: scripts4_120921
One-line: bug fix in testreporter to handle non-existent test status files. 

M       ccsm_utils/Tools/testreporter.pl

================================================================================
Originator: mlevy
Date: 19 Sep 2012
Model: scripts
Version: scripts4_120919
Two-line: Bugfix in pop requires update to cesm_prestage (hybrid / branch runs
          now create rpointer.ocn.init for pop2.buildnml.csh)

M       ccsm_utils/Tools/cesm_prestage
================================================================================
Originator: tcraig
Date: 18 Sep 2012
Model: scripts
Version: scripts4_120918
One-line: runoff separate component

  svn merge $SVNREPO/scripts/trunk_tags/scripts4_120917 $SVNREPO/scripts/branch_tags/rtmcomp_tags/rtmcomp09_scripts4_120917
- add runoff component
- add _L confopts option for run length
- fix timing file for non "daily" coupling runs
- fix grep -a issue on bluefire
- remove "LONG" compsets, now handled by _L
- rof pes set to mach land pes initially, need to work on config_pes.xml file for rof
- update tests for rof (not fully tested)
	
M       create_test
M       sample_grid_file.xml
M       ccsm_utils/Tools/lt_archive.sh
M       ccsm_utils/Tools/st_archive.sh
M       ccsm_utils/Tools/timing/getTiming2.pl
MM      ccsm_utils/Tools/taskmaker.pl
MM      ccsm_utils/Tools/ccsm_getenv
M       ccsm_utils/Tools/setup
M       ccsm_utils/Tools/testcase_end
M       ccsm_utils/Tools/cesm_buildexe
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_compsets.xsl
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/ConfigCase.pm
MM      ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Testcases/PEM_auto_pes_file
M       ccsm_utils/Testcases/PMT_script
M       ccsm_utils/Testcases/SEQ_auto_pes_file
M       ccsm_utils/Testcases/PST_script
M       ccsm_utils/Testcases/PMT_auto_pes_file
M       ccsm_utils/Testcases/NCK_script
M       ccsm_utils/Testcases/PST_auto_pes_file
M       ccsm_utils/Testcases/PEM_script
M       ccsm_utils/Testcases/config_tests.xml
M       ccsm_utils/Testcases/SEQ_script
M       sample_pes_file.xml
M       create_newcase
M       sample_compset_file.xml
================================================================================
Originator: jshollen
Date: 17 Sep 2012
Model: scripts
Version: scripts4_120917
One-line: Fix compare bug in create_test_suite

M       create_test_suite

================================================================================
Originator: erik
Date: 15 Sep 2012
Model: scripts
Version: scripts4_120915
One-line: Fix user_nl_dir option in create_test

      M       create_test

================================================================================
Originator: jedwards
Date: 14 Sep 2012
Model: scripts
Version: scripts4_120914a
One-line: improve parsing of pecount in create_newcase

     M       ccsm_utils/Case.template/ConfigCase.pm

================================================================================
Originator: sacks
Date: 14 Sep 2012
Model: scripts
Version: scripts4_120914
One-line: more flexible match of CCSM_LCOMPSET in set_pes

In determining the PE layout, allow matching CCSM_LCOMPSET if the
attribute value matches any part of the compset name, rather than
requiring a match at the beginning, as is done for most attributes
(similar to GEN_COMPSET_MATCH in config_compsets.xml).

For example, CCSM_LCOMPSET="_10yLONG|_TEST" matches any compsets that
have either "_10yLONG" or "_TEST" anywhere in their long compset name.

(also, more doc changes from Mariana, not documented here)

M       ccsm_utils/Case.template/ConfigCase.pm

================================================================================

Originator: erik
Date: 13 Sep 2012
Model: scripts
Version: scripts4_120913
One-line: Settings for single-point

M       ccsm_utils/Case.template/config_compsets.xml - add more settings for 
              single-point urban resolutions
M       ccsm_utils/Case.template/ConfigCase.pm ------- use exists in pesetup over defined
M       ccsm_utils/Testcases/ERS211d_script ---------- fix set syntax

M       create_newcase -- Set MPILIB=mpi-serial when 1 task set

================================================================================

Originator: mvertens
Date: 12 Sep 2012
Model: scripts
Version: scripts4_120912
One-line: new userdefined machine and Macros now called from setup 

Depends on Machines_120912	
	
M       create_newcase
	- no longer calls generation of Macros
	- no longer need din_loc_root, scratchroot and max_tasks_per_node 
	- generic machine now replaced by userdefined machine name
M       create_test
	- no longer need din_loc_root, scratchroot and max_tasks_per_node 
	- generic machine now replaced by userdefined machine name
M       create_test_suite
	- no longer need din_loc_root, scratchroot and max_tasks_per_node 
	- generic machine now replaced by userdefined machine name
M       ccsm_utils/Tools/preview_namelists
	- fixed bug for seeing multi-instance namelists in CaseDocs
M       ccsm_utils/Tools/setup
	- calls Macros
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/ConfigCase.pm
	- new grouping for auto-documentaiton
M       ccsm_utils/Case.template/config_definition.xml
        - extensive documentation added in ldesc for auto-generated documentation
	
================================================================================
Originator: jshollen
Date: 11 Sep 2012
Model: scripts
Version: scripts4_120911a
One-line: Modified B compsets to 6 days, F compsets to 2 days, E compsets to 6 days, 
I compsets to 2 days, added SMS2,SMS6 test scripts 

A       ccsm_utils/Testcases/SMS2_script
A       ccsm_utils/Testcases/SMS6_script
M       ccsm_utils/Testlists/B41
M       ccsm_utils/Testlists/A01
M       ccsm_utils/Testlists/A02
M       ccsm_utils/Testlists/A03
M       ccsm_utils/Testlists/C01
M       ccsm_utils/Testlists/A04
M       ccsm_utils/Testlists/C02
M       ccsm_utils/Testlists/A05
M       ccsm_utils/Testlists/C03
M       ccsm_utils/Testlists/C04
M       ccsm_utils/Testlists/C41
M       ccsm_utils/Testlists/C05
M       ccsm_utils/Testlists/C42
M       ccsm_utils/Testlists/C43
M       ccsm_utils/Testlists/C44
M       ccsm_utils/Testlists/C46

================================================================================
Originator: fvitt
Date: 11 Sep 2012
Model: scripts
Version: scripts4_120911
One-line: added STRATMAM7 compsets and changed STRATMAM compsets to STRATMAM3 and 
adjustment to MOZSOA compsets
	
M       ccsm_utils/Case.template/config_compsets.xml
 - change *MOZSOA* compsets to use new cam build-namelist use case

 - new compsets:
   B_2000_STRATMAM7_CN
   F_2000_STRATMAM7_CN
  - rename compsets:
   B_2000_STRATMAM_CN -> B_2000_STRATMAM3_CN
   F_2000_STRATMAM_CN -> F_2000_STRATMAM3_CN

M       ccsm_utils/Testlists/C04
M       ccsm_utils/Testlists/C05
 - change test lists for the compset name change "FSTRATMAM" to "FSTRATMAM3"

================================================================================
Originator: erik
Date: 10 Sept 2012
Model: scripts
Version: scripts4_120910
One-line: Convert create_test_suite to perl and extend it's functionality, add I1PT compset, tweak testlists
  - add more CLM tests in especially for single-point/regional
  - add -user_nl_dir option to create_test
  - add some more lengths for restart tests
  - add ability for create_test_suite to process tests with -user_nl_dir options
  - change create_test_suite to perl
  - allow new format for xmlchange for lists of options without needing filename
  - add new xmlquery function to query env_*.xml variables or list for all.
  - add I1PT compset
  - turn mpi-serial on by default for CLM regional grids (requires Machines tag)

>>>>>>>>> Script to query XML variables. Gets you the value (expanded or not), and
>>>>>>>>> filename. "use xmlquery list" to lists ALL variables in ALL files.
A  +    ccsm_utils/Tools/xmlquery -- perl script to query env_*.xml files.

>>>>>>>>> Add new lengths of restart tests and a build-namelist test
A  +    ccsm_utils/Testcases/ERS3y_script --- 3 years
A  +    ccsm_utils/Testcases/ERS211d_script - 211 days
A  +    ccsm_utils/Testcases/ERS48s_script -- 48 steps
A  +    ccsm_utils/Testcases/SBN_script ----- Smoke build-namelist test
          (This tests that a case can be setup and input files are available)
          (Good for fast checking that the XML database is correct)
          (normally use with the -nobatch and -nobuild options to create_test_suite)

>>>>>>>>> Add new test lists
A  +    ccsm_utils/Testlists/lynx.clm.interactive.auxtest
A  +    ccsm_utils/Testlists/allIcompsetsRes.clm.auxtest
A  +    ccsm_utils/Testlists/bluefire.clm.interactive.auxtest
A  +    ccsm_utils/Testlists/shortlist.clm.auxtest

>>>>>>>>> Some changes to get create_test_suite working, new UI to xmlchange
M       create_test ------ Add new options
	      -user_nl_dir 
              -create_suite option and options for it: -reruntests, -nobatch, -nobuild
M       create_newcase --- Correct some documentation, copy xmlquery, may set MPILIB.
M       ccsm_utils/Tools/check_input_data - Return error status if file(s) are missing
M       ccsm_utils/Tools/testcase_end ----- Add COMPARE_NAMELISTS option
M       ccsm_utils/Tools/xmlchange -------- Add new input option where file isn't required
           and you enter a list of variables, and values as:
             xmlchange var1=value1,var2=value2,var3=value3


>>>>>>>>> New I1PT compset, SUPPORTED_BY and new testcases.
M       ccsm_utils/Case.template/config_compsets.xml ---- Add I1PT compset
M       ccsm_utils/Case.template/config_compsets.xsl ---- Add WRF compsets to listings
MM      ccsm_utils/Case.template/config_definition.xml -- Add SUPPORTED_BY
M       ccsm_utils/Testcases/config_tests.xml ----------- Add new testcases:
             SBN, ERS48s, ERS211d, ERS3y

>>>>>>>>> Changes to test lists
M       ccsm_utils/Testlists/bluefire.clm.trans.auxtest -- Add comments
M       ccsm_utils/Testlists/lynx.clm.auxtest ------------ Add more tests/compilers
M       ccsm_utils/Testlists/bluefire.clm.auxtest -------- Add comments
M       ccsm_utils/Testlists/B01 --- add SMS_D.5amazon.ITEST
M       ccsm_utils/Testlists/C01 --- Add single-point spinupCN cases
M       ccsm_utils/Testlists/C02 --- Add single-point spinupCN cases
M       ccsm_utils/Testlists/C03 --- Add single-point spinupCN cases

M       create_test_suite ----- COMPLETELY CHANGED! And rewritten in Perl!
            New options:
     -compset_file <file>    Specifies the compset file to be used.
     -debug                  If you want to run create_test_suite to show what would happen
                             when it runs without actually running tests.
                             (also turns verbose on and autosubmit to off) [Mostly useful for debugging testing]
     -nobatch [on | off]     Run tests interactively rather than submit to queue in cs.submit.*.
     -nobuild [on | off]     Don't run the *.build in cs.submit.*.
     -reruntests [on | off]  Rerun tests when the cs.submit.* script is run again even if they 
                             previously PASSed. (you can also change this in the cs.submit.* file created)

    (-nobatch on -nobuild on are especially useful with lists of SBN tests that just verify
     you can setup a case, the namelist will be built, and you have the files for it)

     Can now also read test lists with: 
        comments (ignores everything after # on a line)
        options to create_test

   Environment variables were removed from create_test_suite, so control of create_test was
   all done by explicit command-line arguments.

Testing:

    bluefire: prebeta tests done All PASS except...
SFAIL SMS.1x1_mexicocityMEX.I1PT.bluefire_ibm.110521  <-- new compset
FAIL  SMS.1x1_numaIA.ICN.bluefire_ibm  <----------------- CLM issue, new test
FAIL  CME.f19_g16.S.bluefire_ibm <----------------------- existing issue
CHECK ERT_PT.T62_g16.GIAF.bluefire_ibm.perf npes=128 tput=15.108 memh=931.674 memr=-0.001 tag=cesm1_1_alpha16e
CHECK ERT.T62_g16.CIAF.bluefire_ibm.perf npes=128 tput=19.774 memh=594.756 memr=-0.001 tag=cesm1_1_alpha16e
IOP: RUN   ERS_IOP.f09_g16.B20TRBPRP.bluefire_ibm.110521 
IOP: RUN   ERS_IOP.f19_f19.FWX.bluefire_ibm.110521 
IOP: RUN   SMS_IOPACGI.ne30_f19_g16.A.bluefire_ibm.110521
    lynx: prebeta tests done All PASS except...
IOP: "FAIL <----------------------------------------- existing issue
IOP: "FAIL
IOP: "FAIL
IOP: "RUN
IOP: "RUN
IOP: "RUN
IOP: "RUN
IOP: "RUN
IOP: "RUN
IOP: "RUN
RUN   ERS_IOP.f19_f19.FWX.lynx_pgi.233949
SFAIL SMS.1x1_mexicocityMEX.I1PT.lynx_pgi.233949  <-- new compset
FAIL  SMS.1x1_numaIA.ICN.lynx_pgi  <----------------- CLM issue, new test
FAIL  CME.f19_g16.S.lynx_pgi  <---------------------- existing issue
CHECK ERT_PT.T62_g16.GIAF.lynx_pgi.perf npes=128 tput=10.654 memh=0.000 memr=0.000 tag=cesm1_1_alpha16e
CHECK ERT.T62_g16.CIAF.lynx_pgi.perf npes=64 tput=8.243 memh=0.000 memr=0.000 tag=cesm1_1_alpha16e


================================================================================
Originator: sacks
Date: 06 Sept 2012
Model: scripts
Version: scripts4_120906
One-line: change some NCK tests from B to BG; add comments in scripts

(also, more doc changes from Mariana, not documented here)

M       ccsm_utils/Testlists/C01
M       ccsm_utils/Testlists/C02
M       ccsm_utils/Testlists/C04
M       ccsm_utils/Testlists/C05
        - change NCK B1850CN to BG1850CN

M       ccsm_utils/Tools/component_compare.sh
M       ccsm_utils/Tools/component_gen_comp
        - add comments

================================================================================
Originator: jedwards
Date: 05 Sept 2012
Model: scripts
Version: scripts4_120905
One-line: Adresses issues in resolving CASEROOT on intrepid, various
changes in doc files are being commited without tags

	
	M ccsm_utils/Tools/setup
	M ccsm_utils/Tools/taskmaker.pl


================================================================================
Originator: mlevy
Date: 29 Aug 2012
Model: scripts
Version: scripts4_120829
One-line: Another fix to documentation generation script because last tag had 
          issues toggling the help on variables appearing in multiple
          namelists.

Note:     There were changes committed to the trunk between tags 120828 and
          120828b that were not tagged, and I didn't make them so I don't
          know what changed... the same thing happened between tag 120828b and
          this tag.

M       doc/nldef2html_pop2

================================================================================
Originator: mlevy
Date: 28 Aug 2012
Model: scripts
Version: scripts4_120828b
One-line: updated documentation generation script to handle pop2 variables with
          the same name appearing in multiple namelists.

M       doc/nldef2html_pop2

================================================================================
Originator: sacks
Date: 28 Aug 2012
Model: scripts
Version: scripts4_120828
One-line: add multi-instance TG tests; add comment to component_compare

M       ccsm_utils/Testlists/bluefire.glc.auxtest
        - add multi-instance tests (add one test, change another)
M       ccsm_utils/Tools/component_compare.sh
        - add comment

================================================================================
Originator: sacks
Date: 27 Aug 2012
Model: scripts
Version: scripts4_120827d
One-line: fixed NCK_script bug for real

M       ccsm_utils/Testcases/NCK_script

================================================================================
Originator: mvertens
Date: 27 Aug 2012
Model: scripts
Version: scripts4_120827c
One-line: fixed NCK_script bug

M       ccsm_utils/Testcases/NCK_script
	
================================================================================
Originator: sacks
Date: 27 Aug 2012
Model: scripts
Version: scripts4_120827b
One-line: modify CME10y_script to stay consistent with CME_script

Removed two instances of './setup -clean; ./setup' that were intentionally
removed from CME_script, but accidentally were not removed from the 10-year
version

M       ccsm_utils/Testcases/CME10y_script

================================================================================
Originator: jshollen
Date: 27 Aug 2012
Model: scripts
Version: scripts4_120827
One-line: modified testreporter to report on IOP tests separately.  

M       ccsm_utils/Tools/testreporter.pl

================================================================================
Originator: sacks
Date: 26 Aug 2012
Model: scripts
Version: scripts4_120826
One-line: add/modify _GLC compsets, modify GLC auxtest lists, add scripts to facilitate testing

This changes answers for BGCN, BG1850CN and FGCN: the change to
config_compsets.xml fixes bugs in these compsets

A       ccsm_utils/Tools/component_gen_comp
A       ccsm_utils/Tools/component_generate.sh
A       ccsm_utils/Tools/component_compare.sh
        - New scripts for doing history file generation & comparison on
          component history files (in contrast to cpl hist files) after
          running CESM tests
        - The main script is component_gen_comp; it uses the other two

M       ccsm_utils/Case.template/config_compsets.xml
        - Add compsets for testing _GLC: IGLONG (10-year IG) and TGG10
          (TG with gland10)
        - Fix BGCN, BG1850CN and FGCN compsets: add CLM_NML_USE_CASE
          - For BGCN & BG1850CN, this meant fixing a typo in the xml; for
            FGCN, had to add missing clause
          - The typo fix also meant that BGCN & BG1850CN now have
            CLM_FORCE_COLDSTART="on", and BG1850CN has a CAM_NML_USE_CASE
        - Refactored definition of TG compsets (adding
          BEG_COMPSET_MATCH clauses, and explicit setting of RUN_STARTDATE)

R  +    ccsm_utils/Testlists/bluefire.glc.auxtest
        - Add tests of TGG10 and IGLONG compsets for more robust testing
        - Change SMS_D.f19_g16.IG.bluefire_ibm to T31_g37 to reduce
          queue wait time
        - Delete one F and one B compset test to reduce test turnaround time
        - Change remaining B test from ERI to SMS to reduce test time
        - Change an I test from ERS to ERI so we still have one non-TG ERI test

A  +    ccsm_utils/Testlists/shortlist.glc.auxtest
        - New lightweight GLC test suite (subset of bluefire.glc1.auxtest)

================================================================================

Originator: mvertens
Date: 25 Aug 2012
Model: scripts
Version: scripts4_120825b
One-line: bugfix - resolved $CASEROOT in $CASE.run

M      ccsm_utils/Tools/setup
	
================================================================================
Originator: mvertens
Date: 25 Aug 2012
Model: scripts
Version: scripts4_120825
One-line: changed executable permissions on setup script

M      ccsm_utils/Tools/setup
       - also added bug fix for calling user_nl_cism.csh
	
================================================================================
Originator: mvertens
Date: 24 Aug 2012
Model: scripts
Version: scripts4_120824
One-line: new streamlined case setup and updates to scripts infrastructure  

1) $case.build, $case.clean_build, $case.submit, $case.l_archive created by create_newcase
   - the above scripts have almost everything inline now - the only csh files these
     scripts now call are
	$CASEROOT/Buildconf/$comp.buildnml.csh           (by $CASE.build and $CASE.run)
        $CASEROOT/Buildconf/$comp.buildexe.csh           (by $CASE.build)	
        $CCSMROOT/scripts/ccsm_utils/Tools/getTiming.csh (by $CASE.run)
        $CCSMROOT/scripts/ccsm_utils/Tools/hist_compare.csh (by $CASE.test)	
2) $case.run, user_nl_xxx create by setup (replacement to configure)
   - the tests below were changed to replace calls to configure with setup and setup -clean 
3) env_build.xml now has all grid information
4) ccsm script names changed to cesm script names	
	
Deleted Files: 	
D       ccsm_utils/Tools/configure
D       ccsm_utils/Tools/ccsm_buildexe.csh
D       ccsm_utils/Tools/ccsm_postrun.csh
D       ccsm_utils/Tools/listfilesin_streams
D       ccsm_utils/Tools/ccsm_buildnml.csh
D       ccsm_utils/Tools/timing/getTiming.pl
D       ccsm_utils/Tools/timing/README.getTiming
D       ccsm_utils/Tools/ccsm_prestage.csh
D       ccsm_utils/Tools/build_streams
D       ccsm_utils/Tools/generate_batch.csh
	
Added Fiels:
A       ccsm_utils/Tools/cesm_prerun_setup
A       ccsm_utils/Tools/cesm_submit
A       ccsm_utils/Tools/cesm_buildnml
A       ccsm_utils/Tools/cesm_prestage
A       ccsm_utils/Tools/cesm_buildstart
A       ccsm_utils/Tools/setup
A       ccsm_utils/Tools/cesm_postrun_setup
A       ccsm_utils/Tools/cesm_buildexe
A       ccsm_utils/Tools/cesm_clean_build
	
Modified Files	
M       create_clone
M       create_newcase
M       ccsm_utils/Tools/check_input_data
M       ccsm_utils/Tools/ccsm_check_lockedfiles
M       ccsm_utils/Tools/archive_metadata.sh
M       ccsm_utils/Tools/preview_namelists
M       ccsm_utils/Tools/clean_build
M       ccsm_utils/Tools/xml2env
M       ccsm_utils/Tools/concat_daily_hist.csh
M       ccsm_utils/Tools/taskmaker.pl
M       ccsm_utils/Tools/ccsm_getenv
M       ccsm_utils/Tools/xmlchange
M       ccsm_utils/Tools/testcase_setup.csh
M       ccsm_utils/Tools/create_production_test
M       ccsm_utils/Case.template/ConfigCase.pm
M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Testcases/PMT_script
M       ccsm_utils/Testcases/PEA_script
M       ccsm_utils/Testcases/PST_script
M       ccsm_utils/Testcases/ERB_script
M       ccsm_utils/Testcases/NCK_script
M       ccsm_utils/Testcases/ERH_script
M       ccsm_utils/Testcases/PEM_script
M       ccsm_utils/Testcases/APT_script
M       ccsm_utils/Testcases/CME10y_script
M       ccsm_utils/Testcases/CME_script
M       ccsm_utils/Testcases/LAR_script
M       ccsm_utils/Testcases/ERI_script
M       ccsm_utils/Testcases/ERI44y_script
M       ccsm_utils/Testcases/PET_script
M       ccsm_utils/Testcases/CME4_script
M       ccsm_utils/Testcases/SEQ_script
	
================================================================================
Originator: sacks
Date: 20 Aug 2012
Model: scripts
Version: scripts4_120820
One-line: add tests for TG compsets, generalize other tests to work with TG

Rationale: Since a reasonable test of TG needs to run for a few years
(which only takes a few minutes on a single processor), I needed to
generalize many of the CESM tests. 

For many tests, this could be done by generalizing the test script to
simply use the compset's default run length. However, for the restart
tests and a few others (e.g., the CME test), the tests set run lengths
that differ from the compset's default. Thus, for some of these tests,
I have added new versions that run for multiple years.

For the CME tests, I have also reworked the existing tests so that the
different CME tests are identical except for a few lines near the top
that set variables. I have not done this for the ERS & ERI tests,
because this will require some more thought about how these tests
should be generalized.

A       ccsm_utils/Tools/year_string.sh
        - helper utility for some of the test scripts (convert integer
          into yyyy string)

A       ccsm_utils/Testcases/CME10y_script
A       ccsm_utils/Testcases/ERI44y_script
A       ccsm_utils/Testcases/ERS20y_script
        - add multi-year versions of tests that set their own run
          length (10-year CME, 44-year ERI, 20-year ERS)

M       ccsm_utils/Testcases/CME_script
M       ccsm_utils/Testcases/CME4_script
        - pulled variables to the top so that these are nearly
          identical to each other and to the new CME10y_script

M       ccsm_utils/Testcases/PMT_script
M       ccsm_utils/Testcases/SMS_script
M       ccsm_utils/Testcases/PEA_script
M       ccsm_utils/Testcases/PST_script
M       ccsm_utils/Testcases/NCK_script
M       ccsm_utils/Testcases/PEM_script
M       ccsm_utils/Testcases/APT_script
M       ccsm_utils/Testcases/PET_script
        - modified to remove hard-coded 5-day run length: now take run
          length from the compset's default, rather than explicitly
          setting STOP_OPTION/STOP_N (for all compsets other than TG,
          this gives the same behavior as before, because all other
          compsets default to a 5-day run already)

M       ccsm_utils/Testcases/config_tests.xml
        - modified test descriptions to be consistent with the changed
          tests: run length now given as "default length" rather than
          "5 day"; STOP_OPTION/STOP_N removed for these tests
        - added descriptions for CME10y, ERI44y, ERS20y

M       ccsm_utils/Testlists/A01
        - added an ERS20y TG test

M       ccsm_utils/Testlists/bluefire.glc.auxtest
        - add TG tests
        - replace non-working tests with working versions: switched
          some resolutions so that all tests have datasets
        - add _ibm in the test names

================================================================================
Originator: sacks
Date: 14 Aug 2012
Model: scripts
Version: scripts4_120814
One-line: add new _GLC compsets, change some default settings for _GLC compsets

- turn on annual history files for all _GLC compsets (this is done by
  setting HIST_OPTION & HIST_N, because cism's history output is
  currently linked to the cpl history output)
- set default run length to 10 years for all T_ compsets
- change forcing data location for TG compset
- add TG1850, TG20TR, TGRCP85 compsets
- add IG4804 & IG4804CN compsets
- add tests for new IG compsets

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Testlists/titan.clm.auxtest
M       ccsm_utils/Testlists/lynx.clm.auxtest
M       ccsm_utils/Testlists/bluefire.clm.auxtest
================================================================================
Originator: tcraig
Date: 10 Aug 2012
Model: scripts
Version: scripts4_120810
One-line: update ICP test

M       ccsm_utils/Testcases/ICP_script
================================================================================
Originator: muszala
Date: 8 Aug 2012
Model: scripts
Version: scripts4_120808
One-line: add support for r01 rtm 

M       scripts/ccsm_utils/Case.template/config_grid.xml
M       scripts/ccsm_utils/Case.template/config_definition.xml
================================================================================
Originator: tcraig
Date: 7 Aug 2012
Model: scripts
Version: scripts4_120807
One-line: update ICP test

M       ccsm_utils/Testcases/config_tests.xml
M       ccsm_utils/Testcases/ICP_script
================================================================================

Originator: jedwards
Date: 1 Aug 2012
Model: scripts
Version: scripts4_120801a
One-line: fix PEA test and intrepid test lists

M            ccsm_utils/Testcases/PEA_script
M            ccsm_utils/Testlists/C05


================================================================================

Originator: sacks
Date: 1 Aug 2012
Model: scripts
Version: scripts4_120801
One-line: fix short-term archiver for multi-instance

M       ccsm_utils/Tools/st_archive.sh

================================================================================

Originator: jshollen
Date: 26 July 2012
Model: scripts
Version: scripts4_120726
One-line: Updated testreporter to use TestStatus.IOP for _IOP test results. 

M       ccsm_utils/Tools/testreporter.pl

================================================================================

Originator: jshollen
Date: 19 July 2012
Model: scripts
Version: scripts4_120719b
One-line: Updated testreporter script to report compiler. 

M       ccsm_utils/Tools/testreporter.pl

================================================================================

Originator: jedwards
Date: 19 July 2012
Model: scripts
Version: scripts4_120719a
One-line: Remove cpl from build script since it is a no-op anyway

M     ccsm_utils/Tools/ccsm_buildexe.csh

================================================================================
Originator: jshollen
Date: 19 July 2012
Model: scripts
Version: scripts4_120719
One-line: added test lists for lynx_intel, janus_intel prealpha, prebeta

A       ccsm_utils/Testlists/lynx_intel.prealpha
A       ccsm_utils/Testlists/janus_intel.prebeta
A       ccsm_utils/Testlists/janus_intel.prealpha
A       ccsm_utils/Testlists/lynx_intel.prebeta

================================================================================
Originator: mlevy
Date: 18 July 2012
Model: scripts
Version: scripts4_120718a
One-line: Include disclaimer than tx1v1 is not for science in description of 
          f19_s11 and T62_s11 resolutions.

M            38773   ccsm_utils/Case.template/config_grid.xml

================================================================================	
Originator: jedwards
Date: 18 July 2012
Model: scripts
Version: scripts4_120718
One-line: fix issues in T85 and T340 grid setups

M            38768   ccsm_utils/Tools/testcase_begin
M            38768   ccsm_utils/Case.template/config_grid.xml

================================================================================	
Originator: tcraig
Date: 17 July 2012
Model: scripts
Version: scripts4_120717
One-line: update for new cice version

- add CICE_DECOMPSETTING, COMP_RUN_BARRIERS
- remove ICE_DOMAIN_PRIMEFACS
- add ne240_t12 grid 
- add ICP testcase

M       ccsm_utils/Tools/configure
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Testcases/config_tests.xml
A  +    ccsm_utils/Testcases/ICP_script
	
================================================================================	
Originator: jshollen
Date: 12 July 2012
Model: scripts
Version: scripts4_120712a
One-line: bug fix in testreporter script. 

M       ccsm_utils/Tools/testreporter.pl

================================================================================	
Originator: jshollen
Date: 12 July 2012
Model: scripts
Version: scripts4_120712
One-line: Remove problematic tests from lynx,intrepid prebeta, add testreporter script. 

A       ccsm_utils/Tools/testreporter.pl
M       ccsm_utils/Testlists/intrepid.prebeta
A       ccsm_utils/Testlists/C05
A       ccsm_utils/Testlists/C46
M       ccsm_utils/Testlists/lynx_pgi.prebeta
M       create_test_suite

================================================================================	
Originator: jedwards
Date: 09 July 2012
Model: scripts
Version: scripts4_120709
One-line:  Refactor IOP test, build gptl before PIO 

	 create_newcase
	 create_test
	 ccsm_utils/Tools/testcase_begin
	 ccsm_utils/Tools/testcase_env.csh
	 ccsm_utils/Tools/testcase_end
             Refactored IOP test to be a repeat of the original test using pnetcdf

	
	 ccsm_utils/Tools/configure
	 ccsm_utils/Tools/ccsm_buildexe.csh
	     Build gptl prior to PIO
	
	 ccsm_utils/Tools/ccsm_postrun.csh
	     Allow for COMP_RUN_BARRIERS variable to not exist
	
	 ccsm_utils/Testcases/ERS6_script
	     Correct messages printed to stdout
	
	 ccsm_utils/Testlists/B01
             Modify tests to include IOP
	
================================================================================	
Originator: jedwards
Date: 02 July 2012
Model: scripts
Version: scripts4_120702
One-line: Corrected parsing of pecount variable
         M ccsm_utils/Case.template/ConfigCase.pm

================================================================================	
Originator: jedwards
Date: 29 June 2012
Model: scripts
Version: scripts4_120629a
One-line: moved single proc test from intrepid, refactored config_pes.xml
parsing
	M create_newcase
	M ccsm_utils/Case.template/ConfigCase.pm
	M ccsm_utils/Testlists/C02
	M ccsm_utils/Testlists/C04

================================================================================
Originator: jshollen
Date: 29 June 2012
Model: scripts
Version: scripts4_120629
One-line: Changed bluefire prealpha ERS2 B compset tests to ERS6, POP will not restart 
		  on day 2.  

M       ccsm_utils/Testlists/A01

================================================================================

Originator: mvertens
Date: 26 June 2012
Model: scripts
Version: scripts4_120626
One-line: new auto-documentaton of namelists, executable now goes into 
	  $EXEROOT/bin/ccsm.exe, rather than $EXEROOT - so $RUNDIR and $EXEROOT
	  are now independent, removed DIN_LOC_ROOT_CSMDATA, prestaging now occurs
	  via linking to $DIN_LOC_ROOT/ccsm4_init rather than on copying, 
	  removed $SHAREROOT, added several new tests that need to be added to
	  auto-test suite and that are shorter

M       SVN_EXTERNAL_DIRECTORIES
M       ccsm_utils/Tools/ccsm_buildexe.csh
M       ccsm_utils/Tools/clean_build
M       ccsm_utils/Tools/ccsm_buildnml.csh
M       ccsm_utils/Tools/ccsm_prestage.csh
M       ccsm_utils/Tools/taskmaker.pl
X       ccsm_utils/Tools/lnd/clm/PTCLM
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/ConfigCase.pm
M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Testcases/config_tests.xml
M       ccsm_utils/Testcases/ERB_script
A       ccsm_utils/Testcases/ERS2_script
M       ccsm_utils/Testcases/ERH_script
M       ccsm_utils/Testcases/ERI_script
A       ccsm_utils/Testcases/ERS6_script
A       ccsm_utils/Testcases/CME4_script
D       ccsm_utils/Testcases/ER3_script
M       ccsm_utils/Testlists/A01
A       doc/nldef2html_cice
A       doc/nldef2html_clm
A       doc/xmldef2html_grid
A       doc/nldef2html_dlnd
D       doc/env_build_list.xml
D       doc/grids_list.xml
A       doc/nldef2html_datm
A       doc/images
A       doc/images/arrow_down.gif
A       doc/images/arrow_right.gif
A       doc/nldef2html_cism
A       doc/nldef2html_drv
A       doc/nldef2html_cam
A       doc/nldef2html_dice
A       doc/xmldef2html_caseroot
D       doc/env_conf_list.xml
A       doc/create_tables
A       doc/xmldef2html_machines
A       doc/nldef2html_pop2
A       doc/nldef2html_docn
D       doc/env_run_list.xml
A       doc/showinfo.js
D       doc/env_case_list.xml
M       create_newcase
	
================================================================================
Originator: fvitt
Date: 25 June 2012
Model: scripts
Version: scripts4_120625
One-line: Correction to *_STRATMAM compsets CAM namelist use case
	  and correction to msmkdir message

M       ccsm_utils/Tools/lt_archive.sh
M       ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: mlevy
Date: June 8, 2012
Model: scripts
Version: scripts4_120608a
One-line: Added env variable $PREVIEW_NML to preview_namelists;
          pop2.buildnml.csh will see this variable and call build-namelist with
          the "-preview" flag, printing out changes made in user_nl_pop2.

M       ccsm_utils/Tools/preview_namelists
================================================================================
Originator: mlevy
Date: June 8, 2012
Model: scripts
Version: scripts4_120608a
One-line: Forgot to update Changelog when I made last tag (here and in perl5lib)

 M      .
M       SVN_EXTERNAL_DIRECTORIES
M       ChangeLog

================================================================================
Originator: mlevy
Date: June 8, 2012
Model: scripts
Version: scripts4_120608
One-line: Updated perl5lib

 M      .
M       SVN_EXTERNAL_DIRECTORIES

================================================================================
Originator: mvertens
Date: June 4, 2012
Model: scripts
Version: scripts4_120604
One-line: grid sizes introduced back into create_newcase - needed for X compsets
	
M       ccsm_utils/Case.template/config_definition.xml
	
================================================================================
Originator: mvertens
Date: June 3, 2012
Model: scripts
Version: scripts4_120603
One-line: grid sizes introduced back into config_definition.xml for X compsets
	
M       ccsm_utils/Case.template/config_definition.xml
	
================================================================================
Originator: mvetens
Date: May 31, 2012
Model: scripts
Version: scripts4_120531
One-line: Bug fix for domain.lnd.ne240np4_gx1v6

M       ccsm_utils/Case.template/config_grid.xml

================================================================================
Originator: jshollen
Date: May 30, 2012
Model: scripts
Version: scripts4_120530
One-line: Add threaded CAM-SE test to bluefire's prebeta suite. 

M       ccsm_utils/Testlists/C45

================================================================================
Originator: mvertens
Date: May 20, 2018
Model: scripts
Version: scripts4_120528
One-line: Simplification of configure, rename of env_mach_pes.xml and
	  removal of many variables from config_definition.xml
	
- configure now only has two commands -case and -clean 
- Buildconf/ no longer has resolved variables
- all config_definition.xml variables that are pure driver namelist variables
  and that are not shared between components have been moved out of config_grid.xml
  and config_definition.xml and into drv/bld/namelist_files/namelist_defaults_drv.xml

M       create_newcase
	- direct copy of xxx.buildnml.csh and xxx.buildexe.csh into Buildconf/
	  no longer need to invoke xxx.cpl7.template
	- remove direct reference to model components such as clm, cam, etc 
	
D       ccsm_utils/Tools/generate_resolved.csh
	- no longer needed
M       create_test
M       sample_grid_file.xml
M       ccsm_utils/Tools/configure
	- only has two options now -case and -clean  
M       ccsm_utils/Tools/archive_metadata.sh
	- moved reference from env_mach_pes.xml to env_configure.xml
	  to reflect that env_configure.xml is locked when configure -case
	  is invoked
M       ccsm_utils/Tools/preview_namelists
	- added error flag if configure has not been invoked first
M       ccsm_utils/Tools/xml2env
	- moved reference from env_mach_pes.xml to env_configure.xml 
M       ccsm_utils/Tools/ccsm_prestage.csh
	- removed reference to DIN_LOC_ROOT_CSMDATA
M       ccsm_utils/Tools/taskmaker.pl
M       ccsm_utils/Tools/ccsm_getenv
M       ccsm_utils/Tools/xmlchange
M       ccsm_utils/Tools/generate_batch.csh
M       ccsm_utils/Case.template/ConfigCase.pm
	- no longer reference model specific components (e.g. cam, clm) 
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/config_definition.xml
        - all config_definition.xml variables that are pure 
	  driver namelist variables and that are not shared 
	  between components have been moved out of config_grid.xml
          and config_definition.xml and into 
	  drv/bld/namelist_files/namelist_defaults_drv.xml
	
M       ccsm_utils/Testcases/PMT_script
M       ccsm_utils/Testcases/CME_script
M       ccsm_utils/Testcases/PEA_script
M       ccsm_utils/Testcases/PST_script
M       ccsm_utils/Testcases/ERB_script
M       ccsm_utils/Testcases/NCK_script
M       ccsm_utils/Testcases/LAR_script
M       ccsm_utils/Testcases/ERH_script
M       ccsm_utils/Testcases/ERI_script
M       ccsm_utils/Testcases/PEM_script
M       ccsm_utils/Testcases/APT_script
M       ccsm_utils/Testcases/PET_script
M       ccsm_utils/Testcases/SEQ_script
M       ccsm_utils/Testcases/ERS_script
	- changed env_mach_pes.xml to env_configure.xml
	- changed configure -cleanmach to -configure -clean
	
================================================================================
Originator: mvertens
Date: May 20, 2012
Model: scripts
Version: scripts4_120520
One-line:  Removal of env_conf.xml
	
remmoval of env_conf.xml is now doable since cam, clm and cice configure called from
cam.buildnml.csh, clm.buildnml.csh and cice.buildnml.csh. configure variables in 
env_conf.xml moved to env_build.xml, namelist variables in env_conf.xml moved to
	
M       create_clone
M       ccsm_utils/Tools/configure
M       ccsm_utils/Tools/ccsm_check_lockedfiles
M       ccsm_utils/Tools/archive_metadata.sh
M       ccsm_utils/Tools/xml2env
M       ccsm_utils/Tools/ccsm_prestage.csh
M       ccsm_utils/Tools/ccsm_getenv
M       ccsm_utils/Tools/xmlchange
M       ccsm_utils/Case.template/config_definition.xsl
M       ccsm_utils/Case.template/ConfigCase.pm
M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Testcases/ERB_script
M       ccsm_utils/Testcases/ERH_script
M       ccsm_utils/Testcases/LAR_script
M       ccsm_utils/Testcases/ERI_script
M       create_newcase
	
================================================================================
Originator: jshollen
Date: May 18, 2012
Model: scripts
Version: scripts4_120518
One-line:  Modify cs.submit script in create_test to recreate and/or rebuild tests in test suites. 

M       create_test

================================================================================
Originator: mvertens
Date: May 17, 2012
Model: scripts
Version: scripts4_120517
One-line: new perl5lib streams functionality

M       create_clone 
	- changes (jim edwards) put in for changes pointed
	  out by Chuck Bardeen
M       SVN_EXTERNAL_DIRECTORIES
        - updated to perl5lib_120517	
	
================================================================================
Originator: jedwards
Date: May 15, 2012
Model: scripts
Version: scripts4_120515a
One-line: Correct some files and paths for ne240 compsets

M         ccsm_utils/Case.template/config_grid.xml	

================================================================================
Originator: dfeddema 
Date: May 15, 2012
Model: scripts
Version: scripts4_120515
One-line: add memory leak and throughput tests 

A       Tools/check_memory.pl
A       Tools/compare_throughput.pl
M       Tools/testcase_end
M       Tools/testcase_setup.csh
M       Testcases/SMS_script
M       Testcases/ERI_script
M       Testcases/ERS_script
M       Testcases/ERT_script

================================================================================
Originator: fvitt
Date: May 14, 2012
Model: scripts
Version: scripts4_120514
One-line: remove default use of ESMF lib for *_SD compsets

M       ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: mvertens
Date: May 10, 2012
Model: scripts
Version: scripts4_120509b
One-line: add new support build-namelist additions for cism and cpl

M       SVN_EXTERNAL_DIRECTORIES
        updated to perl5lib/trunk_tags/perl5lib_120509
M       ccsm_utils/Tools/ccsm_buildnml.csh
M       ccsm_utils/Tools/ccsm_prestage.csh
A       ccsm_utils/Tools/user_nlcreate
	
================================================================================
Originator: tcraig
Date: May 09, 2012
Model: scripts
Version: scripts4_120509a
One-line: add support for wrf component

M       ccsm_utils/Tools/configure
M       ccsm_utils/Tools/st_archive.sh
M       ccsm_utils/Tools/ccsm_buildnml.csh
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/ConfigCase.pm
M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Testcases/config_tests.xml
A  +    ccsm_utils/Testcases/ER3_script
M       ccsm_utils/Testlists/bluefire.prerelease
A  +    ccsm_utils/Testlists/eastwind.rasm.auxtest
A  +    ccsm_utils/Testlists/garnet.rasm.auxtest
A  +    ccsm_utils/Testlists/olympus.rasm.auxtest
A  +    ccsm_utils/Testlists/evergreen.rasm.auxtest
A  +    ccsm_utils/Testlists/bluefire.rasm.auxtest
A  +    ccsm_utils/Testlists/chugach.rasm.auxtest
A  +    ccsm_utils/Testlists/raptor.rasm.auxtest
A       ccsm_utils/Testlists/W01
M       ccsm_utils/Testlists/bluefire.prebeta
M       create_newcase
================================================================================
Originator: jedwards
Date: May 09, 2012
Model: scripts
Version: scripts4_120509
One-line: change case4test to case in all tests, reduce compile time for NCK test.
 
M            36926   ccsm_utils/Testcases/PMT_script
M            36926   ccsm_utils/Testcases/CME_script
M            36926   ccsm_utils/Testcases/PEA_script
M            36926   ccsm_utils/Testcases/PST_script
M            36926   ccsm_utils/Testcases/ERB_script
M            36926   ccsm_utils/Testcases/NCK_script
M            36926   ccsm_utils/Testcases/LAR_script
M            36926   ccsm_utils/Testcases/ERH_script
M            36926   ccsm_utils/Testcases/ERI_script
M            36926   ccsm_utils/Testcases/PEM_script
M            36926   ccsm_utils/Testcases/APT_script

================================================================================
Originator: jshollen
Date: May 01, 2012
Model: scripts
Version: scripts4_120501
One-line: Modify create_newcase so that check_known_problems doesn't call svn. 

M      create_newcase
================================================================================
Originator: jedwards
Date: Apr 23, 2012
Model: scripts
Version: scripts4_120423
One-line: Remove CCSMBUILDONLY env variable

M       create_test
================================================================================

Originator: jshollen
Date: Apr 20, 2012
Model: scripts
Version: scripts4_120420
One-line: Commented out 'known_problems' functionality in create_newcase

M       create_newcase

================================================================================
Origiator: tcraig
Date: Apr 18, 2012
Model: scripts
Version: scripts4_120418
One-line: Va2o mapping, TG compset, getTiming2 tool
	
      - add support for Va2o mapping file
      - migrate source env_mach_specific into ccsm_getenv
      - update to getTiming2.pl for timing output
      - update TG compset for better out of the box setup
      - update S compset to use satm
	
M       ccsm_utils/Tools/ccsm_buildexe.csh
M       ccsm_utils/Tools/ccsm_postrun.csh
D       ccsm_utils/Tools/getTiming2.pl
D       ccsm_utils/Tools/getTiming2.csh
M       ccsm_utils/Tools/timing/getTiming.csh
A  +    ccsm_utils/Tools/timing/getTiming2.pl
M       ccsm_utils/Tools/ccsm_getenv
M       ccsm_utils/Tools/generate_batch.csh
X       ccsm_utils/Tools/lnd/clm/PTCLM
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/config_definition.xml
M       create_newcase
	
================================================================================
Origiator: jshollen
Date: Apr 16, 2012
Model: scripts
Version: scripts4_120416
One-line: Proof-of-concept: create_newcase checks the repository for known problems for
beta tags. 

M       create_newcase

================================================================================

Originator: mlevy
Date: Apr 13, 2012
Model: scripts
Version: scripts4_120413a
One-line:  moved ne16 and ne60 mapping files to follow new convention
	
M       ccsm_utils/Case.template/config_grid.xml
	
================================================================================

Originator: mlevy
Date: Apr 13, 2012
Model: scripts
Version: scripts4_120413
One-line:  added ne16_ne16, ne16_g37, ne60_ne60, and ne60_g16
	
M       ccsm_utils/Case.template/config_grid.xml
	
================================================================================
	
Originator: mvertens
Date: Apr 10, 2012
Model: scripts
Version: scripts4_120410
One-line: changed bug in PRIME_FACTORS for ne120_ne120 and ne240_ne240
	
M       ccsm_utils/Case.template/config_grid.xml
	
================================================================================
	
Originator: jedwards
Date: Apr 09, 2012
Model: scripts
Version: scripts4_120409
One-line: remove check_rundb.  modify some cam-se compsets.

M       create_newcase
D       check_rundb
M       create_test 
M       ccsm_utils/Case.template/config_grid.xml

================================================================================

Originator: mvertens
Date: Apr 05, 2012
Model: scripts
Version: scripts4_120405b
One-line: Script changes necessary to have DARWIN compsets to work

M       ccsm_utils/Tools/configure
M       ccsm_utils/Tools/ccsm_buildnml.csh
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_definition.xml

================================================================================
	
Originator: dfeddema 
Date: Apr 05, 2012
Model: scripts
Version: scripts4_120405
One-line: Modified PST, PEM and PET tests so that they generate cpl.hist files 


M       PST_script
M       PEM_script
M       PET_script

================================================================================
	
Originator: mvertens
Date: Apr 04, 2012
Model: scripts
Version: scripts4_120404
One-line: Changed DATM_CPL_ to DATM_CPLHIST_ to be consistent with new datm build-namelist

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_definition.xml
	
================================================================================
	
Originator: jedwards
Date: Apr 02, 2012
Model: scripts
Version: scripts4_120402
One-line: Corrected ne240 mapping file names, improved generic machine support, 
	removed debug print in create_test, added pio async support in taskmaker.pl
	
  M ccsm_utils/Case.template/ConfigCase.pm
  M ccsm_utils/Case.template/config_grid.xml
  M ccsm_utils/Tools/taskmaker.pl
  M create_test

================================================================================

Originator: mvertens
Date: Mar 29, 2012
Model: scripts
Version: scripts4_120329d
One-line: Fixed ERS test with removal of memory leak capablity

M       ccsm_utils/Testcases/ERS_script
	
================================================================================
	
Originator: mvertens
Date: Mar 29, 2012
Model: scripts
Version: scripts4_120329c
One-line: Added value of gland5UM go GLC_GRID
          (these should all be backwards compatible)

M       ccsm_utils/Case.template/config_definition.xml
	
================================================================================
	
Originator: mvertens
Date: Mar 29, 2012
Model: scripts
Version: scripts4_120329b
One-line: More changes needed to work with new data model build-namelists
          (these should all be backwards compatible)

M       ccsm_utils/Tools/ccsm_buildnml.csh
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_definition.xml
	
================================================================================
	
Originator: mvertens
Date: Mar 29, 2012
Model: scripts
Version: scripts4_120329
One-line: Removed call to check_memory_leak.py for now -and put in hooks for new
	  data model build-namelist functionality

M       SVN_EXTERNAL_DIRECTORIES
        - use of new tag that has generic stream capability (needed for data model
	  build-namelist)
	  https://svn-ccsm-models.cgd.ucar.edu/perl5lib/trunk_tags/perl5lib_120327	 
M       ccsm_utils/Tools/testcase_setup.csh
	- removal of call to check_memory_leak.py   
M       ccsm_utils/Case.template/config_definition.xml
	- new definitions needed for upcoming dlnd with build-namelist
M       ccsm_utils/Testcases/ERS_script
	- removal of call to check_memory_leak.py 
M       create_newcase
        - hooks to not copy data model xxx.cpl7.template files if they are 
	  not there
	
	
================================================================================
	
Originator: jshollen
Date: Mar 28, 2012
Model: scripts
Version: scripts4_120328
One-line: Fixing nonexistent compsets for bluefire, lynx prealpha tests.  

M       ccsm_utils/Testlists/A01
M       ccsm_utils/Testlists/A02

================================================================================

Originator: santos
Date: Mar 27, 2012
Model: scripts
Version: scripts4_120327b
One-line: Add WACCM-X compset (FWX), test pre-beta for bluefire/lynx.

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Testlists/B01
M       doc/compsets_list.xml
M       README

================================================================================
	
Originator: tcraig
Date: Mar 27, 2012
Model: scripts
Version: scripts4_120327a
One-line:  update preview namelist tool to leverage ccsm_buildnml

M       ccsm_utils/Tools/preview_namelists
================================================================================
	
Originator: erik
Date: Mar 27, 2012
Model: scripts
Version: scripts4_120327
One-line:  Fix T42_T42 and T85_T85 grids and set mpi-serial if only one task

Update perl5lib to allow empty strings.

M       ccsm_utils/Case.template/config_grid.xml - Fix T42_T42/T85_T85
M       create_newcase --- If mpilib NOT set, and only task --> set MPILIB to mpi-serial

Note that I think 4x5_T42_gx3v7 is messed up.

================================================================================
	
Originator: mvertens
Date: Mar 26, 2012
Model: scripts
Version: scripts4_120326b
One-line:  Put in new primefacs needed to run cice in prescribed mode on se grid

M       ccsm_utils/Case.template/config_grid.xml
	
================================================================================
Originator: jshollen
Date: Mar 26, 2012
Model: scripts
Version: scripts4_120326
One-line:  Removed ERS T31_T31 and ERS.ne30_g16.F, adding new ne30 tests, and adding MEGAN test for francis. 

M       ccsm_utils/Testlists/A01
M       ccsm_utils/Testlists/A02

================================================================================
Originator: dfeddema 
Date: Mar 25, 2012
Model: scripts
Version: scripts4_120325
One-line:  Modified tests to check for max memory per pe compared to previous baseline.   

M               ccsm_utils/Tools/check_memory_leak.py
M               ccsm_utils/Testcases/ERS_script

================================================================================
Originator: jedwards
Date: Mar 23, 2012
Model: scripts
Version: scripts4_120323
One-line: Corrected behavior of create_newcase when machine_compiler is used as mach flag

M      create_newcase	


================================================================================	
Originator: dfeddema
Date: Mar 16, 2012
Model: scripts
Version: scripts4_120316
One-line: Added create_newcase option for pnetcdf (IOP*). Added pnetcdf tests to
tests lists.  Added check for memory leaks to exact restart (ERS) tests.   

M      create_newcase
M      ccsm_utils/Testlists/bluefire.prebeta
A      ccsm_utils/Testlists/C45
       list of pnetcdf tests
A      ccsm_utils/Tools/check_memory_leak.py
       utility to check coupler logs for indication of memory leaks
M      ccsm_utils/Testcases/ERS_script
M      ccsm_utils/Tools/testcase_setup.csh


================================================================================

Originator: jedwards
Date: Mar 14, 2012
Model: scripts
Version: scripts4_120314
One-line: ifixed issue with aprun when tasks per node exceeds total tasks

M   ccsm_utils/Tools/taskmaker.pl

================================================================================

Originator: fvitt
Date: Mar 12, 2012
Model: scripts
Version: scripts4_120312
One-line: Changes for new MEGAN VOC emission fluxes, new compsets added,
  and turned on CLM/CN for FSTC

  New compsets:
    B_2000_STRATMAM_CN (BSTRATMAM) 
    F_2000_STRATMAM_CN (FSTRATMAM) 
    B_2000_MOZSOA_CN (BMOZSOA) 
    B_2000_MOZMAM_CN (BMOZMAM) 
    F_2000_MOZSOA_CN (FMOZSOA) 
    F_2000_MOZMAM_CN (FMOZMAM) 

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Testlists/C03
M       ccsm_utils/Testlists/C04
M       doc/env_conf_list.xml
	
================================================================================

Originator: jedwards
Date: Mar 09, 2012
Model: scripts
Version: scripts4_120309
One-line: Fix OS missing in config_compilers

M create_newcase
M ccsm_utils/Case.template/ConfigCase.pm 

================================================================================

Originator: erik
Date: Mar 08, 2012
Model: scripts
Version: scripts4_120308
One-line: Fix CLM test lists, move set_compiler to ConfigCase

M       ccsm_utils/Case.template/config_grid.xsl - Formatting
M       ccsm_utils/Case.template/config_grid.xml - Edit descriptions for CLM reg grids
M       ccsm_utils/Case.template/ConfigCase.pm --- Add set_compiler method from
            create_newcase (and parse_hash private method)

>>>>>>>>>> Fix names of regional grids so can run
M       ccsm_utils/Testlists/titan.clm.auxtest
M       ccsm_utils/Testlists/lynx.clm.auxtest
M       ccsm_utils/Testlists/bluefire.clm.auxtest
M       ccsm_utils/Testlists/edinburgh.clm.auxtest

M       create_newcase - Move set_compiler to ConfigCase.pm

================================================================================

Originator: jedwards
Date: Mar 07, 2012
Model: scripts
Version: scripts4_120307
One-line: new namelist variables for pio 

M ccsm_utils/Case.template/config_definition.xml

================================================================================

Originator: jedwards
Date: Mar 06, 2012
Model: scripts
Version: scripts4_120306
One-line: fix issue in create_test_suite

M create_test_suite

================================================================================

Originator: mvertens
Date: Mar 05, 2012
Model: scripts
Version: scripts4_120305
One-line: updated user_nl_add

M       ccsm_utils/Tools/user_nl_add
	
================================================================================
Originator: mvertens
Date: Mar 01, 2012
Model: scripts
Version: scripts4_120301
One-line: remove tabs from namelist generated by build-namelist and
          add routine user_nl_add to add only namelist lines to cesm_namelist

M       SVN_EXTERNAL_DIRECTORIES
	point to new perl5lib - perl5lib_120205
	that removes tabs from namelist generation 
        needed for upcoming POP2 build-namelist addition
A       ccsm_utils/Tools/user_nl_add
        utility to add only namelist lines to cesm_namelist
	
================================================================================
Originator: erik
Date: Feb 29, 2012
Model: scripts
Version: scripts4_120229
One-line: Remove namelist name and ending to user_nl_cam/clm files, 
          allow char config_def variables to have string length

Fix bug: 1427 -- will now abort if case-name length is too long

>>>>>>>>>>>> Remove namelist name and ending to user_nl_cam/clm files
M       create_newcase
>>>>>>>>>>>> Fix syntax errors, allow char to have string-length
>>>>>>>>>>>> will abort if exceeds string-length
M       ccsm_utils/Case.template/config_compsets.xml --- Fix syntax errors
M       ccsm_utils/Case.template/config_grid.xml ------- Fix syntax errors
M       ccsm_utils/Case.template/ConfigCase.pm --------- Allow char to have string
             length char*<len>
M       ccsm_utils/Case.template/config_definition.xml - Set string length for 
             CASE variables (CASE, RUN_REFCASE, DATM_CPL_CASE) and RUN_REFDATE
A       ccsm_utils/Case.template/config_grid.xsl ------- Stylesheet to view config_grid
             file

>>>>>>>>>>>> Add compiler to clm test lists
M       ccsm_utils/Testlists/titan.clm.auxtest
M       ccsm_utils/Testlists/intrepid.clm.auxtest
M       ccsm_utils/Testlists/bluefire.clm.trans.auxtest
A       ccsm_utils/Testlists/lynx.clm.auxtest
M       ccsm_utils/Testlists/bluefire.clm.auxtest
M       ccsm_utils/Testlists/edinburgh.clm.auxtest

================================================================================
Originator: fischer
Date: Feb 22, 2012
Model: scripts
Version: scripts4_120222
One-line: preview_namelist fix for cam5_1_21

M       ccsm_utils/Tools/preview_namelists

================================================================================
Originator: jedwards
Date: Feb 17, 2012
Model: scripts
Version: scripts4_120217
One-line: fix command line parsing in create_test
	M create_test
	
================================================================================
Originator: jedwards
Date: Feb 16, 2012
Model: scripts
Version: scripts4_120216
One-line: changes for compatibility with Machines_120216 which moves compiler options
	  to config_compilers.xml, fixes for test issues 

M            34797   ccsm_utils/Tools/testcase_env.csh
M            34797   ccsm_utils/Tools/generate_batch.csh
M            34797   ccsm_utils/Case.template/config_grid.xml
M            34797   create_test
M            34797   create_newcase
M            34797   create_test_suite



================================================================================	
Originator: mvertens
Date: Feb 13, 2012
Model: scripts
Version: scripts4_120213
One-line: creation of user_nl_comp and preview_namelists in caseroot 

A       sample_grid_file.xml
A       ccsm_utils/Tools/preview_namelists
M       create_newcase
	
================================================================================
Originator: jedwards
Date: Feb 16, 2012
Model: scripts
Version: scripts4_120216
One-line: 

- extend create_test and create_test_suite to handle machine and compiler
  setttings via testname and command line.
- add autosubmit flag to create_test_suite

M       create_test
M       create_test_suite
	
================================================================================
Originator: jedwards
Date: Feb 02, 2012
Model: scripts
Version: scripts4_120202
One-line: port to llnl sierra, update machines_list.xml

================================================================================
Originator: jedwards
Date: Feb 01, 2012
Model: scripts
Version: scripts4_120201
One-line: correct improperly formated fields in config_grid.xml, add more verbose prints to create_newcase
M ccsm_utils/Case.template/config_grid.xml
M create_newcase 

================================================================================        
Originator: fischer
Date: Jan 31, 2012
Model: scripts
Version: scripts4_120131
One-line:  add ESMF test to alpha test list 

M       ccsm_utils/Testlists/A01

================================================================================        
Originator: mvertens, jedwards
Date: Jan 25, 2012
Model: scripts
Version: scripts4_120125
One-line:  Restructure of Machines directory 



M            34081   create_test
M            34081   create_clone
M            34081   SVN_EXTERNAL_DIRECTORIES
D            34081   ccsm_utils/Build
D            34081   ccsm_utils/Build/mkSrcfiles
D            34081   ccsm_utils/Build/mkDepends
D            34081   ccsm_utils/Build/Makefile
X                    ccsm_utils/Tools/perl5lib
D            34081   ccsm_utils/Tools/ccsm_msread
MM           34081   ccsm_utils/Tools/configure
M            34081   ccsm_utils/Tools/ccsm_buildexe.csh
M            34081   ccsm_utils/Tools/ccsm_check_lockedfiles
A  +             -   ccsm_utils/Tools/lt_archive.sh
M            34081   ccsm_utils/Tools/archive_metadata.sh
D            34081   ccsm_utils/Tools/lt_archive.csh
D            34081   ccsm_utils/Tools/ccsm_mswrite
D            34081   ccsm_utils/Tools/ccsm_cpdata
M            34081   ccsm_utils/Tools/clean_build
M            34081   ccsm_utils/Tools/generate_resolved.csh
D            34081   ccsm_utils/Tools/ccsm_sedfile
D            34081   ccsm_utils/Tools/ccsm_msls
D            34081   ccsm_utils/Tools/ccsm_getfile
M           34081   ccsm_utils/Tools/timing/getTiming.pl
D            34081   ccsm_utils/Tools/ccsm_auto.csh
D            34081   ccsm_utils/Tools/lt_archive.pl
MM           34081   ccsm_utils/Tools/taskmaker.pl
M            34081   ccsm_utils/Tools/build_streams
MM           34081   ccsm_utils/Tools/ccsm_getenv
D            34081   ccsm_utils/Tools/ccsm_getinput
M            34081   ccsm_utils/Tools/generate_batch.csh
D            34081   ccsm_utils/Tools/ccsm_splitdf
D            34081   ccsm_utils/Tools/ccsm_msmkdir
X                    ccsm_utils/Tools/lnd/clm/PTCLM
D            34081   ccsm_utils/Tools/ccsm_l_archive.csh
MM           34081   ccsm_utils/Case.template/config_definition.xml
D            34081   ccsm_utils/Components
D            34081   ccsm_utils/Components/mct.buildlib
D            34081   ccsm_utils/Components/csm_share.buildlib
D            34081   ccsm_utils/Components/gptl.buildlib
D            34081   ccsm_utils/Components/pio.buildlib
M            34081   ccsm_utils/Testcases/config_tests.xml
M            34081   ccsm_utils/Testcases/PEA_script
M            34081   ccsm_utils/Testcases/ERS_script
D            34081   ccsm_utils/Testlists/hopp2.prebeta
D            34081   ccsm_utils/Testlists/hopp2.prerelease
A  +             -   ccsm_utils/Testlists/hopper.prebeta
A  +             -   ccsm_utils/Testlists/hopper.prerelease
M            34081   ccsm_utils/Testlists/C02
M            34081   ccsm_utils/Testlists/C03
A  +             -   ccsm_utils/Testlists/S01
M            34081   ccsm_utils/Testlists/bluefire.prebeta
M            34081   ChangeLog
M            34081   create_newcase
M            34081   .

================================================================================        
Originator: tcraig
Date: Jan 23, 2012
Model: scripts
Version: scripts4_120123
One-line:  updates for esmf, fix getTiming bug, extend create_test_suite
       - add E01 testlist for esmf and update bluefire.esmf.auxtest
       - modify csm_share build to include esmf dir whenever USE_ESMF_LIB is true
       - modify create_test_suite to extend parsing of testlist to allow
         both generic testlists (ie. A01) and custom lists (ie. ERS.*) mixed
       - fix getTiming bug introduced in scripts4_111221 that cause test summaries
	 to fail

M       ccsm_utils/Tools/timing/getTiming.csh
M       ccsm_utils/Components/csm_share.buildlib
M       ccsm_utils/Testlists/bluefire.esmf.auxtest
A  +    ccsm_utils/Testlists/E01
M       create_test_suite

================================================================================        
Originator: dfeddema 
Date: Jan 20 2012
Model: scripts
Version: scripts4_120120
One-line:  Modified ERS tests to exit after first failure. 

M            33962   ERS_script
	

================================================================================        
Originator: jedwards
Date: Jan 13, 2012
Model: scripts
Version: scripts4_120113
One-line:  Port to titan, multiinstance bug fix

	M            33577   ccsm_utils/Tools/generate_batch.csh
X                    ccsm_utils/Tools/lnd/clm/PTCLM
A  +             -   ccsm_utils/Testlists/titan.clm.auxtest
D            33577   ccsm_utils/Testlists/jaguarpf.cam.auxtest
A  +             -   ccsm_utils/Testlists/titan.cam.auxtest
D            33577   ccsm_utils/Testlists/jaguarpf.prebeta
A  +             -   ccsm_utils/Testlists/titan.prebeta
D            33577   ccsm_utils/Testlists/jaguarpf.prerelease
D            33577   ccsm_utils/Testlists/jaguarpf.clm.auxtest
A  +             -   ccsm_utils/Testlists/titan.prerelease
M            33577   doc/machines_list.xml
M            33577   doc/chap5.xml


================================================================================        
Originator: mvertens
Date: Dec 27 2011
Model: scripts
Version: scripts4_111227
One-line:  Reintroduced RUN_REFTOD (inadvertently removed in scripts4_111226)

M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/config_definition.xml
	
================================================================================        
Originator: mvertens
Date: Dec 26 2011
Model: scripts
Version: scripts4_111226
One-line:  Made grid specification more general 

M       ccsm_utils/Tools/configure
	- simplified call to generate_cice_decomp.pl  
M       ccsm_utils/Case.template/config_grid.xml
        - put in new default settings and the new variables listed
 	  in config_grid.xml
M       ccsm_utils/Case.template/config_definition.xml
	- the following new xml variables 
            MAP_A2O_PATH 
            MAP_O2A_PATH 
            MAP_L2A_PATH 
            MAP_A2L_PATH 
            MAP_R2O_PATH 
	    - note that now all paths for mapping files are specified
            ICE_DOMAIN_PRIMEFACS 
             - currently only used by cice
             - used by cice utility generate_cice_decomp.pl to 
	       obtain ice block sizes
            - if a cartesian or space curve decomposition cannot 
	      be found based on the icd grid and number of mpi 
	      tasks - then...
              -if ICE_DOMAIN_PRIMEFACS is not zero, the list 
	       will be used to obtain a decomposition where 
	       CICE_DECOMPTYPE=roundrobin and
               CICE_BLCKX>0, CICE_BLCKY=1, CICE_MXBLCKS=1
             -if ICE_DOMAIN_PRIMEFACS is zero, a decomposition 
	      will be generated where 
	      CICE_DECOMPTYPE=roundrobin and
              CICE_BLCKX=1, CICE_BLCKY=1, CICE_MXBLCKS>1
M       create_newcase
        - added new argument -grid_file that permits user to specify 
	  user-specific grid and paths
	
================================================================================        
Originator: erik
Date: Dec 21 2011
Model: scripts
Version: scripts4_111221b
One-line: Add RUN_REFTOD, fix some single-point domain files

Update PTCLM
M       ccsm_utils/Tools/ccsm_msmkdir -- in case used set +t bit on directory create
M       ccsm_utils/Case.template/config_grid.xml ------- Fix some single-point domain files
M       ccsm_utils/Case.template/config_definition.xml - Add RUN_REFTOD
M       create_newcase -- small change so that you can change machine directory 
            and get the -list option to work.

Add some single-point tests to clm auxtest lists
M       ccsm_utils/Testlists/bluefire.clm.auxtest
M       ccsm_utils/Testlists/edinburgh.clm.auxtest
M       ccsm_utils/Testlists/jaguarpf.clm.auxtest

================================================================================        
Originator: mvertens
Date: Dec 21 2011
Model: scripts
Version: scripts4_111221
One-line:  Cleaned up more machine specific information that did not belong in 
	   Tools and introduced ne120_ne120 and ne30_ne30 (gx1v6 mask)
	
M       ccsm_utils/Tools/ccsm_postrun.csh
M       ccsm_utils/Tools/create_production_test_readme
M       ccsm_utils/Tools/timing/getTiming.csh
M       ccsm_utils/Tools/taskmaker.pl
M       ccsm_utils/Tools/create_production_test
M       ccsm_utils/Case.template/config_grid.xml
	
	
================================================================================        
Originator: jedwards
Date: Dec 20 2011
Model: scripts
Version: scripts4_111220
One-line:  Reverted change to make sure env is only read once - this caused 
	   problems when env files were changed by test scripts corrected
	   cpl history path in testcase_end

	M   ccsm_utils/Tools/ccsm_getenv
	M   ccsm_utils/Tools/testcase_end
	
================================================================================        
Originator: jedwards
Date: Dec 19 2011
Model: scripts
Version: scripts4_111219
One-line: Added code to ccsm_getenv to make sure that env is only loaded once,
	  Moved env_mach_specific to the top of that file.

	M   Tools/ccsm_getenv


	
================================================================================        
Originator: mvertens
Date: Dec 18 2011
Model: scripts
Version: scripts4_111218
One-line: removed machine name from .build,clean_build,.run,.l_archive,
	  .submit and .test scripts  
	
M       ccsm_utils/Tools/ccsm_buildexe.csh
M       ccsm_utils/Tools/ccsm_postrun.csh
M       ccsm_utils/Tools/testcase_env.csh
M       ccsm_utils/Tools/check_case
M       ccsm_utils/Tools/create_train
M       ccsm_utils/Tools/generate_batch.csh
M       ccsm_utils/Tools/testcase_setup.csh
M       ccsm_utils/Tools/create_production_test
	
M       ccsm_utils/Tools/configure
        - no longer need ccsm4test in configure
	
M       ccsm_utils/Tools/archive_metadata.sh
	- machine name is obtained only from Macros file 
	
M       ccsm_utils/Tools/testcase_end
	- tests now only compare history files
	
M       ccsm_utils/Case.template/config_definition.xml
	- moved GLC_NEC xml variable to env_conf
	
M       ccsm_utils/Testcases/ERI_script
	- short term archiving turned off just for last stage in order
	  to permit comparison of history files


	
================================================================================        
Originator: jedwards
Date: Dec 16 2011
Model: scripts
Version: scripts4_111216
One-line: add CESM_MAX_THREADS fix ConfigCase bug introduced in previous tag

M       ccsm_utils/Case.template/ConfigCase.pm
M       create_newcase


================================================================================        
Originator: jedwards
Date: Dec 14 2011
Model: scripts
Version: scripts4_111214
One-line: reformat config_machines.xml 

M       ccsm_utils/Case.template/ConfigCase.pm
M       create_newcase
M       create_clone	
M       Build/Makefile  (actually in scripts4_111213)

================================================================================        
Originator: tcraig
Date: Dec 13 2011
Model: scripts
Version: scripts4_111213
One-line: modify build script to buildnml, prestage, buildnml

M       ccsm_utils/Tools/generate_batch.csh
================================================================================        
Originator: fischer
Date: Dec 12 2011
Model: scripts
Version: scripts4_111212
One-line: add 32x64 to config_grid.xml

M       ccsm_utils/Case.template/config_grid.xml

================================================================================        
Originator: trey
Date: Dec 6 2011
Model: scripts
Version: scripts4_111206
One-line: add sticky bit to lt_archive.csh, add CESM1(CAM5) RCP compsets

M       ccsm_utils/Tools/lt_archive.csh
M       ccsm_utils/Case.template/config_compsets.xml

================================================================================        
Originator: tcraig
Date: Dec 2 2011
Model: scripts
Version: scripts4_111202
One-line: remove pop2 rpointer fix from NCK case, now in pop2 template

M       ccsm_utils/Testcases/NCK_script

================================================================================        
Originator: jshollen
Date: Nov 30 2011
Model: scripts
Version: scripts4_111130
One-line: Added /bin/csh -f to all scripts in Testcases

M       ccsm_utils/Testcases/ERU_script
M       ccsm_utils/Testcases/PMT_script
M       ccsm_utils/Testcases/SMS_script
M       ccsm_utils/Testcases/CME_script
M       ccsm_utils/Testcases/PEA_script
M       ccsm_utils/Testcases/PST_script
M       ccsm_utils/Testcases/ERB_script
M       ccsm_utils/Testcases/NCK_script
M       ccsm_utils/Testcases/ERH_script
M       ccsm_utils/Testcases/LAR_script
M       ccsm_utils/Testcases/ERI_script
M       ccsm_utils/Testcases/PEM_script
M       ccsm_utils/Testcases/P4A_script
M       ccsm_utils/Testcases/APT_script
M       ccsm_utils/Testcases/ERP_script
M       ccsm_utils/Testcases/PET_script
M       ccsm_utils/Testcases/PFS_script
M       ccsm_utils/Testcases/SEQ_script
M       ccsm_utils/Testcases/ERS_script
M       ccsm_utils/Testcases/ERT_script

================================================================================        
Originator: erik
Date: Nov 29 2011
Model: scripts
Version: scripts4_111129
One-line: Switch CLMFORC for CLMQIAN, Update PTCLM/perl5lib

 - remove use of DIN_LOC_ROOT_CLMQIAN add DIN_LOC_ROOT_CLMFORC (requires datm8_111115)
 - Update PTCLM to use NetCDF4
 - Update perl5lib to fix hanging problem

M       ccsm_utils/Case.template/config_definition.xml
M       doc/chap4.xml
M       create_newcase
	
================================================================================        
Originator: mvertens
Date: Nov 17 2011
Model: scripts
Version: scripts4_111127
One-line: fixed domain file problems for ne120 and ne240 
	
M       ccsm_utils/Case.template/config_grid.xml

================================================================================        
Originator: tcraig
Date: Nov 13 2011
Model: scripts
Version: scripts4_111113
One-line: add NCPL_BASE_PERIOD env variable
 - and Timing2 tools (not used or tested yet)

A       ccsm_utils/Tools/getTiming2.pl
A       ccsm_utils/Tools/getTiming2.csh
M       ccsm_utils/Case.template/config_definition.xml
================================================================================        
Originator: trey
Date: Nov 10 2011
Model: scripts
Version: scripts4_111110
One-line: added support for DOUT_L_HPSS_ACCNT to lt_archive.csh

Uses HPSS account specified by $DOUT_L_HPSS_ACCNT if it exists and it
larger then 0. Also modified the hsi command to "mkdir -p", creating
any missing directories in the path.
	
M       ccsm_utils/Tools/lt_archive.csh
	
================================================================================        
Originator: mvertens
Date: Nov 09 2011
Model: scripts
Version: scripts4_111109
One-line: add new grid specification for domain files and domain paths

As a result data models no longer need to have these files in the template scripts	
	
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/ConfigCase.pm
M       ccsm_utils/Case.template/config_definition.xml
M       create_newcase
	
================================================================================        
Originator: tcraig
Date: Nov 08 2011
Model: scripts
Version: scripts4_111108
One-line: add NINST_LAYOUT env vars and fix NCK pop rpointer problem

M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Testcases/NCK_script
================================================================================        
Originator: aliceb
Date: Nov 07 2011
Model: scripts
Version: scripts4_111107
One-line: updates to archive_metadata.sh

commented out prompts for case metadata.
Case metadata is now being entered and stored in 
mysql run db at http://csegweb.cgd.ucar.edu.

M      ccsm_utils/Tools/archive_metadata.sh

================================================================================
Originator: mvertens
Date: Nov 5 2011
Model: scripts
Version: scripts4_111105
One-line: Fix for FAMIPCN and for long term arhiver problem

M       ccsm_utils/Case.template/config_compsets.xml
M       create_newcase
	
================================================================================
Originator: mvertens
Date: Oct 31 2011
Model: scripts
Version: scripts4_111031a
One-line: Introduced a new match variable SSTICE_COMPSET_MATCH 

RES_COMPSET_MATCH must now match exactly to the compset longname
SSTICE_COMPSET_MATCH provides a partial match for SSTICE forcing datasets	
	
M       create_newcase
M       config_compsets.xml
M       ConfigCase.pm
	
================================================================================
Originator: mvertens
Date: Oct 31 2011
Model: scripts
Version: scripts4_111031
One-line: Backed out scripts4_111026 - incorrect implementation

M       create_newcase
	
================================================================================
Originator: jedwards
Date: Oct 27 2011
Model: scripts
Version: scripts4_111027
One-line: Fixed NCK test when component uses only 1 task, force mintask
count on intrepid to 512

	
M       create_newcase
M       ccsm_utils/Testcases/NCK_script

================================================================================
Originator: mvertens
Date: Oct 26 2011
Model: scripts
Version: scripts4_111026
One-line: Fixed matching rule for RES_COMPSET_MATCH

Need to match identically - not just part of the pattern
	
M       create_newcase

================================================================================
Originator: jedwards
Date: Oct 25 2011
Model: scripts
Version: scripts4_111025a
One-line: Fixed incorrect filename, added special case for task count on intrepid 

M       ccsm_utils/Case.template/config_compsets.xml
M       create_newcase

================================================================================        
Originator: mvertens
Date: Oct 25 2011
Model: scripts
Version: scripts4_111025
One-line: Fixed problem for initial date for 1850-2005 compsets

M       ccsm_utils/Case.template/config_compsets.xml	
	
================================================================================        
Originator: mvertens
Date: Oct 11 2011
Model: scripts
Version: scripts4_111011b
One-line: New  compset specification mechanism

Major rewrite of compset specification in config_compsets.xml
This is completely documented at the head of the new config_compsets.xml file

M       ccsm_utils/Tools/st_archive.sh
M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/ConfigCase.pm
M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Testlists/B01
M       ccsm_utils/Testlists/C01
M       ccsm_utils/Testlists/C02
M       ccsm_utils/Testlists/C04
M       create_newcase


	
================================================================================        
Originator: trey
Date: Oct 11 2011
Model: scripts
Version: scripts4_111011
One-line: New simple long-term archiver

Added a new long-term archiver. It cleans up after any failed runs of
"lt_archive.pl" and then puts all the files in the short-term archive
into HPSS and automatically deletes files successfully put. It does not
use temporary directories or htar.

A      ccsm_utils/Tools/lt_archive.csh

================================================================================        
Originator: mvertens
Date: Oct 10 2011
Model: scripts
Version: scripts4_111010

For the following test types,  
  APT, SEQ, PEA, PEM, PET, PST, PMT, ERS, ERP, ERT, LAR, PFS, SMS, NCK
tests that should be hybrid initializations had been changed to have RUN_TYPE 
hardwired to startup. This is now changed to have RUN_TYPE actually be hybrid - 
which will change answers for many baseline results

M       ccsm_utils/Testcases/config_tests.xml

================================================================================        
Originator: jedwards
Date: Sep 30 2011
Model: scripts
Version: scripts4_110930

        Replaced -machines_dir with -mach_dir, moved
	CCSM_MACHDIR from env_conf.xml to env_case.xml
	
	changes: 
	create_newcase
	create_clone
	create_test
	ccsm_utils/Case.template/config_definition.xml
	
================================================================================        
Originator: jedwards
Date: Sep 29 2011
Model: scripts
Version: scripts4_110929

	As per cseg meeting discussion, add command line option -machines_dir
        which puts env variable CCSM_MACHDIR into env_conf.xml
	
	changes: 
	create_newcase
	create_clone
	create_test
	ccsm_utils/Tools/generate_batch.csh
	ccsm_utils/Tools/testcase_setup.csh
	ccsm_utils/Case.template/config_definition.xml
	
================================================================================        
Originator: jedwards
Date: Sep 28 2011
Model: scripts
Version: scripts4_110928
One-line: Add gptl.buildlib and build timing as separate lib, remove Machines

     create_clone
     create_newcase

	Added CCSM_MACHDIR variable which allows the Machines directory to be moved
	and set the default to ccsm_utils/Machines
	
    ccsm_utils/Components/gptl.buildlib - ADDED

    ccsm_utils/Machines  - DELETED

-------

scripts4_110922a Changelog entry missing:
	ccsm_utils/Tools/ccsm_prestage.csh: code added to deal with cam2 -> cam
	restart file renaming and maintaing backward compatability.

	
================================================================================        
Originator: erik
Date: Sep 22 2011
Model: scripts
Version: scripts4_110922
One-line: Do cleanall in NCK test

  M ccsm_utils/Testcases/NCK_script ----------- Do cleanall after changing
          NINST_* values as this is required for clm
  M ccsm_utils/Testlists/bluefire.clm.auxtest - Add NCK test in

================================================================================        
Originator: jedwards
Date: Sep 12 2011
Model: scripts
Version: scripts4_110912a
One-line: Update bluefire netcdf path to 4.1.3_seq

        ccsm_utils/Machines/Macros.bluefire
	ccsm_utils/Machines/env_machopts.bluefire
================================================================================        
Originator: erik
Date: Sep 12 2011
Model: scripts
Version: scripts4_110912
One-line: Add ne120 ERS test for I compset to B41

M       ccsm_utils/Testlists/B41 - Switch SMS test at ne120 for ERS test
M       ccsm_utils/Testlists/jaguarpf.clm.auxtest - Add more tests

================================================================================        
Originator: erik
Date: Sep 8 2011
Model: scripts
Version: scripts4_110908
One-line: Switch some ne30_g16 tests in

>>>>>>>>>>>> Switch some ne30_f19_g16 tests for ne30_g16
M       ccsm_utils/Testlists/B01
M       ccsm_utils/Testlists/C01
M       ccsm_utils/Testlists/C02
M       ccsm_utils/Testlists/C03

M       ccsm_utils/Testlists/B41 - Add ne120_g16 test

================================================================================        
Originator: erik
Date: Sep 6 2011
Model: scripts
Version: scripts4_110906
One-line: Add HOMME grid support

Move nonlatlonclm branch to trunk. Requires clm4_0_35 and datm8_110906.

Update PTCLM to PTCLM1_110902

>>>>>>>>>>>>
M       create_newcase - Make sure NTASKS_CPL is set for pecount (fix bug 1403)

>>>>>>>>>>>>
M       ccsm_utils/Build/Makefile ----------------- Change use of strip to compare to null
M       ccsm_utils/Case.template/config_grid.xml -- Add ne30np4_gx1v6, ne240np4_gx1v6, ne120np4_gx1v6
                                                    update ne120np4_0.9x1.25_gx1v6 files
M       ccsm_utils/Testlists/bluefire.clm.auxtest - switch f09 test to ERB.ne30_g16.I_1948-2004.bluefire

>>>>>>>>>>>> Add ne30np4_gx1v6 and ne120np4_gx1v6 tests
M       ccsm_utils/Testlists/C01 - add SMS.ne30_f19_g16.A
M       ccsm_utils/Testlists/C02 - add ERS.ne30_g16.F1850CN, SMS_D.ne30_f19_g16.A
M       ccsm_utils/Testlists/C03 - add ERS.ne30_g16.B1850CN, ERS.ne30_f19_g16.A
M       ccsm_utils/Testlists/C04 - add SMS_D.ne30_16.FCN,    ERS_D.ne30_16.A
M       ccsm_utils/Testlists/C41 - add ERS.ne120_g16.B
M       ccsm_utils/Testlists/C42 - add SMS.ne120_g16.BCN
M       ccsm_utils/Testlists/C43 - add ERS.ne120_g16.F1850
M       ccsm_utils/Testlists/C44 - add SMS.ne120_g16.FCN

================================================================================        
Originator: fischer
Date: Sep 1 2011
Model: scripts
Version: scripts4_110901
One-line:   removed extra grid match for B20TRC5CN, and fix cam filename matching

M       ccsm_utils/Tools/st_archive.sh
M       ccsm_utils/Case.template/config_compsets.xml

================================================================================        
Originator: tcraig
Date: Aug 30 2011
Model: scripts
Version: scripts4_110830
One-line:   multi-instance fixes and features
  - add _N confopts for multi-ensemble setup
  - add NCK test
  - add NCK and _N* tests to testlists
  - add CPL_DECOMP env var 
  - fix st_archive script for multi-ensemble
	
M       create_test
M       ccsm_utils/Tools/configure
M       ccsm_utils/Tools/st_archive.sh
M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Testcases/config_tests.xml
A       ccsm_utils/Testcases/NCK_script
M       ccsm_utils/Testlists/C01
M       ccsm_utils/Testlists/C02
M       ccsm_utils/Testlists/C03
M       ccsm_utils/Testlists/C04
M       create_newcase

================================================================================        
Originator: jedwards
Date: Thursday Aug 11 2011
Model: scripts
Version: scripts4_110811
One-line:   Update Franklin Macros and add a fix for clm standalone build
	   Machines/Macros.franklin
	   Machines/env_machopts.franklin
	   Machines/Macros.postdefs
	
================================================================================        
Originator: jedwards
Date: Wednesday Aug 10 2011
Model: scripts
Version: scripts4_110810
One-line: Update Makefile to correct parallel build errors
	
	   ccsm_utils/Build/Makefile
================================================================================        
Originator: jedwards
Date: Sunday Aug  7 2011
Model: scripts
Version: scripts4_110807
One-line: Update all Macros files for consistancy.
	
	   ccsm_utils/Build/Makefile
	   ccsm_utils/Case.template/config_definition.xml
	   ccsm_utils/Machines/Macros.brutus_im
	   ccsm_utils/Machines/Macros.brutus_io
	   ccsm_utils/Machines/Macros.pingo
	   ccsm_utils/Machines/Macros.intrepid
	   ccsm_utils/Machines/Macros.generic_linux_pathscale
	   ccsm_utils/Machines/Macros.generic_darwin_pgi
	   ccsm_utils/Machines/Macros.prototype_columbia
	   ccsm_utils/Machines/Macros.hector
	   ccsm_utils/Machines/Macros.prototype_atlas
	   ccsm_utils/Machines/Macros.hadley
	   ccsm_utils/Machines/Macros.generic_darwin_intel
	   ccsm_utils/Machines/Macros.generic_linux_lahey
	   ccsm_utils/Machines/Macros.pleiades
	   ccsm_utils/Machines/Macros.generic_xt
	   ccsm_utils/Machines/Macros.generic_linux_pgi
	   ccsm_utils/Machines/Macros.chester
	   ccsm_utils/Machines/Macros.generic_ibm
	   ccsm_utils/Machines/Macros.midnight
	   ccsm_utils/Machines/Macros.pleiades_wes
	   ccsm_utils/Machines/Macros.brutus_pm
	   ccsm_utils/Machines/Macros.prototype_frost
	   ccsm_utils/Machines/Macros.brutus_po

================================================================================        
Originator: jedwards
Date: Fri Aug  5 2011
Model: scripts
Version: scripts4_110805
One-line: Clean up build issues
	Index: ccsm_utils/Build/Makefile
	Index: ccsm_utils/Components/mct.buildlib
	Index: ccsm_utils/Machines/Macros.generic_linux_intel
	Index: ccsm_utils/Machines/Macros.edinburgh_intel
	Index: ccsm_utils/Machines/Macros.franklin
	Index: ccsm_utils/Machines/Macros.edinburgh_pgi
	Index: ccsm_utils/Machines/Macros.postdefs

================================================================================        
Originator: jedwards
Date: Thu Aug  4 2011
Model: scripts
Version: scripts4_110804b
One-line: updated build for lynx and edinburgh
	M            29628   ccsm_utils/Components/pio.buildlib
	M            29628   ccsm_utils/Machines/Macros.lynx_pgi
	M            29628   ccsm_utils/Machines/Macros.edinburgh_lahey
	M            29628   ccsm_utils/Machines/Macros.postdefs


	
================================================================================        
Originator: jedwards
Date: Thu Aug  4 2011
Model: scripts
Version: scripts4_110804a
One-line: updated build

	This change allows pio and mct to be configured and built outside
	of the source directory and avoids having to copy the entire source
	trees for these libraries.
	
COORDINATION REQUIREMENT: scripts4_110804a + pio1_3_5 + MCT2_7_0_110804

M            29606   ccsm_utils/Build/Makefile
M            29606   ccsm_utils/Components/mct.buildlib
M            29606   ccsm_utils/Components/pio.buildlib
M            29606   ccsm_utils/Machines/Macros.postdefs


================================================================================        
================================================================================        
Originator: kauff
Date: Thu Aug  4 11:51:19 MDT 2011
Model: scripts
Version: scripts4_110804
One-line: new/reworked DLND_RUNOFF_MODE, DATM_MODE, DLND_MODE's

DATM_MODE = CORE2 becomes CORE2_NYF + new CORE2_AIF mode
DLND_RUNOFF_MODE = RX1 becomes DIATREN_ANN_RX1 + new DIATREN_IAF_RX1

COORDINATION REQUIREMENT: drvseq4_0_03 + dlnd8_110803 + datm8_110803 + scripts4_110804

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_definition.xml

================================================================================        
Originator: tcraig
Date: Aug 3, 2011
Model: scripts
Version: scripts4_110803
One-line: Update pe layouts for intrepid

M      ccsm_utils/Machines/config_pes.xml

================================================================================        
Originator: jedwards
Date: July 28, 2011
Model: scripts
Version: scripts4_110728
One-line: Update Macros on hopp2 and intrepid, change A02 test from I to
Itest

	
================================================================================        
Originator: erik
Date: July 26, 2011
Model: scripts
Version: scripts4_110726
One-line: Update externals and allows testing on generic machines

Add ability to set the following in create_test and create_test_suite
and pass the input onto create_newcase. This allows create_test and create_test_suite
to be used on generic machines.

+     -scratchroot scratch-directory 
+     -din_loc_root_csmdata inputdata-directory  
+     -max_tasks_per_node number

Externals: Update perl5lib to perl5lib_110726
           Update PTCLM to PTCLM1_110726

M       ccsm_utils/Case.template/config_compsets.xml  Change MOAR case pointed to
           for I1850SPINUPCN compset so it's an 1850 spinup case
A       ccsm_utils/Testlists/intrepid.clm.auxtest --- Add some tests for intrepid
M       ccsm_utils/Testlists/jaguarpf.clm.auxtest --- Change bluefire test to jaguarpf

M       create_test --- Pass scratchroot, din_loc_root_csmdata and max_tasks_per_node
            to create_newcase only if a generic machine is being used
M       create_test_suite - Pass scratchroot, din_loc_root_csmdata and max_tasks_per_node
            if set to create_test

================================================================================        
Originator: jedwards
Date: June 24, 2011
Model: scripts
Version: scripts4_110724
One-line: improve support for Cray systems - can now configure with 5 compilers on lynx

	Index: create_test
	  Added ability to build only without job submission (Tony)     
	
	Index: ccsm_utils/Tools/clean_build
	   Added code to remove autoconf cache used by mct and pio
	
	Index: ccsm_utils/Components/pio.buildlib
	    Need to copy pio directory recursively so that m4 directory is included
	
	Index: ccsm_utils/Machines/Macros.prototype_hera
	Index: ccsm_utils/Machines/Macros.brutus_im
	Index: ccsm_utils/Machines/Macros.brutus_io
	Index: ccsm_utils/Machines/Macros.pingo
	Index: ccsm_utils/Machines/Macros.generic_linux_pathscale
	Index: ccsm_utils/Machines/Macros.prototype_nyblue
	Index: ccsm_utils/Machines/Macros.prototype_atlas
	Index: ccsm_utils/Machines/Macros.edinburgh_intel
	Index: ccsm_utils/Machines/Macros.lynx_cray
	Index: ccsm_utils/Machines/Macros.generic_linux_lahey
	Index: ccsm_utils/Machines/Macros.pleiades
	Index: ccsm_utils/Machines/Macros.lynx_pgi
	Index: ccsm_utils/Machines/Macros.hopper
	Index: ccsm_utils/Machines/Macros.franklin
	Index: ccsm_utils/Machines/Macros.edinburgh_pgi
	Index: ccsm_utils/Machines/Macros.postdefs
	Index: ccsm_utils/Machines/Macros.chester
	Index: ccsm_utils/Machines/Macros.midnight
	Index: ccsm_utils/Machines/Macros.generic_ibm
	Index: ccsm_utils/Machines/Macros.intrepid
	Index: ccsm_utils/Machines/Macros.generic_darwin_pgi
	Index: ccsm_utils/Machines/Macros.prototype_columbia
	Index: ccsm_utils/Machines/Macros.generic_linux_intel
	Index: ccsm_utils/Machines/Macros.hadley
	Index: ccsm_utils/Machines/Macros.kraken
	Index: ccsm_utils/Machines/Macros.generic_darwin_intel
	Index: ccsm_utils/Machines/Macros.generic_xt
	Index: ccsm_utils/Machines/Macros.hopp2
	Index: ccsm_utils/Machines/Macros.edinburgh_lahey
	Index: ccsm_utils/Machines/Macros.generic_linux_pgi
	Index: ccsm_utils/Machines/Macros.bluefire
	Index: ccsm_utils/Machines/Macros.prototype_ranger
	Index: ccsm_utils/Machines/Macros.jaguar
	Index: ccsm_utils/Machines/Macros.prototype_ubgl
	Index: ccsm_utils/Machines/Macros.lynx_pathscale
	Index: ccsm_utils/Machines/Macros.lynx_intel
	Index: ccsm_utils/Machines/Macros.jaguarpf
	Index: ccsm_utils/Machines/Macros.lynx_gnu
	Index: ccsm_utils/Machines/Macros.pleiades_wes
	Index: ccsm_utils/Machines/Macros.brutus_pm
	Index: ccsm_utils/Machines/Macros.brutus_po
	Index: ccsm_utils/Machines/Macros.prototype_frost

        Changed interface to MCT/PIO configure.  Moved common features of Macros files to Macros.postdefs
	extended support to cray,pathscale,intel and gnu (gfortran) compilers on lynx.   This can be used 
	as a template to extend support on similar cray platforms.

	Removed -gopt from default flags in pgi builds - this flag significantly increases the size of the executable
	and can blow memory with little or no benefit otherwise.  

	
	Index: ccsm_utils/Machines/env_machopts.lynx_pgi
	Index: ccsm_utils/Machines/env_machopts.lynx_pathscale
	Index: ccsm_utils/Machines/env_machopts.lynx_intel
	Index: ccsm_utils/Machines/config_machines.xml

	Added support for cray compilers

	
	Index: ccsm_utils/Machines/env_machopts.hopper
	Index: ccsm_utils/Machines/env_machopts.jaguar
	Index: ccsm_utils/Machines/mkbatch.hopper	
	Index: ccsm_utils/Machines/mkbatch.jaguar

	Removed obsolete entrys


	Index: create_newcase
	
	   Added support for Macros.postdefs
================================================================================        
Originator: tcraig
Date: June 10, 2011
Model: scripts
Version: scripts4_110610
One-line: improve test submission and reporting

M       create_test
M       create_test_suite
	
================================================================================        
Originator: tcraig
Date: June 9, 2011
Model: scripts
Version: scripts4_110609
One-line: minor cleanup issues plus lynx and intrepid updates

- remove duplicate F1850C5 compset
- remove F1850CNMAM3 compset
- add valid values for PIO_TYPENAME env variables
- add setenv DVS_MAXNODES 1 for lynx
- change lynx run dir to /glade/scratch from /ptmp
- switch intrepid test to "PL" for ERT_PL.T62_g16.D.intrepid
- increase intrepid time request for larger jobs to 360

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Machines/env_machopts.lynx_pgi
M       ccsm_utils/Machines/env_machopts.lynx_pathscale
M       ccsm_utils/Machines/env_machopts.lynx_intel
M       ccsm_utils/Machines/mkbatch.intrepid
M       ccsm_utils/Machines/config_machines.xml
M       ccsm_utils/Testlists/intrepid.levelC
	
================================================================================        
Originator: tcraig
Date: June 1, 2011
Model: scripts
Version: scripts4_110601
One-line:  updated testlists, add features
  features include 
    - adding -confopts to create_newcase and having create_test use it
    - saving more cpl.log files to baselines
    - removing reference to gx3v5
    - add _CG, _AOA, _AOE confopts options for gregorian calendar and aofluxgrid
	

M       create_test
M       ccsm_utils/Tools/testcase_env.csh
M       ccsm_utils/Tools/testcase_end
X       ccsm_utils/Tools/lnd/clm/PTCLM
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Machines/config_machines.xml
M       ccsm_utils/Testcases/ERU_script
M       ccsm_utils/Testcases/PMT_script
M       ccsm_utils/Testcases/SMS_script
M       ccsm_utils/Testcases/CME_script
M       ccsm_utils/Testcases/PEA_script
M       ccsm_utils/Testcases/PST_script
M       ccsm_utils/Testcases/ERB_script
M       ccsm_utils/Testcases/ERH_script
M       ccsm_utils/Testcases/ERI_script
M       ccsm_utils/Testcases/PEM_script
M       ccsm_utils/Testcases/P4A_script
M       ccsm_utils/Testcases/APT_script
M       ccsm_utils/Testcases/ERP_script
M       ccsm_utils/Testcases/PET_script
M       ccsm_utils/Testcases/PFS_script
M       ccsm_utils/Testcases/SEQ_script
M       ccsm_utils/Testcases/ERS_script
M       ccsm_utils/Testcases/ERT_script
D       ccsm_utils/Testlists/jaguarpf.posttag
D       ccsm_utils/Testlists/jaguarpf.pretag
A       ccsm_utils/Testlists/edinburgh.levelC
A       ccsm_utils/Testlists/jaguarpf.cam.auxtest
D       ccsm_utils/Testlists/kraken.posttag
D       ccsm_utils/Testlists/jaguar.esmf.auxtest
D       ccsm_utils/Testlists/midnight.auxtest
D       ccsm_utils/Testlists/jaguarpf.esmf.release
A       ccsm_utils/Testlists/jaguarpf.levelB
D       ccsm_utils/Testlists/pingo.auxtest
A       ccsm_utils/Testlists/jaguarpf.levelC
A       ccsm_utils/Testlists/jaguarpf.levelD
D       ccsm_utils/Testlists/lynx.pretag
D       ccsm_utils/Testlists/hopper.posttag
D       ccsm_utils/Testlists/bluefire.posttag
A       ccsm_utils/Testlists/franklin.levelC
D       ccsm_utils/Testlists/chester.pretag
D       ccsm_utils/Testlists/bluefire.pretag
D       ccsm_utils/Testlists/franklin.posttag
D       ccsm_utils/Testlists/jaguarpf.esmf.auxtest
A       ccsm_utils/Testlists/kraken.levelD
D       ccsm_utils/Testlists/edinburgh.posttag
A       ccsm_utils/Testlists/intrepid.levelC
D       ccsm_utils/Testlists/edinburgh.pretag
D       ccsm_utils/Testlists/jaguar.posttag
A       ccsm_utils/Testlists/chester.auxtest
D       ccsm_utils/Testlists/jaguar.pretag
D       ccsm_utils/Testlists/intrepid.posttag
A       ccsm_utils/Testlists/hopp2.levelC
A       ccsm_utils/Testlists/lynx.levelA
A       ccsm_utils/Testlists/lynx.levelB
A       ccsm_utils/Testlists/lynx.levelC
D       ccsm_utils/Testlists/jaguar.cam.auxtest
A       ccsm_utils/Testlists/bluefire.levelA
D       ccsm_utils/Testlists/bluefire.esmf.release
A       ccsm_utils/Testlists/bluefire.levelB
A       ccsm_utils/Testlists/bluefire.levelC
A       ccsm_utils/Testlists/bluefire.levelD
D       ccsm_utils/Testlists/hopp2.posttag
A       ccsm_utils/Testlists/edinburgh.levelA
M       create_newcase
M       create_test_suite

================================================================================        
Originator: fischer
Date: May 31, 2011
Model: scripts
Version: scripts4_110531b
One-line:  updated jaguarpf ESMF library, fixed B1850BDRD CO2 variable name

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Machines/env_machopts.jaguarpf

================================================================================
Originator: jedwards
Date: May 31, 2011
Model: scripts
Version: scripts4_110531a
One-line: Add clm restart history files to short term archiver
   M ccsm_utils/Tools/st_archive.sh
================================================================================

Originator: fischer
Date: May 31, 2011
Model: scripts
Version: scripts4_110531
One-line: Update test list, new RUN_REFDATE for BC5

M       ccsm_utils/Case.template/config_compsets.xml
.  B1850CAM5CN and B20TRC5CN have a new RUN_REFDATE of 059-01-01

M       ccsm_utils/Machines/env_machopts.jaguarpf
.  New ESMF library for PGI11.0.0

M       ccsm_utils/Testlists/jaguarpf.posttag
M       ccsm_utils/Testlists/jaguarpf.pretag
.  Add couple hi-res tests ne250_t12 and f02_f02

M       ccsm_utils/Testlists/bluefire.pretag
.  Add f19_g16 BG1850CN test

M       ccsm_utils/Machines/config_machines.xml
-  Set MAX_TASKS_PER_NODE="16" for edinburgh
  

================================================================================        
Originator: fischer
Date: May 26, 2011
Model: scripts
Version: scripts4_110526
One-line: Update test lists, add machine hector, fix B1850BDRD


A       ccsm_utils/Machines/Macros.hector
A       ccsm_utils/Machines/env_machopts.hector
A       ccsm_utils/Machines/mkbatch.hector
M       ccsm_utils/Machines/config_pes.xml
M       ccsm_utils/Machines/config_machines.xml
.  add hector

M       ccsm_utils/Testlists/lynx.pretag
M       ccsm_utils/Testlists/bluefire.posttag
.  Fix ERs vs ERS typo

M       ccsm_utils/Testlists/jaguarpf.esmf.release
M       ccsm_utils/Testlists/bluefire.esmf.release
.  Remove ESMF tests that are expected to fail

M       ccsm_utils/Case.template/config_compsets.xml
.  Add BGC_CO2_PPMV="284.7" to B1850BDRD compset
.  This changes answers for B1850BDRD

M       ccsm_utils/Case.template/config_grid.xml
.  Remove GLC_GRID

================================================================================        
Originator: jedwards
Date: May 25, 2011
Model: scripts
Version: scripts4_110525
One-line: Make hsi directories world readable by default (from bluefire)

M ccsm_utils/Tools/lt_archive.pl

================================================================================        
Originator: fischer
Date: May 20, 2011
Model: scripts
Version: scripts4_110520
One-line: updated PGI/9.0.4 to PGI/11.0.0 for the crays, except for franklin

M       ccsm_utils/Machines/env_machopts.kraken
M       ccsm_utils/Machines/env_machopts.generic_xt
M       ccsm_utils/Machines/env_machopts.jaguarpf
M       ccsm_utils/Machines/env_machopts.hopp2
.  updated to PGI/11.0.0

================================================================================        
Originator: fischer
Date: May 19, 2011
Model: scripts
Version: scripts4_110519a
One-line: Create ESMF tests list for release, update docs, fix FAMIPC5CN test

A       ccsm_utils/Testlists/jaguarpf.esmf.release
A       ccsm_utils/Testlists/bluefire.esmf.release
.  This is an ESMF test list for release testing.  It is pre and post
.  tag lists combined

M       ccsm_utils/Testlists/bluefire.pretag
.  The FAMIPC5CN test should be f09_f09, it was f19_f19

M       doc/env_build_list.xml
M       doc/env_conf_list.xml
M       doc/compsets_list.xml
.  Update documents by running write_docbook_file

================================================================================        
Originator: jedwards
Date: May 19, 2011
Model: scripts
Version: scripts4_110519
One-line: Clean up longterm archiving scripts, add ne240tx0.1 configuration

 ccsm_utils/Tools/ccsm_msread
 ccsm_utils/Tools/ccsm_mswrite
 ccsm_utils/Tools/ccsm_msls
 ccsm_utils/Tools/ccsm_msmkdir

 Removed old NCAR long term archiving tools
	
 ccsm_utils/Tools/ccsm_l_archive.csh

 Removed workaround on path length that was required for old bluefire system
	
 ccsm_utils/Case.template/config_grid.xml
 ccsm_utils/Machines/config_machines.xml
 
 Added ne240np4_0.23x0.31_tx0.1v2 configuration

 ccsm_utils/Machines/mkbatch.bluefire

 Change so that bluefire uses new lt_archive script

	
 ccsm_utils/Tools/lt_archive.pl
   	
 New tool for hsi archiving on bluefire 


================================================================================        
Originator: erik
Date: May 17, 2011
Model: scripts
Version: scripts4_110517
One-line: Add append mode to xmlchange and DATM_CPL variables

Addresses bugs: 1337, 1336, 1158, and 1108

Add an append mode option to xmlchange so that it will append the input value
to the end of existing setting (provided the value isn't UNSET). Only works for
string data. Adds a blank if there is existing data before appending the new value.

The ALIGN years are changed in the compsets: I20TR, and IG20TR. THIS WILL CHANGE
ANSWERS RELATIVE TO THE PREVIOIS TAG.

Add variables for handling running datm with CPLHIST data

   DATM_CPL_CASE ----- Case name of CPLHIST data to use
   DATM_CPL_YR_ALIGN - Simulation year to align with start year
   DATM_CPL_YR_START - Starting year of years to loop over
   DATM_CPL_YR_END --- Ending year of years to loop over

M       ccsm_utils/Tools/xmlchange --------------------- Add append mode
M       ccsm_utils/Case.template/config_compsets.xml --- Make sure DATM_CLMNCEP_ALIGN
            is set, add DATM_CPL vars to I1850SPINUPCN compset
            fix documentation in I compsets, make ALIGN year consistent
M       ccsm_utils/Case.template/ConfigCase.pm --------- Add is_char method
M       ccsm_utils/Case.template/config_definition.xml - Add DATM_CPL_ vars

Work on coverage of tests for I compsets. Add a few tests to make sure all compsets
are tested. Change a few others to give better coverage.

M       ccsm_utils/Testlists/bluefire.clm.auxtest - Make I1850 test I1850CN
M       ccsm_utils/Testlists/lynx.pretag ---------- Add I20TRCN@T31 test
M       ccsm_utils/Testlists/hopper.posttag ------- Change I to ITEST compset
M       ccsm_utils/Testlists/bluefire.posttag ----- Add I20TRCN@f19 test
M       ccsm_utils/Testlists/edinburgh.pretag ----- Add T31 ICN, I20TR tests
M       ccsm_utils/Testlists/edinburgh.clm.auxtest  Add f10 IRCP26CN test
M       ccsm_utils/Testlists/jaguarpf.clm.auxtest - Add I4804CN test, replace
              IRCP85CN test with IRCP60CN test
A  +    ccsm_utils/Testlists/bluefire.clm.trans.auxtest - Rename from rcps
D       ccsm_utils/Testlists/bluefire.clm.rcps.auxtest -- Rename to trans

================================================================================        
Originator: mvertens
Date: May 16, 2011
Model: scripts
Version: scripts4_110516
One-line: fixed T31_g37 B1850CN out of the box functionality

M       README
	- started updating new scientifically supported compsets
M       ccsm_utils/Case.template/config_compsets.xml
	- fixed T31_g37 B1850CN out of the box functionaltiy 
M       ccsm_utils/Machines/config_machines.xml 
        - New path to DIN_LOC_ROOT_CLMQIAN for lynx (from erik)

================================================================================        
Originator: jedwards
Date: May 13, 2011
Model: scripts
Version: scripts4_110513
One-line: Update pretag tests, mkbatch for intrepid
M ccsm_utils/Machines/mkbatch.intrepid
M ccsm_utils/Testlists/lynx.pretag
M ccsm_utils/Testlists/bluefire.posttag

================================================================================
Originator: erik
Date: May 11, 2011
Model: scripts
Version: scripts4_110511
One-line: Fix resolution in testlist

M       ccsm_utils/Testlists/bluefire.clm.auxtest -- Fix typo in resolution for test

================================================================================        
Originator: erik
Date: May 10, 2011
Model: scripts
Version: scripts4_110510
One-line: Add back in changes to I2000 compsets and testlists with ITEST compset

Add a new testlist that runs the rest of the RCP's for CLM

A       ccsm_utils/Testlists/bluefire.clm.rcps.auxtest

Reverts these to scripts4_110505.

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Testlists/kraken.posttag
M       ccsm_utils/Testlists/hadley.auxtest
M       ccsm_utils/Testlists/franklin.posttag
M       ccsm_utils/Testlists/intrepid.posttag
M       ccsm_utils/Testlists/hopp2.posttag

================================================================================        
Originator: erik
Date: May 6, 2011
Model: scripts
Version: scripts4_110506
One-line: Backout changes to I2000 compsets, and testlists with new ITEST compset

This reverts the following files to scripts4_110502

M       ccsm_utils/Case.template/config_compsets.xml - The exception here is that
              I left the new ITEST compsets in.
M       ccsm_utils/Testlists/kraken.posttag
M       ccsm_utils/Testlists/hadley.auxtest
M       ccsm_utils/Testlists/franklin.posttag
M       ccsm_utils/Testlists/intrepid.posttag
M       ccsm_utils/Testlists/hopp2.posttag

================================================================================        
Originator: tcraig
Date: May 5, 2011
Model: scripts
Version: scripts4_110505a
One-line:  Add 2 tests to bluefire.drv.auxtest

M       ccsm_utils/Testlists/bluefire.drv.auxtest

================================================================================        
Originator: erik
Date: May 5, 2011
Model: scripts
Version: scripts4_110505
One-line:  Move PTCLM to external, add doc on MPISERIAL, add ITEST compsets

Make PTCLM an external. Add more documentation about enabling MPISERIAL_SUPPORT.
Add I_TEST_2003, I_TEST_2003_CN compsets, make I2000 compsets run over 1972 to 2004.

D       ccsm_utils/Tools/lnd/clm/PTCLM/*... -- Remove PTCLM make external
...

M       ccsm_utils/Tools/generate_batch.csh -- add line about MPISERIAL_SUPPORT
M       ccsm_utils/Case.template/config_compsets.xml --- New compsets, I2000 run over 1972-2004
M       ccsm_utils/Case.template/config_definition.xml - More doc on MPISERIAL_SUPPORT

M       ccsm_utils/Testlists/bluefire.clm.auxtest --- Switch some tests to ones that
            can work

>>>>>>>>>>>> Use ITEST compset in place of I tests
M       ccsm_utils/Testlists/kraken.posttag
M       ccsm_utils/Testlists/hadley.auxtest
M       ccsm_utils/Testlists/franklin.posttag
M       ccsm_utils/Testlists/newmachine.port.auxtest
M       ccsm_utils/Testlists/intrepid.posttag
M       ccsm_utils/Testlists/hopp2.posttag

>>>>>>>>>>>> Change jaguar list to jaguarpf list
D       ccsm_utils/Testlists/jaguar.clm.auxtest
A  +    ccsm_utils/Testlists/jaguarpf.clm.auxtest

M       doc/chap7.xml --- Add more on MPISERIAL_SUPPORT and USE_MPISERIAL.
>>>>>>>>>>>> Update lists since last time
M       doc/env_build_list.xml
M       doc/grids_list.xml
M       doc/machines_list.xml
M       doc/env_conf_list.xml
M       doc/compsets_list.xml
M       doc/env_run_list.xml

================================================================================        
Originator: fischer
Date: May 2, 2011
Model: scripts
Version: scripts4_110502
One-line: fixed DOC_SSTDATA_YEAR_END=0, added intrepid tests 

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Testlists/intrepid.posttag

================================================================================        
Originator: tcraig
Date: Apr 28, 2011
Model: scripts
Version: scripts4_110428a
One-line: update ESMF lib on jaguarpf

M       ccsm_utils/Machines/env_machopts.jaguarpf
	-  module load esmf/4.0.0rp2_O
	+  module load esmf/5.2.0_O 

================================================================================        
Originator: fischer
Date: Apr 28, 2011
Model: scripts
Version: scripts4_110428
One-line: added new tests, update COSP build dependencies

M       ccsm_utils/Machines/Macros.hopp2
M       ccsm_utils/Machines/Macros.bluefire
M       ccsm_utils/Machines/Macros.franklin
M       ccsm_utils/Machines/Macros.jaguar
M       ccsm_utils/Machines/Macros.jaguarpf
 update COSP build dependencies

M       ccsm_utils/Testlists/bluefire.pretag
 added tests for new compsets

================================================================================        
Originator: Francis Vitt
Date: Apr 27, 2011
Model: scripts
Version: scripts4_110427
One-line: added new compset and adjusted compiler options for pleiades

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Machines/config_pes.xml
M       ccsm_utils/Machines/Macros.pleiades
M       ccsm_utils/Machines/Macros.pleiades_wes

================================================================================        
Originator: hannay
Date: Apr 21, 2011
Model: scripts
Version: scripts4_110421
One-line: update CAM5 compsets:  create E_1850_CAM5_CN

M      ccsm_utils/Case.template/config_compsets.xml
================================================================================        
Originator: fischer,fvitt
Date: Apr 20, 2011
Model: scripts
Version: scripts4_110420
One-line: update WACCM compsets F_2000_WACCM_GHG F_SD_WACCM

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Machines/config_pes.xml
================================================================================        
Originator: tcraig
Date: Apr 19, 2011
Model: scripts
Version: scripts4_110419
One-line: update ESMF lib version

M       ccsm_utils/Machines/config_machines.xml
        ESMF_LIBDIR="/ptmp/svasquez/esmf_install/ESMF_5_2_0-O/lib/"

================================================================================        
Originator: Francis Vitt
Date: Apr 6, 2011
Model: scripts
Version: scripts4_110406
One-line: add/modify CAM-Chem and WACCM compsets and changes for pleiades

M       ccsm_utils/Tools/ccsm_mswrite
 - added remover file command to pleiades write to mass store operation
	
M       ccsm_utils/Tools/st_archive.sh
 - copy *.hs.* (cam history sampling output file) to st_archive

M       ccsm_utils/Tools/ccsm_buildnml.csh
 - copy chemistry mechanism documents to CaseDocs sub-directory

M       ccsm_utils/Tools/testcase_setup.csh
 - corrected spelling of "invocation"

M       ccsm_utils/Case.template/config_compsets.xml
 - added F_2000_TROP_MOZART, F_2000_WACCM_GHG, F_SD_CAMCHEM, and F_SD_BAM compsets
 - added cam build-namelist use case to B_2000_TROP_MOZART 
 - changed "waccm_mozart" to "waccm_mozart_v1" for the waccm AR5 compsets
 - added B_RCP2.6_WACCM_CN and  B_RCP8.5_WACCM_CN
 - removed "hybrid" run type from F_1955-2005_WACCM_CN

M       ccsm_utils/Machines/Macros.pleiades
M       ccsm_utils/Machines/Macros.pleiades_wes
 - compiler options change

M       ccsm_utils/Machines/env_machopts.pleiades
M       ccsm_utils/Machines/env_machopts.pleiades_wes
 - use intel 11 compiler

M       ccsm_utils/Machines/config_pes.xml
 - provide default pe configurations for F_ CAM-Chem/WACCM compsets on bluefire

================================================================================        
Originator: hannay
Date: Mar 28, 2011
Model: scripts
Version: scripts4_110328a
One-line: add/modify CAM5 compsets

M       ccsm_utils/Case.template/config_compsets.xml
	- new compsets for B_1850-2000_CAM5_CN,  F_AMIP_CAM5_CN
        - update initial conditions for B_1850_CAM5, B_1850_CAM5_CN, B_1850-2000_CAM5 
	  F_AMIP_CAM5
	
================================================================================
Originator: fischer
Date: Mar 28, 2011
Model: scripts
Version: scripts4_110328
One-line: update Machines and testlists

M      ccsm_utils/Tools/st_archive.sh
      removed extra lines of code

A      ccsm_utils/Machines/mkbatch.chester
A      ccsm_utils/Machines/env_machopts.chester
A      ccsm_utils/Machines/Macros.chester
A      ccsm_utils/Testlists/chester.pretag
M      ccsm_utils/Machines/config_pes.xml
      added chester for testing

M      ccsm_utils/Machines/config_machines.xml
      changed inputdata dir on lynx to glade
      set up hopper to display error message when selected

M      ccsm_utils/Machines/env_machopts.hopp2
      updated modules

M      ccsm_utils/Testlists/jaguarpf.posttag
      updated tests lists

M      create_newcase
      added code to display error message when invalid machine is selected (hopper)

================================================================================        
Originator: tcraig
Date: Mar 24, 2011
Model: scripts
Version: scripts4_110324a
One-line: modify timing dir management to address issues with long file lists

M       ccsm_utils/Tools/generate_batch.csh

================================================================================        
Originator: tcraig
Date: Mar 24, 2011
Model: scripts
Version: scripts4_110324
One-line: add SAVE_TIMING and change VECT_MAP homme grid default to cart3d

M       ccsm_utils/Tools/ccsm_postrun.csh
M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/config_definition.xml

================================================================================        
Originator: jpe
Date: Mar 22, 2011
Model: scripts
Version: scripts4_110322
One-line: Removed HIRES macro, modified ocean coupling time for tx0.1v2 grids
   doc/env_build_list.xml
   ccsm_utils/Tools/st_archive.sh
   ccsm_utils/Case.template/config_definition.xml
   ccsm_utils/Machines/Macros.cppdefs
   ccsm_utils/Case.template/config_grid.xml
     Removed references to HIRES, changed OCN_NCPL=4 for tx0.1v2 ocean grid
	
   ccsm_utils/Machines/Macros.bluefire
	added -qsclk=micro to build flags, addresses bug 477



================================================================================        
Originator: mvr
Date: Mar 14, 2011
Model: scripts
Version: scripts4_110314
One-line: added grid info for 2.5x3.33

M       ccsm_utils/Case.template/config_grid.xml

================================================================================        
Originator: tcraig
Date: Mar 08, 2011
Model: scripts
Version: scripts4_110308
One-line: update timing tool for drvseq3_1_50

- add bluefire.drv.auxtest
- modify automatic reset of CONTINUE_RUN on RESUBMIT if doing timing tests

M       ccsm_utils/Tools/ccsm_postrun.csh
M       ccsm_utils/Tools/timing/getTiming.pl
A       ccsm_utils/Testlists/bluefire.drv.auxtest
================================================================================        
Originator: mvr
Date: Mar 02, 2011
Model: scripts
Version: scripts4_110302
One-line: adding ability to specify charge account with hsi write command

M       ccsm_utils/Tools/ccsm_msread
M       ccsm_utils/Tools/ccsm_mswrite
M       ccsm_utils/Tools/ccsm_msmkdir
M       ccsm_utils/Case.template/config_definition.xml

================================================================================        
Originator: fischer
Date: Feb 17, 2011
Model: scripts
Version: scripts4_110217
One-line: moved jaguar pretag testing to jaguarpf, fixed FG test names 

A      ccsm_utils/Testlists/jaguarpf.pretag
M      ccsm_utils/Testlists/bluefire.pretag
M      ccsm_utils/Testlists/edinburgh.pretag
M      ccsm_utils/Testlists/jaguar.posttag
================================================================================        
Originator: fischer
Date: Feb 4, 2011
Model: scripts
Version: scripts4_110204
One-line:  fixed typo in ESMF, change BG to BGCN in compsets

M      ccsm_utils/Machines/config_machines.xml
M      ccsm_utils/Testlists/bluefire.pretag
M      ccsm_utils/Testlists/jaguar.pretag

================================================================================        
Originator: fischer
Date: Feb 3, 2011
Model: scripts
Version: scripts4_110203a
One-line:  updated create_train, change modules on jaguar, new pes for B hopp2, new ESMF lib

M      ccsm_utils/Tools/create_train
   updated to be machine independent

M      ccsm_utils/Machines/env_machopts.jaguar
M      ccsm_utils/Machines/env_machopts.jaguarpf
   module swap xt-asyncpe xt-asyncpe/3.7
M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Machines/config_machines.xml
  use ESMF_5_2_0_beta_snapshot_12_64-O on bluefire

================================================================================	
Originator: erik
Date: Feb 3, 2011
Model: scripts
Version: scripts4_110203
One-line:  Add new glc_mec compsets for 1850, transient and transient rcps

  Change I8520 compsets to I20TR to make similar to other compset names.

>>>>>>>>>> Add new glc_mec compsets
   M  ccsm_utils/Case.template/config_compsets.xml

>>>>>>>>>> Change tests based on new compsets
   M  ccsm_utils/Testlists/bluefire.esmf.auxtest - Change BG for BGCN test, IG for IGCN
   M  ccsm_utils/Testlists/bluefire.glc.auxtest -- Change BGCN/FGCN for BG/FG and
             add IG20TR, IGRCP85CN, FG20TRCN, and TG tests
   M  ccsm_utils/Testlists/bluefire.clm.auxtest -- Switch IG1850CN/I8520CN for IG/I20TRCN
             add IGRCP26CN, IGRCP60CN, IGRCP45CN
   M  ccsm_utils/Testlists/lynx.pretag ----------- Switch BG/FG for BGCN/FG1850CN
   M  ccsm_utils/Testlists/jaguar.clm.auxtest ---- Switch I85020CN for I20TRCN
   M  ccsm_utils/Testlists/edinburgh.pretag ------ Switch FG/IG for FG1850CN/IG1850
             add FGCN tests
   M  ccsm_utils/Testlists/jaguar.posttag -------- Switch FG/IG for FGCN1850/IG20TRCN

   M  ccsm_utils/Case.template/config_compsets.xsl - Update list to have list of _GLC
           compsets last

>>>>>>>>>> Some work on PTCLM, fix a few bugs,
>>>>>>>>>> Add ChangeLog/KnownBugs files, add CESM_ROOT and CLM_SOFF env var options to test
>>>>>>>>>> Add --scratchroot option, allow caseprefix to contain directory
>>>>>>>>>> Fix --list option, update for new clm4_0_23 use-case convention
   M  ccsm_utils/Tools/lnd/clm/PTCLM/PTCLM.py
   M  ccsm_utils/Tools/lnd/clm/PTCLM/testcases.csh
   A  ccsm_utils/Tools/lnd/clm/PTCLM/ChangeLog
   A  ccsm_utils/Tools/lnd/clm/PTCLM/KnownBugs

================================================================================	
Originator: tcraig
Date: Feb 1, 2011
Model: scripts
Version: scripts4_110201
One-line:  Add NPFIX and AOFLUX_GRID options (requires drvseq3_1_48)

M       ccsm_utils/Case.template/config_grid.xml
M       ccsm_utils/Case.template/config_definition.xml
M       ccsm_utils/Testlists/bluefire.pretag
	
================================================================================	
Originator: hannay (mvertens)
Date: Jan 11, 2011
Model: scripts
Version: scripts4_110121
One-line:  Update 1850 SSTs for F case and create compsets for CAM5_1850_CN

M      ccsm_utils/Case.template/config_compsets.xml

================================================================================	
Originator: fischer
Date: Jan 11, 2011
Model: scripts
Version: scripts4_110111
One-line:  Changed pgi version of lynx, change queue name on franklin, added IG and FG tests

M      ccsm_utils/Machines/mkbatch.franklin
M      ccsm_utils/Machines/env_machopts.lynx_pgi
a      ccsm_utils/Testlists/hopp2.posttag
M      ccsm_utils/Testlists/lynx.pretag
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/bluefire.pretag
M      ccsm_utils/Testlists/edinburgh.pretag
M      ccsm_utils/Testlists/jaguar.posttag

================================================================================	
Originator: tcraig
Date: Jan 08, 2011
Model: scripts
Version: scripts4_110108
One-line:  update esmf version on bluefire for new metadata capability

- this is just for short-term testing.  a proper release version will be
  available in the next month or two.

-         ESMF_LIBDIR="/contrib/esmf-4.0.0rp2-64/lib"
+         ESMF_LIBDIR="/ptmp/svasquez/esmf_install/ESMF_5_1_0_beta_snapshot_14_64-O/lib"
	
M       ccsm_utils/Machines/config_machines.xml
	
================================================================================	
Originator: Edwards
Date: Jan 06, 2011
Model: scripts
Version: scripts4_110106
One-line:  add flags required for BGP build
	ccsm_utils/Machines/Macros.intrepid

	
================================================================================	
Originator: fischer
Date: Dec 14 2010
Model: scripts
Version: scripts4_101214
One-line:  add support for hopp2, and increase CCSM_CCOST for WACCM

M      ccsm_utils/Case.template/config_compsets.xml

A      ccsm_utils/Machines/mkbatch.hopp2
A      ccsm_utils/Machines/Macros.hopp2
A      ccsm_utils/Machines/env_machopts.hopp2
M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Machines/config_machines.xml

================================================================================	
Originator: Jim
Date: Dec 10 2010
Model: scripts
Version: scripts4_101210
One-line:  Fix timing build, remove AIX flag from bluegene build
	        ccsm_utils/Build/Makefile
	        ccsm_utils/Machines/Macros.cppdefs
	        ccsm_utils/Machines/Macros.intrepid
	        ccsm_utils/Machines/Macros.prototype_frost
	        ccsm_utils/Machines/Macros.prototype_nyblue
	        ccsm_utils/Machines/Macros.prototype_ubgl

================================================================================
Originator: erik
Date: Dec 06 2010
Model: scripts
Version: scripts4_101206
One-line:  Fix run section for lynx from last tag

Also add a few comments about some of the changes made in the last tag.
Move setting of myaprun to run of template script instead of .run file.

M       ccsm_utils/Tools/lnd/clm/PTCLM/testcases.csh  Add USER_CC setting to yong
M       ccsm_utils/Tools/lnd/clm/PTCLM/README ------- Add note about aerdepgrid/ndepgrid
A       ccsm_utils/Tools/lnd/clm/PTCLM/KnownBugs ---- Add file with list of bugs

M       ccsm_utils/Machines/Macros.generic_linux_pathscale - comments
M       ccsm_utils/Machines/Macros.generic_darwin_intel ---- comments
M       ccsm_utils/Machines/mkbatch.lynx_pgi -------- mv setting of myaprun
M       ccsm_utils/Machines/mkbatch.lynx_pathscale -- mv setting of myaprun
M       ccsm_utils/Machines/mkbatch.lynx_intel ------ mv setting of myaprun

================================================================================
Originator: erik
Date: Dec 02 2010
Model: scripts
Version: scripts4_101202
One-line:  USE_MPISERIAL changes and rip out USER_ stuff

      - Add in changes for mpiserial making it easier to support MPISERIAL
      - Add simple run for USE_MPISERIAL for all .run scripts created
      - In all Macro's files if USE_MPISERIAL is TRUE set path to INC_MPI to
        mct/mpi-serial version.
      - Turn MPISERIAL_SUPPORT to TRUE for:
          generic machines, edinburgh*, intrepid, jaguar, jaguarpf, and lynx*
      - Fix bugs: 1209 (pt1 on edinburgh), 1187 (more USE_MPISERIAL), 1184 (pass FFLAGS
        to MCT configure)
      - Rip out USER_* and F_OPTIMIZATION stuff (other than USER_CPPDEFS)
      - Add support for darwin (pgi and intel) and lynx_pathscale
      - Add $GMAKE as a env_mach variable for name of GNU make (needed for darwin)
      - Change F2000_chem compset CCSM_CO2PPMV to 367.0 from Chris Fischer
      - PTCLM updates

>>>>>>>>>>>> Add support for darwin (pgi and intel) (intel working)
  A   ccsm_utils/Machines/Macros.generic_darwin_pgi
  A   ccsm_utils/Machines/env_machopts.generic_darwin_pgi 
  A   ccsm_utils/Machines/mkbatch.generic_darwin_pgi 
  A   ccsm_utils/Machines/mkbatch.generic_darwin_intel
  A   ccsm_utils/Machines/Macros.generic_darwin_intel
  A   ccsm_utils/Machines/env_machopts.generic_darwin_intel

>>>>>>>>>>>> Add support for lynx_pathscale
  A   ccsm_utils/Machines/Macros.lynx_pathscale
  A   ccsm_utils/Machines/env_machopts.lynx_pathscale
  A   ccsm_utils/Machines/mkbatch.lynx_pathscale

  M   ccsm_utils/Build/Makefile -- Remove USER_* stuff, move include MACFILE
          to before ESMF stuff

  M   ccsm_utils/Case.template/config_definition.xml -- Add GMAKE (needed for Darwin)
  M   ccsm_utils/Case.template/config_compsets.xml ---- Change F2000_chem compset
          CCSM_CO2PPMV to 367.0 from Chris Fischer

>>>>>>>>>>>> Change Macros files removing USER_* stuff and
>>>>>>>>>>>> setting MPICH_PATH depending on USE_MPISERIAL
  M   ccsm_utils/Machines/Macros.cppdefs
  M   ccsm_utils/Machines/Macros.prototype_hera
  M   ccsm_utils/Machines/Macros.prototype_atlas
  M   ccsm_utils/Machines/Macros.edinburgh_intel
  M   ccsm_utils/Machines/Macros.kraken
  M   ccsm_utils/Machines/Macros.generic_xt
  M   ccsm_utils/Machines/Macros.lynx_pgi
  M   ccsm_utils/Machines/Macros.edinburgh_lahey
  M   ccsm_utils/Machines/Macros.generic_linux_pgi
  M   ccsm_utils/Machines/Macros.bluefire
  M   ccsm_utils/Machines/Macros.franklin
  M   ccsm_utils/Machines/Macros.edinburgh_pgi
  M   ccsm_utils/Machines/Macros.jaguar
  M   ccsm_utils/Machines/Macros.jaguarpf
  M   ccsm_utils/Machines/Macros.brutus_pm
  M   ccsm_utils/Machines/Macros.brutus_po
  M   ccsm_utils/Machines/Macros.brutus_im
  M   ccsm_utils/Machines/Macros.brutus_io
  M   ccsm_utils/Machines/Macros.pingo
  M   ccsm_utils/Machines/Macros.intrepid
  M   ccsm_utils/Machines/Macros.generic_linux_pathscale
  M   ccsm_utils/Machines/Macros.prototype_columbia
  M   ccsm_utils/Machines/Macros.prototype_nyblue
  M   ccsm_utils/Machines/Macros.generic_linux_intel
  M   ccsm_utils/Machines/Macros.hadley
  M   ccsm_utils/Machines/Macros.generic_linux_lahey
  M   ccsm_utils/Machines/Macros.pleiades
  M   ccsm_utils/Machines/Macros.hopper
  M   ccsm_utils/Machines/Macros.prototype_ranger
  M   ccsm_utils/Machines/Macros.prototype_ubgl
  M   ccsm_utils/Machines/Macros.lynx_intel
  M   ccsm_utils/Machines/Macros.generic_ibm
  M   ccsm_utils/Machines/Macros.midnight
  M   ccsm_utils/Machines/Macros.pleiades_wes
  M   ccsm_utils/Machines/Macros.prototype_frost

  M   ccsm_utils/Machines/config_pes.xml ----- Change "--" to "-" so can view the file
         explicitly set "pt1_pt1" case for edinburgh to fix bug 1209
  M   ccsm_utils/Machines/config_machines.xml  Turn MPISERIAL_SUPPORT to TRUE for
         generic machines, edinburgh*, intrepid, jaguar, jaguarpf, and lynx*
         (previously only bluefire was supported)
  M   ccsm_utils/Tools/ccsm_buildexe.csh ----- Document how to build with ESMF and USE_MPISERIAL

>>>>>>>>>>>> Don't set MPICH_PATH if USE_MPISERIAL is set
  M   ccsm_utils/Machines/env_machopts.generic_linux_pgi
  M   ccsm_utils/Machines/env_machopts.edinburgh_pgi
  M   ccsm_utils/Machines/env_machopts.generic_linux_pathscale
  M   ccsm_utils/Machines/env_machopts.generic_linux_intel
  M   ccsm_utils/Machines/env_machopts.edinburgh_intel
  M   ccsm_utils/Machines/env_machopts.generic_linux_lahey
  M   ccsm_utils/Machines/env_machopts.edinburgh_lahey

>>>>>>>>>>>> Change run files so that if USE_MPISERIAL is FALSE
>>>>>>>>>>>> will simply run without a mpirun
  M   ccsm_utils/Machines/mkbatch.prototype_hera
  M   ccsm_utils/Machines/mkbatch.prototype_atlas
  M   ccsm_utils/Machines/mkbatch.intrepid
  M   ccsm_utils/Machines/mkbatch.brutus_im
  M   ccsm_utils/Machines/mkbatch.brutus_io
  M   ccsm_utils/Machines/mkbatch.pingo
  M   ccsm_utils/Machines/mkbatch.pleiades
  M   ccsm_utils/Machines/mkbatch.lynx_pgi
  M   ccsm_utils/Machines/mkbatch.generic_ibm
  M   ccsm_utils/Machines/mkbatch.prototype_ranger
  M   ccsm_utils/Machines/mkbatch.franklin
  M   ccsm_utils/Machines/mkbatch.prototype_columbia
  M   ccsm_utils/Machines/mkbatch.prototype_frost
  M   ccsm_utils/Machines/mkbatch.generic_linux_pgi
  M   ccsm_utils/Machines/mkbatch.jaguarpf
  M   ccsm_utils/Machines/mkbatch.midnight
  M   ccsm_utils/Machines/mkbatch.pleiades_wes
  M   ccsm_utils/Machines/mkbatch.edinburgh_pgi
  M   ccsm_utils/Machines/mkbatch.hadley
  M   ccsm_utils/Machines/mkbatch.kraken
  M   ccsm_utils/Machines/mkbatch.generic_linux_pathscale
  M   ccsm_utils/Machines/mkbatch.generic_linux_intel
  M   ccsm_utils/Machines/mkbatch.generic_xt
  M   ccsm_utils/Machines/mkbatch.hopper
  M   ccsm_utils/Machines/mkbatch.edinburgh_intel
  M   ccsm_utils/Machines/mkbatch.jaguar
  M   ccsm_utils/Machines/mkbatch.generic_linux_lahey
  M   ccsm_utils/Machines/mkbatch.prototype_ubgl
  M   ccsm_utils/Machines/mkbatch.brutus_pm
  M   ccsm_utils/Machines/mkbatch.brutus_po
  M   ccsm_utils/Machines/mkbatch.prototype_nyblue
  M   ccsm_utils/Machines/mkbatch.lynx_intel
  M   ccsm_utils/Machines/mkbatch.edinburgh_lahey

>>>>>>>>>>>> Build with $GMAKE rather than hardcoded gmake command (needed for Darwin)
  M   ccsm_utils/Components/mct.buildlib
  M   ccsm_utils/Components/csm_share.buildlib
  M   ccsm_utils/Components/pio.buildlib

  M   doc/chap6.xml -- Add a little more about MPICH_PATH and NETCDF_PATH.

>>>>>>>>>>>> Get PTCLM working with changes
>>>>>>>>>>>> PTCLM updates from mpiserial branch
>>>>>>>>>>>> Add PTCLM tests for yong (Mac OS-X laptop)
  M   ccsm_utils/Tools/lnd/clm/PTCLM/PTCLM.py
  M   ccsm_utils/Tools/lnd/clm/PTCLM/testcases.csh
  M   ccsm_utils/Tools/lnd/clm/PTCLM/README

================================================================================
Originator: Chris Fischer
Date: Tue Nov 29 2010
Model: scripts
Version: scripts4_101129
One-line: fix logic in create_newcase for IG and BG pes layouts, fix new trigrid def. 

M      ccsm_utils/Case.template/config_grid.xml
M      create_newcase

================================================================================
Originator: Chris Fischer
Date: Tue Nov 23 2010
Model: scripts
Version: scripts4_101123
One-line: fix edinburgh intel/lahey PIO NO_MPI2 build 

M      ccsm_utils/Machines/Macros.edinburgh_intel
M      ccsm_utils/Machines/Macros.edinburgh_lahey
================================================================================
Originator: Chris Fischer
Date: Mon Nov 22 2010
Model: scripts
Version: scripts4_101122
One-line: pes_file bug fix in create_newcase (bug #1238), B55TRWCN pes update, 
          and CME cost update

M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Testcases/config_tests.xml
M      create_newcase

================================================================================
Originator: Jim Edwards
Date: Wed Nov 17 2010
Model: scripts
Version: scripts4_101117
One-line: Force POP default pio_root=0 improve pio namelist handling
	create_newcase
	ccsm_utils/Machines/config_pes.xml
	ccsm_utls/Case.template/config_definition.xml
	

================================================================================
Originator: Jim Edwards
Date: Thu Nov 12 2010
Model: scripts
Version: scripts4_101112
One-line: Fixed incorrect commit in config_machines.xml, move if def
	  for mpiserial to Macros.cppdefs
	  ccsm_utils/Machines/Macros.cppdefs
	  ccsm_utils/Machines/Macros.hadley
	  ccsm_utils/Machines/Macros.bluefire
	  ccsm_utils/Machines/config_machines.xml
	

================================================================================
Originator: Jim Edwards
Date: Thu Nov 11 2010
Model: scripts
Version: scripts4_101111
One-line: Added initial support for configuration dependent defaults for pio namelist
	  Removed -DNEED_MPI_ROOT in edinburgh Macros
	  
M            25495   ccsm_utils/Case.template/config_definition.xml
M            25495   ccsm_utils/Machines/Macros.edinburgh_intel
M            25495   ccsm_utils/Machines/Macros.edinburgh_lahey
M            25495   ccsm_utils/Machines/Macros.edinburgh_pgi
M            25495   ccsm_utils/Machines/config_pes.xml
M            25495   ccsm_utils/Machines/config_machines.xml
M            25495   create_newcase
	

================================================================================

Originator: Jim Edwards
Date: Mon Nov 8 2010
Model: scripts
Version: scripts4_101108
One-line:  -m "moved pio init from components to driver and consolidated namelist"
	
      ccsm_utils/Case.template/config_compsets.xml
      ccsm_utils/Case.template/ConfigCase.pm
      ccsm_utils/Case.template/config_definition.xml
================================================================================
Originator: jet
Date: Mon Nov 1 2010
Model: scripts
Version: scripts4_101102
One-line:  -m "add support for high res (Eul,HOMME) and flow control mods"
	
M      scripts/ccsm_utils/Case.template/config_compsets.xml
M      scripts/ccsm_utils/Case.template/config_grid.xml
M      scripts/ccsm_utils/Machines/env_machopts.jaguar
M      scripts/ccsm_utils/Machines/env_machopts.jaguarpf
M      scripts/ccsm_utils/Machines/config_pes.xml
================================================================================
Originator: tcraig
Date: Mon Nov 1 2010
Model: scripts
Version: scripts4_101101
One-line: add orbital env variables
	
M      ccsm_utils/Case.template/config_definition.xml
================================================================================
Originator: fischer
Date: Wed Oct 29 2010
Model: scripts
Version: scripts4_101029
One-line: B55TRWCN fix 

M      ccsm_utils/Case.template/config_compsets.xml
================================================================================
Originator: fischer
Date: Wed Oct 27 2010
Model: scripts
Version: scripts4_101027
One-line: load pnetcdf/1.1.1 add WACCM tests, update lynx_intel 

M      ccsm_utils/Machines/env_machopts.kraken
M      ccsm_utils/Machines/env_machopts.franklin
M      ccsm_utils/Machines/env_machopts.jaguar
M      ccsm_utils/Machines/env_machopts.lynx_intel
M      ccsm_utils/Machines/env_machopts.jaguarpf
M      ccsm_utils/Testlists/bluefire.pretag

================================================================================
Originator: fischer
Date: Mon Oct 25 2010
Model: scripts
Version: scripts4_101025
One-line: pes_file bug(#1229) fix and WACCM compsets

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Machines/config_pes.xml
M      ChangeLog
M      create_newcase

================================================================================
Originator: tcraig
Date: Fri Oct 2 2010
Model: scripts
Version: scripts4_101022
One-line: Add MPISERIAL/ESMF check and abort

M      ccsm_utils/Tools/ccsm_buildexe.csh
	
================================================================================
Originator: jedwards, fischer
Date: Mon Oct 18 2010
Model: scripts
Version: scripts4_101018
One-line: update pio config args

M      ccsm_utils/Machines/Macros.bluefire
M      ccsm_utils/Machines/Macros.brutus_pm
M      ccsm_utils/Machines/Macros.brutus_po
M      ccsm_utils/Machines/Macros.edinburgh_intel
M      ccsm_utils/Machines/Macros.edinburgh_lahey
M      ccsm_utils/Machines/Macros.edinburgh_pgi
M      ccsm_utils/Machines/Macros.franklin
M      ccsm_utils/Machines/Mcros.generic_linux_pgi
M      ccsm_utils/Machines/Macros.generic_xt
M      ccsm_utils/Machines/Macros.jaguar
M      ccsm_utils/Machines/Macros.jaguarpf
M      ccsm_utils/Machines/Macros.kraken
M      ccsm_utils/Machines/Macros.lynx_pgi

================================================================================
Originator: fischer
Date: Thu Oct 15 2010
Model: scripts
Version: scripts4_101015b
One-line: made changes to create_train and moved all ${COMP}_NCPL to env_conf

M      ccsm_utils/Tools/create_train
M      ccsm_utils/Case.template/config_definition.xml
M      doc/env_conf_list.xml
M      doc/env_run_list.xml

================================================================================
Originator: fischer
Date: Thu Oct 15 2010
Model: scripts
Version: scripts4_101015a
One-line: update doc for OCN_NCPL move from env_run to env_conf

M      doc/env_conf_list.xml
M      doc/env_run_list.xml

================================================================================

Originator: fischer
Date: Thu Oct 15 2010
Model: scripts
Version: scripts4_101015
One-line: Add T31_g37, train script, mods for intrepid, and bug fixes 

A      ccsm_utils/Tools/load.awk
A      ccsm_utils/Tools/create_train
  added train scripts for jaguar

M      ccsm_utils/Tools/configure
  fixed help text for configure -cleanmach
M      ccsm_utils/Case.template/config_compsets.xml
  added support out of the box runs for B_1850_CN at T31_g37
M      ccsm_utils/Case.template/config_definition.xml
  moved OCN_NCPL from env_run to env_conf  (bug #1216)
M      ccsm_utils/Machines/mkbatch.intrepid
  updated queueing for intreped
M      ccsm_utils/Testlists/bluefire.pretag
  added T31_g37 B_1850_CN tests

================================================================================
Originator: fischer
Date: Thu Oct 14 2010
Model: scripts
Version: scripts4_101014
One-line: update tests lists and machine configs

A      ccsm_utils/Machines/Macros.prototype_hera
A      ccsm_utils/Machines/env_machopts.prototype_hera
M      ccsm_utils/Machines/Macros.prototype_atlas
A      ccsm_utils/Machines/mkbatch.prototype_hera
M      ccsm_utils/Machines/env_machopts.kraken
M      ccsm_utils/Machines/mkbatch.prototype_atlas
M      ccsm_utils/Machines/env_machopts.prototype_atlas
M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Machines/config_machines.xml
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/bluefire.pretag
M      ccsm_utils/Testlists/jaguar.posttag
M      ccsm_utils/Testlists/jaguar.pretag
================================================================================
Originator: tcraig
Date: Tue Oct 5 2010
Model: scripts
Version: scripts4_101005
One-line: add COMP_RUN_BARRIER env variable

M      ccsm_utils/Case.template/config_definition.xml
	
================================================================================
Originator: mvr
Date: Wed Sep 29 2010
Model: scripts
Version: scripts4_100929
One-line: removed obsolete archive scripts

D      ccsm_utils/Tools/archiving
D      ccsm_utils/Tools/archiving/st_archive.sh
D      ccsm_utils/Tools/archiving/lt_archive.sh
D      ccsm_utils/Tools/archiving/README

================================================================================
Originator: eaton
Date: Mon Sep 27 2010
Model: scripts
Version: scripts4_100927a
One-line: mods to Macros files to build COSP simulator as part of atm library

. Add rules to build the COSP files.
M       ccsm_utils/Machines/Macros.bluefire
M       ccsm_utils/Machines/Macros.franklin
M       ccsm_utils/Machines/Macros.jaguar
M       ccsm_utils/Machines/Macros.jaguarpf

================================================================================
Originator: fischer
Date: Mon Sep 27 2010
Model: scripts
Version: scripts4_100927
One-line: changed xt-mpt/3.5.0 to xt-mpt/5.0.0   

M      ccsm_utils/Machines/env_machopts.franklin

================================================================================
Originator: mvr
Date: Thu Sep 16 2010
Model: scripts
Version: scripts4_100916
One-line: added the BATCHQUERY setting for jaguarpf - was UNSET

-fix will help prevent multiple long-term archive jobs from running simultaneously

M      ccsm_utils/Machines/config_machines.xml

================================================================================
Originator: fischer
Date: Wed Sep 01 2010
Model: scripts
Version: scripts4_100901a
One-line: added RCP6.0 compset and tests, coupld of minor changes 

M      ccsm_utils/Case.template/config_compsets.xml

>>>>>>>>>>..  added -DnoI8 flag
M      ccsm_utils/Machines/Macros.intrepid

>>>>>>>>> change lynx configuration and modified testlists
M      ccsm_utils/Machines/mkbatch.lynx_pgi
M      ccsm_utils/Machines/mkbatch.lynx_intel
M      ccsm_utils/Machines/config_pes.xml
A      ccsm_utils/Testlists/lynx.pretag
D      ccsm_utils/Testlists/lynx_pgi.pretag
M      ccsm_utils/Testlists/bluefire.pretag
M      ccsm_utils/Testlists/jaguar.pretag

================================================================================
Originator: mvr
Date: Wed Sep 01 2010
Model: scripts
Version: scripts4_100901
One-line: parallelized the long-term archive script

M      ccsm_utils/Tools/ccsm_l_archive.csh
M      ccsm_utils/Machines/mkbatch.pingo
M      ccsm_utils/Machines/mkbatch.pleiades
M      ccsm_utils/Machines/mkbatch.lynx_pgi
M      ccsm_utils/Machines/mkbatch.bluefire
M      ccsm_utils/Machines/mkbatch.franklin
M      ccsm_utils/Machines/mkbatch.prototype_columbia
M      ccsm_utils/Machines/mkbatch.prototype_frost
M      ccsm_utils/Machines/mkbatch.jaguarpf
M      ccsm_utils/Machines/mkbatch.midnight
M      ccsm_utils/Machines/mkbatch.pleiades_wes
M      ccsm_utils/Machines/mkbatch.edinburgh_pgi
M      ccsm_utils/Machines/mkbatch.hadley
M      ccsm_utils/Machines/mkbatch.kraken
M      ccsm_utils/Machines/mkbatch.hopper
M      ccsm_utils/Machines/mkbatch.edinburgh_intel
M      ccsm_utils/Machines/mkbatch.jaguar
M      ccsm_utils/Machines/mkbatch.intrepid
M      ccsm_utils/Machines/mkbatch.lynx_intel
M      ccsm_utils/Machines/mkbatch.edinburgh_lahey

================================================================================
Originator: erik
Date: Mon Aug 30 2010
Model: scripts
Version: scripts4_100830
One-line: Add rcp6.0 to DATM_PRESAERO options, add in PTCLM directory

>>>>>>>>> Add in I_RCP6.0_CN compset for RCP6.0 and rcp6.0 for DATM_PRESAERO
M       ccsm_utils/Case.template/config_compsets.xml ---- Add I_RCP6.0_CN compset
M       ccsm_utils/Case.template/config_definition.xml -- Add rcp6.0 to DATM_PRESAERO

>>>>>>>>> Add in generic component tree for Tools
A       ccsm_utils/Tools/atm
A       ccsm_utils/Tools/ice
A       ccsm_utils/Tools/ocn
A       ccsm_utils/Tools/glc
A       ccsm_utils/Tools/lnd

>>>>>>>>> Add in PTCLM script for clm
A  +    ccsm_utils/Tools/lnd/clm
A  +    ccsm_utils/Tools/lnd/clm/PTCLM
A  +    ccsm_utils/Tools/lnd/clm/PTCLM/PTCLM.py
A  +    ccsm_utils/Tools/lnd/clm/PTCLM/testcases.csh
A  +    ccsm_utils/Tools/lnd/clm/PTCLM/usr_files
A  +    ccsm_utils/Tools/lnd/clm/PTCLM/usr_files/mkgriddata.TEMPLATE
A  +    ccsm_utils/Tools/lnd/clm/PTCLM/usr_files/mkdatadomain.TEMPLATE
A  +    ccsm_utils/Tools/lnd/clm/PTCLM/PTCLM_sitedata
A  +    ccsm_utils/Tools/lnd/clm/PTCLM/PTCLM_sitedata/cnvrt_trnsyrs2_pftdyntxtfile.pl
A  +    ccsm_utils/Tools/lnd/clm/PTCLM/PTCLM_sitedata/Fluxnet-Canada_sitedata.txt
A  +    ccsm_utils/Tools/lnd/clm/PTCLM/PTCLM_sitedata/EXAMPLE_sitedata.txt
A  +    ccsm_utils/Tools/lnd/clm/PTCLM/PTCLM_sitedata/AmeriFlux_sitedata.txt
A  +    ccsm_utils/Tools/lnd/clm/PTCLM/PTCLM_sitedata/Fluxnet-Canada_soildata.txt
A  +    ccsm_utils/Tools/lnd/clm/PTCLM/PTCLM_sitedata/EXAMPLE_soildata.txt
A  +    ccsm_utils/Tools/lnd/clm/PTCLM/PTCLM_sitedata/AmeriFlux_soildata.txt
A  +    ccsm_utils/Tools/lnd/clm/PTCLM/PTCLM_sitedata/US-Ha1_dynpftdata.txt
A  +    ccsm_utils/Tools/lnd/clm/PTCLM/PTCLM_sitedata/Fluxnet-Canada_pftdata.txt
A  +    ccsm_utils/Tools/lnd/clm/PTCLM/PTCLM_sitedata/EXAMPLE_pftdata.txt
A  +    ccsm_utils/Tools/lnd/clm/PTCLM/PTCLM_sitedata/AmeriFlux_pftdata.txt
A  +    ccsm_utils/Tools/lnd/clm/PTCLM/README

================================================================================
Originator: erik
Date: Tues Aug 24 2010
Model: scripts
Version: scripts4_100824
One-line: Fix bug 1172, bad argument for building threaded with Lahey compiler

Remove -fstack_check argument for compiling with Lahey.

M       ccsm_utils/Machines/Macros.generic_linux_lahey
M       ccsm_utils/Machines/Macros.edinburgh_lahey

================================================================================
Originator: fischer
Date: Thurs Aug 20 2010
Model: scripts
Version: scripts4_100820
One-line: set CCSM_CO2_PPMV=367 for 2000 compsets, added edinburgh Intel infiniband support 
          moved doc and comment changes from release tag into beta tag. Fixed generic intel (see bug 1200) 

M      create_test
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Machines/Macros.generic_linux_intel
M      ccsm_utils/Machines/Macros.edinburgh_intel
M      ccsm_utils/Machines/Macros.lynx_intel
M      ccsm_utils/Machines/mkbatch.edinburgh_intel
M      ccsm_utils/Machines/env_machopts.edinburgh_intel
M      ccsm_utils/Machines/config_pes.xml
M      doc/896pe_layout.jpg
M      doc/grids_list.xml
M      doc/ug.xml
M      doc/machines_list.xml
M      doc/env_conf_list.xml
M      doc/compsets_list.xml
M      doc/chap1.xml
M      doc/chap2.xml
M      doc/chap3.xml
M      doc/chap4.xml
M      doc/chap5.xml
M      doc/chap6.xml
M      doc/app1.xml
M      doc/chap7.xml
M      doc/app2.xml
M      doc/glossary.xml
M      doc/chap8.xml
M      doc/app3.xml
M      doc/chap9.xml
M      doc/app4.xml
M      doc/app5.xml
M      doc/app6.xml
M      doc/bookinfo.xml
M      doc/app7.xml
M      ChangeLog
M      create_newcase
M      create_test_suite

================================================================================
Originator: erik
Date: Thurs Aug 12 2010
Model: scripts
Version: scripts4_100812
One-line: Add START_TOD to env_run

M       ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: fischer
Date: Mon Aug 9 2010
Model: scripts
Version: scripts4_100809a
One-line: updated WACCM_1955-2005 and RCP compsets, added edinburgh pgi infiniband support

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Machines/env_machopts.edinburgh_pgi
M      ccsm_utils/Machines/mkbatch.edinburgh_pgi
M      ccsm_utils/Machines/Macros.edinburgh_pgi

================================================================================
Originator: mvr
Date: Mon Aug 9 2010
Model: scripts
Version: scripts4_100809
One-line: backed out one change from previous tag

backed out the following mod after realizing it would break down if a user was 
re-doing a portion of the model run and older restart files were not cleaned 
out of their run directory...
	
- no longer relies on the time listing of a given restart file type to determine 
  the latest file, but instead uses the timestamp included in the filename (see bug 1192)

M      ccsm_utils/Tools/st_archive.sh

================================================================================
Originator: mvr
Date: Fri Aug 6 2010
Model: scripts
Version: scripts4_100806
One-line: fixes/enhancements to archive scripts

- added source delete option (-d) to hsi put command
- no longer relies on the time listing of a given restart file type to determine 
  the latest file, but instead uses the timestamp included in the filename (see bug 1192)
- no longer exiting if a bad tarfile is encountered, will attempt to process what's there
- bug fix for the occasion where new files are moved into the st archive dir while the 
  lt archiver is processing files there

M      ccsm_utils/Tools/ccsm_mswrite
M      ccsm_utils/Tools/st_archive.sh
M      ccsm_utils/Tools/ccsm_l_archive.csh

================================================================================
Originator: fischer
Date: Wed Aug 4 2010
Model: scripts
Version: scripts4_100804
One-line: add RCP tests
 
M      ccsm_utils/Testlists/lynx_pgi.pretag
M      ccsm_utils/Testlists/bluefire.pretag
M      ccsm_utils/Testlists/jaguar.pretag
================================================================================
Originator: fischer
Date: Tue Jul 30 2010
Model: scripts
Version: scripts4_100730
One-line: updated machines and compsets

M      ccsm_utils/Tools/ccsm_msread
M      ccsm_utils/Tools/ccsm_mswrite
M      ccsm_utils/Tools/ccsm_msls
M      ccsm_utils/Tools/taskmaker.pl
M      ccsm_utils/Tools/ccsm_msmkdir
   add support for pleiades

M      ccsm_utils/Case.template/config_compsets.xml
   add WACCM compsets and fixed F_AMIP_CAM5, F_1850_CN_CHEM

D      ccsm_utils/Machines/mkbatch.prototype_schirra
D      ccsm_utils/Machines/Macros.prototype_schirra
D      ccsm_utils/Machines/env_machopts.prototype_schirra
   schirra decommissioned

M      ccsm_utils/Machines/Macros.intrepid
   add CONFIG_ARGS for pio and mct

A      ccsm_utils/Machines/mkbatch.pleiades
A      ccsm_utils/Machines/Macros.pleiades
A      ccsm_utils/Machines/env_machopts.pleiades
A      ccsm_utils/Machines/mkbatch.pleiades_wes
A      ccsm_utils/Machines/Macros.pleiades_wes
A      ccsm_utils/Machines/env_machopts.pleiades_wes
   add support for pleiades

M      ccsm_utils/Machines/config_pes.xml
    add pleiades, remove schirra, update intrepid
M      ccsm_utils/Machines/config_machines.xml
    add pleiades, remove schirra

================================================================================
Originator: fischer
Date: Tue Jul 27 2010
Model: scripts
Version: scripts4_100727
One-line: fix co2_ppmv values for FW and F1850W

M      ccsm_utils/Case.template/config_compsets.xml
================================================================================
Originator: fischer
Date: Wed Jul 21 2010
Model: scripts
Version: scripts4_100721
One-line: added lynx pgi support and hooks for lynx intel support, updated pes.

A      ccsm_utils/Machines/mkbatch.lynx_pgi
A      ccsm_utils/Machines/Macros.lynx_pgi
A      ccsm_utils/Machines/env_machopts.lynx_pgi
A      ccsm_utils/Machines/Macros.lynx_intel
A      ccsm_utils/Machines/env_machopts.lynx_intel
A      ccsm_utils/Machines/mkbatch.lynx_intel
M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Machines/config_machines.xml
A      ccsm_utils/Testlists/lynx_pgi.pretag
================================================================================
Originator: fischer
Date: Thu Jun 17 2010
Model: scripts
Version: scripts4_100617
One-line: fix env_machopts.franklin bug and added CAM5 threaded tests

M      ccsm_utils/Machines/env_machopts.franklin
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/jaguar.posttag

================================================================================
Originator: fischer
Date: Wed Jun 16 2010
Model: scripts
Version: scripts4_100616
One-line: remove P4A tests

M      ccsm_utils/Testlists/bluefire.esmf.auxtest
M      ccsm_utils/Testlists/bluefire.pretag
	
================================================================================
Originator: mvertens
Date: Tue Jun 15 2010
Model: scripts
Version: scripts4_100615
One-line: fixed create_newcase compset listing to account for missing compsets

M      create_newcase
	        
===============================================================================
Originator: tcraig
Date: Mon Jun 14 2010
Model: scripts
Version: scripts4_100614a
One-line: remove documentation of account and queue options from create_test_suite

M      create_test_suite
	        
===============================================================================
Originator: mvertens
Date: Mon Jun 14 2010
Model: scripts
Version: scripts4_100614
One-line: removed tests for B1850WCN at f45_g37

in testlists, removed ERS.f45_g37.B1850WCN.bluefire, added ERS.f45_g37.FW.bluefire
	
M      ccsm_utils/Testlists/bluefire.esmf.auxtest
M      ccsm_utils/Testlists/bluefire.pretag
	
================================================================================
Originator: mvertens
Date: Sat Jun 12 2010
Model: scripts
Version: scripts4_100612
One-line: changed modules for jaguar and franklin (bfb)

M      ccsm_utils/Machines/env_machopts.franklin
M      ccsm_utils/Machines/env_machopts.jaguar
	
================================================================================
Originator: kauff
Date: Thu Jun 10 2010
Model: scripts
Version: scripts4_100610
One-line: add refcase for B_1850_CAM5, all science runs are tested on two machines

* add <compset NAME="B_1850_CAM5" SHORTNAME="B1850C5" GRID_MATCH="f19_g16"
* change SHORTNAME "F2000W" to "FW" for consistency
* add bluefire.science.auxtest ~ has ERS for all validated science runs

M      ccsm_utils/Case.template/config_compsets.xml
A      ccsm_utils/Testlists/bluefire.science.auxtest
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/jaguar.posttag
M      ccsm_utils/Testlists/jaguar.pretag
M      README ~ expand tab to spaces

================================================================================
Originator: mvertens
Date: Wed Jun 9 2010
Model: scripts
Version: scripts4_100609
One-line:  updated compsets and README for new scientifically supported compsets

Also updated create_newcase to have a more condensed listing and list the
valid resolutions if an invalid resolution is entered. The new README has been
updated to reflect all of the new scientifically supported compsets and 
resolutions.	
	
M      ccsm_utils/Case.template/config_compsets.xml
M      create_newcase
M      README
	
================================================================================
Originator: tcraig
Date: Tue Jun  8 2010
Model: scripts
Version: scripts4_100608a
One-line:  fix typos in documentation

M      ccsm_utils/Case.template/config_definition.xml
M      create_newcase
	
================================================================================
Originator: jwolfe
Date: Tue Jun  8 15:36:02 MDT 2010
Model: scripts
Version: scripts4_100608
One-line:  update cism resolution support, fix to INVALID outputs

- additional resolution support for CISM compsets

M      create_newcase
M      README
M      ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: eaton
Date: Mon Jun  7 18:42:02 MDT 2010
Model: scripts
Version: scripts4_100607
One-line:  update cam use case in 20thC BGC compsets

. change cam use case from 1850-2005_cam4 to 1850-2005_cam4_bgc in the 
  B_1850-2000_BGC-* compsets.

M      ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: mvertens
Date: Sun Jun  6 MDT 2010
Model: scripts
Version: scripts4_100606
One-line:  Removed reference to CAM_DTIME from config_grid.xml for T85

M      ccsm_utils/Case.template/config_grid.xml

================================================================================
Originator: fischer
Date: Fri Jun  4 09:15:44 MDT 2010
Model: scripts
Version: scripts4_100604
One-line:  Change queue name on hopper from debug to regular

M      ccsm_utils/Machines/mkbatch.hopper

================================================================================
Originator: jwolfe
Date: Thu Jun  3 16:04:01 MDT 2010
Model: scripts
Version: scripts4_100603a
One-line:  Update for GLC and a bug fix

- remove -C from CICE debug build on bluefire
- add GRID_MATCH logic for active GLC compsets
- add more pre- and post-tag tests using active GLC compsets

M      ccsm_utils/Machines/Macros.bluefire
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/bluefire.pretag
M      ccsm_utils/Testlists/jaguar.pretag

================================================================================
Originator: eaton
Date: Thu Jun  3 13:31:40 MDT 2010
Model: scripts
Version: scripts4_100603
One-line:  Modify BGC compsets to allow CAM to use data CO2 for rad calcs

M      ccsm_utils/Case.template/config_compsets.xml

<compset NAME="B_1850_BGC-BDRD" SHORTNAME="B1850BDRD"
<compset NAME="B_1850-2000_BGC-BDRD" SHORTNAME="B20TRBDRD"
. add CAM_NAMELIST_OPTS="co2_cycle_rad_passive=.true."

================================================================================
Originator: kauff
Date: Tue Jun 1 2010
Model: scripts
Version: scripts4_100601a
One-line:  New BGC compsets & bluefire posttag tests

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Testlists/bluefire.posttag

<compset NAME="B_1850_BGC-BPRP" SHORTNAME="B1850BPRP"
<compset NAME="B_1850_BGC-BDRD" SHORTNAME="B1850BDRD"
<compset NAME="B_1850-2000_BGC-BPRP" SHORTNAME="B20TRBPRP"
<compset NAME="B_1850-2000_BGC-BDRD" SHORTNAME="B20TRBDRD"

================================================================================
Originator: tcraig
Date: Tue Jun 1 2010
Model: scripts
Version: scripts4_100601
One-line:  Add ERU test and update create_production_test

- add ERU test
- modify create_production_test to use ERU test
- update config_pes.xml for a few recent layouts
- clean up sample_compset_file and remove old variables

M      create_clone
M      ccsm_utils/Tools/create_production_test
M      ccsm_utils/Machines/config_pes.xml
A      ccsm_utils/Testcases/ERU_script
M      ccsm_utils/Testcases/config_tests.xml
M      sample_compset_file.xml
	
================================================================================
Originator: tcraig
Date: Thu May 27 2010
Model: scripts
Version: scripts4_100527
One-line:  Change ESMF lib from rp1 to rp2, add esmf auxtest lists
	
M      ccsm_utils/Machines/env_machopts.jaguar
M      ccsm_utils/Machines/env_machopts.jaguarpf
M      ccsm_utils/Machines/config_machines.xml
A      ccsm_utils/Testlists/bluefire.esmf.auxtest
A      ccsm_utils/Testlists/jaguar.esmf.auxtest
A      ccsm_utils/Testlists/jaguarpf.esmf.auxtest
	
================================================================================
Originator: jwolfe
Date: Tue May 25 2010
Model: scripts
Version: scripts4_100525
One-line:  Change compsets to point to cism rather than gglc

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/ConfigCase.pm
M      ccsm_utils/Case.template/config_definition.xml
M      create_newcase

================================================================================
Originator: tcraig
Date: Mon May 24 2010
Model: scripts
Version: scripts4_100524d
One-line:  Add CaseStatus, check_case, and .submit script capabilities

M      create_test
A      ccsm_utils/Tools/check_case
M      ccsm_utils/Tools/testcase_begin
M      ccsm_utils/Tools/configure
M      ccsm_utils/Tools/ccsm_postrun.csh
M      ccsm_utils/Tools/testcase_end
M      ccsm_utils/Tools/generate_batch.csh
M      create_newcase

================================================================================
Originator: erik
Date: Mon May 24 2010
Model: scripts
Version: scripts4_100524b
One-line:  Change clm testlist a bit, update CLM I-case finidat files

 - Update finidat files for I cases
 - Update clm aux test to include ESMF (_E test) and change IG test to f19

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Case.template/config_compsets.xsl - Remove SCIENCE/STATUS from display
M       ccsm_utils/Testlists/bluefire.clm.auxtest

================================================================================
Originator: kauff
Date: Mon May 24 2010
Model: scripts
Version: scripts4_100524
One-line: rm SCIENCE,STATUS in config_compsets.xml, new: README, G intrepid layout, OCN_CO2_* var

create_newcase cat's scripts/README
rm SCIENCE,STATUS from config_compsets.xml
M      ccsm_utils/Case.template/config_compsets.xml
M      create_newcase
M      README

new intrepid layout, use in intrepid.postag
M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Testlists/intrepid.posttag

add an allowable value to OCN_CO2_TYPE (in env_conv.xml)
add a new OCN_ variable, OCN_CO2_FLUX_OCMIP_BUG_FIX, to env_conf.xml

M      ccsm_utils/Case.template/config_definition.xml
M      doc/compsets_list.xml
M      doc/chap3.xml
	
================================================================================

Originator: mvertens
Date: Sun May 23 2010
Model: scripts
Version: scripts4_100523
One-line: Modified xmlchange to check whether id is in requested file (abort if not)

M      ccsm_utils/Tools/xmlchange
	
================================================================================
Originator: mvertens
Date: Wed May 19 2010
Model: scripts
Version: scripts4_100519b
One-line: Fix compsets to use FAMIPCN correctly 

Note that GRID_MATCH must have an original compset defined with the same name

M      ChangeLog
M      create_newcase
	
================================================================================
Originator: fischer
Date: Wed May 19 2010
Model: scripts
Version: scripts4_100519
One-line: Use default compilers on bluefire, and new PES layout for edinburgh

M      ccsm_utils/Machines/Macros.bluefire
M      ccsm_utils/Machines/config_pes.xml
================================================================================

Originator: kauff
Date: Tue May 18 2010
Model: scripts
Version: scripts4_100518
One-line: bug fix for create_newcase

As of scripts4_100513, create_newcase was causing OCN_NX and OCN_NY to be set
to zero in env_case.xml. Only the X compset uses this value. 
	
M      create_newcase
	
================================================================================
Originator: kauff
Date: Mon May 17 2010
Model: scripts
Version: scripts4_100517
One-line: new ne240_f02_g16, f02_g16; rm CICE_PRESAERO_TYPE; more pop2 archiving


M      ccsm_utils/Case.template/config_grid.xml  <-- new hi-res grids, rm g14,g15
M      ccsm_utils/Machines/config_pes.xml        <-- new hi-res layouts
   +<pes grid_match="0.47x0.63_gx1" compset_match="B">
   +<pes grid_match="0.23x0.31_gx1" compset_match="B">
   +<pes grid_match="ne240np4_0.23x0.31_gx1">

M      ccsm_utils/Testlists/jaguar.posttag       <-- new hi-res tests
M      ccsm_utils/Testlists/bluefire.posttag     <-- new hi-res tests

M      ccsm_utils/Case.template/config_definition.xml <-- rm CICE_PRESAERO_TYPE
M      ccsm_utils/Case.template/config_compsets.xml   <-- rm CICE_PRESAERO_TYPE

M      ccsm_utils/Tools/st_archive.sh   <-- recognizes new pop output

================================================================================
Originator: mvertens
Date: Sun May 16 2010
Model: scripts
Version: scripts4_100516
One-line:  Changed all I,D,G compsets so that DATM_PRESAERO is on by default

M      ccsm_utils/Case.template/config_compsets.xml
	
================================================================================
Originator: jwolfe
Date: Thu May 13 2010
Model: scripts
Version: scripts4_100513
One-line: put GLC_NX and GLC_NY definitions back in create_newcase, necessary for
          dead and data models

M      create_newcase

================================================================================
Originator: kauff
Date: Wed May 10 2010
Model: scripts
Version: scripts4_100512
One-line: bug fix: CAM5 testlist names on bluefire, T62_gx3 pe layouts

M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Testlists/bluefire.pretag
M      ccsm_utils/Testlists/bluefire.posttag

================================================================================
Originator: kauff
Date: Mon May 10 2010
Model: scripts
Version: scripts4_100510a
One-line: add "CAM5" compsets, change "track5" to "cam5"

+ add CM5 compsets, WACCM name cleanup BW1850 -> B1850W
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_definition.xml

+ add CM5 tests, WACCM name cleanup BW1850 -> B1850W
M      ccsm_utils/Testlists/bluefire.pretag
M      ccsm_utils/Testlists/jaguar.pretag

M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/jaguar.posttag
M      ccsm_utils/Testlists/jaguarpf.posttag
M      ccsm_utils/Testlists/franklin.posttag
M      ccsm_utils/Testlists/intrepid.posttag
M      ccsm_utils/Testlists/kraken.posttag
M      ccsm_utils/Testlists/edinburgh.posttag

+ change compset name "TR5" to "C5"
M      ccsm_utiles/Testlists/bluefire.cam.auxtest
M      ccsm_utiles/Testlists/bluefire.cice1.auxtest
M      ccsm_utiles/Testlists/bluefire.cice2.auxtest
M      ccsm_utiles/Testlists/jaguar.cam.auxtest

+ WACCM name cleanup BW1850 -> B1850W
M      ccsm_utils/Machines/config_pes.xml

+ add missing endif
M      create_test

================================================================================
Originator: jwolfe
Date: Mon May 10 2010
Model: scripts
Version: scripts4_100510
One-line: Merge scripts trunk with glc branch to bring glimmer support to the trunk

bug fix to create_test, changing an "endif" to an "end"
M      create_test

add cism files to the short-term archiving script
M      ccsm_utils/Tools/st_archive.sh

remove GRID_MATCHes from IG compsets, because the initial conditions files are
incompatible
M      ccsm_utils/Case.template/config_compsets.xml

remove GLC_GRID from grid configurations, moved to env_conf.xml
M      ccsm_utils/Case.template/config_grid.xml
M      ccsm_utils/Case.template/config_definition.xml

add auxilliary testlist for GLC on bluefire
A      ccsm_utils/Testlists/bluefire.glc.auxtest

add "cism" as new optional GLC land ice component
M      doc/env_case_list.xml
M      create_newcase

================================================================================
Originator: erik
Date: Wed May 5 2010
Model: scripts
Version: scripts4_100505b
One-line: Add CLM_BLDNML_OPTS to env_conf, add IG test to bluefire clm test list

M       ccsm_utils/Case.template/config_definition.xsl - Separate into lists
            to make easier to read, add type and default in table as well.
M       ccsm_utils/Case.template/config_compsets.xsl --- Put more details and 
            easier to view
M       ccsm_utils/Case.template/config_definition.xml - Add CLM_BLDNML_OPTS
            so CLM build-namelist commandline options can be set from your case
            without having to edit the template file. Remove extra " that was
            in CCSM_VOC definition.
M       ccsm_utils/Testlists/bluefire.clm.auxtest ------ Add an IG ERI test in

================================================================================
Originator: kauff
Date: Wed May 5 2010
Model: scripts
Version: scripts4_100505
One-line: fix intrepid layout bugs, fix B_2000_CN_CHEM  CCSM_CO2_PPMV="368.9"

fix bug wrt
   B_2000_CN_CHEM  CCSM_CO2_PPMV="368.9", was 379
fix bugs wrt
   <pes grid_match="1.9x2.5_gx1"  compset_match="B" mach_match="intrepid" pecount="L">
   <pes grid_match="0.9x1.25_gx1" compset_match="B" mach_match="intrepid">

M       ccsm_utils/Machines/config_pes.xml
M       ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: kauff
Date: Tue May 4 2010
Model: scripts
Version: scripts4_100504
One-line: fix intrepid layout bug, add intrepid & bluefire layouts

fix bug wrt
   <pes grid_match="0.47x0.63_tx0.1v2" compset_match="B" mach_match="intrepid">
add layouts for
   <pes grid_match="1.9x2.5_gx1"  compset_match="B" mach_match="intrepid" pecount="S">
   <pes grid_match="1.9x2.5_gx1"  compset_match="B" mach_match="intrepid" pecount="M">
   <pes grid_match="1.9x2.5_gx1"  compset_match="B" mach_match="intrepid" pecount="L">
   <pes grid_match="0.9x1.25_gx1" compset_match="B" mach_match="intrepid">
   <pes grid_match="0.9x1.25_gx1" compset_match="B" mach_match="bluefire">
and group layouts by grid_match

M       ccsm_utils/Machines/config_pes.xml

================================================================================
Originator: tcraig
Date: Mon May 3 2010
Model: scripts
Version: scripts4_100503
One-line: fix create_test bug, fix GET_REFCASE message error, fix mkbatch for intrepid

M      create_test
M      ccsm_utils/Tools/ccsm_prestage.csh
M      ccsm_utils/Machines/mkbatch.intrepid

================================================================================
Originator: tcraig
Date: Tue Apr 27 2010
Model: scripts
Version: scripts4_100427a
One-line: add hopper, fix timing resolved RUNDIR issue, remove test report scripts


test reporting mods:
M       create_test
M       create_test_suite

timing RUNDIR mods:
M       ccsm_utils/Tools/generate_batch.csh

hopper mods:
A       ccsm_utils/Machines/Macros.hopper
A       ccsm_utils/Machines/env_machopts.hopper
A       ccsm_utils/Machines/mkbatch.hopper
M       ccsm_utils/Machines/config_machines.xml
M       ccsm_utils/Machines/config_pes.xml
A       ccsm_utils/Testlists/hopper.posttag
	
================================================================================
Originator: kauff
Date: Tue Apr 27 2010
Model: scripts
Version: scripts4_100427
One-line: new WACCM/CHEM compsets, new T62_t12 layouts

- new WACCM/CHEM compsets from fvitt
- new T62_t12 layouts for C,D,G 
- cleanup of f05_t12 B compset layouts

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Machines/config_pes.xml

================================================================================
Originator: tcraig
Date: Fri Apr 23 2010
Model: scripts
Version: scripts4_100423
One-line: add brutus, update cost computation, add MLIBS to build

- add MLIBS to build to support link step order control.  needed
  for frost where certain Macros defined libs needed to appear last 
  to build properly.
- update frost (testing TBD)
- add brutus machine (testing TBD)
- update ranger (not yet working out of the box)
- update cost computation in timing tool, add BATCH_PES env variable.
- fix multiple configure/build bug, #952

M      ccsm_utils/Build/Makefile
M      ccsm_utils/Tools/configure
M      ccsm_utils/Tools/timing/getTiming.pl
M      ccsm_utils/Tools/generate_batch.csh
M      ccsm_utils/Case.template/config_definition.xml
A      ccsm_utils/Machines/mkbatch.brutus_im
A      ccsm_utils/Machines/mkbatch.brutus_io
A      ccsm_utils/Machines/Macros.brutus_im
A      ccsm_utils/Machines/Macros.brutus_io
M      ccsm_utils/Machines/mkbatch.bluefire
A      ccsm_utils/Machines/env_machopts.brutus_pm
A      ccsm_utils/Machines/env_machopts.brutus_po
M      ccsm_utils/Machines/Macros.prototype_ranger
A      ccsm_utils/Machines/env_machopts.brutus_im
A      ccsm_utils/Machines/env_machopts.brutus_io
M      ccsm_utils/Machines/mkbatch.intrepid
A      ccsm_utils/Machines/mkbatch.brutus_pm
A      ccsm_utils/Machines/mkbatch.brutus_po
A      ccsm_utils/Machines/Macros.brutus_pm
M      ccsm_utils/Machines/Macros.prototype_frost
A      ccsm_utils/Machines/Macros.brutus_po
M      ccsm_utils/Machines/config_machines.xml
A      ccsm_utils/Testlists/brutus_im.auxtest
A      ccsm_utils/Testlists/brutus_io.auxtest
A      ccsm_utils/Testlists/brutus_pm.auxtest
A      ccsm_utils/Testlists/brutus_po.auxtest
	
================================================================================
Originator: jwolfe
Date: Wed Apr 21 2010
Model: scripts
Version: scripts4_100421
One-line: fix default T31_gx3v7 pe layout for bluefire

M      ccsm_utils/Machines/config_pes.xml

================================================================================
Originator: kauff
Date: Tue Apr 20 2010
Model: scripts
Version: scripts4_100420a
One-line: add  OCN_TAVG run-time env vars, jaguar.posttag ERS.T62_t12.C.jaguar

M      ccsm_utils/ccsm_utils/ConfigCase.pm
M      ccsm_utils/ccsm_utils/config_definition.xml
M      ccsm_utils/Testlists/jaguarpf.posttag

================================================================================
Originator: mvertens
Date: Tue Apr 20 2010
Model: scripts
Version: scripts4_100420
One-line: added new DATM_PRESAERO variable to env_conf.xml

M      ccsm_utils/Case.template/config_definition.xml

***This tag is needed to work with datm8_100420***

================================================================================
Originator: tcraig
Date: Tue Apr 06 2010
Model: scripts
Version: scripts4_100406a
One-line: Add ARSC XT5 pingo

A      ccsm_utils/Machines/mkbatch.pingo
A      ccsm_utils/Machines/Macros.pingo
A      ccsm_utils/Machines/env_machopts.pingo
M      ccsm_utils/Machines/config_machines.xml
A      ccsm_utils/Testlists/pingo.auxtest

================================================================================
Originator: kauff/erik
Date: Tue Apr 06 2010
Model: scripts
Version: scripts4_100406
One-line: Add CCSM_VOC, rcp4.5, rcp8.5, AMIP compsets, AMIP tests, increase cray walltime 

M       ccsm_utils/Case.template/config_definition.xml - Add CCSM_VOC to env_conf
M       ccsm_utils/Machines/Macros.cppdefs ------------- If CCSM_VOC is TRUE turn
              VOC CPP define on.

M       ccsm_utils/Case.template/config_compsets.xml   <- rcp, AMIP
M       ccsm_utils/Case.template/config_definition.xml <- rcp, VOC
M       ccsm_utils/Machines/Macros.cppdefs    <- VOC
M       ccsm_utils/Machines/mkbatch.franklin  <- walltime increase
M       ccsm_utils/Machines/mkbatch.jaguar    <- walltime increase
M       ccsm_utils/Machines/mkbatch.jaguarpf  <- walltime increase
M       ccsm_utils/Machines/mkbatch.kraken    <- walltime increase
M       ccsm_utils/Testlists/bluefire.posttag <- AMIP
M       ccsm_utils/Tools/ccsm_prestage.csh    <- fix echo message wrt svn usage

================================================================================
Originator: tcraig
Date: Mon Mar 29 2010
Model: scripts
Version: scripts4_100329
One-line: update esmf lib, fix check_input_data

- fixes to check_input_data when running interactively
- updates to the esmf lib to point to esmf 4.0.0rp1_O on bluefire
  and jaguars.

M      ccsm_utils/Tools/check_input_data
M      ccsm_utils/Machines/env_machopts.jaguar
M      ccsm_utils/Machines/env_machopts.jaguarpf
M      ccsm_utils/Machines/config_machines.xml
	
================================================================================
Originator: tcraig
Date: Fri Mar 26 2010
Model: scripts
Version: scripts4_100326
One-line: fixes for cs.status.xml script

- fix output file permissions
- fix to assume . is not in path

M      create_test
M      create_test_suite
	
================================================================================
Originator: erik
Date: Mon Mar 22 22:34:56 MDT 2010
Model: scripts
Version: scripts4_100322b
One-line: Have transient I cases GRAD_MATCH compsets NOT set START_DATE

Start_date for transient 1850-2000 cases was 1-01-01, should be 1850-01-01
Same with IRCP85CN compset while should be 2005-01-01. Should use the generic
value rather than the value set on the GRID_MATCH.

M      ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: jwolfe
Date: Mon Mar 22 2010
Model: scripts
Version: scripts4_100322
One-line: change POP pe count from 2464 to 2356 for X1 high-res configurations
          on the XTs, small Testlist change

M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Testlists/bluefire.pretag
	
================================================================================
Originator: mvertens
Date: Fri Mar 19 2010
Model: scripts
Version: scripts4_100319
One-line: set CO2_PPMV to 368.9 for present day 

M      ccsm_utils/Case.template/config_compsets.xml
       - set CO2 to new value for present day  
	 CCSM_CO2_PPMV="379.000" => CCSM_CO2_PPMV="368.9"
       - set start date for ramp and 20th century for f19_g16 to 501
       - set start date for ramp for f99_g16 to 863
       - the above start dates wil be used in the ccsm4.0 release	 
	
================================================================================
Originator: jwolfe
Date: Thu Mar 18 2010
Model: scripts
Version: scripts4_100318
One-line: set CLM_FORCE_COLDSTART="on" for all GLC compsets

M      ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: kauff
Date: Wed Mar 17 2010
Model: scripts
Version: scripts4_100317
One-line: remove PEM tests, remove compset G_CPLHIST

M      ccsm_utils/Case.template/config_compsets.xml
M      bluefire.posttag
M      bluefire.cice2.auxtes

================================================================================
Originator: kauff
Date: Tue Mar 16 2010
Model: scripts
Version: scripts4_100316
One-line: intrepid compilers NOT hardcoded, changes to intrepid batch submit

Mods come from Ray Loy at ANL

M      ccsm_utils/Machines/Macros.intrepid
M      ccsm_utils/Machines/mkbatch.intrepid

================================================================================
Originator: kauff
Date: Mon Mar 15 2010
Model: scripts
Version: scripts4_100315
One-line: all xt's: don't set MPICH_UNEX_BUFFER_SIZE, -D_USE_FLOW_CONTROL

Also add jaguarpf testlist, initally same as franklin

M      ccsm_utils/Machines/Macros.kraken
M      ccsm_utils/Machines/Macros.generic_xt
M      ccsm_utils/Machines/Macros.franklin
M      ccsm_utils/Machines/Macros.jaguarpf
M      ccsm_utils/Machines/env_machopts.kraken
M      ccsm_utils/Machines/env_machopts.generic_xt
M      ccsm_utils/Machines/env_machopts.franklin
M      ccsm_utils/Machines/env_machopts.jaguar
M      ccsm_utils/Machines/env_machopts.jaguarpf
A      ccsm_utils/Testlists/jaguarpf.posttag

================================================================================
Originator: jwolfe
Date: Fri Mar 12 2010
Model: scripts
Version: scripts4_100312d
One-line: update XT5 modules, tested to be bfb vs previous set on jaguarpf

M      ccsm_utils/Machines/env_machopts.kraken
M      ccsm_utils/Machines/env_machopts.jaguarpf

================================================================================
Originator: kauff
Date: Fri Mar 12 2010
Model: scripts
Version: scripts4_100312c
One-line: remove CIAF test from jaguar.pretag

M      ccsm_utils/Testlists/jaguar.pretag

================================================================================
Originator: fischer
Date: Fri Mar 12 2010
Model: scripts
Version: scripts4_100312b
One-line: fix MPICH buffer error on franklin   

M      ccsm_utils/Machines/env_machopts.franklin
	
================================================================================
Originator: mvr
Date: Fri Mar 12 2010
Model: scripts
Version: scripts4_100312
One-line: fix to check_input_data regarding its interaction with the svn server 

M      ccsm_utils/Tools/check_input_data
	
================================================================================
Originator: tcraig
Date: Thu Mar 11 2010
Model: scripts
Version: scripts4_100311

M      create_test
One line bug fix (addition of endif)
	
================================================================================
Originator: tcraig
Date: Wed Mar 10 2010
Model: scripts
Version: scripts4_100310d
One-line: update esmf.mk handling, 

- remove esmf logic from Macros files, add it to Makefile
- add support for system env variable ESMFMKFILE
- move USE_MPISERIAL and MPISERIAL_SUPPORT from env_build to env_conf
  due to cam configure interaction and fix impacts of this change.
- add module load esmf/4.0.0r_O to jaguar and jaguarpf
- remove local ESMF_LIBDIR path from config_machines.xml for jaguar and jaguarpf
- add some esmf tests to bluefire and jaguar posttag lists
- add two tests to bluefire.pop2.auxtest
- remove esmf.mk existance test in create_newcase
- fix compset shortname name error for F_1850_CN_CHEM

M      ccsm_utils/Build/Makefile
M      ccsm_utils/Tools/generate_batch.csh
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Machines/Macros.intrepid
M      ccsm_utils/Machines/Macros.generic_linux_pathscale
M      ccsm_utils/Machines/Macros.prototype_schirra
M      ccsm_utils/Machines/Macros.prototype_columbia
M      ccsm_utils/Machines/Macros.prototype_nyblue
M      ccsm_utils/Machines/Macros.cppdefs
M      ccsm_utils/Machines/Macros.generic_linux_intel
M      ccsm_utils/Machines/Macros.prototype_atlas
M      ccsm_utils/Machines/Macros.edinburgh_intel
M      ccsm_utils/Machines/Macros.hadley
M      ccsm_utils/Machines/Macros.kraken
M      ccsm_utils/Machines/Macros.generic_linux_lahey
M      ccsm_utils/Machines/Macros.generic_xt
M      ccsm_utils/Machines/Macros.edinburgh_lahey
M      ccsm_utils/Machines/Macros.generic_linux_pgi
M      ccsm_utils/Machines/Macros.bluefire
M      ccsm_utils/Machines/Macros.prototype_ranger
M      ccsm_utils/Machines/Macros.franklin
M      ccsm_utils/Machines/Macros.edinburgh_pgi
M      ccsm_utils/Machines/Macros.jaguar
M      ccsm_utils/Machines/env_machopts.jaguar
M      ccsm_utils/Machines/Macros.prototype_ubgl
M      ccsm_utils/Machines/Macros.jaguarpf
M      ccsm_utils/Machines/env_machopts.jaguarpf
M      ccsm_utils/Machines/Macros.generic_ibm
M      ccsm_utils/Machines/Macros.midnight
M      ccsm_utils/Machines/Macros.prototype_frost
M      ccsm_utils/Machines/config_machines.xml
M      ccsm_utils/Testcases/PEA_script
M      ccsm_utils/Testlists/jaguar.posttag
M      ccsm_utils/Testlists/bluefire.pop2.auxtest
M      ccsm_utils/Testlists/bluefire.posttag
M      create_newcase
	
================================================================================
Originator: erik
Date: Wed Mar 10 15:34:17 MST 2010
Model: scripts
Version: scripts4_100310c
One-line: All I compsets for 1 and 2 degree use a reference case to startup with

>>>>>>>>> Move xlmtestentry to proper place
D      xmltestentry
A  +   ccsm_utils/Tools/xmltestentry

>>>>>>>>> Add GRIDMATCH for 1-deg and 2-deg for all I compsets to use different datasets
M      ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: kauff
Date: Wed Mar 10 13:55:25 MST 2010
Model: scripts
Version: scripts4_100310
One-line: new PET_PT tests, intrepid B layout, B1850CNCHM compset, test_suite results in xml

- sync test lists wrt rel branch, in particular PET_PT tests
- improved intrepid f05_t12 B compset pe layout
- new compset: <compset NAME="F_1850_CN_CHEM" SHORTNAME="B1850CNCHM"
- new support for test_suite results in xml format

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Testlists/bluefire.pretag
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/jaguar.pretag
M      ccsm_utils/Testlists/jaguar.posttag
M      create_test
A      xmltestentry
M      create_test_suite

================================================================================
Originator: erik
Date: Sat Mar 6 2010
Model: scripts
Version: scripts4_100306
One-line: Add stylesheets for compsets/definition, make CO2 in I compsets 368.9, fix name in clm test

M      ccsm_utils/Case.template/config_compsets.xml -- I2000 compsets have CO2 as
                 368.9 rather than 379.000
M      ccsm_utils/Testlists/bluefire.clm.auxtest ----- Fix name of RCP85 compset

A      ccsm_utils/Case.template/config_definition.xsl - Stylesheet to display 
            config_definition.xml file in a browser
A      ccsm_utils/Case.template/config_compsets.xsl --- Stylesheet to display 
            config_compsets.xml file in a browser

================================================================================
Originator: mvr
Date: Fri Mar 5 2010
Model: scripts
Version: scripts4_100305
One-line: fix to archiving scripts to avoid any hsi commands on ncar machines - hsi commands will be added when hpss is fully functional

M      ccsm_utils/Tools/ccsm_msmkdir
	
================================================================================
Originator: mvertens
Date: Thu Mar 4 2010
Model: scripts
Version: scripts4_100304b
One-line: changes for PEA_P1 test to now work

M      ccsm_utils/Tools/configure
M      ccsm_utils/Testcases/PEA_script
	
Note this also needs a change for cam.cpl7.template that should go into
cam4_1_06	
	
================================================================================
Originator: tcraig
Date: Thu Mar 4 2010
Model: scripts
Version: scripts4_100304
One-line: migrate untested machines to prototype_ name convention, delete calgary

- machines migrated to prototype_ are
    atlas, columbia, nyblue, ranger, schirra, ubgl, frost
- delete calgary_lahey and calgary_pgi

D      ccsm_utils/Machines/env_machopts.columbia
D      ccsm_utils/Machines/env_machopts.frost
D      ccsm_utils/Machines/mkbatch.ranger
A  +   ccsm_utils/Machines/mkbatch.prototype_schirra
D      ccsm_utils/Machines/mkbatch.ubgl
A  +   ccsm_utils/Machines/Macros.prototype_schirra
A  +   ccsm_utils/Machines/Macros.prototype_columbia
A  +   ccsm_utils/Machines/env_machopts.prototype_columbia
A  +   ccsm_utils/Machines/Macros.prototype_nyblue
A  +   ccsm_utils/Machines/env_machopts.prototype_nyblue
A  +   ccsm_utils/Machines/mkbatch.prototype_ranger
D      ccsm_utils/Machines/mkbatch.calgary_lahey
A  +   ccsm_utils/Machines/Macros.prototype_atlas
D      ccsm_utils/Machines/Macros.schirra
D      ccsm_utils/Machines/Macros.calgary_lahey
A  +   ccsm_utils/Machines/mkbatch.prototype_columbia
D      ccsm_utils/Machines/mkbatch.calgary_pgi
D      ccsm_utils/Machines/env_machopts.atlas
A  +   ccsm_utils/Machines/mkbatch.prototype_frost
D      ccsm_utils/Machines/mkbatch.frost
D      ccsm_utils/Machines/env_machopts.calgary_pgi
A  +   ccsm_utils/Machines/env_machopts.prototype_frost
D      ccsm_utils/Machines/Macros.ubgl
D      ccsm_utils/Machines/env_machopts.ubgl
D      ccsm_utils/Machines/Macros.frost
D      ccsm_utils/Machines/Macros.nyblue
D      ccsm_utils/Machines/env_machopts.nyblue
A  +   ccsm_utils/Machines/Macros.prototype_ranger
A  +   ccsm_utils/Machines/env_machopts.prototype_ranger
A  +   ccsm_utils/Machines/mkbatch.prototype_atlas
A  +   ccsm_utils/Machines/Macros.prototype_ubgl
D      ccsm_utils/Machines/mkbatch.schirra
D      ccsm_utils/Machines/mkbatch.columbia
A  +   ccsm_utils/Machines/env_machopts.prototype_ubgl
A  +   ccsm_utils/Machines/env_machopts.prototype_schirra
D      ccsm_utils/Machines/mkbatch.nyblue
D      ccsm_utils/Machines/mkbatch.atlas
A  +   ccsm_utils/Machines/env_machopts.prototype_atlas
D      ccsm_utils/Machines/env_machopts.schirra
A  +   ccsm_utils/Machines/mkbatch.prototype_ubgl
D      ccsm_utils/Machines/Macros.atlas
D      ccsm_utils/Machines/env_machopts.calgary_lahey
A  +   ccsm_utils/Machines/mkbatch.prototype_nyblue
D      ccsm_utils/Machines/Macros.ranger
D      ccsm_utils/Machines/env_machopts.ranger
D      ccsm_utils/Machines/Macros.calgary_pgi
A  +   ccsm_utils/Machines/Macros.prototype_frost
M      ccsm_utils/Machines/config_machines.xml
M      ccsm_utils/Machines/config_pes.xml
D      ccsm_utils/Machines/Macros.columbia
	
================================================================================
Originator: kauff
Date: Mon Mar 1 2010
Model: scripts
Version: scripts4_100301b
One-line: tweak bluefire.posttag PET_PT test

M      ccsm_utils/Testlists/bluefire.posttag

================================================================================
Originator: kauff
Date: Mon Mar 1 2010
Model: scripts
Version: scripts4_100301
One-line: add ref case for f09_g16 BRCP85CN, change "cam3_5_1" to "cam4"

also change "track1" & "cam3_5_1" to "cam4" in cam nml & use cases

M      ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: tcraig,kauff
Date: Thu Feb 25, 2010
Model: scripts
Version: scripts4_100225
One-line: test, pe layout, esmf lib updates

- add bluefire PET_PT,PMT B case tests
- new intrepid hires layout
- bug fix of $found in create_newcase
- add esmf lib path to jaguar and jaguarpf

M      create_newcase
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/bluefire.pretag
M      ccsm_utils/Machines/config_pes.xml

================================================================================
Originator: kauff
Date: Wed Feb 24 2010
Model: scripts
Version: scripts4_100224
One-line: add B_RCP8.5_CN compset

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: tcraig
Date: Mon Feb 22 2010
Model: scripts
Version: scripts4_100222
One-line: updates for pio and mct build

- Add Macros CONFIG_SHELL variable for pio and mct builds
- Update CC option in mct build

M      ccsm_utils/Build/Makefile
M      ccsm_utils/Components/pio.buildlib
M      ccsm_utils/Machines/Macros.intrepid
M      ccsm_utils/Machines/Macros.generic_linux_pathscale
M      ccsm_utils/Machines/Macros.generic_linux_intel
M      ccsm_utils/Machines/Macros.schirra
M      ccsm_utils/Machines/Macros.edinburgh_intel
M      ccsm_utils/Machines/Macros.calgary_lahey
M      ccsm_utils/Machines/Macros.hadley
M      ccsm_utils/Machines/Macros.kraken
M      ccsm_utils/Machines/Macros.generic_linux_lahey
M      ccsm_utils/Machines/Macros.generic_xt
M      ccsm_utils/Machines/Macros.edinburgh_lahey
M      ccsm_utils/Machines/Macros.ubgl
M      ccsm_utils/Machines/Macros.generic_linux_pgi
M      ccsm_utils/Machines/Macros.frost
M      ccsm_utils/Machines/Macros.bluefire
M      ccsm_utils/Machines/Macros.nyblue
M      ccsm_utils/Machines/Macros.franklin
M      ccsm_utils/Machines/Macros.edinburgh_pgi
M      ccsm_utils/Machines/Macros.jaguar
M      ccsm_utils/Machines/Macros.jaguarpf
M      ccsm_utils/Machines/Macros.generic_ibm
M      ccsm_utils/Machines/Macros.midnight
M      ccsm_utils/Machines/Macros.atlas
M      ccsm_utils/Machines/Macros.ranger
M      ccsm_utils/Machines/Macros.calgary_pgi
M      ccsm_utils/Machines/Macros.columbia

================================================================================
Originator: tcraig
Date: Fri Feb 19 2010
Model: scripts
Version: scripts4_100219b
One-line: update cam nml use cases for consistency with cam4_1_02
	
M      ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: tcraig
Date: Fri Feb 19 2010
Model: scripts
Version: scripts4_100219a
One-line: add PST, PMT, and PET tests

M      create_test
M      ccsm_utils/Machines/config_pes.xml
A      ccsm_utils/Testcases/PMT_script
M      ccsm_utils/Testcases/config_tests.xml
A      ccsm_utils/Testcases/PST_script
A      ccsm_utils/Testcases/PMT_auto_pes_file
A      ccsm_utils/Testcases/PST_auto_pes_file
M      ccsm_utils/Testcases/PET_script
M      ccsm_utils/Testlists/jaguar.clm.auxtest
M      ccsm_utils/Testlists/bluefire.clm.auxtest
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/bluefire.cice2.auxtest
M      ccsm_utils/Testlists/bluefire.pretag
================================================================================
Originator: kauff
Date: Fri Feb 19 14:21:34 MST 201
Model: scripts
Version: scripts4_100219
One-line: bluefire xlf compiler patch, rm dublin, mod edinburgh pe layout & tests

M      ccsm_utils/Machines/Macros.bluefire ~ xlf compiler patch
M      ccsm_utils/Testlists/edinburgh.pretag  ~ move some posttag test here
M      ccsm_utils/Testlists/edinburgh.posttag ~ replace f45_g37.B20TRCN with f19_g16.B20TRCN
M      Machines/config_machines.xml
M      Machines/config_pes.xml

D      Machines/Macros.dublin_lahey ~ dublin is gone baby gone
D      Machines/env_machopts.dublin_lahey
D      Machines/Macros.dublin_pgi
D      Machines/env_machopts.dublin_pgi
D      Machines/mkbatch.dublin_pgi
D      Machines/mkbatch.dublin_intel
D      Machines/mkbatch.dublin_lahey
D      Machines/Macros.dublin_intel
D      Machines/env_machopts.dublin_intel
D      Testlists/dublin.posttag

================================================================================
Originator: tcraig
Date: Fri Feb 19 2010
Model: scripts
Version: scripts4_100218
One-line: scripts clean up, fix PET

 - update clean_build script to remove source file for pop
 - update ATM_NCPL and OCN_NCPL defaults from 0
 - remove F and I error message for samegrid
 - add copy of seq_maps.rc to CaseDocs
 - change PET to PST and PMT tests

M      create_test
M      ccsm_utils/Tools/clean_build
M      ccsm_utils/Tools/ccsm_buildnml.csh
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Testcases/config_tests.xml
M      ccsm_utils/Testcases/CME_script
M      ccsm_utils/Testlists/jaguar.clm.auxtest
M      ccsm_utils/Testlists/bluefire.clm.auxtest
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/bluefire.cice2.auxtest
M      ccsm_utils/Testlists/bluefire.pretag
M      create_newcase
	
================================================================================
Originator: tcraig
Date: Mon Feb 15 2010
Model: scripts
Version: scripts4_100216
One-line: updates for jaguar

- change gmake -j to 4 for jaguarpf (was 12)
- add -D_USE_FLOW_CONTROL to jaguar and jaguarpf defaults

M      ccsm_utils/Machines/Macros.jaguar
M      ccsm_utils/Machines/Macros.jaguarpf
M      ccsm_utils/Machines/config_machines.xml
	
================================================================================
Originator: tcraig
Date: Mon Feb 15 2010
Model: scripts
Version: scripts4_100215
One-line: scripts updates

 - set default pes_per_node to max_tasks_per_node, remove pes_per_node from create_newcase
   and from most machines.
 - remove HIST_64BIT env variable
 - remove compset FH
 - correct nx and ny values for homme grids in config_grid.xml
 - add descriptions where they were missing in config_grid.xml
 - add NPFIX variable attribute to homme grids in config_grid.xml
 - fix grid name for T85_T85
 - modify ConfigCase.pm for docbook table generation
 - update descriptions for env variables in config_definition.xml
 - update PNETCDF documentation in generic machines
 - update cost of PET and PEM tests
	
M      ccsm_utils/Tools/testcase_setup.csh
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_grid.xml
M      ccsm_utils/Case.template/ConfigCase.pm
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Machines/Macros.generic_linux_pathscale
M      ccsm_utils/Machines/Macros.generic_linux_intel
M      ccsm_utils/Machines/Macros.generic_linux_lahey
M      ccsm_utils/Machines/Macros.generic_xt
M      ccsm_utils/Machines/Macros.generic_linux_pgi
M      ccsm_utils/Machines/Macros.generic_ibm
M      ccsm_utils/Machines/config_machines.xml
M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Testcases/config_tests.xml
M      ccsm_utils/Testlists/jaguar.posttag
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/bluefire.cice2.auxtest
M      doc/env_build_list.xml
M      doc/env_mach_pes_list.xml
M      doc/env_run_list.xml
M      doc/env_case_list.xml
M      create_newcase

================================================================================
Originator: kauff
Date: 
Model: scripts
Version: scripts4_1002??
One-line: accumulating changes for next tag

- rm run_test_suite - it doesn't do anything
D      run_test_suite

================================================================================
Originator: erik
Date: Tue Feb 9, 2010
Model: scripts
Version: scripts4_100209
One-line: Make sure Fortran mangle definitions are in ALL Macros files

This fixes bug 1091

>>>>>> Add FORTRAN_SAME definition to IBM machines that didn't define it.
M      ccsm_utils/Machines/Macros.schirra
M      ccsm_utils/Machines/Macros.ubgl
M      ccsm_utils/Machines/Macros.nyblue

================================================================================
Originator: erik
Date: Tue Feb 4, 2010
Model: scripts
Version: scripts4_100204
One-line: Remove CLM_DEMAND, add in CLM_CO2_TYPE and CLM_NML_USE_CASE

This change requires clm3_7_02 to work.

Change several CLM_ env conf variables, and add a new "I" compset for rcp-8.5.
Add new clm test list for edinburgh, update documentation with these changes.

>>>>> Remove CLM_DEMAND, add in CLM_CO2_TYPE and CLM_NML_USE_CASE
M      ccsm_utils/Case.template/config_compsets.xml ----- Add new compset for I8521CNR85 
                                                          (1850-2200, rcp=8.5 with CN, starts in 2005)
M      ccsm_utils/Case.template/config_definition.xml --- Add new vars delete CLM_DEMAND
M      ccsm_utils/Testlists/bluefire.clm.auxtest -------- Add a test for new compset: I8521CNR85
A      ccsm_utils/Testlists/edinburgh.clm.auxtest ------- Add new clm test list for edinburgh
>>>>> Update documentation with above changes
M      doc/env_conf_list.xml -- Updated from running write_docbook_file
M      doc/compsets_list.xml -- Updated from running write_docbook_file
M      doc/chap3.xml ---------- Document new variables, remove CLM_DEMAND

================================================================================
Originator: kauff
Date: Tue Feb 2, 2010
Model: scripts
Version: scripts4_100202
One-line: sync config_pe.xml & test lists wrt new tri-grid naming convention

- sync wrt config_grid.xml:<horiz_grid GRID="ne30np4_1.9x2.5_gx1v6" SHORTNAME="ne30_f19_g16"

M    ccsm_utils/Machines/config_pes.xml
M    ccsm_utils/Testlists/jaguar.posttag
M    ccsm_utils/Testlists/bluefire.posttag

================================================================================
Originator: tcraig
Date: Tue Feb 2, 2010
Model: scripts
Version: scripts4_100201d
One-line: fixed ERI problem introduced with move of GET_REFCASE to new env file
	
M      ccsm_utils/Testcases/ERI_script
	
================================================================================
Originator: mvertens
Date: Mon Feb 1, 2010
Model: scripts
Version: scripts4_100201c
One-line: fixed problem in ConfigCase.pm
	
M      ccsm_utils/Case.template/ConfigCase.pm
	
================================================================================
Originator: mvertens
Date: Mon Feb 1, 2010
Model: scripts
Version: scripts4_100201b
One-line: changes needed for clearer display of tables in documentation

M      ccsm_utils/Case.template/ConfigCase.pm
M      doc/machines_list.xml
M      doc/compsets_list.xml
M      doc/write_docbook_file
M      doc/env_run_list.xml
M      create_newcase

================================================================================
Originator: tcraig
Date: Mon Feb 1, 2010
Model: scripts
Version: scripts4_100201a
One-line: fix prestage bug introduced in 100201
	
M      ccsm_utils/Tools/ccsm_prestage.csh

================================================================================
Originator: kauff, tcraig
Date: Mon Feb 1, 2010
Model: scripts
Version: scripts4_100201
One-line: update env vars, add interpid pe layout and other pecounts

- remove ATM_CDF64 LND_CDF64 OCN_CDF64 ICE_CDF64 CPL_CDF64 GLC_CDF64
- clean up grid naming convention for tri-grids...
- new pe layout for: <pes grid_match="0.47x0.63_tx0.1v2" compset_match="B" mach_match="intrepid">
- add support for -h option in create_test
- move GET_REFCASE to env_conf
- move PIO_CONFIG_ARGS to env_build
- add BRNCH_RETAIN_CASENAME to env_conf
- add more pecount options to config_pes.xml including additional threaded layouts
  and support for exactly 1,16,32,64,128,4x4,8x4,16x4 layouts
- turn off automatic export of REFCASE data if it's not available locally
- add module load subversion to jaguars, kraken, franklin
- change "dublin" to "edinburgh" in test list

M      create_test
M      ccsm_utils/Machines/env_machopts.kraken
M      ccsm_utils/Machines/env_machopts.franklin
M      ccsm_utils/Machines/env_machopts.jaguar
M      ccsm_utils/Machines/env_machopts.jaguarpf
M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Tools/ccsm_prestage.csh
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Case.template/config_grid.xml
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Testlists/edinburgh.posttag

================================================================================
Originator: mvr
Date: Jan 27 2010
Model: scripts
Version: scripts4_100127
One-line: change to have lt archiver run as size=0 jobs on jaguar, jaguarpf

M      ccsm_utils/Machines/mkbatch.jaguarpf
M      ccsm_utils/Machines/mkbatch.jaguar

================================================================================
Originator: fischer
Date: Jan 25 2010
Model: scripts
Version: scripts4_100125
One-line: added edinburgh, fixed kraken test failures

- added edinburgh to Machines and Testlists
- commented out setenv MPICH_UNEX_BUFFER_SIZE  .. in
  env_machopts.kraken 

A      ccsm_utils/Machines/env_machopts.edinburgh_pgi
A      ccsm_utils/Machines/Macros.edinburgh_intel
A      ccsm_utils/Machines/Macros.edinburgh_lahey
A      ccsm_utils/Machines/mkbatch.edinburgh_pgi
A      ccsm_utils/Machines/Macros.edinburgh_pgi
A      ccsm_utils/Machines/mkbatch.edinburgh_intel
A      ccsm_utils/Machines/env_machopts.edinburgh_intel
A      ccsm_utils/Machines/mkbatch.edinburgh_lahey
A      ccsm_utils/Machines/env_machopts.edinburgh_lahey
M      ccsm_utils/Machines/env_machopts.kraken
M      ccsm_utils/Machines/config_machines.xml
A      ccsm_utils/Testlists/edinburgh.posttag

================================================================================
Originator: tcraig
Date: Jan 22 2010
Model: scripts
Version: scripts4_100122
One-line: update prestaging and tests

 - update configuration warning message
 - add ERI test
 - update PEA test to be only a 1 pe test with MPISERIAL and MPI
 - update listfilesin_streams and check_input_data for better support for
   DIN_LOC_ROOT and DIN_LOC_ROOT_QIAN
 - update prestage script for addition checks, automatic export, and prestaging
   of REFCASE data.
 - remove PRESTAGE_DATA and DIN_LOC_REFCASE env variables.  leverage existing
   env variables for feature support.  remove CCSM_PECOUNT references
 - update ESMFLIB to official v4 release on bluefire
 - improve CME test to support rerun capability
 - copy listfilesin_streams to CASE Tools directory for future use
 - update test lists, add ERI, update PEA

M      ccsm_utils/Tools/check_input_data
M      ccsm_utils/Tools/listfilesin_streams
M      ccsm_utils/Tools/ccsm_prestage.csh
M      ccsm_utils/Tools/ccsm_getenv
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Machines/config_machines.xml
M      ccsm_utils/Testcases/config_tests.xml
M      ccsm_utils/Testcases/CME_script
M      ccsm_utils/Testcases/PEA_script
M      ccsm_utils/Testcases/ERB_script
M      ccsm_utils/Testcases/ERH_script
A      ccsm_utils/Testcases/ERI_script
M      ccsm_utils/Testcases/P4A_script
M      ccsm_utils/Testlists/jaguar.pretag
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/bluefire.pretag
M      create_newcase

================================================================================
Originator: tcraig
Date: 1/14/10
Model: scripts
Version: scripts4_100114
One-line: updates to generic machines

- add inline documentation with tag GENERIC_USER
- add more examples of possible setups in generic machine scripts

M      ccsm_utils/Machines/env_machopts.generic_linux_pgi
M      ccsm_utils/Machines/Macros.generic_linux_pathscale
M      ccsm_utils/Machines/Macros.generic_linux_intel
M      ccsm_utils/Machines/mkbatch.generic_ibm
M      ccsm_utils/Machines/env_machopts.generic_ibm
M      ccsm_utils/Machines/Macros.generic_linux_lahey
M      ccsm_utils/Machines/mkbatch.generic_linux_pgi
M      ccsm_utils/Machines/Macros.generic_xt
M      ccsm_utils/Machines/env_machopts.generic_xt
M      ccsm_utils/Machines/Macros.generic_linux_pgi
M      ccsm_utils/Machines/mkbatch.generic_linux_pathscale
M      ccsm_utils/Machines/mkbatch.generic_linux_intel
M      ccsm_utils/Machines/mkbatch.generic_xt
M      ccsm_utils/Machines/env_machopts.generic_linux_pathscale
M      ccsm_utils/Machines/env_machopts.generic_linux_intel
M      ccsm_utils/Machines/mkbatch.generic_linux_lahey
M      ccsm_utils/Machines/Macros.generic_ibm
M      ccsm_utils/Machines/env_machopts.generic_linux_lahey

================================================================================
Originator: mvertens
Date: 2010 Jan 13
Model: scripts
Version: scripts4_100113d 
One-line: fixes for documentation
Originator:  kauff

Removed reduntant listing of compsets and removed unsupported machines 
in favor of supported machines that are either routinely or not routinely tested
	
M      ccsm_utils/Machines/config_machines.xml
M      create_newcase
	
================================================================================
Originator: kauff
Date: 2010 Jan 13
Model: scripts
Version: scripts4_100113c 
One-line: fix executable file permission in Machines
Originator:  kauff

svn propset svn:executable ON
M      ccsm_utils/Machines/mkbatch.generic_linux_intel
M      ccsm_utils/Machines/mkbatch.generic_linux_lahey
M      ccsm_utils/Machines/env_machopts.generic_linux_pathscale
M      ccsm_utils/Machines/env_machopts.generic_linux_lahey

================================================================================
Date: 2010 Jan 13
Model: scripts
Version: scripts4_100113b 
One-line: add export to suite cs.submit script.  update test lists, hadley, midnight.

- changes to B20TRCNCHM B1850CNCHM coverage on bluefire, franklin
- add check_input_data export call to cs.submit script for suites.  this
  will help automate downloading of new data when testing.
- update hadley job times
- update midnight env to use pgi9.0.4 instead of pgi7.2.2
- add default pes for T62_gx1 cases.

M      create_test
M      ccsm_utils/Machines/mkbatch.hadley
M      ccsm_utils/Machines/env_machopts.midnight
M      ccsm_utils/Machines/config_pes.xml
M      Testlists/bluefire.posttag
M      Testlists/franklin.posttag

================================================================================
Originator: jwolfes
Date: 2010 Jan 13
Model: scripts
Version: scripts4_100113
One-line: remove obsolete setting from kraken batch script

M      ccsm_utils/Machines/mkbatch.kraken

================================================================================

Originator: kauff, mvertens
Date: 2010 Jan 12
Model: scripts
Version: scripts4_100112c
One-line: add B tests, add B20TRCNCHM , error check create_newcase -pecount option

error check valid pecount options
M      create_newcase

Add new B1850CNCHM compsets
M      ccsm_utils/Case.template/config_compsets.xml

Add new B, B1850CNCHM tests
M      Testlists/hadley.auxtest
M      Testlists/intrepid.posttag
M      Testlists/kraken.posttag
M      Testlists/jauguar.posttag
M      Testlists/bluefire.posttag
M      Testlists/bluefire.cam.auxtest

================================================================================
Originator: jwolfes
Date: 2010 Jan 12
Model: scripts
Version: scripts4_100112b
One-line: fix PIO build for XT's, update kraken machine files
	
remove FFLAGS_PIO variable, remove pnetcdf from PIO build on XTs
M      ccsm_utils/Machines/Macros.kraken
M      ccsm_utils/Machines/Macros.generic_xt
M      ccsm_utils/Machines/Macros.franklin
M      ccsm_utils/Machines/Macros.jaguar
M      ccsm_utils/Machines/Macros.jaguarpf

change kraken to 12 pes per node, use same default pe configurations as jaguarpf
M      ccsm_utils/Machines/config_machines.xml
M      ccsm_utils/Machines/config_pes.xml
	
================================================================================

Originator: mvertens
Date: 2010 Jan 12
Model: scripts
Version: scripts4_100112
One-line: fixed ERB, ERH and P4A tests that were effected by DIN_LOC_REFCASE 
	
fixed hadley DIN_LOC_ROOT_CLMQIAN
M      ccsm_utils/Machines/config_machines.xml
	
fixed problems encountered with running tests that used DIN_LOC_REFCASE with
a non-blank value
M      ccsm_utils/Testcases/ERB_script
M      ccsm_utils/Testcases/ERH_script
M      ccsm_utils/Testcases/P4A_script
	
================================================================================
Originator: mvr
Date: 2010 Jan 11
Model: scripts
Version: scripts4_100111b
One-line: new perl5lib tag - new filenames where a ':' was in the name - needed 
on platforms (windows) that can not handle such filenames
	
 M     .
M      SVN_EXTERNAL_DIRECTORIES

================================================================================
Originator: kauff, tcraig, jedwards
Date: 2010 Jan 11
Model: scripts
Version: scripts4_100111
One-line: Testlist changes, mkbatch.generic_ibm mods, DIN_LOC_REFCASE mods
	
Replace B compsets with B1850CN (B is non-standard)
M      Testlists/bluefire.pop2.auxtest
M      Testlists/franklin.posttag
M      Testlists/hadley.auxtest
M      Testlists/intrepid.posttag
M      Testlists/kraken.posttag

updates to loadleveler configuration
M      ccsm_utils/Machines/mkbatch.generic_ibm

- move DIN_LOC_REFCASE checks to ccsm_prestage to be consistent
  with env_run variable and so checks are done for every run.
- add logic to turn on DIN_LOC_REFCASE checks only when CONTINUE_RUN is FALSE.

M      ccsm_utils/Tools/configure
M      ccsm_utils/Tools/ccsm_buildexe.csh
M      ccsm_utils/Tools/ccsm_prestage.csh

================================================================================
Originator: mvertens
Date: 1/10/2010
Model: scripts
Version: scripts4_100110
One-line: Put in capability for data to be prestaged from ccsm4_init

Introduced new variable in env_run.xml, DIN_LOC_REFCASE This is the
path relative to DIN_LOC_ROOT where standard REFCASE data is kept
(e.g. for B_1850_CN at f09_g16). This variable is automiatically
filled in by create_newcase.

$CASE.$MACH.build will automatically bring this data over into
$RUNDIR. configure will check if it is in prompt the user for
obtaining it form the inputdata repository.
	
M      ccsm_utils/Tools/configure
M      ccsm_utils/Tools/ccsm_buildexe.csh
M      ccsm_utils/Tools/ccsm_prestage.csh
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_definition.xml
	
D      ccsm_utils/Case.template/config_definition.xsl
D      ccsm_utils/Case.template/config_compsets.xsl

================================================================================
Originator: mvertens
Date: 1/8/2010
Model: scripts
Version: scripts4_100108b
One-line: Changed pnetcdf paths for bluefire and frost

These changes are needed to use the new pio code which is now in 
http://parallelio.googlecode.com/svn/trunk_tags/pio1_0_6

M      ccsm_utils/Machines/Macros.bluefire
M      ccsm_utils/Machines/Macros.frost
	
================================================================================
Originator: mvertens,kauff
Date: 1/8/2010
Model: scripts
Version: scripts4_100108
One-line: Put in new comments for science support or lack thereof 
	
Also remove more g35,g14,g15 grid combos (obsolete pop grids)
Change compset name: G1850 -> G1850ECO (uses pop ecosystem)

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/ConfigCase.pm
M      create_newcase
	
================================================================================
Originator: jwolfe
Date: 1/7/2010
Model: scripts
Version: scripts4_100107b
One-line: Update machine files to avoid communication problems

- comment out setting MP_EAGER_LIMIT on bluefire, since the default setting changes
  with processor count and we have had to set this to 0 to run consistently
- add FFLAGS_PIO to all XT Macros files, to avoid a PIO problem that occurs writing
  files when one or more components uses OpenMP threads.  Right now, these flags
  build PIO in debug mode, which has been tested to work on jaguarpf.

M      ccsm_utils/Machines/env_machopts.bluefire
M      ccsm_utils/Machines/Macros.franklin
M      ccsm_utils/Machines/Macros.generic_xt
M      ccsm_utils/Machines/Macros.jaguar
M      ccsm_utils/Machines/Macros.jaguarpf
M      ccsm_utils/Machines/Macros.kraken

================================================================================
Originator: erik
Date: 1/7/2010
Model: scripts
Version: scripts4_100107
One-line: Remove f45_g35 tests from test lists

Remove all of the f45_g35 and replace with f45_g37 in the test lists, as f45_g35
has been removed. Also remove DATM_CO2_TSERIES from config_compsets as was removed
from definition file. And fix I1850CN compset name in bluefire.clm.auxtest test list.

M       ccsm_utils/Case.template/config_compsets.xml
M       ccsm_utils/Testlists/bluefire.cam.auxtest
M       ccsm_utils/Testlists/midnight.auxtest
M       ccsm_utils/Testlists/atlas.auxtest
M       ccsm_utils/Testlists/hadley.auxtest
M       ccsm_utils/Testlists/bluefire.clm.auxtest
M       ccsm_utils/Testlists/bluefire.cice1.auxtest
M       ccsm_utils/Testlists/bluefire.cice2.auxtest


================================================================================
Originator: mvertens
Date: 1/6/2010
Model: scripts
Version: scripts4_100106c
One-line: atm/ocn grids allowed to be identical only for F compsets

M      create_newcase
	
================================================================================
Originator: mvertens
Date: 
Model: scripts
Version: scripts4_100106b
One-line: put in grid optional attributes into compsets

- special input data handling for 1850 control for f09_g16 and f19_g16 starts 
and 20th century f09_g16 - these now start hybrid runs from the 
supported spun-up data in $DIN_LOC_ROOT_CSMDATA/ccsm4_init

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/ConfigCase.pm
M      create_newcase

================================================================================
Originator: kauff, jacob
Date: Wed Jan  6 16:59:25 MST 2010
Model: scripts
Version: scripts4_100106
One-line: fix casename typo, new homme mapping files, intrepid mods

fix typo in B1850RMCN compset name:
M      ccsm_utils/Case.template/config_compsets.xml

new homme mapping files, pe layout, add homme posttag test
M      ccsm_utils/Testlists/jaguar.posttag
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Case.template/config_grid.xml
M      ccsm_utils/Machines/config_pes.xml

intrepid changes from ANL:
M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Machines/mkbatch.intrepid*
M      ccsm_utils/Machines/Macros.intrepid
M      ccsm_utils/Machines/config_machines.xml

================================================================================
Originator: kauff, erik, jwolfe, mvr
Date: Tue Jan  5 16:30:24 MST 2010
Model: scripts
Version: scripts4_100105b
One-line: compset name cleanup; modify handling of initial files in archive scripts

config_compsets.xml cleanup...
   Capitalize compset descriptions.
   Remove "TRACK1"/TR1 out of compset names/shortnames.
   Remove DEV out of TRACK5DEV in CAM compset names/shortnames.
   Remove compsets with CAM that don't have track1/track5, or convert to track1
   Remove WACCM 1995 compsets (B & F).
   Improve descriptions for I compsets.
      I_1850 compsets are now over 1948-1972 datm forcing.
      I_GLC => I_2000_GLC, also set CO2 to 379
      I_1948-2004* also set CO2 to 379
      I_1850-2000* set CO2 to 379, datm over 1948-1972, rather than 1948-2000 aligned on 1895
      I_CN  -> I_2000_CN
   Add transient F case with CN "F_1850-2000_CN".
   Add stylesheets so can display the config_compsets and config_definition
      xml files in a browser (can also help find xml formatting issues).

Change comments in config_definition.xml so can be displayed.
Fix a couple xml errors in config_definition.xml.
Remove DATM_CO2_TSERIES as not used in datm yet.

Remove "-d" flag from POP build on bluefire, as requested by Nancy Norton
Add default configurations for T31_gx3v7 B on bluefire for S, M, and L pe-counts
Add default configuration for 1.9x2.5_gx1v6 B on jaguarpf

initial files will now be treated like restart files during archiving - latest 
files will be kept as part of restart bundle and interim files removed, unless 
specifically requested to be saved via DOUT_S_SAVE_INT_REST_FILES

fix dublin tools path: CCSM_CPRNC="/fs/cgd/csm/tools/cprnc_64/cprnc

M      ccsm_utils/Case.template/config_compsets.xml
A      ccsm_utils/Case.template/config_definition.xsl
A      ccsm_utils/Case.template/config_compsets.xsl
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Testlists/*
M      ccsm_utils/Machines/Macros.bluefire
M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Machines/config_machines.xml
M      ccsm_utils/Tools/st_archive.sh


Update HOMME grids
   - rm     <horiz_grid GRID="ne16np4_T42_gx3v5" SHORTNAME="ne16T42_g35"
   - rm     <horiz_grid GRID="ne30np4_T42_gx3v5" SHORTNAME="ne30T42_g35"
   - add    <horiz_grid GRID="ne30np4_fv1.9x2.5_gx1v6" SHORTNAME="ne30f19_g16"
   - change <horiz_grid GRID="ne240np4_fv23x31_gx1v5" SHORTNAME="ne240f23_g15"

M      ccsm_utils/Case.template/config_grids.xml

================================================================================
Originator: mvertens
Date: Mon Jan 3 2010
Model: scripts
Version: scripts4_100103b
One-line: more updates for User's Guide and release

D      switch_machine
M      README
D      README_production
M      ccsm_utils/Case.template/ConfigCase.pm
M      ccsm_utils/Case.template/config_definition.xml
M      doc/env_build_list.xml
M      doc/env_mach_pes_list.xml
M      doc/env_conf_list.xml
M      doc/env_run_list.xml
D      README_quickstart

================================================================================
Originator: mvertens
Date: Mon Jan 3 2010
Model: scripts
Version: scripts4_100103
One-line: update documentation (for docbook) and config_definition.xml

M      ccsm_utils/Case.template/ConfigCase.pm
       - added capability to write docbook tables for all the xml files
M      ccsm_utils/Case.template/config_definition.xml
       - changed variables CCSM_SSTDATA_ to DOCN_SSTDATA_   
       - added a lot more long description
M      ccsm_utils/Machines/config_pes.xml
M      create_newcase
D      README_bgc

       Note - the doc/ directory is completely different. All the original 
       documentation was removed and the docbook files and xml tables were
       added	
	
================================================================================
Originator: kauff/erik
Date: Thu Dec 31 15:14:00 MST 2009
Model: scripts
Version: scripts4_091231
One-line: rm old grids, B,I,F compset cleanup, F now means docn, new dublin testlist

- rm all grid combos using gx3v5, gx1v4, gx1v5 (obsolete pop grids)
  EXCEPT f02_g15 and all tri-grid compsets because there currently
  is no gx3v7 or gx1v6 (current pop grid) alternative for them
  (they don't have the necessary mapping/grid/domain files to run)
- rm all F compsets using camdom, rename "FD" compsets as "F"
  ("F" used to mean camdom and "FD" meant docn, now "F" means docn)
- rm C_PIO compset
- rm B_HIRES, B_1850-2000 (cam head) compset
- rm USEPIO from config_definition.xml (not used anywhere)
- rm all CASA compsets as CASA will NOT be supported (erik)
- Change name of I compsets to put CN at the end (erik)
- Change name of transient I compsets to use a dash "-" instead of an underscore "_" (erik)
- Add longer descriptions for CLM variables in config_definition file (erik)
- Add longer descriptions for compsets with CLM make use of either Satellite phenology
  (SP) or CN biogeochemistry phenology. (erik)

M      ccsm_utils/Case.template/config_grids.xml
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Testlists/dublin.posttag
M      ccsm_utils/Testlists/jaguar.clm.auxtest   -- Make consistent with new compset names
M      ccsm_utils/Testlists/bluefire.clm.auxtest -- Make consistent with new compset names
M      README_bgc --- Remove documentation about CLM_CO2_TYPE no longer used (erik)
M      ccsm_utils/Build/Makefile -- add SVN header

================================================================================
Originator: kauff
Date: Tue Dec 22 13:34:33 MST 2009
Model: scripts
Version: scripts4_091222
One-line: restore 1850/20th-cent track5 compsets, add T62_t12.G test case

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Testlists/jaguar.posttag

================================================================================
Originator: mvertens
Date: Mon Dec 21 2009
Model: scripts
Version: scripts4_091221
One-line: removed surveyor and added intrepid mods"

D      Macros.surveyor
D      env_machopts.surveyor
M      Macros.intrepid
D      mkbatch.surveyor
M      mkbatch.intrepid
M      config_machines.xml

M      ccsm_utils/Tools/check_input_data
       - corrected error report for input argument

================================================================================
Originator: kauff
Date: Fri Dec 18 17:11:20 MST 2009
Model: scripts
Version: scripts4_091218b
One-line: add 4x5_gx3v7 support, rm track5 compsets

- add 4x5_gx3v7 grid support
M      ccsm_utils/Case.template/config_grid.xml

- T31_g37.B  and f45_g37.BW are part of ccsm4.x -- add tests for these
M      ccsm_utils/jaguar.pretag
M      ccsm_utils/bluefire.pretag

- track5 cam is not part of ccsm4.x  --  rm TR5 compsets
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Testlists/newmachine.port.auxtest  replace BTR5 w/ BTR1CN
M      ccsm_utils/Testlists/midnight.auxtest
M      ccsm_utils/Testlists/atlas.auxtest
M      ccsm_utils/Testlists/hadley.auxtest

================================================================================
Originator: fischer
Date: Fri Dec 18 2009
Model: scripts
Version: scripts4_091218
One-line: add dublin intel and generic machines

M      ccsm_utils/Machines/Macros.generic_linux_intel
M      ccsm_utils/Machines/Macros.dublin_intel
M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Testlists/jaguar.posttag
 
================================================================================
Originator: tcraig
Date: Tue Dec 15 2009
Model: scripts
Version: scripts4_091215
One-line: update check_input_data for commments in input list

M      ccsm_utils/Tools/check_input_data

================================================================================
Originator: tcraig
Date: Mon Dec 14 2009
Model: scripts
Version: scripts4_091214
One-line: Add FDTR1H compset

M      ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: mvertens
Date: Sun Dec 13 2009
Model: scripts
Version: scripts4_091213
One-line: added new generic machine capability

- removed DOUT_L_RCP_ROOT, DOUT_L_RCP, ARCH, SITE xml variables 

- the following new arguments were added to create_newcase that must be specified for a generic machine
   -exedir <name>                ccsm executable directory (EXEROOT will be EXEDIR/CASE) (char)
   -din_loc_root_csmdata <name>  ccsm input data root directory (char)
   -max_tasks_per_node <value>   maximum mpi tasks per machine node (integer)
   -pes_per_node <value>         processors per machine node (integer)

- generic machine settings that must be user editable appear primarily in the env_mach_specific file
    with the entry USER_MUST_FILL_THIS_IN
 
- with these settings new generic machine settings in config_machines.xml are simple - as an example:
    <machine MACH="generic_linux_pgi"
             DESC="generic linux (pgi), os is Linux, batch system is PBS"
             OS="Linux" />
	
- changes are completely backwards compatible
  - verified that ERT.f09_g16.B1850TR1CN.bluefire worked and was bfb with beta36
  - verified that long term archiving and resubmit functionality still worked
	
- tested the generic setting for both generic_ibm (bluefire) and generic_linux_pgi(dublin)
  (note that could easily build on dublin by changing env_mach_specific - but there are
   currently problems running on this platform)
	
M      ccsm_utils/Tools/ccsm_l_archive.csh
M      ccsm_utils/Case.template/config_definition.xml
A      ccsm_utils/Machines/env_machopts.generic_linux_pgi
D      ccsm_utils/Machines/Macros.lightning_path
D      ccsm_utils/Machines/env_machopts.lightning_path
M      ccsm_utils/Machines/mkbatch.ranger
D      ccsm_utils/Machines/mkbatch.lightning_intel
M      ccsm_utils/Machines/mkbatch.ubgl
A      ccsm_utils/Machines/Macros.generic_linux_pathscale
D      ccsm_utils/Machines/env_machopts.lightning_intel
D      ccsm_utils/Machines/mkbatch.lightning_path
A      ccsm_utils/Machines/Macros.generic_linux_intel
A      ccsm_utils/Machines/mkbatch.generic_ibm
A      ccsm_utils/Machines/env_machopts.generic_ibm
A      ccsm_utils/Machines/Macros.generic_linux_lahey
A      ccsm_utils/Machines/mkbatch.generic_linux_pgi
A      ccsm_utils/Machines/Macros.generic_xt
A      ccsm_utils/Machines/env_machopts.generic_xt
A      ccsm_utils/Machines/Macros.generic_linux_pgi
A      ccsm_utils/Machines/mkbatch.dublin_intel
A      ccsm_utils/Machines/mkbatch.generic_linux_pathscale
A      ccsm_utils/Machines/mkbatch.generic_linux_intel
A      ccsm_utils/Machines/mkbatch.generic_xt
D      ccsm_utils/Machines/Macros.lightning_intel
M      ccsm_utils/Machines/mkbatch.columbia
A      ccsm_utils/Machines/env_machopts.generic_linux_pathscale
A      ccsm_utils/Machines/env_machopts.generic_linux_intel
M      ccsm_utils/Machines/mkbatch.atlas
A      ccsm_utils/Machines/mkbatch.generic_linux_lahey
A      ccsm_utils/Machines/Macros.generic_ibm
A      ccsm_utils/Machines/Macros.dublin_intel
A      ccsm_utils/Machines/env_machopts.dublin_intel
A      ccsm_utils/Machines/env_machopts.generic_linux_lahey
M      ccsm_utils/Machines/config_machines.xml
M      ccsm_utils/Machines/config_pes.xml
M      create_newcase
	
================================================================================
Originator: tcraig
Date: Tue Dec 8 2009
Model: scripts
Version: scripts4_091208
One-line: add gregorian calendar and single column support

 - add CALENDAR env variable in env_run
 - add ESMF_INTERFACE and MCT_INTERFACE cpp tokens, needed for esmf register routines
 - remove restriction on identical grids for single column mode
 - update CCSM_SSTDATA_FILE to reflect latest version with calendar attribute in time units

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Machines/Macros.cppdefs
M      create_newcase

================================================================================

Originator: tcraig
Date: Wed Dec 2 2009
Model: scripts
Version: scripts4_091202
One-line: update check_input_data, add "FD" temporary compsets, add ne grids, CCSM_SSTDATA

M      ccsm_utils/Tools/check_input_data
M      ccsm_utils/Tools/ccsm_buildnml.csh
M      ccsm_utils/Tools/ccsm_prestage.csh
M      ccsm_utils/Tools/generate_batch.csh
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_grid.xml
M      ccsm_utils/Case.template/config_definition.xml
	
================================================================================
Originator: kauff
Date: Wed Nov 18 17:01:20 MST 2009
Model: scripts
Version: scripts4_091118
One-line: mods to Testlists in support of release testing

- files ending in .pretag & .posttag are the official ccsm4.0 pre & post tag tests
  all other files get the .auxtest extension

M   atlas.auxtest
M   bluefire.cam.auxtest
M   bluefire.cice1.auxtest
M   bluefire.cice2.auxtest
M   bluefire.clm.auxtest
M   bluefire.pop2.auxtest
M   bluefire.posttag
M   bluefire.pretag
M   dublin.posttag
M   franklin.posttag
M   hadley.auxtest
M   intrepid.posttag
M   jaguar.cam.auxtest
M   jaguar.clm.auxtest
M   jaguar.posttag
M   jaguar.pretag
M   kraken.posttag
M   midnight.auxtest
M   newmachine.port.auxtest

- add LM_BLDNML_OPTS="-use_case glacier_mec" when using glimmer
- expand tabs into spaces
M      config_compsets.xml

================================================================================
Originator: tcraig
Date: Fri Nov 13 2009
Model: scripts
Version: scripts4_091113
One-line: new directory structure for data and drv

- add hadley (UCB Linux with intel compiler)
- change Filepath in csm_share build to be consistent with driver and csm_share changes
- change st archiver to reference to datm instead of datm7 (same)
- change references to DATM7 to DATM (same for other data models)
- change data model sandbox reference from datm7 to datm (same)
	
M      ccsm_utils/Tools/st_archive.sh
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Components/csm_share.buildlib
M      create_newcase
A      ccsm_utils/Machines/Macros.hadley
A      ccsm_utils/Machines/env_machopts.hadley
A      ccsm_utils/Machines/mkbatch.hadley
M      ccsm_utils/Machines/config_machines.xml
A      ccsm_utils/Testlists/hadley.posttag
A      ccsm_utils/Testlists/hadley.pretag
	
================================================================================

Originator: jwolfe
Date: Fri Nov 06 2009
Model: scripts
Version: scripts4_091106b
One-line: Correct p-netcdf module on franklin and remove some lightning files

 - remove lightning_lahey and lightning_pgi files, since lightning is set to be
   decommissioned soon.  Leaving lightning_path and lightning_intel as examples
   until we have other machine files using those compilers
 - fix problem with name of franklin p-netcdf module

M      ccsm_utils/Machines/env_machopts.franklin
D      ccsm_utils/Machines/mkbatch.lightning_lahey
D      ccsm_utils/Machines/env_machopts.lightning_lahey
D      ccsm_utils/Machines/env_machopts.lightning_pgi
D      ccsm_utils/Machines/Macros.lightning_lahey
D      ccsm_utils/Machines/mkbatch.lightning_pgi
D      ccsm_utils/Machines/Macros.lightning_pgi


================================================================================

Originator: tcraig
Date: Fri Nov 06 2009
Model: scripts
Version: scripts4_091106
One-line: fix ERH script for reproducibility (bug #1067)

M      ccsm_utils/Testcases/ERH_script

================================================================================

Originator: jwolfe
Date: Tue Nov 03 2009
Model: scripts
Version: scripts4_091103
One-line: Updated location of inputdata and shared disk space on franklin

M      ccsm_utils/Machines/config_machines.xml

================================================================================

Originator: mvertens
Date: Mon Nov 02 2009
Model: scripts
Version: scripts4_091102
One-line: Misc added compset F_AMIP_1DEG_TRACK1_CMIP5	

M   ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: jwolfe
Date: Fri Oct 30 2009
Model: scripts
Version: scripts4_091030
One-line: Misc machine file clean-up

 - in env_machopts.bluefire, add MP_BUFFER_MEM setting to non-default value of
   64M to help fix Keith Lindsay's mpi hang issue
 - in Macros.bluefire, add "-d" flag to POP compilation to save intermediate
   preprocessing .f90 files
 - change SPEC_FFLAGS to FFLAGS_NOOPT to make usage consistent in all Macros
 - add FFLAGS_NOOPT settings to Macros logic for debug and threading to make
   sure special files are compiled consistently with the rest of CCSM
 - for PGI compilers, change "-fast" to "-O2" due to optimization issues and
   potential changes to the "-fast" definition
 - removed settings for xlf12 on bluefire, since it is now the default compiler
 - add new posttag Testlist  for intrepid

M      ccsm_utils/Machines/Macros.lightning_path
M      ccsm_utils/Machines/Macros.dublin_lahey
M      ccsm_utils/Machines/Macros.dublin_pgi
M      ccsm_utils/Machines/Macros.intrepid
M      ccsm_utils/Machines/Macros.schirra
M      ccsm_utils/Machines/Macros.calgary_lahey
M      ccsm_utils/Machines/Macros.kraken
M      ccsm_utils/Machines/Macros.bluefire
M      ccsm_utils/Machines/env_machopts.bluefire
M      ccsm_utils/Machines/Macros.franklin
M      ccsm_utils/Machines/Macros.jaguar
M      ccsm_utils/Machines/Macros.lightning_intel
M      ccsm_utils/Machines/Macros.jaguarpf
M      ccsm_utils/Machines/Macros.midnight
M      ccsm_utils/Machines/Macros.atlas
M      ccsm_utils/Machines/Macros.lightning_lahey
M      ccsm_utils/Machines/Macros.ranger
M      ccsm_utils/Machines/Macros.calgary_pgi
M      ccsm_utils/Machines/Macros.lightning_pgi
M      ccsm_utils/Machines/Macros.columbia
A      ccsm_utils/Testlists/intrepid.posttag

================================================================================
Originator: erik
Date: Tue Oct 27 2009
Model: scripts
Version: scripts4_091027b
One-line: Update perl5lib to perl5lib_091027

Update the version of perl5lib to perl5lib_091027. This enables the CLM_USRDAT_NAME
feature in scripts.

================================================================================
Originator: fischer
Date: Tue Oct 27 2009
Model: scripts
Version: scripts4_091027
One-line: Added compset_file to create_newcase and create_test.  Added DESC to config_grid.xml

M      create_test                                  Added compset_file option
M      create_newcase                               Added compset_file_option
A      compset_file_sample.xml                      Sample compset file
M      ccsm_utils/Case.template/config_grid.xml     Added DESC variable to grids

================================================================================
Originator: erik
Date: Wed Oct 21 2009
Model: scripts
Version: scripts4_091021
One-line: Some changes to Makefile so clm can use it as source, add in CLM_USRDAT_NAME and DATM_CO2_TSERIES

(moved from branches/co2tseries branch)

M  ccsm_utils/Build/Makefile ---------------------- Add check if: ULIBS, CLIBS, and ULIBDEF defined before being set
                                                    Add comment to MACFILE include, change RM EXEC to EXEC_SE as EXEC 
                                                    no longer exists
M  ccsm_utils/Case.template/config_compsets.xml --- Add DATM_CO2_TSERIES=TRUE for I 1850-2000 compsets
M  ccsm_utils/Case.template/config_definition.xml - Add CLM_USRDAT_NAME and DATM_CO2_TSERIES
M  ccsm_utils/Testlists/bluefire.clm.pretag ------- Add some RLA, RLB, and RLO tests

Testing: Check that A and B cases on bluefire and jaguar work

SMS.f19_g16.BTR1.jaguar ---- Builds
SMS.f19_g16.BTR1.bluefire -- Builds
SMS.f19_g16.A.jaguar ------- Builds and runs
ERS.f19_g16.A.bluefire ----- Builds and runs

Also I cases build and run

================================================================================
Originator: jwolfe
Date: Thu Oct 15 2009
Model: scripts
Version: scripts4_091015
One-line: make PEM and PET tests in testlists more relevant

M      ccsm_utils/Testlists/bluefire.posttag

================================================================================
Originator: tcraig
Date: Tue Oct 13 2009
Model: scripts
Version: scripts4_091013
One-line: general improvements

 - update to kraken machine files
 - updates for intrepid
 - fix ERP and ERT test, rest_option setting proposed by dbailey
 - turn off ERT compare cpl log test (#1020)
 - change cpl log message for missing comm diag output to UNDEF
 - fix test output throughput and memory diagnostics
 - fix create_suite scriptsroot bug  (#1054)
 - fix PFS test output
 - modify test generate behaviour when pass fails, data is not longer
   saved in generate mode
 - add _R* single point option to create_test  (#1047)
 - add _R single point tests to bluefire, dublin, and jaguar posttag lists

M      create_test
M      ccsm_utils/Tools/check_exactrestart.pl
M      ccsm_utils/Tools/testcase_end
M      ccsm_utils/Machines/Macros.intrepid
M      ccsm_utils/Machines/Macros.kraken
M      ccsm_utils/Machines/mkbatch.intrepid
M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Testcases/config_tests.xml
M      ccsm_utils/Testcases/ERP_script
M      ccsm_utils/Testcases/ERT_script
M      ccsm_utils/Testlists/jaguar.posttag
M      ccsm_utils/Testlists/dublin.posttag
M      ccsm_utils/Testlists/bluefire.posttag
M      create_suite
	
================================================================================
Originator: tcraig
Date: Mon Oct 12 2009
Model: scripts
Version: scripts4_091012
One-line: add P4A test and B1850TR1CNRM compset

 - and add 4 tests to bluefire.pretag

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Testcases/config_tests.xml
A      ccsm_utils/Testcases/P4A_script
M      ccsm_utils/Testlists/bluefire.pretag
	
================================================================================
Originator: tcraig
Date: Tue Oct 06 2009
Model: scripts
Version: scripts4_091006
One-line: fix qstat and qsub paths, need full paths right now

- the path /usr/local/torque/bin is not in the default path so qstat and qsub
  need to be defined with the full path.

M      ccsm_utils/Machines/config_machines.xml
	
================================================================================
Originator: tcraig
Date: Thu Oct 01 2009
Model: scripts
Version: scripts4_091001a
One-line: fix MPI_SERIAL mode

- fix MPI_SERIAL mode by changing if test in all Macros files
- add PEA test, compares standard run with 1 pe run in cpl BFB mode
- add PEA and esmf tests to bluefire posttag lists
- add pecount = 1 to config_pes.xml and use that instead of special
  ptsmode_match in create_newcase
- add HAVE_MPI cpp to Macros.cppdefs for timing lib

M      ccsm_utils/Machines/Macros.surveyor
M      ccsm_utils/Machines/Macros.lightning_path
M      ccsm_utils/Machines/Macros.dublin_lahey
M      ccsm_utils/Machines/Macros.dublin_pgi
M      ccsm_utils/Machines/Macros.intrepid
M      ccsm_utils/Machines/Macros.cppdefs
M      ccsm_utils/Machines/Macros.schirra
M      ccsm_utils/Machines/Macros.calgary_lahey
M      ccsm_utils/Machines/Macros.kraken
M      ccsm_utils/Machines/Macros.ubgl
M      ccsm_utils/Machines/Macros.bluefire
M      ccsm_utils/Machines/Macros.frost
M      ccsm_utils/Machines/Macros.nyblue
M      ccsm_utils/Machines/Macros.franklin
M      ccsm_utils/Machines/Macros.jaguar
M      ccsm_utils/Machines/Macros.lightning_intel
M      ccsm_utils/Machines/Macros.jaguarpf
M      ccsm_utils/Machines/Macros.midnight
M      ccsm_utils/Machines/Macros.atlas
M      ccsm_utils/Machines/Macros.lightning_lahey
M      ccsm_utils/Machines/Macros.ranger
M      ccsm_utils/Machines/Macros.calgary_pgi
M      ccsm_utils/Machines/Macros.lightning_pgi
M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Machines/Macros.columbia
M      ccsm_utils/Testcases/config_tests.xml
A      ccsm_utils/Testcases/PEA_script
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/bluefire.pretag
M      create_newcase
		
================================================================================
Originator: jwolfe
Date: Thu Oct 01 2009
Model: scripts
Version: scripts4_091001
One-line: Remove redundant test from bluefire pre-tag testlist 

M      ccsm_utils/Testlists/bluefire.pretag

================================================================================
Originator: mvr
Date: Mon Sep 28 2009
Model: scripts
Version: scripts4_090928
One-line: lt archive fix to enable saving of interim restart files; added handling 
of dart output to st archive script

- fix is for bug #1039
	
M      ccsm_utils/Tools/st_archive.sh
M      ccsm_utils/Tools/ccsm_l_archive.csh
	
================================================================================
Originator: tcraig
Date: Sat Sep 26 2009
Model: scripts
Version: scripts4_090926
One-line: Updates for ESMF

- update perf_summary.pl for cleaner output
- bug fix for PET, PEM, SEQ pe setup in create_test
- update valid values for DLND_RUNOFF_MODE
- add CME test (Compare Mct Esmf implementations)
	
M      create_test
M      ccsm_utils/Tools/timing/perf_summary.pl
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Testcases/config_tests.xml
A      ccsm_utils/Testcases/CME_script
	
================================================================================
Originator: jwolfe
Date: Fri Sep 25 2009
Model: scripts
Version: scripts4_090925
One-line: Update machine files for Cray XT's

- update franklin, jaguar, jaguarpf, and kraken machopts files to use pgi 9.0.2
  and a consistent set of other modules
- update support for jaguarpf to use the right number of PES_PER_NODE
- move jaguarpf pe configurations to be separate from other XT's, since it has
  a different number of PES_PER_NODE

M      ccsm_utils/Machines/env_machopts.kraken
M      ccsm_utils/Machines/env_machopts.franklin
M      ccsm_utils/Machines/env_machopts.jaguar
M      ccsm_utils/Machines/env_machopts.jaguarpf
M      ccsm_utils/Machines/config_machines.xml
M      ccsm_utils/Machines/config_pes.xml

================================================================================
Originator: njn01
Date: Tue Sep 22 2009
Model: scripts
Version: scripts4_090922
One-line: Modify pretag and posttag gx3v7 test lists

M           18413   ccsm_utils/Testlists/bluefire.cice1.pretag
M           18413   ccsm_utils/Testlists/bluefire.pretag

================================================================================
Originator: dbailey
Date: Mon Sep 21 2009
Model: scripts
Version: scripts4_090921
One-line: Additional IPCC/CMIP sea ice tracers.

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: fischer
Date: Wed Sep 16 2009
Model: scripts
Version: scripts4_090916
One-line: new dublin test list

A      ccsm_utils/Testlists/dublin.posttag

================================================================================
Originator: mvr
Date: Tue Sep 15 2009
Model: scripts
Version: scripts4_090915b
One-line: lt archiver enabled to tar daily pop output

M      ccsm_utils/Tools/ccsm_l_archive.csh

================================================================================
Originator: kauff
Date: Tue Sep 15 14:49:08 MDT 2009
Model: scripts
Version: scripts4_090915
One-line: add T62_gx3v7 support, T31/T62_gx3v5 tests changed to gx3v7

M      ccsm_utils/Case.template/config_grid.xml

M      ccsm_utils/Testlists/lightning_path.pretag
M      ccsm_utils/Testlists/franklin.posttag
M      ccsm_utils/Testlists/franklin.pretag
M      ccsm_utils/Testlists/calgary_lahey.pretag
M      ccsm_utils/Testlists/lightning_lahey.pretag
M      ccsm_utils/Testlists/jaguar.posttag
M      ccsm_utils/Testlists/jaguar.pretag
M      ccsm_utils/Testlists/kraken.posttag
M      ccsm_utils/Testlists/kraken.pretag
M      ccsm_utils/Testlists/dublin.pretag
M      ccsm_utils/Testlists/testlist.port
M      ccsm_utils/Testlists/lightning_intel.pretag
M      ccsm_utils/Testlists/lightning.posttag
M      ccsm_utils/Testlists/lightning_pgi.pretag
M      ccsm_utils/Testlists/lightning.pretag
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/bluefire.pretag

================================================================================
Originator: mvr
Date: Sep 10 2009
Model: scripts
Version: scripts4_090910
One-line: remove noopt flag from -qsmp setting due to performance hit

M      ccsm_utils/Machines/Macros.bluefire

================================================================================
Originator: tcraig
Date: Tue Sep  9 2009
Model: scripts
Version: scripts4_090908
One-line: updates to scripts, fixes

 - fix bug #1015, extra restarts in ERT and ERP
 - fix bug in env_build checking wrt cleanmach
 - fix SEQ script for new env_build checking
 - increase pecount="T" pe counts for faster turnaround
 - set CCSM_CCOST to >1 for several slow compsets (BW, BTR5, etc)
 - fix bug #1012, automatic creation of testroot in create_test

M      create_test
M      ccsm_utils/Tools/configure
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Testcases/ERP_script
M      ccsm_utils/Testcases/SEQ_script
M      ccsm_utils/Testcases/ERT_script
	
================================================================================

Originator: tcraig
Date: Mon Sep  7 2009
Model: scripts
Version: scripts4_090907a
One-line: updates to dead/stub models for esmf

 - move stub and dead template scripts to each component bld directory
   and update create_newcase to reflect those changes
 - update env_build and Macros checking and error messages (bug #1018)
 - add S_PRESENT_DAY (S) compset which is all stub + datm
 - update franklin to remove module load xtpe-quadcore and PBS -l feature=quad
   at request of NERSC.  these are obsolete, unnecessary, or problematic with
   next upgrade.
 - add support for _DE and _ED (same thing) test options
 - add ERS_DE.f19_g16.A.bluefire to bluefire.pretag testlist
	

M      ccsm_utils/Tools/configure
M      ccsm_utils/Tools/ccsm_buildexe.csh
M      ccsm_utils/Tools/ccsm_check_lockedfiles
M      ccsm_utils/Tools/clean_build
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_definition.xml
D      ccsm_utils/Components/xlnd.template
D      ccsm_utils/Components/xocn.template
D      ccsm_utils/Components/sglc.template
D      ccsm_utils/Components/xice.template
D      ccsm_utils/Components/xatm.template
D      ccsm_utils/Components/xglc.template
M      ccsm_utils/Components/csm_share.buildlib
D      ccsm_utils/Components/slnd.template
D      ccsm_utils/Components/socn.template
D      ccsm_utils/Components/sice.template
D      ccsm_utils/Components/satm.template
M      ccsm_utils/Machines/mkbatch.franklin
M      ccsm_utils/Machines/env_machopts.franklin
M      ccsm_utils/Testlists/bluefire.pretag
M      create_newcase
	
================================================================================

Originator: kauff
Date: Fri Sep  4 11:10:30 MDT 2009
Model: scripts
Version: scripts4_090904
One-line: new/corrected T31_gx3v7 mapping files

M      Case.template/config_grid.xml

================================================================================
Originator: jwolfe
Date: Sep 03 2009
Model: scripts
Version: scripts4_090903a
One-line: Change CAM thread optimization to -O2 on bluefire

M      ccsm_utils/Machines/Macros.bluefire

================================================================================
Originator: jwolfe
Date: Sep 03 2009
Model: scripts
Version: scripts4_090903
One-line: Change CAM optimization to -O2 on bluefire and update pnetcdf module on kraken

M      ccsm_utils/Machines/env_machopts.kraken
M      ccsm_utils/Machines/Macros.bluefire

================================================================================
Originator: jwolfe
Date: Sep 02 2009
Model: scripts
Version: scripts4_090902a
One-line: add missed files necessary to run on dublin

A      ccsm_utils/Machines/mkbatch.dublin_pgi
A      ccsm_utils/Machines/mkbatch.dublin_lahey

================================================================================
Originator: feiliu
Date: Sep 02 2009
Model: scripts
Version: scripts4_090902
One-line: added capability to automate running ESMF enabled experiments

M      ccsm_utils/Machines/config_machines.xml
M      create_newcase

================================================================================
Originator: kauff
Date: Tue Sep  1 16:34:45 MDT 2009
Model: scripts
Version: scripts4_090901c
One-line: T31_gx3v7 support, new tx0.1 mapping files

- new grid support: T31_gx3v7
- new mapping file: T31_gx3v5 patch map, atm->ocn
- new fv0.5, fv0.25  <-> tx0.1 mapping files (fixes pop domain bug in mapping files).
- document epbal & albav in config_definition.xml

M      Case.template/config_grid.xml
M      ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: fischer
Date: Sep 01 2009
Model: scripts
Version: scripts4_090901a
One-line: added non CN pretag tests for bluefire and jaguar

M      ccsm_utils/Testlists/jaguar.pretag
M      ccsm_utils/Testlists/bluefire.pretag

================================================================================
Originator: jwolfe
Date: Sep 01 2009
Model: scripts
Version: scripts4_090901
One-line: Make XLF12 the default compiler on bluefire

M      ccsm_utils/Machines/Macros.bluefire
M      ccsm_utils/Machines/env_machopts.bluefire

================================================================================
Originator: erik
Date: Aug 29 2009
Model: scripts
Version: scripts4_090829
One-line: Autoflip PES=1, USE_MPISERIAL=TRUE for PTSMODE

Update externals to perl5lib_090829

M      ccsm_utils/Machines/config_pes.xml -- Add setting when ptsmode=TRUE
M      create_newcase ---------------------- Turn USE_MPISERIAL=TRUE and
               pass ptsmode to set_pes.

================================================================================
Originator: jwolfe
Date: Aug 28 2009
Model: scripts
Version: scripts4_090828
One-line: add support for machine dublin with pgi and lahey compilers

A      ccsm_utils/Machines/Macros.dublin_lahey
A      ccsm_utils/Machines/env_machopts.dublin_lahey
A      ccsm_utils/Machines/Macros.dublin_pgi
A      ccsm_utils/Machines/env_machopts.dublin_pgi
M      ccsm_utils/Machines/config_machines.xml
M      ccsm_utils/Machines/config_pes.xml
A      ccsm_utils/Testlists/dublin.pretag

================================================================================
Originator: fischer
Date: Aug 14 2009
Model: scripts
Version: scripts4_090814
One-line: new E_1850_TRACK1_CN compset, and more changes to Testlists

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Testlists/jaguar.posttag
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/bluefire.pretag

================================================================================
Originator: fischer
Date: Aug 11 2009
Model: scripts
Version: scripts4_090811
One-line: Added new WACCM: compset, pes config, and Testlist

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Testlists/bluefire.pretag
M      ccsm_utils/Testlists/bluefire.posttag

================================================================================
Originator: erik
Date: Aug 06 2009
Model: scripts
Version: scripts4_090806
One-line: Update perl5lib, move check for camdom grids to create_newcase

Update perl5lib to perl5lib_090704 -- adds get_usr_file method, needed for new
   clm build-namelist feature.

Move check for camdom atm and ocn grids to be identical to create_newcase rather
than configure. The problem is that you don't know the problem until the configure
step -- you have to rerun create_newcase. This resolves bug 925.

M      ccsm_utils/Tools/configure -- Remove the check for camdom grids move 
                                     it to create_newcase
M      create_newcase -------------- Put the check for camdom grids in set_compset

================================================================================
Originator: mvr,hannay
Date: Aug 05 2009
Model: scripts
Version: scripts4_090805
One-line: Create a new compset for Track5 transient

M      ccsm_utils/Case.template/config_compsets.xml

================================================================================
Originator: Erik Kluzek
Date: Aug 01 2009
Model: scripts
Version: scripts4_090801
One-line: Add some additional clm tests

M      ccsm_utils/Testlists/bluefire.clm.pretag
A      ccsm_utils/Testlists/jaguar.clm.posttag

================================================================================
Originator: Chris Fischer
Date: Thu Jul 30 10:20:50 MDT 2009
Model: scripts
Version: scripts4_090730
One-line: Change to Testlists, pe layout for C compsets on jaguar, and WACCM changes

Changed B20TR1 and B1850TR1 tests to B20TR1CN and B1850TR1CN
Changed pe layout for C compsets on jaguar from the default to the same as on 
bluefire.
Made changes to config_compsets.xml for WACCM 1850 runs

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Machines/config_pes.xml
A      ccsm_utils/Testlists/franklin.posttag
M      ccsm_utils/Testlists/franklin.pretag
M      ccsm_utils/Testlists/jaguar.pretag
M      ccsm_utils/Testlists/bluefire.pretag

================================================================================
Originator: Chris Fischer
Date: Tue Jul 28 15:40:55 MDT 2009
Model: scripts
Version: scripts4_090728c
One-line: Change ERP to ERT tests
M      ccsm_utils/Testlists/franklin.pretag
M      ccsm_utils/Testlists/jaguar.posttag
M      ccsm_utils/Testlists/jaguar.pretag
M      ccsm_utils/Testlists/kraken.pretag
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/bluefire.pretag

================================================================================
Originator: Fei Liu
Date: Tue Jul 28 13:15:56 MDT 2009
Model: scripts
Version: scripts4_090728b
One-line: changes make scripts build with ESMF library turned on or off
M      ccsm_utils/Components/csm_share.buildlib
M      ccsm_utils/Machines/Macros.cppdefs

================================================================================
Originator: erik
Date: Jul 28 2009
Model: scripts
Version: scripts4_090728
One-line:  Add PTS_MODE, PTS_LAT, PTS_LON and command line option to set these

M      ccsm_utils/Case.template/ConfigCase.pm --------- Only add pts_mode stuff if turned on
M      ccsm_utils/Case.template/config_definition.xml - Add PTS_MODE, PTS_LAT, PTS_LON definition
M      create_newcase --------------------------------- Add -pts_lat/-pts_lon command line option
             This will then set PTS_MODE in env_case and PTS_LAT/PTS_LON in env_run
             also ensures that grids are identical, and don't set ???_NX/???_NY is pts-mode

================================================================================
Originator: jwolfe
Date: Wed Jul 22 15:44:19 MDT 2009
Model: scripts
Version: scripts4_090722b
One-line: Change jaguar module to fix build problem with old version

M      ccsm_utils/Machines/env_machopts.jaguar

================================================================================
Originator: kauff
Date: Wed Jul 22 14:04:01 MDT 2009
Model: scripts
Version: scripts4_090722
One-line: add support for jaguarpf, ORNL XT5

Mods provided by Pat Worley

A    Machines/mkbatch.jaguarpf
A    Machines/Macros.jaguarpf
A    Machines/env_machopts.jaguarpf
M    Machines/config_pes.xml
M    Machines/config_machines.xml

================================================================================
Originator: mvr, hannay
Date: Tue Jul 21 2009
Model: scripts
Version: scripts4_090721
One-line: enable environment variable for saving all restart files generated, not 
just those from end-of-run; change the cam use case for track5 1850 compset

M      ccsm_utils/Tools/st_archive.sh
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: erik
Date: Mon Jul 20 2009
Model: scripts
Version: scripts4_090720
One-line: Add CLM1PT option for DATM_MODE, remove CLM_ARB_IC, document CLM_* vars

Add a new CLM1PT option for DATM_MODE. Remove the CLM_ARB_IC variable as it is 
confusing. And make sure CLM_* variables have some documentation on them.

M      ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: erik
Date: Fri Jul 17 2009
Model: scripts
Version: scripts4_090717
One-line: Update perl5lib version, add 2000_control use_cases for CLM 

Update version of perl5lib so that will expand env variable DIN_LOC_ROOT
is chng_din_lc option is on. This allows CLM to use the datm7 version of the
streams template rather than keeping it's own copy.

M      ccsm_utils/Case.template/config_compsets.xml -- Set 2000 control cases
                   for CLM to the 2000_control use_case rather than setting
                   sim_yr to 2000 (in CLM_BLDNML_OPTS).

================================================================================
Originator: mvertens
Date: Thu Jul 16 2009
Model: scripts
Version: scripts4_090716
One-line: put in fix for using csm_share build with ESMF library

M      ccsm_utils/Components/csm_share.buildlib
	
================================================================================

Originator: kauff
Date: Wed Jul 15 21:26:51 MDT 2009
Model: scripts
Version: scripts4_090715
One-line: fv1.9_tx1v1: add support; T62_tx0.1v2: CDF64 true by default, new pe layout

M      ccsm_utils/Case.template/config_grid.xml
M      ccsm_utils/Machines/config_pes.xml

================================================================================
Originator: erik
Date: Tue Jul 14 2009
Model: scripts
Version: scripts4_090714
One-line: Update per5lib, and change CLM_BLDNML_OPTS for clm3_6_46

Update perl5lib version to perl5lib_090714 so has option of converting DIN_LOC_ROOT
to $DIN_LOC_ROOT.

M      ccsm_utils/Case.template/config_compsets.xml -- Change CLM_BLDNML_OPTS
             to use -use_case 1850_control and 20thC_transientinstead of specific
             simulation year. This is needed for clm3_6_46 and later.

================================================================================
Originator: erik/mvertens
Date: Tue Jul 7 2009
Model: scripts
Version: scripts4_090707b
One-line: add BGC spinup compset

Add a compset to use CPLHIST forcing data to run clm. Add a tool to help create files for that usage.
Remove -maxpft for CN with clm, as this is already inside of clm. Change name of B_1850_TRACK5 to include
"DEV" at end to indicate it's a development version.

A      ccsm_utils/Tools/concat_daily_hist.csh --------- Script from Forrest to concatenate data on monthly files
M      ccsm_utils/Case.template/config_compsets.xml --- Remove -maxpft for CN, add I1850SPINUPCN compset, add DEV to TR5 compsets
M      ccsm_utils/Case.template/config_definition.xml - Add options to DATM_MODE: CPLHIST3HrWx,CPLHIST3HrWxHfHrSol
                                                        remove CAMHIST,CPLHIST
M      ccsm_utils/Testlists/bluefire.clm.pretag ------- Add testing of f09_g16 I1850SPINUPCN

================================================================================
Originator: fischer
Date: Tue Jul 7 2009
Model: scripts
Version: scripts4_090707
One-line: added WACCM compsets

M      ccsm_utils/Case.template/config_compsets.xml
	
================================================================================
Originator: tcraig
Date: Wed Jul 1 2009
Model: scripts
Version: scripts4_090701b
One-line: point share filepath to new driver cpl_shr, cpl_mct and cpl_esmf dirs

M      ccsm_utils/Components/csm_share.buildlib
	
================================================================================
	
Originator: fischer
Date: Wed July 1 2009
Model: scripts
Version: scripts4_090701
One-line:  Updated READMEs and comments for env_mach_pes.xml, added new doc from Erika Marcum.

A      doc/styles
A      doc/section1
A      doc/section2
A      doc/section3
A      doc/section4
A      doc/section5
A      doc/section6
A      doc/section7
A      doc/section8
A      doc/images
A      doc/includes
A      doc/index.html
M      README
       Changed references to newer env_mach_pes.xml
M      ccsm_utils/Tools/configure
       Changed references to newer env_mach_pes.xml
M      ccsm_utils/Case.template/ConfigCase.pm
       Added comments about variables in env_mach_pes.xml
M      ccsm_utils/Testlists/bluefire.posttag
       -SMS.f45_g35.FTR1.bluefire
       +SMS.f45_f45.FTR1.bluefire
M      README_quickstart
       Changed references to newer env_mach_pes.xml

================================================================================

Originator: tcraig
Date: Sun Jun 28 2009
Model: scripts
Version: scripts4_090628
One-line: add a few new namelist inputs

 - add SHR_MAP_DOPOLE, NPFIX, AOFLUX_GRID env_run variables
 - add new timers to getTiming.pl for aoflux_grid

M      ccsm_utils/Tools/timing/getTiming.pl
M      ccsm_utils/Case.template/config_definition.xml
	
================================================================================

Originator: kauff
Date: Thu Jun 25 17:29:37 MDT 2009
Model: scripts
Version: scripts4_090625
One-line: add scripts support for T62_tx0.1v2 grid combo

M      ccsm_utils/Case.template/config_grid.xml
M      ccsm_utils/Machines/config_pes.xml

Note: this is for code developement purposes only at this point.
You can configure a case with this resolution, but it won't run...
eg. pop input files for partial-coupling mode (the C case default)
are missing so this code won't run out-of-the-box, but it's a good 
starting point for tx01v2 code development.

================================================================================

Originator: mvr
Date: Wed Jun 24 2009
Model: scripts
Version: scripts4_090624
One-line:  removed -Mrecursive from list of flags for pgi debug compilations

M      ccsm_utils/Machines/Macros.kraken
M      ccsm_utils/Machines/Macros.franklin
M      ccsm_utils/Machines/Macros.jaguar
M      ccsm_utils/Machines/Macros.midnight
M      ccsm_utils/Machines/Macros.atlas
M      ccsm_utils/Machines/Macros.ranger
M      ccsm_utils/Machines/Macros.calgary_pgi
M      ccsm_utils/Machines/Macros.lightning_pgi

================================================================================
Originator: erik
Date: Fri Jun 19 2009
Model: scripts
Version: scripts4_090619
One-line:  Add BCN compsets, Update perl5lib (adding %ymd option), and modify clm test list

Update perl5lib so that you can use the %ymd indicator to generate files with
year-month-day in them.

M      SVN_EXTERNAL_DIRECTORIES
  -ccsm_utils/Tools/perl5lib     https://svn-ccsm-models.cgd.ucar.edu/perl5lib/trunk_tags/perl5lib_090609
  +ccsm_utils/Tools/perl5lib     https://svn-ccsm-models.cgd.ucar.edu/perl5lib/trunk_tags/perl5lib_090613
M      ccsm_utils/Case.template/config_compsets.xml -- Add BCN compsets for TRACK1
M      ccsm_utils/Testlists/bluefire.clm.pretag ------ Add a PET test and make 2 tests debug


================================================================================

Originator: mvertens
Date: Wed Jun 10 2009
Model: scripts
Version: scripts4_090610b
One-line:  added COMP_INTERFACE and USE_ESMF_LIB to config_definitions

M      ccsm_utils/Case.template/config_definition.xml
       added COMP_INTERFACE (default value is MCT) and 
       USE_EMSF_LIB (default value is FALSE)
	    
M      ccsm_utils/Machines/env_machopts.bluefire
       removed MP_SINGLE_THREAD since this was causing a problem for mpiio
	
================================================================================
Originator: fischer
Date: Wed Jun 10 2009
Model: scripts
Version: scripts4_090610
One-line:  Added to the how to.  Added more user instructions to create_production_test, renamed faq.html

  Added instructions on how to transfer a production run to another user to the howto.html.
  Added more user instructions to create_production_test script
  Renamed faq.html to howto.html


M ccsm_utils/Tools/create_production_test
D doc/faq.html
A doc/howto.html

================================================================================
Originator: tcraig
Date: Tue Jun 09 2009
Model: scripts
Version: scripts4_090609b
One-line:  update to perl5lib_090609 to add offset support for streams

M      SVN_EXTERNAL_DIRECTORIES
  -ccsm_utils/Tools/perl5lib     https://svn-ccsm-models.cgd.ucar.edu/perl5lib/trunk_tags/perl5lib_090424
  +ccsm_utils/Tools/perl5lib     https://svn-ccsm-models.cgd.ucar.edu/perl5lib/trunk_tags/perl5lib_090609
	
================================================================================

Originator: mvr
Date: Tue Jun 09 2009
Model: scripts
Version: scripts4_090609
One-line:  made the run script execute tcsh rather than csh to avoid a word/line 
	limit which prevented large processor requests on bluefire

note, there is also a limit with tcsh, but this effectively doubles the allowable 
node request from 8 to 16 (running mpi-only)

M      ccsm_utils/Machines/mkbatch.bluefire

	
================================================================================
Originator: mvr
Date: Fri Jun 05 2009
Model: scripts
Version: scripts4_090605b
One-line:  lt archiver fix to do all tarring in consistent manner so that 
cmp will validate transfers to the archival system

note, problem appeared isolated to the xt's (jaguar, franklin)

M      ccsm_utils/Tools/ccsm_l_archive.csh

	
================================================================================
Originator: njn01
Date: Fri Jun 05 2009
Model: scripts
Version: scripts4_090605
One-line:  Replace "B" compset with "BTRK5DEV"; replace "B1850" compset B1850TR5; update tests

This change was made to protect unsuspecting users from creating a "B" case without
understanding exactly what cam version was included. By removing "B" option, the user 
must select the desired cam "track" when creating a new case.

M           16490   ChangeLog
M           16490   ccsm_utils/Case.template/config_compsets.xml
M           16490   ccsm_utils/Testlists/lightning_path.pretag
M           16490   ccsm_utils/Testlists/midnight.pretag
M           16490   ccsm_utils/Testlists/franklin.pretag
M           16490   ccsm_utils/Testlists/calgary_lahey.pretag
M           16490   ccsm_utils/Testlists/atlas.pretag
M           16490   ccsm_utils/Testlists/lightning_lahey.pretag
M           16490   ccsm_utils/Testlists/jaguar.posttag
M           16490   ccsm_utils/Testlists/jaguar.pretag
M           16490   ccsm_utils/Testlists/kraken.posttag
M           16490   ccsm_utils/Testlists/kraken.pretag
M           16490   ccsm_utils/Testlists/bluefire.pop2.pretag
M           16490   ccsm_utils/Testlists/testlist.port
M           16490   ccsm_utils/Testlists/lightning_intel.pretag
M           16490   ccsm_utils/Testlists/bluefire.cice1.pretag
M           16490   ccsm_utils/Testlists/bluefire.cice2.pretag
M           16490   ccsm_utils/Testlists/lightning_pgi.pretag
M           16490   ccsm_utils/Testlists/lightning.pretag
M           16490   ccsm_utils/Testlists/bluefire.posttag
M           16490   ccsm_utils/Testlists/bluefire.pretag
M           16490   create_newcase
	
================================================================================

Originator: tcraig
Date: Thu Jun 04 2009
Model: scripts
Version: scripts4_090604
One-line: Add capability to run tests with preset pe layouts

- use _P* in the testname to specify preset pe layout.  such as ERS_PS
  will set the pecount to S and ERS_PX1 will set the pecount to X1.
- add ERS_PT.f19_g16.BTR1.bluefire, ERS_PT.f19_f19.F.bluefire to
  bluefire.pretag.
	
M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Testlists/bluefire.pretag
M      create_newcase
	
================================================================================

Originator: erik
Date: Wed Jun 03 2009
Model: scripts
Version: scripts4_090603
One-line: Back out the move of MACFILE include to top of Makefile from scripts4_090527

M      ccsm_utils/Build/Makefile --- Move include of MACFILE to after compile_threaded
                                     section.

================================================================================

Originator: mvertens
Date: Mon Jun 1 2009
Model: scripts
Version: scripts4_090601b
One-line: added new low-resolution tests

M      ccsm_utils/Testlists/bluefire.posttag
	
================================================================================
Originator: mvertens
Date: Mon Jun 1 2009
Model: scripts
Version: scripts4_090601
One-line: fixed ice grid spec for 0.47x0.63_gx1v6

M      ccsm_utils/Case.template/config_grid.xml

================================================================================
Originator: tcraig
Date: Fri May 29 2009
Model: scripts
Version: scripts4_090529
One-line: update env_machopts.bluefire

- update env_machops.bluefire from jim edwards

M      ccsm_utils/Machines/env_machopts.bluefire
	
================================================================================
Originator: mvr,erik
Date: Wed May 27 2009
Model: scripts
Version: scripts4_090527b
One-line: mods to 'how-to' document; bug fix to lt archiver; create_newcase bug 
fix; Add pretag test list for clm for bluefire

M      ccsm_utils/Tools/ccsm_l_archive.csh
- bug fix to only tar history output that was generated monthly
M      doc/faq.html
- updated sections for output data and creating a case
M      create_newcase
- bug fix for invokation of check_rundb
A      ccsm_utils/Testlists/bluefire.clm.pretag 
Add pretag test list for clm for bluefire
Includes 4x5_g35, f19_g16, f09_g16, and f10_f10 (needs to be added to datm to work)
resolution cases, as well as 1850, 2000, 1850-2000 and 1850_CN cases. 


================================================================================
Originator: erik
Date: Wed May 27 2009
Model: scripts
Version: scripts4_090527
One-line: Small updates to allow the clm test-suite to use the ccsm build

Some small updates coming from clm. Minor fix to mkSrcfiles from the clm version.
Small changes to the makefile so that clm will be able to use the Makefile for it's 
own build used by it's test-suite.

M      ccsm_utils/Build/mkSrcfiles -- Make sure to close files when done.
M      ccsm_utils/Build/Makefile ---- Check if MACFILE set before setting generic value
                                      Move include of MACFILE to beginning. Always
                                      include MACFILE even for db_files command goal

================================================================================
Originator: tcraig
Date: Tue May 26 2009
Model: scripts
Version: scripts4_090526
One-line: modify pio cp and updates for intrepid

- updates from mark taylor for intrepid
- modify cp of pio code to the build area to exclude .F90
  pregenerated code with .F90.in files (bug #956)
	
M      ccsm_utils/Tools/timing/perf_summary.pl
M      ccsm_utils/Components/pio.buildlib
M      ccsm_utils/Machines/Macros.surveyor
M      ccsm_utils/Machines/Macros.intrepid
	
================================================================================

Originator: kauff
Date: Fri May 22 15:51:17 MDT 2009
Model: scripts
Version: scripts4_090522
One-line: add support for MPISERIAL (bluefire & calgary_lahey only)

* limited support for MPISERIAL build & run,  requires: pio41_prod
* calgary_lahey build works again
* activate BUDGETS in all C and G compsets (bugz #963)

M      ccsm_utils/Build/Makefile
M      ccsm_utils/Tools/configure
M      ccsm_utils/Tools/ccsm_buildexe.csh
M      ccsm_utils/Tools/generate_batch.csh
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_grid.xml
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Components/mct.buildlib
M      ccsm_utils/Components/pio.buildlib
M      ccsm_utils/Machines/mkbatch.bluefire
M      ccsm_utils/Machines/Macros.cppdefs
M      ccsm_utils/Machines/mkbatch.calgary_lahey
M      ccsm_utils/Machines/Macros.schirra
M      ccsm_utils/Machines/Macros.calgary_lahey
M      ccsm_utils/Machines/Macros.bluefire
M      ccsm_utils/Machines/env_machopts.calgary_lahey
M      ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Machines/config_machines.xml

Notes: 
* MPISERIAL is not part of the test suite, these compset were tested manually...
  -compset A -mach bluefire/calgary_lahey -res f19_g16  with & without MPISERIAL
  -compset I -mach bluefire/calgary_lahey -res pt1_pt1  with & without MPISERIAL
* When compiled for MPISERIAL, job can be run as a batch job or interactively

================================================================================

Originator: tcraig
Date: Mon May 18 2009
Model: scripts
Version: scripts4_090518a
One-line: add BUDGETS=TRUE to all B compsets

M      ccsm_utils/Case.template/config_compsets.xml
	
================================================================================

Originator: tcraig
Date: Mon May 18 2009
Model: scripts
Version: scripts4_090518
One-line: update handling of SMP setting, update faq

M      ccsm_utils/Tools/configure
M      ccsm_utils/Tools/ccsm_buildexe.csh
M      doc/faq.html

- moved SMP settings to ccsm_buildexe from configure.  added BUILD_THREADED
  support in SMP_VALUE settings. (bug #956)
- updated the faq.
	
================================================================================

Originator: mvertens
Date: Thu May 14 2009
Model: scripts
Version: scripts4_090514e
One-line: changes to scripts to removed 1850 specific tests
	
M      ccsm_utils/Case.template/ConfigCase.pm
M      ccsm_utils/Machines/config_machines.xml
M      ccsm_utils/Case.template/config_grid.xml
       Removed typos in config	
	
D      ccsm_utils/Testcases/ERH1850_script
M      ccsm_utils/Testcases/config_tests.xml
D      ccsm_utils/Testcases/ERP1850_script
D      ccsm_utils/Testcases/ERT1850_script
M      ccsm_utils/Testcases/ERB_script
M      ccsm_utils/Testcases/LAR_script
M      ccsm_utils/Testcases/ERH_script
M      ccsm_utils/Testcases/APT_script
M      ccsm_utils/Testcases/ERP_script
D      ccsm_utils/Testcases/ERB1850_script
M      ccsm_utils/Testcases/ERT_script
M      ccsm_utils/Testlists/bluefire.pretag
       Removed any hard-wiring in the tests of start_dates
	
================================================================================
Originator: njn01
Date: Thu May 14 2009
Model: scripts
Version: scripts4_090514d
One-line: changes to faq.html "Other"

M           16068   doc/faq.html
	
================================================================================
Originator: dbailey
Date: Thu May 14 2009
Model: scripts
Version: scripts4_090514c
One-line: changes to faq.html section 4

M      doc/faq.html
	
================================================================================
Originator: fischer
Date: Thu May 14 2009
Model: scripts
Version: scripts4_090514b
One-line: changes to faq.html sections 2 & 3

M      doc/faq.html
	
================================================================================
Originator: fischer
Date: Thu May 14 2009
Model: scripts
Version: scripts4_090514
One-line: changed ERS.f19_g16.F1850 to ERS.f19_f19.F1850 in Testlists

M      ccsm_utils/Testlists/jaguar.posttag
M      ccsm_utils/Testlists/kraken.posttag
M      ccsm_utils/Testlists/bluefire.cice2.pretag
M      ccsm_utils/Testlists/bluefire.posttag
	
================================================================================
Originator: tcraig
Date: Wed May 13 2009
Model: scripts
Version: scripts4_090513
One-line: add G pe layout (bug #949)

M      ccsm_utils/Machines/config_pes.xml
	
================================================================================

Originator: erik
Date: Tue May 6 2009
Model: scripts
Version: scripts4_090506
One-line: fix bug with CLM_DEMAND

Change default CLM_DEMAND value to "null", as that is what is used in CLM build-namelist
to trigger that a value isn't given.

M   ccsm_utils/Case.template/config_definition.xml

================================================================================

Originator: mvertens
Date: Tue May 5 2009
Model: scripts
Version: scripts4_090505c
One-line: resolved xml comment errors

M      ConfigCase.pm

================================================================================
Originator: fischer
Date: Tue May 5 2009
Model: scripts
Version: scripts4_090505b
One-line: added FAQ in doc directory 

A     doc/faq.html
A     doc/config_compsets.xml  
A     doc/config_grid.xml  
A     doc/config_machines.xml  
A     doc/config_compsets.xsl  
A     doc/config_grid.xsl  
A     doc/config_machines.xsl
	
================================================================================
Originator: mvertens
Date: Tue May 5 2009
Model: scripts
Version: scripts4_090505
One-line: changed F compsets to all use camdom and fixed xml related problems 

M      config_compsets.xml
M      ConfigCase.pm
	
================================================================================
Originator: mvertens
Date: Tue Apr 28 2009
Model: scripts
Version: scripts4_090428
One-line summary: numerous fixes to testing scripts

M      ccsm_utils/Testcases/config_tests.xml
A      ccsm_utils/Testcases/ERB1850_script
A      ccsm_utils/Testcases/ERH1850_script
A      ccsm_utils/Testcases/ERP1850_script
A      ccsm_utils/Testcases/ERT1850_script
       - added 1850 tests to enable 20th century tests to work
	 these fixes are not robust - but they were added quickly in order
	 to get the beta15 tag made
	
M      ccsm_utils/Testcases/ERB_script
M      ccsm_utils/Testcases/ERH_script
       - fixes to new name for short term archiving	
	
M      ccsm_utils/Testlists/franklin.pretag
M      ccsm_utils/Testlists/atlas.pretag
M      ccsm_utils/Testlists/jaguar.posttag
M      ccsm_utils/Testlists/jaguar.pretag
M      ccsm_utils/Testlists/kraken.posttag
M      ccsm_utils/Testlists/kraken.pretag
M      ccsm_utils/Testlists/bluefire.pop2.pretag
M      ccsm_utils/Testlists/testlist.port
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/bluefire.pretag
       - removed all references to g15 and replaced them with g16
 	 g15 tests should no longer be done and if there are any problems with g16 
	 they should be noted and resolved
	
================================================================================
Originator: mvertens
Date: Mon Apr 27 2009
Model: scripts
Version: scripts4_090427b
One-line summary: fixes to get transient compset working for CICE

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_definition.xml
	
================================================================================
	
Originator: mvr
Date: Mon Apr 27 2009
Model: scripts
Version: scripts4_090427
One-line summary: lt archive fix for more efficienct tarring; new compset B20TRTR1 
	and additional tests added to pretag lists

M      ccsm_utils/Tools/ccsm_l_archive.csh
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Testlists/franklin.pretag
M      ccsm_utils/Testlists/jaguar.pretag
M      ccsm_utils/Testlists/kraken.pretag
M      ccsm_utils/Testlists/bluefire.pretag

================================================================================

Originator: erik
Date: Fri Apr 24 2009
Model: scripts
Version: scripts4_090424
One-line summary: Update perl5lib, CCSM_BGC=CO2 and CLM:co2_type=diag for all B,E,F cases
    Tweak I-cases a bit

M      SVN_EXTERNAL_DIRECTORIES ------------------------ update perl5lib version
M      ccsm_utils/Case.template/config_compsets.xml ---- modify compsets

================================================================================
	
Originator: erik
Date: Fri Apr 23 2009
Model: scripts
Version: scripts4_090423
One-line summary: Fix bug with checking of tests id's (bug #936)

M      ccsm_utils/Case.template/ConfigCase.pm -- Add CCSM_PECOUNT to ignore_name list
	
================================================================================
Originator: mvertens
Date: Fri Apr 17 2009
Model: scripts
Version: scripts4_090417b
One-line summary: consistency of F compsets with cice settings
from being tarred; other archiver cleanup

added CICE_CONFIG_OPTS="-ntr_aero 0 -ntr_pond 1 -ntr_iage 0" to all F cases

M      ccsm_utils/Case.template/config_compsets.xml
	
================================================================================
	
Originator: mvr
Date: Fri Apr 17 2009
Model: scripts
Version: scripts4_090417
One-line summary: fix to lt archiver to prevent history files within restart directories 
from being tarred; other archiver cleanup

- bug was introduced at scripts4_090107 
	
D      ccsm_utils/Tools/ccsm_s_archive.csh
M      ccsm_utils/Tools/ccsm_postrun.csh
M      ccsm_utils/Tools/st_archive.sh
M      ccsm_utils/Tools/ccsm_l_archive.csh
M      create_newcase

================================================================================
Originator: erik
Date: Thu Apr 16 15:04:25 MDT 2009
Model: scripts
Version: scripts4_090416b
One-line summary: Remove CLM_BGC in compsets, add checks that correct names
    are in the compsets, add I_1850_2000 and I_CN_1850_2000 compsets

M      ccsm_utils/Case.template/config_compsets.xml - Tweak clm compsets (make
                                                      sure sim_year set correctly,
                                                      and fpftdyn demanded for transient
                                                      cases, remove EPS_AGRID setting 
                                                      for I cases.  Add new I transient 
                                                      compsets.
M      ccsm_utils/Case.template/ConfigCase.pm ------- Add is_ignore_name method, tweak
                                                      check for valid names. 
M      create_newcase ------------------------------- Check names.

================================================================================
	
Originator: dbailey
Date: Thu Apr 16 13:42:05 MDT 2009
Model: scripts
Version: scripts4_090416
One-line summary: Modify the CICE pretag tests.
	
M      ccsm_utils/Testlists/bluefire.cice1.pretag
M      ccsm_utils/Testlists/bluefire.cice2.pretag

================================================================================
Originator: mvertens
Date: Apr 14 2009
Model: scripts
Version: scripts4_090414b
One-line summary: various cleanup issues - new compsets, new definition
	
M      ccsm_utils/Tools/ccsm_getenv
       - moved env_run.xml to before env_conf.xml - needed for DIN_LOC_ROOT
	 setting in calls to build-namelists
M      ccsm_utils/Case.template/config_compsets.xml
       - introduced F_CLIMOTEST and G_CPLHIST (both will be modified but they
	 are the place holders for upcoming requirements)
M      ccsm_utils/Case.template/config_grid.xml
       - removed all references to map_r2o_file_r18	
M      ccsm_utils/Case.template/config_definition.xml
       - added CICE_NAMELIST_OPTS
	 
================================================================================
	
Originator: mvr
Date: Apr 14 2009
Model: scripts
Version: scripts4_090414
One-line summary: mods to archiving scripts to better handle mss path length limits

M      ccsm_utils/Tools/st_archive.sh
M      ccsm_utils/Tools/ccsm_l_archive.csh
M      ccsm_utils/Machines/Macros.intrepid
- minor change to machine setting
	
================================================================================
Originator: mvertens
Date: Apr 08 2009
Model: scripts
Version: scripts4_090408
One-line summary: added more changes for CORE2 forcing 

- updated tests lists
- added more documentation to create_clone
- added new CORE2 D compset and C ecosystem functionality
	
M      create_clone
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Components/csm_share.buildlib
M      ccsm_utils/Testlists/franklin.pretag
M      ccsm_utils/Testlists/jaguar.pretag
M      ccsm_utils/Testlists/kraken.pretag
M      ccsm_utils/Testlists/bluefire.pretag
	
================================================================================
Originator: mvertens
Date: Apr 06 2009
Model: scripts
Version: scripts4_090406
One-line summary: added CORE2 and CPLHIST to datm mode

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_definition.xml
	
================================================================================
Originator: tcraig
Date: Apr 05 2009
Model: scripts
Version: scripts4_090405
One-line summary: update testlist

 - update Test lists
 - update TCOST variable for ERP and ERT
 - update memory/tput test tag info

+++ ccsm_utils/Tools/testcase_end      
+++ ccsm_utils/Testcases/config_tests.xml
+++ ccsm_utils/Testlists/franklin.pretag
+++ ccsm_utils/Testlists/jaguar.posttag 
+++ ccsm_utils/Testlists/jaguar.pretag  
+++ ccsm_utils/Testlists/kraken.pretag  
+++ ccsm_utils/Testlists/bluefire.posttag
+++ ccsm_utils/Testlists/bluefire.pretag 

================================================================================
Originator: tcraig
Date: Apr 03 2009
Model: scripts
Version: scripts4_090403
One-line summary: merge in new_env scripts, pio changes, etc

- merge pioscr04_newenv08_scripts4_090402
  - updates for pio usage in driver (pioscr branch)
    - add PIO_CONFIG_OPTS env variable
    - update Macros files
  - bring new scripts branch with new xml files for machine pes on trunk (new_env branch)
- fix jaguar showproj again
- remove .rp files from driver st archive
- update franklin inputdata path
- add link_dirtree scripts
- update atlas scripts from art
- add clean_build scripts
- update Macros and SMP checking

M      create_clone
M      switch_machine
A  +   pes_file_sample.xml
D      pes_file_sample
M      create_test
M      ccsm_utils/Tools/configure
M      ccsm_utils/Tools/ccsm_buildexe.csh
M      ccsm_utils/Tools/check_input_data
M      ccsm_utils/Tools/ccsm_check_lockedfiles
M      ccsm_utils/Tools/archive_metadata.sh
M      ccsm_utils/Tools/st_archive.sh
A      ccsm_utils/Tools/clean_build
D      ccsm_utils/Tools/ccsm_setpes
M      ccsm_utils/Tools/xml2env
M      ccsm_utils/Tools/ccsm_getenv
M      ccsm_utils/Tools/xmlchange
M      ccsm_utils/Tools/generate_batch.csh
M      ccsm_utils/Case.template/ConfigCase.pm
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Components/pio.buildlib
M      ccsm_utils/Machines/Macros.surveyor
M      ccsm_utils/Machines/Macros.lightning_path
D      ccsm_utils/Machines/pes_setup.lightning_pgi
M      ccsm_utils/Machines/Macros.intrepid
D      ccsm_utils/Machines/pes_setup.midnight
M      ccsm_utils/Machines/Macros.schirra
D      ccsm_utils/Machines/pes_setup.kraken
D      ccsm_utils/Machines/pes_setup.lightning_intel
M      ccsm_utils/Machines/Macros.calgary_lahey
M      ccsm_utils/Machines/env_machopts.atlas
M      ccsm_utils/Machines/Macros.kraken
D      ccsm_utils/Machines/pes_setup.columbia
D      ccsm_utils/Machines/pes_setup.nyblue
D      ccsm_utils/Machines/pes_setup.surveyor
D      ccsm_utils/Machines/pes_setup.calgary_lahey
D      ccsm_utils/Machines/pes_setup.lightning_lahey
M      ccsm_utils/Machines/Macros.ubgl
D      ccsm_utils/Machines/pes_setup.jaguar
M      ccsm_utils/Machines/Macros.frost
M      ccsm_utils/Machines/Macros.bluefire
M      ccsm_utils/Machines/Macros.nyblue
D      ccsm_utils/Machines/pes_setup.intrepid
M      ccsm_utils/Machines/Macros.franklin
D      ccsm_utils/Machines/pes_setup.calgary_pgi
M      ccsm_utils/Machines/Macros.jaguar
M      ccsm_utils/Machines/env_machopts.jaguar
M      ccsm_utils/Machines/Macros.lightning_intel
D      ccsm_utils/Machines/pes_setup.frost
D      ccsm_utils/Machines/pes_setup
M      ccsm_utils/Machines/mkbatch.jaguar
D      ccsm_utils/Machines/pes_setup.ranger
M      ccsm_utils/Machines/Macros.midnight
M      ccsm_utils/Machines/Macros.atlas
M      ccsm_utils/Machines/Macros.lightning_lahey
D      ccsm_utils/Machines/pes_setup.ubgl
M      ccsm_utils/Machines/Macros.ranger
D      ccsm_utils/Machines/pes_setup.schirra
D      ccsm_utils/Machines/pes_setup.bluefire
M      ccsm_utils/Machines/Macros.calgary_pgi
M      ccsm_utils/Machines/Macros.lightning_pgi
D      ccsm_utils/Machines/pes_setup.lightning_path
D      ccsm_utils/Machines/pes_setup.atlas
D      ccsm_utils/Machines/pes_setup.franklin
M      ccsm_utils/Machines/config_machines.xml
A  +   ccsm_utils/Machines/config_pes.xml
M      ccsm_utils/Machines/Macros.columbia
A  +   ccsm_utils/Testcases/PEM_auto_pes_file
A  +   ccsm_utils/Testcases/PET_auto_pes_file
A  +   ccsm_utils/Testcases/SEQ_auto_pes_file
M      ccsm_utils/Testcases/PEM_script
M      ccsm_utils/Testcases/PET_script
M      ccsm_utils/Testcases/SEQ_script
M      create_newcase
A      link_dirtree

	
================================================================================
Originator: mvertens
Date: Apr 02 2009
Model: scripts
Version: scripts4_090402
One-line summary: bug fix for parsing of cam,clm namelist options

M      ccsm_utils/Case.template/ConfigCase.pm
	
================================================================================
Originator: mvr
Date: Mar 31 2009
Model: scripts
Version: scripts4_090331
One-line summary: lt archiver bug fix; pre/post tag test lists modified for 
	compset changes of previous tag

M      ccsm_utils/Tools/ccsm_l_archive.csh
- bug fix where data in partial files could be lost when mss was unreachable 
M      ccsm_utils/Testlists/jaguar.posttag
M      ccsm_utils/Testlists/kraken.posttag
M      ccsm_utils/Testlists/bluefire.cice2.pretag
M      ccsm_utils/Testlists/bluefire.posttag

	
================================================================================
Originator: mvertens
Date: Mar 30 2009
Model: scripts
Version: scripts4_090330
One-line summary: Updates for 1850 control and gx1v6 resolution

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_grid.xml
M      ccsm_utils/Case.template/config_definition.xml
	
================================================================================
Originator: erik
Date: Mar 25 2009
Model: scripts
Version: scripts4_090325
One-line summary: Remove CLM_ vars, add replacements, remove CLMNCEP mode, change
                  start-year for I cases, add 10x15_10x15 grid

Remove: CLM_DYNNDEP, CLM_BGC, CLM_CO2_TYPE, DIN_LOC_ROOT_CLMNCEP
Add:    CLM_NAMELIST_OPTS, CLM_FORCE_COLDSTART

Remove CLMNCEP mode for datm7 -- leaving CLM_QIAN mode.

Change start years for I cases to 2003 rather than 1948.

Add a 10x15 to 10x15 grid.

M      ccsm_utils/Case.template/config_compsets.xml ---- Remove deprecated CLM_ vars
          replace with CLM_CONFIG_OPTS, CLM_BLDNML_OPTS, and CLM_NAMELIST_OPTS.
M      ccsm_utils/Case.template/config_definition.xml -- Add/remove vars
M      ccsm_utils/Case.template/config_grid.xml -------- Add 10x15_10x15 grid
M      ccsm_utils/Machines/config_machines.xml --------- Remove DIN_LOC_ROOT_CLMNCEP
M      create_newcase ---- Fix bug 897 so comparisions are string rather than numeric

================================================================================
Originator: mvr
Date: Mar 23 2009
Model: scripts
Version: scripts4_090323
One-line summary: new location of showproj on jaguar, fix to pretag test lists

M      ccsm_utils/Machines/mkbatch.jaguar
M      ccsm_utils/Testlists/franklin.pretag
M      ccsm_utils/Testlists/jaguar.pretag
M      ccsm_utils/Testlists/kraken.pretag

================================================================================
Originator: mvertens
Date: Mar 20 2009
Model: scripts
Version: scripts4_090320
One-line summary: added gx1v6 support for data models 

M      ccsm_utils/Case.template/config_compsets.xml
       - extended compsets to make 1850 and 2000 explicit	
	
M      ccsm_utils/Case.template/config_grid.xml
       - added gx1v6 support thorughout
	 
M      ccsm_utils/Case.template/config_definition.xml
       - added CORE2 to DATM_MODE
       - added RX1 to DLND_RUNOFF_MODE
 	 
M      ccsm_utils/Testlists/franklin.pretag
M      ccsm_utils/Testlists/jaguar.pretag
M      ccsm_utils/Testlists/kraken.pretag
M      ccsm_utils/Testlists/bluefire.pretag
       - added more support for gx1v6 and 1850 and removed gx1v5 tests	
	
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Testlists/bluefire.cice1.pretag
       - Reduce aerosols (ntr_aero) to 3 and change one of the default
         CICE pretag tests to include gx1v6 grid. (dbailey)

================================================================================
Originator: mvr
Date: Mar 17 2009
Model: scripts
Version: scripts4_090317
One-line summary: removed '-overwrite sizemismatch' from msrcp arg list; tweaked 
	pretag test lists for jaguar, franklin, kraken

arg list to msrcp command included '-overwrite sizemismatch' which was not 
allowing model output to overwrite files already written to the mass store
	
M      ccsm_utils/Tools/ccsm_mswrite
M      ccsm_utils/Testlists/franklin.pretag
M      ccsm_utils/Testlists/jaguar.pretag
M      ccsm_utils/Testlists/kraken.pretag

================================================================================
Originator: mvertens
Date: Mar 12 2009
Model: scripts
Version: scripts4_090312
One-line summary: Added new F,B Track1 compsets and updated pretag tests

Updated pre-tag tests to do debug smoke tests for target production runs on 
target machines	
	
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Testlists/franklin.pretag
M      ccsm_utils/Testlists/jaguar.pretag
M      ccsm_utils/Testlists/kraken.pretag
M      ccsm_utils/Testlists/bluefire.pretag

================================================================================
Originator: erik
Date: Mar 10 2009
Model: scripts
Version: scripts4_090310
One-line summary: Add FCN compset and CLM_CONFIG_OPTS to scripts

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_definition.xml

================================================================================
Originator: mvr,mvertens,jwolfe
Date: Mar 03 2009
Model: scripts
Version: scripts4_090303
One-line summary: kraken machine mods, new compset for ramping co2

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/ConfigCase.pm
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Machines/env_machopts.kraken
M      ccsm_utils/Machines/mkbatch.kraken
M      ccsm_utils/Machines/config_machines.xml

================================================================================
Originator: tcraig
Date: Feb 26 2009
Model: scripts
Version: scripts4_090226
One-line summary: Add cpl.rp files to st archiving

M      ccsm_utils/Tools/st_archive.sh

================================================================================
Originator: dbailey
Date: Feb 24 2009
Model: scripts
Version: scripts4_090224
One-line summary: Change default for CICE aerosol deposition and
                  make g15 the default for E compset pretag tests.

M      ccsm_utils/Tools/st_archive.sh
M      ccsm_utils/Tools/archiving/st_archive.sh
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Testlists/bluefire.cice1.pretag
M      ccsm_utils/Testlists/bluefire.pretag

================================================================================

Originator: njn01
Date: Feb 20 2009
Model: scripts
Version: scripts4_090220
One-line summary: Add support for gx1v6 B compset

M           14500   ccsm_utils/Case.template/config_grid.xml
M           14500   ccsm_utils/Testlists/franklin.pretag
M           14500   ccsm_utils/Testlists/jaguar.posttag
M           14500   ccsm_utils/Testlists/jaguar.pretag
M           14500   ccsm_utils/Testlists/bluefire.pop2.pretag
M           14500   ccsm_utils/Testlists/bluefire.posttag
M           14500   ccsm_utils/Testlists/bluefire.pretag
M           14500   ChangeLog
Status against revision:  14500

================================================================================
Originator: jwolfe
Date: Feb 10 2009
Model: scripts
Version: scripts4_090210
One-line summary: update mapping files between T62 and tx1v1

M      ccsm_utils/Case.template/config_grid.xml

================================================================================
Originator: mvertens,mvr
Date: Jan 30 2009
Model: scripts
Version: scripts4_090130
One-line summary: various bug fixes

M      ccsm_utils/Case.template/config_compsets.xml
- fix to several compsets wrt cice_mode, cice_prestype
M      ccsm_utils/Case.template/config_grid.xml
- fix to several grid specifications wrt nx,ny settings
M      ccsm_utils/Case.template/config_definition.xml
- fix to list of valid values wrt to cice_mode, cice_prestype
M      ccsm_utils/Testlists/franklin.pretag
- fix to typo in test specification


	
================================================================================
Originator: mvertens
Date: Jan 26 2009
Model: scripts
Version: scripts4_090126
One-line summary:  Introduced changes necessary to interact with new cice configure

M      ccsm_utils/Case.template/config_grid.xml
       - put in necessary grid info for all cam grids (even if they are not
	 used in the compsets) to have only one place holder for these grids
M      ccsm_utils/Case.template/config_definition.xml
       - introduced CICE_CONFIG_OPTS and CICE_PRESTYPE that will be leveraged
	 by the new cice confgure utility
M      create_newcase
       - copied config_grid.xml into $caseroot Tools/ directory	

Note that this scripts tag is only compatible with cice4_0_20090126 or later	
	
================================================================================
Originator: njn01
Date: Jan 23 2009
Model: scripts
Version: scripts4_090123
One-line summary:  Add pop2 pre-tag testlist for use on bluefire

A               0   ccsm_utils/Testlists/bluefire.pop2.pretag
M           14033   ChangeLog

================================================================================
Originator: njn01
Date: Jan 22 2009
Model: scripts
Version: scripts4_090122
One-line summary:  Add default pe-layout for gx1v5 C compset on bluefire

M           14032   ccsm_utils/Machines/pes_setup.bluefire
M           14032   ChangeLog

================================================================================
Originator: tcraig
Date: Jan 20 2009
Model: scripts
Version: scripts4_090120
One-line summary:  gglc, add new drv namelist (EPS*), change I compsets

- add support for gglc in scripts
- add new drv namelist for domain checking tolerances (EPS*)
- change all I*Q compsets to I* in compsets and test lists
- add AG,BG,EG,FG,IG,XG compsets to support glc development
- decrease default eps_agrid tolerance for domain checking for I cases

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Machines/pes_setup
M      ccsm_utils/Machines/pes_setup.bluefire
M      ccsm_utils/Testlists/jaguar.posttag
M      ccsm_utils/Testlists/kraken.posttag
M      ccsm_utils/Testlists/testlist.port
M      ccsm_utils/Testlists/lightning.posttag
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/bluefire.pretag
M      create_newcase
	
================================================================================
Originator: mvr, tcraig
Date: Jan 12 2009
Model: scripts
Version: scripts4_090112
One-line summary:  archiving bugfix, update timing tool for drvseq3_0_04

M      ccsm_utils/Tools/timing/getTiming.pl
  - update timing tool for runoff resequencing in driver
M      ccsm_utils/Tools/st_archive.sh
  - fix for bug introduced at scripts4_081216b where a copy of the latest cam & 
    clm history files should get to the restart dir and the respective hist dirs 
	
================================================================================
Originator: tcraig
Date: Jan 7 2009
Model: scripts
Version: scripts4_090108
One-line summary:  fix setenv SMP code

M   ccsm_utils/Tools/ccsm_buildexe.csh
	
================================================================================
Originator: tcraig
Date: Jan 7 2009
Model: scripts
Version: scripts4_090107
One-line summary:  merge b07br02 tag to trunk, templates update, bug fixes, features

b07br02_scripts4_081226	
- Buildconf changes
- use RUNDIR everywhere, set it to $EXEROOT/run by default
- separate buildexe, buildnml, prestage.  no more run of buildexe from .run
  script, get rid of BLDTYPE, SETBLD
- build log files only created with buildexe
- SMP status checking, add SMP_BUILD, SMP_VALUE env variables
- switch_machine fixes (may still be broken)\
- change Makefile, mkDepends, mkSrcfiles to be from CASETOOLS
- in all templates, use CASEBUILD, RUNDIR, remove NTASK and NTHRD checks,
  remove lib existance check, remove copy of mkDepends and mkSrcfiles, remove
  copy of .mod files, remove SMP flag sent to gmake
- remove special pop rules for obj/compile location from Macros
- fix testcases for new scripts implementation
- update cice testlists
- update configure -clean*
- move camdom grid check to configure
- minor refactor of locking files
- remove generation of env_derived (was temporary)

b07br01_scripts4_081226	
- add $CASETOOLS and $CASEBUILD to env variables, replace $CASEROOT/Tools
  everywhere.
- move long and short term archive scripts to CASETOOLS and use them there
- rename check-input-data to check_input_data and refactor options, add
  svn check (from chris), improve performance, and other changes.
- fix RESUBMIT problem (bug #867)
- update exact restart compare tool, now checks both history and log files
- add save of timing file in ccsm_baselines for ERT test
- add env_derived file (temporary)
- fix permission issues on baseline saves
- add DOUT_L_HTAR env variable, long term archiver uses this now
- fix coupling frequency of G cases (bug #863)
- add tx1v1 grid, shortname s11
- remove DIN_LOC_ROOT_USER, CCSM_RUNLEN env varaibles
- update valid values info in STOP_OPTION variable (bug #790)
- make INC_NETCDF variable optional in all Macros files (bug #872)
- comment out module list in all env_machopts files
- add cost estimate env varaibles
- update bluefire and jaguar batch time limits based on ESTCOST variable so tests
  don't timeout.
- fix Macros.atlas for mct and pio flags
- add cice test list for bluefire
- fix tmp run file problem in create_suite when running in collections
- fix f02_g15 config_grid.xml typo
- fix kraken testlist C_ECOSYS to CECO and add f05_t12 test to franklin
- rename restart_compare to check_exactrestart

M      ccsm_utils/Build/Makefile
M      ccsm_utils/Tools/configure
A  +   ccsm_utils/Tools/ccsm_buildexe.csh
M      ccsm_utils/Tools/ccsm_s_archive.csh
A  +   ccsm_utils/Tools/check_input_data
M      ccsm_utils/Tools/ccsm_check_lockedfiles
M      ccsm_utils/Tools/archive_metadata.sh
M      ccsm_utils/Tools/ccsm_postrun.csh
M      ccsm_utils/Tools/testcase_env.csh
D      ccsm_utils/Tools/restart_compare.pl
M      ccsm_utils/Tools/ccsm_setpes
M      ccsm_utils/Tools/generate_resolved.csh
A  +   ccsm_utils/Tools/check_exactrestart.pl
D      ccsm_utils/Tools/ccsm_build.csh
D      ccsm_utils/Tools/check-input-data
A  +   ccsm_utils/Tools/ccsm_buildnml.csh
A  +   ccsm_utils/Tools/ccsm_prestage.csh
M      ccsm_utils/Tools/ccsm_getenv
M      ccsm_utils/Tools/testcase_end
M      ccsm_utils/Tools/generate_batch.csh
M      ccsm_utils/Tools/testcase_setup.csh
M      ccsm_utils/Tools/ccsm_l_archive.csh
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_grid.xml
M      ccsm_utils/Case.template/ConfigCase.pm
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Components/xlnd.template
M      ccsm_utils/Components/xocn.template
M      ccsm_utils/Components/sglc.template
M      ccsm_utils/Components/xice.template
M      ccsm_utils/Components/xatm.template
M      ccsm_utils/Components/xglc.template
M      ccsm_utils/Components/mct.buildlib
M      ccsm_utils/Components/csm_share.buildlib
M      ccsm_utils/Components/slnd.template
M      ccsm_utils/Components/pio.buildlib
M      ccsm_utils/Components/socn.template
M      ccsm_utils/Components/sice.template
M      ccsm_utils/Components/satm.template
M      ccsm_utils/Machines/env_machopts.columbia
M      ccsm_utils/Machines/Macros.surveyor
M      ccsm_utils/Machines/Macros.lightning_path
M      ccsm_utils/Machines/env_machopts.lightning_path
M      ccsm_utils/Machines/mkbatch.lightning_intel
M      ccsm_utils/Machines/Macros.intrepid
M      ccsm_utils/Machines/mkbatch.ubgl
M      ccsm_utils/Machines/env_machopts.lightning_intel
M      ccsm_utils/Machines/mkbatch.bluefire
M      ccsm_utils/Machines/mkbatch.lightning_path
M      ccsm_utils/Machines/mkbatch.franklin
M      ccsm_utils/Machines/mkbatch.calgary_lahey
M      ccsm_utils/Machines/mkbatch.lightning_lahey
M      ccsm_utils/Machines/Macros.schirra
M      ccsm_utils/Machines/Macros.calgary_lahey
M      ccsm_utils/Machines/mkbatch.calgary_pgi
M      ccsm_utils/Machines/env_machopts.lightning_lahey
M      ccsm_utils/Machines/Macros.kraken
M      ccsm_utils/Machines/mkbatch.frost
M      ccsm_utils/Machines/env_machopts.lightning_pgi
M      ccsm_utils/Machines/mkbatch.midnight
M      ccsm_utils/Machines/Macros.ubgl
M      ccsm_utils/Machines/pes_setup.jaguar
M      ccsm_utils/Machines/Macros.bluefire
M      ccsm_utils/Machines/Macros.frost
M      ccsm_utils/Machines/Macros.nyblue
M      ccsm_utils/Machines/mkbatch.kraken
M      ccsm_utils/Machines/Macros.franklin
M      ccsm_utils/Machines/Macros.jaguar
M      ccsm_utils/Machines/Macros.lightning_intel
M      ccsm_utils/Machines/mkbatch.columbia
M      ccsm_utils/Machines/mkbatch.schirra
M      ccsm_utils/Machines/env_machopts.schirra
M      ccsm_utils/Machines/mkbatch.jaguar
M      ccsm_utils/Machines/Macros.midnight
M      ccsm_utils/Machines/Macros.atlas
M      ccsm_utils/Machines/env_machopts.midnight
M      ccsm_utils/Machines/Macros.lightning_lahey
M      ccsm_utils/Machines/mkbatch.lightning_pgi
M      ccsm_utils/Machines/Macros.ranger
M      ccsm_utils/Machines/Macros.calgary_pgi
M      ccsm_utils/Machines/Macros.lightning_pgi
M      ccsm_utils/Machines/config_machines.xml
M      ccsm_utils/Machines/Macros.columbia
M      ccsm_utils/Testcases/ERH_script
M      ccsm_utils/Testcases/PEM_script
M      ccsm_utils/Testcases/config_tests.xml
M      ccsm_utils/Testcases/APT_script
M      ccsm_utils/Testcases/ERP_script
M      ccsm_utils/Testcases/ERB_script
M      ccsm_utils/Testcases/PET_script
M      ccsm_utils/Testcases/SEQ_script
M      ccsm_utils/Testcases/ERS_script
M      ccsm_utils/Testcases/ERT_script
M      ccsm_utils/Testlists/franklin.pretag
M      ccsm_utils/Testlists/kraken.posttag
A  +   ccsm_utils/Testlists/bluefire.cice1.pretag
A  +   ccsm_utils/Testlists/bluefire.cice2.pretag
M      switch_machine
M      check_rundb
M      create_newcase
M      create_suite
	
================================================================================
Originator: tcraig
Date: Fri Dec 26 2008
Model: scripts
Version: scripts4_081226
One-line summary:  glc updates

- add glc component
- update timing for glc
- update pe layouts for glc
- set default glc component to sglc in all compsets except X
- add glc grid to resolutions
- add GLC_NEC env variable to set glc elevation classes

M      pes_file_sample
M      README
M      create_test
M      ccsm_utils/Build/Makefile
X      ccsm_utils/Tools/perl5lib
M      ccsm_utils/Tools/st_archive.sh
M      ccsm_utils/Tools/ccsm_setpes
M      ccsm_utils/Tools/generate_resolved.csh
M      ccsm_utils/Tools/timing/getTiming.pl
M      ccsm_utils/Tools/taskmaker.pl
M      ccsm_utils/Tools/ccsm_l_archive.csh
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_grid.xml
M      ccsm_utils/Case.template/ConfigCase.pm
M      ccsm_utils/Case.template/config_definition.xml
A  +   ccsm_utils/Components/sglc.template
A  +   ccsm_utils/Components/xglc.template
M      ccsm_utils/Components/csm_share.buildlib
M      ccsm_utils/Machines/pes_setup.lightning_pgi
M      ccsm_utils/Machines/Macros.cppdefs
M      ccsm_utils/Machines/pes_setup.midnight
M      ccsm_utils/Machines/pes_setup.kraken
M      ccsm_utils/Machines/pes_setup.lightning_intel
M      ccsm_utils/Machines/Macros.calgary_lahey
M      ccsm_utils/Machines/Macros.kraken
M      ccsm_utils/Machines/pes_setup.columbia
M      ccsm_utils/Machines/pes_setup.surveyor
M      ccsm_utils/Machines/pes_setup.nyblue
M      ccsm_utils/Machines/pes_setup.calgary_lahey
M      ccsm_utils/Machines/pes_setup.lightning_lahey
M      ccsm_utils/Machines/pes_setup.jaguar
M      ccsm_utils/Machines/pes_setup.intrepid
M      ccsm_utils/Machines/pes_setup.calgary_pgi
M      ccsm_utils/Machines/Macros.franklin
M      ccsm_utils/Machines/Macros.jaguar
M      ccsm_utils/Machines/pes_setup.frost
M      ccsm_utils/Machines/mkbatch.schirra
M      ccsm_utils/Machines/pes_setup
M      ccsm_utils/Machines/pes_setup.ranger
M      ccsm_utils/Machines/Macros.midnight
M      ccsm_utils/Machines/pes_setup.ubgl
M      ccsm_utils/Machines/Macros.ranger
M      ccsm_utils/Machines/pes_setup.schirra
M      ccsm_utils/Machines/pes_setup.bluefire
M      ccsm_utils/Machines/pes_setup.lightning_path
M      ccsm_utils/Machines/pes_setup.atlas
M      ccsm_utils/Machines/pes_setup.franklin
M      ccsm_utils/Testcases/PEM_script
M      ccsm_utils/Testcases/PET_script
M      ccsm_utils/Testcases/SEQ_script
M      create_newcase

================================================================================
Originator: mvr
Date: Thr Dec 18 2008
Model : scripts
Version: scripts4_081218b
One-line summary: fix to branch and hybrid test scripts to account for changes
in how restart files are bundled

M      ccsm_utils/Testcases/ERH_script
M      ccsm_utils/Testcases/ERB_script

================================================================================
Originator: mvr
Date: Thr Dec 18 2008
Model : scripts
Version: scripts4_081218
One-line summary: fix to jaguar and franklin lt archiver submittal

M      mkbatch.franklin
M      mkbatch.jaguar
M      config_machines.xml
- will now disallow multiple submittals of lt archiver

================================================================================
Originator: mvertens
Date: Tue Dec 17 2008
Model : scripts
Version: scripts4_081217
One-line summary: new mapping files using patch interpolation

The following two new mapping files were introduced
M      ccsm_utils/Case.template/config_grid.xml
       -  MAP_A2OS_FILE="map_fv1.9x2.5_to_gx1v5_patch.081210.nc"
       -  MAP_A2OS_FILE="map_T62_to_gx1v5_patch.081125.nc"

================================================================================

Originator: mvr
Date: Tue Dec 16 2008
Model : scripts
Version: scripts4_081216b
One-line summary: lt and st archive mods and other bug fixes

M      ccsm_utils/Tools/archive_metadata.sh
- updated to work with changes to case directory layout
M      ccsm_utils/Tools/st_archive.sh
- restart files now contained in subdirectory rather than tarball (needed for hires)
M      ccsm_utils/Tools/ccsm_l_archive.csh
- bug fixes and improved handling of partial tarring of history files
M      ccsm_utils/Machines/mkbatch.franklin
M      ccsm_utils/Machines/mkbatch.jaguar
M      ccsm_utils/Machines/config_machines.xml
- small machine specific bug fixes


================================================================================
Originator: tcraig
Date: Tue Dec 16 2008
Model : scripts
Version: scripts4_081216
One-line summary: pio/mct build update

 - update scripts related to new mct/pio build (scripts4_081210)

M      ccsm_utils/Tools/ccsm_build.csh
M      ccsm_utils/Machines/Macros.surveyor
M      ccsm_utils/Machines/Macros.intrepid
M      ccsm_utils/Machines/env_machopts.kraken
M      ccsm_utils/Machines/Macros.ubgl
M      ccsm_utils/Machines/Macros.frost
M      ccsm_utils/Machines/Macros.midnight
M      ccsm_utils/Machines/Macros.calgary_pgi

================================================================================

Originator: tcraig
Date: Mon Dec 15 2008
Model : scripts
Version: scripts4_081215
One-line summary: pio/mct build update

 - update scripts related to new mct/pio build (scripts4_081210)

M    ccsm_utils/Machines/Macros.lightning_path
D    ccsm_utils/Machines/Macros.bangkok.pgf90
M    ccsm_utils/Machines/Macros.calgary_lahey
D    ccsm_utils/Machines/Macros.calgary.pgf90
D    ccsm_utils/Machines/Macros.bangkok
M    ccsm_utils/Machines/Macros.lightning_intel
D    ccsm_utils/Machines/Macros.calgary
M    ccsm_utils/Machines/Macros.lightning_lahey
M    ccsm_utils/Machines/Macros.ranger
M    ccsm_utils/Machines/Macros.calgary_pgi
M    ccsm_utils/Machines/Macros.lightning_pgi

================================================================================

Originator: tcraig
Date: Mon Dec 15 2008
Model : scripts
Version: scripts4_081214
One-line summary: pio/mct build update

 - update scripts related to new mct/pio build (scripts4_081210)

M      ccsm_utils/Components/mct.buildlib
M      ccsm_utils/Components/pio.buildlib
M      ccsm_utils/Machines/Macros.surveyor
M      ccsm_utils/Machines/Macros.lightning_path
M      ccsm_utils/Machines/Macros.intrepid
D      ccsm_utils/Machines/env_machopts.lightning
M      ccsm_utils/Machines/Macros.schirra
M      ccsm_utils/Machines/Macros.calgary_lahey
M      ccsm_utils/Machines/Macros.kraken
M      ccsm_utils/Machines/Macros.bangkok
M      ccsm_utils/Machines/env_machopts.kraken
M      ccsm_utils/Machines/Macros.frost
M      ccsm_utils/Machines/Macros.nyblue
M      ccsm_utils/Machines/Macros.franklin
M      ccsm_utils/Machines/env_machopts.franklin
D      ccsm_utils/Machines/mkbatch.lightning
M      ccsm_utils/Machines/Macros.jaguar
M      ccsm_utils/Machines/env_machopts.jaguar
M      ccsm_utils/Machines/Macros.lightning_intel
D      ccsm_utils/Machines/Macros.lightning
M      ccsm_utils/Machines/Macros.midnight
M      ccsm_utils/Machines/Macros.atlas
M      ccsm_utils/Machines/Macros.lightning_lahey
D      ccsm_utils/Machines/pes_setup.lightning
M      ccsm_utils/Machines/Macros.ranger
M      ccsm_utils/Machines/Macros.calgary_pgi
M      ccsm_utils/Machines/Macros.lightning_pgi
M      ccsm_utils/Machines/Macros.columbia
	
================================================================================
Originator: tcraig
Date: Fri Dec 12 2008
Model : scripts
Version: scripts4_081212
One-line summary: machine updates

 - update intrepid scripts
 - update calgary scripts
 - change dead and stub make file locations, use case versions
 - improve hist_compare tool
 - increase MPICH_UNEX_BUFFER_SIZE on kraken for higher resolution cases
 - move dead_share build out of multiple dead models and into csm_share build
   to eliminate multiple instances of the same build code.  this means the
   dead_share code will be built in all cases (like dshr).

M      ccsm_utils/Tools/hist_compare.csh
M      ccsm_utils/Components/xlnd.template
M      ccsm_utils/Components/xocn.template
M      ccsm_utils/Components/xice.template
M      ccsm_utils/Components/xatm.template
M      ccsm_utils/Components/csm_share.buildlib
M      ccsm_utils/Components/slnd.template
M      ccsm_utils/Components/socn.template
M      ccsm_utils/Components/sice.template
M      ccsm_utils/Components/satm.template
M      ccsm_utils/Machines/Macros.calgary_lahey
M      ccsm_utils/Machines/mkbatch.calgary_pgi
M      ccsm_utils/Machines/env_machopts.kraken
M      ccsm_utils/Machines/pes_setup.intrepid
M      ccsm_utils/Machines/mkbatch.intrepid
M      ccsm_utils/Machines/config_machines.xml
A      ccsm_utils/Testlists/calgary_lahey.pretag

================================================================================
Originator: Jedwards
Date: Wed Dec 10 2008
Model : scripts
Version: scripts4_081210
One-line summary: Refactored mct and pio build process to be more consistant

Note - this tag WILL require a cam tag of cam3_6_24 or later to have a compatible
	pio build.   I have tried to make the correct modifications in Macros files
	but some builds may need further modifications in the Macros.* file.

        ccsm_utils/Build/Makefile
        ccsm_utils/Components/mct.buildlib
        ccsm_utils/Components/pio.buildlib
        ccsm_utils/Machines/Macros.bangkok
        ccsm_utils/Machines/Macros.bluefire
        ccsm_utils/Machines/Macros.calgary_lahey
        ccsm_utils/Machines/Macros.calgary_pgi
        ccsm_utils/Machines/Macros.frost
        ccsm_utils/Machines/Macros.ubgl
        ccsm_utils/Machines/Macros.atlas
        ccsm_utils/Machines/Macros.jaguar
        ccsm_utils/Tools/ccsm_build.csh



	
	
================================================================================
Originator: mvertens
Date: Sun Dec 7 2008
Model: scripts
Version: scripts4_081207
One-line summary:  modified scripts to interact with cam configure more seamlessly
	
Note that this tag WILL require a cam tag of cam3_6_23 or later to have a compatible
cam.cpl7.template	
	
M      create_test
       - fixed bug when test id was not specified

M      ccsm_utils/Case.template/config_compsets.xml
	- significantly cleaned up the compsets to reflect the only new clm
	 supported forcing data (CLMNCEP forcing will not longer be used) 
	 the following compsets were removed
	    -<compset NAME="I_1948_2004" SHORTNAME="I4804"
            -<compset NAME="I_CASA" SHORTNAME="ICASA"
            -<compset NAME="I_CASA_19482004" SHORTNAME="ICASA4804"
            -<compset NAME="I_CN" SHORTNAME="ICN"
            -<compset NAME="I_CN_1984_2004" SHORTNAME="ICN4804"
       - removed following bgc-related compsets that are not being tested and that are broken in cam.
	 when this functionality is needed again, cam will create use cases to set these up automatically. 
	 Discussed this with Keith Lindsay and Brian Eaton on 12/5/2008. 
            -<compset NAME="B_CO2A" SHORTNAME="BADT"
            -<compset NAME="B_CO2B_CN" SHORTNAME="BBN"
            -<compset NAME="B_CO2C_CN_dd" SHORTNAME="BCND"
            -<compset NAME="B_CO2C_CN_dp" SHORTNAME="BCNM"
            -<compset NAME="B_CO2C_CN_pp" SHORTNAME="BCNP"
            -<compset NAME="B_RAMP_CO2" SHORTNAME="BCR"
            -<compset NAME="F_CO2A" SHORTNAME = "FADT"
            -<compset NAME="F_CO2B_CASA" SHORTNAME="FBCA"
            -<compset NAME="F_CO2B_CN" SHORTNAME="FBN"
	
M      ccsm_utils/Case.template/config_grid.xml
	- CAM_DTIME and RTM_NSTEPS were removed
	
M      ccsm_utils/Case.template/config_definition.xml
	- the following variables were removed from config_definition.xml
	  note that cam_dtime is determined in cam.cpl7.template from ATM_NCPL
          note that clm_dtime and rtm_nsteps is determined in clm.cpl7.template from LND_NCPL
             -<entry id="CLM_DTIME" 
             -<entry id="CICE_DTIME" 
             -<entry id="RTM_NSTEPS" 
             -<entry id="iradsw" 
        - the following cam-related variables were removed from config_definition.xml and replaced with
	  CAM_CONFIG_OPTS
             -<entry id="CAM_PLEV" 
             -<entry id="CAM_NADV" 
             -<entry id="CAM_PROG_AERO" 
             -<entry id="CAM_CHEM" 
        - the followin cam-related variables were removed from config_defintion.xml since most of them
	  were no longer functional or working correctly - they correspond to the removal of the related compsets
	  above and will be replaced by cam use cases 
             -<entry id="CAM_CO2_TYPE" 
             -<entry id="CAM_CO2_READ" 
             -<entry id="CAM_CO2FLUX_READ_FUEL" 
             -<entry id="CAM_CO2FLUX_READ_OCN" 
             -<entry id="CAM_RAMP_CO2_START_YMD" 

M      ccsm_utils/Testlists/jaguar.posttag
M      ccsm_utils/Testlists/kraken.posttag
M      ccsm_utils/Testlists/testlist.port
M      ccsm_utils/Testlists/lightning.posttag
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/bluefire.pretag
       - the above tests were modified to reflect the removal of the above compsets and renaming of the I compset	  
	
================================================================================
Originator: jwolfe
Date: Fri Dec 5 2008
Model: scripts
Version: scripts4_081205
One-line summary:  added env variable to increase stack size on XT4s

Added an environment variable setting to increase stack size on XT4 machines,
necessary to run CICE in a threaded context.

M      ccsm_utils/Machines/env_machopts.kraken
M      ccsm_utils/Machines/env_machopts.franklin
M      ccsm_utils/Machines/env_machopts.jaguar

================================================================================
Originator: mvertens
Date: Wed Dec 3 2008
Model: scripts
Version: scripts4_081203
One-line summary:  added status of compset to create_newcase 

Also modified test list for test that are currently failing and will not be
changed soon - this is reflected in the STATUS variable in config_compsets.xml
	
M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Testlists/jaguar.posttag
M      ccsm_utils/Testlists/kraken.posttag
M      ccsm_utils/Testlists/lightning.posttag
M      ccsm_utils/Testlists/bluefire.posttag
M      create_newcase
	
================================================================================
Originator: tcraig
Date: Mon Nov 24 2008
Model: scripts
Version: scripts4_081124
One-line summary:  update to frost and cam options

- remove use of trop_mozart_ghg_paero in CAM_CHEM compset settings
- update CAM_CHEM valid options
- updates to frost scripts from jdennis

M      ccsm_utils/Case.template/config_compsets.xml
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Machines/mkbatch.frost
M      ccsm_utils/Machines/Macros.frost
	
================================================================================
Originator: tcraig
Date: Wed Nov 19 2008
Model: scripts
Version: scripts4_081119b
One-line summary:  merge beta4 branch tag

svn merge $SVNREPO/scripts/trunk_tags/scripts4_081009b $SVNREPO/scripts/branch_tags/beta_scripts4_081009b_tags/beta5_scripts4_081009b

- change CCSM_CO2_PPMV default to 379.0
- add memory tool on bluefire binary launch
- remove DUST cppdef
- update archiver for pop restarts

U    ccsm_utils/Case.template/config_compsets.xml
U    ccsm_utils/Case.template/config_definition.xml
U    ccsm_utils/Machines/mkbatch.bluefire
U    ccsm_utils/Machines/Macros.cppdefs
U    ccsm_utils/Tools/st_archive.sh

================================================================================
Originator: tcraig
Date: Wed Nov 19 2008
Model: scripts
Version: scripts4_081119
One-line summary:  update machines and features

- Add support for lightning_path, lightning_pgi, lightning_lahey,
  lightning_intel, calgary_lahey, calgary_pgi, intrepid, and ubgl
- Update scripts for surveyor and atlas
- Add tests PEM (compare changing component pes, all mpi), PET (compare changing
  component pes, mixed mpi/openmp), and SEQ (compare all concurrent vs all seq)
- Delete old generate_cam and generate_pop_cice tools
- remove regress option in testname in create_test
- add support for SEQ, PEM, and PET tests in create_test (needed to set initial pes)
  and configure (do not clean .test file ever)
- add support for env_build file to be locked and force rebuild.  modify
  behavior of env_mach_decomp so it forces rebuild, not reconfigure.
- add call to ccsm_getenv at the start of the ccsm_postrun script so values
  are updated based on changes a user might have made.
- change formatting of env_mach_pes file to allow ccsm_sedfile to work
- refactor ccsm_build.csh for cleanup and separate buildexe and buildnml
  control
- update ccsm_sedfile logic and grep for more flexibility of env file format
- update lightning pretag tests
- update bluefire.post tests to add SEQ, PEM, and PET tests

Add lots of new files names (not documented here) plus
M      create_test
M      ccsm_utils/Tools/configure
D      ccsm_utils/Tools/generate_cice_pop_decomp
M      ccsm_utils/Tools/ccsm_check_lockedfiles
M      ccsm_utils/Tools/ccsm_postrun.csh
M      ccsm_utils/Tools/ccsm_setpes
M      ccsm_utils/Tools/ccsm_sedfile
D      ccsm_utils/Tools/ConfigInfo.xml
M      ccsm_utils/Tools/ccsm_build.csh
D      ccsm_utils/Tools/generate_cam_decomp
M      ccsm_utils/Tools/generate_batch.csh
M      ccsm_utils/Components/mct.buildlib
M      ccsm_utils/Machines/Macros.calgary.pgf90
M      ccsm_utils/Machines/pes_setup.surveyor
M      ccsm_utils/Machines/Macros.atlas
M      ccsm_utils/Machines/pes_setup.lightning
M      ccsm_utils/Machines/config_machines.xml
M      ccsm_utils/Testcases/config_tests.xml
M      ccsm_utils/Testlists/lightning.pretag
M      ccsm_utils/Testlists/bluefire.posttag
	
================================================================================
Originator: tcraig
Date: Wed Nov 5 2008
Model: scripts
Version: scripts4_081105
One-line summary:  update jaguar modules for upgrade problems

- add
  +  module swap xt-asyncpe xt-asyncpe/1.0c
  +  module swap xt-binutils-quadcore xt-binutils-quadcore/2.0.1
	
M    ccsm_utils/Machines/env_machopts.jaguar
	
================================================================================
Originator: tcraig
Date: Fri Oct 24 2008
Model: scripts
Version: scripts4_081024
One-line summary:  updates for revised timing lib

- fix bug in kraken, jaguar, franklin mkbatch script, divide syntax
- update perf_summary.pl for new timing output
- add TPROF alarm controls for intermediate timing files
- fix scriptsroot bug in create_test
- add creation of $RUNDIR/timing/checkpoints directory

M      create_test
M      ccsm_utils/Tools/timing/perf_summary.pl
M      ccsm_utils/Tools/generate_batch.csh
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Machines/mkbatch.kraken
M      ccsm_utils/Machines/mkbatch.franklin
M      ccsm_utils/Machines/mkbatch.jaguar
M      ccsm_utils/Machines/config_machines.xml
	
================================================================================
Originator: tcraig
Date: Wed Oct 22 2008
Model: scripts
Version: scripts4_081022
One-line summary:  update timing lib, other features

- updates for new timing lib; add timing directory under $RUNDIR for
  timing files, add perf_summary.pl tool to create summary timing file,
  update post-processing getTiming tools for new output, timing files
  written on a per mpi task basis.
- add CCSM_BASELINE to "mach_dout" group, modify generate and compare 
  implementation to support optional setting of env variable CCSM_BASELINE.
- add CCSM_CPRNC to "mach_dout" group,  points to a local copy of cprnc for 
  testing, update hist_compare to use CCSM_CPRNC if it's valid and available,
  otherwise, continue to use ncdump | split.
- modify location of cs.status and cs.submit scripts.  they will now
  go into the testroot directory set in create_suite instead of in the
  local scripts directory.
- fix usage typo in taskmaker.pl
- add copy of TestStatus.out to ccsm_baselines area when generating tests.
  provides history of test result plus provides tput and memory usage in
  ERT mode for use when comparing with the next tag.
- update ERT test to add summary of tput and memory usage as well as
  exact restart using the cpl history file and cprnc if available.
- add CHECK test status option for tput and memory
- new files moved to $CASEROOT/Tools; restart_compare.pl, hist_compare.pl, 
  getTiming.csh, getTiming.pl, perf_summary.pl
- add missing env variable back into config_definition.xml, PES_PER_NODE
- remove use of showproj on kraken and remove BSUB -A option to batch jobs.
- remove -DTHREADED_PTHREADS use in all Macros, add in -DTHREADED_OMP
- update tests for location of restart_compare and hist_compare tools
  (moved to $CASEROOT/Tools, was $CASEROOT).
- remove redundant ERT test in bluefire.posttag
- update r05_to_tx0.1v2 mapping file to map_r05_to_tx0.1v2_r500e1000_080620.nc

M      create_test
M      ccsm_utils/Tools/ccsm_postrun.csh
M      ccsm_utils/Tools/testcase_env.csh
M      ccsm_utils/Tools/ccsm_build.csh
M      ccsm_utils/Tools/timing/getTiming.pl
M      ccsm_utils/Tools/timing/getTiming.csh
A      ccsm_utils/Tools/timing/perf_summary.pl
M      ccsm_utils/Tools/hist_compare.csh
M      ccsm_utils/Tools/taskmaker.pl
M      ccsm_utils/Tools/testcase_end
M      ccsm_utils/Tools/generate_batch.csh
M      ccsm_utils/Tools/testcase_setup.csh
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Case.template/config_grid.xml
M      ccsm_utils/Machines/mkbatch.kraken
M      ccsm_utils/Machines/Macros.calgary
M      ccsm_utils/Machines/Macros.lightning
M      ccsm_utils/Machines/Macros.atlas
M      ccsm_utils/Machines/Macros.bangkok
M      ccsm_utils/Machines/Macros.ranger
M      ccsm_utils/Machines/config_machines.xml
M      ccsm_utils/Machines/Macros.columbia
M      ccsm_utils/Testcases/ERH_script
M      ccsm_utils/Testcases/config_tests.xml
M      ccsm_utils/Testcases/APT_script
M      ccsm_utils/Testcases/ERP_script
M      ccsm_utils/Testcases/ERB_script
M      ccsm_utils/Testcases/ERS_script
M      ccsm_utils/Testcases/ERT_script
M      ccsm_utils/Testlists/bluefire.posttag
M      create_newcase
M      create_suite

================================================================================
Originator: mvr,jwolfe
Date: Thu Oct 09 2008
Model: scripts
Version: scripts4_081009b
One-line summary:  fix to lt archive script for several machines where sourcing 
of the environment files was not happening; addition of grid info for 0.9x1.25 

M      ccsm_utils/Case.template/config_grid.xml
M      ccsm_utils/Machines/mkbatch.kraken
M      ccsm_utils/Machines/mkbatch.franklin
M      ccsm_utils/Machines/mkbatch.jaguar
M      ccsm_utils/Machines/mkbatch.midnight

	
	
================================================================================
Originator: erik
Date: Thu Oct 09 2008
Model: scripts
Version: scripts4_081009

Update version of perl5lib used. Use one with CAM updates in it, and fixed / expanded
unit tests.
	
================================================================================
Originator: mvertens
Date: Wed Oct 08 2008
Model: scripts
Version: scripts4_081008a
	
M  ccsm_utils/Case.template/config_compsets.xml
    - removed B_1990_TYPE1 and F_CAMDOM_TGP
    - removed CAM_CHEM and CAM_USECASE settings from B_PRESENT_DAY and F_CAMDOM
      (UNSET will now be defaulted for those values) 	
M  ccsm_utils/Testlists/lightning.posttag
    - removed F_AICE_SOM and added E	
	
================================================================================
Originator: tcraig
Date: Wed Oct 08 2008
Model: scripts
Version: scripts4_081008
One-line summary:  General scripts fixes

- fix Makefile for INCROOT dependency and -I
- remove INCROOT reference from all Macros references for -I
- add Deppath feature to Makefile
- modify ccsm_check_lockedfiles for Macros change.  Instead of requiring
  a reconfigure, now it resets BUILD_COMPLETE to FALSE, updates the
  locked Macros file, and continues. (bug 827)
- rename F_AICE_SOM to E_PRESENT_DAY, change tests from FAIS to E
- fix SMP flag bug in csm_share.buildlib
- update f05_g15 B pes settings on jaguar, franklin, and kraken for 
  mixed mpi/openmp on 2000 pes.
- remove special seq_timemgr_mod and histflds build rules in all Macros
- update kraken and franklin l_archive batch settings
- redirect stdout/err to a ccsm.log.$LID file for all machines except bluegene's
  (want to test first) that didn't have that implemented.
- fix typo in config_machines.xml for frost, nyblue, and surveyor
- update pretag test lists.  add a debug test to jaguar and add ERP to
  both jaguar and bluefire pretag tests.
- update create_newcase to write out both long and short compset name
  to CCSM_COMPSET variable.

	M      ccsm_utils/Build/Makefile
	M      ccsm_utils/Tools/ccsm_check_lockedfiles
	M      ccsm_utils/Case.template/config_compsets.xml
	M      ccsm_utils/Components/csm_share.buildlib
	M      ccsm_utils/Machines/pes_setup.jaguar
	M      ccsm_utils/Machines/Macros.bluefire
	M      ccsm_utils/Machines/Macros.frost
	M      ccsm_utils/Machines/Macros.surveyor
	M      ccsm_utils/Machines/Macros.nyblue
	M      ccsm_utils/Machines/Macros.bangkok.pgf90
	M      ccsm_utils/Machines/mkbatch.kraken
	M      ccsm_utils/Machines/Macros.franklin
	M      ccsm_utils/Machines/mkbatch.lightning
	M      ccsm_utils/Machines/Macros.jaguar
	M      ccsm_utils/Machines/Macros.calgary
	M      ccsm_utils/Machines/Macros.lightning
	M      ccsm_utils/Machines/mkbatch.franklin
	M      ccsm_utils/Machines/mkbatch.atlas
	M      ccsm_utils/Machines/Macros.schirra
	M      ccsm_utils/Machines/mkbatch.jaguar
	M      ccsm_utils/Machines/pes_setup.kraken
	M      ccsm_utils/Machines/Macros.midnight
	M      ccsm_utils/Machines/Macros.atlas
	M      ccsm_utils/Machines/Macros.calgary.pgf90
	M      ccsm_utils/Machines/Macros.kraken
	M      ccsm_utils/Machines/Macros.bangkok
	M      ccsm_utils/Machines/Macros.ranger
	M      ccsm_utils/Machines/pes_setup.franklin
	M      ccsm_utils/Machines/config_machines.xml
	M      ccsm_utils/Machines/Macros.columbia
	M      ccsm_utils/Testlists/franklin.pretag
	M      ccsm_utils/Testlists/jaguar.pretag
	M      ccsm_utils/Testlists/kraken.pretag
	M      ccsm_utils/Testlists/bluefire.posttag
	M      ccsm_utils/Testlists/bluefire.pretag
	M      create_newcase
	
================================================================================
Originator: erik
Date: Fri Oct 03 2008
Model: scripts
Version: scripts4_081003
One-line summary:  Update perl5lib to fix validate problem for generic namelists

Update the version of perl5lib to fix a validate problem for generic namelists.
This fixes bug 833 for cam in ccsm4_0_alpha37.

================================================================================
Originator: tcraig
Date: Wed Oct 01 2008
Model: scripts
Version: scripts4_081001
One-line summary:  Add midnight, Add H

- add support for midnight (SUN cluster at ARSC)
- fix f45_f45 setup in config_grid.xml
- mvr bug fixes for check-input-data
- add H compset, datm/slnd/cice/pop2
- update to franklin scripts, redirect stdout to ccsm log file

M      ccsm_utils/Tools/check-input-data
M      ccsm_utils/Case.template/config_grid.xml
M      ccsm_utils/Case.template/config_compsets.xml
A      ccsm_utils/Machines/pes_setup.midnight
A      ccsm_utils/Machines/Macros.midnight
A      ccsm_utils/Machines/env_machopts.midnight
M      ccsm_utils/Machines/config_machines.xml
A      ccsm_utils/Machines/mkbatch.midnight
M      ccsm_utils/Machines/mkbatch.franklin
A      ccsm_utils/Testlists/midnight.pretag

================================================================================
Originator: erik
Date: Tue Sep 30 09:58:38 MDT 2008
Model: scripts
Version: scripts4_080930
One-line summary:  Add some I cases, fix typos in I-cases, rename Config.pm to ConfigCase.pm

Move clmIfix02_scripts4_080928 to trunk.

Add some CASA and CLM_QIAN I cases. The CLM_QIAN I cases require a version of datm7 that
has the CLM_QIAN mode for DATM_MODE. Fix some minor typos in other I-cases. Add
DIN_LOC_ROOT_CLM_QIAN path for the CLM_QIAN cases on all machines. Add CLM_DEMAND
and CLM_BLD_NL_OPTIONS so you can pass options to the clm build-namelist script.

Get the pt1 grids working for CLM -- it requires a version of CLM and datm7 that have
this feature enabled.

Rename Config.pm to ConfigCase.pm as on some platforms perl would get confused and use
the wrong Config.pm. Also remove some perl warnings such as use $array[$index] rather
than @array[$index]. Also add some documentation to this module at the top, and move
some of the methods so they are considered private (and rename to prepend name with
an underscore).

Add an option to include the previous December when matching %ym for build_streams.
This required using a newer version of perl5lib in SVN_EXTERNAL_DIRECTORIES
(perl5lib_080924).

Add executable bit to scripts as on some platforms they wouldn't even run correctly
without it.

D      ccsm_utils/Case.template/Config.pm -------------- Rename to ConfigCase.pm
A  +   ccsm_utils/Case.template/ConfigCase.pm ---------- Add documentation, make some
                                                         methods private. Remove some
                                                         perl warnings.
M      ccsm_utils/Case.template/config_compsets.xml ---- Fix I4804 compset, add ICASA4804,
                                                         add CLM_QIAN compsets
M      ccsm_utils/Case.template/config_grid.xml -------- Fix pt1 grid.
M      ccsm_utils/Case.template/config_definition.xml -- Add CLM_DEMAND, CLM_BLD_NL_OPTIONS, 
                                                         and DIN_LOC_ROOT_CLMQIAN. Add
                                                         CLM_QIAN as option to DATM_MODE.
M      create_clone ------------------------------------ Change to ConfigCase.pm
M      switch_machine ---------------------------------- Change to ConfigCase.pm
M      ccsm_utils/Tools/xml2env ------------------------ Change to ConfigCase.pm
M      ccsm_utils/Tools/build_streams ------------------ Add option to include Dec. of previous year
M      ccsm_utils/Tools/xmlchange ---------------------- Change to ConfigCase.pm
M      create_newcase ---------------------------------- Change to ConfigCase.pm
M      ccsm_utils/Machines/config_machines.xml --------- Add DIN_LOC_CLM_QIAN
>>>>>>>>> Add executable bit to scripts that will be run
 M     ccsm_utils/Tools/ccsm_msread
 M     ccsm_utils/Tools/configure
 M     ccsm_utils/Tools/ccsm_s_archive.csh
 M     ccsm_utils/Tools/testcase_env.csh
 M     ccsm_utils/Tools/ccsm_mswrite
 M     ccsm_utils/Tools/ccsm_cpdata
 M     ccsm_utils/Tools/restart_compare.pl
 M     ccsm_utils/Tools/ccsm_setpes
 M     ccsm_utils/Tools/generate_resolved.csh
 M     ccsm_utils/Tools/ccsm_getfile
 M     ccsm_utils/Tools/create_production_test_readme
 M     ccsm_utils/Tools/ccsm_auto.csh
 M     ccsm_utils/Tools/taskmaker.pl
 M     ccsm_utils/Tools/ccsm_getinput
 M     ccsm_utils/Tools/generate_batch.csh
 M     ccsm_utils/Tools/testcase_setup.csh
 M     ccsm_utils/Tools/ccsm_splitdf
 M     ccsm_utils/Tools/ccsm_msmkdir
 M     ccsm_utils/Tools/ccsm_l_archive.csh

================================================================================
Originator: mvertens
Date: Fri Sep 26 15:40:25 MDT 2008
Model: scripts
Version: scripts4_080928
One-line summary: first step in unification of ccsm/cam/clm builds

A new Makefile and machine specific Macros were defined that will permit
cam and clm to use the CCSM build directly. 

1) All generic cpp directives are now set in Machines/Macros.cppdefs
2) All compnent template scripts were modified to pass the following argument to gmake
   gmake complib -j \$GMAKE_J MODEL=cam \
         COMPLIB=\$LIBROOT/libxxx.a \
         SMP=\$SMP 
         MACFILE=\$CASEROOT/Macros.\$MACH 
         USER_CPPDEFS="yyy" 
         -f \$BLDROOT/Makefile || exit 2 
   note that xxx above can be atm,lnd.ocn,ice,cpl
   note that yyy above corresponds to the particular cppdefs needed by that component
3) SMP is no longer hard-wired to TRUE in cam and ccsm templates
4) SLIBS, ULIBS, CLIBS are now defined in the Makefile and Macros - not in ccsm_build.csh
5) all libraries are now buit in $LIBROOT
6) Compile Option changes
   Several machine specific changes were made that will result in non-bfb 
   compare with ccsm4_0_alpha36 results
   - on bluefire, the cam optimization was set back to -O3 (also .f90 files are
     no longer created on bluefire for now - this will be put back in)
   - all -r8 options were removed from Linux systems
7) Testing
   Bluefire pretag and posttag tests were run in the sandbox - all tests pass except the following
    ------------pretag
    BFAIL ERS.f19_g15.X.bluefire.compare.ccsm4_0_alpha36 
    FAIL  ERS_D.f45_g35.B.bluefire.compare.ccsm4_0_alpha36 
    BFAIL ERS.f45_g35.B.bluefire.compare.ccsm4_0_alpha36 
    FAIL  ERB.f45_g35.B.bluefire.compare.ccsm4_0_alpha36 
    FAIL  ERH.f45_g35.B.bluefire.compare.ccsm4_0_alpha36 
    FAIL  ERS_D.T31_g35.B.bluefire.compare.ccsm4_0_alpha36 
    FAIL  ERS_D.f19_f19.F.bluefire.compare.ccsm4_0_alpha36 
    FAIL  ERB.f19_f19.F.bluefire.compare.ccsm4_0_alpha36 
    FAIL  ERH.f19_f19.F.bluefire.compare.ccsm4_0_alpha36 
    BFAIL ERB.f45_g35.D.bluefire.compare.ccsm4_0_alpha36 
    FAIL  ERS.f45_g35.FAIS.bluefire.compare.ccsm4_0_alpha36 
    FAIL  ERB.f45_g35.FAIS.bluefire.compare.ccsm4_0_alpha36 
    FAIL  ERT.f19_g15.B.bluefire 
    ------------posttag
    FAIL  ERT.f19_g15.B.bluefire 
    FAIL  ERB.f19_g15.B.bluefire.compare.ccsm4_0_alpha36 
    FAIL  ERP.f19_g15.B.bluefire 
    FAIL  ERP.f19_g15.B.bluefire.compare.ccsm4_0_alpha36 
    FAIL  ERS.f45T42_g35.B.bluefire.compare.ccsm4_0_alpha36 
    FAIL  ERS.f19_f19.F.bluefire.compare.ccsm4_0_alpha36 
    SFAIL ERS.T31_g35.FBN.bluefire.C.newbld02 
    SFAIL ERB.T31_g35.BCNP.bluefire.C.newbld02 
    FAIL  ERH.f45_g35.BMOZ.bluefire.compare.ccsm4_0_alpha36 
    FAIL  ERS.f19_g15.B18C.bluefire.compare.ccsm4_0_alpha36 
    FAIL  ERS.f19_g15.B20TR.bluefire.compare.ccsm4_0_alpha36 
    FAIL  ERS.f19_g15.F18C.bluefire.compare.ccsm4_0_alpha36 
   Jaguar pretag and posttag tests were run in sandbox - all tests pass except for the following
    ------------pretag
    FAIL  ERS.f19_f19.F.jaguar.compare.ccsm4_0_alpha36 
    FAIL  ERS.f19_g15.B.jaguar.compare.ccsm4_0_alpha36 
    FAIL  ERS.f05_g15.B.jaguar.compare.ccsm4_0_alpha36 
    -------------post
    FAIL  ERT.f19_g15.B.jaguar 
    FAIL  ERB.f19_g15.B.jaguar.compare.ccsm4_0_alpha36 
    FAIL  ERB.T62_g35.CECO.jaguar 
    BFAIL ERB.T62_g35.CECO.jaguar.compare.ccsm4_0_alpha36 
    FAIL  ERB.T31_g35.ICN.jaguar.compare.ccsm4_0_alpha36 
    FAIL  ERH.f19_f19.F.jaguar.compare.ccsm4_0_alpha36 
    FAIL  ERS.f19_g15.B18C.jaguar.compare.ccsm4_0_alpha36 
    FAIL  ERS.f19_g15.F18C.jaguar.compare.ccsm4_0_alpha36 
    FAIL  ERS.f19_g15.G18C.jaguar 
    FAIL  ERS.f19_g15.G18C.jaguar.compare.ccsm4_0_alpha36 

Modified files:	
D      ccsm_utils/Build/Macros.CNL
D      ccsm_utils/Build/Macros.Linux.pgi
D      ccsm_utils/Build/Macros.Linux
M      ccsm_utils/Build/mkSrcfiles
D      ccsm_utils/Build/Macros.AIX
D      ccsm_utils/Build/Macros.Linux.ia64
D      ccsm_utils/Build/Macros.BGL
D      ccsm_utils/Build/Macros.BGP
M      ccsm_utils/Build/Makefile
X      ccsm_utils/Tools/perl5lib
M      ccsm_utils/Tools/ccsm_build.csh
M      ccsm_utils/Tools/ccsm_getenv
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Components/xlnd.template
M      ccsm_utils/Components/xocn.template
M      ccsm_utils/Components/xice.template
M      ccsm_utils/Components/xatm.template
M      ccsm_utils/Components/csm_share.buildlib
M      ccsm_utils/Components/slnd.template
M      ccsm_utils/Components/socn.template
M      ccsm_utils/Components/sice.template
M      ccsm_utils/Components/satm.template
A      ccsm_utils/Machines/Macros.surveyor
A      ccsm_utils/Machines/Macros.bangkok.pgf90
D      ccsm_utils/Machines/Macros.BGL.nyblue
D      ccsm_utils/Machines/Macros.BGP.surveyor
D      ccsm_utils/Machines/pes_setup.bluevista
D      ccsm_utils/Machines/Macros.Linux.lightning
A      ccsm_utils/Machines/Macros.cppdefs
M      ccsm_utils/Machines/env_machopts.lightning
D      ccsm_utils/Machines/Macros.CNL.kraken
A      ccsm_utils/Machines/Macros.schirra
D      ccsm_utils/Machines/Macros.BGL.frost
A      ccsm_utils/Machines/Macros.calgary.pgf90
M      ccsm_utils/Machines/env_machopts.atlas
A      ccsm_utils/Machines/Macros.kraken
A      ccsm_utils/Machines/Macros.bangkok
D      ccsm_utils/Machines/Macros.AIX.bluevista
D      ccsm_utils/Machines/Macros.CNL.jaguar
D      ccsm_utils/Machines/Macros.Linux.pgi.ranger
A      ccsm_utils/Machines/Macros.bluefire
A      ccsm_utils/Machines/Macros.frost
A      ccsm_utils/Machines/Macros.nyblue
D      ccsm_utils/Machines/env_machopts.bluevista
A      ccsm_utils/Machines/Macros.franklin
D      ccsm_utils/Machines/Macros.Linux.ia64.atlas
A      ccsm_utils/Machines/Macros.jaguar
D      ccsm_utils/Machines/Macros.Linux.ifort.columbia
D      ccsm_utils/Machines/Macros.AIX.schirra
D      ccsm_utils/Machines/Macros.Linux.ifort
A      ccsm_utils/Machines/Macros.calgary
A      ccsm_utils/Machines/Macros.lightning
A      ccsm_utils/Machines/Macros.atlas
D      ccsm_utils/Machines/Macros.AIX.bluefire
A      ccsm_utils/Machines/Macros.ranger
D      ccsm_utils/Machines/Macros.CNL.franklin
M      ccsm_utils/Machines/env_machopts.ranger
D      ccsm_utils/Machines/mkbatch.bluevista
A      ccsm_utils/Machines/Macros.columbia
D      ccsm_utils/Testlists/bluevista.pretag
M      ccsm_utils/Testlists/bluefire.pretag
D      ccsm_utils/Testlists/bluevista.posttag
M      create_newcase
	
================================================================================
Originator: tcraig
Date: Sep 23 2008
Model: scripts
Version: scripts4_080923
One-line summary: updates for xt4, minor feature fixes

- update xt4 (kraken, franklin, jaguar) to support mixed
  mpi/openmp modes using aprun command line
- change "die" to "print" in check-input-data tool
- update taskmaker.pl for xt4 requirements, add -aprun option
- add -document option to taskmaker.pl
- change taskmaker.pl implementation to keep idle pes 
  in task list.
- add automatic pes documentation to .run script

M      ccsm_utils/Tools/check-input-data
M      ccsm_utils/Tools/taskmaker.pl
M      ccsm_utils/Tools/generate_batch.csh
M      ccsm_utils/Machines/mkbatch.kraken
M      ccsm_utils/Machines/mkbatch.franklin
M      ccsm_utils/Machines/mkbatch.jaguar
================================================================================
Originator:  mvertens
Date: Sun Sep 14 16:25:51 MDT 2008
Model: scripts
Version: scripts4_080914
One-line summary: ERT bug fix for short term archiving

M ccsm_utils/Testcases/config_tests.xml	
	
================================================================================
Originator: tcraig, mvertens
Date: Thu Sep 11 20:29:03 MDT 2008
Model: scripts
Version: scripts4_080911
One-line summary: updates for prestaging, cleanup, trigrid

- delete input-data-check and get-input-data scripts
- update check-input-data for prestaging
- fix configuration settings for 4x5_T42_gx3v5
- review and cleanup sequencing of timing tool output
- chmod +x some of the script scripts
- fix showproj on jaguar
- move cpl.template and ccsm.template to drv/seq_mct/bld/
  and modify create_newcase to reflect the new locations
- make ccsm_getfile and ccsm_getinput just shells (untill all references
  are deleted - at which point they will be removed) 	
- removed the following variables from config_defintion.xml and mad them
  derived variables in ccsm_getenv
     setenv CASETOOLS   $CASEROOT/Tools       # tools used in case scripts
     setenv RUNDIR      $EXEROOT/ccsm_se      # directory of executable
     setenv OBJROOT     $EXEROOT              # location where code is built
     setenv LIBROOT     $EXEROOT/lib          # location of supplemental libraries 
     setenv INCROOT     $EXEROOT/lib/include  # location of supplemental includes/modfiles
     setenv SCRIPTSROOT $CCSMROOT/scripts                  # root of scripts
     setenv UTILROOT    $CCSMROOT/scripts/ccsm_utils       # root of ccsm script utilities
     setenv BLDROOT     $CCSMROOT/scripts/ccsm_utils/Build # makefiles are here
     setenv CODEROOT    $CCSMROOT/models                   # root of model source code 
     setenv SHAREROOT   $CCSMROOT/models/csm_share         # root of share code 
- removed DIN_LOC_MSROOT from config_definition.xml	
- set DOUT_S to TRUE as default value (short term archiving on by default)
- cleanup format of xml files (made them easier to read)
- cleanup help command on create_suite	
	
D      ccsm_utils/Tools/ccsm_getrestart
D      ccsm_utils/Tools/input-data-check
D      ccsm_utils/Tools/get-input-data
D      ccsm_utils/Components/ccsm.template
D      ccsm_utils/Components/cpl.template
	
M      ccsm_utils/Tools/ccsm_build.csh
M      ccsm_utils/Tools/check-input-data
M      ccsm_utils/Tools/ccsm_getfile
M      ccsm_utils/Tools/ccsm_getenv
M      ccsm_utils/Tools/ccsm_getinput
M      ccsm_utils/Tools/timing/getTiming.pl
M      ccsm_utils/Tools/timing/getTiming.csh
M      ccsm_utils/Case.template/Config.pm
M      ccsm_utils/Case.template/config_definition.xml
M      ccsm_utils/Machines/config_machines.xml
M      ccsm_utils/Testlists/jaguar.posttag
	
M      create_clone
M      create_suite

M      ccsm_utils/Case.template/config_grid.xml

M      ccsm_utils/Machines/pes_setup
M      ccsm_utils/Machines/mkbatch.jaguar
M      ccsm_utils/Machines/pes_setup.bluefire
M      ccsm_utils/Machines/Macros.Linux.pgi.ranger
M      ccsm_utils/Machines/config_machines.xml
	
M      ccsm_utils/Testlists/bluefire.posttag
M      ccsm_utils/Testlists/jaguar.posttag
	
================================================================================

Originator: tcraig
Date: Sun Sep 07 2008
Model: scripts
Version: scripts4_080907
One-line summary: modify drv inputs for drvseq2_0_27 update

update drv namelist drv_threading, samegrid, and budget flags.
reorder timing output to match driver better.

M      scripts/ccsm_utils/Tools/timing/gettiming.pl
M      scripts/ccsm_utils/Case.template/config_definition.xml
M      scripts/ccsm_utils/Components/cpl.template
	
================================================================================
Originator: mvertens,mvr
Date: Fri Sep  5 2008
Model: scripts
Version: scripts4_080905
One-line summary: Update scripts trunk to cpl7refactor branch.

The updated scripts will have the following functionality
  (***please see the readme files in the $caseroot/README directories for more detail***).
  - 5 new xml definition files that will be utilized to create the case and relevant tests
      - ccsm_utils/Case.template/config_definition.xml
      - ccsm_utils/Case.template/ config_compsets.xml
      - ccsm_utils/Case.template/config_grid.xml
      - ccsm_utils/Machines/config_machines.xml
      - ccsm_utils/Testcases/config_tests.xml
  - new create_newcase, create_clone and switch_machine scripts that provide a more consistent
       and robust system
  - *** note that these can be invoked outside of the $scriptsroot directory        
  - *** note that switch_machine should be invoked from the new machine
  - xml files for all variables in config_definition.xml
  - cam and clm will be able to start hooking their template scripts directly into these xml files
  - checking for valid values for all short description where valid_values are specified
         (***note that currently unlimited integers and characters are not checked - just valid values lists)
  - auto documentation of variables in env xml files and in $caseroot/README directory
         (no longer have to manually update a README file)
  - SourceMods for only the Compset that is asked for
  - Component templates in Tools/Templates
  - create_newcase can be invoked from any directory (no longer have the be in the $scriptsroot directory)
  - new configure options
        configure -case
        configure -cleannamelist
        configure -cleanmach
        configure -cleanall
  - cleanup of create_newtest and create_suite arguments (only one set of arguments are now supported)


================================================================================
Originator: kauff
Date: Wed Aug 27 15:12:57 EDT 2008
Model: scripts
Version: scripts4_080827
One-line summary: add threading support for jaguar/kraken/franklin, add intrepid support

- Note: this tag works with ccsm4_0_alpha34
- add support for ALCF BG/P "Intrepid
- franklin is now assumed to be new quad-core XT4.
- add threading support for quad-core XT4 jaguar/kraken/franklin 
  and coordinate/align their env_mach and .run scripts
  Note: pes_setup files have not changed, so there is no threading by default.
M     ccsm_utils/Tools/taskmaker.pl
M     ccsm_utils/Machines/run.cnl.jaguar
M     ccsm_utils/Machines/run.cnl.kraken
M     ccsm_utils/Machines/run.cnl.franklin
M     ccsm_utils/Machines/batch.cnl.kraken
M     ccsm_utils/Machines/batch.cnl.franklin
M     ccsm_utils/Machines/batch.cnl.jaguar
M     ccsm_utils/Machines/env.cnl.franklin
M     ccsm_utils/Machines/env.cnl.kraken
M     ccsm_utils/Machines/env.cnl.jaguar

================================================================================

Originator: kauff
Date: Wed Aug 20 15:24:49 MDT 2008
Model: scripts
Version: scripts4_080820b
One-line summary: add support for kraken

A      ccsm_utils/Machines/batch.cnl.kraken
A      ccsm_utils/Machines/Macros.CNL.kraken
A      ccsm_utils/Machines/pes_setup.kraken
A      ccsm_utils/Machines/env.cnl.kraken
A      ccsm_utils/Machines/l_archive.cnl.kraken
A      ccsm_utils/Machines/run.cnl.kraken

A      ccsm_utils/Testlists/kraken.posttag
A      ccsm_utils/Testlists/kraken.preta
M      ccsm_utils/Testlists/jaguar.pretag (add f05_g15.B test)

================================================================================

Originator: mvr
Date: Wed Aug 20 2008
Model: scripts
Version: scripts4_080820
One-line summary: compile cam (AIX) at O2; fix to generate Depends file with 
each make; fix to st archiver for hybrid runs; fix to always link with qsmp

M      ccsm_utils/Build/Macros.AIX
- now compiling cam with O2 optimization rather than O3
M      ccsm_utils/Build/Makefile
- added fix to generate Depends file with every make
M      ccsm_utils/Tools/st_archive.sh
- fix to avoid initial files of previous runs being archived, as was the 
  case with hybrid runs
M      ccsm_utils/Components/ccsm.template
- now always linking with qsmp, even when threads=1, as workaround to 
  compiler problem in cam with mpi-only and optimization at O2
	
================================================================================
Originator: kauff
Date: Fri Aug  8 18:53:11 MDT 2008
Model: scripts
Version: scripts4_080808
One-line summary: add process-binding for threaded code on bluefire

- add process-binding for threaded code on bluefire
  define ibm thread geometry and use Jim Edward's "ccsm_launch" utility...
  setenv NTHRDS $BIND_THRD_GEOMETRY
  mpirun.lsf /contrib/bin/ccsm_launch ./ccsm.exe >&! ccsm.log.$LID
- simplified env vars in env.ibm.bluefire & env.ibm.bluevista
- new location for showproj on jaguar
- move blueice files to Machines_obsolete

M      ccsm_utils/Tools/taskmaker.pl
M      ccsm_utils/Machines/batch.ibm.bluefire
M      ccsm_utils/Machines/run.ibm.bluefire
M      ccsm_utils/Machines/env.ibm.bluefire
M      ccsm_utils/Machines/env.ibm.bluevista
M      ccsm_utils/Machines/batch.cnl.jaguar
M      ccsm_utils/Machines/l_archive.cnl.jaguar
D      ccsm_utils/Machines/*.blueice
A  +   ccsm_utils/Machines_obsolete/*.blueice

================================================================================

Originator: santos,mvr
Date: Thur Jul 31 2008
Model: scripts
Version: scripts4_080731
One-line summary: added inputdata retrieval scripts to ~Tools

A      ccsm_utils/Tools/input-data-check
A      ccsm_utils/Tools/get-input-data

================================================================================
Originator: erik
Date: Wed Jul 30 10:30:36 MDT 2008
Mode: scripts
Version: scripts4_080730
One-line Summary: Add -input_data_list (-l) option to listfilesin_streams script

-Add an option to listfilesin_streams script to output data files in format needed
 by inputdata file retreiving script.

M      ccsm_utils/Tools/listfilesin_streams

================================================================================

Originator: tcraig,kauff
Date: Tue Jul 22 12:00:51 EDT 2008
Model: scripts
Version: scripts4_080722
One-line Summary: build csm_share threaded, fix sed bug with CONTINUE_RUN

- add capability to build csm_share with THREAD set to TRUE if any component
  has more than 1 thread.  same logic existed in driver.
- fix bug where sed command that changes CONTINUE_RUN to TRUE clobbers env_run file on jaguar

M      ccsm_utils/Components/csm_share.buildlib
M      ccsm_utils/Tools/ccsm_run.csh

================================================================================

Originator: mvr
Date: Mon July 14 2008
Model: scripts
Version: scripts4_080714
One-line Summary: modified the new waccm compset to use cice rather than csim

M      ccsm_utils/Compsets/F_WACCM_1995_SMIN

================================================================================
	
Originator: mvr
Date: Thu July 10 2008
Model: scripts
Version: scripts4_080710
One-line Summary: new compset for waccm; new env var for cam's use_case; machine 
files for schirra, columbia (nasa ames)

D      ccsm_utils/Compsets/B_WACCM_1995
A  +   ccsm_utils/Compsets/B_WACCM_1995_CLIM
A      ccsm_utils/Compsets/F_WACCM_1995_SMIN
M      ccsm_utils/Compsets/F_CAMDOM_TGP
M      ccsm_utils/Compsets/F_CAMDOM_CSIM
M      ccsm_utils/Compsets/F_CAMSOM
M      ccsm_utils/Compsets/F_AICE_SOM
M      ccsm_utils/Compsets/K_PRESENT_DAY
M      ccsm_utils/Compsets/B_CO2C_CN_dp
M      ccsm_utils/Compsets/B_1870_CONTROL
M      ccsm_utils/Compsets/B_RAMP_CO2
M      ccsm_utils/Compsets/F_1870_CONTROL
M      ccsm_utils/Compsets/B_CO2C_CN_pp
M      ccsm_utils/Compsets/B_1870-2000_CONTROL
M      ccsm_utils/Compsets/F_CO2B_CASA
M      ccsm_utils/Compsets/B_1990_TYPE1
M      ccsm_utils/Compsets/F_CAMSOM_CSIM
M      ccsm_utils/Compsets/F_CAMDOM
M      ccsm_utils/Compsets/B_CO2B_CN
M      ccsm_utils/Compsets/F_CO2B_CN
M      ccsm_utils/Compsets/B_HIRES
M      ccsm_utils/Compsets/B_CO2C_CN_dd
M      ccsm_utils/Compsets/B_PRESENT_DAY
A      ccsm_utils/Build/Macros.Linux.ifort
M      ccsm_utils/Case.template/env_conf
A      ccsm_utils/Machines/batch.ibm.schirra
A      ccsm_utils/Machines/env.ibm.schirra
A      ccsm_utils/Machines/l_archive.ibm.schirra
A      ccsm_utils/Machines/run.linux.columbia
A      ccsm_utils/Machines/pes_setup.columbia
A      ccsm_utils/Machines/run.ibm.schirra
A      ccsm_utils/Machines/batch.linux.columbia
A      ccsm_utils/Machines/Macros.AIX.schirra
A      ccsm_utils/Machines/Macros.Linux.ifort.columbia
A      ccsm_utils/Machines/env.linux.columbia
A      ccsm_utils/Machines/l_archive.linux.columbia
A      ccsm_utils/Machines/pes_setup.schirra

note: no testing done with machine files for schirra and columbia
	
================================================================================
Originator: mvr
Date: Thu July 3 2008
Model: scripts
Version: scripts4_080703
One-line Summary: add processor binding to bluefire mpirun command

M      scripts/ccsm_utils/Machines/batch.ibm.bluefire
M      scripts/ccsm_utils/Machines/batch.ibm.bluevista
  increased wall clock time for tests
M      scripts/ccsm_utils/Machines/run.ibm.bluefire
  added processor binding to bluefire calls to mpirun to improve performance
M      scripts/ccsm_utils/Machines/pes_setup.bluefire
  added default setting for f19_f19 resolution
M      scripts/ccsm_utils/Testlists/lightning.posttag
M      scripts/ccsm_utils/Testlists/lightning.pretag
D      scripts/ccsm_utils/Testlists/blueice.posttag
D      scripts/ccsm_utils/Testlists/blueice.pretag
M      scripts/ccsm_utils/Testlists/jaguar.posttag
M      scripts/ccsm_utils/Testlists/phoenix.posttag
M      scripts/ccsm_utils/Testlists/bluefire.posttag
M      scripts/ccsm_utils/Testlists/bluevista.posttag
  cleanup of default test lists
	
================================================================================

Originator: mvr
Date: Wed June 11 2008
Model: scripts
Version: scripts4_080611b
One-line Summary: mods to enable create_newcase to be run from cam; cam,clm lt 
archive script updated for bluefire

M      ccsm_utils/Tools/archiving/lt_archive.sh
M      create_newcase
	
================================================================================

Originator: tcraig
Date: Wed June 11 2008
Model: scripts
Version: scripts4_080611a
One-line Summary: add bluefire + 2 bug fixes

Add bluefire scripts
Fix configure bug with cleanall and env_pes
bug #770 change ccsm_mswrite msrcp -a to msrcp -srcdelete -overwrite sizemismatch
  (untested)

M      ccsm_utils/Tools/ccsm_mswrite
M      ccsm_utils/Case.template/configure
A      ccsm_utils/Machines/batch.ibm.bluefire
A      ccsm_utils/Machines/run.ibm.bluefire
A      ccsm_utils/Machines/env.ibm.bluefire
A      ccsm_utils/Machines/l_archive.ibm.bluefire
A      ccsm_utils/Machines/Macros.AIX.bluefire
A      ccsm_utils/Machines/pes_setup.bluefire
A      ccsm_utils/Testlists/bluefire.posttag
A      ccsm_utils/Testlists/bluefire.pretag
	
================================================================================

Originator: erik
Date: Wed June 11 2008
Model: scripts
Version: scripts4_080611
One-line Summary: Update perl5lib version

Just update perl5lib version to latest.

================================================================================

Originator: mvr
Date: Wed June 4 2008
Model: scripts
Version: scripts4_080604
One-line Summary: changed out-of-the-box pe layout for jaguar 1.9x2.5_gx1 B

M      ccsm_utils/Machines/pes_setup.jaguar
M      ccsm_utils/Machines/pes_setup.franklin
-commented out pe setups that were being overwritten and not used

================================================================================
Originator: erik
Date: Thu May 29 14:10:45 MDT 2008
Model: scripts
Version: scripts4_080529
One-line Summary: Fix CLMNCEP directory name for jaguar, and update perl5lib

Update perl5lib to perl5lib_080522
M      ccsm_utils/Machines/env.cnl.jaguar

================================================================================

Originator: mvr
Date: Tue May 27 2008
Model: scripts
Version: scripts4_080527
One-line Summary: bug fix to lt archive script for cam,clm

M      ccsm_utils/Tools/archiving/lt_archive.sh
arg to msrcp command for specifying project number should be '-proj' not '-p'
	

================================================================================

Originator: mvr
Date: Fri May 23 2008
Model: scripts
Version: scripts4_080523
One-line Summary: atm and ocn grids now explicitly set in grid files rather than 
	parsed from $GRID in env_conf

M      ccsm_utils/Case.template/env_conf
M      ccsm_utils/Grids/T85_T85
M      ccsm_utils/Grids/4x5_4x5
M      ccsm_utils/Grids/1.9x2.5_1.9x2.5
M      ccsm_utils/Grids/0.47x0.63_tx0.1v2
M      ccsm_utils/Grids/0.47x0.63_0.47x0.63
M      ccsm_utils/Grids/0.9x1.25_gx1v5
M      ccsm_utils/Grids/T31_gx3v5
M      ccsm_utils/Grids/4x5_gx3v5
M      ccsm_utils/Grids/1.9x2.5_gx1v3
M      ccsm_utils/Grids/1.9x2.5_gx1v4
M      ccsm_utils/Grids/0.23x0.31_tx0.1v2
M      ccsm_utils/Grids/1.9x2.5_gx1v5
M      ccsm_utils/Grids/0.47x0.63_gx1v5
M      ccsm_utils/Grids/1x1.25_gx1v3
M      ccsm_utils/Grids/0.9x1.25_0.9x1.25
M      ccsm_utils/Grids/1x1.25_gx1v5
M      ccsm_utils/Grids/T31_T31
M      ccsm_utils/Grids/T42_gx1v3
M      ccsm_utils/Grids/T42_gx1v4
M      ccsm_utils/Grids/T42_T42
M      ccsm_utils/Grids/T62_gx1v3
M      ccsm_utils/Grids/T42_gx1v5
M      ccsm_utils/Grids/T62_gx1v4
M      ccsm_utils/Grids/gx1v5_gx1v5
M      ccsm_utils/Grids/T62_gx1v5
M      ccsm_utils/Grids/T42_gx3v5
M      ccsm_utils/Grids/0.23x0.31_gx1v5
M      ccsm_utils/Grids/0.9x1.25_tx0.1v2
M      ccsm_utils/Grids/T62_gx3v5
M      ccsm_utils/Grids/T85_gx1v3
M      ccsm_utils/Grids/T62_gx3v6
M      ccsm_utils/Grids/0.23x0.31_0.23x0.31
M      ccsm_utils/Grids/T85_gx1v4


================================================================================

Originator: erik
Date: Mon May 19 11:14:48 MDT 2008
Model: scripts
Version: scripts4_080519
One-line Summary: add ability to set year range for datm7/CLMNCEP option

Add new variables that apply to CLMNCEP option for datm7.

+setenv DATM_CLMNCEP_YR_ALIGN  "1"      # integer; year align (what simulation year corresponds to starting year below)
+setenv DATM_CLMNCEP_YR_START  "1948"   # integer; starting year to loop data over
+setenv DATM_CLMNCEP_YR_END    "2004"   # integer; ending year to loop data over

setenv DIN_LOC_ROOT_CLMNCEP  

>>>>>>>>>>>>>> Add DATM_CLMNCEP vars and document STOP_DATE usage to give the final date to run simulation to.
M      ccsm_utils/Case.template/env.readme
M      ccsm_utils/Case.template/env_run
M      ccsm_utils/Case.template/env_conf
>>>>>>>>>>>>>> Add DIN_LOC_ROOT_CLMNCEP variable.
M      ccsm_utils/Machines/env.linux.atlas
M      ccsm_utils/Machines/env.linux.ranger
M      ccsm_utils/Machines/env.cnl.franklin
M      ccsm_utils/Machines/env.cnl.jaguar
M      ccsm_utils/Machines/env.ibm.frost
M      ccsm_utils/Machines/env.ibm.surveyor
M      ccsm_utils/Machines/env.ibm.nyblue
M      ccsm_utils/Machines/env.ibm.bluevista
M      ccsm_utils/Machines/env.ibm.blueice
M      ccsm_utils/Machines/env.linux.lightning
>>>>>>>>>>>>>> Add year range to "I" cases, and add new compsets with the 1948-2004 date range.
M      ccsm_utils/Compsets/I_PRESENT_DAY
M      ccsm_utils/Compsets/I_CASA
A      ccsm_utils/Compsets/I_CN_1948-2004
A      ccsm_utils/Compsets/I_1948-2004
M      ccsm_utils/Compsets/I_CN

NOTE: REQUIRES datm7_080519 datm7 VERSION TO WORK!

================================================================================

Originator: tcraig
Date: Fri May 16 04:40:43 MDT 2008
Model: scripts
Version: scripts4_080515b
One-line Summary: add SFAIL feature to create_test

M      create_test

- reorganize a bit of logic in create_test to support
  trapping script/case generation failures better.

===============================================================================
Originator: kauff
Date: Thu May 15 18:44:38 EDT 2008
Model: scripts
Version: scripts4_080515
One-line Summary: add jaguar & surveyor Machines files, archive mods wrt cice rest files

  - add Machines files for surveyor and new jaguar, 
  - add tuned pe counts for jaguar and franklin: fv19 fv05 (kauff)
    M /scripts/trunk/ccsm_utils/Machines/*jaguar*
    M /scripts/trunk/ccsm_utils/Machines/*surveyor*
  - rm Testlist for jaguarcnl (kauff,jacob)
    M /scripts/trunk/ccsm_utils/Testlists/*jaguar*
  - mods to archive scripts including the accounting for additional cice restart files (mvr)
    M /scripts/trunk/ccsm_utils/Tools/archiving/lt_archive.sh
    M /scripts/trunk/ccsm_utils/Tools/archiving/st_archive.sh
    M /scripts/trunk/ccsm_utils/Tools/st_archive.sh

===============================================================================
Originator: mvr
Date: Thurs May 08 2008
Model: scripts
Version: scripts4_080508
One-line Summary: bug fixes to archive scripts

  - fix to avoid tarring history files if done annually
  - fix to no longer validate mss writes by file size only
  - fix to no longer archive clm initial files (avoids potential large tarballs)
M      ccsm_utils/Tools/st_archive.sh
M      ccsm_utils/Tools/ccsm_l_archive.csh


===============================================================================
Originator: tcraig
Date: Wed Apr 30 15:21:01 MDT 2008
Model: scripts
Version: scripts4_080430b
One-line Summary: fix configure bug, update ranger module

  - fix new bug in configure cleanall
  - update module load netcdf for ranger

M      ccsm_utils/Case.template/configure
M      ccsm_utils/Machines/env.linux.ranger

===============================================================================
Originator: mvertens
Date: Wed Apr 30 12:03:17 MDT 2008
Model: scripts
Version: scripts4_080430a
One-line Summary: add compset B_1990_TYPE1
M     	ccsm_utils/Case.template/env_conf
A       ccsm_utils/Compsets/B_1990_TYPE1

The new compset permits CAM to be run in a less expensive configuration with 
just three tracers.
	
===============================================================================
Originator: tcraig
Date: Wed Apr 30 00:03:02 MDT 2008
Model: scripts
Version: scripts4_080430
One-line Summary: add ranger, add new decomp features

  - add ranger
  - upgrade taskmaker.pl
  - change case_docs to CaseDocs
  - bring decomp tools into case directory from cice/pop2 models
  - modify configure for setting blocking in env_pes
  - ability to run hybrid settings for cice in decomp tool

M      ccsm_utils/Tools/ccsm_listcache
M      ccsm_utils/Tools/ccsm_build.csh
M      ccsm_utils/Tools/taskmaker.pl
M      ccsm_utils/Case.template/configure
A      ccsm_utils/Machines/run.linux.ranger
A      ccsm_utils/Machines/env.linux.ranger
A      ccsm_utils/Machines/pes_setup.ranger
A      ccsm_utils/Machines/batch.linux.ranger
A      ccsm_utils/Machines/Macros.Linux.pgi.ranger
M      create_newcase
	
===============================================================================
Originator: tcraig
Date: Tue Apr 15 20:48:15 MDT 2008
Model: scripts
Version: scripts4_080416a
One-line Summary: change blueice test to T62_g15.D to T62_g35.D

   - change blueice pretag test suite D test

M      ccsm_utils/Testlists/blueice.pretag

===============================================================================
Originator: tcraig
Date: Tue Apr 15 20:48:15 MDT 2008
Model: scripts
Version: scripts4_080416
One-line Summary: add tx0.1 runoff and minor fix to case_docs

   - add support for tx0.1 runoff (jd)
   - minor update to case_docs "copy"

M    ccsm_utils/Tools/ccsm_build.csh
M    ccsm_utils/Components/cpl.template

===============================================================================
Originator: tcraig
Date: Mon Apr 14 19:24:43 MDT 2008
Model: scripts
Version: scripts4_080415
One-line Summary: general updates, atlas, decomp tool, dead

   - update decomp tool, remove "pop" part (requires new pop tag)
   - add resolved namelists to case_docs directory
   - update timing tool (best with new driver)
   - update dead templates for new decomps (requires new deads)
   - update atlas scripts
   - add 0.23_0.31_gx1v5 resolution
   - update shortnames for consistency
   - set ATM_CDF64 to TRUE for 0.23x0.31 resolutions

M    ccsm_utils/Tools/ccsm_build.csh
M    ccsm_utils/Tools/ConfigInfo.xml
M    ccsm_utils/Tools/timing/getTiming.pl
M    ccsm_utils/Components/xlnd.template
M    ccsm_utils/Components/xocn.template
M    ccsm_utils/Components/xatm.template
M    ccsm_utils/Components/cpl.template
M    ccsm_utils/Components/xice.template
M    ccsm_utils/Machines/Macros.Linux.ia64.atlas
M    ccsm_utils/Machines/env.linux.atlas
M    ccsm_utils/Grids/0.47x0.63_tx0.1v2
M    ccsm_utils/Grids/0.47x0.63_0.47x0.63
M    ccsm_utils/Grids/0.9x1.25_gx1v5
M    ccsm_utils/Grids/0.23x0.31_tx0.1v2
M    ccsm_utils/Grids/0.47x0.63_gx1v5
M    ccsm_utils/Grids/0.9x1.25_tx0.1v2
A    ccsm_utils/Grids/0.23x0.31_gx1v5
M    ccsm_utils/Grids/0.23x0.31_0.23x0.31
	
===============================================================================
Originator: tcraig
Date: Tue Mar 25 17:24:43 MDT 2008
Model: scripts
Version: scripts4_080325
One-line Summary: minor changes to testlist.port and configure -rmmach

M      ccsm_utils/Case.template/configure
M      ccsm_utils/Testlists/testlist.port

===============================================================================
Originator: tcraig
Date: Mon Mar 24 19:20:07 MDT 2008
Model: scripts
Version: scripts4_080324
One-line Summary: updates for atlas, franklin, blueice

 - change batch.ibm.* from csh to tcsh due to task_geom line length
 - remove "i am here" debugs from franklin scripts
 - shorten PBS lines in batch.cnl.franklin
 - update timing tool for new driver timers, including budget diag, run loop barriers
 - update timing files, add component years/day and compset
 - change C and D compset tests from f19_g15 and f45_g35 to T62_g15 in testlists.
 - change pop 80 pe gx1v* decomp in Config.xml to even decomp
 - change cice 160 pe gx1v* decomp in config.xml to bsize_y 384 for perf
 - fix Macros.Linux Build files, CPPFLAGS vs CPPDEFS
 - add tab to histFldsMod logic in Macros.CNL
 - turn on floating point trapping in Macros.CNL for all compilations
 - updates to atlas scripts, update batch parameters for concurrency, change default pes for low resolutions,
 - rename bluevista/ice.ccsm4.pretag testlists to bluevista/ice.pretag
 - add testlist for atlas and franklin.  update testlists
 - minor mod to ccsm_cpdata to make it more robust
 - add testlist.port summarizing porting tests
 - add CCSM_PECOUNT and CCSM_RUNLEN variables to env_conf, partially supported
 - add some CCSM_PECOUNT and CCSM_RUNLEN parameters to some tests
 - add rmmach option to configure
 - update pes_setup.blueice for some f19_g15 B case settings.
 - add CCSM_RUNLEN feature to batch.ibm.blueice
	
M      ccsm_utils/Build/Macros.CNL
M      ccsm_utils/Build/Macros.Linux.pgi
M      ccsm_utils/Build/Macros.Linux.ia64
M      ccsm_utils/Tools/ccsm_cpdata
M      ccsm_utils/Tools/ccsm_listcache
M      ccsm_utils/Tools/ConfigInfo.xml
M      ccsm_utils/Tools/timing/getTiming.pl
M      ccsm_utils/Tools/testcase_setup.csh
M      ccsm_utils/Case.template/configure
M      ccsm_utils/Case.template/env_conf
M      ccsm_utils/Machines/batch.linux.atlas
M      ccsm_utils/Machines/env.linux.atlas
M      ccsm_utils/Machines/pes_setup.blueice
M      ccsm_utils/Machines/batch.ibm.bluevista
M      ccsm_utils/Machines/batch.ibm.blueice
M      ccsm_utils/Machines/batch.cnl.franklin
M      ccsm_utils/Machines/pes_setup.atlas
M      ccsm_utils/Machines/pes_setup.franklin
M      ccsm_utils/Testcases/ERP
M      ccsm_utils/Testcases/ERS
M      ccsm_utils/Testcases/ERT
M      ccsm_utils/Testcases/SMS
M      ccsm_utils/Testcases/ERB
M      ccsm_utils/Testcases/ERH
M      ccsm_utils/Testlists/bluevista.pretag
A      ccsm_utils/Testlists/atlas.pretag
A      ccsm_utils/Testlists/franklin.pretag
D      ccsm_utils/Testlists/bluevista.ccsm4.pretag
M      ccsm_utils/Testlists/blueice.pretag
A      ccsm_utils/Testlists/testlist.port
D      ccsm_utils/Testlists/blueice.ccsm4.pretag

===============================================================================

Originator: tcraig
Date: Tue Mar 11 23:38:38 EDT 2008
Model: scripts
Version: scripts4_080311
One-line Summary: remove env variable for franklin

M      ccsm_utils/Machines/env.cnl.franklin
	+#setenv MPICH_PTL_MATCH_OFF 1
	
===============================================================================
Originator: tcraig
Date: Tue Mar 11 03:56:34 MDT 2008
Model: scripts
Version: scripts4_080309
One-line Summary: CCSM4 first version!  Merge seqmct48_script_080108 tag to trunk

A shocking number of files are different!  This is for CCSM4.
scripts_080229 has been branched in the repo to scripts/branches/ccsm39br
and scripts/branch_tags/ccsm39br_tags.  ccsm39br_tags/ccsm39br00_scripts_080229
is identical to scripts_080229
	
===============================================================================
===============================================================================
	CCSM4 scripts start from here ^
===============================================================================
===============================================================================

Originator: erik
Date: Fri Feb 29 17:21:59 MST 2008
Model: scripts
Version: scripts_080229
One-line Summary: Fix st_archive names for datm, dice, and dlnd

M      ccsm_utils/Tools/archiving/st_archive.sh --- Correct archive directory names for datm?, dice?, and dlnd? restart files.

===============================================================================

Originator: erik
Date: Mon Feb 25 15:46:44 MST 2008
Model: scripts
Version: scripts_080225
One-line Summary: Short-term archive also archives datm, dice, and dlnd restart files

M      ccsm_utils/Tools/archiving/st_archive.sh --- Archive datm?, dice?, dlnd? restart files.
M      ccsm_utils/Components/clm.template --- Let version number of clm for directories be wildcarded.

===============================================================================

Originator: jwolfe
Date: Thu Jan 31 15:48:32 MST 2008
Model: scripts
Version: scripts_080131a
One-line Summary: add CNL compiler info to mct.buildlib 

- add CNL compiler info to mct.buildlib

M      ccsm_utils/Components/mct.buildlib

===============================================================================

Originator: jwolfe
Date: Thu Jan 31 15:27:50 MST 2008
Model: scripts
Version: scripts_080131
One-line Summary: change machine files for franklin, jaguar, and jaguarcnl from linux to CNL

- have all cnl machines -- franklin, jaguar, and jaguarcnl -- use the new
  Macros.CNL file.  This means deleting the machine files with linux as
  the OS for those machines and replacing them with corresponding CNL-named
  files
- point at new set of modules on jaguar and jaguarcnl to use mpt 3.0, which
  is supposed to fix outstanding mpi issues running under CNL
- point at new location of inputdata on franklin

A      ccsm_utils/Machines/Macros.CNL.franklin
A      ccsm_utils/Machines/Macros.CNL.jaguar
A      ccsm_utils/Machines/Macros.CNL.jaguarcnl
A      ccsm_utils/Machines/batch.cnl.franklin
A      ccsm_utils/Machines/batch.cnl.jaguar
A      ccsm_utils/Machines/batch.cnl.jaguarcnl
A      ccsm_utils/Machines/env.cnl.franklin
A      ccsm_utils/Machines/env.cnl.jaguar
A      ccsm_utils/Machines/env.cnl.jaguarcnl
A      ccsm_utils/Machines/l_archive.cnl.jaguar
A      ccsm_utils/Machines/l_archive.cnl.jaguarcnl
A      ccsm_utils/Machines/run.cnl.franklin
A      ccsm_utils/Machines/run.cnl.jaguar
A      ccsm_utils/Machines/run.cnl.jaguarcnl
D      ccsm_utils/Machines/Macros.Linux.franklin
D      ccsm_utils/Machines/Macros.Linux.jaguar
D      ccsm_utils/Machines/Macros.Linux.jaguarcnl
D      ccsm_utils/Machines/batch.linux.franklin
D      ccsm_utils/Machines/batch.linux.jaguarcnl
D      ccsm_utils/Machines/batch.linux.jaguar
D      ccsm_utils/Machines/env.linux.franklin
D      ccsm_utils/Machines/env.linux.jaguar
D      ccsm_utils/Machines/env.linux.jaguarcnl
D      ccsm_utils/Machines/l_archive.linux.jaguar
D      ccsm_utils/Machines/l_archive.linux.jaguarcnl
D      ccsm_utils/Machines/run.linux.franklin
D      ccsm_utils/Machines/run.linux.jaguar
D      ccsm_utils/Machines/run.linux.jaguarcnl

===============================================================================
Originator: kauff, erik, mvr 
Date: Tue Jan 08 2008
Model: scripts
Version: scripts_080108
One-line Summary: archive scripts now save DGVM files; lt_archive.sh now 
returns to initial directory before launching; removed verbose option 
in d*.template csh scripts

- removed verbose option in d*.template csh scripts

M ccsm_utils/Components/dice.template
M ccsm_utils/Components/dlnd.template
M ccsm_utils/Components/docn.template
M ccsm_utils/Components/datm.template

Allow archiving scripts to also save CLM *.hv.* DGVM history files.

M      ccsm_utils/Tools/ccsm_s_archive.csh ------ save hv files
M      ccsm_utils/Tools/archiving/st_archive.sh - save hv files

bug fix to return to initial dir before lt_archiver is launched
M      ccsm_utils/Tools/archiving/lt_archive.sh

	
===============================================================================

Originator: mvertens
Date: Thu Nov 29 21:24:09 MST 2007
Model: scripts
Version: scripts_071129
One-line Summary: Moved remaining active and data template contents to component
  bld directory	 

D ccsm_utils/Components/dice.template.streams.xml
D ccsm_utils/Components/dlnd.template.streams.xml
D ccsm_utils/Components/docn.template.streams.xml
D ccsm_utils/Components/datm.template.streams.xml
  - moved above xml files to data model bld directory
M ccsm_utils/Components/clm.template
M ccsm_utils/Components/pop2.template
M ccsm_utils/Components/cam.template
M ccsm_utils/Components/dice.template
M ccsm_utils/Components/dlnd.template
M ccsm_utils/Components/docn.template
M ccsm_utils/Components/datm.template
   - moved contents of template to component model
     xxx.cpl6.template in the bld directoyr	 
	
===============================================================================
	
Originator: njn01
Date: Wed Nov 28 15:20 MST 2007
Model: scripts
Version: scripts_071128a
One-line Summary: changed CECO compset

M     ccsm_utils/Compsets/C_ECOSYS  (make this consistent with C_PRESENT_DAY)
   - (bugz #513)	

===============================================================================
Originator: dbailey
Date: Wed Nov 28 10:11:52 MST 2007
Model: scripts
Version: scripts_071128
One-line Summary: Moved cice.template to models/ice/cice/bld

M      ccsm_utils/Components/cice.template (simply points to cice bld)
D      ccsm_utils/Components/cice.template.new
D      ccsm_utils/Components/cice.template.streams.xml
M      ccsm_utils/Components/cam.template (removed mss_irt, rest_pfile, and restart_pfile namelist flags)

===============================================================================

Originator: jwolfe
Date: Thu Nov 15 17:03:37 MST 2007
Model: scripts
Version: scripts_071115
One-line Summary: update templates to include newer resolutions plus machine stuff

- update templates to include recently added higher resolutions
- update machine files for LLNL platform atlas, from Art Mirin

M      ccsm_utils/Components/xlnd.template
M      ccsm_utils/Components/cam_clm_shr.template
M      ccsm_utils/Components/xocn.template
M      ccsm_utils/Components/xatm.template
M      ccsm_utils/Components/cpl.template
M      ccsm_utils/Components/xice.template
M      ccsm_utils/Machines/batch.linux.atlas
M      ccsm_utils/Machines/env.linux.atlas

===============================================================================

Originator: jwolfe
Date: Fri Nov  9 11:36:03 MST 2007
Model: scripts
Version: scripts_071109
One-line Summary: slight updates to match new environment on jaguar under CNL

- match new netcdf module on jaguar
- use absolute path to showproj on jaguar since it is not now in the default path

M      ccsm_utils/Machines/env.linux.jaguar
M      ccsm_utils/Machines/l_archive.linux.jaguar
M      ccsm_utils/Machines/batch.linux.jaguar

===============================================================================
Originator: erik
Date: Wed Nov  7 14:06:16 MST 2007
Model: scripts
Version: scripts_071107b
One-line Summary: changes needed to use latest clm branch tag

M      ccsm_utils/Tools/archiving/st_archive.sh --- Don't abort if <5 rpointer files exist
                                                    only abort if NO rpointer files exist.
M      ccsm_utils/Components/clm.template --------- Remove archiving and rpointer file on namelist
                                                    use build-namelist explicitly.

===============================================================================

Originator: jwolfe, mvr, eaton
Date: Wed Nov  7 12:47:06 MST 2007
Model: scripts
Version: scripts_071107
One-line Summary: bug fix for problem with data models on CNL machines

- modifications to CNL platform Macros files to lower the vectorization level
  for data model compilation (bugz #649)

M      ccsm_utils/Machines/Macros.Linux.jaguarcnl
M      ccsm_utils/Machines/Macros.Linux.franklin
M      ccsm_utils/Machines/Macros.Linux.jaguar

===============================================================================

Originator: jwolfe, mvr, eaton
Date: Tue Nov  6 2007
Model: scripts
Version: scripts_071106
One-line Summary: bug fixes and inclusion of new dependency generator

- modifications to jaguar machine files to add support for CNL
- addition of newer machines franklin and jaguarcnl to subList.txt
- put in cam use cases in compsets and env_conf (bugz #648)
- fix to csm_share bld dependencies (bugz #609)

M      ccsm_utils/Machines/env.linux.jaguar
M      ccsm_utils/Machines/subList.txt
M      xlnd.template
M      xocn.template
M      xatm.template
M      clm.template
M      ccsm_se.template
M      cpl.template
M      csm_share.buildlib
M      pop2.template
M      dlnd.template
M      docn.template
M      datm.template
M      xice.template
M      cice.template.new
M      pop.template
M      cam.template
M      cice.template
M      dice.template

	
===============================================================================

Originator: jwolfe
Date: Fri Nov  2 11:42:42 MDT 2007
Model: scripts
Version: scripts_071102 
One-line Summary: updates for CNL machines

- change default runtime for long-term archiving scripts on jaguar
  and jaguarcnl (bugz #643)
- update jaguarcnl machine files
- update franklin machine files
- update jaguar machine files to be consistent with new CNL operating system

M      ccsm_utils/Machines/env.linux.jaguarcnl
M      ccsm_utils/Machines/Macros.Linux.jaguarcnl
M      ccsm_utils/Machines/env.linux.franklin
M      ccsm_utils/Machines/env.linux.jaguar
M      ccsm_utils/Machines/Macros.Linux.franklin
M      ccsm_utils/Machines/batch.linux.jaguar
M      ccsm_utils/Machines/Macros.Linux.jaguar
M      ccsm_utils/Machines/run.linux.jaguar
M      ccsm_utils/Machines/l_archive.linux.jaguarcnl
M      ccsm_utils/Machines/l_archive.linux.jaguar

===============================================================================

Originator: kauff
Date: Tue Oct 30 11:27:35 MDT 2007
Model: scripts
Version: scripts_071030 
One-line Summary: fix bug wrt ERS_OS.T62_g35.C.* failure (.xml typo)

- fix dice.*.xml typo wrt ssmi RESOLUTION=gx3v5

===============================================================================

Originator: kauff
Date: Fri Oct 26 17:16:48 MDT 2007
Model: scripts
Version: scripts_071026b
One-line Summary: introduce new gx1v5 ssmi climatology data file

- introduce new gx1v5 ssmi climatology data file with corrected time axis
- fixes bug causing ERS_OS.T62_g15.C.* failure

===============================================================================

Originator: kauff,erik
Date: Fri Oct 26 14:20:32 MDT 2007
Model: scripts
Version: scripts_071026
One-line Summary: add new build_streams/stream.txt to cice.template

- add new build_streams/stream.txt functionality to cice.template
- fix bug wrt leaving DIN_LOC_ROOT unresolved in d*.buildnml_prestage.csh 
  this should be resolved at *.build time (not configure time)

M      SVN_EXTERNAL_DIRECTORIES (points to new perl5lib)
.      propset to use new perl5lib changes
M      ccsm_utils/Components/cice.template
M      ccsm_utils/Components/cice.template.streams.xml
M      ccsm_utils/Components/datm.template
M      ccsm_utils/Components/dice.template
M      ccsm_utils/Components/dlnd.template
M      ccsm_utils/Components/docn.template
M      ChangeLog

===============================================================================

Originator: kauff
Date: Fri Oct 19 17:40:36 MDT 2007
Model: scripts
Version: scripts_071019
One-line Summary: stream.txt files are created by templates during configure

- new stream.txt file format requires share3_071009 or later
- stream.txt files are created by templates during configure
- stream.txt files are *not* acquired from inputdata dir
- stream.txt file creation via new build_streams tool
- d*7 models don't prestage anything, data files are read in place from inputdata dir
- cpl6 & data models support new grids
- tx01.v2 grids NOT PRODUCTION READY -- runoff map files are missing,
  but this partial support aids code development only

modified files
M      SVN_EXTERNAL_DIRECTORIES (points to new perl5lib)
X      ccsm_utils/Tools/perl5lib
M      ccsm_utils/Components/cice.template.streams.xml
M      ccsm_utils/Components/datm.template
M      ccsm_utils/Components/datm.template.streams.xml
M      ccsm_utils/Components/dice.template
M      ccsm_utils/Components/dice.template.streams.xml
M      ccsm_utils/Components/dlnd.template
M      ccsm_utils/Components/dlnd.template.streams.xml
M      ccsm_utils/Components/docn.template
M      ccsm_utils/Components/docn.template.streams.xml
M      ChangeLog

New grids, supported by cpl data models only (others to follow)
 T42_gx1v5
  0.9x1.25_gx1v5
 0.47x0.63_gx1v5
 0.23x0.31_tx0.1v2
 0.47x0.63_tx0.1v2
  0.9x1.25_tx0.1v2

===============================================================================

Originator: jwolfe
Date: Fri Oct 12 10:49:36 MDT 2007
Model: scripts
Version: scripts_071012
One-line Summary: Change DOCN_MODE from 'SOM' to 'som' to match docn.template (bugz #641)

M      ccsm_utils/Compsets/D_PRESENT_DAY
M      ccsm_utils/Compsets/F_TICE_SOM

* change DOCN_MODE from 'SOM' to 'som' to match docn.template (bugz #641)

===============================================================================

Originator: jwolfe
Date: Thu Oct 11 14:40:18 MDT 2007
Model: scripts
Version: scripts_071011
One-line Summary: Fix clm template file to work on non-NCAR platforms plus misc bugs

M      ccsm_utils/Components/clm.template
M      ccsm_utils/Machines/batch.linux.atlas
M      ccsm_utils/Machines/env.linux.atlas
M      ccsm_utils/Machines/Macros.Linux.ia64.atlas
A      ccsm_utils/Testlists/jaguar.pretag
M      ccsm_utils/Testlists/jaguar.posttag

* fix clm template to work on non-NCAR machines (bugz #634)
* update llnl atlas machine files
* add jaguar pre-tag testlist

===============================================================================

Originator: erik
Date: Fri Oct  5 13:08:59 MDT 2007
Model: scripts
Version: scripts_071005
One-line Summary: Add tool to build streams text files on the fly

M      ccsm_utils/Components/clm.template 
                Fix problem with 1870_to_2000 -- make consistent with 1870-2000.
		Backout changes to use build-namelist, that require a newer version
                of clm.
M      ccsm_utils/Compsets/B_1870-2000_CONTROL
		Add change from Nancy Norton to turn on CFC's for 1870-2000 simulations.

===============================================================================

Originator: erik
Date: Wed Oct  3 16:16:59 MDT 2007
Model: scripts
Version: scripts_071003
One-line Summary: Add tool to build streams text files on the fly

A      SVN_EXTERNAL_DIRECTORIES ------------- Add perl5lib as an external under ccsm_utils/Tools

A      ccsm_utils/Tools/listfilesin_streams -------- List filenames in a streams text file
A      ccsm_utils/Tools/build_streams -------------- Build a streams text file

Readme file put into header of output streams text files.

A      ccsm_utils/Components/template.streams.xml.readme

Streams templates with information needed to build streams
text files for each model that needs them.

A      ccsm_utils/Components/streams.txt.readme ---------- readme info added to output streams.txt files
M      ccsm_utils/Components/clm.template ---------------- Changes to start using clm build-namelist
         Also has files use $DIN_LOC_ROOT in Buildnml_prestage rather than having it resolved at that point.
         (requires clm version bstrms1_clm3_5_11 or later)
A      ccsm_utils/Components/cice.template.streams.xml
A      ccsm_utils/Components/dice.template.streams.xml
A      ccsm_utils/Components/dlnd.template.streams.xml
A      ccsm_utils/Components/docn.template.streams.xml
A      ccsm_utils/Components/datm.template.streams.xml

To use the build_streams stuff in code -- a version of csm_share that has streams such as csm_share3_070918 
(with the restart bug in bugzilla bug #627 fixed) is required.

===============================================================================

Originator: eaton,mvr
Date: Wed Oct  3
Model: scripts
Version: scripts_071003
One-line Summary: mods for new waccm compset

Changes made:
M      ccsm_utils/Components/cam.template
A      ccsm_utils/Compsets/B_WACCM_1995


===============================================================================

Originator: jwolfe
Date: Mon Oct  1 10:27:24 MDT 2007
Model: scripts
Version: scripts_071001
One-line Summary: minor bug fix to new compset.

Changes made:
M      ccsm_utils/Compsets/B_1870-2000_CONTROL

* added run startdate to compset B_1870-2000_CONTROL

===============================================================================

Originator: jwolfe
Date: Fri Sep 28 12:30:43 MDT 2007
Model: scripts
Version: scripts_070928
One-line Summary: new machines and misc bug fixes.

Changes made:
A      ccsm_utils/Compsets/B_1870-2000_CONTROL
M      ccsm_utils/Testlists/blueice.posttag
M      ccsm_utils/Tools/ccsm_listcache
M      ccsm_utils/Components/clm.template
M      ccsm_utils/Components/cam.template
M      ccsm_utils/Components/pop2.template
M      ccsm_utils/Components/mct.buildlib
A      ccsm_utils/Machines/batch.linux.franklin
A      ccsm_utils/Machines/batch.linux.atlas
A      ccsm_utils/Machines/batch.linux.jaguarcnl
A      ccsm_utils/Machines/env.linux.atlas
A      ccsm_utils/Machines/env.linux.jaguarcnl
A      ccsm_utils/Machines/l_archive.linux.jaguarcnl
M      ccsm_utils/Machines/run.ibm.bluevista
M      ccsm_utils/Machines/pes_setup.blueice
M      ccsm_utils/Machines/Macros.BGL.frost
A      ccsm_utils/Machines/Macros.Linux.jaguarcnl
A      ccsm_utils/Machines/env.linux.franklin
M      ccsm_utils/Machines/env.ibm.frost
A      ccsm_utils/Machines/run.linux.atlas
A      ccsm_utils/Machines/run.linux.jaguarcnl
A      ccsm_utils/Machines/Macros.Linux.franklin
A      ccsm_utils/Machines/Macros.Linux.ia64.atlas
A      ccsm_utils/Machines/run.linux.franklin
M      ccsm_utils/Machines/run.ibm.blueice
A      ccsm_utils/Machines/pes_setup.atlas
A      ccsm_utils/Machines/pes_setup.franklin
A      ccsm_utils/Machines/pes_setup.jaguarcnl
M      ccsm_utils/Case.template/env_conf
M      ccsm_utils/Case.template/SourceMods/src.cice/README

* addition of new machine files for LLNL platform atlas
* addition of new machine files for NERSC platform franklin
* addition of new machine files for ORNL platform jaguarcnl (jaguar using CNL OS)
* update of frost files to move system-specific information to machine files
* improved blueice pe configurations to avoid using partial nodes
* fix in run scripts to avoid long term archiving problems on blueice and
  bluevista (bugz #621)
* fix typo in ccsm_listcache script (bugz #607)
* addition of several new POP2 decompositions
* addition of 1870_2000_control compset and corresponding changes to cam.template
  and clm.template files
* addition of Linux.ia64 logic to mct.buildlib

===============================================================================

Originator: dbailey
Date: Fri Sep 14 10:30:17 MDT 2007
Model: scripts
Version: scripts_070914
One-line Summary: Remove CSIM as an ice model option.

M            6377   ccsm_utils/Tools/ccsm_s_archive.csh
M            6377   ccsm_utils/Tools/ccsm_getrestart
M            6377   ccsm_utils/Case.template/env_setcomp
D            6377   ccsm_utils/Case.template/SourceMods/src.csim
D            6377   ccsm_utils/Case.template/SourceMods/src.csim/README
M            6377   ccsm_utils/Case.template/env.readme
M            6377   ccsm_utils/Case.template/env_conf
D            6377   ccsm_utils/Components/csim.template
M            6377   ccsm_utils/Machines/run.es.moon
M            6377   ccsm_utils/Machines/Macros.Linux.jaguar
M            6377   create_newcase

===============================================================================

Originator: jwolfe
Date: Fri Sep  7 09:33:07 MDT 2007
Model: scripts
Version: scripts_070907
One-line Summary: Small changes in cam.template and clm.template

Changes made:
M      ccsm_utils/Components/cam.template
M      ccsm_utils/Components/clm.template

* small bug fixes to the logic in cam.template for hydrid tests
* minor bug fix in clm.template for finidat in hybrid case only
	
===============================================================================

Originator: jwolfe
Date: Thu Sep  6 09:45:47 MDT 2007
Model: scripts
Version: scripts_070906
One-line Summary: Small change in cam.template

Changes made:
M      ccsm_utils/Components/cam.template

* additional logic in cam.template about when to ignore initial condition year
  and date, which fixes a problem with the hybrid tests
	
===============================================================================

Originator: jwolfe
Date: Fri Aug 31 10:56:07 MDT 2007
Model: scripts
Version: scripts_070831
One-line Summary: Support to make csm_share a library instead of source code
                  in each component.

Changes made:
M      ccsm_utils/Tools/generate_resolved.csh
M      ccsm_utils/Tools/ccsm_build.csh
?      ccsm_utils/Components/csm_share.buildlib
M      ccsm_utils/Components/xlnd.template
M      ccsm_utils/Components/xocn.template
M      ccsm_utils/Components/xatm.template
M      ccsm_utils/Components/clm.template
M      ccsm_utils/Components/cpl.template
M      ccsm_utils/Components/pop2.template
M      ccsm_utils/Components/dlnd.template
M      ccsm_utils/Components/docn.template
M      ccsm_utils/Components/datm.template
M      ccsm_utils/Components/xice.template
M      ccsm_utils/Components/csim.template
M      ccsm_utils/Components/cice.template.new
M      ccsm_utils/Components/pop.template
M      ccsm_utils/Components/cice.template
M      ccsm_utils/Components/dice.template

* remove share code from Filepaths in component templates
* modify generate_resolved.csh and ccsm_build.csh to support building
  csm_share library
* addition of csm_share.buildlib
	
===============================================================================

Originator: erik
Date: Thu Aug 23 09:50:14 MDT 2007
Model: scripts
Version: scripts_070823
One-line Summary: Fix clm build and create resolved namelist

Changes made:
A      ccsm_utils/Components/clm.template

Fix bug 574 -- where land has to rebuild everything each time it does a gmake.
Change so that query's are done in the original configure -- so you end up with
a resolved namelist that you can then freely hand edit. Also fix a minor problem
in the transient scenario. Some of this came from Mariana (mvertens).
	
===============================================================================

Originator: mvertens
Date: Thu Aug 16 13:53:24 MDT 2007
Model: scripts
Version: scripts_070816
One-line Summary: added cice.template.new

Changes made:
A      ccsm_utils/Components/cice.template.new

This will need to be changed to cice.template when cice4_0_20070816 and later
is brought into the ccsm tags	
	
===============================================================================
	
Originator: mvr,eaton
Date: Mon Aug  13 2007
Model: scripts
Version: scripts_070813
One-line Summary:  Modify the call to cam's build-namelist to use default 
	use_case when waccm_mozart or trop_mozart chemistry are active.

Changes made:
M      ccsm_utils/Components/cam.template

===============================================================================
	
Originator: jwolfe
Date: Tue Aug  7 10:06:29 MDT 2007
Model: scripts
Version: scripts_070807b
One-line Summary: remove duplicate test

Changes made:
M      ccsm_utils/Testlists/blueice.pretag

Remove duplicate test

===============================================================================

Originator: jwolfe
Date: Tue Aug  7 10:02:30 MDT 2007
Model: scripts
Version: scripts_070807
One-line Summary: miscellaneous fixes

Changes made:
M      F_CO2B_CASA
M      F_PRESENT_DAY
M      F_AICE_SOM
M      B_CO2C_CN_dp
M      B_CO2C_CN_pp
M      B_CO2B_CN
M      F_TICE_SOM
M      F_CO2B_CN
M      B_CO2C_CN_dd
For all of the above, the following was added:
       -CAM_CHEM trop_mozart_ghg_paero

M      ccsm_utils/Tools/ccsm_l_archive.csh
Minor bug fix to long-term archive script to correctly handle history
files staged for restarts

M      ccsm_utils/Machines/pes_setup.blueice
M      ccsm_utils/Machines/pes_setup.bluevista
	
Change comporder on IBMs to make CAM first, in order to avoid the
possibility of having nodes partly threaded and partly not (Bugz #561)

===============================================================================

Originator: mvertens
Date: Mon Jul 30 12:32:47 MDT 2007
Model: scripts
Version: scripts_070730
One-line Summary: added changes to get 1990_control use case to run correctly

Changes made:
M      ccsm_utils/Components/cam.template
M      ccsm_utils/Compsets/B_PRESENT_DAY
	
	 
Put in fix into cam.template and B compset to have -CAM_CHEM trop_mozart_ghg_paero
turned on in the B present day run.	

===============================================================================
Originator: mvertens
Date: Tue Jul 24 10:32:26 MDT 2007
Model: scripts
Version: scripts_070724
One-line Summary: added new 1990_control use case option to cam.template

Changes made:

M      ccsm_utils/Case.template/env_conf
M      ccsm_utils/Components/cam.template

Added the 1990_control use case option available in cam3_5_05 to cam.template
At this time, it is assumed that an CCSM_IPCC setting of "OFF" implies a 
1990_control run. This is now also documented in env_conf.

===============================================================================
Originator: erik
Date: Mon Jul 23 14:49:12 MDT 2007
Model: scripts
Version: scripts_070723
One-line Summary: Turn CLM_ARB_IC to on by default

Changes made:

M      ccsm_utils/Case.template/env_conf


Allow CLM finidat files to be set to null, so that if a resolution doesn't have
finidat files -- it will run with arbitrary IC.

	
===============================================================================

Originator: mvr
Date: Thu Jul 20 
Model: scripts
Version: scripts_070720
One-line Summary: minor bug fixes to long-term archive script so that partial 
	year tar files are not sent to mass storage

Changes made:

M      ccsm_utils/Tools/ccsm_l_archive.csh

	
===============================================================================
Originator: erik
Date: Thu Jul 19 10:08:34 MDT 2007
Model: scripts
Version: scripts_070719
One-line Summary: Demand that input clm files exist and add option required to change
                  to use non-spunup CLM data

Bugs fixed: 549 (when used with the right clm version)

Changes made:

M      ccsm_utils/Case.template/env_conf -- Add CLM_ARB_IC option
M      ccsm_utils/Components/clm.template -- Add demand option to
           file querys so that something must be returned or prestage fails
           if files are not found.

For using with CLM requires:

	pftintdat07_clm3_5_05

or later...
	
===============================================================================

Originator: erik
Date: Fri Jul  6 13:15:43 MDT 2007
Model: scripts
Version: scripts_070706
One-line Summary: Use new clm configure to prepare build and query for datasets

Changes made:

M      ccsm_utils/Tools/ccsm_cpdata 
        - Correct typo
M      ccsm_utils/Case.template/env_conf 
        - Add CLM_DYNNDEP and CLM_DYNPFT
M      ccsm_utils/Components/clm.template 
        - Use clm configure to prepare build, use query to get datasets

For using with CLM requires:

	bld_070702
	pftintdat06_clm3_5_05

or later...
	
===============================================================================

Originator: jwolfe
Date: Mon Jul  2 10:52:38 MDT 2007
Model: scripts
Version: scripts_070702
One-line Summary: mods to support new intercomponent field transfer specification

Changes made:

M      ccsm_utils/Components/dice.template
       - fixed a typo where the env var $DICE_MODE needs a slash in front of it
M      ccsm_utils/Case.template/configure
       - removed coding that stages a specific cpl_fields_mod file based on the
         value of CCSM_BGC
	
===============================================================================

Originator: jwolfe
Date: Tue Jun 12 11:25:34 MDT 2007
Model: scripts
Version: scripts_070613
One-line Summary: changes to include 1x1.25_gx1v5 

Changes made:

M      ccsm_utils/Components/clm.template
M      ccsm_utils/Components/cpl.template
A      ccsm_utils/Grids/1x1.25_gx1v5
       - added capability to do a 1x1.25_gx1v5 run	 
M      ccsm_utils/Compsets/B_1870_CONTROL
       - modified B_1870_CONTROl to match settings in F_1870_CONTROL
	
===============================================================================
	
Originator: jwolfe
Date: Tue Jun 12 11:25:34 MDT 2007
Model: scripts
Version: scripts_070612b
One-line Summary: changes to jaguar specific files to fix bugs

Changes made:

M      ccsm_utils/Machines/l_archive.linux.jaguar
M      ccsm_utils/Machines/Macros.Linux.jaguar

- removed "mppe=1" from PBS commands because jaguar was complaining about it
- added special compilation control for histFldsMod.F90 in clm (Bugzilla 368)

===============================================================================
Originator: tcraig
Date: Mon Jun 11 19:59:41 MDT 2007
Model: scripts
Version: scripts_070612
One-line Summary: update timing tools for new cam/clm timing output

Changes made:

M    scripts/ccsm_utils/Tools/timing/getTiming.pl
M    scripts/ccsm_utils/Tools/timing/getTiming.csh

- update timing tools for new cam/clm timing output, tested
  on bluevista B case.

===============================================================================
Originator: mvr, dbailey
Date: 070611
Model: scripts
Version: scripts_070611b
One-line Summary: mods to lt archive script for optimizing communications with 
	long-term storage; new ice albedos

Changes made:

M   /scripts/trunk/ccsm_utils/Tools/ccsm_l_archive.csh
	now bundling history files into yearly tar files where appropriate;
	lt archive validation now done merely by listed file size
M   /scripts/trunk/ccsm_utils/Components/cice.template
        new ice albedos based on observation - snow albedos (over ice) now 
	represent slightly older snow

===============================================================================
Originator: mvertens
Date: Mon Jun 11 11:14:19 MDT 2007
Model: scripts
Version: scripts_070611
One-line Summary: rework scripts for 1870 IPCC control
	
Changes made:
The following provides a broad overview followed by a detailed list of files changed:
    - Compsets
    - AMIP case has been removed
    - 1870 F, G and B cases have been added
    - B2 compsets removed (pop is now pop2 in compsets
    - pop is now only present in C_PRESENT_DAY1
    - B_IPCC_A1 also removed
 - Machines                 
    - bluevista16 and bluevista8 have now been moved to obsolete machines
    - pes settings for all machines have CSIM_MODE removed and updated for F case settings
 - Testlists
    - testlists have been changed to account for the above
    - all gx1v4 replaced by gx1v5
    - add ERS [B,F,G]_1870_CONTROL tests to posttag lists on blueice and jaguar (untested)
 - Components
    clm.template 
     - pft physiology file for CN mode
     - added co2_ppmv to namelist
       set ndepdat, finidat and pft-physiology file)
    cam.tempate  
     - invoked build-nameslist for CAM_USECASE setting of 1870_CONTROL from env_conf 
    docn7.template
     - setting of 1870 sst dataset in response to IPCC_MODE setting of 1870_CONTROL 
       setting from env_conf
    cice.template        
     - setting of 1870 sst dataset for prescribed cice in response to 
       CCSM_IPCC_ setting of 1870_CONTROL setting from env_conf
 - env_conf has been changed to reflect new settings
	
M      ccsm_utils/Tools/ccsm_listcache
M      ccsm_utils/Tools/create_production_test_readme
M      ccsm_utils/Machines_obsolete/run.linux.anchorage
M      ccsm_utils/Machines_obsolete/run.cpq.generic_cpq
M      ccsm_utils/Machines_obsolete/env.sx.generic_sx
M      ccsm_utils/Machines_obsolete/run.ibm.blackforest
A  +   ccsm_utils/Machines_obsolete/run.ibm.bluevista16
M      ccsm_utils/Machines_obsolete/run.ibm.cheetah32
M      ccsm_utils/Machines_obsolete/env.ibm.bluesky32
A  +   ccsm_utils/Machines_obsolete/env.ibm.bluevista8
M      ccsm_utils/Machines_obsolete/env.sgi.generic_sgi
A  +   ccsm_utils/Machines_obsolete/l_archive.ibm.bluevista8
A  +   ccsm_utils/Machines_obsolete/Macros.AIX.bluevista8
M      ccsm_utils/Machines_obsolete/batch.ibm.cheetah32
A  +   ccsm_utils/Machines_obsolete/Macros.AIX.bluevista16
M      ccsm_utils/Machines_obsolete/run.sx.generic_sx
M      ccsm_utils/Machines_obsolete/run.ibm.bluesky32
M      ccsm_utils/Machines_obsolete/batch.sx.generic_sx
A  +   ccsm_utils/Machines_obsolete/run.ibm.bluevista8
M      ccsm_utils/Machines_obsolete/env.ibm.generic_ibm
M      ccsm_utils/Machines_obsolete/run.sgi.generic_sgi
M      ccsm_utils/Machines_obsolete/batch.ibm.bluesky
M      ccsm_utils/Machines_obsolete/batch.linux64.ram
M      ccsm_utils/Machines_obsolete/batch.ibm.bluesky32
A  +   ccsm_utils/Machines_obsolete/pes_setup.bluevista16
M      ccsm_utils/Machines_obsolete/env.linux.anchorage
M      ccsm_utils/Machines_obsolete/l_archive.sx.rime
A  +   ccsm_utils/Machines_obsolete/batch.ibm.bluevista16
A  +   ccsm_utils/Machines_obsolete/pes_setup.bluevista8
M      ccsm_utils/Machines_obsolete/env.ibm.blackforest
A  +   ccsm_utils/Machines_obsolete/env.ibm.bluevista16
M      ccsm_utils/Machines_obsolete/env.ibm.cheetah32
A  +   ccsm_utils/Machines_obsolete/l_archive.ibm.bluevista16
M      ccsm_utils/Machines_obsolete/batch.cpq.lemieux
M      ccsm_utils/Machines_obsolete/run.ibm.generic_ibm
A  +   ccsm_utils/Machines_obsolete/batch.ibm.bluevista8
	
M      ccsm_utils/Case.template/env.readme
M      ccsm_utils/Case.template/env_conf
	
M      ccsm_utils/Components/clm.template
M      ccsm_utils/Components/cpl.template
M      ccsm_utils/Components/dlnd.template
M      ccsm_utils/Components/docn.template
M      ccsm_utils/Components/csim.template
M      ccsm_utils/Components/cam.template
M      ccsm_utils/Components/cice.template
M      ccsm_utils/Components/dice.template
	
M      ccsm_utils/Machines/pes_setup.bangkok
M      ccsm_utils/Machines/pes_setup.tempest
M      ccsm_utils/Machines/pes_setup.bluevista
M      ccsm_utils/Machines/pes_setup.blueice
D      ccsm_utils/Machines/run.ibm.bluevista16
D      ccsm_utils/Machines/env.ibm.bluevista8
D      ccsm_utils/Machines/l_archive.ibm.bluevista8
D      ccsm_utils/Machines/Macros.AIX.bluevista8
M      ccsm_utils/Machines/env.linux.jaguar
D      ccsm_utils/Machines/Macros.AIX.bluevista16
M      ccsm_utils/Machines/pes_setup.moon
M      ccsm_utils/Machines/pes_setup.jaguar
M      ccsm_utils/Machines/pes_setup.phoenix
M      ccsm_utils/Machines/pes_setup.seaborg
D      ccsm_utils/Machines/run.ibm.bluevista8
M      ccsm_utils/Machines/pes_setup.frost
M      ccsm_utils/Machines/pes_setup.bassi
M      ccsm_utils/Machines/pes_setup
M      ccsm_utils/Machines/run.linux.jaguar
D      ccsm_utils/Machines/pes_setup.bluevista16
M      ccsm_utils/Machines/pes_setup.lightning
D      ccsm_utils/Machines/batch.ibm.bluevista16
D      ccsm_utils/Machines/pes_setup.bluevista8
D      ccsm_utils/Machines/env.ibm.bluevista16
D      ccsm_utils/Machines/l_archive.ibm.bluevista16
D      ccsm_utils/Machines/batch.ibm.bluevista8
	
M      ccsm_utils/Compsets/C_PRESENT_DAY
M      ccsm_utils/Compsets/F_PRESENT_DAY
M      ccsm_utils/Compsets/G_PRESENT_DAY
M      ccsm_utils/Compsets/F_AICE_SOM
M      ccsm_utils/Compsets/B_CO2A
M      ccsm_utils/Compsets/F_CO2A
M      ccsm_utils/Compsets/B_CO2C_CN_dp
M      ccsm_utils/Compsets/B_1870_CONTROL
D      ccsm_utils/Compsets/B_PRESENT_DAY2
A  +   ccsm_utils/Compsets/B_RAMP_CO2
M      ccsm_utils/Compsets/F_1870_CONTROL
M      ccsm_utils/Compsets/B_CO2C_CN_pp
M      ccsm_utils/Compsets/F_CO2B_CASA
M      ccsm_utils/Compsets/B_TROP_MOZART
D      ccsm_utils/Compsets/B_IPCC_20TH_CENTURY
A  +   ccsm_utils/Compsets/C_PRESENT_DAY1
D      ccsm_utils/Compsets/C_PRESENT_DAY2
D      ccsm_utils/Compsets/B_IPCC_RAMPCO2
D      ccsm_utils/Compsets/G_PRESENT_DAY2
M      ccsm_utils/Compsets/G_1870_CONTROL
M      ccsm_utils/Compsets/B_CO2B_CN
M      ccsm_utils/Compsets/F_TICE_SOM
M      ccsm_utils/Compsets/C_ECOSYS
M      ccsm_utils/Compsets/F_CO2B_CN
D      ccsm_utils/Compsets/F_ICE_AMIP
M      ccsm_utils/Compsets/B_CO2C_CN_dd
M      ccsm_utils/Compsets/B_PRESENT_DAY
	
A  +   ccsm_utils/Testlists/bluevista.pretag
M      ccsm_utils/Testlists/lightning.posttag
M      ccsm_utils/Testlists/lightning.pretag
D      ccsm_utils/Testlists/bluevista16.posttag
D      ccsm_utils/Testlists/bluevista16.pretag
M      ccsm_utils/Testlists/blueice.posttag
M      ccsm_utils/Testlists/blueice.pretag
M      ccsm_utils/Testlists/jaguar.posttag
M      ccsm_utils/Testlists/phoenix.posttag
A  +   ccsm_utils/Testlists/bluevista.posttag

===============================================================================
	
Originator: njn01
Date: Wed May 23 15:53 MDT 2007
Model: scripts
Version: scripts_070523
One-line Summary: rework pop2 build namelist procedure; add support for ocn pes
	
Changes made:

M  4835   ccsm_utils/Components/pop2.template

- rework pop2 build namelist procedure (works with ccsm_pop_2_1_20070523 and higher).
  This simplifies pop2.template and creates fully resolved pop2_in in the case dir
- add support for 240, 384, 480, and 768 ocn pes
	

===============================================================================

Originator: tcraig
Date: Mon May 21 23:13:05 MDT 2007
Model: scripts
Version: scripts_070521
One-line Summary: add compsets, fix bugs
	
Changes made:

A      ccsm_utils/Compsets/B_1870_CONTROL
A      ccsm_utils/Compsets/F_1870_CONTROL
A      ccsm_utils/Compsets/G_1870_CONTROL
D      ccsm_utils/Compsets/B_IPCC_1870
D      ccsm_utils/Compsets/B_IPCC_A1
M      ccsm_utils/Testcases/ERB
M      ccsm_utils/Testcases/ERT_script
M      ccsm_utils/Testlists/blueice.posttag
M      ccsm_utils/Testlists/jaguar.posttag
M      ccsm_utils/Tools/ccsm_s_archive.csh
M      ccsm_utils/Tools/testcase_end

- add [B,F,G]_1870_CONTROL compsets (untested, needs review)
- delete B_IPCC_1870, B_IPCC_A1
- add ERS [B,F,G]_1870_CONTROL tests to posttag lists on blueice and jaguar (untested)
- bug #473 - fix ERT test so history files in diff now properly set.  was
	     a file against itself leading to an "easy" pass.
- bug #476 - update pop short term archiving file list for netcdf option
- bug #489 - fix archiving of generated ERB case

===============================================================================
Originator: njn01
Date: Tue May  8 15:54 MDT 2007
Model: scripts
Version: scripts_070508
One-line Summary: add CHL_OPTION support to pop2.template; remove IPCC and RAMP 
   references from pop.template and pop2.template
	
Changes made:

M 4348   ccsm_utils/Components/pop2.template
           add CHL_OPTION support; remove IPCC and RAMP references
M 4348   ccsm_utils/Components/pop.template
           remove IPCC and RAMP references
M 4348   ChangeLog
           add documentation


===============================================================================
Originator: erik
Date: Wed May  2 15:03:04 MDT 2007
Model: scripts
Version: scripts_070502b
One-line Summary: Update surface datasets used by CLM
	
Changes made:

    ccsm_utils/Components/clm.template
 
    Update the surface and ndep datasets used by CLM with new versions
    that fix the LAI and FMAX problems by Nan R. on May 1st.

===============================================================================
Originator: njn01
Date: Wed May  2 11:37 MDT 2007
Model: scripts
Version: scripts_070502
One-line Summary: change default ocean resolution to gx1v5
	
Changes made:

    ccsm_utils/Case.template/env_conf
 
       -setenv GRID           1.9x2.5_gx1v4   # unlimited, see help
       +setenv GRID           1.9x2.5_gx1v5   # unlimited, see help

===============================================================================

Originator: dbailey
Date: Tue May  1 14:05:14 MDT 2007
Model: scripts
Version: scripts_070501
One-line Summary: New default albedos.
	
Changes made:

in ccsm_utils/Components/cice.template change to...

reset default albedos based on observations.

===============================================================================

Originator: erik
Date: Thu Apr 26 14:34:51 MDT 2007
Model: scripts
Version: scripts_070426
One-line Summary: Update pft-physiology file
	
Changes made:

in ccsm_utils/Components/clm.template change to...

set datapft = pft-physiology.c070207

the previous file used had modifications (c070321) needed for CN.


===============================================================================

Originator: tcraig
Date: Wed Apr 25 13:27:38 MDT 2007
Model: scripts
Version: scripts_070425
One-line Summary: Jaguar mods and minor bug fixes
	
Changes made:

- update jaguar testlist, change g14 to g15
- modify restart_compare.pl to produce better error messages
- modify testcase_end to use restart_compare.pl output for test status

===============================================================================

Originator: tcraig
Date: Tue Apr 24 11:49:09 MDT 2007
Model: scripts
Version: scripts_070424
One-line Summary: Update bug fixes and features
	
Changes made:

- modify ERB test to use reference restart.tar files instead of generic
  ccsm_getfile for prestaging.
- bug #435 add module load netcdf on jaguar
- bug #438 update G2 defaults tasks/threads
- bug #438 fix dlnd.template use of $UTILROOT vs \$UTILROOT
- fix some bugs in getTiming.pl tool
- set bluevista == bluevista16 (SMT mode) and add bluevista8
- bug #163 add long term archive test, LAR
- bug #386 modify long term archiving to "trash" restart files, now
  relying on restart.tar for archiving exclusively.
- add ccsm_msls script to ccsm_utils/Tools.  needed for LAR test
- modify test cleanup to add short term archive area
- add more info on task/thread "levels", use m (machine), r (resolution),
  c (compset)
- set default ocean component from pop to pop2
- change default CCSM_COMPSET variable to B2

===============================================================================
Originator: kauff
Date: Fri Apr 20 11:16:26 MDT 2007
Model: scripts
Version: scripts_070420
One-line Summary: support for gx1v5_gx1v5 & parallel atm7, docn7, dice7
	
Changes made:

- add support for gx1v5_gx1v5 grid combination (eg. datm7 + pop)
- add support for parallel datm7, docn7, dice7
  svn status:
  M      ccsm_utils/Components/cpl.template
  M      ccsm_utils/Components/dlnd.template
  M      ccsm_utils/Components/dice.template
  M      ccsm_utils/Components/docn.template
  M      ccsm_utils/Components/datm.template
  A      ccsm_utils/Grids/gx1v5_gx1v5

===============================================================================

Originator: dbailey
Date: Fri Apr 13 10:37:18 MDT 2007
Model: scripts
Version: scripts_070413
One-line Summary: Prestage prescribed ice netcdf file.
	
Changes made:

- fixes bug #436

cice.template:
@@ -39,10 +39,13 @@
 # set ice datasets
 #-----------------------------------------------------------------------
 set pice_stream = " "
+set pice_stream_nc = " "
 if ($CICE_MODE == 'prescribed_ice_clim') then
-   set pice_stream = sst_clim_hurrell_1x1.070323.stream.txt
+   set pice_stream    = sst_clim_hurrell_1x1.070323.stream.txt
+   set pice_stream_nc = sst_clim_hurrell_1x1_050606.nc
 else if ($CICE_MODE == 'prescribed_ice_amip') then
-   set pice_stream = AMIP_bc1x1_1976-1996.070323.stream.txt
+   set pice_stream    = AMIP_bc1x1_1976-1996.070323.stream.txt
+   set pice_stream_nc = AMIP_bc1x1_1976-1996_050706.nc
 endif

 # set local variables needed to prestage data and create resolved namelist
@@ -99,6 +102,7 @@
 set no_ice_ic = 'default'
 if (\$runtype == startup)  set restart = .false.
 set pice_stream = $pice_stream
+set pice_stream_nc = $pice_stream_nc
 set pice_stream_root = ice/cice
 EOF1

@@ -156,6 +160,7 @@
 if ($prescribed_ice == .true.) then
 cat >> $CASEROOT/Buildnml_Prestage/cice.buildnml_prestage.csh << EOF1
   \$UTILROOT/Tools/ccsm_getinput \$pice_stream_root/\$pice_stream . || exit 99
+  \$UTILROOT/Tools/ccsm_getinput \$pice_stream_root/\$pice_stream_nc . || exit
99
 EOF1
 endif

===============================================================================

Originator: tcraig
Date: Tue Apr 10 20:55:50 MDT 2007
Model: scripts
Version: scripts_070410
One-line Summary: Update bug fixes
	
Changes made:

- move bluesky scripts to Machines_obsolete
- bug #164, fix timing reporting for non-hourly coupling
- bug #356, cloned cases now more robust.  an env_conf cache file is
  created when a cloned case is setup and some error checks were
  added to configure to make sure configure wasn't invoked again
  after env_conf is changed or after a cloned case is created
- bug #373, remove ERB.f45_g35.BMOZ.lightning from lightning.pretag
- bug #390, change ORNL queue from debug to batch for jaguar and phoenix
- bug #391, add ERB_OS.f19_g14.B2 and ERH_OS.f19_g14.B2 to blueice and
  bluevista pretag lists, remove from posttag list.
- bug #422, modify task geometry generation to not fail at high
  processor counts
- bug #431, new pes_setup files relecting more flexiblity

===============================================================================

Originator: erik
Date: Thu Apr  5 13:42:13 EDT 2007
Model: scripts
Version: scripts_070405
One-line Summary: Update clm pft dataset
	
Changes made:

in clm.template:
-set datapft = pft-physiology.c061129
+set datapft = pft-physiology.c070321

================================================================================
Originator: erik
Date: Tue Apr  3 13:42:13 EDT 2007
Model: scripts
Version: scripts_070403b
One-line Summary: Fix clm initial dataset, change tests adding gx1v5 tests
	
Changes made:

-Use new clm initial dataset for 1.9x2.5 gx1v4 resolution (from year 95 of b31.020ws case).
-Test ocean resolution for clm initial file.
-Add ocean gx1v5 tests from Nancy Norton just to blueice for both post and pretag
 (This makes the lists between blueice and bluevista different)
-Move gx1v3 test to gx1v4 and correct test name as per Jeff Lee

================================================================================
Originator: dbailey
Date: Tue Apr  3 10:38 EDT 2007
Model: scripts
Version: scripts_070403
One-line Summary: CICE template support for inputdata conventions.
	
Changes made:

ccsm_utils/Components/cice.template

================================================================================
Originator: njn01
Date: Fri Mar 30 11:24 MDT 2007
Model: scripts
Version: scripts_070330
One-line Summary: Added gx1v5 support; modified G compsets
	
Modifications:
 -- add references to gx1v5 (bug #314)
      ccsm_utils/Components/clm.template
      ccsm_utils/Components/cpl.template
      ccsm_utils/Components/pop2.template
      ccsm_utils/Components/dlnd.template
      ccsm_utils/Components/docn.template
      ccsm_utils/Components/datm.template
      ccsm_utils/Components/dice.template
      ccsm_utils/Machines/pes_setup.blueice
 -- modify G-compset settings (bug #423)
      ccsm_utils/Compsets/G_PRESENT_DAY
      ccsm_utils/Compsets/G_PRESENT_DAY2
 -- point to new surface dataset for CLM with corrected LAI 
    at 1.9x2.5 resolution (EK; bug #419)
      ccsm_utils/Components/clm.template
Additions:
 -- new files: support gx1v5  (bug #314)
      ccsm_utils/Grids/T62_gx1v5
      ccsm_utils/Grids/1.9x2.5_gx1v5
Tests/reviews:
 -- A , 1.9x2.5_gx1v5, 5-day smoke test
 -- B2, 1.9x2.5_gx1v5, 2-month exact restart
 -- C2, T62_gx1v5, 5-day smoke test
 -- G2, T62_gx1v5, 5-day smoke test
 -- G2 case setup inspection by SY & it's ok with DB to 
    modify G compset dfns (3/29/07)
 -- it's ok with EK to include previously untagged clm.template 
    surface dataset changes in this tag (3/30/07)
 -- scripts walkthrough with BK (3/29/07)

================================================================================
Originator: dbailey
Date: Wed Feb 28 09:46:02 MST 2007
Model: scripts
Version: scripts_070228
One-line Summary: Added CICE as the default ice component
	
Changes made:
create_newcase
ccsm_utils/Case.template/SourceMods/src.cice
ccsm_utils/Case.template/env_conf
ccsm_utils/Case.template/env.readme
ccsm_utils/Case.template/env_setcomp
ccsm_utils/Components/cice.template
ccsm_utils/Compsets/* <- csim changed to cice in all active ice compsets
ccsm_utils/Machines/pes_setup.* 
ccsm_utils/Machines/Macros.Linux.jaguar
ccsm_utils/Machines/run.es.moon
ccsm_utils/Tools/ccsm_getrestart
ccsm_utils/Tools/ccsm_listcache
ccsm_utils/Tools/ccsm_s_archive.csh
ccsm_utils/Tools/timing/getTiming.pl
Also some cam namelist changes from Mat:
ccsm_utils/Components/cam.template

================================================================================
Originator: erik
Date: Mon Feb 19 21:53:23 MST 2007
Model: scripts
Version: scripts_070219
One-line Summary: Fix CLM surface dataset needed for T31
	
Changes made:

-Correct CLM surface dataset needed for T31 to Mon Feb 19 21:53:23 MST 2007

================================================================================

Originator: jwolfe
Date: Mon Feb 12 14:46:46 MST 2007
Model: scripts
Version: scripts_070212
One-line Summary: clean-up and bug fixes
	
Changes made:

-added blurb to README about how to run the metadata archive script (mvr)
-added default configuration for 1x1.25_gx1* to pes_setup.blueice (jwolfe)
-remove finidat for start-up in clm template (erik)
-removed -DFORTRANUNDERSCORE from CFLAGS in Macors.linux.jaguar because it is now
 included in CPPDEFS in the Macros.Linux file (Bugzilla #360) (jwolfe)

================================================================================

Originator: erik
Date: Fri Feb  9 21:19:16 MST 2007
Model: scripts
Version: scripts_070209
One-line Summary: new datasets needed for the new Hydrology clm3_expa_89 tag
	
Changes made:

-Move changes from hydro070110_scripts_061213 to trunk
-Change clm.template to use new surface datasets
-Change from Tony on jaguar for csim to flush -Mflushz

================================================================================

Originator: kauff
Date: Thu Feb  8 15:04:13 MST 2007
Model: scripts
Version: scripts_070208
One-line Summary: add support for T340x_gx1v3 for dead models only
	
Changes made:

- add support for T340x_gx1v3 for dead models only (bugz # 351)
- intended for scaling studies only

================================================================================

Originator: mvr
Date: 2007-02-06
Model: scripts
Version: scripts_070206
One-line Summary: improved communications between scripts and ccsm run database
	
Changes made:

M      create_newcase
- now validates unique casenames and allows casename reservations via run database
M      ccsm_utils/Tools/archive_metadata.sh
- some bug fixes and enhanced communications for loading cases into run database
M      create_test
M      ccsm_utils/Tools/create_production_test
M      ccsm_utils/Testcases/ERB_script
M      ccsm_utils/Testcases/ERH_script
M      ccsm_utils/Testcases/APT_script
- added arg to create_newcase calls to skip communications with run database
M      ccsm_utils/Tools/ccsm_build.csh
- removed diagnostic output statement

================================================================================

Originator: Jeff
Date: 2007-01-23
Model: scripts
Version: scripts_070123
One-line Summary: update posttag testlists for bluevista16 and blueice

Changes made:
 
ccsm_utils/Testlists/bluevista16.posttag
ccsm_utils/Testlists/blueice.posttag
- bug #324, change from T41_g13 to T31_g35 for test *ICN*
- bug #315, replace ERB and ERH tests with ERS for test *FBN* 

================================================================================

Originator: jwolfe, tcraig
Date:
Model: scripts
Version: scripts_
One-line Summary: bug fixes, minor improvements

Changes made:
- bug #347, add yod -sz to single exec launch
- bug #343, ERT_script bug
- bug #341/348, single exec build dependency improvement
  change dead templates so they don't always build (preproc.h)
  change mct and mph template, so they cp -p to handle build 
    dependcy better
- bug #334, fix mct build on robin
- update phoenix/jaguar pbs settings due to changes at ornl


================================================================================
Originator: tcraig
Date: Thu Jan 11 04:35:36 MST 2007
Model: scripts
Version: scripts_070110
One-line Summary: add blueice support, minor bug fixes

Changes made:
- add blueice support via blueice Machine files
- add blueice test lists, same as bluevista for now
- update archive_metadata.sh from mvr
- add -Mflushz to Macros.Linux.jaguar (bug #342)
- turn on hsi in ccsm_msread.  requires users at ornl have
  password-free access to the hpss. (bug #335)
- migrate from debug to batch queue as default on jaguar
- add comporder capability to pes_setup.mach files and configure

================================================================================
	
Originator: tcraig
Date: Wed Dec 13 23:38:53 EST 2006
Model: scripts
Version: scripts_061213
One-line Summary: minor bug and feature fixes

Changes made:

- add ERP_OS.T62_g35.C2.bluevista16 to bluevista16.pretag suite (bug #306)
- set default cleanup = off in create_test.  create_suite sets 
  this to on by default when running suites. (bug #245)
- fix csh -f on lightning (bug #317)
- set single executable true by default
- fix jaguar task calculation for single_exec mode and 1 pe components (bug #311)
  
================================================================================
	
Originator: jeff
Date: Mon Nov 27 16:05:25 EST 2006
Model: scripts
Version: scripts_061127
One-line Summary: fixed bug in ccsm_utils/Testlists

Changes made:

- ccsm_utils/Testlists/bluevista16.pretag:
  fixed a bug that was causing run_test_suite to fail to submit all tests.
  
================================================================================
Originator: jwolfe
Date: Wed Nov 22 19:05:25 EST 2006
Model: scripts
Version: scripts_061122
One-line Summary: fixed bug in ccsm_se.template

Changes made:

- Components/ccsm_se.template: fixed a bug that was causing uname to be
  interpreted during the instantiation of this script instead of at runtime.
  This was causing all phoenix tests to fail.

================================================================================
Originator: jwolfe
Date: Fri Nov 17 11:56:05 MST 2006
Model: scripts
Version: scripts_061117
One-line Summary: add bluegene support for machine "frost"

Changes made:

- Components/mct.buildlib: added a configure line for BGL OS
- Machines/Macros.BGL.frost: new file
- Machines/batch.ibm.frost: new file
- Machines/env.ibm.frost: new file
- Machines/pes_setup.frost: new file
- Machines/run.ibm.frost: new file

================================================================================
Originator: njn01
Date: Thu Nov 16 14:51 MST 2006
Model: scripts
Version: scripts_061116a
One-line Summary: mods to correct and activate the long-term archiver 
 at ornl, plus restore ocn env vars to env_conf

Changes made:

- env_conf: restore OCN_CHL_TYPE & OCN_CO2_TYPE to env_conf 
  (missing from scripts_061116)
- Tools/ccsm_msmkdir: modify hsi mkdir command to force the creation 
  of the path (mkdir -p)
- Machines/run.cray.phoenix: activate the submission of the long-term archiver
- Machines/env.linux.jaguar & Machines/env.cray.phoenix: 
  setenv DOUT_L_MSROOT csm/$CASE
- Machines/l_archive.linux.jaguar & Machines/l_archive.cray.phoenix:
  modify #PBS commands
  

================================================================================
Originator: mvr
Date: Thu Nov 16 11:10 MST 2006
Model: scripts
Version: scripts_061116
One-line Summary: mods for new run database, mods to fix bug in Filepath for phoenix

Changes made:

- added file archive_metadata.sh to Tools
- added CCSM_COMPSET and CCSM_REPOTAG to env_conf 
  (create_newcase and ccsm_utils/Case.template/env_conf were modified)
- removed models/csm_share/shr and models/csm_share/cpl from Filepath in 
  ccsm_se.template

================================================================================
Originator: njn01
Date: Mon Nov 13 16:25 MST 2006
Model: scripts
Version: scripts_061113
One-line Summary: bugfixes and add new ocean script options

Changes made:

- add OCN_CHL_TYPE, OCN_CO2_TYPE to env_conf (not needed in listcache)
  (#279,#280)
- add OCN_CHL_TYPE, OCN_CO2_TYPE to  B_CO2C_CN_dp,B_CO2C_CN_pp,B_CO2C_CN_dd
  and C_ECOSYS compsets
- correct errors to ocn_tracers in pop.template and pop2.template (resolved vs unresolved)
  (bug #287)
- add ocean-tracer-specification error detection in pop.template and pop2.template


================================================================================

Originator: tcraig
Date: Thu Nov  9 23:55:00 MST 2006
Model: scripts
Version: scripts_061109
One-line Summary: update testlists, add compsets

Changes made:

- add I_CASA, I_CN
- update new Testlists
- minor update to ERB and ERH to restore env_run after run is complete
	
================================================================================
Originator: tcraig
Date: Tue Nov  7 23:55:00 MST 2006
Model: scripts
Version: scripts_061107
One-line Summary: update testlists and fix other things.

Changes made:

- add G2, F_CO2B_CASA
- remove B_ICE_CLIM, B_ICE_AMIP
- add new pretag and posttag test lists
- remove all *.level* test lists 
- merge run_test_suite capability into create_suite script
- update single executable build in scripts (bug #272)
- fix cs.submit and cs.status (bug #270)
- remove .rc files from cam.template, scripts (bug #267)
- fix minor bug in template files (\$)
- add bluevista16, jaguar to subList.txt
- add ccsm_proj logic for account_name to l_archive scripts
- change default in create_suite to clean on
	
================================================================================
Originator: tcraig
Date: Wed Nov  1 23:55:00 MST 2006
Model: scripts
Version: scripts_061101
One-line Summary: add ERT, update create_production_test

Changes made:

- add ERT test, 2 month exact restart diffing cpl history file
- update ERP test, 2 month exact restart diffing cpl log files
- update create_production_test, works with testcase_setup now,
  minor mods required in create_test, testcase_setup
- add hist_compare.csh script
- update ccsm_sedfile to support negative numbers as arguments
- add *_DATE env variables to env_run and updated cpl.template
- added support for optional ~/.ccsm_proj file with user specified
  default project name

================================================================================
Originator: tcraig
Date: Tue Oct 24 23:55:00 MDT 2006
Model: scripts
Version: scripts_061024
One-line Summary: 

Changes made:

- fix single executable on lightning and jaguar by adding
  -Wl,--allow-multiple-definition
- change brnch_retain_casename default to false for cam and clm namelist

================================================================================
Originator: tcraig
Date: Mon Oct 23 16:55:00 MDT 2006
Model: scripts
Version: scripts_061023
One-line Summary: 

Changes made:

- add jaguar as supported machine in scripts/ccsm_utils/Machines
- add l_archive script for jaguar and phoenix
- fix bug in getTiming.pl which caused error in ocn recv time
- add fv 2d decomp automation to cam.template for cam namelist (bug #256)
- update automatic account name and job name setting in bluevista and
  lightning batch and l_archive scripts (bug #266)
- modify phoenix tmp dir name from /scratch/scr101 to /tmp/work
- add -executable [single/multiple] option to create_test command line and 
  _S _M a supported options for test names, just like -debug [on/off] and
  _D _O.  Examples of valid test names are ERS, ERS_D, ERS_S, ERS_SO, ERS_MD
- modify test cleanup strategy, in ccsm_utils/Tools/testcase_end
- update scripts README

================================================================================

Originator: jwolfe
Date: Thu Oct 19 10:55:00 MDT 2006
Model: scripts
Version: scripts_061019
One-line Summary: bug fix in pop2.template

Changes made:

- removed a redundant and incorrect line from the single executable logic in the
  pop2.template file, reported as Bugzilla bug #260.

================================================================================

Originator: mvertens
Date: Sun Oct 15 12:11:00 MDT 2006
Model: scripts
Version: scripts_061015
One-line Summary: addition of cam_clm_share.template

Changes made:
	
- added cam_clm_share.template for dtime and irad settings
  modified cam.template and clm.template to account for new shared template
- put in bug fix to ccsm_utils for new single executable change
	
================================================================================
	
Originator: tcraig
Date: Thu Oct 12 23:22:17 MDT 2006
Model: scripts
Version: scripts_061012
One-line Summary: various upgrades, bluevista16, PFS, account_name, exe tar

Changes made:

- add automatic account_name setting for bluevista and lightning
  in batch.*.$mach
- add automated setting of ptile in batch.ibm.bluevista based
  on the setting of MAX_TASKS_PER_NODE (for bluevista only) (bug #187)
- add bluevista16, update bluevista16 default pes after load balance testing (bug #173)
- add PFS test, 20 day performance/load balance test
- fix exe tar in ccsm_build.csh for multiple executable case (bug #250,#253)
- change hardwired path in ccsm_build.csh on ORNL platform related to
  build of makdep in single executable mode, fix untested
- add SINGLE_EXEC variable to ccsm_listcache (bug #259)
- change default resolution to 1.9x2.5_gx1v4 (bug #229)

================================================================================

Originator: jwolfe
Date: Mon Oct 02 11:36 MDT 2006
Model: scripts
Version: scripts_061002
One-line Summary: several bug fixes to support ccsm3_1_beta38

Changes made:

- /scripts/ccsm_utils/Testcases/ERH_script
   a fix to the POP restart option
- /scripts/ccsm_utils/Components/mct.buildlib
   removes the optimization override on phoenix

================================================================================

Originator: njn01 
Date: Thu Sep 28 17:36 MDT 2006
Model: scripts
Version: scripts_060928
One-line Summary: modify pop2.template to support 80 pes

Changes made: added option for 80 pes and removed a couple of comment lines

================================================================================

Originator: jwolfe 
Date: Mon Sep 18 11:10:07 MDT 2006
Model: scripts
Version: scripts_060921
One-line Summary: some enhancements for the single executable version and some
bug fixes for problems uncovered during pre-tag testing.

Changes made:

- /scripts/ccsm_utils/Components/clm.template
  fixes a bug in setting DTIME and RTM_NSTEPS in the namelist
  adds $CODEROOT/csm_share/eshr to the Filepath
- /scripts/ccsm_utils/Components/csim.template
  fixes a bug in setting the variable KDYN for all CSIM_MODE's
- /scripts/ccsm_utils/Tools/ccsm_build.csh
  enhances the gmake command to build the single executable, including an
  OS-dependent makdep generation (as is done for the component executables)
  and a THREAD flag based on the use of threading in any of the components
- /scripts/ccsm_utils/Testlists/bluevista.level2
  removes four invalid tests
- /scripts/ccsm_utils/Testlists/phoenix.level2
  removes four invalid tests
 
================================================================================

Originator: jeff   
Date: Mon Sep 18 11:10:07 MDT 2006
Model: scripts
Version: scripts_060918
One-line Summary: set the default time step of clm at T31 to 1200 and add
F_AICE_SOM and F_TICE_SOM.

Changes made:

- /scripts/ccsm_utils/Components/clm.template
  set the default time step at T31 to 1200.
- add /sripts/ccsm_utils/Compsets/F_AICE_SOM
- add /sripts/ccsm_utils/Compsets/F_TICE_SOM
 
================================================================================

Originator: jwolfe
Date: Fri Sep 15 11:09:07 MDT 2006
Model: scripts
Version: scripts_060915
One-line Summary: changes to support the single executable version of CCSM

Changes made:

- all component templates, including data and dead models, have logic added to
  build a library instead of an executable if SINGLE_EXEC is TRUE.
- all component templates, including data and dead models, have logic added to
  check for the existence of a library if SINGLE_EXEC is TRUE, much like they
  already do for executables
- the env.ibm.XXX files are changed for all platforms to set mpmd-specific
  environment variables only if SINGE_EXEC is FALSE
- alternative run methods and some run set-up are added, based on the value
  of SINGLE_EXEC, to all platforms in the Machines directory with the exception
  of moon.  However, some of these methods are first guesses and have not been
  tested, but will not affect anythink while SINGLE_EXEC is FALSE (the current
  default value)
- the environment variable SINGLE_EXEC is added to env_conf with a current
  default value of FALSE
- the script ccsm_build.csh has been modified to support the single executable
  version of CCSM, and slightly rearranged.  All single executable support is
  included by conditional blocks based on the value of SINGLE_EXEC, so the changes
  should not affect the multiple executable version
 
 
================================================================================

Originator: tcraig
Date: Thu Aug 31 14:09:07 MDT 2006
Model: scripts
Version: scripts_060831
One-line Summary: set grid checking in cpl6, update clm datasets, other bug fixes

Changes made:

- add eps_algrid_setting logic to cpl.template to set domain
  compare checking for lat/lons to 1e-12 (was 1e-2) if cam and
  clm are the atm/lnd components. (bug #191)
- update clm datasets, updates lats to be consistent with cam.
  done mostly for T* datasets (bug #181).
- change default time limit for bluevista to 1:30 (was 0:50)
- set default for create_test to cleanup=on (was off) (bug #184)
- change cleanup logic so it only cleans up if the test passes (bug #184)
- fix testcase "generate" PASS logic, did not work properly (bug #208)
- remove testing on lightning for multiple of 2 pes, no longer
  a constraint on that machine.
- modify ERS_script and ERP_script to reset env_run variables after
  2nd run is complete.
 
 
================================================================================
Originator: tcraig
Date: Fri Sep  1 14:43:29 MDT 2006
Model: scripts
Version: scripts_060901
One-line Summary: fix logic for testcase generate status

Changes made:

- Fix logic for testcase generate status.
 
 
================================================================================

Originator: njn01  
Date: Mon Aug 21 14:40:21 MDT 2006
Model: scripts
Version: scripts_060821
One-line Summary: update pop2.template for branch and hybrid starts

Changes made:

- make $TAVG_TS_FILE and $TAVG_TS_FILE.hdr optional for branch and hybrid
  starts
- eliminate *.hdr files from the creation of the rpointer.ocn.restart and
  rpointer.ocn.tavg files (this is done only in branch and hybrid starts)
 
 
================================================================================
Originator: jeff  
Date: Aug 17 10:34:59 MDT 2006
Model: scripts
Version: scripts_0608117
One-line Summary: new create_production_test and new testcase ERP

Changes made:

- add create_production_test                    
  to run restart test before production run
- add testcase ERP 
  to run 60-day restart test
 
================================================================================
Originator: kauff
Date: Aug 16 12:54:59 MDT 2006
Model: scripts
Version: scripts_060816
One-line Summary: new domain.lnd.1.9x2.5_gx1v4 w/ "FV pole fix"

Changes made:

- add B2 & C2 (pop2) tests to pre-tag test lists
- datm/dlnd.template ~ new domain.lnd.1.9x2.5_gx1v4.060810.nc has FV "pole fix"

================================================================================

Originator: njn01
Date: Tue Aug 15 MDT 2006
Model: scripts
Version: scripts_060815
One-line Summary: small modifications for pop2 support

Changes made:
- ccsm_utils/Components/pop2.template: modify supported decompositions
- ccsm_utils/Tools/ccsm_s_archive.csh: mv all pop and pop2 *.pop.d?* files
 
==============================================================================

Originator: kauff
Date: Wed Aug  9 17:26:07 MDT 2006
Model: scripts
Version: scripts_060809
One-line Summary: gx1v4 support, new C_PRESENT_DAY functionality

Changes made:

- env_conf, d*7.template, cpl.template ~
  new DATM_NCPL, CPL_EPBAL, CPL_ALBAV to support compset C_PRESENT_DAY
  (see bugzilla #154,155,156,157)
- env_conf   ~ add SSMI as option to setenv DICE_MODE  (#157)
- Grids ~ add 1.9x2.5_gx1v4, remove 2x2.5_gx1v3 (#177)
- datm7/dlnd7.template ~ now use identical list of domain files
- Components/d*.template ~ add SourceMods/src.share to gmake Filepath (#158)
- env.readme ~ fix comment wrt: If SETBLD is set to AUTO... (#186)
- clm.template ~ set datafrac = fracdata_1.9x2.5_gx1v4_060808.nc
- mct.buildlib ~ for linux compile: OPT="-Oscalar3 -Ovector1 -Ostream1
 
==============================================================================

scripts_060803
Originator: Jeff   
Date: Tue Aug 03 10:30:10 MDT 2006
Model: scripts
Version: scripts_060803
One-line Summary: bug fix and update Testlists

Changes made:

- ccsm_utils/Components/cam.template                                        
  bug fix for CAM_CO2_TYPE = 'prognostic'            
 
- ccsm_utils/Tools/restart_compare.pl                                       
  bug fix       
 
- ccsm_utils/Testlists                                                      
  remove invalid tests
 
==============================================================================

scripts_060801
Originator: Erik   
Date: Tue Aug 01 18:30:10 MDT 2006
Model: scripts
Version: scripts_060801
One-line Summary: 

Changes made:

- Add ability to work with new versions of cam and clm
- Minor bug fixes for timing utility (#164)          
- revert COMPORDER so COMP_ATM is last (#165)
- make BSUB -W active for bluevista long term archiver (#162)
- Add new fracdata sets to clm.template (#143)
- Update a couple pop datasets in pop.template from NN.
- For bluevista, change back to ptiles=8 and MAX_TASKS_PER_NODE=8

==============================================================================

scripts_060720
Originator: tcraig
Date: Wed Jul 19 18:30:10 MDT 2006
Model: scripts
Version: scripts_060720
One-line Summary: 
Platforms used for testing: bluesky, bluevista
  - very limited testing in ccsm3_1_beta34

Changes made:

Minor bug fixes for timing utility (#164), revert COMPORDER so
COMP_ATM is last (#165), make BSUB -W active for bluevista
long term archiver (#162).

Add new fracdata sets to clm.template (#143)
Update a couple pop datasets in pop.template from NN.
	
Modified files:

ccsm_utils/Tools/timing/getTiming.pl
ccsm_utils/Case.template/configure
ccsm_utils/Components/clm.template
ccsm_utils/Components/pop.template
ccsm_utils/Machines/l_archive.ibm.bluevista

==============================================================================

scripts_060719
Originator: Erik
Date: Wed Jul 19 11:26:45 MDT 2006
Model: scripts
Version: scripts_060719
One-line Summary: 
Platforms used for testing: bluesky

Changes made:

Change scripts for the ability to run with cam3_3_15 and clm3_expa_65.
These versions change the namelist name and items in both models, and
instead of redirecting unit 5 -- open a file explicitly.

Modified files:
- ccsm_utils/Tools/ccsm_build.csh -- add option to send namelist filename rather
           than redirect stdin from the input namelist.
- ccsm_utils/Components/cam.template
- ccsm_utils/Components/clm.template

So the ???_stdio.nml namelists created for each component (atm, cpl, ice, lnd, ocn)
either set stdin to the input namelist or set the namelist item nlfile to the input
namelist. The $model.buildexe.csh script determines that nlfile will be set if the
output log file has the following string in it...

"Use NLFILE rather than stdin"
 
===============================================================

scripts_060711
Originator: Jeff 
Date: Tue Jul 11 14:10:44 MDT 2006
Model: scripts
Version: scripts_060711
One-line Summary: add run_test_suite and README_production
Platforms used for testing: bluesky, bluevista, lightning, phoenix

Changes made:

add run_test_suite -- associated changes in:
- add ccsm_utils/Testlists                
- remove scripts/regression_test_list               
- remove scripts/regression_test_run
- modify scripts/create_suite
- modify scripts/create_test 
 
add README_production                         
 
===============================================================

scripts_060710
Originator: Tony 
Date: Mon Jul 10 14:10:44 MDT 2006
Model: scripts
Version: scripts_060710
One-line Summary: add documentation of possible test outcomes to cs.status
Platforms used for testing: 

Changes made:

-add documentation of possible test outcomes to cs.status (bug #137)
-hack in a costscale parameter in timing tool, largely to handle      
    ptiles, could be extended further for other issues, add
    PES_PER_NODE env variable to bluevista env.mach file to support
    feature (see bug #136)    
-add source of env_pes file to timing tool (bug #135)
-add removal of old timing and spmd files from EXE area (bug #134)
 
===============================================================

scripts_060706
Originator: njn01 (with some mods provided by tcraig)
Date: Thu Jul  6 MDT 2006
Model: scripts
Version: scripts_060706
One-line Summary: add support for out-of-the-box "C" and "G" compsets, 
         simplify create_newcase, make misc. changes, & minor error corrections. 
          See below for details
Testing: "C" on bluesky,bluevista 
	 "G" on bluesky
          see http://www.cgd.ucar.edu/~njn01/ccsm_pop/Info/20060630_pop2_res_forcing_tests
        
Platforms used for testing: bluesky, bluevista

Changes made:
 - create_newcase: enable user write permission on files acted upon by
                   ccsm_sedfile; replace block of "mkdir -p" and 
                   "cp $fromcase" statements with a loop construct
 - ccsm_utils/Case.template:
      o configure: add "setenv COMPORDER"
      o configure: error correction: insert "echo" in error message
      o env_conf: add OCN_COUPLING and OCN_ICE_FORCING variables
      o env_setcomp: remove "setenv COMPORDER" line
 - ccsm_utils/Components:
      o pop.template: update/change comments, add info to gx1v4 block,
                      correct error introduced scripts_060530, add support
                      for resolution & forcing-independent pop_in.
      o pop2.template: update/change comments, add support for resolution-
                      and forcing-independent pop2_in
      o clm.template: add error checking that prevents encountering an
                      untrapped error in ccsm_splitdf (tcriag)
 - ccsm_utils/Compsets;
      o C_ECOSYS, C_PRESENT_DAY, G_PRESENT_DAY: add OCN_COUPLING and
        OCN_ICE_FORCING settings
 - ccsm_utils/Machines:
      o change ptil=8 --> ptile=16 in batch.ibm.bluevista (tcraig)
      o change MAX_TASKS_PER_NODE = 8 --> MAX_TASKS_PER_NODE = 16 in
        env.ibm.bluevista (tcraig)
 - ccsm_utils/Tools:
      o add OCN_COUPLING and OCN_ICE_FORCING to ccsm_listcache
 - ccsm_utils/Tools/timing:
      o modify getTiming.pl (tcraig)
 
===============================================================

scripts_060613
Originator: kauff
Date: Tue Jun 13 16:10:44 MDT 2006
Model: scripts
Version: scripts_060613
One-line Summary: dlnd7 replaces dlnd6
Testing: A   on bluesky,bluevista
	 C   on bluesky
Platforms used for testing: bluesky, bluevista

Changes made:

dice7 replaces dice6 -- associated changes in:
- ccsm_utils/Components/dice.template
- ccsm_utils/Case.template/env_conf (wrt DICE_MODE)

===============================================================

scripts_060607
Originator: tcraig
Date: Tue Jun  6 22:56:54 MDT 2006
Model: scripts
Version: scripts_060607
One-line Summary: scripts bug fixes
Testing: standard regression tests
Platforms used for testing: lightning, bluevista

Changes made:

- fix problem in F_CO2B_CN compset, ice model not defined
- fix bug, compset and grid are now optional in create_newcase
- fix bug, change order so env_conf command line options have
  precedent over compset, grid, etc.

===============================================================

scripts_060606
Originator: tcraig
Date: Mon Jun  5 22:58:51 MDT 2006
Model: scripts
Version: scripts_060606
One-line Summary: scripts bug fixes
Testing: standard regression tests
Platforms used for testing: bluesky, lightning, bluevista

Changes made:

- add SCRIPTSROOT env variable to env_run, fixed configure -addmach
  problem.
- fix ability to run scripts remotely (from a different directory)
- add some error checking
- increase default bluevista time limit to 50 minutes from 30 due
  to test cost

===============================================================

scripts_060605
Originator: kauff
Date: Mon Jun  5 18:27:03 MDT 2006
Model: scripts
Version: scripts_060605
One-line Summary: dlnd7 replaces dlnd6
Testing: A   on bluesky,bluevista,lightning
	 C,K on bluesky
Platforms used for testing: bluesky, lightning, bluevista

Changes made:

dlnd7 replaces dlnd6 -- associated changes in:
- ccsm_utils/Components/dlnd.template
- ccsm_utils/Components/cpl.template (wrt R19 runoff)
- ccsm_utils/Complset/C_PRESENT_DAY (wrt R19 runoff)
- ccsm_utils/Case.template/env_conf (wrt DLND_MODE & R19)

===============================================================

scripts_060530
Originator: tcraig
Date: Tue May 30 20:40:07 MDT 2006
Model: scripts
Version: scripts_060530
One-line Summary: Merge tstmods branch, major mods to scripts
Testing: all compsets were run at least once
	 all testcases were run at least once
	 tests associate with regression_test_list were run
	 tested in ccsm3_1_beta27
Platforms used for testing: bluesky, bluevista, lightning, phoenix
Changes made:

- add -env_conf and -env_run options to create_newcase
- add ccsm_utils/Compsets directory and compset files
- add ccsm_utils/Grids directory and grid files
- change check_* and get_* tools in ccsm_utils/Tools
- add ccsm_utils/Tools/ccsm_sedfile tool for modifying
  env_conf and env_run
- modify ccsm_utils/Case.template/env_[conf,run] to support
  error checking
- add shortname capability to scripts for grids and compsets
- rewrite testcases
- add debug option to create_test
- change testcase scripts to reflect changes

===============================================================

scripts_060525
Originator: Nancy Norton
Date: 
Model: scripts
Version: scripts_060525
One-line Summary: Add support to scripts for pop2
Testing: Tested in the ccsm3_1_beta28 context: build and run pop2 on bluesky,
         bluesky32; build and run pop on bluesky32; exact restart, 10-day runs,
         pop2, bluesky.
Platforms used for testing: bluesky, bluesky32, bluevista
Changes made:
     -- create_newcase -- add creation of SourceMods/src.pop2
     -- ccsm_utils/Case.template/configure: add pop2 support
     -- ccsm_utils/Case.template/env.readme: add pop2 support
     -- ccsm_utils/Case.template/env_conf: add pop2 support
     -- ccsm_utils/Case.template/env_setcomp: add pop2 support
     -- ccsm_utils/Components/pop2.template: new script
     -- ccsm_utils/Machines: add pop2 support to all pes_setup.$mach scripts
     -- ccsm_utils/Tools/ccsm_s_archive.csh: add pop2 support
     -- ccsm_utils/Tools/check_res: add T62_gx3v6 to supported resolutions

===============================================================

scripts_060513
Originator: jeff   
Date: Sat May 13 13:17:58 MDT 2006
Model: scripts
Version: scripts_060513
One-line Summary: modify to make regression_test suite work on phoenix
Testing: ccsm3_1_beta28 basecode plus scripts changes
	 A,B cases on bluevista, bluesky, lighting
	 A,C cases on phoenix 
Platforms used for testing: bluesky, bluevista, lightning, phoenix
Changes made:

- update files regression_test_view, regression_test_run, 
  regression_test_README and run_BL_test          
	
===============================================================
scripts_060511
Originator: jwolfe, tcraig
Date: Thu May 11 01:17:58 MDT 2006
Model: scripts
Version: scripts_060510, scripts_060511
One-line Summary: phoenix mods, add MPI_TYPE_MAX to all env scripts
Testing: TER.01a.4x5_gx3v5.B.lightning bfb
	 TER.01a.1.9x2.5_gx1v3.B.bluevista bfb
	 TER.01a.T31_gx3v5.B.tempest bld pass
Platforms used for testing: phoenix, bluevista, lightning, tempest
Changes made:
scripts_060510:
- modify cam.template to handle fv transposes (not general)
- modify Macros.phoenix for change to sw_core
scripts_060511:
- add "setenv MPI_TYPE_MAX 100000" to all env.machine files.  this
  is needed for cam under certain configurations.

===============================================================

scripts_060506
Originator: tcraig
Date: Sat May  6 23:17:58 MDT 2006
Model: scripts
Version: scripts_060506
One-line Summary: bug reports and general fixes
Testing: ccsm3_1_beta27 basecode plus scripts changes
	 A,B,X cases on bluevista, bluesky, lighting
	 B bld test on tempest
	 NO testing on phoenix (beta27 doesn't run out of the box)
Platforms used for testing: bluesky, bluevista, lightning, tempest
Changes made:

- remove latm and L, M, and O compsets (bug #33)
- remove special cam build env settings (bug #54)
- added some hooks for long term harvest fails due to ls (bug #85)
  this bug not fully addressed yet.
- fix cam reference date, add namelist ref_ymd, ref_tod (bug #87)
- comment out module list in env.mach, add to env.run files.
  this cleans up script output on machines with modules.
- remove calls to covert lower to upper case for some env
  variables that were TRUE or FALSE.  this was not done 
  consistently and was not needed.  this is just cleanup.
	
===============================================================

scripts_060428
Originator: njn01   
Date: Fri Apr 28 15:49:14 MDT 2006
Model: scripts
Version: scripts_060428
One-line Summary: update scripts/ccsm_utils/Tools/check_res to support T85_gx1v4
Testing: create_newcase ... -res T85_gx1v4 -compset B and configure  
Platforms used for testing: bluesky
Changes made:

scripts/ccsm_utils/Tools/check_res:
- where resok is set, add the line: T85_gx1v4 \
 
===============================================================

scripts_060427
Originator: jeff   
Date: Thu Apr  27 15:41:34 MDT 2006
Model: scripts
Version: scripts_060427
One-line Summary: update cam.template and csim.template
Testing: see tag summary   
Platforms used for testing: bluevista, phoenix, lightning
Changes made:

cam.template:
- cleanup BGC co2 settings
- add consistency checks for BGC co2 settings
                                    
csim.template:
- change kmt data file name for T31_gx3v5
 
===============================================================

scripts_060419
Originator: tcraig
Date: Wed Apr  19 04:41:34 MDT 2006
Model: scripts
Version: scripts_060419
One-line Summary: update timing tool, phoenix Macros
Testing: see tag summary   
Platforms used for testing: bluevista, phoenix, lightning
Changes made:

update timing tool for dead and data models, for new cpl6 timers,
  improved information, fix a few problems. (bug #81)
fix Macros.UNICOS.phoenix (bug #76)
make debug default queue on phoenix
	
===============================================================

scripts_060405
Originator: tcraig
Date: Wed Apr  5 19:41:34 MDT 2006
Model: scripts
Version: scripts_060405
One-line Summary: update - clm datasets, gx1v4, create_suite
Testing: see tag summary   
Platforms used for testing: bluevista, bluesky, lightning
Changes made:

added T42_gx1v4 (tested on bluevista, B and X)
added T62_gx1v4 (untested)
changed clm.template for new grid, frac, and surf datasets
fix bugs #62, #66, #67
add create_suite file (beta version)
mods to create_test to support create_suite
add new FV domain files to dlnd and datm templates (untested)

===============================================================

scripts_060326
Originator: tcraig
Date: Sun Mar 26 23:21:03 MST 2006
Model: scripts
Version: scripts_060326
One-line Summary: update - test scripts, clone, 
Testing: see tag summary   
Platforms used for testing: bluevista, bluesky, lightning
Changes made:

fix -clone option for create_newcase
add/update casebaseid for TestStatus output
modify pes_setup.bluesky32 file to set task/threads for K cases in 
  cam test suite
modify pes_setup file to improve csim pes setting for prescribed mode
change default batch time limit for lightning, bluesky, phoenix
add -R "span" to batch.bluevista to handle "bad" nodes better
minor mods to getTiming.pl
	
===============================================================

scripts_060321
Originator: tcraig
Date: Tue Mar 21 21:01:18 MST 2006
Model: scripts
Version: scripts_060321
One-line Summary: update - test scripts, phoenix, clone, BR.01a
Testing: see tag summary   
Platforms used for testing: bluevista, bluesky, lightning, phoenix
Changes made:

Fix lots of stuff related to testing including ability to run
  test script from another directory, fix generate and compare,
  improve TestStatus output, add ability to generate and compare
  on same test, fix -clean option, update regression_test_view
  for new TestStatus output, add some extra error checking
update phoenix netcdf module and paths
modify mct build to add -r m option, remove selectrealkind option
add #BSUB -N to bluevista and lightning
turn on MP_LABELIO on bluevista
add -clone option to create_newcase
update README.getTiming
modify cam.template and clm.template to fix BR.01a, restart
  file changes
update ccsm_s_archive script
reduce default batch time limit on some platforms
add default account_name for phoenix #BSUB
Fixes bugs #51, #53, #55
	
===============================================================

scripts_060319
Originator: Jeff
Date: Sun Mar 19 15:37:16 MST 2006
Model: scripts
Version: scripts_060319
One-line Summary: update module netcdf on phoenix 
Testing:                                                              
Platforms used for testing: phoenix   
Changes made:

- ccsm_utils/Machines/env.cray.phoenix            
  module load module load netcdf/3.5.1_r4                      
 
- ccsm_utils/Components/cam.template              
  bug fix
 
===============================================================

scripts_060316
Originator: Jeff
Date: Thu Mar 16 11:37:16 MST 2006
Model: scripts
Version: scripts_060316
One-line Summary: add data files from cam to the restart data set
Testing:                                                              
Platforms used for testing: bluevista
Changes made:

- ccsm_utils/Tools/ccsm_s_archive.csh             
  add the .rccsm and .rh* files in cam to the restart data set
 
- ccsm_utils/Components/cam.template              
  bug fix
 
===============================================================

scripts_060217
Originator: kauff
Date: Fri Feb 17 14:37:16 MST 2006
Model: scripts
Version: scripts_0601217
One-line Summary: fixed typo in docn T42 SOM stream.txt file name
Testing:                                                              
Platforms used for testing: bluevista
Changes made:

- fixed typo in docn T42 SOM stream.txt file name

===============================================================

scripts_060216
Originator: kauff
Date: Thu Feb 16 15:11:52 MST 2006
Model: SCRIPTS
Version: scripts_0601216
One-line Summary: lightning bug workarounds, new dx7 inputdata
Testing:                                                              
Platforms used for testing: lightning, bluesky, bluevista
Changes made:

- kludgy workaround in run.linux.lightning wrt cp rpointer
- new dx7 inputdata dirs
- new support for branch runs in dx7 template files
- new datm7 mode: CLMNCEP
- more useful default CASESTR

===============================================================

scripts_060113
Originator: jeff ( Yen-Huei Lee) and kauff
Date: Fri Jan 13 11:29:34 MST 2006
Model: SCRIPTS
Version: scripts_060113
One-line Summary: fixed lightning bug in scripts_060112
Testing:                                                              
Platforms used for testing: lightning and bluesky
Changes made:

- fixed bug in scripts_060112 -- don't link to esmf lib 
  in env.linux.lightning

===============================================================

scripts_060112
Originator: jeff ( Yen-Huei Lee)
Date: Thu Jan 12 14:26:58 MST 2006
Model: SCRIPTS
Version: scripts_060112
One-line Summary: added support for lightning
Testing:                                                              
Platforms used for testing: lightning
Changes made:

- added files in Machines & Tools as necessary to support lightning
  
===============================================================

scripts_b051207
Originator: jeff ( Yen-Huei Lee)
Date: Wed Dec  7 11:06:57 MST 2005
Model: SCRIPTS
Version: SCRIPTSb051207
One-line Summary: remove esmf module
Testing:                                                              
Platforms used for testing: bluesky                                          
Changes made:

- clm.template: add esmf_wrf_timemgr
- remove esmf.buildlib
- ccsm_build.csh: remove esmf
- generate_resolved.csh: remove esmf
- env.*: remove esmf
  
===============================================================
===============================================================

scripts_b051205
Originator: kauff ( Brian Kauffman)
Date: Mon Dec  5 18:14:58 MST 2005
Model: SCRIPTS
Version: SCRIPTSb051205
One-line Summary: docn7 replaces docn6
Testing: A & D cases work at T31/gx3 T42/gx1
Platforms used for testing: bluesky
Changes made:

scripts support docn7 and not docn6

===============================================================
===============================================================

scripts_b051130
Originator: tcraig ( Anthony Craig)
Date: Wed Nov 30 15:50:01 MST 2005
Model: SCRIPTS
Version: SCRIPTSb051130
One-line Summary: add support for bassi
Testing:    B case                                     
Platforms used for testing:   bassi
Changes made:

Add support for bassi, add 5 files, update check_machines.

===============================================================
===============================================================

scripts_b051121
Originator: jeff ( Yen-Huei Lee)
Date: Mon Nov 21 13:30:16 MST 2005
Model: SCRIPTS
Version: SCRIPTSb051121
One-line Summary: bug fix esmf.buildlib for Linux
Testing:                                                              
Platforms used for testing:                                           
Changes made:

- esmf.buildlib 

===============================================================
===============================================================

scripts_b051118
Originator: jeff ( Yen-Huei Lee)
Date: Fri Nov 18 16:13:50 MST 2005
Model: SCRIPTS
Version: SCRIPTSb051118
One-line Summary: add setting external lib to env.*, bug fix cam.template, add CO2C to CCSM_BGC in env_conf
Testing:                                                              
Platforms used for testing:                                           
Changes made:

- all files of /ccsm_utils/Machines/env.*
- /ccsm_utils/Components/cam.template
- enc_conf

===============================================================
===============================================================

scripts_b051114
Originator: kauff ( Brian Kauffman)
Date: Mon Nov 14 13:26:57 MST 2005
Model: SCRIPTS
Version: SCRIPTSb051114
One-line Summary: G-case uses datm7 (not latm6), full set of domain files for dlnd
Testing:                                                              
Platforms used for testing:                                           
Changes made:

o G-case uses datm7 (not latm6) 
o full set of domain files listed in dlnd.template
o TN460 mode replaces TN460NYF mode

===============================================================
===============================================================

scripts_b051102
Originator: jeff ( Yen-Huei Lee)
Date: Wed Nov  2 11:27:27 MST 2005
Model: SCRIPTS
Version: SCRIPTSb051102
One-line Summary: modify README.new, configure and pop.template
Testing:                                                              
Platforms used for testing:                                           
Changes made:

 README.new
 ccsm_util/Components/pop.template
 ccsm_util/Case.template/configure
 
===============================================================
===============================================================

scripts_b051011
Originator: tcraig ( Anthony Craig)
Date: Tue Oct 11 06:46:43 MDT 2005
Model: SCRIPTS
Version: SCRIPTSb051011
One-line Summary: Add bluevista scripts
Testing:                                                              
Platforms used for testing:                                           
Changes made:

Add bluevista scripts (beta version)

===============================================================
===============================================================

scripts_b051006
Originator: tcraig ( Anthony Craig)
Date: Thu Oct  6 20:59:23 MDT 2005
Model: SCRIPTS
Version: SCRIPTSb051006
One-line Summary: Update scripts to add grid initialization in dlnd via domain file
Testing:                                                              
Platforms used for testing:                                           
Changes made:
 - add domain files in dlnd.template for new domain files with frac
 - modify env.sgi.tempest for NETCDF_MOD env variable

===============================================================
===============================================================

scripts_b051005
Originator: kauff ( Brian Kauffman)
Date: Wed Oct  5 13:23:22 MDT 2005
Model: SCRIPTS
Version: SCRIPTSb051005
One-line Summary: major datm.template change: now supports datm7 and not datm6
Testing: works on bluesky
Platforms used for testing:                                           
Changes made:

o datm6 is out, datm7 is in: new datm.template
o one new setenv DATM_MODE in env_conf also needed to support datm7

===============================================================
===============================================================

scripts_b051004b
Originator: jeff ( Yen-Huei Lee)
Date: Tue Oct  4 13:41:52 MDT 2005
Model: SCRIPTS
Version: SCRIPTSb051004b
One-line Summary: update new csim.template
Testing:                                                              
Platforms used for testing:                                           
Changes made:

ccsm_utils/Components/csim.template       

===============================================================
===============================================================

scripts_b051004
Originator: jeff ( Yen-Huei Lee)
Date: Tue Oct  4 10:41:31 MDT 2005
Model: SCRIPTS
Version: SCRIPTSb051004
One-line Summary: merging scripts_b050107_brnchT_bgc24
Testing:                                                              
Platforms used for testing:                                           
Changes made:

start from scripts_b051002.
===============================================================
===============================================================

===============================================================

scripts_b050107
Originator: weiyu ( Wei Yu)
Date: Fri Jan  7 11:08:52 MST 2005
Model: SCRIPTS
Version: SCRIPTSb050107
One-line Summary: scripts for combine CRAY-x1 and SGI-ram
Testing:                                                              
Platforms used for testing:                                           
Changes made:
scripts by combining  CRAY-x1 and SGI-ram
      
===============================================================
===============================================================

scripts_b041210
Originator: weiyu ( Wei Yu)
Date: Fri Dec 10 12:57:55 MST 2004
Model: SCRIPTS
Version: SCRIPTSb041210
One-line Summary: fix a small bug in script.
Testing:                                                              
Platforms used for testing:                                           
Changes made:
Platforms used for testing: No.
Changes made:
1) "configure" file.
2) rename all of the modules files.
===============================================================
===============================================================

scripts_b041025
Originator: weiyu ( Wei Yu)
Date: Mon Oct 25 09:48:26 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb041025
One-line Summary: Fix ice model bug, add and drop some machines.
Testing:                                                              
Platforms used for testing:                                           
Changes made:
1) add machine tempest.
2) drop machines: cheetah, chinook.
3) fix the bug of the scripts for ice model.
      
===============================================================
===============================================================

scripts_b041014
Originator: weiyu ( Wei Yu)
Date: Thu Oct 14 09:09:23 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb041014
One-line Summary: scripts changed for F case
Testing:                                                              
Platforms used for testing:                                           
Changes made:
1) The scripts changed for adding F case.
2) Disable the function of automatic cleaning archive log files.
===============================================================
===============================================================

scripts_b041005
Originator: weiyu ( Wei Yu)
Date: Tue Oct  5 11:24:38 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb041005
One-line Summary: add path for netcdf.mod
Testing:                                                              
Platforms used for testing:                                           
Changes made:
Define the path for netcdf.mod. User can define in "environment" or use default path. Fix the minor bug for running on lemieux with loading a new module file.

Tested machine: lemieux, blackforest, bluesky, seaborg, cheetah, eagle, Jazz.
===============================================================
===============================================================

scripts_b041005
Originator: weiyu ( Wei Yu)
Date: Tue Oct  5 11:18:53 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb041005
One-line Summary: add path for netcdf.mod
Testing:                                                              
Platforms used for testing:                                           
Changes made:

Define the path for netcdf.mod. User can define in "environment" or use default path. Fix minor bug for running on lemieux with loading a new module file.

Tested machine: blackforest, bluesky, eagle, cheetah, jazz, seaborg.
===============================================================
===============================================================

scripts_b041001
Originator: weiyu ( Wei Yu)
Date: Thu Sep 30 11:24:16 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb041001
One-line Summary: add netcdf.mod for ccsm.
Testing:   Build and ER test                                                           
Platforms used for testing:lemieux,seaborg,blackforest,bluesky,eagle,cheetah,jazz                                           
Changes made:

   1)  Define NETCDF_MOD for the path of netcdf.mod. User can define in "environment" or use default path.

   2)  Fix minor bug for running on lemieux with loading new "module" file.  
     

===============================================================
===============================================================
	
scripts_b040610
Originator: mvertens ( Mariana Vertenstein)
Date: Thu Jun 10 16:24:19 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040610
One-line Summary: fixed minor bugs
Testing:                                                              
Platforms used for testing:                                           
Changes made:
1) fixed minor bugs

===============================================================
===============================================================

scripts_b040602
Originator: hender ( Tom Henderson)
Date: Wed Jun  2 16:37:17 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040602
One-line Summary: New baselineroot option for create_test, extended MSS retention for restart.tars
Testing:    THY.02a.T42_gx1v3.B.bluesky.compare.ccsm3_0_beta22  PASS
Tested with and without new baselineroot option.  
Platforms used for testing:  bluesky
Changes made:
  create_test:  new -baselinroot option
  ccsm_l_archive.csh:  Files in restart.tars now have the same MSS retention 
                       period as other files.  Code review by Mariana 
                       Vertenstein.  
  pop.template:  Added support for one or two MPI processes.  Reviewed by 
                 Nancy Norton and tested by George Carr on bangkok.  

===============================================================

===============================================================
scripts_b040526
Originator: hender ( Tom Henderson)
Date: Wed May 26 11:03:46 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040526
One-line Summary: Equivalent to scripts_b040521_brnchT_ccsm3rel02.  Refined regression tests.
Testing:  
  TBR.01a.T42_gx1v3.B.bluesky.compare.ccsm3_0_beta22  PASS
  TBR.02a.T42_gx1v3.B.bluesky.compare.ccsm3_0_beta22  PASS
  TER.01a.T42_gx1v3.B.bluesky.compare.ccsm3_0_beta22  PASS
  TER.01b.T42_gx1v3.B.bluesky.compare.ccsm3_0_beta22  PASS
  THY.01a.T42_gx1v3.B.bluesky.compare.ccsm3_0_beta22  PASS
  THY.02a.T42_gx1v3.B.bluesky.compare.ccsm3_0_beta22  PASS
  TSM.01a.T42_gx1v3.B.bluesky.compare.ccsm3_0_beta22  PASS
Platforms used for testing:  IBM
Changes made:
Refined regression testing. Added -inputdataroot option.  Improved 
log output and state machine output (they need a lot more...).  
===============================================================

===============================================================
scripts_b040521
Originator: mvertens ( Mariana Vertenstein)
Date: Fri May 21 16:13:14 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040521
One-line Summary: pre-release branch tag
Testing: done by Tom Henderson
Platforms used for testing:IBM, SGI
Changes made:
1) New regression testing scripts 
2) Numerous changes for release branch
===============================================================
	
===============================================================
scripts_b040517
Originator: mvertens ( Mariana Vertenstein)
Date: Mon May 17 13:37:03 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040517
One-line Summary: put in error check into configure
Testing:                                                              
Platforms used for testing:                                           
Changes made:
1) Put in significantly more error checks into 
   ccsm_utils/Case.template/configure	
2) Updated README comments in all SourceMods files except 
   from dead components
3) Removed gx3v4 supported resolutions from Tools/check_res
4) Made CASESTR non-resolved dice.template, dlnd.template,
   docn.template,datm.template, cam.template, clm.template, 
   csim.template, latm.template, cpl.template
5) Updated to pop.template provided by Nancy Norton on 040516
6) Removed recursive copy from create_newcase
===============================================================
	
===============================================================

scripts_b040507
Originator: mvertens ( Mariana Vertenstein)
Date: Fri May  7 12:29:54 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040507
One-line Summary: removed references and tests related to 2000_CONTROL
Testing:                                                              
Platforms used for testing:                                           
Changes made:
1) removed references and tests related to 2000_CONTROL - this
   was requested by Bill Collins	  
===============================================================
	
===============================================================

scripts_b040504
Originator: mvertens ( Mariana Vertenstein)
Date: Tue May  4 16:33:35 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040504
One-line Summary: put in consistent values of inputdata/xxx/dxxx6
Testing:                                                              
Platforms used for testing:                                           
Changes made:
1) made sure that all data models, latm and cpl are
   using inputdata/model/dxxx6 not inputdata/model/dxxx5
===============================================================
	
===============================================================
	
scripts_b040503
Originator: mvertens ( Mariana Vertenstein)
Date: Mon May  3 10:04:43 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040503
One-line Summary: went back to original version of ccsm_cpdata
Testing:                                                              
Platforms used for testing:                                           
Changes made:
1) Due to pop restart problems, need to sometimes overwrite
 existing file - so went back to original version of ccsm_cpdata.

===============================================================
===============================================================

scripts_b040430a
Originator: tcraig ( Anthony Craig)
Date: Fri Apr 30 15:00:11 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040430a
One-line Summary: Update scripts
Testing: seaborg, TER.01a.T42_gx1v3.A.seaborg
Platforms used for testing: seaborg                                          
Changes made:

Add comment in buildnml_prestage file for each component, warns users
  about creating unique input data filenames
Add DIN_LOC_ROOT_USER env variable in all env.* machine files
Change default SETBLD to TRUE in env_run
Add test for exe tar file to check for uniqueness versus previous version,
  remove if not unique, save if unique, in ccsm_build
Modify run.ibm* files for poe.cmdfile, add full path to executable
Check for single inputdata filename in DIN_LOC_ROOT and DIN_LOC_ROOT_USER,
  if duplicate name, error and exit
Clean up ccsm_cpdata error messages, add info about where file was found,
  remove info about all search locations
Add DIN_LOC_ROOT_USER search in ccsm_getinput, secondary to DIN_LOC_ROOT

===============================================================
===============================================================

scripts_b040430
Originator: mvertens ( Mariana Vertenstein)
Date: Fri Apr 30 11:29:44 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040430
One-line Summary: fixed csim template bug
Testing:                                                              
Platforms used for testing:                                           
Changes made:
1) fixed csim.template bug introduced in scripts_b040429b

===============================================================
===============================================================

scripts_b040429b
Originator: mvertens ( Mariana Vertenstein)
Date: Thu Apr 29 21:50:15 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040429b
One-line Summary: updated csim.template 
Testing:                                                              
Platforms used for testing:                                           
Changes made:
1) updated csim.template to always set no_ice_ic to false 
   for continue run

===============================================================
===============================================================

scripts_b040429
Originator: mvertens ( Mariana Vertenstein)
Date: Thu Apr 29 13:23:30 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040429
One-line Summary: put help functionality back into create_newcase, 
	create_test and configure
Testing: Verified that functionality works
Platforms used for testing:                                           
Changes made:
1) put put help functionality back into create_newcase, create_test and configure
2) put in fix for bangkok script (Wei Yu)
3) put in new debug tests (Wei Yu)	

===============================================================
===============================================================

scripts_b040428
Originator: mvertens ( Mariana Vertenstein)
Date: Wed Apr 28 09:43:51 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040428
One-line Summary: updated component templates
Testing:                                                              
Platforms used for testing:                                           
Changes made:
1) updated cpl.template, dlnd.template
2) updated esmf.buildlib to add Cray X1 change

===============================================================
===============================================================

scripts_b040426b
Originator: mvertens ( Mariana Vertenstein)
Date: Mon Apr 26 11:27:40 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040426b
One-line Summary: put raw surface datasets back in clm namelist
Testing:                                                              
Platforms used for testing:                                           
Changes made:
 Put raw surface datasets back in clm prestaging and namelist

===============================================================
===============================================================

scripts_b040426
Originator: mvertens ( Mariana Vertenstein)
Date: Mon Apr 26 10:48:37 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040426
One-line Summary: updated cam and clm tempates
Testing:                                                              
Platforms used for testing:                                           
Changes made:
1) removed raw surface datasets from generated clm resolved
   namelist - still in template file
2) fixed bug generated when prestage directive was removed
3) removed camclm_share from cam and clm namelist

===============================================================
===============================================================

scripts_b040423
Originator: mvertens ( Mariana Vertenstein)
Date: Fri Apr 23 15:50:12 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040423
One-line Summary: added gx3v5 albedos to csim.template
Testing: created t31_gx3v5 and t42_gx1v3 scritps
  and compared csim.buildnml_prestage.csh to verify
  that new albedos were generated correctly	 
Platforms used for testing: sun
Changes made:
	
Added gx3v5 ice albedos to csim.tempate
	
===============================================================
===============================================================

scripts_b040422
Originator: tcraig ( Anthony Craig)
Date: Thu Apr 22 12:57:47 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040422
One-line Summary: Set LOGDIR default to ""
Testing:  TER.01a.T42_gx1v3.X.seaborg, tested with create_newcase
Platforms used for testing:      seaborg                                     
Changes made:

Set LOGDIR default to "" in env_run
Modify testcase_setup to set LOGDIR to $CASEROOT/logs for all tests
Improve error message in ccsm_build.csh when building esmf, mct, 
  mph libs fails

===============================================================
===============================================================

scripts_b040421
Originator: tcraig ( Anthony Craig)
Date: Wed Apr 21 16:56:51 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040421
One-line Summary: Check for existance of binaries, Remove PRESTAGE_DATA parameter
Testing: TER.01a.T42_gx1v3.A.seaborg
         TER.01a.T42_gx1v3.X.seaborg
         TER.01a.T42_gx1v3.B.seaborg
Platforms used for testing:  seaborg                                      
Changes made:

Check for existance of binaries in each component's buildexe script
Remove PRESTAGE_DATA parameter from env_run and all other locations,
  mostly component templates
Fix bug in saving build log files, copy (and run) failed if log/bld 
  did not exist.  ccsm_build now check for directory existance and
  creates it if it doesn't exist.

===============================================================
===============================================================
scripts_b040416
Originator: mvertens ( Mariana Vertenstein)
Date: Fri Apr 16 12:47:51 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040416
One-line Summary: put in fixes to not need . in path
Testing:                                                              
Platforms used for testing: sun (built the scripts)
Changes made:
1) Put in fixes in Machines directory in the run* scripts to
   not need . in user paths
===============================================================
	
===============================================================

scripts_b040412
Originator: tcraig ( Anthony Craig)
Date: Mon Apr 12 12:31:19 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040412
One-line Summary: Remove DOUT_L_MSDIR, add machine to DOUT_L_RCP_ROOT, add restart.tar long term archiving
Testing: Unit testing only on seaborg                                                             
Platforms used for testing:  seaborg
Changes made:

Remove DOUT_L_MSDIR
add mylogin@remotesite.edu: to DOUT_L_RCP_ROOT
add restart.tars files to long term archiving
  will save them to .../$CASE/restart.tars/
  retention period will be 30 
  retention period will be default for files with years ending in 0
Change default runtime on seaborg to 1800 from 7200

===============================================================
===============================================================

scripts_b040408
Originator: tcraig ( Anthony Craig)
Date: Thu Apr  8 13:38:59 MDT 2004
Model: SCRIPTS
Version: SCRIPTSb040408
One-line Summary: Update scripts
Testing:                                                              
Platforms used for testing:                                           
Changes made:

1) fix create_newcase, configure -> ./configure
2) fix configure, mach -> MACH
3) update llq test for harvester submission for all ibms (Buja)
4) remove F case from resok list in check_compset.  still found
   in get_compset.

===============================================================
===============================================================
scripts_b040402
Originator: mvertens ( Mariana Vertenstein)
Date: Fri Apr  2 12:54:37 MST 2004
Model: SCRIPTS
Version: SCRIPTSb040402
One-line Summary: updated to changes needed for  MCT1_9_0a040316
Testing:                                                              
Platforms used for testing:                                           
Changes made:
1) Updated to changes needed for MCT1_9_0a040316
===============================================================
	
===============================================================
scripts_b040330c
Originator: mvertens ( Mariana Vertenstein)
Date: Tue Mar 30 20:24:18 MST 2004
Model: SCRIPTS
Version: SCRIPTSb040330c
One-line Summary: fixed create_test bug
Testing: verified that create_test worked                                                              
Platforms used for testing:                                           
Changes made:
1) removed ./addmach call in Tools/testcase_setup.csh
===============================================================
	
===============================================================
scripts_b040330b
Originator: mvertens ( Mariana Vertenstein)
Date: Tue Mar 30 20:22:46 MST 2004
Model: SCRIPTS
Version: SCRIPTSb040330b
One-line Summary: fixed create_test bug
Testing:                                                              
Platforms used for testing:                                           
Changes made:
1) removed call to ./addmach in Tools/testcase_setup.csh
2) removed output comments added in Case.Template/configure
===============================================================
	
===============================================================
scripts_b040330
Originator: mvertens ( Mariana Vertenstein)
Date: Tue Mar 30 16:16:16 MST 2004
Model: SCRIPTS
Version: SCRIPTSb040330
One-line Summary: removed addmach from script directory, use configure -addmach machine instead
Testing: validated that setting LOGDIR to "" worked correctly
         created numerous test cases with new commands to validate
Platforms used for testing: bluesky                                           
Changes made:
1) removed addmach from Case.template
   added configure -addmach machine instead
   modified create_new case to invoke 
      configure -addmach machine
   in scripts directory instead of 
      addmach -mach
2) added more diagnostic output to configure 
3) modified env_run so that LOGDIR is no longer a derived variable
   although its default value is set as a derived variable
   can now set LOGDIR to "" if do not want standard out to be copied
   back to $LOGDIR
4) DOUT_S is now set to TRUE by default
5) modified all Machines/run.* files to have the following logic:
   # -------------------------------------------------------------------------
   # Save model output stdout and stderr 
   # -------------------------------------------------------------------------

   cd $EXEROOT/cpl 
   set CplLogFile = `ls -1t cpl.log* | head -1` 
   grep 'end of main program' $CplLogFile   || echo "Model did not complete - see $CplLogFile" && exit -1

   cd $EXEROOT
   gzip */*.$LID

   if ($LOGDIR != "") then
     if (! -d $LOGDIR/bld) mkdir -p $LOGDIR/bld || echo " problem in creating $LOGDIR/bld" && exit -1
     cp -p */*buildexe*$LID.* $LOGDIR/bld   || echo "Error in copy of logs " && exit -1
     cp -p */*log*$LID.*      $LOGDIR       || echo "Error in copy of logs " && exit -1
   endif
6) removed references to LOGDIR from ccsm_build.csh
===============================================================
	
===============================================================
scripts_b040326
Originator: mvertens ( Mariana Vertenstein)
Date: Fri Mar 26 14:00:32 MST 2004
Model: SCRIPTS
Version: SCRIPTSb040326
One-line Summary: updated ipcc cam history fields and datasets
Testing:                                                              
Platforms used for testing:                                           
Changes made:
1) updated cam.template for new ipcc mode switches
2) updated cheetah32 task-threads	
===============================================================
	
===============================================================
scripts_b040325
Originator: mvertens ( Mariana Vertenstein)
Date: Thu Mar 25 10:53:38 MST 2004
Model: SCRIPTS
Version: SCRIPTSb040325
One-line Summary: updated topography file in pop
Testing:created scripts and examined output                                                              
Platforms used for testing:ibm
Changes made:
1) added component set O [latm  dlnd dice docn cpl]
2) fixed comments in README and env.README
3) fixed batch script for branch and hybrid to only build 
   reference case interacatively (weiyu)
4) fixed bug in clm.template and cpl.template to
   only obtain restart file if runtype is not 'continue'
5) put in new topography file in pop.template (njn01)
6) updated cheetah, seaborg and bluesky scripts
===============================================================
	
===============================================================
scripts_b040310
Originator: mvertens ( Mariana Vertenstein)
Date: Wed Mar 10 21:59:32 MST 2004
Model: SCRIPTS
Version: SCRIPTSb040310
One-line Summary: fixed archive bug, env.xxx.xxx bugs
Testing: Verified that new xxx.run scripts produced right
task geometry and tasks and threads. Compared all the new
resolved xxx.run files with original resolved xxx.run files 
to validate new env files.	
Platforms used for testing: blackforest, bluesky
Changes made:
1) Fixed bugs in ccsm_s_archive.csh
2) Fixed bugs in env.xxx.machine
3) Added long term archiving for cheetah and cheetah32
   Fixed other problems for cheetah and cheetah32
===============================================================
	
===============================================================
scripts_b040309
Originator: mvertens ( Mariana Vertenstein)
Date: Tue Mar  9 15:33:00 MST 2004
Model: SCRIPTS
Version: SCRIPTSb040309
One-line Summary: removed T200 from cam namelist, fixed branch and hybrid run tests
Testing: Wei Yu verified that BR.02a and HY.02a now work
Platforms used for testing:                                           
Changes made:
1) removed T200 from cam namelist
2) bug fix in create_test for access to subLits
3) bug fix for BR.02a and HY.02a (Wei Yu)
4) added l_archive.cheetah and l_archive.cheetah32 in Machines/	
===============================================================
	
===============================================================
scripts_b040305b
Originator: mvertens ( Mariana Vertenstein)
Date: Thu Mar  4 13:24:27 MST 2004
Model: SCRIPTS
Version: SCRIPTSb040305b
One-line Summary: fixed date bug in co2 ramping test
Testing:                                                              
Platforms used for testing:                                           
Changes made:
1) fixed date bug in co2 ramping test - need date in form
	yyyy-mm-dd (even if yyyy contains a 0xxx-yy-dd)
===============================================================
	
===============================================================
scripts_b040305
Originator: mvertens ( Mariana Vertenstein)
Date: Thu Mar  4 10:25:03 MST 2004
Model: SCRIPTS
Version: SCRIPTSb040305
One-line Summary: added restart tests for ramp_co2 and 2000 ipcc and
	removed old scipt related files
Testing: TER.01e.T42_gx1v3.K.blackforest  
Platforms used for testing: IBM                                          
Changes made:
1) Added restart tests for ramp_co2 and 2000 ipcc. Verified
   that ramping code restarted correctly.
2) Removed the following directories from scripts:
   gui/ gui_run/ test.a1/ test.a2/ tests/ tools/
===============================================================
	
===============================================================
scripts_b040304a
Originator: mvertens ( Mariana Vertenstein)
Date: Wed Mar  3 11:32:21 MST 2004
Model: SCRIPTS
Version: SCRIPTSb040304a
One-line Summary: addition of  restart test for cam ramp co2 code
Testing:                                                              
Platforms used for testing:                                           
Changes made:
fixed bugs in cam.template for ramp co2
added restart test for cam ramp co2 code
===============================================================
	
===============================================================
scripts_b040304
Originator: mvertens ( Mariana Vertenstein)
Date: Wed Mar  3 09:40:21 MST 2004
Model: SCRIPTS
Version: SCRIPTSb040304
One-line Summary: fixed minor bug in blackforest T31 load balance
Testing:                                                              
Platforms used for testing:                                           
Changes made:
fixed minor bug in blackforest T31 load balance
===============================================================
	
===============================================================
scripts_b040302
Originator: mvertens ( Mariana Vertenstein)
Date: Tue Mar  2 08:47:08 MST 2004
Model: SCRIPTS
Version: SCRIPTSb040302
One-line Summary: bug fixes
Testing:                                                              
Platforms used for testing:                                           
Changes made:
1) Bug fixes in csim.template, cam.template and env_mach.blackforest
===============================================================
	
===============================================================
scripts_b040301
Originator: mvertens ( Mariana Vertenstein)
Date: Mon Mar  1 14:21:34 MST 2004
Model: SCRIPTS
Version: SCRIPTSb040301
One-line Summary: updated template files for t31_gx3v5
Testing:                                                              
Platforms used for testing:                                           
Changes made:
1) Updated cam.template, pop.template, csim.template and 
   cpl.template for T31_gx3v5
2) Updated env.ibm.blackforest to have load balanced T31_gx3v5 
   configuration - NOTE only blackforest has been load balanced
   at this time.	
3) Updated cam.template for IPCC configuration
   removed all references to kmxhdc - default will always be used
   removed all references to ozncyc=.true. - this is the default
4) Modified env_conf and cam.template so that co2 ramping is
   an IPCC option	 
5) Put in fix in ccsm_s_archive so that all cam and clm history
   files (not just h0 and h1) will have short timer archiving
===============================================================
	
===============================================================

scripts_b040211
Originator: tcraig ( Anthony Craig)
Date: Wed Feb 11 15:47:58 MST 2004
Model: SCRIPTS
Version: SCRIPTSb040211
One-line Summary: Put gx3v4 back into scripts
Testing:  test scripts generated but not run
Platforms used for testing:  None
Changes made:

Added gx3v4 back into
  Tools/check_res
  Components/clm.template
  Components/pop.template (reviewed only, was still supported)
  Components/csim.template
  Components/cpl.template
  Components/dice.template
  Components/docn.template

Fixed surface dataset bug in clm.template, changed

<    set datasurf = surface-data.128x064_atm.gx1v3_ocn.040209.nc

>    set datasurf = surface-data.096x048_atm.gx3v5_ocn.040209.nc

for T31_gx3v5


===============================================================
===============================================================

scripts_b040210
Originator: tcraig ( Anthony Craig)
Date: Tue Feb 10 14:00:41 MST 2004
Model: SCRIPTS
Version: SCRIPTSb040210
One-line Summary: Several bug fixes and improvements
Testing:                                                              
TER.01a.T31_gx3v5.A.blackforest 
TER.01a.T31_gx3v5.B.blackforest (fails out of box, need pop_in fix)
TER.01a.T31_gx3v5.X.blackforest
TER.01a.T42_gx1v3.B.blackforest
TER.01a.T85_gx1v3.B.bluesky32
TER.01b.T85_gx1v3.B.bluesky32
TER.01c.T85_gx1v3.B.bluesky32
TER.01d.T85_gx1v3.B.bluesky32
T42_gx1v3.B.blackforest, 10 day, info_debug, old/new scripts
T85_gx1v3.B.bluesky, 10 day, info_debug, old/new scripts
Many cases generated and reviewed (not run)

Platforms used for testing:  blackforest, bluesky                                         
Changes made:

- modify test.a1/atm.setup.csh to remove IPCC namelist
- remove "startdate" input in testcase_startup, was dead code,
  "basedate" handles that env_conf input.
- fix M case set, now uses docn
- add source of modules in build scripts, change to generate_batch
- add ER.01c, ER.01d and tested them
- update run.es.moon from Yoshi
- add capability for ocean to be run on 10 or 12 tasks, sets
  decomposition in pop template
- modifications for T31_gx3v5 in docn, dice, csim templates
- modification to cam template for updated datasets for T31 and IPCC cases
- modify copy of modules file in addmach, no more "copy error" if
  modules doesn't exist.
- add optional -testid argument to create_testcase, allows user to
  optionally set the testid.
- added error handling for create_newcase and create_test for unknown
  arguments.  used to just ignore bad inputs
- added modules.cpq.lemieux, required for lemieux
- added build line to batch.$MACH file between cd and llsubmit
- fix src.ocn in pop.template -> src.pop
- updated taskmaker.pl
- added surface forcing file for clm template, T31_gx3v5, B, 
  surface-data.128x064_atm.gx1v3_ocn.040209.nc
- fix datasets for IPCC runs, T42, historical, future
- commented out futureghg dataset for A1, so ER.01d can be tested.  it
  currently uses the wrong dataset, dataset is not yet generated

Verified old and new scripts are bfb for present day T42_gx1v3 and
  T85_gx1v3, present day.  Caveat on divdampn differences in two
  scripts results in bfb.
T31_gx3v5 will run out of the box fine with appropriate changes in pop,
  pop_in
TER.01c, TER.01d work out of the box
pes_file option reviewed, tested (no run) and works in create_test script

===============================================================
===============================================================
===============================================================
