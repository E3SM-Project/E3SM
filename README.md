# cime
Common Infrastructure for Modeling the Earth

**cime** (pronounced: seem) represents the infrastructure code for the 
<a href="http://www2.cesm.ucar.edu">Community Earth System Model </a>
     (CESM) and the Accelerated Climate Model for Energy (ACME).
*cime* includes scripts for configuration, build, and testing of
models, as well as data and stub models for climate components, 
and a code for coupling different climate components together.
     
# Developers

*cime* was initially extracted from CESM as a stand-alone capability in 2015
by members of the CSEG group at NCAR, the software engineering team of
the CESM model. The CSEG group had been developing this infrastrucure within
CESM for a number of years using NSF and DOE funding. After version 4 was released,
a joint development partnership was started with the software engineering group of
ACME, a DOE-funded project, which had branched from CESM in 2014. 
Starting with development of version 5, *cime* is cooperative effort with contributions
and ownership from members of both teams. 

The following table documents the developers who have contributed to *cime*,
showing what versions of they contributed to, and with what source(s) of support.

Name     |     Institution  |  Versions  | Funding Source (for versions)
---------|------------------|----------|----------------------
Mariana Vertenstein | NCAR  |  1 - D   |   NSF, DOE
Jim Edwards |         NCAR  |  1 - D   |   NSF (1-D), DOE(1-2)
Jim Foucar  |         SNL   |  5 - D   |   DOE
Rob Jacob |      ANL  |  5 - D   |   DOE
Andreas Wilke |  ANL  | 5 - D | DOE
Jason Sarich | ANL | 5 - D | DOE
??Sean Santos |         NCAR  |  1 - 4   |   NSF

_Key: D = Current development version (i.e. still active on project)_


# Documentation

See the CESM web site for documentation and information.

The cime directory structure (other than the externals/ directory)
was generated from the following subversion trunk tags:

* cime/driver_cpl		  	        
   https://svn-ccsm-models.cgd.ucar.edu/drv/seq_mct/trunk_tags/drvseq5_1_15

* cime/components/data_comps/datm	        
   https://svn-ccsm-models.cgd.ucar.edu/datm7/trunk_tags/datm8_150310
   
* cime/components/data_comps/dice	        
   https://svn-ccsm-models.cgd.ucar.edu/dice7/trunk_tags/dice8_150310
   
* cime/components/data_comps/dlnd	        
   https://svn-ccsm-models.cgd.ucar.edu/dlnd7/trunk_tags/dlnd8_150310
   
* cime/components/data_comps/docn       	
   https://svn-ccsm-models.cgd.ucar.edu/docn7/trunk_tags/docn8_150310
   
* cime/components/data_comps/drof  	
   https://svn-ccsm-models.cgd.ucar.edu/drof/trunk_tags/drof_150310
   
* cime/components/stub_comps		
   https://svn-ccsm-models.cgd.ucar.edu/stubs/trunk_tags/stubs1_4_08
   
* cime/components/xcpl_comps	
   https://svn-ccsm-models.cgd.ucar.edu/dead7/trunk_tags/dead7_8_04

* cime/machines				
   https://svn-ccsm-models.cgd.ucar.edu/Machines/trunk_tags/Machines_150309

* cime/scripts                            
   https://svn-ccsm-models.cgd.ucar.edu/scripts/trunk_tags/scripts_150309

* cime/share/csm_share	  	        
   https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_150116
   
* cime/share/esmf_wrf_timemgr	        
   https://svn-ccsm-models.cgd.ucar.edu/esmf_wrf_timemgr/trunk_tags/esmf_wrf_timemgr_141217
   
* cime/share/timing                       
   https://svn-ccsm-models.cgd.ucar.edu/timing/trunk_tags/timing_150302

* cime/utils/pythonlib    
   https://svn-ccsm-models.cgd.ucar.edu/scripts/trunk_tags/scripts4_150204/scripts/ccsm_utils/Tools/pythonlib
   
* cime/utils/perl5lib	                
   https://svn-ccsm-models.cgd.ucar.edu/perl5lib/trunk_tags/perl5lib_150302

* cime/tools/load_balancing_tool	
   https://svn-ccsm-models.cgd.ucar.edu/tools/load_balancing_tool/trunk_tags/load_balancing_tool_141008
   
* cime/tools/unit_testing                 
   https://svn-ccsm-models.cgd.ucar.edu/unit_testing/trunk_tags/unit_testing_0_16
   
* cime/tools/statistical_ensemble_test  
   https://svn-ccsm-models.cgd.ucar.edu/validation_testing/trunk_tags/validation_20140708/run_CESM/
   
* cime/tools/mapping                      
   https://svn-ccsm-models.cgd.ucar.edu/tools/mapping/trunk_tags/mapping_141106
   
* cime/tools/cprnc                        
   https://svn-ccsm-models.cgd.ucar.edu/tools/cprnc/trunk_tags/cprnc_150301



