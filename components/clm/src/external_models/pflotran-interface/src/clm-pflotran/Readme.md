# INSTRUCTIONS #

This is the interface, including minimal PFLOTRAN codes, for coupling PFLOTRAN into CLM in DOE sponsored NGEE-Arctic Project. 

The model coupling aims to provide a full alternative solution for CLM-CN's surface-subsurface C/N biogeochemistry and thermal-hydrology, i.e. PFLOTRAN.

***This interface is specifically released for ```PETSC v.3.8x```, and ```E3SM v1.1```. ***


## How do I get set up? ##

**(1)** *git clone the repository*, if not yet properly updating codes. (OPTIONAL)
```
git clone https://github.com/fmyuan/pflotran-elm-interface/src/clm-pflotran
```

*NOTE: when git clone, checkout branch 'default-release-v3.8'*


**(2)** *if coupling CLM with PFLOTRAN, need to build a library named as **libpflotran.a**.*
```
cd ./src/clm-pflotran
```

***FIRST***, run makefile to copy needed PFLOTRAN source code files (*.F90), if PFLOTRAN codes missing. (OPTIONAL)
```
make PETSC_DIR=$PETSC_DIR copy_common_src

(OR, make PETSC_DIR=$PETSC_DIR clean_common_src 
to clean PFLOTRAN source code files, and only leave CLM-PFLOTRAN interface codes)
```

***SECONDLY***, build the library
```
make PETSC_DIR=$PETSC_DIR smoothing2=TRUE libpflotran.a

(OR, make PETSC_DIR=$PETSC_DIR smoothing2=TRUE debugbuild=TRUE libpflotran.a
for a library with '-g -O0')

```

***FINALLY***, build ELM v1.1 with this library, as usual.




## A FEW Specific Notes on How to modify ELM Macro.make and other machine files##
*I.* **Macro (CLM)** or **Macro.make (ALM)** modified for coupling build -
```
ifeq ($(MODEL), clm) 
  ifeq ($(CLM_PFLOTRAN_COLMODE), TRUE) 
    ifeq ($(CLM_PFLOTRAN_COUPLED), TRUE) 
       CPPDEFS +=  -DCOLUMN_MODE 
    endif
  endif

  ifeq ($(CLM_PFLOTRAN_COUPLED), TRUE) 
     FFLAGS +=  -I$(CLM_PFLOTRAN_SOURCE_DIR)
     CPPDEFS +=  -DCLM_PFLOTRAN 
  endif
endif

......

ifeq ($(MODEL), driver) 
  ifeq ($(CLM_PFLOTRAN_COUPLED), TRUE) 
     LDFLAGS +=  -L$(CLM_PFLOTRAN_SOURCE_DIR) -lpflotran $(PETSC_LIB)
  endif
endif

```

*NOTE*: Modified Macro above requires ***4*** alias, setted as following.

*II.* **Makefile**
```
# Set PETSc info if it is being used
ifeq ($(strip $(USE_PETSC)), TRUE)
  ifdef PETSC_PATH
    ifndef INC_PETSC
      INC_PETSC:=$(PETSC_PATH)/include
    endif
    ifndef LIB_PETSC
      LIB_PETSC:=$(PETSC_PATH)/lib
    endif
  else
    $(error PETSC_PATH must be defined when USE_PETSC is TRUE)
  endif

  # Get the "PETSC_LIB" list an env var (UPDATED: 2017-May-19)
  # include $(PETSC_PATH)/conf/variables
  # (1) petsc-git-version: 1a9d3c3c50abf60098813fdf7291fe3540415115
  # (2) in this petsc package, "PETSC_LIB" contains "-L$LIB_PETSC"
  include $(PETSC_PATH)/lib/petsc/conf/variables
  
endif

```

*III.* **config_compilers.xml** editing for all or specific compilers. *NOTE*: if this added, No need to modify 'Macro' or 'Macro.make' under case directory. 

```
<!-- hacking of mach/compiler generated 'Macros.make' for coupling with pflotran -->
<!-- ideally it should go with CLM configuration -->
<compiler>
  <ADD_FFLAGS MODEL="clm" CLM_PFLOTRAN_COUPLED="TRUE"> -I$(CLM_PFLOTRAN_SOURCE_DIR)</ADD_FFLAGS>
  <ADD_CPPDEFS MODEL="clm" CLM_PFLOTRAN_COUPLED="TRUE"> -DCLM_PFLOTRAN </ADD_CPPDEFS>
  <ADD_CPPDEFS MODEL="clm" CLM_PFLOTRAN_COUPLED="TRUE" CLM_PFLOTRAN_COLMODE="TRUE"> -DCOLUMN_MODE </ADD_CPPDEFS>
  <ADD_LDFLAGS MODEL="driver" CLM_PFLOTRAN_COUPLED="TRUE"> -L$(CLM_PFLOTRAN_SOURCE_DIR) -lpflotran $(PETSC_LIB)</ADD_LDFLAGS>
</compiler>
<!-- end of hacking 'Macros.make' for coupling with pflotran -->

```


*IV.* **config_machines.xml** editing for each supported machine. *NOTE*: after './case.setup', edit 'env_mach_specific.xml' to turn on options.

```
      <!-- for CLM-PFLOTRAN coupling, the PETSC_PATH must be defined specifically upon machines -->
      <environment_variables>
        <env name="PETSC_PATH" compiler="gnu" mpilib="openmpi">/software/user_tools/current/cades-ccsi/petsc4pf/openmpi-1.10-gcc-5.3</env>      
        <!-- hack for PFLOTRAN coupling (this is a temporary solution, and user must manually edit env_mach_specific.xml after case.setup)-->
        <env name="CLM_PFLOTRAN_COUPLED">FALSE</env>
        <env name="CLM_PFLOTRAN_COLMODE">FALSE</env>
        <env name="CLM_PFLOTRAN_SOURCE_DIR">/lustre/or-hydra/cades-ccsi/proj-shared/models/pflotran-interface/src/clm-pflotran</env>
      </environment_variables>       

```

***UPDATED: 2018-10-05***