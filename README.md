Running the HSW++ Volcanic Aerosol Injection
================================================================================

This branch maintains the current version of CLDERA's idealized "HSW++" model for
tiered V&V. To run E3SM in this configuration, follow the steps below. 

### 1. 
Clone the repo
```
git clone --recursive git@github.com:sandialabs/CLDERA-E3SM.git -b jhollowed/eam/cldera-sai-module
```

### 2. 
Copy the template case creation script and user namelist settings
```
cp -r /project/projectdirs/m4014/data/HSW/run_scripts location/of/your/directory/
```

### 3. 
Edit the case creation script to point to custom directories for case setup and output directories, and location of the cloned model, namely `MY_E3SM_ROOT`, `CASE_ROOT`, `OUT_ROOT`. Also make any desired edits to the horizontal grid definition, etc.

### 4. 
Edit `user_nl_eam` to provide any desired namelist changes. Alternative namelist files can be provded via an argument to the case creation script. See inline comments for usage.

### 5. 
Run the model (by default runs for 30 days, but run length can be controlled via an argument to the case creation script. see sinline comments for usage.)
```
./create_FIDEAL_ne16pg2L72_SAI.sh
```


### Additional notes
The `--recrusive` flag to `git clone` also clones all sub-modules (CIME, etc.), which may take several minutes. Once complete, cases are ready to be created. This occurs in the standard way, and users familiar with the CIME framework may continue in standard practice. Or, a template case creation script is located on NERSC at 

`/project/projectdirs/m4014/data/HSW/run_scripts/create_FIDEAL_ne16pg2L72_SAI.sh`

Inline comments explain the various operations in detail, as well as a `README` file present at this location. The only line here that will differ from typical case creation is 

`./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS --val "-cldera_sai_trcs -nlev 72"`

where the `cldera_sai_trcs` configure option to EAM activates the aerosol injection. The call to `create_case` also recieves the option `--compset FIDEAL` which selects the idealized Held-Suarez physics package. This compset is the namesake of HSW++, and the heating parameters of the active aerosols have been tuned against the resulting climatology, but the injection can still be activated with any other more complex compset. 

The script also copies a template namelist setting file to the newly created case directory, from 

`/project/projectdirs/m4014/data/HSW/run_scripts/user_nl_eam`

which contains namelist settings controlling the parameters of the aerosol injection, and the subsequent heating processes. Namely, this is where the various pathway nodes are turned on/off via `cldera_sai_formSulfate`, `cldera_sai_stratHeating`, `cldera_sai_surfCooling`. Inline comments explain further. 

-----------------


[![E3SM Logo](https://e3sm.org/wp-content/themes/e3sm/assets/images/e3sm-logo.png)](https://e3sm.org)

Energy Exascale Earth System Model (E3SM)
================================================================================

E3SM is a state-of-the-art fully coupled model of the Earth's climate including
important biogeochemical and cryospheric processes. It is intended to address
the most challenging and demanding climate-change research problems and
Department of Energy mission needs while efficiently using DOE Leadership
Computing Facilities.  

DOI: [10.11578/E3SM/dc.20210927.1](http://dx.doi.org/10.11578/E3SM/dc.20210927.1)

Please visit the [project website](https://e3sm.org) or our [Confluence site](https://acme-climate.atlassian.net/wiki/spaces/DOC/overview)
for further details.

For questions about the model, use [Github Discussions](https://github.com/E3SM-Project/E3SM/discussions)

Table of Contents 
--------------------------------------------------------------------------------
- [Quick Start](#quickstart)
- [Supported Machines](#supportedmachines)
- [Running](#running)
- [Contributing](#contributing)
- [Acknowledge](#acknowledge)
- [License](#license)

Quick Start
--------------------------------------------------------------------------------
The [Quick Start](https://e3sm.org/model/running-e3sm/e3sm-quick-start/) page 
includes instructions on obtaining the necessary code and input data for model 
setup and execution on a supported machine.

Supported Machines 
--------------------------------------------------------------------------------
E3SM is a high-performance computing application and generally requires a
capable compute cluster to run a scientifically validated case at a useful
simulation speed.

To run E3SM, it is recommended that you obtain time on a 
[Supported Machine](https://e3sm.org/model/running-e3sm/supported-machines/).

Running
--------------------------------------------------------------------------------
Please refer to [Running E3SM](https://e3sm.org/model/running-e3sm/) 
 for instructions on running the model. 

Contributing
--------------------------------------------------------------------------------
Please refer to [Contributing](CONTRIBUTING.md) for details on our code development
process.

Acknowledgement
--------------------------------------------------------------------------------
The Energy Exascale Earth System Model (E3SM) Project should be acknowledged in
publications as the origin of the model using
[these guidelines](https://e3sm.org/resources/policies/acknowledge-e3sm/).

In addition, the software should be cited.  For your convenience,
the following BibTeX entry is provided.
```TeX
@misc{e3sm-model,
	title = {{Energy Exascale Earth System Model (E3SM)}},
	author = {{E3SM Project}},
	abstractNote = {{E3SM} is a state-of-the-art fully coupled model of the {E}arth's 
		climate including important biogeochemical and cryospheric processes.},
	howpublished = {[Computer Software] \url{https://dx.doi.org/10.11578/E3SM/dc.20210927.1}},
	url = {https://dx.doi.org/10.11578/E3SM/dc.20210927.1},
	doi = {10.11578/E3SM/dc.20210927.1},
	year = 2021,
	month = sep,
}
```

License
--------------------------------------------------------------------------------
The E3SM model is available under a BSD 3-clause license.
Please see [LICENSE](LICENSE) for details.

