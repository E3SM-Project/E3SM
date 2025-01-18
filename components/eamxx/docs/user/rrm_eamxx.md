# Running EAMxx with a Regionally Refined Mesh (RRM)

Running EAMxx with a RRM allows you run a select region of the globe at high resolution (i.e. 3 km) with the remainder of the globe at a lower resolution (i.e. 25 or 100 km).  This document will point you to the steps required and resources available to assist in developing and running a new RRM.

## Choose Your RRM

What region of the globe do you want to refine?  Your first step should be to check [library of RRM grids/cases](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/3690397775/Library+of+Regionally-Refined+Model+RRM+Grids) that have already been developed to potentially avoid duplicate work.  If you found a RRM that suits your needs, you can skip the next step ("Generate Your RRM").

## Generate Your RRM

Please refer to the offical [E3SM guide for developing new atmosphere grids](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/872579110/Running+E3SM+on+New+Atmosphere+Grids), which provides detailed guidance for developing your RRM.

After you have made all the necessary files for your RRM, you will need to configure your code branch so that it knows about your new grid.  The steps required to do this are documented at the top of the [library of RRM grids/cases page](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/3690397775/Library+of+Regionally-Refined+Model+RRM+Grids).

## Make Your Initial Condition File

The easiest way to generate an initial condition is to use the [HICCUP tool](https://github.com/E3SM-Project/HICCUP), which is a set of flexible and robust python routines to streamline the task of generating a new atmospheric initial condition for E3SM.  Otherwise, please see the step-by-step instuctions if you prefer to [manually generate your initial condition](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/1002373272/Generate+atm+initial+condition+from+analysis+data).

## Assemble Nudging Data (Optional)

If you wish to nudge your simulation, assemble your nudging data in the format required by EAMxx.  Please refer to the [nudging documentation](nudging.md).

In the event that you only want to nudge a portion of your domain, then you will need to generate a nudging weights file.  A common use case for this is when you want the high-resolution region to remain free-running (unnudged) while nudging the coarse domain towards reanalysis or model data.  Please use [this script](https://github.com/E3SM-Project/eamxx-scripts/blob/master/run_scripts/RRM_example_scripts/SCREAMv1_create_nudging_weights.py) as an example of how to generate your nudging weights file.

## Run your RRM

Congratulations, you are now ready to run your EAMxx RRM.  If you are running your RRM in free running mode (not using any nudging) then you simply need to modify an existing EAMxx script and change the resolution to match the one you created for your RRM.

If you are using nudging, then please see this [example script of how to run a nudged EAMxx RRM run](https://github.com/E3SM-Project/eamxx-scripts/blob/master/run_scripts/RRM_example_scripts/SCREAMv1-nudging.CAx32v1pg2.pm-gpu.template.sh).  This example script uses the [California 3-km RRM](https://gmd.copernicus.org/articles/17/3687/2024/), which is on the master branch.
