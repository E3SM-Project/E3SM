#!/bin/csh

# For now, manually build the namelist. Soon this will call the standard CESM
# build-namelist script.

mkdir -p $CASEROOT/CaseDocs

if !(-d $EXEROOT/ocn/obj   ) mkdir -p $EXEROOT/ocn/obj    || exit 2
if !(-d $EXEROOT/ocn/source) mkdir -p $EXEROOT/ocn/source || exit 3 

if !(-d $CASEBUILD/mpas-oconf) mkdir -p $CASEBUILD/mpas-oconf || exit 1
cd $CASEBUILD/mpas-oconf || exit -1

set inst_string = ""
set STREAM_NAME = "streams.ocean"
set NML_NAME = "mpaso_in"

if (-e $CASEROOT/user_nl_mpaso${inst_string}) then
	$UTILROOT/Tools/user_nlcreate                                           \
		-user_nl_file $CASEROOT/user_nl_mpaso${inst_string}                 \
		-namelist_name mpaso_inparm >! $CASEBUILD/mpas-oconf/cesm_namelist 
endif

# Check to see if "-preview" flag should be passed
if ( $?PREVIEW_NML ) then
	set PREVIEW_FLAG = "-preview"
else
	set PREVIEW_FLAG = ""
endif

# Check to see if build-namelist exists in SourceMods
if (-e $CASEROOT/SourceMods/src.mpas-o/build-namelist) then
	set BLD_NML_DIR = $CASEROOT/SourceMods/src.mpas-o
	set CFG_FLAG = "-cfg_dir $CODEROOT/ocn/mpas-o/bld"
else
	set BLD_NML_DIR = $CODEROOT/ocn/mpas-o/bld
	set CFG_FLAG = ""
endif

# Define input_mesh file and graph prefix for mesh
if ( $OCN_GRID == 'oEC60to30' ) then
	set date_stamp = 150107
	set input_mesh = $DIN_LOC_ROOT/ocn/mpas-o/$OCN_GRID/ocean.EC.60-30km.${date_stamp}.nc
	set graph_prefix = $DIN_LOC_ROOT/ocn/mpas-o/$OCN_GRID/mpas-o.graph.info.${date_stamp}
else if ( $OCN_GRID == 'mpas120' ) then
	set date_stamp = 121116
	set input_mesh = $DIN_LOC_ROOT/ocn/mpas-o/$OCN_GRID/ocean120km.${date_stamp}.nc
	set graph_prefix = $DIN_LOC_ROOT/ocn/mpas-o/$OCN_GRID/mpas-o.graph.info.${date_stamp}
endif

# Write mpas-o.input_data_list file
echo "mesh = $input_mesh" > $CASEBUILD/mpas-o.input_data_list
#echo "graph1 = $graph_prefix" >> $CASEBUILD/mpas-o.input_data_list
echo "graph$NTASKS_OCN = $graph_prefix.part.$NTASKS_OCN" >> $CASEBUILD/mpas-o.input_data_list

set MPAS_STREAMS = $CASEBUILD/mpas-oconf/$STREAM_NAME
touch $MPAS_STREAMS
chmod 644 $MPAS_STREAMS

# Write streams file, if there isn't one in SourceMods
if (-e $CASEROOT/SourceMods/src.mpas-o/$STREAM_NAME) then
	cp $CASEROOT/SourceMods/src.mpas-o/$STREAM_NAME $MPAS_STREAMS
else
	cat >! $MPAS_STREAMS << 'EOF'
	<streams>

	<immutable_stream name="mesh"
					  type="none"
					  filename_template="mesh_variables.nc"
	/>

	<immutable_stream name="input"
					  type="input"
'EOF'

	# Breaking file to insert input file location.
	cat >>! $MPAS_STREAMS << EOF
					  filename_template="$input_mesh"
EOF

	cat >>! $MPAS_STREAMS << 'EOF'
					  input_interval="initial_only"/>

	<!--
	The restart stream is actually controlled via the coupler.
	Changing output_interval here will not have any affect on
	the frequency restart files are written.

	Changing the output_interval could cause loss of data.

	The output_interval is set to 1 second to ensure each restart frame has a
	unique file.
	-->
	<immutable_stream name="restart"
					  type="input;output"
					  filename_template="rst.ocn.$Y-$M-$D_$h.$m.$s.nc"
					  filename_interval="output_interval"
					  reference_time="0000-01-01_00:00:00"
					  clobber_mode="truncate"
					  input_interval="initial_only"
					  output_interval="00-00-00_00:00:01"/>

	<!--
	output is the main history output stream. You can add auxiliary streams to
	this stream to include more fields.
	-->

	<stream name="output"
			type="output"
			filename_template="hist.ocn.$Y-$M-$D_$h.$m.$s.nc"
			filename_interval="01-00-00_00:00:00"
			reference_time="0000-01-01_00:00:00"
			clobber_mode="truncate"
			output_interval="00-01-00_00:00:00">

		<stream name="mesh"/>
		<stream name="real_world"/>
		<var_array name="tracers"/>
		<var name="layerThickness"/>
		<var name="ssh"/>
		<var name="maxLevelEdgeTop"/>
		<var name="vertCoordMovementWeights"/>
		<var name="edgeMask"/>
		<var name="vertexMask"/>
		<var name="cellMask"/>
		<var name="refZMid"/>
		<var name="refLayerThickness"/>
		<var name="xtime"/>
		<var name="zMid"/>
		<var name="zTop"/>
		<var name="kineticEnergyCell"/>
		<var name="relativeVorticityCell"/>
		<var name="areaCellGlobal"/>
		<var name="areaEdgeGlobal"/>
		<var name="areaTriangleGlobal"/>
		<var name="volumeCellGlobal"/>
		<var name="volumeEdgeGlobal"/>
		<var name="CFLNumberGlobal"/>

	</stream>

	<!--
	Streams between this line and the auxiliary stream line below are analysis member streams.
	They can be used to perform online analysis of the simulation and control the output of
	the analysis data.
	-->

	<stream name="globalStatsOutput"
			type="output"
			filename_template="analysis_members/ocn.globalStats.$Y-$M-$D_$h.$m.$s.nc"
			filename_interval="01-00-00_00:00:00"
			packages="amGlobalStats"
			clobber_mode="truncate"
			output_interval="0010_00:00:00">

		<var_array name="minGlobalStats"/>
		<var_array name="maxGlobalStats"/>
		<var_array name="sumGlobalStats"/>
		<var_array name="rmsGlobalStats"/>
		<var_array name="avgGlobalStats"/>
		<var_array name="vertSumMinGlobalStats"/>
		<var_array name="vertSumMaxGlobalStats"/>
		<var name="xtime"/>

	</stream>

	<stream name="zonalMeanOutput"
			type="output"
			filename_template="analysis_members/ocn.zonalMeans.$Y-$M-$D_$h.$m.$s.nc"
			filename_interval="01-00-00_00:00:00"
			packages="amZonalMean"
			clobber_mode="truncate"
			output_interval="0000_12:00:00">

		<var_array name="tracersZonalMean"/>
		<var name="xtime"/>
		<var name="binCenterZonalMean"/>
		<var name="binBoundaryZonalMean"/>
		<var name="velocityZonalZonalMean"/>
		<var name="velocityMeridionalZonalMean"/>

	</stream>

	<stream name="sfcAreaOutput"
	        type="output"
	        filename_template="analysis_members/surface_area_weighted_averages.$Y-$M-$D_$h.$m.$s.nc"
	        filename_interval="01-00-00_00:00:00"
	        clobber_mode="truncate"
	        packages="amSfcAreaWeightedAvgPkg"
	        output_interval="00-00-05_00:00:00" >
	
	    <var_array name="minValueWithinOceanRegion"/>
	    <var_array name="maxValueWithinOceanRegion"/>
	    <var_array name="avgValueWithinOceanRegion"/>
	    <var name="xtime"/>
	</stream>
	
	<stream name="waterMassOutput"
	        type="output"
	        filename_template="analysis_members/water_mass_census.$Y-$M-$D_$h.$m.$s.nc"
	        filename_interval="01-00-00_00:00:00"
	        clobber_mode="truncate"
	        packages="amWaterMassCensusPkg"
	        output_interval="00-00-05_00:00:00" >
	
	    <var_array name="waterMassCensusTemperatureValues"/>
	    <var_array name="waterMassCensusSalinityValues"/>
	    <var_array name="waterMassFractionalDistribution"/>
	    <var_array name="potentialDensityOfTSDiagram"/>
	    <var_array name="zPositionOfTSDiagram"/>
	    <var name="xtime"/>
	</stream>
	
	<stream name="layerVolOutput"
	        type="output"
	        filename_template="analysis_members/layer_volume_weighted_averages.$Y-$M-$D_$h.$m.$s.nc"
	        filename_interval="01-00-00_00:00:00"
	        clobber_mode="truncate"
	        packages="amLayerVolWeightedAvgPkg"
	        output_interval="00-00-05_00:00:00" >
	
	    <var_array name="minValueWithinOceanLayerRegion"/>
	    <var_array name="maxValueWithinOceanLayerRegion"/>
	    <var_array name="avgValueWithinOceanLayerRegion"/>
	    <var_array name="minValueWithinOceanVolumeRegion"/>
	    <var_array name="maxValueWithinOceanVolumeRegion"/>
	    <var_array name="avgValueWithinOceanVolumeRegion"/>
	    <var name="xtime"/>
	</stream>

	<stream name="okuboWeissOutput"
	        type="output"
	        filename_template="analysis_members/okuboWeiss.$Y-$M-$D_$h.$m.$s.nc"
	        filename_interval="01-00-00_00:00:00"
	        clobber_mode="truncate"
	        packages="amOkuboWeiss"
	        output_interval="00-00-05_00:00:00" >
	
	    <stream name="mesh"/>
	    <var name="xtime"/>
	    <var name="okuboWeiss"/>
	    <var name="vorticity"/>
	    <var name="eddyID"/>
	</stream>
	
	<stream name="MerHeatTransOutput"
	        type="output"
	        filename_template="analysis_members/meridional_heat_transport.$Y-$M-$D_$h.$m.$s.nc"
	        filename_interval="01-00-00_00:00:00"
	        reference_time="0000-01-01_00:00:00"
	        clobber_mode="truncate"
	        packages="MerHeatTransPkg"
	        output_interval="0001_00:00:00" >
	
	    <var name="xtime"/>
	    <var name="binBoundaryMerHeatTrans"/>
	    <var name="meridionalHeatTransportLatZ"/>
	    <var name="meridionalHeatTransportLat"/>
	    <var name="refZMid"/>
	    <var name="refBottomDepth"/>
	</stream>



	<!--
	All streams below this line are auxiliary streams. They are provided as
	groupings of fields that one might be interested in. You can either enable the
	stream to write a file for the fileds, or add the stream to another stream that
	will already be written.  
	-->

	<stream name="additional_output"
			type="none"
			filename_template="ocn.additional_output.$Y-$M-$D_$h.$m.$s.nc"
			filename_interval="01-00-00_00:00:00"
			reference_time="0000-01-01_00:00:00"
			clobber_mode="truncate"
	>

		<var name="normalVelocity"/>
		<var name="density"/>
		<var name="pressure"/>
		<var name="divergence"/>
		<var name="viscosity"/>
		<var name="vertViscTopOfEdge"/>
		<var name="vertViscTopOfCell"/>
		<var name="vertDiffTopOfCell"/>
		<var name="BruntVaisalaFreqTop"/>
		<var name="RiTopOfCell"/>
		<var name="bulkRichardsonNumber"/>
		<var name="vertAleTransportTop"/>
		<var name="vertVelocityTop"/>

	</stream>

	<stream name="real_world"
			type="none"
			filename_template="ocn.real_world_variables.$Y-$M-$D_$h.$m.$s.nc"
			filename_interval="01-00-00_00:00:00"
			reference_time="0000-01-01_00:00:00"
			clobber_mode="truncate"
	>

		<stream name="mesh"/>
		<var name="normalVelocity"/>
		<var name="velocityZonal"/>
		<var name="velocityMeridional"/>
		<var name="displacedDensity"/>
		<var name="potentialDensity"/>
		<var name="boundaryLayerDepth"/>
		<var name="boundaryLayerDepthEdge"/>
		<var name="indexBoundaryLayerDepth"/>
		<var name="indexSurfaceLayerDepth"/>
		<var name="surfaceFrictionVelocity"/>
		<var name="windStressZonalDiag"/>
		<var name="windStressMeridionalDiag"/>
		<var name="surfaceBuoyancyForcing"/>
		<var name="seaSurfacePressure"/>

	</stream>

	<stream name="averages"
			type="none"
			filename_template="ocn.average_variables.$Y-$M-$D_$h.$m.$s.nc"
			filename_interval="01-00-00_00:00:00"
			reference_time="0000-01-01_00:00:00"
			clobber_mode="truncate"
	>

		<stream name="mesh"/>
		<var_array name="avgTracersSurfaceValue"/>
		<var_array name="avgSurfaceVelocity"/>
		<var_array name="avgSSHGradient"/>
		<var name="nAverage"/>
		<var name="avgSSH"/>
		<var name="avgNormalVelocity"/>
		<var name="avgVelocityZonal"/>
		<var name="avgVelocityMeridional"/>
		<var name="avgVertVelocityTop"/>
		<var name="avgNormalTransportVelocity"/>
		<var name="avgTransportVelocityZonal"/>
		<var name="avgTransportVelocityMeridional"/>
		<var name="avgVertTransportVelocityTop"/>
		<var name="varSSH"/>
		<var name="varNormalVelocity"/>
		<var name="varVelocityZonal"/>
		<var name="varVelocityMeridional"/>

	</stream>

	<stream name="Cartesian"
			type="none"
			filename_template="ocn.Cartesian_variables.$Y-$M-$D_$h.$m.$s.nc"
			filename_interval="01-00-00_00:00:00"
			reference_time="0000-01-01_00:00:00"
			clobber_mode="truncate"
	>

		<stream name="mesh"/>
		<var name="velocityX"/>
		<var name="velocityY"/>

	</stream>

	<stream name="forcing"
			type="none"
			filename_template="ocn.forcing_variables.$Y-$M-$D_$h.$m.$s.nc"
			filename_interval="01-00-00_00:00:00"
			reference_time="0000-01-01_00:00:00"
			clobber_mode="truncate"
	>

		<stream name="mesh"/>
		<var_array name="tracersSurfaceValue"/>
		<var_array name="surfaceVelocity"/>
		<var_array name="SSHGradient"/>
		<var_array name="surfaceTracerFlux"/>
		<var_array name="vertNonLocalFlux"/>
		<var name="surfaceWindStressMagnitude"/>
		<var name="surfaceWindStress"/>
		<var name="surfaceThicknessFlux"/>
		<var name="seaIceEnergy"/>
		<var name="penetrativeTemperatureFlux"/>
		<var name="fractionAbsorbed"/>
		<var name="windStressZonal"/>
		<var name="windStressMeridional"/>
		<var name="latentHeatFlux"/>
		<var name="sensibleHeatFlux"/>
		<var name="longWaveHeatFluxUp"/>
		<var name="longWaveHeatFluxDown"/>
		<var name="seaIceHeatFlux"/>
		<var name="shortWaveHeatFlux"/>
		<var name="evaporationFlux"/>
		<var name="seaIceSalinityFlux"/>
		<var name="seaIceFreshWaterFlux"/>
		<var name="riverRunoffFlux"/>
		<var name="iceRunoffFlux"/>
		<var name="rainFlux"/>
		<var name="snowFlux"/>
		<var name="iceFraction"/>
		<var name="prognosticCO2"/>
		<var name="diagnosticCO2"/>
		<var name="squaredWindSpeed10Meter"/>
		<var name="CO2Flux"/>
		<var name="DMSFlux"/>
		<var name="nAccumulatedCoupled"/>
		<var name="thermalExpansionCoeff"/>
		<var name="salineContractionCoeff"/>

	</stream>

	<stream name="Gent_McWilliams_spherical"
			type="none"
			filename_template="ocn.Gent_McWilliams_spherical_variables.$Y-$M-$D_$h.$m.$s.nc"
			filename_interval="01-00-00_00:00:00"
			reference_time="0000-01-01_00:00:00"
			clobber_mode="truncate"
	>

		<stream name="mesh"/>
		<var name="relativeSlopeTopOfCell"/>
		<var name="relativeSlopeTaperingCell"/>
		<var name="relativeSlopeTopOfCellZonal"/>
		<var name="relativeSlopeTopOfCellMeridional"/>
		<var name="k33"/>
		<var name="GMBolusVelocityZonal"/>
		<var name="GMBolusVelocityMeridional"/>
		<var name="normalGMBolusVelocity"/>
		<var name="vertGMBolusVelocityTop"/>
		<var name="gmStreamFuncTopOfEdge"/>
		<var name="avgNormalGMBolusVelocity"/>
		<var name="avgGMBolusVelocityZonal"/>
		<var name="avgGMBolusVelocityMeridional"/>
		<var name="avgVertGMBolusVelocityTop"/>

	</stream>

	<stream name="Gent_McWilliams_Cartesian"
			type="none"
			filename_template="ocn.Gent_McWilliams_Cartesian_variables.$Y-$M-$D_$h.$m.$s.nc"
			filename_interval="01-00-00_00:00:00"
			reference_time="0000-01-01_00:00:00"
			clobber_mode="truncate"
	>

		<stream name="mesh"/>
		<var name="relativeSlopeTopOfCell"/>
		<var name="relativeSlopeTaperingCell"/>
		<var name="relativeSlopeTopOfCellX"/>
		<var name="relativeSlopeTopOfCellY"/>
		<var name="relativeSlopeTopOfCellZ"/>
		<var name="k33"/>
		<var name="GMBolusVelocityX"/>
		<var name="GMBolusVelocityY"/>
		<var name="normalGMBolusVelocity"/>
		<var name="vertGMBolusVelocityTop"/>
		<var name="gmStreamFuncTopOfEdge"/>
		<var name="avgNormalGMBolusVelocity"/>
		<var name="GMStreamFuncX"/>
		<var name="GMStreamFuncY"/>

	</stream>



	</streams>
'EOF'
endif # Writing streams file

$BLD_NML_DIR/build-namelist $CFG_FLAG $PREVIEW_FLAG                        \
		-infile $CASEBUILD/mpas-oconf/cesm_namelist                        \
		-caseroot $CASEROOT                                                \
		-casebuild $CASEBUILD                                              \
		-scriptsroot $SCRIPTSROOT                                          \
		-inst_string "$inst_string"                                        \
		-date_stamp "$date_stamp"                                          \
		-ocn_grid "$OCN_GRID" || exit -1  

if ( -d ${RUNDIR} ) then
	cp $CASEBUILD/mpas-oconf/mpaso_in ${RUNDIR}/mpaso_in
endif

/bin/cp $MPAS_STREAMS $RUNDIR

if ( -e $CASEROOT/CaseDocs/$STREAM_NAME ) then
	chmod 644 $CASEROOT/CaseDocs/$STREAM_NAME
endif

if ( -e $CASEROOT/CaseDocs/$NML_NAME ) then
	chmod 644 $CASEROOT/CaseDocs/$NML_NAME
endif

cp $MPAS_STREAMS $CASEROOT/CaseDocs/.
cp $RUNDIR/$NML_NAME $CASEROOT/CaseDocs/.
