#!/bin/csh

# For now, manually build the namelist. Soon this will call the standard ACME
# build-namelist script.

mkdir -p $CASEROOT/CaseDocs

if !(-d $EXEROOT/glc/obj ) mkdir -p $EXEROOT/glc/obj || exit 2
if !(-d $EXEROOT/glc/source) mkdir -p $EXEROOT/glc/source || exit 3

if !(-d $CASEBUILD/mpas-liconf) mkdir -p $CASEBUILD/mpas-liconf || exit 1
cd $CASEBUILD/mpas-liconf || exit -1

set inst_string = ""
set STREAM_NAME = "streams.landice"
set NML_NAME = "mpasli_in"

if (-e $CASEROOT/user_nl_mpasli${inst_string}) then
	$UTILROOT/Tools/user_nlcreate                                           \
		-user_nl_file $CASEROOT/user_nl_mpasli${inst_string}                 \
		-namelist_name mpasli_inparm >! $CASEBUILD/mpas-liconf/cesm_namelist
endif

# Check to see if "-preview" flag should be passed
if ( $?PREVIEW_NML ) then
	set PREVIEW_FLAG = "-preview"
else
	set PREVIEW_FLAG = ""
endif

# Check to see if build-namelist exists in SourceMods
if (-e $CASEROOT/SourceMods/src.mpas-li/build-namelist) then
	set BLD_NML_DIR = $CASEROOT/SourceMods/src.mpas-li
	set CFG_FLAG = "-cfg_dir $CODEROOT/glc/mpas-li/bld"
else
	set BLD_NML_DIR = $CODEROOT/glc/mpas-li/bld
	set CFG_FLAG = ""
endif

# Define input_mesh file and graph prefix for each defined GLC/MPAS-LI mesh
if ( $GLC_GRID == 'mpas.gis20km' ) then
        set date_stamp = 150505
	set input_mesh = $DIN_LOC_ROOT/glc/mpas-li/$GLC_GRID/gis20km.${date_stamp}.nc
	set graph_prefix = $DIN_LOC_ROOT/glc/mpas-li/$GLC_GRID/mpas-li.graph.info.${date_stamp}
endif

# Write mpas-li.input_data_list file
echo "mesh = $input_mesh" > $CASEBUILD/mpas-li.input_data_list
#echo "graph1 = $graph_prefix" >> $CASEBUILD/mpas-li.input_data_list
echo "graph$NTASKS_GLC = $graph_prefix.part.$NTASKS_GLC" >> $CASEBUILD/mpas-li.input_data_list



# --------------------------------------
# Setup streams file
# --------------------------------------

set MPAS_STREAMS = $RUNDIR/$STREAM_NAME
touch $MPAS_STREAMS
chmod 644 $MPAS_STREAMS

# Write streams file, if there isn't one in SourceMods
if (-e $CASEROOT/SourceMods/src.mpas-li/$STREAM_NAME) then
	cp $CASEROOT/SourceMods/src.mpas-li/$STREAM_NAME $MPAS_STREAMS
else
	cat >! $MPAS_STREAMS << 'EOF'
	<streams>

	<immutable_stream name="basicmesh"
					  type="none"
					  filename_template="not-to-be-used.nc"
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
					  filename_template="rst.glc.$Y-$M-$D_$h.$m.$s.nc"
					  filename_interval="output_interval"
					  reference_time="0000-01-01_00:00:00"
					  clobber_mode="truncate"
					  input_interval="initial_only"
					  output_interval="10-00-00_00:00:00"/>

	<!--
	output is the main history output stream. You can add auxiliary streams to
	this stream to include more fields.
	-->

	<stream name="output"
			type="output"
			filename_template="hist.glc.$Y-$M-$D_$h.$m.$s.nc"
			filename_interval="01-00-00_00:00:00"
			reference_time="0000-01-01_00:00:00"
			clobber_mode="truncate"
			output_interval="00-01-00_00:00:00">

    <stream name="basicmesh"/>
    <var_array name="tracers"/>
    <var name="xtime"/>
    <var name="thickness"/>
    <var name="layerThickness"/>
    <var name="lowerSurface"/>
    <var name="upperSurface"/>
    <var name="cellMask"/>
    <var name="edgeMask"/>
    <var name="vertexMask"/>
    <var name="normalVelocity"/>
    <var name="uReconstructX"/>
    <var name="uReconstructY"/>

</stream>

</streams>

'EOF'
endif # Writing streams file

$BLD_NML_DIR/build-namelist $CFG_FLAG $PREVIEW_FLAG                        \
		-infile $CASEBUILD/mpas-liconf/cesm_namelist                        \
		-caseroot $CASEROOT                                                \
		-casebuild $CASEBUILD                                              \
		-scriptsroot $SCRIPTSROOT                                          \
		-inst_string "$inst_string"                                        \
                -date_stamp "$date_stamp"                                          \
		-glc_grid "$GLC_GRID" || exit -1

if ( -d ${RUNDIR} ) then
	cp $CASEBUILD/mpas-liconf/mpasli_in ${RUNDIR}/mpasli_in
endif

#cp $MPAS_STREAMS $RUNDIR

if ( -e $CASEROOT/CaseDocs/$STREAM_NAME ) then
	chmod 644 $CASEROOT/CaseDocs/$STREAM_NAME
endif

if ( -e $CASEROOT/CaseDocs/$NML_NAME ) then
	chmod 644 $CASEROOT/CaseDocs/$NML_NAME
endif

cp $MPAS_STREAMS $CASEROOT/CaseDocs/.
cp $RUNDIR/$NML_NAME $CASEROOT/CaseDocs/.

