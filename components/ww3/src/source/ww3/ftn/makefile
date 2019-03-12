 
# WAVEWATCH III executables
# -------------------------
# Grid preprocessor
 
$(aPe)/ww3_grid : $(aPo)/constants.o
$(aPe)/ww3_grid : $(aPo)/w3servmd.o
$(aPe)/ww3_grid : $(aPo)/w3arrymd.o
$(aPe)/ww3_grid : $(aPo)/w3dispmd.o
$(aPe)/ww3_grid : $(aPo)/w3gsrumd.o
$(aPe)/ww3_grid : $(aPo)/w3timemd.o
$(aPe)/ww3_grid : $(aPo)/w3wdatmd.o
$(aPe)/ww3_grid : $(aPo)/w3gdatmd.o
$(aPe)/ww3_grid : $(aPo)/w3adatmd.o
$(aPe)/ww3_grid : $(aPo)/w3idatmd.o
$(aPe)/ww3_grid : $(aPo)/w3odatmd.o
$(aPe)/ww3_grid : $(aPo)/w3triamd.o
$(aPe)/ww3_grid : $(aPo)/w3src4md.o
$(aPe)/ww3_grid : $(aPo)/w3snl1md.o
$(aPe)/ww3_grid : $(aPo)/w3iogrmd.o
$(aPe)/ww3_grid : $(aPo)/ww3_grid.o
	@$(aPb)/link ww3_grid w3wdatmd w3gdatmd w3adatmd w3idatmd w3odatmd   w3triamd w3src4md w3snl1md   w3iogrmd constants w3servmd w3arrymd w3dispmd w3gsrumd w3timemd
 
# Initial conditions program
 
$(aPe)/ww3_strt : $(aPo)/constants.o
$(aPe)/ww3_strt : $(aPo)/w3triamd.o
$(aPe)/ww3_strt : $(aPo)/w3servmd.o
$(aPe)/ww3_strt : $(aPo)/w3arrymd.o
$(aPe)/ww3_strt : $(aPo)/w3dispmd.o
$(aPe)/ww3_strt : $(aPo)/w3gsrumd.o
$(aPe)/ww3_strt : $(aPo)/w3timemd.o
$(aPe)/ww3_strt : $(aPo)/w3gdatmd.o
$(aPe)/ww3_strt : $(aPo)/w3wdatmd.o
$(aPe)/ww3_strt : $(aPo)/w3adatmd.o
$(aPe)/ww3_strt : $(aPo)/w3idatmd.o
$(aPe)/ww3_strt : $(aPo)/w3odatmd.o
$(aPe)/ww3_strt : $(aPo)/w3src4md.o
$(aPe)/ww3_strt : $(aPo)/w3snl1md.o
$(aPe)/ww3_strt : $(aPo)/w3iogrmd.o
$(aPe)/ww3_strt : $(aPo)/w3iorsmd.o
$(aPe)/ww3_strt : $(aPo)/ww3_strt.o
	@$(aPb)/link ww3_strt w3gdatmd w3wdatmd w3adatmd w3idatmd w3odatmd   w3src4md w3snl1md   w3iogrmd w3iorsmd constants w3triamd w3servmd w3arrymd w3dispmd w3gsrumd w3timemd
 
# Field preprocessor
 
$(aPe)/ww3_prep : $(aPo)/constants.o
$(aPe)/ww3_prep : $(aPo)/w3servmd.o
$(aPe)/ww3_prep : $(aPo)/w3timemd.o
$(aPe)/ww3_prep : $(aPo)/w3arrymd.o
$(aPe)/ww3_prep : $(aPo)/w3dispmd.o
$(aPe)/ww3_prep : $(aPo)/w3gsrumd.o
$(aPe)/ww3_prep : $(aPo)/w3fldsmd.o
$(aPe)/ww3_prep : $(aPo)/w3gdatmd.o
$(aPe)/ww3_prep : $(aPo)/w3adatmd.o
$(aPe)/ww3_prep : $(aPo)/w3idatmd.o
$(aPe)/ww3_prep : $(aPo)/w3odatmd.o
$(aPe)/ww3_prep : $(aPo)/w3wdatmd.o
$(aPe)/ww3_prep : $(aPo)/w3triamd.o
$(aPe)/ww3_prep : $(aPo)/w3src4md.o
$(aPe)/ww3_prep : $(aPo)/w3snl1md.o
$(aPe)/ww3_prep : $(aPo)/w3iogrmd.o
$(aPe)/ww3_prep : $(aPo)/ww3_prep.o
	@$(aPb)/link ww3_prep w3gdatmd w3adatmd w3idatmd w3odatmd w3wdatmd w3fldsmd  w3triamd w3src4md w3snl1md    w3iogrmd    constants w3servmd w3timemd  w3arrymd w3dispmd w3gsrumd
 
# NetCDF field preprocessor
 
$(aPe)/ww3_prnc : $(aPo)/constants.o
$(aPe)/ww3_prnc : $(aPo)/w3servmd.o
$(aPe)/ww3_prnc : $(aPo)/w3timemd.o
$(aPe)/ww3_prnc : $(aPo)/w3arrymd.o
$(aPe)/ww3_prnc : $(aPo)/w3dispmd.o
$(aPe)/ww3_prnc : $(aPo)/w3gsrumd.o
$(aPe)/ww3_prnc : $(aPo)/w3fldsmd.o
$(aPe)/ww3_prnc : $(aPo)/w3gdatmd.o
$(aPe)/ww3_prnc : $(aPo)/w3adatmd.o
$(aPe)/ww3_prnc : $(aPo)/w3idatmd.o
$(aPe)/ww3_prnc : $(aPo)/w3odatmd.o
$(aPe)/ww3_prnc : $(aPo)/w3wdatmd.o
$(aPe)/ww3_prnc : $(aPo)/w3triamd.o
$(aPe)/ww3_prnc : $(aPo)/w3src4md.o
$(aPe)/ww3_prnc : $(aPo)/w3snl1md.o
$(aPe)/ww3_prnc : $(aPo)/w3iogrmd.o
$(aPe)/ww3_prnc : $(aPo)/ww3_prnc.o
	@$(aPb)/link ww3_prnc w3gdatmd w3adatmd w3idatmd w3odatmd w3wdatmd w3fldsmd  w3triamd w3src4md w3snl1md   w3iogrmd    constants w3servmd w3timemd w3arrymd w3dispmd w3gsrumd 
 
# Generic shell
 
$(aPe)/ww3_shel : $(aPo)/constants.o
$(aPe)/ww3_shel : $(aPo)/w3servmd.o
$(aPe)/ww3_shel : $(aPo)/w3timemd.o
$(aPe)/ww3_shel : $(aPo)/w3arrymd.o
$(aPe)/ww3_shel : $(aPo)/w3dispmd.o
$(aPe)/ww3_shel : $(aPo)/w3cspcmd.o
$(aPe)/ww3_shel : $(aPo)/w3gsrumd.o
$(aPe)/ww3_shel : $(aPo)/w3namlmd.o
$(aPe)/ww3_shel : $(aPo)/w3fldsmd.o
$(aPe)/ww3_shel : $(aPo)/w3initmd.o
$(aPe)/ww3_shel : $(aPo)/w3wavemd.o
$(aPe)/ww3_shel : $(aPo)/w3wdasmd.o
$(aPe)/ww3_shel : $(aPo)/w3updtmd.o
$(aPe)/ww3_shel : $(aPo)/wmmdatmd.o
$(aPe)/ww3_shel : $(aPo)/w3gdatmd.o
$(aPe)/ww3_shel : $(aPo)/w3wdatmd.o
$(aPe)/ww3_shel : $(aPo)/w3adatmd.o
$(aPe)/ww3_shel : $(aPo)/w3idatmd.o
$(aPe)/ww3_shel : $(aPo)/w3odatmd.o
$(aPe)/ww3_shel : $(aPo)/w3profsmd.o
$(aPe)/ww3_shel : $(aPo)/w3pro3md.o
$(aPe)/ww3_shel : $(aPo)/w3uqckmd.o
$(aPe)/ww3_shel : $(aPo)/w3triamd.o
$(aPe)/ww3_shel : $(aPo)/w3srcemd.o
$(aPe)/ww3_shel : $(aPo)/w3flx1md.o
$(aPe)/ww3_shel : $(aPo)/w3sln1md.o
$(aPe)/ww3_shel : $(aPo)/w3src4md.o
$(aPe)/ww3_shel : $(aPo)/w3snl1md.o
$(aPe)/ww3_shel : $(aPo)/w3sbt1md.o
$(aPe)/ww3_shel : $(aPo)/w3sdb1md.o
$(aPe)/ww3_shel : $(aPo)/w3iogrmd.o
$(aPe)/ww3_shel : $(aPo)/w3iogomd.o
$(aPe)/ww3_shel : $(aPo)/w3iopomd.o
$(aPe)/ww3_shel : $(aPo)/w3iotrmd.o
$(aPe)/ww3_shel : $(aPo)/w3iorsmd.o
$(aPe)/ww3_shel : $(aPo)/w3iobcmd.o
$(aPe)/ww3_shel : $(aPo)/w3iosfmd.o
$(aPe)/ww3_shel : $(aPo)/w3partmd.o
$(aPe)/ww3_shel : $(aPo)/ww3_shel.o
	@$(aPb)/link ww3_shel wmmdatmd w3gdatmd w3wdatmd w3adatmd w3idatmd w3odatmd w3fldsmd w3initmd w3wavemd w3wdasmd w3updtmd w3profsmd w3pro3md w3uqckmd w3triamd w3srcemd  w3flx1md w3sln1md w3src4md w3snl1md w3sbt1md   w3sdb1md      w3iogrmd w3iogomd w3iopomd w3iotrmd w3iorsmd w3iobcmd    w3iosfmd w3partmd constants w3servmd w3timemd  w3arrymd w3dispmd w3cspcmd w3gsrumd  w3namlmd
 
# Multi-grid shell
 
$(aPe)/ww3_multi : $(aPo)/constants.o
$(aPe)/ww3_multi : $(aPo)/w3servmd.o
$(aPe)/ww3_multi : $(aPo)/w3timemd.o
$(aPe)/ww3_multi : $(aPo)/w3arrymd.o
$(aPe)/ww3_multi : $(aPo)/w3dispmd.o
$(aPe)/ww3_multi : $(aPo)/w3cspcmd.o
$(aPe)/ww3_multi : $(aPo)/w3gsrumd.o
$(aPe)/ww3_multi : $(aPo)/wmunitmd.o
$(aPe)/ww3_multi : $(aPo)/w3namlmd.o
$(aPe)/ww3_multi : $(aPo)/wminitmd.o
$(aPe)/ww3_multi : $(aPo)/wmwavemd.o
$(aPe)/ww3_multi : $(aPo)/wmfinlmd.o
$(aPe)/ww3_multi : $(aPo)/wmgridmd.o
$(aPe)/ww3_multi : $(aPo)/wmupdtmd.o
$(aPe)/ww3_multi : $(aPo)/wminiomd.o
$(aPe)/ww3_multi : $(aPo)/w3fldsmd.o
$(aPe)/ww3_multi : $(aPo)/w3initmd.o
$(aPe)/ww3_multi : $(aPo)/w3wavemd.o
$(aPe)/ww3_multi : $(aPo)/w3wdasmd.o
$(aPe)/ww3_multi : $(aPo)/w3updtmd.o
$(aPe)/ww3_multi : $(aPo)/wmmdatmd.o
$(aPe)/ww3_multi : $(aPo)/w3gdatmd.o
$(aPe)/ww3_multi : $(aPo)/w3wdatmd.o
$(aPe)/ww3_multi : $(aPo)/w3adatmd.o
$(aPe)/ww3_multi : $(aPo)/w3idatmd.o
$(aPe)/ww3_multi : $(aPo)/w3odatmd.o
$(aPe)/ww3_multi : $(aPo)/w3profsmd.o
$(aPe)/ww3_multi : $(aPo)/w3pro3md.o
$(aPe)/ww3_multi : $(aPo)/w3uqckmd.o
$(aPe)/ww3_multi : $(aPo)/w3triamd.o
$(aPe)/ww3_multi : $(aPo)/w3srcemd.o
$(aPe)/ww3_multi : $(aPo)/w3flx1md.o
$(aPe)/ww3_multi : $(aPo)/w3sln1md.o
$(aPe)/ww3_multi : $(aPo)/w3src4md.o
$(aPe)/ww3_multi : $(aPo)/w3snl1md.o
$(aPe)/ww3_multi : $(aPo)/w3sbt1md.o
$(aPe)/ww3_multi : $(aPo)/w3sdb1md.o
$(aPe)/ww3_multi : $(aPo)/w3iogrmd.o
$(aPe)/ww3_multi : $(aPo)/w3iogomd.o
$(aPe)/ww3_multi : $(aPo)/w3iopomd.o
$(aPe)/ww3_multi : $(aPo)/wmiopomd.o
$(aPe)/ww3_multi : $(aPo)/w3iotrmd.o
$(aPe)/ww3_multi : $(aPo)/w3iorsmd.o
$(aPe)/ww3_multi : $(aPo)/w3iobcmd.o
$(aPe)/ww3_multi : $(aPo)/w3iosfmd.o
$(aPe)/ww3_multi : $(aPo)/w3partmd.o
$(aPe)/ww3_multi : $(aPo)/ww3_multi.o
	@$(aPb)/link ww3_multi wmmdatmd w3gdatmd w3wdatmd w3adatmd w3idatmd w3odatmd wminitmd wmwavemd wmfinlmd wmgridmd wmupdtmd wminiomd w3fldsmd w3initmd w3wavemd w3wdasmd w3updtmd w3profsmd w3pro3md w3uqckmd w3triamd w3srcemd  w3flx1md w3sln1md w3src4md w3snl1md w3sbt1md   w3sdb1md      w3iogrmd w3iogomd w3iopomd wmiopomd w3iotrmd w3iorsmd w3iobcmd w3iosfmd w3partmd    constants  w3servmd w3timemd w3arrymd w3dispmd w3cspcmd w3gsrumd   wmunitmd w3namlmd
 
# Multi-grid shell sbs version
 
$(aPe)/ww3_sbs1 : $(aPo)/constants.o
$(aPe)/ww3_sbs1 : $(aPo)/w3servmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3timemd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3arrymd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3dispmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3cspcmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3gsrumd.o
$(aPe)/ww3_sbs1 : $(aPo)/wmunitmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3namlmd.o
$(aPe)/ww3_sbs1 : $(aPo)/wminitmd.o
$(aPe)/ww3_sbs1 : $(aPo)/wmwavemd.o
$(aPe)/ww3_sbs1 : $(aPo)/wmfinlmd.o
$(aPe)/ww3_sbs1 : $(aPo)/wmgridmd.o
$(aPe)/ww3_sbs1 : $(aPo)/wmupdtmd.o
$(aPe)/ww3_sbs1 : $(aPo)/wminiomd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3fldsmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3initmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3wavemd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3wdasmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3updtmd.o
$(aPe)/ww3_sbs1 : $(aPo)/wmmdatmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3gdatmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3wdatmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3adatmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3idatmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3odatmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3profsmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3pro3md.o
$(aPe)/ww3_sbs1 : $(aPo)/w3uqckmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3triamd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3srcemd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3flx1md.o
$(aPe)/ww3_sbs1 : $(aPo)/w3sln1md.o
$(aPe)/ww3_sbs1 : $(aPo)/w3src4md.o
$(aPe)/ww3_sbs1 : $(aPo)/w3snl1md.o
$(aPe)/ww3_sbs1 : $(aPo)/w3sbt1md.o
$(aPe)/ww3_sbs1 : $(aPo)/w3sdb1md.o
$(aPe)/ww3_sbs1 : $(aPo)/w3iogrmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3iogomd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3iopomd.o
$(aPe)/ww3_sbs1 : $(aPo)/wmiopomd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3iotrmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3iorsmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3iobcmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3iosfmd.o
$(aPe)/ww3_sbs1 : $(aPo)/w3partmd.o
$(aPe)/ww3_sbs1 : $(aPo)/ww3_sbs1.o
	@$(aPb)/link ww3_sbs1 wmmdatmd w3gdatmd w3wdatmd w3adatmd w3idatmd w3odatmd wminitmd wmwavemd wmfinlmd wmgridmd wmupdtmd wminiomd w3fldsmd w3initmd w3wavemd w3wdasmd w3updtmd w3profsmd w3pro3md w3uqckmd w3triamd w3srcemd  w3flx1md w3sln1md w3src4md w3snl1md w3sbt1md w3sdb1md        w3iogrmd w3iogomd w3iopomd wmiopomd w3iotrmd w3iorsmd w3iobcmd w3iosfmd w3partmd    constants w3servmd w3timemd w3arrymd w3dispmd w3cspcmd w3gsrumd    wmunitmd w3namlmd
 
# Gridded output
 
$(aPe)/ww3_outf : $(aPo)/constants.o
$(aPe)/ww3_outf : $(aPo)/w3servmd.o
$(aPe)/ww3_outf : $(aPo)/w3timemd.o
$(aPe)/ww3_outf : $(aPo)/w3arrymd.o
$(aPe)/ww3_outf : $(aPo)/w3dispmd.o
$(aPe)/ww3_outf : $(aPo)/w3gsrumd.o
$(aPe)/ww3_outf : $(aPo)/w3namlmd.o
$(aPe)/ww3_outf : $(aPo)/wmmdatmd.o
$(aPe)/ww3_outf : $(aPo)/w3gdatmd.o
$(aPe)/ww3_outf : $(aPo)/w3wdatmd.o
$(aPe)/ww3_outf : $(aPo)/w3adatmd.o
$(aPe)/ww3_outf : $(aPo)/w3idatmd.o
$(aPe)/ww3_outf : $(aPo)/w3odatmd.o
$(aPe)/ww3_outf : $(aPo)/w3src4md.o
$(aPe)/ww3_outf : $(aPo)/w3snl1md.o
$(aPe)/ww3_outf : $(aPo)/w3iogrmd.o
$(aPe)/ww3_outf : $(aPo)/w3iogomd.o
$(aPe)/ww3_outf : $(aPo)/ww3_outf.o
	@$(aPb)/link ww3_outf wmmdatmd w3gdatmd w3wdatmd w3adatmd w3idatmd w3odatmd   w3src4md w3snl1md   w3iogrmd w3iogomd constants w3servmd w3timemd w3arrymd w3dispmd w3gsrumd w3namlmd
 
# Point output
 
$(aPe)/ww3_outp : $(aPo)/constants.o
$(aPe)/ww3_outp : $(aPo)/w3servmd.o
$(aPe)/ww3_outp : $(aPo)/w3timemd.o
$(aPe)/ww3_outp : $(aPo)/w3arrymd.o
$(aPe)/ww3_outp : $(aPo)/w3dispmd.o
$(aPe)/ww3_outp : $(aPo)/w3gsrumd.o
$(aPe)/ww3_outp : $(aPo)/w3triamd.o
$(aPe)/ww3_outp : $(aPo)/w3gdatmd.o
$(aPe)/ww3_outp : $(aPo)/w3wdatmd.o
$(aPe)/ww3_outp : $(aPo)/w3adatmd.o
$(aPe)/ww3_outp : $(aPo)/w3idatmd.o
$(aPe)/ww3_outp : $(aPo)/w3odatmd.o
$(aPe)/ww3_outp : $(aPo)/w3flx1md.o
$(aPe)/ww3_outp : $(aPo)/w3sln1md.o
$(aPe)/ww3_outp : $(aPo)/w3src4md.o
$(aPe)/ww3_outp : $(aPo)/w3snl1md.o
$(aPe)/ww3_outp : $(aPo)/w3sbt1md.o
$(aPe)/ww3_outp : $(aPo)/w3sdb1md.o
$(aPe)/ww3_outp : $(aPo)/w3bullmd.o
$(aPe)/ww3_outp : $(aPo)/w3iogrmd.o
$(aPe)/ww3_outp : $(aPo)/w3iopomd.o
$(aPe)/ww3_outp : $(aPo)/w3partmd.o
$(aPe)/ww3_outp : $(aPo)/ww3_outp.o
	@$(aPb)/link ww3_outp w3triamd w3gdatmd w3wdatmd w3adatmd w3idatmd w3odatmd   w3flx1md w3sln1md w3src4md w3snl1md w3sbt1md   w3sdb1md     w3bullmd w3iogrmd w3iopomd w3partmd constants w3servmd w3timemd w3arrymd w3dispmd w3gsrumd
 
# Track output post
 
$(aPe)/ww3_trck : $(aPo)/constants.o
$(aPe)/ww3_trck : $(aPo)/w3servmd.o
$(aPe)/ww3_trck : $(aPo)/w3timemd.o
$(aPe)/ww3_trck : $(aPo)/w3gsrumd.o
$(aPe)/ww3_trck : $(aPo)/w3gdatmd.o
$(aPe)/ww3_trck : $(aPo)/w3odatmd.o
$(aPe)/ww3_trck : $(aPo)/ww3_trck.o
	@$(aPb)/link ww3_trck w3gdatmd w3odatmd     constants w3servmd w3timemd w3gsrumd
 
# Track NetCDF output post
 
$(aPe)/ww3_trnc : $(aPo)/constants.o
$(aPe)/ww3_trnc : $(aPo)/w3servmd.o
$(aPe)/ww3_trnc : $(aPo)/w3timemd.o
$(aPe)/ww3_trnc : $(aPo)/w3gsrumd.o
$(aPe)/ww3_trnc : $(aPo)/w3gdatmd.o
$(aPe)/ww3_trnc : $(aPo)/w3odatmd.o
$(aPe)/ww3_trnc : $(aPo)/ww3_trnc.o
	@$(aPb)/link ww3_trnc w3gdatmd w3odatmd     constants w3servmd w3timemd w3gsrumd
 
# Gridded output (GRIB)
 
$(aPe)/ww3_grib : $(aPo)/constants.o
$(aPe)/ww3_grib : $(aPo)/w3servmd.o
$(aPe)/ww3_grib : $(aPo)/w3timemd.o
$(aPe)/ww3_grib : $(aPo)/w3arrymd.o
$(aPe)/ww3_grib : $(aPo)/w3dispmd.o
$(aPe)/ww3_grib : $(aPo)/w3gsrumd.o
$(aPe)/ww3_grib : $(aPo)/w3namlmd.o
$(aPe)/ww3_grib : $(aPo)/wmmdatmd.o
$(aPe)/ww3_grib : $(aPo)/w3triamd.o
$(aPe)/ww3_grib : $(aPo)/w3gdatmd.o
$(aPe)/ww3_grib : $(aPo)/w3wdatmd.o
$(aPe)/ww3_grib : $(aPo)/w3adatmd.o
$(aPe)/ww3_grib : $(aPo)/w3idatmd.o
$(aPe)/ww3_grib : $(aPo)/w3odatmd.o
$(aPe)/ww3_grib : $(aPo)/w3src4md.o
$(aPe)/ww3_grib : $(aPo)/w3snl1md.o
$(aPe)/ww3_grib : $(aPo)/w3iogrmd.o
$(aPe)/ww3_grib : $(aPo)/w3iogomd.o
$(aPe)/ww3_grib : $(aPo)/ww3_grib.o
	@$(aPb)/link ww3_grib wmmdatmd w3triamd w3gdatmd w3wdatmd w3adatmd w3idatmd w3odatmd   w3src4md w3snl1md   w3iogrmd w3iogomd constants w3servmd w3timemd w3arrymd w3dispmd w3gsrumd w3namlmd
 
# GrADS input file generation (gridded fields)
 
$(aPe)/gx_outf : $(aPo)/constants.o
$(aPe)/gx_outf : $(aPo)/w3servmd.o
$(aPe)/gx_outf : $(aPo)/w3timemd.o
$(aPe)/gx_outf : $(aPo)/w3arrymd.o
$(aPe)/gx_outf : $(aPo)/w3dispmd.o
$(aPe)/gx_outf : $(aPo)/w3gsrumd.o
$(aPe)/gx_outf : $(aPo)/w3namlmd.o
$(aPe)/gx_outf : $(aPo)/wmmdatmd.o
$(aPe)/gx_outf : $(aPo)/w3gdatmd.o
$(aPe)/gx_outf : $(aPo)/w3wdatmd.o
$(aPe)/gx_outf : $(aPo)/w3adatmd.o
$(aPe)/gx_outf : $(aPo)/w3idatmd.o
$(aPe)/gx_outf : $(aPo)/w3odatmd.o
$(aPe)/gx_outf : $(aPo)/w3triamd.o
$(aPe)/gx_outf : $(aPo)/w3src4md.o
$(aPe)/gx_outf : $(aPo)/w3snl1md.o
$(aPe)/gx_outf : $(aPo)/w3sdb1md.o
$(aPe)/gx_outf : $(aPo)/w3iogrmd.o
$(aPe)/gx_outf : $(aPo)/w3iogomd.o
$(aPe)/gx_outf : $(aPo)/gx_outf.o
	@$(aPb)/link gx_outf wmmdatmd w3gdatmd w3wdatmd w3adatmd w3idatmd w3odatmd   w3triamd w3src4md w3snl1md  w3sdb1md     w3iogrmd w3iogomd constants w3servmd w3timemd w3arrymd w3dispmd w3gsrumd w3namlmd
 
# GrADS input file generation for point output
 
$(aPe)/gx_outp : $(aPo)/constants.o
$(aPe)/gx_outp : $(aPo)/w3servmd.o
$(aPe)/gx_outp : $(aPo)/w3timemd.o
$(aPe)/gx_outp : $(aPo)/w3arrymd.o
$(aPe)/gx_outp : $(aPo)/w3dispmd.o
$(aPe)/gx_outp : $(aPo)/w3gsrumd.o
$(aPe)/gx_outp : $(aPo)/w3triamd.o
$(aPe)/gx_outp : $(aPo)/w3gdatmd.o
$(aPe)/gx_outp : $(aPo)/w3wdatmd.o
$(aPe)/gx_outp : $(aPo)/w3adatmd.o
$(aPe)/gx_outp : $(aPo)/w3idatmd.o
$(aPe)/gx_outp : $(aPo)/w3odatmd.o
$(aPe)/gx_outp : $(aPo)/w3flx1md.o
$(aPe)/gx_outp : $(aPo)/w3sln1md.o
$(aPe)/gx_outp : $(aPo)/w3src4md.o
$(aPe)/gx_outp : $(aPo)/w3snl1md.o
$(aPe)/gx_outp : $(aPo)/w3sbt1md.o
$(aPe)/gx_outp : $(aPo)/w3sdb1md.o
$(aPe)/gx_outp : $(aPo)/w3iogrmd.o
$(aPe)/gx_outp : $(aPo)/w3iopomd.o
$(aPe)/gx_outp : $(aPo)/gx_outp.o
	@$(aPb)/link gx_outp w3triamd w3gdatmd w3wdatmd w3adatmd w3idatmd w3odatmd   w3flx1md w3sln1md w3src4md w3snl1md w3sbt1md   w3sdb1md    w3iogrmd w3iopomd constants w3servmd w3timemd w3arrymd w3dispmd w3gsrumd
 
# Gridded NetCDF output
 
$(aPe)/ww3_ounf : $(aPo)/constants.o
$(aPe)/ww3_ounf : $(aPo)/w3servmd.o
$(aPe)/ww3_ounf : $(aPo)/w3timemd.o
$(aPe)/ww3_ounf : $(aPo)/w3arrymd.o
$(aPe)/ww3_ounf : $(aPo)/w3dispmd.o
$(aPe)/ww3_ounf : $(aPo)/w3gsrumd.o
$(aPe)/ww3_ounf : $(aPo)/w3namlmd.o
$(aPe)/ww3_ounf : $(aPo)/w3initmd.o
$(aPe)/ww3_ounf : $(aPo)/wmmdatmd.o
$(aPe)/ww3_ounf : $(aPo)/w3gdatmd.o
$(aPe)/ww3_ounf : $(aPo)/w3wdatmd.o
$(aPe)/ww3_ounf : $(aPo)/w3adatmd.o
$(aPe)/ww3_ounf : $(aPo)/w3idatmd.o
$(aPe)/ww3_ounf : $(aPo)/w3odatmd.o
$(aPe)/ww3_ounf : $(aPo)/w3triamd.o
$(aPe)/ww3_ounf : $(aPo)/w3src4md.o
$(aPe)/ww3_ounf : $(aPo)/w3snl1md.o
$(aPe)/ww3_ounf : $(aPo)/w3iogrmd.o
$(aPe)/ww3_ounf : $(aPo)/w3iogomd.o
$(aPe)/ww3_ounf : $(aPo)/w3iorsmd.o
$(aPe)/ww3_ounf : $(aPo)/w3iopomd.o
$(aPe)/ww3_ounf : $(aPo)/ww3_ounf.o
	@$(aPb)/link ww3_ounf wmmdatmd w3gdatmd w3wdatmd w3adatmd w3idatmd w3odatmd w3initmd  w3triamd w3src4md w3snl1md    w3iogrmd w3iogomd w3iorsmd w3iopomd constants w3servmd w3timemd w3arrymd w3dispmd w3gsrumd w3namlmd
 
# Point NetCDF output
 
$(aPe)/ww3_ounp : $(aPo)/constants.o
$(aPe)/ww3_ounp : $(aPo)/w3servmd.o
$(aPe)/ww3_ounp : $(aPo)/w3timemd.o
$(aPe)/ww3_ounp : $(aPo)/w3arrymd.o
$(aPe)/ww3_ounp : $(aPo)/w3dispmd.o
$(aPe)/ww3_ounp : $(aPo)/w3gsrumd.o
$(aPe)/ww3_ounp : $(aPo)/w3triamd.o
$(aPe)/ww3_ounp : $(aPo)/w3gdatmd.o
$(aPe)/ww3_ounp : $(aPo)/w3wdatmd.o
$(aPe)/ww3_ounp : $(aPo)/w3adatmd.o
$(aPe)/ww3_ounp : $(aPo)/w3idatmd.o
$(aPe)/ww3_ounp : $(aPo)/w3odatmd.o
$(aPe)/ww3_ounp : $(aPo)/w3flx1md.o
$(aPe)/ww3_ounp : $(aPo)/w3sln1md.o
$(aPe)/ww3_ounp : $(aPo)/w3src4md.o
$(aPe)/ww3_ounp : $(aPo)/w3snl1md.o
$(aPe)/ww3_ounp : $(aPo)/w3sbt1md.o
$(aPe)/ww3_ounp : $(aPo)/w3sdb1md.o
$(aPe)/ww3_ounp : $(aPo)/w3bullmd.o
$(aPe)/ww3_ounp : $(aPo)/w3iogrmd.o
$(aPe)/ww3_ounp : $(aPo)/w3iopomd.o
$(aPe)/ww3_ounp : $(aPo)/w3partmd.o
$(aPe)/ww3_ounp : $(aPo)/ww3_ounp.o
	@$(aPb)/link ww3_ounp w3triamd w3gdatmd w3wdatmd w3adatmd w3idatmd w3odatmd   w3flx1md w3sln1md w3src4md w3snl1md w3sbt1md   w3sdb1md     w3bullmd w3iogrmd w3iopomd w3partmd constants w3servmd w3timemd w3arrymd w3dispmd w3gsrumd
 
# Grid splitting
 
$(aPe)/ww3_gspl : $(aPo)/constants.o
$(aPe)/ww3_gspl : $(aPo)/w3servmd.o
$(aPe)/ww3_gspl : $(aPo)/w3timemd.o
$(aPe)/ww3_gspl : $(aPo)/w3arrymd.o
$(aPe)/ww3_gspl : $(aPo)/w3dispmd.o
$(aPe)/ww3_gspl : $(aPo)/w3gsrumd.o
$(aPe)/ww3_gspl : $(aPo)/w3fldsmd.o
$(aPe)/ww3_gspl : $(aPo)/w3gdatmd.o
$(aPe)/ww3_gspl : $(aPo)/w3wdatmd.o
$(aPe)/ww3_gspl : $(aPo)/w3adatmd.o
$(aPe)/ww3_gspl : $(aPo)/w3idatmd.o
$(aPe)/ww3_gspl : $(aPo)/w3odatmd.o
$(aPe)/ww3_gspl : $(aPo)/w3triamd.o
$(aPe)/ww3_gspl : $(aPo)/w3src4md.o
$(aPe)/ww3_gspl : $(aPo)/w3snl1md.o
$(aPe)/ww3_gspl : $(aPo)/w3iogrmd.o
$(aPe)/ww3_gspl : $(aPo)/ww3_gspl.o
	@$(aPb)/link ww3_gspl w3gdatmd w3wdatmd w3adatmd w3idatmd w3odatmd w3fldsmd  w3triamd w3src4md w3snl1md   w3iogrmd     constants w3servmd w3timemd w3arrymd w3dispmd w3gsrumd 
 
# Grid Interpolation
 
$(aPe)/ww3_gint : $(aPo)/constants.o
$(aPe)/ww3_gint : $(aPo)/w3triamd.o
$(aPe)/ww3_gint : $(aPo)/w3servmd.o
$(aPe)/ww3_gint : $(aPo)/w3arrymd.o
$(aPe)/ww3_gint : $(aPo)/w3dispmd.o
$(aPe)/ww3_gint : $(aPo)/w3timemd.o
$(aPe)/ww3_gint : $(aPo)/w3gsrumd.o
$(aPe)/ww3_gint : $(aPo)/w3namlmd.o
$(aPe)/ww3_gint : $(aPo)/wmmdatmd.o
$(aPe)/ww3_gint : $(aPo)/w3gdatmd.o
$(aPe)/ww3_gint : $(aPo)/w3wdatmd.o
$(aPe)/ww3_gint : $(aPo)/w3adatmd.o
$(aPe)/ww3_gint : $(aPo)/w3idatmd.o
$(aPe)/ww3_gint : $(aPo)/w3odatmd.o
$(aPe)/ww3_gint : $(aPo)/w3src4md.o
$(aPe)/ww3_gint : $(aPo)/w3snl1md.o
$(aPe)/ww3_gint : $(aPo)/w3iogrmd.o
$(aPe)/ww3_gint : $(aPo)/w3iogomd.o
$(aPe)/ww3_gint : $(aPo)/ww3_gint.o
	@$(aPb)/link ww3_gint wmmdatmd w3gdatmd w3wdatmd w3adatmd w3idatmd w3odatmd   w3src4md w3snl1md  w3iogrmd w3iogomd constants w3triamd w3servmd  w3arrymd w3dispmd w3timemd w3gsrumd w3namlmd
 
# boundary conditions program
 
$(aPe)/ww3_bound : $(aPo)/constants.o
$(aPe)/ww3_bound : $(aPo)/w3servmd.o
$(aPe)/ww3_bound : $(aPo)/w3timemd.o
$(aPe)/ww3_bound : $(aPo)/w3cspcmd.o
$(aPe)/ww3_bound : $(aPo)/w3adatmd.o
$(aPe)/ww3_bound : $(aPo)/w3gdatmd.o
$(aPe)/ww3_bound : $(aPo)/w3wdatmd.o
$(aPe)/ww3_bound : $(aPo)/w3idatmd.o
$(aPe)/ww3_bound : $(aPo)/w3odatmd.o
$(aPe)/ww3_bound : $(aPo)/w3src4md.o
$(aPe)/ww3_bound : $(aPo)/w3snl1md.o
$(aPe)/ww3_bound : $(aPo)/w3triamd.o
$(aPe)/ww3_bound : $(aPo)/w3iobcmd.o
$(aPe)/ww3_bound : $(aPo)/w3iogrmd.o
$(aPe)/ww3_bound : $(aPo)/w3dispmd.o
$(aPe)/ww3_bound : $(aPo)/w3gsrumd.o
$(aPe)/ww3_bound : $(aPo)/ww3_bound.o
	@$(aPb)/link ww3_bound w3adatmd w3gdatmd w3wdatmd w3idatmd w3odatmd   w3src4md w3snl1md    w3triamd w3iobcmd w3iogrmd w3dispmd w3gsrumd constants w3servmd w3timemd w3cspcmd
 
# NetCDF boundary conditions program
 
$(aPe)/ww3_bounc : $(aPo)/constants.o
$(aPe)/ww3_bounc : $(aPo)/w3servmd.o
$(aPe)/ww3_bounc : $(aPo)/w3timemd.o
$(aPe)/ww3_bounc : $(aPo)/w3cspcmd.o
$(aPe)/ww3_bounc : $(aPo)/w3adatmd.o
$(aPe)/ww3_bounc : $(aPo)/w3gdatmd.o
$(aPe)/ww3_bounc : $(aPo)/w3wdatmd.o
$(aPe)/ww3_bounc : $(aPo)/w3idatmd.o
$(aPe)/ww3_bounc : $(aPo)/w3odatmd.o
$(aPe)/ww3_bounc : $(aPo)/w3src4md.o
$(aPe)/ww3_bounc : $(aPo)/w3snl1md.o
$(aPe)/ww3_bounc : $(aPo)/w3triamd.o
$(aPe)/ww3_bounc : $(aPo)/w3iobcmd.o
$(aPe)/ww3_bounc : $(aPo)/w3iogrmd.o
$(aPe)/ww3_bounc : $(aPo)/w3dispmd.o
$(aPe)/ww3_bounc : $(aPo)/w3gsrumd.o
$(aPe)/ww3_bounc : $(aPo)/ww3_bounc.o
	@$(aPb)/link ww3_bounc w3adatmd w3gdatmd w3wdatmd w3idatmd w3odatmd   w3src4md w3snl1md    w3triamd w3iobcmd w3iogrmd w3dispmd w3gsrumd constants w3servmd w3timemd w3cspcmd
 
# Wave system tracking postprocessor
 
$(aPe)/ww3_systrk : $(aPo)/constants.o
$(aPe)/ww3_systrk : $(aPo)/w3servmd.o
$(aPe)/ww3_systrk : $(aPo)/w3timemd.o
$(aPe)/ww3_systrk : $(aPo)/w3strkmd.o
$(aPe)/ww3_systrk : $(aPo)/ww3_systrk.o
	@$(aPb)/link ww3_systrk  w3strkmd    constants w3servmd w3timemd
 
 
# WAVEWATCH III subroutines
# -------------------------
 
$(aPo)/constants.o : constants.ftn w3macros.h 
	@$(aPb)/ad3 constants
 
$(aPo)/gx_outf.o : gx_outf.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3iogrmd.o $(aPo)/w3iogomd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o
	@$(aPb)/ad3 gx_outf
 
$(aPo)/gx_outp.o : gx_outp.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3iogrmd.o $(aPo)/w3iopomd.o $(aPo)/w3sln1md.o $(aPo)/w3src4md.o $(aPo)/w3snl1md.o $(aPo)/w3sbt1md.o $(aPo)/w3sdb1md.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o $(aPo)/w3dispmd.o
	@$(aPb)/ad3 gx_outp
 
$(aPo)/w3adatmd.o : w3adatmd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3odatmd.o $(aPo)/w3idatmd.o $(aPo)/w3servmd.o
	@$(aPb)/ad3 w3adatmd
 
$(aPo)/w3arrymd.o : w3arrymd.ftn w3macros.h 
	@$(aPb)/ad3 w3arrymd
 
$(aPo)/w3bullmd.o : w3bullmd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3odatmd.o $(aPo)/constants.o $(aPo)/w3timemd.o
	@$(aPb)/ad3 w3bullmd
 
$(aPo)/w3cspcmd.o : w3cspcmd.ftn w3macros.h  $(aPo)/constants.o $(aPo)/w3servmd.o
	@$(aPb)/ad3 w3cspcmd
 
$(aPo)/w3dispmd.o : w3dispmd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3odatmd.o $(aPo)/constants.o $(aPo)/w3servmd.o
	@$(aPb)/ad3 w3dispmd
 
$(aPo)/w3fldsmd.o : w3fldsmd.ftn w3macros.h  $(aPo)/w3idatmd.o $(aPo)/w3timemd.o $(aPo)/w3gsrumd.o
	@$(aPb)/ad3 w3fldsmd
 
$(aPo)/w3flx1md.o : w3flx1md.ftn w3macros.h  $(aPo)/w3odatmd.o $(aPo)/w3servmd.o
	@$(aPb)/ad3 w3flx1md
 
$(aPo)/w3gdatmd.o : w3gdatmd.ftn w3macros.h  $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3gsrumd.o
	@$(aPb)/ad3 w3gdatmd
 
$(aPo)/w3gsrumd.o : w3gsrumd.ftn w3macros.h  $(aPo)/w3servmd.o
	@$(aPb)/ad3 w3gsrumd
 
$(aPo)/w3idatmd.o : w3idatmd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3servmd.o
	@$(aPb)/ad3 w3idatmd
 
$(aPo)/w3initmd.o : w3initmd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3idatmd.o $(aPo)/w3iogrmd.o $(aPo)/w3iogomd.o $(aPo)/w3iopomd.o $(aPo)/w3iorsmd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o $(aPo)/w3arrymd.o $(aPo)/w3dispmd.o
	@$(aPb)/ad3 w3initmd
 
$(aPo)/w3iobcmd.o : w3iobcmd.ftn w3macros.h  $(aPo)/w3cspcmd.o $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3gsrumd.o $(aPo)/w3triamd.o
	@$(aPb)/ad3 w3iobcmd
 
$(aPo)/w3iogomd.o : w3iogomd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3namlmd.o
	@$(aPb)/ad3 w3iogomd
 
$(aPo)/w3iogrmd.o : w3iogrmd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3src4md.o $(aPo)/w3snl1md.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3dispmd.o
	@$(aPb)/ad3 w3iogrmd
 
$(aPo)/w3iopomd.o : w3iopomd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3gsrumd.o $(aPo)/w3triamd.o
	@$(aPb)/ad3 w3iopomd
 
$(aPo)/w3iorsmd.o : w3iorsmd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3servmd.o $(aPo)/w3triamd.o
	@$(aPb)/ad3 w3iorsmd
 
$(aPo)/w3iosfmd.o : w3iosfmd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3partmd.o $(aPo)/constants.o $(aPo)/w3servmd.o
	@$(aPb)/ad3 w3iosfmd
 
$(aPo)/w3iotrmd.o : w3iotrmd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/constants.o $(aPo)/w3timemd.o $(aPo)/w3gsrumd.o
	@$(aPb)/ad3 w3iotrmd
 
$(aPo)/w3namlmd.o : w3namlmd.ftn w3macros.h  $(aPo)/wmmdatmd.o
	@$(aPb)/ad3 w3namlmd
 
$(aPo)/w3odatmd.o : w3odatmd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3servmd.o
	@$(aPb)/ad3 w3odatmd
 
$(aPo)/w3partmd.o : w3partmd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3odatmd.o $(aPo)/constants.o $(aPo)/w3dispmd.o
	@$(aPb)/ad3 w3partmd
 
$(aPo)/w3pro3md.o : w3pro3md.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3idatmd.o $(aPo)/w3uqckmd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o
	@$(aPb)/ad3 w3pro3md
 
$(aPo)/w3profsmd.o : w3profsmd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3idatmd.o $(aPo)/constants.o $(aPo)/w3timemd.o
	@$(aPb)/ad3 w3profsmd
 
$(aPo)/w3sbt1md.o : w3sbt1md.ftn w3macros.h  $(aPo)/w3gdatmd.o
	@$(aPb)/ad3 w3sbt1md
 
$(aPo)/w3sdb1md.o : w3sdb1md.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3odatmd.o $(aPo)/constants.o
	@$(aPb)/ad3 w3sdb1md
 
$(aPo)/w3servmd.o : w3servmd.ftn w3macros.h  $(aPo)/constants.o
	@$(aPb)/ad3 w3servmd
 
$(aPo)/w3sln1md.o : w3sln1md.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3odatmd.o $(aPo)/constants.o $(aPo)/w3servmd.o
	@$(aPb)/ad3 w3sln1md
 
$(aPo)/w3snl1md.o : w3snl1md.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/constants.o
	@$(aPb)/ad3 w3snl1md
 
$(aPo)/w3src4md.o : w3src4md.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3odatmd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3dispmd.o
	@$(aPb)/ad3 w3src4md
 
$(aPo)/w3srcemd.o : w3srcemd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3odatmd.o $(aPo)/w3idatmd.o $(aPo)/w3sln1md.o $(aPo)/w3src4md.o $(aPo)/w3snl1md.o $(aPo)/w3sbt1md.o $(aPo)/w3sdb1md.o $(aPo)/constants.o $(aPo)/w3dispmd.o
	@$(aPb)/ad3 w3srcemd
 
$(aPo)/w3strkmd.o : w3strkmd.ftn w3macros.h 
	@$(aPb)/ad3 w3strkmd
 
$(aPo)/w3timemd.o : w3timemd.ftn w3macros.h 
	@$(aPb)/ad3 w3timemd
 
$(aPo)/w3triamd.o : w3triamd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3idatmd.o $(aPo)/w3iogrmd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3arrymd.o $(aPo)/w3dispmd.o
	@$(aPb)/ad3 w3triamd
 
$(aPo)/w3updtmd.o : w3updtmd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3idatmd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o $(aPo)/w3dispmd.o $(aPo)/w3triamd.o
	@$(aPb)/ad3 w3updtmd
 
$(aPo)/w3uqckmd.o : w3uqckmd.ftn w3macros.h 
	@$(aPb)/ad3 w3uqckmd
 
$(aPo)/w3wavemd.o : w3wavemd.ftn w3macros.h  $(aPo)/w3updtmd.o $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3idatmd.o $(aPo)/w3iogrmd.o $(aPo)/w3iogomd.o $(aPo)/w3iopomd.o $(aPo)/w3iotrmd.o $(aPo)/w3iorsmd.o $(aPo)/w3iobcmd.o $(aPo)/w3iosfmd.o $(aPo)/w3pro3md.o $(aPo)/w3profsmd.o $(aPo)/w3srcemd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o $(aPo)/w3triamd.o
	@$(aPb)/ad3 w3wavemd
 
$(aPo)/w3wdasmd.o : w3wdasmd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o
	@$(aPb)/ad3 w3wdasmd
 
$(aPo)/w3wdatmd.o : w3wdatmd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3odatmd.o $(aPo)/w3servmd.o
	@$(aPb)/ad3 w3wdatmd
 
$(aPo)/wmfinlmd.o : wmfinlmd.ftn w3macros.h  $(aPo)/w3timemd.o $(aPo)/wmmdatmd.o
	@$(aPb)/ad3 wmfinlmd
 
$(aPo)/wmgridmd.o : wmgridmd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3gsrumd.o $(aPo)/w3triamd.o $(aPo)/wmmdatmd.o
	@$(aPb)/ad3 wmgridmd
 
$(aPo)/wminiomd.o : wminiomd.ftn w3macros.h  $(aPo)/w3updtmd.o $(aPo)/w3cspcmd.o $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3iobcmd.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o $(aPo)/wmmdatmd.o
	@$(aPb)/ad3 wminiomd
 
$(aPo)/wminitmd.o : wminitmd.ftn w3macros.h  $(aPo)/w3initmd.o $(aPo)/w3fldsmd.o $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3idatmd.o $(aPo)/w3iogrmd.o $(aPo)/w3iogomd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o $(aPo)/w3gsrumd.o $(aPo)/wmmdatmd.o $(aPo)/wmgridmd.o $(aPo)/wmunitmd.o $(aPo)/wminiomd.o $(aPo)/wmiopomd.o $(aPo)/w3namlmd.o
	@$(aPb)/ad3 wminitmd
 
$(aPo)/wmiopomd.o : wmiopomd.ftn w3macros.h  $(aPo)/w3initmd.o $(aPo)/w3cspcmd.o $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3iopomd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3gsrumd.o $(aPo)/w3triamd.o $(aPo)/wmmdatmd.o
	@$(aPb)/ad3 wmiopomd
 
$(aPo)/wmmdatmd.o : wmmdatmd.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3odatmd.o $(aPo)/w3servmd.o
	@$(aPb)/ad3 wmmdatmd
 
$(aPo)/wmunitmd.o : wmunitmd.ftn w3macros.h  $(aPo)/w3servmd.o
	@$(aPb)/ad3 wmunitmd
 
$(aPo)/wmupdtmd.o : wmupdtmd.ftn w3macros.h  $(aPo)/w3fldsmd.o $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3idatmd.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o $(aPo)/wmmdatmd.o
	@$(aPb)/ad3 wmupdtmd
 
$(aPo)/wmwavemd.o : wmwavemd.ftn w3macros.h  $(aPo)/w3wavemd.o $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3idatmd.o $(aPo)/w3iopomd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o $(aPo)/wmmdatmd.o $(aPo)/wmupdtmd.o $(aPo)/wminiomd.o $(aPo)/wmiopomd.o
	@$(aPb)/ad3 wmwavemd
 
$(aPo)/ww3_bounc.o : ww3_bounc.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3iogrmd.o $(aPo)/w3iobcmd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o
	@$(aPb)/ad3 ww3_bounc
 
$(aPo)/ww3_bound.o : ww3_bound.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3iogrmd.o $(aPo)/w3iobcmd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o
	@$(aPb)/ad3 ww3_bound
 
$(aPo)/ww3_gint.o : ww3_gint.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3idatmd.o $(aPo)/w3iogrmd.o $(aPo)/w3iogomd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o $(aPo)/w3arrymd.o $(aPo)/w3gsrumd.o $(aPo)/w3triamd.o
	@$(aPb)/ad3 ww3_gint
 
$(aPo)/ww3_grib.o : ww3_grib.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3iogrmd.o $(aPo)/w3iogomd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o
	@$(aPb)/ad3 ww3_grib
 
$(aPo)/ww3_grid.o : ww3_grid.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3odatmd.o $(aPo)/w3iogrmd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3arrymd.o $(aPo)/w3dispmd.o $(aPo)/w3gsrumd.o $(aPo)/w3triamd.o
	@$(aPb)/ad3 ww3_grid
 
$(aPo)/ww3_gspl.o : ww3_gspl.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3iogrmd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3arrymd.o
	@$(aPb)/ad3 ww3_gspl
 
$(aPo)/ww3_multi.o : ww3_multi.ftn w3macros.h  $(aPo)/wminitmd.o $(aPo)/wmwavemd.o $(aPo)/wmfinlmd.o $(aPo)/wmmdatmd.o
	@$(aPb)/ad3 ww3_multi
 
$(aPo)/ww3_ounf.o : ww3_ounf.ftn w3macros.h  $(aPo)/w3initmd.o $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3iogrmd.o $(aPo)/w3iogomd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o $(aPo)/w3arrymd.o
	@$(aPb)/ad3 ww3_ounf
 
$(aPo)/ww3_ounp.o : ww3_ounp.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3iogrmd.o $(aPo)/w3iopomd.o $(aPo)/w3partmd.o $(aPo)/w3sln1md.o $(aPo)/w3src4md.o $(aPo)/w3snl1md.o $(aPo)/w3sbt1md.o $(aPo)/w3sdb1md.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o $(aPo)/w3arrymd.o $(aPo)/w3dispmd.o
	@$(aPb)/ad3 ww3_ounp
 
$(aPo)/ww3_outf.o : ww3_outf.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3iogrmd.o $(aPo)/w3iogomd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o $(aPo)/w3arrymd.o
	@$(aPb)/ad3 ww3_outf
 
$(aPo)/ww3_outp.o : ww3_outp.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3iogrmd.o $(aPo)/w3iopomd.o $(aPo)/w3partmd.o $(aPo)/w3bullmd.o $(aPo)/w3sln1md.o $(aPo)/w3src4md.o $(aPo)/w3snl1md.o $(aPo)/w3sbt1md.o $(aPo)/w3sdb1md.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o $(aPo)/w3arrymd.o $(aPo)/w3dispmd.o
	@$(aPb)/ad3 ww3_outp
 
$(aPo)/ww3_prep.o : ww3_prep.ftn w3macros.h  $(aPo)/w3fldsmd.o $(aPo)/w3gdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3iogrmd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o $(aPo)/w3arrymd.o $(aPo)/w3gsrumd.o
	@$(aPb)/ad3 ww3_prep
 
$(aPo)/ww3_prnc.o : ww3_prnc.ftn w3macros.h  $(aPo)/w3fldsmd.o $(aPo)/w3gdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3iogrmd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o $(aPo)/w3arrymd.o $(aPo)/w3gsrumd.o
	@$(aPb)/ad3 ww3_prnc
 
$(aPo)/ww3_sbs1.o : ww3_sbs1.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3idatmd.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o $(aPo)/wminitmd.o $(aPo)/wmwavemd.o $(aPo)/wmfinlmd.o $(aPo)/wmmdatmd.o $(aPo)/wmunitmd.o
	@$(aPb)/ad3 ww3_sbs1
 
$(aPo)/ww3_shel.o : ww3_shel.ftn w3macros.h  $(aPo)/w3initmd.o $(aPo)/w3wavemd.o $(aPo)/w3wdasmd.o $(aPo)/w3fldsmd.o $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3idatmd.o $(aPo)/w3iogrmd.o $(aPo)/w3iogomd.o $(aPo)/w3iopomd.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o
	@$(aPb)/ad3 ww3_shel
 
$(aPo)/ww3_strt.o : ww3_strt.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3wdatmd.o $(aPo)/w3adatmd.o $(aPo)/w3odatmd.o $(aPo)/w3iogrmd.o $(aPo)/w3iorsmd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3arrymd.o $(aPo)/w3dispmd.o $(aPo)/w3gsrumd.o
	@$(aPb)/ad3 ww3_strt
 
$(aPo)/ww3_systrk.o : ww3_systrk.ftn w3macros.h  $(aPo)/w3strkmd.o
	@$(aPb)/ad3 ww3_systrk
 
$(aPo)/ww3_trck.o : ww3_trck.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3odatmd.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o
	@$(aPb)/ad3 ww3_trck
 
$(aPo)/ww3_trnc.o : ww3_trnc.ftn w3macros.h  $(aPo)/w3gdatmd.o $(aPo)/w3odatmd.o $(aPo)/constants.o $(aPo)/w3servmd.o $(aPo)/w3timemd.o
	@$(aPb)/ad3 ww3_trnc
 
# end of WAVEWATCH III subroutines
