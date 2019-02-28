#!/bin/sh
# --------------------------------------------------------------------------- #
# ww3_gspl.sh: Shell wrapper tot split grid into subgrids, including making   #
#              all mod_def files and the relvant part of ww3_multi.inp.       #
#                                                                             #
# use        : ww3_gspl.sh [options] gridname nr_grids                        #
#              See usage function for options.                                #
#                                                                             #
# error codes : Program ends if error occurs in input.                        #
#                                                                             #
# programs used :                                                             #
#       ww3_grid : Grid preprocessor.                                         #
#       ww3_gspl : Grdi splitter.                                             #
#                                                                             #
#                                                      Hendrik L. Tolman      #
#                                                      January 2014           #
#                                                                             #
#    Copyright 2013-2014 National Weather Service (NWS),                      #
#       National Oceanic and Atmospheric Administration.  All rights          #
#       reserved.  WAVEWATCH III is a trademark of the NWS.                   #
#       No unauthorized use without permission.                               #
#                                                                             #
# --------------------------------------------------------------------------- #
# 1. Preparations                                                             #
# --------------------------------------------------------------------------- #
# 1.a Internal variables

  set -e

# 1.a.1 Setup file

# The following line must not be removed: it is a switch for local install
# so that all bin scripts point to the local wwatch3.env
  export ww3_env=/users/sbrus/climate/ww3/wwatch3.env
# For manual install (without install_ww3_tar or install_ww3_svn) make sure to
# either use the generic ww3_env or to add your own ww3_env="${my_directory}"

  if [ ${WWATCH3_ENV} ]; then ww3_env="${WWATCH3_ENV}"; fi # alternate setup file

  home_dir=`pwd`

   info_bot="3  2  1.0  '(12F11.3)'"
  info_obst="3  2  1.0  '(26F5.2)'"
  info_mask="3  2  1    '(66I2)'"  

# 1.a.2 Usage function

  scriptname="`basename $0`"
  optstr="ad:e:f:hil:n:o:rs:t:v"

  nr_it=350
  target='0.75'
  halo_ext=2
  comm_first='0.'
  comm_last='1.'
  frflag='T'

  usage ()
{
cat 2>&1 << EOF

Usage: $scriptname [options]  gridID nr_grid
Required:
  gridID     : name of master grid to be split up
  nr_grid    : number of sub-grids to be generated
Options:
  -a               : use entire assigned cummunicator for each grid
  -h               : help, print this.
  -i               : create template file ww3_gint.inp_tmpl for
                     later integration of output into single grid.
  -d data_dir      : directory with ww3_grid.inp and ancilary data
                      * default is working directory
                      * relative unless starting with '/'
  -e halo_ext      : set halo extension, default is 2
  -o output_dir    : directory for std out redirects
                      * default is working directory
                      * relative unless starting with '/'
  -n n_iter        : maximum number of  interations in ww3_gspl
                      * default = $nr_it
  -t target        : target accuracy in ww3_gspl (%)
                      * default = $target
  -f comm_first    : communicator fraction (first).
                      * default = $comm_first
  -l comm_last     : communicator fraction (last).
                      * default = $comm_last
  -s ww3_multi.inp : name of input file to be modified.
                      * Not set as default.
  -r               : replace file defined under -s, otherwise add .new
  -v               : verbose, show program output
EOF
}
 
# 1.a.3 Process input (part 1)

  args=`getopt $optstr $*`

  if [ $? != 0 ]
  then
    usage
    exit 1
  fi

  set -- $args

  while :
  do
    case "$1" in
    -a) frflag='F' ;;
    -c) shift; info_comm="$1" ;;
    -d) shift; data_dir="$1" ;;
    -e) shift; halo_ext="$1" ;;
    -f) shift; comm_first="$1" ;;
    -h) help=1 ;;
    -i) gint=1 ;;
    -l) shift; comm_last="$1" ;;
    -n) shift; nr_it="$1" ;;
    -o) shift; outp_dir="$1" ;;
    -r) replace=1 ;;
    -s) shift; inp_file="$1" ;;
    -t) shift; target="$1" ;;
    -v) verbose=1 ;;
    --) break ;;
    esac
    shift
  done
  shift

  if [ $help ]
  then
    usage
    exit 1
  fi

  if [ ! $# = 0 ]
  then
    modID="$1" ; shift
  else
    usage
    exit 2
  fi

  if [ ! $# = 0 ]
  then
    nr_grid="$1" ; shift
    if [ "$nr_grid" -lt '1' ] ; then
      nr_grid='1' ; fi
  else
    usage
    exit 3
  fi

# 1.a.4 Setup file processing

  if test -f $ww3_env
  then
    set `grep WWATCH3_DIR $ww3_env` ; shift
    main_dir="$*"
    set `grep WWATCH3_TMP $ww3_env` ; shift
    temp_dir="$*"
    set `grep WWATCH3_SOURCE $ww3_env` ; shift
    source="$*"
    set `grep WWATCH3_LIST $ww3_env` ; shift
    list="$*"
  else
    echo "*** Set-up file $ww3_env not found ***"; echo ' '
    exit 1
  fi

  exe_dir="$main_dir/exe"

# private temp dir to allow for parallel batch processing

  temp_dir=$home_dir/ww3_gspl.temp
 
# 1.a.5 Process input (part 2)


  if [ -z "$data_dir" ] ; then
    data_dir=$home_dir ; fi

  if [ "`echo $data_dir | cut -c1-2`" = './' ] 
  then
    len=`echo $data_dir | wc -c | awk '{print $1}'`
    data_dir="$home_dir/`echo $data_dir | cut -c3-$len`"
  fi

  if [ "`echo $data_dir | cut -c1-1`" != '/' ] ; then
    data_dir="$home_dir/$data_dir" ; fi

  if [ -z "$outp_dir" ] ; then
    outp_dir=$home_dir ; fi

  if [ "`echo $outp_dir | cut -c1-2`" = './' ] 
  then
    len=`echo $outp_dir | wc -c | awk '{print $1}'`
    outp_dir="$home_dir/`echo $outp_dir | cut -c3-$len`"
  fi

  if [ "`echo $outp_dir | cut -c1-1`" != '/' ] ; then
    outp_dir="$home_dir/$outp_dir" ; fi

  mkdir -p $outp_dir

# 1.b ID header  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  echo ' '
  echo '                 *********************************'
  echo '               *** splitting WAVEWATCH III grids ***'
  echo '                 *********************************'
  echo ' '
  echo " Running from      : $home_dir"
  echo " Main directory    : $main_dir"
  echo " Scratch directory : $temp_dir"
  echo " Data directory    : $data_dir"
  echo " Output directory  : $outp_dir"
  echo " Model ID          : $modID"
  echo " Number of grids   : $nr_grid"
  echo ' '

  if [ "$frflag" = 'F' ] 
  then
    echo " Each grid uses maximum communicator (-a) *** ABNORMAL OPERATION ***"
    echo ' '
  fi

# 1.d Set up work space  - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 1.d.1 Genmeral preparations

#  ${main_dir}/bin/w3_clean

  cd $home_dir
  rm -f $inp_file.new
  rm -f ww3_mask.$modID.*
  rm -f mod_def.$modID
  rm -f mod_def.$modID.*
  rm -f ww3_multi.$modID.*
  rm -f part_*

  cd $data_dir
  filelist=`ls`

  cd $outp_dir
  rm -f *.out

  rm -rf $temp_dir
  mkdir -p $temp_dir
  cd $temp_dir
  rm -f *.$modID $modID.* *.$modID.*
  rm -f ww3_grid.*
  rm -f *.bot
  rm -f *.obst
  rm -f *.mask

# 1.d.2 Get files from data directory

  source1=ww3_grid.$modID
  source2=$modID.inp

  if [ -f $data_dir/$source1 ] 
  then
    cp $data_dir/$source1 . 
    echo " Found $source1 ..."
  fi

  if [ -f $data_dir/$source2 ] 
  then
    cp $data_dir/$source2 $source1 
    echo " Found $source2 [renamed $source1] ..."
  fi

  if [ ! -f $source1 ] ; then
    echo "    *** file $source1 not found ; ABORT ***" ; echo ' ' ; exit 1 ; fi

  ifiles="$source1"

  for file in $filelist
  do
    occur=`grep $file $source1 | wc -l | awk '{ print $1}'`
    if [ "$occur" != '0' ] ; then
      cp $data_dir/$file . 
      echo " Found $file ..."
      ifiles="$ifiles $file"
    fi
  done

  sed -n '/^\$.*/!p' $source1 | sed -n '5,/^END OF.*/p' > namelist.data
  echo " Made namelist.data ..."


  if [ -z "$inp_file" ]
  then
    echo " Source for ww3_multi.inp not defined ..."
  else
    echo " Using ww3_multi.inp file [$inp_file] ..."
    if [ $replace ]
    then
       echo "    Replacing this file."
       out_file=$inp_file
    else
       out_file=$inp_file.new
       echo "    New file will be $out_file."
    fi
  fi

  if [ -f $data_dir/$inp_file ] ; then
    cp $data_dir/$inp_file . ; fi

  echo ' '

# --------------------------------------------------------------------------- #
# 2. Make master mod_def file                                                 #
# --------------------------------------------------------------------------- #

  ln -s $source1 ww3_grid.inp

# if [ $verbose ]
# then
#   echo " Making master mod_def file ..."
#   $exe_dir/ww3_grid
# else
    echo " Making master mod_def file, output routed to "
    echo "        $outp_dir/ww3_grid.$modID.out ..."
    $exe_dir/ww3_grid > $outp_dir/ww3_grid.$modID.out
# fi

  rm -f ww3_grid.inp
  mv mod_def.ww3 mod_def.$modID

  for file in $ifiles
  do
    rm -f $file
  done

  echo ' '

# --------------------------------------------------------------------------- #
# 2. Run ww3_gspl                                                             #
# --------------------------------------------------------------------------- #

  if [ "$nr_grid" = '1' ] 
  then

    mv mod_def.$modID $home_dir
    echo ' No splitting of grid required'
    echo ' '
    echo '                     *************************'
    echo '                   *** end of grid splitting ***'
    echo '                     *************************'
    echo ' '
    exit 0
  fi

  if [ $verbose ]
  then
    echo " Run ww3_gspl to get information for sub-grids ..."
  else
    echo " Run ww3_gspl to get information for sub-grids, output routed to"
    echo "        $outp_dir/ww3_gspl.out ..."
  fi

cat > ww3_gspl.inp << EOF
$ -------------------------------------------------------------------- $
$ WAVEWATCH III Grid splitting input file                              $
$ -------------------------------------------------------------------- $
  '$modID'
$
  $nr_grid  $nr_it  $target $halo_ext
$
  $info_bot
  $info_obst
  $info_mask
$
  $comm_first  $comm_last  $frflag
$ -------------------------------------------------------------------- $
$ End of input file                                                    $
$ -------------------------------------------------------------------- $
EOF

  if [ $verbose ]
  then
    $exe_dir/ww3_gspl
    OK="$?"
  else
    $exe_dir/ww3_gspl > $outp_dir/ww3_gspl.out 2>&1
    OK="$?"
  fi
  rm -f ww3_gspl.inp

  if [ "$OK" != 0 ] ; then
    echo ' ' ; echo '   *** ERROR ww3_gspl (ABORT) ***', echo ' ' 
    cat $outp_dir/ww3_gspl.out
    exit 2 ; fi

# --------------------------------------------------------------------------- #
# 3. Process individual grids                                                 #
# --------------------------------------------------------------------------- #

  echo ' '
  echo " Process individual grids ..."

  for file in `ls *.tmpl`
  do
    grdID=`echo $file | sed 's/\./ /g' | awk '{ print $1 }'`
    echo "   Grid $grdID (ww3_grid.$grdID.out) ..."

    sed -e "/NAMELISTS/r namelist.data" \
        -e 's/ NAMELISTS//g' $file > ww3_grid.inp

    $exe_dir/ww3_grid > $outp_dir/ww3_grid.$grdID.out

    mv mod_def.ww3 mod_def.$grdID
    rm -f ww3_grid.inp
    rm -f $grdID.*
  done

  rm -f namelist.data

# --------------------------------------------------------------------------- #
# 4. Process ww3_gint.inp_tmpl                                                #
# --------------------------------------------------------------------------- #

  if [ "$gint" = '1' ]
  then
    echo ' '
    echo " Making ww3_gint.inp_tmpl ..."
    nr_gint=`expr $nr_grid + 1`

cat > ww3_gint.inp_tmpl << EOF
$ -------------------------------------------------------------------- $
$ WAVEWATCH III Grid integration input file                            $
$ -------------------------------------------------------------------- $
  DATE TIME TSTEP NR_STEPS
$
  $nr_gint
$
EOF

  cat ww3_multi.$modID.$nr_grid | awk '{ print $1 }' >> ww3_gint.inp_tmpl
  echo "'$modID'"                                    >> ww3_gint.inp_tmpl

cat >> ww3_gint.inp_tmpl << EOF
$ -------------------------------------------------------------------- $
$ End of input file                                                    $
$ -------------------------------------------------------------------- $
EOF

  fi

# --------------------------------------------------------------------------- #
# 5. Process ww3_multi.inp file                                               #
# --------------------------------------------------------------------------- #

  if [ -z $inp_file ]
  then
    echo ' '
    echo " No file ww3_multi.inp to process ..."
  else
    echo ' '
    echo " Processing file ww3_multi.inp [$inp_file] ..."

    sed -n '/^\$.*/!p' $inp_file > tempfile
    dataold="`head -1 tempfile`"
    rm -f tempfile
    ngrold=`echo $dataold | awk '{ print $1}'`
    ngrnew=`expr $ngrold - 1 + $nr_grid`
    datanew=`echo "$dataold" | sed "s/$ngrold/$ngrnew/"`

    data=`grep $modID $inp_file | head -1 | sed "s/'//g"`
    lev=`echo $data | awk '{ print $2 }'`
    cur=`echo $data | awk '{ print $3 }'`
    wnd=`echo $data | awk '{ print $4 }'`
    ice=`echo $data | awk '{ print $5 }'`
    dt1=`echo $data | awk '{ print $6 }'`
    dt2=`echo $data | awk '{ print $7 }'`
    dt3=`echo $data | awk '{ print $8 }'`
    rnk=`echo $data | awk '{ print $9 }'`
    grp=`echo $data | awk '{ print $10 }'`
    bfl=`echo $data | awk '{ print $13 }'`

    sed -e "s/LEV/$lev/g" \
        -e "s/CUR/$cur/g" \
        -e "s/WND/$wnd/g" \
        -e "s/ICE/$ice/g" \
        -e "s/D1/$dt1/g" \
        -e "s/D2/$dt2/g" \
        -e "s/D3/$dt3/g" \
        -e "s/RANK/$rnk/g" \
        -e "s/GROUP/$grp/g" \
        -e "s/BFLAG/$bfl/g" \
       ww3_multi.$modID.$nr_grid  > ww3_multi.part

    sed -e "s/$dataold/$datanew/g" \
        -e "/$modID/r ww3_multi.part" \
        -ne "/'$modID'/!p" $inp_file > $inp_file.temp

    rm -f ww3_multi.part
    rm -f ww3_multi.$modID.$nr_grid
    if [ "$inp_file" != "$out_file" ] ; then
      rm -f $inp_ifle ; fi
    mv $inp_file.temp $out_file

  fi

# --------------------------------------------------------------------------- #
# 6. saving files                                                             #
# --------------------------------------------------------------------------- #

  echo ' '
  echo " Saving files in $home_dir ..."

  echo "    mod_def files ..."
  mv mod_def.* $home_dir

  for file in ww3_mask.$modID.$nr_grid $out_file ww3_multi.$modID.$nr_grid \
              ww3_gint.inp_tmpl
  do
    if [ -f $file ]
    then
      echo "    file $file ..."
      mv $file $home_dir
    fi
  done

  if [ -f ww3.ctl ]
  then
      echo "    file ww3.ctl ..."
      mv ww3.ctl $home_dir
      echo "    file ww3.ww3_gspl ..."
      mv ww3.ww3_gspl $home_dir
  fi


# --------------------------------------------------------------------------- #
# 6. End of program ID / clean up                                             #
# --------------------------------------------------------------------------- #

  cd $home_dir
  rm -rf $temp_dir
  echo ' '
  echo '                     *************************'
  echo '                   *** end of grid splitting ***'
  echo '                     *************************'
  echo ' '


# End of ww3_gspl.sh -------------------------------------------------------- #
