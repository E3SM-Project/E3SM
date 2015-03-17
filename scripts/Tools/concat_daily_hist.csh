#!/bin/csh -f
#
# This script concatenates daily history files output by the coupler into
# monthly history files.  For each type of history file, a domain file is
# created containing only the domain/grid metadata.  Next, if the correct
# number of daily files exist, the data are concatenated into a single monthly
# file, but the domain metadata are not stored in the new monthly files.
# If the concatenation is successful, the original daily files are removed.
# This script is designed to be run in the $RUNDIR directory prior to
# running the short term archiving script.
#
# Requirements:
#   Environment variables: RUNDIR, CASE
#   Packages: csh, Unix tools, ncdump, NCO Toolkit (ncks and ncrcat)
#
# Adapted by Forrest M. Hoffman from a script by Keith Lindsay and
# Mariana Vertenstein.
# Created: Sat Jul  4 21:43:11 MDT 2009
# Updated: Sun Jul  5 15:44:16 MDT 2009

set dpm = (31 28 31 30 31 30 31 31 30 31 30 31)

set case_fields = `echo ${CASE} | awk -F. '{print NF}'`
@ hist_field = ${case_fields} + 2
@ date_field = ${case_fields} + 3
set comp = cpl

cd ${RUNDIR}

foreach hist_type (`ls | grep -E "^${CASE}\.${comp}\.h.+\.[0-9][0-9][0-9][0-9]-[0-1][0-9]-[0-3][0-9]" | grep '\.nc$' | awk -F. '{print $'${hist_field}'}' | sort | uniq`)
   echo "Processing files for ${comp} history type ${hist_type}"
   # Only the first history file may contain domain data.  If it does, the
   # domain data should be extracted into a separate domain file.
   set dom_src = `ls | grep "^${CASE}\.${comp}\.${hist_type}\.[0-9][0-9][0-9][0-9]-[0-1][0-9]-[0-3][0-9]" | grep '\.nc$' | head -1`
   set dom_cnt = `ncdump -h ${dom_src} | awk '/double/ || /float/ {num = split($2,a, "("); var = a[1]; if (var ~ "dom._.+") print var;}' | wc -l`
   if (${dom_cnt} > 0) then
      set dom_file = "${CASE}.${comp}.${hist_type}.domain.nc"
      rm -f ${dom_file}
      echo "   Creating domain file ${dom_file} from history file ${dom_src}"
      ncks -v '^dom._.+' ${dom_src} ${dom_file}
   else
      echo "   No domain data found in the first ${hist_type} file"
   endif
   foreach year (`ls | grep "^${CASE}\.${comp}\.${hist_type}\.[0-9][0-9][0-9][0-9]-[0-1][0-9]-[0-3][0-9]" | grep '\.nc$' | awk -F. '{print $'${date_field}'}' | awk -F- '{print $1}' | sort | uniq`)
      echo "   Processing daily history files for year ${year}"
      foreach month (`ls | grep "^${CASE}\.${comp}\.${hist_type}\.${year}-[0-1][0-9]-[0-3][0-9]" | grep '\.nc$' |  awk -F. '{print $'${date_field}'}' | awk -F- '{print $2}' | sort | uniq`)
         @ dom_cnt = 0
         @ cnt = 0
         foreach file (`ls | grep "^${CASE}\.${comp}\.${hist_type}\.${year}-${month}-[0-3][0-9]" | grep '\.nc$'`)
            @ dom_cnt = ${dom_cnt} + `ncdump -h ${file} | awk '/double/ || /float/ {num = split($2,a, "("); var = a[1]; if (var ~ "dom._.+") print var;}' | wc -l`
            @ cnt = ${cnt} + 1
         end
         # Do not check to see if we have enough files for a month because some
         # files (like that for runoff) may not gerenate a file until the end
         # of the first day.
         #if (${cnt} == ${dpm[${month}]}) then
            echo "      Concatenating data from files for month ${month}"
            if (${dom_cnt} > 0) then
               echo "         (stripping out domain data)"
               ls | grep "^${CASE}\.${comp}\.${hist_type}\.${year}-${month}-[0-3][0-9]" | grep '\.nc$' | ncrcat -h -x -v '^dom._.+' -o ${CASE}.${comp}.${hist_type}.${year}-${month}.nc && ls | grep "^${CASE}\.${comp}\.${hist_type}\.${year}-${month}-[0-3][0-9]" | grep '\.nc$' | xargs rm -f
            else
               echo "         (contains no domain data)"
               ls | grep "^${CASE}\.${comp}\.${hist_type}\.${year}-${month}-[0-3][0-9]" | grep '\.nc$' | ncrcat -h -o ${CASE}.${comp}.${hist_type}.${year}-${month}.nc && ls | grep "^${CASE}\.${comp}\.${hist_type}\.${year}-${month}-[0-3][0-9]" | grep '\.nc$' | xargs rm -f
            endif
         #else
         #   echo "      Incorrect number of daily files (${cnt}) for month ${month} . . . skipping"
         #endif
      end
   end
   echo ""
end
