#!/bin/csh
#################################################################
# Csh Script to retrieve 58 online Data files of 'ds260.2',
# total 1.03G. This script uses 'wget' to download data.
#
# Highlight this script by Select All, Copy and Paste it into a file;
# make the file executable and run it on command line.
#
# You need pass in your password as a parameter to execute
# this script; or you can set an environment variable RDAPSWD
# if your Operating System supports it.
#
# Contact zji@ucar.edu (Zaihua Ji) for further assistance.
#################################################################

set pswd = $1
if(x$pswd == x && `env | grep RDAPSWD` != '') then
 set pswd = $RDAPSWD
endif
if(x$pswd == x) then
 echo
 echo Usage: $0 YourPassword
 echo
 exit 1
endif
set v = `wget -V |grep 'GNU Wget ' | cut -d ' ' -f 3`
set a = `echo $v | cut -d '.' -f 1`
set b = `echo $v | cut -d '.' -f 2`
if(100 * $a + $b > 109) then
 set opt = 'wget --no-check-certificate'
else
 set opt = 'wget'
endif
set opt1 = '-O Authentication.log --save-cookies auth.rda_ucar_edu --post-data'
set opt2 = "email=mahajans@ornl.gov&passwd=$pswd&action=login"
$opt $opt1="$opt2" https://rda.ucar.edu/cgi-bin/login
set opt1 = "-N --load-cookies auth.rda_ucar_edu"
set opt2 = "$opt $opt1 http://rda.ucar.edu/data/ds260.2/"
set filelist = ( \
  1949.nc.gz \
  1950.nc.gz \
  1951.nc.gz \
  1952.nc.gz \
  1953.nc.gz \
  1954.nc.gz \
  1955.nc.gz \
  1956.nc.gz \
  1957.nc.gz \
  1958.nc.gz \
  1959.nc.gz \
  1960.nc.gz \
  1961.nc.gz \
  1962.nc.gz \
  1963.nc.gz \
  1964.nc.gz \
  1965.nc.gz \
  1966.nc.gz \
  1967.nc.gz \
  1968.nc.gz \
  1969.nc.gz \
  1970.nc.gz \
  1971.nc.gz \
  1972.nc.gz \
  1973.nc.gz \
  1974.nc.gz \
  1975.nc.gz \
  1976.nc.gz \
  1977.nc.gz \
  1978.nc.gz \
  1979.nc.gz \
  1980.nc.gz \
  1981.nc.gz \
  1982.nc.gz \
  1983.nc.gz \
  1984.nc.gz \
  1985.nc.gz \
  1986.nc.gz \
  1987.nc.gz \
  1988.nc.gz \
  1989.nc.gz \
  1990.nc.gz \
  1991.nc.gz \
  1992.nc.gz \
  1993.nc.gz \
  1994.nc.gz \
  1995.nc.gz \
  1996.nc.gz \
  1997.nc.gz \
  1998.nc.gz \
  1999.nc.gz \
  2000.nc.gz \
  2001.nc.gz \
  2002.nc.gz \
  2003.nc.gz \
  2004.nc.gz \
  2005.nc.gz \
  2006.nc.gz \
)
while($#filelist > 0)
 set syscmd = "$opt2$filelist[1]"
 echo "$syscmd ..."
 $syscmd
 shift filelist
end

rm -f auth.rda_ucar_edu Authentication.log
exit 0
