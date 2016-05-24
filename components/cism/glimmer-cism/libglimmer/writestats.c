/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   +                                                              +
   +  writestats.c - part of the Glimmer-CISM ice model           +
   +                                                              +
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   Copyright (C) 2009, 2010
   Glimmer-CISM contributors - see AUTHORS file for list of contributors

   This file is part of Glimmer-CISM.

   Glimmer-CISM is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or (at
   your option) any later version.

   Glimmer-CISM is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Glimmer-CISM.  If not, see <http://www.gnu.org/licenses/>.

   Glimmer-CISM is hosted on BerliOS.de:
   https://developer.berlios.de/projects/glimmer-cism/

 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "writestats.h"
#include "config.inc"

#include <stdio.h>
#include <errno.h>
#include <time.h>
#include <sys/times.h>
#include <sys/utsname.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>



#define CFG_LEN 35
#define BUFFER_LEN 400
#define PERM_FILE (S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)

void FC_FUNC(gf_writestats,GF_WRITESTATS) (const char *resname, const char *cfgname, double wallTime)
{
  struct tms runtime; 
  clock_t clck;
  double userTime,sysTime;
  time_t now;
  struct tm * timeAndDate;
  char dateStr[20];
  struct utsname unameData;
  char outBuffer[BUFFER_LEN+1];
  char hdrBuffer[BUFFER_LEN+1];
  int fd;
  int i,haveLock;
  struct flock fl = { F_WRLCK, SEEK_SET, 0,       0,     0 };
  off_t fileLength;

  /* Jeff: Trim the file name on Jaguar, because spaces are embedded otherwise. */
  char resfname[CFG_LEN+1];
  strncpy(resfname, resname, CFG_LEN);
  i = 0;
  while (isalnum(resfname[i]) && i < CFG_LEN) i++;
  resfname[i] = '\0';

  /* get user and system time */
  clck = times(&runtime);
  userTime = ((double) runtime.tms_utime)/((double) sysconf(_SC_CLK_TCK));
  sysTime = ((double) runtime.tms_stime)/((double) sysconf(_SC_CLK_TCK));

  /* get the current data */
  now = time(NULL);
  timeAndDate = localtime(&now);
  snprintf(dateStr,20,"%4d-%02d-%02d_%02d:%02d",timeAndDate->tm_year+1900, timeAndDate->tm_mon+1, timeAndDate->tm_mday, timeAndDate->tm_hour, timeAndDate->tm_min);

  /* get host name and architecture */
  if ((uname(&unameData))!=0) {
    unameData.nodename[0]='\0';
    unameData.machine[0]='\0';
  }

  /* construct output line */
  snprintf(outBuffer,BUFFER_LEN,"%*s %9.2f %9.2f %8.2f %s %-10s %-6s %-10s \"%s\"\n",-CFG_LEN,cfgname, wallTime, userTime, sysTime, dateStr, \
	   unameData.nodename, unameData.machine, VERSION, GLIMMER_FCFLAGS);
  snprintf(hdrBuffer,BUFFER_LEN,"%*s %9s %9s %-8s %-16s %-10s %-6s %-10s %s\n",-CFG_LEN,"#cfg_file","wall_time","usr_time","sys_time","date","host","arch","version","FCFLAGS");

  /* open output file */
  if ((fd = open(resfname, O_CREAT|O_WRONLY|O_SYNC,PERM_FILE)) == -1) {
    perror("opening result file");
    printf("%s\n",outBuffer);
    return;
  }

  /* get a lock on the file */
  i=0;
  while ((haveLock=fcntl(fd, F_SETLK, &fl))==-1 && i<100000)  i++;
  if (haveLock==-1) {
    close(fd);
    perror("getting lock");
    printf("%s\n",outBuffer);
    return;
  }

  /* go to the end of the file */
  fileLength = lseek(fd,0,SEEK_END);

  /* write data */
  if (fileLength == 0)
    write(fd,hdrBuffer,strlen(hdrBuffer));
  write(fd,outBuffer,strlen(outBuffer));

  /* release the lock */
  fl.l_type = F_UNLCK;
  if (fcntl(fd, F_SETLK, &fl) == -1) {
    perror("unlocking file");
    return;
  }
  
  close(fd);

  return;
}
