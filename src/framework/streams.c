// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>

#ifdef UNDERSCORE
#define open_streams open_streams_
#define close_streams close_streams_
#else
#ifdef DOUBLEUNDERSCORE
#define open_streams open_streams__
#define close_streams close_streams__
#endif
#endif

int fd_out, fd_err;

void open_streams(int * id)
{
   char fname[128];

#ifndef MPAS_DEBUG
   if(*id == 0){
	   sprintf(fname, "log.%4.4i.err", *id);
	   fd_err = open(fname,O_CREAT|O_WRONLY|O_TRUNC,0644);
	   if (dup2(fd_err, 2) < 0) {
		   printf("Error duplicating STDERR\n");
		   return;
	   }

	   sprintf(fname, "log.%4.4i.out", *id);
	   fd_out = open(fname,O_CREAT|O_WRONLY|O_TRUNC,0644);
	   if (dup2(fd_out, 1) < 0) {
		   printf("Error duplicating STDOUT\n");
		   return;
	   }
   } else {
	   sprintf(fname, "/dev/null");
	   fd_err = open(fname,O_CREAT|O_WRONLY|O_TRUNC,0644);
	   if (dup2(fd_err, 2) < 0) {
		   printf("Error duplicating STDERR\n");
		   return;
	   }

	   sprintf(fname, "/dev/null");
	   fd_out = open(fname,O_CREAT|O_WRONLY|O_TRUNC,0644);
	   if (dup2(fd_out, 1) < 0) {
		   printf("Error duplicating STDOUT\n");
		   return;
	   }
   }
#else
   sprintf(fname, "log.%4.4i.err", *id);
   fd_err = open(fname,O_CREAT|O_WRONLY|O_TRUNC,0644);
   if (dup2(fd_err, 2) < 0) {
      printf("Error duplicating STDERR\n");
      return;
   }

   sprintf(fname, "log.%4.4i.out", *id);
   fd_out = open(fname,O_CREAT|O_WRONLY|O_TRUNC,0644);
   if (dup2(fd_out, 1) < 0) {
      printf("Error duplicating STDOUT\n");
      return;
   }
#endif
}

void close_streams()
{

}
