// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//
#include <stdio.h>
#include <stdarg.h>

#define MAX_LINE_LEN 132

char printbuf[MAX_LINE_LEN+2];
char fbuffer[1024];
int nbuf = 0;

void fortprintf(FILE * fd, char * str, ...)
{
   int i, nl, sp, inquotes;
   va_list ap;

   va_start(ap, str);
   i = vsnprintf(fbuffer+nbuf, 1024, str, ap);
   va_end(ap);

   nbuf = nbuf + i;

   inquotes = 0;
   do {
      nl = sp = -1;
      for(i=0; i<MAX_LINE_LEN-1 && i<nbuf; i++) {
         if (fbuffer[i] == '\'' && (fbuffer[i+1] != '\'' || i == nbuf-1)) inquotes = (inquotes + 1) % 2;
         if (fbuffer[i] == '\n') nl = i;
         if (fbuffer[i] == ' ' && i != nbuf-1 && fbuffer[i+1] != '&') sp = i;
      }
      if (nbuf <= MAX_LINE_LEN) sp = -1;

      if (nl > 0) {
         snprintf(printbuf, nl+2, "%s", fbuffer);
         fprintf(fd, "%s", printbuf);
         nl++;
         for(i=0; nl<nbuf; i++, nl++)
            fbuffer[i] = fbuffer[nl];
         nbuf = i;
      }
      else if (sp > 0) {
         snprintf(printbuf, sp+2, "%s", fbuffer);
         i = sp+1;
         if (inquotes) printbuf[i++] = '\'';
         printbuf[i++] = '&';
         printbuf[i++] = '\n';
         printbuf[i++] = '\0';
         fprintf(fd, "%s", printbuf);
         sp++;
         i = 0;
         if (inquotes) {
            inquotes = (inquotes + 1) % 2;
            fbuffer[i++] = '/';
            fbuffer[i++] = '/';
            fbuffer[i++] = '\'';
         }
         for( ; sp<nbuf; i++, sp++)
            fbuffer[i] = fbuffer[sp];
         nbuf = i;
      }
   } while (nl > 0 || sp > 0);

}

void fortprint_flush(FILE * fd)
{
   snprintf(printbuf, nbuf+1, "%s", fbuffer);
   fprintf(fd, "%s", printbuf);
   nbuf = 0;
}
