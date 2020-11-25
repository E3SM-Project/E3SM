// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//
/* File: read_geogrid.c

   Sample subroutine to read an array from the geogrid binary format.

   Michael G. Duda, NCAR/MMM
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#define GEOG_BIG_ENDIAN    0
#define GEOG_LITTLE_ENDIAN 1

/*  In Fortran, use the following as an interface for read_geogrid:

 use iso_c_binding, only : c_char, c_int, c_float, c_ptr, c_loc

 interface
    subroutine read_geogrid(fname, rarray, nx, ny, nz, isigned, endian, &
                            wordsize, status) bind(C)
       use iso_c_binding, only : c_char, c_int, c_float, c_ptr
       character (c_char), dimension(*), intent(in) :: fname
       type (c_ptr), value :: rarray
       integer (c_int), intent(in), value :: nx
       integer (c_int), intent(in), value :: ny
       integer (c_int), intent(in), value :: nz
       integer (c_int), intent(in), value :: isigned
       integer (c_int), intent(in), value :: endian
       integer (c_int), intent(in), value :: wordsize
       integer (c_int), intent(inout) :: status
    end subroutine read_geogrid
 end interface

*/

int read_geogrid(
      char * fname,            /* The name of the file to read from */
      float * rarray,          /* The array to be filled */
      int nx,                /* x-dimension of the array */
      int ny,                /* y-dimension of the array */
      int nz,                /* z-dimension of the array */
      int isigned,           /* 0=unsigned data, 1=signed data */
      int endian,            /* 0=big endian, 1=little endian */
      int wordsize,          /* number of bytes to use for each array element */
      int * status)
{
   int i, ival, cnt, narray;
   int A2, B2;
   int A3, B3, C3;
   int A4, B4, C4, D4;
   unsigned char * c;
   FILE * bfile;

   *status = 0;

   narray = (nx) * (ny) * (nz);

   /* Attempt to open file for reading */
   if (!(bfile = fopen(fname,"rb")))
   {
      *status = 1;
      return 1;
   }

   /* Allocate memory to hold bytes from file and read data */ 
   c = (unsigned char *)malloc(sizeof(unsigned char)* wordsize * narray);
   cnt = fread((void *)c, sizeof(unsigned char), narray * wordsize, bfile);
 
   fclose(bfile);

   if (cnt == 0) 
   {
      *status = 1;
      return 1;
   }

   /* 
      Set up byte offsets for each wordsize depending on byte order.
      A, B, C, D give the offsets of the LSB through MSB (i.e., for 
      word ABCD, A=MSB, D=LSB) in the array from the beginning of a word 
   */
   if (endian == GEOG_BIG_ENDIAN) {
      A2 = 0; B2 = 1;
      A3 = 0; B3 = 1; C3 = 2;
      A4 = 0; B4 = 1; C4 = 2; D4 = 3;
   }
   else {
      B2 = 0; A2 = 1;
      C3 = 0; B3 = 1; A3 = 2;
      D4 = 0; C4 = 1; B4 = 2; A4 = 3;
   }

   /* Convert words from native byte order */
   switch(wordsize) {
      case 1:
         for(i=0; i<narray; i++)
         {
            ival = (int)(c[i]);      
            if ((isigned) && (ival > (1 << 7))) ival -= (1 << 8);
            rarray[i] = (float)ival;
         }
         break;

      case 2:
         for(i=0; i<narray; i++)
         {
            ival = (int)((c[2*i+A2]<<8) | (c[2*i+B2]));      
            if ((isigned) && (ival > (1 << 15))) ival -= (1 << 16);
            rarray[i] = (float)ival;
         }
         break;

      case 3:
         for(i=0; i<narray; i++)
         {
            ival = (int)((c[3*i+A3]<<16) | (c[3*i+B3]<<8) | c[3*i+C3]);      
            if ((isigned) * (ival > (1 << 23))) ival -= (1 << 24);
            rarray[i] = (float)ival;
         }
         break;

      case 4:
         for(i=0; i<narray; i++)
         {
            ival = (int)((c[4*i+A4]<<24) | (c[4*i+B4]<<16) | (c[4*i+C4]<<8) | c[4*i+D4]);      
            if ((isigned) && (ival > (1 << 31))) ival = -(~ival + 1);
            rarray[i] = (float)ival;
         }
         break;
   }

   free(c);

   return 0;
}
