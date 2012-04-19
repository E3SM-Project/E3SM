#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h> // posix io 
#include <unistd.h>    // posix io 
#include <sys/stat.h>  // posix io
#include <fcntl.h>     // posix io

#include <errno.h>

#define MAX_FILES 99
static off_t rl[MAX_FILES];         // record length
// map io handle to fd
static int handlemap[MAX_FILES] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
				   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
				   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
				   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
				   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
				   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
				   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
				   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
				   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
				   -1,-1,-1,-1,-1,-1,-1,-1,-1};
static char *fnmap[MAX_FILES];      // file names
static int nfiles = 0;              // number of files
static int verbose = 0;

#if ( defined FORTRANUNDERSCORE )
#define jropen_direct jropen_direct_
#define jrclose_direct jrclose_direct_
#define jrwrite_direct jrwrite_direct_
#define jrread_direct jrread_direct_
#endif

static inline int gethandle (void);
static int jrioerr (const char *, ...);

// Mimic Fortran direct access open
int jropen_direct (char *fn,         // file to open
		   char *faccess,    // "old" (R/O) or "replace" (WRONLY)
		   int *reclength,   // record length (bytes)
		   int *handle,      // output handle
		   int *ret,         // output return code
		   int nc1, int nc2) // number of chars in Fortran variables
{
  char *access;    // Fortran access

  int flags;       // flags to open file: read-only or write-only
  int fd;          // file descriptor from open
  int lochandle;

  mode_t mode;

  *ret = -1;
  if (*reclength < 1)
    return jrioerr ("jropen_direct: bad reclength=%d bytes\n", *reclength);

  if ( ! (access = malloc (nc2+1)))
    return jrioerr ("jropen_direct: malloc failure %d bytes\n", nc2+1);

  strncpy (access, faccess, nc2);
  access[nc2] = '\0';

  if (nfiles > MAX_FILES-1)
    return jrioerr ("jropen_direct: %d files is too many\n", nfiles);

  // O_DIRECT bypasses system buffering and goes directly to disk

  if (nc2 > 6 && strncmp (access, "replace", 7) == 0) {
    flags = O_CREAT | O_WRONLY;
    mode = S_IRUSR | S_IWUSR;
  } else if (nc2 > 2 && strncmp (access, "old", 3) == 0) {
    flags = O_RDONLY;
    mode = S_IRUSR;
  } else {
    return jrioerr ("jropen_direct: access capability %s unsupported\n", access);
  }

  if ((lochandle = gethandle ()) < 0)
    return jrioerr ("jropen_direct: empty handle not found\n");

  if ( ! (fnmap[lochandle] = malloc (nc1+1)))
    return jrioerr ("jropen_direct: malloc fnmap failure %d bytes\n", nc1+1);

  strncpy (fnmap[lochandle], fn, nc1);
  fnmap[lochandle][nc1] = '\0';

  // Open the file

  if ((fd = open (fnmap[lochandle], flags, mode)) < 0)
    return jrioerr ("jropen_direct: failure to open file %s:%s\n", fnmap[lochandle], strerror(errno));

  if (*reclength < 1)
    return jrioerr ("jropen_direct: bad reclength=%d\n", *reclength);

  handlemap[lochandle] = fd;       // map array index to file descriptor
  rl[lochandle] = *reclength;      // set map value
  *handle = lochandle;             // make handle available to user
  ++nfiles;                        // number of open files

  if (verbose) {
    printf ("jropen_direct: opened direct-access file %s access=%s\n", fnmap[lochandle], access);
    printf ("handle is %d\n", *handle);
    printf ("File descriptor is %d\n", fd);
    printf ("Record length is %d\n", *reclength);
  }

  free (access);
  return (*ret = 0);
}

// Mimic Fortran direct access close
int jrclose_direct (const int *handle)
{
  int i;

  if (*handle < 0 || *handle > MAX_FILES-1)
    return jrioerr ("jrclose_direct: bad handle=%d\n", *handle);

  if (handlemap[*handle] < 0)
    return jrioerr ("jrclose_direct: handle %d not currently assigned\n", *handle);

  if (close (handlemap[*handle]) < 0)
    return jrioerr ("jrclose_direct: fclose error for handle %d:%s\n", *handle, strerror (errno));

  if (verbose) {
    printf ("jrclose_direct: closed direct-access file %s\n", fnmap[*handle]);
    printf ("handle returned to pool is %d\n", *handle);
    printf ("File descriptor returned to pool is %d\n", handlemap[*handle]);
  }
  handlemap[*handle] = -1;  // mark the handle available for future use
  free (fnmap[*handle]);
  --nfiles;
	    
  return 0;
}

// Mimic Fortran direct access write
int jrwrite_direct (const int *handle, const int *recnum, const char *arr, int *bswap, int *ret)
{
  int rc;           // return from "write"
  int fd;           // file descriptor to write to
  int byteswritten; // counts bytes written up to rl[*handle]
  int bytestowrite; // number of bytes from arr to write in the next system call
  int n;

  char *locarr;     // local array in case byte swapping needed

  *ret = -1; // Set return code in case jrioerr gets called

  if (*handle < 0 || *handle > MAX_FILES-1)
    return jrioerr ("jrwrite_direct: bad handle=%d\n", *handle);

  if ((fd = handlemap[*handle]) < 0)
    return jrioerr ("jrwrite_direct: handle %d not currently assigned\n", *handle);

  if (*recnum < 1)
    return jrioerr ("jrwrite_direct: non-positive recnum=%d encountered\n", *recnum);

  if (*bswap != 0 && *bswap != 4 &&*bswap != 8)
    return jrioerr ("jrwrite_direct: bswap must be 0, 4, or 8. Got %d\n", *bswap);

  if (lseek (fd, (*recnum - 1)*rl[*handle], SEEK_SET) < 0)
    return jrioerr ("jrwrite_direct: handle=%d seek failure recnum=%d bytepos=%d\n", 
		    *handle, *recnum, (*recnum-1)*rl[*handle]);

  if (*bswap == 0) {
    for (byteswritten = 0; byteswritten < rl[*handle]; byteswritten += rc) {
      bytestowrite = rl[*handle] - byteswritten;
      if ((rc = write (fd, &arr[byteswritten], bytestowrite)) < 0)
	return jrioerr ("jrwrite_direct: failure writing rec %d of file %s (handle %d fd %d):%s\n", 
			*recnum, fnmap[*handle], *handle, fd, strerror (errno));
      
      if (rc == 0 && byteswritten < rl[*handle])
	return jrioerr ("jrwrite_direct: only able to write %d bytes of rec %d file %s\n", 
			byteswritten, *recnum, fnmap[*handle]);
    }

  } else {

    // Swap bytes before writing

    if ( ! (locarr = malloc (rl[*handle])))
      return jrioerr ("jrwrite_direct: malloc error\n");

    if (*bswap == 4) {
      for (n = 0; n < rl[*handle]; n+=4) {
	locarr[n+3] = arr[n];
	locarr[n+2] = arr[n+1];
	locarr[n+1] = arr[n+2];
	locarr[n]   = arr[n+3];
      }
    } else if (*bswap == 8) {
      for (n = 0; n < rl[*handle]; n+=8) {
	locarr[n+7] = arr[n];
	locarr[n+6] = arr[n+1];
	locarr[n+5] = arr[n+2];
	locarr[n+4] = arr[n+3];
	locarr[n+3] = arr[n+4];
	locarr[n+2] = arr[n+5];
	locarr[n+1] = arr[n+6];
	locarr[n]   = arr[n+7];
      }
    }

    for (byteswritten = 0; byteswritten < rl[*handle]; byteswritten += rc) {
      bytestowrite = rl[*handle] - byteswritten;
      if ((rc = write (fd, &locarr[byteswritten], bytestowrite)) < 0) {
	free (locarr);
	return jrioerr ("jrwrite_direct: failure writing rec %d of file %s (handle %d fd %d):%s\n", 
			*recnum, fnmap[*handle], *handle, fd, strerror (errno));
      }
      
      if (rc == 0 && byteswritten < rl[*handle]) {
	free (locarr);
	return jrioerr ("jrwrite_direct: only able to write %d bytes of rec %d file %s\n", 
			byteswritten, *recnum, fnmap[*handle]);
      }
    }
    free (locarr);
  }

  // *ret = byteswritten;
  return (*ret = 0);
}

// Mimic Fortran direct access read
int jrread_direct (const int *handle, const int *recnum, char *arr, int *bswap, int *ret)
{
  int rc;           // return from "read"
  int fd;           // file descriptor to read to
  int bytesread;    // counts bytes read up to rl[*handle]
  int bytestoread;  // number of bytes from arr to read in the next system call
  int n;

  char locarr[8];   // for byte swapping

  *ret = -1;

  if (verbose) {
    printf ("starting jrread_direct: handle=%d, file descriptor=%d, recnum=%d\n", 
	    *handle, handlemap[*handle], *recnum);
  }

  if (*handle < 0 || *handle > MAX_FILES-1)
    return jrioerr ("jrread_direct: bad handle=%d\n", *handle);

  if (*bswap != 0 && *bswap != 4 &&*bswap != 8)
    return jrioerr ("jrread_direct: bswap must be 0, 4, or 8. Got %d\n", *bswap);

  if ((fd = handlemap[*handle]) < 0) {
    printf ("jrread_direct: fd=%d\n", fd);
    return jrioerr ("jrread_direct: handle %d not currently assigned\n", *handle);
  }

  if (*recnum < 1)
    return jrioerr ("jrread_direct: non-positive recnum=%d encountered\n", *recnum);

  if (lseek (fd, (*recnum - 1)*rl[*handle], SEEK_SET) < 0)
    return jrioerr ("jrread_direct: handle=%d seek failure recnum=%d bytepos=%d\n", 
		    *handle, *recnum, (*recnum-1)*rl[*handle]);
    
  for (bytesread = 0; bytesread < rl[*handle]; bytesread += rc) {
    bytestoread = rl[*handle] - bytesread;
    if ((rc = read (fd, &arr[bytesread], bytestoread)) < 0)
      return jrioerr ("jrread_direct: failure reading rec %d of file %s (handle %d fd %d):%s\n", 
		      *recnum, fnmap[*handle], *handle, fd, strerror (errno));

    if (rc == 0 && bytesread < rl[*handle])
      return jrioerr ("jrread_direct: only able to read %d bytes of rec %d file %s\n", 
		      bytesread, *recnum, fnmap[*handle]);
  }

  // Swap bytes if necessary

  if (*bswap == 4) {
    for (n = 0; n < rl[*handle]; n+=4) {
      locarr[3] = arr[n];
      locarr[2] = arr[n+1];
      locarr[1] = arr[n+2];
      locarr[0] = arr[n+3];

      arr[n  ] = locarr[0];
      arr[n+1] = locarr[1];
      arr[n+2] = locarr[2];
      arr[n+3] = locarr[3];
    }
  } else if (*bswap == 8) {
    for (n = 0; n < rl[*handle]; n+=8) {
      locarr[7] = arr[n];
      locarr[6] = arr[n+1];
      locarr[5] = arr[n+2];
      locarr[4] = arr[n+3];
      locarr[3] = arr[n+4];
      locarr[2] = arr[n+5];
      locarr[1] = arr[n+6];
      locarr[0] = arr[n+7];

      arr[n]   = locarr[0];
      arr[n+1] = locarr[1];
      arr[n+2] = locarr[2];
      arr[n+3] = locarr[3];
      arr[n+4] = locarr[4];
      arr[n+5] = locarr[5];
      arr[n+6] = locarr[6];
      arr[n+7] = locarr[7];
    }
  }

  //  *ret = bytesread;
  return (*ret = 0);
}

static int jrioerr (const char *fmt, ...)
{
  va_list args;
  
  va_start (args, fmt);
  
  if (fmt != NULL)
    (void) vfprintf (stderr, fmt, args);
  
  va_end (args);
  
  return (-1);
}

static inline int gethandle (void)
{
  int i;

  for (i = 0; i < MAX_FILES; ++i)
    if (handlemap[i] < 0)
      return i;

  return jrioerr ("gethandle: no empty handles found\n");
}
