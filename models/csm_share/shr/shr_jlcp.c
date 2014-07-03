#if (defined CATAMOUNT)
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/stat.h>
#include <string.h>

/*
** This emulates the cp command being called from system()
*/
int shr_jlcp(char *from, char *to)
{
  char buf[1000];
  int hasread, fd1, fd2, wrote, written=0, totalw=0,totalr=0;

  fd1 = open(from, O_RDONLY);
  fd2 = open(to, O_CREAT|O_TRUNC|O_WRONLY, S_IREAD|S_IWRITE);

  /*DEBUG printf("JL: Copying %s to %s.\n", from, to);*/
  if ( fd1 < 0 || fd2 < 0 )
  {
    perror("JLCP: unable to open file for copy\n");
    if (fd1 > 0)
      close(fd1);
    if (fd2 > 0)
      close(fd2);
    return -1;
  }

  while ( (hasread = read(fd1,buf,1000)) != 0 )
  {
    if ( hasread == EINTR || hasread == EAGAIN )
      continue;
    else if ( hasread < 0 )
    {
      perror("JLCP: unable to read from file");
      close(fd1);
      close(fd2);
      return hasread;
    } else
    {
      written = 0;
      /*
      ** write() does not guarantee that all will be written.
      ** This loop handles the rare case that write writes
      ** less than the requested number of bytes.
      */
      while ( written < hasread )
      {
        wrote = write(fd2,&(buf[written]),hasread);
        if ( wrote == EINTR || wrote == EAGAIN )
          continue;
        else if ( wrote < 0 )
        {
          perror("JLCP: unable to write to file");
          close(fd1);
          close(fd2);
          return -1;
        }
        else
          written += wrote;
      }
      totalw += written;
      totalr += hasread;
    }
  }

  /*DEBUG printf("JL: read %d wrote: %d\n", totalr, totalw);*/
  close(fd1);
  close(fd2);
  return 0;
} 

/*
** Because Fortan and C do not represent strings the same way, this
** function requires string lengths so that C strings can be made
** and passed to jlcp.  
*/
void shr_jlcp_(char *from, int *fromlen, char *to, int *tolen, int *iret)
{
  char *from_,*to_;

  /* create C strings */
  from_ = malloc(*fromlen+1);
  if ( from_ == NULL )
  {
    *iret = -1;
    perror("JLCP: Could not allocate from_ string");
    return;
  }
  to_ = malloc(*tolen+1);
  if ( to_ == NULL )
  {
    *iret = -1;
    perror("JLCP: Could not allocate to_ string");
    free(from_);
    return;
  }
  
  /* copy strings and add NULL char to the end. */
  bcopy(from,from_,*fromlen);
  from_[*fromlen] = '\0';
  bcopy(to,to_,*tolen);
  to_[*tolen] = '\0';
  
  *iret = shr_jlcp(from_,to_);
  
  /* free the temporary strings */
  free(from_);
  free(to_);
}
#else
/*
  The following is only here since empty file won't compile
*/
#include <stdio.h>
void shr_jlcp_stub_()
{
    printf("JLCP: This stub should NOT be called");
    return;
}
#endif
