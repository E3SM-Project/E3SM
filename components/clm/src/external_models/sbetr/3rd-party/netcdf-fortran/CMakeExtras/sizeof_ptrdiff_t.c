#include <stddef.h>
#include <stdio.h>
int main ()
{
  switch(sizeof(ptrdiff_t))
  {
    case 4:
      printf("%1d",4);
      break;
    case 8:
      printf("%1d",8);
      break;
    default:
      return 1;
  }
  return 0;
}
