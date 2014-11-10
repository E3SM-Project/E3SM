#include <stdio.h>
#include <stdlib.h>

// Convert an index into a list of dimensions. E.g., for index 4 into a
// array defined as a[3][2], will return 1 1.
void idx_to_dim_list(const int ndims, const int gdims[], const int idx,
                     int dim_list[])
{
  int i, curr_idx, next_idx;
  curr_idx = idx;
  // Easiest to start from the right and move left.
  for (i = ndims-1; i >= 0; --i) {
    // This way of doing div/mod is slightly faster than using "/" and "%".
    next_idx = curr_idx / gdims[i];
    dim_list[i] = curr_idx - (next_idx*gdims[i]);
    curr_idx = next_idx;
  }
}

// Expand a region along dimension dim, by incrementing count[i] as much as
// possible, consistent with the map.
//
// Once max_size is reached, the map is exhausted, or the next entries fail
// to match, expand_region updates the count and calls itself with the next
// outermost dimension, until the region has been expanded as much as
// possible along all dimensions.
void expand_region(const int dim, const int gdims[], const int maplen,
                   const int map[], const int region_size,
                   const int region_stride, const int max_size[],
                   int count[])
{
  int i, j, test_idx, expansion_done;
  // Precondition: maplen >= region_size (thus loop runs at least once).

  // Flag used to signal that we can no longer expand the region along
  // dimension dim.
  expansion_done = 0;

  // Expand no greater than max_size along this dimension.
  for (i = 1; i <= max_size[dim]; ++i) {
    // Count so far is at least i.
    count[dim] = i;

    // Now see if we can expand to i+1 by checking that the next
    // region_size elements are ahead by exactly region_stride.
    // Assuming monotonicity in the map, we could skip this for the
    // innermost dimension, but it's necessary past that because the
    // region does not necessarily comprise contiguous values.
    for (j = 0; j < region_size; ++j) {
      test_idx = j + i*region_size;
      // If we have exhausted the map, or the map no longer matches,
      // we are done, break out of both loops.
      if (test_idx >= maplen || map[test_idx] != map[j] + i*region_stride) {
        expansion_done = 1;
        break;
      }
    }
    if (expansion_done) break;
  }

  // Move on to next outermost dimension if there are more left, else return.
  if (dim > 0) {
    expand_region(dim-1, gdims, maplen, map, region_size*count[dim],
                  region_stride*gdims[dim], max_size, count);
  }
}

// Set start and count so that they describe the first region in map.
void find_first_region(const int ndims, const int gdims[],
                       const int maplen, const int map[],
                       int start[], int count[])
{
  int dim;
  int max_size[ndims];
  // Preconditions (which might be useful to check/assert):
  //   ndims is > 0
  //   maplen is > 0
  //   all elements of map are inside the bounds specified by gdims

  idx_to_dim_list(ndims, gdims, map[0], start);

  for (dim = 0; dim < ndims; ++dim) {
    // Can't expand beyond the array edge.
    max_size[dim] = gdims[dim] - start[dim];
  }

  // For each dimension, figure out how far we can expand in that dimension
  // while staying contiguous in the input array.
  //
  // Start with the innermost dimension (ndims-1), and it will recurse
  // through to the outermost dimensions.
  expand_region(ndims-1, gdims, maplen, map, 1, 1, max_size, count);
}

// Find and print all regions. To do something other than print, this will
// have to output a linked list or other dynamic memory, since we don't
// know how many regions there will be on entry...
void get_start_and_count_regions(const int ndims, const int gdims[],
                                 const int maplen, const int map[])
{
  int dim;
  int start[ndims];
  int count[ndims];
  int nmaplen;
  int regionlen;

  nmaplen = 0;

  while(nmaplen < maplen){
    find_first_region(ndims, gdims, maplen-nmaplen, map+nmaplen, start, count);
    printf("start %d %d %d\n", start[0], start[1], start[2]);
    printf("count %d %d %d\n", count[0], count[1], count[2]);
    regionlen = 1;
    for(dim=0; dim<ndims; dim++){
      regionlen *= count[dim];
    }
    nmaplen = nmaplen+regionlen;
  }

}




int main()
{
    int *iomap;
    int ndims=3;
    int gdims[ndims];
    int maplen=38;
    iomap = (int *) calloc(maplen,sizeof(int));
    
    gdims[0]=2;
    gdims[1]=15;
    gdims[2]=12;

    // Old data has four regions:
    //  start 0 0 9 count 1 4 3
    //  start 0 4 0 count 1 4 3
    //  start 0 8 3 count 1 4 3
    //  start 0 12 2 count 2 1 1
    //
    // New data has:
    //  start 0 0 9 count 1 3 3
    //  start 0 3 0 count 1 1 3
    //  start 0 3 9 count 1 1 3
    //  start 0 4 0 count 1 3 3
    //  start 0 8 3 count 1 4 3
    //  start 0 12 2 count 2 1 1

    iomap[0] = 9 ; // 0 0 9
    iomap[1] = 10 ; // 0 0 10
    iomap[2] = 11 ; // 0 0 11
    iomap[3] = 21 ; // 0 1 9
    iomap[4] = 22 ; // 0 1 10
    iomap[5] = 23 ; // 0 1 11
    iomap[6] = 33 ; // 0 2 9
    iomap[7] = 34 ; // 0 2 10
    iomap[8] = 35 ; // 0 2 11
#ifdef NEW_DATA
    iomap[9] = 36 ; // 0 3 0
    iomap[10] = 37 ; // 0 3 1
    iomap[11] = 38 ; // 0 3 2
    iomap[12] = 45 ; // 0 3 9
    iomap[13] = 46 ; // 0 3 10
    iomap[14] = 47 ; // 0 3 11
    iomap[15] = 48 ; // 0 4 0
    iomap[16] = 49 ; // 0 4 1
    iomap[17] = 50 ; // 0 4 2
    iomap[18] = 60 ; // 0 5 0
    iomap[19] = 61 ; // 0 5 1
    iomap[20] = 62 ; // 0 5 2
    iomap[21] = 72 ; // 0 6 0
    iomap[22] = 73 ; // 0 6 1
    iomap[23] = 74 ; // 0 6 2
#else
    iomap[9] = 45 ; // 0 3 9
    iomap[10] = 46 ; // 0 3 10
    iomap[11] = 47 ; // 0 3 11
    iomap[12] = 48 ; // 0 4 0
    iomap[13] = 49 ; // 0 4 1
    iomap[14] = 50 ; // 0 4 2
    iomap[15] = 60 ; // 0 5 0
    iomap[16] = 61 ; // 0 5 1
    iomap[17] = 62 ; // 0 5 2
    iomap[18] = 72 ; // 0 6 0
    iomap[19] = 73 ; // 0 6 1
    iomap[20] = 74 ; // 0 6 2
    iomap[21] = 84 ; // 0 7 0
    iomap[22] = 85 ; // 0 7 1
    iomap[23] = 86 ; // 0 7 2
#endif
    iomap[24] = 99 ; // 0 8 3
    iomap[25] = 100 ; // 0 8 4
    iomap[26] = 101 ; // 0 8 5
    iomap[27] = 111 ; // 0 9 3
    iomap[28] = 112 ; // 0 9 4
    iomap[29] = 113 ; // 0 9 5
    iomap[30] = 123 ; // 0 10 3
    iomap[31] = 124 ; // 0 10 4
    iomap[32] = 125 ; // 0 10 5
    iomap[33] = 135 ; // 0 11 3
    iomap[34] = 136 ; // 0 11 4
    iomap[35] = 137 ; // 0 11 5 
    iomap[36] = 146 ; // 0 12 2
    iomap[37] = 326 ; // 1 12 2


    get_start_and_count_regions(ndims, gdims, maplen, iomap);

    free(iomap);

    return 0;
}
