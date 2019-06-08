#define NULL_CHARACTER '\0'

/*
 use iso_c_binding, only : c_int, c_char

 interface
   subroutine c_pool_hash(hash, key) bind(c)
      use iso_c_binding, only : c_int, c_char
      integer (c_int), intent(inout) :: hash
      character (c_char), dimension(*), intent(in) :: key
   end subroutine c_pool_hash
 end interface
*/

void c_pool_hash(int* hash, char* key)
{
	int i;
	unsigned int whash;

	whash = 0;

	for (i=0; key[i] != NULL_CHARACTER; i++) {
		whash += (unsigned int)key[i];
	} 

	*hash = (int)(whash & 0x7fffffff);
}
