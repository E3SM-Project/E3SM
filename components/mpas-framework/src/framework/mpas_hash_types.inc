   integer, parameter :: TABLESIZE=27183     !< Number of spaces in the table (the number of linked lists)

   type hashnode
      integer :: key
      type (hashnode), pointer :: next
   end type hashnode
 
   type hashnode_ptr
      type (hashnode), pointer :: p        !< Pointer to a list of entries
   end type hashnode_ptr
 
   type hashtable
      integer :: size
      type (hashnode_ptr), dimension(TABLESIZE) :: table !< The hashtable array
   end type hashtable
