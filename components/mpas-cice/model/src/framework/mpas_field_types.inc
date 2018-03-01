   integer, parameter :: MPAS_MISSING_DIM = -999


   integer, parameter :: MPAS_DECOMP_NONDECOMP = 1013, &
                         MPAS_DECOMP_CELLS     = 1014, &
                         MPAS_DECOMP_EDGES     = 1015, &
                         MPAS_DECOMP_VERTICES  = 1016

   ! Derived type for storing fields
   type field5DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:,:,:,:,:), pointer :: array

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(5) :: dimNames
      integer, dimension(5) :: dimSizes
      real (kind=RKIND) :: defaultValue
      real (kind=RKIND) :: missingValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_lists_type), dimension(:), pointer :: attLists => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field5DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field5DReal


   ! Derived type for storing fields
   type field4DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:,:,:,:), pointer :: array

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(4) :: dimNames
      integer, dimension(4) :: dimSizes
      real (kind=RKIND) :: defaultValue
      real (kind=RKIND) :: missingValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_lists_type), dimension(:), pointer :: attLists => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field4DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field4DReal



   ! Derived type for storing fields
   type field3DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:,:,:), pointer :: array

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(3) :: dimNames
      integer, dimension(3) :: dimSizes
      real (kind=RKIND) :: defaultValue
      real (kind=RKIND) :: missingValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_lists_type), dimension(:), pointer :: attLists => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field3DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field3DReal


   ! Derived type for storing fields
   type field2DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:,:), pointer :: array

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(2) :: dimNames
      integer, dimension(2) :: dimSizes
      real (kind=RKIND) :: defaultValue
      real (kind=RKIND) :: missingValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_lists_type), dimension(:), pointer :: attLists => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field2DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field2DReal


   ! Derived type for storing fields
   type field1DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:), pointer :: array

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(1) :: dimNames
      integer, dimension(1) :: dimSizes
      real (kind=RKIND) :: defaultValue
      real (kind=RKIND) :: missingValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_lists_type), dimension(:), pointer :: attLists => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field1DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field1DReal


   ! Derived type for storing fields
   type field0DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND) :: scalar

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      real (kind=RKIND) :: defaultValue
      real (kind=RKIND) :: missingValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      type (att_lists_type), dimension(:), pointer :: attLists => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field0DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field0DReal


   ! Derived type for storing fields
   type field3DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer, dimension(:,:,:), pointer :: array

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(3) :: dimNames
      integer :: defaultValue
      integer :: missingValue
      integer, dimension(3) :: dimSizes
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_lists_type), dimension(:), pointer :: attLists => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field3DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field3DInteger


   ! Derived type for storing fields
   type field2DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer, dimension(:,:), pointer :: array

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(2) :: dimNames
      integer :: defaultValue
      integer :: missingValue
      integer, dimension(2) :: dimSizes
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_lists_type), dimension(:), pointer :: attLists => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field2DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field2DInteger


   ! Derived type for storing fields
   type field1DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer, dimension(:), pointer :: array

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(1) :: dimNames
      integer :: defaultValue
      integer :: missingValue
      integer, dimension(1) :: dimSizes
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_lists_type), dimension(:), pointer :: attLists => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field1DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field1DInteger


   ! Derived type for storing fields
   type field0DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer :: scalar

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      integer :: defaultValue
      integer :: missingValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      type (att_lists_type), dimension(:), pointer :: attLists => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field0DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field0DInteger


   ! Derived type for storing fields
   type field1DChar
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      character (len=StrKIND), dimension(:), pointer :: array

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(1) :: dimNames
      integer, dimension(1) :: dimSizes
      character (len=StrKIND) :: defaultValue
      character (len=StrKIND) :: missingValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_lists_type), dimension(:), pointer :: attLists => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field1DChar), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field1DChar


   ! Derived type for storing fields
   type field0DChar
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      character (len=StrKIND) :: scalar

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND) :: defaultValue
      character (len=StrKIND) :: missingValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      type (att_lists_type), dimension(:), pointer :: attLists => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field0DChar), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field0DChar


   ! Derived type for storing fields
   type field0DLogical
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      logical :: scalar

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      logical :: defaultValue
      logical :: missingValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      type (att_lists_type), dimension(:), pointer :: attLists => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field0DLogical), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field0DLogical

