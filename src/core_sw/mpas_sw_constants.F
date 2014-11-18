! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module sw_constants

   use mpas_grid_types
   use mpas_constants

   public
   save

   type (mpas_pool_type), pointer :: swConfigs, swPackages


   contains


   subroutine sw_constants_init(configPool, packagePool)
      type (mpas_pool_type), pointer :: configPool
      type (mpas_pool_type), pointer :: packagePool

      swConfigs => configPool
      swPackages => swPackages

   end subroutine sw_constants_init



end module sw_constants
