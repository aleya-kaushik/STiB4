module module_poolinfo

!----------------------------------------------------------------------
!
!   SiB4 Pool Informational Module
!
!----------------------------------------------------------------------

use kinds
implicit none


  character(len=20), dimension(:), allocatable :: pool_name_long ! pool names
  character(len=6), dimension(:), allocatable ::  pool_name      ! short pool names
  character(len=4), dimension(:), allocatable ::  pool_type      ! pool types (live, dead)
  character(len=7), dimension(:), allocatable ::  pool_loc       ! pool locations (soil, surface, canopy)
  integer(byte), dimension(:), allocatable ::     pool_indx_lay  ! index for lowest soil layer with carbon

  !...Pool indices
  integer(byte), dimension(:), allocatable :: &
       pool_indx_can,  &  ! pool index numbers for all canopy pools
       pool_indx_sfc,  &  ! pool index numbers for all surface pools
       pool_indx_soil     ! pool index numbers for all soil pools

  integer(byte) pool_indx_leaf    ! pool index number for leaf pool
  integer(byte) pool_indx_froot   ! pool index number for fine root pool
  integer(byte) pool_indx_croot   ! pool index number for coarse root pool
  integer(byte) pool_indx_stwd    ! pool index number for stem/wood pool
  integer(byte) pool_indx_prod    ! pool index number for product pool
  integer(byte) pool_indx_cdb     ! pool index number for coarse dead biomass pool
  integer(byte) pool_indx_metl    ! pool index number for metabolic litter pool
  integer(byte) pool_indx_strl    ! pool index number for structural litter pool
  integer(byte) pool_indx_slit    ! pool index number for soil litter pool
  integer(byte) pool_indx_slow    ! pool index number for soil slow pool
  integer(byte) pool_indx_arm     ! pool index number for soil armored pool

end module module_poolinfo

