
! Reads in first year of LUH2
subroutine luc_init(year)
use module_io, only: luh2_path
use module_luc, only: luc
use module_sibconst, only: luc_movec_switch
!---------------------------------------------------
!...input variables
integer :: year   ! year
!---------------------------------------------------

if (luc_movec_switch) then
    print*, 'Reading in first year of LUH2:',year
endif



end subroutine luc_init