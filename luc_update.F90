
! Reads in specific year of LUH2
subroutine luc_update(year)

use module_io, only: luh2_path
use module_luc, only: luc
use module_sibconst, only: luc_movec_switch

!---------------------------------------------------
!...input variables
integer :: year   ! year
!---------------------------------------------------
character(4) :: syear
logical :: exists

write(syear,1000) year
1000 format (I4)
close(1000)

if (luc_movec_switch) then
    inquire(file=trim(luh2_path)//trim(syear)//'.nc',exist=exists)
    if (exists) then
        print*, 'Reading in ',year,' of LUH2'
    else 
        print*, 'No such file',trim(luh2_path)//trim(syear)//'.nc'
    endif 

endif
end subroutine luc_update