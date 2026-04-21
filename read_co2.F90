
! opens and reads in CO2 timeseries
subroutine read_co2()
use module_io, only: co2_file, spatco2_file
use module_co2_timeseries, only: years,ppm_global,delta_co2
use module_sibconst, only : nsib,spatco2_switch
implicit none

logical :: exists
integer :: i,n_lines,value,year_ix,trendy,deltafile

n_lines = 0

open(newunit=trendy, file = co2_file, action='read')
do
    read(trendy,*,iostat=value) years,ppm_global
    if (value/=0) exit
    n_lines = n_lines+1
end do 
allocate(ppm_global(n_lines), years(n_lines))
rewind(trendy)
do i=1,n_lines
    read(trendy, *) years(i),ppm_global(i)
end do
close(trendy)

!n_lines = 0
!open(newunit=deltafile, file = spatco2_file, action='read')
!do
!    read(deltafile,*,iostat=value) delta_co2
!    if (value/=0) exit
!    n_lines = n_lines+1
!end do 
!
!if (n_lines .ne. nsib) then
!    if (spatco2_switch) then 
!        print*, 'FATAL ERROR: Invalid amount of nsibs in spatco2_file'
!        stop
!    endif
!endif
!
!allocate(delta_co2(12,n_lines))
!
!rewind(deltafile)
!do i=1,n_lines
!    read(deltafile, *) delta_co2(:,i)
!end do
!close(deltafile)



end subroutine read_co2
