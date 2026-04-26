subroutine grid_update()

use module_sibconst, only: nsib,sublarea,subpref,subset,single_pt,subcount, rzero
use module_sibvs
use module_io, only: pbp_larea,pbp_pref, pbp_dtsib, npbp
use module_sib, only: sib

implicit none
integer(i4) :: i,l
if (single_pt) then
       sublarea(1,:) = sibvs(subset(1))%larea(:)
       subpref(1,:) = sibvs(subset(1))%pftref(:)
       if (pbp_dtsib .ne. 0) then
              pbp_larea(1,:) = sibvs(1)%larea(:)
       endif 
else
       if (pbp_dtsib .ne. 0 .AND. subcount .ne. npbp) then
              print *, 'ERROR'
              print *, 'Amount of Points simulated:',subcount,'does not match points written out to pbp:',npbp
              print *, 'For grid_update to work please select npbp = -1'
              stop
       endif

       do i=1, subcount
              sublarea(i,:) = sibvs(subset(i))%larea(:)
              subpref(i,:) = sibvs(subset(i))%pftref(:)
              sib%g(i)%l(:)%larea = rzero
              do l=1,sibvs(subset(i))%gnlu
                     sib%g(i)%l(l)%larea = sibvs(subset(i))%larea(l)
              enddo
              if (pbp_dtsib .ne. 0 .AND. npbp == subcount) then
                     pbp_larea(i,:) = rzero
                     do l=1,sibvs(subset(i))%gnlu
                            pbp_larea(i,l) = sibvs(subset(i))%larea(l)
                            pbp_pref(i,l) = sibvs(subset(i))%pftref(l)
                     enddo
              endif
       enddo
endif 
end subroutine grid_update
