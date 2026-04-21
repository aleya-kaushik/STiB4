!==========================================================================
subroutine pool_graze( &
   poolpft_min, grz_transfer, &
   clim_lai, lai, &
   poollt, pooldt)
!==========================================================================
! ----------------------------------------------------------------
! Calculates grazing, removing carbon from live above-ground pools 
!    and transferring it to dead pools
! ------------------------------------------------------------

use kinds
use module_oparams, only: &
   graze_cfracp, graze_cfracd, &
   graze_minlai, graze_climlai
use module_pparams, only: mol_to_umol
use module_poolinfo, only: &
   pool_indx_can, pool_indx_lay
use module_param, only: pool_param
use module_sib, only: &
   poold_type, pooll_type
use module_sibconst, only: &
   npoolcan, npoolpft, npoollu, &
   ntpool
use module_time, only: &
   dtsib, dtisib, wt_daily

implicit none

! Input variables
real(r8), intent(in) :: clim_lai, lai
real(r4), dimension(npoolpft), intent(in) :: poolpft_min
real(r4), dimension(npoollu+2), intent(in) :: grz_transfer
type(poold_type), intent(inout) :: pooldt
type(pooll_type), intent(inout) :: poollt

! Local variables
integer(i4) :: cp, p
integer(i4) :: m, mref, s
real(r8) :: cpooltot, pool_avail
real(r8) :: grzf, grzc


!----------------------------------------------
!Reset variables
poollt%resp_grz = dzero
poollt%loss_grz = dzero
pooldt%gain_grz_lay = dzero

!Check to make sure canopy/leaf pools are large
!enough to support grazing
cpooltot = sum(poollt%poolpft(pool_indx_can))
pool_avail = sum(poolpft_min(pool_indx_can))
IF ((lai .le. graze_minlai) .or. &
    (cpooltot .lt. pool_avail)) RETURN

!Set grazing rate
IF (clim_lai .gt. graze_climlai) THEN
   grzf = graze_cfracp * wt_daily * dtisib
ELSE
   grzf = graze_cfracd * wt_daily * dtisib
ENDIF

!Increment number of grazing days
poollt%nd_grz = poollt%nd_grz + wt_daily

!----------------------------------------------
!Calculate grazing loss and live removal
DO p=1, npoolcan
   cp = pool_indx_can(p)
   pool_avail = poollt%poolpft(cp) - poolpft_min(cp)
   poollt%loss_grz(p) = MIN(pool_avail, &
        poollt%poolpft(cp) * grzf)
   poollt%poolpft_dloss(cp,1) = &
       poollt%poolpft_dloss(cp,1) + poollt%loss_grz(p)*dtsib
ENDDO

!Determine where grazed carbon loss goes
!   (respired, removed, transferred to dead)
grzc = sum(poollt%loss_grz)
poollt%resp_grz = grzc * grz_transfer(1)
do m=npoolpft+1,ntpool
   mref=m-npoolpft
   if (grz_transfer(mref+2) .gt. dzero) then
       do s=1,pool_indx_lay(m)
           pooldt%gain_grz_lay(mref,s) = grzc &
                 * grz_transfer(mref+2) &
                 * pooldt%poollu_flay(mref,s)
       enddo
       pooldt%poollu_dgain(mref,:) = pooldt%poollu_dgain(mref,:) &
                 + pooldt%gain_grz_lay(mref,:)*dtsib
    endif
enddo


end subroutine pool_graze
