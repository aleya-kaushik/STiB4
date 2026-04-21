!!Subroutine to update the carbon pools
!!
!!===========================================
subroutine pool_update( &
    sibpt, siblon, siblat, pref, &
    assimin, laiin, fparin, &
    equibdt, pooldt, equiblt, poollt)
!!==========================================

use kinds
use module_sib, only: &
   equibd_type, poold_type, &
   equibl_type, pooll_type
use module_sibconst, only: &
   npoollu, npoolpft

implicit none

!...input values
integer(i4), intent(in) :: sibpt, pref
real(r4), intent(in) :: siblon, siblat
real(r8), intent(in) :: assimin, laiin, fparin

type(equibd_type), intent(inout) :: equibdt
type(equibl_type), intent(inout) :: equiblt
type(poold_type), intent(inout) :: pooldt
type(pooll_type), intent(inout) :: poollt

!...local variables
integer(i4) :: p

!------------------------------------------------
!Update pools and equilibrium total gains/losses
!Check carbon balance
!Reset variables

!...dead pools
pooldt%poollu_lay = pooldt%poollu_lay &
   + pooldt%poollu_dgain - pooldt%poollu_dloss

do p=1,npoollu
   equibdt%poollu_totgain(p) = equibdt%poollu_totgain(p) &
       + sum(pooldt%poollu_dgain(p,:))
   equibdt%poollu_totloss(p) = equibdt%poollu_totloss(p) &
       + sum(pooldt%poollu_dloss(p,:))

   pooldt%poollu(p) = sum(pooldt%poollu_lay(p,:))
   equibdt%poollu_max(p) = MAX(equibdt%poollu_max(p), &
       pooldt%poollu(p))
   equibdt%poollu_min(p) = MIN(equibdt%poollu_min(p), &
       pooldt%poollu(p))
enddo

!...live pools
poollt%poolpft_lay = poollt%poolpft_lay &
   + poollt%poolpft_dgain - poollt%poolpft_dloss

do p=1,npoolpft
   equiblt%poolpft_totgain(p) = equiblt%poolpft_totgain(p) &
       + sum(poollt%poolpft_dgain(p,:))
   equiblt%poolpft_totloss(p) = equiblt%poolpft_totloss(p) &
       + sum(poollt%poolpft_dloss(p,:))

   poollt%poolpft(p) = sum(poollt%poolpft_lay(p,:))
   equiblt%poolpft_max(p) = MAX(equiblt%poolpft_max(p), &
       poollt%poolpft(p))
   equiblt%poolpft_min(p) = MIN(equiblt%poolpft_min(p), &
       poollt%poolpft(p))
enddo

!...check carbon balance
call balan_carbon(sibpt, siblon, siblat, pref, &
     assimin, laiin, fparin, pooldt, poollt)

!...reset
pooldt%poollu_dgain = dzero
pooldt%poollu_dloss = dzero
poollt%poolpft_dgain = dzero
poollt%poolpft_dloss = dzero


end subroutine pool_update
