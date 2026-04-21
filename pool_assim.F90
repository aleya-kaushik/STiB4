!==========================================================================
subroutine pool_assim( &
    sibpt, lonsib, latsib, pref, &
    hhti, hlti, shti, slti, gr_frac, &
    assim, rstfac2, tm, poollt)
!==========================================================================
! -----------
! Updates live carbon pools from photosynthetic gains
! -----------

use kinds
use module_poolinfo, only: pool_indx_lay
use module_sib, only: pooll_type
use module_sibconst, only: npoolpft
use module_time, only: dtsib

implicit none

!...input variables
integer(i4), intent(in) :: sibpt, pref
real(r4), intent(in) :: lonsib, latsib
real(r4), intent(in) :: hhti, hlti, shti, slti
real(r4), dimension(npoolpft), intent(in) :: gr_frac
real(r8), intent(in) :: assim, rstfac2, tm
type(pooll_type), intent(inout) :: poollt

!...local variables
integer(i4) :: k, p
integer(i4) :: ialloc
real(r8) :: deltac

!--------------------------------------------------------------------

!Reset variables
poollt%resp_grow = dzero
poollt%resp_nveg = dzero
poollt%gain_assim = dzero
poollt%loss_gresp = dzero

!Check for assimilation to allocate
IF (assim .le. dzero) RETURN

!Check phenology stage allocation fractions,
!...if not one, then in dormancy and respire
!...back out the 'assimilated' carbon
ialloc = nint(sum(poollt%alloc_phen))
IF (ialloc .ne. ione) THEN
    poollt%resp_nveg = assim
    RETURN
ENDIF

!Calculate allocation fraction adjustments
IF ((poollt%aadj_moist) .or. (poollt%aadj_temp)) THEN
    call pool_alloc( &
         sibpt, lonsib, latsib, pref, &
         hhti, hlti, shti, slti, &
         rstfac2, tm, &
         poollt%aadj_moist, poollt%aadj_temp, &
         poollt%alloc_phen, poollt%alloc_moist, &
         poollt%alloc_temp, poollt%alloc)
ELSE
    poollt%alloc_moist(:) = dzero
    poollt%alloc_temp(:) = dzero
    poollt%alloc(:) = poollt%alloc_phen
ENDIF

!Assign photosynthate to live pools
do p=1, npoolpft
   deltac = assim*poollt%alloc(p)
   poollt%gain_assim(p) = deltac
   poollt%loss_gresp(p) = deltac*gr_frac(p)
   do k=1,pool_indx_lay(p)
      poollt%poolpft_dgain(p,k) = poollt%poolpft_dgain(p,k) &
             + poollt%gain_assim(p)*poollt%poolpft_flay(p,k)*dtsib
      poollt%poolpft_dloss(p,k) = poollt%poolpft_dloss(p,k) &
             + poollt%loss_gresp(p)*poollt%poolpft_flay(p,k)*dtsib
   enddo
enddo
poollt%resp_grow = sum(poollt%loss_gresp)

end subroutine pool_assim


!==========================================================================
subroutine pool_alloc( &
    sibpt, lonsib, latsib, pref, &
    hhti, hlti, shti, slti,  &
    rstfac2, tm, &
    adj_moist, adj_temp, &
    alloc_phen, alloc_moist, alloc_temp, &
    alloc)
!==========================================================================
! -----------
!Calculates allocation fractions
!  -Phenology-based allocation from phenology stage
!  -Climate-based adjustments from Friedlingstein et al. (1999)
! -----------

use kinds
use module_oparams, only: &
    aadjustmin, &
    lftit, lftif, lgrw_min, &
    moistmult, wftit, wftif
use module_poolinfo
use module_sibconst, only: &
    npoolpft

implicit none

!...input variables
integer(i4), intent(in) :: sibpt, pref
real(r4), intent(in) :: lonsib, latsib
real(r4), intent(in) :: hhti, hlti, shti, slti
real(r8), intent(in) :: rstfac2, tm
logical, intent(in) :: adj_moist, adj_temp
real(r8), dimension(npoolpft), intent(in) :: alloc_phen
real(r8), dimension(npoolpft), intent(inout) :: &
    alloc_moist, alloc_temp, alloc

!...local variables
integer(i4) :: lp,frp,crp,wp,pp
integer(i4) :: npallow
real(r8) :: atot, aadjust, aadd
real(r8), dimension(npoolpft) :: alloc_now
logical, dimension(npoolpft) :: alloc_allow

real(r8) :: leaf_grw_frz  !scaling factor for leaf growth from temperature
                          !..leaf growth decreases when leaves cold/frozen
real(r8) :: wood_grw_frz  !scaling factor for wood growth from temperature
                          !..no wood growth for freezing temperatures
real(r8) :: wood_grw_temp !scaling factor for wood growth from temperature
                          !..wood growth decreases for low temperatures
real(r8) :: wood_grw_tot  !total wood growth temp scaling factor
real(r8) :: wood_grw_moist !scaling factor for wood growth from moisture
                           !..wood growth decreases with moisture stress

!...misc variables
integer(i4) :: p

!--------------------------------------------------------------------

!Set local variables
lp=pool_indx_leaf
frp=pool_indx_froot
crp=pool_indx_croot
wp=pool_indx_stwd
pp=pool_indx_prod


!Reset allocations
alloc_moist(:) = dzero
alloc_temp(:) = dzero
alloc(:) = dzero
alloc_now(:) = alloc_phen(:)
alloc_allow(:) = .false.

!--------------------------------
!-----ALLOCATION ADJUSTMENTS-----
!--------------------------------

!...scale the factors for temperature stress
IF ((adj_temp) .AND. &
    ((alloc_now(lp) .GT. 0.) .OR. (alloc_now(wp) .GT. 0.))) THEN
   !.....leaf growth decline for cold temperatures
   leaf_grw_frz = 1./(1.+EXP(lftif*(lftit-tm)))
   leaf_grw_frz = MAX(lgrw_min,leaf_grw_frz)

   !.....wood growth decline for cool temperatures
   wood_grw_temp = 1./(1.+EXP(slti*(hlti-tm))) / &
                   (1.+EXP(shti*(tm-hhti)))
   wood_grw_frz = 1./(1.+EXP(wftif*(wftit-tm)))
   wood_grw_tot = wood_grw_temp * wood_grw_frz

   aadjust = MAX(0.,alloc_now(wp)*(1.-wood_grw_tot))  + &
                  MAX(0.,alloc_now(crp)*(1.-wood_grw_tot)) + &
                  MAX(0.,alloc_now(lp)*(1.-leaf_grw_frz))

   if (aadjust >= aadjustmin) then
      alloc_temp(lp) = -1. * MAX(0.,alloc_now(lp)*(1.-leaf_grw_frz))
      alloc_now(lp) = alloc_now(lp) + alloc_temp(lp)
      alloc_temp(wp) = -1. * MAX(0.,alloc_now(wp)*(1-wood_grw_tot))
      alloc_now(wp) = alloc_now(wp) + alloc_temp(wp)
      alloc_temp(crp)= -1. * MAX(0.,alloc_now(crp)*(1-wood_grw_tot))
      alloc_now(crp) = alloc_now(crp) + alloc_temp(crp)

      npallow=0
      do p=1,npoolpft
         IF ((p .ne. lp) .and. (p .ne. wp) .and. &
             (alloc_now(p) .gt. 0.)) THEN
             alloc_allow(p) = .true.
             npallow=npallow+1
         ELSE
             alloc_allow(p) = .false.
         ENDIF
      enddo

      if (npallow == 0) then
         !adjustment fine roots
         alloc_temp(frp) = aadjust
         alloc_now(frp) = aadjust
      else
         aadd = aadjust/dble(npallow)
         do p=1,npoolpft
            IF (alloc_allow(p)) THEN
               alloc_temp(p) = &
                  alloc_temp(p) + aadd
              alloc_now(p) = alloc_now(p) + aadd
            ENDIF
        enddo
      endif

      !check totals
      atot = sum(alloc_temp)
      if ((atot .lt. -0.01) .or. (atot .gt. 0.01)) then
           print*,'---Error with temperature adjustment allocation factors---'
           print*,'   SiB Point/Lon/Lat/PFT: ', sibpt, lonsib, latsib, pref
           print*,'   Temperature Adjustment Total Allocation Factor: ',atot
           print*,'   Temperature Adjustment Contributions: '
           print*, alloc_temp
           stop
       endif

       atot = sum(alloc_now)
       if ((atot .lt. 0.99) .or. (atot .gt. 1.01)) then
           print*,'---Error with total allocation following temp adjustments---'
           print*,'   SiB Point/Lon/Lat/PFT: ',sibpt, lonsib, latsib, pref
           print*,'   Total Allocation Factor: ',atot
           print*,'   Total Allocation Contributions: '
           print*, alloc_now(:)
           stop
       endif

   endif !(aadjust > adjustmin)
ENDIF !use_temp

!----------------
IF ((adj_moist) .AND. &
    ((alloc_now(lp) .GT. 0.) .OR. &
     (alloc_now(wp) .GT. 0.) .OR. & 
     (alloc_now(pp) .GT. 0.))) THEN

   !...scale the factors for moisture stress
   wood_grw_moist = MIN(1.0,MAX(0.,1. - rstfac2*moistmult))
   aadjust = MIN(1.0,MAX(0., &
       alloc_now(lp)*wood_grw_moist + &
       alloc_now(wp)*wood_grw_moist + &
       alloc_now(pp)*wood_grw_moist))

   if (aadjust >= aadjustmin) then
      alloc_moist(lp) = -1. * alloc_now(lp)*wood_grw_moist
      alloc_now(lp) = alloc_now(lp) + alloc_moist(lp)

      alloc_moist(wp) = -1. * alloc_now(wp)*wood_grw_moist
      alloc_now(wp) = alloc_now(wp) + alloc_moist(wp)

      alloc_moist(pp) = -1. * alloc_now(pp)*wood_grw_moist
      alloc_now(pp) = alloc_now(pp) + alloc_moist(pp)

      !put in fine roots
      alloc_moist(frp) = aadjust
       alloc_now(frp) = alloc_now(frp) + alloc_moist(frp)

       !check totals
       atot = sum(alloc_moist)
       if ((atot .lt. -0.01) .or. (atot .gt. 0.01)) then
            print*,'---Error with moisutre adjustment allocation factors---'
            print*,'   SiB Point/Lon/Lat/PFT: ', sibpt, lonsib, latsib, pref
            print*,'   Moisture Adjustment Total Allocation Factor: ',atot
            print*,'   Moisture Adjustment Contributions: '
            print*, alloc_moist
            stop
        endif

        atot = sum(alloc_now)
        if ((atot .lt. 0.99) .or. (atot .gt. 1.01)) then
            print*,'---Error with total allocation following moisture adjustments---'
            print*,'   SiB Point/Lon/Lat/PFT: ', sibpt, lonsib, latsib, pref
            print*,'   Total Allocation Factor: ',atot
            print*,'   Total Allocation Contributions: '
            print*, alloc_now(:)
            stop
        endif

   endif  !aadjust > aadjustmin

ENDIF !use_moist

!-----------------------------
!---SET ALLOCATION FACTORS---
alloc(:) = alloc_now(:)

!...check total allocation
atot = sum(alloc)
if ((atot .lt. 0.99) .or. (atot .gt. 1.01)) then
   print*,'---Error with final allocation factors---'
   print*,'   SiB Point/Lon/Lat/PFT: ', sibpt, lonsib, latsib, pref
   print*,'   Total Allocation Factor: ',atot
   print*,'   Allocation Contributions: '
   print*, alloc(:)
   stop
endif

!...check for negative allocations
do p=1,npoolpft
   if (alloc(p) < 0.) then
      if (alloc(p) > -.001) then
         alloc(p) = 0.
      else
          print*,'---Negative final allocation factors---'
          print*,'   SiB Point/Lon/Lat/PFT: ', sibpt, lonsib, latsib, pref
          print*,'   Total Allocation Factor: ',atot
          print*,'   Allocation Contributions: '
          print*, alloc(:)
          stop
      endif
   endif
enddo

end subroutine pool_alloc




