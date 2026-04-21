!======================================================================
subroutine pool_auto_tran( poolcont, &
     daylen, daylendt, daylenmax, &
     tc, pawfrw, rootf_lay, poollt, &
     gain_transl_lay, poollu_dgain)
!======================================================================

! Description
! ------------
! Calculates the autotrophic transfer from the live pools. 
!
use kinds
use module_param, only: pool_param
use module_poolinfo, only: & 
   pool_indx_lay, pool_indx_leaf
use module_sib, only: pooll_type
use module_sibconst, only: &
   npoolpft, npoollu, &
   ntpool, nsoil
use module_time, only: &
   dtisib, dtsib, steps_per_day

implicit none

!...input variables
type(pool_param), intent(in) :: poolcont
real(r4), intent(in) :: daylenmax
real(r8), intent(in) :: daylen, daylendt
real(r8), intent(in) :: tc, pawfrw
real(r8), dimension(nsoil), intent(in) :: rootf_lay
type(pooll_type), intent(inout) :: poollt
real(r8), dimension(npoollu,nsoil), intent(inout) :: &
    gain_transl_lay, poollu_dgain

!...local variables
integer(byte) :: lp
real(r8) :: qt, tfraclp
real(r8), dimension(npoolpft) :: tfrac !transfer fraction
real(r8), dimension(npoolpft,nsoil) :: &
   poolpft_avail_lay, & !available pool carbon (mol C/m2)
   tloss_lay   !pool transfer loss (mol C/m2/timestep)
real(r8), dimension(npoollu,nsoil) :: lutemp_gainl_lay

!...misc values
integer(i4) :: n,s
integer(i4) :: m,mref

!-----------------------------------------
!...Reset transfer variables
lp = pool_indx_leaf
tfrac(:) = dzero
poolpft_avail_lay(:,:) = dzero
tloss_lay(:,:) = dzero
lutemp_gainl_lay(:,:) = dzero
poollt%tfl_daylen = dzero
poollt%tfl_freeze = dzero
poollt%tfl_dry = dzero
poollt%tfl_total = dzero
poollt%tf_turnover(:) = dzero
poollt%loss_trans_lay(:,:) = dzero
gain_transl_lay(:,:) = dzero

!...Only transfer if there is pool carbon available
IF (sum(poollt%poolpft) .gt. sum(poolcont%poolpft_min)) THEN

do n=1,npoolpft
   do s=1,pool_indx_lay(n)
       poolpft_avail_lay(n,s) = poollt%poolpft_lay(n,s) &
         - poolcont%poolpft_min(n)  &
         + poollt%poolpft_dgain(n,s) &
         - poollt%poolpft_dloss(n,s)
   enddo
enddo

!-----------------------------------------------------------------
!-----Turnover Fractions------------------------------------------
!For leaves, this is supplemented by transfer from shortening days,
! freezing temperatures, and lack of water.
!For all other live pools, this is the only source of transfer.
do n=1, npoolpft
   do s=1,pool_indx_lay(n)
      tloss_lay(n,s) = (1 - poolcont%lresp_eff(n))  &
           * poollt%krater_lay(n,s) * poolpft_avail_lay(n,s) * dtsib

      if (poolpft_avail_lay(n,s) .gt. dzero) then
           tfrac(n) = tfrac(n) &
               + poollt%poolpft_flay(n,s) * tloss_lay(n,s) &
               / poolpft_avail_lay(n,s)
      endif
   enddo
   poollt%tf_turnover(n) = tfrac(n)*steps_per_day
enddo


!----Additional Leaf Pool Loss Fractions-----
!Daylength
if (daylendt .lt. dzero) then
    tfraclp = MIN(poolcont%lt_dmax, MAX(dzero, &
             poolcont%lt_dcoef * &
             (daylenmax - daylen)* &
             (daylenmax-poolcont%lt_dref)))
    tfrac(lp) = tfrac(lp) + tfraclp / steps_per_day
    poollt%tfl_daylen = tfraclp
else
   poollt%tfl_daylen = dzero    
endif

!Freezing
if (tc .lt. poolcont%lt_fref) then
    qt = 0.01 * (poolcont%lt_fref - tc)
    tfraclp = MIN(poolcont%lt_fmax, MAX(dzero, &
               poolcont%lt_fq10**qt - done))
    tfrac(lp) = tfrac(lp) + tfraclp / steps_per_day
    poollt%tfl_freeze = tfraclp
else
    poollt%tfl_freeze = dzero
endif

!Phenology stage
tfrac(lp) = tfrac(lp) + poollt%tfl_pstage / steps_per_day

!Water Deficiency
if (pawfrw .lt. poolcont%lt_wref) then
   tfraclp = MAX(dzero, MIN(poolcont%lt_wmax, &
        poolcont%lt_wcoef*(poolcont%lt_wbase &
        ** (10.*(pawfrw - poolcont%lt_wref)) - 1.)))
   tfrac(lp) = tfrac(lp) + tfraclp / steps_per_day
   poollt%tfl_dry = tfraclp
else
   poollt%tfl_dry = dzero
endif

!Combined factors
tfraclp = MIN(done, MAX(dzero, tfrac(lp)))
poollt%tfl_total = tfraclp * steps_per_day
tloss_lay(lp,1) = poolpft_avail_lay(lp,1) * tfraclp

!-----Live Pool Transfers------
do n=1,npoolpft
    do s=1,pool_indx_lay(n)
       poollt%loss_trans_lay(n,s) = tloss_lay(n,s) * dtisib
       poollt%poolpft_dloss(n,s) = &
              poollt%poolpft_dloss(n,s) + tloss_lay(n,s)
    enddo !s=1,pool_index_lay
enddo !n=1,npoolpft


!----Transfer To Dead Pools---
!...n is the sending/from pool
!...m is the receiving/to pool
tloss_lay = tloss_lay * dtisib !convert to mol C/m2/s
mref=1
do n=1,npoolpft
   do m=npoolpft+1,ntpool
      mref=m-npoolpft
      if (poolcont%pool_trans_frac(n,m) > dzero) then
          !.....transfer from single-layer canopy/soil to surface
          if ((pool_indx_lay(n) .eq. 1) .and. &
              (pool_indx_lay(m) .eq. 1)) then
                 lutemp_gainl_lay(mref,1) = &
                   lutemp_gainl_lay(mref,1) + &
                   tloss_lay(n,1) * poolcont%pool_trans_frac(n,m)

          !.....transfer from single-lay canopy/soil to soil
          elseif ((pool_indx_lay(n) .eq. 1) .and. &
                   (pool_indx_lay(m) .eq. nsoil)) then
                   do s=1,nsoil
                     lutemp_gainl_lay(mref,s) = &
                         lutemp_gainl_lay(mref,s) + &
                         tloss_lay(n,1) * rootf_lay(s) * &
                         poolcont%pool_trans_frac(n,m)
                  enddo

          !.....transfer from soil to soil
           elseif ((pool_indx_lay(n) .eq. nsoil) .and. &
                    (pool_indx_lay(m) .eq. nsoil)) then
                    do s=1,nsoil
                        lutemp_gainl_lay(mref,s) = &
                           lutemp_gainl_lay(mref,s) + &
                           tloss_lay(n,s) * &
                           poolcont%pool_trans_frac(n,m)
                    enddo
           else
                print*, 'Mismatching levels between pool transfers.'
                print*, 'Stopping in pool_trans_auto.'
                stop
            endif  !dead pool transfer
        endif !trans_frac > 0.
   enddo  !m=npoolpft+1,ntpool
enddo !n=1,npoolpft
ENDIF !pools > 0

!...Calculate dead pool gains
gain_transl_lay = lutemp_gainl_lay
poollu_dgain = poollu_dgain + lutemp_gainl_lay*dtsib

end subroutine pool_auto_tran
