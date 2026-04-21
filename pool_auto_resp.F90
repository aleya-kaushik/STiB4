!======================================================================
subroutine pool_auto_resp( poolcont, &
     clim_assim, clim_lai, &
     assimd, lai, tc, td_lay, &
     rootf_lay, poollt, resp_soil, resp_soil_lay)
!======================================================================

! Description
! ------------
! Calculates the autotrophic respiration from the live pools. 
!
! Notes
! -------
! Studies suggest auto resp may scale predominantly 
! with assimilation (e.g. Flexas et al., 2006; Kirschbaum 1988;
!    Meir, 2008; Molchanov, 2009)
!
use kinds
use module_param, only: &
   pool_param
use module_poolinfo, only: & 
   pool_indx_froot, pool_indx_croot, &
   pool_indx_leaf, pool_indx_lay
use module_sib, only: &
   soil_type, pooll_type
use module_sibconst, only: &
   npoolpft, nsoil
use module_time, only: dtsib

implicit none


!...input variables
type(pool_param), intent(in) :: poolcont
real(r8), intent(in) :: clim_assim, clim_lai
real(r8), intent(in) :: assimd, lai, tc
real(r8), dimension(nsoil), intent(in) :: &
      rootf_lay, td_lay
type(pooll_type), intent(inout) :: poollt
real(r8), intent(inout) :: resp_soil
real(r8), dimension(nsoil), intent(inout) :: resp_soil_lay

!...local variables
real(r8) :: slope, yint
real(r8) :: qt, mtemp
real(r8) :: pool_valid, pool_updated
real(r8) :: temp_mrespr, temp_mrespl

!...misc values
integer(i4) :: n,s

!-----------------------------------------
!...Reset respiration variables
poollt%mcr_assim = done
poollt%mcr_freeze = done
poollt%mcr_hot = done
poollt%mcr_scale = done
poollt%mrr_freeze_lay(:) = done
poollt%mrr_hot_lay(:) = done
poollt%mrr_scale_lay(:) = done
poollt%mrr_assim = done
poollt%mrr_freeze = dzero
poollt%mrr_hot = dzero
poollt%mrr_lai = done
poollt%mrr_scale = dzero
poollt%krater_lay(:,:) = dzero
poollt%loss_mresp_lay(:,:) = dzero

!...Only respire if pools are greater than 
!...required minimum
IF (sum(poollt%poolpft) .GT. sum(poolcont%poolpft_min)) THEN

!-----------------------------------------------------------------
!-----Canopy Autotrophic Respiration Scaling Factors-----
!Assimilation Rate Scalar
if (assimd .lt. clim_assim*poolcont%cr_aml) then
   poollt%mcr_assim = poolcont%cr_amin
elseif (assimd .lt. clim_assim*poolcont%cr_amh) then
   slope = (poolcont%cr_amax - poolcont%cr_amin) &
       /(clim_assim*poolcont%cr_amh - clim_assim*poolcont%cr_aml)
   yint = poolcont%cr_amin - slope*clim_assim*poolcont%cr_aml
   poollt%mcr_assim = assimd*slope + yint
else
   poollt%mcr_assim = poolcont%cr_amax
endif

!Freeze Inhibition
poollt%mcr_freeze = MAX(poolcont%cr_fmin, MIN(1.0, &
           EXP(poolcont%cr_fmul*(tc - poolcont%cr_fref))))

!High Temperature Exponential
qt = 0.1 * (tc - poolcont%cr_href)
mtemp = poolcont%cr_hq10**qt
poollt%mcr_hot = MAX(done, &
    MIN(poolcont%cr_hmax, mtemp))

!Combined factors
poollt%mcr_scale = poollt%mcr_assim * poollt%mcr_freeze &
                 * poollt%mcr_hot

!----Soil Autotrophic Respiration Scaling Factors----
!Assimilation Rate Scalar
if (assimd .lt. clim_assim*poolcont%rrt_aml) then
   poollt%mrr_assim = poolcont%rrt_amin
elseif (assimd .lt. clim_assim*poolcont%rrt_amh) then
   slope = (poolcont%rrt_amax - poolcont%rrt_amin) &
       /(clim_assim*poolcont%rrt_amh - clim_assim*poolcont%rrt_aml)
   yint = poolcont%rrt_amin - slope*clim_assim*poolcont%rrt_aml
   poollt%mrr_assim = assimd*slope + yint
else
   poollt%mrr_assim = poolcont%rrt_amax
endif

!LAI
poollt%mrr_lai = MIN(poolcont%rrt_laimax, &
     MAX(poolcont%rrt_laimin, lai / clim_lai))

do s=1,nsoil
   !...Soil Freeze Inhibition
   poollt%mrr_freeze_lay = MAX(poolcont%rrt_fmin, MIN(1.0, &
        EXP(poolcont%rrt_fmul*(td_lay(s) - poolcont%rrt_fref))))
   poollt%mrr_freeze = poollt%mrr_freeze &
       + poollt%mrr_freeze_lay(s) * rootf_lay(s)

   !...Soil High Temp Exponential
    qt = 0.1 * (td_lay(s) - poolcont%rrt_href)
    mtemp = poolcont%rrt_hq10**qt
    poollt%mrr_hot_lay(s) = MAX(done, &
        MIN(poolcont%rrt_hmax, mtemp))
    poollt%mrr_hot = poollt%mrr_hot &
       + poollt%mrr_hot_lay(s) * rootf_lay(s)
enddo

!Combined Factors
poollt%mrr_scale_lay(:) = poollt%mrr_freeze_lay(:) &
   * poollt%mrr_hot_lay(:) * poollt%mrr_lai * poollt%mrr_assim

DO s=1,nsoil
   poollt%mrr_scale = poollt%mrr_scale  &
       + poollt%mrr_scale_lay(s) * rootf_lay(s)  
ENDDO


!-----------------------------------------------------------------
!-----Autotrophic Respiration Rate----
do n=1,npoolpft
    !...only resp if above min pool value
    pool_valid = poolcont%poolpft_min(n)
    pool_updated = poollt%poolpft(n) &
        + sum(poollt%poolpft_dgain(n,:)) &
        - sum(poollt%poolpft_dloss(n,:))
    IF (pool_updated .LE. pool_valid) CYCLE

    !...calculate loss rate
    if (pool_indx_lay(n) .eq. 1) then
        poollt%krater_lay(n,1) = poollt%mcr_scale &
              * poolcont%k_rate(n)
    else
        do s=1,pool_indx_lay(n)
           poollt%krater_lay(n,s) = poollt%mrr_scale_lay(s) &
               * poolcont%k_rate(n)
        enddo
    endif

    !.....calculate/check maintenance resp
    do s=1,pool_indx_lay(n)
       temp_mrespr = poollt%poolpft_lay(n,s) * &
                     poollt%krater_lay(n,s) * &
                     poolcont%lresp_eff(n)
       temp_mrespl = temp_mrespr*dtsib
       pool_updated = pool_updated - temp_mrespl
       IF (pool_updated .LT. pool_valid) CYCLE

       poollt%loss_mresp_lay(n,s) = temp_mrespr
       poollt%poolpft_dloss(n,s) = temp_mrespl &
           + poollt%poolpft_dloss(n,s)
    enddo !s=1,pool_index_lay
    
enddo !n=1,npoolpft
ENDIF !pools are greater than required minimum


!----Save respirations
poollt%resp_auto = sum(poollt%loss_gresp) + sum(poollt%loss_mresp_lay)
poollt%resp_leaf = poollt%loss_gresp(pool_indx_leaf) &
                   + sum(poollt%loss_mresp_lay(pool_indx_leaf,:))
poollt%resp_mntn = sum(poollt%loss_mresp_lay)
poollt%resp_root = &
    sum(poollt%loss_mresp_lay(pool_indx_froot,:)) + &
    sum(poollt%loss_mresp_lay(pool_indx_croot,:)) + &
    poollt%loss_gresp(pool_indx_froot) + &
    poollt%loss_gresp(pool_indx_croot)
resp_soil = poollt%resp_root
resp_soil_lay(:) = &
    poollt%loss_mresp_lay(pool_indx_froot,:) + &
    poollt%loss_mresp_lay(pool_indx_croot,:) + &
    poollt%loss_gresp(pool_indx_froot)*rootf_lay(:) + &
    poollt%loss_gresp(pool_indx_croot)*rootf_lay(:)

end subroutine pool_auto_resp
