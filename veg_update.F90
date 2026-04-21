subroutine veg_update( doy, &
    gref, glon, glat, &
    pnum, pref, &
    iscrop, isgrass, &
    physcont, &
    poollu, snow_cvfc, node_z, &
    poollt, vegt)
!==========================================================================
! ------------------------------------------------------------
! Updates vegetation state
! ------------------------------------------------------------

use kinds
use module_io, only: &
    ztemp, zwind
use module_param, only: &
    phys_param
use module_poolinfo
use module_pparams, only: &
    drytoc, mwc
use module_sib, only: &
    pooll_type, veg_type
use module_sibconst, only: &
    npoolpft, npoollu, nsoil, &
    green_switch, &
    print_stop, print_veg
use module_time, only: wtd_clim
implicit none

!...input variables
integer(i4), intent(in) :: doy, gref, pnum, pref
real(r4), intent(in) :: glon, glat
logical, intent(in) :: iscrop, isgrass
type(phys_param), intent(in) :: physcont

real(r8), dimension(npoollu), intent(in) :: poollu
real(r8), intent(in) :: snow_cvfc
real(r8), dimension(nsoil), intent(in) :: node_z
type(pooll_type), intent(in) :: poollt
type(veg_type), intent(inout) :: vegt

!...local variables
integer(i4) :: s
real(r8) :: lai_pool, lait_pool
real(r8) :: totagb, totgreen

! ------------------------------------------------------------
! Calculate LAI
IF (isgrass) THEN
   lai_pool = poollt%poolpft(pool_indx_leaf) &
              + poollt%poolpft(pool_indx_stwd)
   lait_pool = sum(poollt%poolpft(pool_indx_can)) &
                + poollu(pool_indx_cdb-npoolpft)
ELSE
   lai_pool = poollt%poolpft(pool_indx_leaf)
   lait_pool = lai_pool
ENDIF

vegt%lai = MAX(physcont%laimin, lai_pool &
           * mwc * drytoc * physcont%sla)

vegt%lait = MAX(physcont%laimin, lait_pool &
           * mwc * drytoc * physcont%sla)

vegt%clim_lai = (1.-wtd_clim)*vegt%clim_lai &
           + wtd_clim*vegt%lai

! Calculate FPAR
vegt%fpar = (1.0 - exp(max(min(vegt%lai, &
        max(physcont%laisat, 0.001)),0.0) * &
        log(max(1.0-physcont%fparsat,0.001)) / &
        max(physcont%laisat, 0.001)))

vegt%vcover = (1.0 - exp(max(min(vegt%lait, &
        max(physcont%laisat, 0.001)),0.0) * &
        log(max(1.0-physcont%fparsat,0.001)) / &
        max(physcont%laisat, 0.001)))

!...Original FPAR formula from Sellers et al. (1996)
!vegt%fpar = MAX(.001, & 
!          (1. - exp(-vegt%park * vegt%lai)))
!vegt%vcover = vegt%fpar

! Update PFT aerodynamic properties
call AeroInterpolate( &
        gref, glon, glat, pnum, pref, &
        vegt%lai, vegt%vcover, &
        vegt%z0d, vegt%zp_dispd, vegt%cc1, vegt%cc2)

! Calculate gmudmu
call gmuder(glat, doy, physcont%chil, vegt%gmudmu)

! Update above-ground greenness fraction
if (.not. green_switch) then
    vegt%green = done
elseif (iscrop .or. isgrass) then
    totagb = sum(poollt%poolpft(pool_indx_can))
    totgreen = totagb

    totagb = totagb &
        + poollu(pool_indx_cdb-npoolpft) &
        + poollu(pool_indx_metl-npoolpft) &
        + poollu(pool_indx_strl-npoolpft)
    if (totagb .gt. 0.) then
        vegt%green = totgreen / totagb
    else
        vegt%green = dzero
    endif
else
    vegt%green = done
endif
vegt%green = MAX(dzero, MIN(done, vegt%green))

! Update canopy characteristics
vegt%zpd_adj = vegt%zp_dispd + &
        (physcont%z2 - vegt%zp_dispd) * snow_cvfc
vegt%z0 = max(0.1, vegt%z0d / (physcont%z2 - vegt%zp_dispd) * &
        (physcont%z2 - vegt%zpd_adj))

vegt%zztemp = physcont%z2 - vegt%zpd_adj + ztemp
vegt%zzwind = physcont%z2 - vegt%zpd_adj + zwind

! -------------------
! Print out info
if (print_veg) then
    print*,''
    print('(a,3f10.3)'), &
        '      LAI/FPAR/Green: ', &
        vegt%lai,vegt%fpar,vegt%green
    print('(a,2f10.3)')   , &
        '      z0/zpd_adj:     ', &
        vegt%z0, vegt%zpd_adj
    print('(a,2f10.3)')   , &
        '      cc1/cc2:        ', &
        vegt%cc1, vegt%cc2
    print('(a,g12.4)')    , &
        '      vmax:           ', &
        vegt%vmax
    print'(2a)','        lev   layer (m)    rootf    ', &
                ' frp (mol C/m2)   crp (mol C/m2)'
    do s=1,nsoil
       print'(a,i6,1f10.5,f12.5,a,f12.5,a,f12.5)','     ', &
            s, node_z(s), &
            vegt%rootf(s), '  ', &
            poollt%poolpft_lay(pool_indx_froot,s),'    ', &
            poollt%poolpft_lay(pool_indx_croot,s)
    enddo

    if (print_stop) stop
endif


end subroutine veg_update
