!================SUBROUTINE COS=========================
!
!    Module to set pcosm using either bcosm (500 ppt) or
!    the TM5 mole fraction obtained through inversion that
!    covers seasonal, diurnal and spatial variability.
!
!    Choice for either one is based on cosvar_switch in namel_sibdrv
!    cosvar_switch = .false.  = 500 ppt
!    cosvar_switch = .true. = uses tm5 mixing ratio.
!
!
!-------------------------------------------------------

subroutine set_cos(gprogt)

    use module_pparams, only: p0_sfc, bcosm
    use module_sib, only: sib, gprog_type
    use module_sibconst, only: varcos_switch, subcount
    use module_phosib, only: pressure

    implicit none

    type(gprog_type), intent(inout) :: gprogt

    integer :: i, n ! iteration

    !pressure = dble(gprogt%ps) * 100.0

    if (varcos_switch) then
         do n = 1, subcount
         ! pcosm in pa, cosm_tm5 in ppt
         !gprogt%pcosm = (gprogt%cosm_tm5*p0_sfc)/1.E12
           sib%g(n)%gprogt%pcosm = dble(sib%g(n)%gprogt%cosm_tm5*pressure)/1.E12
         enddo
    else
         do n = 1, subcount
         ! pcosm in pa, bcosm in ppt
         !gprogt%pcosm = (bcosm*p0_sfc)/1.E12
           sib%g(n)%gprogt%pcosm = dble(bcosm*pressure)/1.E12
         enddo
    endif

end subroutine set_cos
