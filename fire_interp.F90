!--------------------------------------------------------------
subroutine fire_interp(indx, lon, lat, &
           sibg)

!--------------------------------------------------------------
!
! This subroutine interpolates the sibdrv fire emissions
!     between their read times
!
!--------------------------------------------------------------

use kinds
use module_io, only: &
   fire_step, fire_seccur, fire_secnext
use module_oparams, only: &
    fire_leaff, fire_stwdf, &
    fire_cdbf, fire_metlf, fire_strlf
use module_param, only: poolcon
use module_pparams, only: &
     mwc, mol_to_umol, &
     month_names!, pdb
use module_pftinfo
use module_poolinfo
use module_sib, only: gridcell_type
use module_sibconst, only: &
   npoolpft, npoollu, nsoil, &
   fireb_print, fireb_stop, fireb_thresh
use module_time, only: &
   month, day, year, &
   dtisib, dtsib, sec_tot, wt_daily

implicit none

!...input variables
integer(i4), intent(in) :: indx
real(r4), intent(in) :: lon, lat
type (gridcell_type), intent(inout) :: sibg

!...parameters
integer(i4), parameter :: isave=5
real(r8) :: dnzero=1.E-10

!...interpolation variables
real(r8) :: facsibdrv  ! scaling factor between driver data points
real(r8) :: totemis, curemis, pcemis

!...distribution variables
integer(i4) :: ntpft
integer(i4), dimension(:), allocatable :: tpref, tpnum
real(r4), dimension(:), allocatable :: tparea
real(r8), dimension(:), allocatable :: tpagb, tpagbtemp
real(r8), dimension(:), allocatable :: tpbiomass
character(len=clen), dimension(:), allocatable :: tpname

real(r8) :: tfarea
integer(i4), dimension(:), allocatable :: tsortref
real(r8), dimension(:), allocatable :: flossb
real(r8), dimension(:,:), allocatable :: flosspft, flosslu

!...net balance variables
logical :: fb_err

!...misc variables
integer :: ido, idoc, l, p, s, myl
integer :: loc_c3g, loc_c4g
real(r8) :: myemis, tempemis
integer(i4), dimension(1) :: tempstore
integer(i4) :: lp,frp,crp,wp,pp
integer(i4) :: cdbp, metlp, strlp, slitp, slowp, armp

!-------Calculatae quasi-equlibrium pools------------
lp =  pool_indx_leaf
frp = pool_indx_froot
crp = pool_indx_croot
wp =  pool_indx_stwd
pp =  pool_indx_prod

cdbp  = pool_indx_cdb-npoolpft
metlp = pool_indx_metl-npoolpft
strlp = pool_indx_strl-npoolpft
slitp = pool_indx_slit-npoolpft
slowp = pool_indx_slow-npoolpft
armp  = pool_indx_arm-npoolpft

!-----------------------------------------------
! reset values
do l=1, sibg%g_nlu
   sibg%l(l)%poollt%resp_fire = dzero
   sibg%l(l)%poollt%rmmd_fire = dzero
   sibg%l(l)%poollt%loss_fire_lay(:,:) = dzero
   sibg%l(l)%pooldt%loss_fire_lay(:,:) = dzero
enddo

! only continue if fire emissions are being used
IF (fire_step .le. izero) RETURN

! get scaling factors
facsibdrv = dble(fire_seccur-sec_tot) / dble(fire_step)

! only continue if fire emissions are valid at this time
IF ((facsibdrv .GT. 1) .OR. (facsibdrv .LT. 0)) THEN
   sibg%gprogt%firec = dzero
   sibg%gprogt%fireco2 = dzero
   RETURN
ENDIF

! interpolate fire C
sibg%gprogt%firec = facsibdrv*sibg%gprogt%firec1 + &
                    (1.-facsibdrv) * sibg%gprogt%firec2

! interpolate fire CO2
sibg%gprogt%fireco2 = facsibdrv*sibg%gprogt%fireco21 &
         + (1.-facsibdrv) * sibg%gprogt%fireco22

! distribute fire emissions per PFT/land unit
if (sibg%gprogt%firec .gt. dzero) then
   totemis = sibg%gprogt%firec * dtsib
   ntpft = sibg%g_nlu
   ! number of land units/PFTs per cell -> =1 for Hyy when run as 1.0ENF
   allocate(tpref(ntpft),tpnum(ntpft),tparea(ntpft))
   allocate(tpagb(ntpft),tpbiomass(ntpft))
   allocate(tpname(ntpft))
   tpref(:) = sibg%l(1:ntpft)%ipft
   tpnum(:) = pft_num(tpref)
   do l=1,ntpft
     tpname(l) = pft_name(tpnum(l))
   enddo
   tparea(:) = sibg%l(1:ntpft)%larea
   tpagb(:) = dzero
   tpbiomass(:) = dzero

   do l=1, ntpft
     !...calculate total above-ground biomass (m-2)
      do p=1, npoolpft
         if (pool_indx_lay(p) .eq. 1) then
            tpagb(l) = tpagb(l) + tparea(l) &
                 * (sibg%l(l)%poollt%poolpft(p) &
                 - sum(sibg%l(l)%poollt%poolpft_dloss(p,:)) &
                 - poolcon(tpnum(l))%poolpft_min(p))
         endif
      enddo

      !...calculate total biomass (m-2)
      tpbiomass(l) = tpbiomass(l) + sum(sibg%l(l)%poollt%poolpft(1:5)) &
           - sum(sibg%l(l)%poollt%poolpft_dloss(1:5,:)) &
           - sum(poolcon(tpnum(l))%poolpft_min(1:5))
      tpbiomass(l) = tpbiomass(l) + sum(sibg%l(l)%pooldt%poollu(1:6)) &
           - sum(sibg%l(l)%pooldt%poollu_dloss(1:6,:))
   enddo

   !...rank PFTs by above ground biomass
   allocate(tpagbtemp(ntpft),tsortref(ntpft))
   tpagbtemp = tpagb

   ! case: both c3g and c4g present
   if ( (ANY(tpname .eq. "c3g")) .and. &
        (ANY(tpname .eq. "c4g")) ) then
     do l=1,ntpft
       if (tpname(l) .eq. "c3g") loc_c3g = l
       if (tpname(l) .eq. "c4g") loc_c4g = l
     enddo

     ! sort grasses to the top, c4g first
     tsortref(1) = loc_c4g
     tsortref(2) = loc_c3g
     tpagbtemp(loc_c4g) = -1
     tpagbtemp(loc_c3g) = -1
     do l=3,ntpft
         tempstore = maxloc(tpagbtemp)
         tsortref(l) = tempstore(1)
         tpagbtemp(tsortref(l)) = -1
     enddo
   endif ! end case selection: c3g & c4g


   ! case: c3g only
   if ( (ANY(tpname .eq. "c3g")) .and. &
        (ALL(tpname .ne. "c4g")) ) then
     do l=1,ntpft
       if (tpname(l) .eq. "c3g") loc_c3g = l
     enddo
     tsortref(1) = loc_c3g
     tpagbtemp(loc_c3g) = -1
     do l=2,ntpft
         tempstore = maxloc(tpagbtemp)
         tsortref(l) = tempstore(1)
         tpagbtemp(tsortref(l)) = -1
     enddo
   endif! end case selection: c3g

   ! case: c4g only
   if ( (ALL(tpname .ne. "c3g")) .and. &
        (ANY(tpname .eq. "c4g")) ) then
     do l=1,ntpft
       if (tpname(l) .eq. "c4g") loc_c4g = l
     enddo
     tsortref(1) = loc_c4g
     tpagbtemp(loc_c4g) = -1
     do l=2,ntpft
         tempstore = maxloc(tpagbtemp)
         tsortref(l) = tempstore(1)
         tpagbtemp(tsortref(l)) = -1
     enddo
   endif! end case selection: c4g

   ! case: no grasses, default to original scheme
   if ( (ALL(tpname .ne. "c3g")) .and. &
        (ALL(tpname .ne. "c4g")) ) then
     do l=1,ntpft
         tempstore = maxloc(tpagbtemp)
         tsortref(l) = tempstore(1)
         tpagbtemp(tsortref(l)) = -1
     enddo
   endif! end case selection: no grasses


   !...check area of top PFTs
   ido=MIN(isave,ntpft)
   idoc=ido
   do l=1,ido
      myl = tsortref(l)
      if (tpagb(myl) .lt. totemis) then
         idoc=MAX(ione, idoc-1)
      endif
   enddo
   ido=idoc
   tfarea = SUM(tparea(tsortref(1:ido)))

   !...remove carbon from top PFTs
   allocate(flosspft(ntpft,npoolpft))
   allocate(flosslu(ntpft,npoollu))
   flosspft(:,:) = dzero
   flosslu(:,:) = dzero
   allocate(flossb(ntpft))
   flossb(:) = dzero

   do l=ido,1,-1
      myl = tsortref(l)
      sibg%l(myl)%poollt%nd_fire = sibg%l(myl)%poollt%nd_fire + wt_daily

      tempemis = totemis*(tparea(myl)/tfarea)
      curemis = tempemis

      !....remove C from leaf pool
      myemis = MIN(fire_leaff*tempemis, MAX(dzero, &
           sibg%l(myl)%poollt%poolpft(lp) &
             - sibg%l(myl)%poollt%poolpft_dloss(lp,1) &
             - poolcon(tpnum(myl))%poolpft_min(lp)))
      flosspft(myl,lp) = myemis
      sibg%l(myl)%poollt%loss_fire_lay(lp,1) = myemis * dtisib
      sibg%l(myl)%poollt%poolpft_dloss(lp,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(lp,1) &
           + myemis
      curemis = curemis - myemis

      !....remove C from wood pool
      myemis = MIN(fire_stwdf*tempemis, MAX(dzero, &
           sibg%l(myl)%poollt%poolpft(wp) &
             - sibg%l(myl)%poollt%poolpft_dloss(wp,1) &
             - poolcon(tpnum(myl))%poolpft_min(wp)))
      sibg%l(myl)%poollt%loss_fire_lay(wp,1) = myemis * dtisib
      sibg%l(myl)%poollt%poolpft_dloss(wp,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(wp,1) &
           + myemis
      flosspft(myl,wp) = myemis
      curemis = curemis - myemis

      !....remove C from metabolic litter
      myemis = MIN(fire_metlf*tempemis, MAX(dzero, &
           sibg%l(myl)%pooldt%poollu(metlp) &
           - sibg%l(myl)%pooldt%poollu_dloss(metlp,1)))
      sibg%l(myl)%pooldt%loss_fire_lay(metlp,1) = myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(metlp,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(metlp,1) &
           + myemis
      flosslu(myl,metlp) = myemis
      curemis = curemis - myemis

      !....remove C from structural litter
      myemis = MIN(fire_strlf*tempemis, MAX(dzero, &
           sibg%l(myl)%pooldt%poollu(strlp) &
           - sibg%l(myl)%pooldt%poollu_dloss(strlp,1)))
      sibg%l(myl)%pooldt%loss_fire_lay(strlp,1) = myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(strlp,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(strlp,1) &
           + myemis
      flosslu(myl,strlp) = myemis
      curemis = curemis - myemis

      !....remove C from coarse dead biomass
      myemis = MIN(fire_cdbf*tempemis, MAX(dzero, &
           sibg%l(myl)%pooldt%poollu(cdbp) &
           - sibg%l(myl)%pooldt%poollu_dloss(cdbp,1)))
      sibg%l(myl)%pooldt%loss_fire_lay(cdbp,1) = myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(cdbp,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(cdbp,1) &
           + myemis
      flosslu(myl,cdbp) = myemis
      curemis = curemis - myemis

      !...remove C from product
      ! IF (curemis .gt. dnzero) THEN
      ! 
      !    !tmppooltot = tmppooltot + sibg%l(myl)%poollt%curpoolpft(pp)
      !
      !    myemis = MIN(curemis, MAX(dzero, &
      !         sibg%l(myl)%poollt%poolpft(pp) &
      !         - sibg%l(myl)%poollt%poolpft_dloss(pp,1) &
      !         - poolcon(tpnum(myl))%poolpft_min(pp)))
      !    sibg%l(myl)%poollt%loss_fire_lay(pp,1) = myemis * dtisib
      !    sibg%l(myl)%poollt%poolpft_dloss(pp,1) = &
      !         sibg%l(myl)%poollt%poolpft_dloss(pp,1) &
      !         + myemis
      !    flosspft(myl,pp) = myemis
      !    curemis = curemis - myemis
      !
      ! ENDIF
      !
      ! !...remove C from roots
      ! IF (curemis .gt. dnzero) THEN
      !
      !    !tmppooltot = tmppooltot + sibg%l(myl)%poollt%curpoolpft(frp)
      !
      !     myemis = MIN(curemis, MAX(dzero, &
      !          sibg%l(myl)%poollt%poolpft(frp) &
      !          - sum(sibg%l(myl)%poollt%poolpft_dloss(frp,:)) &
      !          - poolcon(tpnum(myl))%poolpft_min(frp)))
      !
      !     DO s=1,nsoil
      !        sibg%l(myl)%poollt%loss_fire_lay(frp,s) = &
      !             myemis * dtisib * sibg%l(myl)%poollt%poolpft_flay(frp,s)
      !        sibg%l(myl)%poollt%poolpft_dloss(frp,s) = &
      !             sibg%l(myl)%poollt%poolpft_dloss(frp,s) &
      !             + myemis * sibg%l(myl)%poollt%poolpft_flay(frp,s)
      !     ENDDO
      !     flosspft(myl,frp) = myemis
      !     curemis = curemis - myemis
      !  ENDIF

      !  IF (curemis .gt. dnzero) THEN
      !
      !    !tmppooltot = tmppooltot + sibg%l(myl)%poollt%curpoolpft(crp)
      !
      !      myemis = MIN(curemis, MAX(dzero, &
      !          sibg%l(myl)%poollt%poolpft(crp) &
      !          - sum(sibg%l(myl)%poollt%poolpft_dloss(crp,:)) &
      !          - poolcon(tpnum(myl))%poolpft_min(crp)))
      !
      !     DO s=1,nsoil
      !        sibg%l(myl)%poollt%loss_fire_lay(crp,s) = &
      !             myemis * dtisib * sibg%l(myl)%poollt%poolpft_flay(crp,s)
      !        sibg%l(myl)%poollt%poolpft_dloss(crp,s) = &
      !             sibg%l(myl)%poollt%poolpft_dloss(crp,s) &
      !             + myemis * sibg%l(myl)%poollt%poolpft_flay(crp,s)
      !     ENDDO
      !     flosspft(myl,crp) = myemis
      !     curemis = curemis - myemis
      !  ENDIF
      !
      !  !...remove C from soil litter
      !  IF (curemis .gt. dnzero) THEN
      !
      !    !tmppooltot = tmppooltot + sibg%l(myl)%pooldt%curpoollu(slitp)
      !
      !      myemis = MIN(curemis, MAX(dzero, &
      !           sibg%l(myl)%pooldt%poollu(slitp) &
      !           - sum(sibg%l(myl)%pooldt%poollu_dloss(slitp,:))))
      !
      !      DO s=1,nsoil
      !        sibg%l(myl)%pooldt%loss_fire_lay(slitp,s) = &
      !             myemis * dtisib * sibg%l(myl)%pooldt%poollu_flay(slitp,s)
      !        sibg%l(myl)%pooldt%poollu_dloss(slitp,s) = &
      !             sibg%l(myl)%pooldt%poollu_dloss(slitp,s) &
      !             + myemis * sibg%l(myl)%pooldt%poollu_flay(slitp,s)
      !     ENDDO
      !     flosslu(myl,slitp) = myemis
      !     curemis = curemis - myemis
      ! ENDIF
      !
      ! !...remove C from soil slow
      ! IF (curemis .gt. dnzero) THEN
      !
      !    !tmppooltot = tmppooltot + sibg%l(myl)%pooldt%curpoollu(slowp)
      !
      !      myemis = MIN(curemis, MAX(dzero, &
      !          sibg%l(myl)%pooldt%poollu(slowp) &
      !          - sum(sibg%l(myl)%pooldt%poollu_dloss(slowp,:))))
      !
      !     DO s=1,nsoil
      !        sibg%l(myl)%pooldt%loss_fire_lay(slowp,s) = &
      !             myemis * dtisib * sibg%l(myl)%pooldt%poollu_flay(slowp,s)
      !        sibg%l(myl)%pooldt%poollu_dloss(slowp,s) = &
      !             sibg%l(myl)%pooldt%poollu_dloss(slowp,s) &
      !             + myemis * sibg%l(myl)%pooldt%poollu_flay(slowp,s)
      !     ENDDO
      !     flosslu(myl,slowp) = myemis
      !     curemis = curemis - myemis
      ! ENDIF
      !
      ! !...remove C from soil passive
      ! IF (curemis .gt. dnzero) THEN
      !
      !    !tmppooltot = tmppooltot + sibg%l(myl)%pooldt%curpoollu(armp)
      !
      !     myemis = MIN(curemis, MAX(dzero, &
      !          sibg%l(myl)%pooldt%poollu(armp) &
      !          - sum(sibg%l(myl)%pooldt%poollu_dloss(armp,:))))
      !
      !    DO s=1,nsoil
      !        sibg%l(myl)%pooldt%loss_fire_lay(armp,s) = &
      !             myemis * dtisib * sibg%l(myl)%pooldt%poollu_flay(armp,s)
      !        sibg%l(myl)%pooldt%poollu_dloss(armp,s) = &
      !             sibg%l(myl)%pooldt%poollu_dloss(armp,s) &
      !             + myemis * sibg%l(myl)%pooldt%poollu_flay(armp,s)
      !     ENDDO
      !     flosslu(myl,armp) = myemis
      !     curemis = curemis - myemis
      ! ENDIF

      !...save what C was emitted but not removed
      flossb(myl) = curemis
      sibg%l(myl)%poollt%rmmd_fire = flossb(myl)*dtisib
      sibg%l(myl)%poollt%resp_fire = tempemis * dtisib

   enddo !cycling through top PFTs

   !...check to make sure all emissions have been taken from somewhere
   curemis = totemis - sum(flosspft(:,1:5)) - sum(flosslu(:,1:6)) - sum(flossb)

   IF (curemis .gt. fireb_thresh) THEN
      fb_err = .true.
   else
      fb_err = .false.
   endif
   
   !...Print Results
   IF ((fb_err) .OR. (fireb_print)) THEN
      print*,''
      print*,'---FIRE CARBON---'
      IF (fb_err) THEN
         print('(a)'),'!!Fire Carbon Imbalance!!'
         print*,'Fire Emissions Mismatch (mol C/m2): ', curemis
      ENDIF

      print('(a,a,i3,a,i4)'), &
          '      Date: ', trim(month_names(month)), day, ', ', year
      print('(a,i6,2f8.2)'),  '      Point/Lon/Lat: ', indx, lon, lat
      print('(a,i14)'),       '      Current Model Second: ', sec_tot
      print('(a,2i12)'),      '      Current/Next Fire Second: ', &
              fire_seccur, fire_secnext

      print*,''
      print('(a,i4)'),     '                                ntpft: ', &
           ntpft
      print*,''
      print('(a,f18.8)'),     '                            tfarea: ', &
           tfarea

      print*,''
      print('(a,f18.8)'),     '             sumtotal biomass(m-2): ', &
           sum(tpbiomass)
!      print*,''
!      print('(a,f6.2,a,f6.2,a,f6.2,a,f6.2,a,f6.2)'),&
!           'biomass(m-2): ', &
!           tpbiomass(1),' ',tpbiomass(2),' ',tpbiomass(3),' ',&
!           tpbiomass(4),' ',tpbiomass(5)
      print*,''
      print('(a,f18.8)'),     'sumtotal above-ground biomass(m-2): ', &
           sum(tpagb)
!      print*,''
!      print('(a,f6.2,a,f6.2,a,f6.2,a,f6.2,a,f6.2)'),&
!           'above-ground biomass(m-2): ', &
!           tpagb(1),' ',tpagb(2),' ',tpagb(3),' ',&
!           tpagb(4),' ',tpagb(5)

      print*,''
      print('(a,f18.8)'),     '      Fire C Emissions (umol/m2/s): ', &
           sibg%gprogt%firec*mol_to_umol
      print('(a,f18.8)'),     '      Time-Step C Losses (mol/m2):', totemis

      tempemis = sum(flosspft(:,1:5)) + sum(flosslu(:,1:6))
      print('(a,f18.8)'),     '      SiB4 C Removal (mol/m2):    ', tempemis
      print('(a)'),           '         PFT    Loss          %-BioBurned   Fraction'
         do l=1,ido
            myl = tsortref(l)
            curemis = sum(flosslu(myl,1:6)) + sum(flosspft(myl,1:5))
            if (tpbiomass(myl) .gt. dzero) then
               pcemis = curemis/tpbiomass(myl)*100.
            else
               pcemis = dzero
            endif
            print('(a,i2,2f14.8,a,f6.2)'), '          ', tpref(myl),  &
               curemis, pcemis, '  ',tparea(myl)/tfarea
         enddo
      print('(a,f12.8)'),     '      Non-Matched C Respired: ', sum(flossb)

      IF (fb_err .AND. (fireb_stop)) STOP
   ENDIF  !print

   deallocate(tpref,tpnum,tparea)
   deallocate(tpagb,tpagbtemp)
   deallocate(flosspft,flosslu,flossb)

endif  !firec > 0

end subroutine fire_interp

