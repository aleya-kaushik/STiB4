!...calculation to determine if carbon
!...balance is maintained
subroutine balan_carbon( &
    sibpt, siblon, siblat, pref, &
    assimin, laiin, fparin, pooldt, poollt)

use kinds

use module_poolinfo
use module_pparams, only: &
    mol_to_mg, mol_to_umol, &
    month_names
use module_sib, only: &
    poold_type, pooll_type
use module_sibconst, only: &
    npoolcan, npoollu, npoolpft, &
    carbonb_print, carbonb_stop, carbonb_thresh
use module_time, only: &
    dtsib, day, month, year

implicit none

!...input variables
integer(i4), intent(in) :: sibpt, pref
real(r4), intent(in) :: siblon, siblat
real(r8), intent(in) :: assimin
real(r8), intent(in) :: fparin, laiin
type(poold_type), intent(inout) :: pooldt
type(pooll_type), intent(inout) :: poollt

!...switches
logical, parameter :: checkmol=.true.
real(r8) :: myconvert

!...assimilation variables
real(r8) :: assimtot

!...dead pool variables
real(r8), dimension(npoollu) :: dgain, dloss
real(r8), dimension(npoolpft) :: lgain, lloss

real(r8) :: dcarbonb
real(r8) :: dcarbonin, dcarbonout, dcarbons
real(r8) :: dgaingrz, dgainhrv, &
            dgaintd, dgaintl
real(r8) :: dlossr, dlosstd, dlossf
real(r8) :: dpoolinit, dpoolend, dpoolendc, dpoolchange

!...flux variables
real(r8) :: fassim, fresp
real(r8) :: respleaf, resproot
real(r8) :: respauto, respgrow, respmntn
real(r8) :: resphet, respsoil

!...grazing variables
real(r8) :: grzcarbonb

!...harvest variables
real(r8) :: hrvcarbonb
real(r8) :: hrvresp, hrvrmvd

!...live carbon pool variables
real(r8) :: lcarbonb
real(r8) :: lcarbonin, lcarbonout, lcarbons
real(r8) :: lgaina, lgains
real(r8) :: llossfire
real(r8) :: llossgrz, llossrgrz, llosstgrz
real(r8) :: llosshrv
real(r8) :: llossrg, llossrm, llossrnveg
real(r8) :: llosstd
real(r8) :: lpoolinit, lpoolend, lpoolendc, lpoolchange

real(r8) :: netcarbonb

real(r8) :: tcarbonb
real(r8) :: tlivedead, gdeadlive

!...net balance variables
logical :: cb_err, loc_printout

!...reference indices
integer(i4) :: lp, wp, crp, frp, pp
integer(i4) :: cdb, strl, metl, slit, slow, arm

!...misc variables
integer(i4) :: cref, p
real(r8) :: tempb
logical :: seed_printout

!---------------------------------------------
!...set local variables
cb_err = .false.
lp = pool_indx_leaf
wp = pool_indx_stwd
crp = pool_indx_croot
frp = pool_indx_froot
pp = pool_indx_prod
cdb = pool_indx_cdb - npoolpft
strl = pool_indx_strl - npoolpft
metl = pool_indx_metl - npoolpft
slit = pool_indx_slit - npoolpft
slow = pool_indx_slow - npoolpft
arm = pool_indx_arm - npoolpft

if (checkmol) then
   myconvert = done/mol_to_umol
else
   myconvert = done
endif
seed_printout = .false.


!----- DAILY CARBON BALANCE-----
!...dead pools
dpoolendc = (sum(pooldt%poollup) + sum(pooldt%poollu_dgain) &
     - sum(pooldt%poollu_dloss))*mol_to_mg
tempb = sum(pooldt%poollu)*mol_to_mg - dpoolendc
IF (tempb .GT. carbonb_thresh) THEN
   print*,''
   print('(a)'),'!!Dead Pool Daily Error!!'
   print('(a,a,i3,a,i4)'), &
        'Date: ', trim(month_names(month)), day, ', ', year
   print('(a,i6,2f10.2,i4)'),  &
        'Point/Lon/Lat/PFT: ',sibpt, siblon,siblat,pref
   print('(a,e14.6)'),'Mismatch (Mg C/ha): ',tempb
   print('(a,e14.6)'),'Previous Day Pool: ',sum(pooldt%poollup)*mol_to_mg
   print('(a,2e14.6)'),'Previous Day Loss/Gain: ', &
        sum(pooldt%poollu_dloss)*mol_to_mg,sum(pooldt%poollu_dgain)*mol_to_mg
   print('(a,2e14.6)'),'New Day Calculated/Saved: ',dpoolendc,sum(pooldt%poollu)*mol_to_mg
   if (carbonb_stop) stop
ENDIF

!...live pools
lpoolendc = (sum(poollt%poolpftp) + sum(poollt%poolpft_dgain) &
     - sum(poollt%poolpft_dloss))*mol_to_mg
tempb = sum(poollt%poolpft)*mol_to_mg - lpoolendc
IF (tempb .GT. carbonb_thresh) THEN
   print*,''
   print('(a)'),'!!Live Pool Daily Error!!'
   print('(a,a,i3,a,i4)'), &
        'Date: ', trim(month_names(month)), day, ', ', year
   print('(a,i6,2f10.2,i4)'),  &
        'Point/Lon/Lat/PFT: ',sibpt, siblon,siblat,pref
   print('(a,e14.6)'),'Mismatch (Mg C/ha): ',tempb
   print('(a,e14.6)'),'Previous Day Pool: ',sum(poollt%poolpftp)*mol_to_mg
   print('(a,2e14.6)'),'Previous Day Loss/Gain: ',sum(poollt%poolpft_dloss)*mol_to_mg,sum(poollt%poolpft_dgain)*mol_to_mg
   print('(a,2e14.6)'),'New Day Calculated/Saved: ',lpoolendc,sum(poollt%poolpft)*mol_to_mg
   if (carbonb_stop) stop
ENDIF

!----- TIME-STEP CARBON BALANCE-----
!...set balance variables
!.....assimilation and fluxes
fassim = assimin*mol_to_umol
assimtot = sum(poollt%gain_assim)*dtsib &
             + poollt%resp_nveg*dtsib
tempb = abs(assimin*dtsib - assimtot)
IF (tempb .GT. carbonb_thresh) THEN
   print*,'!!Assimilation Error!!'
   print('(a,i6,2f10.2,i4)'),  &
        ' Point/Lon/Lat/PFT: ',sibpt, siblon,siblat,pref
   print*,'Mismatch (mol C/m2): ',tempb, &
        assimin*dtsib - poollt%resp_nveg*dtsib
   print*,'Assimilation In: ',assimin*dtsib
   print*,'Live Pool Assim Gain: ',sum(poollt%gain_assim)*dtsib
   print*,'Non-Veg Resp Loss:',poollt%resp_nveg*dtsib
   if (carbonb_stop) stop
ENDIF

respgrow = sum(poollt%loss_gresp)*mol_to_umol
tempb = abs(poollt%resp_grow*mol_to_umol - respgrow)
IF (tempb .GT. carbonb_thresh) THEN
   print*,'!!Growth Respiration Error!!'
   print*,'  Resp_Grow: ',poollt%resp_grow*mol_to_umol
   print*,'  Loss_Gresp:',respgrow
   if (carbonb_stop) stop
ENDIF

respmntn = sum(poollt%loss_mresp_lay)*mol_to_umol
tempb = abs(poollt%resp_mntn*mol_to_umol - respmntn)
IF (tempb .GT. carbonb_thresh) THEN
   print*,'!!Maintenance Respiration Error!!'
   print*,'  Resp_Mntn: ',poollt%resp_mntn*mol_to_umol
   print*,'  Loss_Mresp:',respmntn
   if (carbonb_stop) stop
ENDIF

respleaf = (poollt%loss_gresp(lp) &
    + poollt%loss_mresp_lay(lp,1)) * mol_to_umol
tempb = abs(poollt%resp_leaf*mol_to_umol - respleaf)
IF (tempb .GT. carbonb_thresh) THEN
    print*,'!!Leaf Respiration Error!!'
    print*,'  Resp_Leaf: ',poollt%resp_leaf*mol_to_umol
    print*,'  Loss_Leaf:',respleaf
    if (carbonb_stop) stop
ENDIF

resproot = (poollt%loss_gresp(frp) + poollt%loss_gresp(crp) &
    + sum(poollt%loss_mresp_lay(frp,:)) &
    + sum(poollt%loss_mresp_lay(crp,:))) * mol_to_umol
tempb = abs(poollt%resp_root*mol_to_umol - resproot)
IF (tempb .GT. carbonb_thresh) THEN
   print*,'!!Root Respiration Error!!'
   print*,'  Resp_Root: ',poollt%resp_root*mol_to_umol
   print*,'  Loss_Root: ',resproot
   if (carbonb_stop) stop
ENDIF

respauto = sum(poollt%loss_gresp)*mol_to_umol &
     + sum(poollt%loss_mresp_lay)*mol_to_umol
tempb = abs(poollt%resp_auto*mol_to_umol - respauto)
IF (tempb .GT. carbonb_thresh) THEN
   print*,'!!Autotrophic Respiration Error!!'
   if (carbonb_stop) stop
ENDIF

resphet = sum(pooldt%loss_resp_lay)*mol_to_umol
tempb = abs(pooldt%resp_het*mol_to_umol - resphet)
IF (tempb .GT. carbonb_thresh) THEN
    print*,'!!Heterotrophic Respiration Error!!'
    print*,'Loss/Resp/Diff: ', &
        resphet, pooldt%resp_het*mol_to_umol, tempb
    if (carbonb_stop) stop
ENDIF

respsoil = dzero
do p=1,npoolpft
   if (pool_indx_lay(p) .GT. 1) then
       respsoil = respsoil &
          + poollt%loss_gresp(p)*mol_to_umol &
          + sum(poollt%loss_mresp_lay(p,:))*mol_to_umol
   endif
enddo
do p=1,npoollu
   if (pool_indx_lay(p+npoolpft) .GT. 1) then
       respsoil = respsoil &
         + sum(pooldt%loss_resp_lay(p,:))*mol_to_umol
    endif
enddo
tempb = abs(pooldt%resp_soil*mol_to_umol - respsoil)
IF (tempb .GT. carbonb_thresh) THEN
   print*,'!!Soil Respiration Error!!'
   print*,'   Soil_Resp: ',pooldt%resp_soil*mol_to_umol
   print*,'   Soil_RLoss:',respsoil
   if (carbonb_stop) stop
ENDIF

fresp = (pooldt%resp_het + poollt%resp_auto &
        + poollt%resp_nveg + poollt%resp_fire  &
        + poollt%resp_grz + poollt%resp_hrvst) &
        * mol_to_umol

!.....dead pools
do p=1,npoollu
   dgain(p) = sum(pooldt%gain_grz_lay(p,:))*mol_to_mg &
            + sum(pooldt%gain_hrvst_lay(p,:))*mol_to_mg &
            + sum(pooldt%gain_transd_lay(p,:))*mol_to_mg &
            + sum(pooldt%gain_transl_lay(p,:))*mol_to_mg
   dloss(p) = & 
            + sum(pooldt%loss_resp_lay(p,:))*mol_to_mg &
            + sum(pooldt%loss_trans_lay(p,:))*mol_to_mg &
            + sum(pooldt%loss_fire_lay(p,:))*mol_to_mg
enddo
dgain = dgain*dtsib
dloss = dloss*dtsib
dpoolinit = (sum(pooldt%poollu)*mol_to_mg - sum(dgain) + sum(dloss))
dpoolend = sum(pooldt%poollu)*mol_to_mg
dpoolchange = dpoolend - dpoolinit

dgaingrz = sum(pooldt%gain_grz_lay)*mol_to_umol*dtsib
dgainhrv = sum(pooldt%gain_hrvst_lay)*mol_to_umol*dble(dtsib)
dgaintd = sum(pooldt%gain_transd_lay)*mol_to_umol*dtsib
dgaintl = sum(pooldt%gain_transl_lay)*mol_to_umol*dtsib
dcarbonin = dgaingrz + dgainhrv &
              + dgaintd + dgaintl

dlossr = sum(pooldt%loss_resp_lay)*mol_to_umol*dtsib
dlosstd = sum(pooldt%loss_trans_lay)*mol_to_umol*dtsib
dlossf = sum(pooldt%loss_fire_lay)*mol_to_umol*dtsib
dcarbonout = dlossr + dlosstd + dlossf

dcarbons = dpoolchange/mol_to_mg*mol_to_umol
dcarbonb = dcarbonin - dcarbonout - dcarbons
tempb = abs(dcarbonb)/mol_to_umol
if (tempb .GT. carbonb_thresh) cb_err = .true.

!.....live pools
do p=1,npoolpft
   lgain(p) = poollt%gain_assim(p)*mol_to_mg*dtsib &
            + poollt%gain_seed(p)*mol_to_mg
   lloss(p) = &
            sum(poollt%loss_fire_lay(p,:))*mol_to_mg*dtsib &
            + sum(poollt%loss_hrvst_lay(p,:))*mol_to_mg*dtsib &
            + poollt%loss_gresp(p)*mol_to_mg*dtsib &
            + sum(poollt%loss_mresp_lay(p,:))*mol_to_mg*dtsib &
            !+ poollt%resp_nveg*mol_to_mg*dtsib &
            + sum(poollt%loss_trans_lay(p,:))*mol_to_mg*dtsib
enddo
do p=1,npoolcan
   cref = pool_indx_can(p)
   lloss(cref) = lloss(cref) + poollt%loss_grz(p)*mol_to_mg*dtsib
enddo

lpoolinit = (sum(poollt%poolpft)*mol_to_mg - sum(lgain) + sum(lloss))
lpoolend = sum(poollt%poolpft)*mol_to_mg
lpoolchange = lpoolend - lpoolinit

!lgaina = sum(poollt%gain_assim)*mol_to_umol*dtsib
lgaina = assimin*mol_to_umol*dtsib
lgains = sum(poollt%gain_seed)*mol_to_umol
lcarbonin = lgaina + lgains

llossfire = sum(poollt%loss_fire_lay)*mol_to_umol*dtsib
llossgrz = sum(poollt%loss_grz)*mol_to_umol*dtsib
llossrgrz = poollt%resp_grz*mol_to_umol*dtsib
llosstgrz = sum(pooldt%gain_grz_lay)*mol_to_umol*dtsib

llosshrv = sum(poollt%loss_hrvst_lay)*mol_to_umol*dtsib
llossrg = sum(poollt%loss_gresp)*mol_to_umol*dtsib
llossrm = sum(poollt%loss_mresp_lay)*mol_to_umol*dtsib
llossrnveg = poollt%resp_nveg*mol_to_umol*dtsib
llosstd = sum(poollt%loss_trans_lay)*mol_to_umol*dtsib

lcarbonout = llossgrz + llosshrv + llossfire &
     + llossrg + llossrm + llossrnveg + llosstd
lcarbons = lpoolchange/mol_to_mg*mol_to_umol

lcarbonb = lcarbonin - lcarbonout - lcarbons
tempb = abs(lcarbonb)*myconvert
if (tempb .GT. carbonb_thresh) then
    cb_err = .true.
endif

!.....grazing carbon balance
grzcarbonb = llossgrz - llossrgrz - llosstgrz
tempb = abs(grzcarbonb)*myconvert
if (tempb .GT. carbonb_thresh) cb_err = .true.

!.....harvest carbon balance
hrvresp = poollt%resp_hrvst*dtsib*mol_to_umol
hrvrmvd = poollt%rmvd_hrvst*mol_to_umol
hrvcarbonb = llosshrv &
                 - dgainhrv - hrvresp - hrvrmvd
tempb = abs(hrvcarbonb)*myconvert
if (tempb .GT. carbonb_thresh) cb_err = .true.

!.....live-to-dead transfers
tlivedead = llosstd + llosstgrz + (llosshrv - hrvresp - hrvrmvd)
gdeadlive = dgaintl + dgaingrz + dgainhrv
tcarbonb = (tlivedead - gdeadlive) &
         + (llosstd - dgaintl) + (dlosstd - dgaintd)  
tempb = abs(tcarbonb)*myconvert
if (tempb .GT. carbonb_thresh) cb_err = .true.

!.....net carbon balance
netcarbonb = dcarbonb + lcarbonb + tcarbonb &
           + grzcarbonb + hrvcarbonb
tempb = abs(netcarbonb)*myconvert
IF ((tempb .GT. carbonb_thresh) .and. &
     (.not. seed_printout)) cb_err = .true.

!.....print results
loc_printout = cb_err .or. carbonb_print

if (loc_printout) then
    print('(a)'),'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
    print('(a)'), 'TIME-STEP CARBON CYCLE BALANCE'
    print('(a,i6,2f10.2,i4)'),  &
        'Point/Lon/Lat/PFT: ',sibpt, siblon,siblat,pref
    print('(a,a,i3,a,i4)'), &
        'Date: ', trim(month_names(month)), day, ', ', year
    print('(a,f12.6,a,f12.6)'),'LAI/FPAR: ', laiin, '   ', fparin

    print('(a)'),''
    print('(a)'),'--CARBON FLUXES (micromoles C/m2/s)--'
    print('(a,f12.6)'),'  In:  ',fassim
    print('(a,f12.6)'),'  Out: ',fresp
    print('(a)'),      '       Resp        Saved        Calculated'
    print('(a,f12.6,a,f12.6)'), &
                       '       Auto    ', &
          poollt%resp_auto*mol_to_umol, '  ', respauto
    print('(a,f12.6,a,f12.6)'), &
                       '       Het     ', &
          pooldt%resp_het*mol_to_umol, '  ', resphet
    print('(a,f12.6,a,f12.6)'), &
                       '       NVeg    ', &
          poollt%resp_nveg*mol_to_umol, '  ', poollt%resp_nveg*mol_to_umol
    print('(a,f12.6,a,f12.6)'), &
                       '       Graze   ', &
          poollt%resp_grz*mol_to_umol, '  ', poollt%resp_grz*mol_to_umol
    print('(a,f12.6,a,f12.6)'), &
                       '       Harvest ', &
          poollt%resp_hrvst*mol_to_umol, '  ', poollt%resp_hrvst*mol_to_umol
    print('(a,f12.6,a,f12.6)'), &
                       '       Growth  ', &
          poollt%resp_grow*mol_to_umol, '  ', respgrow
    print('(a,f12.6,a,f12.6)'), &
                       '       Leaf    ', &
          poollt%resp_leaf*mol_to_umol, '  ', respleaf
    print('(a,f12.6,a,f12.6)'), &
                       '       Mntn    ', &
          poollt%resp_mntn*mol_to_umol, '  ', respmntn
    print('(a,f12.6,a,f12.6)'), &
                       '       Roots   ', &
          poollt%resp_root*mol_to_umol, '  ', resproot
    print('(a,f12.6,a,f12.6)'), &
                       '       Soil    ', &
          pooldt%resp_soil*mol_to_umol, '  ', respsoil
   print('(a,f12.6,a)'),'  NEE: ',fresp-fassim, '  (Out-In)'


    print('(a)'),''
    print('(a)'),'--CARBON POOLS (Mg C/ha)--'
    print('(a,f14.8,a,e14.6,f14.8)'),'  Live Carbon Init/Change/End: ', &
            lpoolinit, '   ', lpoolchange, lpoolend
    print('(a,f14.8,a,e14.6,f14.8)'),'  Dead Carbon Init/Change/End: ', &
            dpoolinit, '   ', dpoolchange, dpoolend

    print('(a)'),''
    print('(a)'),'  Pool     Orig            Gain            Loss              New'
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  Leaf  ', &
            (poollt%poolpft(lp) - sum(poollt%poolpft_dgain(lp,:)) &
               + sum(poollt%poolpft_dloss(lp,:)))*mol_to_mg, &
               '  ', sum(poollt%poolpft_dgain(lp,:))*mol_to_mg, &
               '  ', sum(poollt%poolpft_dloss(lp,:))*mol_to_mg, &
               '  ', poollt%poolpft(lp)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  FRoot ', &
            (poollt%poolpft(frp)-sum(poollt%poolpft_dgain(frp,:)) &
              + sum(poollt%poolpft_dloss(frp,:)))*mol_to_mg, &
              '  ', sum(poollt%poolpft_dgain(frp,:))*mol_to_mg, &
              '  ', sum(poollt%poolpft_dloss(frp,:))*mol_to_mg, &
              '  ', poollt%poolpft(frp)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  CRoot ', &
            (poollt%poolpft(crp)-sum(poollt%poolpft_dgain(crp,:)) &
               + sum(poollt%poolpft_dloss(crp,:)))*mol_to_mg, &
               '  ', sum(poollt%poolpft_dgain(crp,:))*mol_to_mg, &
               '  ', sum(poollt%poolpft_dloss(crp,:))*mol_to_mg, &
               '  ', poollt%poolpft(crp)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  StWd  ', &
            (poollt%poolpft(wp)-sum(poollt%poolpft_dgain(wp,:)) &
              + sum(poollt%poolpft_dloss(wp,:)))*mol_to_mg, &
              '  ', sum(poollt%poolpft_dgain(wp,:))*mol_to_mg, &
              '  ', sum(poollt%poolpft_dloss(wp,:))*mol_to_mg, &
              '  ', poollt%poolpft(wp)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  Prod  ', &
            (poollt%poolpft(pp)-sum(poollt%poolpft_dgain(pp,:)) &
              + sum(poollt%poolpft_dloss(pp,:)))*mol_to_mg, &
              '  ', sum(poollt%poolpft_dgain(pp,:))*mol_to_mg, &
              '  ', sum(poollt%poolpft_dloss(pp,:))*mol_to_mg, &
              '  ', poollt%poolpft(pp)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  CDB   ', &
            (pooldt%poollu(cdb)-sum(pooldt%poollu_dgain(cdb,:)) &
               + sum(pooldt%poollu_dloss(cdb,:)))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dgain(cdb,:))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dloss(cdb,:))*mol_to_mg, &
              '  ', pooldt%poollu(cdb)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  MetL  ', &
            (pooldt%poollu(metl)-sum(pooldt%poollu_dgain(metl,:)) &
              + sum(pooldt%poollu_dloss(metl,:)))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dgain(metl,:))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dloss(metl,:))*mol_to_mg, &
              '  ', pooldt%poollu(metl)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  StrL  ', &
            (pooldt%poollu(strl)-sum(pooldt%poollu_dgain(strl,:)) &
               + sum(pooldt%poollu_dloss(strl,:)))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dgain(strl,:))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dloss(strl,:))*mol_to_mg, &
              '  ', pooldt%poollu(strl)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  SLit  ', &
            (pooldt%poollu(slit)-sum(pooldt%poollu_dgain(slit,:)) &
               + sum(pooldt%poollu_dloss(slit,:)))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dgain(slit,:))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dloss(slit,:))*mol_to_mg, &
              '  ', pooldt%poollu(slit)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  Slow  ', &
            (pooldt%poollu(slow)-sum(pooldt%poollu_dgain(slow,:)) &
               + sum(pooldt%poollu_dloss(slow,:)))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dgain(slow,:))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dloss(slow,:))*mol_to_mg, &
              '  ', pooldt%poollu(slow)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  Arm   ', &
            (pooldt%poollu(arm)-sum(pooldt%poollu_dgain(arm,:)) &
               + sum(pooldt%poollu_dloss(arm,:)))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dgain(arm,:))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dloss(arm,:))*mol_to_mg, &
              '  ', pooldt%poollu(arm)*mol_to_mg

    print('(a)'),''
    print('(a)'),'--LIVE CARBON BALANCE (micromoles C/m2)--'
    print('(a,e14.6)'),'  Balance Error:  ',lcarbonb
    if (seed_printout) print('(a)'),'  Non-Minimal Error Due To Crop Seed Release'
    print('(a)'),       '  Gain:      Assim          Seed_Release'
    print('(a,e14.6,a,e14.6)'),         '           ', &
         lgaina, ' ', lgains
    print('(a)'),       '  Loss:      Fire           Harvest        Tran_Dead'
    print('(a,e14.6,a,e14.6,a,e14.6)'), '           ', &
         llossfire, ' ', llosshrv, ' ', llosstd
    print('(a)'),       '  Resp:      Growth         Maintenance    Non-Veg'
    print('(a,e14.6,a,e14.6,a,e14.6)'), '           ', &
         llossrg, ' ', llossrm, ' ', llossrnveg
    IF (abs(grzcarbonb) .GT. dzero) THEN
         print('(a)'), ' Graze:      Loss           Resp           Tran_Dead'
         print('(a,e14.6,a,e14.6,a,e14.6)'), '           ', &
               llossgrz, '  ', llossrgrz, '  ', llosstgrz
    ENDIF
    print('(a)'),       '  Net:       In             Out            Stored'
    print('(a,e14.6,a,e14.6,a,e14.6)'), '           ', &
          lcarbonin, '   ', lcarbonout, '  ', lcarbons
    
    print*,''
    print('(a)'),'--DEAD CARBON BALANCE (micromoles C/m2)--'
    print('(a,e14.6)'),'  Balance Error:  ',dcarbonb
    print('(a)'),       '  In:   Grazing      Harvest      Trans_Dead   Trans_Live'
    print('(a,e14.6,a,e14.6,a,e14.6,a,e14.6)'),'      ', &
          dgaingrz, '  ', dgainhrv, ' ', dgaintd, ' ', dgaintl
    print('(a)'),       '  Out:   Loss_Fire   Loss_Resp    Tran_Dead'
    print('(a,3(e14.6,a))'),        '      ', &
          dlossf, ' ', dlossr, ' ', dlosstd
    print('(a)'),       '  Net:   In          Out          Stored'
    print('(a,e14.6,a,e14.6,a,e14.6)'),        '      ', &
          dcarbonin, '  ', dcarbonout, '  ', dcarbons

    print*,''
    print('(a)'),'--TRANSFER CARBON BALANCE (micromoles C/m2)--'
    tcarbonb = (llosstd - dgaintl) + (dlosstd - dgaintd)  

    print('(a,e14.6)'),'  Balance Error:  ', tcarbonb
    print('(a)'),'  Live->Dead:  Diff          LTran_Dead     Gain_DeadL'
    print('(a,e14.6,2(a,e14.6))'),    '              ', &
              llosstd-dgaintl, '  ', llosstd, '  ', dgaintl
    print('(a)'),'  TLive->Dead: Tot Diff      Tot_LLoss      Tot_DGain'
    print('(a,e14.6,2(a,e14.6))'),     '              ', &
              tlivedead - gdeadlive, '', tlivedead, '', gdeadlive
    print('(a)'),'  Dead->Dead:  Diff          DTran_Dead     Gain_DeadD'
    print('(a,e14.6,2(a,e14.6))'),    '              ', &
              dlosstd-dgaintd, '  ', dlosstd, '  ', dgaintd


    IF (abs(hrvcarbonb) .GT. dzero) THEN
       print*,''
       print('(a)'),'--HARVEST CARBON BALANCE (micromoles C/m2)--'
       print('(a,e14.6)'),'  Balance Error:  ', hrvcarbonb
       print('(a)'),'  Live_Loss     Dead_Gain   Resp_Hrvst   Rmvd_Hrvst'
       print('(a,4(e14.6,a))'),'       ', &
             llosshrv, dgainhrv, hrvresp, hrvrmvd
    ENDIF

    print('(a)'),''
    print('(a,e14.6)'),'==>Maximum Carbon Imbalance: ', netcarbonb
             

    if (cb_err) then
       print('(a)'),''
       print('(a)'),'  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       print('(a)'),'  !!!Error in Carbon Balance!!!'
       print('(a)'),'  !!!       Stopping.       !!!'
       print('(a)'),'  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       print('(a)'),''
    endif

    print('(a)'),'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
    if (cb_err .and. carbonb_stop) stop
endif

!...update variables
pooldt%poollup(:) = pooldt%poollu(:)
poollt%poolpftp(:) = poollt%poolpft(:)


end subroutine balan_carbon
