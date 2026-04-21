!==========================================================================
subroutine set_continue_spinup()
!==========================================================================
! Resets the pools for continuing spinup runs.

use kinds
use module_sibconst, only: &
    subcount, spinup_done, &
    npoolpft, npoollu, &
    spinup_continue
use module_io, only: requib_writef
use module_pftinfo, only: pft_num
use module_poolinfo, only: pool_indx_lay
use module_sib, only: sib

implicit none

!...local variables
integer :: i,l,n,s
integer(i4) :: pref, pnum

!... call and run continue version of equipools_calc()
call equipools_calc_continue()

!...Reset equilibrium switch
requib_writef = .false.

!...Set the pools for the next spin-up iteration
if (.not. spinup_done) then
   do i=1,subcount
      if (.not. sib%g(i)%gdiagt%gridcell_spunup) then
          do l=1,sib%g(i)%g_nlu

             if (.not. sib%g(i)%l(l)%equibdt%lupft_spunup) then
                pref = sib%g(i)%l(l)%ipft
                pnum = pft_num(pref)

                call equipools_restart_continue(pnum, &
                    sib%g(i)%l(l)%equiblt%poolpft_end,   &
                    sib%g(i)%l(l)%equiblt%poolpft_equib, &
                    sib%g(i)%l(l)%poollt%poolpft_flay,   &
                    sib%g(i)%l(l)%poollt%poolpft_lay)

                !...reset live pools
                do n=1,npoolpft
                   sib%g(i)%l(l)%equiblt%poolpft_init(n) = &
                       sum(sib%g(i)%l(l)%poollt%poolpft_lay(n,:))
                   sib%g(i)%l(l)%poollt%poolpft(n) = &
                       sib%g(i)%l(l)%equiblt%poolpft_init(n)
                   sib%g(i)%l(l)%equiblt%poolpft_min(n) = &
                       sib%g(i)%l(l)%equiblt%poolpft_init(n)
                   sib%g(i)%l(l)%equiblt%poolpft_max(n) = &
                        sib%g(i)%l(l)%equiblt%poolpft_init(n)

                   sib%g(i)%l(l)%equiblt%poolpft_totgain(n) = dzero
                   sib%g(i)%l(l)%equiblt%poolpft_totloss(n) = dzero

               enddo !n=1,npoolpft
               sib%g(i)%l(l)%poollt%poolpftp(:) = sib%g(i)%l(l)%poollt%poolpft(:)

               !...reset dead pools
               do n=1,npoollu
                  sib%g(i)%l(l)%equibdt%poollu_init(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)
                  sib%g(i)%l(l)%pooldt%poollu(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)
                  sib%g(i)%l(l)%equibdt%poollu_min(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)
                  sib%g(i)%l(l)%equibdt%poollu_max(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)

                  do s=1,pool_indx_lay(n+npoolpft) !(6,11)->pool_indx_lay either 1 or 10
                  !pool_indx_lay(ntpool), poollu_lay(npoollu,nsoil)
                      sib%g(i)%l(l)%pooldt%poollu_lay(n,s) = &
                              sib%g(i)%l(l)%equibdt%poollu_equib(n) * &
                              sib%g(i)%l(l)%pooldt%poollu_flay(n,s)
                   enddo
                   sib%g(i)%l(l)%equibdt%poollu_totgain(n) = dzero
                   sib%g(i)%l(l)%equibdt%poollu_totloss(n) = dzero
               enddo !n=1,npoollu
               sib%g(i)%l(l)%pooldt%poollup(:) = sib%g(i)%l(l)%pooldt%poollu(:)

             endif  !.not. lu_spunup

          enddo  !l=1,g_nlu
      endif  !.not. gridcell_spunup
   enddo  !subcount
endif  !spinup_done

end subroutine set_continue_spinup



!==========================================================================
subroutine equipools_calc_continue()
!==========================================================================
! Calculates the quasi-equilibrium carbon pools.

use kinds
use module_param, only: poolcon
use module_pftinfo, only: &
    pft_num, pft_type, type_bare
use module_poolinfo
use module_sib, only: sib
use module_sibconst, only: &
    single_pt, subcount,   &
    npoolpft, npoollu, &
    spinup_threshold, spinup_done, &
    spinup_continue

implicit none

!...local variables
real(r8) :: ave_gain    !(mol/m2/s) average external inputs per pool
real(r8) :: ave_loss    !(mol/m2/s) average external outputs per pool
real(r8) :: ave_k_rate  !(1/s) average scaled decay rate constant
real(r8) :: pool_init, pool_end, init_ratio, end_ratio
real(r8) :: pdiffr, pdiffi, pdiffe

!...misc variables
integer(byte) :: ptype
integer(i4) :: i,l,n
integer(i4) :: lp,frp,crp,wp,pp
integer(i4) :: cdbp, metlp, strlp, slitp, slowp, armp
integer(i4) :: pref, pnum

!-------Calculatae quasi-equlibrium pools------------
lp =  pool_indx_leaf !ntpool index 1
frp = pool_indx_froot !ntpool index 2
crp = pool_indx_croot !ntpool index 3
wp =  pool_indx_stwd !ntpool index 4
pp =  pool_indx_prod !ntpool index 5

cdbp  = pool_indx_cdb-npoolpft !ntpool index 6, npoollu index 1
metlp = pool_indx_metl-npoolpft !ntpool index 7, npoollu index 2
strlp = pool_indx_strl-npoolpft !ntpool index 8, npoollu index 3
slitp = pool_indx_slit-npoolpft !ntpool index 9, npoollu index 4
slowp = pool_indx_slow-npoolpft !ntpool index 10, npoollu index 5
armp  = pool_indx_arm-npoolpft !ntpool index 11, npoollu index 6

spinup_done = .true.

!...loop through gridcell points
do i=1,subcount
   if (.not. sib%g(i)%gdiagt%gridcell_spunup) then
      sib%g(i)%gdiagt%gridcell_spunup = .true.
      !...loop through landunits
      do l=1,sib%g(i)%g_nlu
         if (.not. sib%g(i)%l(l)%equibdt%lupft_spunup) then

             sib%g(i)%l(l)%equibdt%lupft_spunup = .true.
             pref=sib%g(i)%l(l)%ipft
             pnum=pft_num(pref)
             ptype=pft_type(pnum)

             !--------------
             !--DEAD POOLS--
             do n=1,npoollu
                pool_init = sib%g(i)%l(l)%equibdt%poollu_init(n)
                pool_end = sib%g(i)%l(l)%pooldt%poollu(n)

                !...calculate time average decay rates

                !...for continuing spinup, this avoids pool_end
                !...being set as pool_equib below?? (doesn't work)
                !if (spinup_continue) then
                !   sib%g(i)%l(l)%equibdt%poollu_totgain(n) = dzero
                !   sib%g(i)%l(l)%equibdt%poollu_totloss(n) = dzero
                !endif

                ave_gain = sib%g(i)%l(l)%equibdt%poollu_totgain(n)
                ave_loss = sib%g(i)%l(l)%equibdt%poollu_totloss(n)

                if (pool_end .gt. dzero) then
                     ave_k_rate = ave_loss / pool_end
                else
                     ave_k_rate = dzero
                endif

                !...save the net gains/losses/end
                sib%g(i)%l(l)%equibdt%poollu_gain(n) = ave_gain
                sib%g(i)%l(l)%equibdt%poollu_loss(n) = ave_loss
                sib%g(i)%l(l)%equibdt%poollu_end(n) = pool_end

                !...solve for ratio of input/output
                if (ave_loss .gt. dzero) then
                    sib%g(i)%l(l)%equibdt%poollu_ratio(n) = ave_gain / ave_loss
                endif
                pdiffr=abs(sib%g(i)%l(l)%equibdt%poollu_ratio(n) - done)

                !...solve for equilibrium pool sizes
                if (ave_k_rate .gt. dzero) then
                     sib%g(i)%l(l)%equibdt%poollu_equib(n) = &
                          ave_gain / ave_k_rate
                elseif ((ave_gain > dzero) .or. (ave_loss > dzero)) then
                     sib%g(i)%l(l)%equibdt%poollu_equib(n) = pool_end
                else
                     sib%g(i)%l(l)%equibdt%poollu_equib(n) = dzero
                     sib%g(i)%l(l)%equibdt%poollu_ratio(n) = done
                endif !zero input/output

                !...calculate starting/ending ratios for spinup constraints
                if (sib%g(i)%l(l)%equibdt%poollu_equib(n) .gt. dzero) then
                    init_ratio = pool_init/sib%g(i)%l(l)%equibdt%poollu_equib(n)
                    end_ratio = pool_end/sib%g(i)%l(l)%equibdt%poollu_equib(n)
                else
                    init_ratio = done
                    end_ratio = done
                endif
                pdiffi = abs(init_ratio - done)
                pdiffe = abs(end_ratio - done)

                !...determine if the dead pools are spunup
                if ((pdiffr <= spinup_threshold) .and. &
                    ((pdiffi <= spinup_threshold) .or. &
                     (pdiffe <= spinup_threshold))) then
                     sib%g(i)%l(l)%equibdt%poollu_notdone(n) = .false.
                else
                     sib%g(i)%l(l)%equibdt%poollu_notdone(n) = .true.
                endif
             enddo  !n=1,npoollu

             !...Determine if the surface pools are spun-up
             sib%g(i)%l(l)%equibdt%deadsfc_init  = &
                    sib%g(i)%l(l)%equibdt%poollu_init(cdbp) + &
                    sib%g(i)%l(l)%equibdt%poollu_init(metlp) + &
                    sib%g(i)%l(l)%equibdt%poollu_init(strlp)
             sib%g(i)%l(l)%equibdt%deadsfc_end = &
                    sib%g(i)%l(l)%equibdt%poollu_end(cdbp) + &
                    sib%g(i)%l(l)%equibdt%poollu_end(metlp) + &
                    sib%g(i)%l(l)%equibdt%poollu_end(strlp)
             sib%g(i)%l(l)%equibdt%deadsfc_gain = &
                    sib%g(i)%l(l)%equibdt%poollu_gain(cdbp) + &
                    sib%g(i)%l(l)%equibdt%poollu_gain(metlp) + &
                    sib%g(i)%l(l)%equibdt%poollu_gain(strlp)
             sib%g(i)%l(l)%equibdt%deadsfc_loss = &
                    sib%g(i)%l(l)%equibdt%poollu_loss(cdbp) + &
                    sib%g(i)%l(l)%equibdt%poollu_loss(metlp) + &
                    sib%g(i)%l(l)%equibdt%poollu_loss(strlp)

             if (sib%g(i)%l(l)%equibdt%deadsfc_loss > 0.) then
                    sib%g(i)%l(l)%equibdt%deadsfc_ratio = &
                        sib%g(i)%l(l)%equibdt%deadsfc_gain / &
                        sib%g(i)%l(l)%equibdt%deadsfc_loss
             else
                   sib%g(i)%l(l)%equibdt%deadsfc_ratio = 1.0
             endif

             pdiffr = abs(sib%g(i)%l(l)%equibdt%deadsfc_ratio - 1.0)
             if (pdiffr <= spinup_threshold) then
                 sib%g(i)%l(l)%equibdt%deadsfc_notdone = .false.
             else
                 sib%g(i)%l(l)%equibdt%deadsfc_notdone = .true.
             endif

             !...Determine if the soil pools are spunup
             sib%g(i)%l(l)%equibdt%deadsoil_init = &
                    sib%g(i)%l(l)%equibdt%poollu_init(slitp) + &
                    sib%g(i)%l(l)%equibdt%poollu_init(slowp) + &
                    sib%g(i)%l(l)%equibdt%poollu_init(armp)
             sib%g(i)%l(l)%equibdt%deadsoil_end = &
                    sib%g(i)%l(l)%equibdt%poollu_end(slitp) + &
                    sib%g(i)%l(l)%equibdt%poollu_end(slowp) + &
                    sib%g(i)%l(l)%equibdt%poollu_end(armp)
             sib%g(i)%l(l)%equibdt%deadsoil_gain = &
                    sib%g(i)%l(l)%equibdt%poollu_gain(slitp) + &
                    sib%g(i)%l(l)%equibdt%poollu_gain(slowp) + &
                    sib%g(i)%l(l)%equibdt%poollu_gain(armp)
             sib%g(i)%l(l)%equibdt%deadsoil_loss = &
                    sib%g(i)%l(l)%equibdt%poollu_loss(slitp) + &
                    sib%g(i)%l(l)%equibdt%poollu_loss(slowp) + &
                    sib%g(i)%l(l)%equibdt%poollu_loss(armp)

             if (sib%g(i)%l(l)%equibdt%deadsoil_loss > dzero) then
                    sib%g(i)%l(l)%equibdt%deadsoil_ratio = &
                        sib%g(i)%l(l)%equibdt%deadsoil_gain / &
                        sib%g(i)%l(l)%equibdt%deadsoil_loss
             else
                   sib%g(i)%l(l)%equibdt%deadsoil_ratio = 1.0
             endif

             pdiffr = abs(sib%g(i)%l(l)%equibdt%deadsoil_ratio - 1.0)
             if (pdiffr <= spinup_threshold) then
                 sib%g(i)%l(l)%equibdt%deadsoil_notdone = .false.
             else
                 sib%g(i)%l(l)%equibdt%deadsoil_notdone = .true.
             endif


             !--------------
             !--LIVE POOLS--
             sib%g(i)%l(l)%equiblt%poolpft_end(:) = &
                    sib%g(i)%l(l)%poollt%poolpft(:)

                !...For bare ground do not calculate equilibrium pools,
                !...instead set them to end pools
                if (ptype == type_bare) then
                    sib%g(i)%l(l)%equibdt%poollu_notdone(:) = .false.
                    sib%g(i)%l(l)%equibdt%poollu_equib(:) = &
                          sib%g(i)%l(l)%equibdt%poollu_end(:)

                    sib%g(i)%l(l)%equibdt%deadsfc_notdone = .false.
                    sib%g(i)%l(l)%equibdt%deadsoil_notdone = .false.

                    sib%g(i)%l(l)%equiblt%poolpft_notdone(:) = .false.
                    sib%g(i)%l(l)%equiblt%poolpft_equib(:) = &
                           sib%g(i)%l(l)%equiblt%poolpft_end(:)
                    sib%g(i)%l(l)%equiblt%poolpft_equib(:) = 0.0

                    sib%g(i)%l(l)%equiblt%live_notdone = .false.
                else
                    !...Calculate equilibrium PFT pools
                    do n=1,npoolpft

                        !...calculate time average decay rates
                        ave_gain = sib%g(i)%l(l)%equiblt%poolpft_totgain(n)
                        ave_loss = sib%g(i)%l(l)%equiblt%poolpft_totloss(n)
                        IF (sib%g(i)%l(l)%poollt%poolpft(n) .gt. dzero) THEN
                            ave_k_rate = ave_loss /sib%g(i)%l(l)%poollt%poolpft(n)
                        ELSE
                            ave_k_rate = dzero
                        ENDIF

                        !...save net gains/losses
                        sib%g(i)%l(l)%equiblt%poolpft_gain(n) = ave_gain
                        sib%g(i)%l(l)%equiblt%poolpft_loss(n) = ave_loss
                        if (ave_loss > dzero) then
                           sib%g(i)%l(l)%equiblt%poolpft_ratio(n) = ave_gain/ave_loss
                        else
                           sib%g(i)%l(l)%equiblt%poolpft_ratio(n) = done
                        endif

                        !...solve for equilibrium pool sizes
                        if (ave_k_rate > dzero) then
                            sib%g(i)%l(l)%equiblt%poolpft_equib(n) = &
                                 MAX( poolcon(pnum)%poolpft_min(n), &
                                 ave_gain / ave_k_rate )
                        elseif ((ave_gain > dzero) .or. (ave_loss > dzero)) then
                            sib%g(i)%l(l)%equiblt%poolpft_equib(n) = &
                                MAX( poolcon(pnum)%poolpft_min(n), &
                                sib%g(i)%l(l)%equiblt%poolpft_end(n))
                        else
                           if (((single_pt) .or. (subcount == 1)) .and. &
                                (n .eq. lp)) then
                               print*,'   !!!Error in Equilibrium - Leaf Pool Zero In/Out!!!'
                               print*,'   !!!Setting to Minimum Value!!!'
                               print('(a,2F10.4,a,a)'),'    Input/Output/Pool:', &
                                    ave_gain,ave_k_rate,'  ', pool_name(n)
                                    !stop
                            endif
                            sib%g(i)%l(l)%equiblt%poolpft_equib(n) = poolcon(pnum)%poolpft_min(n)
                            sib%g(i)%l(l)%equiblt%poolpft_ratio(n) = 1.0
                         endif !zero input/output

                         !...determine if the pools are spun-up
                         pdiffr=abs(sib%g(i)%l(l)%equiblt%poolpft_ratio(n) - 1.0)
                         if (pdiffr <= spinup_threshold) then
                             sib%g(i)%l(l)%equiblt%poolpft_notdone(n) = .false.
                         else
                             sib%g(i)%l(l)%equiblt%poolpft_notdone(n) = .true.
                         endif

                    enddo !n=1,npoolpft

                    !...determine if using total live carbon is spun-up
                    sib%g(i)%l(l)%equiblt%live_init = &
                         sib%g(i)%l(l)%equiblt%poolpft_init(lp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_init(wp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_init(frp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_init(crp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_init(pp)
                    sib%g(i)%l(l)%equiblt%live_end  = &
                         sib%g(i)%l(l)%equiblt%poolpft_end(lp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_end(wp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_end(frp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_end(crp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_end(pp)
                    sib%g(i)%l(l)%equiblt%live_gain = &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(lp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(wp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(frp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(crp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(pp)
                    sib%g(i)%l(l)%equiblt%live_loss = &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(lp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(wp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(frp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(crp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(pp)
                    if (sib%g(i)%l(l)%equiblt%live_loss > dzero) then
                          sib%g(i)%l(l)%equiblt%live_ratio = &
                                 sib%g(i)%l(l)%equiblt%live_gain / &
                                 sib%g(i)%l(l)%equiblt%live_loss
                    else
                          sib%g(i)%l(l)%equiblt%live_ratio = 1.
                    endif

                    !...test for spinup determination
                    pdiffr = abs(sib%g(i)%l(l)%equiblt%live_ratio - 1.0)
                    if (pdiffr < spinup_threshold) then
                          sib%g(i)%l(l)%equiblt%live_notdone = .false.
                    else
                          sib%g(i)%l(l)%equiblt%live_notdone = .true.
                    endif

                endif !pft not bare

                !!!set spinup_done using individual pools!!!
                if  (sib%g(i)%l(l)%equiblt%poolpft_notdone(lp)  .or. &
                     sib%g(i)%l(l)%equiblt%poolpft_notdone(frp) .or. &
                     sib%g(i)%l(l)%equiblt%poolpft_notdone(crp) .or. &
                     sib%g(i)%l(l)%equiblt%poolpft_notdone(wp)  .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(cdbp) .or.  &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(metlp) .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(strlp) .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(slitp) .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(slowp) .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(armp)) then
                         sib%g(i)%l(l)%equibdt%lupft_spunup = .false.
                         sib%g(i)%gdiagt%gridcell_spunup = .false.
                         spinup_done = .false.
                 endif

                !!!set spinup_done for PFT pools using live total!!!
                !if (sib%g(i)%l(l)%equibdt%deadsfc_notdone .or. &
                !     sib%g(i)%l(l)%equibdt%deadsoil_notdone .or. &
                !     sib%g(i)%l(l)%equiblt%live_notdone) then
                !        sib%g(i)%l(l)%equibdt%lupft_spunup = .false.
                !        sib%g(i)%gdiagt%gridcell_spunup = .false.
                !        spinup_done = .false.
                !endif

         endif  !.not. lupft_spunup
      enddo !landunit

   endif  !.not. gridcell_spunup
enddo  !subcount

end subroutine equipools_calc_continue


!============================================
subroutine equipools_restart_continue(pnum, &
     poolpft_end, poolpft_equib, poolpft_flay, poolpft_out)
!============================================
! Sets the equilibrium pools to values that can
!   be used to restart a simulation.

use kinds
use module_sibconst, only: &
   npoolpft, nsoil, &
   spinup
use module_pftinfo, only: &
   pft_type, pft_group, &
   type_decid, type_grass, type_crop
use module_poolinfo, only: &
   pool_indx_leaf, pool_indx_stwd, &
   pool_indx_prod, pool_indx_lay

implicit none

!...input variables
integer(i4), intent(in) :: pnum
real(r8), dimension(npoolpft), intent(in) :: poolpft_end, poolpft_equib
real(r8), dimension(npoolpft,nsoil), intent(in) :: poolpft_flay
real(r8), dimension(npoolpft,nsoil), intent(inout) :: poolpft_out

!...local variables
integer(i4) :: n,s
integer(byte) :: ptype, pgroup
real(r8) :: poolrestart
integer(i4) :: lp,wp,pp

!-------Calculatae quasi-equlibrium pools------------
lp =  pool_indx_leaf
wp =  pool_indx_stwd
pp =  pool_indx_prod

!---------------------------
!...Set local variables
ptype = pft_type(pnum)
pgroup = pft_group(pnum)
poolpft_out(:,:) = 0.


!...Set output PFT pools
do n=1,npoolpft !1,5
    poolrestart = poolpft_equib(n)

    !if (.not. spinup) then
      !For deciduous PFTs, set restart to end value
      !for leaf and product pools
      if (ptype == type_decid) then
          if ((n == lp) .or. &
              (n == pp)) then
               poolrestart = MAX(dzero, poolpft_end(n))
          endif
      endif

      !For grass PFTs, set restart to end value
      ! for leaf, stem, and product pools
      if (ptype == type_grass) then
          if ((n == lp) .or. &
              (n == wp) .or. &
              (n == pp)) then
               poolrestart = MAX(dzero, poolpft_end(n))
          endif
      endif

      !For crop PFTs, set restart to end value
      ! for all live pools
      if (ptype == type_crop) then
         poolrestart = MAX(dzero, poolpft_end(n))
      endif
    !endif

    do s=1,pool_indx_lay(n) !for 1,5 ntpool
       poolpft_out(n,s) = poolrestart * &
            poolpft_flay(n,s)
    enddo
enddo


end subroutine equipools_restart_continue


