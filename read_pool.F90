
! Opens and reads in SiB pool parameters.
subroutine read_pool()

use kinds
use module_sibconst, only: &
    cornsoy_switch, &
    npft, ngroup, &
    npoolpft, npoollu, ntpool
use module_io, only: &
    pool_file, poolid
use module_pparams, only: &
    secs_per_day, days_per_year, &
    drytoc, mwc
use module_pftinfo, only: &
    npft_gdd, &
    pft_mze, pft_soy, &
    pft_group, group_grass, group_crop
use module_poolinfo, only: pool_indx_leaf
use module_param, only: &
    physcon, poolcon

implicit none

!...file variables
integer(i4) :: finpft, fingddpft
integer(i4) :: finpool, fingroup
integer(i4) :: num
character(len=3), dimension(npft) :: pftname
character(len=10), dimension(ntpool) :: poolname
character(len=100) :: trash
logical :: iscomment1, iscomment2

real(r4), dimension(npoollu+2) :: graze_trans, harvest_trans

!...misc variables
integer(i4) :: i,j
integer(byte) :: groupref
real(r4) :: poolval

!-------------------
!...Initialize the pool variables
allocate(poolcon(npft))
call init_poolcon(npft, npoolpft, npoollu, poolcon)

!...Open file
print*,'Reading Pool Parameter File: '
print*,'  ',trim(pool_file)
open(unit=poolid,file=trim(pool_file),form='formatted')

read(poolid,*) finpft
if (finpft /= npft) then
    print*,''
    print('(a)'),'!!!Pool File Mistmatch!!!'
    print('(a,i4,a,i4)'),'  Pool file npft: ', &
          finpft,' Sim npft: ',npft
    print*,''
    stop
endif

read(poolid,*) fingddpft
if (fingddpft /= npft_gdd) then
    print*,''
    print('(a)'),'!!!Pool File Mistmatch!!!'
    print('(a,i4,a,i4)'),'  Pool file npft_gdd: ', &
          fingddpft,' Sim npft_gdd: ',npft_gdd
    print*,''
    stop
endif
read(poolid,*) finpool
read(poolid,*) fingroup

if (finpool /= ntpool) then
    print*,''
    print('(a)'),'!!!Pool param file does not match simulation!!!'
    print('(a,i4,a,i4)'),'  Pool file npool: ',finpool,' Sim npool: ',ntpool
    print*,''
    stop
endif

if (fingroup /= ngroup) then
    print*,''
    print('(a)'),'!!!Pool param file does not match simulation!!!'
    print('(a,i4,a,i4)'),'  Pool file ngroup: ',fingroup,' Sim ngroup: ',ngroup
    print*,''
    stop
endif

iscomment1=.true.
iscomment2=.true.
do while ((iscomment1) .or. (iscomment2))
    read(poolid,*) trash
    if (index(trash,'****') .gt. 0) then
        if (iscomment1) then
           iscomment1=.false.
        else 
           iscomment2=.false.
        endif
    endif
enddo

!...Read in the variables
read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
do i=1, npft
   read(poolid,*) num, pftname(i), &
        poolcon(i)%gr_frac(:)
enddo

read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
do i=1, npft
   read(poolid,*) num, pftname(i), &
        poolcon(i)%lresp_eff(:)
enddo

read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
do i=1, npft
   read(poolid,*) num, pftname(i), &
        poolcon(i)%cr_aml,  poolcon(i)%cr_amh,  &
        poolcon(i)%cr_amin, poolcon(i)%cr_amax, &
        poolcon(i)%cr_fmul, poolcon(i)%cr_fref, &
        poolcon(i)%cr_fmin, &
        poolcon(i)%cr_hq10, poolcon(i)%cr_href, &
        poolcon(i)%cr_hmax
enddo


read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
do i=1, npft
   read(poolid,*) num, pftname(i), &
        poolcon(i)%lt_fref, poolcon(i)%lt_fq10, &
        poolcon(i)%lt_fmax, &
        poolcon(i)%lt_dref, poolcon(i)%lt_dcoef, &
        poolcon(i)%lt_dmax, &
        poolcon(i)%lt_wref, poolcon(i)%lt_wbase, &
        poolcon(i)%lt_wcoef, poolcon(i)%lt_wmax
enddo

read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
do i=1, npft
   read(poolid,*) num, pftname(i), &
        poolcon(i)%rrt_aml, poolcon(i)%rrt_amh, &
        poolcon(i)%rrt_amin, poolcon(i)%rrt_amax, &
        poolcon(i)%rrt_fmul, poolcon(i)%rrt_fref, &
        poolcon(i)%rrt_fmin, &
        poolcon(i)%rrt_hq10, poolcon(i)%rrt_href, &
        poolcon(i)%rrt_hmax, poolcon(i)%rrt_laimin, &
        poolcon(i)%rrt_laimax
enddo

read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
do i=1, npft
   read(poolid,*) num, pftname(i), &
        poolcon(i)%hrt_sfc_aml, poolcon(i)%hrt_sfc_amh, &
        poolcon(i)%hrt_sfc_amin, poolcon(i)%hrt_sfc_amax, &
        poolcon(i)%hrt_sfc_fmul, poolcon(i)%hrt_sfc_fref, &
        poolcon(i)%hrt_sfc_fmin, &
        poolcon(i)%hrt_sfc_hq10, poolcon(i)%hrt_sfc_href, &
        poolcon(i)%hrt_sfc_hmax, &
        poolcon(i)%hrt_sfc_pml, poolcon(i)%hrt_sfc_pmin
enddo

read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
do i=1, npft
   read(poolid,*) num, pftname(i), &
        poolcon(i)%hrt_soil_aml, poolcon(i)%hrt_soil_amh, &
        poolcon(i)%hrt_soil_amin, poolcon(i)%hrt_soil_amax, &
        poolcon(i)%hrt_soil_fmul, poolcon(i)%hrt_soil_fref, &
        poolcon(i)%hrt_soil_fmin, &
        poolcon(i)%hrt_soil_hq10, poolcon(i)%hrt_soil_href, &
        poolcon(i)%hrt_soil_hmax, &
        poolcon(i)%hrt_soil_mmin, poolcon(i)%hrt_soil_pawmin
enddo


read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
do i=1, npft
   read(poolid,*) num, pftname(i), &
        poolcon(i)%turnover(:)
enddo
        
read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash

!.....Biological efficiencies and transfers
do i=1,npft
   read(poolid,*) trash
   read(poolid,*) trash
   read(poolid,*) trash
   read(poolid,*) trash
   do j=1, npoollu
      read(poolid,*) num, poolname(j), poolcon(i)%dresp_eff(:,j)
   enddo

   read(poolid,*) trash
   read(poolid,*) trash
   read(poolid,*) trash
   do j=1, ntpool
      read(poolid,*) num, poolname(j), poolcon(i)%pool_trans_frac(:,j)
   enddo
enddo

!.....Specialty transfers
read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
do i=1,npoollu+2
   read(poolid,*) num, poolname(i), graze_trans(i), harvest_trans(i)
enddo   

close(poolid)

!...Save the parameters
do i=1,npft
    groupref = pft_group(i)

   !...Set specialty transfers
   if (groupref .eq. group_grass) then
      poolcon(i)%graze_trans(:) = graze_trans(:)
   endif

   if (groupref .eq. group_crop) then
      if ((sum(harvest_trans) .gt. 0.999) .and. &
          (sum(harvest_trans) .lt. 1.001)) then
          poolcon(i)%harvest_trans(:) = harvest_trans(:)
      else
         print*,''
         print*,'--Incorrect Harvest Transfer Fractions--'
         print*,'Must Sum To 1, Current Sum: ',sum(harvest_trans)
         print*,'Stopping.'
         print*,''
         stop
     endif
   endif

   !...Set calculated parameters
   do j=1,ntpool
      if (poolcon(i)%turnover(j) .gt. 1.e-12) then
         poolcon(i)%k_rate(j) = 1./ &
              (poolcon(i)%turnover(j)*real(secs_per_day)*real(days_per_year))
      endif
   enddo
   
  !...set pool minimum values
  if (physcon(i)%sla .gt. 1.E-12) then
     poolval = physcon(i)%laimin / physcon(i)%sla &
                / real(drytoc) / real(mwc)
      poolcon(i)%poolpft_min(pool_indx_leaf) = poolval
   endif !sla > 0

enddo

if (cornsoy_switch) then
   poolval = max(poolcon(pft_mze)%poolpft_min(1), &
                 poolcon(pft_soy)%poolpft_min(1))
   poolcon(pft_mze)%poolpft_min(pool_indx_leaf) = poolval
   poolcon(pft_soy)%poolpft_min(pool_indx_leaf) = poolval
endif


end subroutine read_pool

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initializes the pool parameters.
subroutine init_poolcon(npft, npoolpft, npoollu, &
     poolcon)

  use kinds
  use module_param, only: pool_param

  !...input variables
  integer(i4), intent(in) :: npft, npoolpft, npoollu
  type(pool_param), dimension(npft), intent(inout) :: poolcon

  !...local variables
  integer(i4) :: i, ntpool

  !...set local variables
  ntpool = npoolpft + npoollu

  !...initialize pool parameters
  do i=1,npft
     allocate(poolcon(i)%lresp_eff(npoolpft))
     poolcon(i)%lresp_eff(:) = rzero
     allocate(poolcon(i)%gr_frac(npoolpft))
     poolcon(i)%gr_frac(:) = rzero

     poolcon(i)%cr_aml = rzero
     poolcon(i)%cr_amh = rzero
     poolcon(i)%cr_amin = rone
     poolcon(i)%cr_amax = rone
     poolcon(i)%cr_fmul = rzero
     poolcon(i)%cr_fref = rzero
     poolcon(i)%cr_fmin = rone
     poolcon(i)%cr_hq10 = rzero
     poolcon(i)%cr_href = rzero
     poolcon(i)%cr_hmax = rone

     poolcon(i)%lt_fq10 = rzero
     poolcon(i)%lt_fref = rzero
     poolcon(i)%lt_fmax = rzero
     poolcon(i)%lt_dcoef = rzero
     poolcon(i)%lt_dref = rzero
     poolcon(i)%lt_dmax = rzero
     poolcon(i)%lt_wref = rzero
     poolcon(i)%lt_wbase = rzero
     poolcon(i)%lt_wcoef = rzero
     poolcon(i)%lt_wmax = rzero

     poolcon(i)%rrt_aml = rzero
     poolcon(i)%rrt_amh = rzero
     poolcon(i)%rrt_amin = rone
     poolcon(i)%rrt_amax = rone
     poolcon(i)%rrt_fmul = rzero
     poolcon(i)%rrt_fref = rzero
     poolcon(i)%rrt_fmin = rone
     poolcon(i)%rrt_hq10 = rzero
     poolcon(i)%rrt_href = rzero
     poolcon(i)%rrt_hmax = rone
     poolcon(i)%rrt_laimin = rone
     poolcon(i)%rrt_laimax = rone
     
     poolcon(i)%hrt_sfc_aml = rzero
     poolcon(i)%hrt_sfc_amh = rzero
     poolcon(i)%hrt_sfc_amin = rone
     poolcon(i)%hrt_sfc_amax = rone
     poolcon(i)%hrt_sfc_fmul = rone
     poolcon(i)%hrt_sfc_fref = rone
     poolcon(i)%hrt_sfc_fmin = rone
     poolcon(i)%hrt_sfc_hq10 = rone
     poolcon(i)%hrt_sfc_href = rone
     poolcon(i)%hrt_sfc_hmax = rone
     poolcon(i)%hrt_sfc_pml = rzero
     poolcon(i)%hrt_sfc_pmin = rone

     poolcon(i)%hrt_soil_aml = rzero
     poolcon(i)%hrt_soil_amh = rzero
     poolcon(i)%hrt_soil_amin = rone
     poolcon(i)%hrt_soil_amax = rone
     poolcon(i)%hrt_soil_fmul = rone
     poolcon(i)%hrt_soil_fref = rone
     poolcon(i)%hrt_soil_fmin = rone
     poolcon(i)%hrt_soil_hq10 = rone
     poolcon(i)%hrt_soil_href = rone
     poolcon(i)%hrt_soil_hmax = rone
     poolcon(i)%hrt_soil_mmin = rone
     poolcon(i)%hrt_soil_pawmin = rone

     allocate(poolcon(i)%dresp_eff(npoollu,npoollu))
     poolcon(i)%dresp_eff(:,:) = rzero
     allocate(poolcon(i)%graze_trans(npoollu+2))
     allocate(poolcon(i)%harvest_trans(npoollu+2))
     poolcon(i)%graze_trans(:) = rzero
     poolcon(i)%harvest_trans(:) = rzero

     allocate(poolcon(i)%turnover(ntpool))
     poolcon(i)%turnover(:) = rzero
     allocate(poolcon(i)%pool_trans_frac(ntpool,ntpool))
     poolcon(i)%pool_trans_frac(:,:) = rzero
     allocate(poolcon(i)%k_rate(ntpool))
     poolcon(i)%k_rate(:) = rzero

     allocate(poolcon(i)%poolpft_min(npoolpft))
     poolcon(i)%poolpft_min(:) = rzero

  enddo !i=1,npft

end subroutine init_poolcon
