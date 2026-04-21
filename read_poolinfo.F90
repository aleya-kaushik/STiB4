
! Opens and reads in pool information.
subroutine read_poolinfo()

    use kinds
    use module_sibconst, only: &
         npoolpft, npoollu, ntpool, &
         npoolcan, npoolsfc, npoolsoil
    use module_io, only: &
         pool_info, piid
    use module_poolinfo

    implicit none

    !...file variables
    integer(i4) :: num
    character(len=30) :: trash

    !...misc variables
    integer(byte) :: iref
    integer(i4) :: i

!---------------------------------
!...Open file
print*,''
print*,'Reading Pool Informational File: '
print*,'  ',trim(pool_info)
open(unit=piid,file=trim(pool_info),form='formatted')

read(piid,*) npoolpft
read(piid,*) npoollu
ntpool = npoolpft + npoollu

read(piid,*) trash
read(piid,*) trash
read(piid,*) trash
read(piid,*) trash
allocate(pool_name_long(ntpool))
allocate(pool_name(ntpool))
allocate(pool_type(ntpool))
allocate(pool_loc(ntpool))
allocate(pool_indx_lay(ntpool))
do i=1,ntpool
   read(piid,'(i6,a20,a2,a6,a4,a4,a4,a7,a3,i2)') &
        num, pool_name_long(i), trash, &
        pool_name(i), trash, &
        pool_type(i), trash, &
        pool_loc(i), trash, &
        pool_indx_lay(i)
enddo

close(piid)

!-----Categorize pools-----------------------------------------
!...Clear out all index variables
npoolcan=0
npoolsfc=0
npoolsoil=0

pool_indx_leaf=0
pool_indx_froot=0
pool_indx_croot=0
pool_indx_stwd=0
pool_indx_prod=0
pool_indx_metl=0
pool_indx_strl=0
pool_indx_slit=0
pool_indx_slow=0
pool_indx_arm=0

!...scan through pool information and count pools of each type
pool_type(:) = adjustl(pool_type(:))
pool_loc(:) = adjustl(pool_loc(:))
pool_name(:) = adjustl(pool_name(:))

do i=1,ntpool
   if(trim(pool_loc(i))=='canopy') &
       npoolcan = npoolcan + bone
   if(trim(pool_loc(i))=='surface') &
       npoolsfc = npoolsfc + bone
   if(trim(pool_loc(i))=='soil') &
        npoolsoil = npoolsoil + bone

    iref=int(i,kind=byte)
    if(trim(pool_name(i))=='leaf') pool_indx_leaf=iref
    if(trim(pool_name(i))=='froot') pool_indx_froot=iref
    if(trim(pool_name(i))=='croot') pool_indx_croot=iref
    if(trim(pool_name(i))=='stwd') pool_indx_stwd=iref
    if(trim(pool_name(i))=='prod') pool_indx_prod=iref
    if(trim(pool_name(i))=='cdb')     pool_indx_cdb=iref
    if(trim(pool_name(i))=='metl')  pool_indx_metl=iref
    if(trim(pool_name(i))=='strl')  pool_indx_strl=iref
    if(trim(pool_name(i))=='slit')    pool_indx_slit=iref
    if(trim(pool_name(i))=='slow')    pool_indx_slow=iref
    if(trim(pool_name(i))=='arm')     pool_indx_arm=iref
enddo

!...set the pool type/location indices
allocate(pool_indx_can(npoolcan))
iref=bone
do i=1, ntpool
   if (trim(pool_loc(i)) == 'canopy') then 
       pool_indx_can(iref)=int(i,kind=byte)
       iref = iref + bone
   endif
enddo

allocate(pool_indx_sfc(npoolsfc))
iref=bone
do i=1, ntpool
   if (trim(pool_loc(i)) == 'surface') then
       pool_indx_sfc(iref)=int(i,kind=byte)
       iref = iref + bone
    endif
enddo

allocate(pool_indx_soil(npoolsoil))
iref=bone
do i=1, ntpool
   if (trim(pool_loc(i)) == 'soil') then
      pool_indx_soil(iref)=int(i,kind=byte)
      iref=iref+bone
    endif
enddo

!----------------------
!...Print Information
print('(a,3i4)'),'   Pool Number (Tot/PFT/LU):     ', ntpool, npoolpft, npoollu
print('(a,3i4)'),'   Pool Location (Can/Sfc/Soil): ', npoolcan, npoolsfc, npoolsoil
print*,''

end subroutine read_poolinfo
