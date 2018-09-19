!File=mod_func.F90
!Author=root
!Created=Tue 11 Sep 2018 03:17:07 PM DST
!Last Modified=Tue 11 Sep 2018 03:17:07 PM DST
module  mod_func
  use mod_constant
  Implicit None

contains 

  function deg2rad(deg)
    use mod_constant
    real(RFREAL) deg,deg2rad
    deg2rad = 2 * PI * deg / 180.0_RFREAL
  end function deg2rad

  subroutine  makedisk(origin,r,nr,ntheta,xyz,conn)
    ! parameters
    real(RFREAL),intent(in):: r,origin(3)
    integer,intent(in):: nr,ntheta
    real(RFREAL),allocatable,intent(out):: xyz(:,:)
    integer,allocatable,intent(out):: conn(:,:)
    ! local variables
    integer:: npnts,ncells,i,j,ct
    real(RFREAL):: x,y,z,dr,dtheta

    ! start
    ! calculate # of points and allocate arrays
    npnts = (nr - 1) * ntheta + 1
    ncells = (nr - 1) * ntheta
    allocate(xyz(3,npnts))
    allocate(conn(5,ncells))
    xyz = 0
    conn = 0
    ! calculate vertex coordinates
    ct = 1
    xyz(:,1) = origin(:)
    do i = 1,nr - 1
      dr = r / (nr - 1) * i
      do j = 1,ntheta
        ct = ct + 1
        ! senity check
        if(ct .gt. npnts)then
          stop "ct larger than npnts"
        end if
        dtheta = 2 * PI / ntheta * (j - 1)
        x = origin(1)
        y = origin(2) + dr * sin(dtheta)
        z = origin(3) + dr * cos(dtheta)
        xyz(1,ct) = x
        xyz(2,ct) = y
        xyz(3,ct) = z
      End do
    End do
    ! connectivity
    ct = 1
    do i = 1,nr - 1
      do j = 1,ntheta
        ! triangles
        if(ct .le. ntheta)then
          conn(1,ct) = TRI
          conn(2,ct) = 1
          conn(3,ct) = j + 1
          if(ct .eq. ntheta) then
            conn(4,ct) = 2
          else 
            conn(4,ct) = j + 2
          endif
          ! quads
        else
          conn(1,ct) = QUA
          if(j .eq. ntheta)then
            conn(2,ct) = (i - 1) * ntheta + 2
            conn(3,ct) = (i - 2)* ntheta + 2
          else 
            conn(2,ct) = (i - 1) * ntheta + j + 2
            conn(3,ct) = (i - 2)* ntheta + j + 2
          endif
          conn(4,ct) = (i - 2) * ntheta + j + 1
          conn(5,ct) = (i - 1)* ntheta + j + 1
        End if
        ct = ct + 1
      End do
    End do

  end subroutine  makedisk


  subroutine  makecylinder(origin,r,l,ang,nx,ntheta,xyz,conn)
    ! parameter
    real(RFREAL),intent(in):: r,origin(3),ang,l
    integer,intent(in):: nx,ntheta
    real(RFREAL),allocatable,intent(out):: xyz(:,:)
    integer,allocatable,intent(out):: conn(:,:)
    ! local variables
    integer:: npnts,ncells,i,j,ct
    real(RFREAL):: x,y,z,dr,dx,dtheta

    ! start
    ! calculate # of points and allocate arrays
    npnts = nx * ntheta
    ncells = (nx - 1) * ntheta
    allocate(xyz(3,npnts))
    allocate(conn(5,ncells))
    xyz = 0
    conn = 0
    ! calculate vertex coordinates
    ct = 1
    xyz(:,1) = origin(:)
    do i = 1,nx
      dx = l / (nx - 1) * (i - 1)
      dr =  r + dx * tan(deg2rad(ang))
      do j = 1,ntheta
        ! senity check
        if(ct .gt. npnts)then
          stop "ct larger than npnts"
        end if
        dtheta = 2 * PI / ntheta * (j - 1)
        x = origin(1) + dx
        y = origin(2) + dr * sin(dtheta)
        z = origin(3) + dr * cos(dtheta)
        xyz(1,ct) = x
        xyz(2,ct) = y
        xyz(3,ct) = z
        ct = ct + 1
      End do
    End do
    ! connectivity
    ct = 1
    do i = 1,nx - 1
      do j = 1,ntheta
        ! quads
        conn(1,ct) = QUA
        if(j .eq. ntheta)then
          conn(2,ct) = i * ntheta + 1
          conn(3,ct) = (i - 1) * ntheta + 1
        else
          conn(2,ct) = i * ntheta + j + 1
          conn(3,ct) = (i - 1) * ntheta + j + 1
        endif
        conn(4,ct) = (i - 1) * ntheta + j
        conn(5,ct) = i * ntheta + j
        ct = ct + 1
      End do
    End do

  end subroutine  makecylinder


  subroutine  write_cgns(xyz,conn)
    use cgns
    ! write cgns mesh file
    ! parameters 
    real(RFREAL),intent(in)::xyz(:,:)
    integer,intent(in)::conn(:,:)
    ! local variables
    integer(cgsize_t) isize(1,3)
    integer(cgsize_t),allocatable::jelem(:,:)
    integer(cgsize_t) nelem_start,nelem_end
    integer index_file,index_section,ielem_no,ier,iset,iphysdim,icelldim,   &
            index_base,index_zone,index_coord,nbdyelem,                     &
            npts,ncells,dimconn,i,j,ct
    character*32 basename,zonename,sectname

    npts = size(xyz,2)
    iphysdim = size(xyz,1)
    ncells = size(conn,2)
    dimconn = size(conn,1)
    if(iphysdim .ne. 3) stop "dim .ne. 3"
    if(dimconn .ne. 5) stop "conn dimension wrong"

    !   WRITE X, Y, Z GRID POINTS TO CGNS FILE
    write(6,'('' start writing cgns file grid.cgns'')')
    !   open CGNS file for write
    call cg_open_f('grid.cgns',CG_MODE_WRITE,index_file,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
    !   create base (user can give any name)
    basename='Base'
    icelldim=2
    iphysdim=2
    call cg_base_write_f(index_file,basename,icelldim,iphysdim,index_base,ier)
    !   define zone name (user can give any name)
    zonename = 'Zone  1'
    !   vertex size
    isize(1,1) = npts
    !   cell size
    isize(1,2) = ncells
    !   boundary vertex size (zero if elements not sorted)
    isize(1,3)=0
    !   create zone
    call cg_zone_write_f(index_file,index_base,zonename,isize,               &
      Unstructured,index_zone,ier)
    !   write grid coordinates (user must use SIDS-standard names here)
    call cg_coord_write_f(index_file,index_base,index_zone,RealDouble,       &
      'CoordinateX',xyz(1,:),index_coord,ier)
    call cg_coord_write_f(index_file,index_base,index_zone,RealDouble,       &
      'CoordinateY',xyz(2,:),index_coord,ier)
    call cg_coord_write_f(index_file,index_base,index_zone,RealDouble,       &
      'CoordinateZ',xyz(3,:),index_coord,ier)

    !  unsorted boundary elements
    nbdyelem=0
    nelem_start = 1
    nelem_end = 1
    ! read and write cell connectivity
    ct = 0
    do i = 2,ncells
      if((conn(1,i) .ne. conn(1,i-1)) .or. (i .eq. ncells))then
        ! need to write a cgns section
        ct = ct + 1
        if(i .eq. ncells) nelem_end = nelem_end + 1
        allocate(jelem(4,nelem_end - nelem_start + 1))
        jelem = conn(2:5,nelem_start:nelem_end)
        select case(conn(1,i-1))
        case(QUA)
          !  write QUAD_4 element connectivity
          write(sectname,'(''Quad'',i4)') ct
          call cg_section_write_f(index_file,index_base,index_zone,            &
            trim(sectname),QUAD_4,nelem_start,nelem_end,nbdyelem,              &
            jelem(1:4,:),index_section,ier)
        case(TRI)
          ! write TRI_3 element connectivity
          write(sectname,'(''Tri'',i4)') ct
          call cg_section_write_f(index_file,index_base,index_zone,            &
            trim(sectname),TRI_3,nelem_start,nelem_end,nbdyelem,               &
            jelem(1:3,:),index_section,ier)
        case default
          write(*,*) "element type error, element type", conn(1,i), "not found"
          stop
        end select
        deallocate(jelem)
        nelem_start = nelem_end + 1
        nelem_end = nelem_start
      else
        nelem_end = nelem_end + 1
      End if
    enddo
    !   close CGNS file
    call cg_close_f(index_file,ier)
    write(6,'('' Successfully wrote unstructured grid to file'',             &
      '' grid.cgns'')')
  end subroutine  write_cgns


  subroutine write_debug(xyz,conn)
    ! parameters 
    real(RFREAL),intent(in)::xyz(:,:)
    integer,intent(in)::conn(:,:)
    ! local variables
    integer ifile,i,j,nlines,ncol
    ! start
    ifile = 100
    open(ifile,file='debug.out',status='replace',form='formatted')
    nlines = size(xyz,2)
    ncol = size(xyz,1)
    write(ifile,'(3f8.4)') ((xyz(i,j),i=1,3),j=1,nlines)
    nlines =  size(conn,2)
    write(ifile,'(5i6)') ((conn(i,j),i=1,5),j=1,nlines)
    close(ifile)

  end subroutine write_debug

end module  mod_func
