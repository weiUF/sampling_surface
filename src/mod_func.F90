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
    deg2rad = PI * deg / 180.0_RFREAL
  end function deg2rad
  
  subroutine makeline(origin,dist,n,xyz,conn)
    ! make line for sampling hydrodynamic quantities for source calculation
    ! arguments
    real(RFREAL),intent(in)::origin(3),dist(3)
    integer,intent(in)::n
    real(RFREAL),allocatable,intent(out)::xyz(:,:)
    integer,allocatable,intent(out):: conn(:,:)
    ! local variables
    integer:: i

    allocate(xyz(3,n))
    allocate(conn(5,n-1))
    
    do i = 1,n
      xyz(1,i) = origin(1) + dble(i - 1) / (n - 1) * dist(1)
      xyz(2,i) = origin(2) + dble(i - 1) / (n - 1) * dist(2)
      xyz(3,i) = origin(3) + dble(i - 1) / (n - 1) * dist(3)
      conn(1,i) = BAR
      conn(2,i) = i
      conn(3,i) = i + 1
    end do



  end subroutine makeline


  subroutine  makedisk(origin,r,nr,ntheta,xyz,conn)
    ! arguments
    real(RFREAL),intent(in):: r,origin(3)
    integer,intent(in):: nr,ntheta
    real(RFREAL),allocatable,intent(out):: xyz(:,:)
    integer,allocatable,intent(out):: conn(:,:)
    ! local variables
    integer:: npnts,ncells,i,j,ct,nintsct,ct_intsct
    real(RFREAL):: x,y,z,dr,dtheta

    ! start
    ! calculate # of points and allocate arrays
    npnts = (nr - 1) * ntheta + 1
    ncells = (nr - 1) * ntheta
    allocate(xyz(3,npnts))
    allocate(conn(5,ncells))
    xyz = 0
    conn = 0
    !! number of intersections
    !! intsct is to ensure grid point on the intersection of planes
    !nintsct = size(intsct)
    !! sanity check: intsct(i+1) > intsct(i) 
    !do i = 1,nintsct - 1
    !  if(intsct(i + 1) .le. intsct(i)) then
    !    write(6,*) 'intersection radius array is not in ascending order.'
    !    stop
    !  End if
    !enddo
    !ct_intsct = 1
    ! calculate vertex coordinates
    ct = 1
    xyz(:,1) = origin(:)
    do i = 1,nr - 1
      dr = r / (nr - 1) * i
      !! check if planes intersect within this layer of cell
      !if(dr .gt. intsct(ct_intsct))then
      !  dr = intsct(ct_intsct)
      !  ct_intsct = ct_intsct + 1
      !  ! sanity check
      !  if(ct_intsct .gt. nintsct)then
      !    write(6,*) 'STOP: largest intersection r need to equal r of disk'
      !    write(6,*) 'r:',dr,'sect:',intsct(nintsct)
      !    stop
      !  End if
      !End if
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
          conn(4,ct) = j + 1
          if(ct .eq. ntheta) then
            conn(3,ct) = 2
          else 
            conn(3,ct) = j + 2
          endif
          ! quads
        else
          conn(1,ct) = QUA
          if(j .eq. ntheta)then
            conn(3,ct) = (i - 1) * ntheta + 2
            conn(2,ct) = (i - 2)* ntheta + 2
          else 
            conn(3,ct) = (i - 1) * ntheta + j + 2
            conn(2,ct) = (i - 2)* ntheta + j + 2
          endif
          conn(5,ct) = (i - 2) * ntheta + j + 1
          conn(4,ct) = (i - 1)* ntheta + j + 1
        End if
        ct = ct + 1
      End do
    End do

  end subroutine  makedisk


  subroutine  makehollowdisk(origin,r1,r2,nr,ntheta,xyz,conn,direction)
    ! arguments
    real(RFREAL),intent(in):: r1,r2,origin(3)
    integer,intent(in):: nr,ntheta,direction
    real(RFREAL),allocatable,intent(out):: xyz(:,:)
    integer,allocatable,intent(out):: conn(:,:)
    ! local variables
    integer:: npnts,ncells,i,j,ct,temp
    real(RFREAL):: x,y,z,dr,dtheta
    ! start
    if(r1 .eq. 0._RFREAL) then
      call makedisk(origin,r2,nr,ntheta,xyz,conn)
    else
      ! calculate # of points and allocate arrays
      npnts = nr * ntheta
      ncells = (nr - 1) * ntheta
      allocate(xyz(3,npnts))
      allocate(conn(5,ncells))
      xyz = 0
      conn = 0
      ! calculate vertex coordinates
      ct = 1
      do i = 1,nr
        dr = r1 + (r2 - r1) / (nr - 1) * (i - 1)
        do j = 1,ntheta
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
          ct = ct + 1
        End do
      End do
      ! connectivity
      ct = 1
      do i = 1,nr - 1
        do j = 1,ntheta
          ! senity check
          if(ct .gt. ncells)then
            stop "ct larger than ncells"
          end if
          conn(1,ct) = QUA
          if(j .eq. ntheta)then
            conn(2,ct) = i * ntheta + 1
            conn(3,ct) = (i - 1)* ntheta + 1
          else 
            conn(2,ct) = i * ntheta + j + 1
            conn(3,ct) = (i - 1) * ntheta + j + 1
          endif
          conn(4,ct) = (i - 1) * ntheta + j
          conn(5,ct) = i * ntheta + j
          if(direction .eq. 1)then
            ! face x+ direction, 
            !invert conn to flow right hand rule
            temp = conn(2,ct)
            conn(2,ct) = conn(3,ct)
            conn(3,ct) = temp
            temp =  conn(4,ct)
            conn(4,ct) = conn(5,ct)
            conn(5,ct) = temp
          elseif(direction .eq. 0)then
            ! face x- dirction, do nothing
          else
            print *, "face direction for hollow disk not recognized: ", direction
            stop
          endif
          ct = ct + 1
        End do
      End do
    endif
  end subroutine  makehollowdisk

  subroutine  makecylinder(origin,r,l,ang,nx,ntheta,xyz,conn)
    ! arguments
    real(RFREAL),intent(in):: r,origin(3),ang,l
    integer,intent(in):: nx,ntheta
    real(RFREAL),allocatable,intent(out):: xyz(:,:)
    integer,allocatable,intent(out):: conn(:,:)
    ! local variables
    integer:: npnts,ncells,i,j,ct
    real(RFREAL):: x,y,z,dr,dx,dtheta,growthrate,ddx,dx1

    ! start
    growthrate=1.01
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
    dx1 = l * (1 - growthrate) / (1 - growthrate ** (nx - 1))
    do i = 1,nx
      dx = dx1 * (1 - growthrate ** (i - 1)) / (1 - growthrate)
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

  subroutine fillin_array(xyz,conn,xyz_part,conn_part,npts,ncells,weight_array,weight,pointrange,subfaceno)
    ! arguments
    real(RFREAL),intent(in)::xyz_part(:,:),weight
    integer,intent(in)::conn_part(:,:)
    real(RFREAL),intent(inout)::xyz(:,:),weight_array(:)
    integer,intent(inout)::conn(:,:),ncells,npts,pointrange(:),subfaceno
    ! local variables
    xyz(:,npts + 1: npts + size(xyz_part,2)) = xyz_part
    conn(2:,ncells + 1: ncells + size(conn_part,2)) = conn_part(2:,:) + npts
    conn(1,ncells + 1: ncells + size(conn_part,2)) = conn_part(1,:) 
    weight_array(ncells + 1: ncells + size(conn_part,2)) = weight

    npts = npts + size(xyz_part,2)
    ncells = ncells + size(conn_part,2)
    pointrange(subfaceno) = ncells
    subfaceno = subfaceno + 1

  end subroutine fillin_array

  subroutine enclose_surface(r,x,ang,nx,nr1,nr2,ntheta,line,NptsperLine)
    ! make enclosed surface consist of disk, hollowdisk and cylinder    
    !-----------------------------------------------------------
    ! usage:
    !
    !-----------------------------------------------------------
    
    ! arguments
    real(RFREAL),intent(in)::r(:),x(:),ang,line(:,:,:)
    integer,intent(in)::nx,nr1,nr2,ntheta,NptsperLine
    ! local variables 
    !real(RFREAL),allocatable::r2(:,:)
    integer,parameter::maxpts = 500000
    integer,parameter::maxcells = 500000
    real(RFREAL):: r1,dr1,r2,dr2,dx,origin(3),weight_func
    integer::nsurfaces,nendcaps,npts,ncells,nr_temp,nx_temp,i,j
    real(RFREAL)::xyz(3,maxpts),weight(maxcells)
    integer:: conn(5,maxcells),pointrange(100),subfaceno
    real(RFREAL),allocatable::xyz_temp(:,:)
    integer,allocatable::conn_temp(:,:)

    ! start
    nsurfaces = size(r) - 1
    nendcaps = size(x) - 1
    npts = 0
    ncells = 0
    subfaceno = 1
    !allocate(r2(nsurfaces))
    do i = 1,nsurfaces  
      origin = (/x(1),0._RFREAL,0._RFREAL/)
      r1 = r(i)
      dr1 = r(i + 1) - r(i)
      nr_temp = nint(dr1 / (r(ubound(r,1)) - r(1)) * nr1)
      call makehollowdisk(origin,r1,r1 + dr1,nr_temp,ntheta,xyz_temp,conn_temp,0)
      call fillin_array(xyz,conn,xyz_temp,conn_temp,npts,ncells,weight,1._RFREAL,pointrange,subfaceno)
      deallocate(xyz_temp)
      deallocate(conn_temp)
      if (i .eq. 1)then
        r1 = 0
        dr1 = r(2)
      End if
      do j = 1,nendcaps
        ! make cylinder between disks (begin cap - end cap or end cap - end cap)
        weight_func = (nendcaps + 1.0_RFREAL - j) / nendcaps
        origin = (/x(j),0._RFREAL,0._RFREAL/)
        dx = x(j + 1) - x(j)
        if(i .eq. 1)then
          r2 = r1
          dr2 = dr1 + tan(deg2rad(ang)) * dx
        else
          r2 = r1 + tan(deg2rad(ang)) * dx
          dr2 = dr1
        End if
        nx_temp = nint(dx / (x(ubound(x,1)) - x(1)) * nx)
        call makecylinder(origin,r1+dr1,dx,ang,nx_temp,ntheta,xyz_temp,conn_temp)
        call fillin_array(xyz,conn,xyz_temp,conn_temp,npts,ncells,weight,weight_func,pointrange,subfaceno)
        deallocate(xyz_temp)
        deallocate(conn_temp)
        ! end cap
        weight_func = 1.0_RFREAL / nendcaps
        origin = (/x(j + 1),0.0_RFREAL,0.0_RFREAL/)
        nr_temp = nint(dr2 / (r2 + dr2 + r(ubound(r,1)) - r(i+1)) * nr2)
        !debug
        !print*, 'surface:',i,'endcap:',j,'nr:', nr_temp, 'rtotal:',r2 + dr2 + r(ubound(r,1)) - r(j)
        call makehollowdisk(origin,r2,r2 + dr2,nr_temp,ntheta,xyz_temp,conn_temp,1)
        call fillin_array(xyz,conn,xyz_temp,conn_temp,npts,ncells,weight,weight_func,pointrange,subfaceno)
        deallocate(xyz_temp)
        deallocate(conn_temp)
        r1 = r2
        dr1 = dr2
      End do
    End do

    ! add sampling points on lines -- lipline and centerline
    ! CenterLine
    open(102,file='CL_LL.out')
    write(102,'(a,i8)') '#nlines',size(line,1)
    write(102,'(a,i8)') '#pnts', NptsperLine
    call makeline(line(1,1,:),line(1,2,:),NptsperLine,xyz_temp,conn_temp)
    call fillin_array(xyz,conn,xyz_temp,conn_temp,npts,ncells,weight,0._RFREAL,pointrange,subfaceno)
    xyz(:,npts + 1: npts + size(xyz_temp,2)) = xyz_temp
    npts = npts + NptsperLine
    write(102,'(2i8)') npts - NptsperLine + 1, npts
    deallocate(xyz_temp)
    deallocate(conn_temp)
    print *, 'Centerline: ', NptsperLine, 'pnts, origin: ', line(1,1,:)
    ! LipLine
    write(102,'(a,i8)') '#pnts', NptsperLine
    call makeline(line(2,1,:),line(2,2,:),NptsperLine,xyz_temp,conn_temp)
    call fillin_array(xyz,conn,xyz_temp,conn_temp,npts,ncells,weight,0._RFREAL,pointrange,subfaceno)
    xyz(:,npts + 1: npts + size(xyz_temp,2)) = xyz_temp
    npts = npts + NptsperLine
    write(102,'(2i8)') npts - NptsperLine + 1, npts
    deallocate(xyz_temp)
    deallocate(conn_temp)
    close(102)
    print *, 'Lipline: ', NptsperLine, 'pnts, origin: ', line(2,1,:)

    ! write info to output, write cngs file
    print *,pointrange(:subfaceno)
    print * , npts, ncells
    call write_cgns(xyz(:,:npts),conn(:,:ncells),weight(:ncells),pointrange(:subfaceno),nsurfaces,nendcaps)

    
    call write_debug(xyz(:,:npts),conn(:,:ncells))
  
  end subroutine enclose_surface


  subroutine  write_cgns(xyz,conn,weight,pointrange_in,nsurfaces,nendcaps)
    use cgns
    implicit none
    ! write cgns mesh file
    ! arguments 
    real(RFREAL),intent(in)::xyz(:,:),weight(:)
    integer,intent(in):: conn(:,:)
    integer,intent(in)::pointrange_in(:)
    integer,intent(in):: nsurfaces,nendcaps
    ! local variables
    integer(cgsize_t) isize(1,3)
    integer(cgsize_t),allocatable::jelem(:,:),bcpointlist(:)
    integer(cgsize_t) nelem_start,nelem_end
    integer index_file,index_section,ielem_no,ier,iset,iphysdim,icelldim,   &
      index_base,index_zone,index_coord,index_bc,nbdyelem,index_flow,       &
      npts,ncells,dimconn,i,j,k,ct,nbcpts
    character*32 basename,zonename,sectname,error_msg

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
        case(BAR)
          ! write BAR_2 element connectivity
          write(sectname,'(''BAR'',i4)') ct
          call cg_section_write_f(index_file,index_base,index_zone,            &
            trim(sectname),BAR_2,nelem_start,nelem_end,nbdyelem,               &
            jelem(1:2,:),index_section,ier)
        case default
          write(*,*) "element type error, element type", conn(1,i-1), "not found"
          stop
        end select
        deallocate(jelem)
        nelem_start = nelem_end + 1
        nelem_end = nelem_start
      else
        nelem_end = nelem_end + 1
      End if
    enddo

    ! write weight function array
    !call cg_goto_f(index_file,index_base,ier,'Zone_t',index_zone,'end')
    !call cg_user_data_write_f('weights',ier)
    !call cg_goto_f(index_file,index_base,ier,'Zone_t',index_zone,'weights',0,'end')
    !call cg_gridlocation_write_f(CellCenter, ier)
    !call cg_ptset_write_f(PointRange,2,(/1,ncells/),ier)
    !call cg_array_write_f('weight',RealDouble,1,ncells,weight,ier)

    call cg_sol_write_f(index_file,index_base,index_zone,'face weight',CellCenter,index_flow,ier)
    call cg_field_write_f(index_file,index_base,index_zone,1,RealDouble,'weight',weight,index_flow,ier)

    ! write BC for each FWH surfaces

    !if(size(pointrange) - 1 .ne. nsurf * (2 * nend + 1))then
    !  print*, "size of pointrange wrong, should be ",nsurf * (2 * nend + 1)
    !  stop
    !End if
    
    do i = 1,nsurfaces
      write(sectname,'(''BC'',i4)') i
      allocate(bcpointlist(pointrange_in(i * (2 * nendcaps + 1))))
      ct = 0
      npts = 0
      do j = 1,(i - 1) * (2 * nendcaps + 1),2
        if(j .eq. 1) then 
          nelem_start = 1
        else
          nelem_start = pointrange_in(j - 1) + 1
        endif
        nelem_end = pointrange_in(j)
        npts = npts + nelem_end - nelem_start + 1
        
        ! temp
        print *, 'i:',i,'j:',j
        print * , 'start:',nelem_start,'end:',nelem_end
        ! end temp

        do k = nelem_start,nelem_end
          ct = ct + 1
          bcpointlist(ct) = k
        enddo
      enddo
      do j = (i - 1) * (2 * nendcaps + 1) + 1,i  * (2 * nendcaps + 1)
        if(j .eq. 1) then 
          nelem_start = 1
        else
          nelem_start = pointrange_in(j - 1) + 1
        endif
        nelem_end = pointrange_in(j)
        npts = npts + nelem_end - nelem_start + 1

        ! temp
        print *, 'i:',i,'j:',j
        print * , 'start:',nelem_start,'end:',nelem_end
        ! end temp

        do k = nelem_start,nelem_end
          ct = ct + 1
          bcpointlist(ct) = k
        enddo
      enddo

      ! temp
      !print *, 'i:',i
      !print * , 'npts:',npts, 'ct:', ct
      !print *, bcpointlist(:npts)
      ! end temp

      call cg_boco_write_f(index_file,index_base,index_zone,trim(sectname),  &
        BCTypeNull,PointList,npts,bcpointlist(:npts),index_bc,ier)
      call cg_goto_f(index_file,index_base,ier,'Zone_t',1,                   &
        'ZoneBC_t',1,'BC_t',index_bc,'end')
      !call cg_gridlocation_write_f(FaceCenter,ier)
      ! ? here write FaceCenter returns error ier=1
      call cg_gridlocation_write_f(CellCenter,ier)
      ! temp
      !open(101,file='debug.out')
      !write(101,'(ai4)') 'bcpointlist',i
      !write(101,'(10i8)') bcpointlist(:npts)
      !close(101)
      ! end temp
      deallocate(bcpointlist)
    enddo


    !   close CGNS file
    call cg_close_f(index_file,ier)
    write(6,'('' Successfully wrote unstructured grid to file'',             &
      '' grid.cgns'')')
  end subroutine  write_cgns


  subroutine write_debug(xyz,conn)
    ! arguments 
    real(RFREAL),intent(in)::xyz(:,:)
    integer,intent(in)::conn(:,:)
    ! local variables
    integer ifile,i,j,nlines,ncol
    ! start
    ifile = 100
    open(ifile,file='debug.out',status='replace',form='formatted')
    nlines = size(xyz,2)
    ncol = size(xyz,1)
    write(ifile, '(i8)') nlines
    write(ifile,'(3f8.4)') ((xyz(i,j),i=1,3),j=1,nlines)
    !nlines =  size(conn,2)
    !write(ifile,'(5i6)') ((conn(i,j),i=1,5),j=1,nlines)
    close(ifile)

  end subroutine write_debug

end module  mod_func
