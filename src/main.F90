!File=main.F90
!Author=wei
!Created=Sun 1* Sep 2018 03:53:02 PM DST
!Last Modified=Sun 1* Sep 2018 03:53:02 PM DST

Program main
  use mod_constant
  use mod_func
  Implicit None

  ! test
  real(RFREAL),allocatable::xyz(:,:),origin(:),sect(:)
  integer,allocatable::conn(:,:)

  real(RFREAL) r(4),x(6),ang,line(2,2,3)
  integer nx,nr1,nr2,ntheta,NptsperLine


  r = (/0.0275,0.0508, 0.1016, 0.1524/)
  x = (/-0.018, 1.276, 1.338, 1.4, 1.462, 1.524/)
  ang = 4.535_RFREAL
  nx = 400
  nr1 = 20
  nr2 = 40
  ntheta = 30

!  r = (/0.0275, 0.1016, 0.2032/)
!!  r = (/0.0275, 0.032, 0.0508/)
!  x = (/-0.018, 1.276, 1.338, 1.4/)
!  ang = 4.535_RFREAL
!  nx = 100
!  nr1 = 10
!  nr2 = 20
!  ntheta = 16


  ! centerline and Lipline
  line(1,1,:) = (/0.0, 0.0, 0.0/)     ! CL: origin
  line(1,2,:) = (/1.0, 0.0, 0.0/)     ! CL: dx,dy,dz
  line(2,1,:) = (/0.0, 0.0254, 0.0/)  ! LL: origin
  line(2,2,:) = (/1.0, 0.0, 0.0/)     ! LL: dx,dy,dz
  NptsperLine =  500  ! Npts per sampling line
  
  call enclose_surface(r,x,ang,nx,nr1,nr2,ntheta,line,NptsperLine)
 
!  call makedisk(origin,10.0_RFREAL,21,21,sect,xyz,conn)
!  
! ! call write_cgns(xyz,conn)
! ! call write_debug(xyz,conn)
!  deallocate(xyz)
!  deallocate(conn)
!
!  call makecylinder(origin,1.0_RFREAL,10.0_RFREAL,5.0_RFREAL,21,21,xyz,conn)
!
!  !call write_debug(xyz,conn)
!  call write_cgns(xyz,conn)
!  deallocate(xyz)
!  deallocate(conn)
!  
!  call makehollowdisk(origin,5.0_RFREAL,10.0_RFREAL,301,301,xyz,conn)
!!  call write_debug(xyz,conn)
!  deallocate(xyz)
!  deallocate(conn)


End Program main
