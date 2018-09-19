!File=main.F90
!Author=wei
!Created=Sun 1* Sep 2018 03:53:02 PM DST
!Last Modified=Sun 1* Sep 2018 03:53:02 PM DST
Program main
  use mod_constant
  use mod_func
  Implicit None

  ! test
  real(RFREAL),allocatable::xyz(:,:),origin(:)
  integer,allocatable::conn(:,:)

  origin = (/0.0,0.0,0.0/)
  
  call makedisk(origin,10.0_RFREAL,201,201,xyz,conn)
!  write(*,*) xyz
!  write(*,*) conn 
!  call write_debug(xyz,conn)
  deallocate(xyz)
  deallocate(conn)

  call makecylinder(origin,10.0_RFREAL,10.0_RFREAL,5.0_RFREAL,101,101,xyz,conn)
!  write(*,*) xyz
!  write(*,*) conn 
!
  call write_debug(xyz,conn)
  call write_cgns(xyz,conn)
  deallocate(xyz)
  deallocate(conn)

End Program main
