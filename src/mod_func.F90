!File=mod_func.F90
!Author=root
!Created=Tue 11 Sep 2018 03:17:07 PM DST
!Last Modified=Tue 11 Sep 2018 03:17:07 PM DST
module  mod_func
    use constants
    Implicit None

contains 

    subroutine  disk(origin,r,nr,ntheta,xyz,conn)
        ! parameters
        real(RFREAL),intent(in):: r,origin(3)
        integer,intent(in):: nr,ntheta
        real(RFREAL),allocatable,intent(out):: xyz(:,:)
        integer,allocatable,intent(out):: conn(:,:)
        ! local variables
        integer:: npnts,ncells,i,j,ct
        real(RFREAL):: x,y,z

        ! start
        ! calculate # of points and allocate arrays
        npnts = (nr - 1) * ntheta + 1
        ncells = (nr - 1) * ntheta
        allocate(xyz(3,npnts))
        allocate(conn(5,ncells))
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
            conn(4,ct) = j + 2
            ! quads
        else
            conn(1,ct) = QUA
            conn(2,ct) = (i - 1) * ntheta + j
            conn(3,ct) = (i - 1) * ntheta + j + 1
            conn(4,ct) = i * ntheta + j
            conn(5,ct) = i * ntheta + j + 1
        End if
        ct = ct + 1
        End do
        End do

    end subroutine  disk


    subroutine  write_vtk(xyz,conn)
      ! parameters
      real(RFREAL),intent(in):: xyz(:,:)   
      integer,intent(in)::conn
      ! variables


    end subroutine  write_vtk

end module  mod_func
