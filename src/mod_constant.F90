module mod_constant
  implicit none
  
  integer,parameter:: SPREAL = selected_real_kind(6,37)
  integer,parameter:: RFREAL = selected_real_kind(15,307)
  real,parameter:: PI = 3.1415926536
  integer,parameter:: TRI = -1,                      &
                      QUA = -2

end module

