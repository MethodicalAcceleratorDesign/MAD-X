!The Full Polymorphic Package
!Copyright (C) Etienne Forest and Frank Schmidt
! See file a_scratch_size
module define_newda
  use precision_constants
  implicit none
  logical(lp) warnda
  integer firsthole
  logical(lp),parameter::checkass=.true.,newdadynamical=.true., PACKING=.false.
  integer,parameter::totaldol=3000

  type  taylorlow
     real(dp) ,  DIMENSION(:), POINTER :: R
     ! integer  id
     integer, POINTER :: id
     integer,  DIMENSION(:), POINTER :: NZ
     logical(lp), dimension(:), pointer ::  yes
     integer, POINTER :: m
     ! integer  m
  end  type taylorlow

end module define_newda
