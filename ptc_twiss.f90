subroutine ptc_twiss(lhc,icav)
  USE S_fitting
  use set_para
  use twisspara
  implicit none
  integer i,no,mynd2,npara,mynpa,nda,icase,icav,ij
  real(dp) x(6),suml
  type(layout) LHC
  type(real_8) y(6)
  real(kind(1d0)) get_value
  type(twiss) tw
  type(fibre), POINTER :: current
!------------------------------------------------------------------------------

  nda=0
  suml=zero

  icase = get_value('ptc ','icase ')
  if(icase==4) then
     default=default+only_4d+NOCAVITY
     mynpa=4
  elseif(icase==5) then
     default=default+delta
     mynpa=5
  else
     mynpa=6
  endif
  if(mynpa==6.and.icav==0) then
     default=default+delta
     mynpa=5
  endif

  default=default+time
  CALL UPDATE_STATES
  call print(default,6)

  x(:)=zero
  call find_orbit(lhc,x,1,default)
  print*,"Closed orbit: ",x
!  current=>LHC%start
!  suml=zero
!  do i=1,LHC%n
!     call track(lhc,x,i,i+1,default)
!     suml=suml+current%MAG%P%ld
!     write(20,'(a,13(1x,1p,e21.14))') current%MAG%name,suml
!     current=>current%next
!  enddo

  no = get_value('ptc ','no ')

  call init(default,no,nda,BERZ,mynd2,npara)
  call alloc(y)
  y=npara
  Y=X
  call track(lhc,y,1,default)
  call alloc(tw)
  tw=y

  y=npara
  Y=X
  current=>LHC%start
  suml=zero
  do i=1,LHC%n
     call track(lhc,y,i,i+1,default)
     suml=suml+current%MAG%P%ld
     tw=y
     write(20,'(a,13(1x,1p,e21.14))') current%MAG%name,suml,tw%mu(1),tw%mu(2),tw%beta(1,1),tw%beta(1,2)
     current=>current%next
  enddo
  call kill(tw)
  CALL kill(y)

END subroutine ptc_twiss
