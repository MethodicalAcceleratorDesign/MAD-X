subroutine ptc_normal(lhc,icav)
  USE S_fitting
  implicit none
  integer i,no,mynd2,npara,mynpa,nda,icase,icav,ij
  real(dp) x(6),deltap0,deltap
  type(layout) LHC
  type(real_8) y(6)
  type(normalform) n
  real(kind(1d0)) get_value
  type(fibre), POINTER :: current
!------------------------------------------------------------------------------

  nda=0

  icase = get_value('ptc ','icase ')
  deltap0 = get_value('ptc ','deltap ')
  deltap = zero
  select case(icase)
  CASE(4)
     default=default+only_4d+NOCAVITY
     mynpa=4
  CASE(5)
     default=default+delta
     deltap=deltap0
     mynpa=5
  CASE(6)
     mynpa=6
  CASE DEFAULT
     default=default+only_4d+NOCAVITY
     mynpa=4
  END SELECT
  if(mynpa==6.and.icav==0) then
     default=default+delta
     mynpa=5
  endif

  default=default+time

  CALL UPDATE_STATES
  call print(default,6)

  x(:)=zero
  if(icase.eq.5) x(5)=deltap
  call find_orbit(lhc,x,1,default)
  print*,"Closed orbit: ",x

  no = get_value('ptc ','no ')

  call init(default,no,nda,BERZ,mynd2,npara)
  call alloc(y)
  y=npara
  Y=X
  call track(lhc,y,1,default)
  call alloc(n)
  n=y
  call daprint(y,18)
  write(19,'(/a/)') 'Disperion, First and Higher Orders'
  call daprint(n%A1,19)
  write(19,'(/a/)') 'Tunes, Chromaticities and Anharmonicities'
!  call daprint(n%A_t,19)
!  call daprint(n%A,19)
  call daprint(n%dhdj,19)
  call kill(n)
  CALL kill(y)

END subroutine ptc_normal
