subroutine ptc_twiss(lhc,icav)
  USE S_fitting
  use set_para
  use twisspara
  implicit none
  integer i,no,mynd2,npara,mynpa,nda,icase,icav,ij
  real(dp) x(6),suml,deltap0,deltap
  type(layout) LHC
  type(real_8) y(6)
  real(kind(1d0)) get_value
  type(twiss) tw
  type(fibre), POINTER :: current
!------------------------------------------------------------------------------

  nda=0
  suml=zero

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
     write(20,'(a,13(1x,1p,e21.14))') current%MAG%name,suml,tw%mu(1),tw%mu(2),tw%beta(1,1),tw%beta(1,2),tw%beta(2,1),tw%beta(2,2),&
!     write(20,'(a,13(1x,1p,e21.14))') current%MAG%name,suml,tw%mu(1),tw%mu(2),tw%mu(3),tw%beta(1,1),tw%beta(1,2),tw%beta(1,3),&
!tw%beta(2,1),tw%beta(2,2),&
     tw%disp(1),tw%disp(3)
     current=>current%next
  enddo
  call kill(tw)
  CALL kill(y)

END subroutine ptc_twiss
