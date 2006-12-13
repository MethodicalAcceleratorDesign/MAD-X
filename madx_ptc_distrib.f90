module madx_ptc_distrib_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! madx_ptc_distrib module
! Piotr K. Skowronski (CERN)
!
! This module contains service for trackin

  use madx_keywords
  use madx_ptc_intstate_module, only : getdebug
  implicit none

  include "madx_ptc_distrib.inc"
  save
  private

  !============================================================================================
  !  PUBLIC INTERFACE
  public                         :: aremomentson
  public                         :: putmoments
  public                         :: initmoments
  public                         :: allocmoments
  public                         :: killmoments
  public                         :: getdistrtype
  public                         :: setemittances
  public                         :: addmoment
  public                         :: getnmoments
  public                         :: getmomentstabcol

  !============================================================================================
  !  PRIVATE
  !    data structures

  real(dp), pointer              ::  normmoments(:,:,:)
  real(dp)                       ::  sigmas(6)
  integer                        ::  distributiontype(3)
  
  integer, parameter             ::  maxnmoments=100
  type (momentdef), private      ::  moments(maxnmoments)
  integer                        ::  nmoments = 0

  type(damap)                    :: damapa  !da map needed to calculation initialized with sigmas
  type(taylor)                   :: function_to_average ! taylor used to calculate averages
  
  !============================================================================================
  !  PRIVATE
  !    routines

  public                         :: makegaus
  public                         :: makeflat5
  public                         :: makeflat56
  private                        :: filter
  
     
contains
  !____________________________________________________________________________________________

  !____________________________________________________________________________________________


  subroutine addmoment(x,px,y,py,dp,t,tableIA, columnIA, parametric )
    implicit none
    include 'twissa.fi'
    integer               :: x,px,y,py,dp,t
    integer               :: columnIA(*)
    integer               :: tableIA(*)
    integer               :: parametric
    integer, parameter    :: zeroasciicode = IACHAR ( '0' )
    integer               :: i
    real(kind(1d0))       :: get_value

    nmoments = nmoments + 1

    do i=1,6
       moments(nmoments)%iarray(1) =  x
       moments(nmoments)%iarray(2) =  px
       moments(nmoments)%iarray(3) =  y
       moments(nmoments)%iarray(4) =  py
       moments(nmoments)%iarray(5) =  dp
       moments(nmoments)%iarray(6) =  t
    enddo
    
    moments(nmoments)%table = charconv(tableIA)
    moments(nmoments)%column = charconv(columnIA)

    moments(nmoments)%table(tableIA(1)+1:tableIA(1)+1)=achar(0)
    moments(nmoments)%column(columnIA(1)+1:columnIA(1)+1)=achar(0)

    if (parametric /= 0) then
      print*,"To be made as parametric variable"
      moments(nmoments)%index = nmoments !for the time being it is the same as the index 
    else
      moments(nmoments)%index = 0   
    endif    
    
    
!    if (getdebug()>3) then
       print  *,"addmoment : <", moments(nmoments)%iarray(1:6) ,&
      &                   ">,<", moments(nmoments)%column,&
      &                   ">,<", moments(nmoments)%table, ">)"
!    endif



  end subroutine addmoment

  !____________________________________________________________________________________________

  integer function getnmoments()
    implicit none
       getnmoments = nmoments
  end function getnmoments
  !____________________________________________________________________________________________

  subroutine  getmomentstabcol(n, tabn, coln)
    implicit none
    integer              :: n
    character(20)        :: tabn
    character(17)        :: coln
    
    if ( (n < 1) .and. (n > nmoments ) ) then
      
      tabn(1:1) = achar(0)
      coln(1:1) = achar(0)
      call fort_warn("getmomentstabcol","index out of range")
      return
    endif
    
   
   tabn = moments(n)%table  
   coln = moments(n)%column  
    
    
  end subroutine 
  !____________________________________________________________________________________________


  logical(lp) function aremomentson()
    implicit none

      if ( associated(normmoments) ) then
         aremomentson = my_true
      else   
         aremomentson = my_false
      endif
    
  end function aremomentson
  !____________________________________________________________________________________________

  subroutine putmoments(n,name,s,y)
    implicit none
    integer              :: n !fibre number
    character(*)         :: name !fibre name
    real(kind(1d0))      :: s  !position along the orbit
    type(real_8),target  :: y(6)!input 6 dimensional function (polynomial) : Full MAP: A*YC*A_1
    real(dp)             :: v
    logical              :: set
    integer              :: i,j,k
    
    
    
    if ( .not. associated(normmoments) ) then
      return
    endif
    
!     print*, "#################################################"
!     print*, "#################################################"
!     print*, "#################################################"
!     print*, "Moments for fibre ", n,name," at ", s,"m"
     
    

    do i=1, nmoments
       
       set = .false.
       do j=1,c_%nd2
         do k = 1, moments(i)%iarray(j)
           
           if (set) then
             function_to_average = function_to_average*y(j)%t
           else
             function_to_average = y(j)%t
             set = .true.
           endif  
         enddo
       enddo

!       function_to_average=y(1)*y(1) ! just a function (taylor series)

       call cfu(function_to_average,filter,function_to_average) !cycling i.e. put the form factor to the function 
       function_to_average=function_to_average.o.damapa ! replaces x px y py ... by sigma1, sigma2, etc

       
       v = function_to_average.sub.0
       
       call double_to_table(moments(i)%table,moments(i)%column,v)

!       print*,moments(i)%iarray
!       call print(function_to_average,6)
!       print*, v
!       print*,"----------------------------------"
     
    enddo    

    call augmentcountmomtabs(s)
    
  end subroutine putmoments

  !_________________________________________________________________________________

  subroutine setemittances(emix,emiy,emiz)
    implicit none
    real(dp)         :: emix
    real(dp)         :: emiy
    real(dp)         :: emiz
    
    if (emix .lt. zero) then
      print *, "X Emittance is less then 0"
      return
    endif

    if (emiy .lt. zero) then
      print *, "Y Emittance is less then 0"
      return
    endif

    if (emiz .lt. zero) then
      print *, "Z Emittance is less then 0"
      return
    endif
    
    print*,"============================="    
    print*, emix,emiy,emiz
    
    sigmas(1) = sqrt(emix)
    sigmas(2) = sigmas(1)
    sigmas(3) = sqrt(emiy)
    sigmas(4) = sigmas(3)

    sigmas(5) = sqrt(emiz)
    sigmas(6) = sigmas(5)

    print*, sigmas
    
    
  end subroutine setemittances

  !_________________________________________________________________________________

  subroutine initmoments()
    implicit none
    integer                         :: no
    integer                         :: i ! dimension
    integer                         :: get_string 
    character(len=48), dimension(3) :: disttypes
    integer, dimension(3)           :: stringlength
    character(48)                   :: tmpstring
    real(kind(1d0))                 :: get_value
    ! This routine must be called before any init in ptc_twiss is performed
    ! since it initialize Bertz for its purpose
    !
    !
    if ( associated(normmoments) ) then
      deallocate(normmoments)
    endif
    
    
    
    i = get_value('ptc_twiss ','moments ')

    if (i .eq. 0 ) then
      call fort_warn("initmoments","Moments are not requested.")
      return
    endif


    stringlength(1) = get_string('ptc_twiss ','xdistr ',disttypes(1))
    stringlength(2) = get_string('ptc_twiss ','ydistr ',disttypes(2))
    stringlength(3) = get_string('ptc_twiss ','zdistr ',disttypes(3))
    
    no = get_value('ptc_twiss ','no ')
    if ( no <= 1 ) then
      call fort_warn('madx_ptc_distrib.f90 <initmoments>:','Order in twiss is smaller then 1')
      return
    endif
    
    allocate(normmoments(3, 0:no, 0:no))
    normmoments = zero
    
    do i=1,3
       tmpstring = disttypes(i)
       select case(tmpstring(1:stringlength(i)))
         case ('gauss')
           print*, "initmoments: Gauss distribution for dimension ", i
           call makegaus(no,i)
           distributiontype(i) = distr_gauss
         case ('flat5')
           print*, "initmoments: Flat distribution for dimension ", i
           call makeflat5(no,i)
           distributiontype(i) = distr_flat5
         case ('flat56')
           print*, "initmoments: Flat distribution for dimension ", i
           call makeflat56(no,i)
           distributiontype(i) = distr_flat56
         case default
           call fort_warn("initmoments","Distribution not recognized")
           print*, "initmoments: Distribution ", tmpstring(1:stringlength(i)), "not recognized"
           print*, "initmoments: Using default Gaussian for dimension ", i
           call makegaus(no,i)
           distributiontype(i) = distr_gauss
       end select
      
    enddo


     
    

  end subroutine initmoments
  !_________________________________________________________________________________
  subroutine allocmoments
    implicit none
    !allocates variables needed for moments calculations
    ! namely damapa and fucntion_to_avarage

      if ( .not. associated(normmoments) ) then
        return
      endif
      
      call alloc(damapa)
      damapa = sigmas
      call alloc(function_to_average)
    
  end subroutine allocmoments

  !_________________________________________________________________________________
  subroutine killmoments
    implicit none
    !cleans everything allocated in in initmoments and allocmoments
      
      if ( .not. associated(normmoments) ) then
        return
      endif
      
      deallocate(normmoments)

      call kill(damapa)
      call kill(function_to_average)

  end subroutine killmoments  
  !_________________________________________________________________________________

  real(dp) function filter(e)  
    implicit none
    integer e(:)
    integer i

    filter=one

    do i=1,c_%nd
     
     filter=filter*normmoments(i,e(2*i-1),e(2*i))
     
    enddo

  end function filter


  !_________________________________________________________________________________
  !_________________________________________________________________________________
  !_________________________________________________________________________________
!!  
   subroutine makegaus(no,d)
   implicit none
   integer no !order
   integer d ! dimension number
   integer i,j,jn(2)
   type(taylor) x,p,f
   type(Taylorresonance) fr

      print*, "Making Gauss distributions"

      call init(no,1,0,0)

      call alloc(x,p,f)
      call alloc(fr)

      x=1.d0.mono.1
      p=1.d0.mono.2

      f=one
      do i=0,no
       do j=i,no

         if(mod(i,2)/=0) cycle
         if(mod(j,2)/=0) cycle
         f=x**i*p**j

         fr=f

         jn(1)=(i+j)/2;
         jn(2)=jn(1);
         normmoments(d,i,j)=(fr%cos.sub.jn)
         normmoments(d,i,j)=normmoments(d,i,j)*singlefac(jn(1))*two**(jn(1))
         normmoments(d,j,i)=normmoments(d,i,j)
         
         print*, "mom(",i,",",j,")=",normmoments(d,i,j)
       enddo
      enddo
      call kill(x,p,f)
      call kill(fr)

   end subroutine makegaus

  !_________________________________________________________________________________

   subroutine makeflat5(no,d)
   implicit none
   integer no
   integer d !dimension number
   real(dp), allocatable :: m(:,:)
   integer i,j 

    print*, "Making flat distribution "

     do i=0,no
       do j=i,no

       if(mod(i,2)/=0) cycle

        normmoments(d,i,j)=(three)**(i/2)/(i+1) !delta assumed flat distribution and
        normmoments(d,j,i)=normmoments(d,i,j)   !and L is the delta function
        print*, "mom(",i,",",j,")=",normmoments(d,i,j)

       enddo
     enddo

   end subroutine makeflat5

  !_________________________________________________________________________________

   subroutine makeflat56(no,d)
   implicit none
   integer no
   integer d !dimension number
   real(dp), allocatable :: m(:,:)
   integer i,j 

    print*, "Making flat in delta and T distributions"

     do i=0,no,2
       do j=i,no,2

        normmoments(d,i,j)=(three)**(i/2)/(i+1) !delta assumed flat distribution and
        normmoments(d,i,j)=normmoments(d,i,j)*(three)**(j/2)/(j+1) !delta assumed flat distribution and
        
        normmoments(d,j,i)=normmoments(d,i,j)   !and L is the delta function

        print*, "mom(",i,",",j,")=",normmoments(d,i,j)

       enddo
     enddo

   end subroutine makeflat56
  !_________________________________________________________________________________
  !_________________________________________________________________________________
  !_________________________________________________________________________________

  
  !_________________________________________________________________________________
   
   real(dp) function singlefac(n)
   implicit none
   integer n,i

      singlefac=one
      do i=1,n
       singlefac=singlefac*i
      enddo

   end function singlefac

  !_________________________________________________________________________________

   integer function getdistrtype(axis)
   implicit none
   integer axis  
   
     getdistrtype =  distributiontype(axis)
   end function getdistrtype 
   

  !_________________________________________________________________________________
   
   
   real(dp) function getsigma(axis)
   implicit none
   integer axis  
     
     if ( (axis > 0) .and. (axis < 7) ) then
       getsigma = sigmas(axis)
     else
       getsigma = zero
     endif 
     
   end function getsigma

end module madx_ptc_distrib_module
!_________________________________________________________________________________
!_________________________________________________________________________________
!_________________________________________________________________________________


integer function w_ptc_getnmoments()
  use madx_ptc_distrib_module 
  implicit none
      w_ptc_getnmoments = getnmoments()
end function w_ptc_getnmoments

!_________________________________________________________________________________

subroutine  w_ptc_getmomentstabcol(n, tabn, coln)
  use madx_ptc_distrib_module 
  implicit none
  integer              :: n
  character(20)        :: tabn
  character(17)        :: coln

  call getmomentstabcol(n, tabn, coln)

end subroutine w_ptc_getmomentstabcol
