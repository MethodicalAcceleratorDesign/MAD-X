module madx_ptc_distrib_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  ! madx_ptc_distrib module
  ! Piotr K. Skowronski (CERN)
  !
  ! This module contains service for trackin

  use madx_keywords
  use madx_ptc_module
  USE madx_ptc_knobs_module
  use madx_ptc_intstate_module, only : getdebug
  implicit none

  include "madx_ptc_distrib.inc"
  save
  private

  !============================================================================================
  !  PUBLIC INTERFACE
  public                         :: momfirstinit
  public                         :: aremomentson
  public                         :: ptc_moments
  public                         :: putmoments
  public                         :: initmoments
  public                         :: allocmoments
  public                         :: killmoments
  public                         :: getdistrtype
  public                         :: setemittances
  public                         :: setsigma
  public                         :: printsigmas
  public                         :: addmoment
  public                         :: getnmoments
  public                         :: getmomentstabcol

  !============================================================================================
  !  PRIVATE
  !    data structures
  integer                        ::  firstinitdone = 0
  real(dp), pointer              ::  normmoments(:,:,:)
  real(dp)                       ::  sigmas(6)
  integer                        ::  distributiontype(3)

  integer, parameter             ::  maxnmoments=100
  type (momentdef), private      ::  moments(maxnmoments)
  integer                        ::  nmoments = 0

  type(gmap)                     :: gmapa  !da map needed to calculation initialized with sigmas
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


  subroutine addmoment(x,px,y,py,t,dp,tableIA, columnIA, parametric )
    use twissafi
    implicit none
    integer               :: x,px,y,py,dp,t
    integer               :: columnIA(*)
    integer               :: tableIA(*)
    integer               :: parametric
    integer, parameter    :: zeroasciicode = IACHAR ( '0' )
    character(48) charconv
    ! name of the column corresponds to MADX nomenclature (5 col is dp/p)
    ! iarray corresponds to

    nmoments = nmoments + 1

    moments(nmoments)%iarray(1) =  x
    moments(nmoments)%iarray(2) =  px
    moments(nmoments)%iarray(3) =  y
    moments(nmoments)%iarray(4) =  py
    moments(nmoments)%iarray(5) =  dp
    moments(nmoments)%iarray(6) =  t

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


    if (getdebug()>0) then
       print  *,"addmoment : <", moments(nmoments)%iarray(1:6) ,&
            &                   ">,<", moments(nmoments)%column,&
            &                   ">,<", moments(nmoments)%table, ">)"
    endif



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
    external             :: fort_warn
    
    if ( (n < 1) .and. (n > nmoments ) ) then

       tabn(1:1) = achar(0)
       coln(1:1) = achar(0)
       call fort_warn("getmomentstabcol","index out of range")
       return
    endif


    tabn = moments(n)%table
    coln = moments(n)%column


  end subroutine getmomentstabcol
  !____________________________________________________________________________________________

  subroutine momfirstinit()
    implicit none
     if (firstinitdone == 1) then
       return
     endif
     
     nullify(normmoments)
     firstinitdone = 1
     
  end subroutine momfirstinit

  !____________________________________________________________________________________________

  logical(lp) function aremomentson()
    implicit none

    if ( associated(normmoments) ) then
       aremomentson = my_true
    else
       aremomentson = my_false
    endif

  end function aremomentson
  !_________________________________________________________________

  subroutine ptc_moments(order)
    implicit none
    integer       :: order,mynd2,npara,nda
    integer       :: i,ii,iii
    type(real_8)  :: y(6)
    integer       :: restart_sequ,advance_node
    external      :: fort_warn, fort_info, seterrorflag, makemomentstables
    
    if (nmoments < 1) then
       call fort_info("ptc_moments","No moments specified for calculation.")
       return
    endif

    if (mapsorder < 1) then
       call seterrorflag(1,"ptc_moments",&
            "Maps are not available. Did you run ptc_twiss with savemaps=true ?")
       return
    endif

    if (.not. associated(maps)) then
       return
    endif

    if (mapsicase == 5) then
       call fort_warn("ptc_moments","For the time being the calculation of moments is not available in 5D case.")
       return
    endif

    if ((mapsicase == 5) .and. (sigmas(5) .le. 0)) then
       call fort_warn("ptc_moments","Spread in dp/p in undefined and it won't be taken taken to the account")
       print*,"In 5D case you have to specify either"
       print*,"       -  SIGE in the BEAM command or"
       print*,"       -  ET in the BEAM command  AND  BETZ with in ptc_twiss"
    endif

    if ((mapsicase == 6) .and. (sigmas(5) .le. 0)) then
       call fort_warn("ptc_moments","Spread in dp/p in undefined and it won't be taken taken to the account")
       print*,"In 6D case you have to specify longitudinal emittance SIGE in the BEAM command"
    endif


    call initmoments()
    call makemomentstables();

    nda = getnknobsall() !defined in madx_ptc_knobs

    !print*, "In moments order ", order
    mynd2 = 0
    npara = 0
    call init(default,order,nda,BERZ,mynd2,npara)

    call allocmoments()

    call alloc(y)

    iii=restart_sequ()

    do i=lbound(maps,1),ubound(maps,1)

       !       if (i == MY_RING%n) then
       !         call ptc_setdebuglevel(1)
       !       endif

       do ii=1,6
          y(ii) = maps(i)%unimap(ii)
       enddo

       call putmoments(i,maps(i)%name,maps(i)%s,y)

       iii=advance_node()

    enddo
100 continue

    call ptc_setdebuglevel(0)

    call kill(y)
    call killmoments()

  end subroutine ptc_moments
  !____________________________________________________________________________________________


  subroutine putmoments(n,name,s,y)
    implicit none
    integer              :: n !fibre number
    character(*)         :: name !fibre name
    real(dp)             :: s  !position along the orbit
    type(real_8),target  :: y(6)!input 6 dimensional function (polynomial) : Full MAP: A*YC*A_1
    real(kind(1d0))      :: v
    logical              :: set
    integer              :: i,j,k,e(6)
    integer              :: debug
    
    !    integer              :: last
    !    integer              :: index !standarf function

    !    last = index(name,'$END')
    !    if ( last == 0) then
    debug = getdebug()
    !    else
    !      debug = 10;
    !    endif

    if ( .not. associated(normmoments) ) then
       return
    endif

    if (debug > 3) then

       print*, "#################################################"
       print*, "#################################################"
       print*, "#################################################"
       print*, "Moments for fibre ", n,name," at ", s,"m"
       print*, "ND2=",c_%nd2
       print*, "NPARA=",c_%npara

       print*, "Function 1"
       call print(y(1),6)
       print*, ""
       print*, "Function 5"
       call print(y(5),6)
       print*, ""
       print*, "Function 6"
       call print(y(6),6)
       print*, ""

       print*, "GMap"
       call print(gmapa,6)
       print*, ""

    endif

    !--moments--!
    do i=1, nmoments

       function_to_average = zero
       set = .false.

       do j=1,c_%npara


          if (debug .gt. 3) write(*,'(a6,i1,a6,i1,a6,i1)',ADVANCE='NO') "nmom=",i," ndim=",j," pow=",moments(i)%iarray(j)
          do k = 1, moments(i)%iarray(j)
             if (set) then
                function_to_average = function_to_average*y(j)%t
                if (debug .gt. 3) write(*,'(a1)',ADVANCE='NO') "*"
             else
                function_to_average = y(j)%t
                set = .true.
                if (debug .gt. 3) write(*,'(a1)',ADVANCE='NO') "|"
             endif
          enddo
          if (debug .gt. 3) write(*,*) "->"

       enddo

       !       function_to_average=y(1)*y(1) ! just a function (taylor series)

       !        if (debug > 3) then
       !          print*, "function_to_average"
       !          call print(function_to_average,6)
       !        endif

       if (debug > 5) then
          print*, "Function to average"
          call print(function_to_average,6)
       endif


       call cfu(function_to_average,filter,function_to_average) !cycling i.e. put the form factor to the function

       if (debug > 4) then
          print*, "After cfu"
          call print(function_to_average,6)
       endif

       function_to_average=function_to_average.o.gmapa ! replaces x px y py ... by sigma1, sigma2, etc

       if (debug > 3) then
          print*, "Averaging (gmapped)"
          call print(function_to_average,6)

       endif

       v = function_to_average.sub.0

       if (c_%npara == 5) then
          if (debug > 3) then
              print*, v
          endif
          e = 0
          do k=1,c_%no
             e(5) = k
             if (debug > 3) then
                print*, "s^",k,"=",(sigmas(5)**(k))
                print*, "f = ", (function_to_average.sub.e)
                print*, (function_to_average.sub.e) * (sigmas(5)**(2*k))
             endif
             v = v + (function_to_average.sub.e) * (sigmas(5)**(k)) !was to 2*k
             if (debug > 9) then
                 print*, v
             endif
          enddo
       endif


       call double_to_table_curr(moments(i)%table,moments(i)%column,v)

       if (debug > 2) then
          print*, "Final  ",i," =  ", v
          print*,moments(i)%iarray
          call print(function_to_average,6)
          print*,"######################################################################################"
       endif


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


    if (getdebug() > 1) then
        print*, "Setting Sigmas (Emittances)"
    endif

    sigmas(1) = sqrt(emix)
    sigmas(2) = sigmas(1)
    sigmas(3) = sqrt(emiy)
    sigmas(4) = sigmas(3)

    sigmas(5) = sqrt(emiz)
    sigmas(6) = sigmas(5)

    if (getdebug() > 1) then
        print *, "Current sigmas setemittances ", sigmas
    endif



  end subroutine setemittances
  !_________________________________________________________________________________


  subroutine setsigma(ndim,sig)
    implicit none
    integer          :: ndim
    real(kind(1d0))  :: sig

    if (sig .lt. zero) then
       print *, "X Emittance is less then 0"
       return
    endif

    if ( (ndim .le. 0) .or. (ndim .gt. 6)) then
       print *, "Unknown dimension code", ndim
       return
    endif


    if (getdebug() > 1) then
       print *, "Setting sigma for ", ndim
       print *, "Current sigmas (setsigma) ", sigmas
    endif

    sigmas(ndim) = sig

  end subroutine setsigma
  !_________________________________________________________________________________

  subroutine printsigmas
    implicit none

    print*,"Sigmas:", sigmas
  end subroutine printsigmas
  !_________________________________________________________________________________

  subroutine initmoments()
    implicit none
    integer                         :: no
    integer                         :: i ! dimension
    integer                         :: get_string
    character(len=48), dimension(3) :: disttypes
    character(len=48)               :: cmdname
    integer, dimension(3)           :: stringlength
    character(48)                   :: tmpstring
    integer                         :: getcurrentcmdname
    external             :: fort_warn
    ! This routine must be called before any init in ptc_twiss is performed
    ! since it initialize Bertz for its purpose
    !
    !
    if ( associated(normmoments) ) then
       deallocate(normmoments)
    endif

    if (nmoments < 1) then
       !      call fort_warn("initmoments","No moments specified for calculation.")
       return
    endif


    i = getcurrentcmdname(cmdname);

    if (i .eq. 0 ) then
       call fort_warn("initmoments","Can not get the current command name.")
       return
    endif

    stringlength(1) = get_string(cmdname,'xdistr ',disttypes(1))
    stringlength(2) = get_string(cmdname,'ydistr ',disttypes(2))
    stringlength(3) = get_string(cmdname,'zdistr ',disttypes(3))

    !we take what was available in the last ptc_twiss and go to maximum order we can go
    no =  mapsorder
    if ( no < 1 ) then
       call fort_warn('madx_ptc_distrib.f90 <initmoments>:','Order in twiss is smaller then 1')
       return
    endif
    no = no*2 !we take what was available in the last ptc_twiss and go to maximum order we can go

    allocate(normmoments(3, 0:no, 0:no))
    normmoments = zero

    do i=1,3
       tmpstring = disttypes(i)
       select case(tmpstring(1:stringlength(i)))
       case ('gauss')
          if (getdebug() > 1) then
              print*, "initmoments: Gauss distribution for dimension ", i
          endif
          call makegaus(no,i)
          distributiontype(i) = distr_gauss
       case ('flat5')
          if (getdebug() > 1) then
              print*, "initmoments: Flat distribution for dimension ", i
          endif
          call makeflat5(no,i)
          distributiontype(i) = distr_flat5
       case ('flat56')
          if (getdebug() > 1) then
              print*, "initmoments: Flat distribution for dimension ", i
          endif
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
    ! namely gmapa and fucntion_to_avarage

    if ( .not. associated(normmoments) ) then
       return
    endif

    call alloc(gmapa,c_%nv)
    gmapa = sigmas
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

    call kill(gmapa)
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
       if (getdebug() > 4) then
          print*, "normmoments(",i, e(2*i-1), e(2*i),")=", normmoments(i,e(2*i-1),e(2*i))
       endif
    enddo

    if (getdebug() > 3) then

       print*,"filter(",e(1:6),")=",filter
       print*,"=================="
    endif

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

    if (getdebug() > 1) then
        print*, "Making Gauss distributions for dimension ", d
    endif

    call init(no,1,0,0)

    call alloc(x,p,f)
    call alloc(fr)

    x=one.mono.1
    p=one.mono.2

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

          if (getdebug() > 1) then
              print*, "mom(",i,",",j,")=",normmoments(d,i,j)
          endif
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
    integer i,j

    if (getdebug() > 1) then
        print*, "Making flat distribution "
    endif

    do i=0,no
       do j=i,no

          if(mod(i,2)/=0) cycle

          normmoments(d,i,j)=(three)**(i/2)/(i+1) !delta assumed flat distribution and
          normmoments(d,j,i)=normmoments(d,i,j)   !and L is the delta function

          if (getdebug() > 1) then
              print*, "mom(",i,",",j,")=",normmoments(d,i,j)
          endif

       enddo
    enddo

  end subroutine makeflat5

  !_________________________________________________________________________________

  subroutine makeflat56(no,d)
    implicit none
    integer no
    integer d !dimension number
    integer i,j

    if (getdebug() > 1) then
        print*, "Making flat in delta and T distributions"
    endif

    do i=0,no,2
       do j=i,no,2

          normmoments(d,i,j)=(three)**(i/2)/(i+1) !delta assumed flat distribution and
          normmoments(d,i,j)=normmoments(d,i,j)*(three)**(j/2)/(j+1) !delta assumed flat distribution and

          normmoments(d,j,i)=normmoments(d,i,j)   !and L is the delta function

          if (getdebug() > 1) then
              print*, "mom(",i,",",j,")=",normmoments(d,i,j)
          endif

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
