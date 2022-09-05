module madx_ptc_knobs_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  ! madx_ptc_knobs module
  ! Piotr K. Skowronski (CERN)
  !
  ! This module contains the routines for PTC knobs manipulation (aka pol_block's)
  ! MAD-X command ptc_knob configures the knobs.
  ! - addknob sets a pol_block(s) in this module
  ! - setknobs sets the configured previously knobs on a layout
  use madx_keywords
  use madx_ptc_intstate_module, only : getdebug
  use twiss0fi
  implicit none

  include "madx_ptc_knobs.inc"
  save
  private
  !============================================================================================
  !  PUBLIC INTERFACE
  public                         :: putusertable
  public                         :: addpush
  public                         :: cleartables

  public                         :: addknob   !adds knobs (called from c)
  public                         :: addknobi
  public                         :: resetknobs ! removes knobs
  public                         :: setknobs ! sets added knobs on a layout
  public                         :: getnknobsm ! returns number of knobs that are magnets properties (set via pol_block)
  public                         :: getnknobis ! returns number if knobs that are initial condition at the beginning of a line
  public                         :: getnknobsall ! returns number of all knobs getnknobsm+getnknobis
  public                         :: resultswithknobs ! called by ptc_twiss to dump the requested results
  public                         :: parametrictwiss  ! returns lattice functions with dependence on knobs
  public                         :: getknobinicond
  public                         :: getnpushes ! returns number of selected parameters to be stored

  public                         :: writeparresults
  public                         :: killparresult
  public                         :: finishknobs
  public                         :: twissfctn
  public                         :: getknobsnames
  public                         :: getfctnsnames
  public                         :: getfunctionat
  public                         :: getfunctionvalueat
  public                         :: getlengthat
  public                         :: setparvalue
  public                         :: setknobvalue
  public                         :: filltables
  public                         :: setnameintbles

  !============================================================================================
  !  PRIVATE
  !    data structures

  integer, parameter                    ::  maxnpushes=100
  type (tablepush_poly), private        ::  pushes(maxnpushes)
  integer                               ::  npushes = 0
  character(20), private                ::  tables(maxnpushes)  !tables names of existing pushes - each is listed only ones
  integer,  private                     ::  ntables = 0 !number of distinctive tables
  real(dp), private, allocatable        ::  Dismom(:,:)    ! <xnormal_(2*i-1)**(2j)>= dismon(i,j)*I_i**j


  integer, parameter                    ::  maxnpolblocks=20
  type (pol_block), target              ::  polblocks(maxnpolblocks)
  integer                               ::  npolblocks=0
  integer                               ::  nknobs=0
  type(pol_block_inicond)               ::  knobi
  integer                               ::  nknobi


  integer, parameter                    ::  maxpar = 100
  integer                               ::  nmapels = 0
  !  type(mapelresult),target              ::  mapels(maxpar)

  real(dp), allocatable                 ::  spos(:)
  real(dp), allocatable                 ::  deltaes(:) !array with energy increase for each element with respect to the beginning
  real(dp), allocatable                 ::  parvals(:) ! temp array with parameter values, to find out a function value for some parameters
  !  type(taylor), allocatable             ::  results(:,:)
  type(universal_taylor), allocatable, target   ::  results(:,:)

  character(48)                         ::  twisstablename

  integer, parameter                    ::  textbufferlength = 100000
  character(textbufferlength)           ::  textbuffer


  integer                               ::  currentrow = 0

  type(real_8), private, dimension(3)   ::  testold
  type(real_8), private, dimension(3)   ::  phase

  integer, private, allocatable         ::  E(:)  !array to pick taylor coefficients with .sub. and .par.
  type(taylor), private                 ::  ave(6,6,3)
  type(real_8), private                 ::  tx, ty
  type(real_8), private                 ::  dph
  type(real_8), private                 ::  test
  integer, private                      ::  taylorsallocated = 0

  !  integer, private, dimension(6)        :: j1 = (/1,0,0,0,0,0/)
  !  integer, private, dimension(6)        :: j2 = (/0,1,0,0,0,0/)
  !  integer, private, dimension(6)        :: j3 = (/0,0,1,0,0,0/)
  !  integer, private, dimension(6)        :: j4 = (/0,0,0,1,0,0/)
  integer, private, dimension(6)        :: j5 = (/0,0,0,0,1,0/)
  integer, private, dimension(6)        :: j6 = (/0,0,0,0,0,1/)
  integer, private, dimension(6,6)       :: fo = &
       reshape(          (/1,0,0,0,0,0,&
       0,1,0,0,0,0,&
       0,0,1,0,0,0,&
       0,0,0,1,0,0,&
       0,0,0,0,1,0,&
       0,0,0,0,0,1 /), &
       (/6,6/) )

  !============================================================================================
  !  PRIVATE
  !    routines
  private                               :: augment_counts
  private                               :: issuchtableexist
  private                               :: putnameintables
  private                               :: allocparresult






contains
  !____________________________________________________________________________________________

  function getknobinicond()
    implicit none
    type(pol_block_inicond) :: getknobinicond

    getknobinicond%beta = knobi%beta
    getknobinicond%alfa = knobi%alfa
    getknobinicond%dispersion = knobi%dispersion

  end function getknobinicond

  !____________________________________________________________________________________________


  subroutine nullify(inpbic)
    implicit none
    type(pol_block_inicond)  :: inpbic

    inpbic%alfa = 0
    inpbic%beta = 0
    inpbic%dispersion = 0

  end subroutine nullify

  subroutine putusertable(n,name,s,denergy,y,A)!here y is transfermap and A is a_ (in twiss it is called y)
    !puts the coefficients in tables as defined in array pushes
    implicit none
    integer              :: n !fibre number
    character(*)         :: name !fibre name
    real(kind(1d0))      :: s  !position along the orbit
    real(dp)             :: denergy
    type(real_8),target  :: y(6)!input 6 dimensional function (polynomial) : Full MAP: theTransferMap 
    type(real_8),target  :: A(6)!input 6 dimensional function (polynomial) : A_script
    type(real_8),pointer :: e !element in array
    real(kind(1d0))      :: coeff
    integer              :: i !iterator
    logical              :: pblockson

    if (getdebug() > 3) then
       print *,"madx_ptc_tablepush :putusertable n ",n," name ", name
    endif


    if ((getnknobsall() > 0) .and. (currentrow > 0)) then
       !if there are any knobs and knobs where initialized (i.e. currentrow > 0)
       ! otherwise the results table is not allocated
       pblockson = .true.  !it means there are parameters on
       spos(currentrow) = s
       deltaes(currentrow) = denergy
       call parametric_coord(A)
       call parametrictwiss(A)
    else
       pblockson = .false.
    endif
    
    if (name /= '$end$ ') then
      call putnameintables()
    else
      call setnameintbles('$end$ ')!twiss gets names automatically from C element iterator
      if (getdebug()>3) then
        print*,"PTC_NORMAL mode"
      endif
    endif  
    
    !if (present(y) == .false.) then
    !  return
    !endif
      
    if (c_%npara == 6) then
       call sixdmode()
    else
       do i=1,npushes

          e => y(pushes(i)%element)

          if (pushes(i)%pushtab) then
             !save value defined by ptc_select in the user TFS table
             coeff  = e.sub.(pushes(i)%monomial)
             if (getdebug()>3) then
                write(6,'(a13, a10, a3, f9.6, a10, i1, 5(a13), i3)') &
                     &        "4or5D putting",pushes(i)%monomial,"=",coeff," arr_row ", pushes(i)%element,&
                     &        " in table ", pushes(i)%tabname," at column ", pushes(i)%colname, &
                     &        " for fibre no ",n
             endif

             call double_to_table_curr(pushes(i)%tabname, pushes(i)%colname, coeff);
          endif

          if ( pblockson .and. (pushes(i)%index > 0) ) then
             !save polynomial in function of knobs in the results array
             results(currentrow,pushes(i)%index) = e.par.(pushes(i)%monomial)

             if (getdebug()>3) then
                write(6,*) &
                     &        "4or5D putting ",pushes(i)%monomial," arr_row ", pushes(i)%element,&
                     &        " at named ", pushes(i)%colname, &
                     &        " for fibre no ",n
                write(6,*) "currentrow is ", currentrow," index ",pushes(i)%index
                call print(results(currentrow,pushes(i)%index),6)
             endif
          endif

       enddo
    endif

    call augment_counts()

    if (currentrow > 0) then
       !if currnet row == 0 it means that knobs were not initialized
       currentrow = currentrow + 1
    endif

    !____________________________________________________________________________________________
  contains
    !____________________________________________________________________________________________
    subroutine sixdmode()
      implicit none
      character bufchar
      character*(6) monstr

      do i=1,npushes

         if (pushes(i)%element < 5) then
            e => y(pushes(i)%element)
         else
            if (pushes(i)%element == 5) then
               e => y(6) !6th coordinate  is d(cT) or cT
            else
               e => y(5) !5th coordinate is dp/p
            endif
         endif

         monstr = pushes(i)%monomial
         bufchar = monstr(5:5)
         monstr(5:5) = monstr(6:6)
         monstr(6:6) = bufchar

         if (pushes(i)%pushtab) then
            coeff = e.sub.monstr
            if ((mod(iachar(monstr(6:6))-iachar('0'),2).eq.1).neqv.(pushes(i)%element.eq.5)) coeff = -coeff

            if (getdebug()>3) then
               write(6,'(a21, a10, a3, f9.6, a10, i1, 5(a13), i3)') &
                    &        "6D mode Put in table ",pushes(i)%monomial,"=",coeff," arr_row ", pushes(i)%element,&
                    &        " in table ", pushes(i)%tabname," at column ", pushes(i)%colname, &
                    &        " for fibre no ",n
            endif

            call double_to_table_curr(pushes(i)%tabname, pushes(i)%colname, coeff);
         endif

         if ( pblockson .and. (pushes(i)%index > 0) ) then

            results(currentrow,pushes(i)%index) = e.par.monstr


            if (getdebug()>3) then
               write(6,*) &
                    &        "6Dmode parametric ",pushes(i)%monomial," arr_row ", pushes(i)%element,&
                    &        " at named ", pushes(i)%colname, &
                   &        " for fibre no ",n
               write(6,*) "currentrow is ", currentrow," index ",pushes(i)%index
               call print(results(currentrow,pushes(i)%index),6)
            endif

         endif

      enddo
    end subroutine sixdmode

    !____________________________________________________________________________________________

  end subroutine putusertable
  !____________________________________________________________________________________________

  subroutine fillusertables()
    implicit none
    integer                         :: i,e
    integer                         :: nelems
    real(kind(1d0))                 :: coeff
    type(universal_taylor), pointer :: t

    call cleartables()

    nelems = ubound(results,1)
    do e=1, nelems
       do i=1,npushes
          if ( (pushes(i)%pushtab) .and. (pushes(i)%index > 0) ) then
             t => results(e,pushes(i)%index)
             coeff = gettaylorvalue(t)
             call double_to_table_curr(pushes(i)%tabname, pushes(i)%colname, coeff);
          endif
       enddo

       call augment_counts()

    enddo
  end subroutine fillusertables
  !____________________________________________________________________________________________

  subroutine filltwisstable()
    !puts the coefficients in tables as defined in array pushes
    implicit none
    integer              :: i,ii !iterator
    integer              :: nelems !iterator
    integer, parameter   :: n = 6 !counts
    integer, parameter   :: fillntwisses  = disp4 - beta11 + 1
    integer, parameter   :: ntwissesde = gama33 - beta11 + 1
    real(kind(1d0))      :: opt_fun(ntwisses) ! opt_fun(fundim) does not allocate enough space
                                              ! -> invoke undefined behavior in few places below.
    type(universal_taylor), pointer :: t

    if (.not. ALLOCATED(results)) then
       return
    endif

    call reset_count(twisstablename)

    nelems = ubound(results,1)

    if (nelems .lt. (currentrow-1)) then
       call fort_warn("filltwisstable",&
            "It seems the last ptc_twiss has failed")
       nelems = currentrow - 1
    endif



    do i=1, nelems

       do ii=beta11,ntwisses
          t => results(i,ii)
          opt_fun(ii)=  gettaylorvalue(t)
       enddo

       do ii=beta11,ntwissesde
          opt_fun(ii) = opt_fun(ii)*deltaes(i)
       enddo

       call vector_to_table_curr(twisstablename, 'beta11 ', opt_fun(beta11), fillntwisses)
       call vector_to_table_curr(twisstablename, 'x ', opt_fun(kn_x), n)

       call augmentcountonly(twisstablename)

    enddo

  end subroutine filltwisstable

  !____________________________________________________________________________________________

  subroutine inittables()
    implicit none

    !     print*, " no=",c_%no
    !     print*, " nv=",c_%nv
    !     print*, " nd=",c_%nd
    !     print*, " nd2=",c_%nd2
    !     print*, " ndpt=",c_%ndpt
    !     print*, " npara_fpp=",c_%npara_fpp
    !     print*, " npara=",c_%npara
    !     print*, " np_pol=",c_%np_pol
    !     print*, " setknob=",c_%setknob
    !

    deallocate(dismom)
    allocate(dismom(c_%nd,0:c_%no/2))

  end subroutine inittables
  !____________________________________________________________________________________________

  subroutine cleartables()
    implicit none
    integer  :: i ! iterator
    ! we do not deal here with the main twiss table, hence we do not clear it
    do i=1,ntables
       if (getdebug()>3) then
           print *,"Clearing ",tables(i)
       end if
       call reset_count(tables(i))
    enddo

  end subroutine cleartables
  !____________________________________________________________________________________________

  subroutine augment_counts()
    implicit none
    integer  :: i ! iterator

    do i=1,ntables
       if (getdebug()>3) then
           print *,"Augmenting ",tables(i)
       end if
       call augmentcountonly(tables(i)) !we need to use special augement,
       !cause the regular one looks for for
       !variables names like columns to fill the table
    enddo

  end subroutine augment_counts
  !____________________________________________________________________________________________
  subroutine setnameintbles(name)
    implicit none
    character(*)         :: name !fibre name
    integer       :: i ! iterator
    
    do i=1,ntables
       if (getdebug()>2) then
           print *,"Putting name in ",tables(i)
       end if
       call string_to_table_curr(tables(i),"name ",name)
    enddo

  end subroutine setnameintbles
  !____________________________________________________________________________________________
  subroutine putnameintables()
    implicit none
    integer       :: i ! iterator
    
    do i=1,ntables
       if (getdebug()>2) then
           print *,"Putting name in ",tables(i)
       end if
       !puts automatically name of the current element
       !does not work for ptc_normal
       call string_to_table_curr(tables(i),"name ","name ")
    enddo

  end subroutine putnameintables
  !____________________________________________________________________________________________


  subroutine addpush(table,column,element,monomial)
    use twissafi
    implicit none
    integer   table(*)
    integer   column(*)
    integer   element
    integer   monomial(*)
    logical   addtable
    logical   parametric
    real(kind(1d0))            :: get_value
    character(48) charconv


    parametric = get_value('ptc_select ','parametric ') .ne. 0
    if ( (parametric .eqv. .false.) .and. (table(1) == 0) ) then
       call fort_warn("addpush","Neither table specified neither parametric value. Ignoring")
       return
    endif

    npushes = npushes + 1

    pushes(npushes)%tabname = charconv(table)
    pushes(npushes)%colname = charconv(column)

    pushes(npushes)%element = element
    pushes(npushes)%monomial = charconv(monomial)
    !imput "int string"  has the length of the string at the first plave
    pushes(npushes)%tabname(table(1)+1:table(1)+1)=achar(0)
    pushes(npushes)%colname(column(1)+1:column(1)+1)=achar(0)
    pushes(npushes)%monomial(monomial(1)+1:monomial(1)+1)=achar(0)

    if (column(1) == 0) then
       write(pushes(npushes)%colname,'(i1,a1,a6)') element,"_",pushes(npushes)%monomial
    endif


    if (table(1) > 0) then
       pushes(npushes)%pushtab = .true.
       addtable = .not. issuchtableexist(pushes(npushes)%tabname) !add to table to the list onl

    else
       pushes(npushes)%pushtab = .false.
       addtable = .false. !no table name
    endif

    if (parametric) then
       nmapels = nmapels + 1
       pushes(npushes)%index = ntwisses+nmapels
    else
       pushes(npushes)%index = 0
    endif


    if (getdebug()>3) then
       print  *,"madx_ptc_tablepush : addpush(<",&
            &          pushes(npushes)%element,">,<",pushes(npushes)%monomial,">)"
       print  *,"madx_ptc_tablepush : colname <",pushes(npushes)%colname,">"
       print  *,"madx_ptc_tablepush : parametric results index ", pushes(npushes)%index
       if (pushes(npushes)%pushtab) then
          print  *,"madx_ptc_tablepush : table <",pushes(npushes)%tabname,">"
       else
          print  *,"madx_ptc_tablepush : not pushing to table"
       endif
    endif


    if ( addtable) then
       ntables = ntables + 1
       tables(ntables) = pushes(npushes)%tabname
       if (getdebug()>3)  then
          print *,"Table has been added to the tables list ", tables(ntables)
       endif

    endif

  end subroutine addpush
  !____________________________________________________________________________________________


  logical(lp) function issuchtableexist(tname)
    implicit none
    character(20) :: tname !name of the table to be checked if already is listed in table names array
    integer       :: i! iterator

    issuchtableexist = .false.

    do i=1, ntables
       if (tables(i) == tname) then
          issuchtableexist = .true.
          return
       endif
    enddo

  end function issuchtableexist


  !____________________________________________________________________________________________
  subroutine addknobi(nameIA)
    use twissafi
    implicit none
    integer     :: nameIA(*)
    character(48)              :: name
    character(48) charconv

    if (nknobi >= maxnpolblocks) then
       call fort_warn("addknob","Can not add more knobs, array with initial knobs if full")
       return
    endif

    name = charconv(nameIA)
    select case(name(1:6))
    case ('beta11')
       knobi%beta(1) = nknobi+1
    case ('beta22')
       knobi%beta(2) = nknobi+1
    case ('beta33')
       knobi%beta(3) = nknobi+1
    case ('alfa11')
       knobi%alfa(1) = nknobi+1
    case ('alfa22')
       knobi%alfa(2) = nknobi+1
    case ('alfa33')
       knobi%alfa(3) = nknobi+1
    case ('disp1')
       knobi%dispersion(1) = nknobi+1
    case ('disp2')
       knobi%dispersion(2) = nknobi+1
    case ('disp3')
       knobi%dispersion(3) = nknobi+1
    case ('disp4')
       knobi%dispersion(4) = nknobi+1

    case default
       print*, "Variable not recognized"
       print*,"parameter ",name,"not recognized"
       return
    end select

    nknobi = nknobi + 1


  end subroutine addknobi
  !____________________________________________________________________________________________
  subroutine addknob(fibrenameIA)
    use twissafi
    implicit none
    integer     :: fibrenameIA(*)
    character(48)              :: fibrename
    integer     :: nint, ndble, k, int_arr(nmax), char_l(nmax), i
    real(kind(1d0))            :: d_arr(nmax)
    character(400)             :: char_a
    type (pol_block), pointer  :: pb
    logical(lp) :: exactmatch
    real(kind(1d0))            :: get_value
    character(48) charconv


    if (npolblocks >= maxnpolblocks) then
       call fort_warn("addknob","Can not add more knobs, array with pol_blocks if full")
       return
    endif

    npolblocks = npolblocks + 1
    pb => polblocks(npolblocks)
    pb = 0

    fibrename = charconv(fibrenameIA)

    k = index(fibrename,':')
    if (k > 0) then
       pb%name = fibrename(1:k-1)
    else
       pb%name = fibrename
    endif

    if (getdebug() > 1) then
        print *,"addknob: pb%name is ", pb%name," npolblocks=",npolblocks
    end if


    exactmatch = get_value('ptc_knob ','exactmatch ') .ne. 0
    if (exactmatch) then
       if (getdebug() > 1) then
           print*,"addknob: Using Exact name match: ", fibrename
       end if
       pb%vorname = fibrename
    else
       pb%n_name =  len_trim(pb%name)
       if (getdebug() > 1) then
          print*,"addknob: Using Not Exact name match:"
          print*,"    all elements starting with ", pb%name
          print*,"    number of first letters    ", pb%n_name
       endif
    endif

    call comm_para('kn ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
    if (getdebug()>2) then
        print*, "there is ",nint, " kn's: ", int_arr(1:nint)
    end if
    do i = 1, nint
       if (int_arr(i) < 0) then
          exit
       endif
       int_arr(i) = int_arr(i) + 1 !madx numerates from 0
       nknobs = nknobs + 1
       pb%ibn(int_arr(i)) = nknobs
       if (getdebug()>0) then
           print*, "Set normal mulitpole component ", int_arr(i),"as ", nknobs, "parameter of PTC"
       end if
    enddo

    call comm_para('ks ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
    if (getdebug()>2) then
        print*, "there is ",nint, " ks's: ", int_arr(1:nint)
    end if
    do i = 1, nint
       if (int_arr(i) < 0) then
          exit
       endif
       int_arr(i) = int_arr(i) + 1 !madx numerates from 0
       nknobs = nknobs + 1
       pb%ian(int_arr(i)) = nknobs
       if (getdebug()>0) then
           print*, "Set skew mulitpole component ", int_arr(i)," as ", nknobs, "parameter of PTC"
       end if
    enddo

  end subroutine addknob
  !____________________________________________________________________________________________

  subroutine setknobs(alayout)
    use twissafi
    implicit none
    type(layout),target      :: alayout
    type (pol_block), pointer  :: pb
    integer     :: i,j,k
    character(48) charconv

!    if (.not. associated(alayout)) then
!       call fort_warn("setknobs","Passed pointer to a layout is not associated")
!       return
!    endif

    twisstablename = table_name !defined in twissa.fi

    if (getdebug() > 2 ) then
       print*, "setknobs: There are ", npolblocks, "pol_blocks"
    endif

    if ((npolblocks == 0) .and. (nknobs == 0) .and. (nknobi == 0)) then
       if (ALLOCATED(results)) call killparresult()
       currentrow = -1
       return
    endif

    do i=1, npolblocks
       pb => polblocks(i)
       if (getdebug() > 2 ) then
          print*, "setknobs: Setting pol_block ", i," for ",pb%name
       endif
       alayout = pb
    enddo

    if (ALLOCATED(results)) call killparresult()

    call allocparresult(alayout%n)

    do i=1,c_%nd2
       do j=1,c_%nd2
          do k=1,c_%nd
             call alloc(ave(i,j,k))
          enddo
       enddo
    enddo

    call alloc(phase)
    call alloc(testold)
    call alloc(test)
    call alloc(dph)
    call alloc(tx, ty)

    taylorsallocated = 1

    currentrow = 1

    !    print*, "setknobs: All pol_blocks set"

  end subroutine setknobs
  !_________________________________________________________________________________
  subroutine allocparresult(n)
    implicit none
    integer       :: n
    integer       :: nr, np
    integer       :: i,j

    np = c_%nv - c_%npara
    if (np <= 0) then
       call fort_warn("addpush","Number of parameters is 0")
       currentrow = -1  ! this is the signal that initialization has failed
       return;
    endif

    allocate(parvals(np))
    parvals(:) = zero

    nr = ntwisses+nmapels
    allocate(results(n,nr))
    do i=1,n
       do j=1,nr
          results(i,j)=0
       enddo
    enddo
    allocate(spos(n))
    allocate(deltaes(n))
    allocate(e(c_%npara_fpp))


  end subroutine allocparresult
  !_________________________________________________________________________________

  subroutine killparresult
    implicit none
    integer       :: i,j

    if (.not. ALLOCATED(results)) then
       return
    endif

    !    remnestant of the taylor type
    if(getdebug() > 2) then
       print*, "killparresult: Shape of the current array: "
       print*, "1D",lbound(results,1),ubound(results,1)
       print*, "2D",lbound(results,2),ubound(results,2)
    endif

    do i=lbound(results,1),ubound(results,1)
       do j=lbound(results,2),ubound(results,2)
          !          print*,i,j
          results(i,j) = -1
          !           call kill(results(i,j))
       enddo
    enddo

    deallocate(spos)
    deallocate(deltaes)
    deallocate(results)
    deallocate(parvals)
    deallocate(e)

    currentrow = 0

  end subroutine killparresult
  !____________________________________________________________________________________________

  subroutine resetknobs()
    implicit none
    integer     :: i

    call nullify(knobi)
    nknobi = 0

    do i=1, npolblocks
       polblocks(i) = 0
    enddo
    nknobs = 0
    npolblocks = 0
  end subroutine resetknobs
  !____________________________________________________________________________________________

  function getnknobsm()
    implicit none
    integer     :: getnknobsm

    getnknobsm = nknobs

  end function getnknobsm
  !____________________________________________________________________________________________

  function getnknobis()
    implicit none
    integer     :: getnknobis

    getnknobis = nknobi

  end function getnknobis
  !____________________________________________________________________________________________

  function getnknobsall()
    implicit none
    integer     :: getnknobsall

    getnknobsall = nknobs + nknobi

  end function getnknobsall
  !____________________________________________________________________________________________

  function getnpushes()
    implicit none
    integer     :: getnpushes

    getnpushes = npushes

  end function getnpushes
  !____________________________________________________________________________________________

  subroutine getfctnsnames
    implicit none
    character*(100) ::name
    character*(8) ::twname
    integer i,j
    integer last

    do i=1, ntwisses
       twname = twissnames(i) ! twissname now known to be of fixed length 7
       twname(8:8) = achar(0)
       call madxv_setfctnname(i,twname)
    enddo

    do i=ntwisses, ntwisses+nmapels
       do j=1,npushes
          if (pushes(j)%index /= i ) cycle

          name = pushes(j)%colname

          last = len_trim(name) + 1
          if (last > 100) last = 100
          name(last:last) = achar(0)

          call madxv_setfctnname(i, name)

       enddo

    enddo


  end subroutine getfctnsnames
  !____________________________________________________________________________________________

  subroutine getknobsnames()
    implicit none
    character*(100) ::name
    character*(48)  ::pname
    integer last
    integer i,j,n

    n = 0

    !magnets properties
    do i=1, npolblocks
       do j=1, nmax
          if (polblocks(i)%ian(j) /= 0 ) then
             n = n + 1

             if (polblocks(i)%vorname == ' ') then
                write(name,'(A1,i2.2,A2,A,A,I2)') "p",polblocks(i)%ian(j),": ",& !convertion to madx counting
                     polblocks(i)%name(1:len_trim(polblocks(i)%name)), "* skew ",j-1
             else
                write(name,'(A1,i2.2,A2,A,A,I2)') "p",polblocks(i)%ian(j),": ",& !convertion to madx counting
                     polblocks(i)%vorname(1:len_trim(polblocks(i)%vorname)), " skew ",j-1
             endif
             print*, "n=",n, name

             last = len_trim(name) + 1
             if (last > 100) last = 100
             name(last:last) = achar(0)

             call madxv_setknobname(n,name)
          endif

          if (polblocks(i)%ibn(j) /= 0 ) then
             n = n + 1
             if (polblocks(i)%vorname == ' ') then
                write(name,'(A1,i2.2,A2,A,A,I2)') "p",polblocks(i)%ibn(j),": ",& !convertion to madx counting
                     polblocks(i)%name(1:len_trim(polblocks(i)%name)) , "* normal ",j-1
             else
                write(name,'(A1,i2.2,A2,A,A,I2)') "p",polblocks(i)%ibn(j),": ",& !convertion to madx counting
                     polblocks(i)%vorname(1:len_trim(polblocks(i)%vorname)) , " normal ",j-1
             endif
             print*, "n=",n, name

             last = len_trim(name) + 1
             if (last > 100) last = 100
             name(last:last) = achar(0)
             call madxv_setknobname(n,name)
          endif

       enddo
    enddo

    !initial conditions
    do i=1, nknobi
       n = n + 1
       pname = 'NULL'
       if (knobi%beta(1) == i) pname='initial beta11'
       if (knobi%beta(2) == i) pname='initial beta22'
       if (knobi%beta(3) == i) pname='initial beta33'
       if (knobi%alfa(1) == i) pname='initial alfa11'
       if (knobi%alfa(2) == i) pname='initial alfa22'
       if (knobi%alfa(3) == i) pname='initial alfa33'
       if (knobi%dispersion(1) == i) pname='initial disp1'
       if (knobi%dispersion(2) == i) pname='initial disp2'
       if (knobi%dispersion(3) == i) pname='initial disp3'
       if (knobi%dispersion(4) == i) pname='initial disp3'

       if (pname(1:4) == 'NULL') then
          print*,"Knob of Initial Condition ", i," can not be recognized"
          cycle
       endif

       write(name,'(A1,i2.2,A2,A15)') "p",n,": ",pname
       last = len_trim(name) + 1
       if (last > 48) last = 48
       name(last:last) = achar(0)
       call madxv_setknobname(n,name)
    enddo


  end subroutine getknobsnames
  !____________________________________________________________________________________________

  subroutine getfunctionat( el, n)
    implicit none
    integer n, el
    integer eqlength

    if (.not. ALLOCATED(results)) then
       return
    endif


    if ( (el < lbound(results,1)) .or. (el > ubound(results,1)) ) then
       return
    endif

    if ( (n < lbound(results,2)) .or. (n > ubound(results,2)) ) then
       return
    endif

    call getpareq(results(el,n),textbuffer)

    eqlength = ( LEN_TRIM(textbuffer) + 1 )

    if (eqlength > textbufferlength) eqlength = textbufferlength
    textbuffer(eqlength:eqlength) = achar(0)
    call madxv_setfunctionat(el, n, textbuffer)


  end subroutine getfunctionat
  !____________________________________________________________________________________________



  !____________________________________________________________________________________________

  function getfunctionvalueat( n, el)
    implicit none
    real(dp)                        :: getfunctionvalueat
    integer                         :: n, el
    type(universal_taylor), pointer :: t

    if (.not. ALLOCATED(results)) then
       return
    endif

    if ( (el < lbound(results,1)) .or. (el > ubound(results,1)) ) then
       return
    endif

    if ( (n < lbound(results,2)) .or. (n > ubound(results,2)) ) then
       return
    endif

    getfunctionvalueat = zero

    t => results(el,n)

    if (.not. associated(t)) then
       print*,"Getting function ",n," at el ",el
       getfunctionvalueat = zero
    else
       getfunctionvalueat = gettaylorvalue(t)
    endif


  end function getfunctionvalueat
  !____________________________________________________________________________________________

  function gettaylorvalue(ut)
    implicit none
    real(dp)              :: gettaylorvalue
    type(universal_taylor), pointer :: ut
    integer i,ii,np
    real(dp)   ::  c,p


    if (.not. associated(ut)) then
       call fort_warn("gettaylorvalue",&
            "provided taylor is not associated")

       gettaylorvalue = zero
       return
    endif

    if (.not. associated(ut%n)) then
       call fort_warn("gettaylorvalue",&
            "provided taylor is not associated")

       gettaylorvalue = zero
       return
    endif

    if (.not. allocated(parvals)) then
       call fort_warn("gettaylorvalue",&
            "array with parameter values is not allocated")

       gettaylorvalue = zero
       return
    endif


    if(ut%n == 0) then
       if (getdebug() > 3) then
          print*,"no coefficients in the taylor"
       endif

       gettaylorvalue = zero
       return
    endif

    gettaylorvalue = zero

    do i = 1,ut%n

       c = ut%c(i)
       !       print*, "coef",i," is ",c

       do ii = c_%npara + 1, ut%nv

          if (ut%j(i,ii) /= 0) then

             np =  ii - c_%npara
             p = parvals(np)
             !           print*, "par",np," is ",p
             if (ut%j(i,ii) /= 1) then
                p = p**ut%j(i,ii)
                !             print*, "par",np," to power",ut%j(i,ii), " is ",p
             endif

             c = c*p
             !           print*,"coef",i, "X par",np," to power",ut%j(i,ii), " is ",c

          endif

       enddo
       !       print*,"COEF",i, " is ",c
       gettaylorvalue = gettaylorvalue + c
       !       print*,"FUNCTION after ", i," is ",gettaylorvalue

    enddo



  end function gettaylorvalue

  !____________________________________________________________________________________________
  !____________________________________________________________________________________________
  !____________________________________________________________________________________________

  subroutine twissfctn(y,ave)   !  Computes <x_i x_j> assuming linearity and no parameters
    implicit none
    type(real_8) y(6)
    real(dp) ave(6,6,3)
    integer e(6)
    integer i,j,k

    print*,"c_%nd2 is ", c_%nd2
    print*,"c_%nd is ", c_%nd

    e=0
    ave=zero
    do i=1,c_%nd2
       do j=i,c_%nd2
          do k=1,c_%nd
             e(k*2-1)=1
             ave(i,j,k)=ave(i,j,k)+(y(i).sub.e)*(y(j).sub.e)
             e(k*2-1)=0
             e(k*2)=1
             ave(i,j,k)=ave(i,j,k)+(y(i).sub.e)*(y(j).sub.e)
             e(2*k)=0
             ave(j,i,k)=ave(i,j,k)

             !              if (ave(i,j,k) /= zero) then
             !                write(6,'(a3,i1,i1,i1,a3,f16.8)') 'ave', i,j,k,' = ',ave(i,j,k)
             !              endif
          enddo
       enddo
    enddo


    !    print*, "Beta11", ave(1,1,1)
    !    print*, "Alfa11", ave(1,1,2)
    !    print*, "Gama11", ave(1,1,3)
    !
    !    print*, "Beta12", ave(3,3,1)
    !    print*, "Beta12", ave(3,3,2)
    !    print*, "Beta12", ave(3,3,3)
    !
    !    print*, "Beta13", ave(5,5,3)
    !
    !    print*, "Beta21", ave(1,3,1)
    !    print*, "Beta22", ave(1,5,1)
    !    print*, "Beta23", ave(3,5,1)
    !
    !    print*, "Beta31", ave(3,1,1)
    !    print*, "Beta32", ave(5,1,1)
    !    print*, "Beta33", ave(5,5,1)


  end subroutine twissfctn

  !_________________________________________________________________________________
  subroutine parametric_coord(y)   !  Computes <x_i x_j> assuming linearity and with parameters
    implicit none
    type(real_8) y(6)

    e(:)=0

    results(currentrow, kn_x)  = y(1).par.e
    results(currentrow, kn_px) = y(2).par.e

    results(currentrow, kn_y)  = y(3).par.e
    results(currentrow, kn_py) = y(4).par.e

    if (c_%npara_fpp == 5) then
       results(currentrow, kn_dp) = y(5).par.e
       results(currentrow, kn_t)  = zero
    else if(c_%npara_fpp == 6) then
       results(currentrow, kn_dp) = y(5).par.e
       results(currentrow, kn_t)  = y(6).par.e
    else
       results(currentrow, kn_dp) = zero
       results(currentrow, kn_t)  = zero
    endif

  end subroutine parametric_coord

  !_________________________________________________________________________________
  subroutine parametrictwiss(y)   !  Computes <x_i x_j> assuming linearity and with parameters
    implicit none
    type(real_8) y(6)
    real(dp) :: epsil=1e-12
    integer i,j,k

    e=0
    do i=1,c_%nd2
       do j=i,c_%nd2
          do k=1,c_%nd
             ave(i,j,k)=zero
             e(k*2-1)=1
             ave(i,j,k)= ave(i,j,k) + (y(i)%t.par.e)*(y(j)%t.par.e)
             e(k*2-1)=0
             e(k*2)=1
             ave(i,j,k)=morph(ave(i,j,k))+ (y(i).par.e)*(y(j).par.e) !line * does the same, here taylor is morphed to polimorph,
             e(2*k)=0                                                !and above taylor component of polimorph is used explicitely
             ave(j,i,k)=ave(i,j,k)

          enddo
       enddo
    enddo
    
    ave(1,2,1) = -ave(1,2,1)
    results(currentrow, alfa11) = ave(1,2,1)
    results(currentrow, beta11) = ave(1,1,1)
    results(currentrow, gama11) = ave(2,2,1)

    ave(3,4,2) = -ave(3,4,2)
    results(currentrow, alfa22) = ave(3,4,2)
    results(currentrow, beta22) = ave(3,3,2)
    results(currentrow, gama22) = ave(4,4,2)

    !----------------

    ave(1,2,2) = -ave(1,2,2)
    results(currentrow, alfa12) = ave(1,2,2)
    results(currentrow, beta12) = ave(1,1,2)
    results(currentrow, gama12) = ave(2,2,2)

    ave(3,4,1) = -ave(3,4,1)
    results(currentrow, alfa21) = ave(3,4,1)
    results(currentrow, beta21) = ave(3,3,1)
    results(currentrow, gama21) = ave(4,4,1)


    !----------------

    if ( c_%nd == 3) then
       ave(1,2,3) = -ave(1,2,3)
       results(currentrow, alfa13) = ave(1,2,3) !-
       results(currentrow, beta13) = ave(1,1,3)  !-
       results(currentrow, gama13) = ave(2,2,3)  !-

       ave(3,4,3) =-ave(3,4,3)
       results(currentrow, alfa23) = ave(3,4,3) !-
       results(currentrow, beta23) = ave(3,3,3)  !-
       results(currentrow, gama23) = ave(4,4,3)  !-

       ave(5,6,3) = -ave(5,6,3)
       results(currentrow, alfa33) = ave(5,6,3) !-
       results(currentrow, beta33) = ave(6,6,3) !-
       results(currentrow, gama33) = ave(5,5,3)

       ave(5,6,2) =-ave(5,6,2)
       results(currentrow, alfa32) = ave(5,6,2) !
       results(currentrow, beta32) = ave(6,6,2)  !-
       results(currentrow, gama32) = ave(5,5,2)  !-

       ave(5,6,1) = -ave(5,6,1)
       results(currentrow, alfa31) = ave(5,6,1) !-
       results(currentrow, beta31) = ave(6,6,1)  !-
       results(currentrow, gama31) = ave(5,5,1)  !-

    else

       results(currentrow, alfa13) = zero
       results(currentrow, beta13) = zero
       results(currentrow, gama13) = zero
       results(currentrow, alfa23) = zero
       results(currentrow, beta23) = zero
       results(currentrow, gama23) = zero
       results(currentrow, alfa33) = zero
       results(currentrow, beta33) = zero
       results(currentrow, gama33) = zero
       results(currentrow, alfa32) = zero
       results(currentrow, beta32) = zero
       results(currentrow, gama32) = zero  !-
       results(currentrow, alfa31) = zero
       results(currentrow, beta31) = zero
       results(currentrow, gama31) = zero

    endif

    if( (c_%npara==5)       .or.  (c_%ndpt/=0) ) then
       !when there is no cavity it gives us dispersions
       do i=1,4
          ave(i,1,1)=(Y(i)%t.par.J5)
       enddo
    elseif (c_%nd2 == 6) then
       do i=1,4
          ave(i,1,1) =              (Y(i)%t.par.J5)*(Y(6)%t.par.J6)
          ave(i,1,1) = ave(i,1,1) + (Y(i)%t.par.J6)*(Y(5)%t.par.J5)
       enddo
    else
       do i=1,4
          ave(i,1,1)=zero
       enddo
    endif

    results(currentrow, disp1) = ave(1,1,1)
    results(currentrow, disp2) = ave(2,1,1)
    results(currentrow, disp3) = ave(3,1,1)
    results(currentrow, disp4) = ave(4,1,1)


!!!!!!!!!!!!!!!!
    ! phase advance!
!!!!!!!!!!!!!!!!


    k = 2
    if(c_%nd2==6.and.c_%ndpt==0) k = 3

    do i=1, k
       tx = Y(2*i -1).PAR.fo(2*i,:)
       ty = Y(2*i-1).PAR.fo(2*i-1,:)
       TEST=ATAN2(tx,ty) /TWOPI

       IF(test<zero.AND.abs(test)>EPSIL)  TEST=TEST+one
       DPH=TEST-TESTOLD(i)

       IF(dph<one.AND.abs(dph)>EPSIL) DPH=DPH+one
       IF(dph>half) DPH=DPH-one

       PHASE(i)=PHASE(i)+DPH
       TESTOLD(i)=TEST

    enddo


    results(currentrow, mu1) = phase(1)
    results(currentrow, mu2) = phase(2)
    results(currentrow, mu3) = phase(3)




    results(currentrow, beta11p) = zero
    results(currentrow, beta12p) = zero
    results(currentrow, beta13p) = zero
    results(currentrow, beta21p) = zero
    results(currentrow, beta22p) = zero
    results(currentrow, beta23p) = zero
    results(currentrow, beta31p) = zero
    results(currentrow, beta32p) = zero
    results(currentrow, beta33p) = zero
    results(currentrow, alfa11p) = zero
    results(currentrow, alfa12p) = zero
    results(currentrow, alfa13p) = zero
    results(currentrow, alfa21p) = zero
    results(currentrow, alfa22p) = zero
    results(currentrow, alfa23p) = zero
    results(currentrow, alfa31p) = zero
    results(currentrow, alfa32p) = zero
    results(currentrow, alfa33p) = zero
    results(currentrow, gama11p) = zero
    results(currentrow, gama12p) = zero
    results(currentrow, gama13p) = zero
    results(currentrow, gama21p) = zero
    results(currentrow, gama22p) = zero
    results(currentrow, gama23p) = zero
    results(currentrow, gama31p) = zero
    results(currentrow, gama32p) = zero
    results(currentrow, gama33p) = zero

    results(currentrow, disp1p) = zero
    results(currentrow, disp2p) = zero
    results(currentrow, disp3p) = zero
    results(currentrow, disp4p) = zero

    results(currentrow, disp1p2) = zero
    results(currentrow, disp2p2) = zero
    results(currentrow, disp3p2) = zero
    results(currentrow, disp4p2) = zero

    results(currentrow, disp1p3) = zero
    results(currentrow, disp2p3) = zero
    results(currentrow, disp3p3) = zero
    results(currentrow, disp4p3) = zero

    
!    print*, "parametrictwiss(",currentrow,",bet11):", beta12
!    call print(results(currentrow,beta12),6)
    
    !     print*, "Beta X"
    !     call daprint(ave(1,1,1),6)
    !
    !     print*, "Beta Y"
    !     call daprint(ave(3,3,2),6)
    !
    !     print*, "Beta Z"
    !     call daprint(ave(5,5,3),6)



  end subroutine parametrictwiss
  !_________________________________________________________________________________

  subroutine finishknobs
    implicit none
    integer i,j,k


    if (taylorsallocated == 0) then
       return
    endif

    call kill(phase)
    call kill(testold)
    call kill(test)
    call kill(dph)
    call kill(tx, ty)

    do i=1,c_%nd2
       do j=1,c_%nd2
          do k=1,c_%nd
             call kill(ave(i,j,k))
          enddo
       enddo
    enddo

    taylorsallocated = 0

  end subroutine finishknobs
  !_________________________________________________________________________________

  subroutine resultswithknobs(n,name,y)
    implicit none
    integer              :: n !fibre number
    character(*)         :: name !fibre name
    type(real_8),target  :: y(6)!input 6 dimensional function (polynomial)
    !    real(kind(1d0))      :: coeff
    !    integer              :: i,ii,j,k !iterator
    !    type(taylor)         :: ave(6,6,3)

    print*," resultswithknobs not yet implemented ",n,name
    call print(y(1),6)

  end subroutine resultswithknobs

  !_________________________________________________________________________________
  subroutine writeparresults(filenameIA)
    use twissafi
    implicit none
    integer   filenameIA(*)
    integer       :: mf !macro file descriptor
    character(48) :: filename
    character(48) :: fmat
    integer       :: i,j, nel, slen
    logical       :: fmt_ptc, fmt_tex
    integer, parameter         :: length=16
    character(length)          :: name
    integer                    :: get_string
    integer                    :: restart_sequ,advance_node
    character(48) charconv
    
    if(getdebug()>1) then 
     print*," writeparresults "
     print*,"nv=",c_%nv
     print*,"nd2=",c_%nd2
     print*,"np=",c_%np
     print*,"ndpt=",c_%ndpt 
     print*,"=>",c_%nv-c_%nd2-c_%np
    endif
    
    if (.not. ALLOCATED(results)) then
       call fort_warn("writeparresults","Array with parametric results is not present.")
       print*, "writeparresults tip: it might have been erased the ptc_end command."
       return
    endif

    if (filenameIA(1) > 0) then
       filename = charconv(filenameIA)
       call kanalnummer(mf)
       open(unit=mf,file=filename)
    else
       mf = 6
    endif

    fmt_ptc = .false.
    fmt_tex = .false.

    i = get_string('ptc_printparametric ','format ',fmat)
    if (i > 0) then
       print*,"ptc_printparametric: format is ", fmat(1:i)

       select case(fmat(1:i))
       case ('ptc')
          print*, "ptc_printparametric: Recognized PTC native format"
          fmt_ptc = .true.
       case ('tex')
          print*, "ptc_printparametric: Recognized LaTeX native format"
          fmt_tex = .true.
       case default
          print*, "ptc_printparametric: Format not recognized - using default PTC"
          fmt_ptc = .true.
       end select
    else
       print*, "Got empty Format - using default PTC"
       fmt_ptc = .true.
    endif


    print*,"ptc_printparametric : currentrow is ", currentrow

    nel = restart_sequ()

    do i=1,currentrow-1

       call node_name(name,length)

       write(mf,*) "Magnet ", i," ",name(1:length)


       do j=1,ntwisses
          !        print*, "Writing i,j",i,j
          write(mf,*) twissnames(j)
          if (fmt_ptc)  call print(results(i,j),mf)
          if (fmt_tex)  call printpareq(results(i,j),mf)
          write(mf,*) " "
       enddo


       do j=1,npushes
          if (pushes(j)%index < 1 ) cycle
          !        print*, "Writing i, j->index",i,j,pushes(j)%index
          slen = len_trim(pushes(j)%colname) - 1
          write(mf,*) pushes(j)%colname(1:slen)
          if (fmt_ptc)  call print(results(i,pushes(j)%index),mf)
          if (fmt_tex)  call printpareq(results(i,pushes(j)%index),mf)
          write(mf,*) " "
       enddo
       write(mf,*) "======================================="

       nel = advance_node()

    enddo

    if (mf /= 6) close(mf)

  end subroutine writeparresults
  !_________________________________________________________________________________
  function getparname(n)
    implicit none
    character*(3)  :: getparname
    integer        :: n

    write(getparname,'(A1,i2.2)') "p",n - c_%npara

  end function getparname
  !_________________________________________________________________________________

  subroutine setparvalue(n,v)
    implicit none
    integer :: n
    real :: v

    if (.not. allocated(parvals)) then
       call fort_warn("setparvalue",&
            "array with parameter values is not allocated")
       return
    endif

    if ( (n < 1) .and. (n > ubound(parvals,1)) ) then
       call fort_warn("setparvalue","Array index out of range")
    endif

    parvals(n) = v

  end subroutine setparvalue
  !____________________________________________________________________________________________

  subroutine setknobvalue(fibrenameIA)
    use twissafi
    implicit none
    integer                    :: fibrenameIA(*)
    character(48)              :: fibrename
    integer                    :: i
    integer                    :: kn, ks, par
    real                       :: v
    real(kind(1d0))            :: get_value
    logical(lp)                :: refreshtables
    character(48) charconv
    par = -1
    fibrename = charconv(fibrenameIA)

    !    print *,"setknobvalue: fibrename is ", fibrename


    kn = get_value('ptc_setknobvalue ','kn ') + 1 !madx numerates from 0
    ks = get_value('ptc_setknobvalue ','ks ') + 1

    if ( (kn>0) .and. (ks>0) ) then
       call fort_warn("setknobvalue","Both kn and ks can not be specified together");
       return
    endif

    if ( (kn<=0) .and. (ks<=0) ) then

       select case(fibrename(1:6))
       case ('BETA11')
          par = knobi%beta(1)
       case ('BETA22')
          par = knobi%beta(2)
       case ('BETA33')
          par = knobi%beta(3)
       case ('ALFA11')
          par = knobi%ALFA(1)
       case ('ALFA22')
          par = knobi%ALFA(2)
       case ('ALFA33')
          par = knobi%alfa(3)
       case ('DISP1 ')
          par = knobi%DISPersion(1)
       case ('DISP2 ')
          par = knobi%DISPersion(2)
       case ('DISP3 ')
          par = knobi%DISPersion(3)
       case ('DISP4 ')
          par = knobi%dispersion(4)

       case default
          print*, "Name of initial condition parameter ",fibrename," not recognized"
          return
       end select

       par = par + getnknobsm() !init cond starts as parameters after field components
    else
       do i=1, npolblocks
          if ( polblocks(i)%name == fibrename(1:nlp)) then

             if ( kn>0 ) then
                par = polblocks(i)%ibn(kn)
             elseif ( ks>0 ) then
                par = polblocks(i)%ian(ks)
             endif

             exit
          endif
       enddo

       if (par < 0) then
          call fort_warn("setknobvalue","There is no knob defined on such element");
          return
       endif

    endif


    v = get_value('ptc_setknobvalue ','value ')

    if (getdebug() > 1) then
       print*, "Setting parameter ",par,"(el=",fibrename(1:16),", kn=",kn,", ks=",ks," ) to ", v
    endif

    if (getdebug() > 1) then
        print*, "Setting par ", par, " to ", v, fibrename(1:16)
    end if
    call setparvalue(par, v)

    refreshtables = get_value('ptc_setknobvalue ','refreshtables ') .ne. 0

    if (refreshtables) then
       call filltables()
    endif


  end subroutine setknobvalue
  !____________________________________________________________________________________________

  subroutine filltables
    implicit none
    call filltwisstable()
    !    call cleartables()
    call fillusertables()
  end subroutine filltables
  !____________________________________________________________________________________________

  function getlengthat(n)
    implicit none
    real(dp)         :: getlengthat
    integer          :: n

    print*, "getlengthat, n is ", n

    if (.not. ALLOCATED(spos)) then
       return
    endif

    if ( (n < 1) .and. (n > ubound(spos,1)) ) then
       call fort_warn("getlengthat","position out of range")
       getlengthat = -one
    endif

    print*, "getlengthat, spos at n is ", spos(n)

    getlengthat = spos(n)

  end function getlengthat
  !_________________________________________________________________________________


  subroutine printpareq(ut,iunit)
    implicit none
    type(universal_taylor)  :: ut
    integer                 :: iunit, eqlen
     
    if( .not. associated (ut%n)) then
       call fort_warn("printpareq", "this universal taylor is void")
       write(iunit,'(A)') "this universal taylor is void"
       return
    endif
    
    if(ut%nv /= c_%nv) then
       call fort_warn("printpareq",&
            "number of variables of this universal taylor is different from currnet TPSA")
       call print(ut,6)
       print*,"nv=",c_%nv
       print*,"nd2=",c_%nd2
       print*,"np=",c_%np
       print*,"ndpt=",c_%ndpt 
       print*,"=>",c_%nv-c_%nd2-c_%np
       return
    endif


    call getpareq(ut,textbuffer)
    eqlen = len_trim(textbuffer) + 1;
    write(iunit,'(A)')  textbuffer(1:eqlen)
    

  end subroutine printpareq
  !_________________________________________________________________________________

  subroutine getpareq(ut,string)
    implicit none
    type(universal_taylor) :: ut
    integer                :: i,ii
    integer                :: cpos, last
    character              :: sign
    character*(20)         :: parname
    character(*)           :: string

    !    print*,"getpareq"

    if (.not. associated(ut%n)) then
       call fort_warn("getpareq",&
            "provided taylor is not allocated")
       write(string,'(A)') ' 0 '
       return
    endif


    if(ut%nv /= c_%nv) then
       call fort_warn("getpareq",&
            "number of variables of this universal taylor is different from currnet TPSA")
       return
    endif


    if(ut%n == 0) then
       write(string,'(A)') ' 0 '
       if (getdebug() > 3) then
          print*,"no coefficients in the taylor"
       endif
       return
    endif

    cpos = 1
    last = len(string)
    sign = ' '

    if (getdebug() > 3) then
       print*,"There is ", ut%n, " coefficients "
    endif

    do i = 1,ut%n

       if ( ut%c(i) < zero ) then
          sign = ' '
       endif

       write(string(cpos:last),'(A1,A1,G21.14)') " ",sign,ut%c(i);
       cpos = len_trim(string) + 1;

       if ( ( cpos + 200) > last) then
          call fort_warn("routine madx_ptc_knobs.f90::getpareq",&
               "Buffer for taylor equation is too small! Bailing out!")
          stop;
       endif

       do ii = 1, ut%nv
          !         print*, "cpos=",cpos," last=",last
          if (ut%j(i,ii) /= 0) then
             write(string(cpos:last),'(A1)') "*"; cpos = len_trim(string) + 1;
             parname = getparname(ii)
             write(string(cpos:last),'(a)') parname(1:LEN_TRIM(parname)); cpos = len_trim(string) + 1;

             if (ut%j(i,ii) /= 1) then
                write(string(cpos:last),'(a1)')  "^"; cpos = len_trim(string) + 1;
                write(parname,'(i3)')  ut%j(i,ii)
                parname = ADJUSTL(parname)
                write(string(cpos:last),'(a)')  parname(1:LEN_TRIM(parname)); cpos = len_trim(string) + 1;
             endif
          endif

       enddo

       write(string(cpos:last),'(A1)') " "; cpos = len_trim(string) + 1;
       sign = '+'
    enddo


    !    print*, string
  end subroutine getpareq


end module madx_ptc_knobs_module



! backup code for the results verification
!     call twissfctn(scv,ave)
!     print*, "==============================================================="
!     print*, "==============================================================="
!     print*, "==============================================================="
!     print*, "==============================================================="
!     do ii=1,c_%nd2
!       do jj=ii,c_%nd2
!         do kk=1,c_%nd
!
!          if (ave(ii,jj,kk) /= zero) then
!            write(6,'(a3,i1,i1,i1,a3,f16.8)') 'ave', ii,jj,kk,' = ',ave(ii,jj,kk)
!            do ll=1,3
!              do mm=1,3
!
!                v = abs(ave(ii,jj,kk) - tw%beta(ll,mm))
!                if ( v < 1e-16_dp) then
!                  print*, "v(",ll,",",mm,")=",v
!                  write(6,'(a3,i1,i1,i1,a7,i1,i1)') 'ave', ii,jj,kk,' = beta',ll,mm
!                endif
!
!                v = abs(ave(ii,jj,kk) + tw%alfa(ll,mm))
!                if ( v < 1e-16_dp) then
!                  print*, "v(",ll,",",mm,")=",v
!                  write(6,'(a3,i1,i1,i1,a7,i1,i1)') 'ave', ii,jj,kk,' = alfa',ll,mm
!                endif
!
!
!                v = abs(ave(ii,jj,kk) - tw%gama(ll,mm))
!                if ( v < 1e-16_dp) then
!                  print*, "v(",ll,",",mm,")=",v
!                  write(6,'(a3,i1,i1,i1,a7,i1,i1)') 'ave', ii,jj,kk,' = gama',ll,mm
!                endif
!
!              enddo
!
!              v = abs(ave(ii,jj,kk) - tw%mu(ll))
!              if ( v < 1e-16_dp) then
!                print*, "v(",ll,",",mm,")=",v
!                write(6,'(a3,i1,i1,i1,a7,i1)') 'ave', ii,jj,kk,' = mu',ll
!              endif
!
!              v = abs(ave(ii,jj,kk) - tw%disp(ll))
!              if ( v < 1e-16_dp) then
!                print*, "v(",ll,",",mm,")=",v
!                write(6,'(a3,i1,i1,i1,a7,i1)') 'ave', ii,jj,kk,' = disp',ll
!              endif
!
!              v = abs(ave(ii,jj,kk) - tw%disp(ll+1))
!              if ( v < 1e-16_dp) then
!                print*, "v(",ll,",",mm,")=",v
!                write(6,'(a3,i1,i1,i1,a7,i1)') 'ave', ii,jj,kk,' = disp',ll+1
!              endif
!
!
!            enddo
!
!          endif
!
!       enddo
!    enddo
! enddo
!
!     print*, "==============================================================="
!     print*, "==============================================================="
!     print*, "==============================================================="
!     print*, "==============================================================="


!____________________________________________________________________________________________
!____________________________________________________________________________________________
!____________________________________________________________________________________________
!____________________________________________________________________________________________


function w_ptc_getnknobs()
  use madx_ptc_knobs_module
  implicit none
  integer w_ptc_getnknobs

  w_ptc_getnknobs = getnknobsall()

end function w_ptc_getnknobs
!____________________________________________________________________________________________


function w_ptc_getlengthat(n)
  use madx_ptc_knobs_module
  implicit none
  real(kind(1d0)) :: w_ptc_getlengthat
  integer         :: n

  w_ptc_getlengthat = getlengthat(n)

end function w_ptc_getlengthat
!____________________________________________________________________________________________



subroutine w_ptc_getfctnsnames()
  use madx_ptc_knobs_module
  implicit none

  call getfctnsnames()

end subroutine w_ptc_getfctnsnames
!____________________________________________________________________________________________

subroutine w_ptc_getknobsnames()
  use madx_ptc_knobs_module
  implicit none

  call getknobsnames()

end subroutine w_ptc_getknobsnames
!____________________________________________________________________________________________


subroutine w_ptc_getfunctionat(e,n)
  use madx_ptc_knobs_module
  implicit none
  integer e,n

  call getfunctionat(e,n)

end subroutine w_ptc_getfunctionat

!____________________________________________________________________________________________


function w_ptc_getfunctionvalueat(n,e)
  use madx_ptc_knobs_module
  implicit none
  real(kind(1d0)) :: w_ptc_getfunctionvalueat
  integer e,n


  w_ptc_getfunctionvalueat = getfunctionvalueat(n,e)

end function w_ptc_getfunctionvalueat

!____________________________________________________________________________________________


subroutine w_ptc_rviewer()
  implicit none
  integer :: rviewer
  integer :: res
  res =  rviewer()
end subroutine w_ptc_rviewer

!____________________________________________________________________________________________


subroutine w_ptc_setparvalue(n,v)
  use madx_ptc_knobs_module
  implicit none
  integer :: n
  real    :: v

  call setparvalue(n,v)

end subroutine w_ptc_setparvalue
