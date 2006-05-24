module madx_ptc_tablepush_module
  !This module enables the user to store any variable of PTC in a MAD-X table
  !what gives access to them from the other modules, which the most important is MATCH
  use madx_keywords
  use madx_ptc_intstate_module, only : getdebug
  implicit none
  save
  private

  !============================================================================================
  !  PUBLIC INTERFACE
  public                                      :: putusertable
  public                                      :: addpush
  public                                      :: cleartables

  !============================================================================================
  !  PRIVATE
  !    data structures
  type tablepush_poly
     character(20) :: tabname  !table name  to put the coefficients
     character(20) :: colname  !column name to put the coefficients
     integer       :: element  !defines the element of 6D function
     character(10) :: monomial !defiens the monomial f.g. 100000 is coeff of x, and 020300 is coeff of px^2*py^3
  end type tablepush_poly

  type (tablepush_poly), private, dimension(20) :: pushes
  integer                                       :: npushes = 0
  character(20), private, dimension(20)         :: tables  !tables names of existing pushes - each is listed only ones
  integer,       private                        :: ntables = 0 !number of distictive tables
  real(dp),      private, allocatable           :: Dismom(:,:)    ! <xnormal_(2*i-1)**(2j)>= dismon(i,j)*I_i**j

  !    routines
  private                                       :: augment_counts
  private                                       :: issuchtableexist
  private                                       :: putnameintables


contains
  !____________________________________________________________________________________________

  subroutine putusertable(n,name,y)
    !puts the coefficients in tables as defined in array pushes
    implicit none
    integer         :: n !fibre number
    character(*)   :: name !fibre name
    type(real_8),target  :: y(6)!input 6 dimensional function (polynomial)
    type(real_8),pointer :: e !element in array
    real(kind(1d0))      :: coeff
    integer              :: i,ii !iterator
    integer              :: at !iterator

    !    print *,"madx_ptc_tablepush :putusertable "
    !    call daprint(y(1),6)

    call putnameintables()

    if (c_%npara == 6) then
       call sixdmode()
    else
       do i=1,npushes

          e => y(pushes(i)%element)
          coeff = e.sub.(pushes(i)%monomial)
          if (getdebug()>3) then
             write(6,'(a13, a10, a3, f9.6, a10, i1, 5(a13), i3)') &
                  &        "Putting coef ",pushes(i)%monomial,"=",coeff," arr_row ", pushes(i)%element,&
                  &        " in table ", pushes(i)%tabname," at column ", pushes(i)%colname, &
                  &        " for fibre no ",n
          endif

          call double_to_table(pushes(i)%tabname, pushes(i)%colname, coeff);

       enddo
    endif

    call augment_counts()

    !____________________________________________________________________________________________
  contains
    !____________________________________________________________________________________________
    subroutine sixdmode()
      implicit none
      integer ele
      character bufchar
      character(10) monstr

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

         coeff = e.sub.(monstr)

         if (getdebug()>3) then
            write(6,'(a13, a10, a3, f9.6, a10, i1, 5(a13), i3)') &
                 &        "Put 6D coef ",pushes(i)%monomial,"=",coeff," arr_row ", pushes(i)%element,&
                 &        " in table ", pushes(i)%tabname," at column ", pushes(i)%colname, &
                 &        " for fibre no ",n
         endif

         call double_to_table(pushes(i)%tabname, pushes(i)%colname, coeff);

      enddo
    end subroutine sixdmode

    !____________________________________________________________________________________________

  end subroutine putusertable
  !____________________________________________________________________________________________

  subroutine inittables()
    implicit none
    integer  :: i ! iterator

    deallocate(dismom)

    allocate(dismom(c_%nd,0:c_%no/2))


  end subroutine inittables

  subroutine cleartables()
    implicit none
    integer  :: i ! iterator
    ! we do not deal here with the main twiss table, hence we do not clear it
    do i=1,ntables
       if (getdebug()>3) print *,"Clearing ",tables(i)
       call reset_count(tables(i))
    enddo

  end subroutine cleartables
  !____________________________________________________________________________________________

  subroutine augment_counts()
    implicit none
    integer  :: i ! iterator

    do i=1,ntables
       if (getdebug()>3) print *,"Augmenting ",tables(i)
       call augmentcountonly(tables(i)) !we need to use special augement,
       !cause the regular one looks for for
       !variables names like columns to fill the table
    enddo

  end subroutine augment_counts
  !____________________________________________________________________________________________

  subroutine putnameintables()
    implicit none
    integer       :: i ! iterator
    do i=1,ntables
       if (getdebug()>2) print *,"Putting name in ",tables(i)
       call string_to_table(tables(i),"name ","name ")
    enddo

  end subroutine putnameintables
  !____________________________________________________________________________________________


  subroutine addpush(table,column,element,monomial)
    implicit none
    include 'twissa.fi'
    integer   table(*)
    integer   column(*)
    integer   element
    integer   monomial(*)

    npushes = npushes + 1
    pushes(npushes)%tabname = charconv(table)
    pushes(npushes)%colname = charconv(column)
    pushes(npushes)%element = element
    pushes(npushes)%monomial = charconv(monomial)
    !imput "int string"  has the length of the string at the first plave
    pushes(npushes)%tabname(table(1)+1:table(1)+1)=achar(0)
    pushes(npushes)%colname(column(1)+1:column(1)+1)=achar(0)
    pushes(npushes)%monomial(monomial(1)+1:monomial(1)+1)=achar(0)

    if (getdebug()>3) then
       print  *,"madx_ptc_tablepush : addpush(",&
            &          pushes(npushes)%tabname,">,<",pushes(npushes)%colname,">,<",&
            &          pushes(npushes)%element,">,<",pushes(npushes)%monomial,">)"
    endif

    if ( issuchtableexist(pushes(npushes)%tabname) .eqv. .false.) then
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

  subroutine average_x_i_x_j(y,ave,i,j)   !  Computes <x_i x_j>
    implicit none
    type(real_8) y(6)
    type(taylor) ave
    type(taylorresonance) tr
    integer i,j

    call alloc(tr)

    if(j/=0) then
       tr=y(i)%t*y(j)%t
    else
       tr=y(i)%t
    endif

    call cfu(tr%cos,filter,ave)

    call kill(tr)
  end subroutine average_x_i_x_j
  !_________________________________________________________________________________

  real(dp) function filter(e)   !  Computes <x_i x_j>
    implicit none
    integer e(:)
    integer i

    filter=one

    do i=1,c_%nd
       if(e(2*i-1)/=e(2*i)) then
          filter=zero
          return
       else
          filter=filter*dismom(i,e(2*i))
       endif
    enddo

  end function filter
  !_________________________________________________________________________________



  subroutine make_gaussian(plane,I0)
    implicit none
    integer plane,i
    real(dp) I0

    dismom(plane,0)=one

    do i=1,c_%no/2
       dismom(plane,i)=i*I0*two*dismom(plane,i-1)
    enddo

  end   subroutine make_gaussian
  !_________________________________________________________________________________

  subroutine make_ring(plane,I0)
    implicit none
    integer plane,i
    real(dp) I0

    dismom(plane,0)=one

    do i=1,c_%no/2
       dismom(plane,i)=I0*two*dismom(plane,i-1)
    enddo

  end   subroutine make_ring
  !_________________________________________________________________________________


end module madx_ptc_tablepush_module
