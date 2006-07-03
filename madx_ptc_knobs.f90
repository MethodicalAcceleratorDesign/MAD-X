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
  implicit none
  save
  private

  !============================================================================================
  !  PUBLIC INTERFACE
  public                                      :: addknob   !adds knobs (called from c)
  public                                      :: resetknobs ! removes knobs
  public                                      :: setknobs ! sets added knobs on a layout
  public                                      :: getnknobs ! returns number
  public                                      :: resultswithknobs ! called by ptc_twiss to dump the requested results
  public                                      :: parametrictwiss  ! returns lattice functions with dependence on knobs
  
  integer, parameter         ::   maxnpolblocks=10
  type (pol_block), target   ::   polblocks(maxnpolblocks) 
  integer                    ::   npolblocks=0
  integer                    ::   nknobs=0
  
  
contains
  !____________________________________________________________________________________________
  subroutine addknob(fibrenameIA)
    implicit none
    include 'twissa.fi'
    integer                    :: fibrenameIA(*)
    character(48)              :: fibrename
    integer                    :: nint, ndble, k, int_arr(nmax), char_l(nmax), i
    real(kind(1d0))            :: d_arr(nmax)
    character(400)             :: char_a
    type (pol_block), pointer  :: pb
    logical(lp)                :: exactmatch
    real(kind(1d0))            :: get_value
    
    
    npolblocks = npolblocks + 1
    pb => polblocks(npolblocks)
    pb = 0
    
    fibrename = charconv(fibrenameIA)
    
    print *,"addknob: fibrename is ", fibrename
    pb%name = fibrename
!    pb%n_name =  fibrenameIA(1)  !firt element of this array contatins length of the string
    
    exactmatch = get_value('ptc_knob ','exactmatch ') .ne. 0
    if (exactmatch) then
      print*,"addknob: Using Exact name match: ", fibrename
      pb%vorname = fibrename
    else
      print*,"addknob: Using Not Exact name match: all elements starting with ", fibrename
    endif

    call comm_para('kn ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
    print*, "there is ",nint, " kn's: ", int_arr(1:nint)
    do i = 1, nint
     if (int_arr(i) < 1) then
       exit
     endif
     nknobs = nknobs + 1
     pb%ibn(int_arr(i)) = nknobs
     print*, "Set normal mulitpole component ", int_arr(i),"as ", nknobs, "parameter of PTC"
    enddo
    
    call comm_para('ks ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
    print*, "there is ",nint, " ks's: ", int_arr(1:nint)
    do i = 1, nint
     if (int_arr(i) < 1) then
       exit
     endif
     nknobs = nknobs + 1
     pb%ian(int_arr(i)) = nknobs
     print*, "Set skew mulitpole component ", int_arr(i)," as ", nknobs, "parameter of PTC"
    enddo

  end subroutine addknob
  !____________________________________________________________________________________________

  subroutine setknobs(alayout)
    implicit none
    type(layout),pointer       :: alayout
    type (pol_block), pointer  :: pb
    integer                    :: i

    print*, "setknobs: There is ", npolblocks, "pol_blocks"
    do i=1, npolblocks
      pb => polblocks(i)
      print*, "setknobs: Setting pol_block ", i," for ",pb%name
      alayout = pb
    enddo

    print*, "setknobs: All pol_blocks set"
     
  end subroutine setknobs
  !____________________________________________________________________________________________

  subroutine resetknobs()
    implicit none
    integer                    :: i
    
    do i=1, npolblocks
      polblocks(i) = 0
    enddo
    nknobs = 0
    npolblocks = 0
  end subroutine resetknobs
  !____________________________________________________________________________________________

  function getnknobs()
    implicit none
    integer                    :: getnknobs
    
    getnknobs = nknobs

  end function getnknobs
  
  
  !_________________________________________________________________________________
  subroutine parametrictwiss(y,ave)   !  Computes <x_i x_j> assuming linearity and with parameters
    implicit none
    type(real_8) y(6)
    type(taylor) ave(6,6,3)
    integer i,j,k
    integer, allocatable :: e(:)

    allocate(e(c_%npara_fpp))

    e=0
    do i=1,c_%nd2
       do j=i,c_%nd2
          do k=1,c_%nd
             ave(i,j,k)=zero
             e(k*2-1)=1
             ave(i,j,k)=      ave(i,j,k) + (y(i)%t.par.e)*(y(j)%t.par.e) !*
             e(k*2-1)=0
             e(k*2)=1
             ave(i,j,k)=morph(ave(i,j,k))+ (y(i).par.e)*(y(j).par.e) !line * does the same, here taylor is morphed to polimorph,
             e(2*k)=0                                                !and above taylor component of polimorph is used explicitely
             ave(j,i,k)=ave(i,j,k)
          enddo
       enddo
    enddo

    deallocate(e)
    
    print*, "Beta X"
    call daprint(ave(1,1,1),6)

    print*, "Beta Y"
    call daprint(ave(1,1,2),6)

    print*, "Beta Z"
    call daprint(ave(1,1,3),6)
    
  end subroutine parametrictwiss
  !_________________________________________________________________________________

  subroutine resultswithknobs(n,name,y)
    implicit none
    integer              :: n !fibre number
    character(*)         :: name !fibre name
    type(real_8),target  :: y(6)!input 6 dimensional function (polynomial)
    real(kind(1d0))      :: coeff
    integer              :: i,ii,j,k !iterator
    type(taylor)         :: ave(6,6,3)
    
    do i=1,6
       do j=1,6
          do k=1,3
             call alloc(ave(i,j,k))
          enddo
       enddo
    enddo
    
    call print(y(1),6)
    call parametrictwiss(y,ave)

    do i=1,6
       do j=1,6
          do k=1,3
             call kill(ave(i,j,k))
          enddo
       enddo
    enddo

  end subroutine resultswithknobs
  
  !____________________________________________________________________________________________
end module madx_ptc_knobs_module
