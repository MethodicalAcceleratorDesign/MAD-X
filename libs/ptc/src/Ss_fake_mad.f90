module madx_ptc_module
  !  use madx_keywords
  use S_fitting_new
  implicit none
  public
  TYPE(INTERNAL_STATE),POINTER :: my_state
  TYPE(layout),POINTER :: my_ring,bmadl
  type(mad_universe), pointer :: m_u,m_t

contains

  subroutine ptc_INI()
    implicit none

    allocate(m_u)
    allocate(m_t)
    allocate(bmadl)
    call set_up_universe(m_u)
    call append_empty_layout(m_u)
    call set_up_universe(m_t)
    call append_empty_layout(m_t)
    call set_up(bmadl)
    bmadl%NAME='BMAD REUSED FIBRE LAYOUT'
  END subroutine ptc_ini

  subroutine ptc_ini_no_append()
    implicit none

    allocate(m_u)
    call set_up_universe(m_u)
    allocate(m_t)
    call set_up_universe(m_t)
    allocate(bmadl)
    call set_up(bmadl)
    bmadl%NAME='BMAD REUSED FIBRE LAYOUT'

  END subroutine ptc_ini_no_append

  subroutine ptc_end()
    implicit none
    integer i

    call kill_universe(m_t)
    call kill_universe(m_u)
    call kill_tpsa
    call kill(bmadl)
    do i=1,size(s_b)
       call nul_coef(s_b(i))
    enddo
    deallocate(s_b)



    firsttime_coef=.true.

  end subroutine ptc_end


  subroutine create_fibre_reuse(key,EXCEPTION,magnet_only)  
    implicit none

    logical(lp), optional :: magnet_only
    type(keywords) key
    INTEGER EXCEPTION  !,NSTD0,METD0
    logical(lp) doneit,append
    type(fibre), pointer :: current


     if(associated(bmadl%end)) then
      IF(ASSOCIATED(bmadl%T)) THEN
         CALL kill_NODE_LAYOUT(bmadl%T)  !  KILLING THIN LAYOUT
         nullify(bmadl%T)
        if(lielib_print(12)==1) WRITE(6,*) " NODE LAYOUT HAS BEEN KILLED "
       ENDIF      
        bmadl%end=-1
       else
        call append_empty(bmadl)
     endif

     call  create_fibre(bmadl%end,key,EXCEPTION,magnet_only)

    bmadl%closed=.true.

    doneit=.true.
    call ring_l(bmadl,doneit)

    call survey(bmadl)
    call MAKE_NODE_LAYOUT( bmadl)
  end subroutine create_fibre_reuse

end module madx_ptc_module

character * 48 function charconv(tint)
  !----------------------------------------------------------------------*
  ! purpose:                                                             *
  !   converts integer array to string (based on ascii)                  *
  ! input:                                                               *
  !   tint  (int array)  1 = length, rest = string                       *
  !----------------------------------------------------------------------*
  implicit none
  integer tint(*)
  integer i, j, m, n
  parameter (m = 128)
  character *(m) letter
  data letter /                                                     &
       &'                                !"#$%&''()*+,-./0123456789:;<=>?@&
       &ABCDEFGHIJKLMNOPQRSTUVWXYZ[ ]^_`abcdefghijklmnopqrstuvwxyz{|}~'/
  charconv = ' '
  n = tint(1)
  do i = 1, n
     j = tint(i+1)
     if (j .lt. m)  charconv(i:i) = letter(j:j)
  enddo
end function charconv
