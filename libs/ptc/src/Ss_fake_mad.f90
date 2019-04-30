module madx_ptc_module
  !  use madx_keywords
  use S_fitting_new
  implicit none
  public
  TYPE(INTERNAL_STATE),POINTER :: my_state
  TYPE(layout),POINTER :: my_ring,bmadl
  type(mad_universe), pointer :: m_u=>null(),m_t=>null();

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

    call point_m_u(m_u,m_t)
    ind_spin(1,1)=1+6;ind_spin(1,2)=2+6;ind_spin(1,3)=3+6;
    ind_spin(2,1)=4+6;ind_spin(2,2)=5+6;ind_spin(2,3)=6+6;
    ind_spin(3,1)=7+6;ind_spin(3,2)=8+6;ind_spin(3,3)=9+6;    
    k1_spin(1)=1;k2_spin(1)=1;
    k1_spin(2)=1;k2_spin(2)=2;
    k1_spin(3)=1;k2_spin(3)=3;
    k1_spin(4)=2;k2_spin(4)=1;
    k1_spin(5)=2;k2_spin(5)=2;
    k1_spin(6)=2;k2_spin(6)=3;
    k1_spin(7)=3;k2_spin(7)=1;
    k1_spin(8)=3;k2_spin(8)=2;
    k1_spin(9)=3;k2_spin(9)=3;
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
    call point_m_u(m_u,m_t)
    ind_spin(1,1)=1+6;ind_spin(1,2)=2+6;ind_spin(1,3)=3+6;
    ind_spin(2,1)=4+6;ind_spin(2,2)=5+6;ind_spin(2,3)=6+6;
    ind_spin(3,1)=7+6;ind_spin(3,2)=8+6;ind_spin(3,3)=9+6;    
    k1_spin(1)=1;k2_spin(1)=1;
    k1_spin(2)=1;k2_spin(2)=2;
    k1_spin(3)=1;k2_spin(3)=3;
    k1_spin(4)=2;k2_spin(4)=1;
    k1_spin(5)=2;k2_spin(5)=2;
    k1_spin(6)=2;k2_spin(6)=3;
    k1_spin(7)=3;k2_spin(7)=1;
    k1_spin(8)=3;k2_spin(8)=2;
    k1_spin(9)=3;k2_spin(9)=3;
  END subroutine ptc_ini_no_append

  subroutine ptc_end(graphics_maybe,flat_file,file_m_u,file_m_t)
    implicit none
    integer i,i_layout
    character(120) filename
    logical, optional :: flat_file
    integer, optional :: graphics_maybe
    type(layout), pointer :: mring
    character(48) command_gino
    character(*), optional :: file_m_u,file_m_t
    if(.not.associated(m_u)) return
    if(present(graphics_maybe)) then
        if(graphics_maybe>=1) call open_gino_graphics
      do i=2,graphics_maybe
        command_gino="MINI" 
        call call_gino(command_gino)  
      enddo
    endif
    if(present(graphics_maybe)) then
     if(graphics_maybe>0) call close_gino_graphics
    endif
    if(present(flat_file)) then
      if(flat_file) then
        if(associated(m_t%start)) then  ! two universes
           if(present(file_m_u).and.present(file_m_t)) then
            write(6,*) "printing the universes in ", file_m_u(1:len_trim(file_m_u)), &
             file_m_t(1:len_trim(file_m_t))
            call print_universe(m_u,file_m_u)
            call print_universe_pointed(m_u,m_t,file_m_t)

           else
           write(6,*) "printing the universes in ",' m_u.txt and m_t.txt '
            call print_universe(m_u,'m_u.txt')
            call print_universe_pointed(m_u,m_t,'m_t.txt')
           endif
        else
             mring=>m_u%start
       do i=1,m_u%n
        write(filename,*) "flat",i,".txt"
        call context(filename)
        write(6,*) "printing flat file ",filename(1:len_trim(filename))
          call print_new_flat(mring,filename)
          mring=>mring%next
       enddo
       endif
      endif
    endif

    call kill_map_cp()
    call kill_universe(m_t)
    call kill_universe(m_u)
    call kill_tpsa
    call kill(bmadl)
!    do i=1,size(s_b)
       call nul_coef(s_E)
       call nul_coef(S_B_FROM_V)
  !  enddo
  !  deallocate(s_b)



    firsttime_coef=.true.

  end subroutine ptc_end




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
