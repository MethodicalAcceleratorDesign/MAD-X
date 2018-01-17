!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

MODULE S_FIBRE_BUNDLE
  USE S_DEF_ELEMENT
  !  USE   S_ELEMENTS
  ! Implementation of abstract data type as a linked layout
  IMPLICIT NONE
  public
  private unify_mad_universe

  PRIVATE kill_layout,kill_info,alloc_info,copy_info
  private dealloc_fibre,append_fibre   !, alloc_fibre public now also as alloc
  !  private null_it0
  private move_to_p,move_to_name_old,move_to_nameS,move_to_name_FIRSTNAME 
  PRIVATE append_EMPTY_FIBRE
  PRIVATE FIND_PATCH_0,FIND_PATCH_bmad_0
  PRIVATE FIND_PATCH_p_new
  PRIVATE INDEX_0,index_1
  private FIND_POS_in_universe,FIND_POS_in_layout,super_dealloc_fibre
  TYPE(LAYOUT), PRIVATE, POINTER:: LC
  logical :: superkill=.false.
  logical(lp),TARGET :: use_info=.true., use_info_m=.false.
  integer, target :: nsize_info = 70
  private zero_fibre
  INTEGER :: INDEX_0=0, INDEX_1=0
  INTEGER :: INDEX_node=0
  logical(lp),PRIVATE,PARAMETER::T=.TRUE.,F=.FALSE.
  real(dp),target :: eps_pos=1e-10_dp
  integer(2),parameter::it0=0,it1=1,it2=2,it3=3,it4=4,it5=5,it6=6,it7=7,it8=8,it9=9
  INTEGER,private,parameter :: IPOS =1000

  INTERFACE kill
     MODULE PROCEDURE kill_layout
     MODULE PROCEDURE dealloc_fibre
     MODULE PROCEDURE kill_info
     MODULE PROCEDURE kill_NODE_LAYOUT
     MODULE PROCEDURE de_Set_Up_ORBIT_LATTICE
     MODULE PROCEDURE kill_BEAM_BEAM_NODE
  END INTERFACE

  INTERFACE super_kill
     MODULE PROCEDURE super_dealloc_fibre
  end INTERFACE

  INTERFACE alloc
     !     MODULE PROCEDURE set_up
     MODULE PROCEDURE alloc_fibre
     MODULE PROCEDURE alloc_info
     MODULE PROCEDURE ALLOC_BEAM_BEAM_NODE
  END INTERFACE

  INTERFACE copy
     MODULE PROCEDURE copy_info
  END INTERFACE

  INTERFACE append
     MODULE PROCEDURE append_fibre
  END INTERFACE

  INTERFACE append_EMPTY
     MODULE PROCEDURE append_EMPTY_FIBRE
  END INTERFACE

  INTERFACE move_to
     MODULE PROCEDURE move_to_p
     MODULE PROCEDURE move_to_name_old
     MODULE PROCEDURE move_to_nameS
     MODULE PROCEDURE move_to_name_FIRSTNAME
   END INTERFACE

  INTERFACE FIND_PATCH
     MODULE PROCEDURE FIND_PATCH_0
  END INTERFACE

  INTERFACE FIND_PATCH_bmad
     MODULE PROCEDURE FIND_PATCH_bmad_0
  END INTERFACE



  INTERFACE FIND_pos
     MODULE PROCEDURE FIND_POS_in_layout
     MODULE PROCEDURE FIND_POS_in_universe
  END INTERFACE




  interface assignment (=)
     ! MODULE PROCEDURE null_it0
     MODULE PROCEDURE zero_fibre
  end interface

CONTAINS

  SUBROUTINE alloc_info( c ) ! Does the full allocation of fibre and initialization of internal variables
    implicit none
    type(info),target, intent(inout):: c

    allocate(c%s) ;c%s=0.0_dp;
    allocate(c%beta(nsize_info));c%beta=0.0_dp;
    allocate(c%fix(6));c%fix=0.0_dp;
    allocate(c%fix0(6));c%fix0=0.0_dp;
    allocate(c%pos(2));c%pos=0.0_dp;
    if(use_info_m) then 
      allocate(c%m(6,6))
      c%m=0
    endif

  end SUBROUTINE alloc_info

  SUBROUTINE copy_info( c,d ) ! Does the full allocation of fibre and initialization of internal variables
    implicit none
    type(info),target, intent(in)::c
    type(info),target,  intent(inout)::d

    !   d%name=c%name
    d%s=c%s
    d%beta=c%beta
    d%fix=c%fix
    d%fix0=c%fix0
    d%pos=c%pos
    if(associated(c%m))d%m=c%m
  end SUBROUTINE copy_info

  SUBROUTINE kill_info( c ) ! Does the full allocation of fibre and initialization of internal variables
    implicit none
    type(info),target, intent(inout):: c

    !   deallocate(c%name)
    deallocate(c%s)
    deallocate(c%fix)
    deallocate(c%fix0)
    deallocate(c%beta)
    deallocate(c%pos)
    if(associated(c%m)) deallocate(c%m)

  end SUBROUTINE kill_info

  SUBROUTINE APPEND_mad_like( L, el )  ! Used in MAD-Like input
    implicit none
    TYPE (fibre),target :: el
    TYPE (fibre), POINTER :: Current
    TYPE (layout), TARGET, intent(inout):: L
    L%N=L%N+1
    CALL ALLOCATE_FIBRE(Current);
    current%mag=>el%mag
    current%magp=>el%magp
    current%CHART=>el%CHART
    current%PATCH=>el%PATCH
    if(use_info) current%i=>el%i
    current%dir=>el%dir
    !  OCTOBER 2007
    !        current%P0C=>el%P0C
    current%BETA0=>el%BETA0
    current%GAMMA0I=>el%GAMMA0I
    current%GAMBET=>el%GAMBET
    current%MASS=>el%MASS
    current%AG=>el%AG
    current%CHARGE=>el%CHARGE

    current%PARENT_LAYOUT=>L
    if(L%N==1) current%next=> L%start
    Current % previous => L % end  ! point it to next fibre
    if(L%N>1)  THEN
       L % end % next => current      !
    ENDIF

    L % end => Current
    if(L%N==1) L%start=> Current

    L%LASTPOS=L%N ;
    L%LAST=>CURRENT;

  END SUBROUTINE APPEND_mad_like


  SUBROUTINE kill_layout( L )  ! Destroys a layout
    implicit none
    TYPE (fibre), POINTER :: Current
    TYPE (layout), TARGET, intent(inout):: L
    logical(lp) doneit
    if(lielib_print(12)==1) write(6,*) "Killing Layout",L%name
    CALL LINE_L(L,doneit)
    nullify(current)
    IF(ASSOCIATED(L%T)) THEN
       CALL kill_NODE_LAYOUT(L%T)  !  KILLING THIN LAYOUT
       nullify(L%T)
       if(lielib_print(12)==1) WRITE(6,*) " NODE LAYOUT HAS BEEN KILLED "
    ENDIF
    IF(ASSOCIATED(L%DNA)) THEN
       DEALLOCATE(L%DNA)
       if(lielib_print(12)==1) WRITE(6,*) " DNA CONTENT HAS BEEN DEALLOCATED "
    ENDIF
    !    IF(ASSOCIATED(L%con)) THEN
    !       DEALLOCATE(L%con)
    !       if(lielib_print(12)==1) WRITE(6,*) " CONNECTOR CONTENT HAS BEEN KILLED "
    !    ENDIF
    !    IF(ASSOCIATED(L%con1)) THEN
    !       DEALLOCATE(L%con1)
    !       if(lielib_print(12)==1) WRITE(6,*) " CONNECTOR CONTENT HAS BEEN DEALLOCATED "
    !    ENDIF
    !    IF(ASSOCIATED(L%con2)) THEN
    !       DEALLOCATE(L%con2)
    !       if(lielib_print(12)==1) WRITE(6,*) " CONNECTOR CONTENT HAS BEEN DEALLOCATED "
    !    ENDIF
    !    IF(ASSOCIATED(L%girder)) THEN
    !       DEALLOCATE(L%girder)
    !       if(lielib_print(12)==1) WRITE(6,*) " GIRDER CONTENT HAS BEEN DEALLOCATED "
    !    ENDIF

    LC=> L  ! USED TO AVOID DNA MEMBERS
    Current => L % end      ! end at the end
    DO WHILE (ASSOCIATED(L % end))
       L % end => Current % previous  ! update the end before disposing
       call dealloc_fibre(Current)
       Current => L % end     ! alias of last fibre again
       L%N=L%N-1
    END DO
    call de_set_up(L)
    if(lielib_print(12)==1) WRITE(6,*) 'Layout killed '
  END SUBROUTINE kill_layout


  SUBROUTINE APPEND_fibre( L, el ) ! Standard append that clones everything
    implicit none
    TYPE (fibre),target, intent(in) :: el
    TYPE (fibre), POINTER :: Current
    TYPE (layout), TARGET,intent(inout):: L
    logical(lp) doneit
    CALL LINE_L(L,doneit)
    L%N=L%N+1
    nullify(current)
    call alloc_fibre(current)
    !   call copy(el%magp,current%mag)       ! 2010 etienne does not understand!
    !   call copy(current%mag,current%magp) ! 2010 etienne does not understand!
    !   call copy(el%mag,current%mag)   ! 2010 etienne does not understand!
    call copy(el%mag,current%mag)       ! 2010 etienne replaces!
    call copy(current%mag,current%magp) ! 2010 etienne replaces!
    !   write(6,*) " used "
    !if(associated(current%CHART))
    call copy(el%CHART,current%CHART)
    !if(associated(current%patch))
    call copy(el%PATCH,current%PATCH)
    if(use_info.and.associated(current%patch)) call copy(el%i,current%i)
    current%dir=el%dir
    !        current%P0C    =el%P0C
    current%BETA0  =el%BETA0
    current%GAMMA0I=el%GAMMA0I
    current%GAMBET =el%GAMBET
    current%MASS  =el%MASS
    current%AG  =el%AG
    current%CHARGE =el%CHARGE

    current%PARENT_LAYOUT=>L
    current%mag%PARENT_FIBRE=>current
    current%magP%PARENT_FIBRE=>current
    !    current%magp%PARENT_FIBRE=>current
    if(L%N==1) then
     current%next=> L%start
      index_0=index_0+1
      L%index=index_0    
    endif

    Current % previous => L % end  ! point it to next fibre
    if(L%N>1)  THEN
       L % end % next => current      !
    ENDIF

    L % end => Current
    if(L%N==1) L%start=> Current
    current%pos=l%n

    L%LASTPOS=L%N ; L%LAST=>CURRENT;
    CALL RING_L(L,doneit)
  END SUBROUTINE APPEND_fibre

  SUBROUTINE APPEND_clone( L, muonfactor,charge ) ! Standard append that clones everything
    implicit none
    TYPE (fibre), POINTER :: Current
    TYPE (layout), TARGET,intent(inout):: L
    logical(lp) doneit
    real(dp),optional :: charge
    real(dp),optional :: muonfactor
    real(dp) mu
    real(dp) ch
    CALL LINE_L(L,doneit)
    L%N=L%N+1
    nullify(current)
    call alloc_fibre(current)
    !    if(use_info.and.associated(current%patch)) call copy(el%i,current%i)
    current%dir=1
    mu=1.0_dp
    ch=1
    if(present(muonfactor)) mu=muonfactor
    if(present(charge)) ch=charge
    ! OCT 2007
    !        current%P0C=ONE
    current%BETA0=1.0_dp
    current%GAMMA0I=1.0_dp
    current%GAMBET=0.0_dp
    current%MASS=mu*pmae
    current%AG=A_particle
    current%CHARGE=ch

    current%pos=l%n

    current%PARENT_LAYOUT=>L
    current%mag%PARENT_FIBRE=>current
    current%magP%PARENT_FIBRE=>current
    !    current%magp%PARENT_FIBRE=>current
    if(L%N==1) current%next=> L%start
    Current % previous => L % end  ! point it to next fibre
    if(L%N>1)  THEN
       L % end % next => current      !
    ENDIF

    L % end => Current
    if(L%N==1) L%start=> Current

    L%LASTPOS=L%N ; L%LAST=>CURRENT;
    CALL RING_L(L,doneit)
  END SUBROUTINE APPEND_clone








  SUBROUTINE move_to_p( L,current,POS ) ! Moves current to the i^th position
    implicit none
    TYPE (fibre), POINTER :: Current
    TYPE (layout), TARGET, intent(inout):: L
    integer i,k,POS

    !    CALL LINE_L(L,doneit)  !TGV
    I=mod_n(POS,L%N)
 
    IF(L%LASTPOS==0) THEN
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,a72,/),(1X,a72))'
       write(6,*) " L%LASTPOS=0 : ABNORMAL UNLESS LINE EMPTY"
       write(*,'(a7,i4)')" L%N = ",L%N
       ! call !write_e(-124)
    ENDIF

    nullify(current);
    Current => L%LAST

    k=L%LASTPOS
    IF(I>=L%LASTPOS) THEN
       DO K=L%LASTPOS,I-1
          !      DO WHILE (ASSOCIATED(Current).and.k<i) !TGV
          !          k=k+1 !TGV
          Current => Current % next
       END DO
    ELSE
       DO K=L%LASTPOS,I+1,-1
          !       DO WHILE (ASSOCIATED(Current).and.k>i) !TGV
          !          k=k-1 !TGV
          Current => Current % PREVIOUS
       END DO
    ENDIF
 

    L%LASTPOS=I; L%LAST => Current;
    !    CALL RING_L(L,doneit) ! TGV
  END SUBROUTINE move_to_p

 SUBROUTINE move_to_p_safe( L,current,POS ) ! Moves current to the i^th position
    implicit none
    TYPE (fibre), POINTER :: Current
    TYPE (layout), TARGET, intent(inout):: L
    integer i,k,POS

    !    CALL LINE_L(L,doneit)  !TGV
    I=mod_n(POS,L%N)
 

    nullify(current);
    Current => L%start

 
 
       DO K=1,I-1
          Current => Current % next
       END DO
  
  END SUBROUTINE move_to_p_safe


  SUBROUTINE move_to_name_old( L,current,name,pos,reset) ! moves to next one in list called name
    implicit none
    logical(lp),optional :: reset    
    TYPE (fibre), POINTER :: Current
    TYPE (layout), TARGET, intent(inout):: L
    integer, optional, intent(inout):: pos
    character(*), intent(in):: name
    CHARACTER(nlp) S1NAME
    integer i,poss

    logical(lp) foundit
    TYPE (fibre), POINTER :: p
    
    if(present(reset)) then
     if(reset) then
       l%lastpos=1
       l%last=>L%start
     endif
    endif
    
    foundit=.false.
    S1NAME=name
    CALL CONTEXT(S1name)
    
    nullify(p)
    p=>l%last%next

    if(.not.associated(p)) goto 100
    do i=1,l%n
       if(p%mag%name==s1name) then
          foundit=.true.
          goto 100
       endif
       p=>p%next
       if(.not.associated(p)) goto 100
    enddo
100 continue
    if(foundit) then
       current=>p
       poss=mod_n(l%lastpos+i,l%n)
       l%lastpos=poss
       l%last=>current
    else
       poss=0
       WRITE(6,*) " Fibre not found in move_to_name_old ",S1name
    endif
    if(present(pos)) pos=poss
  END SUBROUTINE move_to_name_old

  SUBROUTINE move_to_partial( L,current,name,pos) ! moves to next one in list called name
    implicit none
    TYPE (fibre), POINTER :: Current
    TYPE (layout), TARGET, intent(inout):: L
    integer,optional, intent(inout):: pos
    character(*), intent(in):: name
    CHARACTER(nlp) S1NAME
    integer i,poss

    logical(lp) foundit
    TYPE (fibre), POINTER :: p

    foundit=.false.
    S1NAME=name
    CALL CONTEXT(S1name)
 
    nullify(p)
    p=>l%last%next

    if(.not.associated(p)) goto 100
    do i=1,l%n
       if(index(p%mag%name,s1name(1:len_trim(s1name)))/=0) then
          foundit=.true.
 
          goto 100
       endif
       p=>p%next
       if(.not.associated(p)) goto 100
    enddo
100 continue
    if(foundit) then
       current=>p
       poss=mod_n(l%lastpos+i,l%n)
       l%lastpos=poss
       l%last=>current
    else
       poss=0
    endif
 
 if(present(pos)) pos=poss

  END SUBROUTINE move_to_partial

  SUBROUTINE move_to_name_FIRSTNAME( L,current,name,VORNAME,pos) ! moves to next one in list called name
    implicit none
    TYPE (fibre), POINTER :: Current
    TYPE (layout), TARGET, intent(inout):: L
    integer, intent(inout):: pos
    character(*), intent(in):: name,VORNAME
    CHARACTER(nlp) S1NAME,S2NAME
    integer i

    logical(lp) foundit
    TYPE (fibre), POINTER :: p

    foundit=.false.
    S1NAME=name
    S2NAME=VORNAME
    CALL CONTEXT(S1name)
    CALL CONTEXT(S2name)

    nullify(p)
    p=>l%last%next

    if(.not.associated(p)) goto 100
    do i=1,l%n
       if(p%mag%name==s1name.AND.p%mag%VORname==S2NAME) then
          foundit=.true.
          goto 100
       endif
       p=>p%next
       if(.not.associated(p)) goto 100
    enddo
100 continue
    if(foundit) then
       current=>p
       pos=mod_n(l%lastpos+i,l%n)
       l%lastpos=pos
       l%last=>current
    else
       pos=0
       WRITE(6,*) " Did not find in move_to_name_FIRSTNAME"
       WRITE(6,*) s1name,s2name
    endif
  END SUBROUTINE move_to_name_FIRSTNAME

  SUBROUTINE move_to_nameS( L,current,name,posR,POS) ! moves to next one in list called name
    implicit none
    TYPE (fibre), POINTER :: Current
    TYPE (layout), TARGET, intent(inout):: L
    integer, intent(inout):: pos,POSR
    character(*), intent(in):: name
    CHARACTER(nlp) S1NAME
    integer i,IC

    logical(lp) foundit
    TYPE (fibre), POINTER :: p

    foundit=.false.
    S1NAME=name
    CALL CONTEXT(S1name)

    nullify(p)
    p=>l%START
    IC=0
    if(.not.associated(p)) goto 100
    do i=1,l%n
       if(p%mag%name==s1name) then
          IC=IC+1
          IF(IC==POSR) THEN
             foundit=.true.
             goto 100
          ENDIF
       endif
       p=>p%next
       if(.not.associated(p)) goto 100
    enddo
100 continue
    if(foundit) then
       current=>p
       pos=mod_n(i,l%n)
       l%lastpos=pos
       l%last=>current
    else
       pos=0
    endif
  END SUBROUTINE move_to_nameS

 SUBROUTINE move_to_i( L,current,POS) !      move_to_i   ! move to ith fibre
    implicit none
    TYPE (fibre), POINTER :: Current
    TYPE (layout), TARGET, intent(inout):: L
    integer, intent(inout):: pos
    integer i

    logical(lp) foundit
    TYPE (fibre), POINTER :: p

    foundit=.false.

    nullify(p)
    p=>l%START
    if(.not.associated(p)) goto 100

    do i=1,l%n
       if(p%pos==pos) then
         foundit=my_true
         exit
       endif
       if(l%lastpos>=pos) then
        p=>p%previous
       else
        p=>p%next
       endif
       if(.not.associated(p)) goto 100
    enddo
100 continue
    if(foundit) then
       current=>p
       l%lastpos=pos
       l%last=>current
    else
       pos=0
    endif
  END SUBROUTINE move_to_i




  SUBROUTINE Set_Up( L ) ! Sets up a layout: gives a unique negative index
    implicit none
    TYPE (layout),TARGET, INTENT(INOUT):: L
    type(mad_universe), pointer :: madu
!  new 2012.9.7
    nullify(madu)
   if(associated(L%parent_universe) ) madu=>L%parent_universe
!    

    CALL NULLIFY_LAYOUT(L)
    
    if(associated(madu) ) L%parent_universe=>madu
     nullify(madu)   

    ALLOCATE(L%closed);  ALLOCATE(L%lastpos);ALLOCATE(L%NAME);ALLOCATE(L%HARMONIC_NUMBER);
    ALLOCATE(L%NTHIN);ALLOCATE(L%THIN);ALLOCATE(L%INDEX);
    ALLOCATE(L%n);
    L%closed=.false.;
    L%NTHIN=0;L%THIN=0.0_dp;
    L%N=0;
    L%lastpos=0;L%NAME='No name assigned';
!    INDEX_0=INDEX_0+1
!    L%INDEX=INDEX_0
    L%HARMONIC_NUMBER=0
  END SUBROUTINE Set_Up




  SUBROUTINE de_Set_Up( L ) ! deallocates layout content
    implicit none
    TYPE (layout),TARGET, INTENT(INOUT):: L
    deallocate(L%closed);deallocate(L%lastpos);deallocate(L%NAME);deallocate(L%HARMONIC_NUMBER);
    deallocate(L%INDEX);
    deallocate(L%NTHIN);deallocate(L%THIN);
    deallocate(L%n);          !deallocate(L%parent_universe)   left out
    IF(ASSOCIATED(L%T)) deallocate(L%T);
  END SUBROUTINE de_Set_Up


  SUBROUTINE nullIFY_LAYOUT( L ) ! Nullifies layout content,i
    implicit none
    !   integer , intent(in) :: i
    TYPE (layout),TARGET, intent(inout) :: L
    !   if(i==0) then
    nullify(L%T)  ! THIN LAYOUT
    nullify(L%DNA)  ! THIN LAYOUT
    !    nullify(L%CON)  ! THIN LAYOUT
    !    nullify(L%CON1)  ! THIN LAYOUT
    !    nullify(L%CON2)  ! THIN LAYOUT
    !    nullify(L%girder)  ! THIN LAYOUT
    nullify(L%parent_universe)
    nullify(L%INDEX)
    nullify(L%HARMONIC_NUMBER)
    nullify(L%NAME)
    nullify(L%CLOSED,L%N )
    nullify(L%NTHIN      )
    nullify(L%THIN   )
    nullify(L%LASTPOS )  ! POSITION OF LAST VISITED
    nullify(L%LAST )! LAST VISITED
    !
    nullify(L%END )
    nullify(L%START )
    nullify(L%START_GROUND )! STORE THE GROUNDED VALUE OF START DURING CIRCULAR SCANNING
    nullify(L%END_GROUND )! STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING
    !   nullify(L%NEXT )! STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING
    !   nullify(L%PREVIOUS )! STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING
    !  nullify(L%parent_universe ) ! left out
    !  else
    !    !w_p=0
    !    !w_p%nc=1
    !    !w_p%fc='(1((1X,a72)))'
    !    write(6,*) " Only =0 permitted (nullify) "
    !    ! call !write_e(100)
    ! endif
  END SUBROUTINE nullIFY_LAYOUT



  SUBROUTINE LINE_L(L,doneit) ! makes into line temporarily
    implicit none
    TYPE (layout), TARGET, intent(inout):: L
    logical(lp) doneit
    doneit=.false.
    if(L%closed)  then
       if(associated(L%end%next)) then
          L%end%next=>L%start_ground
          doneit=.true.
       endif
       if(associated(L%start%previous)) then
          L%start%previous=>L%end_ground
       endif
    endif
  END SUBROUTINE LINE_L

  SUBROUTINE RING_L(L,doit) ! Brings back to ring if needed
    implicit none
    TYPE (layout), TARGET, intent(inout):: L
    logical(lp) doit
    if(L%closed.and.doit)  then
       if(.NOT.(associated(L%end%next))) then
          L%start_ground=>L%end%next      ! saving grounded pointer
          L%end%next=>L%start
       endif
       if(.NOT.(associated(L%start%previous))) then
          L%end_ground=>L%start%previous  ! saving grounded pointer
          L%start%previous=>L%end
       endif
    endif
  END SUBROUTINE RING_L


  SUBROUTINE APPEND_POINT( L, el )   ! Appoints without cloning
    implicit none
    TYPE (fibre),POINTER :: el
    TYPE (fibre), POINTER :: Current
    TYPE (layout), TARGET, intent(inout):: L
    !    type(fibre), pointer :: p
    logical(lp) doneit
    TYPE(fibre_appearance), POINTER :: D
    !    nullify(p);
    CALL LINE_L(L,doneit)
    L%N=L%N+1
    CALL ALLOCATE_FIBRE(Current);
    ALLOCATE(Current%PATCH);


    CURRENT%PARENT_LAYOUT=>L     !
    current%mag=>el%mag
    current%magp=>el%magp
    current%CHART=>el%CHART
    current%PATCH=0   ! new patches always belong to fibre  ! this was the error Weishi
    IF(EL%PATCH%PATCH/=0) THEN
       IF(.NOT.ASSOCIATED(CURRENT%PATCH)) CURRENT%PATCH=0
       CALL COPY(EL%PATCH,current%PATCH)
    ENDIF
    !    if(use_info) current%i=>el%i
    if(use_info) then
       allocate(current%i)
       call alloc(current%i)
    endif

    ALLOCATE(current%DIR)   !;
    !    ALLOCATE(current%P0C);
    ALLOCATE(current%BETA0);
    ALLOCATE(current%GAMMA0I);
    ALLOCATE(current%GAMBET);
    ALLOCATE(current%MASS);
    ALLOCATE(current%AG);
    ALLOCATE(current%CHARGE);
    current%dir=el%dir
    !    current%P0C=el%P0C
    current%BETA0=el%BETA0
    current%GAMMA0I=el%GAMMA0I
    current%GAMBET=el%GAMBET
    current%MASS=el%MASS
    current%AG=el%AG
    current%CHARGE=el%CHARGE


    ALLOCATE(Current%pos);
    current%pos=l%n
    !    current%P0C=el%P0C
    !    current%BETA0=el%BETA0
    if(L%N==1) then 
      current%next=> L%start
      index_1=index_1+1
      L%index=index_1
    endif
    Current % previous => L % end  ! point it to next fibre
    if(L%N>1)  THEN
       L % end % next => current      !
    ENDIF

    L % end => Current
    if(L%N==1) L%start=> Current
    if(.not.associated(current%pos)) allocate(current%pos)
    current%pos=l%n

    L%LASTPOS=L%N ;
    L%LAST=>CURRENT;
    CALL RING_L(L,doneit)

    IF(.NOT.ASSOCIATED(CURRENT%MAG%DOKO)) THEN
       ALLOCATE(CURRENT%MAG%DOKO)
       NULLIFY(CURRENT%MAG%DOKO%NEXT)
       CURRENT%MAG%DOKO%PARENT_FIBRE=>CURRENT
    ELSE
       D=>CURRENT%MAG%DOKO
       DO WHILE(ASSOCIATED(D%NEXT))
          D=>D%NEXT
       ENDDO
       ALLOCATE(D%NEXT)
       D=>D%NEXT
       D%PARENT_FIBRE=>CURRENT
       NULLIFY(D%NEXT)
    ENDIF

  END SUBROUTINE APPEND_POINT




  SUBROUTINE append_EMPTY_FIBRE( L )  ! Creates an empty fibre to be filled later
    implicit none
    TYPE (fibre), POINTER :: Current
    TYPE (layout), TARGET, intent(inout):: L
    L%N=L%N+1
    CALL ALLOCATE_FIBRE(Current)
    if(L%N==1) then
     current%next=> L%start
    endif
    Current % previous => L % end  ! point it to next fibre
    if(L%N>1)  THEN
       L%end%next => current      !
    ENDIF

    L % end => Current
    if(L%N==1) L%start=> Current
    if(.not.associated(current%pos)) allocate(current%pos)
    current%pos=l%n

    L%LASTPOS=L%N ;
    L%LAST=>CURRENT;
    current%parent_layout=>L
  END SUBROUTINE append_EMPTY_FIBRE

  SUBROUTINE append_NOT_SO_EMPTY_FIBRE( L )  ! Creates an empty fibre to be filled later
    implicit none
    TYPE (fibre), POINTER :: Current
    TYPE (layout), TARGET, intent(inout):: L
    L%N=L%N+1
    CALL ALLOC(Current)
    if(L%N==1) current%next=> L%start
    Current % previous => L % end  ! point it to next fibre
    if(L%N>1)  THEN
       L%end%next => current      !
    ENDIF

    L % end => Current
    if(L%N==1) L%start=> Current
    if(.not.associated(current%pos)) allocate(current%pos)
    current%pos=l%n

    L%LASTPOS=L%N ;
    L%LAST=>CURRENT;
    current%parent_layout=>L
  END SUBROUTINE append_NOT_SO_EMPTY_FIBRE

  SUBROUTINE NULL_FIBRE(CURRENT)  ! nullifies fibre content
    implicit none
    TYPE (fibre), TARGET, intent(inout):: Current
    nullify(Current%dir); !nullify(Current%P0C);nullify(Current%BETA0);
    nullify(Current%magp);nullify(Current%mag);nullify(Current%CHART);nullify(Current%PATCH);
    nullify(current%next);nullify(current%previous);
    nullify(current%PARENT_LAYOUT);
    nullify(current%T1,current%T2,current%TM);
    nullify(current%i,current%pos,current%loc);


    !    nullify(Current%P0C);
    nullify(Current%BETA0);
    nullify(Current%GAMMA0I);
    nullify(Current%GAMBET);
    nullify(Current%MASS);
    nullify(Current%AG);
    nullify(Current%CHARGE);

    nullify(current%P,current%N);

    !    nullify(current%PARENT_CHART);nullify(current%PARENT_MAG);
  END SUBROUTINE NULL_FIBRE

  SUBROUTINE ALLOCATE_FIBRE(CURRENT)   ! allocates and nullifies current's content
    implicit none
    TYPE (fibre), POINTER :: Current
    NULLIFY(CURRENT)
    ALLOCATE(Current)
    CALL NULL_FIBRE(CURRENT)
  END SUBROUTINE ALLOCATE_FIBRE

  SUBROUTINE ALLOCATE_DATA_FIBRE(CURRENT) ! Allocates pointers in fibre
    implicit none
    TYPE (fibre), TARGET, intent(inout):: Current
    ALLOCATE(Current%dir); ! ALLOCATE(Current%P0C);ALLOCATE(Current%BETA0);
    ALLOCATE(Current%magp);ALLOCATE(Current%mag);

    ALLOCATE(Current%CHART);
    ALLOCATE(Current%PATCH);
    ALLOCATE(Current%pos);

    !    ALLOCATE(Current%P0C);
    ALLOCATE(Current%BETA0);
    ALLOCATE(Current%GAMMA0I);
    ALLOCATE(Current%GAMBET);
    ALLOCATE(Current%MASS);
    ALLOCATE(Current%AG);
    ALLOCATE(Current%CHARGE);
    if(use_info) then
       allocate(Current%i)
       call alloc(Current%i)
    endif

  END SUBROUTINE ALLOCATE_DATA_FIBRE

  SUBROUTINE alloc_fibre( c ) ! Does the full allocation of fibre and initialization of internal variables
    implicit none
    type(fibre),pointer:: c
    CALL ALLOCATE_FIBRE(C)
    CALL ALLOCATE_DATA_FIBRE(C)
    c%DIR=1
    !    C%P0C = ONE
    C%BETA0 = 1.0_dp
    C%GAMMA0I = 1.0_dp
    C%GAMBET = 1.0_dp
    C%MASS = 1.0_dp
    C%MASS = A_particle
    C%CHARGE = 1


    !  c%P0C=zero
    !  c%BETA0=zero
    c%mag=0
    c%magp=0
    !if(associated(c%CHART))
    c%CHART=0
    !if(associated(c%PATCH))
    c%PATCH=0
  end SUBROUTINE alloc_fibre

  SUBROUTINE zero_fibre( c,i ) ! Does the full allocation of fibre and initialization of internal variables
    implicit none
    type(fibre),target,intent(inout):: c
    integer, intent(in) :: i
    if(i==0) then
       c%DIR=1
       !    C%P0C = ONE
       C%BETA0 = 1.0_dp
       C%GAMMA0I = 1.0_dp
       C%GAMBET = 1.0_dp
       C%MASS = 1.0_dp
       C%ag = a_particle
       C%CHARGE = 1
       !       c%P0C=zero
       !       c%BETA0=zero
       c%mag=0
       c%magp=0
       if(associated(c%CHART)) c%CHART=0
       if(associated(c%PATCH)) c%PATCH=0
    elseif(i==-1) then
       IF(ASSOCIATED(LC,c%mag%PARENT_FIBRE%PARENT_LAYOUT).or.superkill) THEN    ! ORDINARY
          IF(ASSOCIATED(c%magP)) then  !  2010_1
             c%magp=-1;
             deallocate(c%magP);
          ENDIF
          IF(ASSOCIATED(c%mag)) then   !  2010_1 changed order with above
             c%mag=-1;
             deallocate(c%mag);
          ENDIF
          IF(ASSOCIATED(c%CHART)) then  !.AND.(.NOT.ASSOCIATED(c%PARENT_CHART))) THEN
             C%CHART=-1
             deallocate(c%CHART);
          ENDIF
          IF(ASSOCIATED(c%PATCH)) then  !.AND.(.NOT.ASSOCIATED(c%PARENT_PATCH))) THEN
             C%PATCH=-1
             deallocate(c%PATCH);
          ENDIF
       ELSE   ! POINTED LAYOUT
          IF(.NOT.ASSOCIATED(c%mag%PARENT_FIBRE%CHART,c%CHART)) then
             C%CHART=-1
             deallocate(c%CHART);
          ENDIF
          IF(.NOT.ASSOCIATED(c%mag%PARENT_FIBRE%PATCH,c%PATCH)) then
             C%PATCH=-1
             deallocate(c%PATCH);
          ENDIF
       ENDIF

       IF(ASSOCIATED(c%DIR)) THEN
          deallocate(c%DIR);
       ENDIF
       !       IF(ASSOCIATED(c%P0C)) THEN
       !          deallocate(c%P0C);
       !       ENDIF
       IF(ASSOCIATED(c%BETA0)) THEN
          deallocate(c%BETA0);
       ENDIF
       IF(ASSOCIATED(c%GAMMA0I)) THEN
          deallocate(c%GAMMA0I);
       ENDIF
       IF(ASSOCIATED(c%GAMBET)) THEN
          deallocate(c%GAMBET);
       ENDIF
       IF(ASSOCIATED(c%MASS)) THEN
          deallocate(c%MASS);
       ENDIF
       IF(ASSOCIATED(c%ag)) THEN
          deallocate(c%ag);
       ENDIF
       IF(ASSOCIATED(c%CHARGE)) THEN
          deallocate(c%CHARGE);
       ENDIF

       !       IF(ASSOCIATED(C%N)) nullify(C%N)
       !       IF(ASSOCIATED(C%P)) nullify(C%P)
       nullify(C%N)
       nullify(C%P)

   !!! maybe missing per Sagan 2012.3.18       
 !    IF(ASSOCIATED(C%T1)) THEN
 !         if(associated(C%T1,C%TM)) nullify(C%TM)
 !         deallocate(C%T1);
 !         deallocate(C%T2);
 !    ENDIF
    nullify(C%T1,C%T2,C%Tm)
  !!! maybe missing per Sagan 2012.3.18       
    
       IF(ASSOCIATED(c%pos)) THEN
          deallocate(c%pos);
       ENDIF
       IF(ASSOCIATED(c%loc)) THEN
          deallocate(c%loc);
       ENDIF

       IF(ASSOCIATED(C%TM)) deallocate(C%TM);

       IF(ASSOCIATED(c%i).and.use_info) THEN
          call kill(c%i);
          deallocate(c%i);
       ENDIF

    else
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(1((1X,a72)))'
       write(6,*) "Error in zero_fibre "
       ! call !write_e(100)
    endif
  end SUBROUTINE zero_fibre

  SUBROUTINE SUPER_zero_fibre( c,i ) ! Does the full allocation of fibre and initialization of internal variables
    implicit none
    type(fibre),target,intent(inout):: c
    integer, intent(in) :: i
    if(i==0) then
       c%DIR=1
       !    C%P0C = ONE
       C%BETA0 = 1.0_dp
       C%GAMMA0I = 1.0_dp
       C%GAMBET = 1.0_dp
       C%MASS = 1.0_dp
       C%ag = a_particle
       C%CHARGE = 1
       !       c%P0C=zero
       !       c%BETA0=zero
       c%mag=0
       c%magp=0
       if(associated(c%CHART)) c%CHART=0
       if(associated(c%PATCH)) c%PATCH=0
    elseif(i==-1) then
       !       IF(ASSOCIATED(c%mag)) then  !.AND.(.NOT.ASSOCIATED(c%PARENT_MAG))) THEN
       c%mag=-1;
       deallocate(c%mag);
       !      ENDIF
       !      IF(ASSOCIATED(c%magP)) then  !.AND.(.NOT.ASSOCIATED(c%PARENT_MAG))) THEN
       c%magp=-1;
       deallocate(c%magP);
       !       ENDIF
       !       IF(ASSOCIATED(c%CHART)) then  !.AND.(.NOT.ASSOCIATED(c%PARENT_CHART))) THEN
       C%CHART=-1
       deallocate(c%CHART);
       !       ENDIF
       !       IF(ASSOCIATED(c%PATCH)) then  !.AND.(.NOT.ASSOCIATED(c%PARENT_PATCH))) THEN
       C%PATCH=-1
       deallocate(c%PATCH);
       !       ENDIF


       !      IF(ASSOCIATED(c%DIR)) THEN
       deallocate(c%DIR);
       !      ENDIF
       !      IF(ASSOCIATED(c%BETA0)) THEN
       deallocate(c%BETA0);
       !      ENDIF
       !      IF(ASSOCIATED(c%GAMMA0I)) THEN
       deallocate(c%GAMMA0I);
       !      ENDIF
       !      IF(ASSOCIATED(c%GAMBET)) THEN
       deallocate(c%GAMBET);
       !      ENDIF
       !      IF(ASSOCIATED(c%MASS)) THEN
       deallocate(c%MASS);
       deallocate(c%ag);
       !      ENDIF
       !      IF(ASSOCIATED(c%CHARGE)) THEN
       deallocate(c%CHARGE);
       !      ENDIF
       !       IF(ASSOCIATED(C%N)) nullify(C%N)
       !       IF(ASSOCIATED(C%P)) nullify(C%P)
       nullify(C%N)
       nullify(C%P)

       IF(ASSOCIATED(C%T1)) THEN
          deallocate(C%T1);
          deallocate(C%T2);
          deallocate(C%TM);
       ENDIF
       IF(ASSOCIATED(c%i)) THEN
          call kill(c%i);
          deallocate(c%i);
       ENDIF
             IF(ASSOCIATED(c%pos)) THEN
       deallocate(c%pos);
             ENDIF
       IF(ASSOCIATED(c%loc)) deallocate(c%loc);

    else
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(1((1X,a72)))'
       write(6,*) "Error in zero_fibre "
       ! call !write_e(100)
    endif
  end SUBROUTINE SUPER_zero_fibre




  SUBROUTINE dealloc_fibre( c ) ! destroys internal data  if it is not pointing (i.e. not a parent)
    implicit none
    type(fibre),pointer :: c
    IF(ASSOCIATED(C)) THEN
       CALL zero_fibre(c,-1)
       deallocate(c);
    ENDIF
  end SUBROUTINE dealloc_fibre

  SUBROUTINE super_dealloc_fibre( c ) ! destroys internal data  if it is not pointing (i.e. not a parent)
    implicit none
    type(fibre),pointer :: c
    IF(ASSOCIATED(C)) THEN
       CALL super_zero_fibre(c,-1)
       deallocate(c);
    ENDIF
  end SUBROUTINE super_dealloc_fibre

  !  MORE FUNNY APPENDING
  SUBROUTINE APPEND_FLAT( L, el, NAME )  ! points unless called "name" in which case it clones
    implicit none
    TYPE (layout), TARGET, intent(inout):: L
    TYPE (fibre), POINTER :: el
    CHARACTER(*) NAME
    CHARACTER(nlp) NAME1

    NAME1=NAME
    CALL CONTEXT(NAME1)

    IF(EL%MAG%NAME==NAME1) THEN  !FULL CLONING
       CALL APPEND(L,EL)
    ELSE ! FULL POINTING
       CALL APPEND_POINT(L,EL)
    ENDIF
  END SUBROUTINE APPEND_FLAT


  !  EUCLIDEAN ROUTINES
  SUBROUTINE CHECK_NEED_PATCH(EL1,EL2_NEXT,PREC,PATCH_NEEDED) ! check need of  PATCHES
    IMPLICIT NONE
    TYPE (FIBRE), TARGET,INTENT(IN) :: EL1
    TYPE (FIBRE),TARGET,OPTIONAL, INTENT(INOUT) :: EL2_NEXT
    TYPE (FIBRE),POINTER :: EL2
    REAL(DP)  D(3),ANG(3)
    REAL(DP) ENT(3,3),EXI(3,3),ENT0(3,3),EXI0(3,3)
    REAL(DP), POINTER,DIMENSION(:)::A,B
    INTEGER  DIR
    REAL(DP)   PREC
    INTEGER A_YZ,A_XZ
    LOGICAL(LP)  DISCRETE,ene
    INTEGER I,PATCH_NEEDED
    REAL(DP) NORM,pix(3)

    PATCH_NEEDED=0
    pix=0.0_dp
    pix(1)=pi
    DIR=1
    DISCRETE=.FALSE.
    ANG=0.0_dp
    D=0.0_dp

    IF(PRESENT(EL2_NEXT)) THEN
       EL2=>EL2_NEXT
    ELSE
       EL2=>EL1%NEXT
    ENDIF



    IF(EL1%DIR*EL2%DIR==1) THEN   !   1
       IF(EL1%DIR==1) THEN
          EXI=EL1%CHART%F%EXI
          B=>EL1%CHART%F%B
          ENT=EL2%CHART%F%ENT
          A=>EL2%CHART%F%A
          A_XZ=1;A_YZ=1;
       ELSE
          EXI=EL1%CHART%F%ENT
          exi0=exi
          call geo_rot(exi,pix,1,basis=exi0)
          B=>EL1%CHART%F%A
          ENT=EL2%CHART%F%EXI
          ent0=ent
          call geo_rot(ent,pix,1,basis=ent0)
          A=>EL2%CHART%F%B
          ! A_XZ=1;A_YZ=1;
          A_XZ=-1;A_YZ=-1;
       ENDIF
    ELSE                          !   1
       IF(EL1%DIR==1) THEN
          EXI=EL1%CHART%F%EXI
          B=>EL1%CHART%F%B
          ENT=EL2%CHART%F%EXI
          ent0=ent
          call geo_rot(ent,pix,1,basis=ent0)
          A=>EL2%CHART%F%B
          A_XZ=1;A_YZ=-1;
       ELSE
          EXI=EL1%CHART%F%ENT
          exi0=exi
          call geo_rot(exi,pix,1,basis=exi0)
          B=>EL1%CHART%F%A
          ENT=EL2%CHART%F%ENT
          A=>EL2%CHART%F%A
          A_XZ=-1;A_YZ=1;
       ENDIF
    ENDIF                     !   1

    CALL FIND_PATCH(B,EXI,A,ENT,D,ANG)

  ! if(a_xz==-1.or.a_yz==-1) then
  !   write(6,*) "  discrete ",  a_xz,A_YZ
  !  endif


    NORM=0.0_dp
    DO I=1,3
       NORM=NORM+ABS(D(I))
    ENDDO
 
    IF(NORM>=PREC) THEN
       if(global_verbose)  Write(6,'(a20,1x,3(1x,g12.5))') " Patch Translations ",d
       D=0.0_dp
       PATCH_NEEDED=PATCH_NEEDED+1
    ENDIF
    NORM=0.0_dp
    DO I=1,3
       NORM=NORM+ABS(ANG(I))
    ENDDO
    ene=(NORM<=PREC.and.(A_XZ==1.and.A_YZ==1)).or.(NORM<=PREC.and.(A_XZ==-1.and.A_YZ==-1))
    IF(.not.ene) THEN
        if(global_verbose) Write(6,'(a20,1x,3(1x,g12.5))') "  Patch rotations  ",ang
        ANG=0.0_dp
       PATCH_NEEDED=PATCH_NEEDED+10
    ENDIF


    if(ABS((EL2%MAG%P%P0C-EL1%MAG%P%P0C)/EL1%MAG%P%P0C)>PREC) then
     if(global_verbose) Write(6,'(a20,1x,1(1x,g12.5))') "  Patch Energy      ", &
      ABS((EL2%MAG%P%P0C-EL1%MAG%P%P0C)/EL1%MAG%P%P0C)
     PATCH_NEEDED=PATCH_NEEDED+100
    endif

    DISCRETE=.false.
    IF(ANG(1)/TWOPI<-0.25_dp) THEN
       write(6,*) " #1 Discrete IMPOSSIBLE ANG(1)/TWOPI<-0.25_dp"
       DISCRETE=.TRUE.
    ENDIF
    IF(ANG(1)/TWOPI>0.25_dp) THEN
       write(6,*) " #2 Discrete IMPOSSIBLE  ANG(1)/TWOPI>0.25_dp"
        DISCRETE=.TRUE.
    ENDIF
    IF(ANG(2)/TWOPI<-0.25_dp) THEN
       write(6,*) " #3 Discrete IMPOSSIBLE ANG(2)/TWOPI<-0.25_dp"
        DISCRETE=.TRUE.
    ENDIF
    IF(ANG(1)/TWOPI>0.25_dp) THEN
       write(6,*) " #4 Discrete IMPOSSIBLE  ANG(1)/TWOPI>0.25_dp"
        DISCRETE=.TRUE.
    ENDIF

    !    IF(DISCRETE) THEN
    !       WRITE(6,*)  " NO GEOMETRIC PATCHING POSSIBLE : MORE THAN 90 DEGREES BETWEEN FACES "
    !       STOP 1123
    !    ENDIF

    if(discrete) then
       PATCH_NEEDED=PATCH_NEEDED-1000
    endif

   ! norm=abs(el1%mag%p%p0c-el2%mag%p%p0c)
   ! ene=(norm>prec)

!    if(ene) then
!       PATCH_NEEDED=PATCH_NEEDED+100
!    endif

  END SUBROUTINE CHECK_NEED_PATCH

  SUBROUTINE remove_patch(r,geometry,energy) ! check need of  PATCHES
    IMPLICIT NONE
    TYPE (layout), target ::  r
    TYPE (FIBRE), pointer ::  p
    integer i
    logical(lp), optional :: geometry,energy
    logical(lp) g,e

    g=my_true
    e=my_true

    if(present(energy)) e=energy
    if(present(geometry)) g=geometry

    p=>r%start

    do i=1,r%n
       if(g) p%patch%patch=0
       if(e) p%patch%energy=0
       p=>p%next
    enddo


  end SUBROUTINE remove_patch

  SUBROUTINE FIND_PATCH_P_new(EL1,EL2_NEXT,D,ANG,DIR,ENERGY_PATCH,PREC,patching) ! COMPUTES PATCHES
    IMPLICIT NONE
    TYPE (FIBRE), target, INTENT(INOUT) :: EL1
    TYPE (FIBRE),TARGET,OPTIONAL, INTENT(INOUT) :: EL2_NEXT
    TYPE (FIBRE),POINTER :: EL2
    REAL(DP), INTENT(INOUT) :: D(3),ANG(3)
    REAL(DP) ENT(3,3),EXI(3,3),ENT0(3,3),EXI0(3,3)
    REAL(DP), POINTER,DIMENSION(:)::A,B
    INTEGER, INTENT(IN) ::  DIR
    LOGICAL(LP), OPTIONAL, INTENT(IN) ::  ENERGY_PATCH

    LOGICAL(LP), OPTIONAL, INTENT(out) ::  patching
    REAL(DP), OPTIONAL, INTENT(IN) ::  PREC
    INTEGER A_YZ,A_XZ
    LOGICAL(LP) ENE,DOIT,DISCRETE
    INTEGER LOC,I,PATCH_NEEDED
    REAL(DP) NORM,pix(3)
    PATCH_NEEDED=1
    pix=0.0_dp
    pix(1)=pi

    DISCRETE=.FALSE.
    IF(PRESENT(EL2_NEXT)) THEN
       LOC=-1
       EL2=>EL2_NEXT
    ELSE
       LOC=1
       EL2=>EL1%NEXT
    ENDIF
    ENE=.FALSE.
    IF(PRESENT(ENERGY_PATCH)) ENE=ENERGY_PATCH
    DOIT=ASSOCIATED(EL1%CHART%F).AND.ASSOCIATED(EL2%CHART%F)
    IF(DIR==1) THEN
       DOIT=DOIT.AND.(ASSOCIATED(EL2%PATCH))
    ELSE
       DOIT=DOIT.AND.(ASSOCIATED(EL1%PATCH))
    ENDIF
    IF(DOIT) THEN
       IF(EL1%DIR*EL2%DIR==1) THEN   !   1
          IF(EL1%DIR==1) THEN
             EXI=EL1%CHART%F%EXI
             B=>EL1%CHART%F%B
             ENT=EL2%CHART%F%ENT
             A=>EL2%CHART%F%A
             A_XZ=1;A_YZ=1;
          ELSE
             EXI=EL1%CHART%F%ENT
             exi0=exi
             call geo_rot(exi,pix,1,basis=exi0)
             B=>EL1%CHART%F%A
             ENT=EL2%CHART%F%EXI
             ent0=ent
             call geo_rot(ent,pix,1,basis=ent0)
             A=>EL2%CHART%F%B
             !  A_XZ=1;A_YZ=1;
             A_XZ=-1;A_YZ=-1;
          ENDIF
       ELSE                          !   1
          IF(EL1%DIR==1) THEN
             EXI=EL1%CHART%F%EXI
             B=>EL1%CHART%F%B
             ENT=EL2%CHART%F%EXI
             ent0=ent
             call geo_rot(ent,pix,1,basis=ent0)
             A=>EL2%CHART%F%B
             A_XZ=1;A_YZ=-1;
          ELSE
             EXI=EL1%CHART%F%ENT
             exi0=exi
             call geo_rot(exi,pix,1,basis=exi0)
             B=>EL1%CHART%F%A
             ENT=EL2%CHART%F%ENT
             A=>EL2%CHART%F%A
             A_XZ=-1;A_YZ=1;
          ENDIF
       ENDIF                     !   1

       CALL FIND_PATCH(B,EXI,A,ENT,D,ANG)

       IF(PRESENT(PREC)) THEN
          NORM=0.0_dp
          DO I=1,3
             NORM=NORM+ABS(D(I))
          ENDDO
          IF(NORM<=PREC) THEN
             D=0.0_dp
             PATCH_NEEDED=PATCH_NEEDED+1
          ENDIF
          NORM=0.0_dp
          DO I=1,3
             NORM=NORM+ABS(ANG(I))
          ENDDO
          IF(NORM<=PREC.and.(A_XZ==1.and.A_YZ==1)) THEN
             ANG=0.0_dp
             PATCH_NEEDED=PATCH_NEEDED+1
          ELSEIF(NORM<=PREC.and.(A_XZ==-1.and.A_YZ==-1)) THEN  ! added 2008.6.18
             ANG=0.0_dp
             PATCH_NEEDED=PATCH_NEEDED+1
          ENDIF
          IF(PATCH_NEEDED==3) THEN
             PATCH_NEEDED=0
          ELSE
             PATCH_NEEDED=1
          ENDIF
       ENDIF
       if(PRESENT(PREC)) then
          norm=abs(2*(el1%mag%p%p0c-el2%mag%p%p0c)/(el1%mag%p%p0c+el2%mag%p%p0c))
          ene=ene.and.(norm>prec)
       endif

       IF(DIR==1) THEN

          EL2%PATCH%A_X2=A_YZ
          EL2%PATCH%A_X1=A_XZ
          EL2%PATCH%A_D=D
          EL2%PATCH%A_ANG=ANG
          SELECT CASE(EL2%PATCH%PATCH)
          CASE(it0,it1)
             EL2%PATCH%PATCH=1*PATCH_NEEDED
          CASE(it2,it3)
             EL2%PATCH%PATCH=PATCH_NEEDED + 2     ! etienne 2008.05.29
          END SELECT
          IF(ENE) THEN

             SELECT CASE(EL2%PATCH%ENERGY)
             CASE(it0,it1)
                EL2%PATCH%ENERGY=1
             CASE(it2,it3)
                EL2%PATCH%ENERGY=3
             END SELECT
          ENDIF

       ELSEIF(DIR==-1) THEN

          EL1%PATCH%B_X2=A_YZ    !  BUG WAS EL2
          EL1%PATCH%B_X1=A_XZ    !
          EL1%PATCH%B_D=D
          EL1%PATCH%B_ANG=ANG
          SELECT CASE(EL1%PATCH%PATCH)
          CASE(it0,it2)
             EL1%PATCH%PATCH=2*PATCH_NEEDED
          CASE(it1,it3)
             EL1%PATCH%PATCH=2*PATCH_NEEDED + 1     ! etienne 2008.05.29
          END SELECT
          IF(ENE) THEN
             SELECT CASE(EL2%PATCH%ENERGY)
             CASE(it0,it2)
                EL1%PATCH%ENERGY=2
             CASE(it1,it3)
                EL1%PATCH%ENERGY=3
             END SELECT
          ENDIF
       ENDIF
    ELSE ! NO FRAME

       !w_p=0
       !w_p%NC=3
       !w_p%FC='(2(1X,A72,/),(1X,A72))'
       write(6,*) " NO GEOMETRIC PATCHING POSSIBLE : EITHER NO FRAMES IN PTC OR NO PATCHES "
       write(6,'(A16,1X,L1,1X,L1)')  " CHARTS 1 AND 2 ", ASSOCIATED(EL1%CHART%F), ASSOCIATED(EL2%CHART%F)
       WRITE(6,'(A16,1X,L1,1X,L1)')  "PATCHES 1 AND 2 ", ASSOCIATED(EL1%PATCH), ASSOCIATED(EL2%PATCH)
       ! call ! WRITE_I

       IF(DIR==1) THEN

          IF(ASSOCIATED(EL2%PATCH)) THEN
             IF(ENE) THEN
                SELECT CASE(EL2%PATCH%ENERGY)
                CASE(it0,it1)
                   EL2%PATCH%ENERGY=1
                CASE(it2,it3)
                   EL2%PATCH%ENERGY=3
                END SELECT
             ENDIF
          ELSE
           write(6,*) " NOT EVEN ENERGY PATCH POSSIBLE ON ELEMENT 1 "
          ENDIF

       ELSEIF(DIR==-1) THEN

          IF(ASSOCIATED(EL2%PATCH)) THEN
             IF(ENE) THEN
                SELECT CASE(EL2%PATCH%ENERGY)
                CASE(it0,it2)
                   EL1%PATCH%ENERGY=2
                CASE(it1,it3)
                   EL1%PATCH%ENERGY=3
                END SELECT
             ENDIF
          ELSE
           write(6,*) " NOT EVEN ENERGY PATCH POSSIBLE ON ELEMENT 1 "
          ENDIF
       ENDIF

    ENDIF

    DISCRETE=.false.
    IF(ANG(1)/TWOPI<-0.25_dp) THEN
       DISCRETE=.TRUE.
    ENDIF
    IF(ANG(1)/TWOPI>0.25_dp) THEN
       DISCRETE=.TRUE.
    ENDIF
    IF(ANG(2)/TWOPI<-0.25_dp) THEN
       DISCRETE=.TRUE.
    ENDIF
    IF(ANG(1)/TWOPI>0.25_dp) THEN
       DISCRETE=.TRUE.
    ENDIF

    IF(DISCRETE) THEN
       if(.not.present(patching))  write(6,*) " NO GEOMETRIC PATCHING POSSIBLE : MORE THAN 90 DEGREES BETWEEN FACES "
    ENDIF

    if(present(patching)) patching=.not.discrete

  END SUBROUTINE FIND_PATCH_P_new

  SUBROUTINE FIND_PATCH_0(EL1,EL2_NEXT,NEXT,ENERGY_PATCH,PREC,patching) ! COMPUTES PATCHES
    IMPLICIT NONE
    TYPE (FIBRE),target :: EL1
    TYPE (FIBRE),TARGET,OPTIONAL, INTENT(INOUT) :: EL2_NEXT
    TYPE (FIBRE),POINTER :: EL2
    REAL(DP)  D(3),ANG(3)
    REAL(DP), OPTIONAL :: PREC
    LOGICAL(LP), OPTIONAL, INTENT(IN) ::  NEXT,ENERGY_PATCH
    INTEGER DIR
    LOGICAL(LP) ENE,NEX
    LOGICAL(LP), OPTIONAL, INTENT(out) ::  patching
    IF(PRESENT(EL2_NEXT)) THEN
       EL2=>EL2_NEXT
    ELSE
       EL2=>EL1%NEXT
    ENDIF
    NEX=.FALSE.
    ENE=.FALSE.
    IF(PRESENT(NEXT)) NEX=NEXT
! changed 2012.12.04 fixing for bmad patch element
!    if(associated(el1,el1%parent_layout%start)) then
!     if(.not.nex) then
!      nex=my_true
!     endif
!    endif
!    if(associated(el1%next,el1%parent_layout%start)) then
! Test assoc el2%parent_layout for dummy_fibre used in Bmad
    if (associated(el2%parent_layout)) then
      if(associated(el2, el2%parent_layout%start)) then
        nex=my_false
      endif
    endif


    
    el1%PATCH%B_X1=1
    el1%PATCH%B_X2=1
    el1%PATCH%B_D=0.0_dp
    el1%PATCH%B_ANG=0.0_dp
    el1%PATCH%B_T=0.0_dp
    el1%PATCH%B_L=0.0_dp

    EL2%PATCH%A_X1=1
    EL2%PATCH%A_X2=1
    EL2%PATCH%A_D=0.0_dp
    EL2%PATCH%A_ANG=0.0_dp
    EL2%PATCH%A_T=0.0_dp
    EL2%PATCH%A_L=0.0_dp

    if(el1%PATCH%patch==3) then
       el1%PATCH%patch=1
    elseIF(el1%PATCH%patch==2) then
       el1%PATCH%patch=0
    endif

    if(el1%PATCH%energy==3) then
       el1%PATCH%ENERGY=1
    elseIF(el1%PATCH%energy==2) then
       el1%PATCH%ENERGY=0
    endif

    if(el1%PATCH%time==3) then
       el1%PATCH%time=1
    elseIF(el1%PATCH%time==2) then
       el1%PATCH%time=0
    endif


    if(EL2%PATCH%patch==3) then
       EL2%PATCH%patch=2
    elseIF(EL2%PATCH%patch==1) then
       EL2%PATCH%patch=0
    endif

    if(EL2%PATCH%energy==3) then
       EL2%PATCH%ENERGY=2
    elseIF(EL2%PATCH%energy==1) then
       EL2%PATCH%ENERGY=0
    endif

    if(EL2%PATCH%time==3) then
       EL2%PATCH%time=2
    elseIF(EL2%PATCH%time==1) then
       EL2%PATCH%time=0
    endif

    IF(PRESENT(ENERGY_PATCH)) then
       ENE=ENERGY_PATCH
    else
       if(ABS((EL2%MAG%P%P0C-EL1%MAG%P%P0C)/EL1%MAG%P%P0C)>eps_fitted) ENE=.TRUE.
    endif
    DIR=-1  ; IF(NEX) DIR=1;
    D=0.0_dp;ANG=0.0_dp;

    CALL FIND_PATCH_P_new(EL1,EL2,D,ANG,DIR,ENERGY_PATCH=ENE,prec=PREC,patching=patching)


  END SUBROUTINE FIND_PATCH_0

  SUBROUTINE FIND_PATCH_bmad_0(EL1,EL2,D,ANG) ! COMPUTES PATCHES
    IMPLICIT NONE
    TYPE (FIBRE),TARGET, INTENT(INOUT) :: EL1
    TYPE (FIBRE),TARGET, INTENT(INOUT) :: EL2
    REAL(DP)  D(3),ANG(3),ent(3,3),exi(3,3),a(3),b(3)


    if(el1%dir==1) then
     ent=el1%chart%f%exi
     a=el1%chart%f%b
    else
     ent=el1%chart%f%ent
     a=el1%chart%f%a
    endif

    if(el2%dir==1) then
     exi=el2%chart%f%ent
     b=el2%chart%f%a
    else
     exi=el2%chart%f%exi
     b=el2%chart%f%b
    endif


    call FIND_PATCH_BMAD(A,ENT,B,EXI,D,ANG)

  END SUBROUTINE FIND_PATCH_bmad_0

  ! UNIVERSE STUFF

  SUBROUTINE Set_Up_UNIVERSE( L ) ! Sets up a layout: gives a unique negative index
    implicit none
    TYPE (MAD_UNIVERSE), TARGET, intent(inout):: L
    CALL NULLIFY_UNIVERSE(L)
    ALLOCATE(L%n);
    ALLOCATE(L%SHARED);
    ALLOCATE(L%LASTPOS);
    ALLOCATE(L%NF);
    L%N=0;
    L%SHARED=0;
    L%LASTPOS=0;
    L%NF=0;
  END SUBROUTINE Set_Up_UNIVERSE

  SUBROUTINE kill_last_layout( L )  ! Destroys a layout
    implicit none
    TYPE (LAYOUT), POINTER :: Current,Current1
    TYPE (MAD_UNIVERSE), TARGET, intent(inout):: L
    nullify(current)
    nullify(current1)
    Current => L % end      ! end at the end
    !    DO WHILE (ASSOCIATED(L % end))
    Current1 => L % end      ! end at the end
    L % end => Current % previous  ! update the end before disposing
    call kill_layout(Current)
    Current => L % end     ! alias of last fibre again
    L%N=L%N-1
    deallocate(Current1)
    !   END DO
    !    call de_Set_Up_UNIVERSE(L)
  END SUBROUTINE kill_last_layout

  SUBROUTINE kill_UNIVERSE( L )  ! Destroys a layout
    implicit none
    TYPE (LAYOUT), POINTER :: Current,Current1
    TYPE (MAD_UNIVERSE), TARGET, intent(inout):: L
    nullify(current)
    nullify(current1)
    Current => L % end      ! end at the end
    DO WHILE (ASSOCIATED(L % end))
       Current1 => L % end      ! end at the end
       L % end => Current % previous  ! update the end before disposing
     !  WRITE(6,*) ' killing last layout '
       call kill_layout(Current)
     !  WRITE(6,*) ' killed last layout '
       Current => L % end     ! alias of last fibre again
       L%N=L%N-1
       deallocate(Current1)
    END DO
    call de_Set_Up_UNIVERSE(L)
  END SUBROUTINE kill_UNIVERSE

  SUBROUTINE kill_layout_in_universe( L )  ! Destroys a layout
    implicit none
    TYPE (LAYOUT), POINTER :: L,C1,c2
    TYPE (MAD_UNIVERSE), pointer ::  u

   if(.not.associated(l)) then
      write(6,*) " There is nothing to kill "
      return
    endif


     u=>l%parent_universe    
     
     if(u%nf/=0) then
      write(6,*) " You cannot kill a layout in a tied Universe "
      return
     endif

    if(associated(u%start,u%end)) then
     call kill_layout(u%start)
     call de_Set_Up_UNIVERSE(u)
     call Set_Up_UNIVERSE(u)
   !  write(6,*) " 1 "
     return
    elseif(u%n==2) then
       if(associated(l,u%start)) then
        call kill_layout(L)
        u%start=>u%end
   !  write(6,*) " start 2"
       else
        call kill_layout(L)
        u%end=>u%start
       endif
      u%n=1
   !   write(6,*) " end 2"
      return
    endif

       if(associated(l,u%start)) then
        C1=>l%next
        call kill_layout(L)
        u%start=>c1
     !    write(6,*) " start >2"
       elseif (associated(l,u%end)) then
        C1=>l%previous
        call kill_layout(L)
        u%end=>c1
     !    write(6,*) " end >2"
       else
        C1=>l%previous
        C2=>l%next
        call kill_layout(L)
        c1%next=>c2
        c2%previous=>c1        
     !    write(6,*) " middle >2"

      endif

      u%n=u%n-1
  END SUBROUTINE kill_layout_in_universe

  SUBROUTINE FIND_POS_in_universe(C,i )  ! Finds the location "i" of the fibre C in layout L
    implicit none
    INTEGER, INTENT(INOUT) :: I
    TYPE (layout), POINTER :: C
    TYPE (layout), POINTER :: P
    NULLIFY(P);
    P=>C
    I=0
    DO WHILE(ASSOCIATED(P))
       I=I+1
       P=>P%PREVIOUS
    ENDDO
  END SUBROUTINE FIND_POS_in_universe

  SUBROUTINE MOVE_TO_LAYOUT_I( L,current,i ) ! Moves current to the i^th position
    implicit none
    TYPE (LAYOUT), POINTER :: Current
    TYPE (MAD_UNIVERSE), TARGET, intent(inout):: L
    integer i,k

    nullify(current);
    Current => L%START
    IF(I<=L%N) THEN
       DO K=1,I-1
          CURRENT=>CURRENT%NEXT
       ENDDO
    ELSE
       WRITE(6,*) "FATAL ERROR IN MOVE_TO_LAYOUT_I ",I,L%N
       STOP 900
    ENDIF
  END SUBROUTINE MOVE_TO_LAYOUT_I


  SUBROUTINE MOVE_TO_LAYOUT_name( L,current,name ) ! Moves current to the i^th position
    implicit none
    TYPE (LAYOUT), POINTER :: Current
    TYPE (MAD_UNIVERSE), TARGET, intent(inout):: L
    integer i,k
    character(120) name,name1

    name1=name
    call context(name1)


    nullify(current);
    do i=1,l%n
       Current => L%START
       call context(current%name)
       IF(current%NAME==NAME1) RETURN
       IF(I<=L%N) THEN
          DO K=1,I-1
             CURRENT=>CURRENT%NEXT
             call context(current%name)
             IF(current%NAME==NAME1) RETURN
          ENDDO
       ELSE
          WRITE(6,*) "FATAL ERROR IN MOVE_TO_LAYOUT_I ",I,L%N
          STOP 900
       ENDIF
    enddo
  END SUBROUTINE MOVE_TO_LAYOUT_name

  SUBROUTINE de_Set_Up_UNIVERSE( L ) ! deallocates layout content
    implicit none
    TYPE (MAD_UNIVERSE), TARGET, intent(inout):: L
    deallocate(L%n);
    deallocate(L%SHARED);
    deallocate(L%NF);
    deallocate(L%LASTPOS);
  END SUBROUTINE de_Set_Up_UNIVERSE

  SUBROUTINE nullIFY_UNIVERSE( L ) ! Nullifies layout content,i
    implicit none
    TYPE (MAD_UNIVERSE), TARGET, intent(inout):: L
    nullify(L%N)
    nullify(L%SHARED)

    nullify(L%END )! STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING
    nullify(L%START )! STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING
    nullify(L%NF )  ! POSITION OF LAST VISITED
    nullify(L%LASTPOS )  ! POSITION OF LAST VISITED
    nullify(L%LAST )! LAST VISITED

  END SUBROUTINE nullIFY_UNIVERSE


  SUBROUTINE APPEND_EMPTY_LAYOUT( L )   ! Appoints without cloning
    implicit none
    TYPE (MAD_UNIVERSE), TARGET, intent(inout):: L
    TYPE (LAYOUT),POINTER :: current
    nullify(current);
    L%N=L%N+1

    allocate(current)
    CALL SET_UP(current)
    current%parent_universe=>L

    if(L%N==1) then
       L%start=>current
       L%end=>current
       nullify(current%previous)
       nullify(current%next)
       return
    endif
    Current % previous => L % end  ! point it to next fibre
    L % end % next => current      !

    L % end => Current

  END SUBROUTINE APPEND_EMPTY_LAYOUT


  SUBROUTINE locate_in_universe(F,i,j)
    IMPLICIT NONE
    integer i,j
    TYPE(FIBRE),pointer ::  F


    call FIND_POS(f%mag%PARENT_FIBRE%parent_layout, f%mag%PARENT_FIBRE,j )

    call FIND_POS( f%mag%PARENT_FIBRE%parent_layout,i )


  END SUBROUTINE locate_in_universe

  SUBROUTINE FIND_POS_in_layout(L, C,i )  ! Finds the location "i" of the fibre C in layout L
    implicit none
    INTEGER, INTENT(INOUT) :: I
    TYPE(LAYOUT) L
    TYPE (fibre), POINTER :: C
    TYPE (fibre), POINTER :: P
    NULLIFY(P);

    !    CALL LINE_L(L,doneit)  ! TGV
    I=0
    IF(ASSOCIATED(C,L%START)) THEN
     I=1
     RETURN
    ENDIF
    P=>L%start%NEXT
    I=2
    DO WHILE(.NOT.ASSOCIATED(P,C))
       I=I+1
       P=>P%NEXT
       if(i>1000000) then
        write(6,*) " not found in FIND_POS_in_layout "
        i=0
        exit
       endif
    ENDDO

    !    CALL RING_L(L,doneit)
  END SUBROUTINE FIND_POS_in_layout

  SUBROUTINE unify_mad_universe(M_U,N)
    implicit none
    type(MAD_UNIVERSE),TARGET :: M_U
    type(layout),pointer :: L
    integer i,k,N0
    type(fibre),pointer :: c,c0
    INTEGER, OPTIONAL :: N
    ! used in TIE_MAD_UNIVERSE
    N0=M_U%N
    IF(PRESENT(N)) N0=N

    IF(N0>M_U%N) THEN
       WRITE(6,*) " ERROR IN unify_mad_universe"
    ENDIF

    k=0
    l=>m_u%start
    do i=1,N0-1
       k=k+l%n
       l%end%N=>l%next%start
       l%next%start%P=>l%end
       l=>l%next
    enddo
    l%end%N=>m_u%start%start
    m_u%start%start%P=>l%end
    k=k+l%n

    write(6,*) "universe has ",k," fibres"
    k=0
    l=>m_u%start

    k=0
    c0=>l%start
    c=>l%start
    do while(.true.)
       k=k+1
       c=>c%N
       if(associated(c0,c)) exit
    enddo
    write(6,*) "universe has ",k," fibres"

  end  SUBROUTINE unify_mad_universe

  SUBROUTINE TIE_MAD_UNIVERSE(M_U,N)
    implicit none
    type(layout),pointer :: L
    integer i,j,N0,K
    INTEGER, OPTIONAL :: N
    type(fibre),pointer :: c
    type(MAD_UNIVERSE),TARGET :: M_U
    N0=M_U%N
    ! ties universe from  layout 1 to layout N; otherwise ties it all
    ! with new pointers fibre%N and fibre%P. (Next and previous; circular list)
    ! See move_to_name
    !  m_u%nf  the numbers of fibres tied together
    ! fibre%loc  location in the tied universed

    IF(PRESENT(N)) N0=N

    IF(N0>M_U%N) THEN
       WRITE(6,*) " ERROR IN TIE_MAD_UNIVERSE"
    ENDIF
    K=1
    l=>m_u%start
    do i=1,N0
       C=>L%START
       do j=1,L%N
          C%N=>C%NEXT
          C%P=>C%PREVIOUS
          if(.not.associated(c%loc)) allocate(c%loc)
          c%loc=k
          K=K+1
          C=>C%NEXT
       enddo
       L=>L%NEXT
    enddo
    k=k-1
    WRITE(6,*) K," FIBRES COMPUTED IN TIE_MAD_UNIVERSE"
    CALL unify_mad_universe(M_U,N)
    m_u%nf=k
    m_u%last=>m_u%start%start
    m_u%lastpos=1
  end SUBROUTINE TIE_MAD_UNIVERSE

  subroutine gUniverse_max_n(u,n)
    !use build_lattice
    implicit none
    integer n,i
    type(mad_universe), target :: u
    type(layout), pointer :: L
    n=0

    l=>u%start
    do i=1,u%n
       n=n+l%n
       l=>l%next
    enddo

  end subroutine gUniverse_max_n


  subroutine gUniverse_max_node_n(u,n)
    !use build_lattice
    implicit none
    integer n,i
    type(mad_universe), target :: u
    type(layout), pointer :: L
    n=0

    l=>u%start
    do i=1,u%n
       if(associated(l%t) ) n=n+l%t%n
       l=>l%next
    enddo

  end subroutine gUniverse_max_node_n


  SUBROUTINE move_to_name( m_u,current,name,pos,next)
    ! moves to next one in list called name in tied universe
    implicit none
    TYPE (fibre), POINTER :: Current
    TYPE (mad_universe), target :: m_u
    integer, intent(inout):: pos
    character(*), intent(in):: name
    CHARACTER(nlp) S1NAME
    integer i
    logical(lp), optional :: next
    logical(lp) ne

    logical(lp) foundit,b
    TYPE (fibre), POINTER :: p
    TYPE (fibre), POINTER :: pb
    TYPE (fibre), POINTER :: pa

    !   locates magnet with name "name"
    ! it searches back and forth

    ne=.true.
    if(present(next)) ne=next
    foundit=.false.
    b=.false.
    S1NAME=name
    CALL CONTEXT(S1name)

    nullify(p)
    p=>m_u%last
    pb=>p%p
    pa=>p%n
    if(.not.associated(p)) goto 100
    do i=1,m_u%nf/2+1
       if(i==1.and..not.ne) then
          if(p%mag%name==s1name) then
             foundit=.true.
             b=.true.
             pb=>p
             goto 100
          endif
       endif
       if(pb%mag%name==s1name) then
          foundit=.true.
          b=.true.
          goto 100
       endif
       if(pa%mag%name==s1name) then
          foundit=.true.
          goto 100
       endif
       pa=>pa%n
       pb=>pb%p
    enddo
100 continue
    if(foundit) then
       if(b) then
          current=>pb
          pos=mod_n(m_u%lastpos-i,m_u%nf)
       else
          current=>pa
          pos=mod_n(m_u%lastpos+i,m_u%nf)
       endif
       m_u%lastpos=pos
       m_u%last=>current
    else
       pos=0
       write(6,*) " did not find ",S1name, "in tied universe "
    endif
  END SUBROUTINE move_to_name

  !  THIN LENS STRUCTURE STUFF


  SUBROUTINE NULL_THIN(T)  ! nullifies THIN content
    implicit none
    TYPE (INTEGRATION_NODE), TARGET, intent(inout):: T
    NULLIFY(T%PARENT_NODE_LAYOUT)
    NULLIFY(T%PARENT_FIBRE)
    !    NULLIFY(T%BB)
    NULLIFY(T%S)
    NULLIFY(T%lost)
    NULLIFY(T%delta_rad_out)
    NULLIFY(T%delta_rad_in)
    NULLIFY(T%ref)
    !    NULLIFY(T%ORBIT)
    NULLIFY(T%a,T%ENT)
    NULLIFY(T%B,T%EXI)
    !    NULLIFY(T%BT)
    NULLIFY(T%NEXT)
    NULLIFY(T%PREVIOUS)
    NULLIFY(T%BB)
!    NULLIFY(T%T)
    !    NULLIFY(T%WORK)
    !    NULLIFY(T%USE_TPSA_MAP)
    !    NULLIFY(T%TPSA_MAP)
    !    NULLIFY(T%INTEGRATION_NODE_AFTER_MAP)
  END SUBROUTINE NULL_THIN

  SUBROUTINE ALLOCATE_THIN(CURRENT)   ! allocates and nullifies current's content
    implicit none
    TYPE (INTEGRATION_NODE), POINTER :: Current
    NULLIFY(CURRENT)
    ALLOCATE(Current)
    CALL NULL_THIN(CURRENT)

    ALLOCATE(CURRENT%S(5))
    ALLOCATE(CURRENT%ds_ac)
    ALLOCATE(CURRENT%lost)
    ALLOCATE(CURRENT%delta_rad_in)
    ALLOCATE(CURRENT%delta_rad_out)
    ALLOCATE(CURRENT%ref(4))
    CURRENT%lost=0
    CURRENT%ref=0.0_dp
    CURRENT%delta_rad_in=0.0_dp
    CURRENT%delta_rad_out=0.0_dp
    CURRENT%ds_ac=0.0_dp
    !    ALLOCATE(CURRENT%ORBIT(6))
    ALLOCATE(CURRENT%pos_in_fibre)
    ALLOCATE(CURRENT%pos)
    ALLOCATE(CURRENT%CAS)
    ALLOCATE(CURRENT%TEAPOT_LIKE)
    !    ALLOCATE(CURRENT%USE_TPSA_MAP)

    !    ALLOCATE(CURRENT%A(3),CURRENT%ENT(3,3))
    !    ALLOCATE(CURRENT%B(3),CURRENT%EXI(3,3))
    !    CURRENT%A=ZERO
    !    CURRENT%ENT=GLOBAL_FRAME
    !    CURRENT%B=ZERO
    !    CURRENT%EXI=GLOBAL_FRAME

    CURRENT%pos_in_fibre=-100
    CURRENT%pos=-100
    CURRENT%CAS=-100
    CURRENT%TEAPOT_LIKE=-100
    !    CURRENT%USE_TPSA_MAP=MY_FALSE
  END SUBROUTINE ALLOCATE_THIN

  !  SUBROUTINE ALLOCATE_NODE_MAP(CURRENT)   ! allocates and nullifies current's content
  !    implicit none
  !    TYPE (INTEGRATION_NODE), POINTER :: Current
  !    ALLOCATE(CURRENT%ORBIT(6))
  !    ALLOCATE(CURRENT%TPSA_MAP)
  !    CURRENT%USE_TPSA_MAP=MY_FALSE
  !    CURRENT%ORBIT=ZERO
  !  END SUBROUTINE ALLOCATE_NODE_MAP

  SUBROUTINE nullIFY_NODE_LAYOUT( L ) ! Nullifies layout content,i
    implicit none
    !   integer , intent(in) :: i
    TYPE (NODE_layout), TARGET, intent(inout):: L
    !   if(i==0) then
    nullify(L%INDEX)
    nullify(L%NAME)
    nullify(L%CLOSED,L%N )
    nullify(L%LASTPOS )  ! POSITION OF LAST VISITED
    nullify(L%LAST )! LAST VISITED
    !
    nullify(L%END )
    nullify(L%START )
    nullify(L%START_GROUND )! STORE THE GROUNDED VALUE OF START DURING CIRCULAR SCANNING
    nullify(L%END_GROUND )! STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING
    nullify(L%parent_LAYOUT )!
    nullify(L%ORBIT_LATTICE )!



  END SUBROUTINE nullIFY_NODE_LAYOUT

  SUBROUTINE Set_Up_NODE_LAYOUT( L ) ! Sets up a layout: gives a unique  index
    implicit none
    TYPE (NODE_LAYOUT), TARGET, intent(inout):: L
    CALL NULLIFY_NODE_LAYOUT(L)
    ALLOCATE(L%closed);  ALLOCATE(L%lastpos);ALLOCATE(L%NAME);
    ALLOCATE(L%INDEX);
    ALLOCATE(L%n);
    L%closed=.false.;
    L%N=0;
    L%lastpos=0;L%NAME='NEMO';
    NULLIFY(L%LAST)
    INDEX_node=INDEX_node+1
    L%INDEX=INDEX_node
  END SUBROUTINE Set_Up_NODE_LAYOUT

  SUBROUTINE APPEND_EMPTY_THIN( L )  ! Creates an empty fibre to be filled later
    implicit none
    TYPE (INTEGRATION_NODE), POINTER :: Current
    TYPE (NODE_LAYOUT), TARGET, intent(inout):: L
    !    LOGICAL(LP) doneit

    L%N=L%N+1
    CALL ALLOCATE_THIN(Current)
    if(L%N==1) current%next=> L%start
    Current % previous => L % end  ! point it to next fibre
    if(L%N>1)  THEN
       L%end%next => current      !
    ENDIF

    L % end => Current
    if(L%N==1) L%start=> Current

    L%LASTPOS=L%N ;
    L%LAST=>CURRENT;

  END SUBROUTINE APPEND_EMPTY_THIN


  SUBROUTINE allocate_node_frame( L )  ! Creates an empty fibre to be filled later
    implicit none
    TYPE (INTEGRATION_NODE), POINTER :: Current
    TYPE (LAYOUT), TARGET, intent(inout):: L
    integer i


    Current=>L%T%START
    do i=1,L%T%N
       IF(.NOT.ASSOCIATED(CURRENT%A)) THEN
          ALLOCATE(CURRENT%A(3),CURRENT%ENT(3,3))
          ALLOCATE(CURRENT%B(3),CURRENT%EXI(3,3))
          CURRENT%A=0.0_dp
          CURRENT%ENT=GLOBAL_FRAME
          CURRENT%B=0.0_dp
          CURRENT%EXI=GLOBAL_FRAME
       ENDIF
       Current=>CURRENT%NEXT
    ENDDO
  end SUBROUTINE allocate_node_frame

  SUBROUTINE LINE_L_THIN(L,doneit) ! makes into line temporarily
    implicit none
    TYPE (NODE_LAYOUT), TARGET, intent(inout):: L
    logical(lp) doneit
    doneit=.false.
    if(L%closed)  then
       if(associated(L%end%next)) then
          L%end%next=>L%start_ground
          doneit=.true.
       endif
       if(associated(L%start%previous)) then
          L%start%previous=>L%end_ground
       endif
    endif
  END SUBROUTINE LINE_L_THIN

  SUBROUTINE RING_L_THIN(L,doit) ! Brings back to ring if needed
    implicit none
    TYPE (NODE_LAYOUT), TARGET, intent(inout):: L
    logical(lp) doit
    if(L%closed.and.doit)  then
       if(.NOT.(associated(L%end%next))) then
          L%start_ground=>L%end%next      ! saving grounded pointer
          L%end%next=>L%start
       endif
       if(.NOT.(associated(L%start%previous))) then
          L%end_ground=>L%start%previous  ! saving grounded pointer
          L%start%previous=>L%end
       endif
    endif
  END SUBROUTINE RING_L_THIN

  SUBROUTINE  DEALLOC_INTEGRATION_NODE(T)
    IMPLICIT NONE
 !!! maybe missing per Sagan 2012.3.18       
  ! TYPE(INTEGRATION_NODE), TARGET, INTENT(INOUT) :: T
   TYPE(INTEGRATION_NODE), pointer :: T
!!! maybe missing per Sagan 2012.3.18       
 
    !    IF(ASSOCIATED(T%bb)) then
    !      CALL KILL(t%bb)
    !      DEALLOCATE(T%bb)
    !    endif

    IF(ASSOCIATED(T%TEAPOT_LIKE)) DEALLOCATE(T%TEAPOT_LIKE)
    IF(ASSOCIATED(T%delta_rad_in)) DEALLOCATE(T%delta_rad_in)
    IF(ASSOCIATED(T%delta_rad_out)) DEALLOCATE(T%delta_rad_out)
    IF(ASSOCIATED(T%ref)) DEALLOCATE(T%ref)
    IF(ASSOCIATED(T%a)) DEALLOCATE(T%a)
    IF(ASSOCIATED(T%ent)) DEALLOCATE(T%ent)
    IF(ASSOCIATED(T%b)) DEALLOCATE(T%b)
    IF(ASSOCIATED(T%exi)) DEALLOCATE(T%exi)
    IF(ASSOCIATED(T%S)) DEALLOCATE(T%S)
    IF(ASSOCIATED(T%DS_ac)) DEALLOCATE(T%DS_ac)
    IF(ASSOCIATED(T%lost)) DEALLOCATE(T%lost)
    !    IF(ASSOCIATED(T%ORBIT)) DEALLOCATE(T%ORBIT)
    IF(ASSOCIATED(T%pos_in_fibre)) DEALLOCATE(T%pos_in_fibre)
    IF(ASSOCIATED(T%POS)) DEALLOCATE(T%POS)
    IF(ASSOCIATED(T%CAS)) DEALLOCATE(T%CAS)
    IF(ASSOCIATED(T%BB)) THEN
       CALL KILL(T%BB)
       DEALLOCATE(T%BB)
    ENDIF
!    IF(ASSOCIATED(T%T)) THEN
!       CALL KILL(T%T)
!       DEALLOCATE(T%T)
!    ENDIF
    !    IF(ASSOCIATED(T%TPSA_MAP)) THEN
    !       CALL KILL(T%TPSA_MAP)
    !       DEALLOCATE(T%TPSA_MAP)
    !    ENDIF
    !    IF(ASSOCIATED(T%USE_TPSA_MAP)) DEALLOCATE(T%USE_TPSA_MAP)
    !    IF(ASSOCIATED(T%TPSA_MAP)) THEN
    !       CALL KILL(T%TPSA_MAP)
    !       DEALLOCATE(T%TPSA_MAP)
    !    ENDIF
!!! maybe missing per Sagan 2012.3.18       
          DEALLOCATE(T)
!!! maybe missing per Sagan 2012.3.18       

  END SUBROUTINE  DEALLOC_INTEGRATION_NODE

  SUBROUTINE kill_NODE_LAYOUT( L )  ! Destroys a layout
    implicit none
    TYPE (INTEGRATION_NODE), POINTER :: Current
    TYPE (NODE_LAYOUT), POINTER ::  L
    logical(lp) doneit
    IF(.NOT.ASSOCIATED(L)) RETURN
    CALL LINE_L_THIN(L,doneit)

    IF(ASSOCIATED(L%ORBIT_LATTICE)) THEN
       CALL de_Set_Up_ORBIT_LATTICE(L%ORBIT_LATTICE)  !  KILLING ORBIT LATTICE
       !(NO LINKED LIST DE_SET_UP_... = KILL_... )
       if(lielib_print(12)==1) WRITE(6,*) " ORBIT LATTICE HAS BEEN KILLED "
    ENDIF


    nullify(current)
    Current => L % end      ! end at the end
    DO WHILE (ASSOCIATED(L % end))
       L % end => Current % previous  ! update the end before disposing
       call DEALLOC_INTEGRATION_NODE(Current)
       Current => L % end     ! alias of last fibre again
       L%N=L%N-1
    END DO
    call de_Set_Up_NODE_LAYOUT(L)
    DEALLOCATE(L);
    NULLIFY(L);
  END SUBROUTINE kill_NODE_LAYOUT

  SUBROUTINE de_Set_Up_ORBIT_LATTICE( L ) ! deallocates layout content
    implicit none
    TYPE (ORBIT_LATTICE),POINTER :: L
    INTEGER I

    DO I=1,L%ORBIT_N_NODE+1
       !       CALL KILL_ORBIT_NODE(L%ORBIT_NODES,I)
       CALL KILL_ORBIT_NODE1(L%ORBIT_NODES(I))
    ENDDO
    deallocate(L%ORBIT_NODES)
    deallocate(L%ORBIT_N_NODE)
    deallocate(L%ORBIT_USE_ORBIT_UNITS)
    deallocate(L%ORBIT_WARNING)
    deallocate(L%ORBIT_P0C)
    deallocate(L%ORBIT_BETA0)
    deallocate(L%ORBIT_LMAX)
    deallocate(L%orbit_kinetic)
    deallocate(L%orbit_brho)
    deallocate(L%ORBIT_MAX_PATCH_TZ)
    deallocate(L%ORBIT_mass_in_amu)
    deallocate(L%ORBIT_gammat)
    deallocate(L%ORBIT_harmonic)
    deallocate(L%ORBIT_L)
    deallocate(L%ORBIT_CHARGE)
    deallocate(L%STATE)
    deallocate(L%orbit_energy)
    deallocate(L%ORBIT_OMEGA_after,L%orbit_gamma)
    !    deallocate(L%orbit_dppfac)
    deallocate(L%orbit_deltae)
    deallocate(L%accel)
    if(associated(L%dt)) deallocate(L%dt)
   nullify(L%tp)

    !    deallocate(L%dxs6,L%xs6,L%freqb,L%freqa,L%voltb,L%volta,L%phasa,L%phasb)
    deallocate(L)

  END SUBROUTINE de_Set_Up_ORBIT_LATTICE




  SUBROUTINE KILL_ORBIT_NODE1(ORBIT_LAYOUT_node)
    IMPLICIT NONE
    TYPE(ORBIT_NODE), TARGET, intent(inout):: ORBIT_LAYOUT_node
    DEALLOCATE(ORBIT_LAYOUT_node%LATTICE)
    DEALLOCATE(ORBIT_LAYOUT_node%DPOS)
    DEALLOCATE(ORBIT_LAYOUT_node%ENTERING_TASK)
    DEALLOCATE(ORBIT_LAYOUT_node%PTC_TASK)
    DEALLOCATE(ORBIT_LAYOUT_node%CAVITY)
  END SUBROUTINE KILL_ORBIT_NODE1

  SUBROUTINE ALLOC_ORBIT_NODE1(ORBIT_LAYOUT_node,NL)
    IMPLICIT NONE
    TYPE(ORBIT_NODE), TARGET, intent(inout):: ORBIT_LAYOUT_node
    INTEGER NL

    ALLOCATE(ORBIT_LAYOUT_node%LATTICE(1:NL))
    ALLOCATE(ORBIT_LAYOUT_node%DPOS)
    ALLOCATE(ORBIT_LAYOUT_node%ENTERING_TASK)
    ALLOCATE(ORBIT_LAYOUT_node%PTC_TASK)
    ALLOCATE(ORBIT_LAYOUT_node%CAVITY)

    ORBIT_LAYOUT_node%LATTICE(1:NL)=0.0_dp
    ORBIT_LAYOUT_node%DPOS=0
    ORBIT_LAYOUT_node%ENTERING_TASK=0
    ORBIT_LAYOUT_node%PTC_TASK=0
    ORBIT_LAYOUT_node%CAVITY=MY_FALSE

  END SUBROUTINE ALLOC_ORBIT_NODE1

  SUBROUTINE Set_Up_ORBIT_LATTICE(O,N,U)
    IMPLICIT NONE
    TYPE(ORBIT_LATTICE), TARGET, intent(inout):: O
    INTEGER N
    LOGICAL(lp)  ::  U

    if(N>0) THEN
       ALLOCATE(O%ORBIT_NODES(N))
    ELSE
       ALLOCATE(O%ORBIT_N_NODE);O%ORBIT_N_NODE=N
       ALLOCATE(O%ORBIT_USE_ORBIT_UNITS);O%ORBIT_USE_ORBIT_UNITS=U
       ALLOCATE(O%ORBIT_WARNING);O%ORBIT_WARNING=0
       ALLOCATE(O%ORBIT_OMEGA);O%ORBIT_OMEGA=1.0_dp
       ALLOCATE(O%ORBIT_P0C);O%ORBIT_P0C=1.0_dp
       ALLOCATE(O%ORBIT_BETA0);O%ORBIT_BETA0=1.0_dp
       ALLOCATE(O%ORBIT_LMAX);O%ORBIT_LMAX=0.0_dp
       ALLOCATE(O%orbit_kinetic);O%orbit_kinetic=0.0_dp
       ALLOCATE(O%ORBIT_MAX_PATCH_TZ);O%ORBIT_MAX_PATCH_TZ=0.0_dp
       ALLOCATE(O%ORBIT_mass_in_amu);O%ORBIT_mass_in_amu=0.0_dp
       ALLOCATE(O%ORBIT_gammat);O%ORBIT_gammat=0.0_dp
       ALLOCATE(O%ORBIT_L);O%ORBIT_L=0.0_dp
       ALLOCATE(O%ORBIT_harmonic);O%ORBIT_harmonic=1.0_dp
       ALLOCATE(O%ORBIT_CHARGE);O%ORBIT_CHARGE=1
       ALLOCATE(O%STATE);O%STATE=DEFAULT
       ALLOCATE(O%orbit_brho);O%orbit_brho=1.0_dp
       ALLOCATE(O%orbit_energy);O%orbit_energy=0.0_dp;
       ALLOCATE(O%orbit_gamma);O%orbit_gamma=0.0_dp;
       !    ALLOCATE(O%orbit_dppfac);O%orbit_dppfac=zero;
       ALLOCATE(O%orbit_deltae);O%orbit_deltae=0.0_dp;
       ALLOCATE(O%ORBIT_OMEGA_after);O%ORBIT_OMEGA_after=1.0_dp
       !    ALLOCATE(O%dxs6,O%xs6,O%freqb,O%freqa,O%voltb,O%volta,O%phasa,O%phasb)
       ALLOCATE(O%accel);
       nullify(O%dt);
       nullify(O%tp);
       nullify(O%parent_layout);
       !    O%freqb=zero
       !    O%freqa=zero
       !    O%voltb=zero
       !    O%volta=zero
       !    O%phasa=zero
       !    O%phasb=zero
       !    O%xs6=zero
       !    O%dxs6=zero
       O%accel=my_false
    ENDIF

    !   REAL(DP), pointer  ::  orbit_dppfac ! GET_dppfac
    !   REAL(DP), pointer  ::  orbit_deltae ! GET_deltae
    !   REAL(DP), pointer  ::  ORBIT_OMEGA_after
    !   REAL(DP), pointer  ::  freqb,freqa,voltb,volta,phasa,phasb,xs6,dxs6


  END SUBROUTINE Set_Up_ORBIT_LATTICE


  SUBROUTINE de_Set_Up_NODE_LAYOUT( L ) ! deallocates layout content
    implicit none
    TYPE (NODE_LAYOUT), TARGET, intent(inout):: L
    deallocate(L%closed);deallocate(L%lastpos);deallocate(L%NAME);
    deallocate(L%INDEX);
    deallocate(L%n);          !deallocate(L%parent_universe)   left out
    IF(ASSOCIATED(L%ORBIT_LATTICE)) deallocate(L%ORBIT_LATTICE);
  END SUBROUTINE de_Set_Up_NODE_LAYOUT

  SUBROUTINE move_to_INTEGRATION_NODE( L,current,POS ) ! Moves current to the i^th position
    implicit none
    TYPE (INTEGRATION_NODE), POINTER :: Current
    TYPE (NODE_LAYOUT), TARGET, intent(inout):: L
    integer i,k,POS,nt
    nt=l%n
    I=mod_n(POS,L%N)

    !    CALL LINE_L_THIN(L,doneit)   ! TGV

    IF(L%LASTPOS==0) THEN
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,a72,/),(1X,a72))'
       write(6,*) " L%LASTPOS=0 : ABNORMAL UNLESS LINE EMPTY"
       write(*,'(a7,i4)')" L%N = ",L%N
       ! call !write_e(-124)
    ENDIF

    nullify(current);
    Current => L%LAST

    k=L%LASTPOS

    IF(I>=L%LASTPOS) THEN

       !       DO WHILE (ASSOCIATED(Current).and.k<i)    !TGV
       DO WHILE (k<nt.and.k<i)
          k=k+1
          Current => Current % next
       END DO
    ELSE
       !       DO WHILE (ASSOCIATED(Current).and.k>i)   !TGV
       DO WHILE (k>1.and.k>i)
          k=k-1
          Current => Current % PREVIOUS
       END DO
    ENDIF
    L%LASTPOS=I; L%LAST => Current;
    !    CALL RING_L_THIN(L,doneit)
  END SUBROUTINE move_to_INTEGRATION_NODE    !TGV

  !  Beam beam stuff

  SUBROUTINE ALLOC_BEAM_BEAM_NODE(B)
    IMPLICIT NONE
    TYPE(BEAM_BEAM_NODE),POINTER :: B

    allocate(B)
    !    ALLOCATE(B%DS)
    ALLOCATE(B%S)
    ALLOCATE(B%FK)
    ALLOCATE(B%SX)
    ALLOCATE(B%SY)
    ALLOCATE(B%XM)
    ALLOCATE(B%YM)
    !    ALLOCATE(B%DPOS)
    ALLOCATE(B%bbk(2))
    !    ALLOCATE(B%mid(3,3))
    !    ALLOCATE(B%o(3))
    ALLOCATE(B%A(3))
    ALLOCATE(B%D(3))
    !    ALLOCATE(B%beta0)
    ALLOCATE(B%A_X1)
    ALLOCATE(B%A_X2)
    ALLOCATE(B%PATCH)
    B%PATCH=.FALSE.
    B%A_X1=1
    B%A_X2=1
    !    B%beta0=one
    !    B%mid=global_frame
    !    B%o=zero
    B%A=0.0_dp
    B%D=0.0_dp
    B%bbk=0.0_dp
    B%SX=1.0_dp
    B%Sy=1.0_dp
    B%XM=0.0_dp
    B%YM=0.0_dp
    !    B%DS=ZERO
    B%S=0.0_dp
    !    B%DPOS=0
    B%FK=0.0_dp
  END SUBROUTINE ALLOC_BEAM_BEAM_NODE

  SUBROUTINE KILL_BEAM_BEAM_NODE(B)
    IMPLICIT NONE
    TYPE(BEAM_BEAM_NODE),POINTER :: B

    !    DEALLOCATE(B%DS)
    DEALLOCATE(B%FK)
    DEALLOCATE(B%SX)
    DEALLOCATE(B%SY)
    DEALLOCATE(B%XM)
    DEALLOCATE(B%YM)
    DEALLOCATE(B%s)
    !    DEALLOCATE(B%DPOS)
    DEALLOCATE(B%bbk)
    !    DEALLOCATE(B%mid)
    !    DEALLOCATE(B%O)
    DEALLOCATE(B%A)
    DEALLOCATE(B%D)
    !    DEALLOCATE(B%beta0)
    DEALLOCATE(B%A_X1)
    DEALLOCATE(B%A_X2)
    DEALLOCATE(B%PATCH)

  !  DEALLOCATE(B)

  END SUBROUTINE KILL_BEAM_BEAM_NODE
!!!!  aperture stuff

  SUBROUTINE assign_aperture(p,kindaper,R,X,Y,dx,dy,pos)
    IMPLICIT NONE
    integer, optional :: pos
    integer kindaper
    REAL(DP) R(:),X,Y,dx,dy
    type(fibre), pointer :: P

    if(.NOT.ASSOCIATED(P%MAG%p%aperture)) THEN
       call alloc(P%MAG%p%aperture)
       call alloc(P%MAGP%p%aperture)
    ENDIF
    if(kindaper/=0) then
       P%MAG%p%aperture%kind = kindaper
       P%MAGP%p%aperture%kind = kindaper
       P%MAG%p%aperture%r    = R
       P%MAG%p%aperture%x    = X
       P%MAG%p%aperture%y    = y
       P%MAG%p%aperture%dx    = dX
       P%MAG%p%aperture%dy    = dy
       P%MAGP%p%aperture%r    = R
       P%MAGP%p%aperture%x    = X
       P%MAGP%p%aperture%y    = y
       P%MAGP%p%aperture%dx    = dX
       P%MAGP%p%aperture%dy    = dy
       if(present(pos)) then
          P%MAG%p%aperture%pos=pos
         P%MAGP%p%aperture%pos=pos
       endif
    endif

  end SUBROUTINE assign_aperture

  SUBROUTINE assign_one_aperture(L,mpos,kindaper,R,X,Y,dx,dy,pos)
    IMPLICIT NONE
    TYPE(LAYOUT),TARGET :: L
    integer, optional :: pos
    integer mpos,kindaper
    REAL(DP) R(:),X,Y,dx,dy
    type(fibre), pointer :: P

    call move_to(L,p,mpos)

    call assign_aperture(p,kindaper,R,X,Y,dx,dy,pos)

  end SUBROUTINE assign_one_aperture

 SUBROUTINE toggle_ONE_aperture(p)
    IMPLICIT NONE
    integer pos
    type(fibre), pointer :: P


    if(ASSOCIATED(P%MAG%p%aperture)) THEN
       P%MAG%p%aperture%kind = -P%MAG%p%aperture%kind
       P%MAGP%p%aperture%kind = P%MAG%p%aperture%kind
    ENDIF

  end SUBROUTINE toggle_ONE_aperture

 SUBROUTINE TURN_OFF_ONE_aperture(R,pos)
    IMPLICIT NONE
    TYPE(LAYOUT),TARGET :: R
    integer pos
    type(fibre), pointer :: P

    call toggle_ONE_aperture(p)



  end SUBROUTINE TURN_OFF_ONE_aperture 



 SUBROUTINE toggle_aperture(R,pos)
    IMPLICIT NONE
    TYPE(LAYOUT),TARGET :: R
    integer pos
    type(fibre), pointer :: P

    call toggle_ONE_aperture(p)

  end SUBROUTINE toggle_aperture

 SUBROUTINE turn_off_aperture(p)
    IMPLICIT NONE
    type(fibre), pointer :: P

    if(ASSOCIATED(P%MAG%p%aperture)) THEN
       IF(P%MAG%p%aperture%pos<IPOS/2) THEN
        P%MAG%p%aperture%pos = ipos+P%MAG%p%aperture%pos 
        P%MAGP%p%aperture%pos  = ipos+P%MAG%p%aperture%pos 
       ELSE
        WRITE(6,*) " ERROR APERTURE OFF ALREADY IN ", P%POS,P%MAG%NAME
       ENDIF
    ENDIF

  end SUBROUTINE turn_off_aperture

 SUBROUTINE turn_ON_aperture(p)
    IMPLICIT NONE
    type(fibre), pointer :: P

    if(ASSOCIATED(P%MAG%p%aperture)) THEN
       IF(P%MAG%p%aperture%pos>IPOS/2) THEN
         P%MAG%p%aperture%pos =  P%MAG%p%aperture%pos - IPOS
         P%MAGP%p%aperture%pos  = P%MAG%p%aperture%pos - IPOS
       ELSE
        WRITE(6,*) " ERROR APERTURE ON ALREADY IN ", P%POS,P%MAG%NAME
       ENDIF
    ENDIF
 

  end SUBROUTINE turn_ON_aperture


END MODULE S_FIBRE_BUNDLE
