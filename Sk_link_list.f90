!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file Sa_rotation_mis
MODULE S_FIBRE_BUNDLE
  !USE   S_FRAME
  USE   S_ELEMENTS
  ! Implementation of abstract data type as a linked layout
  IMPLICIT NONE

  PRIVATE kill_layout,kill_info,alloc_info,copy_info
  private dealloc_fibre,append_fibre   !, alloc_fibre public now also as alloc
  private null_it0
  private move_to_p,move_to_name,move_to_name2,move_from_to_name
  private PRINT_P_FIBRE,ANALYSE_FIBRE,PRINT_BASIC,READ_BASIC
  private PRINT_MAG,READ_MAG,PRINT_MAG_CHART,READ_MAG_CHART,PRINT_PATCH,READ_PATCH
  private PRINT_CHART,read_CHART,FIND_PATCH_p,FIND_PATCH_0
  PRIVATE INDEX
  logical(lp),TARGET :: with_chart=.true.
  logical(lp),TARGET :: with_patch=.true.
  logical(lp),TARGET :: use_info=.false.
  private zero_fibre

  INTEGER :: INDEX=0
  logical(lp),PRIVATE,PARAMETER::T=.TRUE.,F=.FALSE.

  TYPE FIBRE_CLONE
     INTEGER INDEX
     logical(lp) PARENT_LAYOUT
     logical(lp) PARENT_PATCH
     logical(lp) PARENT_CHART
     logical(lp) PARENT_MAG
  END TYPE FIBRE_CLONE

  type info
     !     character(nlp),pointer :: name
     real(sp),pointer :: s
     real(sp),pointer ::  beta(:)
     real(sp),pointer ::  fix0(:)
     real(sp),pointer ::  fix(:)
     real(sp), pointer:: pos(:)
  END type info

  TYPE FIBRE
     !  BELOW ARE THE DATA CARRIED BY THE NODE
     INTEGER,POINTER ::DIR
     REAL(DP),POINTER ::P0C,BETA0
     TYPE(PATCH),POINTER ::PATCH
     TYPE(CHART),POINTER ::CHART
     TYPE (ELEMENT), POINTER ::  MAG
     TYPE (ELEMENTP),POINTER ::  MAGP
     !  END OF DATA
     !  POINTER TO THE MAGNETS ON EACH SIDE OF THIS NODE
     TYPE (FIBRE),POINTER :: PREVIOUS
     TYPE (FIBRE),POINTER :: NEXT
     !  POINTING TO PARENT LAYOUT AND PARENT FIBRE DATA
     TYPE (LAYOUT),POINTER :: PARENT_LAYOUT
     TYPE (FIBRE),POINTER ::  PARENT_PATCH
     TYPE (FIBRE),POINTER ::  PARENT_CHART
     TYPE (FIBRE),POINTER ::  PARENT_MAG
     type(info),pointer ::i
  END TYPE FIBRE

  TYPE LAYOUT
     CHARACTER(120), POINTER ::  NAME ! IDENTIFICATION
     INTEGER, POINTER ::  INDEX,CHARGE ! IDENTIFICATION, CHARGE SIGN
     logical(lp),POINTER ::CLOSED
     INTEGER,  POINTER :: N     ! TOTAL ELEMENT IN THE CHAIN
     INTEGER,POINTER ::NTHIN  ! NUMBER IF THIN LENSES IN COLLECTION  (FOR SPEED ESTIMATES)
     REAL(DP),  POINTER :: THIN    ! PARAMETER USED FOR AUTOMATIC CUTTING INTO THIN LENS
     !POINTERS OF LINK LAYOUT
     INTEGER, POINTER :: LASTPOS   ! POSITION OF LAST VISITED
     TYPE (FIBRE), POINTER :: LAST ! LAST VISITED
     !
     TYPE (FIBRE), POINTER :: END
     TYPE (FIBRE), POINTER :: START
     TYPE (FIBRE), POINTER :: START_GROUND ! STORE THE GROUNDED VALUE OF START DURING CIRCULAR SCANNING
     TYPE (FIBRE), POINTER :: END_GROUND ! STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING
  END TYPE LAYOUT


  INTERFACE kill
     MODULE PROCEDURE kill_layout
     MODULE PROCEDURE dealloc_fibre
     MODULE PROCEDURE kill_info
  END INTERFACE

  INTERFACE alloc
     MODULE PROCEDURE set_up
     MODULE PROCEDURE alloc_fibre
     MODULE PROCEDURE alloc_info
  END INTERFACE

  INTERFACE copy
     MODULE PROCEDURE copy_info
  END INTERFACE

  INTERFACE append
     MODULE PROCEDURE append_fibre
  END INTERFACE

  INTERFACE move_to
     MODULE PROCEDURE move_to_p
     MODULE PROCEDURE move_to_name
     MODULE PROCEDURE move_to_name2
     MODULE PROCEDURE move_from_to_name
  END INTERFACE

  INTERFACE FIND_PATCH
     MODULE PROCEDURE FIND_PATCH_p
     MODULE PROCEDURE FIND_PATCH_0
  END INTERFACE





  interface assignment (=)
     MODULE PROCEDURE null_it0
     MODULE PROCEDURE zero_fibre
  end interface

CONTAINS

  SUBROUTINE alloc_info( c ) ! Does the full allocation of fibre and initialization of internal variables
    implicit none
    type(info),pointer:: c

    !     allocate(c%name); c%name=' ';
    allocate(c%s) ;c%s=0.d0;
    allocate(c%beta(40));c%beta=0.d0;
    allocate(c%fix(6));c%fix=0.d0;
    allocate(c%fix0(6));c%fix0=0.d0;
    allocate(c%pos(2));c%pos=0.d0;

  end SUBROUTINE alloc_info

  SUBROUTINE copy_info( c,d ) ! Does the full allocation of fibre and initialization of internal variables
    implicit none
    type(info) c,d

    !     d%name=c%name
    d%s=c%s
    d%beta=c%beta
    d%fix=c%fix
    d%fix0=c%fix0
    d%pos=c%pos

  end SUBROUTINE copy_info

  SUBROUTINE kill_info( c ) ! Does the full allocation of fibre and initialization of internal variables
    implicit none
    type(info),pointer:: c

    !     deallocate(c%name)
    deallocate(c%s)
    deallocate(c%fix)
    deallocate(c%fix0)
    deallocate(c%beta)
    deallocate(c%pos)

  end SUBROUTINE kill_info

  SUBROUTINE kill_mad_like( L )  ! Destroys a layout
    implicit none
    TYPE (fibre), POINTER :: c
    TYPE (layout) L
    logical(lp) doneit
    CALL LINE_L(L,doneit)
    nullify(c)
    c => L % end      ! end at the end
    DO WHILE (ASSOCIATED(L % end))
       L % end => c % previous  ! update the end before disposing
       IF(ASSOCIATED(c%mag).AND.(.NOT.ASSOCIATED(c%PARENT_MAG))) THEN
          c%mag=-1;
          deallocate(c%mag);  ! AIMIN CHANGES FOR MS4.0
       ENDIF
       IF(ASSOCIATED(c%magP).AND.(.NOT.ASSOCIATED(c%PARENT_MAG))) THEN
          c%magp=-1;
          deallocate(c%magP);  ! AIMIN CHANGES FOR MS4.0
       ENDIF
       IF(ASSOCIATED(c%CHART).AND.(.NOT.ASSOCIATED(c%PARENT_CHART))) THEN
          C%CHART=-1
          deallocate(c%CHART);  ! AIMIN CHANGES FOR MS4.0
       ENDIF
       IF(ASSOCIATED(c%PATCH).AND.(.NOT.ASSOCIATED(c%PARENT_PATCH))) THEN
          C%PATCH=-1
          deallocate(c%PATCH);  ! AIMIN CHANGES FOR MS4.0
       ENDIF
       IF(ASSOCIATED(c%DIR)) THEN
          deallocate(c%DIR);  ! AIMIN CHANGES FOR MS4.0
       ENDIF
       IF(ASSOCIATED(c%P0C)) THEN
          deallocate(c%P0C);  ! AIMIN CHANGES FOR MS4.0
       ENDIF
       IF(ASSOCIATED(c%BETA0)) THEN
          deallocate(c%BETA0);  ! AIMIN CHANGES FOR MS4.0
       ENDIF
       deallocate(c);
       c => L % end     ! alias of last fibre again
       L%N=L%N-1
    END DO
    call de_set_up(L)
  END SUBROUTINE kill_mad_like


  SUBROUTINE APPEND_mad_like( L, el )  ! Used in MAD-Like input
    implicit none
    TYPE (fibre),target :: el
    TYPE (fibre), POINTER :: Current
    TYPE (layout), TARGET:: L
    L%N=L%N+1
    CALL ALLOCATE_FIBRE(Current);
    current%mag=>el%mag
    current%magp=>el%magp
    current%CHART=>el%CHART
    current%PATCH=>el%PATCH
    if(use_info) current%i=>el%i
    current%dir=>el%dir
    current%P0C=>el%P0C
    current%BETA0=>el%BETA0

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
    TYPE (layout) L
    logical(lp) doneit
    CALL LINE_L(L,doneit)
    nullify(current)
    Current => L % end      ! end at the end
    DO WHILE (ASSOCIATED(L % end))
       L % end => Current % previous  ! update the end before disposing
       call dealloc_fibre(Current)
       Current => L % end     ! alias of last fibre again
       L%N=L%N-1
    END DO
    call de_set_up(L)
  END SUBROUTINE kill_layout


  SUBROUTINE APPEND_fibre( L, el ) ! Standard append that clones everything
    implicit none
    TYPE (fibre), intent(in) :: el
    TYPE (fibre), POINTER :: Current
    TYPE (layout), TARGET,intent(inout):: L
    logical(lp) doneit
    CALL LINE_L(L,doneit)
    L%N=L%N+1
    nullify(current)
    call alloc_fibre(current)
    call copy(el%magp,current%mag)
    call copy(current%mag,current%magp)
    call copy(el%mag,current%mag)
    if(associated(current%CHART)) call copy(el%CHART,current%CHART)
    if(associated(current%patch))call copy(el%PATCH,current%PATCH)
    if(use_info.and.associated(current%patch)) call copy(el%i,current%i)
    current%dir=el%dir
    current%P0C=el%P0C
    current%BETA0=el%BETA0

    current%PARENT_LAYOUT=>L
    if(L%N==1) current%next=> L%start
    Current % previous => L % end  ! point it to next fibre
    if(L%N>1)  THEN
       L % end % next => current      !
    ENDIF

    L % end => Current
    if(L%N==1) L%start=> Current

    L%LASTPOS=L%N ; L%LAST=>CURRENT;
    CALL RING_L(L,doneit)
  END SUBROUTINE APPEND_fibre







  SUBROUTINE FIND_POS(L, C,i )  ! Finds the location "i" of the fibre C in layout L
    implicit none
    INTEGER, INTENT(INOUT) :: I
    logical(lp) doneit
    TYPE(LAYOUT) L
    TYPE (fibre), POINTER :: C
    TYPE (fibre), POINTER :: P
    NULLIFY(P);
    P=>C
    CALL LINE_L(L,doneit)
    I=0
    DO WHILE(ASSOCIATED(P))
       I=I+1
       P=>P%PREVIOUS
    ENDDO

    CALL RING_L(L,doneit)
  END SUBROUTINE FIND_POS






  SUBROUTINE move_to_p( L,current,i ) ! Moves current to the i^th position
    implicit none
    TYPE (fibre), POINTER :: Current
    TYPE (layout) L
    integer i,k
    logical(lp) doneit

    CALL LINE_L(L,doneit)
    if(i>l%n.AND.(.NOT.L%CLOSED)) THEN
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       w_p%c(1)= " i>l%n PERMITTED ONLY IN RINGS"
       CALL WRITE_E(-123)
    ENDIF
    IF(L%LASTPOS==0) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1X,a72))'
       w_p%c(1)= " L%LASTPOS=0 : ABNORMAL UNLESS LINE EMPTY"
       write(w_p%c(2),'(a7,i4)')" L%N = ",L%N
       CALL WRITE_E(-124)
    ENDIF
    nullify(current);
    Current => L%LAST

    k=L%LASTPOS
    IF(I>=L%LASTPOS) THEN
       DO WHILE (ASSOCIATED(Current).and.k<i)
          k=k+1
          Current => Current % next
       END DO
    ELSE
       DO WHILE (ASSOCIATED(Current).and.k>i)
          k=k-1
          Current => Current % PREVIOUS
       END DO
    ENDIF
    L%LASTPOS=I; L%LAST => Current;
    CALL RING_L(L,doneit)
  END SUBROUTINE move_to_p


  SUBROUTINE move_from_to_name( L,c1,POSC1,current,name,pos)
    ! Moves from (c1,posc1) to current called "name" (posc1<=0 then finds posc1)
    implicit none
    TYPE (fibre), POINTER :: c1
    TYPE (fibre), POINTER :: Current
    TYPE (layout), intent(inout):: L
    integer, intent(inout):: pos,POSC1
    character(*), intent(in):: name
    CHARACTER(nlp) S1NAME
    integer i

    logical(lp) foundit
    TYPE (fibre), POINTER :: p

    foundit=.false.
    S1NAME=name
    CALL CONTEXT(S1name)

    nullify(p)
    IF(POSC1<=0) THEN
       CALL FIND_POS(L, C1,POSC1 )
    ENDIF
    p=>c1%NEXT
    if(.not.associated(p)) goto 100
    do i=1,l%n-1
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
       pos=mod_n(POSC1+i,l%n)
    else
       pos=0
    endif
  END SUBROUTINE move_from_to_name


  SUBROUTINE move_to_name( L,current,name,pos) ! moves to next one in list called name
    implicit none
    TYPE (fibre), POINTER :: Current
    TYPE (layout), intent(inout):: L
    integer, intent(inout):: pos
    character(*), intent(in):: name
    CHARACTER(nlp) S1NAME
    integer i

    logical(lp) foundit
    TYPE (fibre), POINTER :: p

    foundit=.false.
    S1NAME=name
    CALL CONTEXT(S1name)

    nullify(p)
    p=>l%last%next
    if(.not.associated(p)) goto 100
    do i=1,l%n-1
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
       pos=mod_n(l%lastpos+i,l%n)
       l%lastpos=pos
       l%last=>current
    else
       pos=0
    endif
  END SUBROUTINE move_to_name

  SUBROUTINE move_to_FLAT( L,current,name,POS)
    ! find in a simple flat file (Same as move_to_name but starts with oneself and scan completely)
    implicit none
    TYPE (fibre), POINTER :: Current
    TYPE (layout), intent(inout):: L
    integer, intent(inout):: pos
    character(*), intent(in):: name
    CHARACTER(nlp) S1NAME
    integer i

    logical(lp) foundit
    TYPE (fibre), POINTER :: p

    foundit=.false.
    S1NAME=name
    CALL CONTEXT(S1name)

    nullify(p)
    p=>l%last
    if(.not.associated(p)) goto 100
    do i=0,l%n-1
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
       pos=mod_n(l%lastpos+i,l%n)
       l%lastpos=pos
       l%last=>current
    else
       pos=0
    endif
  END SUBROUTINE move_to_FLAT


  SUBROUTINE move_to_name2( L,current,name,vorname,pos) !Same as move_to_name but matches name and vorname
    implicit none
    TYPE (fibre), POINTER :: Current
    TYPE (layout), intent(inout):: L
    integer, intent(inout):: pos
    character(*), intent(in):: name,vorname
    CHARACTER(nlp) S1NAME,s2name
    integer i

    logical(lp) foundit
    TYPE (fibre), POINTER :: p

    foundit=.false.
    S1NAME=name
    S2NAME=vorname
    CALL CONTEXT(S1name)
    CALL CONTEXT(S2name)

    nullify(p)
    p=>l%last%next
    if(.not.associated(p)) goto 100
    do i=1,l%n-1
       if((p%mag%name==s1name).and.(p%mag%vorname==s2name)) then
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
    endif
  END SUBROUTINE move_to_name2


  SUBROUTINE Set_Up( L ) ! Sets up a layout: gives a unique negative index
    implicit none
    TYPE (layout) L
    L=0
    ALLOCATE(L%closed);  ALLOCATE(L%lastpos);ALLOCATE(L%NAME);
    ALLOCATE(L%NTHIN);ALLOCATE(L%THIN);ALLOCATE(L%INDEX);ALLOCATE(L%CHARGE);
    ALLOCATE(L%n);
    L%closed=.false.;
    L%NTHIN=0;L%THIN=zero;
    L%N=0;
    L%lastpos=0;L%NAME='NEMO';
    INDEX=INDEX-1
    L%INDEX=INDEX
    L%CHARGE=1

  END SUBROUTINE Set_Up

  SUBROUTINE de_Set_Up( L ) ! deallocates layout content
    implicit none
    TYPE (layout) L
    deallocate(L%closed);deallocate(L%lastpos);deallocate(L%NAME);
    deallocate(L%INDEX);    deallocate(L%CHARGE);
    deallocate(L%NTHIN);deallocate(L%THIN);
    deallocate(L%n);
  END SUBROUTINE de_Set_Up


  SUBROUTINE null_it0( L,i ) ! Nullifies layout content
    implicit none
    integer , intent(in) :: i
    TYPE (layout), intent(inout) :: L
    if(i==0) then
       nullify(L%INDEX)
       nullify(L%CHARGE)
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
    else
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       w_p%c(1)= " Only =0 permitted (nullify) "
       CALL WRITE_E(100)
    endif
  END SUBROUTINE null_it0


  SUBROUTINE LINE_L(L,doneit) ! makes into line temporarily
    implicit none
    TYPE (layout) L
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
    TYPE (layout) L
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
    TYPE (layout), TARGET:: L
    type(fibre), pointer :: p
    logical(lp) doneit
    nullify(p);
    CALL LINE_L(L,doneit)
    L%N=L%N+1
    CALL ALLOCATE_FIBRE(Current);


    !  FINDING THE VERY ORIGINAL FIBRE  RECURSIVELY
    p=>el
    do while(associated(p))
       CURRENT%PARENT_MAG=>P
       p=>P%PARENT_MAG
    enddo
    p=>el
    do while(associated(p))
       CURRENT%PARENT_PATCH=>P
       p=>P%PARENT_PATCH
    enddo
    p=>el
    do while(associated(p))
       CURRENT%PARENT_CHART=>P
       p=>P%PARENT_CHART
    enddo
    !  END OF FINDING THE VERY ORIGINAL FIBRE
    CURRENT%PARENT_LAYOUT=>EL%PARENT_LAYOUT
    current%mag=>el%mag
    current%magp=>el%magp
    current%CHART=>el%CHART
    current%PATCH=>el%PATCH
    if(use_info) current%i=>el%i
    ALLOCATE(current%DIR);ALLOCATE(current%P0C);ALLOCATE(current%BETA0);
    current%dir=el%dir
    current%P0C=el%P0C
    current%BETA0=el%BETA0
    if(L%N==1) current%next=> L%start
    Current % previous => L % end  ! point it to next fibre
    if(L%N>1)  THEN
       L % end % next => current      !
    ENDIF

    L % end => Current
    if(L%N==1) L%start=> Current

    L%LASTPOS=L%N ;
    L%LAST=>CURRENT;
    CALL RING_L(L,doneit)

  END SUBROUTINE APPEND_POINT

  SUBROUTINE APPEND_EMPTY( L )  ! Creates an empty fibre to be filled later
    implicit none
    TYPE (fibre), POINTER :: Current
    TYPE (layout) L
    L%N=L%N+1
    CALL ALLOCATE_FIBRE(Current)
    if(L%N==1) current%next=> L%start
    Current % previous => L % end  ! point it to next fibre
    if(L%N>1)  THEN
       L % end % next => current      !
    ENDIF

    L % end => Current
    if(L%N==1) L%start=> Current

    L%LASTPOS=L%N ;
    L%LAST=>CURRENT;

  END SUBROUTINE APPEND_EMPTY

  SUBROUTINE NULL_FIBRE(CURRENT)  ! nullifies fibre content
    implicit none
    TYPE (fibre), POINTER :: Current
    nullify(Current%dir); nullify(Current%P0C);nullify(Current%BETA0);
    nullify(Current%magp);nullify(Current%mag);nullify(Current%CHART);nullify(Current%PATCH);
    nullify(current%next);nullify(current%previous);
    nullify(current%PARENT_LAYOUT);nullify(current%PARENT_PATCH);
    nullify(current%PARENT_CHART);nullify(current%PARENT_MAG);
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
    TYPE (fibre), POINTER :: Current
    ALLOCATE(Current%dir); ALLOCATE(Current%P0C);ALLOCATE(Current%BETA0);
    ALLOCATE(Current%magp);ALLOCATE(Current%mag);

    if(with_CHART) ALLOCATE(Current%CHART);
    if(with_patch) ALLOCATE(Current%PATCH);
  END SUBROUTINE ALLOCATE_DATA_FIBRE

  SUBROUTINE alloc_fibre( c ) ! Does the full allocation of fibre and initialization of internal variables
    implicit none
    type(fibre),pointer:: c
    CALL ALLOCATE_FIBRE(C)
    CALL ALLOCATE_DATA_FIBRE(C)
    if(use_info) then
       allocate(c%i)
       call alloc(c%i)
    endif
    c%DIR=1
    c%P0C=zero
    c%BETA0=zero
    c%mag=0
    c%magp=0
    if(associated(c%CHART)) c%CHART=0
    if(associated(c%PATCH)) c%PATCH=0
  end SUBROUTINE alloc_fibre

  SUBROUTINE zero_fibre( c,i ) ! Does the full allocation of fibre and initialization of internal variables
    implicit none
    type(fibre),intent(inout):: c
    integer, intent(in) :: i
    if(i==0) then
       c%DIR=1
       c%P0C=zero
       c%BETA0=zero
       c%mag=0
       c%magp=0
       if(associated(c%CHART)) c%CHART=0
       if(associated(c%PATCH)) c%PATCH=0
    elseif(i==-1) then
       IF(ASSOCIATED(c%i).and.use_info) THEN
          call kill(c%i);
          deallocate(c%i);
       ENDIF
       IF(ASSOCIATED(c%mag).AND.(.NOT.ASSOCIATED(c%PARENT_MAG))) THEN
          c%mag=-1;
          deallocate(c%mag);  ! AIMIN CHANGES FOR MS4.0
       ENDIF
       IF(ASSOCIATED(c%magP).AND.(.NOT.ASSOCIATED(c%PARENT_MAG))) THEN
          c%magp=-1;
          deallocate(c%magP);  ! AIMIN CHANGES FOR MS4.0
       ENDIF
       IF(ASSOCIATED(c%CHART).AND.(.NOT.ASSOCIATED(c%PARENT_CHART))) THEN
          C%CHART=-1
          deallocate(c%CHART);  ! AIMIN CHANGES FOR MS4.0
       ENDIF
       IF(ASSOCIATED(c%PATCH).AND.(.NOT.ASSOCIATED(c%PARENT_PATCH))) THEN
          C%PATCH=-1
          deallocate(c%PATCH);  ! AIMIN CHANGES FOR MS4.0
       ENDIF
       IF(ASSOCIATED(c%DIR)) THEN
          deallocate(c%DIR);  ! AIMIN CHANGES FOR MS4.0
       ENDIF
       IF(ASSOCIATED(c%P0C)) THEN
          deallocate(c%P0C);  ! AIMIN CHANGES FOR MS4.0
       ENDIF
       IF(ASSOCIATED(c%BETA0)) THEN
          deallocate(c%BETA0);  ! AIMIN CHANGES FOR MS4.0
       ENDIF

    else
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       w_p%c(1)= "Error in zero_fibre "
       CALL WRITE_E(100)
    endif
  end SUBROUTINE zero_fibre



  SUBROUTINE dealloc_fibre( c ) ! destroys internal data  if it is not pointing (i.e. not a parent)
    implicit none
    type(fibre),pointer :: c
    IF(ASSOCIATED(C)) THEN
       c=-1
       deallocate(c);
    ENDIF
  end SUBROUTINE dealloc_fibre

  !  MORE FUNNY APPENDING
  SUBROUTINE APPEND_FLAT( L, el, NAME )  ! points unless called "name" in which case it clones
    implicit none
    TYPE (layout), TARGET :: L
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

  ! switching routines
  SUBROUTINE switch_to_kind7( el )  ! switch to kind7
    implicit none
    TYPE (fibre), POINTER :: el
    ! This routines switches to kind7 (not exact) from kind2,10,16
    select case(el%mag%kind)
    case(kind10,kind16,kind2)
       el%magp=-1
       el%mag%L=el%mag%p%ld
       el%mag%p%lc=el%mag%p%ld
       el%mag%p%exact=.false.
       el%magp=0
    end select

    select case(el%mag%kind)
    case(kind10)
       EL%MAG%TP10=-1
       deallocate(EL%MAG%TP10)
       el%mag%kind=KIND7
       CALL SETFAMILY(EL%MAG)
       CALL COPY(EL%MAG,EL%MAGP)
    case(kind16)
       EL%MAG%k16=-1
       deallocate(EL%MAG%k16)
       el%mag%kind=KIND7
       CALL SETFAMILY(EL%MAG)
       CALL COPY(EL%MAG,EL%MAGP)
    case(KIND2)
       el%mag%kind=KIND7
       CALL SETFAMILY(EL%MAG)
       CALL COPY(EL%MAG,EL%MAGP)
    end select

  END SUBROUTINE switch_to_kind7

  !  Euclidean routines
  SUBROUTINE FIND_PATCH_p(EL1,EL2,D,ANG,dir,energy_patch) ! computes patches
    implicit none
    TYPE (fibre), intent(inout) :: EL1,EL2
    real(dp), intent(inout) :: D(3),ANG(3)
    real(dp) ENT(3,3),EXI(3,3)
    real(dp), POINTER,dimension(:)::A,B
    integer, intent(in) ::  dir
    logical(lp), optional, intent(in) ::  energy_patch
    integer i
    logical(lp) ene,DOIT

    ene=.false.
    if(present(energy_patch)) ene=energy_patch
    DOIT=associated(EL1%CHART%f).and.associated(EL2%CHART%f)
    IF(DIR==1) THEN
       DOIT=DOIT.AND.(associated(EL2%PATCH))
    ELSE
       DOIT=DOIT.AND.(associated(EL1%PATCH))
    ENDIF
    if(DOIT) then
       if(el1%DIR*EL2%DIR==1) THEN   !   1
          IF(EL1%DIR==1) THEN
             EXI=EL1%CHART%f%EXI
             B=>EL1%CHART%f%B
             ENT=EL2%CHART%f%ENT
             A=>EL2%CHART%f%A
          ELSE
             EXI=EL1%CHART%f%ENT
             B=>EL1%CHART%f%A
             ENT=EL2%CHART%f%EXI
             A=>EL2%CHART%f%B
          ENDIF
       ELSE                          !   1
          IF(EL1%DIR==1) THEN
             EXI=EL1%CHART%f%EXI
             B=>EL1%CHART%f%B
             ENT=EL2%CHART%f%EXI
             A=>EL2%CHART%f%B
          ELSE
             EXI=EL1%CHART%f%ENT
             B=>EL1%CHART%f%A
             ENT=EL2%CHART%f%ENT
             A=>EL2%CHART%f%A
          ENDIF
          do i=1,3
             EXI(1,i)=-EXI(1,i)
             EXI(3,i)=-EXI(3,i)
          enddo
       ENDIF                     !   1

       CALL FIND_PATCH(B,EXI,A,ENT,D,ANG)

       if(dir==1) then

          EL2%PATCH%A_D=D
          EL2%PATCH%A_ANG=ANG
          EL2%PATCH%PATCH=.TRUE.
          EL2%PATCH%energy=ene
          if(EL1%dir*EL2%dir==-1) then
             EL2%PATCH%A_D(3)=-EL2%PATCH%A_D(3)
             EL2%PATCH%A_ANG(1)=-EL2%PATCH%A_ANG(1)
             EL2%PATCH%A_ANG(2)=-EL2%PATCH%A_ANG(2)
          endif

       elseif(dir==-1) then

          EL1%PATCH%b_D=D
          EL1%PATCH%b_ANG=ANG
          EL1%PATCH%PATCH=.TRUE.
          EL1%PATCH%energy=ene
          if(EL1%dir*EL2%dir==-1) then
             EL1%PATCH%b_D(3)=-EL1%PATCH%b_D(3)
             EL1%PATCH%b_ANG(1)=-EL1%PATCH%b_ANG(1); EL1%PATCH%b_ANG(2)=-EL1%PATCH%b_ANG(2)
          endif
       endif
    else ! no frame

       w_p=0
       w_p%nc=3
       w_p%fc='(2(1X,a72,/),(1X,a72))'
       w_p%c(1)= " No Geometric Patching possible : either no frames in PTC or no patches "
       write(w_p%c(2),'(a16,1x,L1,1x,L1)')  " charts 1 and 2 ", associated(EL1%CHART%f), associated(EL2%CHART%f)
       write(w_p%c(3),'(a16,1x,L1,1x,L1)')  "patches 1 and 2 ", associated(EL1%PATCH), associated(EL2%PATCH)
       CALL WRITE_I

       if(dir==1) then

          IF(ASSOCIATED(EL2%PATCH)) THEN
             EL2%PATCH%energy=ene
          ELSE
             w_p=0
             w_p%nc=1
             w_p%fc='((1X,a72))'
             w_p%c(1)= " Not even Energy patch possible on element 2 "
             CALL WRITE_I
          ENDIF

       elseif(dir==-1) then

          IF(ASSOCIATED(EL2%PATCH)) THEN
             EL1%PATCH%energy=ene
          ELSE
             w_p=0
             w_p%nc=1
             w_p%fc='((1X,a72))'
             w_p%c(1)= " Not even Energy patch possible on element 1 "
             CALL WRITE_I
          ENDIF
       endif

    endif

  end SUBROUTINE FIND_PATCH_p

  SUBROUTINE FIND_PATCH_0(EL1,EL2,next,energy_patch) ! computes patches
    implicit none
    TYPE (fibre), intent(inout) :: EL1,EL2
    real(dp)  D(3),ANG(3)
    logical(lp), optional, intent(in) ::  next,energy_patch
    integer dir
    logical(lp) ene,nex

    nex=.false.
    ene=.false.
    if(present(next)) nex=next
    if(present(energy_patch)) ene=energy_patch
    dir=-1  ; if(nex) dir=1;
    d=zero;ang=zero;
    call FIND_PATCH(EL1,EL2,D,ANG,dir,energy_patch=ene)


  end SUBROUTINE FIND_PATCH_0

  ! Writing and reading routines

  SUBROUTINE PRINT_FLAT( L, filen )
    IMPLICIT NONE
    TYPE(LAYOUT), INTENT(INOUT):: L
    character(*)  filen
    INTEGER mf,i
    logical(lp) DONEIT,TEST
    TYPE(FIBRE_CLONE) A
    type(fibre), pointer :: p
    mf=NEWFILE
    open(unit=mf,file=filen,status='unknown',recl=255)
    CALL PRINT_BASIC(L,MF)
    CALL LINE_L(L,doneit)
    nullify(p);

    P=>L%START
    I=0
    DO WHILE(ASSOCIATED(P))
       I=I+1
       A=ANALYSE_FIBRE(P)
       IF(A%PARENT_LAYOUT) THEN
          IF(A%INDEX/=L%INDEX) THEN
             w_p=0
             w_p%nc=1
             w_p%fc='((1X,a72))'
             w_p%c(1)= " DIFFERENT LAYOUT NOT YET SUPPORTED "
             CALL WRITE_E(-999)
          ELSE
             TEST=A%PARENT_CHART.OR.A%PARENT_MAG.OR.A%PARENT_PATCH
             WRITE(MF,'(1X,A35,1X,i4)') " "
             IF(TEST) THEN
                WRITE(MF,'(1X,A35,1X,i4)') "POINTED FIBRE PARTIALLY OR TOTALLY ",I
                WRITE(MF,'(1X,I4,4(1X,L1))') A%INDEX,A%PARENT_LAYOUT,A%PARENT_PATCH,A%PARENT_CHART,A%PARENT_MAG
                CALL PRINT_P_FIBRE( L,P,A,MF )
             ELSE                         !1234567890123456789012345678901234567890
                WRITE(MF,'(1X,A15,1X,i4)') "ORIGINAL FIBRE ",I
                CALL PRINT_P_FIBRE( L,P,A,MF )
             ENDIF
          ENDIF
       ENDIF
       P=>P%NEXT
    ENDDO

    CALL ring_L(L,doneit)

    mf=CLOSEFILE


  END SUBROUTINE PRINT_FLAT

  SUBROUTINE READ_FLAT( L, filen )
    IMPLICIT NONE
    TYPE(LAYOUT), INTENT(INOUT):: L
    character(*)  filen
    CHARACTER(25) A1
    CHARACTER(1) A0
    INTEGER mf,N_TEMP
    logical(lp) DONEIT
    TYPE(FIBRE_CLONE) A
    mf=NEWFILE
    open(unit=mf,file=filen)

    CALL SET_UP(L)
    CALL READ_BASIC(L,MF)
    N_TEMP=L%N
    L%N=0

1001 continue
    READ(MF,"(a25)",END=1000) A1
    READ(MF,*) A1  ;
    CALL CONTEXT(A1) ;
    A0=A1(1:1)
    SELECT CASE(a0)
    CASE('O')
       CALL APPEND_EMPTY(L)
       CALL ALLOCATE_DATA_FIBRE(L%END)
       A=FIBRE_CLONE(0,F,F,F,F)
       CALL READ_FULL_FIBRE( L,A,MF )
       L%END%MAGP=0
       CALL COPY(L%END%MAG,L%END%MAGP)
    CASE('P')
       CALL APPEND_EMPTY(L)
       ALLOCATE(L%END%DIR);ALLOCATE(L%END%P0C);ALLOCATE(L%END%BETA0);
       READ(MF,'(1X,I4,4(1X,L1))') A%INDEX,A%PARENT_LAYOUT,A%PARENT_PATCH,A%PARENT_CHART,A%PARENT_MAG
       CALL READ_FULL_FIBRE( L,A,MF )


    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       w_p%c(1)=  "PROBLEMS IN READ_FLAT "
       CALL WRITE_E(101)
    END SELECT
    goto 1001
1000 CONTINUE
    w_p=0
    w_p%nc=3
    w_p%fc='((1X,a72,/),(1X,A120,/),(1X,a72))'
    w_p%c(1)=  "FLAT FILE "
    w_p%c(2)=  FILEN
    w_p%c(3)=  " READ "
    CALL WRITE_I
    mf=CLOSEFILE

    DONEIT=.true.
    if(l%closed) then
       call ring_l(l,DONEIT)
    else
       call line_l(l,DONEIT)
    endif

  END SUBROUTINE READ_FLAT

  SUBROUTINE READ_FULL_FIBRE( L,A,MF )
    implicit none
    TYPE(FIBRE_CLONE) A
    TYPE (LAYOUT) L
    TYPE (fibre),POINTER :: el,P
    INTEGER POS(3),MF
    logical(lp) MIS
    CHARACTER(nlp) NAME
    character(vp) VORNAME
    CHARACTER(20) A1,A3

    POS=0
    NULLIFY(EL);EL=>L%END;
    READ(MF,200) A1,EL%DIR,A3,EL%P0C,EL%BETA0,NAME,VORNAME
200 FORMAT(A6,(1X,I4,1X),A15,2(1X,G20.14),2(1X,A16))
    IF(A%PARENT_MAG) THEN
       READ(MF,100) A1,POS(1)
       CALL MOVE_TO(L,P,POS(1))
       EL%MAG=>P%MAG
       EL%MAGP=>P%MAGP
       MIS=EL%MAG%MIS
    ELSE
       EL%CHART=0
       EL%PATCH=0
       CALL READ_MAG( EL%MAG, MF )
       EL%MAG%NAME=NAME
       EL%MAG%VORNAME=VORNAME
       MIS=EL%MAG%MIS
    ENDIF
    IF(A%PARENT_CHART) THEN
       READ(MF,100) A1,POS(2)
       CALL MOVE_TO(L,P,POS(2))
       EL%CHART=>P%CHART
    else
       CALL read_CHART(EL%CHART,MIS,MF)
    ENDIF
    IF(A%PARENT_PATCH) THEN
       READ(MF,100) A1,POS(3)
       CALL MOVE_TO(L,P,POS(3))
       EL%PATCH=>P%PATCH
    else
       CALL READ_PATCH(EL%PATCH,MF)
    ENDIF

    IF(.NOT.A%PARENT_MAG) THEN
       EL%mag%p%f%mid(:,:) = EL%CHART%f%mid(:,:)
       EL%mag%p%f%o(:)     = EL%CHART%f%o(:)
       EL%mag%p%f%ent(:,:) = EL%CHART%f%ent(:,:)
       EL%mag%p%f%a(:)     = EL%CHART%f%a(:)
       EL%mag%p%f%exi(:,:) = EL%CHART%f%exi(:,:)
       EL%mag%p%f%b(:)     = EL%CHART%f%b(:)
    ENDIF
100 FORMAT(1X,A20,1x,i4)
  END SUBROUTINE READ_FULL_FIBRE




  FUNCTION ANALYSE_FIBRE( el )
    implicit none
    TYPE(FIBRE_CLONE) ANALYSE_FIBRE
    TYPE (fibre),POINTER :: el
    ANALYSE_FIBRE=FIBRE_CLONE(0,F,F,F,F)
    IF(ASSOCIATED(EL%PARENT_LAYOUT)) THEN
       ANALYSE_FIBRE%PARENT_LAYOUT=T
       ANALYSE_FIBRE%INDEX=EL%PARENT_LAYOUT%INDEX
    ENDIF
    IF(ASSOCIATED(EL%PARENT_CHART)) ANALYSE_FIBRE%PARENT_CHART=T
    IF(ASSOCIATED(EL%PARENT_PATCH)) ANALYSE_FIBRE%PARENT_PATCH=T
    IF(ASSOCIATED(EL%PARENT_MAG))   ANALYSE_FIBRE%PARENT_MAG  =T
  END FUNCTION ANALYSE_FIBRE

  SUBROUTINE PRINT_P_FIBRE( L,el,A,MF )
    implicit none
    TYPE(FIBRE_CLONE) A
    TYPE (LAYOUT) L
    TYPE (fibre),POINTER :: el
    INTEGER POS(3),MF
    logical(lp) MIS
    POS=0
    WRITE(MF,200) "DIR = ",EL%DIR," P0C & BETA0 = ",EL%P0C,EL%BETA0,EL%MAG%NAME,EL%MAG%VORNAME
200 FORMAT(A6,(1X,I4,1X),A15,2(1X,G20.14),2(1X,A16))
    IF(A%PARENT_MAG) THEN
       CALL FIND_POS(L,EL%PARENT_MAG,POS(1))
       WRITE(MF,100) "ORIGINAL ELEMENT AT ",POS(1)
       MIS=EL%PARENT_MAG%MAG%MIS
    ELSE
       MIS=EL%MAG%MIS
       CALL PRINT_MAG(EL%MAG,MF)
    ENDIF
    IF(A%PARENT_CHART) THEN
       CALL FIND_POS(L,EL%PARENT_CHART,POS(2))
       WRITE(MF,100) "ORIGINAL CHART AT   ",POS(2)
    else
       CALL PRINT_CHART(EL%CHART,MIS,MF)
    ENDIF
    IF(A%PARENT_PATCH) THEN
       CALL FIND_POS(L,EL%PARENT_PATCH,POS(3))
       WRITE(MF,100) "ORIGINAL PATCH AT   ",POS(3)
    else
       CALL PRINT_PATCH(EL%PATCH,MF)
    ENDIF
100 FORMAT(1X,A20,1x,i4)

  END SUBROUTINE PRINT_P_FIBRE

  SUBROUTINE PRINT_BASIC( L, MF )
    IMPLICIT NONE
    TYPE(LAYOUT), INTENT(INOUT):: L
    INTEGER MF
    WRITE(MF,*) " BASIC INFORMATION ABOUT THE LAYOUT "
    WRITE(MF,'(A120)') L%NAME
    WRITE(MF,101) " INDEX = ",L%INDEX," CLOSED = ",L%CLOSED," N = ", L%N," CHARGE = ",L%CHARGE
    WRITE(MF,102) " THIN = ",L%THIN," NTHIN = ",L%NTHIN
    !  WRITE(MF,103) " CIRCUMFERENCE = ",L%CIRCUMFERENCE
    !  WRITE(MF,104) " ENERGY = ",L%ENERGY," KINETIC = ",L%KINETIC," P0C = ",L%P0C," BRHO = ",L%BRHO," BETA0 = ",L%BETA0
    WRITE(MF,105) " USER DEFINED METHOD = ",NEW_METHOD
101 FORMAT(1X,A9,I4,1X,A10,L1,1X,A5,I4,1X,A10,I4)
102 FORMAT(1X,A10,(1X,G20.14),1X,A9,I4)
103 FORMAT(A18,(1X,G20.14))
105 FORMAT(A23,(1X,L1))
104 FORMAT(A10,(1X,G20.14),A11,(1X,G20.14),A7,(1X,G20.14),A8,(1X,G20.14),A9,(1X,G20.14))
    IF(NEW_METHOD) then
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       w_p%c(1)=  " New method not yet supported in PRINT_BASIC "
       CALL WRITE_E(-934)
    endif
  END SUBROUTINE PRINT_BASIC

  SUBROUTINE READ_BASIC( L, MF )
    IMPLICIT NONE
    TYPE(LAYOUT), INTENT(INOUT):: L
    CHARACTER(25) A1,A2,A3,A4
    !    CHARACTER(25) A5
    INTEGER MF

    READ(MF,*) A1    !" BASIC INFORMATION ABOUT THE LAYOUT "
    READ(MF,'(A120)') L%NAME
    READ(MF,101) A1,L%INDEX,A2,L%CLOSED,A3, L%N,A4,L%CHARGE
    READ(MF,102) A1,L%THIN,A2,L%NTHIN
    !  READ(MF,103) A1,L%CIRCUMFERENCE
    !  READ(MF,104) A1,L%ENERGY,A2,L%KINETIC,A3,L%P0C,A4,L%BRHO,A5,L%BETA0
    READ(MF,105) A1,NEW_METHOD
101 FORMAT(1X,A9,I4,1X,A10,L1,1X,A5,I4,1X,A10,I4)
102 FORMAT(1X,A10,(1X,G20.14),1X,A9,I4)
103 FORMAT(A18,(1X,G20.14))
105 FORMAT(A23,(1X,L1))
104 FORMAT(A10,(1X,G20.14),A11,(1X,G20.14),A7,(1X,G20.14),A8,(1X,G20.14),A9,(1X,G20.14))
    IF(NEW_METHOD) STOP 934

  END SUBROUTINE READ_BASIC

  SUBROUTINE PRINT_MAG( EL, MF )
    IMPLICIT NONE
    TYPE(ELEMENT), INTENT(INOUT)::EL
    INTEGER MF
    CALL PRINT_MAG_CHART(EL,MF)
    SELECT CASE(EL%KIND)
    CASE(kind0)
    case(KIND1)
    case(KIND2)
       WRITE(MF,100) " MAD FRINGE GARBAGE ",EL%FINT,EL%HGAP,EL%H1,EL%H2
    case(KIND3)
    case(KIND4)
       WRITE(MF,101) "VOLTAGE= ",el%VOLT," FREQUENCY = ",el%FREQ," PHASE = ",el%PHAS,&
            &" THIN= ",el%THIN," FRINGE = ",el%c4%FRINGE
    case(KIND5)
       WRITE(MF,102) " B_SOL = ",EL%B_SOL
    case(KIND6)
       WRITE(MF,100) " MAD FRINGE GARBAGE ",EL%FINT,EL%HGAP,EL%H1,EL%H2
    case(KIND7)
       WRITE(MF,100) " MAD FRINGE GARBAGE ",EL%FINT,EL%HGAP,EL%H1,EL%H2
    case(KIND8)
    case(KIND9)
    case(KIND10)
       WRITE(MF,100) " MAD FRINGE GARBAGE ",EL%FINT,EL%HGAP,EL%H1,EL%H2
       WRITE(MF,103) " DRIFTKICK = ",EL%TP10%DRIFTKICK
    CASE(KIND11:KIND14)
    CASE(KIND15)
       WRITE(MF,104) " SEPTUM VOLTAGE = ",EL%VOLT
    CASE(KIND16)
       WRITE(MF,100) " MAD FRINGE GARBAGE ",EL%FINT,EL%HGAP,EL%H1,EL%H2
       WRITE(MF,105) " DRIFTKICK = ",EL%K16%DRIFTKICK, " LIKEMAD = ",EL%K16%LIKEMAD
    CASE(KIND17)
       WRITE(MF,102) " B_SOL = ",EL%B_SOL
    case(KINDFITTED)
    case(KINDUSER1)
       CALL PRINT_USER(EL%U1,MF)
    case(KINDUSER2)
       CALL PRINT_USER(EL%U2,MF)
    case default
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       WRITE(w_p%c(1),'((1X,I4,1X,A28))')  el%kind," not supported in print_mag "
       CALL WRITE_E(0)
    END SELECT
102 FORMAT(A9,(1X,G20.14))
101 FORMAT(A9,(1X,G20.14),A14,(1X,G20.14),A10,(1X,G20.14),A7,L1,A10,L1)
100 FORMAT(A20,4(1X,G20.14))
103 FORMAT(A13,(1X,L1))
104 FORMAT(A18,(1X,G20.14))
105 FORMAT(A13,(1X,L1),A11,(1X,L1))

  END SUBROUTINE PRINT_MAG

  SUBROUTINE READ_MAG( EL, MF )
    IMPLICIT NONE
    TYPE(ELEMENT), INTENT(INOUT)::EL
    INTEGER MF
    CHARACTER(20) A1,A2,A3,A4,A5
    logical(lp) FRINGE
    CALL READ_MAG_CHART(EL,MF)
    SELECT CASE(EL%KIND)
    CASE(kind0)
       CALL SETFAMILY(EL)
    case(KIND1)
       CALL SETFAMILY(EL)
    case(KIND2)
       READ(MF,100) A1,EL%FINT,EL%HGAP,EL%H1,EL%H2
       CALL SETFAMILY(EL)
    case(KIND3)
       CALL SETFAMILY(EL)
    case(KIND4)
       ALLOCATE(EL%VOLT,EL%FREQ,EL%PHAS,EL%DELTA_E,EL%THIN)
       READ(MF,101) A1,el%VOLT,A2,el%FREQ,A3,el%PHAS,&
            &A4,el%THIN,A5,FRINGE
       CALL SETFAMILY(EL)
       EL%C4%FRINGE=FRINGE;
    case(KIND5)
       ALLOCATE(EL%B_SOL)
       READ(MF,102) A1,EL%B_SOL
       CALL SETFAMILY(EL)
    case(KIND6)
       READ(MF,100) A1,EL%FINT,EL%HGAP,EL%H1,EL%H2
       CALL SETFAMILY(EL)
    case(KIND7)
       READ(MF,100) A1,EL%FINT,EL%HGAP,EL%H1,EL%H2
       CALL SETFAMILY(EL)
    case(KIND8)
       CALL SETFAMILY(EL)
    case(KIND9)
       CALL SETFAMILY(EL)
    case(KIND10)
       CALL SETFAMILY(EL)
       READ(MF,100) A1,EL%FINT,EL%HGAP,EL%H1,EL%H2
       READ(MF,103) A1,EL%TP10%DRIFTKICK
    CASE(KIND11:KIND14)
       CALL SETFAMILY(EL)
    CASE(KIND15)
       ALLOCATE(EL%VOLT)
       READ(MF,104) A1,EL%VOLT
       CALL SETFAMILY(EL)
    CASE(KIND16)
       CALL SETFAMILY(EL)
       READ(MF,100) A1,EL%FINT,EL%HGAP,EL%H1,EL%H2
       READ(MF,105) A1,EL%K16%DRIFTKICK, A2,EL%K16%LIKEMAD
    CASE(KIND17)
       ALLOCATE(EL%B_SOL)
       READ(MF,102) A1,EL%B_SOL
       CALL SETFAMILY(EL)
    case(KINDFITTED)
       CALL SETFAMILY(EL)
    case(KINDUSER1)
       CALL SETFAMILY(EL)
       CALL READ_USER(EL%U1,MF)
    case(KINDUSER2)
       CALL SETFAMILY(EL)
       CALL READ_USER(EL%U2,MF)
    case default
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       WRITE(w_p%c(1),'((1X,I4,1X,A28))')  el%kind," not supported in print_mag "
       CALL WRITE_E(0)
    END SELECT
102 FORMAT(A9,(1X,G20.14))
101 FORMAT(A9,(1X,G20.14),A14,(1X,G20.14),A10,(1X,G20.14),A7,L1,A10,L1)
100 FORMAT(A20,4(1X,G20.14))
103 FORMAT(A13,(1X,L1))
104 FORMAT(A18,(1X,G20.14))
105 FORMAT(A13,(1X,L1),A11,(1X,L1))
  END SUBROUTINE READ_MAG

  SUBROUTINE PRINT_MAG_CHART( EL, MF )
    IMPLICIT NONE
    TYPE(ELEMENT), INTENT(INOUT)::EL
    INTEGER MF
    WRITE(MF,'(4(1X,A8,1X,I4),1X,A14,L1,1X,A9,L1)') " NMUL = ",EL%P%NMUL," METH = ",EL%P%METHOD,"  NST = ",&
         EL%P%NST," KIND = ",EL%KIND," PERMFRINGE = ",EL%PERMFRINGE," EXACT = ",EL%P%EXACT
    WRITE(MF,'(4(1X,A5,1X,G20.14),(1X,A8,1X,G20.14))') " L = ",EL%L ,"LD = ",EL%P%LD  ,"B0 = ",EL%P%B0    ,&
         "LC = ",EL%P%LC,"TILTD = ",EL%P%TILTD
    WRITE(MF,'(1X,A6,4(1X,G20.14),(1X,A7,1X,G20.14,1X,G20.14))') "P0C = ",EL%P%P0C,EL%P%BETA0,EL%P%GAMMA0I,&
         EL%P%GAMBET  ,"EDGE = ",EL%P%EDGE(1),EL%P%EDGE(2)

    IF(EL%P%NMUL>0) WRITE(MF,110) EL%BN,EL%AN
    IF(EL%MIS) THEN
       WRITE(MF,'(A35,2(1X,L1,1X),2(A7,3(1X,G20.14)) )') " MISALIGNEMENTS AND INTERNAL FRAME ",EL%MIS,&
            EL%EXACTMIS,"EL%D = ",EL%D,"EL%R = ",EL%R
       WRITE(MF,1001) EL%P%f%A, EL%P%f%ENT
       WRITE(MF,1001) EL%P%f%O, EL%P%f%MID
       WRITE(MF,1001) EL%P%f%B, EL%P%f%EXI
    ELSE
       WRITE(MF,'(A35,2(1X,L1,1X))')" MISALIGNEMENTS AND INTERNAL FRAME ",EL%MIS,EL%EXACTMIS
    ENDIF

110 FORMAT(12(1X,G20.14))
1001 FORMAT(12(1X,G20.14))
  END SUBROUTINE PRINT_MAG_CHART

  SUBROUTINE READ_MAG_CHART( EL, MF )
    IMPLICIT NONE
    TYPE(ELEMENT), INTENT(INOUT)::EL
    INTEGER MF,NMUL,METHOD,NST,KIND
    logical(lp) PERMFRINGE,EXACT
    CHARACTER(20) A1,A2,A3,A4,A5
    CHARACTER(35) A6
    real(dp) D(3),R(3),L
    READ(MF,'(4(1X,A8,1X,I4),1X,A14,L1,1X,A9,L1)') A1,NMUL,A2,METHOD,A3,NST,A4,KIND,A5,PERMFRINGE,A6,EXACT
    EL=NMUL;EL%P%METHOD=METHOD;EL%P%NST=NST;EL%KIND=KIND;EL%PERMFRINGE=PERMFRINGE;EL%P%EXACT=EXACT;
    READ(MF,'(4(1X,A5,1X,G20.14),(1X,A8,1X,G20.14))') A5,L,A1,EL%P%LD  ,A2,EL%P%B0    ,A3,EL%P%LC,A4,EL%P%TILTD
    EL%L=L;
    READ(MF,'(1X,A6,4(1X,G20.14),(1X,A7,1X,G20.14,1X,G20.14))') A1,EL%P%P0C,EL%P%BETA0,EL%P%GAMMA0I,EL%P%GAMBET &
         ,A2,EL%P%EDGE(1),EL%P%EDGE(2)
    IF(EL%P%NMUL>0) READ(MF,110) EL%BN,EL%AN

    READ(MF,'(A35,2(1X,L1,1X),2(A7,3(1X,G20.14)) )')A6,EL%MIS,EL%EXACTMIS,A1,D,A2,R
    IF(EL%MIS) THEN
       IF(.NOT.ASSOCIATED(EL%D) ) ALLOCATE(EL%D(3))
       IF(.NOT.ASSOCIATED(EL%R) ) ALLOCATE(EL%R(3))
       EL%D=D;EL%R=R;
       READ(MF,1001) EL%P%f%A, EL%P%f%ENT
       READ(MF,1001) EL%P%f%O, EL%P%f%MID
       READ(MF,1001) EL%P%f%B, EL%P%f%EXI
    ENDIF

110 FORMAT(12(1X,G20.14))
1001 FORMAT(12(1X,G20.14))
  END SUBROUTINE READ_MAG_CHART

  SUBROUTINE PRINT_PATCH( EL, MF )
    IMPLICIT NONE
    TYPE(PATCH), INTENT(INOUT)::EL
    INTEGER MF
    WRITE(MF,*) "************** THE PATCHES **************"
    WRITE(MF,100) " POSITION PATCH = ",EL%PATCH," POSITION PATCH = ",EL%ENERGY," TIME PATCH = ",EL%TIME
    IF(EL%PATCH) THEN
       WRITE(MF,101) " A_D = ",EL%A_D ," A_ANG = ",EL%A_ANG
       WRITE(MF,101) " B_D = ",EL%B_D ," B_ANG = ",EL%B_ANG
    ENDIF
    IF(EL%TIME)  WRITE(MF,102) " A_T = ", EL%A_T," B_T = ",EL%B_T
100 FORMAT(A18,L1,A18,L1,A14,L1)
101 FORMAT(A7,3(1X,G20.14),A9,3(1X,G20.14))
102 FORMAT(A7,(1X,G20.14),A7,(1X,G20.14))
  END SUBROUTINE PRINT_PATCH

  SUBROUTINE READ_PATCH( EL, MF )
    IMPLICIT NONE
    TYPE(PATCH), INTENT(INOUT)::EL
    INTEGER MF
    CHARACTER(20) A1,A2,A3
    READ(MF,*) A1
    READ(MF,100) A1,EL%PATCH,A2,EL%ENERGY,A3,EL%TIME
    IF(EL%PATCH) THEN
       READ(MF,101) A1,EL%A_D ,A2,EL%A_ANG
       READ(MF,101) A1,EL%B_D ,A2,EL%B_ANG
    ENDIF
    IF(EL%TIME)  READ(MF,102) A1, EL%A_T,A2,EL%B_T
100 FORMAT(A18,L1,A18,L1,A14,L1)
101 FORMAT(A7,3(1X,G20.14),A9,3(1X,G20.14))
102 FORMAT(A7,(1X,G20.14),A7,(1X,G20.14))
  END SUBROUTINE READ_PATCH

  SUBROUTINE PRINT_CHART( EL,MIS, MF )
    IMPLICIT NONE
    TYPE(CHART), INTENT(INOUT)::EL
    logical(lp) MIS
    INTEGER MF
    WRITE(MF,*) "************** THE LOCAL FRAME **************"
    WRITE(MF,100) "L = ",EL%L  ,"ALPHA = ",EL%ALPHA    ,"A_XY = ",EL%A_XY
    WRITE(MF,101) EL%f%A, EL%f%ENT
    WRITE(MF,101) EL%f%O, EL%f%MID
    WRITE(MF,101) EL%f%B, EL%f%EXI
    IF(MIS) THEN
       WRITE(MF,103) " D_IN  = ",EL%D_IN ,  " ANG_IN  = ",EL%ANG_IN
       WRITE(MF,103) " D_OUT = ",EL%D_OUT , " ANG_OUT = ",EL%ANG_OUT
    ENDIF
103 FORMAT(A9,3(1X,G20.14),A11,3(1X,G20.14))
100 FORMAT(1X,A4,1X,G20.14,1X,A8,1X,G20.14,1X,A7,1X,G20.14)
101 FORMAT(12(1X,G20.14))
  END SUBROUTINE PRINT_CHART

  SUBROUTINE read_CHART( EL, MIS,MF )
    IMPLICIT NONE
    TYPE(CHART), INTENT(INOUT)::EL
    INTEGER MF
    logical(lp) MIS
    character(11) A1,A2,A3
    read(MF,*) a1
    READ(MF,100) A1,EL%L  ,A2,EL%ALPHA    ,A3,EL%A_XY
    read(MF,101) EL%f%A, EL%f%ENT
    read(MF,101) EL%f%O, EL%f%MID
    read(MF,101) EL%f%B, EL%f%EXI
    IF(MIS) THEN
       read(MF,103) A1,EL%D_IN ,  A1,EL%ANG_IN
       read(MF,103) A1,EL%D_OUT , A1,EL%ANG_OUT
    ENDIF
103 FORMAT(A9,3(1X,G20.14),A11,3(1X,G20.14))
100 FORMAT(1X,A4,1X,G20.14,1X,A8,1X,G20.14,1X,A7,1X,G20.14)
101 FORMAT(12(1X,G20.14))
  END SUBROUTINE read_CHART



END MODULE S_FIBRE_BUNDLE
