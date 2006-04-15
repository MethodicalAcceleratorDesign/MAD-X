!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file A_SCRATCH_SIZE.F90

MODULE S_FIBRE_BUNDLE
  USE S_DEF_ELEMENT
  !  USE   S_ELEMENTS
  ! Implementation of abstract data type as a linked layout
  IMPLICIT NONE
  public

  PRIVATE kill_layout,kill_info,alloc_info,copy_info
  private dealloc_fibre,append_fibre   !, alloc_fibre public now also as alloc
  !  private null_it0
  private move_to_p,move_to_name,move_to_name2,move_from_to_name
  PRIVATE append_EMPTY_FIBRE
  PRIVATE FIND_PATCH_0
  PRIVATE FIND_PATCH_p_new
  PRIVATE INDEX
  private FIND_POS_in_universe, FIND_POS_in_layout

  logical(lp),TARGET :: use_info=.false.
  private zero_fibre
  INTEGER :: INDEX=0
  logical(lp),PRIVATE,PARAMETER::T=.TRUE.,F=.FALSE.


  INTERFACE kill
     MODULE PROCEDURE kill_layout
     MODULE PROCEDURE dealloc_fibre
     MODULE PROCEDURE kill_info
  END INTERFACE

  INTERFACE alloc
     !     MODULE PROCEDURE set_up
     MODULE PROCEDURE alloc_fibre
     MODULE PROCEDURE alloc_info
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
     MODULE PROCEDURE move_to_name
     MODULE PROCEDURE move_to_name2
     MODULE PROCEDURE move_from_to_name
  END INTERFACE

  INTERFACE FIND_PATCH
     MODULE PROCEDURE FIND_PATCH_0
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
    type(info), intent(inout):: c

    allocate(c%s) ;c%s=zero;
    allocate(c%beta(40));c%beta=zero;
    allocate(c%fix(6));c%fix=zero;
    allocate(c%fix0(6));c%fix0=zero;
    allocate(c%pos(2));c%pos=zero;


  end SUBROUTINE alloc_info

  SUBROUTINE copy_info( c,d ) ! Does the full allocation of fibre and initialization of internal variables
    implicit none
    type(info), intent(in)::c
    type(info),  intent(inout)::d

    !   d%name=c%name
    d%s=c%s
    d%beta=c%beta
    d%fix=c%fix
    d%fix0=c%fix0
    d%pos=c%pos

  end SUBROUTINE copy_info

  SUBROUTINE kill_info( c ) ! Does the full allocation of fibre and initialization of internal variables
    implicit none
    type(info), intent(inout):: c

    !   deallocate(c%name)
    deallocate(c%s)
    deallocate(c%fix)
    deallocate(c%fix0)
    deallocate(c%beta)
    deallocate(c%pos)

  end SUBROUTINE kill_info

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
    !    current%P0C=>el%P0C
    !    current%BETA0=>el%BETA0

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
    !    current%P0C=el%P0C
    !    current%BETA0=el%BETA0

    current%PARENT_LAYOUT=>L
    current%mag%PARENT_FIBRE=>current
    current%magp%PARENT_FIBRE=>current
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







  SUBROUTINE FIND_POS_in_layout(L, C,i )  ! Finds the location "i" of the fibre C in layout L
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
  END SUBROUTINE FIND_POS_in_layout





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
    CALL NULLIFY_LAYOUT(L)
    ALLOCATE(L%closed);  ALLOCATE(L%lastpos);ALLOCATE(L%NAME);ALLOCATE(L%HARMONIC_NUMBER);
    ALLOCATE(L%NTHIN);ALLOCATE(L%THIN);ALLOCATE(L%INDEX);ALLOCATE(L%CHARGE);
    ALLOCATE(L%n);
    L%closed=.false.;
    L%NTHIN=0;L%THIN=zero;
    L%N=0;
    L%lastpos=0;L%NAME='NEMO';
    INDEX=INDEX-1
    L%INDEX=INDEX
    L%CHARGE=1
    L%HARMONIC_NUMBER=0
  END SUBROUTINE Set_Up




  SUBROUTINE de_Set_Up( L ) ! deallocates layout content
    implicit none
    TYPE (layout) L
    deallocate(L%closed);deallocate(L%lastpos);deallocate(L%NAME);deallocate(L%HARMONIC_NUMBER);
    deallocate(L%INDEX);    deallocate(L%CHARGE);
    deallocate(L%NTHIN);deallocate(L%THIN);
    deallocate(L%n);          !deallocate(L%parent_universe)   left out
  END SUBROUTINE de_Set_Up


  SUBROUTINE nullIFY_LAYOUT( L ) ! Nullifies layout content,i
    implicit none
    !   integer , intent(in) :: i
    TYPE (layout), intent(inout) :: L
    !   if(i==0) then
    nullify(L%INDEX)
    nullify(L%CHARGE)
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
    !    w_p=0
    !    w_p%nc=1
    !    w_p%fc='(1((1X,a72)))'
    !    w_p%c(1)= " Only =0 permitted (nullify) "
    !    CALL WRITE_E(100)
    ! endif
  END SUBROUTINE nullIFY_LAYOUT



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
    !    type(fibre), pointer :: p
    logical(lp) doneit
    !    nullify(p);
    CALL LINE_L(L,doneit)
    L%N=L%N+1
    CALL ALLOCATE_FIBRE(Current);
    ALLOCATE(Current%PATCH);

    !  FINDING THE VERY ORIGINAL FIBRE  RECURSIVELY
    !    p=>el
    !    do while(associated(p))
    !       CURRENT%PARENT_MAG=>P
    !       p=>P%PARENT_MAG
    !    enddo
    !    p=>el
    !    do while(associated(p))
    !       CURRENT%PARENT_PATCH=>P
    !       p=>P%PARENT_PATCH
    !    enddo
    !    p=>el
    !    do while(associated(p))
    !       CURRENT%PARENT_CHART=>P
    !       p=>P%PARENT_CHART
    !    enddo
    !  END OF FINDING THE VERY ORIGINAL FIBRE
    CURRENT%PARENT_LAYOUT=>EL%PARENT_LAYOUT
    current%mag=>el%mag
    current%magp=>el%magp
    current%CHART=>el%CHART
    current%PATCH=0   ! new patches always belong to fibre  ! this was the error Weishi
    if(use_info) current%i=>el%i
    ALLOCATE(current%DIR)   !;ALLOCATE(current%P0C);ALLOCATE(current%BETA0);
    current%dir=el%dir
    !    current%P0C=el%P0C
    !    current%BETA0=el%BETA0
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




  SUBROUTINE append_EMPTY_FIBRE( L )  ! Creates an empty fibre to be filled later
    implicit none
    TYPE (fibre), POINTER :: Current
    TYPE (layout) L
    L%N=L%N+1
    CALL ALLOCATE_FIBRE(Current)
    if(L%N==1) current%next=> L%start
    Current % previous => L % end  ! point it to next fibre
    if(L%N>1)  THEN
       L%end%next => current      !
    ENDIF

    L % end => Current
    if(L%N==1) L%start=> Current

    L%LASTPOS=L%N ;
    L%LAST=>CURRENT;

  END SUBROUTINE append_EMPTY_FIBRE

  SUBROUTINE NULL_FIBRE(CURRENT)  ! nullifies fibre content
    implicit none
    TYPE (fibre), POINTER :: Current
    nullify(Current%dir); !nullify(Current%P0C);nullify(Current%BETA0);
    nullify(Current%magp);nullify(Current%mag);nullify(Current%CHART);nullify(Current%PATCH);
    nullify(current%next);nullify(current%previous);
    nullify(current%PARENT_LAYOUT);
    !    nullify(current%PARENT_PATCH);
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
    TYPE (fibre), POINTER :: Current
    ALLOCATE(Current%dir); ! ALLOCATE(Current%P0C);ALLOCATE(Current%BETA0);
    ALLOCATE(Current%magp);ALLOCATE(Current%mag);

    ALLOCATE(Current%CHART);
    ALLOCATE(Current%PATCH);
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
    !  c%P0C=zero
    !  c%BETA0=zero
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
       !       c%P0C=zero
       !       c%BETA0=zero
       c%mag=0
       c%magp=0
       if(associated(c%CHART)) c%CHART=0
       if(associated(c%PATCH)) c%PATCH=0
    elseif(i==-1) then
       IF(ASSOCIATED(c%i).and.use_info) THEN
          call kill(c%i);
          deallocate(c%i);
       ENDIF
       IF(ASSOCIATED(c%mag)) then  !.AND.(.NOT.ASSOCIATED(c%PARENT_MAG))) THEN
          c%mag=-1;
          deallocate(c%mag);
       ENDIF
       IF(ASSOCIATED(c%magP)) then  !.AND.(.NOT.ASSOCIATED(c%PARENT_MAG))) THEN
          c%magp=-1;
          deallocate(c%magP);
       ENDIF
       IF(ASSOCIATED(c%CHART)) then  !.AND.(.NOT.ASSOCIATED(c%PARENT_CHART))) THEN
          C%CHART=-1
          deallocate(c%CHART);
       ENDIF
       IF(ASSOCIATED(c%PATCH)) then  !.AND.(.NOT.ASSOCIATED(c%PARENT_PATCH))) THEN
          C%PATCH=-1
          deallocate(c%PATCH);
       ENDIF
       IF(ASSOCIATED(c%DIR)) THEN
          deallocate(c%DIR);
       ENDIF
       !      IF(ASSOCIATED(c%P0C)) THEN
       !         deallocate(c%P0C);
       !      ENDIF
       !      IF(ASSOCIATED(c%BETA0)) THEN
       !         deallocate(c%BETA0);
       !      ENDIF

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
    case(kind10,kind16,kind2,KIND20)
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
    case(kind16,KIND20)
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

  !  EUCLIDEAN ROUTINES

  SUBROUTINE FIND_PATCH_P_new(EL1,EL2_NEXT,D,ANG,DIR,ENERGY_PATCH,PREC) ! COMPUTES PATCHES
    IMPLICIT NONE
    TYPE (FIBRE), INTENT(INOUT) :: EL1
    TYPE (FIBRE),TARGET,OPTIONAL, INTENT(INOUT) :: EL2_NEXT
    TYPE (FIBRE),POINTER :: EL2
    REAL(DP), INTENT(INOUT) :: D(3),ANG(3)
    REAL(DP) ENT(3,3),EXI(3,3)
    REAL(DP), POINTER,DIMENSION(:)::A,B
    INTEGER, INTENT(IN) ::  DIR
    LOGICAL(LP), OPTIONAL, INTENT(IN) ::  ENERGY_PATCH
    REAL(DP), OPTIONAL, INTENT(IN) ::  PREC
    INTEGER A_YZ,A_XZ
    LOGICAL(LP) ENE,DOIT,DISCRETE
    INTEGER LOC,I,PATCH_NEEDED
    REAL(DP) NORM,pix(3)
    PATCH_NEEDED=1
    pix=zero
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
             call geo_rot(exi,pix,1,basis=exi)
             B=>EL1%CHART%F%A
             ENT=EL2%CHART%F%EXI
             call geo_rot(ent,pix,1,basis=ent)
             A=>EL2%CHART%F%B
             A_XZ=-1;A_YZ=-1;
          ENDIF
       ELSE                          !   1
          IF(EL1%DIR==1) THEN
             EXI=EL1%CHART%F%EXI
             B=>EL1%CHART%F%B
             ENT=EL2%CHART%F%EXI
             call geo_rot(ent,pix,1,basis=ent)
             A=>EL2%CHART%F%B
             A_XZ=1;A_YZ=-1;
          ELSE
             EXI=EL1%CHART%F%ENT
             call geo_rot(exi,pix,1,basis=exi)
             B=>EL1%CHART%F%A
             ENT=EL2%CHART%F%ENT
             A=>EL2%CHART%F%A
             A_XZ=-1;A_YZ=1;
          ENDIF
       ENDIF                     !   1

       CALL FIND_PATCH(B,EXI,A,ENT,D,ANG)

       IF(PRESENT(PREC)) THEN
          NORM=ZERO
          DO I=1,3
             NORM=NORM+ABS(D(I))
          ENDDO
          IF(NORM<=PREC) THEN
             D=ZERO
             PATCH_NEEDED=PATCH_NEEDED+1
          ENDIF
          NORM=ZERO
          DO I=1,3
             NORM=NORM+ABS(ANG(I))
          ENDDO
          IF(NORM<=PREC.and.(A_XZ==1.and.A_YZ==1)) THEN
             ANG=ZERO
             PATCH_NEEDED=PATCH_NEEDED+1
          ENDIF
          IF(PATCH_NEEDED==3) THEN
             PATCH_NEEDED=0
          ELSE
             PATCH_NEEDED=1
          ENDIF
       ENDIF

       IF(DIR==1) THEN

          EL2%PATCH%A_X2=A_YZ
          EL2%PATCH%A_X1=A_XZ
          EL2%PATCH%A_D=D
          EL2%PATCH%A_ANG=ANG
          SELECT CASE(EL2%PATCH%PATCH)
          CASE(0,1)
             EL2%PATCH%PATCH=1*PATCH_NEEDED
          CASE(2,3)
             EL2%PATCH%PATCH=3*PATCH_NEEDED
          END SELECT
          IF(ENE) THEN
             SELECT CASE(EL2%PATCH%ENERGY)
             CASE(0,1)
                EL2%PATCH%ENERGY=1
             CASE(2,3)
                EL2%PATCH%ENERGY=3
             END SELECT
          ENDIF

       ELSEIF(DIR==-1) THEN

          EL1%PATCH%B_X2=A_YZ    !  BUG WAS EL2
          EL1%PATCH%B_X1=A_XZ    !
          EL1%PATCH%B_D=D
          EL1%PATCH%B_ANG=ANG
          SELECT CASE(EL1%PATCH%PATCH)
          CASE(0,2)
             EL1%PATCH%PATCH=2*PATCH_NEEDED
          CASE(1,3)
             EL1%PATCH%PATCH=3*PATCH_NEEDED
          END SELECT
          IF(ENE) THEN
             SELECT CASE(EL2%PATCH%ENERGY)
             CASE(0,2)
                EL1%PATCH%ENERGY=2
             CASE(1,3)
                EL1%PATCH%ENERGY=3
             END SELECT
          ENDIF
       ENDIF
    ELSE ! NO FRAME

       W_P=0
       W_P%NC=3
       W_P%FC='(2(1X,A72,/),(1X,A72))'
       W_P%C(1)= " NO GEOMETRIC PATCHING POSSIBLE : EITHER NO FRAMES IN PTC OR NO PATCHES "
       WRITE(W_P%C(2),'(A16,1X,L1,1X,L1)')  " CHARTS 1 AND 2 ", ASSOCIATED(EL1%CHART%F), ASSOCIATED(EL2%CHART%F)
       WRITE(W_P%C(3),'(A16,1X,L1,1X,L1)')  "PATCHES 1 AND 2 ", ASSOCIATED(EL1%PATCH), ASSOCIATED(EL2%PATCH)
       CALL WRITE_I

       IF(DIR==1) THEN

          IF(ASSOCIATED(EL2%PATCH)) THEN
             IF(ENE) THEN
                SELECT CASE(EL2%PATCH%ENERGY)
                CASE(0,1)
                   EL2%PATCH%ENERGY=1
                CASE(2,3)
                   EL2%PATCH%ENERGY=3
                END SELECT
             ENDIF
          ELSE
             W_P=0
             W_P%NC=1
             W_P%FC='((1X,A72))'
             W_P%C(1)= " NOT EVEN ENERGY PATCH POSSIBLE ON ELEMENT 2 "
             CALL WRITE_I
          ENDIF

       ELSEIF(DIR==-1) THEN

          IF(ASSOCIATED(EL2%PATCH)) THEN
             IF(ENE) THEN
                SELECT CASE(EL2%PATCH%ENERGY)
                CASE(0,2)
                   EL1%PATCH%ENERGY=2
                CASE(1,3)
                   EL1%PATCH%ENERGY=3
                END SELECT
             ENDIF
          ELSE
             W_P=0
             W_P%NC=1
             W_P%FC='((1X,A72))'
             W_P%C(1)= " NOT EVEN ENERGY PATCH POSSIBLE ON ELEMENT 1 "
             CALL WRITE_I
          ENDIF
       ENDIF

    ENDIF

    DISCRETE=.false.
    IF(ANG(1)/TWOPI<-c_0_25) THEN
       DISCRETE=.TRUE.
    ENDIF
    IF(ANG(1)/TWOPI>c_0_25) THEN
       DISCRETE=.TRUE.
    ENDIF
    IF(ANG(2)/TWOPI<-c_0_25) THEN
       DISCRETE=.TRUE.
    ENDIF
    IF(ANG(1)/TWOPI>c_0_25) THEN
       DISCRETE=.TRUE.
    ENDIF

    IF(DISCRETE) THEN
       W_P=0
       W_P%NC=1
       W_P%FC='(2(1X,A72,/),(1X,A72))'
       W_P%C(1)= " NO GEOMETRIC PATCHING POSSIBLE : MORE THAN 90 DEGREES BETWEEN FACES "
       CALL WRITE_I
    ENDIF


  END SUBROUTINE FIND_PATCH_P_new

  SUBROUTINE FIND_PATCH_0(EL1,EL2_NEXT,NEXT,ENERGY_PATCH,PREC) ! COMPUTES PATCHES
    IMPLICIT NONE
    TYPE (FIBRE), INTENT(INOUT) :: EL1
    TYPE (FIBRE),TARGET,OPTIONAL, INTENT(INOUT) :: EL2_NEXT
    TYPE (FIBRE),POINTER :: EL2
    REAL(DP)  D(3),ANG(3)
    REAL(DP), OPTIONAL :: PREC
    LOGICAL(LP), OPTIONAL, INTENT(IN) ::  NEXT,ENERGY_PATCH
    INTEGER DIR
    LOGICAL(LP) ENE,NEX

    IF(PRESENT(EL2_NEXT)) THEN
       EL2=>EL2_NEXT
    ELSE
       EL2=>EL1%NEXT
    ENDIF

    NEX=.FALSE.
    ENE=.FALSE.
    IF(PRESENT(NEXT)) NEX=NEXT
    IF(PRESENT(ENERGY_PATCH)) then
       ENE=ENERGY_PATCH
    else
       if(ABS((EL2%MAG%P%P0C-EL1%MAG%P%P0C)/EL1%MAG%P%P0C)>eps_fitted) ENE=.TRUE.
    endif
    DIR=-1  ; IF(NEX) DIR=1;
    D=ZERO;ANG=ZERO;

    CALL FIND_PATCH_P_new(EL1,EL2,D,ANG,DIR,ENERGY_PATCH=ENE,prec=PREC)


  END SUBROUTINE FIND_PATCH_0




  ! UNIVERSE STUFF

  SUBROUTINE Set_Up_UNIVERSE( L ) ! Sets up a layout: gives a unique negative index
    implicit none
    TYPE (MAD_UNIVERSE) L
    CALL NULLIFY_UNIVERSE(L)
    ALLOCATE(L%n);
    ALLOCATE(L%SHARED);
    L%N=0;
    L%SHARED=0;
  END SUBROUTINE Set_Up_UNIVERSE

  SUBROUTINE kill_UNIVERSE( L )  ! Destroys a layout
    implicit none
    TYPE (LAYOUT), POINTER :: Current,Current1
    TYPE (MAD_UNIVERSE) L
    nullify(current)
    nullify(current1)
    Current => L % end      ! end at the end
    DO WHILE (ASSOCIATED(L % end))
       Current1 => L % end      ! end at the end
       L % end => Current % previous  ! update the end before disposing
       call kill_layout(Current)
       Current => L % end     ! alias of last fibre again
       L%N=L%N-1
       deallocate(Current1)
    END DO
    call de_Set_Up_UNIVERSE(L)
  END SUBROUTINE kill_UNIVERSE

  SUBROUTINE FIND_POS_in_universe(C,i )  ! Finds the location "i" of the fibre C in layout L
    implicit none
    INTEGER, INTENT(INOUT) :: I
    logical(lp) doneit
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
    TYPE (MAD_UNIVERSE) L
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

  SUBROUTINE de_Set_Up_UNIVERSE( L ) ! deallocates layout content
    implicit none
    TYPE (MAD_UNIVERSE) L
    deallocate(L%n);
    deallocate(L%SHARED);
  END SUBROUTINE de_Set_Up_UNIVERSE

  SUBROUTINE nullIFY_UNIVERSE( L ) ! Nullifies layout content,i
    implicit none
    TYPE (MAD_UNIVERSE), intent(inout) :: L
    nullify(L%N)
    nullify(L%SHARED)

    nullify(L%END )! STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING
    nullify(L%START )! STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING

  END SUBROUTINE nullIFY_UNIVERSE


  SUBROUTINE APPEND_EMPTY_LAYOUT( L )   ! Appoints without cloning
    implicit none
    TYPE (MAD_UNIVERSE), TARGET:: L
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


  SUBROUTINE locate_in_universe(F,i_tot,i,j)
    IMPLICIT NONE
    integer i_tot,i,j
    integer k
    TYPE(FIBRE),pointer ::  F
    TYPE(layout),pointer ::  L


    call FIND_POS(f%mag%PARENT_FIBRE%parent_layout, f%mag%PARENT_FIBRE,j )

    call FIND_POS( f%mag%PARENT_FIBRE%parent_layout,i )

    i_tot=0

    L=>f%mag%PARENT_FIBRE%parent_layout%parent_universe%START

    do k=1,i-1
       i_tot= L%N+I_TOT
       L=>L%NEXT
    enddo

    I_TOT=I_ToT+ J



  END SUBROUTINE locate_in_universe



END MODULE S_FIBRE_BUNDLE
