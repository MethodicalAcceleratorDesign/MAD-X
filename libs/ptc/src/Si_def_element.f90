!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN


MODULE S_DEF_ELEMENT
 ! USE S_DEF_KIND
  !  USE USER_kind1
  !  USE USER_kind2
  USE sagan_WIGGLER

  IMPLICIT NONE
  public
  logical(lp),PARAMETER::BERZ=.TRUE.,ETIENNE=.NOT.BERZ
  logical(lp) :: USE_TPSAFIT=.TRUE.  ! USE GLOBAL ARRAY INSTEAD OF PERSONAL ARRAY
  logical(lp), target :: set_tpsafit=.false.
  logical(lp), target :: set_ELEMENT=.false.
  real(dp) , target :: scale_tpsafit=1.0_dp
  real(dp), target :: tpsafit(lnv) !   used for fitting with tpsa in conjunction with pol_block
  PRIVATE copy_el_elp,copy_elp_el,copy_el_el
  PRIVATE cop_el_elp,cop_elp_el,cop_el_el
  private ZERO_EL,ZERO_ELP
  !  PRIVATE MAGPSTATE,MAGSTATE
  PRIVATE SETFAMILYR,SETFAMILYP
  PRIVATE ADD_ANBNR,ADD_ANBNP,bL_0,EL_BL,ELp_BL,COPY_BL,UNARYP_BL
  PRIVATE ELp_POL,bLPOL_0
  PRIVATE work_0,work_r,ELp_WORK,EL_WORK,WORK_EL,WORK_ELP,BL_EL,BL_ELP,unaryw_w
  PRIVATE ZERO_ANBN,ZERO_ANBN_R,ZERO_ANBN_P
  private null_EL,null_ELp
  logical(lp), PRIVATE :: VERBOSE = .FALSE.
  logical(lp), PRIVATE :: GEN = .TRUE.
  logical(lp),TARGET :: ALWAYS_EXACTMIS=.TRUE.
  logical(lp),TARGET :: FEED_P0C=.FALSE.
  integer, TARGET :: np_pol=0
  !  logical(lp) :: isomorphism_MIS=.TRUE.  !Not needed anymore always should be true
  private put_aperture_el,put_aperture_elp
  integer :: mfpolbloc=0
  logical(lp),target :: recirculator_cheat=my_false
  PRIVATE TRACKR,TRACKP
  logical(lp), target :: restore_mag=my_true,restore_magp=my_true
  ! Old home for element and elementp, now in sh_def_kind







  INTERFACE EQUAL
     MODULE PROCEDURE copy_el_elp                              ! need upgrade
     MODULE PROCEDURE copy_elp_el                              ! need upgrade
     MODULE PROCEDURE copy_el_el                ! need upgrade
  end  INTERFACE

  INTERFACE COPY
     MODULE PROCEDURE cop_el_elp                              ! need upgrade
     MODULE PROCEDURE cop_elp_el                              ! need upgrade
     MODULE PROCEDURE cop_el_el                ! need upgrade
     MODULE PROCEDURE COPY_BL
  end  INTERFACE

  INTERFACE ADD
     MODULE PROCEDURE ADD_ANBNR
     MODULE PROCEDURE ADD_ANBNP
  end  INTERFACE

  INTERFACE ZERO_ANBN
     MODULE PROCEDURE ZERO_ANBN_R
     MODULE PROCEDURE ZERO_ANBN_P
  end  INTERFACE


  INTERFACE OPERATOR (+)
     MODULE PROCEDURE unaryP_BL
  END INTERFACE

  INTERFACE OPERATOR (+)
     MODULE PROCEDURE unaryw_w
  END INTERFACE


 INTERFACE print
     MODULE PROCEDURE print_work
  end  INTERFACE

  INTERFACE SETFAMILY
     MODULE PROCEDURE SETFAMILYR                              ! need upgrade
     MODULE PROCEDURE SETFAMILYP                              ! need upgrade
  end  INTERFACE

  INTERFACE null_ELEment
     MODULE PROCEDURE null_EL                               ! need upgrade
     MODULE PROCEDURE null_ELp                              ! need upgrade
  end  INTERFACE

  INTERFACE put_aperture
     MODULE PROCEDURE put_aperture_el                               ! need upgrade
     MODULE PROCEDURE put_aperture_elp                              ! need upgrade
  end  INTERFACE


  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE ZERO_EL                 ! NEED UPGRADE
     MODULE PROCEDURE ZERO_ELP                  ! NEED UPGRADE
     !     MODULE PROCEDURE MAGSTATE              ! need upgrade IF STATES EXPANDED
     !     MODULE PROCEDURE MAGPSTATE             ! need upgrade IF STATES EXPANDED
     ! Multipole block setting
     MODULE PROCEDURE BL_0
     MODULE PROCEDURE EL_BL
     MODULE PROCEDURE ELp_BL
     MODULE PROCEDURE BL_EL
     MODULE PROCEDURE BL_ELP
     ! polymorphism
     MODULE PROCEDURE bLPOL_0
     MODULE PROCEDURE ELp_POL
     ! energy/mass retrieving
     MODULE PROCEDURE work_0
     MODULE PROCEDURE work_r
     MODULE PROCEDURE ELp_WORK
     MODULE PROCEDURE EL_WORK
     MODULE PROCEDURE WORK_EL
     MODULE PROCEDURE WORK_ELP
  END INTERFACE


  INTERFACE TRACK
     !  INTERFACE TRACK
     MODULE PROCEDURE TRACKR
     MODULE PROCEDURE TRACKP
     !  END INTERFACE
     ! END old Sj_elements
  END INTERFACE





CONTAINS


  SUBROUTINE TRACKR(EL,X,K,MID)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(ELEMENT),INTENT(INOUT):: EL
    TYPE(WORM),OPTIONAL, INTENT(INOUT):: MID
    TYPE(INTERNAL_STATE) K

    if(associated(el%p%aperture)) then
     if(el%p%dir*el%p%aperture%pos==0.OR.el%p%dir*el%p%aperture%pos==-1) call CHECK_APERTURE(EL%p%aperture,X)
    endif
    !    if(other_program) then
    !       call track_R(x)
    !       return
    !    endif
    SELECT CASE(EL%KIND)
    CASE(KIND0)
       IF(PRESENT(MID)) CALL XMID(MID,X,0)
       IF(PRESENT(MID)) CALL XMID(MID,X,1)   ! ADDED FOR NST=1 IN MARKER FOR THIN_LAYOUT SURVEY
    case(KIND1)
       CALL TRACK(EL%D0,X,K,MID)
    case(KIND2)
       CALL TRACK(EL%K2,X,k,MID)
    case(KIND3)
       CALL TRACK(EL%K3,X,k,MID)
    case(KIND4)
       CALL TRACK(EL%C4,X,k,MID)
    case(KIND5)
       CALL TRACK(EL%S5,X,k,MID)
    case(KIND6)
       CALL TRACK(EL%T6,X,k,MID)
    case(KIND7)
       CALL TRACK(EL%T7,X,k,MID)
    case(KIND8)
       CALL TRACK(EL%S8,X,k,MID)
    case(KIND9)
       CALL TRACK(EL%S9,X,k,MID)
    case(KIND10)
       CALL TRACK(EL%TP10,X,k,MID)
    CASE(KIND11:KIND14)
       call TRACK(EL%MON14,X,k,MID)
    CASE(KIND15)
       call TRACK(EL%SEP15,X,k,MID)
    CASE(KIND16,KIND20)
       call TRACK(EL%K16,X,k,MID)
    CASE(KIND18)
       call TRACK(EL%RCOL18,X,k,MID)
    CASE(KIND19)
       call TRACK(EL%ECOL19,X,k,MID)
    CASE(KIND21)
       call TRACK(EL%CAV21,X,k,MID)
    CASE(KIND22)
       call TRACK(EL%HE22,X,k,MID)
    case(KINDWIGGLER)
       call TRACK(EL%WI,X,k,MID)
    case(KINDPA)
       call TRACK(EL%PA,X,k,MID)
    case(kindsuperdrift)
       call TRACK(EL%SDR,X,k,MID)
    case(KINDABELL)
       call TRACK(EL%AB,X,k,MID)

    case default
 
       write(6,'(1x,i4,a21)') el%kind," not supported TRACKR"
       ! call !write_e(0)
    END SELECT
    if(associated(el%p%aperture)) then
     if(el%p%dir*el%p%aperture%pos==0.OR.el%p%dir*el%p%aperture%pos==1)  call CHECK_APERTURE(EL%p%aperture,X)
    endif
  END SUBROUTINE TRACKR

  SUBROUTINE TRACKP(EL,X,K)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(ELEMENTP),INTENT(INOUT):: EL
    !    TYPE(WORM_8),OPTIONAL, INTENT(INOUT):: MID
    TYPE(INTERNAL_STATE) K

    if(associated(el%p%aperture)) then
     if(el%p%dir*el%p%aperture%pos==0.OR.el%p%dir*el%p%aperture%pos==-1)  call CHECK_APERTURE(EL%p%aperture,X)
    endif
    SELECT CASE(EL%KIND)
    CASE(KIND0)
       !       IF(PRESENT(MID)) CALL XMID(MID,X,0)
    case(KIND1)
       CALL TRACK(EL%D0,X,K)
    case(KIND2)
       CALL TRACK(EL%K2,X,k)
    case(KIND3)
       CALL TRACK(EL%K3,X,k)
    case(KIND4)
       CALL TRACK(EL%C4,X,k)
    case(KIND5)
       CALL TRACK(EL%S5,X,k)
    case(KIND6)
       CALL TRACK(EL%T6,X,k)
    case(KIND7)
       CALL TRACK(EL%T7,X,k)
    case(KIND8)
       CALL TRACK(EL%S8,X,k)
    case(KIND9)
       CALL TRACK(EL%S9,X,k)
    case(KIND10)
       CALL TRACK(EL%TP10,X,k)
    CASE(KIND11:KIND14)
       call TRACK(EL%MON14,X,k)
    CASE(KIND15)
       call TRACK(EL%SEP15,X,k)
    CASE(KIND16,KIND20)
       call TRACK(EL%K16,X,k)
    CASE(KIND18)
       call TRACK(EL%RCOL18,X,k)
    CASE(KIND19)
       call TRACK(EL%ECOL19,X,k)
    CASE(KIND21)
       call TRACK(EL%CAV21,X,k)
    CASE(KIND22)
       call TRACK(EL%HE22,X,k)
    case(KINDWIGGLER)
       call TRACK(EL%WI,X,k)
    case(KINDPA)
       call TRACK(EL%PA,X,k)
    case(kindsuperdrift)
       call TRACK(EL%SDR,X,k)
    case(KINDABELL)
       call TRACK(EL%AB,X,k)
    case default
       write(6,'(1x,i4,a21)') el%kind," not supported TRACKP"
       ! call !write_e(0)
    END SELECT
    if(associated(el%p%aperture)) then
     if(el%p%dir*el%p%aperture%pos==0.OR.el%p%dir*el%p%aperture%pos==1)  call CHECK_APERTURE(EL%p%aperture,X)
    endif
  END SUBROUTINE TRACKP

  !  SUBROUTINE TRACK_R(X)
  !    IMPLICIT NONE
  !    REAL(DP) X(6),x6,xp,yp
  !    INTEGER icharef
  !    COMMON/ptc/ icharef
  !
  !
  !    if(j_global==1) return  ! skipping OBJECT OF ZGOUBI = TRACKING COMMAND INTERNAL TO ZGOUBI
  !    icharef=0

  !
  !    x(1)=x(1)*c_100
  !    x(3)=x(3)*c_100
  !    x6=x(6)*c_100
  !
  !    xp=x(2)/root((one+x(5))**2-x(2)**2-x(4)**2)
  !    yp=x(4)/root((one+x(5))**2-x(2)**2-x(4)**2)
  !    x(2)=atan(xp)*c_1d3
  !    x(4)=atan(yp/root(one+xp**2))*c_1d3
  !
  !    x(6)=x(5)
  !    x(5)=x6
  !
  !    !call track_z(x,j_global,j_global)
  !
  !    x6=x(5)/c_100
  !    x(5)=x(6)
  !    x(6)=x6
  !
  !    x(1)=x(1)/c_100
  !    x(3)=x(3)/c_100
  !    xp=tan(x(2)/c_1d3)
  !    yp=tan(x(4)/c_1d3)*root(one+xp**2)

  !    x(2)=(one+x(5))*xp/root(one+xp**2+yp**2)
  !    x(4)=(one+x(5))*yp/root(one+xp**2+yp**2)

  !   icharef=1

  !  END SUBROUTINE TRACK_R

  !  SUBROUTINE TRACK_P(X)
  !    IMPLICIT NONE
  !    TYPE(REAL_8) X(6)

  ! track_zp is a fortran external routine using numerical differentiation
  !call track_zp(x,j_global,j_global)
  !    WRITE(6,*) " NOT SUPPORTED "
  !    STOP 111
  !  END SUBROUTINE TRACK_P




  SUBROUTINE  work_0(S2,S1)
    implicit none
    type (work),INTENT(inOUT):: S2
    INTEGER,INTENT(IN):: S1

    S2%BETA0=1.0_dp
    S2%energy=0.0_dp
    S2%kinetic=0.0_dp
    S2%p0c=0.0_dp
    S2%brho=0.0_dp
    S2%mass=0.0_dp
    S2%gamma0I=0.0_dp
    S2%gambet=0.0_dp
    if(s1/=0) then
       S2%rescale=.true.
       s2%power=s1
    else
       S2%rescale=.false.
       s2%power=0
    endif
  END SUBROUTINE work_0



  SUBROUTINE  work_r(S2,S1)
    implicit none
    type (work),INTENT(inOUT):: S2
    real(dp),INTENT(IN):: S1

    !    S2%energy=-(S2%energy+s1)
    !  VERBOSE = .FALSE.
    IF(FEED_P0C) THEN
       call find_energy(s2,P0C=S1+S2%P0C)
    ELSE
       call find_energy(s2,ENERGY=S1+S2%energy)
    ENDIF
    !  VERBOSE = .TRUE.
  END SUBROUTINE work_r

  SUBROUTINE  print_work(S2,mf)
    implicit none
    type (work),INTENT(inOUT):: S2
    integer,INTENT(IN):: mf
    
    write(mf,*) "Beta0 = ",s2%beta0
    write(mf,*) "Mass  = ",s2%mass
    write(mf,*) "Energy = ",s2%Energy
    write(mf,*) "Kinetic Energy = ",s2%kinetic
    write(mf,*) "p0c = ",s2%p0c
    write(mf,*) "gamma  = ",1.d0/s2%gamma0i
    write(mf,*) "p0c = ",s2%p0c
     write(mf,*) "brho = ",s2%brho
    write(mf,*) "rescale and power = ",s2%rescale,s2%power


  END SUBROUTINE print_work

  function  unaryw_w(S1)
    implicit none
    type (WORK),INTENT(IN):: S1
    TYPE(WORK) unaryw_w
    unaryw_w=s1
    unaryw_w%rescale=.false.

  end   function  unaryw_w

  SUBROUTINE  ELp_WORK(S2,S1)
    implicit none
    type (WORK),INTENT(IN):: S1
    TYPE(ELEMENTP),INTENT(inOUT):: S2
    integer i

    if(s1%rescale) then
       if(s2%p%nmul/=0) then   ! doing for crab also
          do i=1,s2%P%nmul
             s2%bn(i)=s2%bn(i)*(S2%P%P0C/S1%P0C)**S1%power
             s2%an(i)=s2%an(i)*(S2%P%P0C/S1%P0C)**S1%power
          enddo
          CALL ADD(s2,1,1,0.0_dp)
       endif
       if(associated(s2%B_sol))  s2%B_sol=s2%B_sol*(S2%P%P0C/S1%P0C)**S1%power

       !       if(s2%kind==kinduser1) call scale_user1(s2%u1,S2%P%P0C,S1%P0C,S1%power)
       !       if(s2%kind==kinduser2) call scale_user2(s2%u2,S2%P%P0C,S1%P0C,S1%power)
       if(s2%kind==KINDwiggler) call scale_sagan(s2%wi,S2%P%P0C,S1%P0C,S1%power)

    else   ! ramping
       !       S2%P%BETA0=S1%BETA0
       !       S2%P%GAMMA0I=S1%GAMMA0I
       !       S2%P%GAMBET=S1%GAMBET
       S2%P%P0C=S1%P0C
    endif

  END SUBROUTINE ELp_WORK

  SUBROUTINE  EL_WORK(S2,S1)
    implicit none
    type (WORK),INTENT(IN):: S1
    TYPE(ELEMENT),INTENT(inOUT):: S2
    integer i

    if(s1%rescale) then
       if(s2%p%nmul/=0) then
          do i=1,s2%P%nmul
             s2%bn(i)=s2%bn(i)*(S2%P%P0C/S1%P0C)**S1%power
             s2%an(i)=s2%an(i)*(S2%P%P0C/S1%P0C)**S1%power
          enddo
          CALL ADD(s2,1,1,0.0_dp)
       endif
       if(associated(s2%B_sol))  s2%B_sol=s2%B_sol*(S2%P%P0C/S1%P0C)**S1%power
       !       if(s2%kind==kinduser1) call scale_user1(s2%u1,S2%P%P0C,S1%P0C,S1%power)
       !       if(s2%kind==kinduser2) call scale_user2(s2%u2,S2%P%P0C,S1%P0C,S1%power)
       if(s2%kind==KINDwiggler) call scale_sagan(s2%wi,S2%P%P0C,S1%P0C,S1%power)
    else  ! ramping


       !       S2%P%BETA0=S1%BETA0
       !       S2%P%GAMMA0I=S1%GAMMA0I
       !       S2%P%GAMBET=S1%GAMBET
       S2%P%P0C=S1%P0C
    endif


  END SUBROUTINE EL_WORK


  SUBROUTINE  WORK_EL(S1,S2)
    implicit none
    type (WORK),INTENT(inOUT):: S1
    TYPE(ELEMENT),INTENT(IN):: S2

    S1=S1%POWER

    !    S1%P0C=-S2%P%P0C
    !  VERBOSE = .FALSE.
    call find_energy(s1,P0C=S2%P%P0C)
    !  VERBOSE = .TRUE.

  END SUBROUTINE WORK_EL

  SUBROUTINE  WORK_ELp(S1,S2)
    implicit none
    type (WORK),INTENT(inOUT):: S1
    TYPE(ELEMENTP),INTENT(IN):: S2

    S1=S1%POWER

    !    S1%P0C=-S2%P%P0C
    !  VERBOSE = .FALSE.
    call find_energy(s1,P0C=S2%P%P0C)
    !  VERBOSE = .TRUE.

  END SUBROUTINE WORK_ELp


  integer function mod_n(i,j)
    implicit none
    integer, intent(in) :: i,j
    integer k
 
    k=i
    if(i<1) then
       do while(k<1)
          k=k+j
       enddo
    endif
    mod_n=mod(k,j)
    if(mod_n==0) mod_n=j
  end function  mod_n

  SUBROUTINE  bL_0(S2,S1)
    implicit none
    type (MUL_BLOCK),INTENT(OUT):: S2
    INTEGER,INTENT(IN):: S1
    INTEGER I

    IF(S1>=0.OR.S1<=nmax) THEN
       do i = 1,nmax
          s2%aN(i)=0.0_dp
          s2%bN(i)=0.0_dp
       enddo
       s2%natural=1
       s2%nmul=S1
       s2%ADD=0
    ELSEIF(S1>NMAX) THEN
 
       write(6,'(A38,1X,I4)') " NMAX NOT BIG ENOUGH: PLEASE INCREASE ",NMAX
 
    ELSE
      stop 135
    ENDIF

  END SUBROUTINE bL_0

  SUBROUTINE  bLPOL_0(S2,S1)
    implicit none
    type (POL_BLOCK),INTENT(OUT):: S2
    INTEGER,INTENT(IN):: S1
    INTEGER I

    !    IF(S1>=0.and.S1<=nmax) THEN
    do i = 1,nmax
       s2%SAN(i)=1.0_dp
       s2%SBN(i)=1.0_dp
       s2%IaN(i)=0
       s2%IbN(i)=0
    enddo
    !    S2%user1=0
    !    S2%user2=0
    S2%SAGAN=0
    S2%SVOLT=1.0_dp
    S2%SFREQ=1.0_dp
    S2%SPHAS=1.0_dp
    S2%SB_SOL=1.0_dp
    S2%IVOLT=0
    S2%IFREQ=0
    S2%IPHAS=0
    S2%IB_SOL=0
    s2%npara=S1
    s2%g=0
    s2%np=0
    s2%nb=0
    !     s2%NMUL=0
    s2%NAME=' '
    s2%N_NAME=0
    s2%VORNAME=' '
    !    s2%CHECK_NMUL=.TRUE.
    nullify(s2%tpsafit);nullify(s2%set_tpsafit);
    nullify(s2%set_ELEMENT);

    IF(USE_TPSAFIT) then
       s2%tpsafit=>tpsafit
       s2%set_tpsafit=>set_tpsafit
       s2%set_ELEMENT=>set_ELEMENT
    endif

    if(s1>0) then
       c_%npara_fpp=0   ! backwards compatible
    endif

  END SUBROUTINE bLPOL_0

  SUBROUTINE  BL_EL(S1,S2)
    implicit none
    type (MUL_BLOCK),INTENT(out):: S1
    TYPE(ELEMENT),INTENT(IN):: S2
    INTEGER I

    IF(S2%P%NMUL>NMAX) THEN
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(1((1X,A72)))'
       write(*,'(A21,1X,I4,1X,I4)')  " NMAX NOT BIG ENOUGH ", S2%P%NMUL,NMAX
       ! call !write_e(456)
    ENDIF
    S1=S2%P%NMUL

    DO I=1,S2%P%NMUL
       s1%AN(I)=s2%AN(I)
       s1%BN(I)=s2%BN(I)
    ENDDO

  END SUBROUTINE BL_EL

  SUBROUTINE  BL_ELP(S1,S2)
    implicit none
    type (MUL_BLOCK),INTENT(out):: S1
    TYPE(ELEMENTP),INTENT(IN):: S2
    INTEGER I

    IF(S2%P%NMUL>NMAX) THEN
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(1((1X,A72)))'
       write(*,'(A21,1X,I4,1X,I4)')  " NMAX NOT BIG ENOUGH ", S2%P%NMUL,NMAX
       ! call !write_e(456)
    ENDIF
    S1=S2%P%NMUL

    DO I=1,S2%P%NMUL
       s1%AN(I)=s2%AN(I)
       s1%BN(I)=s2%BN(I)
    ENDDO

  END SUBROUTINE BL_ELP

  SUBROUTINE  EL_BL(S2,S1)
    implicit none
    type (MUL_BLOCK),INTENT(IN):: S1
    TYPE(ELEMENT),INTENT(inOUT):: S2
    INTEGER I

    IF(S2%P%NMUL>NMAX) THEN
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(1((1X,A72)))'
       write(*,'(A21,1X,I4,1X,I4)')  " NMAX NOT BIG ENOUGH ", S2%P%NMUL,NMAX
       ! call !write_e(456)
    ENDIF
    IF(s1%nmul>s2%P%nmul) CALL ADD(s2,s1%nmul,1,0.0_dp)

    DO I=1,S2%P%NMUL
       s2%AN(I)=S1%ADD*s2%AN(I)+s1%AN(I)
       s2%BN(I)=S1%ADD*s2%BN(I)+s1%BN(I)
    ENDDO
    CALL ADD(s2,1,1,0.0_dp)

  END SUBROUTINE EL_BL

  SUBROUTINE  ELp_BL(S2,S1)
    implicit none
    type (MUL_BLOCK),INTENT(IN):: S1
    TYPE(ELEMENTP),INTENT(inOUT):: S2
    INTEGER I

    IF(S2%P%NMUL>NMAX) THEN
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(1((1X,A72)))'
       write(*,'(A21,1X,I4,1X,I4)')  " NMAX NOT BIG ENOUGH ", S2%P%NMUL,NMAX
       ! call !write_e(456)
    ENDIF

    IF(s1%nmul>s2%P%nmul) CALL ADD(s2,s1%nmul,1,0.0_dp)
    DO I=1,S2%P%NMUL
       s2%AN(I)=S1%ADD*s2%AN(I)+s1%AN(I)
       s2%BN(I)=S1%ADD*s2%BN(I)+s1%BN(I)
    ENDDO
    CALL ADD(s2,1,1,0.0_dp)


  END SUBROUTINE ELp_BL

  SUBROUTINE  ELp_POL(S2,S1)
    implicit none
    type (POL_BLOCK),INTENT(IN):: S1
    TYPE(ELEMENTP),INTENT(inOUT):: S2
    logical(lp) DOIT                    !,checkname
    CHARACTER(nlp) S1NAME
    CHARACTER(vp)    S1VORNAME


    IF(S2%P%NMUL>NMAX) THEN
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(1((1X,A72)))'
       write(*,'(A21,1X,I4,1X,I4)')  " NMAX NOT BIG ENOUGH ", S2%P%NMUL,NMAX
       ! call !write_e(456)
    ENDIF

    S1NAME=S1%name
    S1VORNAME=S1%VORname
    CALL CONTEXT(S1name)
    CALL CONTEXT(S1vorname)
    CALL CONTEXT(S2%name)
    CALL CONTEXT(S2%vorname)

    DOIT=.TRUE.
    IF(S1NAME/=' ') THEN
       if(s1%n_name==0) then
          IF(S1NAME/=S2%NAME) DOIT=.FALSE.
       else
          IF(S1NAME(1:s1%n_name)/=S2%NAME(1:s1%n_name)) DOIT=.FALSE.
       endif
    ENDIF


    IF(S1VORNAME/=' ') THEN
       IF(S1VORNAME/=S2%VORNAME.or.S1NAME/=S2%NAME) DOIT=.FALSE.
    ENDIF


    IF(DOIT) THEN
       IF(.not.S1%SET_TPSAFIT.AND.(.NOT.SET_ELEMENT)) THEN
          if(s2%knob) then
             write(6,'(A45,A16)')" BE CAREFUL USING A POL_BLOCK ON SAME MAGNET ",S2%NAME
          ENDIF
       endif
       s2%knob=.TRUE.
       call ELp_POL_force(S2,S1)
    ENDIF


  END SUBROUTINE ELp_POL

  SUBROUTINE  ELp_POL_force(S2,S1)
    implicit none
    type (POL_BLOCK),INTENT(IN):: S1
    TYPE(ELEMENTP),INTENT(inOUT):: S2
    INTEGER I,S1NMUL
    logical(lp) DOIT,DONEIT                    !,checkname

    IF(S2%P%NMUL>NMAX) THEN
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(1((1X,A72)))'
       write(*,'(A21,1X,I4,1X,I4)')  " NMAX NOT BIG ENOUGH ", S2%P%NMUL,NMAX
       ! call !write_e(456)
    ENDIF


    DOIT=.TRUE.



    s2%knob=.TRUE.

    !       IF(S1%NPARA>=4.AND.S1%NPARA<=6) THEN
    DONEIT=.FALSE.

    !        IF(S1%CHECK_NMUL) THEN
    S1NMUL=0
    DO I=NMAX,1,-1
       IF(s1%IAN(I)/=0.OR.s1%IBN(I)/=0)  THEN
          S1NMUL=I
          if(s1%IAN(I)>c_%np_pol) c_%np_pol=s1%IAN(I)
          if(s1%IBN(I)>c_%np_pol) c_%np_pol=s1%IBN(I)
          GOTO 100
       ENDIF
    ENDDO
100 CONTINUE
    !          CALL SET_FALSE(S1%CHECK_NMUL)
    !        ENDIF

    IF(S1NMUL>S2%P%NMUL) then
       CALL ADD(S2,S1NMUL,1,0.0_dp)  !etienne
    endif
    DO I=1,S1NMUL
       IF(S1%IAN(I)>0) THEN
          s2%AN(I)%I=S1%IAN(I)+S1%NPARA
          s2%AN(I)%S=S1%SAN(I)
          s2%AN(I)%KIND=3
!          s2%AN(I)%g=S1%g
!          s2%AN(I)%nb=S1%nb
          DONEIT=.TRUE.
          IF(S1%SET_TPSAFIT) THEN
             s2%aN(I)%R=s2%aN(I)%R+scale_tpsafit*s2%AN(I)%S*s1%TPSAFIT(S1%IAN(I))
          ENDIF
          IF(S1%SET_ELEMENT) THEN
             s2%PARENT_FIBRE%MAG%aN(I)=s2%aN(I)%R
          ENDIF
       ENDIF
       IF(S1%IBN(I)>0) THEN
          s2%BN(I)%I=S1%IBN(I)+S1%NPARA
          s2%BN(I)%S=S1%SBN(I)
          s2%BN(I)%KIND=3
!          s2%BN(I)%g=S1%g
!          s2%BN(I)%nb=S1%nb
          DONEIT=.TRUE.
          IF(S1%SET_TPSAFIT) THEN
             s2%BN(I)%R=s2%BN(I)%R+scale_tpsafit*s2%BN(I)%S*s1%TPSAFIT(S1%IBN(I))
          ENDIF
          IF(S1%SET_ELEMENT) THEN
             s2%PARENT_FIBRE%MAG%BN(I)=s2%BN(I)%R
          ENDIF
       ENDIF
    ENDDO
    IF(DONEIT.AND.(S1%SET_TPSAFIT.OR.S1%SET_ELEMENT)) THEN
       CALL ADD(S2,1,1,0.0_dp)     !etienne
    ENDIF
    IF(S2%KIND==KIND4) THEN    ! CAVITY
       DONEIT=.FALSE.                     ! NOT USED HERE
       IF(S1%IVOLT>0) THEN
          s2%VOLT%I=S1%IVOLT+S1%NPARA
          s2%VOLT%S=S1%SVOLT
          s2%VOLT%KIND=3
!          s2%VOLT%g=S1%g
!          s2%VOLT%nb=S1%nb
          DONEIT=.TRUE.
          if(S1%IVOLT>c_%np_pol) c_%np_pol=S1%IVOLT
          IF(S1%SET_TPSAFIT) THEN
             s2%VOLT%R=s2%VOLT%R+scale_tpsafit*s2%VOLT%S*s1%TPSAFIT(S1%IVOLT)
          ENDIF
          IF(S1%SET_ELEMENT) THEN
             s2%PARENT_FIBRE%MAG%VOLT=s2%VOLT%R
          ENDIF
       ENDIF
       IF(S1%IFREQ>0) THEN
          s2%FREQ%I=S1%IFREQ+S1%NPARA
          s2%FREQ%S=S1%SFREQ
!          s2%FREQ%g=S1%g
!          s2%FREQ%nb=S1%nb
          s2%FREQ%KIND=3
          if(S1%IFREQ>c_%np_pol) c_%np_pol=S1%IFREQ
          IF(S1%SET_TPSAFIT) THEN
             s2%FREQ%R=s2%FREQ%R+scale_tpsafit*s2%FREQ%S*s1%TPSAFIT(S1%IFREQ)
          ENDIF
          IF(S1%SET_ELEMENT) THEN
             s2%PARENT_FIBRE%MAG%FREQ=s2%FREQ%R
          ENDIF
          DONEIT=.TRUE.
       ENDIF
       IF(S1%IPHAS>0) THEN
          s2%PHAS%I=S1%IPHAS+S1%NPARA
          s2%PHAS%S=S1%SPHAS
          s2%PHAS%KIND=3
 !         s2%PHAS%g=S1%g
 !         s2%PHAS%nb=S1%nb
          DONEIT=.TRUE.
          if(S1%IPHAS>c_%np_pol) c_%np_pol=S1%IPHAS
          IF(S1%SET_TPSAFIT) THEN
             s2%PHAS%R=s2%PHAS%R+scale_tpsafit*s2%PHAS%S*s1%TPSAFIT(S1%IPHAS)
          ENDIF
          IF(S1%SET_ELEMENT) THEN
             s2%PARENT_FIBRE%MAG%PHAS=s2%PHAS%R
          ENDIF
       ENDIF
    ENDIF
    IF(S2%KIND==KIND21) THEN    ! CAVITY
       DONEIT=.FALSE.                     ! NOT USED HERE
       IF(S1%IVOLT>0) THEN
          s2%VOLT%I=S1%IVOLT+S1%NPARA
          s2%VOLT%S=S1%SVOLT
!          s2%PHAS%g=S1%g
!          s2%PHAS%nb=S1%nb
          s2%VOLT%KIND=3
          if(S1%IVOLT>c_%np_pol) c_%np_pol=S1%IVOLT
          DONEIT=.TRUE.
          IF(S1%SET_TPSAFIT) THEN
             s2%VOLT%R=s2%VOLT%R+scale_tpsafit*s2%VOLT%S*s1%TPSAFIT(S1%IVOLT)
          ENDIF
          IF(S1%SET_ELEMENT) THEN
             s2%PARENT_FIBRE%MAG%VOLT=s2%VOLT%R
          ENDIF
       ENDIF
       IF(S1%IFREQ>0) THEN
          s2%FREQ%I=S1%IFREQ+S1%NPARA
          s2%FREQ%S=S1%SFREQ
!          s2%FREQ%g=S1%g
!          s2%FREQ%nb=S1%nb
          s2%FREQ%KIND=3
          if(S1%IFREQ>c_%np_pol) c_%np_pol=S1%IFREQ
          IF(S1%SET_TPSAFIT) THEN
             s2%FREQ%R=s2%FREQ%R+scale_tpsafit*s2%FREQ%S*s1%TPSAFIT(S1%IFREQ)
          ENDIF
          IF(S1%SET_ELEMENT) THEN
             s2%PARENT_FIBRE%MAG%FREQ=s2%FREQ%R
          ENDIF
          DONEIT=.TRUE.
       ENDIF
       IF(S1%IPHAS>0) THEN
          s2%PHAS%I=S1%IPHAS+S1%NPARA
          s2%PHAS%S=S1%SPHAS
  !        s2%PHAS%g=S1%g
  !        s2%PHAS%nb=S1%nb
          s2%PHAS%KIND=3
          if(S1%IPHAS>c_%np_pol) c_%np_pol=S1%IPHAS
          DONEIT=.TRUE.
          IF(S1%SET_TPSAFIT) THEN
             s2%PHAS%R=s2%PHAS%R+scale_tpsafit*s2%PHAS%S*s1%TPSAFIT(S1%IPHAS)
          ENDIF
          IF(S1%SET_ELEMENT) THEN
             s2%PARENT_FIBRE%MAG%PHAS=s2%PHAS%R
          ENDIF
       ENDIF
    ENDIF
    IF(S2%KIND==KIND5) THEN    ! SOLENOID
       DONEIT=.FALSE.
       IF(S1%IB_SOL>0) THEN
          s2%B_SOL%I=S1%IB_SOL+S1%NPARA
          s2%B_SOL%S=S1%SB_SOL
    !      s2%B_SOL%g=S1%g
    !      s2%B_SOL%nb=S1%nb
          s2%B_SOL%KIND=3
          DONEIT=.TRUE.
          if(S1%IB_SOL>c_%np_pol) c_%np_pol=S1%IB_SOL
          IF(S1%SET_TPSAFIT) THEN
             s2%B_SOL%R=s2%B_SOL%R+scale_tpsafit*s2%B_SOL%S*s1%TPSAFIT(S1%IB_SOL)
          ENDIF
          IF(S1%SET_ELEMENT) THEN
             s2%PARENT_FIBRE%MAG%B_SOL=s2%B_SOL%R
          ENDIF
       ENDIF
    ENDIF
    !    IF(S2%KIND==kinduser1) THEN    ! new element
    !       DONEIT=.FALSE.                     ! NOT USED HERE
    !       call ELp_POL_user1(S2%u1,S1,DONEIT)
    !    ENDIF
    !    IF(S2%KIND==kinduser2) THEN    ! new element
    !       DONEIT=.FALSE.                     ! NOT USED HERE
    !       call ELp_POL_user2(S2%u2,S1,DONEIT)
    !    ENDIF
    IF(S2%KIND==KINDWIGGLER) THEN    ! new element
       DONEIT=.FALSE.                     ! NOT USED HERE
       call ELp_POL_SAGAN(S2%WI,S2%PARENT_FIBRE%MAG%WI,S1,DONEIT)
    ENDIF



  END SUBROUTINE ELp_POL_force

  SUBROUTINE  ELp_POL_print(S2)
    implicit none
    TYPE(ELEMENTP),INTENT(inOUT):: S2
    INTEGER I
    type(work) w



    !          CALL SET_FALSE(S1%CHECK_NMUL)
    !        ENDIF

    DO I=1,S2%P%NMUL
       IF(s2%AN(I)%KIND==3) THEN
          w=s2
          write(mfpolbloc,'(a16,a8,1x,i4,2(1x,e18.8))') s2%name, ' MAD AN ',i,s2%aN(I)%R*MADFAC(I),s2%aN(I)%R*w%brho*MADFAC(I)
       ENDIF
       IF(s2%bN(I)%KIND==3) THEN
          w=s2
          write(mfpolbloc,'(a16,a8,1x,i4,2(1x,e18.8))') s2%name, ' MAD BN ',i,s2%BN(I)%R*MADFAC(I),s2%BN(I)%R*w%brho*MADFAC(I)
       endif
    ENDDO
    IF(S2%KIND==KIND4.or.S2%KIND==KIND21) THEN    ! CAVITY
       IF(s2%VOLT%KIND==3) THEN
          write(mfpolbloc,*) s2%name, ' VOLT ',s2%VOLT%R
       ENDIF
       IF(s2%FREQ%KIND==3) THEN
          write(mfpolbloc,*) s2%name, ' FREQ ',s2%FREQ%R
       ENDIF
       IF(s2%PHAS%KIND==3) THEN
          write(mfpolbloc,*) s2%name, ' PHAS ',s2%PHAS%R
       ENDIF
    ENDIF
    IF(S2%KIND==KIND5) THEN    ! SOLENOID
       IF(s2%B_SOL%KIND==3) THEN
          write(mfpolbloc,*) s2%name, ' B_SOL ',s2%B_SOL%R
       ENDIF
    ENDIF

    !    IF(S2%KIND==KINDWIGGLER) THEN    ! new element
    !       DONEIT=.FALSE.                     ! NOT USED HERE
    !       call ELp_POL_SAGAN(S2%WI,S1,DONEIT)
    !    ENDIF



  END SUBROUTINE ELp_POL_print





  SUBROUTINE  COPY_BL(S1,S2)
    implicit none
    type (MUL_BLOCK),INTENT(IN):: S1
    TYPE(MUL_BLOCK),INTENT(OUT):: S2
    INTEGER I

    DO I=1,NMAX
       s2%AN(I)=s1%AN(I)
       s2%BN(I)=S1%BN(I)
    ENDDO

    S2%NMUL     =S1%NMUL
    S2%ADD      =S1%ADD
    S2%NATURAL  =S1%NATURAL

  END SUBROUTINE COPY_BL


  FUNCTION  UNARYP_BL(S1)
    implicit none
    type (MUL_BLOCK),INTENT(IN):: S1
    type (MUL_BLOCK) UNARYP_BL

    CALL COPY(S1,UNARYP_BL)
    UNARYP_BL%ADD=1

  END FUNCTION UNARYP_BL




  !  SUBROUTINE SETFAMILYR(EL,T,t_ax,t_ay,NTOT,ntot_rad,NTOT_REV,ntot_rad_REV,ND2)
  SUBROUTINE SETFAMILYR(EL,T)  !,angc,xc,dc,h)  !,NTOT,ntot_rad,NTOT_REV,ntot_rad_REV,ND2)
    IMPLICIT NONE
    TYPE(ELEMENT), INTENT(INOUT) ::EL
    !    INTEGER,OPTIONAL :: NTOT,ntot_rad,NTOT_REV,ntot_rad_REV,ND2
    type(tree_element),OPTIONAL :: T(:) !,t_ax(:),t_ay(:)
 !   real(dp), optional :: angc,xc,dc,h
   ! EL%P%permfringe=>EL%permfringe
    SELECT CASE(EL%KIND)
    CASE(KIND1)
       if(.not.ASSOCIATED(EL%D0))ALLOCATE(EL%D0)
       EL%D0%P=>EL%P
       EL%D0%L=>EL%L
    CASE(KIND2)
       IF(EL%P%EXACT) THEN
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,A72,/,1X,A72))'
          !w_p%c(1)=" ERROR IN SETFAMILYR "
          write(*,'(A43,1x,I4)') " 'NOT EXACT' OPTION NOT SUPPORTED FOR KIND ", EL%KIND
          ! call !write_e(222)
       ENDIF
       if(.not.ASSOCIATED(EL%K2)) THEN
          ALLOCATE(EL%K2)
          EL%K2=0
       ELSE
          EL%K2=-1
          EL%K2=0
       ENDIF
       !       if(.not.ASSOCIATED(EL%K2))ALLOCATE(EL%K2)
       EL%K2%P=>EL%P
       EL%K2%L=>EL%L
       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       EL%K2%AN=>EL%AN
       EL%K2%BN=>EL%BN
       EL%K2%FINT=>EL%FINT
       EL%K2%HGAP=>EL%HGAP
       EL%K2%H1=>EL%H1
       EL%K2%H2=>EL%H2
       EL%K2%VA=>EL%VA
       EL%K2%VS=>EL%VS
       NULLIFY(EL%K2%F);ALLOCATE(EL%K2%F);EL%K2%F=1;
    CASE(KIND3)
       if(.not.ASSOCIATED(EL%K3)) THEN
          ALLOCATE(EL%K3)
          el%K3=0
       ELSE
          el%K3=-1
          el%K3=0
       ENDIF
       EL%K3%P=>EL%P
       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       EL%K3%AN=>EL%AN
       EL%K3%BN=>EL%BN
       ALLOCATE(EL%K3%ls);EL%K3%ls=1.0_dp
       ALLOCATE(EL%K3%hf);EL%K3%hf=0
       ALLOCATE(EL%K3%vf);EL%K3%vf=0
       ALLOCATE(EL%K3%thin_h_foc);EL%K3%thin_h_foc=0
       ALLOCATE(EL%K3%thin_v_foc);EL%K3%thin_v_foc=0
       ALLOCATE(EL%K3%thin_h_angle);EL%K3%thin_h_angle=0
       ALLOCATE(EL%K3%thin_v_angle);EL%K3%thin_v_angle=0
       ALLOCATE(EL%K3%patch);EL%K3%patch=my_false
       EL%K3%B_SOL=>EL%B_SOL
       NULLIFY(EL%k3%dx);ALLOCATE(EL%k3%dx);EL%k3%dx=0.d0;
       NULLIFY(EL%k3%dy);ALLOCATE(EL%k3%dy);EL%k3%dy=0.d0;
       NULLIFY(EL%k3%pitch_x);ALLOCATE(EL%k3%pitch_x);EL%k3%pitch_x=0.d0;
       NULLIFY(EL%k3%pitch_y);ALLOCATE(EL%k3%pitch_y);EL%k3%pitch_y=0.d0;
    CASE(kindsuperdrift)
       if(.not.ASSOCIATED(EL%sdr)) THEN
          ALLOCATE(EL%sdr)
          EL%sdr=0
       ELSE
          EL%sdr=-1
          EL%sdr=0
       ENDIF
       EL%SDR%P=>EL%P
       EL%SDR%L=>EL%L

       ALLOCATE(EL%SDR%D(3));EL%SDR%D=0.0_dp;
       ALLOCATE(EL%SDR%ANG(3));EL%SDR%ANG=0.0_dp;
       ALLOCATE(EL%SDR%a_x1);EL%SDR%a_x1=1.0_dp;
       ALLOCATE(EL%SDR%a_x2);EL%SDR%a_x2=1.0_dp;

    CASE(kind4)
       if(.not.ASSOCIATED(EL%C4)) THEN
          ALLOCATE(EL%C4)
          el%C4=0
       ELSE
          el%C4=-1
          el%C4=0
       ENDIF
       EL%C4%P=>EL%P
       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       EL%C4%AN=>EL%AN
       EL%C4%BN=>EL%BN
       EL%C4%L=>EL%L
       EL%C4%VOLT=>EL%VOLT
       EL%C4%FREQ=>EL%FREQ
       EL%C4%PHAS=>EL%PHAS
       !       EL%C4%P0C=>EL%P0C
       EL%C4%DELTA_E=>EL%DELTA_E
       EL%C4%THIN=>EL%THIN
       ALLOCATE(EL%C4%N_BESSEL);EL%C4%N_BESSEL=0
       ALLOCATE(EL%C4%cavity_totalpath);EL%C4%cavity_totalpath=cavity_totalpath
       ALLOCATE(EL%C4%phase0);EL%C4%phase0=phase0
       ALLOCATE(EL%C4%NF);EL%C4%NF=N_CAV4_F
       ALLOCATE(EL%C4%F(N_CAV4_F));EL%C4%F=0.0_dp;EL%C4%F(1)=1.0_dp;
       ALLOCATE(EL%C4%A);EL%C4%A=0.0_dp;
       ALLOCATE(EL%C4%R);EL%C4%R=1.0_dp;
       ALLOCATE(EL%C4%always_on);EL%C4%always_on=my_false;
       ALLOCATE(EL%C4%PH(N_CAV4_F));EL%C4%PH=0.0_dp;
       ALLOCATE(EL%C4%t);EL%C4%t=0.0_dp;

    CASE(KIND21)
       if(.not.ASSOCIATED(EL%CAV21)) THEN
          ALLOCATE(EL%CAV21)
          el%CAV21=0
       ELSE
          el%CAV21=-1
          el%CAV21=0
       ENDIF
       EL%CAV21%P=>EL%P
       EL%CAV21%L=>EL%L
       EL%CAV21%VOLT=>EL%VOLT
       EL%CAV21%FREQ=>EL%FREQ
       EL%CAV21%PHAS=>EL%PHAS
       !       EL%C4%P0C=>EL%P0C
       EL%CAV21%DELTA_E=>EL%DELTA_E
       EL%CAV21%THIN=>EL%THIN
       ALLOCATE(EL%CAV21%PSI);EL%CAV21%PSI=0.0_dp
       ALLOCATE(EL%CAV21%DVDS);EL%CAV21%DVDS=0.0_dp
       ALLOCATE(EL%CAV21%DPHAS);EL%CAV21%DPHAS=0.0_dp
       ALLOCATE(EL%CAV21%cavity_totalpath);EL%CAV21%cavity_totalpath=cavity_totalpath
       ALLOCATE(EL%CAV21%phase0);EL%CAV21%phase0=phase0
       ALLOCATE(EL%CAV21%always_on);EL%CAV21%always_on=my_false;
    CASE(KIND22)
       if(.not.ASSOCIATED(EL%HE22)) THEN
          ALLOCATE(EL%HE22)
          el%HE22=0
       ELSE
          el%HE22=-1
          el%HE22=0
       ENDIF
       EL%HE22%P=>EL%P
       EL%HE22%L=>EL%L
       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       EL%HE22%AN=>EL%AN
       EL%HE22%BN=>EL%BN
       EL%HE22%FREQ=>EL%FREQ
       EL%HE22%PHAS=>EL%PHAS
       ALLOCATE(EL%HE22%N_BESSEL);EL%HE22%N_BESSEL=0
       ALLOCATE(EL%HE22%fake_shift(6));EL%HE22%fake_shift=0
    CASE(KIND5)
       if(.not.ASSOCIATED(EL%S5)) THEN
          ALLOCATE(EL%S5)
          EL%S5=0
       ELSE
          EL%S5=-1
          EL%S5=0
       ENDIF
       EL%S5%P=>EL%P
       EL%S5%L=>EL%L
       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       EL%S5%AN=>EL%AN
       EL%S5%BN=>EL%BN
       EL%S5%FINT=>EL%FINT      ! added may 31st 2004
       EL%S5%HGAP=>EL%HGAP
       EL%S5%H1=>EL%H1
       EL%S5%H2=>EL%H2
       EL%S5%VA=>EL%VA
       EL%S5%VS=>EL%VS
       EL%S5%B_SOL=>EL%B_SOL
       NULLIFY(EL%s5%dx);ALLOCATE(EL%s5%dx);EL%s5%dx=0.d0;
       NULLIFY(EL%s5%dy);ALLOCATE(EL%s5%dy);EL%s5%dy=0.d0;
       NULLIFY(EL%s5%pitch_x);ALLOCATE(EL%s5%pitch_x);EL%s5%pitch_x=0.d0;
       NULLIFY(EL%s5%pitch_y);ALLOCATE(EL%s5%pitch_y);EL%s5%pitch_y=0.d0;
    CASE(KIND6)
       IF(EL%P%EXACT.AND.EL%P%B0/=0.0_dp) THEN
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,A72,/,1X,A72))'
          !w_p%c(1)=" ERROR IN SETFAMILYR "
          write(*,'(A37,1x,I4)') " EXACT OPTION NOT SUPPORTED FOR KIND ", EL%KIND
          ! call !write_e(777)
       ENDIF
       if(.not.ASSOCIATED(EL%T6)) THEN
          ALLOCATE(EL%T6)
          el%T6=0
       ELSE
          el%T6=-1
          el%T6=0
       ENDIF
       EL%T6%P=>EL%P
       EL%T6%L=>EL%L
       IF(EL%P%NMUL==0)       THEN
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,A72,/,1X,A72))'
          !w_p%c(1)= " ERROR IN SETFAMILYR "
          !w_p%c(2)= " ERROR ON T6: SLOW THICK "
          ! call !write_e(0)
       ENDIF
       EL%T6%AN=>EL%AN
       EL%T6%BN=>EL%BN
       EL%T6%FINT=>EL%FINT
       EL%T6%HGAP=>EL%HGAP
       EL%T6%H1=>EL%H1
       EL%T6%H2=>EL%H2
       EL%T6%VA=>EL%VA
       EL%T6%VS=>EL%VS
       nullify(EL%T6%MATX);ALLOCATE(EL%T6%MATX(2,3));
       nullify(EL%T6%MATY);ALLOCATE(EL%T6%MATY(2,3));
       nullify(EL%T6%LX);ALLOCATE(EL%T6%LX(6));
       nullify(EL%T6%LY);ALLOCATE(EL%T6%LY(3));
    CASE(KIND7)
       IF(EL%P%EXACT.AND.EL%P%B0/=0.0_dp) THEN
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,A72,/,1X,A72))'
          !w_p%c(1)=" ERROR IN SETFAMILYR "
          write(*,'(A37,1x,I4)') " EXACT OPTION NOT SUPPORTED FOR KIND ", EL%KIND
          ! call !write_e(777)
       ENDIF
       !       if(.not.ASSOCIATED(EL%T7))ALLOCATE(EL%T7)
       if(.not.ASSOCIATED(EL%T7)) THEN
          ALLOCATE(EL%T7)
          EL%T7=0
       ELSE
          EL%T7=-1
          EL%T7=0
       ENDIF
       EL%T7%P=>EL%P
       EL%T7%L=>EL%L
       IF(EL%P%NMUL==0)       THEN
          !w_p=0
          !w_p%nc=1
          !w_p%fc='((1X,A72))'
          !w_p%c(1)= "ERROR ON T7: FAST THICK "
          ! call !write_e(0)
       ENDIF
       EL%T7%AN=>EL%AN
       EL%T7%BN=>EL%BN
       EL%T7%FINT=>EL%FINT
       EL%T7%HGAP=>EL%HGAP
       EL%T7%H1=>EL%H1
       EL%T7%H2=>EL%H2
       EL%T7%VA=>EL%VA
       EL%T7%VS=>EL%VS
       NULLIFY(EL%T7%F);ALLOCATE(EL%T7%F);EL%T7%F=1;
       nullify(EL%T7%MATX);ALLOCATE(EL%T7%MATX(2,3));
       nullify(EL%T7%MATY);ALLOCATE(EL%T7%MATY(2,3));
       nullify(EL%T7%LX);ALLOCATE(EL%T7%LX(3));
       nullify(EL%T7%RMATX);ALLOCATE(EL%T7%RMATX(2,3));
       nullify(EL%T7%RMATY);ALLOCATE(EL%T7%RMATY(2,3));
       nullify(EL%T7%RLX);ALLOCATE(EL%T7%RLX(3));
       IF(GEN) call GETMAT7(EL%T7)
    CASE(KIND8)
       if(.not.ASSOCIATED(EL%S8))ALLOCATE(EL%S8)
       EL%S8%P=>EL%P
       IF(EL%P%NMUL==0)       THEN
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,A72,/,1X,A72))'
          !w_p%c(1)= " ERROR IN SETFAMILYR "
          !w_p%c(2)= "ERROR ON S8:  NORMAL SMI "
          ! call !write_e(0)
       ENDIF
       EL%S8%BN=>EL%BN
    CASE(KIND9)
       if(.not.ASSOCIATED(EL%S9))ALLOCATE(EL%S9)
       EL%S9%P=>EL%P
       IF(EL%P%NMUL==0)       THEN
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,A72,/,1X,A72))'
          !w_p%c(1)= " ERROR IN SETFAMILYR "
          !w_p%c(2)= "ERROR ON S9: SKEW SMI "
          ! call !write_e(0)
       ENDIF
       EL%S9%AN=>EL%AN
    CASE(KIND10)
       IF(.not.EL%P%EXACT) THEN
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,A72,/,1X,A72))'
          !w_p%c(1)=" ERROR IN SETFAMILYR "
          write(*,'(A43,1x,I4)') " 'NOT EXACT' OPTION NOT SUPPORTED FOR KIND ", EL%KIND
          ! call !write_e(777)
       ENDIF
       if(.not.ASSOCIATED(EL%TP10)) THEN
          ALLOCATE(EL%TP10)
          EL%TP10=0
       ELSE
          EL%TP10=-1
          EL%TP10=0
       ENDIF
       EL%TP10%P=>EL%P
       EL%TP10%L=>EL%L
       IF(EL%P%NMUL==0.OR.EL%P%NMUL>SECTOR_NMUL)       THEN
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,A72,/,1X,A72))'
          !w_p%c(1)= " ERROR IN SETFAMILYR "
          !w_p%c(2)= "ERROR ON TP10: TEAPOT "
          ! call !write_e(0)
       ENDIF
       EL%TP10%AN=>EL%AN
       EL%TP10%BN=>EL%BN
       EL%TP10%FINT=>EL%FINT
       EL%TP10%HGAP=>EL%HGAP
       EL%TP10%H1=>EL%H1
       EL%TP10%H2=>EL%H2
       EL%TP10%VA=>EL%VA
       EL%TP10%VS=>EL%VS


       NULLIFY(EL%TP10%BF_X);
       NULLIFY(EL%TP10%BF_Y);

       EL%TP10%ELECTRIC=>EL%ELECTRIC

        ALLOCATE(EL%TP10%BF_X(S_E%N_MONO))
        ALLOCATE(EL%TP10%BF_Y(S_E%N_MONO))
 

        NULLIFY(EL%TP10%VM);ALLOCATE(EL%TP10%VM(S_E%N_MONO))

        EL%TP10%VM=0.0_dp
        EL%TP10%BF_X=0.0_dp
        EL%TP10%BF_Y=0.0_dp

       NULLIFY(EL%TP10%DRIFTKICK);ALLOCATE(EL%TP10%DRIFTKICK);EL%TP10%DRIFTKICK=.true.;
!       if(EL%ELECTRIC) then
        NULLIFY(EL%TP10%E_X);ALLOCATE(EL%TP10%E_X(S_E%N_MONO))
        NULLIFY(EL%TP10%E_Y);ALLOCATE(EL%TP10%E_Y(S_E%N_MONO))
        NULLIFY(EL%TP10%PHI);ALLOCATE(EL%TP10%PHI(S_E%N_MONO))

        NULLIFY(EL%TP10%AE);ALLOCATE(EL%TP10%AE(SECTOR_NMUL_max));
        NULLIFY(EL%TP10%BE);ALLOCATE(EL%TP10%BE(SECTOR_NMUL_max));


        EL%TP10%E_X=0.0_dp
        EL%TP10%E_Y=0.0_dp
        EL%TP10%PHI=0.0_dp

        EL%TP10%AE=0.0_DP;
        EL%TP10%BE=0.0_DP;
        call GETAEBE(EL%TP10) ! not efective here because ae=be=0 but need on magnetic field
 !      ELSE
 !       call GETANBN(EL%TP10)  
 !      endif

       NULLIFY(EL%TP10%F);ALLOCATE(EL%TP10%F);EL%TP10%F=1;
    CASE(KIND11:KIND14)
       if(.not.ASSOCIATED(EL%MON14)) THEN
          ALLOCATE(EL%MON14)
          el%MON14=0
       ELSE
          el%MON14=-1
          el%MON14=0
       ENDIF
       EL%MON14%P=>EL%P
       EL%MON14%L=>EL%L
       nullify(EL%MON14%X);ALLOCATE(EL%MON14%X);EL%MON14%X=0.0_dp;
       nullify(EL%MON14%Y);ALLOCATE(EL%MON14%Y);EL%MON14%Y=0.0_dp
    CASE(KIND15)
       if(.not.ASSOCIATED(EL%SEP15))ALLOCATE(EL%SEP15)
       EL%SEP15%P=>EL%P
       EL%SEP15%L=>EL%L
       EL%SEP15%VOLT=>EL%VOLT
       EL%SEP15%PHAS=>EL%PHAS
    CASE(KIND16,KIND20)
       if(.not.ASSOCIATED(EL%K16)) THEN
          ALLOCATE(EL%K16)
          el%K16=0
       ELSE
          el%K16=-1
          el%K16=0
       ENDIF
       EL%K16%P=>EL%P
       EL%K16%L=>EL%L
       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       EL%K16%AN=>EL%AN
       EL%K16%BN=>EL%BN
       EL%K16%FINT=>EL%FINT
       EL%K16%HGAP=>EL%HGAP
       EL%K16%H1=>EL%H1
       EL%K16%H2=>EL%H2
       EL%K16%VA=>EL%VA
       EL%K16%VS=>EL%VS
       NULLIFY(EL%K16%DRIFTKICK);ALLOCATE(EL%K16%DRIFTKICK);EL%K16%DRIFTKICK=.true.;
       NULLIFY(EL%K16%LIKEMAD);ALLOCATE(EL%K16%LIKEMAD);EL%K16%LIKEMAD=.false.;
       NULLIFY(EL%K16%F);ALLOCATE(EL%K16%F);EL%K16%F=1;
    CASE(KIND17)
       if(.not.ASSOCIATED(EL%ENGE17)) THEN
          ALLOCATE(EL%ENGE17)
          el%ENGE17=0
       ELSE
          el%ENGE17=-1
          el%ENGE17=0
       ENDIF
       EL%ENGE17%P=>EL%P
       EL%ENGE17%L=>EL%L
       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       EL%ENGE17%AN=>EL%AN
       EL%ENGE17%BN=>EL%BN
       NULLIFY(EL%ENGE17%F);ALLOCATE(EL%ENGE17%F);EL%ENGE17%F=1.0_dp;
       NULLIFY(EL%ENGE17%D);ALLOCATE(EL%ENGE17%D);EL%ENGE17%D=1.0_dp;
       NULLIFY(EL%ENGE17%A);ALLOCATE(EL%ENGE17%A(0:N_ENGE));EL%ENGE17%A=0.0_dp;
       NULLIFY(EL%ENGE17%nbessel);ALLOCATE(EL%ENGE17%nbessel);EL%ENGE17%nbessel=0;
    CASE(KIND18)
       if(.not.ASSOCIATED(EL%RCOL18)) THEN
          ALLOCATE(EL%RCOL18)
          EL%RCOL18=0
       ELSE
          EL%RCOL18=-1
          EL%RCOL18=0
       ENDIF
       EL%RCOL18%P=>EL%P
       EL%RCOL18%L=>EL%L
!       nullify(EL%RCOL18%A);!ALLOCATE(EL%RCOL18%A);CALL ALLOC(EL%RCOL18%A)
    CASE(KIND19)
       if(.not.ASSOCIATED(EL%ECOL19)) THEN
          ALLOCATE(EL%ECOL19)
          EL%ECOL19=0
       ELSE
          EL%ECOL19=-1
          EL%ECOL19=0
       ENDIF
       EL%ECOL19%P=>EL%P
       EL%ECOL19%L=>EL%L
!       nullify(EL%ECOL19%A);!ALLOCATE(EL%ECOL19%A);CALL ALLOC(EL%ECOL19%A)
    CASE(KINDWIGGLER)
       if(.not.ASSOCIATED(EL%WI)) THEN
          ALLOCATE(EL%WI)
          EL%WI=0
       ELSE
          EL%WI=-1
          EL%WI=0
       ENDIF
       EL%WI%P=>EL%P
       EL%WI%L=>EL%L
       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       EL%WI%AN=>EL%AN
       EL%WI%BN=>EL%BN
       CALL POINTERS_SAGAN(EL%WI)
    CASE(KINDpa)
       if(.not.ASSOCIATED(EL%pa)) THEN
          ALLOCATE(EL%pa)
          EL%PA=0
       ELSE
          EL%pa=-1
          EL%pa=0
       ENDIF
       EL%pa%P=>EL%P
       EL%pa%L=>EL%L
       !       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       !       EL%mu%AN=>EL%AN
       !       EL%mu%BN=>EL%BN
       CALL POINTERS_pancake(EL%pa,T) !,angc,xc,dc,h) !,t_ax,t_ay)
    CASE(KINDabell)
       if(.not.ASSOCIATED(EL%ab)) THEN
          ALLOCATE(EL%ab)
          EL%ab=0
       ELSE
          EL%ab=-1
          EL%ab=0
       ENDIF
       EL%ab%P=>EL%P
       EL%ab%L=>EL%L
       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       EL%ab%AN=>EL%AN
       EL%ab%BN=>el%BN
   !real(dp), POINTER :: DZ => null(), T(:) => null()
   !complex(dp), POINTER :: B(:,:) => null()
   !INTEGER , POINTER :: N,M => null()  

       ALLOCATE(EL%ab%m);EL%ab%m=m_abell;
       ALLOCATE(EL%ab%n);EL%ab%n=n_abell;
       ALLOCATE(EL%ab%dz(0:m_abell));EL%ab%dz=0
       ALLOCATE(EL%ab%t(0:m_abell));EL%ab%t=0.0_dp;
       ALLOCATE(EL%ab%b(0:m_abell,-n_abell/2:n_abell/2-1));EL%ab%b=0.0_dp; 
       ALLOCATE(EL%ab%tE(0:m_abell));EL%ab%tE=0.0_dp;
       ALLOCATE(EL%ab%E(0:m_abell,-n_abell/2:n_abell/2-1));EL%ab%E=0.0_dp; 
       CALL POINTERS_abell(EL%ab) !,angc,xc,dc,h) !,t_ax,t_ay)
    END SELECT
  END SUBROUTINE SETFAMILYR


  SUBROUTINE SETFAMILYP(EL,T)  !,angc,xc,dc,h)  !,NTOT,ntot_rad,NTOT_REV,ntot_rad_REV,ND2)
    !  SUBROUTINE SETFAMILYP(EL,T,t_ax,t_ay,NTOT,ntot_rad,NTOT_REV,ntot_rad_REV,ND2)
    IMPLICIT NONE
    TYPE(ELEMENTP), INTENT(INOUT) ::EL
    !    INTEGER,OPTIONAL :: NTOT,ntot_rad,NTOT_REV,ntot_rad_REV,ND2
    type(tree_element),OPTIONAL :: T(:) !,t_ax(:),t_ay(:)
 !   real(dp), optional :: angc,xc,dc,h
!    EL%P%permfringe=>EL%permfringe
    SELECT CASE(EL%KIND)
    CASE(KIND1)
       if(.not.ASSOCIATED(EL%D0))ALLOCATE(EL%D0)
       EL%D0%P=>EL%P
       EL%D0%L=>EL%L
    CASE(KIND2)
       IF(EL%P%EXACT) THEN
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,A72,/,1X,A72))'
          !w_p%c(1)=" ERROR IN SETFAMILYP "
          write(*,'(A43,1x,I4)') " 'NOT EXACT' OPTION NOT SUPPORTED FOR KIND ", EL%KIND
          ! call !write_e(222)
       ENDIF
       if(.not.ASSOCIATED(EL%K2)) THEN
          ALLOCATE(EL%K2)
          EL%K2=0
       ELSE
          EL%K2=-1
          EL%K2=0
       ENDIF
       !       if(.not.ASSOCIATED(EL%K2))ALLOCATE(EL%K2)
       EL%K2%P=>EL%P
       EL%K2%L=>EL%L
       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       EL%K2%AN=>EL%AN
       EL%K2%BN=>EL%BN
       EL%K2%FINT=>EL%FINT
       EL%K2%HGAP=>EL%HGAP
       EL%K2%H1=>EL%H1
       EL%K2%H2=>EL%H2
       EL%K2%VA=>EL%VA
       EL%K2%VS=>EL%VS
       NULLIFY(EL%K2%F);ALLOCATE(EL%K2%F);EL%K2%F=1;
    CASE(KIND3)
       if(.not.ASSOCIATED(EL%K3)) THEN
          ALLOCATE(EL%K3)
          el%K3=0
       ELSE
          el%K3=-1
          el%K3=0
       ENDIF
       EL%K3%P=>EL%P
       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       EL%K3%AN=>EL%AN
       EL%K3%BN=>EL%BN
       ALLOCATE(EL%K3%ls);EL%K3%ls=1
       ALLOCATE(EL%K3%hf);CALL ALLOC(EL%K3%hf);EL%K3%hf=0.0_dp
       ALLOCATE(EL%K3%vf);CALL ALLOC(EL%K3%vf);EL%K3%vf=0.0_dp
       ALLOCATE(EL%K3%thin_h_foc);CALL ALLOC(EL%K3%thin_h_foc);EL%K3%thin_h_foc=0.0_dp
       ALLOCATE(EL%K3%thin_v_foc);CALL ALLOC(EL%K3%thin_v_foc);EL%K3%thin_v_foc=0.0_dp
       ALLOCATE(EL%K3%thin_h_angle);CALL ALLOC(EL%K3%thin_h_angle);EL%K3%thin_h_angle=0.0_dp
       ALLOCATE(EL%K3%thin_v_angle);CALL ALLOC(EL%K3%thin_v_angle);EL%K3%thin_v_angle=0.0_dp
       ALLOCATE(EL%K3%patch);EL%K3%patch=my_false
       EL%K3%B_SOL=>EL%B_SOL
       NULLIFY(EL%k3%dx);ALLOCATE(EL%k3%dx);EL%k3%dx=0.d0;
       NULLIFY(EL%k3%dy);ALLOCATE(EL%k3%dy);EL%k3%dy=0.d0;
       NULLIFY(EL%k3%pitch_x);ALLOCATE(EL%k3%pitch_x);EL%k3%pitch_x=0.d0;
       NULLIFY(EL%k3%pitch_y);ALLOCATE(EL%k3%pitch_y);EL%k3%pitch_y=0.d0;
    CASE(kindsuperdrift)
       if(.not.ASSOCIATED(EL%sdr)) THEN
          ALLOCATE(EL%sdr)
          EL%sdr=0
       ELSE
          EL%sdr=-1
          EL%sdr=0
       ENDIF
       EL%SDR%P=>EL%P
       EL%SDR%L=>EL%L
 
       ALLOCATE(EL%SDR%a_x1);EL%SDR%a_x1=1.0_dp;
       ALLOCATE(EL%SDR%a_x2);EL%SDR%a_x2=1.0_dp;
       ALLOCATE(EL%SDR%D(3));EL%SDR%D=0.0_dp;
       ALLOCATE(EL%SDR%ANG(3));EL%SDR%ANG=0.0_dp;

    CASE(KIND4)
       if(.not.ASSOCIATED(EL%C4)) THEN
          ALLOCATE(EL%C4)
          el%C4=0
       ELSE
          el%C4=-1
          el%C4=0
       ENDIF
       EL%C4%P=>EL%P
       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       EL%C4%AN=>EL%AN
       EL%C4%BN=>EL%BN
       EL%C4%L=>EL%L
       EL%C4%VOLT=>EL%VOLT
       EL%C4%FREQ=>EL%FREQ
       EL%C4%PHAS=>EL%PHAS
       !       EL%C4%P0C=>EL%P0C
       EL%C4%DELTA_E=>EL%DELTA_E
       EL%C4%THIN=>EL%THIN
       ALLOCATE(EL%C4%N_BESSEL);EL%C4%N_BESSEL=0
       ALLOCATE(EL%C4%cavity_totalpath);EL%C4%cavity_totalpath=cavity_totalpath
       ALLOCATE(EL%C4%phase0);EL%C4%phase0=phase0
       ALLOCATE(EL%C4%NF);EL%C4%NF=N_CAV4_F
       ALLOCATE(EL%C4%F(N_CAV4_F));CALL ALLOC(EL%C4%F,N_CAV4_F);EL%C4%F(1)=1.0_dp;
       ALLOCATE(EL%C4%A);CALL ALLOC(EL%C4%A);EL%C4%A=0.0_dp;
       ALLOCATE(EL%C4%R);CALL ALLOC(EL%C4%R);EL%C4%R=1.0_dp;
       ALLOCATE(EL%C4%always_on);EL%C4%always_on=my_false;
       ALLOCATE(EL%C4%PH(N_CAV4_F));CALL ALLOC(EL%C4%PH,N_CAV4_F);
       ALLOCATE(EL%C4%t);EL%C4%t=0.0_dp;
    CASE(KIND21)
       if(.not.ASSOCIATED(EL%CAV21)) THEN
          ALLOCATE(EL%CAV21)
          el%CAV21=0
       ELSE
          el%CAV21=-1
          el%CAV21=0
       ENDIF
       EL%CAV21%P=>EL%P
       EL%CAV21%L=>EL%L
       EL%CAV21%VOLT=>EL%VOLT
       EL%CAV21%FREQ=>EL%FREQ
       EL%CAV21%PHAS=>EL%PHAS
       !       EL%C4%P0C=>EL%P0C
       EL%CAV21%DELTA_E=>EL%DELTA_E
       EL%CAV21%THIN=>EL%THIN
       ALLOCATE(EL%CAV21%PSI);CALL ALLOC(EL%CAV21%PSI);EL%CAV21%PSI=0.0_dp
       ALLOCATE(EL%CAV21%DVDS);CALL ALLOC(EL%CAV21%DVDS);EL%CAV21%DVDS=0.0_dp
       ALLOCATE(EL%CAV21%DPHAS);CALL ALLOC(EL%CAV21%DPHAS);EL%CAV21%DPHAS=0.0_dp
       ALLOCATE(EL%CAV21%cavity_totalpath);EL%CAV21%cavity_totalpath=cavity_totalpath
       ALLOCATE(EL%CAV21%always_on);EL%CAV21%always_on=my_false;
       ALLOCATE(EL%CAV21%phase0);EL%CAV21%phase0=phase0
    CASE(KIND22)
       if(.not.ASSOCIATED(EL%HE22)) THEN
          ALLOCATE(EL%HE22)
          el%HE22=0
       ELSE
          el%HE22=-1
          el%HE22=0
       ENDIF
       EL%HE22%P=>EL%P
       EL%HE22%L=>EL%L
       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       EL%HE22%AN=>EL%AN
       EL%HE22%BN=>EL%BN
       EL%HE22%FREQ=>EL%FREQ
       EL%HE22%PHAS=>EL%PHAS
       ALLOCATE(EL%HE22%N_BESSEL);EL%HE22%N_BESSEL=0
       ALLOCATE(EL%HE22%fake_shift(6));call alloc(EL%HE22%fake_shift);EL%HE22%fake_shift=0
    CASE(KIND5)
       if(.not.ASSOCIATED(EL%S5)) THEN
          ALLOCATE(EL%S5)
          EL%S5=0
       ELSE
          EL%S5=-1
          EL%S5=0
       ENDIF
       EL%S5%P=>EL%P
       EL%S5%L=>EL%L
       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       EL%S5%AN=>EL%AN
       EL%S5%BN=>EL%BN
       EL%S5%FINT=>EL%FINT      ! added may 31st 2004
       EL%S5%HGAP=>EL%HGAP
       EL%S5%H1=>EL%H1
       EL%S5%H2=>EL%H2
       EL%S5%VA=>EL%VA
       EL%S5%VS=>EL%VS
       EL%S5%B_SOL=>EL%B_SOL
       NULLIFY(EL%s5%dx);ALLOCATE(EL%s5%dx);EL%s5%dx=0.d0;
       NULLIFY(EL%s5%dy);ALLOCATE(EL%s5%dy);EL%s5%dy=0.d0;
       NULLIFY(EL%s5%pitch_x);ALLOCATE(EL%s5%pitch_x);EL%s5%pitch_x=0.d0;
       NULLIFY(EL%s5%pitch_y);ALLOCATE(EL%s5%pitch_y);EL%s5%pitch_y=0.d0;
    CASE(KIND6)
       IF(EL%P%EXACT.AND.EL%P%B0/=0.0_dp) THEN
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,A72,/,1X,A72))'
          !w_p%c(1)=" ERROR IN SETFAMILYP "
          write(*,'(A37,1x,I4)') " EXACT OPTION NOT SUPPORTED FOR KIND ", EL%KIND
          ! call !write_e(777)
       ENDIF
       if(.not.ASSOCIATED(EL%T6)) THEN
          ALLOCATE(EL%T6)
          el%T6=0
       ELSE
          el%T6=-1
          el%T6=0
       ENDIF
       EL%T6%P=>EL%P
       EL%T6%L=>EL%L
       IF(EL%P%NMUL==0)       THEN
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,A72,/,1X,A72))'
          !w_p%c(1)= " ERROR IN SETFAMILYP "
          !w_p%c(2)= "ERROR ON T6: SLOW THICK "
          ! call !write_e(0)
       ENDIF
       EL%T6%AN=>EL%AN
       EL%T6%BN=>EL%BN
       EL%T6%FINT=>EL%FINT
       EL%T6%HGAP=>EL%HGAP
       EL%T6%H1=>EL%H1
       EL%T6%H2=>EL%H2
       EL%T6%VA=>EL%VA
       EL%T6%VS=>EL%VS
       nullify(EL%T6%MATX);ALLOCATE(EL%T6%MATX(2,3));
       nullify(EL%T6%MATY);ALLOCATE(EL%T6%MATY(2,3));
       nullify(EL%T6%LX);ALLOCATE(EL%T6%LX(6));
       nullify(EL%T6%LY);ALLOCATE(EL%T6%LY(3));
    CASE(KIND7)
       IF(EL%P%EXACT.AND.EL%P%B0/=0.0_dp) THEN
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,A72,/,1X,A72))'
          !w_p%c(1)=" ERROR IN SETFAMILYP "
          write(*,'(A37,1x,I4)') " EXACT OPTION NOT SUPPORTED FOR KIND ", EL%KIND
          ! call !write_e(777)
       ENDIF
       if(.not.ASSOCIATED(EL%T7)) THEN
          ALLOCATE(EL%T7)
          EL%T7=0
       ELSE
          EL%T7=-1
          EL%T7=0
       ENDIF
       EL%T7%P=>EL%P
       EL%T7%L=>EL%L
       IF(EL%P%NMUL==0)       THEN
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,A72,/,1X,A72))'
          !w_p%c(1)= " ERROR IN SETFAMILYP "
          !w_p%c(2)= "ERROR ON T7: FAST THICK "
          ! call !write_e(0)
       ENDIF
       EL%T7%AN=>EL%AN
       EL%T7%BN=>EL%BN
       EL%T7%FINT=>EL%FINT
       EL%T7%HGAP=>EL%HGAP
       EL%T7%H1=>EL%H1
       EL%T7%H2=>EL%H2
       EL%T7%VA=>EL%VA
       EL%T7%VS=>EL%VS
       NULLIFY(EL%T7%F);ALLOCATE(EL%T7%F);EL%T7%F=1;
       nullify(EL%T7%MATX);  ALLOCATE(EL%T7%MATX(2,3));
       nullify(EL%T7%MATY);  ALLOCATE(EL%T7%MATY(2,3));
       nullify(EL%T7%LX);    ALLOCATE(EL%T7%LX(3));
       nullify(EL%T7%RMATX); ALLOCATE(EL%T7%RMATX(2,3));
       nullify(EL%T7%RMATY); ALLOCATE(EL%T7%RMATY(2,3));
       nullify(EL%T7%RLX);   ALLOCATE(EL%T7%RLX(3));
       CALL ALLOC(EL%T7)
       IF(GEN) call GETMAT7(EL%T7)
    CASE(KIND8)
       if(.not.ASSOCIATED(EL%S8))ALLOCATE(EL%S8)
       EL%S8%P=>EL%P
       IF(EL%P%NMUL==0)       THEN
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,A72,/,1X,A72))'
          !w_p%c(1)= " ERROR IN SETFAMILYP "
          !w_p%c(2)= "ERROR ON S8:  NORMAL SMI "
          ! call !write_e(0)
       ENDIF
       EL%S8%BN=>EL%BN
    CASE(KIND9)
       if(.not.ASSOCIATED(EL%S9))ALLOCATE(EL%S9)
       EL%S9%P=>EL%P
       IF(EL%P%NMUL==0)       THEN
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,A72,/,1X,A72))'
          !w_p%c(1)= " ERROR IN SETFAMILYP "
          !w_p%c(2)= "ERROR ON S9: SKEW SMI "
          ! call !write_e(0)
       ENDIF
       EL%S9%AN=>EL%AN
    CASE(KIND10)
       IF(.not.EL%P%EXACT) THEN
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,A72,/,1X,A72))'
          !w_p%c(1)=" ERROR IN SETFAMILYP "
          write(*,'(A43,1x,I4)') " 'NOT EXACT' OPTION NOT SUPPORTED FOR KIND ", EL%KIND
          ! call !write_e(777)
       ENDIF
       if(.not.ASSOCIATED(EL%TP10)) THEN
          ALLOCATE(EL%TP10)
          EL%TP10=0
       ELSE
          EL%TP10=-1
          EL%TP10=0
       ENDIF
       EL%TP10%P=>EL%P
       EL%TP10%L=>EL%L
       IF(EL%P%NMUL==0.OR.EL%P%NMUL>SECTOR_NMUL)       THEN
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,A72,/,1X,A72))'
          !w_p%c(1)= " ERROR IN SETFAMILYP "
          !w_p%c(2)= "ERROR ON TP10: TEAPOT "
          ! call !write_e(0)
       ENDIF
       EL%TP10%AN=>EL%AN
       EL%TP10%BN=>EL%BN
       EL%TP10%FINT=>EL%FINT
       EL%TP10%HGAP=>EL%HGAP
       EL%TP10%H1=>EL%H1
       EL%TP10%H2=>EL%H2
       EL%TP10%VA=>EL%VA
       EL%TP10%VS=>EL%VS
       NULLIFY(EL%TP10%BF_X); 
       NULLIFY(EL%TP10%BF_Y);
       EL%TP10%ELECTRIC=>EL%ELECTRIC
        ALLOCATE(EL%TP10%BF_X(S_B_FROM_V%N_MONO))
        ALLOCATE(EL%TP10%BF_Y(S_B_FROM_V%N_MONO))
       NULLIFY(EL%TP10%DRIFTKICK);ALLOCATE(EL%TP10%DRIFTKICK);EL%TP10%DRIFTKICK=.true.;
  !     CALL ALLOC(EL%TP10)

        NULLIFY(EL%TP10%VM);ALLOCATE(EL%TP10%VM(S_E%N_MONO))


 !      if(EL%ELECTRIC) then
        NULLIFY(EL%TP10%E_X);ALLOCATE(EL%TP10%E_X(S_E%N_MONO))
        NULLIFY(EL%TP10%E_Y);ALLOCATE(EL%TP10%E_Y(S_E%N_MONO))
        NULLIFY(EL%TP10%PHI);ALLOCATE(EL%TP10%PHI(S_E%N_MONO))

        NULLIFY(EL%TP10%AE);ALLOCATE(EL%TP10%AE(SECTOR_NMUL_max)); 
        NULLIFY(EL%TP10%BE);ALLOCATE(EL%TP10%BE(SECTOR_NMUL_max)); 

        call alloc(EL%TP10)
        call GETAEBE(EL%TP10) ! not efective here because ae=be=0 but need on magnetic field
 !      ELSE
!        call alloc(EL%TP10)
!        call GETANBN(EL%TP10)
!       endif

       NULLIFY(EL%TP10%F);ALLOCATE(EL%TP10%F);EL%TP10%F=1;
    CASE(KIND11:KIND14)
       if(.not.ASSOCIATED(EL%MON14)) THEN
          ALLOCATE(EL%MON14)
          el%MON14=0
       ELSE
          el%MON14=-1
          el%MON14=0
       ENDIF
       EL%MON14%P=>EL%P
       EL%MON14%L=>EL%L
       nullify(EL%MON14%X);ALLOCATE(EL%MON14%X);EL%MON14%X=0.0_dp;
       nullify(EL%MON14%Y);ALLOCATE(EL%MON14%Y);EL%MON14%Y=0.0_dp
    CASE(KIND15)
       if(.not.ASSOCIATED(EL%SEP15))ALLOCATE(EL%SEP15)
       EL%SEP15%P=>EL%P
       EL%SEP15%L=>EL%L
       EL%SEP15%VOLT=>EL%VOLT
       EL%SEP15%PHAS=>EL%PHAS
    CASE(KIND16,KIND20)
       if(.not.ASSOCIATED(EL%K16)) THEN
          ALLOCATE(EL%K16)
          el%K16=0
       ELSE
          el%K16=-1
          el%K16=0
       ENDIF
       EL%K16%P=>EL%P
       EL%K16%L=>EL%L
       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       EL%K16%AN=>EL%AN
       EL%K16%BN=>EL%BN
       EL%K16%FINT=>EL%FINT
       EL%K16%HGAP=>EL%HGAP
       EL%K16%H1=>EL%H1
       EL%K16%H2=>EL%H2
       EL%K16%VA=>EL%VA
       EL%K16%VS=>EL%VS
       NULLIFY(EL%K16%DRIFTKICK);ALLOCATE(EL%K16%DRIFTKICK);EL%K16%DRIFTKICK=.true.;
       NULLIFY(EL%K16%LIKEMAD);ALLOCATE(EL%K16%LIKEMAD);EL%K16%LIKEMAD=.false.;
       NULLIFY(EL%K16%F);ALLOCATE(EL%K16%F);EL%K16%F=1;
    CASE(KIND17)
       if(.not.ASSOCIATED(EL%ENGE17)) THEN
          ALLOCATE(EL%ENGE17)
          el%ENGE17=0
       ELSE
          el%ENGE17=-1
          el%ENGE17=0
       ENDIF
       EL%ENGE17%P=>EL%P
       EL%ENGE17%L=>EL%L
       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       EL%ENGE17%AN=>EL%AN
       EL%ENGE17%BN=>EL%BN
       NULLIFY(EL%ENGE17%F);ALLOCATE(EL%ENGE17%F);EL%ENGE17%F=1.0_dp;
       NULLIFY(EL%ENGE17%D);ALLOCATE(EL%ENGE17%D);EL%ENGE17%D=1.0_dp;
       NULLIFY(EL%ENGE17%A);ALLOCATE(EL%ENGE17%A(0:N_ENGE));EL%ENGE17%A=0.0_dp;
       NULLIFY(EL%ENGE17%nbessel);ALLOCATE(EL%ENGE17%nbessel);EL%ENGE17%nbessel=0;
    CASE(KIND18)
       if(.not.ASSOCIATED(EL%RCOL18)) THEN
          ALLOCATE(EL%RCOL18)
          EL%RCOL18=0
       ELSE
          EL%RCOL18=-1
          EL%RCOL18=0
       ENDIF
       EL%RCOL18%P=>EL%P
       EL%RCOL18%L=>EL%L
!       nullify(EL%RCOL18%A);!ALLOCATE(EL%RCOL18%A);CALL ALLOC(EL%RCOL18%A)
    CASE(KIND19)
       if(.not.ASSOCIATED(EL%ECOL19)) THEN
          ALLOCATE(EL%ECOL19)
          EL%ECOL19=0
       ELSE
          EL%ECOL19=-1
          EL%ECOL19=0
       ENDIF
       EL%ECOL19%P=>EL%P
       EL%ECOL19%L=>EL%L

    CASE(KINDWIGGLER)
       if(.not.ASSOCIATED(EL%WI)) THEN
          ALLOCATE(EL%WI)
          EL%WI=0
       ELSE
          EL%WI=-1
          EL%WI=0
       ENDIF
       EL%WI%P=>EL%P
       EL%WI%L=>EL%L
       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       EL%WI%AN=>EL%AN
       EL%WI%BN=>EL%BN
       CALL POINTERS_SAGAN(EL%WI)
       CALL ALLOC(EL%WI)
    CASE(KINDpa)
       if(.not.ASSOCIATED(EL%pa)) THEN
          ALLOCATE(EL%pa)
          EL%PA=0
       ELSE
          EL%pa=-1
          EL%pa=0
       ENDIF
       EL%pa%P=>EL%P
       EL%pa%L=>EL%L
       !       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       !       EL%mu%AN=>EL%AN
       !       EL%mu%BN=>EL%BN

       CALL POINTERS_pancake(EL%pa,T) !,angc,xc,dc,h)  !,t_ax,t_ay)

    CASE(KINDabell)
       if(.not.ASSOCIATED(EL%ab)) THEN
          ALLOCATE(EL%ab)
          EL%ab=0
       ELSE
          EL%ab=-1
          EL%ab=0
       ENDIF
       EL%ab%P=>EL%P
       EL%ab%L=>EL%L
       IF(EL%P%NMUL==0) CALL ZERO_ANBN(EL,1)
       EL%ab%AN=>EL%AN
       EL%ab%BN=>el%BN
   !real(dp), POINTER :: DZ => null(), T(:) => null()
   !complex(dp), POINTER :: B(:,:) => null()
   !INTEGER , POINTER :: N,M => null()  


       ALLOCATE(EL%ab%m);EL%ab%m=m_abell;
       ALLOCATE(EL%ab%n);EL%ab%n=n_abell;
       ALLOCATE(EL%ab%dz(0:m_abell));EL%ab%dz=0
       ALLOCATE(EL%ab%t(0:m_abell));EL%ab%t=0.0_dp;
       ALLOCATE(EL%ab%b(0:m_abell,-n_abell/2:n_abell/2-1));EL%ab%b=0.0_dp; 
       ALLOCATE(EL%ab%tE(0:m_abell));EL%ab%tE=0.0_dp;
       ALLOCATE(EL%ab%E(0:m_abell,-n_abell/2:n_abell/2-1));EL%ab%E=0.0_dp; 
       CALL POINTERS_abell(EL%ab) !,angc,xc,dc,h) !,t_ax,t_ay)
    END SELECT

  END SUBROUTINE SETFAMILYP




  SUBROUTINE ZERO_ANBN_R(EL,N)
    IMPLICIT NONE
    TYPE(ELEMENT), INTENT(INOUT) ::EL
    INTEGER, INTENT(IN) ::N
    INTEGER I

    IF(N<=0) RETURN
    IF(ASSOCIATED(EL%AN)) DEALLOCATE(EL%AN)
    IF(ASSOCIATED(EL%BN)) DEALLOCATE(EL%BN)
    EL%p%NMUL=N
    ALLOCATE(EL%AN(EL%p%NMUL),EL%BN(EL%p%NMUL))

    DO I=1,EL%P%NMUL
       EL%AN(I)=0.0_dp
       EL%BN(I)=0.0_dp
    ENDDO

  END SUBROUTINE ZERO_ANBN_R

  SUBROUTINE ZERO_ANBN_P(EL,N)
    IMPLICIT NONE
    TYPE(ELEMENTP), INTENT(INOUT) ::EL
    INTEGER, INTENT(IN) ::N

    IF(N<=0) RETURN
    IF(ASSOCIATED(EL%AN)) DEALLOCATE(EL%AN)
    IF(ASSOCIATED(EL%BN)) DEALLOCATE(EL%BN)
    EL%P%NMUL=N
    ALLOCATE(EL%AN(EL%P%NMUL),EL%BN(EL%P%NMUL))
    CALL ALLOC(EL%AN,EL%P%NMUL);CALL ALLOC(EL%BN,EL%P%NMUL);

  END SUBROUTINE ZERO_ANBN_P

  SUBROUTINE transfer_ANBN(EL,ELP,VR,DVR,VP,DVP,T)
    IMPLICIT NONE
    TYPE(ELEMENT),TARGET, INTENT(INOUT) ::EL
    TYPE(ELEMENTp),TARGET, INTENT(INOUT) ::ELp
    real(dp),OPTIONAL :: VR
    real(dp),OPTIONAL :: DVR,T
    TYPE(REAL_8),OPTIONAL :: VP
    TYPE(REAL_8),OPTIONAL :: DVP
    INTEGER N
    !     DV=(XS%AC(n)%X(1)*COS(EL%theta_ac)-XS%AC(n)%X(2)*SIN(EL%theta_ac))
    !     V=EL%DC_ac+EL%A_ac*DV
    !     DV=el%D_ac*DV
    !CALL transfer_ANBN(EL,ELP,VR=V,DVR=DV)

    if(EL%KIND==kind1) return

    if(associated(EL%ramp)) then
    
      if(EL%KIND/=kind15) then
          do n=1,EL%P%NMUL
             EL%BN(N)= EL%ramp%table(0)%bn(n)
             EL%AN(N)= EL%ramp%table(0)%an(n)
             ELP%BN(N)= ELP%ramp%table(0)%bn(n)
             ELP%AN(N)= ELP%ramp%table(0)%an(n)
          enddo  
      else
            EL%VOLT=EL%ramp%table(0)%bn(1)*COS(twopi*EL%ramp%table(0)%an(1)*T/clight+EL%ramp%table(0)%bn(2))+EL%ramp%table(0)%an(2)
           ELP%VOLT=EL%ramp%table(0)%bn(1)*COS(twopi*EL%ramp%table(0)%an(1)*T/clight+EL%ramp%table(0)%bn(2))+EL%ramp%table(0)%an(2)
         write(6,*) " volt ",el%volt,EL%ramp%table(0)%bn(1)
      endif
      
      if(EL%ramp%table(0)%b_t/=0.0_dp) then
          if(EL%parent_fibre%PATCH%TIME==0) EL%parent_fibre%PATCH%TIME=2
          if(EL%parent_fibre%PATCH%TIME==1) EL%parent_fibre%PATCH%TIME=3
          EL%parent_fibre%PATCH%b_T=EL%ramp%table(0)%b_t
        else
          if(EL%parent_fibre%PATCH%TIME==2) EL%parent_fibre%PATCH%TIME=0
          if(EL%parent_fibre%PATCH%TIME==3) EL%parent_fibre%PATCH%TIME=1
        EL%parent_fibre%PATCH%b_T=0.0_dp
      endif
          
    else

      IF(EL%P%NMUL>=1) THEN
        if(present(VR))then
          do n=1,EL%P%NMUL
             EL%BN(N)= vR*EL%D0_BN(N)+DVR*EL%D_BN(N) 
             EL%AN(N)= vR*EL%D0_AN(N)+DVR*EL%D_AN(N)
             ELP%BN(N)= vR*EL%D0_BN(N)+DVR*EL%D_BN(N)
             ELP%AN(N)= vR*EL%D0_AN(N)+DVR*EL%D_AN(N)
          enddo
        else
          do n=1,EL%P%NMUL
             EL%BN(N)= vp*EL%D0_BN(N)+DVp*EL%D_BN(N)
             EL%AN(N)= vp*EL%D0_AN(N)+DVp*EL%D_AN(N)
             ELP%BN(N)= vp*EL%D0_BN(N)+DVp*EL%D_BN(N)
             ELP%AN(N)= vp*EL%D0_AN(N)+DVp*EL%D_AN(N)
          enddo
       endif
    
    
     endif 
   endif
       if(el%kind==kind10) then
          call GETANBN(EL%TP10)
          call GETANBN(ELP%TP10)
       endif
       if(el%kind==kind7) then
          call GETMAT7(EL%T7)
          call GETMAT7(ELP%T7)
       endif

  END SUBROUTINE transfer_ANBN

  SUBROUTINE restore_ANBN(R)
    IMPLICIT NONE
    type(layout), target :: R
    type(fibre), pointer :: P
    INTEGER N


    p=>r%start

    do N=1,R%N
       IF(P%MAG%SLOW_AC/=0) THEN
          CALL restore_ANBN_SINGLE(P%MAG,P%MAGP)
       ELSE
          CYCLE
       ENDIF
       P=>P%NEXT
    ENDDO

  END SUBROUTINE restore_ANBN

  SUBROUTINE restore_ANBN_SINGLE(EL,ELP)
    IMPLICIT NONE
    TYPE(ELEMENT),TARGET, INTENT(INOUT) ::EL
    TYPE(ELEMENTp),TARGET, INTENT(INOUT) ::ELp
    INTEGER N

    IF(EL%P%NMUL>=1) THEN
       do n=1,EL%P%NMUL
          if(restore_mag) then
             EL%BN(N)= EL%D0_BN(N)
             EL%AN(N)= EL%D0_AN(N)
          endif
          if(restore_magp) then
             ELp%BN(N)= EL%D0_BN(N)
             ELp%AN(N)= EL%D0_AN(N)
          endif
       enddo
       if(el%kind==kind10) then
          if(restore_mag)call GETANBN(EL%TP10)
          if(restore_magp)call GETANBN(ELp%TP10)
       endif
       if(el%kind==kind7) then
          if(restore_mag) call GETMAT7(EL%T7)
          if(restore_magp) call GETMAT7(ELp%T7)
       endif
    ENDIF

  END SUBROUTINE restore_ANBN_SINGLE

  SUBROUTINE force_restore_ANBN_SINGLE(EL,ELP)
    IMPLICIT NONE
    TYPE(ELEMENT),TARGET, INTENT(INOUT) ::EL
    TYPE(ELEMENTp),TARGET, INTENT(INOUT) ::ELp
    logical(lp) rm,rmp

    rm=restore_mag
    rmp=restore_magp
    restore_mag=my_true
    restore_magp=my_true

    call restore_ANBN_SINGLE(EL,ELP)

    restore_mag=rm
    restore_magp=rmp

  END SUBROUTINE force_restore_ANBN_SINGLE

  SUBROUTINE force_restore_ANBN(R)
    IMPLICIT NONE
    type(layout), target :: R
    type(fibre), pointer :: P
    INTEGER N


    p=>r%start

    do N=1,R%N
       IF(P%MAG%SLOW_AC/=0) CALL force_restore_ANBN_SINGLE(P%MAG,P%MAGP)
       P=>P%NEXT
    ENDDO

  END SUBROUTINE force_restore_ANBN


  SUBROUTINE ADD_ANBNR(EL,NM,F,V,electric)
    IMPLICIT NONE
    TYPE(ELEMENT), INTENT(INOUT) ::EL
    real(dp), INTENT(IN) ::V
    INTEGER, INTENT(IN) ::NM,F
    INTEGER I,N
    real(dp), ALLOCATABLE,dimension(:)::AN,BN
    logical(lp), optional :: electric
    logical(lp) elec
    elec=my_false
    if(present(electric)) elec=electric
    if(elec.and.(.not.EL%KIND==kind10)) return

if(elec) then
    N=NM
    IF(NM<0) N=-N
    if(N>SECTOR_NMUL_max) THEN
     WRITE(6,*) " ADD_ANBNR NOT PERMITTED N>SECTOR_NMUL  " ,N,SECTOR_NMUL_max
     STOP
    ENDIF
    ! ALREADY THERE
       IF(NM>0) THEN
          EL%TP10%BE(N)= F*EL%TP10%BE(N)+V*volt_i
       ELSE
          EL%TP10%AE(N)= F*EL%TP10%AE(N)+V*volt_i
       ENDIF
       if(el%kind==kind10) then
          call GETAEBE(EL%TP10)
       endif
else

    if(EL%KIND==kind1) return
    N=NM
    IF(NM<0) N=-N
    ! ALREADY THERE
    IF(EL%P%NMUL>=N) THEN
       IF(NM>0) THEN
          EL%BN(N)= F*EL%BN(N)+V
       ELSE
          EL%AN(N)= F*EL%AN(N)+V
       ENDIF
       if(el%kind==kind10) then
        if(el%electric) then
        call GETAEBE(EL%TP10)
        else
         call GETANBN(EL%TP10)
        endif
       endif
       if(el%kind==kind7) then
          call GETMAT7(EL%T7)
       endif
       RETURN
    ENDIF

    allocate(AN(N),BN(N))
    DO I=1,EL%P%NMUL
       AN(I)=EL%AN(I)
       BN(I)=EL%BN(I)
    ENDDO
    DO I=EL%P%NMUL+1,N
       AN(I)=0.0_dp
       BN(I)=0.0_dp
    ENDDO
    IF(NM<0) THEN
       AN(N)=V
    ELSE
       BN(N)=V
    ENDIF


    IF(ASSOCIATED(EL%AN)) DEALLOCATE(EL%AN)
    IF(ASSOCIATED(EL%BN)) DEALLOCATE(EL%BN)
    EL%P%NMUL=N
    ALLOCATE(EL%AN(EL%P%NMUL),EL%BN(EL%P%NMUL))

    DO I=1,EL%P%NMUL
       EL%AN(I)=AN(I)
       EL%BN(I)=BN(I)
    ENDDO

    DEALLOCATE(AN);DEALLOCATE(BN);

    SELECT CASE(EL%KIND)
       !    CASE(KIND2,KIND3,KIND5,KIND6,KIND17)
       !       select case(EL%KIND)
    case(kind2)
       EL%K2%AN=>EL%AN
       EL%K2%BN=>EL%BN
    case(kind3)
       EL%K3%AN=>EL%AN
       EL%K3%BN=>EL%BN
    case(kind4)
       EL%C4%AN=>EL%AN
       EL%C4%BN=>EL%BN
    case(kind5)
       EL%S5%AN=>EL%AN
       EL%S5%BN=>EL%BN
    case(kind6)
       EL%T6%AN=>EL%AN
       EL%T6%BN=>EL%BN
    CASE(KIND7)
       EL%T7%AN=>EL%AN
       EL%T7%BN=>EL%BN
       call GETMAT7(EL%T7)
    CASE(KIND8)
       EL%S8%BN=>EL%BN
    CASE(KIND9)
       EL%S9%AN=>EL%AN
    CASE(KIND10)
       EL%TP10%AN=>EL%AN
       EL%TP10%BN=>EL%BN
       if(el%electric) then
        call GETAEBE(EL%TP10)
       else
        call GETANBN(EL%TP10)
       endif
    CASE(KIND16,KIND20)
       EL%K16%AN=>EL%AN
       EL%K16%BN=>EL%BN
       !    CASE(KINDuser1)
       !       EL%U1%AN=>EL%AN
       !       EL%U1%BN=>EL%BN
       !    CASE(KINDuser2)
       !       EL%U2%AN=>EL%AN
       !       EL%U2%BN=>EL%BN
    CASE(kind17)
       EL%ENGE17%AN=>EL%AN
       EL%ENGE17%BN=>EL%BN
    CASE(KINDWIGGLER)
       EL%WI%AN=>EL%AN
       EL%WI%BN=>EL%BN
    case(kind22)
       EL%HE22%AN=>EL%AN
       EL%HE22%BN=>EL%BN
    case default
       !w_p=0
       !w_p%nc=1
       !w_p%fc='((1X,A72,/,1X,A72))'
       write(*,'(A13,A24,A27)')" THIS MAGNET ", MYTYPE(EL%KIND), " CANNOT ACCEPT ANs AND BNs "
       ! call !write_e(988)
    END SELECT
endif

    !    if(el%kind==kind10) then
    !    call GETANBN(EL%TP10)
    !    endif
    !    if(el%kind==kind7) then
    !       call GETMAT7(EL%T7)
    !    endif

  END SUBROUTINE ADD_ANBNR

  SUBROUTINE ADD_ANBNP(EL,NM,F,V,electric)
    IMPLICIT NONE
    TYPE(ELEMENTP), INTENT(INOUT) ::EL
    real(dp), INTENT(IN) ::V
    INTEGER, INTENT(IN) ::NM,F
    INTEGER I,N
    TYPE(REAL_8), ALLOCATABLE,dimension(:)::AN,BN
    logical(lp), optional :: electric
    logical(lp) elec
    elec=my_false
    if(present(electric)) elec=electric
    if(elec.and.(.not.EL%KIND==kind10)) return
if(elec) then
    N=NM
    IF(NM<0) N=-N
    if(N>SECTOR_NMUL_max) THEN
     WRITE(6,*) " ADD_ANBNP NOT PERMITTED N>SECTOR_NMUL  " ,N,SECTOR_NMUL_max
     STOP
    ENDIF
    ! ALREADY THERE
       IF(NM>0) THEN
          EL%TP10%BE(N)= F*EL%TP10%BE(N)+V
       ELSE
          EL%TP10%AE(N)= F*EL%TP10%AE(N)+V
       ENDIF
       if(el%kind==kind10) then
          call GETAEBE(EL%TP10)
       endif
else

    if(EL%KIND==kind1) return

    N=NM
    IF(NM<0) N=-N
    ! ALREADY THERE
    IF(EL%P%NMUL>=N) THEN
       IF(NM>0) THEN
          EL%BN(N)= F*EL%BN(N)+V*volt_i
       ELSE
          EL%AN(N)= F*EL%AN(N)+V*volt_i
       ENDIF
       if(el%kind==kind10) then
        if(el%electric) then
        call GETAEBE(EL%TP10)
        else
         call GETANBN(EL%TP10)
        endif
       endif
       if(el%kind==kind7) then
          call GETMAT7(EL%T7)     !etienne
       endif
       RETURN
    ENDIF

    allocate(AN(N),BN(N))
    CALL ALLOC(AN,N);CALL ALLOC(BN,N);
    DO I=1,EL%P%NMUL
       AN(I)=EL%AN(I)
       BN(I)=EL%BN(I)
    ENDDO
    IF(NM<0) THEN
       AN(N)=V
    ELSE
       BN(N)=V
    ENDIF

    CALL KILL(EL%AN,EL%P%NMUL);CALL KILL(EL%BN,EL%P%NMUL);
    IF(ASSOCIATED(EL%AN)) DEALLOCATE(EL%AN)
    IF(ASSOCIATED(EL%BN)) DEALLOCATE(EL%BN)
    EL%P%NMUL=N
    ALLOCATE(EL%AN(EL%P%NMUL),EL%BN(EL%P%NMUL))
    CALL ALLOC(EL%AN,EL%P%NMUL);CALL ALLOC(EL%BN,EL%P%NMUL);  ! BUG there

    DO I=1,EL%P%NMUL
       EL%AN(I)=AN(I)
       EL%BN(I)=BN(I)
    ENDDO

    DEALLOCATE(AN);DEALLOCATE(BN);

    SELECT CASE(EL%KIND)
       !   CASE(KIND2,KIND3,KIND5,KIND6,KIND17)
       !      select case(EL%KIND)
    case(kind2)
       EL%K2%AN=>EL%AN
       EL%K2%BN=>EL%BN
    case(kind3)
       EL%K3%AN=>EL%AN
       EL%K3%BN=>EL%BN
    case(kind4)
       EL%C4%AN=>EL%AN
       EL%C4%BN=>EL%BN
    case(kind5)
       EL%S5%AN=>EL%AN
       EL%S5%BN=>EL%BN
    case(kind6)
       EL%T6%AN=>EL%AN
       EL%T6%BN=>EL%BN
    CASE(KIND7)
       EL%T7%AN=>EL%AN
       EL%T7%BN=>EL%BN
       call GETMAT7(EL%T7)
    CASE(KIND8)
       EL%S8%BN=>EL%BN
    CASE(KIND9)
       EL%S9%AN=>EL%AN
    CASE(KIND10)
       EL%TP10%AN=>EL%AN
       EL%TP10%BN=>EL%BN
        if(el%electric) then
        call GETAEBE(EL%TP10)
        else
         call GETANBN(EL%TP10)
        endif
    CASE(KIND16,KIND20)
       EL%K16%AN=>EL%AN
       EL%K16%BN=>EL%BN
       !    CASE(KINDuser1)
       !       EL%U1%AN=>EL%AN
       !       EL%U1%BN=>EL%BN
       !    CASE(KINDuser2)
       !       EL%U2%AN=>EL%AN
       !       EL%U2%BN=>EL%BN
    CASE(kind17)
       EL%ENGE17%AN=>EL%AN
       EL%ENGE17%BN=>EL%BN
    case(kind22)
       EL%HE22%AN=>EL%AN
       EL%HE22%BN=>EL%BN
    CASE(KINDWIGGLER)
       EL%WI%AN=>EL%AN
       EL%WI%BN=>EL%BN
    case default
       !w_p=0
       !w_p%nc=1
       !w_p%fc='((1X,A72,/,1X,A72))'
       write(*,'(A13,A24,A27)')" THIS MAGNET ", MYTYPE(EL%KIND), " CANNOT ACCEPT ANs AND BNs "
       ! call !write_e(987)
    END SELECT
ENDIF
    !if(el%kind==kind10) then
    !call GETANBN(EL%TP10)
    !endif
    !if(el%kind==kind7) then
    !   call GETMAT7(EL%T7)
    !endif

  END SUBROUTINE ADD_ANBNP



  SUBROUTINE null_EL(EL)
    IMPLICIT NONE
    TYPE(ELEMENT), INTENT(INOUT)::EL
    nullify(EL%KIND);
    nullify(EL%PLOT);
    nullify(EL%NAME);nullify(EL%vorname);nullify(EL%electric);
nullify(EL%filef,el%fileb);
!    nullify(EL%PERMFRINGE);
    nullify(EL%L);
    nullify(EL%AN);nullify(EL%BN);
    nullify(EL%FINT);nullify(EL%HGAP);
    nullify(EL%H1);nullify(EL%H2);
    nullify(EL%VA);nullify(EL%VS);nullify(EL%ene);
    nullify(EL%VOLT);nullify(EL%FREQ);nullify(EL%PHAS);nullify(EL%DELTA_E);
    nullify(EL%lag);
    nullify(EL%B_SOL);
    nullify(EL%slow_ac);
    nullify(EL%a_ac);
    nullify(EL%theta_ac);
    nullify(EL%DC_ac);
    nullify(EL%D_AC);nullify(EL%D_AN);nullify(EL%D_BN);nullify(EL%D0_AN);nullify(EL%D0_BN);
    nullify(EL%THIN);
    nullify(EL%MIS); !nullify(EL%EXACTMIS);
    !    nullify(EL%D);nullify(EL%R);
    nullify(EL%D0);
    nullify(EL%K2);
    nullify(EL%K16);
    nullify(EL%K3);
    nullify(EL%C4);
    nullify(EL%CAV21);
    nullify(EL%HE22);
    nullify(EL%S5);
    nullify(EL%T6);
    !    nullify(EL%M22);
    nullify(EL%T7);
    nullify(EL%S8);
    nullify(EL%S9);
    nullify(EL%TP10);
    nullify(EL%MON14);
    nullify(EL%SEP15);
    nullify(EL%RCOL18);
    nullify(EL%ECOL19);
    !    nullify(EL%U1);
    !    nullify(EL%U2);
    nullify(EL%WI);
    nullify(EL%RAMP);
    nullify(EL%PA);
    nullify(EL%P);
    nullify(EL%siamese);
    nullify(EL%girders);
    nullify(EL%assembly);
    nullify(EL%SIAMESE_FRAME);
    nullify(EL%girder_FRAME);
    nullify(EL%doko);
    nullify(EL%forward,EL%backWARD,el%usef,el%useb,el%skip_ptc_f,el%skip_ptc_b,el%do1mapf,el%do1mapb);   
  end SUBROUTINE null_EL

  SUBROUTINE null_ELp(EL)
    IMPLICIT NONE
    TYPE(ELEMENTP), INTENT(INOUT)::EL

    nullify(EL%KNOB);
    nullify(EL%KIND);
    nullify(EL%NAME);nullify(EL%vorname);nullify(EL%electric);

!    nullify(EL%PERMFRINGE);
    nullify(EL%L);
    nullify(EL%AN);nullify(EL%BN);
    nullify(EL%FINT);nullify(EL%HGAP);
    nullify(EL%H1);nullify(EL%H2);
    nullify(EL%VA);nullify(EL%VS);
    nullify(EL%VOLT);nullify(EL%FREQ);nullify(EL%PHAS);nullify(EL%DELTA_E);
    nullify(EL%B_SOL);
    nullify(EL%slow_ac);
    nullify(EL%a_ac);
    nullify(EL%theta_ac);
    nullify(EL%DC_ac);
    nullify(EL%D_AC);nullify(EL%D_AN);nullify(EL%D_BN);nullify(EL%D0_AN);nullify(EL%D0_BN);
    nullify(EL%THIN);
    nullify(EL%MIS);  !nullify(EL%EXACTMIS);
    !    nullify(EL%D);nullify(EL%R);
    nullify(EL%D0);
    nullify(EL%K2);
    nullify(EL%K16);
    nullify(EL%K3);
    nullify(EL%C4);
    nullify(EL%CAV21);
    nullify(EL%HE22);
    nullify(EL%S5);
    nullify(EL%T6);
    !    nullify(EL%M22);
    nullify(EL%T7);
    nullify(EL%S8);
    nullify(EL%S9);
    nullify(EL%TP10);
    nullify(EL%MON14);
    nullify(EL%SEP15);
    nullify(EL%RCOL18);
    nullify(EL%ECOL19);
    !    nullify(EL%U1);
    !    nullify(EL%U2);
    nullify(EL%WI);
    nullify(EL%RAMP);
    nullify(EL%PA);
    nullify(EL%P);
    nullify(EL%PARENT_FIBRE);
    nullify(EL%forward,EL%backWARD,el%usef,el%useb,el%skip_ptc_f,el%skip_ptc_b,el%do1mapf,el%do1mapb);
  end SUBROUTINE null_ELp



  SUBROUTINE ZERO_EL(EL,I)
    IMPLICIT NONE
    TYPE(ELEMENT),target, INTENT(INOUT)::EL
    INTEGER, INTENT(IN)::I

    IF(I==-1) THEN
       DEALLOCATE(EL%KIND);
       DEALLOCATE(EL%PLOT);
       DEALLOCATE(EL%recut);
       DEALLOCATE(EL%even);
       DEALLOCATE(EL%NAME);DEALLOCATE(EL%VORNAME);DEALLOCATE(EL%electric);
       DEALLOCATE(EL%L);
        DEALLOCATE(EL%filef,el%fileb);
       DEALLOCATE(EL%MIS); !DEALLOCATE(EL%EXACTMIS);
       call kill(EL%P)    ! AIMIN MS 4.0
!       DEALLOCATE(EL%PERMFRINGE);
       !       IF(ASSOCIATED(EL%R)) DEALLOCATE(EL%R)
       !       IF(ASSOCIATED(EL%D)) DEALLOCATE(EL%D)
       IF(ASSOCIATED(EL%AN)) DEALLOCATE(EL%AN)
       IF(ASSOCIATED(EL%BN)) DEALLOCATE(EL%BN)
       IF(ASSOCIATED(EL%FINT)) DEALLOCATE(EL%FINT)
       IF(ASSOCIATED(EL%HGAP)) DEALLOCATE(EL%HGAP)
       IF(ASSOCIATED(EL%H1)) DEALLOCATE(EL%H1)
       IF(ASSOCIATED(EL%H2)) DEALLOCATE(EL%H2)
       IF(ASSOCIATED(EL%VA)) DEALLOCATE(EL%VA)
       IF(ASSOCIATED(EL%VS)) DEALLOCATE(EL%VS)
              IF(ASSOCIATED(EL%ene)) DEALLOCATE(EL%ene)
       IF(ASSOCIATED(EL%VOLT)) DEALLOCATE(EL%VOLT)
       IF(ASSOCIATED(EL%lag)) DEALLOCATE(EL%lag)
       IF(ASSOCIATED(EL%FREQ)) DEALLOCATE(EL%FREQ)
       IF(ASSOCIATED(EL%PHAS)) DEALLOCATE(EL%PHAS)
       IF(ASSOCIATED(EL%DELTA_E)) DEALLOCATE(EL%DELTA_E)
       IF(ASSOCIATED(EL%B_SOL)) DEALLOCATE(EL%B_SOL)
       IF(ASSOCIATED(EL%slow_ac)) DEALLOCATE(EL%slow_ac)
       IF(ASSOCIATED(EL%a_ac)) DEALLOCATE(EL%a_ac)
       IF(ASSOCIATED(EL%theta_ac)) DEALLOCATE(EL%theta_ac)
       IF(ASSOCIATED(EL%DC_ac)) DEALLOCATE(EL%DC_ac)
       IF(ASSOCIATED(EL%D_AC)) DEALLOCATE(EL%D_AC)
       IF(ASSOCIATED(EL%D_AN)) DEALLOCATE(EL%D_AN)
       IF(ASSOCIATED(EL%D_BN)) DEALLOCATE(EL%D_BN)
       IF(ASSOCIATED(EL%D0_AN)) DEALLOCATE(EL%D0_AN)
       IF(ASSOCIATED(EL%D0_BN)) DEALLOCATE(EL%D0_BN)
       IF(ASSOCIATED(EL%THIN)) DEALLOCATE(EL%THIN)
       IF(ASSOCIATED(EL%d0)) DEALLOCATE(EL%d0)       ! drift
!       IF(ASSOCIATED(EL%K2)) DEALLOCATE(EL%K2)       ! INTEGRATOR
       IF(ASSOCIATED(EL%K2)) then
          EL%K2=-1
        DEALLOCATE(EL%K2)       ! SOLENOID
       endif
       !       IF(ASSOCIATED(EL%K16)) DEALLOCATE(EL%K16)       ! INTEGRATOR
       !       IF(ASSOCIATED(EL%K3)) DEALLOCATE(EL%K3)       !  THIN LENS
       IF(ASSOCIATED(EL%K3)) then
          !          IF(ASSOCIATED(EL%K3%hf)) DEALLOCATE(EL%K3%hf)
          !          IF(ASSOCIATED(EL%K3%vf)) DEALLOCATE(EL%K3%vf)
          !          IF(ASSOCIATED(EL%K3%thin_h_foc)) DEALLOCATE(EL%K3%thin_h_foc)
          !          IF(ASSOCIATED(EL%K3%thin_v_foc)) DEALLOCATE(EL%K3%thin_v_foc)
          !          IF(ASSOCIATED(EL%K3%thin_h_angle)) DEALLOCATE(EL%K3%thin_h_angle)
          !          IF(ASSOCIATED(EL%K3%thin_v_angle)) DEALLOCATE(EL%K3%thin_v_angle)
          !          IF(ASSOCIATED(EL%K3%patch)) DEALLOCATE(EL%K3%patch)
          EL%K3=-1
          DEALLOCATE(EL%K3)
       endif

       IF(ASSOCIATED(EL%S5)) then
          EL%S5=-1
        DEALLOCATE(EL%S5)       ! SOLENOID
       endif
       !       IF(ASSOCIATED(EL%T6)) DEALLOCATE(EL%T6)       ! INTEGRATOR
       !       IF(ASSOCIATED(EL%T7)) DEALLOCATE(EL%T7)       ! INTEGRATOR
       IF(ASSOCIATED(EL%S8)) DEALLOCATE(EL%S8)       ! NORMAL SMI
       IF(ASSOCIATED(EL%S9)) DEALLOCATE(EL%S9)       ! SKEW SMI
       !       IF(ASSOCIATED(EL%TP10)) DEALLOCATE(EL%TP10)   ! SECTOR TEAPOT
       IF(ASSOCIATED(EL%T6)) THEN
          EL%T6=-1
          DEALLOCATE(EL%T6)   ! thick sixtrack
       ENDIF
       !       IF(ASSOCIATED(EL%M22)) THEN
       !          EL%M22=-1
       !          DEALLOCATE(EL%M22)   ! thick sixtrack
       !       ENDIF
       IF(ASSOCIATED(EL%T7)) THEN
          EL%T7=-1
          DEALLOCATE(EL%T7)   ! thick
       ENDIF
       IF(ASSOCIATED(EL%C4)) THEN
          EL%C4=-1
          DEALLOCATE(EL%C4)   ! MONITOR
       ENDIF
       IF(ASSOCIATED(EL%CAV21)) THEN
          EL%CAV21=-1
          DEALLOCATE(EL%CAV21)   ! MONITOR
       ENDIF
       IF(ASSOCIATED(EL%HE22)) THEN
          EL%HE22=-1
          DEALLOCATE(EL%HE22)   ! MONITOR
       ENDIF
       IF(ASSOCIATED(EL%TP10)) then
          EL%TP10=-1
          DEALLOCATE(EL%TP10)   ! SECTOR TEAPOT
       ENDIF
       IF(ASSOCIATED(EL%MON14)) THEN
          EL%MON14=-1
          DEALLOCATE(EL%MON14)   ! MONITOR
       ENDIF
       IF(ASSOCIATED(EL%RCOL18)) THEN
          EL%RCOL18=-1
          DEALLOCATE(EL%RCOL18)   ! RCOLLIMATOR
       ENDIF
       IF(ASSOCIATED(EL%ECOL19)) THEN
          EL%ECOL19=-1
          DEALLOCATE(EL%ECOL19)   ! ECOLLIMATOR
       ENDIF
       IF(ASSOCIATED(EL%SEP15)) DEALLOCATE(EL%SEP15)       ! ELSEPARATOR
       IF(ASSOCIATED(EL%K16)) then
          EL%K16=-1
          DEALLOCATE(EL%K16)       ! INTEGRATOR
       endif
       !       IF(ASSOCIATED(EL%U1))        then
       !          el%U1=-1     !USER DEFINED MAGNET
       !          DEALLOCATE(EL%U1)
       !       ENDIF

       !       IF(ASSOCIATED(EL%U2))        then
       !          el%U2=-1     !USER DEFINED MAGNET
       !          DEALLOCATE(EL%U2)
       !       ENDIF

       IF(ASSOCIATED(EL%WI))        then
          el%WI=-1     !USER DEFINED MAGNET
          DEALLOCATE(EL%WI)
       ENDIF

       IF(ASSOCIATED(EL%forward))        then
          call kill(EL%forward)     
          DEALLOCATE(EL%forward)
       ENDIF


       IF(ASSOCIATED(EL%backWARD))        then
          call kill(EL%backWARD)     
          DEALLOCATE(EL%backWARD)
       ENDIF

    IF(ASSOCIATED(EL%skip_ptc_f))DEALLOCATE(EL%skip_ptc_f)
    IF(ASSOCIATED(EL%skip_ptc_b))DEALLOCATE(EL%skip_ptc_b)
    IF(ASSOCIATED(el%do1mapf))DEALLOCATE(el%do1mapf)
    IF(ASSOCIATED(el%do1mapb))DEALLOCATE(el%do1mapb)
    IF(ASSOCIATED(el%usef))DEALLOCATE(el%usef)
    IF(ASSOCIATED(el%useb))DEALLOCATE(el%useb)

       IF(ASSOCIATED(EL%ramp))        then
          el%ramp=-1     !USER DEFINED MAGNET
          DEALLOCATE(EL%ramp)
       ENDIF

       IF(ASSOCIATED(EL%PARENT_FIBRE))        then
          nullify(EL%PARENT_FIBRE)
       ENDIF
       IF(ASSOCIATED(EL%DOKO))        then
          nullify(EL%DOKO)
       ENDIF
       nullify(EL%siamese);
       nullify(EL%girders);
       IF(ASSOCIATED(EL%SIAMESE_FRAME))        then
          call kill_af(EL%SIAMESE_FRAME)
          DEALLOCATE(EL%SIAMESE_FRAME)
       ENDIF
       IF(ASSOCIATED(EL%girder_FRAME))        then
          call kill_af(EL%girder_FRAME)
          DEALLOCATE(EL%girder_FRAME)
       ENDIF


    elseif(I>=0)       then

       !FIRST nullifies

       call null_ELEment(el)

       call alloc(el%P);

       ALLOCATE(EL%KIND);EL%KIND=0;
       ALLOCATE(EL%PLOT);EL%PLOT=MY_TRUE;
       ALLOCATE(EL%RECUT);EL%RECUT=MY_TRUE;
       ALLOCATE(EL%even);EL%even=MY_false;
       ALLOCATE(EL%NAME);ALLOCATE(EL%VORNAME);ALLOCATE(EL%electric);
       ALLOCATE(EL%filef,el%fileb);
       ALLOCATE(EL%skip_ptc_f);   EL%skip_ptc_f=0;   ALLOCATE(EL%skip_ptc_b);EL%skip_ptc_b=0  ;
       ALLOCATE(el%do1mapf);   el%do1mapf=.false. ;   ALLOCATE(el%do1mapb);el%do1mapb=.false.  ;
  ALLOCATE(el%usef);   el%usef=.false. ;   ALLOCATE(el%useb);el%useb=.false.  ;

       EL%NAME=' ';EL%NAME=TRIM(ADJUSTL(EL%NAME));
       EL%VORNAME=' ';EL%VORNAME=TRIM(ADJUSTL(EL%VORNAME));el%filef=' ';el%fileb=' ';
       EL%electric=solve_electric
!       ALLOCATE(EL%PERMFRINGE);EL%PERMFRINGE=.FALSE.;  ! PART OF A STATE INITIALIZED BY EL=DEFAULT
       ALLOCATE(EL%L);EL%L=0.0_dp;
       ALLOCATE(EL%MIS);
       !       ALLOCATE(EL%girder_index);
       !       ALLOCATE(EL%EXACTMIS);
       EL%MIS=.FALSE.;
       !       EL%EXACTMIS=ALWAYS_EXACTMIS;
       !       allocate(el%r(3));allocate(el%d(3));
       !      el%r=zero;el%d=zero;

       !       EL=DEFAULT;
       !   ANBN
       CALL ZERO_ANBN(EL,I)
       ALLOCATE(EL%FINT(2));EL%FINT(1)=0.5_dp;EL%FINT(2)=0.5_dp;
       ALLOCATE(EL%HGAP(2));EL%HGAP(1)=0.0_dp;EL%HGAP(2)=0.0_dp;
       ALLOCATE(EL%H1);EL%H1=0.0_dp;
       ALLOCATE(EL%H2);EL%H2=0.0_dp;
       ALLOCATE(EL%VA);EL%VA=0.0_dp;
       ALLOCATE(EL%VS);EL%VS=0.0_dp;
              ALLOCATE(EL%ene);EL%ene=0.0_dp;
       !       ALLOCATE(EL%theta_ac); EL%theta_ac= zero ;
       !       ALLOCATE(EL%a_ac);  EL%a_ac = zero;
       !       ALLOCATE(EL%DC_ac); EL%DC_ac= zero ;
       ALLOCATE(EL%slow_ac); EL%slow_ac=0 ;
    ENDIF

  END SUBROUTINE ZERO_EL

  SUBROUTINE ZERO_ELP(EL,I)
    IMPLICIT NONE
    TYPE(ELEMENTP),target, INTENT(INOUT)::EL
    INTEGER, INTENT(IN)::I
    INTEGER J

    IF(I==-1) THEN

       IF(ASSOCIATED(EL%P%NMUL))THEN
          IF(EL%P%NMUL>0) THEN
             DO  J=1,EL%P%NMUL
                CALL KILL(EL%AN(J))
                CALL KILL(EL%BN(J))
             ENDDO
             IF(ASSOCIATED(EL%AN)) DEALLOCATE(EL%AN)
             IF(ASSOCIATED(EL%BN)) DEALLOCATE(EL%BN)
          ENDIF
       ENDIF

       IF(ASSOCIATED(EL%d0)) DEALLOCATE(EL%d0)       ! drift
!       IF(ASSOCIATED(EL%K2)) DEALLOCATE(EL%K2)       ! INTEGRATOR
       IF(ASSOCIATED(EL%K2)) then
          EL%K2=-1
        DEALLOCATE(EL%K2)       ! SOLENOID
       endif
       !       IF(ASSOCIATED(EL%K16)) DEALLOCATE(EL%K16)       ! INTEGRATOR
       !       IF(ASSOCIATED(EL%K3)) DEALLOCATE(EL%K3)       !  THIN LENS
       IF(ASSOCIATED(EL%K3)) then
          EL%K3=-1
          DEALLOCATE(EL%K3)
          !          IF(ASSOCIATED(EL%K3%hf)) DEALLOCATE(EL%K3%hf)
          !          IF(ASSOCIATED(EL%K3%vf)) DEALLOCATE(EL%K3%vf)
          !          IF(ASSOCIATED(EL%K3%thin_h_foc)) DEALLOCATE(EL%K3%thin_h_foc)
          !          IF(ASSOCIATED(EL%K3%thin_v_foc)) DEALLOCATE(EL%K3%thin_v_foc)
          !          IF(ASSOCIATED(EL%K3%thin_h_angle)) DEALLOCATE(EL%K3%thin_h_angle)
          !          IF(ASSOCIATED(EL%K3%thin_v_angle)) DEALLOCATE(EL%K3%thin_v_angle)
          !          IF(ASSOCIATED(EL%K3%patch)) DEALLOCATE(EL%K3%patch)
       endif

       IF(ASSOCIATED(EL%C4)) THEN
          EL%C4=-1
          DEALLOCATE(EL%C4)       ! CAVITY
          CALL KILL(EL%VOLT)
          CALL KILL(EL%FREQ)
          CALL KILL(EL%PHAS)
          IF(ASSOCIATED(EL%VOLT)) DEALLOCATE(EL%VOLT)
          IF(ASSOCIATED(EL%FREQ)) DEALLOCATE(EL%FREQ)
          IF(ASSOCIATED(EL%PHAS)) DEALLOCATE(EL%PHAS)
          IF(ASSOCIATED(EL%DELTA_E)) DEALLOCATE(EL%DELTA_E)
          IF(ASSOCIATED(EL%THIN)) DEALLOCATE(EL%THIN)
       ENDIF
       IF(ASSOCIATED(EL%CAV21)) THEN
          EL%CAV21=-1
          DEALLOCATE(EL%CAV21)       ! CAVITY
          CALL KILL(EL%VOLT)
          CALL KILL(EL%FREQ)
          CALL KILL(EL%PHAS)
          IF(ASSOCIATED(EL%VOLT)) DEALLOCATE(EL%VOLT)
          IF(ASSOCIATED(EL%FREQ)) DEALLOCATE(EL%FREQ)
          IF(ASSOCIATED(EL%PHAS)) DEALLOCATE(EL%PHAS)
          IF(ASSOCIATED(EL%DELTA_E)) DEALLOCATE(EL%DELTA_E)
          IF(ASSOCIATED(EL%THIN)) DEALLOCATE(EL%THIN)
       ENDIF

       IF(ASSOCIATED(EL%HE22)) THEN
          EL%HE22=-1
          DEALLOCATE(EL%HE22)       ! CAVITY
          CALL KILL(EL%FREQ)
          CALL KILL(EL%PHAS)
          IF(ASSOCIATED(EL%FREQ)) DEALLOCATE(EL%FREQ)
          IF(ASSOCIATED(EL%PHAS)) DEALLOCATE(EL%PHAS)
       ENDIF

       IF(ASSOCIATED(EL%S5)) THEN
          EL%S5=-1
          DEALLOCATE(EL%S5)       ! solenoid
          !          CALL KILL(EL%B_SOL)    ! sagan
          !         IF(ASSOCIATED(EL%B_SOL)) DEALLOCATE(EL%B_SOL)     ! sagan
       ENDIF
       IF(ASSOCIATED(EL%T6)) then
          EL%T6=-1
          DEALLOCATE(EL%T6)       ! INTEGRATOR
       endif
       !       IF(ASSOCIATED(EL%M22)) then
       !          EL%M22=-1
       !          DEALLOCATE(EL%M22)       ! INTEGRATOR
       !       endif
       IF(ASSOCIATED(EL%T7)) then
          EL%T7=-1
          DEALLOCATE(EL%T7)       ! INTEGRATOR
       ENDIF
       IF(ASSOCIATED(EL%S8)) DEALLOCATE(EL%S8)       ! SMI KICK
       IF(ASSOCIATED(EL%S9)) DEALLOCATE(EL%S9)       ! SKEW SMI KICK
       IF(ASSOCIATED(EL%MON14)) THEN
          EL%MON14=-1
          DEALLOCATE(EL%MON14)   ! MONITOR
       ENDIF
       IF(ASSOCIATED(EL%RCOL18)) THEN
          EL%RCOL18=-1
          DEALLOCATE(EL%RCOL18)   ! RCOLLIMATOR
       ENDIF
       IF(ASSOCIATED(EL%ECOL19)) THEN
          EL%ECOL19=-1
          DEALLOCATE(EL%ECOL19)   ! ECOLLIMATOR
       ENDIF
       IF(ASSOCIATED(EL%K16)) then
          EL%K16=-1
          DEALLOCATE(EL%K16)       ! INTEGRATOR
       endif
       IF(ASSOCIATED(EL%SEP15)) THEN
          DEALLOCATE(EL%SEP15)       ! CAVITY
          CALL KILL(EL%VOLT); CALL KILL(EL%PHAS);
          IF(ASSOCIATED(EL%VOLT)) DEALLOCATE(EL%VOLT)
          IF(ASSOCIATED(EL%PHAS)) DEALLOCATE(EL%PHAS)
       ENDIF
       IF(ASSOCIATED(EL%TP10)) then
          EL%TP10=-1
          DEALLOCATE(EL%TP10)       ! INTEGRATOR SECTOR EXACT
       ENDIF
       !       IF(ASSOCIATED(EL%U1))        then
       !          el%U1=-1
       !          DEALLOCATE(EL%U1)
       !       ENDIF
       !       IF(ASSOCIATED(EL%U2))        then
       !          el%U2=-1
       !          DEALLOCATE(EL%U2)
       !       ENDIF

       IF(ASSOCIATED(EL%PARENT_FIBRE))        then
          nullify(EL%PARENT_FIBRE)
       ENDIF
       IF(ASSOCIATED(EL%WI))        then
          el%WI=-1
          DEALLOCATE(EL%WI)
       ENDIF

       IF(ASSOCIATED(EL%forward))        then
          call kill(EL%forward)     
          DEALLOCATE(EL%forward)
       ENDIF

       IF(ASSOCIATED(EL%backWARD))        then
          call kill(EL%backWARD)     
          DEALLOCATE(EL%backWARD)
       ENDIF

    IF(ASSOCIATED(EL%skip_ptc_f))DEALLOCATE(EL%skip_ptc_f)
    IF(ASSOCIATED(EL%skip_ptc_b))DEALLOCATE(EL%skip_ptc_b)    
    IF(ASSOCIATED(el%do1mapb))DEALLOCATE(el%do1mapb)
    IF(ASSOCIATED(el%do1mapf))DEALLOCATE(el%do1mapf)
    IF(ASSOCIATED(el%usef))DEALLOCATE(el%usef)
    IF(ASSOCIATED(el%useb))DEALLOCATE(el%useb)

       IF(ASSOCIATED(EL%ramp))        then
          el%ramp=-1     !USER DEFINED MAGNET
          DEALLOCATE(EL%ramp)
       ENDIF


       !       IF(ASSOCIATED(EL%PARENT_FIBRE))        then
       !          nullify(EL%PARENT_FIBRE)
       !       ENDIF


       DEALLOCATE(EL%KIND);DEALLOCATE(EL%KNOB);
       DEALLOCATE(EL%NAME);DEALLOCATE(EL%VORNAME);DEALLOCATE(EL%electric);
!       DEALLOCATE(EL%PERMFRINGE);
       CALL KILL(EL%L);DEALLOCATE(EL%L);
       CALL KILL(EL%FINT);DEALLOCATE(EL%FINT);
       CALL KILL(EL%HGAP);DEALLOCATE(EL%HGAP);
       CALL KILL(EL%H1);DEALLOCATE(EL%H1);
       CALL KILL(EL%H2);DEALLOCATE(EL%H2);
       CALL KILL(EL%VA);DEALLOCATE(EL%VA);
       CALL KILL(EL%VS);DEALLOCATE(EL%VS);
       DEALLOCATE(EL%MIS); !DEALLOCATE(EL%EXACTMIS);

       IF(ASSOCIATED(EL%slow_ac))DEALLOCATE(EL%slow_ac)
       IF(ASSOCIATED(EL%a_ac)) then
          call kill(el%a_ac)
          DEALLOCATE(EL%a_ac)
       endif
       IF(ASSOCIATED(EL%theta_ac)) then
          call kill(el%theta_ac)
          DEALLOCATE(EL%theta_ac)
       endif
       IF(ASSOCIATED(EL%DC_ac)) then
          call kill(el%DC_ac)
          DEALLOCATE(EL%DC_ac)
       endif
       IF(ASSOCIATED(EL%D_AC)) then
          call kill(el%D_AC)
          DEALLOCATE(EL%D_AC)
       endif
       IF(ASSOCIATED(EL%D_AN)) then
          call kill(el%D_AN)
          DEALLOCATE(EL%D_AN)
       endif
       IF(ASSOCIATED(EL%D_BN)) then
          call kill(el%D_BN)
          DEALLOCATE(EL%D_BN)
       endif
       IF(ASSOCIATED(EL%D0_AN)) then
          call kill(el%D0_AN)
          DEALLOCATE(EL%D0_AN)
       endif
       IF(ASSOCIATED(EL%D0_BN)) then
          call kill(el%D0_BN)
          DEALLOCATE(EL%D0_BN)
       endif



       call kill(EL%P)        ! call kill(EL%P)    ! AIMIN MS 4.0

       !       IF(ASSOCIATED(EL%R)) DEALLOCATE(EL%R)
       !       IF(ASSOCIATED(EL%D)) DEALLOCATE(EL%D)
       !       IF(ASSOCIATED(EL%B_SOL)) DEALLOCATE(EL%B_SOL)  ! sagan

       IF(ASSOCIATED(EL%B_SOL)) then ! sagan
          CALL KILL(EL%B_SOL) ! sagan
          DEALLOCATE(EL%B_SOL)     ! sagan
       endif   ! sagan

       IF(ASSOCIATED(EL%THIN)) DEALLOCATE(EL%THIN)


    elseif(I>=0)       then

       !FIRST nullifies


       call null_ELEment(el)

       call alloc(el%P)

       ALLOCATE(EL%KIND);EL%KIND=0;ALLOCATE(EL%KNOB);EL%KNOB=.FALSE.;
       ALLOCATE(EL%NAME);ALLOCATE(EL%VORNAME);ALLOCATE(EL%electric);
       ALLOCATE(EL%skip_ptc_f);   EL%skip_ptc_f=0 ;   ALLOCATE(EL%skip_ptc_b);EL%skip_ptc_b=0  ;
       ALLOCATE(el%do1mapb);   el%do1mapb=.false. ;   ALLOCATE(el%do1mapf);el%do1mapf=.false.  ;
  ALLOCATE(el%usef);   el%usef=.false. ;   ALLOCATE(el%useb);el%useb=.false.  ;

       EL%NAME=' ';EL%NAME=TRIM(ADJUSTL(EL%NAME));
       EL%VORNAME=' ';EL%VORNAME=TRIM(ADJUSTL(EL%VORNAME));
       EL%electric=solve_electric
!       ALLOCATE(EL%PERMFRINGE);EL%PERMFRINGE=.FALSE.;  ! PART OF A STATE INITIALIZED BY EL=DEFAULT
       ALLOCATE(EL%L);CALL ALLOC(EL%L);EL%L=0.0_dp;
       ALLOCATE(EL%MIS);
       ! ALLOCATE(EL%EXACTMIS);
       EL%MIS=.FALSE.;
       !  EL%EXACTMIS=ALWAYS_EXACTMIS;
       !       allocate(el%r(3));allocate(el%d(3));
       !       el%r=zero;el%d=zero;
       !      EL=DEFAULT;
       !   ANBN
       CALL ZERO_ANBN(EL,I)
       ALLOCATE(EL%FINT(2));CALL ALLOC(EL%FINT);EL%FINT(1)=0.5_dp;EL%FINT(2)=0.5_dp;
       ALLOCATE(EL%HGAP(2));CALL ALLOC(EL%HGAP);EL%HGAP(1)=0.0_dp;EL%HGAP(2)=0.0_dp;
       ALLOCATE(EL%H1);CALL ALLOC(EL%H1);EL%H1=0.0_dp;
       ALLOCATE(EL%H2);CALL ALLOC(EL%H2);EL%H2=0.0_dp;
       ALLOCATE(EL%VA);CALL ALLOC(EL%VA);EL%VA=0.0_dp;
       ALLOCATE(EL%VS);CALL ALLOC(EL%VS);EL%VS=0.0_dp;
       !       ALLOCATE(EL%theta_ac);CALL ALLOC(EL%theta_ac); EL%theta_ac= zero ;
       !       ALLOCATE(EL%a_ac);CALL ALLOC(EL%a_ac);  EL%a_ac = zero;
       !       ALLOCATE(EL%DC_ac); EL%DC_ac= zero ;
       ALLOCATE(EL%slow_ac); EL%slow_ac=0 ;
    ENDIF

  END SUBROUTINE ZERO_ELP

  SUBROUTINE cop_el_elp(EL,ELP)
    IMPLICIT NONE
    TYPE(ELEMENT),INTENT(IN)::  EL
    TYPE(ELEMENTP),INTENT(inOUT)::  ELP
    CALL EQUAL(ELP,EL)
  END SUBROUTINE cop_el_elp

  SUBROUTINE cop_elp_el(EL,ELP)
    IMPLICIT NONE
    TYPE(ELEMENTP),INTENT(IN)::  EL
    TYPE(ELEMENT),INTENT(inOUT)::  ELP
    CALL EQUAL(ELP,EL)
  END SUBROUTINE       cop_elp_el

  SUBROUTINE cop_el_el(EL,ELP)
    IMPLICIT NONE
    TYPE(ELEMENT),INTENT(IN)::  EL
    TYPE(ELEMENT),INTENT(inOUT)::  ELP
    CALL EQUAL(ELP,EL)
  END SUBROUTINE       cop_el_el



  SUBROUTINE copy_el_elp(ELP,EL)
    IMPLICIT NONE
    TYPE(ELEMENT),INTENT(IN)::  EL
    TYPE(ELEMENTP),INTENT(inOUT)::  ELP
    INTEGER J,i,N

!    ELP%PERMFRINGE=EL%PERMFRINGE
    ELP%NAME=EL%NAME
    ELP%electric=EL%electric
    ELP%vorname=EL%vorname
    ELP%KIND=EL%KIND
    ELP%L=EL%L
    ELP%FINT(1)=EL%FINT(1)
    ELP%FINT(2)=EL%FINT(2)
    ELP%HGAP(1)=EL%HGAP(1)
    ELP%HGAP(2)=EL%HGAP(2)
    ELP%H1=EL%H1
    ELP%H2=EL%H2
    ELP%VA=EL%VA
    ELP%VS=EL%VS
    !    if(associated(el%siamese)) elp%siamese=>el%siamese
    !    if(associated(el%girder)) elp%girder=>el%girder
    ELP%slow_ac=EL%slow_ac

    IF(ASSOCIATED(EL%a_ac)) then
       ELP%a_ac=EL%a_ac
    endif
    IF(ASSOCIATED(EL%theta_ac)) then
       ELP%theta_ac=EL%theta_ac
    endif
    IF(ASSOCIATED(EL%DC_ac)) then
       ELP%DC_ac=EL%DC_ac
    endif



    IF(ASSOCIATED(EL%D_AN)) then

       IF(EL%P%NMUL>0) THEN
          IF(EL%P%NMUL/=ELP%P%NMUL.and.ELP%P%NMUL/=0) THEN
             call kill(ELP%D_AN,ELP%P%NMUL);call kill(ELP%D_bN,ELP%P%NMUL);
             call kill(ELP%D0_AN,ELP%P%NMUL);call kill(ELP%D0_bN,ELP%P%NMUL);
             DEALLOCATE(ELP%D_AN );DEALLOCATE(ELP%D_BN )
             DEALLOCATE(ELP%D0_AN );DEALLOCATE(ELP%D0_BN )
          endif
          if(.not.ASSOCIATED(ELP%D_AN)) THEN
             ALLOCATE(ELP%D_AN(EL%P%NMUL),ELP%D_BN(EL%P%NMUL))
             ALLOCATE(ELP%D0_AN(EL%P%NMUL),ELP%D0_BN(EL%P%NMUL))
          ENDIF


          CALL ALLOC(ELP%D_AN,EL%P%NMUL)
          CALL ALLOC(ELP%D_BN,EL%P%NMUL)
          CALL ALLOC(ELP%D0_AN,EL%P%NMUL)
          CALL ALLOC(ELP%D0_BN,EL%P%NMUL)
          DO I=1,EL%P%NMUL
             ELP%D_AN(I) = EL%D_AN(I)
             ELP%D_BN(I) = EL%D_BN(I)
             ELP%D0_AN(I) = EL%D0_AN(I)
             ELP%D0_BN(I) = EL%D0_BN(I)
          ENDDO

       ENDIF

    endif




    IF(EL%P%NMUL>0) THEN
       IF(EL%P%NMUL/=ELP%P%NMUL.and.ELP%P%NMUL/=0) THEN
          call kill(ELP%AN,ELP%P%NMUL);call kill(ELP%bN,ELP%P%NMUL);
          DEALLOCATE(ELP%AN );DEALLOCATE(ELP%BN )
       endif
       if(.not.ASSOCIATED(ELP%AN)) THEN
          ALLOCATE(ELP%AN(EL%P%NMUL),ELP%BN(EL%P%NMUL))
       ENDIF


       CALL ALLOC(ELP%AN,EL%P%NMUL)
       CALL ALLOC(ELP%BN,EL%P%NMUL)
       DO I=1,EL%P%NMUL
          ELP%AN(I) = EL%AN(I)
          ELP%BN(I) = EL%BN(I)
       ENDDO

    ENDIF
    ELP%P=EL%P

    ! MISALIGNMENTS
    ELP%MIS=EL%MIS
    !    ELP%EXACTMIS=EL%EXACTMIS

    !    IF(ASSOCIATED(EL%R)) THEN
    !       if(.not.ASSOCIATED(ELP%R))  ALLOCATE(ELP%R(3))

    !       DO I=1,3
    !          ELP%R(I)=EL%R(I)
    !       ENDDO
    !    ENDIF
    !    IF(ASSOCIATED(EL%D)) THEN
    !       if(.not.ASSOCIATED(ELP%D))  ALLOCATE(ELP%D(3))

    !       DO I=1,3
    !          ELP%D(I)=EL%D(I)
    !       ENDDO
    !    ENDIF

    IF(EL%KIND==KIND1) CALL SETFAMILY(ELP)
    IF(EL%KIND==KIND2) then
  !     write(6,*) associated(elp%k2)
  !     ELP%K2=0   ! new 2014.8.5
       CALL SETFAMILY(ELP)
       ELP%K2%F=EL%K2%F
    ENDIF
    IF(EL%KIND==KIND16.OR.EL%KIND==KIND20) THEN
       CALL SETFAMILY(ELP)
       ELP%K16%DRIFTKICK=EL%K16%DRIFTKICK
       ELP%K16%LIKEMAD=EL%K16%LIKEMAD
       ELP%K16%F=EL%K16%F
    ENDIF

    IF(EL%KIND==KIND3) THEN
       if(.not.ASSOCIATED(ELP%K3)) ALLOCATE(ELP%K3)
       ELP%K3=0
       if(.not.ASSOCIATED(ELP%B_SOL)) ALLOCATE(ELP%B_SOL       )
       CALL ALLOC( ELP%B_SOL)
       ELP%B_SOL = EL%B_SOL
       CALL SETFAMILY(ELP)
       ELP%K3%hf=EL%K3%hf
       ELP%K3%vf=EL%K3%vf
       ELP%K3%thin_h_foc=EL%K3%thin_h_foc
       ELP%K3%thin_v_foc=EL%K3%thin_v_foc
       ELP%K3%thin_h_angle=EL%K3%thin_h_angle
       ELP%K3%thin_v_angle=EL%K3%thin_v_angle
       ELP%K3%patch=EL%K3%patch
       ELP%K3%ls=EL%K3%ls
       ELP%k3%DX=EL%k3%DX
       ELP%k3%DY=EL%k3%DY
       ELP%k3%PITCH_X=EL%k3%PITCH_X
       ELP%k3%PITCH_Y=EL%k3%PITCH_Y
    ENDIF


    IF(EL%KIND==KIND4) THEN         !
       if(.not.ASSOCIATED(ELP%C4)) ALLOCATE(ELP%C4)
       ELP%C4=0
       if(.not.ASSOCIATED(ELP%VOLT)) ALLOCATE(ELP%VOLT,ELP%FREQ,ELP%PHAS,ELP%DELTA_E       )
       if(.not.ASSOCIATED(ELP%THIN)) ALLOCATE(ELP%THIN       )
       CALL ALLOC( ELP%VOLT)
       CALL ALLOC( ELP%FREQ)
       CALL ALLOC( ELP%PHAS)
       ELP%VOLT = EL%VOLT
       ELP%FREQ = EL%FREQ
       ELP%PHAS = EL%PHAS
       ELP%DELTA_E = EL%DELTA_E               ! DELTA_E IS real(dp)
       ELP%THIN = EL%THIN
       N_CAV4_F=EL%C4%NF
       CALL SETFAMILY(ELP)
       ELP%C4%N_BESSEL = EL%C4%N_BESSEL
       ELP%C4%cavity_totalpath = EL%C4%cavity_totalpath
       ELP%C4%phase0 = EL%C4%phase0
       DO I=1,EL%C4%NF
          ELP%C4%F(I)=EL%C4%F(I)
          ELP%C4%PH(I)=EL%C4%PH(I)
       ENDDO
       ELP%C4%t=EL%C4%t
       ELP%C4%R=EL%C4%R
       ELP%C4%A=EL%C4%A
       ELP%C4%Always_on=EL%C4%Always_on
    ENDIF

    IF(EL%KIND==kindsuperdrift) THEN         !
       if(.not.ASSOCIATED(ELP%SDR)) ALLOCATE(ELP%SDR)
       ELP%SDR=0

       CALL SETFAMILY(ELP)

       ELP%SDR%a_x1=EL%SDR%a_x1
       ELP%SDR%a_x2=EL%SDR%a_x2
       ELP%SDR%D=EL%SDR%D
       ELP%SDR%ANG=EL%SDR%ANG


    ENDIF
    IF(EL%KIND==KIND21) THEN         !
       if(.not.ASSOCIATED(ELP%CAV21)) ALLOCATE(ELP%CAV21)
       ELP%CAV21=0
       if(.not.ASSOCIATED(ELP%VOLT)) ALLOCATE(ELP%VOLT,ELP%FREQ,ELP%PHAS,ELP%DELTA_E       )
       if(.not.ASSOCIATED(ELP%THIN)) ALLOCATE(ELP%THIN       )
       CALL ALLOC( ELP%VOLT)
       CALL ALLOC( ELP%FREQ)
       CALL ALLOC( ELP%PHAS)
       ELP%VOLT = EL%VOLT
       ELP%FREQ = EL%FREQ
       ELP%PHAS = EL%PHAS
       ELP%DELTA_E = EL%DELTA_E               ! DELTA_E IS real(dp)
       ELP%THIN = EL%THIN
       CALL SETFAMILY(ELP)
       ELP%CAV21%PSI = EL%CAV21%PSI
       ELP%CAV21%DVDS = EL%CAV21%DVDS
       ELP%CAV21%DPHAS = EL%CAV21%DPHAS
       ELP%CAV21%cavity_totalpath = EL%CAV21%cavity_totalpath
       ELP%CAV21%phase0 = EL%CAV21%phase0
       ELP%CAV21%Always_on=EL%CAV21%Always_on
    ENDIF

    IF(EL%KIND==KIND22) THEN         !
       if(.not.ASSOCIATED(ELP%HE22)) ALLOCATE(ELP%HE22)
       ELP%HE22=0
       if(.not.ASSOCIATED(ELP%FREQ)) ALLOCATE(ELP%FREQ,ELP%PHAS)
       CALL ALLOC( ELP%FREQ)
       CALL ALLOC( ELP%PHAS)
       ELP%FREQ = EL%FREQ
       ELP%PHAS = EL%PHAS
       CALL SETFAMILY(ELP)
       ELp%HE22%N_BESSEL=EL%HE22%N_BESSEL 
       do i=1,6
        ELp%HE22%fake_shift(i)=EL%HE22%fake_shift(i)
       enddo
    ENDIF

    IF(EL%KIND==KIND5) THEN         !
       if(.not.ASSOCIATED(ELP%B_SOL)) ALLOCATE(ELP%B_SOL       )
       CALL ALLOC( ELP%B_SOL)
       if(.not.ASSOCIATED(ELP%s5)) ALLOCATE(ELP%s5)
       ELP%S5=0 ! 2014.8.5
       ELP%B_SOL = EL%B_SOL
       CALL SETFAMILY(ELP)
       ELP%S5%DX=EL%S5%DX
       ELP%S5%DY=EL%S5%DY
       ELP%S5%PITCH_X=EL%S5%PITCH_X
       ELP%S5%PITCH_Y=EL%S5%PITCH_Y
    ENDIF



    IF(EL%KIND==KIND6) CALL SETFAMILY(ELP)



    IF(EL%KIND==KIND7) THEN         !
       GEN=.FALSE.
       CALL SETFAMILY(ELP)
       IF(.NOT.GEN) THEN !.NOT.GEN
          ELP%T7%F=EL%T7%F
          DO J=1,3
             ELP%T7%LX(J)=EL%T7%LX(J)
             ELP%T7%RLX(J)=EL%T7%RLX(J)
             DO I=1,2
                ELP%T7%MATX(I,J)=EL%T7%MATX(I,J)
                ELP%T7%MATY(I,J)=EL%T7%MATY(I,J)
                ELP%T7%RMATX(I,J)=EL%T7%RMATX(I,J)
                ELP%T7%RMATY(I,J)=EL%T7%RMATY(I,J)
             ENDDO
          ENDDO
       ENDIF !.NOT.GEN
       GEN=.TRUE.
    ENDIF

    IF(EL%KIND==KIND8) CALL SETFAMILY(ELP)

    IF(EL%KIND==KIND9) CALL SETFAMILY(ELP)

    IF(EL%KIND==KIND10) THEN
       CALL SETFAMILY(ELP)
       ELP%TP10%DRIFTKICK=EL%TP10%DRIFTKICK
       ELP%TP10%F=EL%TP10%F
       IF(EL%ELECTRIC) THEN
 !     do i=1,S_E%N_MONO
 !       ELP%TP10%E_X(I)=EL%TP10%E_X(i)
 !       ELP%TP10%E_Y(I)=EL%TP10%E_Y(I)
 !     enddo
 !     do i=1,S_E%N_MONO
 !       ELP%TP10%PHI(I)=EL%TP10%PHI(I)
 !     enddo

        DO I=1,SECTOR_NMUL_max     
         ELP%TP10%AE(I)=EL%TP10%AE(I)     
         ELP%TP10%BE(I)=EL%TP10%BE(I)     
        enddo   
        call GETAEBE(ELP%TP10)     
       ENDIF
    ENDIF

    IF(EL%KIND>=KIND11.AND.EL%KIND<=KIND14) THEN
       CALL SETFAMILY(ELP)
       ELP%MON14%X=EL%MON14%X
       ELP%MON14%Y=EL%MON14%Y
    ENDIF

    IF(EL%KIND==KIND18) THEN
       CALL SETFAMILY(ELP)
       !ELP%RCOL18%A=EL%RCOL18%A
    ENDIF

    IF(EL%KIND==KIND19) THEN
       CALL SETFAMILY(ELP)
     !  ELP%ECOL19%A=EL%ECOL19%A
    ENDIF

    IF(EL%KIND==KIND15) THEN         !
       if(.not.ASSOCIATED(ELP%VOLT)) ALLOCATE(ELP%VOLT)
       if(.not.ASSOCIATED(ELP%PHAS)) ALLOCATE(ELP%PHAS)
       CALL ALLOC( ELP%VOLT)
       CALL ALLOC( ELP%PHAS)
       ELP%VOLT = EL%VOLT
       ELP%PHAS = EL%PHAS
       CALL SETFAMILY(ELP)
    ENDIF

    !    IF(EL%KIND==KINDUSER1) THEN         !
    !       CALL SETFAMILY(ELP)
    !       CALL COPY(EL%U1,ELP%U1)
    !    ENDIF

    !    IF(EL%KIND==KINDUSER2) THEN         !
    !       CALL SETFAMILY(ELP)
    !       CALL COPY(EL%U2,ELP%U2)
    !    ENDIF

    IF(EL%KIND==KINDWIGGLER) THEN         !
       CALL SETFAMILY(ELP)
       CALL COPY(EL%WI,ELP%WI)
    ENDIF

    IF(ASSOCIATED(EL%RAMP)) THEN         !
       CALL COPY_RAMPING(EL%RAMP,ELP%RAMP)
    ENDIF
    
    
    
    IF(EL%KIND==KINDPA) THEN         !
       CALL SETFAMILY(ELP,t=EL%PA%B)  !,EL%PA%angc,EL%PA%xc,EL%PA%dc,EL%PA%h)  !,EL%PA%ax,EL%PA%ay)
       CALL COPY(EL%PA,ELP%PA)
    ENDIF
    IF(EL%KIND==KINDABELL) THEN 
       M_ABELL=EL%ab%M
       N_ABELL=EL%ab%N       !
       CALL SETFAMILY(ELP)  !,EL%PA%angc,EL%PA%xc,EL%PA%dc,EL%PA%h)  !,EL%PA%ax,EL%PA%ay)
       CALL COPY(EL%ab,ELP%ab)
    ENDIF
    !    IF(ASSOCIATED(EL%PARENT_FIBRE))        then
    !       ELP%PARENT_FIBRE=>EL%PARENT_FIBRE
    !    ENDIF

       IF(ASSOCIATED(EL%backWARD))        then
         if(associated(elp%backWARD)) then
          call kill(ELp%backWARD)     
          DEALLOCATE(ELp%backWARD)
         endif
         allocate(ELp%backWARD(3))
         do i=1,3
         call alloc_tree(ELp%backWARD(i),EL%backWARD(i)%n,EL%backWARD(i)%np)
         enddo
         call COPY_TREE_N(EL%backWARD,ELp%backWARD)

       ENDIF


       IF(ASSOCIATED(EL%forward))        then
         if(associated(elp%forward)) then
          call kill(ELp%forward)     
          DEALLOCATE(ELp%forward)
         endif
         allocate(elp%forward(3))
         do i=1,3
         call alloc_tree(elp%forward(i),el%forward(i)%n,el%forward(i)%np)
         enddo
         call COPY_TREE_N(EL%forward,ELp%forward)

       ENDIF
       ELp%skip_ptc_f=EL%skip_ptc_f
       ELp%skip_ptc_b=EL%skip_ptc_b
       elp%do1mapf=el%do1mapf
       elp%do1mapb=el%do1mapb
          ELp%usef= EL%usef
          ELp%useb= EL%useb
  END SUBROUTINE copy_el_elp





  SUBROUTINE copy_elp_el(ELP,EL)
    IMPLICIT NONE
    TYPE(ELEMENTP),INTENT(IN)::  EL
    TYPE(ELEMENT),INTENT(inOUT)::  ELP
    INTEGER I,J,N

    !    if(associated(el%siamese)) elp%siamese=>el%siamese
    !    if(associated(el%girder)) elp%girder=>el%girder
!    ELP%PERMFRINGE=EL%PERMFRINGE
    ELP%electric=EL%electric
    ELP%vorname=EL%vorname
    ELP%KIND=EL%KIND
    ELP%L=EL%L
    ELP%FINT(1)=EL%FINT(1)
    ELP%FINT(2)=EL%FINT(2)
    ELP%HGAP(1)=EL%HGAP(1)
    ELP%HGAP(2)=EL%HGAP(2)
    ELP%H1=EL%H1
    ELP%H2=EL%H2
    ELP%VA=EL%VA
    ELP%VS=EL%VS
    ELP%slow_ac=EL%slow_ac

    IF(ASSOCIATED(EL%a_ac)) then
       ELP%a_ac=EL%a_ac
    endif
    IF(ASSOCIATED(EL%theta_ac)) then
       ELP%theta_ac=EL%theta_ac
    endif
    IF(ASSOCIATED(EL%DC_ac)) then
       ELP%DC_ac=EL%DC_ac
    endif


    IF(ASSOCIATED(EL%D_AN)) then

       IF(EL%P%NMUL>0) THEN
          IF(EL%P%NMUL/=ELP%P%NMUL.and.ELP%P%NMUL/=0) THEN
             DEALLOCATE(ELP%D_AN );DEALLOCATE(ELP%D_BN )
             DEALLOCATE(ELP%D0_AN );DEALLOCATE(ELP%D0_BN )
          endif
          if(.not.ASSOCIATED(ELP%D_AN)) THEN
             ALLOCATE(ELP%D_AN(EL%P%NMUL),ELP%D_BN(EL%P%NMUL))
             ALLOCATE(ELP%D0_AN(EL%P%NMUL),ELP%D0_BN(EL%P%NMUL))
          ENDIF

          DO I=1,EL%P%NMUL
             ELP%D_AN(I) = EL%D_AN(I)
             ELP%D_BN(I) = EL%D_BN(I)
             ELP%D0_AN(I) = EL%D0_AN(I)
             ELP%D0_BN(I) = EL%D0_BN(I)
          ENDDO

       ENDIF

    endif





    IF(EL%P%NMUL>0) THEN
       IF(EL%P%NMUL/=ELP%P%NMUL.and.ELP%P%NMUL/=0) THEN
          DEALLOCATE(ELP%AN );DEALLOCATE(ELP%BN )
       endif
       if(.not.ASSOCIATED(ELP%AN)) THEN
          ALLOCATE(ELP%AN(EL%P%NMUL),ELP%BN(EL%P%NMUL))
       ENDIF

       DO I=1,EL%P%NMUL
          ELP%AN(I) = EL%AN(I)
          ELP%BN(I) = EL%BN(I)
       ENDDO

    ENDIF
    ELP%P=EL%P



    ! MISALIGNMENTS
    ELP%MIS=EL%MIS
    !    ELP%EXACTMIS=EL%EXACTMIS

    !    IF(ASSOCIATED(EL%R)) THEN
    !       if(.not.ASSOCIATED(ELP%R))  ALLOCATE(ELP%R(3))

    !       DO I=1,3
    !          ELP%R(I)=EL%R(I)
    !       ENDDO
    !    ENDIF
    !    IF(ASSOCIATED(EL%D)) THEN
    !       if(.not.ASSOCIATED(ELP%D))  ALLOCATE(ELP%D(3))

    !       DO I=1,3
    !          ELP%D(I)=EL%D(I)
    !       ENDDO
    !    ENDIF

    IF(EL%KIND==KIND1) CALL SETFAMILY(ELP)

    IF(EL%KIND==KIND2) then
 !      ELP%K2=0   ! new 2014.8.5
       CALL SETFAMILY(ELP)
       ELP%K2%F=EL%K2%F
    ENDIF
    IF(EL%KIND==KIND16.OR.EL%KIND==KIND20) THEN
       CALL SETFAMILY(ELP)
       ELP%K16%DRIFTKICK=EL%K16%DRIFTKICK
       ELP%K16%LIKEMAD=EL%K16%LIKEMAD
       ELP%K16%F=EL%K16%F
    ENDIF

    IF(EL%KIND==KIND3) THEN
       if(.not.ASSOCIATED(ELP%K3)) ALLOCATE(ELP%K3)
       ELP%K3=0
       if(.not.ASSOCIATED(ELP%B_SOL)) ALLOCATE(ELP%B_SOL       )
       ELP%B_SOL = EL%B_SOL
       CALL SETFAMILY(ELP)
       ELP%K3%hf=EL%K3%hf
       ELP%K3%vf=EL%K3%vf
       ELP%K3%thin_h_foc=EL%K3%thin_h_foc
       ELP%K3%thin_v_foc=EL%K3%thin_v_foc
       ELP%K3%thin_h_angle=EL%K3%thin_h_angle
       ELP%K3%thin_v_angle=EL%K3%thin_v_angle
       ELP%K3%patch=EL%K3%patch
       ELP%K3%ls=EL%K3%ls
       ELP%k3%DX=EL%k3%DX
       ELP%k3%DY=EL%k3%DY
       ELP%k3%PITCH_X=EL%k3%PITCH_X
       ELP%k3%PITCH_Y=EL%k3%PITCH_Y
    ENDIF


    IF(EL%KIND==KIND4) THEN         !
       if(.not.ASSOCIATED(ELP%C4)) ALLOCATE(ELP%C4)
       ELP%C4=0
       if(.not.ASSOCIATED(ELP%VOLT)) ALLOCATE(ELP%VOLT,ELP%FREQ,ELP%PHAS,ELP%DELTA_E       )
       if(.not.ASSOCIATED(ELP%THIN)) ALLOCATE(ELP%THIN       )
       ELP%VOLT = EL%VOLT
       ELP%FREQ = EL%FREQ
       ELP%PHAS = EL%PHAS
       ELP%DELTA_E = EL%DELTA_E
       ELP%THIN = EL%THIN
       N_CAV4_F=EL%C4%NF
       CALL SETFAMILY(ELP)
       ELP%C4%N_BESSEL = EL%C4%N_BESSEL
       ELP%C4%cavity_totalpath = EL%C4%cavity_totalpath
       ELP%C4%phase0 = EL%C4%phase0
       DO I=1,EL%C4%NF
          ELP%C4%F(I)=EL%C4%F(I)
          ELP%C4%PH(I)=EL%C4%PH(I)
       ENDDO
       ELP%C4%t=EL%C4%t
       ELP%C4%R=EL%C4%R
       ELP%C4%A=EL%C4%A
       ELP%C4%Always_on=EL%C4%Always_on
    ENDIF

    IF(EL%KIND==kindsuperdrift) THEN         !
       if(.not.ASSOCIATED(ELP%SDR)) ALLOCATE(ELP%SDR)
       ELP%SDR=0

       CALL SETFAMILY(ELP)

       ELP%SDR%a_x1=EL%SDR%a_x1
       ELP%SDR%a_x2=EL%SDR%a_x2
       ELP%SDR%D=EL%SDR%D
       ELP%SDR%ANG=EL%SDR%ANG

    ENDIF

    IF(EL%KIND==KIND21) THEN         !
       if(.not.ASSOCIATED(ELP%CAV21)) ALLOCATE(ELP%CAV21)
       ELP%CAV21=0
       if(.not.ASSOCIATED(ELP%VOLT)) ALLOCATE(ELP%VOLT,ELP%FREQ,ELP%PHAS,ELP%DELTA_E       )
       if(.not.ASSOCIATED(ELP%THIN)) ALLOCATE(ELP%THIN       )
       ELP%VOLT = EL%VOLT
       ELP%FREQ = EL%FREQ
       ELP%PHAS = EL%PHAS
       ELP%DELTA_E = EL%DELTA_E
       ELP%THIN = EL%THIN
       CALL SETFAMILY(ELP)
       ELP%CAV21%PSI = EL%CAV21%PSI
       ELP%CAV21%DVDS = EL%CAV21%DVDS
       ELP%CAV21%DPHAS = EL%CAV21%DPHAS
       ELP%CAV21%cavity_totalpath = EL%CAV21%cavity_totalpath
       ELP%CAV21%phase0 = EL%CAV21%phase0
       ELP%CAV21%Always_on=EL%CAV21%Always_on
    ENDIF

    IF(EL%KIND==KIND22) THEN         !
       if(.not.ASSOCIATED(ELP%HE22)) ALLOCATE(ELP%HE22)
       ELP%HE22=0
       if(.not.ASSOCIATED(ELP%FREQ)) ALLOCATE(ELP%FREQ,ELP%PHAS)
       ELP%FREQ = EL%FREQ
       ELP%PHAS = EL%PHAS
       CALL SETFAMILY(ELP)
       ELp%HE22%N_BESSEL=EL%HE22%N_BESSEL
       do i=1,6
        ELp%HE22%fake_shift(i)=EL%HE22%fake_shift(i)
       enddo
    ENDIF

    IF(EL%KIND==KIND5) THEN         !
       if(.not.ASSOCIATED(ELP%S5)) ALLOCATE(ELP%S5)
       if(.not.ASSOCIATED(ELP%B_SOL)) ALLOCATE(ELP%B_SOL       )

       ELP%S5=0 ! 2014.8.5      
       ELP%B_SOL = EL%B_SOL
       CALL SETFAMILY(ELP)
       ELP%S5%DX=EL%S5%DX
       ELP%S5%DY=EL%S5%DY
       ELP%S5%PITCH_X=EL%S5%PITCH_X
       ELP%S5%PITCH_Y=EL%S5%PITCH_Y
    ENDIF

    !    IF(EL%KIND==KIND17) THEN         !
    !       !       if(.not.ASSOCIATED(ELP%S17)) ALLOCATE(ELP%S17)
    !       if(.not.ASSOCIATED(ELP%B_SOL)) ALLOCATE(ELP%B_SOL       )
    !       ELP%B_SOL = EL%B_SOL
    !       CALL SETFAMILY(ELP)
    !    ENDIF

    IF(EL%KIND==KIND6) CALL SETFAMILY(ELP)

    !    IF(EL%KIND==KIND22) THEN
    !       i=0;j=0;k=0;l=0;
    !       if(associated(EL%M22%T_REV)) i=EL%M22%T_REV%N
    !       if(associated(EL%M22%T_rad_REV)) j=EL%M22%T_rad_REV%N
    !       if(associated(EL%M22%T)) k=EL%M22%T%N
    !       if(associated(EL%M22%T_rad)) l=EL%M22%T_rad%N
    !       CALL SETFAMILY(ELP,NTOT=k,ntot_rad=l,NTOT_REV=i,ntot_rad_REV=j,ND2=6)
    !       ELP%M22%DELTAMAP=EL%M22%DELTAMAP
    !
    !       if(associated(EL%M22%T))  CALL COPY_TREE(EL%M22%T,ELP%M22%T)
    !       if(associated(EL%M22%T_rad)) CALL COPY_TREE(EL%M22%T_rad,ELP%M22%T_rad)
    !       if(associated(EL%M22%T_REV)) CALL COPY_TREE(EL%M22%T_REV,ELP%M22%T_REV)
    !       if(associated(EL%M22%T_rad_REV)) CALL COPY_TREE(EL%M22%T_rad_REV,ELP%M22%T_rad_REV)
    !    ENDIF

    IF(EL%KIND==KIND7) THEN         !
       GEN=.FALSE.
       CALL SETFAMILY(ELP)
       IF(.NOT.GEN) THEN !.NOT.GEN
          ELP%T7%F=EL%T7%F
          DO J=1,3
             ELP%T7%LX(J)=EL%T7%LX(J)
             ELP%T7%RLX(J)=EL%T7%RLX(J)
             DO I=1,2
                ELP%T7%MATX(I,J)=EL%T7%MATX(I,J)
                ELP%T7%MATY(I,J)=EL%T7%MATY(I,J)
                ELP%T7%RMATX(I,J)=EL%T7%RMATX(I,J)
                ELP%T7%RMATY(I,J)=EL%T7%RMATY(I,J)
             ENDDO
          ENDDO
       ENDIF !.NOT.GEN
       GEN=.TRUE.

    ENDIF


    IF(EL%KIND==KIND8) CALL SETFAMILY(ELP)

    IF(EL%KIND==KIND9) CALL SETFAMILY(ELP)

    IF(EL%KIND==KIND10) THEN
       CALL SETFAMILY(ELP)
       ELP%TP10%DRIFTKICK=EL%TP10%DRIFTKICK
       ELP%TP10%F=EL%TP10%F
       IF(EL%ELECTRIC) THEN
 !     do i=1,S_E%N_MONO
 !       ELP%TP10%E_X(I)=EL%TP10%E_X(i)
 !       ELP%TP10%E_Y(I)=EL%TP10%E_Y(I)
 !     enddo
 !     do i=1,S_E%N_MONO
 !       ELP%TP10%PHI(I)=EL%TP10%PHI(I)
 !     enddo

        DO I=1,SECTOR_NMUL_max     
         ELP%TP10%AE(I)=EL%TP10%AE(I)     
         ELP%TP10%BE(I)=EL%TP10%BE(I)     
        enddo 
        call GETAEBE(ELP%TP10)         
       ENDIF
       
    ENDIF

    IF(EL%KIND>=KIND11.AND.EL%KIND<=KIND14) THEN
       CALL SETFAMILY(ELP)
       ELP%MON14%X=EL%MON14%X
       ELP%MON14%Y=EL%MON14%Y
    ENDIF

    IF(EL%KIND==KIND18) THEN
       CALL SETFAMILY(ELP)
    !   ELP%RCOL18%A=EL%RCOL18%A
    ENDIF

    IF(EL%KIND==KIND19) THEN
       CALL SETFAMILY(ELP)
     !  ELP%ECOL19%A=EL%ECOL19%A
    ENDIF

    IF(EL%KIND==KIND15) THEN         !
       if(.not.ASSOCIATED(ELP%SEP15)) ALLOCATE(ELP%SEP15)
       if(.not.ASSOCIATED(ELP%VOLT)) ALLOCATE(ELP%VOLT)
       if(.not.ASSOCIATED(ELP%PHAS)) ALLOCATE(ELP%PHAS)
       ELP%VOLT = EL%VOLT
       ELP%PHAS = EL%PHAS
       CALL SETFAMILY(ELP)
    ENDIF

    !    IF(EL%KIND==KINDUSER1) THEN         !
    !       CALL SETFAMILY(ELP)
    !       CALL COPY(EL%U1,ELP%U1)
    !    ENDIF

    !    IF(EL%KIND==KINDUSER2) THEN         !
    !       CALL SETFAMILY(ELP)
    !       CALL COPY(EL%U2,ELP%U2)
    !    ENDIF

    IF(EL%KIND==KINDWIGGLER) THEN         !
       CALL SETFAMILY(ELP)
       CALL COPY(EL%WI,ELP%WI)
    ENDIF
    
       IF(ASSOCIATED(EL%RAMP)) THEN         !
       CALL COPY_RAMPING(EL%RAMP,ELP%RAMP)
    ENDIF
 
    IF(EL%KIND==KINDPA) THEN         !
       CALL SETFAMILY(ELP,t=EL%PA%B) !,EL%PA%angc,EL%PA%xc,EL%PA%dc,EL%PA%h)  !,EL%PA%ax,EL%PA%ay)
       CALL COPY(EL%PA,ELP%PA)
    ENDIF

    IF(EL%KIND==KINDABELL) THEN 
       M_ABELL=EL%ab%M
       N_ABELL=EL%ab%N        !
       CALL SETFAMILY(ELP)  !,EL%PA%angc,EL%PA%xc,EL%PA%dc,EL%PA%h)  !,EL%PA%ax,EL%PA%ay)
       CALL COPY(EL%ab,ELP%ab)
    ENDIF

       IF(ASSOCIATED(EL%backWARD))        then
         if(associated(elp%backWARD)) then
          call kill(ELp%backWARD)     
          DEALLOCATE(ELp%backWARD)
         endif
         allocate(ELp%backWARD(3))
         do i=1,3
         call alloc_tree(ELp%backWARD(i),EL%backWARD(i)%n,EL%backWARD(i)%np)
         enddo
         call COPY_TREE_N(EL%backWARD,ELp%backWARD)
 
       ENDIF
 

       IF(ASSOCIATED(EL%forward))        then
         if(associated(elp%forward)) then
          call kill(ELp%forward)     
          DEALLOCATE(ELp%forward)
         endif
         allocate(ELp%forward(3))
         do i=1,3
         call alloc_tree(ELp%forward(i),EL%forward(i)%n,EL%forward(i)%np)
         enddo
         call COPY_TREE_N(EL%forward,ELp%forward)
 
       ENDIF
       ELp%skip_ptc_f=EL%skip_ptc_f
       ELp%skip_ptc_b=EL%skip_ptc_b
       elp%do1mapf=el%do1mapf
       elp%do1mapb=el%do1mapb
          ELp%usef= EL%usef
          ELp%useb= EL%useb
  END SUBROUTINE copy_elp_el



  SUBROUTINE copy_el_el(ELP,EL)
    IMPLICIT NONE
    TYPE(ELEMENT),INTENT(IN)::  EL
    TYPE(ELEMENT),INTENT(inOUT)::  ELP
    INTEGER I,J,n


    !    if(associated(el%siamese)) elp%siamese=>el%siamese
    !    if(associated(el%girder)) elp%girder=>el%girder
!    ELP%PERMFRINGE=EL%PERMFRINGE
    ELP%NAME=EL%NAME
    ELP%electric=EL%electric
    ELP%vorname=EL%vorname
    ELP%RECUT=EL%RECUT
    ELP%even=EL%even
    ELP%KIND=EL%KIND
    ELP%PLOT=EL%PLOT
    ELP%L=EL%L
    ELP%FINT(1)=EL%FINT(1)
    ELP%FINT(2)=EL%FINT(2)
    ELP%HGAP(1)=EL%HGAP(1)
    ELP%HGAP(2)=EL%HGAP(2)
    ELP%H1=EL%H1
    ELP%H2=EL%H2
    ELP%VA=EL%VA
    ELP%VS=EL%VS
        ELP%ene=EL%ene
    ELP%slow_ac=EL%slow_ac
        ELp%filef=EL%filef
        elp%fileb=EL%fileb

    IF(ASSOCIATED(EL%a_ac)) then
       ELP%a_ac=EL%a_ac
    endif
    IF(ASSOCIATED(EL%theta_ac)) then
       ELP%theta_ac=EL%theta_ac
    endif
    IF(ASSOCIATED(EL%DC_ac)) then
       ELP%DC_ac=EL%DC_ac
    endif

    IF(ASSOCIATED(EL%D_AN)) then

       IF(EL%P%NMUL>0) THEN
          IF(EL%P%NMUL/=ELP%P%NMUL.and.ELP%P%NMUL/=0) THEN
             DEALLOCATE(ELP%D_AN );DEALLOCATE(ELP%D_BN )
             DEALLOCATE(ELP%D0_AN );DEALLOCATE(ELP%D0_BN )
          endif
          if(.not.ASSOCIATED(ELP%D_AN)) THEN
             ALLOCATE(ELP%D_AN(EL%P%NMUL),ELP%D_BN(EL%P%NMUL))
             ALLOCATE(ELP%D0_AN(EL%P%NMUL),ELP%D0_BN(EL%P%NMUL))
          ENDIF

          DO I=1,EL%P%NMUL
             ELP%D_AN(I) = EL%D_AN(I)
             ELP%D_BN(I) = EL%D_BN(I)
             ELP%D0_AN(I) = EL%D0_AN(I)
             ELP%D0_BN(I) = EL%D0_BN(I)
          ENDDO

       ENDIF

    endif




    IF(EL%P%NMUL>0) THEN
       IF(EL%P%NMUL/=ELP%P%NMUL.and.ELP%P%NMUL/=0) THEN
          DEALLOCATE(ELP%AN );DEALLOCATE(ELP%BN )
       endif
       if(.not.ASSOCIATED(ELP%AN)) THEN
          ALLOCATE(ELP%AN(EL%P%NMUL),ELP%BN(EL%P%NMUL))
       ENDIF

       DO I=1,EL%P%NMUL
          ELP%AN(I) = EL%AN(I)
          ELP%BN(I) = EL%BN(I)
       ENDDO

    ENDIF
    ELP%P=EL%P



    ! MISALIGNMENTS
    ELP%MIS=EL%MIS
    !    ELP%EXACTMIS=EL%EXACTMIS

    !    IF(ASSOCIATED(EL%R)) THEN
    !       if(.not.ASSOCIATED(ELP%R))  ALLOCATE(ELP%R(3))
    !       DO I=1,3
    !          ELP%R(I)=EL%R(I)
    !       ENDDO
    !    ENDIF
    !   IF(ASSOCIATED(EL%D)) THEN
    !       if(.not.ASSOCIATED(ELP%D))  ALLOCATE(ELP%D(3))
    !       DO I=1,3
    !          ELP%D(I)=EL%D(I)
    !       ENDDO
    !    ENDIF

    IF(EL%KIND==KIND1) CALL SETFAMILY(ELP)

    IF(EL%KIND==KIND2) then
 !      ELP%K2=0   ! new 2014.8.5
       CALL SETFAMILY(ELP)
       ELP%K2%F=EL%K2%F
    ENDIF
    IF(EL%KIND==KIND16.OR.EL%KIND==KIND20) THEN
       CALL SETFAMILY(ELP)
       ELP%K16%DRIFTKICK=EL%K16%DRIFTKICK
       ELP%K16%LIKEMAD=EL%K16%LIKEMAD
       ELP%K16%F=EL%K16%F
    ENDIF

    IF(EL%KIND==KIND3) THEN
       if(.not.ASSOCIATED(ELP%K3)) ALLOCATE(ELP%K3)
       ELP%K3=0
       if(.not.ASSOCIATED(ELP%B_SOL)) ALLOCATE(ELP%B_SOL       )
       ELP%B_SOL = EL%B_SOL
       CALL SETFAMILY(ELP)
       ELP%K3%hf=EL%K3%hf
       ELP%K3%vf=EL%K3%vf
       ELP%K3%thin_h_foc=EL%K3%thin_h_foc
       ELP%K3%thin_v_foc=EL%K3%thin_v_foc
       ELP%K3%thin_h_angle=EL%K3%thin_h_angle
       ELP%K3%thin_v_angle=EL%K3%thin_v_angle
       ELP%K3%patch=EL%K3%patch
       ELP%K3%ls=EL%K3%ls
       ELP%k3%DX=EL%k3%DX
       ELP%k3%DY=EL%k3%DY
       ELP%k3%PITCH_X=EL%k3%PITCH_X
       ELP%k3%PITCH_Y=EL%k3%PITCH_Y
    ENDIF

    IF(EL%KIND==KIND4) THEN         !
       if(.not.ASSOCIATED(ELP%C4)) ALLOCATE(ELP%C4)
       ELP%C4=0
       if(.not.ASSOCIATED(ELP%VOLT)) ALLOCATE(ELP%VOLT,ELP%FREQ,ELP%PHAS,ELP%DELTA_E       )
       if(.not.ASSOCIATED(ELP%THIN)) ALLOCATE(ELP%THIN       )
       if(.not.ASSOCIATED(ELP%lag)) ALLOCATE(ELP%lag       )
       ELP%lag = EL%lag
       ELP%VOLT = EL%VOLT
       ELP%FREQ = EL%FREQ
       ELP%PHAS = EL%PHAS
       ELP%DELTA_E = EL%DELTA_E
       ELP%THIN = EL%THIN
       N_CAV4_F=EL%C4%NF
       CALL SETFAMILY(ELP)
       ELP%C4%N_BESSEL = EL%C4%N_BESSEL
       ELP%C4%cavity_totalpath = EL%C4%cavity_totalpath
       ELP%C4%phase0 = EL%C4%phase0
       DO I=1,EL%C4%NF
          ELP%C4%F(I)=EL%C4%F(I)
          ELP%C4%PH(I)=EL%C4%PH(I)
       ENDDO
       ELP%C4%t=EL%C4%t
       ELP%C4%R=EL%C4%R
       ELP%C4%A=EL%C4%A
       ELP%C4%Always_on=EL%C4%Always_on
    ENDIF


    IF(EL%KIND==kindsuperdrift) THEN         !
       if(.not.ASSOCIATED(ELP%SDR)) ALLOCATE(ELP%SDR)
       ELP%SDR=0

       CALL SETFAMILY(ELP)

       ELP%SDR%a_x1=EL%SDR%a_x1
       ELP%SDR%a_x2=EL%SDR%a_x2
       ELP%SDR%D=EL%SDR%D
       ELP%SDR%ANG=EL%SDR%ANG


    ENDIF

    IF(EL%KIND==KIND21) THEN         !
       if(.not.ASSOCIATED(ELP%CAV21)) ALLOCATE(ELP%CAV21)
       ELP%CAV21=0
       if(.not.ASSOCIATED(ELP%VOLT)) ALLOCATE(ELP%VOLT,ELP%FREQ,ELP%PHAS,ELP%DELTA_E       )
       if(.not.ASSOCIATED(ELP%THIN)) ALLOCATE(ELP%THIN       )
       if(.not.ASSOCIATED(ELP%lag)) ALLOCATE(ELP%lag       )
       ELP%lag = EL%lag
       ELP%VOLT = EL%VOLT
       ELP%FREQ = EL%FREQ
       ELP%PHAS = EL%PHAS
       ELP%DELTA_E = EL%DELTA_E
       ELP%THIN = EL%THIN
       CALL SETFAMILY(ELP)
       ELP%CAV21%PSI = EL%CAV21%PSI
       ELP%CAV21%DVDS = EL%CAV21%DVDS
       ELP%CAV21%DPHAS = EL%CAV21%DPHAS
       ELP%CAV21%cavity_totalpath = EL%CAV21%cavity_totalpath
       ELP%CAV21%phase0 = EL%CAV21%phase0
       ELP%CAV21%Always_on=EL%CAV21%Always_on
    ENDIF

    IF(EL%KIND==KIND22) THEN         !
       if(.not.ASSOCIATED(ELP%HE22)) ALLOCATE(ELP%HE22)
       ELP%HE22=0
       if(.not.ASSOCIATED(ELP%FREQ)) ALLOCATE(ELP%FREQ,ELP%PHAS)
       ELP%FREQ = EL%FREQ
       ELP%PHAS = EL%PHAS
       CALL SETFAMILY(ELP)
       ELp%HE22%N_BESSEL=EL%HE22%N_BESSEL
       do i=1,6
        ELp%HE22%fake_shift(i)=EL%HE22%fake_shift(i)
       enddo
    ENDIF

    IF(EL%KIND==KIND5) THEN         !
       if(.not.ASSOCIATED(ELP%S5)) ALLOCATE(ELP%S5)
       if(.not.ASSOCIATED(ELP%B_SOL)) ALLOCATE(ELP%B_SOL       )
       ELP%S5=0 ! 2014.8.5
       ELP%B_SOL = EL%B_SOL
       CALL SETFAMILY(ELP)
       ELP%S5%DX=EL%S5%DX
       ELP%S5%DY=EL%S5%DY
       ELP%S5%PITCH_X=EL%S5%PITCH_X
       ELP%S5%PITCH_Y=EL%S5%PITCH_Y
    ENDIF

    !    IF(EL%KIND==KIND17) THEN         !
    !       !      if(.not.ASSOCIATED(ELP%S17)) ALLOCATE(ELP%S17)
    !       if(.not.ASSOCIATED(ELP%B_SOL)) ALLOCATE(ELP%B_SOL       )
    !       ELP%B_SOL = EL%B_SOL
    !       CALL SETFAMILY(ELP)
    !    ENDIF

    IF(EL%KIND==KIND6) CALL SETFAMILY(ELP)

    !    IF(EL%KIND==KIND22) THEN
    !       i=0;j=0;k=0;l=0;
    !       if(associated(EL%M22%T_REV)) i=EL%M22%T_REV%N
    !       if(associated(EL%M22%T_rad_REV)) j=EL%M22%T_rad_REV%N
    !       if(associated(EL%M22%T)) k=EL%M22%T%N
    !       if(associated(EL%M22%T_rad)) l=EL%M22%T_rad%N
    !       CALL SETFAMILY(ELP,NTOT=k,ntot_rad=l,NTOT_REV=i,ntot_rad_REV=j,ND2=6)
    !       ELP%M22%DELTAMAP=EL%M22%DELTAMAP
    !
    !       if(associated(EL%M22%T))  CALL COPY_TREE(EL%M22%T,ELP%M22%T)
    !       if(associated(EL%M22%T_rad)) CALL COPY_TREE(EL%M22%T_rad,ELP%M22%T_rad)
    !       if(associated(EL%M22%T_REV)) CALL COPY_TREE(EL%M22%T_REV,ELP%M22%T_REV)
    !       if(associated(EL%M22%T_rad_REV)) CALL COPY_TREE(EL%M22%T_rad_REV,ELP%M22%T_rad_REV)
    !    ENDIF

    IF(EL%KIND==KIND7) THEN         !
       GEN=.FALSE.
       CALL SETFAMILY(ELP)
       IF(.NOT.GEN) THEN !.NOT.GEN
          ELP%T7%F=EL%T7%F
          DO J=1,3
             ELP%T7%LX(J)=EL%T7%LX(J)
             ELP%T7%RLX(J)=EL%T7%RLX(J)
             DO I=1,2
                ELP%T7%MATX(I,J)=EL%T7%MATX(I,J)
                ELP%T7%MATY(I,J)=EL%T7%MATY(I,J)
                ELP%T7%RMATX(I,J)=EL%T7%RMATX(I,J)
                ELP%T7%RMATY(I,J)=EL%T7%RMATY(I,J)
             ENDDO
          ENDDO
       ENDIF !.NOT.GEN
       GEN=.TRUE.
    ENDIF


    IF(EL%KIND==KIND8) CALL SETFAMILY(ELP)

    IF(EL%KIND==KIND9) CALL SETFAMILY(ELP)

    IF(EL%KIND==KIND10) THEN
       CALL SETFAMILY(ELP)
       ELP%TP10%DRIFTKICK=EL%TP10%DRIFTKICK
       ELP%TP10%F=EL%TP10%F
       IF(EL%ELECTRIC) THEN
 !     do i=1,S_E%N_MONO
 !       ELP%TP10%E_X(I)=EL%TP10%E_X(i)
 !       ELP%TP10%E_Y(I)=EL%TP10%E_Y(I)
 !     enddo
 !     do i=1,S_E%N_MONO
 !       ELP%TP10%PHI(I)=EL%TP10%PHI(I)
 !     enddo

        DO I=1,SECTOR_NMUL_max     
         ELP%TP10%AE(I)=EL%TP10%AE(I)     
         ELP%TP10%BE(I)=EL%TP10%BE(I)     
        enddo  
        call GETAEBE(ELP%TP10)         
       ENDIF
    ENDIF

    IF(EL%KIND>=KIND11.AND.EL%KIND<=KIND14) THEN
       CALL SETFAMILY(ELP)
       ELP%MON14%X=EL%MON14%X
       ELP%MON14%Y=EL%MON14%Y
    ENDIF

    IF(EL%KIND==KIND18) THEN
       CALL SETFAMILY(ELP)
     !  ELP%RCOL18%A=EL%RCOL18%A
    ENDIF

    IF(EL%KIND==KIND19) THEN
       CALL SETFAMILY(ELP)
    !   ELP%ECOL19%A=EL%ECOL19%A
    ENDIF

    IF(EL%KIND==KIND15) THEN         !
       if(.not.ASSOCIATED(ELP%SEP15)) ALLOCATE(ELP%SEP15)
       if(.not.ASSOCIATED(ELP%VOLT)) ALLOCATE(ELP%VOLT)
       if(.not.ASSOCIATED(ELP%PHAS)) ALLOCATE(ELP%PHAS)
       ELP%VOLT = EL%VOLT
       ELP%PHAS = EL%PHAS
       CALL SETFAMILY(ELP)
    ENDIF

    !    IF(EL%KIND==KINDUSER1) THEN         !
    !       CALL SETFAMILY(ELP)
    !       CALL COPY(EL%U1,ELP%U1)
    !    ENDIF

    !    IF(EL%KIND==KINDUSER2) THEN         !
    !       CALL SETFAMILY(ELP)
    !       CALL COPY(EL%U2,ELP%U2)
    !    ENDIF

    IF(EL%KIND==KINDWIGGLER) THEN         !
       CALL SETFAMILY(ELP)
       CALL COPY(EL%WI,ELP%WI)
    ENDIF

    IF(ASSOCIATED(EL%RAMP)) THEN         !
       CALL COPY_RAMPING(EL%RAMP,ELP%RAMP)
    ENDIF    
    
    IF(EL%KIND==KINDPA) THEN         !
       CALL SETFAMILY(ELP,t=EL%PA%B) !,EL%PA%angc,EL%PA%xc,EL%PA%dc,EL%PA%h)   !,EL%PA%ax,EL%PA%ay)
       CALL COPY(EL%PA,ELP%PA)
    ENDIF
    IF(EL%KIND==KINDABELL) THEN 
       M_ABELL=EL%ab%M
       N_ABELL=EL%ab%N        !
       CALL SETFAMILY(ELP)  !,EL%PA%angc,EL%PA%xc,EL%PA%dc,EL%PA%h)  !,EL%PA%ax,EL%PA%ay)
       CALL COPY(EL%ab,ELP%ab)
    ENDIF
    !    IF(ASSOCIATED(EL%PARENT_FIBRE))        then
    !       ELP%PARENT_FIBRE=>EL%PARENT_FIBRE
    !    ENDIF


       IF(ASSOCIATED(EL%backWARD))        then
         if(associated(elp%backWARD)) then
          call kill(ELp%backWARD)     
          DEALLOCATE(ELp%backWARD)
         endif
         allocate(ELp%backWARD(3))
         do i=1,3
         call alloc_tree(ELp%backWARD(i),EL%backWARD(i)%n,EL%backWARD(i)%np)
         enddo
         call COPY_TREE_N(EL%backWARD,ELp%backWARD)

       ENDIF
 

       IF(ASSOCIATED(EL%forward))        then
         if(associated(elp%forward)) then
          call kill(ELp%forward)     
          DEALLOCATE(ELp%forward)
         endif
         allocate(ELp%forward(3))
         do i=1,3
         call alloc_tree(ELp%forward(i),EL%forward(i)%n,EL%forward(i)%np)
         enddo
         call COPY_TREE_N(EL%forward,ELp%forward)

       ENDIF


       ELp%skip_ptc_f=EL%skip_ptc_f
       ELp%skip_ptc_b=EL%skip_ptc_b
       elp%do1mapf=el%do1mapf
       elp%do1mapb=el%do1mapb
          ELp%usef= EL%usef
          ELp%useb= EL%useb
  END SUBROUTINE copy_el_el


  SUBROUTINE reset31(ELP)
    IMPLICIT NONE
    TYPE(ELEMENTP),INTENT(inOUT)::  ELP
    INTEGER I

    ELP%knob=.FALSE.

    CALL resetpoly_R31(ELP%L)         ! SHARED BY EVERYONE
    CALL resetpoly_R31(ELP%FINT(1))         ! SHARED BY EVERYONE
    CALL resetpoly_R31(ELP%FINT(2))         ! SHARED BY EVERYONE
    CALL resetpoly_R31(ELP%HGAP(1))         ! SHARED BY EVERYONE
    CALL resetpoly_R31(ELP%HGAP(2))         ! SHARED BY EVERYONE
    CALL resetpoly_R31(ELP%H1)         ! SHARED BY EVERYONE
    CALL resetpoly_R31(ELP%H2)         ! SHARED BY EVERYONE
    CALL resetpoly_R31(ELP%VA)         ! SHARED BY EVERYONE
    CALL resetpoly_R31(ELP%VS)         ! SHARED BY EVERYONE
    if(associated(ELP%theta_ac)) CALL resetpoly_R31(ELP%theta_ac)         ! SHARED BY EVERYONE
    if(associated(ELP%a_ac)) CALL resetpoly_R31(ELP%a_ac)         ! SHARED BY EVERYONE
    if(associated(ELP%DC_ac)) CALL resetpoly_R31(ELP%DC_ac)         ! SHARED BY EVERYONE
    if(associated(ELP%D_ac)) then
       CALL resetpoly_R31(ELP%D_ac)         ! SHARED BY EVERYONE
       IF(ELP%P%NMUL>0) THEN             ! SHARED BY A LOT
          DO I=1,ELP%P%NMUL
             CALL resetpoly_R31(ELP%d_AN(I))
             CALL resetpoly_R31(ELP%d_BN(I))
             CALL resetpoly_R31(ELP%d0_AN(I))
             CALL resetpoly_R31(ELP%d0_BN(I))
          ENDDO
       ENDIF
    endif
    IF(ELP%P%NMUL>0) THEN             ! SHARED BY A LOT
       DO I=1,ELP%P%NMUL
          CALL resetpoly_R31(ELP%AN(I))
          CALL resetpoly_R31(ELP%BN(I))
       ENDDO
    ENDIF

    IF(ELP%KIND==KIND10) THEN
        IF(ELP%ELECTRIC) THEN
           DO I=1,SIZE(ELP%tp10%BE)
              CALL resetpoly_R31(ELP%tp10%AE(I))
              CALL resetpoly_R31(ELP%tp10%BE(I))
           ENDDO
        ENDIF
    ENDIF

    IF(ELP%KIND==KIND4) THEN
       CALL resetpoly_R31(ELP%VOLT)
       CALL resetpoly_R31(ELP%FREQ )
       CALL resetpoly_R31(ELP%PHAS )
       DO I=1,ELP%C4%NF
          CALL resetpoly_R31(ELP%C4%F(I))
          CALL resetpoly_R31(ELP%C4%PH(I))
       ENDDO
       CALL resetpoly_R31(ELP%C4%A )
       CALL resetpoly_R31(ELP%C4%R )

       !      CALL resetpoly_R31(ELP%P0C )
    ENDIF

    IF(ELP%KIND==KIND3) THEN
       CALL resetpoly_R31(ELP%K3%hf)
       CALL resetpoly_R31(ELP%K3%vf)
       CALL resetpoly_R31(ELP%K3%thin_h_foc)
       CALL resetpoly_R31(ELP%K3%thin_v_foc)
       CALL resetpoly_R31(ELP%K3%thin_h_angle )
       CALL resetpoly_R31(ELP%K3%thin_v_angle)
       CALL resetpoly_R31(ELP%B_SOL)
    ENDIF

    IF(ELP%KIND==KIND21) THEN
       CALL resetpoly_R31(ELP%VOLT)
       CALL resetpoly_R31(ELP%FREQ )
       CALL resetpoly_R31(ELP%PHAS )
       CALL resetpoly_R31(ELP%CAV21%PSI )
       CALL resetpoly_R31(ELP%CAV21%DVDS )
       CALL resetpoly_R31(ELP%CAV21%DPHAS )
    ENDIF

    IF(ELP%KIND==KIND22) THEN
       CALL resetpoly_R31(ELP%FREQ )
       CALL resetpoly_R31(ELP%PHAS )
    ENDIF

    IF(ELP%KIND==KIND15) THEN          ! NEW 2002.11.16
       CALL resetpoly_R31(ELP%VOLT)
       CALL resetpoly_R31(ELP%PHAS )
    ENDIF

    IF(ELP%KIND==KIND5) THEN
       CALL resetpoly_R31(ELP%B_SOL)
    ENDIF



    !    IF(ELP%KIND==KINDUSER1) THEN
    !       CALL reset_U1(ELP%U1)
    !    ENDIF

    !    IF(ELP%KIND==KINDUSER2) THEN
    !       CALL reset_U2(ELP%U2)
    !    ENDIF

    IF(ELP%KIND==KINDWIGGLER) THEN
       CALL reset_WI(ELP%WI)
    ENDIF

    IF(ELP%KIND==KINDPA) THEN
       CALL reset_PA(ELP%PA)
    ENDIF

    IF(ELP%KIND==KINDabell) THEN
       CALL reset_abell(ELP%ab)
    ENDIF


  END SUBROUTINE reset31

  SUBROUTINE  find_energy(t,KINETIC,ENERGY,P0C,BRHO,beta0,gamma)
    implicit none
    type(work) ,INTENT(INout):: t
    real(dp) XMC2,cl,CU,ERG,beta0i,GAMMA0,GAMMA2,CON
    logical(lp) PROTON
    real(dp) KINETIC1,ENERGY1,P0C1,BRHO1,beta01,gamma1   !  private here
    real(dp), optional ::   KINETIC,ENERGY,P0C,BRHO,beta0,gamma   !  private here
    real(dp)  gamma0I,gamBET  ! private here

    gamma1=0.0_dp
    kinetic1=0.0_dp
    ENERGY1=0.0_dp
    beta01=0.0_dp
    brho1=0.0_dp
    p0c1=0.0_dp
    if(present(gamma)) gamma1=-gamma
    if(present(KINETIC)) kinetic1=-kinetic
    if(present(energy))  energy1=-energy
    if(present(BETa0))   BETa01=-BETa0
    if(present(brho) )    brho1=-brho
    if(present(p0c) )    p0c1=-p0c

    PROTON=.NOT.ELECTRON
    cl=(clight/1e8_dp)
    CU=55.0_dp/24.0_dp/SQRT(3.0_dp)

    if(electron) then
       XMC2=muon*pmae
    elseif(proton) then
       XMC2=pmap
    endif
    if(ENERGY1<0) then
       ENERGY1=-ENERGY1
       erg=ENERGY1
       p0c1=SQRT(erg**2-xmc2**2)
    endif
    if(kinetic1<0) then
       kinetic1=-kinetic1
       erg=kinetic1+xmc2
       p0c1=SQRT(erg**2-xmc2**2)
    endif
    if(brho1<0) then
       brho1=-brho1
       p0c1=SQRT(brho1**2*(cl/10.0_dp)**2)
    endif
    if(BETa01<0) then
       BETa01=-BETa01
       p0c1=xmc2*BETa01/SQRT(1.0_dp-BETa01**2)
    endif

    if(p0c1<0) then
       p0c1=-p0c1
    endif

    if(gamma1<0) then
       gamma1=-gamma1
       erg=gamma1*xmc2
       p0c1=sqrt(erg**2-XMC2**2)
    endif

    erg=SQRT(p0c1**2+XMC2**2)
    kinetic1=ERG-xmc2
    BETa01=SQRT(kinetic1**2+2.0_dp*kinetic1*XMC2)/erg
    beta0i=1.0_dp/BETa01
    GAMMA0=erg/XMC2
 
    CON=3.0_dp*CU*CGAM*HBC/2.0_dp*TWOPII/pmae**3
!    CON=3.0_dp*CU*CGAM*HBC/2.0_dp*TWOPII/XMC2**3
    CRAD=CGAM*TWOPII   !*ERG**3
    CFLUC=CON  !*ERG**5
    GAMMA2=erg**2/XMC2**2
    brho1=SQRT(ERG**2-XMC2**2)*10.0_dp/cl
    if(verbose) then
       write(6,*) ' p0c = ',p0c1
       write(6,*)' GAMMA0 = ',SQRT(GAMMA2)
       write(6,*)' BRHO = ',brho1
      write(6,*)"CRAD AND CFLUC ", crad ,CFLUC
    endif
 
    gamma0I=XMC2*BETa01/p0c1
    GAMBET=(XMC2/p0c1)**2

    t%kinetic=kinetic1
    t%energy =ERG
    t%BETa0=BETa01
    t%BRHO=brho1
    t%p0c=p0c1
    t%gamma0I=gamma0I
    t%gambet=gambet
    t%mass=xmc2


  END SUBROUTINE find_energy

  subroutine put_aperture_el(el,kind,r,x,y,dx,dy)
    implicit none
    real(dp),intent(in):: r(2),x,y,dx,dy
    integer,intent(in):: kind
    type(element),intent(inout):: el

    if(.not.associated(el%p%aperture)) call alloc(el%p%aperture)
    el%p%aperture%dx=dx
    el%p%aperture%dy=dy
    el%p%aperture%x=x
    el%p%aperture%y=y
    el%p%aperture%r=r
    el%p%aperture%kind=kind
  end  subroutine put_aperture_el

  subroutine put_aperture_elp(el,kind,r,x,y,dx,dy)
    implicit none
    real(dp),intent(in):: r(2),x,y,dx,dy
    integer,intent(in):: kind
    type(elementp),intent(inout):: el

    if(.not.associated(el%p%aperture)) call alloc(el%p%aperture)
    el%p%aperture%dx=dx
    el%p%aperture%dy=dy
    el%p%aperture%x=x
    el%p%aperture%y=y
    el%p%aperture%r=r
    el%p%aperture%kind=kind
  end  subroutine put_aperture_elp

  subroutine remove_aperture_el(el)
    implicit none
    type(element),intent(inout):: el

    if(associated(el%p%aperture)) then
       CALL kill(el%p%APERTURE)
       DEALLOCATE(el%p%APERTURE);
    endif
  end  subroutine remove_aperture_el

  subroutine remove_aperture_elp(el)
    implicit none
    type(elementp),intent(inout):: el

    if(associated(el%p%aperture)) then
       CALL kill(el%p%APERTURE)
       DEALLOCATE(el%p%APERTURE);
    endif
  end  subroutine remove_aperture_elp


  SUBROUTINE decode_element(EL)
    IMPLICIT NONE
    TYPE(ELEMENT),INTENT(INOUT):: EL

    SELECT CASE(EL%KIND)
    CASE(KIND0)
       print*,"this is KIND0"
    case(KIND1)
       print*,"this is KIND1"
    case(KIND2)
       print*,"this is KIND2"
    case(KIND3)
       print*,"this is KIND3"
    case(KIND4)
       print*,"this is KIND4"
    case(KIND5)
       print*,"this is KIND5"
    case(KIND6)
       print*,"this is KIND6"
    case(KIND7)
       print*,"this is KIND7"
    case(KIND8)
       print*,"this is KIND8"
    case(KIND9)
       print*,"this is KIND9"
    case(KIND10)
       print*,"this is KIND11"
    CASE(KIND11)
       print*,"this is KIND12"
    CASE(KIND12)
       print*,"this is KIND13"
    CASE(KIND13)
       print*,"this is KIND14"
    CASE(KIND14)
       print*,"this is KIND11"
    CASE(KIND15)
       print*,"this is KIND15"
    CASE(KIND16)
       print*,"this is KIND16"
    CASE(KIND18)
       print*,"this is KIND18"
    CASE(KIND19)
       print*,"this is KIND19"
    CASE(KIND20)
       print*,"this is KIND20"
    CASE(KIND21)
       print*,"this is KIND21"
    CASE(KIND22)
       print*,"this is KIND22"
    case(KINDWIGGLER)
       print*,"this is KINDWIGGLER"
    case(KINDPA)
       print*,"this is KINDPA"
    case(kindsuperdrift)
       print*,"this is KINDSUPERDRIFT"
    case(KINDABELL)
       print*,"this is KINDABELL"

    case default
 
       write(6,'(1x,i4,a21)') el%kind," not supported decode_element"
       ! call !write_e(0)
    END SELECT
     
    
  end SUBROUTINE decode_element

END MODULE S_DEF_ELEMENT
