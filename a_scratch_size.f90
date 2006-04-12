!The Full Polymorphic Package
!Copyright (C) Etienne Forest and Frank Schmidt

! This library is free software; you can redistribute it and/or
! modify it under the terms of the The Clarified Artistic License.
! However the source files C_dabnew.f90 and D_Lielib.f90 are derived from
! the LBNL versions of Berz's DA-package and Forest's analysis package.
! Their distribution and commercial usage are governed by the laws of
! the United States of America and more specifically by the USA Department
! of Energy. The above license may be partly or totally void with respect to
! these two files.

! All the other files were created as part of Forest's work as an employee
! of the Japanese Ministry of Culture and Education and are, to best of our
! knowledge, his intellectual property under Japanese Law.

! THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES,
! INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF MERCHANTIBILITY AND
! FITNESS FOR A PARTICULAR PURPOSE.


module precision_constants
  implicit none
  public


  !Precision
  integer,parameter::nlp=16
  integer,parameter::vp=16
  integer,parameter::lp=4
  integer,parameter::sp=kind(1e0)
  integer,parameter::dp=selected_real_kind(2*precision(1e0_sp))
  !  linked list for scratch variables of TPSA
  integer  :: newscheme_max =200
  !Numbers double
  real(dp),parameter::zero=0e0_dp,one=1e0_dp,two=2e0_dp,three=3e0_dp,four=4e0_dp,five=5e0_dp

  real(dp),parameter::six=6e0_dp,seven=7e0_dp,eight=8e0_dp,nine=9e0_dp,ten=10e0_dp
  real(dp),parameter::eleven=11e0_dp,twelve=12e0_dp,c_14=14e0_dp,c_15=15e0_dp
  real(dp),parameter::c_16=16e0_dp,c_20=20e0_dp,c_24=24e0_dp,c_27=27e0_dp,c_30=30e0_dp
  real(dp),parameter::c_32=32e0_dp,c_40=40e0_dp,c_48=48e0_dp,c_50=50e0_dp
  real(dp),parameter::c_55=55e0_dp,c_72=72e0_dp,c_80=80e0_dp,c_85=85e0_dp
  real(dp),parameter::c_90=90e0_dp,c_100=100e0_dp,c_120=120e0_dp,c_160=160e0_dp
  real(dp),parameter::c_180=180e0_dp,c_1024=1024e0_dp

  real(dp),parameter::half=0.5e0_dp,c_0_8=0.8e0_dp,c_0_1=0.1e0_dp,c_0_75=0.75e0_dp
  real(dp),parameter::c_0_4375=0.4375e0_dp,c_4d_1=0.4e0_dp,c_0_125=0.125e0_dp
  real(dp),parameter::c_1d_2=1e-2_dp,c_1d_3=1e-3_dp,c_0_0001=1e-04_dp
  real(dp),parameter::c_0_25d_3=0.25e-03_dp,c_0_25=0.25e0_dp,c_0_5d_3=0.5e-03_dp
  real(dp),parameter::c_1d_8=1e-8_dp,c_1_2=1.2e0_dp,c_0_2=0.2e0_dp
  real(dp),parameter::c_1d3=1e3_dp,c_1d5=1e5_dp,c_1d6=1e6_dp
  real(dp),parameter::c_1d9=1e9_dp,c_111110=1.1111e5_dp,c_2d5=2e5_dp
  real(dp),parameter::c_1d8=1e8_dp,c_1d10=1e10_dp,c_1d4=1e4_dp
  real(dp),parameter::c_1d30=1e30_dp,c_1d36=1e36_dp
  real(dp),parameter::c_221=221e0_dp,c_981=981e0_dp,c_867=867e0_dp,c_102=102e0_dp
  real(dp),parameter::c_183=183e0_dp,c_678=678e0_dp,c_472=472e0_dp,c_66=66e0_dp
  real(dp),parameter::c_716=716e0_dp,c_2079=2079e0_dp,c_1002=1002e0_dp,c_834=834e0_dp
  real(dp),parameter::c_454=454e0_dp,c_82=82e0_dp,c_41=41e0_dp,c_216=216e0_dp,c_272=272e0_dp
  real(dp),parameter::c_840=840e0_dp

  real(dp),parameter::c_1d_5=1e-05_dp

  real(dp),TARGET :: hyperbolic_aperture=ten
  !Numbers single
  real(sp),parameter::c_1e7=1e7_sp,c_0e0=0e0_sp
  logical(LP) :: mytrue=.true.,myfalse=.false.
  !Mathematical Constants
  real(dp),PARAMETER::pi=3.141592653589793238462643383279502e0_dp,twopi=two*pi,pih=pi*half
  real(dp),PARAMETER::twopii=one/twopi,pil=pih-c_4d_1,pim=pih+c_4d_1
  real(dp),PARAMETER::rpi4=1.772453850905516027298167e0_dp/two
  real(dp),parameter::RAD_TO_DEG_=c_180/pi,DEG_TO_RAD_=pi/c_180

  !Physical Constants
  real(dp),PARAMETER::pmae=5.10998902e-04_dp
  ![GeV]
  real(dp),PARAMETER::pmap=0.938271998e0_dp
  ![GeV]
  real(dp),PARAMETER::CLIGHT=2.99792458e8_dp
  ![m/s]
  real(dp),PARAMETER::hbar=6.58211889e-25_dp
  ![GeV*s]
  real(dp),PARAMETER::qelect=1.602176462e-19_dp
  ![A*s]
  real(dp),PARAMETER::eps_0=8.854187817e-12_dp
  ![A*S/V*m]
  real(dp),PARAMETER::class_e_radius=qelect/four/pi/eps_0/pmae/c_1d9
  ![m]
  real(dp),PARAMETER::CGAM=four*pi/three*class_e_radius/pmae**3
  ![m/Gev^3] old: 8.846056192e-05_dp
  real(dp),PARAMETER::HBC=hbar*CLIGHT
  ![GeV*m] old: HBC=1.9732858e-16_dp

  !Smallness Parameters
  real(dp),PARAMETER::epsmac=1e-7_dp,c_1d_20=1e-20_dp,c_1d_37=1e-37_dp
  !c_dabnew.f90
  real(dp)::eps=1e-38_dp
  !c_dabnew.f90
  real(dp),PARAMETER::machep=1e-17_dp
  !d_lielib.f90
  real(dp),parameter::PUNY=1e-38_dp
  !e_define_newda.f90: is used in newdamul: could be zero
  real(dp),PARAMETER::EPSdolmac=1e-7_dp
  !f_newda.f90
  real(dp)::EPSdol=1e-37_dp
  !f_newda.f90
  real(dp),PARAMETER::eps_tpsalie=1e-9_dp
  !j_tpsalie.f90
  real(dp),PARAMETER::eps_real_poly=1e-6_dp
  !m_real_polymorph.f90
  real(dp),PARAMETER::eps_rot_mis1=1e-6_dp,eps_rot_mis2=1e-9_dp,epsdif=1e-6_dp
  !Sa_rotation_mis.f90
  real(dp),PARAMETER::eps_extend_poly=1e-10_dp
  !Sb_extend_poly.f90
  real(dp),PARAMETER::eps_fitted=1e-11_dp
  !Sg_0_fitted.f90
  real(dp),PARAMETER::c_1d_11=1e-11_dp
  !Sg_0_fitted.f90
  real(dp)::eps_def_kind=1e-9_dp
  !Sh_DEF_KIND.f90
  real(dp),PARAMETER::deps_tracking=1e-6_dp
  !Sm_tracking.f90
  real(dp),PARAMETER::c_1d_7=1e-7_dp
  !d_lielib.f90,g_newLielib.f90,k_tpsalie_analysis.f90,So_fitting.f90
  real(dp),PARAMETER::c_1d_38=1e-38_dp
  !e_defineda.f90
  real(dp),parameter::c_1d_40=1e-40_dp
  !Sm_tracking.f90
  real(dp),parameter::c_1d_15=1e-15_dp
  !u1_twiss.f90
  real(dp),parameter::c_1d_6=1e-06_dp
  !k_tpsalie_analysis.f90,m_real_polymorph.f90,Sa_rotation_mis.f90
  real(dp),parameter::c_1d_9=1e-09_dp
  !d_lielib.f90,j_tpsalie.f90,Sa_rotation_mis.f90
  real(dp),parameter::c_1d_10=1e-10_dp
  !d_lielib.f90,Sb_extend_poly.f90,Sm_tracking.f90,Sn_mad_like.f90

  real(dp),parameter::c_2_2d_7=2.2e-07_dp,c_1_35d_8=1.35e-08_dp,c_2_2d_8=2.2e-08_dp
  !Sq_userroutine.f90
  real(dp),parameter::c_2_7d_8=2.7e-08_dp
  !Sq_userroutine.f90

  !Magic Numbers
  real(dp),PARAMETER::c_9999_12345=9999.12345e0_dp
  !c_da.f90
  real(dp),parameter::c_0_284496736=0.284496736e0_dp
  !c_da.f90
  real(dp),parameter::c_0_254829592=0.254829592e0_dp
  !c_da.f90
  real(dp),parameter::c_0_3275911=0.3275911e0_dp
  !c_da.f90
  real(dp),parameter::c_1_061405429=1.061405429e0_dp
  !c_da.f90
  real(dp),parameter::c_1_421413741=1.421413741e0_dp
  !c_da.f90
  real(dp),parameter::c_1_453152027=1.453152027e0_dp
  !c_da.f90
  real(dp),parameter::c_720=720e0_dp,c_30240=30240e0_dp
  !d_lielib.f
  real(dp),parameter::c_1209600=1209600e0_dp,c_21772800=21772800e0_dp
  !d_lielib.f
  real(dp),parameter::c_1_17767998417887=1.17767998417887e0_dp
  !Sf_STATUS.f90
  real(dp),parameter::c_0_235573213359357=0.235573213359357e0_dp
  !Sf_STATUS.f90
  real(dp),parameter::c_0_78451361047756=0.78451361047756e0_dp
  !Sf_STATUS.f90
  real(dp),parameter::c_0_999=0.999e0_dp,c_4_2d_2=4.2e-02_dp
  !Sq_userroutine.f90
  real(dp),parameter::c_4_2d_3=4.2e-03_dp,c_6_6d_8=6.6e-08_dp
  !Sq_userroutine.f90
  real(dp),parameter::c_2_2d_3=2.2e-03_dp,c_2_2d_4=2.2e-04_dp
  !Sq_userroutine.f90
  real(dp),parameter::c_0_9=0.9e0_dp,c_2_5=2.5e0_dp,c_0_148=0.148e0_dp
  !Sq_userroutine.f90
  real(dp),parameter::c_0_005=5e-03_dp,c_0_012=1.2e-02_dp,c_1_5=1.5e0_dp
  !Sq_userroutine.f90

  real(dp),parameter::c_0_002=2e-03_dp,c_0_05=5e-02_dp,c_0_216=0.216e0_dp
  !Sq_userroutine.f90
  real(dp),parameter::c_0_7=0.7e0_dp,c_1_2d_5=1.2e-05_dp,c_1d7=1e7_dp
  !Sq_userroutine.f90
  ! Constant Symplectic integrator schemes
  real(dp) YOSK(0:4), YOSD(4)    ! FIRST 6TH ORDER OF YOSHIDA
  real(dp),PARAMETER::AAA=-0.25992104989487316476721060727823e0_dp  ! fourth order integrator
  real(dp),PARAMETER::FD1=half/(one+AAA),FD2=AAA*FD1,FK1=one/(one+AAA),FK2=(AAA-one)*FK1
  ! end of symplectic integrator coefficients

  integer,private,parameter::n_read_max=20,NCAR=120
  private EQUAL_i,EQUAL_Si,EQUAL_r,EQUAL_c,WRITE_G
  private read_d,read_int,read_int_a,read_d_a
  logical(lp), parameter:: my_true=.true.
  logical(lp), parameter:: my_false=.false.
  logical(lp),target :: global_verbose=.true.

  type info_window
     character(3) adv
     integer nc,nr,ni
     character(NCAR) c(n_read_max)
     character(NCAR) Fc
     real(dp) r(n_read_max)
     character(NCAR) FR
     integer i(n_read_max)
     character(NCAR) FI
  end type info_window

  TYPE CONTROL
     ! Da stuff
     real(dp),POINTER :: total_da_size !  in megabytes
     integer,POINTER :: lda_used       !  maximum number of da variables in Berz's
     logical(lp),POINTER  :: OLD    ! = true  = bERZ
     logical(lp),POINTER  :: real_warning  ! = true
     integer,POINTER :: no       ! order of da
     integer,POINTER :: nv       ! number of variables
     integer,POINTER :: nd       ! degrees of freedom
     integer,POINTER :: nd2      ! phase space dimension
     integer,POINTER :: np       ! number of parameters
     integer,POINTER :: ndpt     ! constant energy variable position is different from zero
     integer,POINTER :: NPARA     ! PARAMETER LOCATION IN PTC
     integer,POINTER :: npara_fpp     ! PARAMETER LOCATION IN FPP
     integer,POINTER :: np_pol     ! parameters produced through pol_block
     logical(lp),pointer :: knob
     logical(lp),pointer :: OTHER_PROGRAM
     integer, pointer :: NDPT_OTHER
     logical(lp), pointer :: setknob
     REAL(dp),POINTER     :: da_absolute_aperture  ! in case one tracks with da.
     !

     LOGICAL(lp),POINTER  :: ROOT_CHECK   !=.TRUE. performs check in roots and hyperbolic if true
     LOGICAL(lp),POINTER  :: CHECK_STABLE !=.TRUE. particle status
     LOGICAL(lp),POINTER  :: CHECK_MADX_APERTURE  !=.TRUE. false means particle lost in aperture
     LOGICAL(lp),POINTER  :: APERTURE_FLAG       !=.TRUE. aperture checks globally done
     LOGICAL(lp),POINTER  :: check_iteration       ! checks iteration in fitted magnet for symplectic tracking (not used)
     LOGICAL(lp),POINTER  :: check_interpolate_x       ! check if lost by being outside the interpolation region
     LOGICAL(lp),POINTER  :: check_interpolate_y      ! check if lost by being outside the interpolation region
     LOGICAL(lp),POINTER  :: check_x_min     ! check if lost by aperture fitted now
     LOGICAL(lp),POINTER  :: check_x_max     ! check if lost by aperture fitted now
     LOGICAL(lp),POINTER  :: check_y_min     ! check if lost by aperture fitted now
     LOGICAL(lp),POINTER  :: check_y_max     ! check if lost by aperture fitted now

     LOGICAL(lp),POINTER  :: WATCH_USER     ! FALSE NORMALLY : WATCHES USER FOR FAILING TO CHECK APERTURES

     REAL(dp),POINTER     :: absolute_aperture     !=1e3_dp generic aperture check
     real(dp),POINTER :: hyperbolic_aperture  ! controls crashes in exponentials

     !
     LOGICAL(lp),POINTER   :: NEW_METHOD ! new integration method exists

     ! influence fibre creation
     integer, pointer :: MADTHICK        !
     integer, pointer :: MADTHIN_NORMAL
     integer, pointer :: MADTHIN_SKEW
     INTEGER, pointer::  NSTD,METD       ! number of steps and integration method
     LOGICAL(lp), pointer :: MADLENGTH   ! =.false. rbend crazy length in mad8 as input
     LOGICAL(lp), pointer :: MAD         !=.false. mad definition of multipole for input only
     LOGICAL(lp), pointer :: EXACT_MODEL != .false. exact model used
     logical(lp),pointer :: ALWAYS_EXACTMIS  !=.TRUE. exact formula in tracking used for that element
     INTEGER, pointer:: CAVITY_TOTALPATH ! REAL PILL B0X =1 , FAKE =0
     INTEGER,pointer :: HIGHEST_FRINGE !=2  quadrupole fringe ON IF FRINGE PRESENT
     ! creates a reverse propagator
     INTEGER,pointer ::FIBRE_DIR         !=1 or -1 for reversed
     ! creates a reverse propagator and a reversed ring in combination with above
     logical(lp),pointer ::FIBRE_flip    !=.true.
     !  x_prime true means noncanonical outside magnets. x(5) variables stays the same.
     LOGICAL(lp), pointer :: x_prime != .false.

     ! fill once and never touch again

     INTEGER, pointer :: SECTOR_NMUL_MAX     != 10 maxwell equations is solved to order 10 in exact sectors
     INTEGER, pointer :: SECTOR_NMUL     != 4  MULTIPOLES IN TEAPOT BEND ALLOWED BY DEFAULT
     real(dp), pointer :: wedge_coeff(:)     ! QUAD_KICK IN WEDGE
     logical(lp), pointer :: MAD8_WEDGE      ! QUAD_KICK + FRINGE IF FRINGE IS OUT.





     LOGICAL(lp), pointer:: electron     !  electron if true otherwise proton
     real(dp), pointer :: massfactor     !=one  sets variable muon and electron must be true
     ! global on the fly
     LOGICAL(lp), pointer :: stoch_in_rec != .false. stochastic radiation in straigth elements (global)
     logical(lp),pointer :: FEED_P0C   !=.FALSE.  work takes p0c instead of energy
     logical(lp),pointer :: ALWAYS_EXACT_PATCHING  !=.TRUE. patching done correctly
     ! used to output horror messages
     logical,pointer :: stable_da !=.true.  interrupts DA if check_da is true
     logical,pointer :: check_da  !=.true.
     logical(lp),pointer :: OLD_IMPLEMENTATION_OF_SIXTRACK  !=.true.
     real(dp),pointer :: phase0 ! default phase in cavity
     logical(lp), pointer :: global_verbose
     character*120 message
  end TYPE CONTROL

  type(control) c_

  type(info_window),TARGET:: w_i
  type(info_window),TARGET:: w_ii
  type(info_window),TARGET::  r_i
  type(info_window),POINTER:: W_P
  type(info_window),POINTER:: R_P

  INTERFACE assignment (=)
     MODULE PROCEDURE EQUAL_Si
     MODULE PROCEDURE EQUAL_i
     MODULE PROCEDURE EQUAL_r
     MODULE PROCEDURE EQUAL_c
  end  INTERFACE

  INTERFACE WRITE_I
     MODULE PROCEDURE WRITE_G  ! not private
  END INTERFACE
  INTERFACE WRITE_E
     MODULE PROCEDURE WRITE_G  ! not private
  END INTERFACE

  INTERFACE read
     MODULE PROCEDURE read_d
     MODULE PROCEDURE read_int
     MODULE PROCEDURE read_int_a
     MODULE PROCEDURE read_d_a
  END INTERFACE


contains

  !  SUBROUTINE check_stability(c)
  !    IMPLICIT NONE

  !  end   SUBROUTINE check_stability

  ! Symplectic integrator routines setting coefficients

  SUBROUTINE MAKE_YOSHIDA
    IMPLICIT NONE
    integer i

    YOSK(4)=zero
    YOSK(3)=c_0_78451361047756
    YOSK(2)=c_0_235573213359357
    YOSK(1)=-c_1_17767998417887
    YOSK(0)=one-two*(YOSK(1)+YOSK(2)+YOSK(3))

    do i=4,1,-1
       YOSD(i)=(YOSK(i)+YOSK(i-1))/two
    enddo

    do i=3,0,-1
       YOSK(i+1)=YOSK(I)
    enddo


  END SUBROUTINE MAKE_YOSHIDA


  !  Printing routines for screen

  subroutine  get_ncar(n)
    implicit none
    integer n
    n=ncar
  end subroutine  get_ncar

  SUBROUTINE  EQUAL_Si(S2,S1)
    implicit none
    type (info_window),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1
    integer i
    if(s1==1)then
       s2%adv='NO'
    else
       s2%adv='YES'
    endif
    s2%ni=0
    s2%nC=0
    s2%nR=0
    s2%Fi=' '
    s2%FC=' '
    s2%FR=' '
    do i=1,n_read_max
       s2%i(i)=0
       s2%R(i)=zero
       s2%C(i)=' '
    enddo
  END SUBROUTINE EQUAL_Si

  SUBROUTINE  EQUAL_i(S2,S1)
    implicit none
    type (info_window),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1(:)
    integer n,i
    n=size(s1)
    s2%ni=n
    do i=1,n
       s2%i(i)=s1(i)
    enddo
  END SUBROUTINE EQUAL_i

  SUBROUTINE  EQUAL_r(S2,S1)
    implicit none
    type (info_window),INTENT(inOUT)::S2
    real(dp),INTENT(IN)::S1(:)
    integer n,i
    n=size(s1)
    s2%nr=n
    do i=1,n
       s2%r(i)=s1(i)
    enddo
  END SUBROUTINE EQUAL_r

  SUBROUTINE  EQUAL_c(S2,S1)
    implicit none
    type (info_window),INTENT(inOUT)::S2
    character(*),INTENT(IN)::S1(*)
    integer i
    do i=1,s2%nc
       s2%c(i)=' '
       s2%c(i)=s1(i)
    enddo
  END SUBROUTINE EQUAL_c

  SUBROUTINE WRITE_G(IEX)
    IMPLICIT NONE
    INTEGER, OPTIONAL :: IEX
    INTEGER I,MYPAUSE,IPAUSE
    if(.not.global_verbose) return
    IF(W_P%NC/=0) THEN
       if(W_P%FC/=' ') then
          WRITE(6,W_P%FC,advance=W_P%ADV) (W_P%C(I), I=1,W_P%NC)
       else
          do i=1,W_P%NC
             WRITE(6,*) W_P%C(I)
          enddo
       endif
    ENDIF
    IF(W_P%NI/=0) THEN
       if(W_P%FI/=' ') then
          WRITE(6,W_P%FI,advance=W_P%ADV) (W_P%I(I), I=1,W_P%NI )
       else
          do i=1,W_P%NI
             WRITE(6,*) W_P%I(I)
          enddo
       endif
    ENDIF
    IF(W_P%NR/=0) THEN
       if(W_P%FR/=' ') then
          WRITE(6,W_P%FR,advance=W_P%ADV) (W_P%R(I), I=1,W_P%NR)
       else
          do i=1,W_P%NR
             WRITE(6,*) W_P%R(I)
          enddo
       endif
    ENDIF
    if(W_P%ADV=='NO') then
       WRITE(6,*) " "
    endif
    IF(PRESENT(IEX)) THEN
       if(iex==-1) stop
       IPAUSE=MYPAUSE(IEX)
    ENDIF
  END SUBROUTINE WRITE_G

  SUBROUTINE read_int(IEX)
    IMPLICIT NONE
    integer, intent(inout):: iex
    read(5,*)  iex
  END SUBROUTINE read_int

  SUBROUTINE read_int_a(IEX,n)
    IMPLICIT NONE
    integer, intent(inout):: iex(:)
    integer, intent(in):: n
    integer i
    read(5,*) (iex(i),i=1,n)
  END SUBROUTINE read_int_a

  SUBROUTINE read_d(IEX)
    IMPLICIT NONE
    real(dp), intent(inout)::   iex
    read(5,*) iex
  END SUBROUTINE read_d

  SUBROUTINE read_d_a(IEX,n)
    IMPLICIT NONE
    real(dp), intent(inout)::   iex(:)
    integer, intent(in):: n
    integer i
    read(5,*) (iex(i),i=1,n)
  END SUBROUTINE read_d_a

end module precision_constants


module scratch_size
  implicit none
  public
  integer,parameter::ndumt=10                   !  Number of scratch level
end module scratch_size

module file_handler
  implicit none
  public
  PRIVATE INTFILE,INTFILE_K
  integer,private,parameter::nf=3,MFI=20,MFO=99
  integer,PRIVATE :: winterfile(nf) =(/40,41,42/)
  logical :: myfile(MFI:MFO)
  logical,private :: doit=.true.,ginofile=.false.

  type file_
     logical MF
     ! AIMIN CHANGES FOR MS4.0
     !   logical :: mf=.false.
  end type file_

  type file_K
     logical MF
     ! AIMIN CHANGES FOR MS4.0
     !  logical :: mf=.false.
  end type file_K

  TYPE(FILE_)  NEWFILE
  TYPE(FILE_K)  closeFILE


  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE INTFILE
     MODULE PROCEDURE INTFILE_K
  END  INTERFACE

CONTAINS


  SUBROUTINE  INTFILE(S2,S1)
    ! gino module here
    implicit none
    INTEGER,INTENT(inout)::S2
    type (FILE_),INTENT(in)::S1
    integer i

    if(ginofile) then
       ! put gino compatible code here

    else
       if(doit) then
          call zerofile
          doit=.false.
       endif
       I=MFI
       DO WHILE(myfile(i).AND.I<=MFO)
          I=I+1
       ENDDO
       IF(I==MFO) THEN
          !       w_p=0
          !       w_p%nc=1
          !       w_p%fc='(1(1X,A120))'
          !       w_p%c(1)=  " NO MORE UNITS AVAILABLE: INTFILE "
          !       call write_e(MFO)
       ENDIF
       S2=I
       myfile(I)=.true.
       if(s1%mf) then
          !      w_p=0
          !      w_p%nc=1
          !      w_p%fc='(1(1X,A120))'
          !      write(w_p%c(1),'(1x,L1,1x)')   s1%mf
          !      call write_i
       endif
    endif
  END SUBROUTINE INTFILE


  SUBROUTINE  INTFILE_K(S2,S1)
    implicit none
    INTEGER,INTENT(inOUT)::S2
    type (FILE_K),INTENT(IN)::S1
    integer IPAUSE,MYPAUSES
    character(120) line


    if(ginofile) then
       ! put gino compatible code here

    else
       IF(S2>=MFI.AND.S2<MFO) THEN
          IF(myfile(S2)) THEN
             myfile(S2)=.false.
          ELSE
             WRITE(line,'(a30)') "PROBLEMS WITH UNITS: INTFILE_K"
             WRITE(line(31:120),'(1x,i4,1x,L1)') S2, myfile(S2)
             IPAUSE=MYPAUSES(1,line)
          ENDIF
       ELSE
          WRITE(line(1:30),'(a30)') "PROBLEMS WITH UNITS: INTFILE_K"
          WRITE(line(31:120),'(1x,i4,1x,L1)') S2, myfile(S2)
          IPAUSE=MYPAUSES(2,line)
       ENDIF
       CLOSE(S2)
       S2=-S2
       if(s1%mf) THEN
          WRITE(line,'(a9,L1)') " s1%mf = ",s1%mf
          IPAUSE=MYPAUSES(3,line)
       ENDIF
    endif
  END SUBROUTINE INTFILE_K


  SUBROUTINE  ZEROFILE
    implicit none
    integer ipause, mypause
    integer i

    DO I=MFI,MFO
       myfile(i)=.FALSE.
    ENDDO
    do i=1,nf
       myfile(winterfile(i))=.true.
    enddo
  END SUBROUTINE ZEROFILE

  SUBROUTINE KanalNummer(iff)
    implicit none
    INTEGER, INTENT(OUT) :: iff

    LOGICAL :: opened, exists
    INTEGER :: i

    DO i= 9999, 7, -1
       INQUIRE(UNIT= i, EXIST= exists, OPENED= opened)
       IF (exists .AND. (.NOT. opened)) GOTO 20
    ENDDO
    WRITE (UNIT= *, FMT= *) ' cannot find free unit within the range 7-9999..'
    CALL ReportOpenFiles
    STOP
20  CONTINUE
    iff= I
  END SUBROUTINE KanalNummer

  SUBROUTINE ReportOpenFiles
    implicit none

    LOGICAL :: opened, exists, named
    CHARACTER(LEN= 400) :: name
    integer :: i

    DO i= 9999, 7, -1
       INQUIRE(UNIT= i, EXIST= exists, OPENED= opened)
       IF (exists .AND. opened) THEN
          INQUIRE(UNIT= i, NAMED= named, NAME= name)
          write (*,7010) i, name(:LEN_TRIM(name))
       ENDIF
    ENDDO
7010 FORMAT(' iUnit:',I3,', name: "',A,'"')

  END SUBROUTINE ReportOpenFiles

end module file_handler


module my_own_1D_TPSA
  USE precision_constants
  implicit none
  public
  private input_real_in_my_1D_taylor
  integer :: n_tpsa_exp = 10
  INTEGER, PARAMETER :: N_my_1D_taylor=9  ! SHOULD BE AS ENGE_N

  PRIVATE mul,dmulsc,dscmul
  PRIVATE div,ddivsc,dscdiv,Idivsc
  PRIVATE add,unaryADD,daddsc,dscadd
  PRIVATE subs,unarySUB,dsubsc,dscsub,POW
  PRIVATE  DEXPT,DLOGT,DSQRTT
  PRIVATE  DCOST,DSINT

  ! Definitions
  TYPE my_1D_taylor
     real(dp) a(0:N_my_1D_taylor)

  END TYPE my_1D_taylor

  INTERFACE assignment (=)
     MODULE PROCEDURE input_my_1D_taylor_in_real   !@1 &nbsp; Taylor=real
     MODULE PROCEDURE input_real_in_my_1D_taylor  !@1 &nbsp; real=Taylor
  end  INTERFACE

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE mul   !@1 &nbsp; Taylor * Taylor
     MODULE PROCEDURE dmulsc    !@1 &nbsp; Taylor * Real(dp)
     MODULE PROCEDURE dscmul     !@1 &nbsp;  Real(dp) * Taylor
  END INTERFACE

  INTERFACE OPERATOR (/)
     MODULE PROCEDURE div    !@1 &nbsp; Taylor / Taylor
     MODULE PROCEDURE ddivsc  !@1 &nbsp; Taylor / Real(dp)
     MODULE PROCEDURE dscdiv  !@1 &nbsp; Real(dp) / Taylor
     MODULE PROCEDURE Idivsc  !@1 &nbsp; Taylor / Integer  <font color="#FF0000">&nbsp; &#8594; &nbsp; added because useful in example code</font>
  END INTERFACE

  INTERFACE OPERATOR (+)
     MODULE PROCEDURE add       !@1 &nbsp; Taylor + Taylor
     MODULE PROCEDURE unaryADD  !@1 &nbsp;  +Taylor
     MODULE PROCEDURE daddsc    !@1 &nbsp;  Taylor + Real(dp)
     MODULE PROCEDURE dscadd    !@1 &nbsp;  Real(dp) + Taylor
  END INTERFACE

  INTERFACE OPERATOR (-)
     MODULE PROCEDURE subs      !@1 &nbsp; Taylor - Taylor
     MODULE PROCEDURE unarySUB     !@1 &nbsp;  -Taylor
     MODULE PROCEDURE dsubsc   !@1 &nbsp;  Taylor - Real(dp)
     MODULE PROCEDURE dscsub    !@1 &nbsp;  Real(dp) - Taylor
  END INTERFACE


  INTERFACE OPERATOR (**)
     MODULE PROCEDURE POW    !@1 &nbsp; Taylor ** Integer
  END INTERFACE


  !@3  <p><i><font size="5">&nbsp;Overloading standard procedures </font></i></p>

  INTERFACE exp
     MODULE PROCEDURE DEXPT   !@1 &nbsp; exp(Taylor)
  END INTERFACE

  INTERFACE LOG
     MODULE PROCEDURE DLOGT   !@1 &nbsp; log(Taylor)
  END INTERFACE

  INTERFACE SQRT
     MODULE PROCEDURE DSQRTT    !@1 &nbsp; sqrt(Taylor)
  END INTERFACE

  INTERFACE COS
     MODULE PROCEDURE DCOST    !@1 &nbsp; cos(Taylor)
  END INTERFACE

  INTERFACE SIN
     MODULE PROCEDURE DSINT   !@1 &nbsp; sin(Taylor)
  END INTERFACE



  !@3  <p><i><font size="5">User defined operator</font></i></p>


  ! Destructors and Constructors for my_1D_taylor


contains




  FUNCTION add( S1, S2 )
    implicit none
    TYPE (my_1D_taylor) add
    TYPE (my_1D_taylor), INTENT (IN) :: S1, S2

    add%a=S1%a + S2%a

  END FUNCTION add

  FUNCTION daddsc( S1, sc )
    implicit none
    TYPE (my_1D_taylor) daddsc
    TYPE (my_1D_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    daddsc=s1
    daddsc%a(0)= s1%a(0) + sc

  END FUNCTION daddsc

  FUNCTION dscadd( sc ,  S1)
    implicit none
    TYPE (my_1D_taylor) dscadd
    TYPE (my_1D_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    dscadd=s1
    dscadd%a(0)= s1%a(0) + sc

  END FUNCTION dscadd

  FUNCTION unaryadd( S1 )
    implicit none
    TYPE (my_1D_taylor) unaryadd
    TYPE (my_1D_taylor), INTENT (IN) :: S1

    unaryadd=s1

  END FUNCTION unaryadd


  FUNCTION subs( S1, S2 )
    implicit none
    TYPE (my_1D_taylor) subs
    TYPE (my_1D_taylor), INTENT (IN) :: S1, S2

    subs%a=S1%a - S2%a

  END FUNCTION subs

  FUNCTION unarySUB( S1 )
    implicit none
    TYPE (my_1D_taylor) unarySUB
    TYPE (my_1D_taylor), INTENT (IN) :: S1

    unarySUB%a=-s1%a

  END FUNCTION unarySUB

  FUNCTION dsubsc( S1, sc )
    implicit none
    TYPE (my_1D_taylor) dsubsc
    TYPE (my_1D_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    dsubsc=s1
    dsubsc%a(0)= s1%a(0) - sc

  END FUNCTION dsubsc


  FUNCTION dscsub( sc , S1 )
    implicit none
    TYPE (my_1D_taylor) dscsub
    TYPE (my_1D_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    dscsub=-s1   ! uses unary sub
    dscsub=sc + dscsub   ! uses add

  END FUNCTION dscsub


  FUNCTION mul( S1, S2 )
    implicit none
    TYPE (my_1D_taylor) mul
    TYPE (my_1D_taylor), INTENT (IN) :: S1, S2
    INTEGER I,J
    mul%a=0.0_DP
    DO I=0,N_my_1D_taylor
       DO J=0,N_my_1D_taylor
          IF(I+J>N_my_1D_taylor) CYCLE
          mul%a(I+J)=S1%a(I)*S2%a(J)+ mul%a(I+J)
       ENDDO
    ENDDO

  END FUNCTION mul

  FUNCTION dmulsc( S1, sc )
    implicit none
    TYPE (my_1D_taylor) dmulsc
    TYPE (my_1D_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    dmulsc%a= s1%a*sc

  END FUNCTION dmulsc

  FUNCTION dscmul( sc ,S1 )
    implicit none
    TYPE (my_1D_taylor) dscmul
    TYPE (my_1D_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    dscmul%a= s1%a*sc

  END FUNCTION dscmul

  FUNCTION ddivsc( S1, sc )
    implicit none
    TYPE (my_1D_taylor) ddivsc
    TYPE (my_1D_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    ddivsc%a= s1%a/sc

  END FUNCTION ddivsc

  FUNCTION Idivsc( S1, sc )
    implicit none
    TYPE (my_1D_taylor) Idivsc
    TYPE (my_1D_taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: sc

    Idivsc%a= s1%a/sc

  END FUNCTION Idivsc

  FUNCTION POW( S1,N )
    implicit none
    TYPE (my_1D_taylor) POW , T
    TYPE (my_1D_taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: N
    INTEGER I

    POW=one

    IF(N<=0) THEN
       T=one/S1
       DO I=1,-N
          POW=POW*T
       ENDDO
    ELSE
       DO I=1,N
          POW=POW*S1
       ENDDO

    ENDIF

  END FUNCTION POW

  FUNCTION inv( S1 )
    implicit none
    TYPE (my_1D_taylor) inv,t,TT
    TYPE (my_1D_taylor), INTENT (IN) :: S1
    INTEGER I
    T=S1/S1%A(0)
    T%A(0)=0.D0

    TT=T
    inv=one

    DO I=1,N_my_1D_taylor
       INV=INV-TT
       TT=-TT*T
    ENDDO

    INV=INV/S1%A(0)

  END FUNCTION inv

  FUNCTION DIV( S1, S2 )
    implicit none
    TYPE (my_1D_taylor) DIV
    TYPE (my_1D_taylor), INTENT (IN) :: S1, S2

    DIV=INV(S2)
    DIV=S1*DIV
  END FUNCTION DIV

  FUNCTION dscdiv(  sc , S1)
    implicit none
    TYPE (my_1D_taylor) dscdiv
    TYPE (my_1D_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    dscdiv=INV(S1)
    dscdiv= SC *  dscdiv

  END FUNCTION dscdiv

  !  DEFININING my_1D_taylor=CONSTANT
  subroutine input_real_in_my_1D_taylor( S2, S1 )
    implicit none
    real(dp), INTENT (IN) :: S1
    TYPE (my_1D_taylor), INTENT (inout) :: S2

    S2%a=0.0_dp
    S2%a(0)=s1

  END subroutine input_real_in_my_1D_taylor

  subroutine input_my_1D_taylor_in_real( S2, S1 )
    implicit none
    TYPE (my_1D_taylor), INTENT (IN) :: S1
    real(dp), INTENT (inout) :: S2

    S2=S1%a(0)

  END subroutine input_my_1D_taylor_in_real

  !  DEFININING EXP

  FUNCTION DEXPT( S1 )
    implicit none
    TYPE (my_1D_taylor) DEXPT,t,tt
    TYPE (my_1D_taylor), INTENT (IN) :: S1
    integer I
    T=S1
    T%A(0)=0.0_DP


    DEXPT=1.0_dp
    tt=1.0_dp

    do i=1,N_my_1D_taylor
       tt=tt*t/i
       DEXPT=DEXPT + tt
    enddo

    DEXPT=DEXPT*EXP(S1%A(0))

  END FUNCTION DEXPT


  FUNCTION DLOGT( S1 )
    implicit none
    TYPE (my_1D_taylor) DLOGT,t,TT
    TYPE (my_1D_taylor), INTENT (IN) :: S1
    INTEGER I

    T=S1/S1%A(0)
    T%A(0)=0.0_DP

    TT=T
    do i=1,N_my_1D_taylor
       DLOGT=TT/I+DLOGT
       TT=-TT*T
    ENDDO

    DLOGT=DLOGT+ LOG(S1%A(0))

  END FUNCTION  DLOGT

  FUNCTION DSQRTT( S1 )
    implicit none
    TYPE (my_1D_taylor) DSQRTT,t
    TYPE (my_1D_taylor), INTENT (IN) :: S1
    WRITE(6,*) " MARDE "
    STOP 666
    T=S1/S1%A(0)
    T%A(0)=0.0_DP

    DSQRTT=1.0_DP + T/2.0_DP - T**2/8.0_DP+T**3/16.0_DP
    DSQRTT=DSQRTT* SQRT(S1%A(0))

  END FUNCTION  DSQRTT

  FUNCTION DCOST(S1)
    implicit none
    TYPE (my_1D_taylor) DCOST,t
    TYPE (my_1D_taylor), INTENT (IN) :: S1
    WRITE(6,*) " MARDE "
    STOP 666

    T=S1
    T%A(0)=0.0_DP

    DCOST=COS(S1%A(0))*(1.0_DP-T**2/2.0_DP)-SIN(S1%A(0))*(T-T**3/6.0_DP)

  END FUNCTION  DCOST

  FUNCTION DSINT(S1)
    implicit none
    TYPE (my_1D_taylor) DSINT,t
    TYPE (my_1D_taylor), INTENT (IN) :: S1

    T=S1
    T%A(0)=0.0_DP

    DSINT=SIN(S1%A(0))*(1.0_DP-T**2/2.0_DP)+COS(S1%A(0))*(T-T**3/6.0_DP)

  END FUNCTION  DSINT

end module my_own_1D_TPSA


integer function mypause(i)
  use precision_constants
  implicit none
  !
  ! Replaces obsolescent feature pause
  !
  integer i
  !
  !  write (*,'(A,i6)') ' PAUSE: ',i
  w_p=1
  w_p=(/i/); w_p%fi='(1x,i4)'
  w_p%nc=1
  w_p%c(1)=' ipause=mypause(0)  ';w_p%fc='((A8,1x))'

  ! call write_i
  ! read(*,*) I
  mypause=i
  ! mypause=sqrt(dble(-i))
end function mypause

integer function mypauses(i,string)
  use precision_constants
  implicit none
  !
  ! Replaces obsolescent feature pause
  !
  integer i,l,n
  character(*) string
  !
  !  write (*,'(A,i6)') ' PAUSE: ',i
  w_p=1
  w_p=(/i/); w_p%fi='(1x,i4)'
  w_p%nc=2
  l=len(string)
  call get_ncar(n)
  if(l>n) l=120
  w_p%c(1)=string(1:l)
  w_p%c(2)=' ipause=mypause(0)  ';w_p%fc='((A120,1x,/,a8,1x,))'

  call write_i
  read(*,*) I
  mypauses=i
  mypauses=sqrt(dble(-i))
end function mypauses
