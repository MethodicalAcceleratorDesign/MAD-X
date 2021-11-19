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
  integer,parameter  :: newscheme_max =200 
  integer,private,parameter::n_read_max=20,NCAR=120
  private read_d,read_int,read_int_a,read_d_a
  !Double precision
  real(kind(1d0)) :: doublenum = 0d0
  !Precision
  integer,parameter::nlp=24
  integer,parameter::vp=nlp
  integer,parameter::lp=4
  ! double precision
  integer,parameter::sp=kind(1.e0)
  integer,parameter::dp=selected_real_kind(2*precision(1.e0))
  ! quadrupole precision
  !  integer,parameter::sp=selected_real_kind(2*precision(1.e0))
  !  integer,parameter::dp=selected_real_kind(4*precision(1.e0))
  !Logicals
  logical(lp),parameter:: my_true=.true.
  logical(lp),parameter:: my_false=.false.
  logical(lp),target :: global_verbose=.false.
  logical(lp),target :: no_hyperbolic_in_normal_form=.true.
  integer,private,parameter::nmax=400
  real(dp),private,parameter::tiny=1e-20_dp
  private ludcmp_nr,lubksb_nr
  private ludcmp_nr0,lubksb_nr0
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
  real(dp),parameter::c_1d_2=1e-2_dp,c_1d_3=1e-3_dp,c_0_0001=1e-4_dp
  real(dp),parameter::c_0_25d_3=0.25e-3_dp,c_0_25=0.25e0_dp,c_0_5d_3=0.5e-3_dp
  real(dp),parameter::c_1d_8=1e-8_dp,c_1_2=1.2e0_dp,c_0_2=0.2e0_dp
  real(dp),parameter::c_1d3=1e3_dp,c_1d5=1e5_dp,c_1d6=1e6_dp
  real(dp),parameter::c_1d9=1e9_dp,c_111110=1.1111e5_dp,c_2d5=2e5_dp
  real(dp),parameter::c_1d8=1e8_dp,c_1d10=1e10_dp,c_1d4=1e4_dp
  real(dp),parameter::c_1d30=1e30_dp,c_1d36=1e36_dp
  real(dp),parameter::c_221=221e0_dp,c_981=981e0_dp,c_867=867e0_dp,c_102=102e0_dp
  real(dp),parameter::c_183=183e0_dp,c_678=678e0_dp,c_472=472e0_dp,c_66=66e0_dp
  real(dp),parameter::c_716=716e0_dp,c_2079=2079e0_dp,c_1002=1002e0_dp,c_834=834e0_dp
  real(dp),parameter::c_454=454e0_dp,c_82=82e0_dp,c_41=41e0_dp,c_216=216e0_dp,c_272=272e0_dp
  real(dp),parameter::c_840=840e0_dp,c_360=360e0_dp,c_300=300e0_dp,c_137=137e0_dp
  real(dp),parameter::c_0_28=0.28e0_dp,c_0_31=0.31e0_dp,c_1_8=1.8e0_dp
  real(dp),parameter::c_1d_5=1e-5_dp, c_0_3079=0.3079e0_dp
  !Mathematical Constants
  real(dp),TARGET :: hyperbolic_aperture=10.0_dp
  real(dp),parameter::pi=3.141592653589793238462643383279502e0_dp,twopi=2.0_dp*pi,pih=pi*0.5_dp
  real(dp),parameter::twopii=1.0_dp/twopi,pil=pih-0.4_dp,pim=pih+0.4_dp
  !  real(dp),parameter::rpi4=1.772453850905516027298167e0_dp/two ! rpi4=sqrt(pi)/2
  real(dp)::rpi4
  real(dp),parameter::RAD_TO_DEG_=180.0_dp/pi,DEG_TO_RAD_=pi/180.0_dp
  !Physical Constants
  real(dp),parameter::A_ELECTRON=1.15965218128e-3_dp ! NIST CODATA 2018
   real(dp),parameter::A_MUON=1.16592089e-3_dp        ! NIST CODATA 2018
   real(dp),parameter::A_PROTON=1.79284735e-0_dp      !
   real(dp),parameter::pmaMUON=105.6583755E-3_dp      ! NIST CODATA 2018
   real(dp) :: e_muon = 0.d0, volt_c=1.0e-3_dp, volt_i=1.0_dp
 !  real(dp),parameter:: pmadt = 1.875612793e0_dp    ! sateesh
  !  real(dp),parameter:: pmah3 = 2.808391e0_dp    ! sateesh
  !  real(dp),parameter:: A_dt = -0.142987272e0_dp    ! sateesh
  !  real(dp),parameter:: a_h3 =-4.183963e0_dp    ! sateesh
  logical(lp),  public :: longprint = my_true
  logical(lp),  public :: madxprint = my_true

   real(dp) :: A_particle=A_ELECTRON
   real(dp),parameter::pmae=5.1099895000E-4_dp     ! NIST CODATA 2018
   real(dp),parameter::pmae_amu=5.446170214876E-4_dp  ! NIST CODATA 2018 [GeV] 5.446170214876E-4
   real(dp),parameter::pmap=0.93827208816e0_dp     ! NIST CODATA 2018 [GeV]
   real(dp),parameter::CLIGHT=2.99792458e8_dp       ! exact [m/s]
   real(dp),parameter::hbar=6.582119569e-25_dp    ! NIST CODATA 2018 [GeV*s]
   real(dp),parameter::dhbar=1.054571817e-34_dp   ! NIST CODATA 2018 [J*s]
   real(dp),parameter::qelect=1.602176634e-19_dp ! NIST CODATA 2018 [A*s]
   real(dp),parameter::eps_0=8.854187817e-12_dp     ! exact [A*S/V*m]
   real(dp),parameter::class_e_radius=qelect/4.0_dp/pi/eps_0/pmae/1e9_dp ![m]
   real(dp),parameter::CGAM=4.0_dp*pi/3.0_dp*class_e_radius/pmae**3 ![m/Gev^3] old: 8.846056192e-5_dp
   real(dp),parameter::HBC=hbar*CLIGHT              ![GeV*m] old: HBC=1.9732858e-16_dp
  real(dp),parameter:: alpha_fine= 137.035999084_dp
   real(dp),parameter::class_c_radius=qelect/4.0_dp/pi/eps_0/1e9_dp ![m]
   real(dp),parameter::CGAM0=4.0_dp*pi/3.0_dp*class_c_radius   !/pmae**4 ![m/Gev^3] old: 8.846056192e-5_dp

  ![GeV*m] old: HBC=1.9732858e-16_dp
  !Smallness Parameters
   !Smallness Parameters
  real(dp),parameter::epsmac=1e-7_dp,c_1d_20=1e-20_dp,c_1d_37=1e-37_dp
  real(dp) :: epsflo=1.e-10_dp
  real(dp),parameter::mybig=1e38_dp
  real(dp),parameter::machep=1e-17_dp
  real(dp),parameter::PUNY=1e-38_dp
  real(dp),parameter::EPSdolmac=1e-7_dp
  real(dp),parameter::eps_tpsalie=1e-9_dp
  real(dp),parameter::eps_real_poly=1e-6_dp
  real(dp),parameter::eps_rot_mis1=1e-6_dp,eps_rot_mis2=1e-9_dp,epsdif=1e-6_dp
  real(dp),parameter::eps_extend_poly=1e-10_dp
  real(dp),parameter::eps_fitted=1e-11_dp
  real(dp),parameter::c_1d_11=1e-11_dp
  real(dp),parameter::c_1d_16=1e-16_dp
  real(dp),parameter::eps_def_kind=1e-9_dp
  real(dp),parameter::deps_tracking=1e-6_dp
  real(dp),parameter::c_1d_7=1e-7_dp
  real(dp),parameter::c_1d_38=1e-38_dp
  real(dp),parameter::c_1d_40=1e-40_dp
  real(dp),parameter::c_1d_15=1e-15_dp
  real(dp),parameter::c_1d_6=1e-6_dp
  real(dp),parameter::c_1d_9=1e-9_dp
  real(dp),parameter::c_1d_10=1e-10_dp
  real(dp),parameter::c_2_2d_7=2.2e-7_dp,c_1_35d_8=1.35e-8_dp,c_2_2d_8=2.2e-8_dp
  real(dp),parameter::c_2_7d_8=2.7e-8_dp
  !Magic Numbers
  real(dp),parameter::c_9999_12345=9999.12345e0_dp
  real(dp),parameter::c_0_284496736=0.284496736e0_dp
  real(dp),parameter::c_0_254829592=0.254829592e0_dp
  real(dp),parameter::c_0_3275911=0.3275911e0_dp
  real(dp),parameter::c_1_061405429=1.061405429e0_dp
  real(dp),parameter::c_1_421413741=1.421413741e0_dp
  real(dp),parameter::c_1_453152027=1.453152027e0_dp
  real(dp),parameter::c_720=720e0_dp,c_30240=30240e0_dp
  real(dp),parameter::c_1209600=1209600e0_dp,c_21772800=21772800e0_dp
  real(dp),parameter::c_1_17767998417887=1.17767998417887e0_dp
  real(dp),parameter::c_0_235573213359357=0.235573213359357e0_dp
  real(dp),parameter::c_0_78451361047756=0.78451361047756e0_dp
  real(dp),parameter::c_0_999=0.999e0_dp,c_4_2d_2=4.2e-2_dp
  real(dp),parameter::c_4_2d_3=4.2e-3_dp,c_6_6d_8=6.6e-8_dp
  real(dp),parameter::c_2_2d_3=2.2e-3_dp,c_2_2d_4=2.2e-4_dp
  real(dp),parameter::c_0_9=0.9e0_dp,c_2_5=2.5e0_dp,c_0_148=0.148e0_dp
  real(dp),parameter::c_0_005=5e-3_dp,c_0_012=1.2e-2_dp,c_1_5=1.5e0_dp
  real(dp),parameter::c_0_002=2e-3_dp,c_0_05=5e-2_dp,c_0_216=0.216e0_dp
  real(dp),parameter::c_0_7=0.7e0_dp,c_1_2d_5=1.2e-5_dp,c_1d7=1e7_dp
  real(dp),parameter:: suntao=1.6021766208e-19_dp/299792458.0_dp/9.10938356e-31_dp
  ! Constant Symplectic integrator schemes
  real(dp) YOSK(0:4), YOSD(4)    ! FIRST 6TH ORDER OF YOSHIDA
  real(dp) wyosh(0:7),wyoshid(0:15),wyoshik(15)    ! FIRST 8TH ORDER OF YOSHIDA
  real(dp) ck8(11),bk8(11),ak8(11,11)
  real(dp),parameter::AAA=-0.25992104989487316476721060727823e0_dp  ! fourth order integrator
  real(dp),parameter::FD1=0.5_dp/(1.0_dp+AAA),FD2=AAA*FD1,FK1=1.0_dp/(1.0_dp+AAA),FK2=(AAA-1.0_dp)*FK1
  ! end of symplectic integrator coefficients
  !Initialized numbers
  real(dp)::eps=1e-38_dp
  real(dp)::EPSdol=1e-37_dp
  LOGICAL(lp),target  :: s_aperture_CHECK=.TRUE.
  LOGICAL(lp),TARGET  :: ROOT_CHECK=.TRUE.
  LOGICAL(lp),TARGET  :: CHECK_STABLE=.TRUE.
  LOGICAL(lp),TARGET  :: WATCH_USER=.FALSE.
  !  LOGICAL(lp) :: ALLOW_TRACKING=.true.
  LOGICAL(lp),TARGET  :: CHECK_MADX_APERTURE=.TRUE.
  LOGICAL(lp),TARGET  :: APERTURE_FLAG=.true.

  REAL(dp),TARGET   :: absolute_aperture=1.0_dp, t_aperture =1.d6
  integer,TARGET :: wherelost=0
  logical(lp),TARGET :: stable_da =.true.
  logical(lp),TARGET :: check_da =.true.
  logical(lp),TARGET :: print_frame =.true.
  logical(lp),TARGET :: sixtrack_compatible =.false.
  integer ,target ::  spin_normal_position=2
  real(dp),target ::  da_absolute_aperture=1e6_dp
  real(dp),pointer :: crash => null()
  INTEGER,  TARGET :: NPARA_original
  logical  :: default_tpsa=.false.,drawframeblack=.false.    ! drawframeblack for gino
  logical, target :: lingyun_yang=.false.
  integer, target :: last_tpsa=0
  integer, target :: c_last_tpsa=0
  integer :: mf_herd=0
  character*255 :: print_herd="PRINT_HERD.TXT"
  character*255 :: initial_setting="FINAL_SETTINGS.TXT"
  character*255 :: final_setting="FINAL_SETTINGS.TXT"
  character*255 :: def_orbit_node="def_orbit_node.txt"
  character*255 :: file_block_name="noprint"
  real(dp) :: lmax=1.e38_dp
  logical(lp) :: printdainfo=my_false
  integer   lielib_print(15)
  DATA lielib_print /0,0,0,1,0,0,0,0,0,0,0,1,0,1,1/
  integer :: SECTOR_NMUL_MAX=22
  INTEGER, target :: SECTOR_NMUL = 11
 
!  integer, parameter :: no_e=5  !  electric 
  logical(lp) :: use_complex_in_ptc=.false.
  logical(lp) :: change_sector=my_true
  real(dp) :: xlost(6)=0.0_dp
  integer :: limit_int0(2) =(/4,18/)
  character(1024) :: messagelost
  integer, target :: ndpt_bmad = 0, only2d =0, addclock=0
  integer,TARGET :: HIGHEST_FRINGE=2
  logical :: use_quaternion = .false.
  logical :: use_tpsa = .false.
  logical :: conversion_xprime_in_abell=.true.
    logical(lp) :: use_bmad_units=.false.,inside_bmad=.false.
  character(18) ::         format1="(1(1x,g23.16,1x))"
  character(18) ::         format2="(2(1x,g23.16,1x))"
  character(18) ::         format3="(3(1x,g23.16,1x))"
  character(18) ::         format4="(4(1x,g23.16,1x))"
  character(18) ::         format5="(5(1x,g23.16,1x))"
  character(18) ::         format6="(6(1x,g23.16,1x))"
  character(18) ::         format7="(7(1x,g23.16,1x))"
  character(18) ::         format8="(8(1x,g23.16,1x))"
  character(18) ::         format9="(9(1x,g23.16,1x))"
  character(19) ::         format10="(10(1x,g23.16,1x))"
  !  logical(lp) :: fixed_found
  !  lielib_print(1)=1   lieinit prints info
  !  lielib_print(2)=1   expflo warning if no convergence
  !  lielib_print(3)=1   Shows details in flofacg
  !  lielib_print(4)=1   prints thin layout information
  !  lielib_print(5)=1  order in orbital normal form
  !  lielib_print(6)=1  symplectic condition
  !  lielib_print(7)=-1  go manual in normal form  (use auto command in fpp)
  !  lielib_print(8)=-1  To use nplane from FPP normalform%plane
  !  lielib_print(9)=1  print in checksymp(s1,norm) in j_tpsalie.f90
  !  lielib_print(10)=1  print radation gain 
  !  lielib_print(11)=1  print warning about Teng-Edwards
  !  lielib_print(12)=1  print info in make_node_layout
  !  lielib_print(13)=1  print info of normal form kernel into file kernel.txt and kernel_spin.txt
  !  lielib_print(14)=1  print info about recutting
  !  lielib_print(15)=1  print info during flat file reading and printing

  INTERFACE read
     MODULE PROCEDURE read_d
     MODULE PROCEDURE read_int
     MODULE PROCEDURE read_int_a
     MODULE PROCEDURE read_d_a
  END INTERFACE

  type CONTROL
     ! Da stuff
     real(dp),pointer :: total_da_size => null() !  in megabytes
     integer,pointer :: lda_used => null()       !  maximum number of da variables in Berz's
     logical(lp),pointer  :: OLD => null()    ! = true  = bERZ
     logical(lp),pointer  :: real_warning => null()  ! = true
     integer,pointer :: no  => null()      ! order of da
     integer,pointer :: nv  => null()      ! number of variables
     integer,pointer :: nd  => null()      ! degrees of freedom
     integer,pointer :: nd2 => null()      ! phase space dimension
     integer,pointer :: np  => null()      ! number of parameters in fpp
     integer,pointer :: nspin => null()       ! number of spin variables (0 or 3)
     integer,pointer :: SPIN_pos => null()       ! position of spin variables (0 or 3)
     integer,pointer :: ndpt     => null() ! constant energy variable position is different from zero
     integer,pointer ::ndptb  => null()  
     integer,pointer :: NPARA    => null() ! PARAMETER LOCATION IN PTC in fpp
     integer,pointer :: npara_fpp=> null()     ! PARAMETER LOCATION IN FPP or PTC
     integer,pointer :: np_pol   => null()  ! parameters produced through pol_block
     integer,pointer :: nd2t   => null()  ! harmonic planes minus clocks
     integer,pointer :: nd2harm   => null()  ! harmonic plane
     integer,pointer :: ndc2t   => null()  ! 0 or 2 : jordan planes     
     integer,pointer :: pos_of_delta   => null()  !  constant delta
     integer,pointer :: rf   => null()  !   # of modulated planes

     logical(lp),pointer :: knob => null()
     logical(lp),pointer :: valishev => null()
     !     integer, pointer :: NDPT_OTHER
     logical(lp), pointer :: setknob => null()
     REAL(dp),pointer     :: da_absolute_aperture => null()  ! in case one tracks with da.
     !

     integer,pointer :: wherelost => null()     ! counting lost particles in integration nodes
     logical(lp),pointer  :: ROOT_CHECK => null()   !=.TRUE. performs check in roots and hyperbolic if true
     logical(lp),pointer  :: CHECK_STABLE => null() !=.TRUE. particle status
     logical(lp),pointer  :: CHECK_MADX_APERTURE => null()  !=.TRUE. false means particle lost in aperture
     logical(lp),pointer  :: APERTURE_FLAG => null()       !=.TRUE. aperture checks globally done (default)
     logical(lp),pointer  :: s_aperture_CHECK => null()       !=.TRUE. aperture checks globally done (default)


     logical(lp),pointer  :: WATCH_USER => null()     ! FALSE NORMALLY : WATCHES USER FOR FAILING TO CHECK APERTURES

     REAL(dp),pointer     :: absolute_aperture => null()     !=1e3_dp generic aperture check
     real(dp),pointer :: hyperbolic_aperture => null()  ! controls crashes in exponentials


     ! influence fibre creation

     integer, pointer :: MADTHICK => null()        !
     integer, pointer :: MADTHIN_NORMAL => null()
     integer, pointer :: MADTHIN_SKEW => null()
     integer, pointer::  NSTD => null(),METD => null()       ! number of steps and integration method
     logical(lp), pointer :: MADLENGTH => null()   ! =.false. rbend crazy length in mad8 as input
     logical(lp), pointer :: MAD       => null()   !=.false. mad definition of multipole for input only
     logical(lp), pointer :: EXACT_MODEL => null() != .false. exact model used
     logical(lp), pointer :: ALWAYS_EXACTMIS => null()  !=.TRUE. exact formula in tracking used for that element
     logical(lp),pointer :: ALWAYS_knobs => null()  !=.false. ptc knob default status
     logical(lp),pointer :: recirculator_cheat => null()  ! =.false.  if true energy patches use the time formula always
     logical(lp),pointer :: sixtrack_compatible => null() !  to insure some sixtrack compatibility default=false
     integer, pointer:: CAVITY_TOTALPATH => null() ! REAL PILL B0X =1 , FAKE =0  default
     integer,pointer :: HIGHEST_FRINGE => null() !=2  quadrupole fringe ON IF FRINGE PRESENT
     logical(lp),pointer :: do_beam_beam => null()   ! obvious meaning: false normally
     ! creates a reverse propagator
     integer,pointer ::FIBRE_DIR => null()         !=1 or -1 for reversed
     real(dp),pointer ::INITIAL_CHARGE => null()         ! =1 or -1 AND  ADJUST THE MASS IS THE PREFERED MODE
     ! creates a reverse propagator and a reversed ring in combination with above
     logical(lp),pointer ::FIBRE_flip => null()    !=.true.
     !  x_prime true means noncanonical outside magnets. x(5) variables stays the same.
     real(dp), pointer :: eps_pos => null()
     ! fill once and never touch again

     integer, pointer :: SECTOR_NMUL_MAX => null()     != 10 maxwell equations is solved to order 10 in exact sectors
     integer, pointer :: SECTOR_NMUL     => null() != 4  MULTIPOLES IN TEAPOT BEND ALLOWED BY DEFAULT
     real(dp), pointer :: wedge_coeff(:) => null()     ! QUAD_KICK IN WEDGE
     logical(lp), pointer :: MAD8_WEDGE  => null()     ! QUAD_KICK + FRINGE IF FRINGE IS OUT.

     logical(lp), pointer:: electron     => null()!  electron if true otherwise proton
     real(dp), pointer :: massfactor     => null()!=one  sets variable muon and electron must be true
     ! global on the fly
     logical(lp), pointer :: compute_stoch_kick => null() != .false. store stochastic kick for stochastic tracking
     logical(lp),pointer :: FEED_P0C => null()   !=.FALSE.  work takes p0c instead of energy
     logical(lp),pointer :: ALWAYS_EXACT_PATCHING => null()  !=.TRUE. patching done correctly
     ! used to output horror messages
     logical(lp),pointer :: stable_da => null() !=.true.  interrupts DA if check_da is true
     logical(lp),pointer :: check_da  => null() !=.true.
     logical(lp),pointer :: OLD_IMPLEMENTATION_OF_SIXTRACK => null()  !=.true.
     real(dp),pointer :: phase0  => null()! default phase in cavity
     logical(lp), pointer :: global_verbose => null()
     logical(lp), pointer :: no_hyperbolic_in_normal_form => null()! unstable produces exception
     integer, pointer :: ndpt_bmad => null() 
  end type CONTROL

  type(control) c_
 

 

contains

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


 real(dp) function bran(xran)
    implicit none
    !     ************************************
    !
    !     VERY SIMPLE RANDOM NUMBER GENERATOR
    !
    !-----------------------------------------------------------------------------
    !
    real(dp) xran
    !
    xran = xran + 10.0_dp
    if(xran.gt.1e4_dp) xran = xran - 9999.12345e0_dp
    bran = abs(sin(xran))
    bran = 10*bran
    bran = bran - int(bran)
    !      IF(BRAN.LT. c_0_1) BRAN = BRAN + c_0_1
    !
    return
  end function bran

  function mat_norm(m)
    implicit none
    real(dp) mat_norm
    real(dp) m(:,:)
    integer i,j

    mat_norm=0.0_dp
    do i=1,size(m,dim=1)
       do j=1,size(m,dim=2)
          mat_norm=mat_norm+abs(m(i,j))
       enddo
    enddo

  end function mat_norm



  ! Symplectic integrator routines setting coefficients

  SUBROUTINE MAKE_YOSHIDA
    IMPLICIT NONE
    integer i ,id 
    real(dp) a,b,b8,s21,ak

    YOSK(4)=0.0_dp
    YOSK(3)=0.78451361047756e0_dp
    YOSK(2)=0.235573213359357e0_dp
    YOSK(1)=-1.17767998417887e0_dp
    YOSK(0)=1.0_dp-2.0_dp*(YOSK(1)+YOSK(2)+YOSK(3))

    do i=4,1,-1
       YOSD(i)=(YOSK(i)+YOSK(i-1))/2.0_dp
    enddo

    do i=3,0,-1
       YOSK(i+1)=YOSK(I)
    enddo
!wyosh(1,1)= -0.161582374150097E1_dp   
!wyosh(1,2)= -0.244699182370524E1_dp   
!wyosh(1,3)= -0.716989419708120E-2_dp  
!wyosh(1,4)= 0.244002732616735E1_dp    
!wyosh(1,5)= 0.157739928123617E0_dp   
!wyosh(1,6)= 0.182020630970714E1_dp  
!wyosh(1,7)= 0.104242620869991E1_dp    

!wyosh(2,1)=  -0.169248587770116E-2_dp  
!wyosh(2,2)=  0.289195744315849E1_dp  
!wyosh(2,3)=  0.378039588360192E-2_dp   
!wyosh(2,4)= -0.289688250328827E1_dp   
!wyosh(2,5)=  0.289105148970595E1_dp   
!wyosh(2,6)=  -0.233864815101035E1_dp   
!wyosh(2,7)=  0.148819229202922E1_dp  

!wyosh(3,1)=  0.311790812418427E0_dp
!wyosh(3,2)=  -0.155946803821447E1_dp
!wyosh(3,3)=   -0.167896928259640E1_dp
!wyosh(3,4)=  0.166335809963315E1_dp
!wyosh(3,5)=  -0.106458714789183E1_dp
!wyosh(3,6)=  0.136934946416871E1_dp
!wyosh(3,7)=  0.629030650210433E0_dp

 ! real(dp) wyosh(0:7),wyoshid(0:15),wyoshik(15)    ! FIRST 8TH ORDER OF YOSHIDA

wyosh(1)= 0.102799849391985e0_dp 
wyosh(2)= -0.196061023297549E1_dp 
wyosh(3)= 0.193813913762276E1_dp  
wyosh(4)= -0.158240635368243E0_dp 
wyosh(5)= -0.144485223686048E1_dp 
wyosh(6)= 0.253693336566229E0_dp  
wyosh(7)= 0.914844246229740E0_dp  

!wyosh(5,1)=0.227738840094906E-1_dp
!wyosh(5,2)=0.252778927322839E1_dp
!wyosh(5,3)=-0.719180053552772E-1_dp
!wyosh(5,4)=0.536018921307285E-2_dp
!wyosh(5,5)=-0.204809795887393E1_dp
!wyosh(5,6)=0.107990467703699E0_dp
!wyosh(5,7)=0.130300165760014E1_dp

 
wyosh(0)=1.0_dp

do i=1,7
wyosh(0)=wyosh(0)-2*wyosh(i)
enddo

!wyosh(0:7),wyoshid(0:15),wyoshik(15) 
id=1
wyoshid(0)=wyosh(7)/2
do i=7,1,-1
 wyoshik(id)=wyosh(i)
 wyoshid(id)=(wyosh(i)+wyosh(i-1))/2
id=id+1
enddo
  wyoshik(id)=wyosh(0)
  wyoshid(id)=(wyosh(0)+wyosh(1))/2
  id=id+1

do i=1,6
 wyoshik(id)=wyosh(i)
 wyoshid(id)=(wyosh(i)+wyosh(i+1))/2
 id=id+1
enddo
 wyoshik(id)=wyosh(7)
 wyoshid(id)=wyosh(7)/2

if(.false.) then
write(6,*) id
!pause 777
 a=wyoshid(0) 
 b=0
 write(6,*) 0,wyoshid(0) 
do id=1,15
  a=a+wyoshid(id)
  b=b+wyoshik(id)
 write(6,*) id,wyoshid(id),wyoshik(id)
enddo
write(6,*) a,b
 stop 999
endif

ck8=0
bk8=0
ak8=0
!!! these are free parameters
ak8(10,5)=1.0_dp/9.0_dp
bk8(8)=49.0_dp/180.0_dp     ! cannot be zero

!!!  
s21=sqrt(21.0_dp)
b8=bk8(8)
ak=ak8(10,5)
bk8(7) = -b8+49.0_dp/180.0_dp


!  real(dp) ck8(11),bk8(11),ak8(11,11)
bk8(1)=1.0_dp/20.0_dp
bk8(9)=16.0_dp/45.0_dp
bk8(10)=49.0_dp/180.0_dp
bk8(11)=1.0_dp/20.0_dp

ck8(2)=1.0_dp/2.0_dp
ck8(3)=1.0_dp/2.0_dp
ck8(4)=(7.0_dp+s21)/14.0_dp
ck8(5)=(7.0_dp+s21)/14.0_dp
ck8(6)=1.0_dp/2.0_dp
ck8(7)=(7.0_dp-s21)/14.0_dp
ck8(8)=(7.0_dp-s21)/14.0_dp
ck8(9)=1.0_dp/2.0_dp
ck8(10)=(7.0_dp+s21)/14.0_dp
ck8(11)=1.0_dp

ak8(2,1)=0.5_dp
ak8(3,1)=0.25_dp
ak8(3,2)=0.25_dp
ak8(4,1)=1.0_dp/7.0_dp
ak8(4,2)= (-7.0_dp-3.0_dp*s21)/98.0_dp
ak8(4,3)= (21.0_dp+5.0_dp*s21)/49.0_dp
ak8(5,1)= (11.0_dp+s21)/84.0_dp
ak8(5,3)= (4.0_dp*s21)/63.0_dp + 2.0_dp/7.0_dp
ak8(5,4)= (21.0_dp-s21)/252.0_dp
ak8(6,1)= (5.0_dp+s21)/48.0_dp
ak8(6,3)= (9.0_dp+s21)/36.0_dp
ak8(6,4)= (-231.0_dp+14*s21)/360.0_dp
ak8(6,5)= (63.0_dp-7*s21)/80.0_dp
ak8(7,1)= (10.0_dp-s21)/42.0_dp
ak8(9,1)=  1.0_dp/32.0_dp
ak8(10,1)=  1.0_dp/14.0_dp
ak8(10,9)=  4.0_dp*s21/35.0_dp+132.0_dp/245.0_dp
ak8(11,9)=  (28.0_dp-28.0_dp*s21)/45.0_dp 
ak8(11,10)=  (49.0_dp-7.0_dp*s21)/18.0_dp 

ak8(7,3) = -(24.0_dp/35.0_dp)*ak-136.0_dp/105.0_dp+(-(12.0_dp/245.0_dp)*ak+ (656.0_dp/2205.0_dp))*s21


ak8(7,4) = 7.0_dp -(3.0_dp/10.0_dp)*ak*s21-(71.0_dp/45.0_dp)*s21 + (3.0_dp/10.0_dp)*ak


ak8(7,5) = -(3.0_dp/10.0_dp)*ak + (3.0_dp/10.0_dp)*ak*s21-43.0_dp/6.0_dp+ (169.0_dp/105.0_dp)*s21
ak8(7,6) = -(277.0_dp/735.0_dp)*s21+181.0_dp/105.0_dp+(12.0_dp/245.0_dp)*ak*s21+24.0_dp/35.0_dp*ak
ak8(8,1) = -(180.0_dp*b8*s21-49.0_dp*s21-1800*b8+343.0_dp)/b8/7560.0_dp
ak8(8,5) = -(441.0_dp*ak*s21-3240.0_dp*ak8(7,5)*b8-28.0_dp*s21+882.0_dp*ak8(7,5)-2205.0_dp*ak &
          +147.0_dp)/b8/3240.0_dp
ak8(8,6) =  (72.0_dp*ak*s21+1620.0_dp*ak8(7,6)*b8-29.0_dp*s21-441.0_dp*ak8(7,6)-252.0_dp*ak+119.0_dp )/1620.0_dp/b8

ak8(8,3) = -(900.0_dp*b8*s21+11340.0_dp*ak8(7,2)*b8+11340.0_dp*ak8(8,6)*b8-98.0_dp*s21-3087.0_dp*ak8(7,2) - 4860.0_dp*b8 &
+686.0_dp )/11340.0_dp/b8

ak8(8,7) = 49.0_dp/1620.0_dp/b8

ak8(8,4) = ( ck8(8)**2/2 -  ak8(8,2)*ck8(2) - ak8(8,3)*ck8(3) - ak8(8,5)*ck8(5)- ak8(8,6)*ck8(6)-ak8(8,7)*ck8(7))/ck8(4)


ak8(9,3) =   (1.0_dp/8.0_dp)*ak*s21-(1.0_dp/8.0_dp)*ak-1.0_dp/72.0_dp*s21 + 1.0_dp/72.0_dp

ak8(9,4) =   -49.0_dp/288.0_dp -(7.0_dp/32.0_dp)*ak*s21 +(7.0_dp/288.0_dp)*s21 + (49.0_dp/32.0_dp)*ak
! etienne
ak8(9,5) =   (7.0_dp/32.0_dp)*ak*s21-(35.0_dp/576.0_dp)*s21-(49.0_dp/32.0_dp)*ak+21.0_dp/64.0_dp


ak8(9,6) =   -(1.0_dp/8.0_dp)*ak*s21+(1.0_dp/8.0_dp)*ak+(1.0_dp/72.0_dp)*s21 + 5.0_dp/36.0_dp

ak8(9,7) =   (91.0_dp/576.0_dp)+(7.0_dp/192.0_dp)*s21-(585.0_dp/1568.0_dp)*b8*s21-(405.0_dp/224.0_dp)*b8

ak8(9,8) = (585.0_dp/1568.0_dp)*b8*s21 +(405.0_dp/224.0_dp)*b8

ak8(10,3) = -(6.0_dp/49.0_dp)*ak*s21 -(2.0_dp/7.0_dp)*ak + (2.0_dp/147.0_dp)*s21 +2.0_dp/63.0_dp

ak8(10,4) = 1.0_dp/9.0_dp - ak

ak8(10,6) = (2.0_dp/7.0_dp)*ak - (803.0_dp/2205.0_dp) + (6.0_dp/49.0_dp)*ak*s21 -(59.0_dp/735.0_dp)*s21

ak8(10,7) = 1.0_dp/9.0_dp + (1.0_dp/42.0_dp)*s21 +(2295.0_dp/686.0_dp)*b8 + (495.0_dp/686.0_dp)*b8*s21

ak8(10,8) = -(2295.0_dp/686.0_dp)*b8 - (495.0_dp/686.0_dp)*b8*s21

ak8(11,3) = (2.0_dp/3.0_dp)*ak*s21 -(2.0_dp/3.0_dp)*ak - (2.0_dp/27.0_dp)*s21 +2.0_dp/27.0_dp

ak8(11,4) = -(7.0_dp/6.0_dp)*ak*s21 +(7.0_dp/54.0_dp)*s21 + (49.0_dp/6.0_dp)*ak -49.0_dp/54.0_dp

ak8(11,5) = (7.0_dp/27.0_dp)*s21 -(77.0_dp/54.0_dp)  - (49.0_dp/6.0_dp)*ak + (7.0_dp/6.0_dp)*ak*s21

ak8(11,6) = (2.0_dp/3.0_dp) * ak - (64.0_dp/135.0_dp) -(2.0_dp/3.0_dp)*ak*s21 +(94.0_dp/135.0_dp)*s21

ak8(11,7) =  7.0_dp/18.0_dp -(265.0_dp/98.0_dp) *b8*s21 - (215.0_dp/14.0_dp)*b8

ak8(11,8) = (295.0_dp/98.0_dp)*b8*s21 + (215.0_dp/14.0_dp)*b8

  END SUBROUTINE MAKE_YOSHIDA

  SUBROUTINE input_sector(se2,se1)
    implicit none
    logical ttt,t1,t2
    integer se1,se2

    !    if(.not.change_sector) return

    t1=(SECTOR_NMUL_MAX/=se1)
    t2=(SECTOR_NMUL/=se2)

    ttt=t1.or.t2

    if(ttt) then
       if(change_sector) then
  !        write(6,*) " SECTOR_NMUL_MAX is changed from ",SECTOR_NMUL_MAX," to ",se1
          write(6,*) " SECTOR_NMUL is changed from ",SECTOR_NMUL," to ",se2
          write(6,*) " GLOBAL VARIABLES that can no longer be changed"
          SECTOR_NMUL_MAX=se1
          SECTOR_NMUL=se2
       else
 !         if(t1) write(6,*) " sector_nmul_max CANNOT be changed from ",SECTOR_NMUL_MAX," to ",se1
          if(t2) write(6,*) " sector_nmul CANNOT be changed from ",SECTOR_NMUL," to ",se2
          write(6,*) " Watch out : The are GLOBAL VARIABLES "
       endif
    endif

  end subroutine input_sector

  !  Printing routines for screen

  subroutine  get_ncar(n)
    implicit none
    integer n
    n=ncar
  end subroutine  get_ncar

 



 

!!!!!!!!!!!!!!!!!! old   special for lielib: others in sa_extend_poly
  REAL(DP) FUNCTION  ARCCOS_lielib(X)  ! REPLACES ACOS(X)
    IMPLICIT NONE
    REAL(DP),INTENT(IN)::X
    IF(.NOT.c_%CHECK_STABLE) then
       ARCCOS_lielib=0.0_dp
       return
    endif
    IF((ABS(X)>1.0_dp).AND.c_%ROOT_CHECK) THEN
       ARCCOS_lielib=0.0_dp
       c_%CHECK_STABLE=.FALSE.
       messagelost="a_scratch_size.f90 ARCCOS_lielib: abs(x)>1"
    ELSEIF(ABS(X)<=1.0_dp) THEN
       ARCCOS_lielib=ACOS(X)
    ELSE      !  IF X IS NOT A NUMBER
       ARCCOS_lielib=0.0_dp
       c_%CHECK_STABLE=.FALSE.
       messagelost="a_scratch_size.f90 ARCCOS_lielib: x is NaN"
    ENDIF

  END FUNCTION ARCCOS_lielib

  REAL(DP) FUNCTION  LOGE_lielib(X)  ! REPLACES ACOS(X)
    IMPLICIT NONE
    REAL(DP),INTENT(IN)::X
    IF(.NOT.c_%CHECK_STABLE) then
       LOGE_lielib=0.0_dp
       return
    endif

    IF(X<=0.0_dp.AND.c_%ROOT_CHECK) THEN
       LOGE_lielib=0.0_dp
       c_%CHECK_STABLE=.FALSE.
       messagelost="a_scratch_size.f90 LOGE_lielib: x<0"
    ELSE
       LOGE_lielib=LOG(X)
    ENDIF

  END FUNCTION LOGE_lielib

 subroutine dofma(nt,dt,xlist,pxlist,q0,qmean,qvar)
! demin's Laskar routine
! nt number of turns (periods) nt=2048
! time winows in units of turns   dt=512   max dt=nt/2
! xlist(1:nt), pxlist(1:nt) one-plane data usually x-px
! q0 closed orbit tune in x-px plane in revolutions
! output qmean average tune in nt window
! qvar = variance of measured tune
!
      implicit none
      
      ! Input parameters
      integer nt, dt
      real(dp) :: xlist(1:nt), pxlist(1:nt)
      real(dp) q0 
      
      ! Output
      real(dp) qmean, qvar
      
      ! Local
      real(dp) qmin, qmax, afind, bfind, stepfind, find, findnu, xnu
      real(dp) sumr, sumi, tw, rr, sumnorm
      real(dp), allocatable :: nulist(:)
      integer np, nwindow, iw, ik, ik1, i
      
      ! Preparations
      qmin = floor(q0*2.0)/2.0+0.001
      qmax = qmin + 0.498
      np = nt/2
      nwindow = (nt-np)/dt + 1
      allocate(nulist(1:nwindow))
      qmean = 0.0
      qvar = 0.0
      
      ! FMA calculation, Ref. D. Shatilov, Phys. Rev. ST Accel. Beams 14, 014001 (2011)
      do iw = 0, nwindow-1
        afind = qmax
        bfind = qmin
        stepfind = 0.002
        do while(stepfind>1.0e-12)
          find = 0.0
          findnu = 0.0
          xnu = afind
          do while(xnu .ge. bfind)
            sumr = 0.0
            sumi = 0.0
            ik = 0
            do while(ik .le. np)
              tw = 2.0 * dble(ik) / np - 1.0
              rr = twopi * ik * xnu
              ik1 = ik + iw * dt + 1
              sumr = sumr + (xlist(ik1) * cos(rr) + pxlist(ik1) * sin(rr)) * (1.0 + cos(pi*tw))
              sumi = sumi + (-xlist(ik1) * sin(rr) + pxlist(ik1) * cos(rr)) * (1.0 + cos(pi*tw))
              ik = ik + 1
            enddo
            sumnorm = sqrt((sumr*sumr+sumi*sumi)/np)
            if(find < sumnorm) then
              find = sumnorm
              findnu = xnu  
            end if
            xnu = xnu - stepfind
          enddo
          afind = findnu + stepfind
          bfind = findnu - stepfind
          stepfind = stepfind / 10.0
        enddo
        nulist(iw+1) = findnu
      enddo
      
      do i = 1, nwindow
        qmean = qmean + nulist(i)
      enddo
      qmean = qmean/nwindow
      
      do i = 1, nwindow
        qvar = qvar + (nulist(i)-qmean) * (nulist(i)-qmean)
      enddo
      qvar = sqrt(qvar/nwindow)

      return
      
   end subroutine dofma


  subroutine c_matinv(a,ai,n,nmx,ier)

    implicit none

    integer i,ier,j,n,nmx
    integer,dimension(nmax)::indx
    complex(dp) d
    complex(dp),dimension(nmx,nmx)::a,ai
    complex(dp),dimension(nmax,nmax)::aw
    !
    !    if((.not.C_%STABLE_DA)) then
    !       if(c_%watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif

    aw(1:n,1:n) = a(1:n,1:n)

    call ludcmp_nr(aw,n,nmax,indx,d,ier)
    if (ier .eq. 132) return

    ai(1:n,1:n) = 0.0_dp
    !    forall (i = 1:n) ai(i,i) = one
    do i=1,n
       ai(i,i) = 1.0_dp
    enddo

    do j=1,n
       call lubksb_nr(aw,n,nmax,indx,ai(1,j),nmx)
    enddo

  end subroutine c_matinv


  !
  subroutine ludcmp_nr(a,n,np,indx,d,ier)
    implicit none
    !     ************************************
    !
    !     THIS SUBROUTINE DECOMPOSES A MATRIX INTO LU FORMAT
    !     INPUT A: NXN MATRIX - WILL BE OVERWRITTEN BY THE LU DECOMP.
    !           NP: PHYSICAL DIMENSION OF A
    !           INDX: ROW PERMUTATION VECTOR
    !           D: EVEN OR ODD ROW INTERCHANGES
    !
    !     REFERENCE: NUMERICAL RECIPIES BY PRESS ET AL (CAMBRIDGE) PG. 35
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ier,imax,j,k,n,np
    integer,dimension(np)::indx
    complex(dp) d,dum,sum
    real(dp) aamax
    complex(dp),dimension(np,np)::a
    complex(dp),dimension(nmax)::vv
    !
    !    if((.not.C_%STABLE_DA)) then
    !       if(c_%watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    ier=0
    d=1.0_dp
    do i=1,n
       aamax=0.0_dp
       do j=1,n
          if(abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
       enddo
       if(aamax.eq.0.0_dp) then
          ier=132
          return
       endif
       vv(i)=1.0_dp/aamax
    enddo
    do j=1,n
       if(j.gt.1) then
          do i=1,j-1
             sum=a(i,j)
             if(i.gt.1) then
                do k=1,i-1
                   sum=sum-a(i,k)*a(k,j)
                enddo
                a(i,j)=sum
             endif
          enddo
       endif
       aamax=0.0_dp
       do i=j,n
          sum=a(i,j)
          if (j.gt.1) then
             do k=1,j-1
                sum=sum-a(i,k)*a(k,j)
             enddo
             a(i,j)=sum
          endif
          dum=vv(i)*abs(sum)
          if(abs(dum).ge.aamax) then
             imax=i
             aamax=dum
          endif
       enddo
       if (j.ne.imax) then
          do k=1,n
             dum=a(imax,k)
             a(imax,k)=a(j,k)
             a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
       endif
       indx(j)=imax
       if(j.ne.n) then
          if(a(j,j).eq.0.0_dp) a(j,j)=tiny
          dum=1.0_dp/a(j,j)
          do i=j+1,n
             a(i,j)=a(i,j)*dum
          enddo
       endif
    enddo
    if(a(n,n).eq.0.0_dp) a(n,n)=tiny
    return
  end subroutine ludcmp_nr
  !
  subroutine lubksb_nr(a,n,np,indx,b,nmx)
    implicit none
    !     ************************************
    !
    !     THIS SUBROUTINE SOLVES SET OF LINEAR EQUATIONS AX=B,
    !     INPUT A: NXN MATRIX IN lu FORM GIVEN BY ludcmp_nr
    !           NP: PHYSICAL DIMENSION OF A
    !           INDX: ROW PERMUTATION VECTOR
    !           D: EVEN OR ODD ROW INTERCHANGES
    !           B: RHS OF LINEAR EQUATION - WILL BE OVERWRITTEN BY X
    !
    !     REFERENCE: NUMERICAL RECIPIES BY PRESS ET AL (CAMBRIDGE) PG. 36
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ii,j,ll,n,nmx,np
    integer,dimension(np)::indx
    complex(dp) sum
    complex(dp),dimension(np,np)::a
    complex(dp),dimension(nmx)::b
    !
    !    if((.not.C_%STABLE_DA)) then
    !       if(c_%watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    ii = 0
    do i=1,n
       ll = indx(i)
       sum = b(ll)
       b(ll) = b(i)
       if(ii.ne.0) then
          do j=ii,i-1
             sum = sum-a(i,j)*b(j)
          enddo
       else if (abs(sum).ne.0.0_dp) then
          ii = i
       endif
       b(i)=sum
    enddo
    do i=n,1,-1
       sum=b(i)
       if(i.lt.n) then
          do j=i+1,n
             sum = sum-a(i,j)*b(j)
          enddo
       endif

       b(i)=sum/a(i,i)

    enddo
    return
  end subroutine lubksb_nr
  !

  subroutine matinv(a,ai,n,nmx,ier)

    implicit none

    integer i,ier,j,n,nmx
    integer,dimension(nmax)::indx
    real(dp) d
    real(dp),dimension(nmx,nmx)::a,ai
    real(dp),dimension(nmax,nmax)::aw
    !
    !    if((.not.C_%STABLE_DA)) then
    !       if(c_%watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif

    aw(1:n,1:n) = a(1:n,1:n)

    call ludcmp_nr0(aw,n,nmax,indx,d,ier)
    if (ier .eq. 132) return

    ai(1:n,1:n) = 0.0_dp
    !    forall (i = 1:n) ai(i,i) = one
    do i=1,n
       ai(i,i) = 1.0_dp
    enddo

    do j=1,n
       call lubksb_nr0(aw,n,nmax,indx,ai(1,j),nmx)
    enddo

  end subroutine matinv


  !
  subroutine ludcmp_nr0(a,n,np,indx,d,ier)
    implicit none
    !     ************************************
    !
    !     THIS SUBROUTINE DECOMPOSES A MATRIX INTO LU FORMAT
    !     INPUT A: NXN MATRIX - WILL BE OVERWRITTEN BY THE LU DECOMP.
    !           NP: PHYSICAL DIMENSION OF A
    !           INDX: ROW PERMUTATION VECTOR
    !           D: EVEN OR ODD ROW INTERCHANGES
    !
    !     REFERENCE: NUMERICAL RECIPIES BY PRESS ET AL (CAMBRIDGE) PG. 35
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ier,imax,j,k,n,np
    integer,dimension(np)::indx
    real(dp) aamax,d,dum,sum
    real(dp),dimension(np,np)::a
    real(dp),dimension(nmax)::vv
    !
    !    if((.not.C_%STABLE_DA)) then
    !       if(c_%watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    ier=0
    d=1.0_dp
    do i=1,n
       aamax=0.0_dp
       do j=1,n
          if(abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
       enddo
       if(aamax.eq.0.0_dp) then
          ier=132
          return
       endif
       vv(i)=1.0_dp/aamax
    enddo
    do j=1,n
       if(j.gt.1) then
          do i=1,j-1
             sum=a(i,j)
             if(i.gt.1) then
                do k=1,i-1
                   sum=sum-a(i,k)*a(k,j)
                enddo
                a(i,j)=sum
             endif
          enddo
       endif
       aamax=0.0_dp
       do i=j,n
          sum=a(i,j)
          if (j.gt.1) then
             do k=1,j-1
                sum=sum-a(i,k)*a(k,j)
             enddo
             a(i,j)=sum
          endif
          dum=vv(i)*abs(sum)
          if(dum.ge.aamax) then
             imax=i
             aamax=dum
          endif
       enddo
       if (j.ne.imax) then
          do k=1,n
             dum=a(imax,k)
             a(imax,k)=a(j,k)
             a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
       endif
       indx(j)=imax
       if(j.ne.n) then
          if(a(j,j).eq.0.0_dp) a(j,j)=tiny
          dum=1.0_dp/a(j,j)
          do i=j+1,n
             a(i,j)=a(i,j)*dum
          enddo
       endif
    enddo
    if(a(n,n).eq.0.0_dp) a(n,n)=tiny
    return
  end subroutine ludcmp_nr0
  !
  subroutine lubksb_nr0(a,n,np,indx,b,nmx)
    implicit none
    !     ************************************
    !
    !     THIS SUBROUTINE SOLVES SET OF LINEAR EQUATIONS AX=B,
    !     INPUT A: NXN MATRIX IN lu FORM GIVEN BY ludcmp_nr
    !           NP: PHYSICAL DIMENSION OF A
    !           INDX: ROW PERMUTATION VECTOR
    !           D: EVEN OR ODD ROW INTERCHANGES
    !           B: RHS OF LINEAR EQUATION - WILL BE OVERWRITTEN BY X
    !
    !     REFERENCE: NUMERICAL RECIPIES BY PRESS ET AL (CAMBRIDGE) PG. 36
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ii,j,ll,n,nmx,np
    integer,dimension(np)::indx
    real(dp) sum
    real(dp),dimension(np,np)::a
    real(dp),dimension(nmx)::b
    !
    !    if((.not.C_%STABLE_DA)) then
    !       if(c_%watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    ii = 0
    do i=1,n
       ll = indx(i)
       sum = b(ll)
       b(ll) = b(i)
       if(ii.ne.0) then
          do j=ii,i-1
             sum = sum-a(i,j)*b(j)
          enddo
       else if (sum.ne.0.0_dp) then
          ii = i
       endif
       b(i)=sum
    enddo
    do i=n,1,-1
       sum=b(i)
       if(i.lt.n) then
          do j=i+1,n
             sum = sum-a(i,j)*b(j)
          enddo
       endif

       b(i)=sum/a(i,i)

    enddo
    return
  end subroutine lubksb_nr0
  !
end module precision_constants


module scratch_size
  implicit none
  public
  integer,parameter::ndumt=10 
  integer,parameter::c_ndumt=10                   !  Number of scratch level
end module scratch_size

module file_handler
  use precision_constants
  implicit none
  public
  private INTFILE,INTFILE_K
  integer,private,parameter::nf=3,MFI=20,MFO=99
  integer,private :: winterfile(nf) =(/40,41,42/)
  logical(lp) :: myfile(MFI:MFO)
  logical(lp),private :: doit=.true.,ginofile=.false.

  type file_
     logical(lp) MF
     ! AIMIN CHANGES FOR MS4.0
     !   logical(lp) :: mf=.false.
  end type file_

  type file_K
     logical(lp) MF
     ! AIMIN CHANGES FOR MS4.0
     !  logical(lp) :: mf=.false.
  end type file_K

  type(FILE_)  NEWFILE
  type(FILE_K)  closeFILE


  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE INTFILE
     MODULE PROCEDURE INTFILE_K
  END  INTERFACE

CONTAINS


  SUBROUTINE  INTFILE(S2,S1)
    ! gino module here
    implicit none
    integer,INTENT(inout)::S2
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
          !       ! call !write_e(MFO)
       ENDIF
       S2=I
       myfile(I)=.true.
       if(s1%mf) then
          !      w_p=0
          !      w_p%nc=1
          !      w_p%fc='(1(1X,A120))'
          !      write(w_p%c(1),'(1x,L1,1x)')   s1%mf
          !      ! call ! WRITE_I
       endif
    endif
  END SUBROUTINE INTFILE


  SUBROUTINE  INTFILE_K(S2,S1)
    implicit none
    integer,INTENT(inOUT)::S2
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
    integer i

    DO I=MFI,MFO
       myfile(i)=.FALSE.
    ENDDO
    do i=1,nf
       myfile(winterfile(i))=.true.
    enddo
  END SUBROUTINE ZEROFILE

  SUBROUTINE KanalNummer(iff,file,old)
    implicit none
    integer, INTENT(OUT) :: iff
    character(*),optional :: file
    logical,optional :: old
    logical :: opened, exists
    integer :: i,ier

    DO i= 9999, 7, -1
       INQUIRE(UNIT= i, EXIST= exists, OPENED= opened)
       IF (exists .AND. (.NOT. opened)) GOTO 20
    ENDDO
    WRITE (UNIT= *, FMT= *) ' cannot find free unit within the range 7-9999..'
    CALL ReportOpenFiles
    STOP
20  CONTINUE
    iff= I

    if(present(file)) then
       if(present(old)) then
        if(old) then
        open(unit=iff,file=file,status='OLD',iostat=ier)
         if(ier/=0) then
          write(6,*) " file ",file," does not exist "
          stop 864
          endif
        else
         open(unit=iff,file=file,status='NEW',iostat=ier)
         if(ier/=0) then
          write(6,*) " file ",file,"  exists already "
          stop 865
         endif
        endif
       else
        open(unit=iff,file=file)
       endif
    endif
  END SUBROUTINE KanalNummer

  SUBROUTINE ReportOpenFiles
    implicit none

    logical :: opened, exists, named
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


  SUBROUTINE CONTEXT( STRING, nb,dollar,maj )
    IMPLICIT NONE
    CHARACTER(*) STRING
    CHARACTER(1) C1
    integer, optional :: nb
    logical(lp), optional :: dollar ,maj
    integer I,J,K,nb0,count
     logical(lp) dol,ma
    nb0=0
    dol=.false.
    ma=.true.
    if(present(nb)) nb0=1
    if(present(dollar)) dol=dollar
    if(present(maj)) ma=maj
    J = 0
    count=0
    DO I = 1, LEN (STRING)
       C1 = STRING(I:I)
       if(dol) then
        if(c1=='$') c1="_"
       endif
       STRING(I:I) = ' '
       IF( C1 .NE. ' ' ) THEN
          if(count/=0.and.nb0==1) then
             J = J + 1
             STRING(J:J) = ' '
             count=0
          endif
          J = J + 1
          K = ICHAR( C1 )
          IF( K .GE. ICHAR('a') .AND. K .LE. ICHAR('z').and.ma ) THEN
             C1 = CHAR( K - ICHAR('a') + ICHAR('A') )
          ENDIF
          STRING(J:J) = C1
       else
          count=count+1
       ENDIF
    ENDDO
    string=string(1:len_trim(string))
    RETURN
  END  SUBROUTINE CONTEXT

  SUBROUTINE create_name(name_root,ind,suffix,filename)
    implicit none
    CHARACTER(*) name_root,suffix,filename
    integer :: ind
    character(20) temp
    call context(name_root)
    call context(suffix)

    write(temp,*) ind
    call context(temp)

    filename=' '

    filename=name_root(1:len_trim(name_root))//temp(1:len_trim(temp))//'.'//suffix(1:len_trim(suffix))

  END SUBROUTINE create_name

end module file_handler


module my_own_1D_TPSA
  USE precision_constants
  implicit none
  !  public
  private input_real_in_my_1D_taylor
  integer :: n_tpsa_exp = 10
  integer, parameter :: N_my_1D_taylor=31  ! SHOULD BE AS ENGE_N
  integer, private :: No_my_1D_taylor=N_my_1D_taylor  ! SHOULD BE AS ENGE_N

  private mul,dmulsc,dscmul,INV
  private div,ddivsc,dscdiv,Idivsc
  private add,unaryADD,daddsc,dscadd
  private subs,unarySUB,dsubsc,dscsub,POW
  private  DEXPT,DLOGT,DSQRTT
  private  DCOST,DSINT

  ! Definitions
  type my_1D_taylor
     real(dp) a(0:N_my_1D_taylor)

  END type my_1D_taylor

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

  subroutine set_my_taylor_no(no)
    implicit none
    integer no

    if(no>N_my_1D_taylor) then
       No_my_1D_taylor=N_my_1D_taylor
       write(6,*) " warning NO too big in set_my_taylor_no: recompile FPP if needed "
    else
       No_my_1D_taylor=no
    endif

  end subroutine set_my_taylor_no

  FUNCTION add( S1, S2 )
    implicit none
    type (my_1D_taylor) add
    type (my_1D_taylor), INTENT (IN) :: S1, S2

    add%a=S1%a + S2%a

  END FUNCTION add

  FUNCTION daddsc( S1, sc )
    implicit none
    type (my_1D_taylor) daddsc
    type (my_1D_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    daddsc=s1
    daddsc%a(0)= s1%a(0) + sc

  END FUNCTION daddsc

  FUNCTION dscadd( sc ,  S1)
    implicit none
    type (my_1D_taylor) dscadd
    type (my_1D_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    dscadd=s1
    dscadd%a(0)= s1%a(0) + sc

  END FUNCTION dscadd

  FUNCTION unaryadd( S1 )
    implicit none
    type (my_1D_taylor) unaryadd
    type (my_1D_taylor), INTENT (IN) :: S1

    unaryadd=s1

  END FUNCTION unaryadd


  FUNCTION subs( S1, S2 )
    implicit none
    type (my_1D_taylor) subs
    type (my_1D_taylor), INTENT (IN) :: S1, S2

    subs%a=S1%a - S2%a

  END FUNCTION subs

  FUNCTION unarySUB( S1 )
    implicit none
    type (my_1D_taylor) unarySUB
    type (my_1D_taylor), INTENT (IN) :: S1

    unarySUB%a=-s1%a

  END FUNCTION unarySUB

  FUNCTION dsubsc( S1, sc )
    implicit none
    type (my_1D_taylor) dsubsc
    type (my_1D_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    dsubsc=s1
    dsubsc%a(0)= s1%a(0) - sc

  END FUNCTION dsubsc


  FUNCTION dscsub( sc , S1 )
    implicit none
    type (my_1D_taylor) dscsub
    type (my_1D_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    dscsub=-s1   ! uses unary sub
    dscsub=sc + dscsub   ! uses add

  END FUNCTION dscsub


  FUNCTION mul( S1, S2 )
    implicit none
    type (my_1D_taylor) mul
    type (my_1D_taylor), INTENT (IN) :: S1, S2
    integer I,J
    mul%a=0.0_dp
    DO I=0,No_my_1D_taylor
       DO J=0,No_my_1D_taylor
          IF(I+J>No_my_1D_taylor) CYCLE
          mul%a(I+J)=S1%a(I)*S2%a(J)+ mul%a(I+J)
       ENDDO
    ENDDO

  END FUNCTION mul

  FUNCTION dmulsc( S1, sc )
    implicit none
    type (my_1D_taylor) dmulsc
    type (my_1D_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    dmulsc%a= s1%a*sc

  END FUNCTION dmulsc

  FUNCTION dscmul( sc ,S1 )
    implicit none
    type (my_1D_taylor) dscmul
    type (my_1D_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    dscmul%a= s1%a*sc

  END FUNCTION dscmul

  FUNCTION ddivsc( S1, sc )
    implicit none
    type (my_1D_taylor) ddivsc
    type (my_1D_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    ddivsc%a= s1%a/sc

  END FUNCTION ddivsc

  FUNCTION Idivsc( S1, sc )
    implicit none
    type (my_1D_taylor) Idivsc
    type (my_1D_taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc

    Idivsc%a= s1%a/sc

  END FUNCTION Idivsc

  FUNCTION POW( S1,N )
    implicit none
    type (my_1D_taylor) POW , T
    type (my_1D_taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: N
    integer I

    POW=1.0_dp

    IF(N<=0) THEN
       T=1.0_dp/S1
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
    type (my_1D_taylor) inv,t,TT
    type (my_1D_taylor), INTENT (IN) :: S1
    integer I
    T=S1/S1%A(0)
    T%A(0)=0.0_dp

    TT=T
    inv=1.0_dp

    DO I=1,No_my_1D_taylor
       INV=INV-TT
       TT=-TT*T
    ENDDO

    INV=INV/S1%A(0)

  END FUNCTION inv

  FUNCTION DIV( S1, S2 )
    implicit none
    type (my_1D_taylor) DIV
    type (my_1D_taylor), INTENT (IN) :: S1, S2

    DIV=INV(S2)
    DIV=S1*DIV
  END FUNCTION DIV

  FUNCTION dscdiv(  sc , S1)
    implicit none
    type (my_1D_taylor) dscdiv
    type (my_1D_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    dscdiv=INV(S1)
    dscdiv= SC *  dscdiv

  END FUNCTION dscdiv

  !  DEFININING my_1D_taylor=CONSTANT
  subroutine input_real_in_my_1D_taylor( S2, S1 )
    implicit none
    real(dp), INTENT (IN) :: S1
    type (my_1D_taylor), INTENT (inout) :: S2

    S2%a=0.0_dp
    S2%a(0)=s1

  END subroutine input_real_in_my_1D_taylor

  subroutine input_my_1D_taylor_in_real( S2, S1 )
    implicit none
    type (my_1D_taylor), INTENT (IN) :: S1
    real(dp), INTENT (inout) :: S2

    S2=S1%a(0)

  END subroutine input_my_1D_taylor_in_real

  !  DEFININING EXP

  FUNCTION DEXPT( S1 )
    implicit none
    type (my_1D_taylor) DEXPT,t,tt
    type (my_1D_taylor), INTENT (IN) :: S1
    integer I
    T=S1
    T%A(0)=0.0_dp


    DEXPT=1.0_dp
    tt=1.0_dp

    do i=1,No_my_1D_taylor
       tt=tt*t/i
       DEXPT=DEXPT + tt
    enddo

    DEXPT=DEXPT*EXP(S1%A(0))

  END FUNCTION DEXPT


  FUNCTION DLOGT( S1 )
    implicit none
    type (my_1D_taylor) DLOGT,t,TT
    type (my_1D_taylor), INTENT (IN) :: S1
    integer I

    T=S1/S1%A(0)
    T%A(0)=0.0_dp
    DLOGT=0.0_dp
    TT=T
    do i=1,No_my_1D_taylor
       DLOGT=TT/I+DLOGT
       TT=-TT*T
    ENDDO

    DLOGT=DLOGT+ LOG(S1%A(0))

  END FUNCTION  DLOGT

  FUNCTION DSQRTT( S1 )
    implicit none
    type (my_1D_taylor) DSQRTT 
    type (my_1D_taylor), INTENT (IN) :: S1
    ! WRITE(6,*) " MARDE "
    ! STOP 666
    ! T=S1/S1%A(0)
    ! T%A(0)=zero

    ! DSQRTT=one + T/two - T**2/eight+T**3/c_16
    ! DSQRTT=DSQRTT* SQRT(S1%A(0))
    DSQRTT=log(s1)/2.0_dp
    DSQRTT=exp(DSQRTT)
  END FUNCTION  DSQRTT

  FUNCTION DCOST(S1)
    implicit none
    type (my_1D_taylor) DCOST,t
    type (my_1D_taylor), INTENT (IN) :: S1
    WRITE(6,*) " MARDE "
    STOP 666

    T=S1
    T%A(0)=0.0_dp

    DCOST=COS(S1%A(0))*(1.0_dp-T**2/2.0_dp)-SIN(S1%A(0))*(T-T**3/6.0_dp)

  END FUNCTION  DCOST

  FUNCTION DSINT(S1)
    implicit none
    type (my_1D_taylor) DSINT,t
    type (my_1D_taylor), INTENT (IN) :: S1

    T=S1
    T%A(0)=0.0_dp

    DSINT=SIN(S1%A(0))*(1.0_dp-T**2/2.0_dp)+COS(S1%A(0))*(T-T**3/6.0_dp)

  END FUNCTION  DSINT

end module my_own_1D_TPSA

module gauss_dis

  use precision_constants
  public
  private alex_ISEED
  integer :: alex_ISEED=1000
  !cccccccccccccccccccccccc     Program 2.13     ccccccccccccccccccccccccc
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !                                                                      c
  ! Please Note:                                                         c
  !                                                                      c
  ! (1) This computer program is part of the book, "An Introduction to   c
  !     Computational Physics," written by Tao Pang and published and    c
  !     copyrighted by Cambridge University Press in 1997.               c
  !                                                                      c
  ! (2) No warranties, express or implied, are made for this program.    c
  !                                                                      c
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
contains

  SUBROUTINE gaussian_seed(i)
    implicit none
    integer i
    alex_ISEED=i
  end SUBROUTINE gaussian_seed

  SUBROUTINE GRNF(X,cut)
    implicit none
    real(dp) r1,r2,x,cut


1   R1 = -LOG(1.0_dp-RANF())
    R2 = 2.0_dp*PI*RANF()
    R1 = SQRT(2.0_dp*R1)
    X  = R1*COS(R2)
    if(abs(x)>cut) goto 1
    ! Y  = R1*SIN(R2)
    RETURN
  END SUBROUTINE GRNF
  !
  real(dp) FUNCTION RANF()
    implicit none
    integer ia,ic,iq,ir,ih,il,it
    DATA IA/16807/,IC/2147483647/,IQ/127773/,IR/2836/
    IH = alex_ISEED/IQ
    IL = MOD(alex_ISEED,IQ)
    IT = IA*IL-IR*IH
    IF(IT.GT.0) THEN
       alex_ISEED = IT
    ELSE
       alex_ISEED = IC+IT
    END IF
    RANF = alex_ISEED/FLOAT(IC)
    RETURN
  END FUNCTION RANF

end module gauss_dis


module my_own_linear_tpsa   
use precision_constants  
use file_handler 
implicit none
!@3 <font color="#FF0000" size="5"><p>Private routines only accessible through an interface</font></p>
 private input_my_linear_taylor_in_real
 private input_real_in_my_linear_taylor
 private input_comp_in_my_linear_taylor
 private input_my_linear_taylor_in_comp
 private input_my_linear_taylor_in_my_linear_taylor
 private mul   
 private dmulsc 
 private dscmul 
 private dmulscc
 private dscmulc
 private div     
 private ddivsc  
 private ddivscc 
 private dscdiv  
 private Idivsc  
 private add      
 private unaryADD 
 private daddsc   
 private dscadd   
 private daddscc  
 private dscaddc  
 private subs     
 private unarySUB 
 private dsubsc   
 private dscsub   
 private dsubscc  
 private dscsubc  
 private POW,powr
 private DEXPT  
 private DLOGT  
 private DSQRTT 
 private DCOST
 private DSINT
 private atan2t,dtant,atant
 private tpsa_expt
 private MONOT
 private alloc_my_linear_taylor  
 private morpht
 private clean_taylor
 !@  <hr>
  complex   (dp), parameter :: I_=(0,1) !@1 &nbsp; complex number <font face="Times New Roman">&#8730;-1</font>
  real(dp) ::  epsclean=1.d-10   !@1 &nbsp; number used to insure that small numbers are set to zero
 integer :: n_tpsa_exp = 10      !@1 &nbsp; used in a "bad" way to do the exponential of a Taylor series
  real(dp) ::  epsprint=1.d-10

integer,parameter :: n_mono =6


! Definitions                
 TYPE my_linear_taylor             
    complex(dp) a(0:n_mono)      !@2  &nbsp; Taylor series テイラー展開（テイラーてんかい) 浮動小数点数（ふどうしょうすうてんすう）
END TYPE my_linear_taylor
 
type(my_linear_taylor) dxl(n_mono)
!type(my_linear_taylor) dx_1,dx_2,dx_3


 

  INTERFACE assignment (=)
     MODULE PROCEDURE input_my_linear_taylor_in_my_linear_taylor   !@1 &nbsp; Taylor=Taylor
     MODULE PROCEDURE input_my_linear_taylor_in_real   !@1 &nbsp; real=Taylor
     MODULE PROCEDURE input_real_in_my_linear_taylor  !@1 &nbsp; Taylor=real 
     MODULE PROCEDURE input_real_in_my_linear_taylors  !@1 &nbsp; Taylor=real 
     MODULE PROCEDURE input_comp_in_my_linear_taylor  !@1 &nbsp; complex=Taylor
     MODULE PROCEDURE input_my_linear_taylor_in_comp   !@1 &nbsp; complex=Taylor
  end  INTERFACE     
!@ <font color="#800000"><b>N.B. </b></font>In this toy DA, Taylor=Taylor is not 
!@ needed because Fortran90 does automatically the correct job. This is <b>
!@ <font color="#FF0000">not</font></b> true when we overload a pointer as in the 
!@ case of Berz's package.

  INTERFACE OPERATOR (*)  !@2  &nbsp; Operator overloading (or ad hoc polymorphism)  オペレータ(演算子)のオーバーロード
!@3  <b><font color="#FF0000">Taylor=Taylor*Taylor is interesting </font></b> click below </br>
     MODULE PROCEDURE mul   !@1 &nbsp; Taylor * Taylor   
     MODULE PROCEDURE dmulsc    !@1 &nbsp; Taylor * Real(dp)
     MODULE PROCEDURE dscmul     !@1 &nbsp;  Real(dp) * Taylor
     MODULE PROCEDURE dmulscc    !@1 &nbsp; Taylor * complex(dp)
     MODULE PROCEDURE dscmulc     !@1 &nbsp;  complex(dp) * Taylor
  END INTERFACE
  
  INTERFACE OPERATOR (/)  !@2  &nbsp; Operator overloading (or ad hoc polymorphism)  オペレータ(演算子)のオーバーロード
!@3  <b><font color="#FF0000">Taylor=Taylor/Taylor is also interesting because first case separating TPSA from DA part to create a nilpotent operator!</font></b> click below </br>
     MODULE PROCEDURE div    !@1 &nbsp; Taylor / Taylor
     MODULE PROCEDURE ddivsc  !@1 &nbsp; Taylor / Real(dp)
     MODULE PROCEDURE ddivscc  !@1 &nbsp; Taylor / complex(dp)
     MODULE PROCEDURE dscdiv  !@1 &nbsp; Real(dp) / Taylor
     MODULE PROCEDURE Idivsc  !@1 &nbsp; Taylor / Integer  <font color="#FF0000">&nbsp; &#8594; &nbsp; added because useful in example code</font>
  END INTERFACE

  INTERFACE OPERATOR (+)  !@2  &nbsp; Operator overloading (or ad hoc polymorphism)  オペレータ(演算子)のオーバーロード
     MODULE PROCEDURE add       !@1 &nbsp; Taylor + Taylor
     MODULE PROCEDURE unaryADD  !@1 &nbsp;  +Taylor
     MODULE PROCEDURE daddsc    !@1 &nbsp;  Taylor + Real(dp)
     MODULE PROCEDURE dscadd    !@1 &nbsp;  Real(dp) + Taylor
     MODULE PROCEDURE daddscc    !@1 &nbsp;  Taylor + complex(dp)
     MODULE PROCEDURE dscaddc    !@1 &nbsp;  complex(dp) + Taylor
  END INTERFACE

  INTERFACE OPERATOR (-)  !@2  &nbsp; Operator overloading (or ad hoc polymorphism)  オペレータ(演算子)のオーバーロード
     MODULE PROCEDURE subs      !@1 &nbsp; Taylor - Taylor
     MODULE PROCEDURE unarySUB     !@1 &nbsp;  -Taylor
     MODULE PROCEDURE dsubsc   !@1 &nbsp;  Taylor - Real(dp)
     MODULE PROCEDURE dscsub    !@1 &nbsp;  Real(dp) - Taylor
     MODULE PROCEDURE dsubscc   !@1 &nbsp;  Taylor - complex(dp)
     MODULE PROCEDURE dscsubc    !@1 &nbsp;  complex(dp) - Taylor
  END INTERFACE


  INTERFACE OPERATOR (**)  !@2  &nbsp; Operator overloading (or ad hoc polymorphism)  オペレータ(演算子)のオーバーロード
     MODULE PROCEDURE POW    !@1 &nbsp; Taylor ** Integer 
     MODULE PROCEDURE POWr    !@1 &nbsp; Taylor ** real(dp) 
  END INTERFACE


!@3  <p><i><font size="5">&nbsp;Overloading standard procedures </font></i></p>

  INTERFACE exp
     MODULE PROCEDURE DEXPT   !@1 &nbsp; exp(Taylor)
  END INTERFACE

  INTERFACE tpsa_exp
     MODULE PROCEDURE tpsa_expt   !@1 &nbsp; exp(Taylor) using a sum up to n_tpsa_exp
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

  INTERFACE tan
     MODULE PROCEDURE dtant   !@1 &nbsp; tan(Taylor)
  END INTERFACE

  INTERFACE atan
     MODULE PROCEDURE atant   !@1 &nbsp; tan(Taylor)
  END INTERFACE

  INTERFACE atan2
     MODULE PROCEDURE atan2t   !@1 &nbsp; Atan2(Taylor,Taylor)
  END INTERFACE




!@3  <p><i><font size="5">User defined operator</font></i></p>

  INTERFACE OPERATOR (.mon.)
     MODULE PROCEDURE MONOT    !@1 <font color="#FF0000">Creates the monomial </font><font color="#0000FF">r x<sub>i</sub></font>
  END INTERFACE
 
! Destructors and Constructors for my_linear_taylor

  INTERFACE alloc
     MODULE PROCEDURE alloc_my_linear_taylor   !@1 &nbsp; for compatibility with FPP
  END INTERFACE

  INTERFACE kill
     MODULE PROCEDURE alloc_my_linear_taylor   !@1 &nbsp; for compatibility with FPP
  END INTERFACE

  INTERFACE morph
     MODULE PROCEDURE morpht   !@1 &nbsp; for compatibility with FPP
  END INTERFACE

  INTERFACE clean
     MODULE PROCEDURE clean_taylor   !@1 &nbsp;  removes numbers smaller than epsclean
  END INTERFACE

contains


  subroutine init_linear_taylor 
  implicit none
  integer i

   do i=1,n_mono
    dxl(i)%a=0.0_dp
    dxl(i)%a(i)=1
   enddo
 end subroutine init_linear_taylor

  subroutine clean_taylor( S2 )
    implicit none
    TYPE (my_linear_taylor), INTENT (INout) :: S2
    integer i
    
    do i=0,n_mono
     if(abs (real(S2%a(i),kind=dp) ) <epsclean) S2%a(i)= i_*aimag(S2%a(i))
     if(abs(aimag(S2%a(i)))  <epsclean)         S2%a(i)=real(S2%a(i),kind=dp)
    enddo     
  END subroutine clean_taylor

  subroutine clean_taylor_PRINT( S2 )
    implicit none
    TYPE (my_linear_taylor), INTENT (INout) :: S2
    integer i
    
    do i=0,n_mono
     if(abs (real(S2%a(i),kind=dp) ) <EPSPRINT) S2%a(i)= i_*aimag(S2%a(i))
     if(abs(aimag(S2%a(i)))  <EPSPRINT)         S2%a(i)=real(S2%a(i),kind=dp)
    enddo     
  END subroutine clean_taylor_PRINT


  subroutine alloc_my_linear_taylor( S2 )
    implicit none
    TYPE (my_linear_taylor), INTENT (INout) :: S2

     S2%a=0.0_dp
     
  END subroutine alloc_my_linear_taylor



  FUNCTION morpht( t )
    implicit none
    TYPE (my_linear_taylor) morpht
    TYPE (my_linear_taylor), INTENT(IN) :: t
    
	morpht= t
  
  END FUNCTION morpht

  FUNCTION MONOT( R, I )
    implicit none
    TYPE (my_linear_taylor) MONOT
    REAL(DP), INTENT(IN) :: R
    INTEGER, INTENT(IN) :: I
 
 

 
     MONOT=0.0_DP
     MONOT%a(i)= R

  END FUNCTION MONOT

 

 
 
 

  FUNCTION add( S1, S2 )
    implicit none
    TYPE (my_linear_taylor) add
    TYPE (my_linear_taylor), INTENT (IN) :: S1, S2

     add%a=S1%a + S2%a     

  END FUNCTION add

  FUNCTION daddsc( S1, sc )
    implicit none
    TYPE (my_linear_taylor) daddsc
    TYPE (my_linear_taylor), INTENT (IN) :: S1 
    real(dp), INTENT (IN) :: sc
     
     daddsc=s1
     daddsc%a(0)= s1%a(0) + sc    
     call clean(daddsc)
  END FUNCTION daddsc

  FUNCTION dscadd( sc ,  S1)
    implicit none
    TYPE (my_linear_taylor) dscadd
    TYPE (my_linear_taylor), INTENT (IN) :: S1 
    real(dp), INTENT (IN) :: sc
     
     dscadd=s1
     dscadd%a(0)= s1%a(0) + sc    
     call clean(dscadd)

  END FUNCTION dscadd

  FUNCTION daddscc( S1, sc )
    implicit none
    TYPE (my_linear_taylor) daddscc
    TYPE (my_linear_taylor), INTENT (IN) :: S1 
    complex(dp), INTENT (IN) :: sc
     
     daddscc=s1
     daddscc%a(0)= s1%a(0) + sc    
     call clean(daddscc)

  END FUNCTION daddscc

  FUNCTION dscaddc( sc ,  S1)
    implicit none
    TYPE (my_linear_taylor) dscaddc
    TYPE (my_linear_taylor), INTENT (IN) :: S1 
    complex(dp), INTENT (IN) :: sc
     
     dscaddc=s1
     dscaddc%a(0)= s1%a(0) + sc    
     call clean(dscaddc)
  END FUNCTION dscaddc

  FUNCTION unaryadd( S1 )
    implicit none
    TYPE (my_linear_taylor) unaryadd
    TYPE (my_linear_taylor), INTENT (IN) :: S1 
     
     unaryadd=s1

  END FUNCTION unaryadd


  FUNCTION subs( S1, S2 )
    implicit none
    TYPE (my_linear_taylor) subs
    TYPE (my_linear_taylor), INTENT (IN) :: S1, S2

     subs%a=S1%a - S2%a     
     call clean(subs)

  END FUNCTION subs

  FUNCTION unarySUB( S1 )
    implicit none
    TYPE (my_linear_taylor) unarySUB
    TYPE (my_linear_taylor), INTENT (IN) :: S1 
     
     unarySUB%a=-s1%a

  END FUNCTION unarySUB

  FUNCTION dsubsc( S1, sc )
    implicit none
    TYPE (my_linear_taylor) dsubsc
    TYPE (my_linear_taylor), INTENT (IN) :: S1 
    real(dp), INTENT (IN) :: sc
     
     dsubsc=s1
     dsubsc%a(0)= s1%a(0) - sc    
     call clean(dsubsc)

  END FUNCTION dsubsc


  FUNCTION dscsub( sc , S1 )
    implicit none
    TYPE (my_linear_taylor) dscsub
    TYPE (my_linear_taylor), INTENT (IN) :: S1 
    real(dp), INTENT (IN) :: sc
     
     dscsub=-s1   ! uses unary sub
     dscsub=sc + dscsub   ! uses add   
     call clean(dscsub)

  END FUNCTION dscsub

  FUNCTION dsubscc( S1, sc )
    implicit none
    TYPE (my_linear_taylor) dsubscc
    TYPE (my_linear_taylor), INTENT (IN) :: S1 
    complex(dp), INTENT (IN) :: sc
     
     dsubscc=s1
     dsubscc%a(0)= s1%a(0) - sc    
     call clean(dsubscc)
  END FUNCTION dsubscc


  FUNCTION dscsubc( sc , S1 )
    implicit none
    TYPE (my_linear_taylor) dscsubc
    TYPE (my_linear_taylor), INTENT (IN) :: S1 
    complex(dp), INTENT (IN) :: sc
     
     dscsubc=-s1   ! uses unary sub
     dscsubc=sc + dscsubc   ! uses add   
     call clean(dscsubc)

  END FUNCTION dscsubc


  FUNCTION mul( S1, S2 )
    implicit none
    TYPE (my_linear_taylor) mul
    TYPE (my_linear_taylor), INTENT (IN) :: S1, S2
    integer i
 
      mul%a=0.0_dp
     mul%a(0)=S1%a(0)*S2%a(0)
     do i=1,n_mono

 
      
      mul%a(i)=S1%a(i)*S2%a(0) +S2%a(i)*S1%a(0) 
     
     
     enddo
     

  END FUNCTION mul


  

  FUNCTION dmulsc( S1, sc )
    implicit none
    TYPE (my_linear_taylor) dmulsc
    TYPE (my_linear_taylor), INTENT (IN) :: S1 
    real(dp), INTENT (IN) :: sc

     dmulsc%a= s1%a*sc    
     call clean(dmulsc)

  END FUNCTION dmulsc

  FUNCTION dscmul( sc ,S1 )
    implicit none
    TYPE (my_linear_taylor) dscmul
    TYPE (my_linear_taylor), INTENT (IN) :: S1 
    real(dp), INTENT (IN) :: sc

     dscmul%a= s1%a*sc    


  END FUNCTION dscmul

  FUNCTION dmulscc( S1, sc )
    implicit none
    TYPE (my_linear_taylor) dmulscc
    TYPE (my_linear_taylor), INTENT (IN) :: S1 
    complex(dp), INTENT (IN) :: sc

     dmulscc%a= s1%a*sc    


  END FUNCTION dmulscc

  FUNCTION dscmulc( sc ,S1 )
    implicit none
    TYPE (my_linear_taylor) dscmulc
    TYPE (my_linear_taylor), INTENT (IN) :: S1 
    complex(dp), INTENT (IN) :: sc

     dscmulc%a= s1%a*sc    

  END FUNCTION dscmulc

  FUNCTION ddivsc( S1, sc )
    implicit none
    TYPE (my_linear_taylor) ddivsc
    TYPE (my_linear_taylor), INTENT (IN) :: S1 
    real(dp), INTENT (IN) :: sc
     
     ddivsc%a= s1%a/sc    


  END FUNCTION ddivsc

  FUNCTION ddivscc( S1, sc )
    implicit none
    TYPE (my_linear_taylor) ddivscc
    TYPE (my_linear_taylor), INTENT (IN) :: S1 
    complex(dp), INTENT (IN) :: sc
     
     ddivscc%a= s1%a/sc    

  END FUNCTION ddivscc



  FUNCTION Idivsc( S1, sc )
    implicit none
    TYPE (my_linear_taylor) Idivsc
    TYPE (my_linear_taylor), INTENT (IN) :: S1 
    INTEGER, INTENT (IN) :: sc
     
     Idivsc%a= s1%a/sc    

  END FUNCTION Idivsc

  FUNCTION POW( S1,N )
    implicit none
    TYPE (my_linear_taylor) POW , T
    TYPE (my_linear_taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: N
INTEGER I

    POW=1.0_dp

    IF(N<0) THEN
     T=1.0_dp/S1
     DO I=1,-N
 POW=POW*T
 ENDDO
ELSE
     DO I=1,N
 POW=POW*S1
 ENDDO

ENDIF

    END FUNCTION POW

  FUNCTION POWr( S1,N )
    implicit none
    TYPE (my_linear_taylor) POWr  
    TYPE (my_linear_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: N
 

       powr=log(s1)*n
       powr=exp(powr)


    END FUNCTION POWr

  FUNCTION inv( S1 )
    implicit none
    TYPE (my_linear_taylor) inv,t,tt
    TYPE (my_linear_taylor), INTENT (IN) :: S1
    integer i
    
T=S1/S1%A(0)
T%A(0)=0.D0
 INV=1.0_DP
 TT=1.0_DP
 TT=-T*TT
 INV=inv+TT

    INV=INV/S1%A(0)

  END FUNCTION inv

  FUNCTION DIV( S1, S2 )
    implicit none
    TYPE (my_linear_taylor) DIV
    TYPE (my_linear_taylor), INTENT (IN) :: S1, S2
!@3<b><font color="#FF0000">Go to function </font></b><a href="#INV">INV</a></br>
   DIV=INV(S2)
   DIV=S1*DIV

  END FUNCTION DIV

  FUNCTION dscdiv(  sc , S1)
    implicit none
    TYPE (my_linear_taylor) dscdiv
    TYPE (my_linear_taylor), INTENT (IN) :: S1 
    real(dp), INTENT (IN) :: sc
     
     dscdiv=INV(S1)
     dscdiv= SC *  dscdiv


  END FUNCTION dscdiv

!  DEFININING my_linear_taylor=CONSTANT   
  subroutine input_real_in_my_linear_taylor( S2, S1 )
    implicit none
    real(dp), INTENT (IN) :: S1
    TYPE (my_linear_taylor), INTENT (inout) :: S2

     S2%a=0.0_dp 
     S2%a(0)=s1

  END subroutine input_real_in_my_linear_taylor

!  DEFININING my_linear_taylor=CONSTANT   
  subroutine input_real_in_my_linear_taylors( S2, S1 )
    implicit none
    real(dp), INTENT (IN) :: S1
    TYPE (my_linear_taylor), INTENT (inout) :: S2(:)
    integer i
     do i=1,size(s2)
      S2(i)%a=0.0_dp 
      S2(i)%a(0)=s1
     enddo
  END subroutine input_real_in_my_linear_taylors

  subroutine input_comp_in_my_linear_taylor( S2, S1 )
    implicit none
    complex(dp), INTENT (IN) :: S1
    TYPE (my_linear_taylor), INTENT (inout) :: S2

     S2%a=0.0_dp 
     S2%a(0)=s1

  END subroutine input_comp_in_my_linear_taylor

  subroutine input_my_linear_taylor_in_comp( S2, S1 )
    implicit none
    TYPE (my_linear_taylor), INTENT (IN) :: S1
    complex(dp), INTENT (inout) :: S2

     S2=S1%a(0) 
     
  END subroutine input_my_linear_taylor_in_comp
  
    subroutine input_my_linear_taylor_in_my_linear_taylor( S2, S1 )
    implicit none
    TYPE (my_linear_taylor), INTENT (IN) :: S1
    TYPE (my_linear_taylor), INTENT (inout) :: S2

     S2%a=S1%a 
     
  END subroutine input_my_linear_taylor_in_my_linear_taylor

    subroutine input_my_linear_taylor_in_real( S2, S1 )
    implicit none
    TYPE (my_linear_taylor), INTENT (IN) :: S1
    real(dp), INTENT (inout) :: S2

     S2=S1%a(0) 
     
  END subroutine input_my_linear_taylor_in_real


!  DEFININING EXP

  FUNCTION DEXPT( S1 )
    implicit none
    TYPE (my_linear_taylor) DEXPT,t,tt
    TYPE (my_linear_taylor), INTENT (IN) :: S1

T=S1
T%A(0)=0.0_DP
    

DEXPT=1.0_dp
tt=1.0_dp


	 tt=tt*t
     DEXPT=DEXPT + tt


    DEXPT=EXP(S1%A(0))*DEXPT


  END FUNCTION DEXPT 


  FUNCTION DLOGT( S1 )
    implicit none
    TYPE (my_linear_taylor) DLOGT,T,TT
    TYPE (my_linear_taylor), INTENT (IN) :: S1

T=S1/S1%A(0); T%A(0)=0.0_DP;
      DLOGT=0.0_dp
     TT=-1.0_DP

      TT=-T*TT
      DLOGT=DLOGT+TT    

     
    DLOGT=DLOGT+ LOG(S1%A(0))
     call clean(DLOGT)

  END FUNCTION  DLOGT 

  FUNCTION DSQRTT( S1 )
    implicit none
    TYPE (my_linear_taylor) DSQRTT
    TYPE (my_linear_taylor), INTENT (IN) :: S1
!@3  <p><i><sup><font size="4" color="#0066FF"> </font></sup></i></p>

    DSQRTT= log(s1)/2.0_dp
    DSQRTT=exp(DSQRTT)



  END FUNCTION  DSQRTT 

  FUNCTION DCOST(S1)
    implicit none
    TYPE (my_linear_taylor) DCOST
    TYPE (my_linear_taylor), INTENT (IN) :: S1
!@3  <p><i><sup><font size="4" color="#0066FF"> </font></sup></i></p>

DCOST=(exp(i_*s1)+exp(-i_*s1))/2.0_dp


  END FUNCTION  DCOST 

  FUNCTION atant(S1)
    implicit none
    TYPE (my_linear_taylor) atant,c
    TYPE (my_linear_taylor), INTENT (IN) :: S1
 
!@3  <p><i><sup><font size="4" color="#0066FF"> </font></sup></i></p>
     c=sqrt(1.0_dp/(1.0_dp+s1**2))
     c=c+i_*s1*c  ! cos(atant)+i * sin(atant)

     atant=log(c)/i_



  END FUNCTION  atant 

  FUNCTION DtanT(S1)
    implicit none
    TYPE (my_linear_taylor) DtanT
    TYPE (my_linear_taylor), INTENT (IN) :: S1
    integer i
!@3  <p><i><sup><font size="4" color="#0066FF"> </font></sup></i></p>

    DtanT=sin(s1)/cos(s1)


  END FUNCTION  DtanT 

  FUNCTION atan2t(S1,s2)
    implicit none
    TYPE (my_linear_taylor) atan2t
    TYPE (my_linear_taylor), INTENT (IN) :: S1,s2


     atan2t=atan(s1/s2)
     atan2t%a(0)=0.0_dp

     atan2t= atan2t + atan2(real(s1%a(0),kind=dp),real(s2%a(0),kind=dp))


  END FUNCTION  atan2t 

  FUNCTION DSINT(S1)
    implicit none
    TYPE (my_linear_taylor) DSINT
    TYPE (my_linear_taylor), INTENT (IN) :: S1
    
!@3  <p><i><sup><font size="4" color="#0066FF">Lazy implementation: should be  generalized, only good to 4th order.</font></sup></i></p>
 
DSINT=(exp(i_*s1)-exp(-i_*s1))/2.0_dp/i_

  END FUNCTION  DSINT 

  FUNCTION tpsa_expt( S1 )
    implicit none
    TYPE (my_linear_taylor) tpsa_expt,t
    TYPE (my_linear_taylor), INTENT (IN) :: S1
    integer i

   
   T=1.0_DP
   tpsa_expt=1.0_DP

   do i=1,n_tpsa_exp
    T=S1*T/I
    tpsa_expt=tpsa_expt+T
   enddo
 
   
  END FUNCTION tpsa_expt 


end module my_own_linear_tpsa



 
 
integer function mypause(i)
  use precision_constants
  implicit none
  !
  ! Replaces obsolescent feature pause
  !
  integer i
  !
 
write(6,*) ' ipause=mypause(0)  '

   read(5,*) I
  mypause=i
 
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
 
 write(6,*) string
write(6,*) ' ipause=mypause(0)  ' 
 read(*,*) I
  mypauses=i
 
end function mypauses
