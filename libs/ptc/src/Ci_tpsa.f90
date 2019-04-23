!The Full Polymorphic Package
!Copyright (C) Etienne Forest
!    ndct=iabs(ndpt-ndptb)  ! 1 if coasting, otherwise 0
!    ndc2t=2*ndct  ! 2 if coasting, otherwise 0
!    nd2t=nd2-2*rf-ndc2t   !  size of harmonic oscillators minus modulated clocks
!    ndt=nd2t/2        ! ndt number of harmonic oscillators minus modulated clocks
!    nd2harm=nd2t+2*rf  !!!!  total dimension of harmonic phase space
!    ndharm=ndt+rf  !!!! total number of harmonic planes
!
MODULE c_TPSA
  !use newda
  use definition
  use file_handler
  use tree_element_MODULE
  IMPLICIT NONE
  public
  integer,private::nd2par,nd2part,nd2partt
  integer,private,target ::pos_of_delta  

  integer,private,dimension(lnv)::jfil,jfilt

  private equal,DAABSEQUAL,Dequaldacon ,equaldacon ,Iequaldacon,derive,DEQUALDACONS  !,AABSEQUAL 2002.10.17
  private pow, GETORDER,CUTORDER,getchar,GETint,GETORDERMAP  !, c_bra_v_spinmatrix
  private getdiff,getdATRA  ,mul,dmulsc,dscmul,GETintmat   !,c_spinor_spinmatrix
  private mulsc,scmul,imulsc,iscmul,DAREADTAYLORS,c_pri_c_ray
  private div,ddivsc,dscdiv,divsc,scdiv,idivsc,iscdiv,equalc_ray_r6r 
  private unaryADD,add,daddsca,dscadd,addsc,scadd,iaddsc,iscadd,print_ql
  private unarySUB,subs,dsubsc,dscsub,subsc,scsub,isubsc,iscsub,c_clean_taylors
  private c_allocda,c_killda,c_a_opt,K_opt,c_,c_allocdas,filter_part
  private dexpt,dcost,dsint,dtant,DAPRINTTAYLORS,c_clean_yu_w,mul_ql_m,mul_ql_cm
  PRIVATE GETCHARnd2,GETintnd2,dputchar,dputint, filter,check_j,c_dputint0,c_dputint0r
  private GETintnd2t,equalc_cspinor_cspinor,c_AIMAG,c_real,equalc_ray_ray,EQUALql_q,EQUALq_ql,EQUALql_i,EQUALql_ql
  PRIVATE DEQUAL,REQUAL,varf,varf001,equalc_spinor_cspinor,EQUALql_r  !,CHARINT,pbbrav,cpbbrav
  !  PUBLIC VAR,ASS
  private pbbra,liebra,full_absT,c_asstaylor,getcharnd2s,GETintnd2s,GETintk
  private shiftda,shift000,cDEQUAL,pri,rea,cfu000,alloc_DA,alloc_c_spinmatrix,cpbbra
  private alloc_c_damap,c_DPEKMAP,c_DPOKMAP,kill_c_damap,kill_c_spinmatrix,c_etcct,c_spinmatrix_mul_cray
  private EQUALspinmatrix,c_trxtaylor,powmap,POWMAPs,alloc_c_vector_field,kill_c_vector_field
  private alloc_c_normal_form,kill_c_normal_form,c_EQUALVEC,c_spinmatrix_spinmatrix,c_IdentityEQUALVEC,qua_ql
  private liebraquaternion,pow_tpsaMAP,c_concat_quaternion_ray
  private EQUALql_cmap,EQUALcmap_ql,EQUAL_complex_quaternion_c_quaternion,EQUAL_c_quaternion_complex_quaternion
  private NO,ND,ND2,NP,NDPT,NV,ndptb,rf
  integer, target :: NP,NO,ND,ND2,NDPT,NV,ndptb,rf
  private nd_used
  integer nd_used
  logical(lp):: do_linear_ac_longitudinal=.true.
  private old
!private map_mul_vec
  logical(lp) :: old  
  logical(lp),target  :: c_real_warning =.true.
  logical(lp) :: c_mess_up_vector=.false. 
  real(dp) :: a_mess=0.d0 , b_mess=1.d0
  integer :: i_piotr(3)= (/0,0,0/)

  PRIVATE null_it,Set_Up,de_Set_Up,LINE_L,RING_L,kill_DALEVEL,dealloc_DASCRATCH,set_up_level
  private insert_da,append_da,GETINTegrate,c_pek000,c_pok000,cDEQUALDACON
  private cdaddsc,cdscadd,cdsubsc,cdscsub,cdmulsc,cdscmul,cddivsc,cdscdiv
 private equalc_t_ct,equalc_ct_c,matrixMAPr,MAPmatrixr,r_MAPmatrixr,c_EQUALMAP,c_IdentityEQUALMAP
 private c_IdentityEQUALSPIN,c_pri_spinmatrix,c_pri_map,r_matrixMAPr
 private equalc_cmap_map,c_bra_v_ct,c_bra_v_dm,equalc_cvec_vec,c_expflo,c_expflo_map
 private alloc_c_factored_lie,kill_c_factored_lie,c_expflo_fac
 private c_trxspinmatrix,c_inv_as,sqrtt,alloc_c_spinor,kill_c_spinor,c_complex_spinmatrix,c_trxspinmatrixda
 private c_spinmatrix_add_spinmatrix,c_exp_spinmatrix,unarySUB_spinor,c_spinor_add_spinor,c_taylor_spinor
 private c_IdentityEQUALSPINOR,c_spinmatrix_spinor,c_logt,c_pri_factored_lie,equalc_map_cmap
 private c_expflo_fac_inv,c_logc,c_complex_spinor,c_real_spinor,GETORDERSPINMATRIX,c_pri_spinor
 private c_spinor_cmap,c_adjoint_vec,c_adjoint,c_trxtaylor_da,c_spinmatrix_sub_spinmatrix,c_spinor_cmap_tpsa
 PRIVATE CUTORDERMAP,CUTORDERspin,CUTORDERspinor,c_concat_tpsa,c_concat_spinor_ray,GETORDERquaternion
  type(C_dalevel) c_scratchda(ndumt)   !scratch levels of DA using linked list
integer, private,target :: nd2t=6,ndt=3,ndc2t=2,ndct=1,nd2harm,ndharm
!integer, private, parameter :: ndim2t=10
logical(lp), private ::   c_similarity=my_false
logical(lp) :: symp =my_false
logical(lp) :: c_normal_auto=my_true,c_verbose=my_true
integer :: spin_def_tune=1   !, private 
integer :: order_gofix=1
logical(lp) :: time_lie_choice=my_false,courant_snyder_teng_edwards=my_true,dosymp=my_false
  private copy_damap_matrix,copy_matrix_matrix,invert_22,ALLOC_33t,kill_33t,matmul_33,print_33t
  private A_OPT_C_damap,K_OPT_c_damap,equalc_t,equalt_c,daddsco,scdaddo,daddsc,scdadd,matmulr_33
private equal_real8_cmap,equal_cmap_real8,EQUAL_c_map_RAY8,EQUAL_RAY8_c_map,c_add_vf,real_mul_vec
private c_sub_vf,c_spinor_sub_spinor,matmult_33,EQUALq_i
private c_IdentityEQUALfactored,c_log_spinmatrix,c_concat_c_ray,equalc_ray_r6,equalc_r6_ray
private dotc_spinor,c_spinor_spinor,c_read_spinmatrix,c_read_map,c_concat_spinmatrix_ray
 
integer,private,parameter::ndd=6
private c_concat_vector_field_ray,CUTORDERVEC,kill_c_vector_field_fourier,alloc_c_vector_field_fourier
private complex_mul_vec,equal_c_vector_field_fourier,c_IdentityEQUALVECfourier
private c_add_map,c_sub_map,c_read_spinor,flatten_c_factored_lie_r,c_EQUALcray,c_read_quaternion
integer :: n_fourier=12,n_extra=0
logical :: remove_tune_shift=.false.
complex(dp) :: n_cai=-2*i_
integer :: complex_extra_order=0
logical :: special_extra_order_1=.true.
real(dp) :: epso_factor =1000.d0 ! for log
logical(lp):: extra_terms_log=.false. 
logical :: add_constant_part_concat=.true.,assume_c_quaternion_normalised=.true.
private EQUAL_c_spinmatrix_probe,EQUAL_c_spinmatrix_3_by_3,EQUAL_3_by_3_probe,EQUAL_probe_c_spinmatrix
private EQUAL_probe_3_by_3,equalc_cspinor_spinor,EQUAL_3_by_3_c_spinmatrix
private EQUALq_r,EQUALq_8_c,EQUALq_c_8,EQUALq,POWq,c_invq,subq,mulq,addq,alloc_c_quaternion,kill_c_quaternion
private c_pri_quaternion,CUTORDERquaternion,c_trxquaternion,EQUALq_c_r,EQUALq_r_c,mulcq,c_exp_quaternion
private equalc_quaternion_c_spinor,equalc_spinor_c_quaternion,unarySUB_q,c_trxquaternion_tpsa
private c_exp_vectorfield_on_quaternion,c_vector_field_quaternion,addql,subql,mulqdiv,powql
!private equal_map_real8,equal_map_complex8,equal_real8_map,equal_complex8_map
real(dp) dts
real(dp), private :: sj(6,6)
logical :: use_new_stochastic_normal_form=.true.
logical :: qphase=.true.

type q_linear
 complex(dp) mat(6,6)
 complex(dp)  q(0:3,0:6) 
end type q_linear



type(q_linear) q_phasor,qi_phasor

  INTERFACE assignment (=)
     MODULE PROCEDURE EQUAL
     MODULE PROCEDURE EQUALq_r
     MODULE PROCEDURE EQUALq_8_c
     MODULE PROCEDURE EQUALq_c_8
     MODULE PROCEDURE EQUALql_q
     MODULE PROCEDURE EQUALq_ql
     MODULE PROCEDURE EQUALql_i
     MODULE PROCEDURE EQUALq_i
     MODULE PROCEDURE EQUALql_r
     MODULE PROCEDURE EQUALq
     MODULE PROCEDURE EQUALq_c_r
     MODULE PROCEDURE EQUALq_r_c
     MODULE PROCEDURE EQUALql_ql
     MODULE PROCEDURE EQUAL_complex_quaternion_c_quaternion
     MODULE PROCEDURE EQUAL_c_quaternion_complex_quaternion
     MODULE PROCEDURE EQUALql_cmap
     MODULE PROCEDURE EQUALcmap_ql
     MODULE PROCEDURE cDEQUAL
     MODULE PROCEDURE DEQUAL  ! added 2002.10.17    ! check2002.10.17
     MODULE PROCEDURE REQUAL   ! added 2002.10.17   ! check2002.10.17
     MODULE PROCEDURE cDEQUALDACON
     MODULE PROCEDURE Dequaldacon
     MODULE PROCEDURE DEQUALDACONS ! C_taylor(:) = real(dp))
     MODULE PROCEDURE Iequaldacon
     MODULE PROCEDURE equalc_t_ct   !# complex Taylor
     MODULE PROCEDURE equalc_ct_c
     MODULE PROCEDURE equalc_t   !# Taylor
     MODULE PROCEDURE equalt_c
     MODULE PROCEDURE equal_real8_cmap  !# put cmap in real_8(6)
     MODULE PROCEDURE equal_cmap_real8  ! put real_8(6) in cmap
     MODULE PROCEDURE c_EQUALMAP  !#  c_damap=c_damap
    MODULE PROCEDURE c_IdentityEQUALMAP
    MODULE PROCEDURE c_IdentityEQUALSPIN
    MODULE PROCEDURE c_IdentityEQUALSPINOR   !# c_spinor = 0,1,2, or 3
    MODULE PROCEDURE c_IdentityEQUALVEC
    MODULE PROCEDURE c_IdentityEQUALfactored
    MODULE PROCEDURE c_IdentityEQUALVECfourier
    MODULE PROCEDURE EQUAL_c_spinmatrix_probe
    MODULE PROCEDURE EQUAL_c_spinmatrix_3_by_3
    module procedure EQUAL_3_by_3_c_spinmatrix
    MODULE PROCEDURE EQUAL_3_by_3_probe
    MODULE PROCEDURE EQUAL_probe_3_by_3
    MODULE PROCEDURE EQUAL_probe_c_spinmatrix

    MODULE PROCEDURE matrixMAPr
    MODULE PROCEDURE r_matrixMAPr
    MODULE PROCEDURE MAPmatrixr
    MODULE PROCEDURE r_MAPmatrixr
   !  MODULE PROCEDURE c_DPEKMAP
  !   MODULE PROCEDURE c_DPOKMAP
     MODULE PROCEDURE EQUALspinmatrix

     MODULE PROCEDURE EQUAL_c_map_RAY8   !#  c_damap=probe_8
     MODULE PROCEDURE EQUAL_RAY8_c_map   !#  probe_8=c_damap
     MODULE PROCEDURE c_EQUALcray       !#  c_ray=integer
      MODULE PROCEDURE equalc_cmap_map 
      MODULE PROCEDURE equalc_map_cmap

      MODULE PROCEDURE equalc_cvec_vec

      MODULE PROCEDURE c_EQUALVEC   !# c_vector_field=c_vector_field
      MODULE PROCEDURE equalc_cspinor_cspinor
      MODULE PROCEDURE equalc_ray_r6
      MODULE PROCEDURE equalc_ray_r6r
      MODULE PROCEDURE equalc_r6_ray
      MODULE PROCEDURE equalc_ray_ray
      MODULE PROCEDURE equal_c_vector_field_fourier
      MODULE PROCEDURE equalc_spinor_cspinor
      MODULE PROCEDURE equalc_cspinor_spinor
      MODULE PROCEDURE flatten_c_factored_lie_r !# same as flatten_c_factored_lie
      MODULE PROCEDURE equalc_quaternion_c_spinor
      MODULE PROCEDURE equalc_spinor_c_quaternion

      MODULE PROCEDURE equal_map_real8
      MODULE PROCEDURE equal_map_complex8    ! replaces c_dpokmap
      MODULE PROCEDURE equal_real8_map
      MODULE PROCEDURE equal_complex8_map   ! replaces c_dpekmap

  end  INTERFACE


  INTERFACE OPERATOR (+)
     MODULE PROCEDURE unaryADD  !@2 This is a unary operation
     MODULE PROCEDURE add
     MODULE PROCEDURE addq
     MODULE PROCEDURE addql

     MODULE PROCEDURE daddsco   !# c_damap + real(6)
     MODULE PROCEDURE scdaddo   !# real(6) + c_damap
     MODULE PROCEDURE daddsc    !# c_damap + probe_8
     MODULE PROCEDURE scdadd    !# probe_8 + c_damap

     MODULE PROCEDURE cdaddsc
     MODULE PROCEDURE daddsca
     MODULE PROCEDURE cdscadd
     MODULE PROCEDURE iscadd

     MODULE PROCEDURE dscadd
     MODULE PROCEDURE addsc
     MODULE PROCEDURE scadd
     MODULE PROCEDURE iaddsc
     MODULE PROCEDURE c_spinmatrix_add_spinmatrix
     MODULE PROCEDURE c_spinor_add_spinor  ! adding c_spinor
     MODULE PROCEDURE c_add_vf  !#  adding vector field for exp(ad VF)
     MODULE PROCEDURE c_add_map
 
  END INTERFACE



  INTERFACE OPERATOR (-)
     MODULE PROCEDURE unarySUB

     MODULE PROCEDURE subq
     MODULE PROCEDURE subql
     MODULE PROCEDURE subs
     MODULE PROCEDURE dsubsc
     MODULE PROCEDURE dscsub
     MODULE PROCEDURE cdsubsc
     MODULE PROCEDURE cdscsub
     MODULE PROCEDURE subsc
     MODULE PROCEDURE scsub
     MODULE PROCEDURE isubsc
     MODULE PROCEDURE iscsub

     MODULE PROCEDURE unarySUB_vec
     MODULE PROCEDURE unarySUB_q
     MODULE PROCEDURE unarySUB_spinor
     MODULE PROCEDURE c_spinor_sub_spinor
     MODULE PROCEDURE c_spinmatrix_sub_spinmatrix
     MODULE PROCEDURE c_sub_vf
     MODULE PROCEDURE c_sub_map
  END INTERFACE

  INTERFACE OPERATOR (*)
 !    MODULE PROCEDURE transform_vector_field_by_map
     MODULE PROCEDURE mul
     MODULE PROCEDURE mulq
     MODULE PROCEDURE mulql
     MODULE PROCEDURE mulcq
     MODULE PROCEDURE mul_ql_m
     MODULE PROCEDURE mul_ql_cm

     MODULE PROCEDURE dmulsc
     MODULE PROCEDURE cdmulsc
     MODULE PROCEDURE dscmul
     MODULE PROCEDURE cdscmul
     MODULE PROCEDURE mulsc
     MODULE PROCEDURE scmul
     MODULE PROCEDURE imulsc
     MODULE PROCEDURE iscmul


     MODULE PROCEDURE c_concat      !# c_damap o  c_damap
     MODULE PROCEDURE c_trxtaylor   !  c_taylor o  c_damap   
     MODULE PROCEDURE c_bra_v_ct     !# F.grad taylor                !v1
     MODULE PROCEDURE c_bra_v_dm     !# c_damap=(F.grad) c_damap     !v2
!     MODULE PROCEDURE c_bra_v_v     !   (exp(F.grad) H) . grad
     MODULE PROCEDURE c_spinmatrix_spinmatrix !#  Spinmatrix*Spinmatrix
     MODULE PROCEDURE c_complex_spinmatrix  !# c*Spinmatrix
     MODULE PROCEDURE c_taylor_spinor    !# taylor * spinor
     MODULE PROCEDURE c_complex_spinor   !#  complex * spinor
     MODULE PROCEDURE c_real_spinor      !#  real(dp) * spinor
     MODULE PROCEDURE c_spinor_spinor    !# spinor x spinor 
     MODULE PROCEDURE c_trxspinmatrix    !# c_spinmatrix=  c_spinmatrix * c_damap
     MODULE PROCEDURE c_trxquaternion    !# c_quaternion=  c_quaternion * c_damap
     MODULE PROCEDURE c_spinmatrix_spinor  !# matrix * spinor
     MODULE PROCEDURE c_spinor_cmap        !# spinor * c_damap
     MODULE PROCEDURE real_mul_vec  ! real(dp)*vf
     MODULE PROCEDURE complex_mul_vec  !# complex(dp)*vf
     MODULE PROCEDURE real_mul_map     !# real(dp)*c_damap
 !    MODULE PROCEDURE c_spinor_spinmatrix    !# spinor.L   spinmatrix
     MODULE PROCEDURE c_vector_field_quaternion  !# (f.grad + q)quaternion
                                                 !  bra_v_q(f,quaternion) = f.grad quaternion
!    uses  map_mul_vec 
     MODULE PROCEDURE map_mul_vec_q  !   c_damap * c_vector_field means "transform field by map" 
  END INTERFACE

  INTERFACE OPERATOR (.o.) 
  !
     module procedure c_concat_c_ray    !# c_taylor .o. c_ray  
     module procedure c_concat_map_ray  !# c_ray= c_damap .o. c_ray
     module procedure c_trxtaylor_da   !# c_taylor= c_taylor .o. c_damap
     module procedure c_concat_tpsa     !# c_damap .o.  c_damap
     module procedure c_concat_spinor_ray   !# c_spinor= c_spinor .o.  c_ray
! not overloaded
   !  module procedure c_concat_spinmatrix_ray !# c_spinmatrix= c_spinmatrix .o.  c_ray
   !  module procedure c_concat_quaternion_ray !# c_quaternion= c_quaternion .o.  c_ray
   !  Instead these overloaded
     MODULE PROCEDURE c_spinmatrix_mul_cray !# c_ray%s = spin_matrix.o.c_ray%x  c_ray%S(1:3)
     MODULE PROCEDURE c_quaternion_mul_cray !# c_ray%q = c_quaternion.o.c_ray%x  c_ray%q  c_quaternion.o.c_ray%x**(-1)

  !   module procedure c_concat_vector_field_ray  !# c_vector_field .o.  cray
     module procedure c_trxspinmatrixda     !# c_spinmatrix=c_spinmatrix .o.  c_damap
     module procedure c_trxquaternion_tpsa   !# c_quaternion=  c_quaternion .o. c_damap
     MODULE PROCEDURE c_spinor_cmap_tpsa        !# spinor * c_damap
   END INTERFACE 
 
  INTERFACE OPERATOR (.oo.) 
     module procedure pow_tpsaMAP
  END INTERFACE 

  INTERFACE OPERATOR (/)
     MODULE PROCEDURE div
     MODULE PROCEDURE cddivsc
     MODULE PROCEDURE ddivsc
     MODULE PROCEDURE cdscdiv
     MODULE PROCEDURE dscdiv
     MODULE PROCEDURE divsc
     MODULE PROCEDURE scdiv
     MODULE PROCEDURE idivsc
     MODULE PROCEDURE iscdiv
     MODULE PROCEDURE mulqdiv
  END INTERFACE

  INTERFACE OPERATOR (**)
     MODULE PROCEDURE POW
     MODULE PROCEDURE POWQ
     MODULE PROCEDURE POWql
     MODULE PROCEDURE powmap
     MODULE PROCEDURE powmaps
  END INTERFACE



  ! New Operators

  INTERFACE OPERATOR (.cmono.)
     MODULE PROCEDURE c_dputint0   !@1 &nbsp; single integer
     MODULE PROCEDURE c_dputint0r   !@1 &nbsp; single integer
     MODULE PROCEDURE dputintr    !@1 &nbsp; Accepts J(nv)
     MODULE PROCEDURE dputcharr   !@1 &nbsp; Accepts String such as '12'
     MODULE PROCEDURE dputint    !@1 &nbsp; Accepts J(nv)
     MODULE PROCEDURE dputchar   !@1 &nbsp; Accepts String such as '12'
  END INTERFACE

  INTERFACE OPERATOR (.var.)
     MODULE PROCEDURE varf        !@1 &nbsp; replaces var (overloads DAVAR)
     MODULE PROCEDURE varf001       !@1 replaces var001
  END INTERFACE

  INTERFACE OPERATOR (.d.)
     MODULE PROCEDURE getdiff    !@1 takes derivatives
  END INTERFACE

  INTERFACE OPERATOR (.i.)
     MODULE PROCEDURE GETINTegrate    !@1 takes anti-derivatives
  END INTERFACE


  INTERFACE OPERATOR (.SUB.)
     MODULE PROCEDURE GETORDER
     MODULE PROCEDURE getchar
     MODULE PROCEDURE GETint
     MODULE PROCEDURE GETORDERMAP  ! with negative integer map.sub.i the spin is handled with iabs(i)-1
     MODULE PROCEDURE GETORDERSPINMATRIX ! FOR SPIN MATRICES
     MODULE PROCEDURE GETORDERquaternion
  END INTERFACE


  INTERFACE OPERATOR (.index.)
    MODULE PROCEDURE GETintmat
  END INTERFACE


 

  INTERFACE OPERATOR (.PAR.)
     MODULE PROCEDURE getcharnd2
     MODULE PROCEDURE GETintnd2
  END INTERFACE

  INTERFACE OPERATOR (.part.)
     MODULE PROCEDURE GETintnd2t
  END INTERFACE


  INTERFACE OPERATOR (<=)
     MODULE PROCEDURE getcharnd2s
     MODULE PROCEDURE GETintnd2s
     MODULE PROCEDURE GETintk
  END INTERFACE

  INTERFACE OPERATOR (.CUT.)
     MODULE PROCEDURE CUTORDER
     MODULE PROCEDURE CUTORDERMAP
     MODULE PROCEDURE CUTORDERvec
     MODULE PROCEDURE CUTORDERspin
     MODULE PROCEDURE CUTORDERspinor
     MODULE PROCEDURE CUTORDERquaternion
  END INTERFACE



  INTERFACE OPERATOR (.dot.)
     MODULE PROCEDURE dotc_spinor    ! Used internally primarily
  END INTERFACE

  INTERFACE OPERATOR (.K.)
     MODULE PROCEDURE getdATRA    ! Used internally primarily
  END INTERFACE

  INTERFACE OPERATOR (.pb.)
     MODULE PROCEDURE pbbra
  END INTERFACE
  
    INTERFACE OPERATOR (.cpb.)
     MODULE PROCEDURE cpbbra
    END INTERFACE
    
  INTERFACE OPERATOR (.lb.)
  !
     MODULE PROCEDURE liebra    !#  F=<G,H> includes spin now 
     MODULE PROCEDURE liebraquaternion !# used by liebra
  END INTERFACE

 ! INTERFACE getvf
 !    MODULE PROCEDURE pbbrav
  !end INTERFACE getvf

  !INTERFACE cgetvf
  !   MODULE PROCEDURE cpbbrav
  !end INTERFACE cgetvf

  ! intrisic functions overloaded

  INTERFACE q_part
     MODULE PROCEDURE qua_ql
  end INTERFACE q_part


  INTERFACE c_phasor
     MODULE PROCEDURE from_phasor
  end INTERFACE c_phasor


  INTERFACE ci_phasor
     MODULE PROCEDURE to_phasor
  end INTERFACE ci_phasor

  ! Exponential of Lie Operators

  INTERFACE clean
!     MODULE PROCEDURE c_clean
     MODULE PROCEDURE c_clean_spinor
     MODULE PROCEDURE c_clean_taylor
     MODULE PROCEDURE c_clean_spinmatrix
     MODULE PROCEDURE c_clean_damap
     MODULE PROCEDURE c_clean_vector_field
     MODULE PROCEDURE c_clean_yu_w
     MODULE PROCEDURE c_clean_quaternion
     MODULE PROCEDURE c_clean_taylors
  end INTERFACE clean
  ! Exponential of Lie Operators

 
  INTERFACE c_simil
     MODULE PROCEDURE c_adjoint
     MODULE PROCEDURE c_adjoint_vec
  end INTERFACE c_simil

  INTERFACE texp_inv
     MODULE PROCEDURE c_expflo_fac_inv
  END INTERFACE

  INTERFACE exp_inv
     MODULE PROCEDURE c_expflo_fac_inv   ! v9
  END INTERFACE

  INTERFACE real
     MODULE PROCEDURE c_real
  END INTERFACE

  INTERFACE aimag
     MODULE PROCEDURE c_AIMAG
  END INTERFACE

  INTERFACE exp
     MODULE PROCEDURE c_expflo    ! flow on c_taylor     !v3
     MODULE PROCEDURE c_expflo_map                       !v4
     MODULE PROCEDURE c_expflo_fac                       !v7
     MODULE PROCEDURE c_exp_spinmatrix
     MODULE PROCEDURE c_exp_vectorfield_on_quaternion    !v6
     MODULE PROCEDURE c_exp_quaternion
     MODULE PROCEDURE exp_ad    ! exp(<F,>)F    F is a vector field  !v5
  END INTERFACE


  INTERFACE texp
     MODULE PROCEDURE c_expflo   
     MODULE PROCEDURE c_expflo_map
     MODULE PROCEDURE c_expflo_fac
     MODULE PROCEDURE c_exp_spinmatrix
     MODULE PROCEDURE c_exp_vectorfield_on_quaternion
     MODULE PROCEDURE c_exp_quaternion
     MODULE PROCEDURE exp_ad    ! exp(<F,>)F    F is a vector field
  END INTERFACE
  
  INTERFACE abs
     MODULE PROCEDURE DAABSEQUAL  ! remove 2002.10.17
  END INTERFACE
  INTERFACE dabs
     MODULE PROCEDURE DAABSEQUAL  ! remove 2002.10.17
  END INTERFACE

  INTERFACE exp
     MODULE PROCEDURE dexpt
  END INTERFACE
  INTERFACE dexp
     MODULE PROCEDURE dexpt
  END INTERFACE
  INTERFACE cexp
     MODULE PROCEDURE dexpt
  END INTERFACE
  INTERFACE cdexp
     MODULE PROCEDURE dexpt
  END INTERFACE

  INTERFACE log
     MODULE PROCEDURE c_logt
     MODULE PROCEDURE c_logc
     MODULE PROCEDURE c_logf  !# log of a map see subroutine c_flofacg
     MODULE PROCEDURE c_log_spinmatrix  !#  spinor=log(s)
! c_logf_spin is not overloaded
  END INTERFACE


  INTERFACE cos
     MODULE PROCEDURE dcost
  END INTERFACE
  INTERFACE cdcos
     MODULE PROCEDURE dcost
  END INTERFACE
  INTERFACE dcos
     MODULE PROCEDURE dcost
  END INTERFACE
  INTERFACE ccos
     MODULE PROCEDURE dcost
  END INTERFACE




  INTERFACE sin
     MODULE PROCEDURE dsint
  END INTERFACE
  INTERFACE cdsin
     MODULE PROCEDURE dsint
  END INTERFACE
  INTERFACE ccsin
     MODULE PROCEDURE dsint
  END INTERFACE
  INTERFACE dsin
     MODULE PROCEDURE dsint
  END INTERFACE

  INTERFACE sqrt
     MODULE PROCEDURE sqrtt
  END INTERFACE


  INTERFACE tan
     MODULE PROCEDURE dtant
  END INTERFACE

  INTERFACE dtan
     MODULE PROCEDURE dtant
  END INTERFACE

  ! Non-intrisic Functions

  INTERFACE c_pek
     MODULE PROCEDURE c_pek000 !  
  END INTERFACE

  INTERFACE c_pok
     MODULE PROCEDURE c_pok000  !  
  END INTERFACE

  INTERFACE shiftda
     MODULE PROCEDURE shift000  ! not private
  END INTERFACE

  !  INTERFACE var
  !     MODULE PROCEDURE var000  ! not private
  !     MODULE PROCEDURE var001  ! not private
  !  END INTERFACE

  INTERFACE cfu
     MODULE PROCEDURE c_cfu000  ! not private
  END INTERFACE

  INTERFACE full_abs
     MODULE PROCEDURE full_absT
  END INTERFACE

    INTERFACE daread
       MODULE PROCEDURE c_rea
       module procedure c_read_spinmatrix
       module procedure c_read_map
       MODULE PROCEDURE c_read_spinor
       MODULE PROCEDURE DAREADTAYLORS
       MODULE PROCEDURE c_read_quaternion
    END INTERFACE

    INTERFACE read
       MODULE PROCEDURE c_rea
       module procedure c_read_spinmatrix
       module procedure c_read_map
       MODULE PROCEDURE c_read_spinor
       MODULE PROCEDURE DAREADTAYLORS
       MODULE PROCEDURE c_read_quaternion
    END INTERFACE

    INTERFACE daprint
       MODULE PROCEDURE c_pri
       MODULE PROCEDURE c_pri_spinmatrix
       MODULE PROCEDURE c_pri_map
       MODULE PROCEDURE c_pri_vec
       MODULE PROCEDURE c_pri_factored_lie
       MODULE PROCEDURE c_pri_spinor
       MODULE PROCEDURE print_33t
       MODULE PROCEDURE c_pri_quaternion
      MODULE PROCEDURE DAPRINTTAYLORS
      MODULE PROCEDURE c_pri_c_ray
    END INTERFACE

    INTERFACE print
       MODULE PROCEDURE c_pri
       MODULE PROCEDURE c_pri_spinmatrix
       MODULE PROCEDURE c_pri_map
       MODULE PROCEDURE c_pri_vec
       MODULE PROCEDURE c_pri_factored_lie
       MODULE PROCEDURE c_pri_spinor
       MODULE PROCEDURE print_33t
       MODULE PROCEDURE DAPRINTTAYLORS
       MODULE PROCEDURE c_pri_quaternion
       MODULE PROCEDURE print_ql
      MODULE PROCEDURE c_pri_c_ray
    END INTERFACE


  ! Constructors and Destructors

  INTERFACE alloc
     MODULE PROCEDURE c_allocda
     MODULE PROCEDURE c_a_opt
     MODULE PROCEDURE A_OPT_c_damap
     MODULE PROCEDURE c_allocdas
     MODULE PROCEDURE alloc_c_spinmatrix
     MODULE PROCEDURE alloc_c_vector_field_fourier
     MODULE PROCEDURE alloc_c_damap
     MODULE PROCEDURE alloc_c_vector_field
     MODULE PROCEDURE alloc_c_factored_lie
     MODULE PROCEDURE alloc_c_normal_form
     MODULE PROCEDURE alloc_c_spinor
     MODULE PROCEDURE alloc_c_yu_w
     MODULE PROCEDURE alloc_c_quaternion
  END INTERFACE

  INTERFACE ALLOC_nn
     MODULE PROCEDURE ALLOC_33t
  END INTERFACE

  INTERFACE kill_nn
     MODULE PROCEDURE kill_33t
  END INTERFACE


  INTERFACE matmul_nn
     MODULE PROCEDURE matmul_33
  END INTERFACE

  INTERFACE matmulr_nn
     MODULE PROCEDURE matmulr_33
     MODULE PROCEDURE matmult_33
  END INTERFACE

  INTERFACE KILL
     MODULE PROCEDURE c_killda
     MODULE PROCEDURE c_killdas
     MODULE PROCEDURE K_opt
     MODULE PROCEDURE K_OPT_c_damap
     MODULE PROCEDURE kill_c_damap
     MODULE PROCEDURE kill_c_spinmatrix
     MODULE PROCEDURE kill_c_vector_field
     MODULE PROCEDURE kill_c_factored_lie
     MODULE PROCEDURE kill_c_vector_field_fourier
     MODULE PROCEDURE kill_c_normal_form
     MODULE PROCEDURE kill_c_spinor
     MODULE PROCEDURE kill_c_yu_w
     MODULE PROCEDURE kill_c_quaternion
  END INTERFACE

  INTERFACE alloctpsa
     MODULE PROCEDURE c_allocda
  END INTERFACE

  INTERFACE KILLtpsa
     MODULE PROCEDURE c_killda
  END INTERFACE

  INTERFACE MAKESO3
     MODULE PROCEDURE quaternion_to_matrix_in_c_damap
     MODULE PROCEDURE q_linear_to_matrix
     MODULE PROCEDURE q_linear_to_3_by_3_by_6
  END INTERFACE
  ! management routines

  INTERFACE ass
     MODULE PROCEDURE c_asstaylor   !2000.12.25
  END INTERFACE

  INTERFACE AVERAGE
     MODULE PROCEDURE c_AVERAGE   !2000.12.25
  END INTERFACE




CONTAINS






 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine c_get_indices(n,mf)
!#general: informationc_%
!# In the arrary n(>11), important parameters of the normal
!# form can be retrieved.
!# If mf/=0, they are printed on file mf.
!# n(1)=NO
!# n(2)=ND
!# n(3)=ND2
!# n(4)=NV
!# n(5)=Ndpt
!# n(6)=ndptb
!# n(7)=np
!# n(8)=rf*2
!# n(9)=ndc2t
!# n(10)=nd2t
!# n(11)=nd2harm
 implicit none
 integer n(:),mf

  if(size(n)<11) then
   write(6,*) " index array to small in c_get_indices "
   stop
  endif
  n(1)=NO
  n(2)=ND
  n(3)=ND2
  n(4)=NV
  n(5)=Ndpt
  n(6)=ndptb
  n(7)=np
  n(8)=rf*2
  n(9)=ndc2t
  n(10)=nd2t
  n(11)=nd2harm

if(mf/=0) then
  write(mf,"(11(a7))") " NO   ","   ND ","   ND2","   NV ","  NDPT "," NDPTB ","   NP ","    RF "," NDC2T ", &
                      "  ND2T ","  HARM "
  write(mf,"(11(3x,i2,2x))") n(1:11)
endif

end subroutine c_get_indices

  subroutine c_count_taylor(n,ns,ne)
!#restricted : information
!# Counts number of c_taylor allocated;
!# it is used for debugging purposes.
    implicit none
    integer n,ns,ne,i
    call c_count_da(n)
    ns=0
    do i=1,ndumt
       ns=c_scratchda(i)%n+ns
    enddo
    ne=n-ns
  end subroutine c_count_taylor

  FUNCTION scdadd( S2,S1  )
    implicit none
    TYPE (probe_8) scdadd
    TYPE (c_damap), INTENT (IN) :: S1
    type(probe) , INTENT (IN) :: S2
    integer localmaster,i,j,nd2h,dc
    type(taylor) d

    call alloc(d)

    !   call ass(scdadd)
    scdadd%u=my_false
    scdadd%E_ij=0.0_dp
    scdadd%nac=s2%nac
    scdadd%use_q=s2%use_q  
    if(doing_ac_modulation_in_ptc) then
       dc=2*rf
    else
       dc=0
    endif

    if(c_%ndpt==0) then    ! 1                    
      nd2h=nd2t
       do i=1,nd2h                 !   from 1-4 or 1-6 (if ndpt=0)
          localmaster=master
          call ass(scdadd%x(i))
          !       scdadd%x(i)=s1%m%v(i)+s2%x(i)

          d=s1%V(i)+s2%x(i)-(s1%V(i).sub.'0')
          scdadd%x(i)=d
          master=localmaster
       enddo
  
       do i=nd2h+1,6
          localmaster=master
          call ass(scdadd%x(i))
!  !       if((c_%npara==5+dc).AND.I==5) then   ! npr
          if(((c_%npara==5+dc).AND.I==5+ndpt_bmad).or.((c_%npara==3+dc).AND.I==5+ndpt_bmad)) then   ! npr
             scdadd%x(i)=s2%x(i)+(1.0_dp.mono.c_%npara)
          else
             scdadd%x(i)=s2%x(i)
          endif
          master=localmaster
       enddo

    else        ! 1
      nd2h=nd2t+2
       do i=1,nd2h                 !   from 1-4 or 1-6 (if ndpt=0)
          localmaster=master
          call ass(scdadd%x(i))
          !       scdadd%x(i)=s1%m%v(i)+s2%x(i)

          d=s1%V(i)+s2%x(i)-(s1%V(i).sub.'0')
          scdadd%x(i)=d
          master=localmaster
       enddo
    endif       ! 1


 do i=1,scdadd%nac
    localmaster=master
    call ass(scdadd%AC(i)%x(1))
    !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
     j=2*scdadd%nac-(2*i-1)
     d=s1%V(C_%ND2-j) +addclock*s2%AC(i)%x(1)
    scdadd%ac(i)%x(1)=d
    master=localmaster
    localmaster=master
     j=2*scdadd%nac-(2*i)
    call ass(scdadd%AC(i)%x(2))
    !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
    d=s1%V(C_%ND2-j) +addclock*s2%AC(i)%x(2)
    scdadd%ac(i)%x(2)=d
    master=localmaster
    localmaster=master
    call ass(scdadd%AC(i)%om)
!    call ass(scdadd%AC%t)
    !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
    scdadd%AC(i)%om=s2%AC(i)%om
    scdadd%AC(i)%t=s2%AC(i)%t
    master=localmaster
enddo
    !    endif


 !   if(use_quaternion)   THEN
       DO J=0,3
          localmaster=master
          call ass(scdadd%q%x(J))
          d=S1%q%x(j)
          scdadd%q%x(J)=d
          master=localmaster
       ENDDO
!else
    DO I=1,3
       !          call ass(scdadd%s%x(i))
       DO J=1,3
          localmaster=master
          call ass(scdadd%s(J)%x(i))
          d=S1%S%s(I,J)
          scdadd%s(J)%x(i)=d
          master=localmaster
       ENDDO

    ENDDO

!endif

    scdadd%e_ij=s1%e_ij
    scdadd%x0(1:6)=s2%x

    call kill(d)

  END FUNCTION scdadd 


  FUNCTION daddsc( S1,S2  )
    implicit none
    TYPE (probe_8) daddsc
    TYPE (c_damap), INTENT (IN) :: S1
    type(probe) , INTENT (IN) :: S2
    integer localmaster,i,j,nd2h,dc
    type(taylor) d

    call alloc(d)

    !   call ass(daddsc)
    daddsc%u=my_false
    daddsc%E_ij=0.0_dp
    daddsc%nac=s2%nac

    daddsc%nac=s2%nac
    daddsc%use_q=s2%use_q
    if(doing_ac_modulation_in_ptc) then
       dc=2*rf
    else
       dc=0
    endif

    if(c_%ndpt==0) then    ! 1                    
      nd2h=nd2t
       do i=1,nd2h                 !   from 1-4 or 1-6 (if ndpt=0)
          localmaster=master
          call ass(daddsc%x(i))
          !       scdadd%x(i)=s1%m%v(i)+s2%x(i)

          d=s1%V(i)+s2%x(i)-(s1%V(i).sub.'0')
          daddsc%x(i)=d
          master=localmaster
       enddo
  
       do i=nd2h+1,6
          localmaster=master
          call ass(daddsc%x(i))
!  !       if((c_%npara==5+dc).AND.I==5) then   ! npr
          if(((c_%npara==5+dc).AND.I==5+ndpt_bmad).or.((c_%npara==3+dc).AND.I==5+ndpt_bmad)) then   ! npr
             daddsc%x(i)=s2%x(i)+(1.0_dp.mono.c_%npara)
          else
             daddsc%x(i)=s2%x(i)
          endif
          master=localmaster
       enddo

    else        ! 1
      nd2h=nd2t+2
       do i=1,nd2h                 !   from 1-4 or 1-6 (if ndpt=0)
          localmaster=master
          call ass(daddsc%x(i))
          !       scdadd%x(i)=s1%m%v(i)+s2%x(i)

          d=s1%V(i)+s2%x(i)-(s1%V(i).sub.'0')
          daddsc%x(i)=d
          master=localmaster
       enddo
    endif       ! 1


do i=1,daddsc%nac
    localmaster=master
    call ass(daddsc%AC(i)%x(1))
    !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
     j=2*daddsc%nac-(2*i-1)
     d=s1%V(C_%ND2-j) +addclock*s2%AC(i)%x(1)
    daddsc%ac(i)%x(1)=d
    master=localmaster
    localmaster=master
     j=2*daddsc%nac-(2*i)
    call ass(daddsc%AC(i)%x(2))
    !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
    d=s1%V(C_%ND2-j) +addclock*s2%AC(i)%x(2)
    daddsc%ac(i)%x(2)=d
    master=localmaster
    localmaster=master
    call ass(daddsc%AC(i)%om)
!    call ass(daddsc%AC%t)
    !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
    daddsc%AC(i)%om=s2%AC(i)%om
    daddsc%AC(i)%t=s2%AC(i)%t
    master=localmaster
    !    endif
enddo
    DO I=1,3
       !          call ass(scdadd%s%x(i))
       DO J=1,3
          localmaster=master
          call ass(daddsc%s(J)%x(i))
          d=S1%S%s(I,J)
          daddsc%s(J)%x(i)=d
          master=localmaster
       ENDDO

    ENDDO

       !          call ass(scdadd%s%x(i))
       DO J=0,3
          localmaster=master
          call ass(daddsc%q%x(J))
          d=S1%q%x(j)
          daddsc%q%x(J)=d
          master=localmaster
       ENDDO

    daddsc%e_ij=s1%e_ij
    daddsc%x0(1:6)=s2%x
    call kill(d)

  END FUNCTION daddsc 

 

  FUNCTION daddsco( S1, S2 )
    implicit none
    TYPE (real_8) daddsco(ndd)
    TYPE (c_damap), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: S2(ndd)
    type(taylor) t
    integer localmaster,nd2h,i
 
     nd2h=nd2-2*rf 
 
 !   do i=1,nd2
 !      localmaster=master
 !      call ass(daddsco(i))
 !      daddsco(i)=s1%v(i)+s2(i)-(s1%v(i).sub.'0')
 !      master=localmaster
 !   enddo
 !   do i=nd2+1,ndd
 !      localmaster=master
 !      call ass(daddsco(i))
 !      if(nd2<=4.and.(c_%npara==3.or.c_%npara==5.or.c_%npara==8).and.i==5+ndpt_bmad) then
 !         If(ndpt_bmad==0) then
 !          if(nd2==4) daddsco(i)=s2(i)+(1.0_dp.mono.'00001')
 !          if(nd2==2) daddsco(i)=s2(i)+(1.0_dp.mono.'001')
 !         endif
 !      else
 !         daddsco(i)=s2(i)
 !      endif
 !      master=localmaster
 !   enddo


    call alloc(t)

    do i=1,nd2h
       localmaster=master
       call ass(daddsco(i))
       t= s1%v(i)+s2(i)-(s1%v(i).sub.'0')
       daddsco(i)=t
       master=localmaster
    enddo
    do i=nd2h+1,ndd
       localmaster=master
       call ass(daddsco(i))
       if(nd2==4.and.(c_%npara==5.or.c_%npara==8).and.i==5+ndpt_bmad) then
          daddsco(i)=s2(i)+(1.0_dp.mono.'00001')
       elseif(nd2==2.and.(c_%npara==3.or.c_%npara==6).and.i==5+ndpt_bmad) then
          daddsco(i)=s2(i)+(1.0_dp.mono.'001')
       else
          daddsco(i)=s2(i)
       endif
       master=localmaster
    enddo


    call kill(t)

  END FUNCTION daddsco

  FUNCTION scdaddo( S2,S1  )
    implicit none
    TYPE (real_8) scdaddo(ndd)
    TYPE (c_damap), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: S2(ndd)
    integer localmaster  ,nd2h,i
    type(taylor) t

      nd2h=nd2-2*rf 

    call alloc(t)
    do i=1,nd2h
       localmaster=master
       call ass(scdaddo(i))
       t=s1%v(i)+s2(i)-(s1%v(i).sub.'0')
       scdaddo(i)=t
       master=localmaster
    enddo
    do i=nd2h+1,ndd
       localmaster=master
       call ass(scdaddo(i))
       if(nd2==4.and.(c_%npara==5.or.c_%npara==8).and.i==5+ndpt_bmad) then
          scdaddo(i)=s2(i)+(1.0_dp.mono.'00001')
       elseif(nd2==2.and.(c_%npara==3.or.c_%npara==6).and.i==5+ndpt_bmad) then
          scdaddo(i)=s2(i)+(1.0_dp.mono.'001')
       else
          scdaddo(i)=s2(i)
       endif
       master=localmaster
    enddo

    call kill(t)

  END FUNCTION scdaddo

  FUNCTION unaryADD( S1 )
    implicit none
    TYPE (c_taylor) unaryADD
    TYPE (c_taylor), INTENT (IN) :: S1
    integer localmaster

    IF(.NOT.C_STABLE_DA) then
     unaryADD%i=0
     RETURN
    endif

    localmaster=c_master

    !    call check(s1)
    call ass(unaryADD)

    unaryADD=s1

    c_master=localmaster

  END FUNCTION unaryADD

  FUNCTION unarySUB( S1 )
    implicit none
    TYPE (c_taylor) unarySUB
    TYPE (c_taylor), INTENT (IN) :: S1
    integer localmaster

    IF(.NOT.C_STABLE_DA) then
     unarySUB%i=0
     RETURN
    endif

    localmaster=c_master

    !    call check(s1)
    call ass(unarySUB)

    ! unarySUB=(-one)*s1
    !    if(old) then
    call c_dacmu(s1%i,(-1.0_dp,0.0_dp),c_temp%i)
    call c_dacop(c_temp%i,unarySUB%i)
    !   else
    !      call newc_dacmu(s1%j,-one,unarySUB%j)
    !      !  call newc_dacmu(s1%j,-one,c_temp%i%il)
    !      !  call newc_dacop(c_temp%i%il,unarySUB%j)
    !   endif
    c_master=localmaster

  END FUNCTION unarySUB


  FUNCTION unarySUB_vec( S1 )
    implicit none
    TYPE (c_vector_field) unarySUB_vec
    TYPE (c_vector_field), INTENT (IN) :: S1
    integer localmaster,i

    IF(.NOT.C_STABLE_DA) then
     unarySUB_vec%v%i=0
     RETURN
    endif
    localmaster=c_master

    !    call check(s1)
    unarySUB_vec%n=s1%n
    call c_ass_vector_field(unarySUB_vec)
    
    do i=1,s1%n
     unarySUB_vec%v(i)=-s1%v(i)
    enddo

    if(use_quaternion)   THEN
     unarySUB_vec%q=-s1%q
    endif
    unarySUB_vec%nrmax=s1%nrmax
    unarySUB_vec%eps=s1%eps

    c_master=localmaster

  END FUNCTION unarySUB_vec



  FUNCTION unarySUB_q( S1 )
    implicit none
    TYPE (c_quaternion) unarySUB_q
    TYPE (c_quaternion), INTENT (IN) :: S1
    integer localmaster,i

    IF(.NOT.C_STABLE_DA) then
     unarySUB_q%x(0)%i=0
     RETURN
    endif
    localmaster=c_master

    call c_ass_quaternion(unarySUB_q)
    
    do i=0,3
     unarySUB_q%x(i)=-s1%x(i)
    enddo


    c_master=localmaster

  END FUNCTION unarySUB_q
  
  FUNCTION unarySUB_spinor( S1 )
    implicit none
    TYPE (c_spinor) unarySUB_spinor
    TYPE (c_spinor), INTENT (IN) :: S1
    integer localmaster,i

    IF(.NOT.C_STABLE_DA) then
     unarySUB_spinor%v%i=0
     RETURN
    endif
    localmaster=c_master


    call c_ass_spinor(unarySUB_spinor)
    
    do i=1,3
     unarySUB_spinor%v(i)=-s1%v(i)
    enddo


    c_master=localmaster

  END FUNCTION unarySUB_spinor

 subroutine normalise_spinor(n1) 
!#general: manipulation
!# The c_spinor n1 is normalised.
    implicit none
    type(c_spinor), intent(inout) :: n1
    type(c_taylor) norm

    call alloc(norm);

    norm=n1.dot.n1
    norm=1.0_dp/sqrt(norm)
     n1=norm*n1

    call kill(norm);

 end subroutine normalise_spinor 

  subroutine orthogonalise_spin_matrix(s) 
!#general: manipulation
!# The c_spinmatrix s is normalised, that is to say
!# forced into O(3).

    implicit none
    TYPE (c_spinmatrix) s
    type(c_spinor) n1,n2,n3
    type(c_taylor) norm
    integer i
    call alloc(norm);
    call alloc(n1);  call alloc(n2);  call alloc(n3);
    do i=1,3
     n1%v(i)=s%s(1,i)
     n2%v(i)=s%s(2,i)
    enddo

    call normalise_spinor(n1)

     norm=-(n1.dot.n2)

     n2=n2+norm*n1

    call normalise_spinor(n2)

    n3=n1*n2


    call normalise_spinor(n3)

     do i=1,3
     s%s(1,i)=n1%v(i)
     s%s(2,i)=n2%v(i)
     s%s(3,i)=n3%v(i)
    enddo


    call kill(norm);
    call kill(n1);  call kill(n2);  call kill(n3);

  end subroutine orthogonalise_spin_matrix

  FUNCTION dotc_spinor( S1, S2 ) ! spin routine
    implicit none
    TYPE (c_taylor) dotc_spinor
    TYPE (c_spinor), INTENT (IN) :: S1,S2

    INTEGER I
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=c_master


    !    call checkdamap(s1)
    call c_asstaylor(dotc_spinor)

    dotc_spinor=0.0_dp

    DO I=1,3
       dotc_spinor=dotc_spinor+s1%v(i)*s2%v(i)
    ENDDO


    c_master=localmaster

  END FUNCTION dotc_spinor

  SUBROUTINE  c_maketree(S1,s2)
    implicit none
    type (c_taylor),INTENT(IN)::S1
    type (c_taylor),INTENT(inOUT):: s2
    IF(.NOT.C_STABLE_DA) RETURN

    !    if(old) then
    call c_mtree((/s1%i/),1,(/s2%i/),1)
    !   else
    !      call newc_dacop(s1%j,s2%j)
    !   endif
  END SUBROUTINE c_maketree

  SUBROUTINE  c_allocda(S1)
!*
    implicit none
    type (c_taylor),INTENT(INOUT)::S1

    !    IF(first_time) THEN
    IF(c_last_tpsa==0) THEN
 
      write(6,*) " No TPSA package ever initialized c_allocda " 
 
       ! call !write_e(111)
    ENDIF
    !    if(old) then
    s1%i=0
    call c_etall1(s1%i)
    !    else
    !       call nullnewda(s1%j)
    !       call allocnewda(s1%j)
    !    endif
  END SUBROUTINE c_allocda

 SUBROUTINE  alloc_c_quaternion(S2)
    implicit none
    type (c_quaternion),INTENT(INOUT)::S2
    integer i
     do i=0,3
      call alloc(s2%x(i))
    enddo
    
  END SUBROUTINE alloc_c_quaternion

 SUBROUTINE  kill_c_quaternion(S2)
    implicit none
    type (c_quaternion),INTENT(INOUT)::S2
    integer i
     do i=0,3
      call kill(s2%x(i))
    enddo

  END SUBROUTINE kill_c_quaternion

  SUBROUTINE  c_a_opt(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
!*
    implicit none
    type (c_taylor),INTENT(INout)::S1,S2
    type (c_taylor),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
    call c_allocda(s1)
    call c_allocda(s2)
    if(present(s3)) call c_allocda(s3)
    if(present(s4)) call c_allocda(s4)
    if(present(s5)) call c_allocda(s5)
    if(present(s6)) call c_allocda(s6)
    if(present(s7)) call c_allocda(s7)
    if(present(s8)) call c_allocda(s8)
    if(present(s9)) call c_allocda(s9)
    if(present(s10))call c_allocda(s10)
  END SUBROUTINE c_a_opt

  SUBROUTINE  K_OPT(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
!*
    implicit none
    type (c_taylor),INTENT(INout)::S1,S2
    type (c_taylor),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
    call c_killda(s1)
    call c_killda(s2)
    if(present(s3)) call c_killda(s3)
    if(present(s4)) call c_killda(s4)
    if(present(s5)) call c_killda(s5)
    if(present(s6)) call c_killda(s6)
    if(present(s7)) call c_killda(s7)
    if(present(s8)) call c_killda(s8)
    if(present(s9)) call c_killda(s9)
    if(present(s10))call c_killda(s10)
  END SUBROUTINE K_opt


  SUBROUTINE  A_opt_c_damap(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
!*
    implicit none
    type (c_damap),INTENT(INout)::S1,S2
    type (c_damap),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
    call alloc(s1)
    call alloc(s2)
    if(present(s3)) call alloc(s3)
    if(present(s4)) call alloc(s4)
    if(present(s5)) call alloc(s5)
    if(present(s6)) call alloc(s6)
    if(present(s7)) call alloc(s7)
    if(present(s8)) call alloc(s8)
    if(present(s9)) call alloc(s9)
    if(present(s10))call alloc(s10)
  END SUBROUTINE A_opt_c_damap


  SUBROUTINE  K_OPT_c_damap(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
!*
    implicit none
    type (c_damap),INTENT(INout)::S1,S2
    type (c_damap),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
    call KILL(s1)
    call KILL(s2)
    if(present(s3)) call KILL(s3)
    if(present(s4)) call KILL(s4)
    if(present(s5)) call KILL(s5)
    if(present(s6)) call KILL(s6)
    if(present(s7)) call KILL(s7)
    if(present(s8)) call KILL(s8)
    if(present(s9)) call KILL(s9)
    if(present(s10))call KILL(s10)
  END SUBROUTINE K_OPT_c_damap

  SUBROUTINE  c_allocdas(S1,k)
!*
    implicit none
    type (c_taylor),INTENT(INOUT),dimension(:)::S1
    INTEGER,optional,INTENT(IN)::k
    INTEGER J,i,N

    if(present(k)) then
       I=LBOUND(S1,DIM=1)
       N=LBOUND(S1,DIM=1)+K-1
    else
       I=LBOUND(S1,DIM=1)
       N=UBOUND(S1,DIM=1)
    endif

    DO   J=I,N
       CALL c_allocda(S1(j))
    ENDDO

  END SUBROUTINE c_allocdas

  SUBROUTINE  c_killda(S1)
!*
    implicit none
    type (c_taylor),INTENT(INOUT)::S1
    !    if(old) then
    call c_DADAL1(s1%i)
    !    else
    !       call KILLNEWDAs(s1%j)
    !    endif

  END SUBROUTINE c_killda

  SUBROUTINE  c_killdaS(S1,k)
!*
    implicit none
    type (c_taylor),INTENT(INOUT),dimension(:)::S1
    INTEGER,optional,INTENT(IN)::k
    INTEGER J,i,N

    if(present(k)) then
       I=LBOUND(S1,DIM=1)
       N=LBOUND(S1,DIM=1)+K-1
    else
       I=LBOUND(S1,DIM=1)
       N=UBOUND(S1,DIM=1)
    endif

    DO   J=I,N
       CALL c_killda(S1(j))
    ENDDO

  END SUBROUTINE c_killdaS

  SUBROUTINE  alloc_c_damap(S1)
!*
    implicit none
    type (c_damap),INTENT(INOUT) :: S1
    INTEGER i,N



     if(s1%n==0) then   
       s1%n=nd2
    endif

    n=s1%n

    call alloc(s1%v,n)

    CALL alloc(s1%s)
    CALL alloc(s1%q)
    s1%q=1.0_dp
    
    !s1%q%x(1)=1.0_dp

    do i=n+1,size(s1%v)
     s1%v(i)%i=0
    enddo
    s1%e_ij=0.0_dp
    s1%x0=0.0_dp
    s1%tpsa=use_tpsa
  END SUBROUTINE alloc_c_damap



  SUBROUTINE  alloc_c_yu_w(S1)
!*
    implicit none
    type (c_yu_w),INTENT(INOUT) :: S1
    INTEGER i,j

     if(s1%n==0) then   
       s1%n=(no-1)/2
    endif
    
    allocate(s1%w(nd2t,0:s1%n))
    
    do i=1,nd2t
    do j=0,s1%n
     call alloc(s1%w(i,j))
    enddo
    enddo
 

  END SUBROUTINE alloc_c_yu_w

  SUBROUTINE  kill_c_yu_w(S1)
!*
    implicit none
    type (c_yu_w),INTENT(INOUT) :: S1
    INTEGER i,j
    

    
    do i=1,nd2t
    do j=0,s1%n
     call kill(s1%w(i,j))
    enddo
    enddo
    s1%n=0
    deallocate(s1%w)

  END SUBROUTINE kill_c_yu_w

  SUBROUTINE  alloc_c_vector_field(S1)
!*
    implicit none
    type (c_vector_field),INTENT(INOUT) :: S1
    INTEGER i,n


    s1%eps   = eps_tpsalie
    s1%nrmax = nrmax
     if(s1%n==0) then   
       s1%n=nd2
    endif

    n=s1%n

    call alloc(s1%v,n)

    call alloc(s1%q)

    do i=n+1,size(s1%v)
     s1%v(i)%i=0
    enddo

  END SUBROUTINE alloc_c_vector_field


  SUBROUTINE  alloc_c_factored_lie(S1,N)
!*
    implicit none
    type (c_factored_lie),INTENT(INOUT) :: S1
    integer,optional :: n
    INTEGER i
 !   logical existed

   ! existed=.false.
   ! if(s1%n/=0) existed=.true.

    if(present(n)) then
     s1%n=n
    else
     s1%n=no
    endif
     
     if(associated(s1%f) ) then
      deallocate(s1%f)
      nullify(s1%f)
     endif

     allocate(s1%f(s1%n))
    do i=1,s1%n
     call alloc(s1%f(i))
    enddo


  END SUBROUTINE alloc_c_factored_lie

  SUBROUTINE  alloc_c_normal_form(S1)
!*
    implicit none
    type (c_normal_form),INTENT(INOUT) :: S1
    
    call alloc(s1%g)
    call alloc(s1%ker)
    call alloc(s1%n)
    call alloc(s1%a2)
    call alloc(s1%a1)
    call alloc(s1%a_t)
    call alloc(s1%AS)
    call alloc(s1%Atot)
    s1%as=1
    s1%nres=0
    s1%positive=.true.
    s1%m=0  ! orbital resonance
    s1%ms=0  ! spin resonance
    s1%g%dir=-1
    s1%ker%dir=1
    s1%s_ij0=0
    s1%s_ijr=0
      s1%tune=0
      s1%damping=0
      s1%spin_tune=0

  END SUBROUTINE alloc_c_normal_form

  SUBROUTINE  kill_c_normal_form(S1)
!*
    implicit none
    type (c_normal_form),INTENT(INOUT) :: S1
    
    call kill(s1%g)
    call kill(s1%ker)
    call kill(s1%n)
    call kill(s1%a_t)
    call kill(s1%a2)
    call kill(s1%a1)
    call kill(s1%AS)
    call kill(s1%Atot)
    s1%positive=.true.
    s1%nres=0
    s1%m=0  ! orbital resonance
    s1%ms=0  ! spin resonance
    s1%g%dir=-1
    s1%ker%dir=1
    s1%s_ij0=0
    s1%s_ijr=0
      s1%tune=0
      s1%damping=0
      s1%spin_tune=0

  END SUBROUTINE kill_c_normal_form

  SUBROUTINE  kill_c_factored_lie(S1)
!*
    implicit none
    type (c_factored_lie),INTENT(INOUT) :: S1
    INTEGER i


    do i=1,s1%n
     call kill(s1%f(i))
    enddo
    s1%n=0
    deallocate(s1%f)

  END SUBROUTINE kill_c_factored_lie

  SUBROUTINE  kill_c_damap(S1)
!*
    implicit none
    type (c_damap),INTENT(INOUT) :: S1

    call kill(s1%v,s1%n)
    CALL kill(s1%s)
    CALL kill(s1%q)
    s1%n=0
    s1%e_ij=0.0_dp
  END SUBROUTINE kill_c_damap

  SUBROUTINE  kill_c_vector_field(S1)
!*
    implicit none
    type (c_vector_field),INTENT(INOUT) :: S1
 


    call kill(s1%v,s1%n)

    call kill(s1%q)

 
    s1%n=0

  END SUBROUTINE kill_c_vector_field

  SUBROUTINE  alloc_c_spinmatrix(S1) ! spin routine
!*
    implicit none
    type (c_spinmatrix),INTENT(INOUT) :: S1
    INTEGER i,J

   do i=1,3
   do j=1,3
       call alloc(s1%s(i,j))
   enddo
   enddo

  END SUBROUTINE alloc_c_spinmatrix


  SUBROUTINE  alloc_c_spinor(S1) ! spin routine
!*
    implicit none
    type (c_spinor),INTENT(INOUT) :: S1
    INTEGER i

   do i=1,3
       call alloc(s1%v(i))
   enddo

  END SUBROUTINE alloc_c_spinor

  SUBROUTINE  kill_c_spinor(S1) ! spin routine
!*
    implicit none
    type (c_spinor),INTENT(INOUT) :: S1
    INTEGER i

   do i=1,3
       call kill(s1%v(i))
   enddo

  END SUBROUTINE kill_c_spinor


  SUBROUTINE  kill_c_spinmatrix(S1) ! spin routine
!*
    implicit none
    type (c_spinmatrix),INTENT(INOUT) :: S1
    INTEGER i,J

   do i=1,3
   do j=1,3
       call kill(s1%s(i,j))
   enddo
   enddo

  END SUBROUTINE kill_c_spinmatrix




  function  c_real(S1)
    implicit none
    type (c_taylor) c_real
    type (c_taylor),INTENT(IN)::S1
    integer i,n,localmaster
    complex(dp) x
    complex(dp) value
    integer, allocatable :: j(:)


    IF(.NOT.C_STABLE_DA) then
     c_real%i=0
     RETURN
    endif

    localmaster=c_master

    !    call check(s1)
    call ass(c_real)

    !    if(old) then
 
    if(s1%i==0) call c_crap1("EQUAL 2") ! call allocw(s1)

    allocate(j(c_%nv))
    c_real=0.0_dp
   
    call c_taylor_cycle(s1,size=n)

    do i=1,n
       call c_taylor_cycle(s1,ii=i,value=value,j=j)
       x=real(value)
       c_real=c_real+(x.cmono.j)
   enddo
   
     c_master=localmaster

   deallocate(j)
  END  function c_real

 function  c_aimag(S1)
    implicit none
    type (c_taylor) c_aimag
    type (c_taylor),INTENT(IN)::S1
    integer i,n,localmaster
    complex(dp) x
    complex(dp) value
    integer, allocatable :: j(:)

    IF(.NOT.C_STABLE_DA) then
     c_aimag%i=0
     RETURN
    endif

    localmaster=c_master

    !    call check(s1)
    call ass(c_aimag)

    !    if(old) then
 
    if(s1%i==0) call c_crap1("EQUAL 2") ! call allocw(s1)

    allocate(j(c_%nv))
    c_aimag=0.0_dp
   
    call c_taylor_cycle(s1,size=n)

    do i=1,n
       call c_taylor_cycle(s1,ii=i,value=value,j=j)
       x=aimag(value)
       c_aimag=c_aimag+(x.cmono.j)
   enddo
   
     c_master=localmaster

   deallocate(j)
  END  function c_aimag


  SUBROUTINE  equalc_t(S2,S1)
!*
    implicit none
    type (c_taylor),INTENT(inOUT)::S2
    type (taylor),INTENT(IN)::S1
    integer i,n
    real(dp) x
    integer, allocatable :: j(:)
    IF(.NOT.C_STABLE_DA) RETURN

!    call c_check_snake
    call check_snake
    !    if(old) then
    if(s2%i==0) then
       call c_crap1("EQUAL 1 in tpsa") !call allocw(s2)
    endif
    if(s1%i==0) call crap1("EQUAL 2") ! call allocw(s1)

    allocate(j(c_%nv))
    s2=0.0_dp
   
    call taylor_cycle(s1,size=n)

    do i=1,n
       call taylor_cycle(s1,ii=i,value=x,j=j)
       s2=s2+(x.cmono.j)
   enddo
   

   deallocate(j)
  END SUBROUTINE equalc_t

 SUBROUTINE  equalt_c(S2,S1)
!*
    implicit none
    type (taylor),INTENT(inOUT)::S2
    type (c_taylor),INTENT(IN)::S1
    integer i,n
    complex(dp) value
    integer, allocatable :: j(:)
    IF(.NOT.C_STABLE_DA) RETURN

    call c_check_snake
 !   call check_snake
    !    if(old) then
    if(s2%i==0) then
       call crap1("EQUAL 1 in tpsa") !call allocw(s2)
    endif
    if(s1%i==0) call c_crap1("EQUAL 2") ! call allocw(s1)

    allocate(j(nv))
    s2=0.0_dp
   
    call c_taylor_cycle(s1,size=n)

    do i=1,n
       call c_taylor_cycle(s1,ii=i,value=value,j=j)
       s2=s2+(value.mono.j)
   enddo
   
   deallocate(j)
  END SUBROUTINE equalt_c

 SUBROUTINE  equalc_ray_ray(S2,S1)
!*
    implicit none
    type (c_ray),INTENT(inOUT)::S2
     type (c_ray),INTENT(IN)::S1 
     
     s2%x  = s1%x
     s2%s1 = s1%s1
     s2%s2 = s1%s2
     s2%s3 = s1%s3
     s2%q = s1%q
     s2%n = s1%n
     s2%x0  = s1%x0
  END SUBROUTINE equalc_ray_ray

 SUBROUTINE  equalc_ray_r6(S2,S1)
!*
    implicit none
    type (c_ray),INTENT(inOUT)::S2
    complex(dp),INTENT(IN)::S1(:)
    integer i


    IF(.NOT.C_STABLE_DA) RETURN
    
     s2%x=0.0_dp
     s2%x0=0.0_dp
     do i=1,size(s1)
       s2%x(i)=s1(i)
     enddo
   !  s2%q=1.0_dp
     s2%n=size(s1)
      s2%s1=0
      s2%s2=0
      s2%s3=0
           
      s2%s1(1)=1
      s2%s2(2)=1
      s2%s3(3)=1
     
  END SUBROUTINE equalc_ray_r6

 SUBROUTINE  equalc_ray_r6r(S2,S1)
!*
    implicit none
    type (c_ray),INTENT(inOUT)::S2
    real(dp),INTENT(IN)::S1(:)
    integer i


    IF(.NOT.C_STABLE_DA) RETURN
    
     s2%x=0.0_dp
     do i=1,size(s1)
       s2%x(i)=s1(i)
     enddo
  !   s2%q=1.0_dp
      s2%s1=0
      s2%s2=0
      s2%s3=0
     s2%n=size(s1)           
      s2%s1(1)=1
      s2%s2(2)=1
      s2%s3(3)=1
     
  END SUBROUTINE equalc_ray_r6r

 SUBROUTINE  equalc_r6_ray(S1,S2)
!*
    implicit none
    type (c_ray),INTENT(in)::S2
    complex(dp),INTENT(inOUT)::S1(:)
    integer i


    IF(.NOT.C_STABLE_DA) RETURN
    

     do i=1,size(s1)
       s1(i)=s2%x(i)
     enddo

     
  END SUBROUTINE equalc_r6_ray


 SUBROUTINE  equalc_t_ct(S2,S1)
!*
    implicit none
    type (c_taylor),INTENT(inOUT)::S2
    type (complextaylor),INTENT(IN)::S1
    integer i,n
    complex(dp) x
    real(dp) value
    integer, allocatable :: j(:)
    IF(.NOT.C_STABLE_DA) RETURN

!    call c_check_snake
    call check_snake
    !    if(old) then
    if(s2%i==0) then
       call c_crap1("EQUAL 1 in tpsa") !call allocw(s2)
    endif
    if(s1%r%i==0) call crap1("EQUAL 2") ! call allocw(s1)

    allocate(j(c_%nv))
    s2=0.0_dp
   
    call taylor_cycle(s1%r,size=n)

    do i=1,n
       call taylor_cycle(s1%r,ii=i,value=value,j=j)
       x=value
       s2=s2+(x.cmono.j)
   enddo
   
    call taylor_cycle(s1%i,size=n)

    do i=1,n
       call taylor_cycle(s1%i,ii=i,value=value,j=j)
       x=value*i_
       s2=s2+(x.cmono.j)
   enddo


   deallocate(j)
  END SUBROUTINE equalc_t_ct

  SUBROUTINE  equalc_ct_c(S2,S1)
!*
    implicit none
    type (complextaylor),INTENT(inOUT)::S2
    type (c_taylor),INTENT(IN)::S1
    integer i,n
    complex(dp) value
    integer, allocatable :: j(:)
    IF(.NOT.C_STABLE_DA) RETURN

    call c_check_snake
 !   call check_snake
    !    if(old) then
    if(s2%r%i==0) then
       call c_crap1("EQUAL 1 in tpsa") !call allocw(s2)
    endif
    if(s1%i==0) call crap1("EQUAL 2") ! call allocw(s1)

    allocate(j(nv))
    s2=0.0_dp
   
    call c_taylor_cycle(s1,size=n)

    do i=1,n
       call c_taylor_cycle(s1,ii=i,value=value,j=j)
       s2=s2+(value.mono.j)
   enddo
   
   deallocate(j)
  END SUBROUTINE equalc_ct_c

  SUBROUTINE  equalc_cmap_map(S2,S1)
!*
    implicit none
    type (c_damap),INTENT(inOUT)::S2
    type (damap),INTENT(IN)::S1
    type(complextaylor) ct
    integer i,i1(4),i2(4)



    call check_snake
    call alloc(ct)
    call liepeek(i1,i2)
     if(i1(4)/=s2%n) then
       write(6,*) "Error in equalc_cmap_map",i1(4),s2%n
     endif
    do i=1,i1(4)
     ct=s1%v(i)
!     ct=ct-(ct.sub.'0')
      s2%v(i)=ct
    enddo

    call kill(ct)

 end SUBROUTINE  equalc_cmap_map

  SUBROUTINE  equalc_map_cmap(S2,S1)
!*
    implicit none
    type (damap),INTENT(inOUT)::S2
    type (c_damap),INTENT(IN)::S1
    type(complextaylor) ct
    integer i,i1(4),i2(4)



    call c_check_snake

    call alloc(ct)
    call liepeek(i1,i2)
     if(i1(4)/=s1%n) then
       write(6,*) "Error in equalc_map_cmap",i1(4),s1%n
     endif
    do i=1,i1(4)
     ct=s1%v(i)
!     ct=ct-(ct.sub.'0')
      s2%v(i)=ct%r
    enddo
 
    call kill(ct)

 end SUBROUTINE  equalc_map_cmap

  SUBROUTINE  equal_real8_cmap(S2,S1)
!*
    implicit none
    type (real_8),INTENT(inOUT)::S2(:)
    type (c_damap),INTENT(IN)::S1
    type(taylor) ct
    integer i


    call c_check_snake

    call alloc(ct)

    do i=1,nd2
     ct=s1%v(i)

      s2(i)=ct 
    enddo
 
    call kill(ct)

 end SUBROUTINE  equal_real8_cmap

  SUBROUTINE  equal_cmap_real8(S1,S2)
!*
    implicit none
    type (real_8),INTENT(in)::S2(:)
    type (c_damap),INTENT(inOUT)::S1
    type(taylor) ct
    integer i


    call check_snake

    call alloc(ct)

    do i=1,nd2
     ct=s2(i)

      s1%v(i)=ct 

    enddo
 
    call kill(ct)

 end SUBROUTINE  equal_cmap_real8

  subroutine EQUAL_c_spinmatrix_probe(S,R)
!*
    implicit none
    TYPE(probe), INTENT(IN) :: R
    TYPE(c_spinmatrix), INTENT(INOUT) :: S

    s%s(1,1)=r%s(1)%x(1);    s%s(2,1)=r%s(1)%x(2);    s%s(3,1)=r%s(1)%x(3);
    s%s(1,2)=r%s(2)%x(1);    s%s(2,2)=r%s(2)%x(2);    s%s(3,2)=r%s(2)%x(3);
    s%s(1,3)=r%s(3)%x(1);    s%s(2,3)=r%s(3)%x(2);    s%s(3,3)=r%s(3)%x(3);

  END subroutine EQUAL_c_spinmatrix_probe

  subroutine EQUAL_probe_c_spinmatrix(R,S)
!*
    implicit none
    TYPE(probe), INTENT(INout) :: R
    TYPE(c_spinmatrix), INTENT(IN) :: S

    r%s(1)%x(1)=s%s(1,1);    r%s(1)%x(2)=s%s(2,1);    r%s(1)%x(3)=s%s(3,1);
    r%s(2)%x(1)=s%s(1,2);    r%s(2)%x(2)=s%s(2,2);    r%s(2)%x(3)=s%s(3,2);
    r%s(3)%x(1)=s%s(1,3);    r%s(3)%x(2)=s%s(2,3);    r%s(3)%x(3)=s%s(3,3);

  END subroutine EQUAL_probe_c_spinmatrix

  subroutine EQUAL_c_spinmatrix_3_by_3(S,R)
!*
    implicit none
    real(dp), INTENT(IN) :: R(3,3)
    TYPE(c_spinmatrix), INTENT(INOUT) :: S
    integer i,j
    do i=1,3
    do j=1,3
    s%s(i,j)=r(i,j) 
    enddo
    enddo

  END subroutine EQUAL_c_spinmatrix_3_by_3


  subroutine EQUAL_3_by_3_c_spinmatrix(R,S)
!*
    implicit none
    real(dp), INTENT(INOUT) :: R(3,3)
    TYPE(c_spinmatrix), INTENT(IN) :: S
    integer i,j
    do i=1,3
    do j=1,3
    r(i,j)=s%s(i,j)
    enddo
    enddo

  END subroutine EQUAL_3_by_3_c_spinmatrix


  subroutine EQUAL_3_by_3_probe(R,S)
!*
    implicit none
    real(dp), INTENT(INout) :: R(3,3)
    TYPE(probe), INTENT(IN) :: S
    integer i,j

    r(1,1)=s%s(1)%x(1);    r(2,1)=s%s(1)%x(2);    r(3,1)=s%s(1)%x(3);
    r(1,2)=s%s(2)%x(1);    r(2,2)=s%s(2)%x(2);    r(3,2)=s%s(2)%x(3);
    r(1,3)=s%s(3)%x(1);    r(2,3)=s%s(3)%x(2);    r(3,3)=s%s(3)%x(3);

  END subroutine EQUAL_3_by_3_probe

  subroutine EQUAL_probe_3_by_3(S,R)
!*
    implicit none
    real(dp), INTENT(IN) :: R(3,3)
    TYPE(probe), INTENT(INout) :: S
    integer i,j

    s%s(1)%x(1)=r(1,1);    s%s(1)%x(2)=r(2,1);    s%s(1)%x(3)=r(3,1);
    s%s(2)%x(1)=r(1,2);    s%s(2)%x(2)=r(2,2);    s%s(2)%x(3)=r(3,2);
    s%s(3)%x(1)=r(1,3);    s%s(3)%x(2)=r(2,3);    s%s(3)%x(3)=r(3,3);

  END subroutine EQUAL_probe_3_by_3

  subroutine EQUAL_c_map_RAY8(DS,R)
!*
    implicit none
    TYPE(probe_8), INTENT(IN) :: R
    TYPE(c_damap), INTENT(INOUT) :: DS
    real(dp) m(6,6)
    type(taylor) t

    INTEGER I,J,nd2t1
    logical(lp) rad_in

    call alloc(t)

     nd2t1=c_%nd2-2*rf 
 !   nd2t1=C_%ND2
 !   if(doing_ac_modulation_in_ptc) then
 !      nd2t1=C_%ND2-2
 !   endif



    DO I=1,nd2t1
       t=R%X(I)
       DS%V(I)=t
    ENDDO

 !   DO I=nd2t1+1,C_%ND2
 !      t=R%ac%x(i-nd2t1)
 !      DS%V(I)=t
 !   ENDDO

    j=1
    DO I=nd2t1+1,C_%ND2,2
       t=R%ac(j)%x(1)
       DS%V(I)=t
       t=R%ac(j)%x(2)
       DS%V(I+1)=t
       j=j+1
    ENDDO

! quaternion
    if(use_quaternion)   THEN
    DO I=0,3

          t=r%q%x(i)
          DS%q%x(i)=t

    ENDDO
else
    DO I=1,3
       DO J=1,3
          t=R%S(J)%X(I)
          DS%S%s(I,J)=t
       ENDDO
    ENDDO
endif

    call check_rad(r%e_ij,rad_in)
    ds%e_ij=0.0_dp
    if(rad_in) then
       m=ds
       ds%e_ij=matmul(matmul(m,r%e_ij),transpose(m))
    endif
DS%x0=0
ds%x0(1:6)=r%x0
    call kill(t)

  END subroutine EQUAL_c_map_RAY8

  subroutine EQUAL_RAY8_c_map(R,DS)
!*
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: R
    TYPE(c_damap), INTENT(IN) :: DS
    real(dp) m(6,6)
    logical(lp) rad_in
    INTEGER I,J,nd2t1
    type(taylor) t

    call alloc(t)

    nd2t1=c_%nd2-2*rf

    DO I=1,nd2t1
       t=DS%V(I)
       R%X(I)= t  
    ENDDO
 
    
 !     do i=nd2t1+1,c_%nd2
 !      t=DS%V(I)
 !      r%ac%x(i-nd2t1) =  t
 !     enddo

    j=1
    DO I=nd2t1+1,C_%ND2,2
       t=DS%V(I)
       R%ac(j)%x(1)=t
       t=DS%V(I+1) 
       R%ac(j)%x(2)=t
       j=j+1
    ENDDO

! quaternion
if(use_quaternion)   THEN
    DO I=0,3
          t=DS%q%x(i)
          r%q%x(i)=t
    ENDDO
else
    DO J=1,3
       DO I=1,3
          t=DS%S%s(I,J) 
          R%S(J)%X(I)=t     
       ENDDO
    ENDDO
endif
    call c_check_rad(ds%e_ij,rad_in)
        r%e_ij=0.0_dp
    if(rad_in) then
       m=ds**(-1)
       r%e_ij=matmul(matmul(m,ds%e_ij),transpose(m))
    endif
 
   call kill(t)

  END subroutine EQUAL_RAY8_c_map

 

  SUBROUTINE  equalc_cvec_vec(S2,S1)
!*
    implicit none
    type (c_vector_field),INTENT(inOUT)::S2
    type (vecfield),INTENT(IN)::S1
    type(complextaylor) ct
    integer i,i1(4),i2(4)



    call check_snake

    call alloc(ct)

    call liepeek(i1,i2)
     if(i1(4)/=s2%n) then
       write(6,*) "Error in equalc_cmap_map",i1(4),s2%n
     endif
    
    do i=1,i1(4)
     ct=s1%v(i)
      s2%v(i)=ct
    enddo
 
    call kill(ct)

 end SUBROUTINE  equalc_cvec_vec

  SUBROUTINE  equalc_cspinor_cspinor(S2,S1) ! spin routine
!*
    implicit none
    type (c_spinor),INTENT(inOUT)::S2
    type (c_spinor),INTENT(IN)::S1

    integer i 

    call check_snake

    do i=1,3
      s2%v(i)=s1%v(i)
    enddo


 end SUBROUTINE  equalc_cspinor_cspinor

  SUBROUTINE  equalc_spinor_c_quaternion(S2,S1) ! spin routine
!*
    implicit none
    type (c_spinor),INTENT(inOUT)::S2
    type (c_quaternion),INTENT(IN)::S1

    integer i 

    call check_snake

    do i=1,3
      s2%v(i)=s1%x(i)
    enddo


 end SUBROUTINE  equalc_spinor_c_quaternion

  SUBROUTINE  equalc_quaternion_c_spinor(S2,S1) ! spin routine
!*
    implicit none
    type (c_quaternion),INTENT(inOUT)::S2
    type (c_spinor),INTENT(IN)::S1

    integer i 

    call check_snake
    s2%x(1)=0.0_dp
    do i=1,3
      s2%x(i)=s1%v(i)
    enddo


 end SUBROUTINE  equalc_quaternion_c_spinor

  SUBROUTINE  equalc_spinor_cspinor(S2,S1) ! spin routine
!*
    implicit none
    type (spinor),INTENT(inOUT)::S2
    type (c_spinor),INTENT(IN)::S1

    integer i 

    call check_snake

    do i=1,3
      s2%x(i)=s1%v(i)
    enddo


 end SUBROUTINE  equalc_spinor_cspinor

  SUBROUTINE  equalc_cspinor_spinor(S1,S2) ! spin routine
!*
    implicit none
    type (spinor),INTENT(in)::S2
    type (c_spinor),INTENT(inOUT)::S1

    integer i 

    call check_snake

    do i=1,3
      s1%v(i)=s2%x(i)
    enddo


 end SUBROUTINE  equalc_cspinor_spinor

  SUBROUTINE  c_DPEKMAP(S2,S1)
    implicit none
    complex(dp),INTENT(inOUT),dimension(:)::S2
    type (c_damap),INTENT(IN)::S1
    IF(.NOT.C_STABLE_DA) RETURN
    call c_check_snake
    ! if(old) then
    CALL c_DAPEK0(S1%V%I,S2,s1%n)
    !    else
    !       CALL newDAPEK0(S1%V%J,S2,nd2)
    !    endif
  END SUBROUTINE c_DPEKMAP

  SUBROUTINE  c_DPOKMAP(S1,S2)
    implicit none
    complex(dp),INTENT(IN),dimension(:)::S2
    type (c_damap),INTENT(inOUT)::S1
    IF(.NOT.C_STABLE_DA) RETURN
 
    ! if(old) then
    if(s1%V(1)%i==0) call crap1("DPOKMAP 1") !call allocw_old(s1%V(1))   !call etall(s1%V%i,ND2)
    CALL c_DAPOK0(S1%V%I,S2,s1%n)
    !    else
    !       if(.NOT. ASSOCIATED(s1%V(1)%j%r))  call crap1("DPOKMAP 2") !  !call allocw_old(s1%V(1))  !call newetall(s1%V%j,ND2)
    !       CALL NEWDAPOK0(S1%V%J,S2,nd2)
    !    endif
  END SUBROUTINE c_DPOKMAP

  SUBROUTINE  EQUAL(S2,S1)
!*
    implicit none
    type (c_taylor),INTENT(inOUT)::S2
    type (c_taylor),INTENT(IN)::S1
    IF(.NOT.C_STABLE_DA) RETURN

    call c_check_snake
    !    if(old) then
    if(s2%i==0) then
       call c_crap1("EQUAL 1 in tpsa") !call allocw(s2)
    endif
    if(s1%i==0) call c_crap1("EQUAL 2") ! call allocw(s1)
    CALL c_dacop(S1%I,S2%I)
 
  END SUBROUTINE EQUAL

  SUBROUTINE  equal_map_real8(S2,S1)
!*
    implicit none
    type (c_damap),INTENT(inOUT)::S2
    real(dp),INTENT(IN),dimension(:)::S1     

    integer i
    IF(.NOT.C_STABLE_DA) RETURN

    call c_check_snake
    !    if(old) then
    if(s2%v(1)%i==0) then
       call c_crap1("EQUAL 1 in tpsa") !call allocw(s2)
    endif

    do i=1, min(s2%n,size(s1,1)) 
     S2%v(i)=s1(i)
    enddo
 
  END SUBROUTINE equal_map_real8

  SUBROUTINE  equal_map_complex8(S2,S1)
!*
    implicit none
    type (c_damap),INTENT(inOUT)::S2
    complex(dp),INTENT(IN),dimension(:)::S1   

    integer i
    IF(.NOT.C_STABLE_DA) RETURN

    call c_check_snake
    !    if(old) then
    if(s2%v(1)%i==0) then
       call c_crap1("EQUAL 1 in tpsa") !call allocw(s2)
    endif
 
    do i=1, min(s2%n,size(s1,1)) 
     S2%v(i)=s1(i)
    enddo
 
  END SUBROUTINE equal_map_complex8

  SUBROUTINE  equal_real8_map(S2,S1)
!*
    implicit none
    real(dp),INTENT(inOUT),dimension(:)::S2
    type (c_damap),INTENT(IN)::S1 
    integer i
    IF(.NOT.C_STABLE_DA) RETURN

    call c_check_snake
    !    if(old) then
    if(s1%v(1)%i==0) then
       call c_crap1("EQUAL 1 in tpsa") !call allocw(s2)
    endif
 
    do i=1, min((s1%n),size(s2,1)) 
     s2(i)=S1%v(i)
    enddo
 
  END SUBROUTINE equal_real8_map

  SUBROUTINE  equal_complex8_map(S2,S1)
!*
    implicit none
    complex(dp),INTENT(inOUT),dimension(:)::S2
    type (c_damap),INTENT(IN)::S1
    integer i
    IF(.NOT.C_STABLE_DA) RETURN

    call c_check_snake
    !    if(old) then
    if(s1%v(1)%i==0) then
       call c_crap1("EQUAL 1 in tpsa") !call allocw(s2)
    endif
    
    do i=1, min(s1%n,size(s2,1)) 
     s2(i)=S1%v(i)
    enddo
 
  END SUBROUTINE equal_complex8_map

!skowron to bypass strange gfortran error when using s2=s1     
 SUBROUTINE  equal_c_tayls(S2,S1)
    implicit none
    type (c_taylor),INTENT(inOUT)::S2
    type (c_taylor),INTENT(IN)::S1

    call equal(s2,s1)

  end SUBROUTINE  equal_c_tayls 

  SUBROUTINE  EQUALspinmatrix(S2,S1) ! spin routine
!*
    implicit none
    type (c_spinmatrix),INTENT(inOUT)::S2
    type (c_spinmatrix),INTENT(IN)::S1
    integer i,j
    IF(.NOT.C_STABLE_DA) RETURN

     call c_check_snake
 
 
     do i=1,3
     do j=1,3
      s2%s(i,j)=s1%s(i,j)
     enddo
     enddo

  END SUBROUTINE EQUALspinmatrix

  SUBROUTINE  cDEQUAL(R1,S2)
!*
    implicit none
    type (c_taylor),INTENT(IN)::S2
    complex(dp), INTENT(inOUT)::R1
    IF(.NOT.C_STABLE_DA) RETURN
    call c_check_snake

    R1=S2.SUB.'0'
  END SUBROUTINE cDEQUAL

  SUBROUTINE  DEQUAL(R1,S2)
!*
    implicit none
    type (c_taylor),INTENT(IN)::S2
    real(dp), INTENT(inOUT)::R1
    IF(.NOT.C_STABLE_DA) RETURN
    call c_check_snake

    R1=S2.SUB.'0'
  END SUBROUTINE DEQUAL

  SUBROUTINE  REQUAL(R1,S2)
!*
    implicit none
    type (c_taylor),INTENT(IN)::S2
    REAL(SP), INTENT(inOUT)::R1
    IF(.NOT.C_STABLE_DA) RETURN

    if(c_real_warning) call c_real_stop
    call c_check_snake

    R1=S2.SUB.'0'

  END SUBROUTINE REQUAL

  function  DAABSEQUAL(S2)
    implicit none
    type (c_taylor),INTENT(IN)::S2
    real(dp) DAABSEQUAL
    IF(.NOT.C_STABLE_DA) RETURN

    call c_check_snake
    DAABSEQUAL=abs(S2.sub.'0')

  END function DAABSEQUAL


  SUBROUTINE  cDEQUALDACON(S2,R1)
!*
    implicit none
    type (c_taylor),INTENT(inOUT)::S2
    complex(dp), INTENT(IN)::R1
    IF(.NOT.C_STABLE_DA) RETURN

    !    if(old) then
    if(s2%i==0)  call c_crap1("DEQUALDACON 1") !call allocw(s2)
    CALL c_dacon(S2%I,R1)
    !    else
    !       IF (.NOT. ASSOCIATED(s2%j%r)) call c_crap1("DEQUALDACON 2") !call allocw(s2)
    !       CALL newDACON(S2%j,R1)
    !    endif
  END SUBROUTINE cDEQUALDACON

  SUBROUTINE  DEQUALDACON(S2,R1)
!*
    implicit none
    type (c_taylor),INTENT(inOUT)::S2
    real(dp), INTENT(IN)::R1
    complex(dp) rt
    IF(.NOT.C_STABLE_DA) RETURN

    !    if(old) then
    if(s2%i==0)  call c_crap1("DEQUALDACON 1") !call allocw(s2)
    rt=r1
    CALL c_dacon(S2%I,rt)
    !    else
    !       IF (.NOT. ASSOCIATED(s2%j%r)) call c_crap1("DEQUALDACON 2") !call allocw(s2)
    !       CALL newDACON(S2%j,R1)
    !    endif
  END SUBROUTINE DEQUALDACON

  SUBROUTINE  DEQUALDACONS(S2,R1)
!*
    implicit none
    type (c_taylor),INTENT(inOUT)::S2(:)
    real(dp), INTENT(IN)::R1
    complex(dp) rt
    integer k
    IF(.NOT.C_STABLE_DA) RETURN

    !    if(old) then
    if(s2(1)%i==0)  call c_crap1("DEQUALDACON 1") !call allocw(s2)
    do k=lbound(s2,1),ubound(s2,1)
    rt=r1
    CALL c_dacon(S2(k)%I,rt)
    enddo
    !    else
    !       IF (.NOT. ASSOCIATED(s2%j%r)) call c_crap1("DEQUALDACON 2") !call allocw(s2)
    !       CALL newDACON(S2%j,R1)
    !    endif
  END SUBROUTINE DEQUALDACONS


  SUBROUTINE  EQUALDACON(S2,R1)
!*
    implicit none
    type (c_taylor),INTENT(inOUT)::S2
    REAL(SP), INTENT(IN)::R1
    real(dp) R2
    IF(.NOT.C_STABLE_DA) RETURN
    if(c_real_warning) call c_real_stop
  !  call c_check_snake

    if(c_real_warning) call c_real_stop
    !    if(old) then
    if(s2%i==0) call c_crap1("EQUALDACON 1") !call allocw(s2)
    !    else
    !       IF (.NOT. ASSOCIATED(s2%j%r)) call c_crap1("EQUALDACON 2") !call allocw(s2)
    !    endif
    r2=REAL(r1,kind=DP)
    s2=r2
  END SUBROUTINE EQUALDACON

  SUBROUTINE  IEQUALDACON(S2,R1)
!*
    implicit none
    type (c_taylor),INTENT(inOUT)::S2
    INTEGER, INTENT(IN)::R1
    real(dp) r2
    IF(.NOT.C_STABLE_DA) RETURN
 !   call c_check_snake


    ! if(old) then
    if(s2%i==0) call c_crap1("IEQUALDACON 1") !call allocw(s2)
    !    else
    !       IF (.NOT. ASSOCIATED(s2%j%r)) call c_crap1("IEQUALDACON 2") !call allocw(s2)
    !    endif
    r2=REAL(r1,kind=DP)
    s2=r2
  END SUBROUTINE IEQUALDACON

  FUNCTION dexpt( S1 )
    implicit none
    TYPE (c_taylor) dexpt
    TYPE (c_taylor), INTENT (IN) :: S1
    integer localmaster

    IF(.NOT.C_STABLE_DA) then
     dexpt%i=0
     RETURN
    endif
    localmaster=c_master

    !    call check(s1)
    call ass(dexpt)

    ! if(old) then
    call c_dafun('EXP ',s1%i,c_temp%i)
    call c_dacop(c_temp%i,dexpt%i)
    !    else
    !       call newdafun('EXP ',s1%j,dexpt%j)
    !    endif

    c_master=localmaster

  END FUNCTION dexpt



  FUNCTION c_logt( S1 )
    implicit none
    TYPE (c_taylor) c_logt
    TYPE (c_taylor), INTENT (IN) :: S1
    integer localmaster

    IF(.NOT.C_STABLE_DA) then
     c_logt%i=0
     RETURN
    endif
    localmaster=c_master

    !    call check(s1)
    call ass(c_logt)

    ! if(old) then
    call c_dafun('LOG ',s1%i,c_temp%i)
    call c_dacop(c_temp%i,c_logt%i)
    !    else
    !       call newdafun('EXP ',s1%j,dexpt%j)
    !    endif

    c_master=localmaster

  END FUNCTION c_logt

  SUBROUTINE  flatten_c_factored_lie_r(St,S1)
  implicit none
  type (c_factored_lie),INTENT(IN) :: S1
  type (c_vector_field),INTENT(inOUT) :: ST
  call flatten_c_factored_lie(S1,ST)
end  SUBROUTINE  flatten_c_factored_lie_r

  SUBROUTINE  flatten_c_factored_lie(S1,ST)
!#general: manipulation
!# Type c_factored_lie is a product of Lie exponents.
!# Lie map=exp(s1%f(1).grad)...exp(s1%f(s1%n).grad) if s1%dir=1.
!# or Lie map=exp(s1%f(s1%n).grad)...exp(s1%f(1).grad) if s1%dir=-1.
!# In some cases, these exponents commute and can be added; for example
!# if the map is a rotation.
!# s2 = s1%f(1)+...+s1%f(s1%n)
  implicit none
  type (c_factored_lie),INTENT(IN) :: S1
  type (c_vector_field),INTENT(inOUT) :: ST
  integer i,j
  type (c_vector_field) S2
!  This routine assumes commutation of all the factored Lie exponents
    
    s2%n = s1%f(1)%n
    call alloc(s2)

    s2%eps   = s1%f(1)%eps
    s2%nrmax = s1%f(1)%nrmax
    
 
    s2=0
    do j=1,s1%n
if(use_quaternion)   THEN
       s2%q=s2%q+s1%f(j)%q
 
endif


       do i=1,s2%n
       s2%v(i)=s2%v(i)+s1%f(j)%v(i)
      enddo      
    enddo
    st=s2
  if(complex_extra_order==1.and.special_extra_order_1) s2=s2.cut.no
   call kill(s2)
  END SUBROUTINE flatten_c_factored_lie
! etienne1
  FUNCTION c_logf_spin( s1,h,epso,n,tpsa )
    implicit none
    TYPE (c_vector_field) c_logf_spin
    TYPE (c_vector_field) , optional :: h
    TYPE (c_damap), INTENT (INOUT) :: S1
     real(dp), optional :: epso
     integer, optional :: n
     logical, optional :: tpsa
     logical da
    integer localmaster,k,i
    TYPE (c_vector_field) t,t2
    type(c_damap) mt,t_1,e2,s1t
    real(dp) xnorm1,epsone,xn,xnorm2,r
    IF(.NOT.C_STABLE_DA) then
     c_logf_spin%v%i=0
     RETURN
    endif
    da=.not.s1%tpsa
    
    if(present(tpsa)) da=.not.tpsa
    localmaster=c_master
    t%n=s1%n;t2%n=s1%n;mt%n=s1%n;t_1%n=s1%n;e2%n=s1%n;s1t%n=s1%n;
    call alloc(t);call alloc(t2);call alloc(mt,t_1,e2,s1t);

    c_logf_spin%n=s1%n
    call c_ass_vector_field(c_logf_spin)
    c_logf_spin=0
     if(present(h)) c_logf_spin=h
    s1t=s1
    xnorm1=0.0_dp
if(use_quaternion) then
   call c_full_norm_quaternion(s1t%q,k,xnorm1)
else
  write(6,*) "log no longer available for SO(3) "
 stop 1959
endif
    do i=1,s1t%n
       if(da) then
         s1t%v(i)=s1t%v(i)-(s1t%v(i).sub.'0')
        endif
       r=full_abs(s1t%v(i))
       xnorm1=xnorm1+r
    enddo

    if(present(epso)) then
     epsone= epso
    else
     epsone= xnorm1/epso_factor
    endif
       if(lielib_print(3)==1) write(6,*) epsone,xnorm1
!     eps=abs(epsone)/2  !1.d5
    xn=1e36_dp

    k=1000
  if(present(n)) k=n
    do i=1,k

      mt=exp(-c_logf_spin,s1t)

      t=c_1_vf_q(mt)   !  extracts mt-1 as a vector field including spin
 
      t_1=1
      t_1=mt-t_1     !  same as a map
      e2=(t*t_1)
 
      t2=(-0.5d0)*c_1_vf_q(e2,0)
      t=t+t2
 
      call c_full_norm_vector_field(t,xnorm1)
      if(lielib_print(3)==1) write(6,*) i,xn,xnorm1

 
    if(xnorm1.lt.epsone) then !.and.i>=10) then
     call c_full_norm_vector_field(t2,xnorm2)
 
      c_logf_spin=c_logf_spin+t+0.5e0_dp*(c_logf_spin.lb.t)
      if(xnorm1>=xn.and.(.not.present(n))) exit
      xn=xnorm1
    else
      c_logf_spin=c_logf_spin+t
      xn=xnorm1
    endif

   enddo
   if(i>k-10.and.(.not.present(n))) then
     write(6,*) " no convergence in c_logf_spin "
   !  stop 1984
   endif
    call kill(t);call kill(t2);call kill(mt,t_1,e2,s1t);
!     call c_flofacg(s1,c_logf_spin,epso,n)
  if(complex_extra_order==1.and.special_extra_order_1) c_logf_spin=c_logf_spin.cut.no
    c_master=localmaster

  END FUNCTION c_logf_spin

! etienne2
  FUNCTION c_logf( s1,h,epso,n,tpsa )
!#internal: manipulation
!# Accessable with the interface log if desired.
!# Takes the logarithm of a nonlinear map s1
!# s1=exp(log(s1).grad) I
!# The map must be near the identity.
!# H, epso and n are optional.
!# epso is a small positive number estimated automatically.
!# The algorithm goes nonlinear when the norm of the correction
!# is smaller than epso. See Chap.11.
!# H is a guess for the logarithm: optional.

    implicit none
    TYPE (c_vector_field) c_logf
    TYPE (c_vector_field) , optional :: h
    TYPE (c_damap), INTENT (INOUT) :: S1
     real(dp), optional :: epso
     integer, optional :: n
     logical, optional :: tpsa
    type(c_damap) s1t
     logical da
    integer localmaster,i
    IF(.NOT.C_STABLE_DA) then
     c_logf%v%i=0
     RETURN
    endif
    localmaster=c_master
     s1t%n=s1%n
     call alloc(s1t)
    s1t=s1
    da=.true.
    if(present(tpsa)) da=.not.tpsa
    if(da) then
     do i=1,s1t%n
         s1t%v(i)=s1t%v(i)-(s1t%v(i).sub.'0')
     enddo
    endif

    c_logf%n=s1%n
    call c_ass_vector_field(c_logf)
    c_logf=0
     if(present(h)) c_logf=h
     call c_flofacg(s1t,c_logf,epso,n)
    call kill(s1t)
    c_master=localmaster
  if(complex_extra_order==1.and.special_extra_order_1) c_logf=c_logf.cut.no
  END FUNCTION c_logf

  FUNCTION c_logc( S1 )
    implicit none
    complex(dp) c_logc
    complex(dp), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     c_logc=0
     RETURN
    endif
    localmaster=c_master

    c_logc  = LOG(abs(S1)) + (0.0_dp,1.0_dp)*atan2(aimag(S1),real(S1))


  END FUNCTION c_logc

  FUNCTION FULL_ABST( S1 )
!#general: manipulation
!# This routine computes the norm of the c_taylor s1.
!# This routine is called with the interface
!# full_abs(s1) where s1 is a c_taylor. 
    implicit none
    real(dp) FULL_ABST
    TYPE (c_taylor), INTENT (IN) :: S1

    IF(.NOT.C_STABLE_DA) then
     full_abst=0
     RETURN
    endif
    !    call check(s1)

    ! if(old) then
    CALL c_DAABS(S1%I,FULL_ABST)
    !    else
    !       CALL newDAABS(S1%j,FULL_ABST)
    !    endif

  END FUNCTION FULL_ABST




  FUNCTION dtant( S1 )
    implicit none
    TYPE (c_taylor) dtant
    TYPE (c_taylor), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     dtant%i=0
     RETURN
    endif
    localmaster=c_master

    !    call check(s1)
    call ass(dtant)

    ! if(old) then
    call c_dafun('SIN ',s1%i,c_temp%i)
    call c_dacop(c_temp%i,dtant%i)
    call c_dafun('COS ',s1%i,c_temp%i)
    call c_dadiv(dtant%i,c_temp%i,dtant%i)
    !    else
    !       call newdafun('SIN ',s1%j,c_temp%il)
    !       call newc_dacop(c_temp%il,dtant%j)
    !       call newdafun('COS ',s1%j,c_temp%il)
    !       call newdadiv(dtant%j,c_temp%il,dtant%j)
    !    endif

    c_master=localmaster

  END FUNCTION dtant


  FUNCTION dcost( S1 )
    implicit none
    TYPE (c_taylor) dcost
    TYPE (c_taylor), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     dcost%i=0
     RETURN
    endif
    localmaster=c_master



    !    call check(s1)
    call ass(dcost)

    ! if(old) then
    call c_dafun('COS ',s1%i,c_temp%i)
    call c_dacop(c_temp%i,dcost%i)
    !    else
    !       call newdafun('COS ',s1%j,dcost%j)
    !    endif

    c_master=localmaster

  END FUNCTION dcost

  FUNCTION dsint( S1 )
    implicit none
    TYPE (c_taylor) dsint
    TYPE (c_taylor), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     dsint%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(dsint)
    ! if(old) then
    call c_dafun('SIN ',s1%i,c_temp%i)
    call c_dacop(c_temp%i,dsint%i)


    c_master=localmaster

  END FUNCTION dsint

  FUNCTION sqrtt( S1 )
    implicit none
    TYPE (c_taylor) sqrtt
    TYPE (c_taylor), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     sqrtt%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(sqrtt)
    ! if(old) then
    call c_dafun('SQRT',s1%i,c_temp%i)
    call c_dacop(c_temp%i,sqrtt%i)


    c_master=localmaster

  END FUNCTION sqrtt


  FUNCTION mul( S1, S2 )
    implicit none
    TYPE (c_taylor) mul
    TYPE (c_taylor), INTENT (IN) :: S1, S2
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     mul%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    !    call check(s2)
    call ass(mul)

    ! if(old) then
    call c_damul(s1%i,s2%i,c_temp%i)
    call c_dacop(c_temp%i,mul%i)
    !    else
    !       call newdamul(s1%j,s2%j,mul%j)
    !    endif

    c_master=localmaster

  END FUNCTION mul

  FUNCTION pbbra( S1, S2 )
    implicit none
    TYPE (c_taylor) pbbra
    TYPE (c_taylor), INTENT (IN) :: S1, S2
    integer localmaster
    integer i
    IF(.NOT.C_STABLE_DA) then
     pbbra%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    !    call check(s2)
    call ass(pbbra)

    ! if(old) then
    pbbra=0.0_dp
    do i=1,nd-rf
       pbbra=(s1.d.(2*i-1))*(s2.d.(2*i))-(s2.d.(2*i-1))*(s1.d.(2*i))+pbbra
    enddo
    !    call DAPOI(s1%i,s2%i,c_temp%i,nd)
    !    call c_dacop(c_temp%i,pbbra%i)
    !    else
    !       call newDAPOI(s1%j,s2%j,c_temp%il,nd)
    !       call newc_dacop(c_temp%il,pbbra%j)
    !    endif

    c_master=localmaster

  END FUNCTION pbbra


!  FUNCTION pbbrav( S1 )
!    implicit none
!    TYPE (c_vector_field) pbbrav
!    TYPE (c_taylor), INTENT (IN) :: S1
!    type(c_damap) s2
!    integer localmaster
!    integer i!
!
!    IF(.NOT.C_STABLE_DA) then
!     pbbrav%v%i=0
!     RETURN
!    endif
!    localmaster=c_master
!    call alloc(s2)
!     pbbrav%n=nd2
!     call c_ass_vector_field(pbbrav)
!    pbbrav=0
!s2=1
!    do i=1,nd2
!     pbbrav%v(i)=s1.pb.s2%v(i)
!    enddo
!    c_master=localmaster
!    call kill(s2)
!  END FUNCTION pbbrav

!  FUNCTION cpbbrav( S1, S2 )
!    implicit none
!    TYPE (c_vector_field) cpbbrav
!    TYPE (c_taylor), INTENT (IN) :: S1
!    type(c_damap) s2
!    integer localmaster
!    integer i
!
!    IF(.NOT.C_STABLE_DA) then
!     cpbbrav%v%i=0
!     RETURN
!    endif
!    localmaster=c_master
!    call alloc(s2)
!    s2=1
!     cpbbrav%n=nd2
!     call c_ass_vector_field(cpbbrav)
!    cpbbrav=0
!
!    do i=1,nd2
!     cpbbrav%v(i)=s1.cpb.s2%v(i)
!    enddo
!    c_master=localmaster
!    call kill(s2)
!  END FUNCTION cpbbrav


FUNCTION cpbbra( S1, S2 )
    implicit none
    TYPE (c_taylor) cpbbra
    TYPE (c_taylor), INTENT (IN) :: S1, S2
    integer localmaster

    IF(.NOT.C_STABLE_DA) then
     cpbbra%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    !    call check(s2)
    call ass(cpbbra)

    ! if(old) then
    cpbbra=n_cai*(s1.pb.s2) 
 
    c_master=localmaster

  END FUNCTION cpbbra

 

  FUNCTION liebraquaternion( S1, S2 )
    implicit none
    TYPE (c_quaternion) liebraquaternion
    TYPE (c_quaternion), INTENT (IN) :: S1, S2
    integer localmaster

    IF(.NOT.C_STABLE_DA) then
     liebraquaternion%x(1)%i=0
     RETURN
    endif

    localmaster=c_master


 
     call c_ass_quaternion(liebraquaternion)
 

    liebraquaternion=s1*s2-s2*s1


    c_master=localmaster

  END FUNCTION liebraquaternion

  FUNCTION liebra( S1, S2 )
    implicit none
    TYPE (c_vector_field) liebra
    TYPE (c_vector_field), INTENT (IN) :: S1, S2
    integer localmaster,i

    IF(.NOT.C_STABLE_DA) then
     liebra%v%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    !    call check(s2)
     liebra%n=s1%n
     call c_ass_vector_field(liebra)
    liebra=0
     call c_etcom(s1%v%i,s2%v%i,liebra%v%i,s1%n)

   


    if(use_quaternion)   THEN
     liebra%q=s2%q*s1%q-s1%q*s2%q
     liebra%q=liebra%q+c_bra_v_q(s1,s2%q)-c_bra_v_q(s2,s1%q)
   ! c_bra_v_q
    endif
         if(complex_extra_order==1.and.special_extra_order_1) liebra=liebra.cut.no
    c_master=localmaster

  END FUNCTION liebra

 

 
  FUNCTION GETORDER( S1, S2 )
    implicit none
    TYPE (c_taylor) GETORDER
    TYPE (c_taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     GETORDER%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(GETORDER)

    ! if(old) then
    CALL c_TAKE(S1%I,S2,c_temp%i)
    call c_dacop(c_temp%i,GETORDER%i)
    !    else
    !       CALL NEWTAKE(S1%J,S2,c_temp%iL)
    !       call NEWc_dacop(c_temp%iL,GETORDER%J)
    !    endif
    c_master=localmaster

  END FUNCTION GETORDER


  FUNCTION GETORDERMAP( S1, S2 )
    implicit none
    TYPE (c_damap) GETORDERMAP
    TYPE (c_damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    INTEGER I,s22
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     GETORDERMAP%v%i=0
     RETURN
    endif
    localmaster=c_master

    s22=iabs(s2)
    GETORDERMAP%n=s1%n
    call c_assmap(GETORDERMAP)
    GETORDERMAP=s1
    DO I=1,s1%n
       GETORDERMAP%V(I)=(GETORDERMAP%V(I)).SUB.S22
    ENDDO
    
    if(s2<0) s22=s22-1

if(use_quaternion)   THEN
    GETORDERMAP%q=GETORDERMAP%q.SUB.S22
else
    GETORDERMAP%s=GETORDERMAP%s.SUB.S22
endif


!    DO I=1,3
!    DO j=1,3
!       GETORDERMAP%s%s(I,j)=(GETORDERMAP%s%s(I,j)).SUB.S22
!    ENDDO
!    ENDDO

    c_master=localmaster

  END FUNCTION GETORDERMAP


  FUNCTION GETORDERquaternion( S1, S2 ) ! spin routine function
    implicit none
    TYPE (c_quaternion) GETORDERquaternion
    TYPE (c_quaternion), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    INTEGER I,J
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     GETORDERquaternion%x(1)%i=0
     RETURN
    endif
    localmaster=c_master
    call c_ass_quaternion(GETORDERquaternion)

    DO I=0,3
       GETORDERquaternion%x(I)=(S1%x(I)).SUB.S2
    ENDDO


    c_master=localmaster

  END FUNCTION GETORDERquaternion


  FUNCTION GETORDERSPINMATRIX( S1, S2 ) ! spin routine function
    implicit none
    TYPE (c_spinmatrix) GETORDERSPINMATRIX
    TYPE (c_spinmatrix), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    INTEGER I,J
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     GETORDERSPINMATRIX%s%i=0
     RETURN
    endif
    localmaster=c_master
    call c_ass_spinmatrix(GETORDERSPINMATRIX)



 
 
    DO I=1,3
    DO j=1,3
       GETORDERSPINMATRIX%s(I,j)=(S1%s(I,j)).SUB.S2
    ENDDO
    ENDDO

    c_master=localmaster

  END FUNCTION GETORDERSPINMATRIX

!!!!!!!!!!!!!!   programming extraction of nth order with parameters   !!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION from_phasor(k)
    implicit none
    TYPE (c_damap) from_phasor
    integer, optional :: k
    INTEGER I,k1
    integer localmaster
    TYPE (c_damap) from_phasori
    complex(dp) N    
    real(dp) xn
    IF(.NOT.C_STABLE_DA) then
     from_phasor%v%i=0
     RETURN
    endif
    localmaster=c_master

    
    call alloc(from_phasori)
    
    
    k1=1
    from_phasor%n=nd2
    call c_assmap(from_phasor)
    from_phasor=1
    
    if( present(k) ) k1=k

    if(symp) then
        n=sqrt(2.e0_dp)
        n_cai=1
          do i=1,ndt
             from_phasori%v(2*i-1)=((1.0_dp.cmono.(2*i-1))+i_*(1.0_dp.cmono.(2*i)))/N
             from_phasori%v(2*i)=((1.0_dp.cmono.(2*i-1))-i_*(1.0_dp.cmono.(2*i)))/N/(-i_)
          enddo
          do i=nd,nd-rf+1,-1
             from_phasori%v(2*i-1)=((1.0_dp.cmono.(2*i-1))+i_*(1.0_dp.cmono.(2*i)))/N
             from_phasori%v(2*i)=((1.0_dp.cmono.(2*i-1))-i_*(1.0_dp.cmono.(2*i)))/N/(-i_)
          enddo
          from_phasor=from_phasori**(-1)
    else   
     xn=abs(n_cai)
     if(xn>1.5_dp) then
      n=1.e0_dp
      else
      if(aimag(n_cai)/=-1) then
       Write(6,*) "n_cai can only be -2i or -i "
       stop
      endif
      n=sqrt(2.e0_dp)
     endif
          do i=1,ndt
             from_phasor%v(2*i-1)=n*((0.5_dp.cmono.(2*i-1))+(0.5_dp.cmono.(2*i)))
             from_phasor%v(2*i)=n*((0.5_dp.cmono.(2*i-1))-(0.5_dp.cmono.(2*i)))/i_
          enddo

          do i=nd,nd-rf+1,-1
             from_phasor%v(2*i-1)=n*((0.5_dp.cmono.(2*i-1))+(0.5_dp.cmono.(2*i)))
             from_phasor%v(2*i)=n*((0.5_dp.cmono.(2*i-1))-(0.5_dp.cmono.(2*i)))/i_
          enddo
     endif


     from_phasor=from_phasor**(k1)

    call kill(from_phasori)
    !!! spin later 

    c_master=localmaster

  END FUNCTION from_phasor

  FUNCTION to_phasor(k)
    implicit none
    TYPE (c_damap) to_phasor
    integer, optional :: k
    integer localmaster,k1
 
    IF(.NOT.C_STABLE_DA) then
     to_phasor%v%i=0
     RETURN
    endif
    localmaster=c_master

    k1=-1
    if(Present(k)) k1=-k
    
    
    to_phasor%n=nd2
    call c_assmap(to_phasor)

    to_phasor=from_phasor(k1)
    

    c_master=localmaster

  END FUNCTION to_phasor

  FUNCTION CUTORDER( S1, S2 )
    implicit none
    TYPE (c_taylor) CUTORDER
    TYPE (c_taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     CUTORDER%i=0
     RETURN
    endif
    localmaster=c_master

    !    call check(s1)
    call ass(CUTORDER)

    ! if(old) then
    call c_datrunc(S1%I,s2,CUTORDER%i)
    !    call c_dacop(S1%I,CUTORDER%i)

    !    DO I=S2,NO
    !       CALL TAKE(CUTORDER%I,I,c_temp%i)
    !       CALL DASUB(CUTORDER%I,c_temp%i,CUTORDER%I)
    !    ENDDO
    !    else
    !       call NEWc_dacop(S1%J,CUTORDER%J)
    !       DO I=S2,NO
    !          CALL NEWTAKE(CUTORDER%J,I,c_temp%iL)
    !          CALL NEWDASUB(CUTORDER%J,TEMPL,CUTORDER%J)
    !       ENDDO
    !    endif
    c_master=localmaster

  END FUNCTION CUTORDER

  FUNCTION CUTORDERMAP( S1, S2 )
    implicit none
    TYPE (c_DAMAP) CUTORDERMAP
    TYPE (c_DAMAP), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2
    integer localmaster,I,s22
    IF(.NOT.C_STABLE_DA) then
     CUTORDERMAP%v%i=0
     RETURN
    endif
    localmaster=c_master
    s22=iabs(s2)
     CUTORDERMAP%N=S1%N
    call C_assMAP(CUTORDERMAP)
      CUTORDERMAP=S1
    
     DO I=1,CUTORDERMAP%N
      CUTORDERMAP%V(I)=CUTORDERMAP%V(I).CUT.S22
     ENDDO
     if(s2<0) s22=s22-1
if(use_quaternion)   THEN
      CUTORDERMAP%q=CUTORDERMAP%q.cut.s22
else
      CUTORDERMAP%s=CUTORDERMAP%s.cut.s22
endif


    c_master=localmaster

  END FUNCTION CUTORDERMAP



  FUNCTION CUTORDERVEC( S1, S2 )
    implicit none
    TYPE (c_vector_field) CUTORDERVEC
    TYPE (c_vector_field), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2
    integer localmaster,I
    IF(.NOT.C_STABLE_DA) then
     CUTORDERVEC%v%i=0
     RETURN
    endif
    localmaster=c_master

     CUTORDERVEC%N=S1%N
     CUTORDERVEC%nrmax=S1%nrmax

    call C_ass_vector_field(CUTORDERVEC)

      CUTORDERVEC=S1

     DO I=1,CUTORDERVEC%N
      CUTORDERVEC%V(I)=CUTORDERVEC%V(I).CUT.S2
     ENDDO

     DO I=0,3
      CUTORDERVEC%q%x(I)= CUTORDERVEC%q%x(I).CUT.s2
     ENDDO

    c_master=localmaster

  END FUNCTION CUTORDERVEC

  FUNCTION CUTORDERspin( S1, S2 ) ! spin routine function
    implicit none
    TYPE (c_spinmatrix) CUTORDERspin
    TYPE (c_spinmatrix), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2
    integer localmaster,I,j
    IF(.NOT.C_STABLE_DA) then
     CUTORDERspin%s%i=0
     RETURN
    endif
    localmaster=c_master


    call c_ass_spinmatrix(CUTORDERspin)
      CUTORDERspin=S1

     DO I=1,3
     do j=1,3
      CUTORDERspin%s(i,j)=CUTORDERspin%s(i,j).CUT.S2
     ENDDO
    enddo

    c_master=localmaster

  END FUNCTION CUTORDERspin

  FUNCTION CUTORDERquaternion( S1, S2 ) ! spin routine function
    implicit none
    TYPE (c_quaternion) CUTORDERquaternion
    TYPE (c_quaternion), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2
    integer localmaster,I,j
    IF(.NOT.C_STABLE_DA) then
     CUTORDERquaternion%x(1)=0
     RETURN
    endif
    localmaster=c_master


    call c_ass_quaternion(CUTORDERquaternion)
      CUTORDERquaternion=S1

     DO I=0,3
      CUTORDERquaternion%x(i)=CUTORDERquaternion%x(i) 
    enddo

    c_master=localmaster

  END FUNCTION CUTORDERquaternion

  FUNCTION CUTORDERspinor( S1, S2 ) ! spin routine function
    implicit none
    TYPE (c_spinor) CUTORDERspinor
    TYPE (c_spinor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2
    integer localmaster,I
    IF(.NOT.C_STABLE_DA) then
     CUTORDERspinor%v(1)%i=0
     RETURN
    endif
    localmaster=c_master


    call c_ass_spinor(CUTORDERspinor)
      CUTORDERspinor=S1

     DO I=1,3
      CUTORDERspinor%v(i)=CUTORDERspinor%v(i).CUT.S2
     ENDDO

    c_master=localmaster

  END FUNCTION CUTORDERspinor

  FUNCTION dputchar( S1, S2 )
    implicit none
    TYPE (c_taylor) dputchar
    complex(dp), INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2
    CHARACTER (LEN = LNV)  resul
    integer j(lnv),i,io
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     dputchar%i=0
     RETURN
    endif
    localmaster=c_master


    call ass(dputchar)


    resul = trim(ADJUSTL (s2))

    do i=1,lnv
       j(i)=0
    enddo

    !frs get around compiler problem
    nd2par= len(trim(ADJUSTL (s2)))
    !frs    do i=1,len(trim(ADJUSTL (s2)))
io=0
    do i=1,nd2par
       CALL  CHARINT(RESUL(I:I),J(I))
         io=io+j(i)
       if(i>nv) then    
          if(j(i)>0) then
             dputchar=0.0_dp
             !             call var(dputchar,zero,0)
    c_master=localmaster
             return
          endif
       endif
    enddo
 
          if(io>no) then
             dputchar=0.0_dp
             !             call var(dputchar,zero,0)
    c_master=localmaster
             return
          endif    

    dputchar=0.0_dp
    !    call var(dputchar,zero,0)
    
    CALL c_pok(dputchar,j,s1)
    c_master=localmaster

  END FUNCTION dputchar

  FUNCTION dputcharr( S1r, S2 )
    implicit none
    TYPE (c_taylor) dputcharr
    real(dp), INTENT (IN) :: S1r
    CHARACTER(*)  , INTENT (IN) ::  S2
    CHARACTER (LEN = LNV)  resul
    integer j(lnv),i,io
    integer localmaster
    complex(dp) s1
    IF(.NOT.C_STABLE_DA) then
     dputcharr%i=0
     RETURN
    endif
    localmaster=c_master

    s1=s1r
    call ass(dputcharr)


    resul = trim(ADJUSTL (s2))

    do i=1,lnv
       j(i)=0
    enddo

    !frs get around compiler problem
    nd2par= len(trim(ADJUSTL (s2)))
    !frs    do i=1,len(trim(ADJUSTL (s2)))
io=0
    do i=1,nd2par
       CALL  CHARINT(RESUL(I:I),J(I))
         io=io+j(i)
       if(i>nv) then
          if(j(i)>0) then
             dputcharr=0.0_dp
             !             call var(dputchar,zero,0)
    c_master=localmaster
             return
          endif
       endif
    enddo


          if(io>no) then
             dputcharr=0.0_dp
             !             call var(dputchar,zero,0)
    c_master=localmaster
             return
          endif    
    dputcharr=0.0_dp
    !    call var(dputchar,zero,0)
    CALL c_pok(dputcharr,j,s1)
    c_master=localmaster

  END FUNCTION dputcharr

  FUNCTION dputint( S1, S2 )
    implicit none
    TYPE (c_taylor) dputint
    complex(dp), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2(:)
    integer j(lnv),i,io
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     dputint%i=0
     RETURN
    endif
    localmaster=c_master


    call ass(dputint)



    do i=1,lnv
       j(i)=0
    enddo

    !frs get around compiler problem
    nd2par= size(s2)
    !frs    do i=1,len(trim(ADJUSTL (s2)))
    do i=1,nd2par
       j(i)=s2(i)
    enddo
io=0
    do i=1,nd2par
        io=io+j(i)
       if(i>nv) then
          if(j(i)>0) then
             !             call var(dputint,zero,0)
             dputint=0.0_dp
    c_master=localmaster
             return
          endif
       endif
    enddo

 
          if(io>no) then
             dputint=0.0_dp
             !             call var(dputchar,zero,0)
    c_master=localmaster
             return
          endif    


    dputint=0.0_dp
    !    call var(dputint,zero,0)
    CALL c_pok(dputint,j,s1)
    c_master=localmaster

  END FUNCTION dputint

  FUNCTION dputintr( S1r, S2 )
    implicit none
    TYPE (c_taylor) dputintr
    real(dp), INTENT (IN) :: S1r
    integer  , INTENT (IN) ::  S2(:)
    integer j(lnv),i,io
    integer localmaster
    complex(dp) s1
    IF(.NOT.C_STABLE_DA) then
     dputintr%i=0
     RETURN
    endif
    localmaster=c_master

    s1=s1r
    call ass(dputintr)



    do i=1,lnv
       j(i)=0
    enddo

    !frs get around compiler problem
    nd2par= size(s2)
    !frs    do i=1,len(trim(ADJUSTL (s2)))
    do i=1,nd2par
       j(i)=s2(i)
    enddo
io=0
    do i=1,nd2par
        io=io+j(i)
       if(i>nv) then
          if(j(i)>0) then
             !             call var(dputint,zero,0)
             dputintr=0.0_dp
    c_master=localmaster
             return
          endif
       endif
    enddo

          if(io>no) then
             dputintr=0.0_dp
             !             call var(dputchar,zero,0)
    c_master=localmaster
             return
          endif    

    dputintr=0.0_dp
    !    call var(dputint,zero,0)
    CALL c_pok(dputintr,j,s1)
    c_master=localmaster

  END FUNCTION dputintr

   FUNCTION c_dputint0( S1, S2 )
    implicit none
    TYPE (c_taylor) c_dputint0
    complex(dp), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2
    integer j(lnv)
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     c_dputint0%i=0
     RETURN
    endif
    localmaster=c_master


    call ass(c_dputint0)
    j=0
    if(s2>nv) then
       c_dputint0=0.d0
        c_master=localmaster
       return
    endif
    if(s2==0) then
       c_dputint0=s1
       c_master=localmaster
       return
    endif

    c_dputint0=0.0_dp
    !    call var(c_dputint0,zero,s2)

    j(s2)=1
    CALL c_pok(c_dputint0,j,s1)

    c_master=localmaster

  END FUNCTION c_dputint0

   FUNCTION c_dputint0r( S1, S2 )
    implicit none
    TYPE (c_taylor) c_dputint0r
    real(dp), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2
    integer j(lnv)
    integer localmaster
    complex(dp) ss
    IF(.NOT.C_STABLE_DA) then
     c_dputint0r%i=0
     RETURN
    endif
    localmaster=c_master
   
    ss=s1

    call ass(c_dputint0r)
    j=0
    if(s2>nv) then
       c_dputint0r=0.0_dp
    c_master=localmaster
       return
    endif
    if(s2==0) then
       c_dputint0r=ss
    c_master=localmaster
       return
    endif

    c_dputint0r=0.0_dp
    !    call var(c_dputint0,zero,s2)

    j(s2)=1
    CALL c_pok(c_dputint0r,j,ss)

    c_master=localmaster

  END FUNCTION c_dputint0r


  FUNCTION GETCHARnd2s( S1, S2 )
    implicit none
    TYPE (c_taylor) GETCHARnd2s
    TYPE (c_taylor), INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2

    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     GETCHARnd2s%i=0
     RETURN
    endif
    localmaster=c_master


    call ass(GETCHARnd2s)


    GETCHARnd2s=s1.par.s2
    call  shiftda(GETCHARnd2s,GETCHARnd2s, len(trim(ADJUSTR (s2) )))

    c_master=localmaster


  END FUNCTION GETCHARnd2s

  FUNCTION GETintnd2s( S1, S2 )
    implicit none
    TYPE (c_taylor) GETintnd2s
    TYPE (c_taylor), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2(:)

    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     GETintnd2s%i=0
     RETURN
    endif
    localmaster=c_master


    call ass(GETintnd2s)


    GETintnd2s=s1.par.s2

    call  shiftda(GETintnd2s,GETintnd2s, size(s2) )

    c_master=localmaster


  END FUNCTION GETintnd2s

  FUNCTION GETintk( S1, S2 )
    implicit none
    TYPE (c_taylor) GETintk
    TYPE (c_taylor), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2

    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     GETintk%i=0
     RETURN
    endif
    localmaster=c_master


    call ass(GETintk)



    call  shiftda(s1,GETintk, s2 )

    c_master=localmaster


  END FUNCTION GETintk



  FUNCTION GETchar( S1, S2 ) 

    implicit none
    complex(dp) GETchar,r1
    TYPE (c_taylor), INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2
    CHARACTER (LEN = LNV)  resul
    integer j(lnv),i,c,cm
    IF(.NOT.C_STABLE_DA) then
     GETchar=0
     RETURN
    endif

 

    resul = s2    
    call context(resul) 

    do i=1,lnv
       j(i)=0
    enddo

 
    nd2par= len_trim(resul)  




 
    do i=1,nd2par
       CALL  CHARINT(RESUL(I:I),J(I))
    enddo

c=0
    do i=nv+1,lnv
       c=j(i)+c
    enddo
cm=0
    do i=1,nv
       cm=j(i)+cm
    enddo

!if(c>0.or.cm>nv) then  ! 2017.1.16
if(c>0.or.cm>no) then
r1=0.0_dp
else
    CALL c_dapek(S1%I,j,r1)
endif


    GETchar=r1

  END FUNCTION GETchar


 

  FUNCTION GETint( S1, S2 )
    implicit none
    complex(dp) GETint,r1
    TYPE (c_taylor), INTENT (IN) :: S1
    integer , INTENT (IN) ::  S2(:)
    integer j(lnv),i,c,cm
    IF(.NOT.C_STABLE_DA) then
     GETint=0
     RETURN
    endif

 

    do i=1,lnv
       j(i)=0
    enddo

   
    nd2par= size(s2)
   
    do i=1,nd2par
       J(I)=s2(i)
    enddo

    c=0
    do i=nv+1,lnv
       c=j(i)+c
    enddo
cm=0
    do i=1,nv
       cm=j(i)+cm
    enddo


!if(c>0.or.cm>nv) then  ! 2017.1.16
if(c>0.or.cm>no) then

r1=0.0_dp
else
    CALL c_dapek(S1%I,j,r1)
endif
 
    GETint=r1

  END FUNCTION GETint

  FUNCTION GETintmat( S1, S2 )
    implicit none
    complex(dp) GETintmat,r1
    TYPE (c_taylor), INTENT (IN) :: S1
    integer , INTENT (IN) ::  S2
    integer j(lnv),i,c,cm
    IF(.NOT.C_STABLE_DA) then
     GETintmat=0
    endif

 
    if(s2==0) then
     GETintmat=S1
     return
    endif
    do i=1,lnv
       j(i)=0
    enddo
    j(s2)=1

    c=0
    do i=nv+1,lnv
       c=j(i)+c
    enddo
cm=0
    do i=1,nv
       cm=j(i)+cm
    enddo


!if(c>0.or.cm>nv) then  ! 2017.1.16
if(c>0.or.cm>no) then

r1=0.0_dp
else
    CALL c_dapek(S1%I,j,r1)
endif
 
    GETintmat=r1

  END FUNCTION GETintmat



  FUNCTION GETdiff( S1, S2 )
    implicit none
    TYPE (c_taylor) GETdiff
    TYPE (c_taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     GETdiff%i=0
     RETURN
    endif


     localmaster=c_master
    

    !    call check(s1)
    call ass(GETdiff)
 
    ! if(old) then
    CALL c_dader(S2,S1%I,c_temp%i)
!    call c_dacop(c_temp%i,GETdiff%i)

    getdiff=c_temp
    !    else
    !       CALL NEWdader(S2,S1%J,TEMPL)
    !       call NEWc_dacop(tempL,GETdiff%J)
    !    endif
 
 
    c_master=localmaster

  END FUNCTION GETdiff

  FUNCTION GETINTegrate( S1, S2 )
    implicit none
    TYPE (c_taylor) GETINTegrate
    TYPE (c_taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2
    integer localmaster,n,i
    type(c_taylor) t,x
    complex(dp) value
    integer, allocatable :: jc(:)
    IF(.NOT.C_STABLE_DA) then
     GETINTegrate%i=0
     RETURN
    endif

    localmaster=c_master

    allocate(jc(nv))
    jc=0
    !    call check(s1)
    call ass(GETINTegrate)
    call alloc(t,x)
    t=s1
    x=0
    call c_taylor_cycle(t,size=n)

    do i=1,n
       call c_taylor_cycle(t,ii=i,value=value,j=jc)
 
         x=((value/(jc(s2)+1)).cmono.jc)*((1.0_dp,0.0_dp).cmono.s2)+x

    enddo

    GETINTegrate=x

    call kill(t,x)
    deallocate(jc)
    c_master=localmaster

  END FUNCTION GETINTegrate

! computes Lie polynomial in relevant planes : s2=1 generating function, ss=0 potential
! Does not included the nonsymplectic modulated planes
  FUNCTION getpb( S1,S1p, S2 )  
    implicit none
    TYPE (c_taylor) getpb
    TYPE (c_vector_field),optional, INTENT (IN) :: S1
    TYPE (c_damap),optional, INTENT (IN) :: S1p
    INTEGER,optional, INTENT (IN) ::  S2
    integer localmaster,n,i,j,l,fac,nd2here,ss
    type(c_taylor) t,x
    complex(dp) value
    integer, allocatable :: jc(:)
    IF(.NOT.C_STABLE_DA) then
     getpb%i=0
     RETURN
    endif

    localmaster=c_master

    nd2here=nd2-2*rf
    ss=-1
    if(present(s2)) ss=s2
    allocate(jc(c_%nv))
    jc=0
    !    call check(s1)
    call ass(getpb)
    getpb=0.0_dp
    call alloc(t,x)

    do j=1,nd2here
    if(present(s1)) then
     t=s1%v(j)
    elseif(present(s1p)) then
     t=s1p%v(j)
     if(.not.present(s2))ss=1
     else
     write(6,*) " error in getpb "
     stop
    endif
    x=0
    call c_taylor_cycle(t,size=n)

    do i=1,n
       call c_taylor_cycle(t,ii=i,value=value,j=jc)
        fac=0
        do l=1,nd2here
         fac=jc(l)+fac
        enddo
        fac=fac+1
       if(ss<0) then  !  fixed bug 2017 jan 9
        if(mod(j,2)==0) then  
          x=((value/fac).cmono.jc)*(1.0_dp.cmono.(j-1))+x
        else
          x=ss*((value/fac).cmono.jc)*(1.0_dp.cmono.(j+1))+x        
        endif
       else
          x=((value/fac).cmono.jc)*(1.0_dp.cmono.(j))+x
       endif
    enddo

    getpb=x+getpb


    enddo ! j

    call kill(t,x)
    deallocate(jc)
    c_master=localmaster

  END FUNCTION getpb

    FUNCTION cgetpb( S1, s1p,S2 )  
    implicit none
    TYPE (c_taylor) cgetpb
    TYPE (c_vector_field),optional, INTENT (IN) :: S1
    INTEGER, optional,INTENT (IN) ::  S2
    TYPE (c_damap),optional, INTENT (IN) :: S1p
    integer localmaster 
    localmaster=master

 
 
    call ass(cgetpb)
     cgetpb=getpb( s1=S1, s1p=s1p,s2=S2 )/n_cai
    c_master=localmaster
    END FUNCTION cgetpb

  FUNCTION getpb_from_transverse( S1,f,S1p, S2 )  
    implicit none
    TYPE (c_taylor) getpb_from_transverse
    TYPE (c_vector_field),optional, INTENT (IN) :: S1
    TYPE (c_damap),optional, INTENT (IN) :: S1p
    INTEGER,optional, INTENT (IN) ::  S2
    integer localmaster,n,i,j,l,fac,nd2here,ss,k
    type(c_taylor) t,x
    complex(dp) value
    real(dp) c
    integer, allocatable :: jc(:)
    type(c_vector_field), INTENT (INout) :: f
    IF(.NOT.C_STABLE_DA) then
     getpb_from_transverse%i=0
     RETURN
    endif

    localmaster=c_master

    nd2here=4 !nd2-2*rf
    ss=-1
    if(present(s2)) ss=s2
    allocate(jc(c_%nv))
    jc=0
    !    call check(s1)
    call ass(getpb_from_transverse)
    getpb_from_transverse=0.0_dp
    call alloc(t,x)
    f=0

    do j=1,nd2here
    if(present(s1)) then
     t=s1%v(j)
    elseif(present(s1p)) then
     t=s1p%v(j)
     if(.not.present(s2))ss=1
     else
     write(6,*) " error in getpb "
     stop
    endif
    x=0 
    call c_taylor_cycle(t,size=n)

    do i=1,n
       call c_taylor_cycle(t,ii=i,value=value,j=jc)
        fac=0
        do l=1,nd2here
         fac=jc(l)+fac
        enddo
        fac=fac+1
       if(ss<0) then  !  fixed bug 2017 jan 9
        if(mod(j,2)==0) then  
          x=((value/fac).cmono.jc)*(1.0_dp.cmono.(j-1))+x
          jc(j-1)=jc(j-1)+1
          value=value/fac
        else
          x=ss*((value/fac).cmono.jc)*(1.0_dp.cmono.(j+1))+x  
          jc(j+1)=jc(j+1)+1  
          value=ss*value/fac    
        endif
        do k=1,nd
         l=2*k
          call derive(jc,l,c)
          if(jc(l)>=0) f%v(2*k-1)=f%v(2*k-1)-((value*c).cmono.jc)
          jc(2*k)=jc(2*k)+1
         l=2*k-1
          call derive(jc,l,c)
          if(jc(l)>=0) f%v(2*k)=f%v(2*k)+((value*c).cmono.jc)
          jc(2*k-1)=jc(2*k-1)+1
        enddo
       else
          x=((value/fac).cmono.jc)*(1.0_dp.cmono.(j))+x
       endif
    enddo
    getpb_from_transverse=x+getpb_from_transverse
    

    enddo ! j

    call kill(t,x)
    deallocate(jc)
    c_master=localmaster

  END FUNCTION getpb_from_transverse

  subroutine derive(j,k,c)
  implicit none
   integer, intent(inout):: J(:)
   integer, intent(in):: k
   real(dp), intent(out) :: c

   if(j(k)==0) then 
     c=0.d0
    j(k)=j(k)-1
   else
    c=j(k)
    j(k)=j(k)-1
   endif
  end subroutine derive



      FUNCTION getvectorfield( S1,s2 )  
!#internal: getvectorfield
!# produce vector field S1.grad= :s2:
!
    implicit none
    TYPE (c_vector_field) getvectorfield
    TYPE (c_taylor), INTENT (IN) :: S1
    integer, optional, INTENT (IN) :: S2
    integer localmaster ,i, ss
    localmaster=c_master
    
        getvectorfield%n=nd2
    call c_ass_vector_field(getvectorfield)
    getvectorfield=0
    ss=-1
    if(present(s2)) ss=s2
    if(ss==-1) then
     do i=1,nd2-2*rf
      getvectorfield%v(i)= s1.pb.(1.0_dp.cmono.i)
     enddo
    else
     do i=1,nd2-2*rf
      getvectorfield%v(i)= s1.d.i
     enddo
    endif
    c_master=localmaster
     if(complex_extra_order==1.and.special_extra_order_1) getvectorfield=getvectorfield.cut.no

      END FUNCTION getvectorfield 
    
    FUNCTION cgetvectorfield( S1 )  
    implicit none
    TYPE (c_vector_field) cgetvectorfield
    TYPE (c_taylor), INTENT (IN) :: S1
 
    integer localmaster ,i 
    localmaster=master
    
    cgetvectorfield%n=nd2
    call c_ass_vector_field(cgetvectorfield)

cgetvectorfield=0
 
     do i=1,nd2-2*rf
      cgetvectorfield%v(i)= s1.cpb.(1.0_dp.cmono.i)
     enddo
 
     c_master=localmaster
     if(complex_extra_order==1.and.special_extra_order_1) cgetvectorfield=cgetvectorfield.cut.no
    END FUNCTION cgetvectorfield 
    

  FUNCTION GETdatra( S1, S2 )
    implicit none
    TYPE (c_taylor) GETdatra
    TYPE (c_taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     GETdatra%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(GETdatra)

    ! if(old) then
    CALL c_datra(S2,S1%I,c_temp%i)
    call c_dacop(c_temp%i,GETdatra%i)
    !    else
    !       CALL NEWdatra(S2,S1%J,TEMPL)
    !       call NEWc_dacop(tempL,GETdatra%J)
    !    endif
    c_master=localmaster

  END FUNCTION GETdatra

  FUNCTION POW( S1, R2 )
    implicit none
    TYPE (c_taylor) POW
    TYPE (c_taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    INTEGER I,R22
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     POW%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(POW)

    ! if(old) then
    CALL c_dacon(c_temp%i,(1.0_dp,0.0_dp))

    R22=IABS(R2)
    DO I=1,R22
       CALL c_damul(c_temp%i,S1%I,c_temp%i)
    ENDDO
    IF(R2.LT.0) THEN
       CALL c_dadic(c_temp%i,(1.0_dp,0.0_dp),c_temp%i)
    ENDIF
    call c_dacop(c_temp%i,POW%i)

    c_master=localmaster
  END FUNCTION POW

  
  FUNCTION cdmulsc( S1, sc )
    implicit none
    TYPE (c_taylor) cdmulsc
    TYPE (c_taylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     cdmulsc%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(cdmulsc)

    ! if(old) then
    call c_dacmu(s1%i,sc,c_temp%i)
    call c_dacop(c_temp%i,cdmulsc%i)
    !    else
    !       call newc_dacmu(s1%j,sc,dmulsc%j)
    !    endif

    c_master=localmaster
  END FUNCTION cdmulsc

  FUNCTION dmulsc( S1, sc )
    implicit none
    TYPE (c_taylor) dmulsc
    TYPE (c_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    complex(dp) sct
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     dmulsc%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(dmulsc)
    sct=sc
    ! if(old) then
    call c_dacmu(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,dmulsc%i)
    !    else
    !       call newc_dacmu(s1%j,sc,dmulsc%j)
    !    endif

    c_master=localmaster
  END FUNCTION dmulsc

  FUNCTION mulsc( S1, sc )
    implicit none
    TYPE (c_taylor) mulsc
    TYPE (c_taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    complex(dp) sct
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     mulsc%i=0
     RETURN
    endif
    localmaster=c_master


    if(c_real_warning) call c_real_stop
    !    call check(s1)
    call ass(mulsc)
     sct=sc
    ! if(old) then
    call c_dacmu(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,mulsc%i)
    !    else
    !       call newc_dacmu(s1%j,REAL(sc,kind=DP),mulsc%j)
    !    endif
    c_master=localmaster
  END FUNCTION mulsc

  FUNCTION imulsc( S1, sc )
    implicit none
    TYPE (c_taylor) imulsc
    TYPE (c_taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    integer localmaster
    complex(dp) sct
    IF(.NOT.C_STABLE_DA) then
     imulsc%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(imulsc)
    sct=sc

    ! if(old) then
    call c_dacmu(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,imulsc%i)
    !    else
    !       call newc_dacmu(s1%j,REAL(sc,kind=DP),imulsc%j)
    !    endif

    c_master=localmaster
  END FUNCTION imulsc

  FUNCTION cdscmul( sc,S1 )
    implicit none
    TYPE (c_taylor) cdscmul
    TYPE (c_taylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     cdscmul%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(cdscmul)

    ! if(old) then
    call c_dacmu(s1%i,sc,c_temp%i)
    call c_dacop(c_temp%i,cdscmul%i)
    !    else
    !       call newc_dacmu(s1%j,sc,dscmul%j)
    !    endif

    c_master=localmaster

  END FUNCTION cdscmul

  FUNCTION dscmul( sc,S1 )
    implicit none
    TYPE (c_taylor) dscmul
    TYPE (c_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    complex(dp) sct
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     dscmul%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(dscmul)
    sct=sc
    ! if(old) then
    call c_dacmu(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,dscmul%i)
    !    else
    !       call newc_dacmu(s1%j,sc,dscmul%j)
    !    endif

    c_master=localmaster

  END FUNCTION dscmul

  FUNCTION scmul( sc,S1 )
    implicit none
    TYPE (c_taylor) scmul
    TYPE (c_taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    complex(dp) sct
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     scmul%i=0
     RETURN
    endif
    localmaster=c_master


    if(c_real_warning) call c_real_stop
    !    call check(s1)
    call ass(scmul)
    
    sct=sc
   
    ! if(old) then
    call c_dacmu(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,scmul%i)
    !    else
    !       call newc_dacmu(s1%j,REAL(sc,kind=DP),scmul%j)
    !    endif

    c_master=localmaster

  END FUNCTION scmul

  FUNCTION iscmul( sc,S1 )
    implicit none
    TYPE (c_taylor) iscmul
    TYPE (c_taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    complex(dp) sct
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     iscmul%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(iscmul)
    sct=sc
    ! if(old) then
    call c_dacmu(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,iscmul%i)
    !    else
    !       call newc_dacmu(s1%j,REAL(sc,kind=DP),iscmul%j)
    !    endif

    c_master=localmaster

  END FUNCTION iscmul

  FUNCTION div( S1, S2 )
    implicit none
    TYPE (c_taylor) div
    TYPE (c_taylor), INTENT (IN) :: S1, S2
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     div%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    !    call check(s2)
    call ass(div)

    ! if(old) then
    call c_dadiv(s1%i,s2%i,c_temp%i)
    call c_dacop(c_temp%i,div%i)
    !    else
    !       call newdadiv(s1%j,s2%j,templ)
    !       call newc_dacop(templ,div%j)
    !    endif

    c_master=localmaster
  END FUNCTION div

  FUNCTION cdscdiv( sc,S1 )
    implicit none
    TYPE (c_taylor) cdscdiv
    TYPE (c_taylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     cdscdiv%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(cdscdiv)

    ! if(old) then
    call c_dadic(s1%i,sc,c_temp%i)
    call c_dacop(c_temp%i,cdscdiv%i)
    !    else
    !       call newc_dadic(s1%j,sc,dscdiv%j)
    !    endif

    c_master=localmaster

  END FUNCTION cdscdiv

  FUNCTION dscdiv( sc,S1 )
    implicit none
    TYPE (c_taylor) dscdiv
    TYPE (c_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    complex(dp) sct
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     dscdiv%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(dscdiv)
    sct=sc

    ! if(old) then
    call c_dadic(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,dscdiv%i)
    !    else
    !       call newc_dadic(s1%j,sc,dscdiv%j)
    !    endif

    c_master=localmaster

  END FUNCTION dscdiv

  FUNCTION scdiv( sc,S1 )
    implicit none
    TYPE (c_taylor) scdiv
    TYPE (c_taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    complex(dp) sct
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     scdiv%i=0
     RETURN
    endif
    localmaster=c_master


    if(c_real_warning) call c_real_stop
    !    call check(s1)
    call ass(scdiv)
    
     sct=sc

    ! if(old) then
    call c_dadic(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,scdiv%i)
    !    else
    !       call newc_dadic(s1%j,REAL(sc,kind=DP),scdiv%j)
    !    endif

    c_master=localmaster
  END FUNCTION scdiv

  FUNCTION iscdiv( sc,S1 )
    implicit none
    TYPE (c_taylor) iscdiv
    TYPE (c_taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    complex(dp) sct
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     iscdiv%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(iscdiv)
    sct=sc
    ! if(old) then
    call c_dadic(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,iscdiv%i)
    !    else
    !       call newc_dadic(s1%j,REAL(sc,kind=DP),iscdiv%j)
    !    endif


    c_master=localmaster
  END FUNCTION iscdiv

  FUNCTION cddivsc( S1, sc )
    implicit none
    TYPE (c_taylor) cddivsc
    TYPE (c_taylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     cddivsc%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(cddivsc)


    ! if(old) then
    call c_dacdi(s1%i,sc,c_temp%i)
    call c_dacop(c_temp%i,cddivsc%i)
    !    else
    !       call newc_dacdi(s1%j,sc,ddivsc%j)
    !    endif
    c_master=localmaster

  END FUNCTION cddivsc


  FUNCTION ddivsc( S1, sc )
    implicit none
    TYPE (c_taylor) ddivsc
    TYPE (c_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    complex(dp) sct
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     ddivsc%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(ddivsc)

    sct=sc
    ! if(old) then
    call c_dacdi(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,ddivsc%i)
    !    else
    !       call newc_dacdi(s1%j,sc,ddivsc%j)
    !    endif
    c_master=localmaster

  END FUNCTION ddivsc

  FUNCTION divsc( S1, sc )
    implicit none
    TYPE (c_taylor) divsc
    TYPE (c_taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    complex(dp) sct
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     divsc%i=0
     RETURN
    endif
    localmaster=c_master


    if(c_real_warning) call c_real_stop
    !    call check(s1)
    call ass(divsc)
    sct=sc
    ! if(old) then
    call c_dacdi(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,divsc%i)
    !    else
    !       call newc_dacdi(s1%j,REAL(sc,kind=DP),divsc%j)
    !    endif
    c_master=localmaster

  END FUNCTION divsc


  FUNCTION idivsc( S1, sc )
    implicit none
    TYPE (c_taylor) idivsc
    TYPE (c_taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    integer localmaster
    complex(dp) sct
    IF(.NOT.C_STABLE_DA) then
     idivsc%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(idivsc)
    
    sct=sc

    ! if(old) then
    call c_dacdi(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,idivsc%i)
    !    else
    !       call newc_dacdi(s1%j,REAL(sc,kind=DP),idivsc%j)
    !    endif
    c_master=localmaster

  END FUNCTION idivsc


  FUNCTION add( S1, S2 )
    implicit none
    TYPE (c_taylor) add
    TYPE (c_taylor), INTENT (IN) :: S1, S2
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     add%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    !    call check(s2)
    call ass(add)


    ! if(old) then
    call c_daadd(s1%i,s2%i,add%i)
    !  call c_dacop(c_temp%i,add%i)
    !  call daadd(s1%i,s2%i,c_temp%i)
    !  call c_dacop(c_temp%i,add%i)
    !    else
    !       call newdaadd(s1%j,s2%j,add%j)
    !    endif

    c_master=localmaster

  END FUNCTION add

  FUNCTION addq( S1, S2 )
    implicit none
    TYPE (c_quaternion) addq
    TYPE (c_quaternion), INTENT (IN) :: S1, S2
    integer i,localmaster
              localmaster=c_master
              call c_ass_quaternion(addq)
       do i=0,3
        addq%x(i)=s1%x(i)+s2%x(i)
       enddo
          c_master=localmaster
 
  END FUNCTION addq

  FUNCTION mulq( S1, S2 )
    implicit none
    TYPE (c_quaternion) mulq
    TYPE (c_quaternion), INTENT (IN) :: S1, S2
    type(c_taylor) temp(0:3)
    integer i,localmaster

              localmaster=c_master
              call c_ass_quaternion(mulq)
 
       call alloc(temp)
 

          temp(0)=s1%x(0)*s2%x(0)-s1%x(1)*s2%x(1)-s1%x(2)*s2%x(2)-s1%x(3)*s2%x(3)
          
         temp(1)= s1%x(2)*s2%x(3)-s1%x(3)*s2%x(2)
         temp(2)= s1%x(3)*s2%x(1)-s1%x(1)*s2%x(3)
         temp(3)= s1%x(1)*s2%x(2)-s1%x(2)*s2%x(1)

        do i=1,3
         temp(i)= temp(i) + s1%x(0)*s2%x(i)+ s1%x(i)*s2%x(0)
        enddo

        do i=0,3
              mulq%x(i)=temp(i)
       enddo
               c_master=localmaster    
       call kill(temp)
  END FUNCTION mulq

  FUNCTION mulcq( S1, S2 )
    implicit none
    TYPE (c_quaternion) mulcq
    complex(dp), INTENT (IN) :: S1
    TYPE (c_quaternion), INTENT (IN) ::  S2
    integer i,localmaster
              localmaster=c_master
              call c_ass_quaternion(mulcq)

        do i=0,3
              mulcq%x(i)=s1*s2%x(i)
        enddo
               c_master=localmaster   
  END FUNCTION mulcq



  FUNCTION subq( S1, S2 )
    implicit none
    TYPE (c_quaternion) subq
    TYPE (c_quaternion), INTENT (IN) :: S1, S2
    integer i,localmaster
              localmaster=c_master
              call c_ass_quaternion(subq)
       do i=0,3
                 subq%x(i)=s1%x(i)-s2%x(i)
       enddo
               c_master=localmaster 

  END FUNCTION subq

  FUNCTION c_invq( S1 )
    implicit none
    TYPE (c_quaternion) c_invq
    TYPE (c_quaternion), INTENT (IN) :: S1
    type(c_taylor) norm
    type(c_taylor) temp(0:3)
    integer i,localmaster

    IF(.NOT.C_%STABLE_DA) then
     c_invq%x(1)=0
     RETURN
    endif
         localmaster=c_master
         call c_ass_quaternion(c_invq)
       call alloc(temp)
      call alloc(norm)
 if(assume_c_quaternion_normalised) then
              norm=1.0_dp
else
              norm=s1%x(0)**2+s1%x(1)**2+s1%x(2)**2+s1%x(3)**2
endif
            temp(0)=s1%x(0)
              do i=1,3
                temp(i)=-s1%x(i)
              enddo
              do i=0,3
                temp(i)=temp(i)/norm
              enddo

    do i=0,3
         c_invq%x(i)=temp(i)
    enddo

          c_master=localmaster
          call kill(norm)
          call kill(temp)

 
  END FUNCTION c_invq

 

  FUNCTION POWq( S1, R2 )
    implicit none
    TYPE (c_quaternion) POWq,qtemp
    TYPE (c_quaternion), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    INTEGER I,R22
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       POWq=0.0_dp

      RETURN
    endif
         localmaster=c_master
    call c_ass_quaternion(powq)

     do i=0,3
     call alloc(qtemp%x(i))  
    enddo  
     qtemp=1.0_dp

    R22=IABS(R2)
    DO I=1,R22
       qtemp=qtemp*s1
    ENDDO
    IF(R2.LT.0) THEN
       qtemp=c_invq(qtemp)
    ENDIF

     do i=0,3
         powq%x(i)=qtemp%x(i)
     enddo
     do i=0,3
     call kill(qtemp%x(i))  
    enddo  
          c_master=localmaster
  END FUNCTION POWq

  FUNCTION POWql( S1, R2 )
    implicit none
    TYPE (q_linear) POWql,qtemp
    TYPE (q_linear), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    INTEGER I,R22
    integer localmaster
 
 

      qtemp=0
!     qtemp%q=0
!     qtemp%q(0,0)=1.0_dp

    R22=IABS(R2)
    DO I=1,R22
       qtemp=qtemp*s1
    ENDDO
    IF(R2.LT.0) THEN
       qtemp=inv_q_linear(qtemp)
    ENDIF
      powql=qtemp
 

   
  END FUNCTION POWql

  SUBROUTINE  EQUALq(S2,S1)
    implicit none
    integer ipause, mypauses
    type (c_quaternion),INTENT(inOUT)::S2
    type (c_quaternion),INTENT(IN)::S1
    integer i

      do i=0,3
        s2%x(i)=s1%x(i)
      enddo

 end   SUBROUTINE  EQUALq

  SUBROUTINE  EQUALq_c_r(S2,S1)
    implicit none
    integer ipause, mypauses
    type (c_quaternion),INTENT(inOUT)::S2
    type (quaternion),INTENT(IN)::S1
    integer i
 
 
      do i=0,3
        s2%x(i)=s1%x(i) 
      enddo
  

 end   SUBROUTINE  EQUALq_c_r

  SUBROUTINE  EQUALq_r_c(S2,S1)
    implicit none
    integer ipause, mypauses
    type (quaternion),INTENT(inOUT)::S2
    type (c_quaternion),INTENT(IN)::S1
    integer i


      do i=0,3

        s2%x(i)=s1%x(i)
      enddo

 end   SUBROUTINE  EQUALq_r_c

  SUBROUTINE  EQUALql_i(S2,S1)
    implicit none
    integer ipause, mypauses
    type (q_linear),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1
    integer i 

      S2%q(0:3,0:6)=0
      s2%mat=0.0_dp
       do i=1,6
      s2%mat(i,i)=1.0_dp
       enddo
      S2%q(s1,0)= 1.0_dp

 end   SUBROUTINE  EQUALql_i

  SUBROUTINE  EQUALql_r(S2,S1)
    implicit none
    integer ipause, mypauses
    type (q_linear),INTENT(inOUT)::S2
    real(dp) ,INTENT(IN)::S1
    integer i,j

      S2%q(0:3,0:6)=0
 
      S2%q(0,0)= s1

 end   SUBROUTINE  EQUALql_r

  function  qua_ql(S2)
    implicit none
    integer ipause, mypauses
    type (q_linear) qua_ql
    type (q_linear),INTENT(inOUT)::S2
 
      qua_ql=0
      qua_ql%q = s2%q 


 end   function  qua_ql



  SUBROUTINE  EQUALql_q(S2,S1)
    implicit none
    integer ipause, mypauses
    type (q_linear),INTENT(inOUT)::S2
    type (c_quaternion),INTENT(IN)::S1
    integer i,j

      S2%q=0
      do i=0,3
       do j=0,min(6,nd2)
        s2%q(i,j)=s1%x(i).index.j
       enddo
      enddo

 end   SUBROUTINE  EQUALql_q

  SUBROUTINE  EQUALql_cmap(S2,S1)
    implicit none
    integer ipause, mypauses
    type (q_linear),INTENT(inOUT)::S2
    type (c_damap),INTENT(IN)::S1
  
       s2%mat=s1
       s2=s1%q
      
 end   SUBROUTINE  EQUALql_cmap

  SUBROUTINE  EQUALcmap_ql(S1,S2)
    implicit none
    integer ipause, mypauses
    type (q_linear),INTENT(in)::S2
    type (c_damap),INTENT(INOUT)::S1
  
       s1=s2%mat
       s1%q=s2 
      
 end   SUBROUTINE  EQUALcmap_ql


  SUBROUTINE  EQUALq_ql(S1,S2)
    implicit none
    integer ipause, mypauses
    type (q_linear),INTENT(in)::S2
    type (c_quaternion),INTENT(INout)::S1
    integer i,j

      s1=0.0_dp

      do i=0,3
        s1%x(i) = s2%q(i,0)   +s1%x(i)
       do j=1,min(6,nd2)
        s1%x(i)= s2%q(i,j)*dz_c(j)+s1%x(i)
       enddo
      enddo

 end   SUBROUTINE  EQUALq_ql


  SUBROUTINE  EQUAL_c_quaternion_complex_quaternion(S1,S2)
    implicit none
    type (complex_quaternion),INTENT(in)::S2
    type (c_quaternion),INTENT(INout)::S1
    integer i 

 

      do i=0,3
        s1%x(i) = s2%x(i)  
      enddo

 end   SUBROUTINE  EQUAL_c_quaternion_complex_quaternion

  SUBROUTINE  EQUAL_complex_quaternion_c_quaternion(S1,S2)
    implicit none
    type (c_quaternion),INTENT(in)::S2
    type (complex_quaternion),INTENT(INout)::S1
    integer i 

 

      do i=0,3
        s1%x(i) = s2%x(i)  
      enddo


 end   SUBROUTINE  EQUAL_complex_quaternion_c_quaternion



  SUBROUTINE  EQUALql_ql(S1,S2)
    implicit none
    type (q_linear),INTENT(in)::S2
    type (q_linear),INTENT(INout)::S1
 

   s1%mat=s2%mat
   S1%q(0:3,0:6)=S2%q(0:3,0:6)
 

 end   SUBROUTINE  EQUALql_ql

  SUBROUTINE  print_ql(S2,imaginary,mf)
    implicit none
    type (q_linear),INTENT(inOUT)::S2
    integer, optional :: mf
    logical, optional :: imaginary
    integer i,mff

      mff=6
     if(present(mf)) mff=mf 

write(mff,*) " Orbital Matrix "
if(present(imaginary) )write(mff,*) "Real part "
      do i=1,6
 
         write(mff,'(6(1x,G21.14))') real(s2%mat(i,1:min(6,nd2)))
       
      enddo
if(present(imaginary) ) then
if(imaginary) then
write(mff,*) "Imaginary part "
      do i=1,6
 
         write(mff,'(6(1x,G21.14))') aimag(s2%mat(i,1:min(6,nd2)))
       
      enddo
endif
endif



write(mff,*) " Quaternion Matrix "
if(present(imaginary) )write(mff,*) "Real part "

      do i=0,3
 
         write(mff,'(7(1x,G21.14))') real(s2%q(i,0:min(6,nd2)))
       
      enddo

if(present(imaginary) ) then
if(imaginary) then
write(mff,*) "Imaginary part "

      do i=0,3
 
         write(mff,'(7(1x,G21.14))') aimag(s2%q(i,0:min(6,nd2)))
       
      enddo
endif
endif

 end   SUBROUTINE  print_ql

  function   inv_symplectic66(S1)
    implicit none
    integer ipause, mypauses
    real(dp) inv_symplectic66(6,6)
    real(dp),intent(in) :: s1(6,6)

      inv_symplectic66=-matmul(matmul(sj,transpose(s1)),sj)

 end   function  inv_symplectic66

  function   inv_q_linear(S1)
    implicit none
    integer ipause, mypauses
    type(q_linear) inv_q_linear 
    type(q_linear),intent (IN) :: s1 
    integer i
       call c_matinv(s1%mat,inv_q_linear%mat,6,6,i)
      inv_q_linear%q(0,0:6)=s1%q(0,0:6)
      inv_q_linear%q(1:3,0:6)=-s1%q(1:3,0:6)

      inv_q_linear=inv_q_linear*inv_q_linear%mat

 end   function  inv_q_linear

  function   mulqdiv(S1,s2)
    implicit none
    integer ipause, mypauses
    type(q_linear) mulqdiv 
    type(q_linear),intent (IN) :: s1 ,s2

     mulqdiv=s1*inv_q_linear(s2)

 end   function  mulqdiv


  function   mul_ql_m(S1,S2)
    implicit none
    integer ipause, mypauses
    type (q_linear) mul_ql_m
    type (q_linear),INTENT(IN)::S1
    real(dp), INTENT(IN)::S2(6,6)
    integer i,j

      mul_ql_m=s1

 !     mul_ql_m%q(0:3,0)=S1%q(0:3,0)

      do i=0,3
       mul_ql_m%q(i,1:6)= matmul(s1%q(i,1:6),s2)
      enddo

 end   function  mul_ql_m


  function   mul_ql_cm(S1,S2)
    implicit none
    integer ipause, mypauses
    type (q_linear) mul_ql_cm
    type (q_linear),INTENT(IN)::S1
    complex(dp), INTENT(IN)::S2(6,6)
    integer i,j

      mul_ql_cm=s1

 !     mul_ql_m%q(0:3,0)=S1%q(0:3,0)

      do i=0,3
       mul_ql_cm%q(i,1:6)= matmul(s1%q(i,1:6),s2)
      enddo

 end   function  mul_ql_cm

  function   addql(S1,S2)
    implicit none
    integer ipause, mypauses
    type (q_linear) addql
    type (q_linear),INTENT(IN)::S1
    type (q_linear),INTENT(IN)::S2
    integer i,j

      addql%q=s1%q+s2%q
      addql%mat=s1%mat+s2%mat


 end   function  addql


  function   mulql(S1,S2)
    implicit none
    integer ipause, mypauses
    type (q_linear) mulql
    type (q_linear),INTENT(IN)::S1
    type (q_linear),INTENT(IN)::S2
    complex(dp) temp(0:3),q0(0:3),p0(0:3),p1(0:3,6),q1(0:3,6)
    integer i,j
    type(q_linear) qt
        
 
        qt=s1
        qt=qt*s2%mat

        mulql%q=0
        mulql%mat=matmul(s1%mat,s2%mat)
        
        q0(0:3)=qt%q(0:3,0)
        p0(0:3)=s2%q(0:3,0)
        q1(0:3,1:6)=qt%q(0:3,1:6)
        p1(0:3,1:6)=s2%q(0:3,1:6)        

         temp(0)= q0(0)*p0(0)-q0(1)*p0(1)-q0(2)*p0(2)-q0(3)*p0(3)
         temp(1)= q0(2)*p0(3)-q0(3)*p0(2)
         temp(2)= q0(3)*p0(1)-q0(1)*p0(3)
         temp(3)= q0(1)*p0(2)-q0(2)*p0(1)

        do i=1,3
         temp(i)= temp(i) + q0(0)*p0(i)+ q0(i)*p0(0)
        enddo
              mulql%q(0:3,0)=temp(0:3)
       
do j=1,6
           mulql%q(0,j)= q0(0)*p1(0,j)-q0(1)*p1(1,j)-q0(2)*p1(2,j)-q0(3)*p1(3,j)+mulql%q(0,j)
           mulql%q(1,j)= q0(2)*p1(3,j)-q0(3)*p1(2,j)+mulql%q(1,j)
           mulql%q(2,j)= q0(3)*p1(1,j)-q0(1)*p1(3,j)+mulql%q(2,j)
           mulql%q(3,j)= q0(1)*p1(2,j)-q0(2)*p1(1,j)+mulql%q(3,j)

        do i=1,3
         mulql%q(i,j)= mulql%q(i,j) + q0(0)*p1(i,j)+ q0(i)*p1(0,j)
        enddo

           mulql%q(0,j)=  p0(0)*q1(0,j)-p0(1)*q1(1,j)-p0(2)*q1(2,j)-p0(3)*q1(3,j)+mulql%q(0,j)
           mulql%q(1,j)= -p0(2)*q1(3,j)+p0(3)*q1(2,j)+mulql%q(1,j)
           mulql%q(2,j)= -p0(3)*q1(1,j)+p0(1)*q1(3,j)+mulql%q(2,j)
           mulql%q(3,j)= -p0(1)*q1(2,j)+p0(2)*q1(1,j)+mulql%q(3,j)

        do i=1,3
         mulql%q(i,j)= mulql%q(i,j) + p0(0)*q1(i,j)+ p0(i)*q1(0,j)
        enddo
        
enddo

 end   function  mulql



  function   subql(S1,S2)
    implicit none
    integer ipause, mypauses
    type (q_linear) subql
    type (q_linear),INTENT(IN)::S1
    type (q_linear),INTENT(IN)::S2
    integer i,j

      subql%q=s1%q-s2%q
      subql%mat=s1%mat-s2%mat


 end   function  subql



  SUBROUTINE  EQUALq_c_8(S2,S1)
    implicit none
    integer ipause, mypauses
    type (c_quaternion),INTENT(inOUT)::S2
    type (quaternion_8),INTENT(IN)::S1
    integer i
    type(complextaylor) ct
    call alloc(ct)
      do i=0,3
        ct=s1%x(i)%t
        s2%x(i)=ct
      enddo
    call kill(ct)

 end   SUBROUTINE  EQUALq_c_8

  SUBROUTINE  EQUALq_8_c(S2,S1)
    implicit none
    integer ipause, mypauses
    type (quaternion_8),INTENT(inOUT)::S2
    type (c_quaternion),INTENT(IN)::S1
    integer i
    type(complextaylor) ct
    call alloc(ct)
      do i=0,3
        ct=s1%x(i)
        s2%x(i)=morph(ct)  !morph(s1%x(i))
      enddo
    call kill(ct)
 end   SUBROUTINE  EQUALq_8_c

  SUBROUTINE  EQUALq_r(S2,S1)
    implicit none
    integer ipause, mypauses
    type (c_quaternion),INTENT(inOUT)::S2
    real(dp),INTENT(IN)::S1
    integer i

      do i=0,3
        s2%x(i)=0.0_dp
      enddo
        s2%x(0)=s1

 end   SUBROUTINE  EQUALq_r


  SUBROUTINE  EQUALq_i(S2,S1)
    implicit none
    integer ipause, mypauses
    type (c_quaternion),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1
    integer i

      do i=0,3
        s2%x(i)=0.0_dp
      enddo
        s2%x(s1)=1

 end   SUBROUTINE  EQUALq_i

  subroutine  quaternion_to_matrix_in_c_damap(p)
    implicit none
    TYPE(c_damap), INTENT(INOUT) :: p
    type(c_quaternion) s,sf
    integer i,j

     call alloc(s)
     call alloc(sf)

    do i=1,3
     s=0.0_dp
     s%x(i)=1.0_dp
     sf=p%q*s*p%q**(-1)
     do j=1,3
      p%s%s(j,i)=sf%x(j)
     enddo
    enddo
     call kill(s)
     call kill(sf)

    end subroutine  quaternion_to_matrix_in_c_damap

  subroutine  q_linear_to_matrix (q_lin,m)
    implicit none
    TYPE(c_spinmatrix), INTENT(INOUT) :: m
    TYPE(q_linear), INTENT(IN) :: q_lin
    type(c_quaternion) s,sf
    type(c_damap) p
    integer i,j

     call alloc(s)
     call alloc(sf)
     call alloc(p)
     p%q=q_lin
    do i=1,3
     s=0.0_dp
     s%x(i)=1.0_dp
     sf=p%q*s*p%q**(-1)
     do j=1,3
      m%s(j,i)=sf%x(j)
     enddo
    enddo
     call kill(s)
     call kill(sf)
     call kill(p)

    end subroutine  q_linear_to_matrix 

  subroutine  q_linear_to_3_by_3_by_6 (q_lin,m)
    implicit none
    real(dp), INTENT(INOUT) :: m(3,3,0:6)
    TYPE(q_linear), INTENT(IN) :: q_lin
    type(q_linear) sf,q,s
    integer i,j

!type q_linear
! complex(dp) mat(6,6)
! complex(dp)  q(0:3,0:6) 
!end type q_linear

     q=1
     q%q=q_lin%q
     s%mat=0
     s%q=0
     m=0
    do i=1,3
     s=i

     sf=q*s*q**(-1)

     do j=1,3
      m(j,i,0:6)=sf%q(j,0:6)
     enddo
    enddo
 

    end subroutine  q_linear_to_3_by_3_by_6 


  FUNCTION cdaddsc( S1, sc )
    implicit none
    TYPE (c_taylor) cdaddsc
    TYPE (c_taylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     cdaddsc%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(cdaddsc)

    ! if(old) then
    call c_dacad(s1%i,sc,c_temp%i)
    call c_dacop(c_temp%i,cdaddsc%i)
    !    else
    !       call newc_dacad(s1%j,sc,daddsc%j)
    !    endif
    c_master=localmaster

  END FUNCTION cdaddsc

  FUNCTION daddsca( S1, sc )
    implicit none
    TYPE (c_taylor) daddsca
    TYPE (c_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    complex(dp) sct
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     daddsca%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(daddsca)
    sct=sc
    ! if(old) then
    call c_dacad(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,daddsca%i)
    !    else
    !       call newc_dacad(s1%j,sc,daddsc%j)
    !    endif
    c_master=localmaster

  END FUNCTION daddsca


  FUNCTION addsc( S1, sc )
    implicit none
    TYPE (c_taylor) addsc
    TYPE (c_taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    integer localmaster
    complex(dp) sct
    IF(.NOT.C_STABLE_DA) then
     addsc%i=0
     RETURN
    endif
    localmaster=c_master

    sct=sc
    if(c_real_warning) call c_real_stop
    !    call check(s1)
    call ass(addsc)


    ! if(old) then
    call c_dacad(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,addsc%i)
    !    else
    !       call newc_dacad(s1%j,REAL(sc,kind=DP),addsc%j)
    !    endif
    c_master=localmaster

  END FUNCTION addsc

  FUNCTION iaddsc( S1, sc )
    implicit none
    TYPE (c_taylor) iaddsc
    TYPE (c_taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    integer localmaster

    complex(dp) sct
    IF(.NOT.C_STABLE_DA) then
     iaddsc%i=0
     RETURN
    endif
    localmaster=c_master

    sct=sc
    !    call check(s1)
    call ass(iaddsc)

    ! if(old) then
    call c_dacad(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,iaddsc%i)
    !    else
    !       call newc_dacad(s1%j,REAL(sc,kind=DP),iaddsc%j)
    !    endif
    c_master=localmaster

  END FUNCTION iaddsc

  FUNCTION cdscadd( sc,S1)
    implicit none
    TYPE (c_taylor) cdscadd
    TYPE (c_taylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     cdscadd%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(cdscadd)

    ! if(old) then

    call c_dacad(s1%i,sc,c_temp%i)
    call c_dacop(c_temp%i,cdscadd%i)
    !    else
    !       call newc_dacad(s1%j,sc,dscadd%j)
    !    endif
    c_master=localmaster

  END FUNCTION cdscadd

  FUNCTION dscadd( sc,S1)
    implicit none
    TYPE (c_taylor) dscadd
    TYPE (c_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    complex(dp) sct
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     dscadd%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(dscadd)

    ! if(old) then
    sct=sc
    call c_dacad(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,dscadd%i)
    !    else
    !       call newc_dacad(s1%j,sc,dscadd%j)
    !    endif
    c_master=localmaster

  END FUNCTION dscadd

  FUNCTION scadd( sc,S1)
    implicit none
    TYPE (c_taylor) scadd
    TYPE (c_taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    complex(dp) sct
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     scadd%i=0
     RETURN
    endif
    localmaster=c_master

    sct=sc
    if(c_real_warning) call c_real_stop
    !    call check(s1)
    call ass(scadd)

    ! if(old) then
    call c_dacad(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,scadd%i)
    !    else
    !       call newc_dacad(s1%j,REAL(sc,kind=DP),scadd%j)
    !    endif

    c_master=localmaster

  END FUNCTION scadd

  FUNCTION iscadd( sc,S1)
    implicit none
    TYPE (c_taylor) iscadd
    TYPE (c_taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    complex(dp) sct
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     iscadd%i=0
     RETURN
    endif
    localmaster=c_master

    sct=sc
    !    call check(s1)
    call ass(iscadd)


    ! if(old) then
    call c_dacad(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,iscadd%i)
    !    else
    !       call newc_dacad(s1%j,REAL(sc,kind=DP),iscadd%j)
    !    endif
    c_master=localmaster

  END FUNCTION iscadd

  FUNCTION subs( S1, S2 )
    implicit none
    TYPE (c_taylor) subs
    TYPE (c_taylor), INTENT (IN) :: S1, S2
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     subs%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    !    call check(s2)
    call ass(subs)



    ! if(old) then
    call c_dasub(s1%i,s2%i,c_temp%i)
    call c_dacop(c_temp%i,subs%i)
    !    else
    !       call newdasub(s1%j,s2%j,subs%j)
    !    endif
    c_master=localmaster

  END FUNCTION subs

  FUNCTION cdsubsc( S1, sc )
    implicit none
    TYPE (c_taylor) cdsubsc
    TYPE (c_taylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     cdsubsc%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(cdsubsc)
    ! if(old) then
    call c_dacsu(s1%i,sc,c_temp%i)
    call c_dacop(c_temp%i,cdsubsc%i)
    !    else
    !       call newdacsu(s1%j,sc,dsubsc%j)
    !    endif

    c_master=localmaster


  END FUNCTION cdsubsc

  FUNCTION dsubsc( S1, sc )
    implicit none
    TYPE (c_taylor) dsubsc
    TYPE (c_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    complex(dp) sct
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     dsubsc%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(dsubsc)
    sct=sc
    ! if(old) then
    call c_dacsu(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,dsubsc%i)
    !    else
    !       call newdacsu(s1%j,sc,dsubsc%j)
    !    endif

    c_master=localmaster


  END FUNCTION dsubsc


  FUNCTION subsc( S1, sc )
    implicit none
    TYPE (c_taylor) subsc
    TYPE (c_taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    complex(dp) sct
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     subsc%i=0
     RETURN
    endif
    localmaster=c_master


    if(c_real_warning) call c_real_stop
    !    call check(s1)
    call ass(subsc)
    sct=sc
    ! if(old) then
    call c_dacsu(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,subsc%i)
    !    else
    !       call newdacsu(s1%j,REAL(sc,kind=DP),subsc%j)
    !    endif
    c_master=localmaster

  END FUNCTION subsc

  FUNCTION isubsc( S1, sc )
    implicit none
    TYPE (c_taylor) isubsc
    TYPE (c_taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    integer localmaster
    complex(dp) x
     IF(.NOT.C_STABLE_DA) then
     isubsc%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(isubsc)
    x=sc
    ! if(old) then
    call c_dacsu(s1%i,x,c_temp%i)
    call c_dacop(c_temp%i,isubsc%i)
    !    else
    !       call newdacsu(s1%j,REAL(sc,kind=DP),isubsc%j)
    !    endif
    c_master=localmaster

  END FUNCTION isubsc

 FUNCTION cdscsub( sc,S1)
    implicit none
    TYPE (c_taylor) cdscsub
    TYPE (c_taylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     cdscsub%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(cdscsub)

    ! if(old) then
    call c_dasuc(s1%i,sc,c_temp%i)
    call c_dacop(c_temp%i,cdscsub%i)
    !    else
    !       call newdasuc(s1%j,sc,dscsub%j)
    !    endif
    c_master=localmaster

  END FUNCTION cdscsub

  FUNCTION dscsub( sc,S1)
    implicit none
    TYPE (c_taylor) dscsub
    TYPE (c_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    complex(dp) sct
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     dscsub%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(dscsub)
    sct=sc
    ! if(old) then
    call c_dasuc(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,dscsub%i)
    !    else
    !       call newdasuc(s1%j,sc,dscsub%j)
    !    endif
    c_master=localmaster

  END FUNCTION dscsub

  FUNCTION scsub( sc,S1)
    implicit none
    TYPE (c_taylor) scsub
    TYPE (c_taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    complex(dp)sct
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     scsub%i=0
     RETURN
    endif

    localmaster=c_master


    if(c_real_warning) call c_real_stop
    !    call check(s1)
    call ass(scsub)
    sct=sc
    ! if(old) then
    call c_dasuc(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,scsub%i)
    !    else
    !       call newdasuc(s1%j,REAL(sc,kind=DP),scsub%j)
    !    endif
    c_master=localmaster

  END FUNCTION scsub

  FUNCTION iscsub( sc,S1)
    implicit none
    TYPE (c_taylor) iscsub
    TYPE (c_taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    integer localmaster
    complex(dp) sct
     IF(.NOT.C_STABLE_DA) then
     iscsub%i=0
     RETURN
    endif
    localmaster=c_master


    !    call check(s1)
    call ass(iscsub)
    sct=sc
    ! if(old) then
    call c_dasuc(s1%i,sct,c_temp%i)
    call c_dacop(c_temp%i,iscsub%i)

    c_master=localmaster

  END FUNCTION iscsub

  !  These are new general TPSA-Routines

  FUNCTION varf( S1, S2 )
    implicit none
    TYPE (c_taylor) varf
    complex(dp), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     varf%i=0
     RETURN
    endif
    localmaster=c_master


    call ass(varf)

    varf=S1 + ((1.0_dp,0.0_dp).cmono.S2)

    c_master=localmaster

  END FUNCTION varf

  FUNCTION varf001( S1, S2 )
    implicit none
    TYPE (c_taylor) varf001
    complex(dp), INTENT (IN) :: S1(2)
    integer  , INTENT (IN) ::  S2
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     varf001%i=0
     RETURN
     endif
    localmaster=c_master


    call ass(varf001)

    varf001=S1(1) + (s1(2).cmono.S2)

    c_master=localmaster

  END FUNCTION varf001



  SUBROUTINE  shift000(S1,S2,s)
    implicit none
    INTEGER,INTENT(IN)::s
    type (c_taylor),INTENT(IN)::S1
    type (c_taylor),INTENT(inout)::S2
    IF(.NOT.C_STABLE_DA) RETURN

    ! if(old) then
    if(s2%i==0) call c_crap1("shift000  1" )  !call etall1(s2%i)
    CALL c_DAshift(s1%i,s2%i,s)
    !   else
    !      if(.NOT. ASSOCIATED(s2%j%r))call c_crap1("shift000  2" )   ! call newetall(s2%j,1)
    !
    !      CALL NEWDAshift(s1%j,s2%j,s)
    !   endif
    !
  END SUBROUTINE shift000


  SUBROUTINE  c_pek000(S1,J,R1)
    implicit none
    INTEGER,INTENT(IN),dimension(:)::j
    complex(dp),INTENT(inOUT)::R1
    type (c_taylor),INTENT(IN)::S1
 !   integer k
    IF(.NOT.C_STABLE_DA) RETURN
    ! if(old) then
    if(s1%i==0) call c_crap1("c_pek000  1" )  !call etall1(s1%i)
 !   k=s1%i
!    write(6,*) r1,k

    CALL c_DApek(s1%i,j,r1)
    !    else
    !       if(.NOT. ASSOCIATED(s1%j%r)) call c_crap1("c_pek000  2" ) ! newetall(s1%j,1)
    !
    !      CALL newDApek(s1%j,j,r1)
    !    endif
    !
  END SUBROUTINE c_pek000

  SUBROUTINE  c_pok000(S1,J,R1)
    implicit none
    INTEGER,INTENT(in),dimension(:)::j
    complex(dp),INTENT(in)::R1
    type (c_taylor),INTENT(inout)::S1
    IF(.NOT.C_STABLE_DA) RETURN

    if(check_j(j)/=0) return
    ! if(old) then
    if(s1%i==0) call c_crap1("c_pok000 1" )  ! call etall1(s1%i)
    CALL c_DApok(s1%i,j,r1)
    !    else
    !       if(.NOT. ASSOCIATED(s1%j%r)) call c_crap1("c_pok000  2" )  ! call newetall(s1%j,1)
    !
    !       CALL newDApok(s1%j,j,r1)
    !    endif
    !
  END SUBROUTINE c_pok000

  SUBROUTINE  c_TAYLOR_ran(S1,r1,R2,f)
    implicit none
    real(dp),INTENT(in)::R1
    real(dp),INTENT(inout)::R2
    type (c_taylor),INTENT(inout)::S1
    logical,optional :: f(:)
    integer j,i
    integer, allocatable:: je(:)
    complex(dp) v
    type (c_taylor) t
    IF(.NOT.C_STABLE_DA) RETURN

    !
    !     THIS SUBROUTINE FILLS THE DA VECTOR A WITH RANDOM ENTRIES.
    !     FOR R1 > 0, THE VECTOR IS FILLED WITH REALS,
    !     FOR R1 < 0, THE VECTOR IS FILLED WITH SINGLE DIGIT INTEGERS
    !     ABS(R1) IS THE FILLING FACTOR
    ! if(old) then
    if(s1%i==0) call c_crap1("c_TAYLOR_ran  1" )  ! call etall1(s1%i)
    call c_daran(s1%i,r1,R2)
        
       if(present(f)) then
       call alloc(t)
        allocate(je(nv))
        j=1
        do while(.true.) 
          call  c_cycle(s1,j,v ,je); if(j==0) exit;
          do i=1,nv
           if(je(i)/=0.and.(.not.f(i))) then
             v=0.0_dp
             exit
           endif
          enddo
         t=t+(v.cmono.je)
        enddo
        s1=t
        deallocate(je)
       call kill(t)
       endif
    !    else
    !       if(.NOT. ASSOCIATED(s1%j%r))call c_crap1("c_TAYLOR_ran  2" )  !  call newetall(s1%j,1)
    !
    !       call newdaran(s1%j,r1,R2)
    !    endif
    !
  END SUBROUTINE c_TAYLOR_ran

  
  SUBROUTINE  c_cfu000(S2,FUN,S1)
    implicit none
    type (c_taylor),INTENT(INOUT)::S1
    type (c_taylor),INTENT(IN)::S2
    complex(dp) FUN
    EXTERNAL FUN
    IF(.NOT.C_STABLE_DA) RETURN
    ! if(old) then
    if(s1%i==0) call c_crap1("c_cfu000  1" )  !  call etall1(s1%i)
    CALL c_DACFU(s2%i,FUN,s1%i)
    !    else
    !       if(.NOT. ASSOCIATED(s1%j%r))call c_crap1("c_cfu000  2" )  !  call newetall(s1%j,1)
    !       CALL NEWDACFU(s2%J,FUN,s1%J)
    !    endif

  END SUBROUTINE c_cfu000

 
  SUBROUTINE  c_taylor_eps(r1)
    implicit none
    real(dp),INTENT(INOUT)::r1
    IF(.NOT.C_STABLE_DA) RETURN
    ! if(old) then
    CALL c_DAeps(r1)
    !   else
    !      CALL newDAeps(r1)
    !   endif

  END SUBROUTINE c_taylor_eps



  FUNCTION GETCHARnd2( S1, S2 )
    implicit none
    TYPE (c_taylor) GETCHARnd2,junk
    TYPE (c_taylor), INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2
    CHARACTER (LEN = LNV)  resul
    integer i,k
    integer localmaster

     IF(.NOT.C_STABLE_DA) then
     GETCHARnd2%i=0
     RETURN
     endif
    localmaster=c_master

   ! ndel=0
    !    call check(s1)
    call ass(GETCHARnd2)

    call alloc(junk)
    resul = trim(ADJUSTR (s2))

    do i=1,lnv
       jfil(i)=0
    enddo

    nd2par= len(trim(ADJUSTR (s2)))

    !frs get around compiler problem
    !frs    do i=1,len(trim(ADJUSTR (s2)))
    do i=1,nd2par
       CALL  CHARINT(RESUL(I:I),Jfil(I))
       if(i>nv) then
          if(Jfil(i)>0) then
             GETCHARnd2=0.0_dp
             return
          endif
       endif

    enddo

    !do i=nd2+ndel+1,nv
    do i=nd2par+1,nv
       if(jfil(i)/=0) then
 
            write(6,*)" error in getchar for .para. "
          ! call !write_e(0)
          stop
       endif
    enddo

    call cfu(s1,filter,junk)

    !DO I=1,ND2+ndel
    DO I=1,ND2par
       DO K=1,JFIL(I)
          JUNK=JUNK.K.I
       ENDDO
    ENDDO

    GETCHARnd2=junk

    call kill(junk)
    c_master=localmaster

  END FUNCTION GETCHARnd2

  FUNCTION GETintnd2( S1, S2 )
    implicit none
    TYPE (c_taylor) GETintnd2,junk
    TYPE (c_taylor), INTENT (IN) :: S1
    integer , INTENT (IN) ::  S2(:)
    integer i,k
    integer localmaster

     IF(.NOT.C_STABLE_DA) then
     GETintnd2%i=0
     RETURN
     endif
    localmaster=c_master

    !    call check(s1)
    call ass(GETintnd2)

    call alloc(junk)

    do i=1,lnv
       jfil(i)=0
    enddo
    nd2par=size(s2)
  !  ndel=0

    !frs get around compiler problem
    !frs    do i=1,len(trim(ADJUSTR (s2)))
    do i=1,nd2par
       Jfil(I)=s2(i)
       if(i>nv) then
          if(Jfil(i)>0) then
             GETintnd2=0.0_dp
             return
          endif
       endif

    enddo

    !do i=nd2+ndel+1,nv
    do i=nd2par+1,nv
       if(jfil(i)/=0) then
 
            write(6,*)" error in GETintnd2 for .para. "
          ! call !write_e(0)
          stop
       endif
    enddo

    call cfu(s1,filter,junk)

    !DO I=1,ND2+ndel
    DO I=1,ND2par
       DO K=1,JFIL(I)
          JUNK=JUNK.K.I
       ENDDO
    ENDDO

    GETintnd2=junk

    call kill(junk)
    c_master=localmaster

  END FUNCTION GETintnd2

  FUNCTION GETintnd2t( S1, S22 )
    implicit none
    TYPE (c_taylor) GETintnd2t,junk
    TYPE (c_taylor), INTENT (IN) :: S1
    type(sub_taylor), INTENT (IN) :: S22
    integer s2(lnv)
    integer i
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     GETintnd2t%i=0
     RETURN
     endif
    localmaster=c_master

    !    call check(s1)
    call ass(GETintnd2t)

    call alloc(junk)

    do i=1,lnv
       jfilt(i)=0
    enddo
    s2=s22%j
    nd2part=s22%min
    nd2partt=s22%max
  !  ndel=0
    !frs get around compiler problem
    !frs    do i=1,len(trim(ADJUSTR (s2)))
    do i=nd2part,nd2partt
       jfilt(I)=s2(i)
       if(i>nv) then
          if(jfilt(i)>0) then
             GETintnd2t=0.0_dp
             return
          endif
       endif

    enddo

    !do i=nd2+ndel+1,nv
    do i=nd2partt+1,nv
       if(jfilt(i)/=0) then
 
            write(6,*)" error in GETintnd2t for .part_taylor. "
          ! call !write_e(0)
          stop
       endif
    enddo

    call cfu(s1,c_filter_part,junk)

    !DO I=1,ND2+ndel
    !    DO I=1,ND2par
    !       DO K=1,jfilt(I)
    !          JUNK=JUNK.K.I
    !       ENDDO
    !    ENDDO

    GETintnd2t=junk

    call kill(junk)
    c_master=localmaster

  END FUNCTION GETintnd2t


  SUBROUTINE  c_taylor_cycle(S1,size,ii,VALUE,J)
    implicit none
    type (c_taylor),INTENT(IN)::S1
    integer,optional, intent(inout):: size
    integer,optional, intent(in):: ii
    integer,optional, intent(inout)::J(:)
    complex(dp), OPTIONAL, intent(inout):: value
    INTEGER ipresent,ILLA
    complex(dp) VALUE0
    IF(.NOT.C_STABLE_DA) RETURN
    ! if(old) THEN
    IF(PRESENT(J).AND.PRESENT(VALUE).and.present(ii)) THEN
       call c_dacycle(S1%i,ii,value,illa,J)
    ELSEif(present(size)) then
       call c_dacycle(S1%i,ipresent,value0,size)
    else
       write(6,*) "error in taylor_cycle"
       stop 888
    ENDIF

  END SUBROUTINE c_taylor_cycle

  SUBROUTINE  c_cycle(S1,i,value,j)
    implicit none
    type (c_taylor),INTENT(in)::S1
     integer i,n
    complex(dp) value
    integer  j(:)
!!!  always initialize i=1  
!!! check for i=0 and then exit your loop    
 
    IF(.NOT.C_STABLE_DA) RETURN

       call c_taylor_cycle(s1,size=n)

       if(i>n) then
        i=0
        return
       endif


       call c_taylor_cycle(s1,ii=i,value=value,j=j)

       i=i+1
 
  END SUBROUTINE c_cycle

  subroutine c_check_snake()
!*
    implicit none
    c_master=c_master+1
    select case (c_master)
    case(1:c_ndumt)
       if(c_iass0user(c_master)>c_scratchda(c_master)%n.or.c_scratchda(c_master)%n>newscheme_max) then
 
            write(6,*) "c_iass0user(c_master),c_scratchda(c_master)%n,newscheme_max"
          write(6,*) c_iass0user(c_master),c_scratchda(c_master)%n,newscheme_max 
          ! call !write_e
          call c_ndum_warning_user
       endif
       c_iass0user(c_master)=0
    case(c_ndumt+1:)
 
       write(6,*) "Should not be here" 
 
       ! call !write_e(101)
    end select
    C_master=C_master-1
  end subroutine c_check_snake

  ! functions used inside other routines

 


  function check_j(j)
    implicit none
    integer check_j
    INTEGER,INTENT(in),dimension(:)::j
    integer i,no
    check_j=0
    IF(.NOT.C_STABLE_DA) RETURN



    no=0
    do i=1,size(j)
       no=j(i)+no
    enddo

    if(no>no) then
       check_j=no
       return
    endif

    do i=nv+1,size(j)
       if(j(i)/=0) then
          check_j=-i
       endif
    enddo
  end function check_j

  function check_harmonic_order(j)
    implicit none
    integer check_harmonic_order
    integer i
    integer,dimension(:)::j

    check_harmonic_order=0
    !do i=1,nd2+ndel
    do i=1,nd2t
       check_harmonic_order=check_harmonic_order+j(i)
    enddo

  end  function check_harmonic_order

  function filter(j)
    implicit none
    complex(dp) filter
    integer i
    integer,dimension(:)::j

    filter=1.0_dp
    !do i=1,nd2+ndel
    do i=1,nd2par
       if(jfil(i)/=j(i)) filter=0.0_dp
    enddo

  end  function filter

  function c_filter_part(j)
    implicit none
    complex(dp) c_filter_part
    integer i
    integer,dimension(:)::j
    !    WRITE(6,*) jfilt(1:4)
    !    WRITE(6,*)nd2part,nd2partt
    c_filter_part=1.0_dp
    !do i=1,nd2+ndel
    do i=nd2part,nd2partt
       if(jfilt(i)/=j(i)) c_filter_part=0.0_dp
    enddo

  end  function c_filter_part

  !  i/o routines

  SUBROUTINE  c_pri_c_ray(S1,MFILE,prec,dospin)
    implicit none
        INTEGER,OPTIONAL,INTENT(IN)::MFILE
    REAL(DP),OPTIONAL,INTENT(INOUT)::prec
    type (c_ray),INTENT(INout)::S1
    logical, optional :: dospin
    integer i,j,k,mfi
    logical(lp) dos
    real(dp) norm

     dos=.true.
     mfi=6
     if(present(mfile)) mfi=mfile

     if(present(dospin)) dos=dospin

    write(mfi,*) "  "
    write(mfi,*) s1%n, " Dimensional ray "
    do i=1,s1%n
     write(mfi,*) s1%x(i)
    enddo
    
if(dos) then
  if(use_quaternion) then
    write(mfi,*) " Quaternion "
    do i=0,3
     write(mfi,*) s1%q%x(i)
    enddo  
    write(mfi,*) " Three spin directions from quaternion "
 write(mfi,*) "S1"
    do i=1,3
     write(mfi,*) s1%S1(i)
    enddo    
 write(mfi,*) "S2"
    do i=1,3
     write(mfi,*) s1%S2(i)
    enddo 
 write(mfi,*) "S3"
    do i=1,3
     write(mfi,*) s1%S3(i)
    enddo
  else
    write(mfi,*) " Three spin directions "
write(mfi,*) "S1"
    do i=1,3
     write(mfi,*) s1%S1(i)
    enddo    
 write(mfi,*) "S2"
    do i=1,3
     write(mfi,*) s1%S2(i)
    enddo 
 write(mfi,*) "S3"
    do i=1,3
     write(mfi,*) s1%S3(i)
    enddo
  endif
else
         write(mfi,*) " Spin results not printed per user's request "
endif

  END SUBROUTINE c_pri_c_ray


  SUBROUTINE  c_pri_map(S1,MFILE,prec,dospin)
    implicit none
        INTEGER,OPTIONAL,INTENT(IN)::MFILE
    REAL(DP),OPTIONAL,INTENT(INOUT)::prec
    type (c_damap),INTENT(INout)::S1
    logical, optional :: dospin
    integer i,j,k,mfi
    logical(lp) rad_in,dos
    real(dp) norm
     dos=.true.
     mfi=6
     if(present(mfile)) mfi=mfile

     if(present(dospin)) dos=dospin

    write(mfi,*) "  "
    if(s1%tpsa) then
     write(mfi,*) s1%n, " Dimensional TPSA map around z=0 "
    else
     write(mfi,*) s1%n, " Dimensional DA map (around chosen orbit in map%x0) "
    endif
    do i=1,s1%n
     call c_pri(s1%v(i),mfile,prec)
    enddo
    
if(dos) then
        call c_full_norm_spin(s1%s,k,norm)
        if(k==-1) then
          write(mfi,*) " Spin Matrix "
          call c_pri_spinmatrix(S1%s,MFILE,prec)  
         endif
        if(k==0) then
         write(mfi,*) " No Spin Matrix "
        endif
        if(k==1) then
         write(mfi,*) " Spin Matrix is identity "
        endif
else
         write(mfi,*) " Spin Matrix is not printed per user's request "
endif
if(dos) then
        call c_full_norm_quaternion(s1%q,k,norm)
        if(k==-1) then
          write(mfi,*) " Quaternion  "
          call c_pri_quaternion(S1%q,MFILE,prec)  
         endif
        if(k==0) then
         write(mfi,*) " No c_quaternion "
        endif
        if(k==1) then
         write(mfi,*) " c_quaternion is identity "
        endif
else
         write(mfi,*) " c_quaternion is not printed per user's request "
endif
            call c_check_rad(s1%e_ij,rad_in)
        if(rad_in) then
         write(mfi,*) "Stochastic Radiation "
          do i=1,6
          do j=1,6
           write(mfi,*) i,j,s1%e_ij(i,j)
          enddo
          enddo
        else
         write(mfi,*) "No Stochastic Radiation "
        endif   

  END SUBROUTINE c_pri_map

  SUBROUTINE  c_pri_quaternion(S2,mf,prec)
    implicit none
    integer ipause, mypauses,i,k
    type (c_quaternion),INTENT(INOUT)::S2
    integer,optional :: mf
    real(dp), optional :: prec
    i=6
    if(present(mf)) i=mf
      write(i,*) " c_quaternion "
    do k=0,3
      call print(s2%x(k),i,prec)
    enddo

  END SUBROUTINE c_pri_quaternion

  SUBROUTINE  c_read_quaternion(S2,mfile)
    implicit none
    integer mfile,i
    type (c_quaternion),INTENT(INOUT)::S2
      character(255) line   
 
     read(mfile,'(a255)') line
    do i=0,3
     call c_rea(s2%x(i),mfile)    
    enddo

  END SUBROUTINE c_read_quaternion


  SUBROUTINE  c_read_map(S1,MFILE)
    implicit none
        INTEGER,OPTIONAL,INTENT(IN)::MFILE
    type (c_damap),INTENT(inout)::S1
    integer i,j,i1,j2
    character(255) line 

    read(mfile,'(a255)') line
    read(mfile,*) s1%n 
    do i=1,s1%n
     call c_rea(s1%v(i),mfile)
    enddo
        read(mfile,'(a255)') line

        if(index(line,"No")/=0) then
           s1%s=0
        elseif(index(line,"id")/=0) then
           s1%s=1
        else
           call c_read_spinmatrix(S1%s,MFILE) 
        endif
  
        read(mfile,'(a255)') line

        if(index(line,"No")/=0) then
           s1%q=0.0_dp
        elseif(index(line,"id")/=0) then
           s1%q=1.0_dp
        else
           call c_read_quaternion(S1%q,MFILE) 
        endif

        read(mfile,'(a255)') line
        if(index(line,"No")/=0) then
         s1%e_ij=0.0_dp
        else
          do i=1,6
          do j=1,6
           read(mfile,*) i1,j2,s1%e_ij(i1,j2)
          enddo
          enddo
        endif   
      
  END SUBROUTINE c_read_map


  SUBROUTINE  c_pri_vec(S1,MFILE,DEPS,dospin)
    implicit none
        INTEGER,OPTIONAL,INTENT(IN)::MFILE
    REAL(DP),OPTIONAL,INTENT(INOUT)::DEPS
    type (c_vector_field),INTENT(INout)::S1
    integer i,mfi,k
    real(dp) norm
    logical, optional :: dospin
    logical(lp) dos
 
     dos=.true.    
    if(present(dospin)) dos=dospin
    
     mfi=6
     if(present(mfile)) mfi=mfile

    write(mfi,*) "  "
    write(mfi,*) s1%n, " Dimensional Vector Field "
    do i=1,s1%n
     call c_pri(s1%v(i),mfile,deps)
    enddo
    
 

    if(dos) then
 
        call c_full_norm_quaternion(s1%q,k,norm)
        if(k==-1) then
          write(mfi,*) " Quaternion  "
          call c_pri_quaternion(S1%q,MFILE,prec=deps)  
         endif
        if(k==0) then
         write(mfi,*) " No c_quaternion "
        endif
        if(k==1) then
         write(mfi,*) " c_quaternion is identity "
        endif
 
    else
         write(mfi,*) " c_quaternion is not printed per user's request "
    endif

 


  
  END SUBROUTINE c_pri_vec

  SUBROUTINE  c_pri_factored_lie(S1,MFILE,DEPS)
    implicit none
        INTEGER,OPTIONAL,INTENT(IN)::MFILE
    REAL(DP),OPTIONAL,INTENT(INOUT)::DEPS
    type (c_factored_lie),INTENT(IN)::S1
    integer i,mfi
     mfi=6
     if(present(mfile)) mfi=mfile

    write(mfi,*) "  "
    if(s1%dir==1) then
     write(mfi,*) " Dragt-Finn Representation ", s1%dir
    else
     write(mfi,*) " Inverse Dragt-Finn Representation ", s1%dir
    endif
    write(mfi,*) s1%n, "  Vector Fields "
    do i=1,s1%n
       write(mfi,*) " Vector field number ",i
     call c_pri_vec(s1%f(i),mfile,deps)
    enddo
    
  END SUBROUTINE c_pri_factored_lie


  SUBROUTINE  c_pri_spinmatrix(S1,MFILE,prec) ! spin routine
    implicit none
        INTEGER,OPTIONAL,INTENT(IN)::MFILE
    REAL(DP),OPTIONAL,INTENT(INOUT)::prec
    type (c_spinmatrix),INTENT(IN)::S1
    integer i,j,mfi
     mfi=6
     if(present(mfile)) mfi=mfile

    do i=1,3
    do j=1,3
     write(mfi,*) " "
     write(mfi,*) i,j
     write(mfi,*) " "
     call c_pri(s1%s(i,j),mfile,prec)
    enddo
    enddo

 
  END SUBROUTINE c_pri_spinmatrix


  SUBROUTINE  c_read_spinmatrix(S1,MFILE) ! spin routine
    implicit none
    INTEGER,INTENT(IN)::MFILE

    type (c_spinmatrix),INTENT(INout)::S1
    integer i,j,i1,j1
    character(255) line

    do i=1,3
    do j=1,3
     read(mfile,'(a255)') line
     read(mfile,*) i1,j1
     read(mfile,'(a255)') line
     call c_rea(s1%s(i1,j1),mfile)
    enddo
    enddo

 
  END SUBROUTINE c_read_spinmatrix


  SUBROUTINE  c_full_norm_spin(S1,k,EPS) ! spin routine
    implicit none
    REAL(DP),OPTIONAL,INTENT(INOUT)::EPS
    real(dp) deps
    type (c_spinmatrix),INTENT(IN)::S1
    integer i,k

      deps=0.0_dp
     k=-1
     call c_full_norm_spinmatrix(s1,deps)

    if(DEPS==0.0_dp) k=0
    if(DEPS==3.0_dp) then
    DEPS=0.0_dp
      do i=1,3
!        DEPS=DEPS+full_abs(s1%s(i,i))
        DEPS=DEPS+(s1%s(i,i))   ! in case of minus 1 
       enddo
     if(DEPS==3.0_dp)  k=1
    endif
    if(present(eps)) eps=deps
  END SUBROUTINE c_full_norm_spin

  SUBROUTINE  c_norm_spin(S1,k,EPS) ! spin routine
    implicit none
    REAL(DP),OPTIONAL,INTENT(INOUT)::EPS
    real(dp) deps
    type (c_spinmatrix),INTENT(IN)::S1
    integer i,k

      deps=0.0_dp
     k=-1
     call c_norm_spinmatrix(s1,deps)

    if(DEPS==0.0_dp) k=0

    if(abs(deps-3.0_dp)<=1.d-7) then
     DEPS=0.0_dp
      do i=1,3
        DEPS=DEPS+abs(s1%s(i,i))
       enddo
     if(abs(deps-3.0_dp)<=1.d-7) then
        k=1
        deps=abs(deps-3.0_dp)
      endif
    endif
    if(present(eps)) eps=deps
  END SUBROUTINE c_norm_spin


  SUBROUTINE  c_pri_spinor(S1,MFILE,prec) ! spin routine
    implicit none
        INTEGER,OPTIONAL,INTENT(IN)::MFILE
    REAL(DP),OPTIONAL,INTENT(INOUT)::prec
    type (c_spinor),INTENT(IN)::S1
    integer i,mfi
     mfi=6
     if(present(mfile)) mfi=mfile

      write(mfi,*) " Complex Spinor "
    do i=1,3
 
     write(mfi,*) " "
     write(mfi,*) i
     write(mfi,*) " "
     call c_pri(s1%v(i),mfile,prec)

    enddo

 
  END SUBROUTINE c_pri_spinor


  SUBROUTINE  c_read_spinor(S1,MFILE) ! spin routine
    implicit none
    INTEGER,INTENT(IN)::MFILE
    type (c_spinor),INTENT(IN)::S1
    integer i
    character(120) line

     ! write(mfile,*) " Complex Spinor "
      read(mfile,'(a120)') line
    do i=1,3
 
      read(mfile,'(a120)') line
      read(mfile,'(a120)') line
      read(mfile,'(a120)') line

     call c_rea(s1%v(i),mfile)

    enddo

 
  END SUBROUTINE c_read_spinor

  SUBROUTINE  c_pri(S1,MFILE,prec)
    implicit none
    INTEGER,OPTIONAL,INTENT(IN)::MFILE
    REAL(DP),OPTIONAL,INTENT(INOUT)::prec
    type (c_taylor),INTENT(IN)::S1
    REAL(DP) deps
    integer mfi
mfi=6
if(present(mfile)) mfi=mfile
    IF(PRESENT(prec)) THEN
       deps=-1.0_dp
       CALL c_taylor_eps(DEPS)
       CALL c_taylor_eps(PREC)
    ENDIF
   
    ! if(old) then
    if(print77) then
       CALL c_DAPRI77(s1%i,mfi)
    else
       CALL c_DAPRI(s1%i,mfi)
    endif

    IF(PRESENT(prec))  CALL c_taylor_eps(deps)

  END SUBROUTINE c_pri

  SUBROUTINE  DAPRINTTAYLORS(S1,MFILE,PREC)
    implicit none
    INTEGER,OPTIONAL,INTENT(IN)::MFILE
    type (C_TAYLOR),INTENT(IN)::S1(:)
    REAL(DP),OPTIONAL,INTENT(INOUT)::PREC
    INTEGER I,mfi
mfi=6
if(present(mfile)) mfi=mfile

    DO I=1,size(S1)
       if(s1(i)%i>0) then
          if(size(S1)>1) write(mfi,*) "Taylor #",i
          CALL C_PRI(s1(i),MFILE,PREC)
       endif
    ENDDO
  END SUBROUTINE DAPRINTTAYLORS

  SUBROUTINE  c_REA(S1,MFILE)
    implicit none
    INTEGER,INTENT(in)::MFILE
    type (c_taylor),INTENT(IN)::S1

    if(s1%i==0)call c_crap1("REA  1" ) !  call etall1(s1%i)

    if(read77) then
       CALL c_darea77(s1%i,MFILE)
    else
       CALL c_DAREA(s1%i,MFILE)
    endif


  END SUBROUTINE c_REA

  SUBROUTINE  DAREADTAYLORS(S1,MFILE)
    implicit none
    INTEGER,INTENT(in)::MFILE
    type (C_TAYLOR),INTENT(INOUT)::S1(NDIM2)
    INTEGER I

    DO I=1,ND2
       CALL C_REA(s1(I),MFILE)
    ENDDO

  END SUBROUTINE DAREADTAYLORS

  ! Universal Taylor Routines   (Sagan's Stuff)



 
 
  ! End of Universal Taylor Routines




  ! Warning Routines

  subroutine c_crap1(STRING)
!*
    implicit none
    CHARACTER(*) STRING

    
     write(6,*) "ERROR IN :"
    write(6,*)  STRING
    ! call !write_e(3478)

  end subroutine c_crap1

  SUBROUTINE c_real_stop()
!*
    implicit none
    integer i(1),j

 
    write(6,*) " You are using a kind(1.0_dp) "
    write(6,*)" set c_real_warning to false to permit this "
    stop 
    i(j)=0
    c_real_warning=.false.

  END   SUBROUTINE c_real_stop


  SUBROUTINE  c_ndum_warning_user()
!*
    implicit none
    integer ipause,II(0:1)


 
    write(6,*) " *  Should never be here in New Linked List Scheme               *"
 
      write(6,*) " do you want a crash? "
    ! call !write_e
     stop
    ii(2000*ipause)=0

  end SUBROUTINE  c_ndum_warning_user

  ! End of  Warning Routines

  ! linked list of da for scratch levels

  SUBROUTINE Set_Up( L ) ! Sets up a layout: gives a unique negative index
!*
    implicit none
    TYPE (c_dalevel) L
    call null_it(L)
    ALLOCATE(L%n);
    ALLOCATE(L%CLOSED);
    L%closed=.FALSE.
    L%N=0
  END SUBROUTINE Set_Up

  SUBROUTINE de_Set_Up( L ) ! deallocates layout content
!*
    implicit none
    TYPE (c_dalevel) L
    deallocate(L%closed);
    deallocate(L%n);
  END SUBROUTINE de_Set_Up



  SUBROUTINE null_it( L ) ! Nullifies layout content
!*
    implicit none
    TYPE (c_dalevel), intent(inout) :: L
    nullify(L%N )
    nullify(L%CLOSED )
    nullify(L%PRESENT )
    !
    nullify(L%END )
    nullify(L%START )
    nullify(L%START_GROUND )! STORE THE GROUNDED VALUE OF START DURING CIRCULAR SCANNING
    nullify(L%END_GROUND )! STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING
  END SUBROUTINE null_it

  SUBROUTINE LINE_L(L,doneit) ! makes into line temporarily
!*
    implicit none
    TYPE (c_DALEVEL) L
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
!*
    implicit none
    TYPE (c_DALEVEL) L
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

  SUBROUTINE APPEND_DA( L ) ! Standard append that clones everything
!*
    implicit none
    TYPE (c_dascratch), POINTER :: Current
    TYPE (c_DALEVEL), TARGET,intent(inout):: L
    logical(lp) doneit
    CALL LINE_L(L,doneit)
    L%N=L%N+1
    nullify(current);ALLOCATE(current);

    call c_alloc_DA(current)

    if(L%N==1) current%next=> L%start
    Current % previous => L % end  ! point it to next fibre
    if(L%N>1)  THEN
       L % end % next => current      !
    ENDIF

    L % end => Current
    if(L%N==1) L%start=> Current
    L%PRESENT=>CURRENT    ! ALWAYS IF APPENDING
    CALL RING_L(L,doneit)
  END SUBROUTINE APPEND_DA

  SUBROUTINE INSERT_DA( L ) ! Standard append that clones everything
!*
    implicit none
    logical(lp) :: doneitt=.true.
    TYPE (c_dascratch), POINTER :: Current
    TYPE (c_DALEVEL), TARGET,intent(inout):: L
    IF(L%N>1.AND.(.NOT.ASSOCIATED(L%PRESENT,L%END))) THEN

       L%N=L%N+1
       nullify(current);ALLOCATE(current);

       call c_alloc_DA(current)

       Current % previous => L % PRESENT   ! 2P -> 2
       Current % NEXT => L % PRESENT%NEXT  ! 2P -> 3
       L%PRESENT%NEXT=> CURRENT            ! 2  -> 2P
       Current % NEXT%PREVIOUS => CURRENT  ! 3  -> 2P
       L%PRESENT=>CURRENT                  ! 2P BECOMES 3

    ELSE

       CALL APPEND_DA( L )
       if(L%N==1) THEN
          L%CLOSED=.TRUE.
          CALL RING_L(L,doneitt)
       ENDIF

    ENDIF
  END SUBROUTINE INSERT_DA

  SUBROUTINE c_alloc_DA( c ) ! Does the full allocation of fibre and initialization of internal variables
!*
    implicit none
    type(c_dascratch),pointer:: c
    ALLOCATE(C%T)
    CALL ALLOC(C%T)
    NULLIFY(C%NEXT)
    NULLIFY(C%PREVIOUS)

  end SUBROUTINE c_alloc_DA

  SUBROUTINE kill_DALEVEL( L )  ! Destroys a layout
!*
    implicit none
    TYPE (c_DASCRATCH), POINTER :: Current
    TYPE (c_DALEVEL) L
    logical(lp) doneit
    CALL LINE_L(L,doneit)
    nullify(current)
    Current => L % end      ! end at the end
    DO WHILE (ASSOCIATED(L % end))
       L % end => Current % previous  ! update the end before disposing
       call dealloc_DASCRATCH(Current)
       Current => L % end     ! alias of last fibre again
       L%N=L%N-1
    END DO
    call de_set_up(L)
  END SUBROUTINE kill_DALEVEL

  SUBROUTINE dealloc_DASCRATCH( c ) ! destroys internal data  if it is not pointing (i.e. not a parent)
!*
    implicit none
    type(c_DASCRATCH),pointer :: c
    IF(ASSOCIATED(C)) THEN
       CALL KILL(C%T)
       IF(ASSOCIATED(C%T)) DEALLOCATE(C%T)
       !       IF(ASSOCIATED(C%NEXT)) DEALLOCATE(C%NEXT)
       !       IF(ASSOCIATED(C%PREVIOUS)) DEALLOCATE(C%PREVIOUS)
       deallocate(c);
    ENDIF
  end SUBROUTINE dealloc_DASCRATCH

  SUBROUTINE set_up_level()
!*
    implicit none
    integer i
    do i=1,c_ndumt
       call set_up(c_scratchda(i))
       !    do j=1,n
       !      call INSERT_da(scratchda(i))
       !    enddo
       !    scratchda(i)%CLOSED=.TRUE.
       !    CALL RING_L(scratchda(i),.TRUE.)
    enddo

  end   SUBROUTINE set_up_level

  SUBROUTINE c_report_level()
!*
    implicit none
    integer i
    if(associated(c_scratchda(1)%n)) then
       do i=1,c_ndumt
 
          write(6,'(a6,1x,i4,a5,1x,i4,1x,a7)') "Level ",i, " has ",c_scratchda(i)%n, "Taylors"
          !          write(w_p%c(1),'(a6,1x,i4,a5,1x,i4,1x,a7)') "Level ",i, " has ",scratchda(i)%n, "Taylors"
          !          ! call !write_e
       enddo
    endif
  END   SUBROUTINE c_report_level

  ! end linked list of da for scratch levels

  ! Assignments Routines

  subroutine c_ASSIGN()
!*
    implicit none
    integer i
    do i=1,c_ndumt
       c_iassdoluser(i)=0
       c_iass0user(i)=0
    enddo
    ! if(old) then
    CALL c_etall1(c_dummy)
    call alloc(c_temp)
    !    else
    !       CALL allocnewda(c_dummyl)
    !       call allocnewda(templ)
    !    endif
    CALL set_up_level
  end subroutine c_ASSIGN

  subroutine c_deassign()
!*
    implicit none
    integer i
    do i=1,c_ndumt
       c_iassdoluser(i)=0
       c_iass0user(i)=0
    enddo
    ! if(old) then
    CALL c_dADAL1(c_dummy)
    call kill(c_temp)
    !    else
    !       CALL KILLnewdaS(c_dummyl)
    !       call KILLnewdaS(templ)
    !    endif
    do i=1,c_ndumt
       CALL kill_DALEVEL(c_scratchda(I))
    ENDDO
  end subroutine c_deassign

  subroutine c_asstaylor(s1)
!*
    implicit none
    TYPE (c_taylor) s1
    !  lastmaster=master  ! 2002.12.13

    select case(c_master)
    case(0:c_ndumt-1)
       c_master=c_master+1
    case(c_ndumt)
       write(6,*) " cannot indent anymore in c_asstaylor ",ndumt
       read(5,*) c_master
       stop 444
    end select
    !    write(26,*) "   taylor ",master
    call c_ass0(s1)
    s1=0.0_dp

  end subroutine c_asstaylor

  subroutine c_ass0(s1)
!*
    implicit none
    integer ipause, mypause
    TYPE (c_taylor) s1

    IF(c_MASTER>c_NDUMT.or.c_master==0) THEN
       WRITE(6,*) "more scratch level needed ",master,NDUMT
       ipause=mypause(123)
       write(6,*) 1/sqrt(-dble(1000+c_master))
       stop 123
    ENDIF

    if(.not.no_ndum_check) c_iass0user(c_master)=c_iass0user(c_master)+1
    if(c_iass0user(c_master)>c_scratchda(c_master)%n) then
       call INSERT_DA( c_scratchda(c_master) )
    ELSE
       c_scratchda(c_master)%PRESENT=>c_scratchda(c_master)%PRESENT%NEXT
    ENDIF
    ! if(old) then
    s1%i=c_scratchda(c_master)%PRESENT%T%i
    !    else
    !       s1%j=scratchda(c_master)%PRESENT%T%j
    !    endif


  end subroutine c_ASS0

!!!!!  assign routines !!!
  subroutine c_assmap(s1)
!*
    implicit none
    TYPE (c_damap) s1
    integer i,j
   
    select case(c_master)
    case(0:c_ndumt-1)
       c_master=c_master+1
    case(c_ndumt)
       write(6,*) " cannot indent anymore in c_assmap ",ndumt
       read(5,*) c_master
       stop 444
    end select

    do i=1,s1%n
       call c_ass0(s1%v(i))
       s1%v(i)=0.0_dp
    enddo

    do i=1,3
    do j=1,3
       call c_ass0(s1%s%s(i,j))
       s1%s%s(i,j)=0.0_dp
    enddo
    enddo

    do i=0,3
       call c_ass0(s1%q%x(i))
       s1%q%x(i)=0.0_dp
    enddo
       s1%x0=0.0_dp
       s1%tpsa=use_tpsa
  end subroutine c_assmap

  subroutine c_ass_quaternion(s1) ! spin routine
!*
    implicit none
    TYPE (c_quaternion) s1
    integer i
   
    select case(c_master)
    case(0:c_ndumt-1)
       c_master=c_master+1
    case(c_ndumt)
       write(6,*) " cannot indent anymore in c_ass_spinmatrix ",c_ndumt
       read(5,*) c_master
       stop 444
    end select



    do i=0,3

       call c_ass0(s1%x(i))
       s1%x(i)=0.0_dp
    enddo
    s1%x(0)=1.0_dp

  end subroutine c_ass_quaternion

  subroutine c_ass_spinmatrix(s1) ! spin routine
!*
    implicit none
    TYPE (c_spinmatrix) s1
    integer i,j
   
    select case(c_master)
    case(0:c_ndumt-1)
       c_master=c_master+1
    case(c_ndumt)
       write(6,*) " cannot indent anymore in c_ass_spinmatrix ",c_ndumt
       read(5,*) c_master
       stop 444
    end select



    do i=1,3
    do j=1,3
       call c_ass0(s1%s(i,j))
       s1%s(i,j)=0.0_dp
    enddo
    enddo

  end subroutine c_ass_spinmatrix


  subroutine c_ass_spinor(s1) ! spin routine
!*
    implicit none
    TYPE (c_spinor) s1
    integer i
   
    select case(c_master)
    case(0:c_ndumt-1)
       c_master=c_master+1
    case(c_ndumt)
       write(6,*) " cannot indent anymore in c_ass_spinmatrix ",c_ndumt
       read(5,*) c_master
       stop 444
    end select



    do i=1,3
       call c_ass0(s1%v(i))
       s1%v(i)=0.0_dp
    enddo

  end subroutine c_ass_spinor

  subroutine c_ass_vector_field(s1)
!*
    implicit none
    TYPE (c_vector_field) s1
    integer i

    select case(c_master)
    case(0:c_ndumt-1)
       c_master=c_master+1
    case(c_ndumt)
       write(6,*) " cannot indent anymore in c_ass_vector_field ",ndumt
       read(5,*) c_master
       stop 444
    end select

    do i=1,s1%n
       call c_ass0(s1%v(i))
       s1%v(i)=0.0_dp
    enddo
 
      do i=0,3
       call c_ass0(s1%q%x(i))
       s1%q%x(i)=0.0_dp
      enddo
    s1%eps   = eps_tpsalie
    s1%nrmax = nrmax

  end subroutine c_ass_vector_field

 SUBROUTINE  c_norm(S1,S2,prec)
    implicit none
    type (c_taylor),INTENT(INOUT)::S2
    type (c_taylor), intent(INOUT):: s1
    real(dp) prec
    INTEGER ipresent,n,I,illa
    complex(dp) value,v
    INTEGER, allocatable :: j(:)
    type (c_taylor) t

    call alloc(t)
    t=0.0_dp
    ipresent=1
    call c_dacycle(S1%I,ipresent,value,n)

    allocate(j(nv))

    do i=1,N
       call c_dacycle(S1%I,i,value,illa,j)
       v=0.0_dp
       if(abs(value)>prec) v=abs(value)
          t=t+(v.cmono.j)
!       endif
    ENDDO
    s2=t
    deallocate(j)
    call kill(t)

  END SUBROUTINE c_norm

  ! remove small numbers

 complex(dp) function  c_clean(s1,prec)
 implicit none
 complex(dp) v,s1
 real(dp) x,prec

  x=s1
  v=0
       if(abs(x)>prec) v=x
       x=aimag(s1)
       if(abs(x)>prec) v=v+i_*x
  c_clean=v

 end function c_clean


  SUBROUTINE  c_clean_yu_w(S1,S2,prec)
    implicit none
    type (c_yu_w),INTENT(INOUT)::S2
    type (c_yu_w), intent(INOUT):: s1
    real(dp) prec
    integer i,j

     do i=1,size(s1%w,1)
     do j=0,s1%n
      call c_clean_taylor(s1%w(i,j),s2%w(i,j),prec)   
     enddo
     enddo
    
end   SUBROUTINE  c_clean_yu_w

  SUBROUTINE  c_clean_taylor(S1,S2,prec,r)
    implicit none
    type (c_taylor),INTENT(INOUT)::S2
    type (c_taylor), intent(INOUT):: s1
    real(dp) prec
    INTEGER ipresent,n,I,illa,k
    complex(dp) value,v,y
    real(dp) x,xx
    INTEGER, allocatable :: j(:)
    type(c_ray),optional :: r
    type (c_taylor) t
    type(c_ray) s
    
    s%x=0
    s%s1=0
    s%s2=0
    s%s3=0
    s%x=1.0_dp
    if(present(r)) s=r
    call alloc(t)
    t=0.0_dp
    ipresent=1
    call c_dacycle(S1%I,ipresent,value,n)

    allocate(j(nv))

    do i=1,N
       call c_dacycle(S1%I,i,value,illa,j)
       v=0.0_dp
    if(present(r)) then
       y=value
      do k=1,nv
       y=y*r%x(k)**j(k)
      enddo
     else
      y=value
    endif
       xx=y
       x=value
       if(abs(xx)>prec) v=x
       xx=aimag(y)
       x=aimag(value)      
       if(abs(xx)>prec) v=v+i_*x
!       if(abs(value)>prec) then
          t=t+(v.cmono.j)
!       endif
    ENDDO
    s2=t
    deallocate(j)
    call kill(t)

  END SUBROUTINE c_clean_taylor


  SUBROUTINE  c_clean_spinmatrix(S1,S2,prec,r) ! spin routine
    implicit none
    type (c_spinmatrix),INTENT(INOUT)::S2
    type (c_spinmatrix), intent(INOUT):: s1
    real(dp) prec
    integer i,j
    type(c_ray),optional :: r
    do i=1,3
    do j=1,3
       call clean(s1%s(i,j),s2%s(i,j),prec,r)
    enddo
    enddo


  END SUBROUTINE c_clean_spinmatrix

  SUBROUTINE  c_clean_quaternion(S1,S2,prec,r) ! spin routine
    implicit none
    type (c_quaternion),INTENT(INOUT)::S2
    type (c_quaternion), intent(INOUT):: s1
    real(dp) prec
    integer i,j
    type(c_ray),optional :: r

    do i=0,3
    
       call clean(s1%x(i),s2%x(i),prec,r)
    
    enddo


  END SUBROUTINE c_clean_quaternion

  SUBROUTINE  c_clean_spinor(S1,S2,prec,r) ! spin routine
    implicit none
    type (c_spinor),INTENT(INOUT)::S2
    type (c_spinor), intent(INOUT):: s1
    real(dp) prec
    type(c_ray),optional :: r
    integer i

    do i=1,3
       call c_clean_taylor(s1%v(i),s2%v(i),prec,r)
    enddo

  END SUBROUTINE c_clean_spinor

  SUBROUTINE  c_clean_damap(S1,S2,prec,r)
    implicit none
    type (c_damap),INTENT(INOUT)::S2
    type (c_damap), intent(INOUT):: s1
    real(dp) prec
    integer i
    type(c_ray),optional :: r

    do i=1,nd2
       call c_clean_taylor(s1%v(i),s2%v(i),prec,r)
    enddo
    
    call c_clean_spinmatrix(s1%s,s2%s,prec)    
    call c_clean_quaternion(s1%q,s2%q,prec)

  END SUBROUTINE c_clean_damap


  SUBROUTINE  c_clean_vector_field(S1,S2,prec,r)
    implicit none
    type (c_vector_field),INTENT(INOUT)::S2
    type (c_vector_field), intent(INOUT):: s1
    real(dp) prec
    integer i
    type(c_ray),optional :: r

    do i=1,nd2
       call c_clean_taylor(s1%v(i),s2%v(i),prec)
    enddo

    do i=0,3
       call c_clean_taylor(s1%q%x(i),s2%q%x(i),prec)
    enddo

  END SUBROUTINE c_clean_vector_field



! typical input if code PTC is calling c_init
!call c_init(c_%NO,c_%nd,c_%np,c_%ndpt,number_of_ac_plane,ptc=my_true)  

  subroutine c_init(NO1,NV1,np1,ndpt1,AC_rf,ptc)  !,spin
    implicit none
    integer, intent(in) :: NO1,NV1
    integer, optional :: np1,ndpt1,AC_RF
    logical(lp), optional :: ptc  !spin,
    integer ndpt_ptc,i

   ! order_gofix=no1
     if(associated(dz_c)) then
      call kill(dz_c)
      deallocate(dz_c)
      nullify(dz_c)
     endif
     call set_da_pointers()

     C_STABLE_DA=.true.
     C_watch_user=.true.
     read77=.true.
     print77=.true.
   if(c_last_tpsa/=0) then
       call c_DEASSIGN
       c_last_tpsa=1
   endif

 ndpt_ptc=0

 RF=0
IF(PRESENT(AC_RF)) RF=AC_RF

!!!! some junk for PTC's bad handling !!!!
if(present(ptc)) then
 if(ptc) then

   if(ndpt1/=0) then
     ndpt_ptc=2*rf
   endif
 endif
endif
!!!! end of  some junk for PTC's bad handling !!!!

  NDPT=0
  ndptb=0
!spin_on=MY_FALSE
if(present(np1)) then  ! This true for any map calculation
 if(present(ndpt1) ) then
  if(ndpt1/=0) then
!REAL FPP moves ndpt to 7 when there is one modulation plane
! now put back to 5 if ndpt_ptc=2 
! This is due to horrible gymnastic in the real FPP used by PTC
   NDPT= NDPT1-ndpt_ptc  
! location of time, ndptb is computed here ndptb=6
   ndptb=ndpt+1
  if(mod(ndpt,2)==0) ndptb=ndpt-1  
  endif
 endif
 NP=np1
 NO=no1
 ND=nv1  ! nv1 is just nd if map used 
 ND2=2*nd !!!!  total dimension of phase space
 nv=nd2+np !!!!  total number of Taylor variables

    ndct=iabs(ndpt-ndptb)  ! 1 if coasting, otherwise 0
    ndc2t=2*ndct  ! 2 if coasting, otherwise 0
    nd2t=nd2-2*rf-ndc2t   !  size of harmonic oscillators minus modulated clocks
    ndt=nd2t/2        ! ndt number of harmonic oscillators minus modulated clocks
    nd2harm=nd2t+2*rf  !!!!  total dimension of harmonic phase space
    ndharm=ndt+rf  !!!! total number of harmonic planes

else
 if(present(ndpt1).or.present(AC_RF).or.present(ptc)) then
   write(6,*) " error : nonsensical input in c_init"
   read(5,*) 
   stop
 endif
 NP=0
 NO=no1
 ND=0
 ND2=0
 NDPT=0
 NV=nv1
endif
!if(present(spin)) spin_on=spin
!write(6,*) "ndc2t,nd2t,nd2harm,nd2"
!write(6,*) ndc2t,nd2t,nd2harm,nd2

      call c_daini(no,nv,0)
    c_master=0  !  master=1   2002.12.2

    CALL c_ASSIGN
    allocate(dz_c(nv))
    call alloc(dz_c)

    do i=1,nv
     dz_c(i)=1.0_dp.cmono.i   
    enddo
! for fast inversion in 
    sj=0
    do i=1,3
     sj(2*i-1,2*i)=1
     sj(2*i,2*i-1)=-1
    enddo 
q_phasor=c_phasor()
qi_phasor=ci_phasor()

c_%rf=>rf
c_%nd2t=>nd2t
c_%nd2harm=>nd2harm
c_%ndc2t=>ndc2t
c_%no=>NO
c_%NDPT=>NDPT
c_%ND=>ND
c_%ND2=>ND2
c_%ndptb=>ndptb
c_%ndpt=>ndpt

c_%pos_of_delta=>pos_of_delta
c_%pos_of_delta=0
if(ndpt/=0) then
c_%pos_of_delta=ndpt
else
 i=nv1-nd2harm-2*rf-np
  if(i/=0) c_%pos_of_delta=nd2harm+1 
endif




!    ndct=iabs(ndpt-ndptb)  ! 1 if coasting, otherwise 0
!    ndc2t=2*ndct  ! 2 if coasting, otherwise 0
!    nd2t=nd2-2*rf-ndc2t   !  size of harmonic oscillators minus modulated clocks
!    ndt=nd2t/2        ! ndt number of harmonic oscillators minus modulated clocks
!    nd2harm=nd2t+2*rf  !!!!  total dimension of harmonic phase space
!    ndharm=ndt+rf  !!!! total number of harmonic planes
!

  end subroutine c_init

  subroutine c_init_all(NO1,NV1,np1,ndpt1,AC_rf,ptc)  !,spin
    implicit none
    integer, intent(in) :: NO1,NV1
    integer, optional :: np1,ndpt1,AC_RF
    logical(lp), optional :: ptc  
    call c_init(NO1,NV1,np1,ndpt1,AC_rf,ptc)
     call init(NO,nd,np,ndpt) 

c_%nd2t=>nd2t
c_%nd2harm=>nd2harm
c_%ndc2t=>ndc2t
c_%no=>NO
c_%NDPT=>NDPT
c_%ND=>ND
c_%ND2=>ND2
c_%ndptb=>ndptb
c_%ndpt=>ndpt
write(6,*) ndpt

 end   subroutine c_init_all

    subroutine c_etcct(x,n1,y,n2,z)
!*
    implicit none
    !  Z=XoY
    integer i,nt,n1,k,n2
    integer,dimension(lnv)::ie,iv
    integer,dimension(:)::x,y,z
    if(.not.c_stable_da) return

    nt=nv-n2
    if(nt.gt.0) then
       do k=1,nt
         call c_etall1(ie(k)) 
       enddo
       do i=n2+1,nv
          call c_davar(ie(i-n2),(0.0_dp,0.0_dp),i)
       enddo
       do i=n2+1,nv
          iv(i)=ie(i-n2)
       enddo
    endif
    do i=1,n2
       iv(i)=y(i)
    enddo
    call c_dacct(x,n1,iv,nv,z,n1)
    if(nt.gt.0) then
       do k=1,nt
         call c_DADAL1(ie(k)) 
       enddo
    endif
    return
  end subroutine c_etcct

 subroutine c_etinv(x,y)
!*
    implicit none
    ! Y=X^-1
    integer i
    type(c_damap) ie1,ie2 
    type(c_damap), intent(inout):: x,y
    if(.not.c_stable_da) return
    
    ie1%n=nv;ie2%n=nv;
    call alloc(ie1)
    call alloc(ie2)
      do i=1,x%n
       ie1%v(i)=x%v(i)
      enddo
      do i=x%n+1,nv
         ie1%v(i)=1.0_dp.cmono.i
      enddo

    call c_dainv(ie1%v(1:nv)%i,nv,ie2%v(1:nv)%i,nv)
       do i=1,x%n
       y%v(i)=ie2%v(i)
      enddo

if(use_quaternion)   THEN
            y%q=x%q*y
      y%q=y%q**(-1)
else
       y%s=x%s*y
      call c_inv_as(y%s,y%s)
endif






    call kill(ie1)

    call kill(ie2)

  end subroutine c_etinv


 FUNCTION transform_vector_field_by_map(S1,S2)
    implicit none
    TYPE (c_vector_field) transform_vector_field_by_map 
    TYPE (c_vector_field), INTENT (IN) :: S1
    TYPE (c_damap), INTENT (IN) :: S2
    TYPE (c_damap) s2i
    integer localmaster
    integer i,k
     IF(.NOT.C_STABLE_DA) then
     transform_vector_field_by_map%v%i=0
     RETURN
     endif

    localmaster=c_master


    transform_vector_field_by_map%n=S1%n

    call c_ass_vector_field(transform_vector_field_by_map)

    call alloc(s2i)
 
     
     s2i=s2**(-1) 

     do k=1,s2%n
      transform_vector_field_by_map%v(k)=0.d0
     do i=1,s2%n
       transform_vector_field_by_map%v(k)=s1%v(i)*(s2i%v(k).d.i)+transform_vector_field_by_map%v(k)
     enddo      
     enddo
 
     do k=1,s2%n
      transform_vector_field_by_map%v(k)=transform_vector_field_by_map%v(k)*s2
     enddo

    transform_vector_field_by_map%nrmax=s1%nrmax
    transform_vector_field_by_map%eps=s1%eps

    call kill(s2i)
    c_master=localmaster
    if(complex_extra_order==1.and.special_extra_order_1) transform_vector_field_by_map=transform_vector_field_by_map.cut.no

  END FUNCTION transform_vector_field_by_map

  FUNCTION c_concat(S1,S2)
    implicit none
    TYPE (c_damap) c_concat,t1,t2,tempnew
    TYPE (c_damap), INTENT (IN) :: S1, S2
    complex(dp) f2it(6,6),f2i(6,6)
    integer i
    logical(lp) rad1,rad2
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_concat%v%i=0
     RETURN
     endif
    localmaster=c_master
    c_concat%n=s1%n
    call c_assmap(c_concat)
    t1%n=s1%n
    t2%n=s2%n
    tempnew%n=t1%n
    call alloc(t1);call alloc(t2);call alloc(tempnew);



    t1=s1;t2=s2;

    do i=1,t1%n
     t1%v(i)=t1%v(i)-(t1%v(i).sub.'0')
    enddo
    do i=1,t2%n
     t2%v(i)=t2%v(i)-(t2%v(i).sub.'0')
    enddo

!    v1=s1     ! change oct 2004.10

    ! if(old) then

    call c_etcct(t1%v%i,t1%n,t2%v%i,t2%n,tempnew%v%i)
    
    if(add_constant_part_concat) then
     do i=1,t1%n
      tempnew%v(i)=tempnew%v(i)+(s1%v(i).sub.'0')
     enddo
    endif

     t1%s = t1%s*t2
     t1%q = t1%q*t2
    
     tempnew%s=t1%s*t2%s
      tempnew%q=t1%q*t2%q

 if(.not.c_similarity) then   
    call c_check_rad(t1%e_ij,rad1)
    call c_check_rad(t2%e_ij,rad2)
   
    if((rad1.or.rad2).and.nd2==6) then
  !   write(6,*) " stochastic "
      t1=t1.sub.1
      f2i=t1    
      f2it=transpose(f2i)
      tempnew%e_ij=t1%e_ij + matmul(matmul(f2i,t2%e_ij),f2it)
    endif
endif
    c_concat=tempnew
    c_concat%x0=s2%x0
    c_concat%tpsa=s2%tpsa
    if(complex_extra_order==1.and.special_extra_order_1) c_concat=c_concat.cut.no


    call kill(t1);call kill(t2);call kill(tempnew);
    c_master=localmaster

  END FUNCTION c_concat



  FUNCTION c_concat_tpsa(s2,s1)
    implicit none
    TYPE (c_damap) c_concat_tpsa,t2,t1,tempnew,t0,t0i
    TYPE (c_damap), INTENT (IN) :: s2, s1
    complex(dp) f2it(6,6),f2i(6,6)

    logical(lp) rad1,rad2
    integer localmaster,i
     IF(.NOT.C_STABLE_DA) then
     c_concat_tpsa%v%i=0
     RETURN
     endif
    localmaster=c_master
    c_concat_tpsa%n=s2%n
    call c_assmap(c_concat_tpsa)
    t0%n=s2%n
    t2%n=s2%n
    t1%n=s1%n
    tempnew%n=t2%n
    call alloc(t2,t1,tempnew);  



    t2=s2;t1=s1;






    call c_etcct(t2%v%i,t2%n,t1%v%i,t1%n,tempnew%v%i)

 
     tempnew%s = t2%s.o.t1
     tempnew%q = t2%q.o.t1
 
    
     tempnew%s=t2%s*t1%s
     tempnew%q=t2%q*t1%q
 
 
 if(.not.c_similarity) then   
    call c_check_rad(t2%e_ij,rad1)
    call c_check_rad(t1%e_ij,rad2)
   
    if((rad1.or.rad2).and.nd2==6) then
  !   write(6,*) " stochastic "
      t2=t2.sub.1
      f2i=t2    
      f2it=transpose(f2i)
      tempnew%e_ij=t2%e_ij + matmul(matmul(f2i,t1%e_ij),f2it)
    endif
endif
    c_concat_tpsa=tempnew
    c_concat_tpsa%x0=s1%x0   
    c_concat_tpsa%tpsa=s1%tpsa

   if(complex_extra_order==1.and.special_extra_order_1) c_concat_tpsa=c_concat_tpsa.cut.no


     call kill(t2,t1,tempnew);  
    c_master=localmaster

  END FUNCTION c_concat_tpsa


  FUNCTION MAKETPSA(s1)
    implicit none
    TYPE (c_damap) MAKETPSA,t0,t1
    TYPE (c_damap), INTENT (IN) :: s1

    integer localmaster,i
     IF(.NOT.C_STABLE_DA) then
     MAKETPSA%v%i=0
     RETURN
     endif
    localmaster=c_master
    MAKETPSA%n=s1%n
    call c_assmap(MAKETPSA)

    t0%n=s1%n
    t1%n=s1%n
    call alloc(t0,t1);  
     t1=s1
     t0=1
      do i=1,s1%n
       t0%v(i)=t0%v(i)-s1%x0(i)
     enddo
     call c_etcct(t1%v%i,t1%n,t0%v%i,t0%n,t1%v%i)
    
 
     t1%s = t1%s.o.t0
     t1%q = t1%q.o.t0
 

 

    MAKETPSA=t1
    MAKETPSA%tpsa=.true.
    MAKETPSA%x0=s1%x0
 
     call kill(t0,t1);  
    c_master=localmaster
   end   FUNCTION MAKETPSA

  FUNCTION MAKEDA(s1)
    implicit none
    TYPE (c_damap) MAKEDA,t0,t1
    TYPE (c_damap), INTENT (IN) :: s1

    integer localmaster,i
     IF(.NOT.C_STABLE_DA) then
     MAKEDA%v%i=0
     RETURN
     endif
    localmaster=c_master
    MAKEDA%n=s1%n
    call c_assmap(MAKEDA)

    t0%n=s1%n
    t1%n=s1%n
    call alloc(t0,t1);  
     t1=s1
     t0=1
      do i=1,s1%n
       t0%v(i)=t0%v(i)+s1%x0(i)
     enddo
     call c_etcct(t1%v%i,t1%n,t0%v%i,t0%n,t1%v%i)
    
     t1%s = t1%s.o.t0
     t1%q = t1%q.o.t0
 
    MAKEDA=t1
    MAKEDA%tpsa=.false.
    MAKEDA%x0=s1%x0
 
     call kill(t0,t1);  
    c_master=localmaster
   end   FUNCTION MAKEDA


  FUNCTION c_adjoint(S1,S2,i)
    implicit none
    TYPE (c_damap) c_adjoint 
    TYPE (c_damap), INTENT (IN) :: S1, S2
    integer , INTENT (IN) :: i
    logical(lp) rad1
    integer localmaster
    complex(dp) f2it(6,6),f2i(6,6)
     IF(.NOT.C_STABLE_DA) then
     c_adjoint%v%i=0
     RETURN
     endif
     
    localmaster=c_master

    c_similarity=my_true
    c_adjoint%n=s1%n
    call c_assmap(c_adjoint)
 
    if(i==1) then
     c_adjoint=s1*s2*s1**(-1)    
    else
     c_adjoint=s1**(-1)*s2*s1
    endif
    
     if(c_similarity) then   
    call c_check_rad(s2%e_ij,rad1)
 
   
    if(rad1.and.nd2==6) then
   !  write(6,*) " stochastic "
      if(i==1) then
        f2i=s1   
      else
        f2i=s1**(-1)
      endif
      f2it=transpose(f2i)
      c_adjoint%e_ij= matmul(matmul(f2i,s2%e_ij),f2it)
    endif

    endif
 
       if(complex_extra_order==1.and.special_extra_order_1) c_adjoint=c_adjoint.cut.no
    c_similarity=my_false
  
    c_master=localmaster

  END FUNCTION c_adjoint

  FUNCTION c_adjoint_vec(S1,S2,i)
    implicit none
    TYPE (c_damap) c_adjoint_vec 
    TYPE (c_vector_field), INTENT (IN) :: S1
    TYPE (c_damap), INTENT (IN) :: S2
    integer , INTENT (IN) :: i
    logical(lp) rad1
    integer localmaster
    complex(dp) f2it(6,6),f2i(6,6)
     IF(.NOT.C_STABLE_DA) then
     c_adjoint_vec%v%i=0
     RETURN
     endif
     
    localmaster=c_master

    c_similarity=my_true
    
    c_adjoint_vec%n=s1%n
    call c_assmap(c_adjoint_vec)

    if(i==1) then
     c_adjoint_vec=exp(s1)*s2
     c_adjoint_vec=exp(-s1,c_adjoint_vec)
    else
     c_adjoint_vec=exp(-s1)*s2
     c_adjoint_vec=exp(s1,c_adjoint_vec)
    endif
    
     if(c_similarity) then   
    call c_check_rad(s2%e_ij,rad1)
 
   
    if(rad1.and.nd2==6) then
  !   write(6,*) " stochastic "
      if(i==1) then
        f2i=exp(s1)   
      else
        f2i=exp(-s1)
      endif
      f2it=transpose(f2i)
      c_adjoint_vec%e_ij= matmul(matmul(f2i,s2%e_ij),f2it)
    endif

    endif
    
       if(complex_extra_order==1.and.special_extra_order_1) c_adjoint_vec=c_adjoint_vec.cut.no

    c_similarity=my_false
    
    c_master=localmaster

  END FUNCTION c_adjoint_vec




  FUNCTION c_spinmatrix_spinmatrix(S1,S2) ! spin routine function
    implicit none
    TYPE (c_spinmatrix) c_spinmatrix_spinmatrix 
    TYPE (c_spinmatrix), INTENT (IN) :: S1, S2
    integer i,j,k

    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_spinmatrix_spinmatrix%s%i=0
     RETURN
     endif
    localmaster=c_master

    call c_ass_spinmatrix(c_spinmatrix_spinmatrix)

 
    c_spinmatrix_spinmatrix=0

     do i=1,3
      do j=1,3
       do k=1,3
          c_spinmatrix_spinmatrix%s(i,k)=s1%s(i,j)*s2%s(j,k)+c_spinmatrix_spinmatrix%s(i,k)
       enddo
      enddo
     enddo
     
        if(complex_extra_order==1.and.special_extra_order_1) c_spinmatrix_spinmatrix=c_spinmatrix_spinmatrix.cut.no

    c_master=localmaster

  END FUNCTION c_spinmatrix_spinmatrix

  FUNCTION c_spinmatrix_mul_cray(S1,S2) ! spin routine function
    implicit none
    TYPE (c_ray) c_spinmatrix_mul_cray 
    TYPE (c_spinmatrix), INTENT (IN) :: S1
    TYPE (c_ray), INTENT (IN) :: S2 
    TYPE (c_spinmatrix) temp
    integer i,j

    call alloc(temp)
    temp=s1
    c_spinmatrix_mul_cray=s2
    temp=c_concat_spinmatrix_ray(temp,s2) !.o.s2
    c_spinmatrix_mul_cray%s1=0.0_dp
    c_spinmatrix_mul_cray%s2=0.0_dp
    c_spinmatrix_mul_cray%s3=0.0_dp

     do i=1,3
      do j=1,3
 
          c_spinmatrix_mul_cray%s1(i)=temp%s(i,j)*s2%s1(j)+c_spinmatrix_mul_cray%s1(i)
          c_spinmatrix_mul_cray%s2(i)=temp%s(i,j)*s2%s2(j)+c_spinmatrix_mul_cray%s2(i)
          c_spinmatrix_mul_cray%s3(i)=temp%s(i,j)*s2%s3(j)+c_spinmatrix_mul_cray%s3(i)
 
      enddo
     enddo
!     c_spinmatrix_mul_cray%x=s2%x
 
     call kill(temp)
  END FUNCTION c_spinmatrix_mul_cray

  FUNCTION c_quaternion_mul_cray(S1,S2) ! spin routine function
    implicit none
    TYPE (c_ray) c_quaternion_mul_cray 
    TYPE (c_quaternion), INTENT (IN) :: S1
    TYPE (c_ray), INTENT (IN) :: S2 
    TYPE (c_quaternion) temp
    TYPE (complex_quaternion) tempc,mc
    integer i,j

      
   
    call alloc(temp)
    mc%x=0.0_dp
    temp=s1
    c_quaternion_mul_cray=s2
    tempc=c_concat_quaternion_ray(temp,s2)  !.o.s2
    c_quaternion_mul_cray%q=tempc
         mc%x=0.0_dp
         mc%x(1:3)=s2%s1
         mc=tempc*mc*tempc**(-1)
         c_quaternion_mul_cray%s1=mc%x(1:3) 
         mc%x=0.0_dp
         mc%x(1:3)=s2%s2
         mc=tempc*mc*tempc**(-1)
         c_quaternion_mul_cray%s2=mc%x(1:3) 
         mc%x=0.0_dp
         mc%x(1:3)=s2%s3
         mc=tempc*mc*tempc**(-1)
         c_quaternion_mul_cray%s3=mc%x(1:3) 
 !   c_quaternion_mul_cray%q=tempc*s2%q*tempc**(-1)
!    c_quaternion_mul_cray%x=s2%x
 
     call kill(temp)
  END FUNCTION c_quaternion_mul_cray

 FUNCTION c_spinmatrix_spinor(S1,S2) ! spin routine function
    implicit none
    TYPE (c_spinor) c_spinmatrix_spinor 
    TYPE (c_spinmatrix), INTENT (IN) :: S1
    TYPE (c_spinor), INTENT (IN) :: S2
    integer localmaster
    integer i,j
     IF(.NOT.C_STABLE_DA) then
     c_spinmatrix_spinor%v%i=0
     RETURN
     endif
    localmaster=c_master

    call c_ass_spinor(c_spinmatrix_spinor)

 
    c_spinmatrix_spinor=0

     do i=1,3
      do j=1,3
          c_spinmatrix_spinor%v(i)=s1%s(i,j)*s2%v(j)+c_spinmatrix_spinor%v(i)
      enddo
     enddo
     
         if(complex_extra_order==1.and.special_extra_order_1) c_spinmatrix_spinor=c_spinmatrix_spinor.cut.no

    c_master=localmaster

  END FUNCTION c_spinmatrix_spinor

  FUNCTION c_transpose(S1)
    implicit none
    TYPE (c_damap) c_transpose,t1
    TYPE (c_damap), INTENT (IN) :: S1
    complex(dp) f2t(ndim2t,ndim2t)

    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_transpose%v%i=0
     RETURN
     endif

    localmaster=c_master

    c_transpose%n=s1%n
    call c_assmap(c_transpose)
 
    t1%n=s1%n

    call alloc(t1);

    t1=s1.sub.1
    
    f2t=t1
    f2t=transpose(f2t)

    c_transpose=f2t


    call kill(t1);

    c_master=localmaster

  END FUNCTION c_transpose




 FUNCTION c_spinor_cmap(S1,S2) ! spin routine function
    implicit none
    TYPE (c_spinor) c_spinor_cmap 
    TYPE (c_damap), INTENT (IN) :: S2
    TYPE (c_spinor), INTENT (IN) :: S1
    integer localmaster
    integer i
     IF(.NOT.C_STABLE_DA) then
     c_spinor_cmap%v%i=0
     RETURN
     endif

    localmaster=c_master

    call c_ass_spinor(c_spinor_cmap)

 
    c_spinor_cmap=0

     do i=1,3
          c_spinor_cmap%v(i)=s1%v(i)*s2
     enddo
     
          if(complex_extra_order==1.and.special_extra_order_1) c_spinor_cmap=c_spinor_cmap.cut.no
    c_master=localmaster

  END FUNCTION c_spinor_cmap

 FUNCTION c_spinor_cmap_tpsa(S1,S2) ! spin routine function
    implicit none
    TYPE (c_spinor) c_spinor_cmap_tpsa 
    TYPE (c_damap), INTENT (IN) :: S2
    TYPE (c_spinor), INTENT (IN) :: S1
    integer localmaster
    integer i
     IF(.NOT.C_STABLE_DA) then
     c_spinor_cmap_tpsa%v%i=0
     RETURN
     endif

    localmaster=c_master

    call c_ass_spinor(c_spinor_cmap_tpsa)

 
    c_spinor_cmap_tpsa=0

     do i=1,3
          c_spinor_cmap_tpsa%v(i)=s1%v(i)*s2
     enddo
     
          if(complex_extra_order==1.and.special_extra_order_1) c_spinor_cmap_tpsa=c_spinor_cmap_tpsa.cut.no
    c_master=localmaster

  END FUNCTION c_spinor_cmap_tpsa

  FUNCTION c_complex_spinmatrix(S1,S2) ! spin routine function
    implicit none
    TYPE (c_spinmatrix) c_complex_spinmatrix 
    TYPE (c_spinmatrix), INTENT (IN) :: S2
    complex (dp), INTENT (IN) :: S1
    integer i,j

    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_complex_spinmatrix%s%i=0
     RETURN
     endif
    localmaster=c_master

    call c_ass_spinmatrix(c_complex_spinmatrix)

 
    c_complex_spinmatrix=0

     do i=1,3
      do j=1,3
          c_complex_spinmatrix%s(i,j)=s1*s2%s(i,j) 
      enddo
     enddo
     
 
    c_master=localmaster

  END FUNCTION c_complex_spinmatrix

  FUNCTION c_spinmatrix_add_spinmatrix(S1,S2) ! spin routine function
    implicit none
    TYPE (c_spinmatrix) c_spinmatrix_add_spinmatrix 
    TYPE (c_spinmatrix), INTENT (IN) :: S1, S2
    integer i,j

    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_spinmatrix_add_spinmatrix%s%i=0
     RETURN
     endif
    localmaster=c_master

    call c_ass_spinmatrix(c_spinmatrix_add_spinmatrix)

 
    c_spinmatrix_add_spinmatrix=0

     do i=1,3
      do j=1,3

          c_spinmatrix_add_spinmatrix%s(i,j)=s1%s(i,j)+s2%s(i,j) 

      enddo
     enddo
     
 
    c_master=localmaster

  END FUNCTION c_spinmatrix_add_spinmatrix

 FUNCTION c_spinmatrix_sub_spinmatrix(S1,S2) ! spin routine function
    implicit none
    TYPE (c_spinmatrix) c_spinmatrix_sub_spinmatrix 
    TYPE (c_spinmatrix), INTENT (IN) :: S1, S2
    integer i,j

    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_spinmatrix_sub_spinmatrix%s%i=0
     RETURN
     endif
    localmaster=c_master

    call c_ass_spinmatrix(c_spinmatrix_sub_spinmatrix)

 
    c_spinmatrix_sub_spinmatrix=0

     do i=1,3
      do j=1,3

          c_spinmatrix_sub_spinmatrix%s(i,j)=s1%s(i,j)-s2%s(i,j) 

      enddo
     enddo
     
 
    c_master=localmaster

  END FUNCTION c_spinmatrix_sub_spinmatrix



  FUNCTION c_spinor_add_spinor(S1,S2) ! spin routine function
    implicit none
    TYPE (c_spinor) c_spinor_add_spinor 
    TYPE (c_spinor), INTENT (IN) :: S1, S2
    integer i

    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_spinor_add_spinor%v%i=0
     RETURN
     endif
    localmaster=c_master

    call c_ass_spinor(c_spinor_add_spinor)

     do i=1,3
          c_spinor_add_spinor%v(i)=s1%v(i)+s2%v(i) 
     enddo
     
    c_master=localmaster

  END FUNCTION c_spinor_add_spinor


  FUNCTION c_spinor_sub_spinor(S1,S2) ! spin routine function
    implicit none
    TYPE (c_spinor) c_spinor_sub_spinor 
    TYPE (c_spinor), INTENT (IN) :: S1, S2
    integer i

    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_spinor_sub_spinor%v%i=0
     RETURN
     endif
    localmaster=c_master

    call c_ass_spinor(c_spinor_sub_spinor)

     do i=1,3
          c_spinor_sub_spinor%v(i)=s1%v(i)-s2%v(i) 
     enddo
     
    c_master=localmaster

  END FUNCTION c_spinor_sub_spinor



 FUNCTION c_taylor_spinor(S1,S2) ! spin routine function
    implicit none
    TYPE (c_spinor) c_taylor_spinor 
    TYPE (c_spinor), INTENT (IN) ::  S2
    TYPE (c_taylor), INTENT (IN) :: S1
    integer i

    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_taylor_spinor%v%i=0
     RETURN
     endif
    localmaster=c_master

    call c_ass_spinor(c_taylor_spinor)

     do i=1,3
          c_taylor_spinor%v(i)=s1*s2%v(i) 
     enddo
     
 
    c_master=localmaster

  END FUNCTION c_taylor_spinor

 FUNCTION c_complex_spinor(S1,S2) ! spin routine function
    implicit none
    TYPE (c_spinor) c_complex_spinor 
    TYPE (c_spinor), INTENT (IN) ::  S2
    complex(dp), INTENT (IN) :: S1
    integer i

    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_complex_spinor%v%i=0
     RETURN
     endif
    localmaster=c_master

    call c_ass_spinor(c_complex_spinor)

     do i=1,3
          c_complex_spinor%v(i)=s1*s2%v(i) 
     enddo
     
 
    c_master=localmaster

  END FUNCTION c_complex_spinor

 FUNCTION c_real_spinor(S1,S2) ! spin routine function
    implicit none
    TYPE (c_spinor) c_real_spinor 
    TYPE (c_spinor), INTENT (IN) ::  S2
    real(dp), INTENT (IN) :: S1
    integer i

    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_real_spinor%v%i=0
     RETURN
     endif
    localmaster=c_master

    call c_ass_spinor(c_real_spinor)

     do i=1,3
          c_real_spinor%v(i)=s1*s2%v(i) 
     enddo
    
    c_master=localmaster

  END FUNCTION c_real_spinor

 FUNCTION c_spinor_spinor(S1,S2) ! spin routine function
    implicit none
    TYPE (c_spinor) c_spinor_spinor 
    TYPE (c_spinor), INTENT (IN) ::  S2
    TYPE (c_spinor), INTENT (IN) :: S1


    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_spinor_spinor%v%i=0
     RETURN
     endif
    localmaster=c_master

    call c_ass_spinor(c_spinor_spinor)

          c_spinor_spinor%v(1)=s1%v(2)*s2%v(3)-s1%v(3)*s2%v(2)
          c_spinor_spinor%v(2)=s1%v(3)*s2%v(1)-s1%v(1)*s2%v(3)
          c_spinor_spinor%v(3)=s1%v(1)*s2%v(2)-s1%v(2)*s2%v(1)
           if(complex_extra_order==1.and.special_extra_order_1) c_spinor_spinor=c_spinor_spinor.cut.no

    c_master=localmaster

  END FUNCTION c_spinor_spinor

  FUNCTION c_trxspinmatrix( S1, S2 ) ! spin routine function
    implicit none
    TYPE (c_spinmatrix) c_trxspinmatrix
    TYPE (c_spinmatrix), INTENT (IN) :: S1
    TYPE (c_damap), INTENT (IN) ::  S2
    integer i,j
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_trxspinmatrix%s%i=0
     RETURN
     endif
    localmaster=c_master

    call c_ass_spinmatrix(c_trxspinmatrix)
     
      do i=1,3
       do j=1,3
        c_trxspinmatrix%s(i,j)=s1%s(i,j)*s2
       enddo
      enddo
     if(complex_extra_order==1.and.special_extra_order_1) c_trxspinmatrix=c_trxspinmatrix.cut.no
    c_master=localmaster

  END FUNCTION c_trxspinmatrix

  FUNCTION c_trxquaternion( S1, S2 ) ! spin routine function
    implicit none
    TYPE (c_quaternion) c_trxquaternion
    TYPE (c_quaternion), INTENT (IN) :: S1
    TYPE (c_damap), INTENT (IN) ::  S2
    integer i
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_trxquaternion%x(1)%i=0
     RETURN
     endif
    localmaster=c_master

    call c_ass_quaternion(c_trxquaternion)
     
      do i=0,3

        c_trxquaternion%x(i)=s1%x(i)*s2
       enddo
      
     if(complex_extra_order==1.and.special_extra_order_1) c_trxquaternion=c_trxquaternion.cut.no
    c_master=localmaster

  END FUNCTION c_trxquaternion

  FUNCTION c_trxquaternion_tpsa( S1, S2 ) ! spin routine function
    implicit none
    TYPE (c_quaternion) c_trxquaternion_tpsa
    TYPE (c_quaternion), INTENT (IN) :: S1
    TYPE (c_damap), INTENT (IN) ::  S2
    integer i
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_trxquaternion_tpsa%x(1)%i=0
     RETURN
     endif
    localmaster=c_master

    call c_ass_quaternion(c_trxquaternion_tpsa)
     
      do i=0,3

        c_trxquaternion_tpsa%x(i)=s1%x(i).o.s2
       enddo
      
     if(complex_extra_order==1.and.special_extra_order_1) c_trxquaternion_tpsa=c_trxquaternion_tpsa.cut.no
    c_master=localmaster

  END FUNCTION c_trxquaternion_tpsa


  FUNCTION c_trxspinmatrixda( S1, S2 ) ! spin routine function
    implicit none
    TYPE (c_spinmatrix) c_trxspinmatrixda
    TYPE (c_spinmatrix), INTENT (IN) :: S1
    TYPE (c_damap), INTENT (IN) ::  S2
    integer i,j
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_trxspinmatrixda%s%i=0
     RETURN
     endif
    localmaster=c_master

    call c_ass_spinmatrix(c_trxspinmatrixda)
     
      do i=1,3
       do j=1,3
        c_trxspinmatrixda%s(i,j)=s1%s(i,j).o.s2
       enddo
      enddo
         if(complex_extra_order==1.and.special_extra_order_1) c_trxspinmatrixda=c_trxspinmatrixda.cut.no
    c_master=localmaster

  END FUNCTION c_trxspinmatrixda



  FUNCTION c_trxtaylor( S1, S2 )
    implicit none
    TYPE (c_taylor) c_trxtaylor
    TYPE (c_taylor), INTENT (IN) :: S1
    TYPE (c_damap), INTENT (IN) ::  S2
    TYPE (c_damap)  S22,temp
    integer i
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_trxtaylor%i=0
     RETURN
     endif
    localmaster=c_master



 
    call c_ASStaylor(c_trxtaylor)

    s22%n=s2%n
    call alloc(s22)
    temp%n=s2%n
    call alloc(temp)

    s22=s2
    do i=1,s22%n
     s22%v(i)=s22%v(i)-(s22%v(i).sub.'0')
    enddo

    temp%v(1)=s1
        call c_etcct(temp%v%i,temp%n,s22%v%i,s22%n,temp%v%i)
   ! temp=temp*s22  !(to prevent a recursive call in c_concat)

    c_trxtaylor=temp%v(1)

    call kill(temp)
    call kill(s22)
    c_master=localmaster

  END FUNCTION c_trxtaylor


  FUNCTION c_trxtaylor_da( S1, S2 )
    implicit none
    TYPE (c_taylor) c_trxtaylor_da
    TYPE (c_taylor), INTENT (IN) :: S1
    TYPE (c_damap), INTENT (IN) ::  S2
    TYPE (c_damap)  S22,temp

    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_trxtaylor_da%i=0
     RETURN
     endif
    localmaster=c_master



 
    call c_ASStaylor(c_trxtaylor_da)

    s22%n=s2%n
    call alloc(s22)
    temp%n=s2%n
    call alloc(temp)

    s22=s2
   ! do i=1,s22%n
   !  s22%v(i)=s22%v(i)-(s22%v(i).sub.'0')
   ! enddo

    temp%v(1)=s1
        call c_etcct(temp%v%i,temp%n,s22%v%i,s22%n,temp%v%i)
   ! temp=temp*s22  !(to prevent a recursive call in c_concat)

    c_trxtaylor_da=temp%v(1)

    call kill(temp)
    call kill(s22)
    c_master=localmaster

  END FUNCTION c_trxtaylor_da


 FUNCTION c_concat_spinor_ray( S1, S2 )
    implicit none
    TYPE (c_spinor) c_concat_spinor_ray
    TYPE (c_spinor), INTENT (IN) :: S1
    TYPE (c_ray), INTENT (IN) ::  S2
    integer i
    integer localmaster
    
     IF(.NOT.C_STABLE_DA) then
     c_concat_spinor_ray=0
     RETURN
     endif

    localmaster=c_master
     
    call c_ass_spinor(c_concat_spinor_ray)

    do i=1,3
     c_concat_spinor_ray%v(i)=s1%v(i).o.s2
    enddo
     
    c_master=localmaster

  END FUNCTION c_concat_spinor_ray


 FUNCTION c_concat_spinmatrix_ray( S1, S2 )
    implicit none
    TYPE (c_spinmatrix) c_concat_spinmatrix_ray
    TYPE (c_spinmatrix), INTENT (IN) :: S1
    TYPE (c_ray), INTENT (IN) ::  S2
    integer i,j
    integer localmaster
    
     IF(.NOT.C_STABLE_DA) then
     c_concat_spinmatrix_ray=1
     RETURN
     endif

    localmaster=c_master
     
    call c_ass_spinmatrix(c_concat_spinmatrix_ray)

    do i=1,3
    do j=1,3
     c_concat_spinmatrix_ray%s(i,j)=s1%s(i,j).o.s2
    enddo
    enddo
     
    c_master=localmaster

  END FUNCTION c_concat_spinmatrix_ray

 FUNCTION c_concat_quaternion_ray( S1, S2 )
    implicit none
    TYPE (c_quaternion) c_concat_quaternion_ray
    TYPE (c_quaternion), INTENT (IN) :: S1
    TYPE (c_ray), INTENT (IN) ::  S2
    integer i,j
    integer localmaster
    
     IF(.NOT.C_STABLE_DA) then
     c_concat_quaternion_ray=1.0_dp
     RETURN
     endif

    localmaster=c_master
     
    call c_ass_quaternion(c_concat_quaternion_ray)

    do i=0,3
      c_concat_quaternion_ray%x(i)=s1%x(i).o.s2
     enddo
    c_master=localmaster

  END FUNCTION c_concat_quaternion_ray



  FUNCTION c_concat_c_ray( S1, S2 )
    implicit none
    complex(dp) c_concat_c_ray
    TYPE (c_taylor), INTENT (IN) :: S1
    TYPE (c_ray), INTENT (IN) ::  S2
    TYPE (c_damap)  temp
    integer i
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_concat_c_ray=0
     RETURN
     endif
    localmaster=c_master

    c_concat_c_ray=0.0_dp

 
    temp%n=nv
    call alloc(temp)

    do i=1,nv
     temp%v(i)=s2%x(i)  !-s2%x0(i)
    enddo
     

      c_concat_c_ray=s1.o.temp


    call kill(temp)
    c_master=localmaster

  END FUNCTION c_concat_c_ray


  FUNCTION c_concat_map_ray( S1, S2 )
    implicit none
    TYPE (c_ray) c_concat_map_ray
    TYPE (c_damap), INTENT (IN) :: S1
    TYPE (c_ray), INTENT (IN) ::  S2
    TYPE (c_quaternion) q

    type(c_ray) temp2
    real(dp) norm
    integer i
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_concat_map_ray%x=0
     RETURN
     endif
    localmaster=c_master
    c_concat_map_ray=s2
  !  c_concat_map_ray%x=0.0_dp
  !  c_concat_map_ray%n=s2%n

    temp2=0
    temp2%n=s2%n
    temp2=s2
    temp2%x0=0.0_dp
     norm=0
    do i=1,s1%n
     norm=norm+abs(s2%x0(i))
    enddo

   ! if(norm>eps_tpsalie.and.S1%tpsa) then
   !  write(6,*) "Both c_ray and c_damap are tpsa: not allowed "
   !  stop 997
   ! endif

     if(S1%tpsa) then
      do i=1,s1%n
        temp2%x(i)=s2%x(i)  !-s1%x0(i)
       enddo
     else     
      do i=1,s1%n
        temp2%x(i)=s2%x(i)-s2%x0(i)
       enddo
    endif
 
      if(use_quaternion) then
       c_concat_map_ray=s1%q.o.temp2    !s2
      else
       c_concat_map_ray=s1%s.o.temp2    !s2
      endif
 ! order important because the first one puts s2%x into c_concat_map_ray%x
     do i=1,s1%n
      c_concat_map_ray%x(i)=s1%v(i).o.temp2
     enddo

    c_master=localmaster

  END FUNCTION c_concat_map_ray

  FUNCTION c_concat_vector_field_ray( S1, S2 )
    implicit none
    TYPE (c_ray) c_concat_vector_field_ray
    TYPE (c_vector_field), INTENT (IN) :: S1
    TYPE (c_ray), INTENT (IN) ::  S2
 
 
    integer i 
 
     IF(.NOT.C_STABLE_DA) then
     c_concat_vector_field_ray%x=0
     RETURN
     endif
 

    c_concat_vector_field_ray%x=0.0_dp

     do i=1,s1%n
      c_concat_vector_field_ray%x(i)=s1%v(i).o.s2
     enddo

 
  END FUNCTION c_concat_vector_field_ray

  FUNCTION c_bra_v_ct( S1, S2 )
    implicit none
    TYPE (c_taylor) c_bra_v_ct
    TYPE (c_vector_field), INTENT (IN) :: S1
    TYPE (c_taylor), INTENT (IN) ::  S2
    TYPE (c_taylor)  S22 
 
    integer i
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_bra_v_ct%i=0 
     RETURN
     endif
    localmaster=c_master

 
    call c_ASStaylor(c_bra_v_ct)
 
    call alloc(s22 )
    

    s22=(0.0_dp,0.0_dp)

     do i=1,s1%n
 
     s22=s22 + s1%v(i)*(s2.d.i)
 
    enddo
 
    c_bra_v_ct=s22
    
    call kill(s22  )
    c_master=localmaster

  END FUNCTION c_bra_v_ct

  
    FUNCTION c_bra_v_q( S1, S2 )
    implicit none
    TYPE (c_quaternion) c_bra_v_q
    TYPE (c_vector_field), INTENT (IN) :: S1
    TYPE (c_quaternion), INTENT (IN) ::  S2
    TYPE (c_quaternion)  S22 
    integer i,j
    integer localmaster
    
     IF(.NOT.C_STABLE_DA) then
     c_bra_v_q%x(0)%i=0 
     RETURN
     endif
    localmaster=c_master

     localmaster=c_master
    call c_ass_quaternion(c_bra_v_q)

 
    call alloc(s22 )
    

      s22=0.0_dp

     do i=1,s1%n
      do j=0,3
        s22%x(j)=s22%x(j) + s1%v(i)*(s2%x(j).d.i)
      enddo
     enddo
 
    c_bra_v_q=s22
    
    call kill(s22  )
    c_master=localmaster

    END FUNCTION c_bra_v_q
    
 FUNCTION c_bra_v_dm( S1, S2 )
    implicit none
    TYPE (c_damap) c_bra_v_dm
    TYPE (c_vector_field), INTENT (IN) :: S1
    TYPE (c_damap), INTENT (IN) ::  S2
    TYPE (c_damap)  S22 

    integer i,j,k
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     c_bra_v_dm%v%i=0
     RETURN
     endif
    localmaster=c_master

     c_bra_v_dm%n=s2%n
    call c_ASSmap(c_bra_v_dm)
     
    s22%n=s2%n
    call alloc(s22)
 

    s22=0

    do i=1,s1%n
     do j=1,s2%n
      s22%v(j)=s22%v(j)+ s1%v(i)*(s2%v(j).d.i)
     enddo

!     do j=1,3
!      do k=1,3
!       s22%s%s(j,k)=s22%s%s(j,k)+ s1%v(i)*(s2%s%s(j,k).d.i)
!      enddo
!     enddo
    enddo
!!  from the spin itself in the vector field replacing the above and adding s1%om
!etienne
    if(use_quaternion) then
     s22%q=s1*s2%q
    endif

    if(complex_extra_order==1.and.special_extra_order_1) c_bra_v_dm=c_bra_v_dm.cut.no
     c_bra_v_dm=s22

    call kill(s22)
    c_master=localmaster

  END FUNCTION c_bra_v_dm



  FUNCTION POWMAP( S1, R2 )
    implicit none
    TYPE (c_damap) POWMAP
    TYPE (c_damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    TYPE (c_damap) S11
    INTEGER I,R22
    integer localmaster
     IF(.NOT.C_STABLE_DA) then
     POWMAP%v%i=0
     RETURN
     endif
    localmaster=c_master

     POWMAP%N=s1%n
    call c_assmap(POWMAP)
    
    s11%n=s1%n
    call alloc(s11)

    s11=1

    R22=IABS(R2)
    DO I=1,R22
       s11=s1*s11
   ENDDO

    IF(R2.LT.0) THEN

 !     CALL c_etinv1(S11%v%i,S11%v%i,s11%n)
        CALL c_etinv(S11,S11)

    ENDIF

  

    powmap=s11
         if(complex_extra_order==1.and.special_extra_order_1) powmap=powmap.cut.no

    call kill(s11)

    c_master=localmaster

  END FUNCTION POWMAP


  FUNCTION pow_tpsaMAP( S1, R2 )
    implicit none
    TYPE (c_damap) pow_tpsaMAP
    TYPE (c_damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    TYPE (c_damap) S11,s0
    INTEGER I,R22
    integer localmaster
    complex(dp) v(lnv)
     v=0
     IF(.NOT.C_STABLE_DA) then
     pow_tpsaMAP%v%i=0
     RETURN
     endif
    localmaster=c_master

     pow_tpsaMAP%N=s1%n
    call c_assmap(pow_tpsaMAP)
    
    s11%n=s1%n
    s0%n=s1%n
    call alloc(s11,s0)

    s11=1

    R22=IABS(R2)
    DO I=1,R22
       s11=s1.o.s11
   ENDDO
    do i=1,s11%n
       v(i)=s11%v(i).sub.'0'
     enddo
    IF(R2.LT.0) THEN

 !     CALL c_etinv1(S11%v%i,S11%v%i,s11%n)
        CALL c_etinv(S11,S11)

    ENDIF

    do i=1,s1%n
       s0%v(i)=(1.0_dp.cmono.i)-v(i)
       !s11%v(i)=s11%v(i)-(s11%v(i).sub.'0')
    enddo
    s11=s11.o.s0
    pow_tpsaMAP=s11
         if(complex_extra_order==1.and.special_extra_order_1) pow_tpsaMAP=pow_tpsaMAP.cut.no

    call kill(s11,s0)

    c_master=localmaster

  END FUNCTION pow_tpsaMAP

  FUNCTION POWMAPs( SS1, R2 )
    implicit none
    TYPE (c_spinmatrix) POWMAPs
    TYPE (c_spinmatrix), INTENT (IN) :: SS1
    INTEGER, INTENT (IN) :: R2
    TYPE (c_damap) s1,S11
    INTEGER I,R22
    integer localmaster
    IF(.NOT.C_STABLE_DA) then
     POWMAPs%s%i=0
     RETURN
     endif
    localmaster=c_master


    call c_ass_spinmatrix(POWMAPs)
    
    s11%n=nv
    call alloc(s11)
    s1%n=nv
    call alloc(s1)

    s11=1
    s1=1
    s1%s=ss1

    R22=IABS(R2)
    DO I=1,R22
       s11=s1*s11
   ENDDO

    IF(R2.LT.0) THEN

 !     CALL c_etinv1(S11%v%i,S11%v%i,s11%n)
        CALL c_etinv(S11,S11)

    ENDIF

  

    powmaps=s11%s
         if(complex_extra_order==1.and.special_extra_order_1) powmaps=powmaps.cut.no

    call kill(s11)
    call kill(s1)

    c_master=localmaster

  END FUNCTION POWMAPs

 SUBROUTINE  c_EQUALMAP(S2,S1)
!*
    implicit none
    type (c_damap),INTENT(inOUT)::S2
    type (c_damap),INTENT(IN)::S1
    integer i
    IF(.NOT.C_STABLE_DA) RETURN
    
    call c_check_snake

 
    do i=1,min(s1%n,s2%n)
       s2%v(i)=s1%v(i)
    enddo
     if(use_quaternion) then
     s2%q=s1%q
     else
     s2%s=s1%s
    endif
     s2%e_ij=s1%e_ij
     s2%x0=s1%x0
     s2%tpsa=s1%tpsa
  END SUBROUTINE c_EQUALMAP

 SUBROUTINE  c_EQUALVEC(S2,S1)
!*
    implicit none
    type (c_vector_field),INTENT(inOUT)::S2
    type (c_vector_field),INTENT(IN)::S1
    integer i
    IF(.NOT.C_STABLE_DA) RETURN
    
    call c_check_snake
 
    do i=1,s1%n 
       s2%v(i)=s1%v(i)
    enddo

 

    do i=0,3
       s2%q%x(i)=s1%q%x(i)
    enddo

     s2%n=s1%n 
     s2%nrmax=s1%nrmax 
     s2%eps=s1%eps

  END SUBROUTINE c_EQUALVEC

SUBROUTINE  c_EQUALcray(S2,S1)
!*
    implicit none
    type (c_ray),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1

     s2%x=0.0_dp
     s2%n=0 
     s2%x0=0.0_dp
     s2%s1=0.0_dp
     s2%s2=0.0_dp
     s2%s3=0.0_dp
     s2%s1(1)=1
     s2%s2(2)=1
     s2%s3(3)=1
     s2%q=0.0_dp
     s2%q%x(s1)=1.0_dp
  END SUBROUTINE c_EQUALcray

  SUBROUTINE  c_IdentityEQUALMAP(S2,S1)
!*
    implicit none
    type (c_damap),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1
    integer i
    real(dp) s1r
    IF(.NOT.C_STABLE_DA) RETURN

    ! if(old) then
        s2%s=s1
        s1r=s1
        s1r=s1
        s2%q=s1r
    IF(S1.EQ.1) then
     do i=1,s2%n
      s2%v(i)=1.0_dp.cmono.i
     enddo
!     s2%s=1
    elseIF(S1.EQ.0)  then
     do i=1,s2%n
      s2%v(i)=(0.0_dp,0.0_dp)
     enddo
!     s2%s=0
   endif

     s2%e_ij=0.0_dp
  END SUBROUTINE c_IdentityEQUALMAP


  SUBROUTINE  c_IdentityEQUALSPIN(S2,S1)
!*
    implicit none
    type (c_spinmatrix),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1
    integer i,j
    IF(.NOT.C_STABLE_DA) RETURN



    IF(S1.EQ.1) then
     do i=1,3
     do j=1,3
      if(i/=j) then 
        s2%s(i,j)=(0.0_dp,0.0_dp)
      else
        s2%s(i,j)=(1.0_dp,0.0_dp)
      endif
     enddo
     enddo
    elseIF(S1.EQ.0) then
     do i=1,3
     do j=1,3
        s2%s(i,j)=(0.0_dp,0.0_dp)
     enddo
     enddo
   endif

  END SUBROUTINE c_IdentityEQUALSPIN

  SUBROUTINE  c_IdentityEQUALSPINOR(S2,S1)
!*
    implicit none
    type (c_spinor),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1
    integer i
    IF(.NOT.C_STABLE_DA) RETURN

     do i=1,3
        s2%v(i)=(0.0_dp,0.0_dp)
     enddo

    IF(S1>0.and.s1<=3) then
      s2%v(s1)=(1.0_dp,0.0_dp)
    endif


  END SUBROUTINE c_IdentityEQUALSPINOR


  SUBROUTINE  c_IdentityEQUALVEC(S2,S1)
!*
    implicit none
    type (c_vector_field),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1
    integer i
    IF(.NOT.C_STABLE_DA) RETURN
     if(s1/=0) then 
       write(6,*) "c_IdentityEQUALVEC"
      stop
     endif
     do i=1,s2%n
        s2%v(i)=(0.0_dp,0.0_dp)
     enddo

 
        s2%q=0.0_dp
     
  END SUBROUTINE c_IdentityEQUALVEC


  SUBROUTINE  c_IdentityEQUALfactored(S2,S1)
!*
    implicit none
    type (c_factored_lie),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1
    integer i
    IF(.NOT.C_STABLE_DA) RETURN
    
    if(s1/=0) then
     s2%n=s1
    else
     s2%n=no
    endif
     do i=1,size(s2%f)
        s2%f(i)=0   !s1
     enddo

  END SUBROUTINE c_IdentityEQUALfactored
  
  SUBROUTINE  matrixMAPr(S2,S1)
!*
    implicit none
    complex(dp),INTENT(inOUT)::S2(:,:)            !(ndim2,ndim2)
    type (c_damap),INTENT(IN)::S1
    integer i,j,JL(lnv)
    IF(.NOT.C_STABLE_DA) RETURN
    call c_check_snake

    do i=1,lnv
       JL(i)=0
    enddo
    ! if(old) then
    do i=1,min(S1%n,size(s2,1))
       do j=1,min(S1%n,size(s2,2))
          JL(j)=1
          call c_dapek(S1%v(i)%i,JL,s2(i,j))
          JL(j)=0
       enddo
    enddo

  END SUBROUTINE matrixMAPr

  SUBROUTINE  r_matrixMAPr(S2,S1)
!*
    implicit none
    real(dp),INTENT(inOUT)::S2(:,:)            !(ndim2,ndim2)
    type (c_damap),INTENT(IN)::S1
    integer i,j,JL(lnv)
    complex(dp) x
 
    IF(.NOT.C_STABLE_DA) RETURN
    call c_check_snake

    do i=1,lnv
       JL(i)=0
    enddo

    ! if(old) then
    do i=1,S1%n
       do j=1,S1%n
          JL(j)=1
          call c_dapek(S1%v(i)%i,JL,x)
          s2(i,j)=x
          JL(j)=0
       enddo
    enddo
 
 
  END SUBROUTINE r_matrixMAPr

  SUBROUTINE  MAPmatrixr(S1,S2)
!*
    implicit none
    complex(dp),INTENT(in)::S2(:,:)       !    (ndim2,ndim2)
    type (c_damap),INTENT(inout)::S1
    integer i,j,JL(lnv)
    IF(.NOT.C_STABLE_DA) RETURN
    do i=1,lnv
       JL(i)=0
    enddo

    do i=1,s1%n
     s1%v(i)=(0.0_dp,0.0_dp)
    enddo

    ! if(old) then
    do i=1,s1%n  !size(s2,1)
       do j=1,s1%n  !,size(s2,2)
          JL(j)=1
          call c_dapok(S1%v(i)%i,JL,s2(i,j))  
          JL(j)=0
       enddo
    enddo

  END SUBROUTINE MAPmatrixr

  SUBROUTINE  r_MAPmatrixr(S1,S2)
!*
    implicit none
    real(dp),INTENT(in)::S2(:,:)       !    (ndim2,ndim2)
    type (c_damap),INTENT(inout)::S1
    integer i,j,JL(lnv)
    complex(dp) x
    IF(.NOT.C_STABLE_DA) RETURN

    do i=1,lnv
       JL(i)=0
    enddo

    do i=1,s1%n
     s1%v(i)=(0.0_dp,0.0_dp)
    enddo

    ! if(old) then
    do i=1,s1%n  !size(s2,1)
       do j=1,s1%n   !,size(s2,2)
          JL(j)=1
          x=s2(i,j)
          call c_dapok(S1%v(i)%i,JL,x)
          JL(j)=0
       enddo
    enddo

  END SUBROUTINE r_MAPmatrixr

!!!!!!!!!!  TPSA LIE PART !!!!!!!


subroutine c_linear_a(xy,a1)
!#internal: normal
!# This routine linearises the linear part of the map ONLY.
!# For a full harmonic system the call c_linear_a(xy,a1) will result in
!# R=a1**(-1)*xy*a1.
!# The map R can be an amplitude dependent rotation, or a rotation followed by a drift
!# in the energy plane, or even a rotation sink if radiation is present. 
!# R can also have rotations for clocks concerning AC modulation. (See Chap.4 of my Springer book)
    implicit none
    integer i,j
    type(c_damap), intent(inout) ::  xy,a1 
    real(dp) reval(ndim2t),imval(ndim2t),vr(ndim2t,ndim2t),vi(ndim2t,ndim2t),vrt(ndim2t,ndim2t),vit(ndim2t,ndim2t)
    real(dp) fm0(ndim2t,ndim2t),x(ndim2t/2),xx(ndim2t/2)
    integer idef(ndim2t/2)
    real(dp), allocatable :: fm(:,:),fmi(:,:),fmii(:,:)
    type(c_damap) s1


    if(.not.c_stable_da) return

 
    idef=0
 
    allocate(fm(xy%n,xy%n),fmi(xy%n,xy%n),fmii(xy%n,xy%n))
 

     fm=xy
     fm0=0

     !!! To understand this, it is better to first understand ndpt=0
     !!! in PTC ndpt=0 includes :    nd2=4 and nd2=6 (cavity on)
     !!! Ndpt =5 is PTC's coasting beam normalisation
     !!! if there are no magnet modulation,then the matrix massaging below will do very little
     !!! see  explanation below at the "else"

     !!! We are diagonalising a Lie map which is just the transposed
     !!! of the matrix in the linear case
     !!! The usual nonlinear map acts on rays and so does its matrix.
     !!! Lie maps, Hamiltonian operators, act on functions.
     !!! In the world of Linear maps, this involves taking the transposed of the matrix
     !!! which represents the linear part of the map
    if(ndpt==0) then
      fm0(1:nd2,1:nd2)=transpose(fm(1:nd2,1:nd2))
    else
      if(c_verbose) then
        write(6,*) " nd2t,nd2,nd2harm "
        write(6,*) nd2t,nd2,nd2harm
      endif
      ! Consider the following example
      ! ndc2t,nd2t,nd2harm,nd2
      !           2           4           6           8
      ! nd2= total size of phase space
      !
      !  ndc2t=2 dimension of coasting phase space (2 or 0), here it is 2
      !
      ! nd2t=4  dimension of pseudo-harmonic  phase space due to orbit (not magnet modulations) 
      ! here it is 4, which is a normal PTC run. 
      ! in accelerator physics, with mid-plane symmetry, it could be 2. But PTC only permits 4 and 6.
      !
      ! nd2harm=6  Total size of pseudo-harmonic  planes including magnet modulation. In the above case,
      ! nd2-nd2harm=2  So there is a coasting plane. 
      !  

      !! The purpose of all the gymnastic below is to exchange the order of the planes
      ! IF AND ONLY IF the longitudinal is coasting
      ! The issue is that  call c_eig6(fm0,reval,imval,vr,vi) should be called only on
      ! matrix which is pseudo-harmonic. 
      fmi=0
      fmii=0
      do i=1,nd2t
        fmi(i,i)=1
        fmii(i,i)=1
      enddo
      fmi(nd2-1,nd2t+1)=1
      fmi(nd2,nd2t+2)=1
      fmii(nd2t+1,nd2-1)=1
      fmii(nd2t+2,nd2)=1

      do i=1,2*rf
        fmi(nd2t+i,nd2t+2+i)=1
        fmii(nd2t+2+i,nd2t+i)=1
      enddo

  !      fmi(1:nd2,1:nd2)=matmul(fmi(1:nd2,1:nd2),fmii(1:nd2,1:nd2))

      fm0(1:nd2,1:nd2)=matmul( fm(1:nd2,1:nd2),fmii(1:nd2,1:nd2))
      fm0(1:nd2,1:nd2)=matmul( (fmi(1:nd2,1:nd2)),fm0(1:nd2,1:nd2))


      fm0(1:nd2harm,1:nd2harm)=transpose(fm0(1:nd2harm,1:nd2harm))
! The diagonalisation is done only on the planes
! 1...nd2harm which are now in the front of the matrix
! Thus the Coasting beam is moved in the last plane  (nd2-1,nd2)
    
    endif
 

    call c_eig6(fm0,reval,imval,vr,vi)

    !! This routine will locate the modulated plane
    !! It assumes that the modulated plane are rotations
    !! It will fail if some orbital planes are exactly rotations, say beta=1.00000000 and no coupling
    call  c_locate_modulated_magnet_planes(fm0,idef,reval)

    call alloc(s1)
 
    if(.not.c_normal_auto) then


      if(c_verbose) then
        write(6,*) " "
        do i=1,nd2harm  !t
          write(6,'(i2,2(1x,G21.14))') i, reval(i),imval(i)
        enddo
        write(6,*) " "
        do i=1,nd2harm  !t
           write(6,*) i
          write(6,'(6(1x,G21.14))') vr(i,1:6)
          write(6,'(6(1x,G21.14))') vi(i,1:6)
        enddo
      endif
      
      call c_locate_planes(vr,vi,idef)

      if(i_piotr(1)==0) then   !.and.nd<=ndharm) then
        write(6,'(a82)') " The order of the planes has been guessed using the algorithm in c_locate_planes "
        write(6,'(a41)') " Hopefully it is correct! Please check! "
        write(6,'(a18,100(i2,1x))') " Order guessed -> ",idef(1:ndharm)

        write(6,*) " Manually identify the location of distinct tunes in the order you like "
        read(5,*) idef(1:ndharm)
      else
        
        do j=1,size(i_piotr)
           idef(j)=i_piotr(j)
        enddo
        
        if(c_verbose) then
          write(6,'(a82)') " The order of the planes has been defined by the user with i_piotr array "
          write(6,'(a20,100(i2,1x))') " Order defuined -> ",idef(1:ndharm)
        endif
        
      endif
    else
      !!  c_locate_planes tries to locate the planes: important if there is little coupling
      !! I suppose that this routine can be refined
      call c_locate_planes(vr,vi,idef)
      if(c_verbose) then
        do i=1,nd2harm  !t
          write(6,'(i2,2(1x,G21.14))') i, reval(i),imval(i)
        enddo
    
        write(6,'(a82)') " The order of the planes has been guessed using the algorithm in c_locate_planes "
        write(6,'(a41)') " Hopefully it is correct! Please check! "
        write(6,'(a18,100(i2,1x))') " Order guessed -> ",idef(1:ndharm)
      endif
    endif 



    if(c_mess_up_vector) then
      vrt=a_mess*vr-b_mess*vi
      vit=a_mess*vi+b_mess*vr
      vr=vrt
      vi=vit
    endif

    x=0.0_dp
    xx=1.0_dp

    

    do i=1,ndharm  !ndt                

      x(i)=0.0_dp
      xx(i)=1.0_dp

      !!! Here the Eigenvectors of the transposed matric (Lie map)
      !!! are normalised so that the Poisson bracket is one.
      !!! If the map is not Hamiltonian (radiation), no harm is done.
      !!! If it is Hamiltonian, that alone will insure that the linear
      !!! canonical transformation is symplectic.

      do j=1,ndharm  !ndt
           x(i)=vr(2*j-1,idef(i))*vi(2*j,idef(i))-vr(2*j,idef(i))*vi(2*j-1,idef(i))+x(i)
      enddo

    enddo

    do i=1,ndharm  !ndt
      if(x(i).lt.0.0_dp) xx(i)=-1.0_dp
      x(i)=SQRT(abs(x(i)))
    enddo


    fm=0
    do i=1,xy%n
      fm(i,i)=1.0_dp
    enddo

    do i=1,nd2harm  !nd2t
      do j=1,ndharm  !ndt
        fm(2*j-1,i)=vr(i,idef(j))*xx(j)/x(j)
        fm(2*j,i)=vi(i,idef(j))/x(j)
      enddo
    enddo
    
   
    a1 = fm

    call alloc(s1)
    !if(courant_snyder_teng_edwards) then
    !!!! Here starts the gymnastics that insures a Courant-Snyder/Teng-Edwards 
    !!! choice A_12=A_34==0
    !
    !    s1=1
    !
    !       do i=1,ndharm  !ndt 
    !          p=ATAN(-fm(2*i-1,2*i)/fm(2*i,2*i))
    !          s1%v(2*i-1) =(COS(p).cmono.(2*i-1))+(sin(p).cmono.(2*i))
    !          s1%v(2*i)   =(COS(p).cmono.(2*i))-(sin(p).cmono.(2*i-1))
    !       enddo
    !       a1=s1*a1
    !       s1=1
    !       fm=a1
    !       ! adjust sa to have sa(1,1)>0 and sa(3,3)>0 rotate by pi if necessary.
    !       do i=1,ndharm  !ndt
    !          p=1.0_dp
    !          if(fm(2*i-1,2*i-1).lt.0.0_dp) p=-1.0_dp
    !           s1%v(2*i-1) = p.cmono.(2*i-1)
    !           s1%v(2*i)   = p.cmono.(2*i)
    !       enddo
    !       a1=s1*a1
    !endif
    !!!!!  In the case of magnet modulations, planes are put back in their places.

    
    if(ndpt/=0) then
      fm=a1
      fm0(1:nd2,1:nd2)=matmul( fm(1:nd2,1:nd2),fmi(1:nd2,1:nd2))
      fm(1:nd2,1:nd2)=matmul( fmii(1:nd2,1:nd2),fm0(1:nd2,1:nd2))
      a1=fm
    endif



    !!! In the case of a symplectic map, without magnet modulation,
    !!! everything is done.
    !!! However magnet modulations do not produce a symplectic matrix
    !!! if there is a longitudinal coasting beam.
    !!! Therefore time must be adjusted to make sure that the map
    !!! is truly diagonalised.
    a1%s=1
    a1%q=1.0_dp
    a1=a1**(-1)
    

    if(rf>0.and.ndpt>0.and.do_linear_ac_longitudinal) then 
    !  if(ndpt>0) then 
      call c_linear_ac_longitudinal(xy,a1,s1) 
    s1%s=1
    s1%q=1.0_dp
      a1=a1*s1
    endif
!    a1%s=1  ! make sure spin is non-zero
    call kill(s1)    
    deallocate(fm,fmi,fmii)   


    return
  end subroutine c_linear_a

subroutine c_linear_a_stoch(xy,a1)
!#internal: normal
!# This routine linearises the linear part of the map ONLY.
!# For a full harmonic system the call c_linear_a(xy,a1) will result in
!# R=a1**(-1)*xy*a1.
!# The map R can be an amplitude dependent rotation, or a rotation followed by a drift
!# in the energy plane, or even a rotation sink if radiation is present. 
!# R can also have rotations for clocks concerning AC modulation. (See Chap.4 of my Springer book)
    implicit none
    integer i,j
    type(c_damap), intent(inout) ::  xy,a1 
    real(dp) reval(ndim2t),imval(ndim2t),vr(ndim2t,ndim2t),vi(ndim2t,ndim2t),vrt(ndim2t,ndim2t),vit(ndim2t,ndim2t)
    real(dp) fm0(ndim2t,ndim2t),x(ndim2t/2),xx(ndim2t/2)
    integer idef(ndim2t/2)
    real(dp), allocatable :: fm(:,:),fmi(:,:),fmii(:,:)
    type(c_damap) s1
    real(dp) :: eps_eigen = 1.e-12_dp,norm

    if(.not.c_stable_da) return

 
    idef=0
 
    allocate(fm(xy%n,xy%n),fmi(xy%n,xy%n),fmii(xy%n,xy%n))
 

     fm=xy
     fm0=0

     !!! To understand this, it is better to first understand ndpt=0
     !!! in PTC ndpt=0 includes :    nd2=4 and nd2=6 (cavity on)
     !!! Ndpt =5 is PTC's coasting beam normalisation
     !!! if there are no magnet modulation,then the matrix massaging below will do very little
     !!! see  explanation below at the "else"

     !!! We are diagonalising a Lie map which is just the transposed
     !!! of the matrix in the linear case
     !!! The usual nonlinear map acts on rays and so does its matrix.
     !!! Lie maps, Hamiltonian operators, act on functions.
     !!! In the world of Linear maps, this involves taking the transposed of the matrix
     !!! which represents the linear part of the map
    if(ndpt==0) then
      fm0(1:nd2,1:nd2)=transpose(fm(1:nd2,1:nd2))
    else
      if(c_verbose) then
        write(6,*) " nd2t,nd2,nd2harm "
        write(6,*) nd2t,nd2,nd2harm
      endif
      ! Consider the following example
      ! ndc2t,nd2t,nd2harm,nd2
      !           2           4           6           8
      ! nd2= total size of phase space
      !
      !  ndc2t=2 dimension of coasting phase space (2 or 0), here it is 2
      !
      ! nd2t=4  dimension of pseudo-harmonic  phase space due to orbit (not magnet modulations) 
      ! here it is 4, which is a normal PTC run. 
      ! in accelerator physics, with mid-plane symmetry, it could be 2. But PTC only permits 4 and 6.
      !
      ! nd2harm=6  Total size of pseudo-harmonic  planes including magnet modulation. In the above case,
      ! nd2-nd2harm=2  So there is a coasting plane. 
      !  

      !! The purpose of all the gymnastic below is to exchange the order of the planes
      ! IF AND ONLY IF the longitudinal is coasting
      ! The issue is that  call c_eig6(fm0,reval,imval,vr,vi) should be called only on
      ! matrix which is pseudo-harmonic. 
      fmi=0
      fmii=0
      do i=1,nd2t
        fmi(i,i)=1
        fmii(i,i)=1
      enddo
      fmi(nd2-1,nd2t+1)=1
      fmi(nd2,nd2t+2)=1
      fmii(nd2t+1,nd2-1)=1
      fmii(nd2t+2,nd2)=1

      do i=1,2*rf
        fmi(nd2t+i,nd2t+2+i)=1
        fmii(nd2t+2+i,nd2t+i)=1
      enddo

  !      fmi(1:nd2,1:nd2)=matmul(fmi(1:nd2,1:nd2),fmii(1:nd2,1:nd2))

      fm0(1:nd2,1:nd2)=matmul( fm(1:nd2,1:nd2),fmii(1:nd2,1:nd2))
      fm0(1:nd2,1:nd2)=matmul( (fmi(1:nd2,1:nd2)),fm0(1:nd2,1:nd2))


      fm0(1:nd2harm,1:nd2harm)=transpose(fm0(1:nd2harm,1:nd2harm))
! The diagonalisation is done only on the planes
! 1...nd2harm which are now in the front of the matrix
! Thus the Coasting beam is moved in the last plane  (nd2-1,nd2)
    
    endif
 

    call c_eig6(fm0,reval,imval,vr,vi)

norm=0.0_dp
do i=1,6
 norm=norm+ abs(reval(i))+abs(imval(i))-1.0_dp
enddo

if(c_verbose) then
do i=1,6
write(6,*) i,reval(i),imval(i)
enddo
write(6,*) " norm ",norm
do i=1,6
write(6,'(i4,6(1x,g13.6))') i,vr(1:6,i)
write(6,'(i4,6(1x,g13.6))') i,vi(1:6,i)
enddo
endif


if(norm<eps_eigen) then
 a1=1
 deallocate(fm,fmi,fmii)   
return
endif
    !! This routine will locate the modulated plane
    !! It assumes that the modulated plane are rotations
    !! It will fail if some orbital planes are exactly rotations, say beta=1.00000000 and no coupling
  
    x=0.0_dp
    xx=1.0_dp

    idef(1)=1
    idef(2)=3
    idef(3)=5

    do i=1,ndharm  !ndt                

      x(i)=0.0_dp
      xx(i)=1.0_dp

      !!! Here the Eigenvectors of the transposed matric (Lie map)
      !!! are normalised so that the Poisson bracket is one.
      !!! If the map is not Hamiltonian (radiation), no harm is done.
      !!! If it is Hamiltonian, that alone will insure that the linear
      !!! canonical transformation is symplectic.

      do j=1,ndharm  !ndt
           x(i)=vr(2*j-1,idef(i))*vi(2*j,idef(i))-vr(2*j,idef(i))*vi(2*j-1,idef(i))+x(i)
      enddo

    enddo

    do i=1,ndharm  !ndt
      if(x(i).lt.0.0_dp) xx(i)=-1.0_dp
      x(i)=SQRT(abs(x(i)))
    enddo


    fm=0
    do i=1,xy%n
      fm(i,i)=1.0_dp
    enddo

    do i=1,nd2harm  !nd2t
      do j=1,ndharm  !ndt
        fm(2*j-1,i)=vr(i,idef(j))*xx(j)/x(j)
        fm(2*j,i)=vi(i,idef(j))/x(j)
      enddo
    enddo
    
   
    a1 = fm

    call alloc(s1)
 
    
 



    !!! In the case of a symplectic map, without magnet modulation,
    !!! everything is done.
    !!! However magnet modulations do not produce a symplectic matrix
    !!! if there is a longitudinal coasting beam.
    !!! Therefore time must be adjusted to make sure that the map
    !!! is truly diagonalised.
    a1%s=1
    a1%q=1.0_dp



 

    a1=a1**(-1)
    

 
!    a1%s=1  ! make sure spin is non-zero
    call kill(s1)    
    deallocate(fm,fmi,fmii)   


    return
  end subroutine c_linear_a_stoch

subroutine c_locate_planes(vr,vi,idef)
!#restricted: normal
!# Tries to locate the planes in c_linear_a, 
!# so that if decoupled Qx correspounds to Q1 it is properly identified.
!# If c_normal_auto.false. it is not called and the user is 
!# asked to choose planes manually. 
    implicit none
    real(dp), intent(in) ::  vr(ndim2t,ndim2t),vi(ndim2t,ndim2t)
    real(dp) t(ndim2t,ndim2t)
    integer idef(ndim2t/2)
    real(dp) r,rmax
    logical doit
    integer j,k,jmax,i,irf
!    idef=0
! write(6,*) " nd2,ndc2t ",nd2,ndc2t 
!pause 
irf=0
!      do j=2,nd2-ndc2t,2
     do j=1,nd-ndct -irf*rf
 !      rmax=0
 !      kmax=0
       do k=1,(nd2-ndc2t)/2
         r=abs(vr(2*k-1,2*j))+abs(vr(2*k,2*j))+abs(vi(2*k-1,2*j))+abs(vi(2*k,2*j))
         t(j,k)=r


       enddo
  !    idef(j/2)=2*kmax-1
    enddo
   

do k=1,(nd2-ndc2t)/2
       rmax=0
       jmax=0
     do j=1,nd-ndct-irf*rf


        r=t(j,k) 

          if(r>rmax) then
            rmax=r
            jmax=j
          endif

       enddo

       idef(k)=2*jmax-1
     enddo


 !!!! checking   maybe equal tunes....
      doit=.false.
      do j=1,(nd2-ndc2t-irf*2*rf )/2
      do k=1,(nd2-ndc2t-irf*2*rf )/2
 
         if(j==k) cycle
        if(idef(j)==idef(k)) doit=.true.
     
        if(doit) exit
      enddo
      enddo

     if(doit) then
       if(c_verbose) write(6,*) "warning : trouble locating planes, so chosen arbitrarily "
       do k=1,(nd2-ndc2t-irf*2*rf )/2
        idef(k)=2*k-1
        enddo
     endif
!       if(abs(reval(j)-r)<epsflo) then
!        idef(i/2)=j


end subroutine c_locate_planes

subroutine c_locate_modulated_magnet_planes(fm0,idef,reval)
!#restricted: normal
! It locates modulated magnet planes in c_linear_a.
    implicit none
    real(dp), intent(in) ::  fm0(ndim2t,ndim2t)
    integer idef(ndim2t/2)
    real(dp) reval(ndim2t),r 

    integer i,j
    idef=0


    do i=nd2-ndc2t,nd2-ndc2t-2*rf+1,-2

     r=fm0(i,i)
     do j=1,nd2
       if(abs(reval(j)-r)<epsflo) then
        idef(i/2)=j
        exit
       endif
     enddo
    enddo

end subroutine c_locate_modulated_magnet_planes

subroutine c_linear_ac_longitudinal(xy,a1,ac)
!#internal:  normal
!# This routine fixes the linear normal form when a Jordan normal form is mixed with 
!# an ac-modulated clock.
!# This routines makes the adjustements described in Sec.4.7.4 of my Springer book.
    implicit none
    integer i
    type(c_damap), intent(inout) ::  xy,a1,ac 
    complex(dp) m(ndim2t,ndim2t) , cc

    type(c_damap) s1
    integer jc(ndim2t)
    complex  r(ndim2t)

    if(.not.c_stable_da) return



    call alloc(s1)
     
     !!! s1 is diagonal and a drift in the longitudinal plane
     !!! if no AC present
     s1=from_phasor(-1)*a1**(-1)*(xy.sub.1)*a1*from_phasor()
      

     ac=1
     m=s1
     do i=1,nd2
      jc=0;jc(i)=1;
      r(i)=s1%v(i).sub.jc  ! storing the complex tunes
     enddo

     do i=nd2,nd2-2*rf+1,-1 ! Do looping over the AC frequencies
      jc=0; jc(i)=1;
      cc=s1%v(ndptb).sub.jc  ! these are non zero if AC present
      cc=cc/(r(i)-1.0_dp)
      ac%v(ndptb)=ac%v(ndptb)+(cc.cmono.jc) ! AC map will remove these terms
     enddo

 !    do i=1,nd2! Do looping over all frequencies
 !     if(i/=ndptb.and.i/=ndpt) then
 !      jc=0; jc(i)=1;
 !      cc=s1%v(ndptb).sub.jc  ! these are non zero if AC present or nonsymplectic coordinates used
 !      cc=cc/(r(i)-1.0_dp)
 !      ac%v(ndptb)=ac%v(ndptb)+(cc.cmono.jc) ! AC map will remove these terms
 !     endif
 !    enddo
 
!!!!!!!!!     s1=ac**(-1)*s1*ac  ! this should be fully normalised

     ac=from_phasor()*ac*from_phasor(-1)  ! puting AC in cartesian
 
    call kill(s1)

  end subroutine c_linear_ac_longitudinal


 subroutine c_gofix(xy,a1)
!#restricted:  normal 
!# The routine c_gofix computes a1 such that 
!# "call c_gofix(m,a1)" produces
!# m0=a1**(-1)*xy*a1 to order "order_gofix" around the parameter dependent fixed point.
!# order_gofix is defaulted to 1.
!# dosymp is defaulted to false. 
!# dosymp should be true if order_gofix > 1 (if c_gofix is used on its own).
!# c_canonise will take care of a1 if needed to higher order.

    implicit none
    integer i,j,ndloc
    type(c_damap), intent(inout) ::  xy,a1 
    type(c_damap) w,v,rel,x
    type(c_taylor) t1
    real(dp), allocatable :: mt(:,:),mv(:)
    integer, allocatable :: je(:)

    
  !  ndct=iabs(ndpt-ndptb)  ! 1 if coasting, otherwise 0
  !  ndc2t=2*ndct  ! 2 if coasting, otherwise 0
  !  nd2t=nd2-2*rf-ndc2t   !  size of harmonic oscillators minus modulated clocks
  !  ndt=nd2t/2        ! ndt number of harmonic oscillators minus modulated clocks
  !  nd2harm=nd2t+2*rf  !!!!  total dimension of harmonic phase space
  !  ndharm=ndt+rf  !!!! total number of harmonic planes
!   write(6,*) " nd2,ndct,ndc2t,nd2t,ndt,nd2harm,rf,ndpt "
!   write(6,'(8(1x,i4))') nd2,ndct,ndc2t,nd2t,ndt,nd2harm,rf,ndpt

  if(.not.c_stable_da) return

     ndloc=0
                 

    ! COMPUTATION OF A1 AND A1I USING DAINV
    rel%n=nv; w%n=nv; v%n=nv;x%n=nv
    call alloc(rel);call alloc(w);call alloc(v);call alloc(x);
    call alloc(t1);

    rel=1
    v=1
 

    x=xy

    do i=1,nd2  
     if(i/=ndpt.and.i/=ndptb) then
        v%v(i)=x%v(i)-rel%v(i)    !   V= X-1  where X is the map
       else
        if(mod(i,2)==0) ndloc=i/2  ! degree of freedom of coasting plane
      endif
    enddo

    v=v.cut.(order_gofix+1)
 
    w=v**(-1)    !  W= (Map-1)**-1   

    x=0
    x%s=1    ! spin part is identity
    x%q=1.0_dp 
    do i=nd2+1,nv
     x%v(i)=1.d0.cmono.i  !  Identity in all the parameter planes
    enddo
    if(ndpt/=0) x%v(ndpt)=1.d0.cmono.ndpt !  If coasting, then energy is a parameter

    v=w*x  ! v contains the fixed point, for example v(1)= eta_x * x_pt + ...
         

    a1=v

  ! however a1 is not a  transformation, we must add the identity (done at the end) 
  ! also we must add some stuff to time to make it symplectic if coasting
  ! because the Lie operator which produces v(1)= eta_x * x_pt + ...
  ! is, within a sign, eta_x * x_pt * p_x-eta_x_prime * x_pt * x-...  which affects time
if(ndpt/=0) then
if(dosymp) then
      ! if(ndpt/=0) then 
         t1=0
         do i=1,nd
           if(i/=ndloc) then
            x%v(2*i)  =(-1)**(2*i-1)*(a1%v(2*i-1))
            x%v(2*i-1)=(-1)**(2*i  )*(a1%v(2*i))
            v%v(2*i)=   x%v(2*i).d.ndpt
            v%v(2*i-1)= x%v(2*i-1).d.ndpt

          endif
         enddo
         do i=1,nd
           if(i/=ndloc) then
             t1=-(1.0_dp.cmono.(2*i-1))*v%v(2*i-1)+t1  ! first order
             t1=-(1.0_dp.cmono.(2*i))*v%v(2*i)+t1      ! first order
             t1=-0.5_dp*(x%v(2*i-1)*v%v(2*i)-x%v(2*i)*v%v(2*i-1))+t1  ! second order
! because  eta_x * x_pt * p_x-eta_x_prime * x_pt * x-...  is linear in the transverse
!          there are NO terms higher than second order
           endif
         enddo
         t1=(-1)**ndpt*t1
         a1%v(ndptb)=(1.0_dp.cmono.ndptb)+t1 !!! effect on  time added to identity map in the time-energy plane
         a1%v(ndpt)=1.0_dp.cmono.ndpt
      !  endif
!!!! end of the coasting beam gymnastic !!!!!!!!!!!!!!

!!! identity is added to a1 except coasting plane (already there) !!! 
    do i=1,nd2 
     if(i/=ndptb.and.i/=ndpt) a1%v(i)=(1.d0.cmono.i)+a1%v(i)
    enddo

!call print(a1,6)
!pause 333

else

!  v=1
  v=xy

 allocate(mt(nd2t,nd2t),mv(nd2t)); allocate(je(nv));

  mt=0.0_dp
 do i=1,nd2t
  do j=1,nd2t
   je=0
   je(j)=1
   mt(i,j)=v%v(i).sub.je
  enddo
 enddo

  mt=-transpose(mt)

  do i=1,nd2t
     mt(i,i)=1.0_dp+mt(i,i)
  enddo
  call matinv(mt,mt,nd2t,nd2t,i)



 mv=0.0_dp
  do i=1,nd2t
   je=0
    je(i)=1
    mv(i)=-xy%v(ndptb).sub.je
  enddo

 mv=matmul(mt,mv)

    do i=1,nd2 
     if(i/=ndpt) then 
       a1%v(i)=(1.d0.cmono.i)+a1%v(i)  ! ndpt is already identity
     endif
    enddo 

    do i=1,nd2t 
     a1%v(ndptb)=a1%v(ndptb)+mv(i)*(1.d0.cmono.i)
    enddo
 deallocate(mt,mv);deallocate(je);
! endif ! npdt/=0
 endif       

else ! npdt=0

    do i=1,nd2 
       a1%v(i)=(1.d0.cmono.i)+a1%v(i)  ! ndpt is already identity
    enddo 

endif

 
 
    call kill(t1);
    call kill(v);call kill(w);call kill(rel);call kill(x);
    return
  end subroutine c_gofix

 
subroutine c_factor_map(m,l,f,dir)  
!#general: manipulation
!#  This routine factors a map m as:
!#  m= l exp(f.grad) if dir=1
!#  m= exp(f.grad) l if dir =-1
!#  l is assumed to be a linear composition operator https://en.wikipedia.org/wiki/Composition_operator .
implicit none
    type(c_damap) , intent(inout) :: m,l
    type(c_vector_field), intent(inout) :: f
    integer dir
    type(c_damap) t
    real(dp) epso
 
    call alloc(t)

    t=m

 
     l=t.sub.(-1)
 

    if(dir==1) then  ! Dragt-Finn direction 
     t=t*l**(-1)
    else
     t=l**(-1)*t
    endif
     !epso=-no
     f=c_logf_spin(t)  !,epso=epso)
     !call c_flofacg(t,f,epso)
     
    call kill(t)
    
end subroutine c_factor_map
!!  no spin here
subroutine c_canonise(at,a_cs,a0,a1,a2,phase,irot) 
!#general: manipulation
!# This routine is of great pedagogical importance.
!# It is restricted to the orbital motion,
!# c_full_canonise also includes spin.
!# It factors a canonical transformation as
!# at=a_cs o rotation(phase) as explaned in Chap.7 of my Springer book.
!# a_cs = a_0 o a_1 o a_2
!# The flag courant_snyder_teng_edwards, defaulted to true, controls the linear part of a_cs.
    implicit none                             
    type(c_damap) , intent(inout) :: at,a_cs 
    type(c_damap) , optional, intent(inout) :: a2,a1,a0
    type(c_taylor) ,optional, intent(inout) :: phase(:)
     integer ,optional, intent(in) :: irot
    call c_full_canonise(at,a_cs,a0=a0,a1=a1,a2=a2,phase=phase,irot=irot)

end subroutine c_canonise


subroutine c_full_factor_map(U,Q,U_0,U_1,U_2) 
!#general: manipulation
!# U = Q o U_0 o U_1 o U_2 
     implicit none
    type(c_damap) , intent(inout) :: U,Q,U_0,U_1,U_2
    type(c_damap) a
    call alloc(a)
    qphase=.false.
    call c_full_canonise(U,a,Q,U_0,U_1,U_2,irot=0) 
    qphase=.true.
    call kill(a)
end subroutine c_full_factor_map   

subroutine c_full_canonise(at,a_cs,as,a0,a1,a2,rotation,phase,nu_spin,irot) 
!#general: manipulation
!# This routine is of great pedagogical importance.
!# It factors a canonical transformation as
!# at=a_cs o rotation(phase,nu_spin) as explained in Chap.7 of my Springer book.
!# a_cs = a_s o a_0 o a_1 o a_2
    implicit none
    type(c_damap) , intent(inout) :: at,a_cs 
    type(c_damap) , optional, intent(inout) :: as,a2,a1,a0,rotation
    type(c_taylor) ,optional, intent(inout) :: phase(:),nu_spin
    integer ,optional, intent(in) :: irot
    type(c_damap) ar,att,phi,a0t,a1t,a2t,ast
    type(c_taylor) pha,tune_spin
    integer i,kspin,ir
    real(dp) norm
    ir=1
    if(present(irot)) ir=irot
    call alloc(ar)
    call alloc(att)
    call alloc(a0t)
    call alloc(a1t)
    call alloc(a2t)    
    call alloc(ast)
    call alloc(phi)
    call alloc(pha,tune_spin)   
    
  !  at= (a,s) =  (a,I) o  (I,s)
  !  call c_full_norm_spin(at%s,kspin,norm)  
    call c_full_norm_spin_map(at,kspin,norm)
 
 ! storing the spin because the code is careless (sets it to zero)   
      if(kspin==-1) then
         att=at
         att%s=1
         att%q=1.0_dp
         ast=1
         ast%s=at%s
         ast%q=at%q 
      else
       att=at
       ast=1
      endif
 
!  at= (a,s) =  (att,I) o  (I,ats)

      att%s=0
      att%q=0.0_dp

      ar=1
 
    call extract_a0(att,a0t)
 
    phi=1

!call print(phi%v(6),6)
!pause 1
    if(ir/=1) then
    call extract_a1(att,a1t)
   else
    call extract_a1(att,a1t,phi)
   endif
!    if(present(phase))     ar=ar*phi
!call print(phi%v(6),6)
!pause 2
    if(ir/=1) then
        call extract_a2(att)
     else
    call extract_a2(att,ar)
      phi=ar*phi
    endif
 
    a2t=att
    if(kspin==-1) then 
    if(use_quaternion)   THEN
      ast%q=ast%q*phi**(-1)
     else
      ast%s=ast%s*phi**(-1)
     endif

    endif
    if(present(phase).or.present(rotation))     then 
         if(present(rotation)) then
           rotation=phi
           rotation%s=1
           rotation%q=1.0_dp
         endif
         phi=c_simil(from_phasor(-1),phi,1)
    endif

    if(use_quaternion)   THEN
      a0t%q=1.0_dp
      a1t%q=1.0_dp
      a2t%q=1.0_dp
     else
      a0t%s=1
      a1t%s=1
      a2t%s=1
    endif
    if(kspin==-1.and.ir==1) then

     call c_remove_y_rot(ast,ar,tune_spin)
     if(present(nu_spin) ) nu_spin=-tune_spin/twopi+nu_spin      !  changed 2018.11.01
         if(present(rotation)) then
           rotation=rotation*ar
         endif     
    endif
    a_cs=a0t*a1t*a2t

      if(kspin==-1) then
        if(use_quaternion) then
         ast%q=ast%q*a_cs**(-1)
        else
         ast%s=ast%s*a_cs**(-1)    !! at= ast*a_cs * rotation 
        endif
      endif

    a_cs=ast*a_cs

    if(present(phase).and.ir==1) then
     do i=1,nd2t/2
       pha=(phi%v(2*i).k.(2*i))
       pha=(-i_*log(pha)/twopi).cut.no
       phase(i)=pha+phase(i)
     enddo
      if(ndpt/=0) then
        if(mod(ndpt,2)==0) then
         i=ndpt/2
        else
         i=ndptb/2
        endif
       phase(i)=phase(i)+phi%v(ndptb)-(1.0_dp.cmono.ndptb)
      endif
    endif
 

    if(present(a0)) a0=a0t
 
    if(present(a1)) a1=a1t

    if(present(a2)) a2=a2t

    if(present(as)) as=ast 

    call kill(ar)
    call kill(att)
    call kill(a0t)
    call kill(a1t)
    call kill(a2t)
    call kill(ast)
    call kill(phi)
    call kill(pha,tune_spin)
end subroutine c_full_canonise

subroutine c_identify_resonance(j,n,c) 
    implicit none
    integer, intent(inout) :: J(:)
    complex(dp), intent(out) :: c
    integer, intent(out) :: n
    integer i

    do i=1,ndt*2
       if(j(i)/=0) exit
    enddo
        n=i
        if(mod(n,2)==0) then
          n=n-1
         else
          n=n+1
        endif
        c=1.0_dp/(j(i)*n_cai)
        j(i)=j(i)-1
        

end subroutine c_identify_resonance

subroutine c_full_factorise(at,as,a0,a1,a2,dir) 
!#general: manipulation
!# a_t = a_0 o a_1 o a_2 o a_s  for dir=1
    implicit none
    type(c_damap) , intent(inout) :: at 
    type(c_damap) , optional, intent(inout) :: as,a2,a1,a0
    integer,optional :: dir

    type(c_damap) att,a0t,a1t,a2t,ast
    type(c_taylor) p
    integer i,kspin,ii
    real(dp) norm


    call alloc(att)
    call alloc(a0t)
    call alloc(a1t)
    call alloc(a2t)    
    call alloc(ast)
 
 
    ii=1
    if(present(dir)) ii=dir
  !  at= (a,s) =  (a,I) o  (I,s)
    call c_full_norm_spin(at%s,kspin,norm)  
 ! storing the spin because the code is careless (sets it to zero)   
      if(kspin==-1) then
         att=at
         att%s=1
         att%q=1.0_dp
         ast=1
         ast%s=at%s 
      else
       att=at
       ast=1
      endif
 
!  at= (a,s) =  (att,I) o  (I,ats)

      att%s=0


 
    call extract_only_a0(att,a0t)

!call print(phi%v(6),6)
!pause 1
    
    call extract_only_a1(att,a1t)
!    if(present(phase))     ar=ar*phi
!call print(phi%v(6),6)
!pause 2
!    if(no>1) 
 !   call extract_a2(att,ar)

!call print(phi%v(6),6)
!pause 3

    a2t=att



    a0t%s=1
    a1t%s=1
    a2t%s=1
    a0t%q=1.0_dp
    a1t%q=1.0_dp
    a2t%q=1.0_dp


    if(ii==-1) then
     att=a0t*a1t*a2t
     ast=att*ast*att**(-1)
     a1t=a0t*a1t*a0t**(-1)
     a2t=a0t*a2t*a0t**(-1)
     a2t=a1t*a2t*a1t**(-1)
    endif


    if(present(a0)) a0=a0t
 
    if(present(a1)) a1=a1t

    if(present(a2)) a2=a2t

    if(present(as)) as=ast 


    call kill(att)
    call kill(a0t)
    call kill(a1t)
    call kill(a2t)
    call kill(ast)
 
 
end subroutine c_full_factorise

 
 subroutine c_normal(xy,n,dospin,no_used,rot,phase,nu_spin)
!#general:  normal
!# This routine normalises the map xy
!# xy = n%a_t**(-1)*r*n%a_t 
!# The linear part of r is described in Chap.4 for the orbital part
!# and in Chap.6 for the spin. The nonlinear parts are in Chap.5 and 6.
!# Dospin must be set to .true. if spin is to be normalised.
!# Resonances can be left in the map. Thir number is in n%nres.
!# They are nres rosnances The kth resonance is n%m(i,k).Q_i+n%ms(k)=integer

    implicit none
    type(c_damap) , intent(inout) :: xy
    type(c_damap) m1,ri,nonl,a1,a2,mt,AS 
    type(c_normal_form), intent(inout) ::  n
    type(c_damap), optional :: rot
    type(c_taylor), optional :: phase(:),nu_spin
    type(taylor) c1,s1
    integer,optional :: no_used
    integer i,j,k,l,kr,not
    integer, allocatable :: je(:)
    logical(lp) removeit,rad_in
    complex(dp) v,lam,egspin(3)
    complex(dp), allocatable :: eg(:)
    real(dp) norm,alpha,prec !,cx,sx
    logical(lp), optional :: dospin
    logical dospinr
    type(c_spinor) n0,nr
    type(c_quaternion) qn0,qnr
    integer mker, mkers,mdiss,mdis
    if(lielib_print(13)/=0) then
     call kanalnummer(mker,"kernel.txt")
     call kanalnummer(mdis,"distortion.txt")
     call kanalnummer(mkers,"kernel_spin.txt")
     call kanalnummer(mdiss,"distortion_spin.txt")
    endif

    not=no
    if(present(no_used)) then
      not=no_used  ! sometimes only linear stuff is needed
    else
       if(complex_extra_order==1.and.special_extra_order_1) not=not-1
    endif

    dospinr=.false.
    if(present(dospin)) dospinr=dospin
    call alloc(m1);call alloc(nonl);call alloc(a1);call alloc(a2);call alloc(ri);
 
    allocate(je(nv))    
    allocate(eg(xy%n))

    
    m1=xy


    ! Brings the map to the parameter dependent fixed point
    ! including the coasting beam gymnastic: time-energy is canonical
    ! but energy is constant. (Momentum compaction, phase slip etc.. falls from there)
 
    call  c_gofix(m1,a1) 

     m1=c_simil(a1,m1,-1)
    ! Does the the diagonalisation into a rotation
    call c_linear_a(m1,a2)  

    !!! Now the linear map is normalised
    m1=c_simil(a2,m1,-1)
    !!! We go into the phasors' basis
    m1=c_simil(from_phasor(-1),m1,1)

 
    ri=(m1.sub.-1)**(-1) 
    ri%s=1  ! make spin identity
    ri%q=1.0_dp  ! make spin identity


    !!! The tunes are stored for the nonlinear normal form recursive algorithm
    do k=1,xy%n
      if(coast(k)) then
        eg(k)=1
       else
        je=0
        je(k)=1
        eg(k)=ri%v(k).sub.je
       endif  
    enddo
 
    n%ker=0  ! In case reusing normal form

    do i=2,not
      if(lielib_print(13)/=0) then
        write(mdis,*) " **************************************** " 
        write(mdis,*) "Order ",i
        write(mker,*) " **************************************** " 
        write(mker,*) "Order ",i
      endif

      nonl=(m1*ri)
      nonl= exp_inv(n%ker,nonl)
      nonl=nonl.sub.i


      do k=1,xy%n
        if(lielib_print(13)/=0) then
          write(mdis,*) " **************************************** " 
          write(mdis,*) "field component ",k
          write(mker,*) " **************************************** " 
          write(mker,*) "field component ",k
        endif

        n%g%f(i)%v(k)=0.0_dp
        n%ker%f(i)%v(k)=0.0_dp


        j=1

        do while(.true.) 

           call  c_cycle(nonl%v(k),j,v ,je); if(j==0) exit;
           call check_kernel(k,xy%n,je,removeit)

           if(n%nres>0.and.removeit) then 
             do kr=1,n%nres
               if(n%ms(kr)/=0) cycle  ! a spin resonance
               call check_resonance(k,xy%n,je,kr,n%m,removeit)
               if(.not.removeit) then
                 exit
               endif
             enddo
           endif

          if(removeit) then

            lam=1.0_dp
            je(k)=je(k)-1
            do l=1,xy%n 
              if(coast(l)) cycle 
              lam=lam*eg(l)**je(l)
            enddo

            if(lielib_print(13)/=0) then
                 write(mdis,*) k
                 write(mdis,'(6(1x,i4))') je(1:nd2)
                 write(mdis,*) v
                 write(mdis,*) abs(v/(1-lam))
            endif

            je(k)=je(k)+1

            n%g%f(i)%v(k)=n%g%f(i)%v(k)+(v.cmono.je)/(1.0_dp-lam)

          else ! Put in the kernel

            if(lielib_print(13)/=0) then
               je(k)=je(k)-1
               write(mker,*) k
               write(mker,'(6(1x,i4))') je(1:nd2)
               write(mker,*) v
               write(mker,*) abs(v/(1-lam))
               je(k)=je(k)+1
            endif
               n%ker%f(i)%v(k)=n%ker%f(i)%v(k)+(v.cmono.je)
            endif

        enddo  ! over monomial
      enddo  ! over vector index

      m1=c_simil(n%g%f(i),m1,-1)

    enddo

    n%a_t=a1*a2*from_phasor()*texp(n%g)*from_phasor(-1)

    n%a1=a1
    n%a2=a2

!!!!!   here we put the normalised linear part into the factored vector field
!!!!!   not necessary but useful
    do k=1,xy%n
       if(.not.coast(k)) then    
        je=0
        je(k)=1     
        n%ker%f(1)%v(k)=n%ker%f(1)%v(k)-(log(eg(k)).cmono.je)

        if(mod(k,2)==1) then
            n%tune((k+1)/2)=aimag(log(eg(k)))/twopi
            n%damping((k+1)/2)=real(log(eg(k)))
            if(n%tune((k+1)/2)<0.and.n%positive) n%tune((k+1)/2)=n%tune((k+1)/2)+1.0_dp
        endif
       endif 
      enddo

        if(nd2t==6) then
           if(n%tune(3)>0.5d0) n%tune(3)=n%tune(3)-1.0_dp
        endif 

       if(ndpt/=0) then
        je=0
        je(ndpt)=1
        lam=(ri%v(ndptb).sub.je) 
        n%ker%f(1)%v(ndptb)=n%ker%f(1)%v(ndptb)-(lam.cmono.je)
            if(mod(ndpt,2)==0) then
              n%tune(ndpt/2)=-lam
            else
              n%tune(ndptb/2)=-lam
            endif
       endif


    if(dospinr) then

if(use_quaternion)then
      call c_full_norm_quaternion(m1%q,k,norm) 
else
      call c_full_norm_spin(m1%s,k,norm)   
endif
      if(k>=0) then
        dospinr=.false.
         if(use_quaternion)  then
           write(6,*) " no quaternion spin in map: dospin command ignored "
         else
            write(6,*) " no spin matrix in map: dospin command ignored "
        endif
     endif
    endif


    if(dospinr) then
      call alloc(n0) 
      call alloc(nr)
      call alloc(mt) 
      call alloc(AS) 
      call alloc(qnr)
      n%AS=1
 

if(use_quaternion)then
      call c_normal_spin_linear_quaternion(m1,m1,n%AS,alpha)

      ri=1 ; ri%q=m1%q.sub.0 ; ! exp(theta_0 L_y)   (2)
!      sx=sqrt(ri%q%x(1)**2+ri%q%x(2)**2+ri%q%x(3)**2)
!      cx=ri%q%x(0)
!write(6,*) alpha
!      alpha=-(-2*atan2(sx,cx))
!write(6,*) alpha
!pause 723
      egspin(3)=cos(alpha)-i_*sin(alpha)
      egspin(2)=1.0_dp
      egspin(1)=cos(alpha)+i_*sin(alpha) 
else
       call c_normal_spin_linear(m1,m1,n%AS,n0)  ! (1)
       ri=1 ; ri%s=m1%s.sub.0 ; ! exp(theta_0 L_y)   (2)
      egspin(3)=ri%s%s(1,1)-i_*ri%s%s(1,3)
      egspin(2)=1.0_dp
      egspin(1)=ri%s%s(1,1)+i_*ri%s%s(1,3)
endif
 
 



 
      if(lielib_print(13)/=0) then
        write(mdiss,*) " eg(1:4),spin_def_tune" ,spin_def_tune
        write(mdiss,*)eg(1)
        write(mdiss,*)eg(2)
        write(mdiss,*)eg(3)
        write(mdiss,*)eg(4)
        write(mdiss,*) " egspin(1:3)" 
        write(mdiss,*)egspin(1)
        write(mdiss,*)egspin(2)
        write(mdiss,*)egspin(3)
      endif
      if(lielib_print(13)/=0) then
        write(mkers,*) " eg(1:4),spin_def_tune" ,spin_def_tune
        write(mkers,*)eg(1)
        write(mkers,*)eg(2)
        write(mkers,*)eg(3)
        write(mkers,*)eg(4)
        write(mkers,*) " egspin(1:3)" 
        write(mkers,*)egspin(1)
        write(mkers,*)egspin(2)
        write(mkers,*)egspin(3)
      endif
      
      !!! tune is taken from egspin(1) or egspin(3)   spin_def_tune= +/- 1
       n%spin_tune=aimag(log(egspin(2-spin_def_tune))/twopi)  
 
      ! because  exp(a L_y) x = x- a z + O(a**2)
       ri=ri**(-1) ! exp(-alpha_0 L_y)   (3)


if(use_quaternion)then
       nonl=m1.sub.1 ; nonl%q=1.0_dp ;nonl=nonl**(-1)  ! R_0^-1      (4)  
else
     nonl=m1.sub.1 ; nonl%s=1 ;nonl=nonl**(-1)  ! R_0^-1      (4)  
endif
!       nonl=m1.sub.1 ; nonl%s=1 ;nonl=nonl**(-1)  ! R_0^-1      (4)          
        

       do i=1,no    !+2
          if(lielib_print(13)/=0) then
            write(mdiss,*) " **************************************** " 
            write(mdiss,*) "Order ",i
            write(mkers,*) " **************************************** " 
            write(mkers,*) "Order ",i
          endif
  
        
          mt=m1*ri !  S*exp(-theta_0 L_y)    (5)
 
 
if(use_quaternion)then
       n0=mt%q
else
      call c_find_om_da(mt%s,n0)   ! (4)  
endif

           call c_n0_to_nr(n0,n0)   ! n0 = > eigen-operator of spin   (7)
          n0=n0*nonl               !  no * R^-1      (8)

          nr=0
          
          do k=1,3
            if(lielib_print(13)/=0) then 
              write(mdiss,*) " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ " 
              write(mdiss,*) "Spin component ",k
              write(mdiss,*) " "
              write(mkers,*) " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ " 
              write(mkers,*) "Spin component ",k
              write(mkers,*) " "
            endif
            
            j=1
            do while(.true.) 
              call  c_cycle(n0%v(k),j,v ,je); if(j==0) exit;
              call check_kernel_spin(k,xy%n,je,removeit)
              if(n%nres>0.and.removeit) then 
                do kr=1,n%nres
                  call check_resonance_spin(k,xy%n,je,kr,n%ms,n%m,removeit)
                  if(.not.removeit) then
                    exit
                  endif
                enddo
              endif
                
              if(removeit) then
                lam=egspin(k) 
                do l=1,xy%n 
                  if(coast(l)) cycle 
                  lam=lam*eg(l)**je(l)
                enddo
               
                if(lielib_print(13)/=0) then 
                  !do kr=1,nd2
	!je(kr)=-(-1)**kr*je(kr)
                  !enddo
                  write(mdiss,'(6(1x,i4))') je(1:nd2)
                  write(mdiss,*)lam
                  write(mdiss,*) v
                  write(mdiss,*) abs(v/(1-lam))
                  !do kr=1,nd2
	!je(kr)=-(-1)**kr*je(kr)
                  !enddo
                endif
              nr%v(k)=nr%v(k) +(v.cmono.je)/(1.0_dp-lam)   ! (9)
            else
              if(lielib_print(13)/=0) then 
                do kr=1,nd2
                  je(kr)=-(-1)**kr*je(kr)
                enddo      
                write(mkers,'(6(1x,i4))') je(1:nd2)
                write(mkers,*) v
                do kr=1,nd2
                  je(kr)=-(-1)**kr*je(kr)
                enddo
              endif
            endif
          enddo ! cycle
        enddo ! k
        
        call c_nr_to_n0(nr,nr)  !   (10)


if(use_quaternion)then
qnr=nr
 AS=1 ; AS%q=exp(qnr)
else
      AS=1 ; AS%s=exp(nr)*AS%s 
endif


 
        n%AS=n%AS*AS             ! (12)
 

 
        m1=c_simil(AS,m1,-1) 
  

       enddo

      n%AS=from_phasor()*n%AS*from_phasor(-1)
 
      n%AS=n%A_t*n%AS*n%a_t**(-1)
 
      call kill(AS) 
      call kill(mt) 
      call kill(n0) 
      call kill(nr) 
      call kill(qnr) 
    endif
      

    n%n=c_simil(from_phasor(),m1,1)
    n%Atot=n%as*n%a_t

 


    if(present(rot)) then
      rot=n%Atot**(-1)*xy*n%Atot
    endif
    
    if(present(nu_spin)) nu_spin=0.0_dp
      
    if(present(phase)) then
      
      do i=1,size(phase)
         phase(i)=0.0_dp
      enddo
        
      if(present(rot)) then
        m1=rot
      else
        m1=n%Atot**(-1)*xy*n%Atot
      endif
          qphase=.false.
      call c_full_canonise(m1,a1,phase=phase,nu_spin=nu_spin)
       if(dospinr.and.present(nu_spin)) then
        if(real(nu_spin.sub.'0')<0) nu_spin=-nu_spin   ! 2018.11.01  to match phase advance
       endif
          qphase=.true.
    endif


    call c_check_rad(m1%e_ij,rad_in)
    if(rad_in) call c_normal_radiation(m1,n)

    call kill(m1);call kill(nonl);call kill(a1);call kill(a2);call kill(ri);

      deallocate(eg)
      deallocate(je)

    if(lielib_print(13)/=0) then
     close(mker)
     close(mdis)
     close(mdiss)
     close(mkers)
    endif


   
 end subroutine c_normal


 subroutine c_normal_spin_linear_quaternion(m_in,m_out,as,alpha) 
!#restricted: normal
!# This routine normalises the constant part of the spin matrix. 
!# m_out=as**(-1)*m_in*as
  implicit none
  type(c_damap), intent(inout) :: m_in,m_out,as
  type(quaternion) q0,q1,e_y,q3,qs
  real(dp) alpha,cosalpha,sinalpha

q0=m_in%q.sub.0

         as=1

q1=q0
q1%x(0)=0.0_dp
qs=1.0_dp/sqrt(q1%x(1)**2+q1%x(2)**2+q1%x(3)**2)
q1=q1*qs   ! q1=n

e_y=0.0_dp
e_y%x(2)=1.0_dp
 

q3=q1*e_y

 ! q3 =-n.j + n x j . l

cosalpha=-q3%x(0)

sinalpha=sqrt(q3%x(1)**2+q3%x(2)**2+q3%x(3)**2)



alpha= atan2(sinalpha,cosalpha)

 

if(alpha==0.and.cosalpha/=-1.0_dp) then
! write(6,*)sinalpha,cosalpha
! write(6,*) "weird in c_normal_spin_linear_quaternion "
! pause 123 
 q3=1.0_dp


else

if(abs(cosalpha+1.0_dp)<=1.e-16_dp)  then
 q3=-1.0_dp
else 
 q3%x(0)=cos(alpha/2)
 q3%x(1:3)=-sin(alpha/2)*q3%x(1:3)/sinalpha 

endif


endif

  

as%q=q3   
        m_out=c_simil(AS,m_in,-1)
q0=m_out%q

alpha=2*atan2(q0%x(2),q0%x(0))
 
 end  subroutine c_normal_spin_linear_quaternion

 



 
  subroutine c_normal_spin_linear(m_in,m_out,as,n0,as_ext) 
!#restricted: normal
!# This routine normalises the constant part of the spin matrix. 
!# m_out=as**(-1)*m_in*as
!# n0 is the invariant spin direction of type c_spinor.
!# Optional map as_ext is for a crazy Stern-Gerlach extension (not used).
  implicit none
  type(c_damap), intent(inout) :: m_in,m_out,as
  type(c_damap), optional :: as_ext
  type(c_spinor), intent (inout) :: n0
  type(c_taylor) s(3)
  integer i
!  as_ext is for crazy Stern-Gerlach extension
        as=1
        call c_find_n0(m_in%s,n0,linear=my_true)
        call c_find_as(n0,AS%s)

        m_out=c_simil(AS,m_in,-1)

    if(present(as_ext)) then
        as_ext=1
        as_ext%s=AS%s**(-1)
     do i=1,3
      call alloc(s(i))
     enddo
  
     s(1)=(1.0_dp.cmono.(nd2-1))
     s(3)=(1.0_dp.cmono.nd2)
     s(2)= sqrt(1.0_dp-s(1)**2-s(3)**2)

     as_ext%v(nd2-1)=0.0_dp
     as_ext%v(nd2)=0.0_dp
     do i=1,3
      as_ext%v(nd2-1)=as_ext%v(nd2-1)+ as_ext%s%s(1,i)*s(i)
      as_ext%v(nd2)  =as_ext%v(nd2)  + as_ext%s%s(3,i)*s(i)
     enddo      
     as_ext%s=0
     
     do i=1,3
      call kill(s(i))
     enddo

    endif 

 end  subroutine c_normal_spin_linear

 subroutine c_convert_spin(xy,m1) 
!#general: Stern-Gerlach
!# Bizarre routine to convert a spectator spin
!# into a Stern-Gerlach type map.
    implicit none
    type(c_damap) , intent(inout) :: xy,m1
    type(c_taylor) s(3)
    integer i
    
    do i=1,3
     call alloc(s(i))
    enddo
    m1=xy    
   
    s(1)=(1.0_dp.cmono.(nd2-1))
    s(3)=(1.0_dp.cmono.nd2)
    s(2)= sqrt(1.0_dp-s(1)**2-s(3)**2)
     
     m1%v(nd2-1)=0.0_dp
     m1%v(nd2)=0.0_dp
     do i=1,3
      m1%v(nd2-1)=m1%s%s(1,i)*s(i)+m1%v(nd2-1)
      m1%v(nd2)=m1%s%s(3,i)*s(i)+m1%v(nd2)
     enddo
     m1%s=0

    do i=1,3
     call kill(s(i))
    enddo

 end subroutine c_convert_spin


    function coast(l)
    implicit none 
    logical(lp) coast
    integer l

    coast=my_false
    if(ndpt/=0) then
     coast=l>nd2t.and.l<nd2-2*rf+1 
    endif
      
    end function coast

 subroutine c_normal_radiation(m1,n)
!#general:  normal
!# This routine normalises a beam envelope.
!# The map m1 enters already in phasors basis
!# therefore c_normal_radiation is unlikely to be used on its own.
    implicit none
    type(c_damap) , intent(inout) :: m1
    type(c_normal_form), intent(inout) ::  n
    complex(dp) r(6,6)
    integer i,j
   


    r=m1

    do i=1,6 ; do j=1,6;
     n%s_ijr(i,j)= 1.0_dp/(1.d0- r(i,i)*r(j,j))*m1%e_ij(i,j)
    enddo ;enddo;

     do i=1,3
      n%emittance(i)=abs(n%s_ijr(2*i-1,2*i))/2.0_dp
     enddo

   m1%e_ij= n%s_ijr  !using m1 to transform equilibrium beam sizes
   m1=c_simil(n%a_t,c_simil(from_phasor(),m1,1),1)
   
   n%s_ij0=m1%e_ij

   end  subroutine  c_normal_radiation 

 subroutine c_stochastic_kick(m,ait,ki,eps)
!#general:  stochastic kick
!# This routine creates a random kick 
!# It diagonalises the beam envelope kick by noting
!# that this object is dual to Lie operators
    implicit none
    type(c_damap) , intent(inout) :: m
    real(dp), intent(out) :: ait(6,6),ki(6)
    real(dp), intent(in) :: eps
    integer i,j
    real(dp) norm,normy,f(6,6),s(6,6), at(6,6),ai(6,6),m1(6,6)
    real(dp) b(6,6),a(6,6),d(6,6)
    type(c_vector_field) vf
    type(c_damap) id,as
 
    logical yhere


!!!! put stochastic kick in front per Sagan 
    m1=m**(-1)     !matmul(matmul(mat,ma%e_ij),transpose(mat))  not necessary I think
    m%e_ij=matmul(matmul(m1,m%e_ij),transpose(m1))
!!!!!!!!!!!!!!!!!!!

    norm=0.0_dp
    do i=1,6
    do j=1,6
     norm=abs(m%e_ij(i,j)) + norm
    enddo
    enddo
       if(norm==0) then
         ki=0
         ait=0
         return
      endif



     vf%n=m%n
    call alloc(vf)
    call alloc(id,as)
 


    norm=norm/10

    f=0
    s=0
    do i=1,3
     s(2*i-1,2*i)=1
     s(2*i,2*i-1)=-1
    enddo
    do i=1,6
    do j=1,6
     b(i,j)=m%e_ij(i,j)
     f(i,j)=m%e_ij(i,j)/norm
    enddo
    enddo

if(.not.use_new_stochastic_normal_form) then    
!!!! In rings without errors and middle plane symmetry
!!!  the y-part of the envelope is exactly zero
!!! So I add a y-part to make sure that my normal form
!!! does not crap out
    normy=0
    do i=1,6
    do j=1,6
     yhere=i==3.or.i==4.or.j==3.or.j==4
     if(yhere) normy=normy +abs(f(i,j))
    enddo
    enddo
     if (c_verbose) write(6,*) " norm y",norm

     if((normy/10)/norm<eps) then
       f(3,3)=0.234567_dp*twopi 
       f(4,4)=0.234567_dp*twopi 
     endif
endif
     f=matmul(f,s)
    do i=1,6
    do j=1,6
     vf%v(i)=vf%v(i)+ f(j,i) * (1.0_dp.cmono.j)
    enddo
    enddo

    id=exp(vf)

if(use_new_stochastic_normal_form) then
call  c_linear_a_stoch(id,as)
else
call c_linear_a(id,as)
endif 
    a=as
    at=transpose(a)
 if (c_verbose) then
    do i=1,6
          write(6,'(6(1x,G21.14))') a(i,1:6)
    enddo
 endif
!!!!  initially  !!!!
!sigma_f= M sigma M^t + B
!
!
!  D=  A^t B A 
    d=matmul(at,b)
    d=matmul(d,a)

! D = is diagonal
!  consider  a random variable r_i (i=1,6)  where <r_i>=0 and <r_i r_j>=delta_ij
! then  
  
    do i=1,6
     ki(i)=d(i,i)
    enddo
     do i=1,3
      if(ki(2*i-1)<0) then
       write(6,*) "fluctuations ill defined in plane ",i
        write(6,*) i,ki(2*i-1:2*i)
         ki(2*i-1)=0
         ki(2*i)=0
       endif
      if(ki(2*i)<0) then
       write(6,*) "fluctuations ill defined in plane ",i
        write(6,*) i,ki(2*i-1:2*i)
         ki(2*i-1)=0
         ki(2*i)=0
       endif
      if(ki(2*i-1)/=0.and.ki(2*i)/=0) then
        if(ki(2*i-1)-ki(2*i)/=0) then
          if(abs( (ki(2*i-1)-ki(2*i))/(abs(ki(2*i-1))+abs(ki(2*i))))>1.d-4) then
          write(6,*) "fluctuations ill defined in plane ",i
          write(6,*) i,ki(2*i-1:2*i)
           ki(2*i-1)=0
           ki(2*i)=0
          endif
        endif
      endif
     enddo 
    do i=1,6
     ki(i)=sqrt(ki(i))
    enddo
!  construct    z_i =ki(i) * r_i
! then x= Ait z  is the appropriate kick
! in the original space. 
! 

 if(global_verbose)  then
    write(6,*)" Stochastic kick"
    do i=1,6
    do j=1,6
    if(abs(d(i,j))>eps*norm) then
     write(6,*) i,j,d(i,j)
    endif
    enddo
    enddo
endif

    ai=-matmul(matmul(s,at),s)
    ait=transpose(ai)

!  B= ait*d*ai

    call kill(vf)
    call kill(id,as)
 

end subroutine c_stochastic_kick

    subroutine check_kernel(k,n,je,removeit)
!#internal: normal
!# This routine identifies terms in an orbital vector field that
!# are not in the kernel of a complete normalisation.
!# Namely it preserves generators of rotation.
!# Due to the vagueries of accelerator physics,
!# it leaves non-linear generators of rotations even in the
!# case of a damped oscillator.
    implicit none
    logical(lp) removeit
    integer i,k,n,je(:),t,o

    if(remove_tune_shift) then
      removeit=my_true
      t=0
      if(k/=0) je(k)=je(k)-1
      o=0
      do i=1,n,2
         if(coast(i)) cycle 
       o=o+je(i)+je(i+1)
       t=t+abs(je(i)-je(i+1))
      enddo
        if(t==0.and.o==0) removeit=my_false
      if(k/=0) je(k)=je(k)+1
     else
      removeit=my_true
      t=0
      if(k/=0) je(k)=je(k)-1
      do i=1,n,2
         if(coast(i)) cycle 
       t=t+abs(je(i)-je(i+1))
      enddo
        if(t==0) removeit=my_false
      if(k/=0) je(k)=je(k)+1
     endif
    end subroutine check_kernel


    subroutine check_kernel_spin(k,n,je,removeit)
!#internal: normal
!# This routine identifies terms in spin-orbital vector field that
!# are not in the kernel of a complete normalisation.

    implicit none
    logical(lp) removeit
    integer i,k,n,je(:),t
     
    removeit=my_true
    if(mod(k,2)/=0) return
    t=0
    do i=1,n,2
       if(coast(i)) cycle 
     t=t+abs(je(i)-je(i+1))
    enddo
  !    if(k==1) then
  !     t=t+iabs(1)   !  k is +/- 1 for spin s1 and s3    
  !    elseif(k==3) then
  !      t=t+iabs(-1)   !  k is +/- 1 for spin s1 and s3    
  !   endif
      if(t==0) removeit=my_false

    end subroutine check_kernel_spin

    subroutine check_resonance(k,n,je,kr,m,removeit)
!#internal: normal
!# This routine identifies terms in an orbital vector field that
!# are left per user's request.
!# This is used if a resonance family is to be left in the map.
!# See Sec.5.4 of Springer book.

    implicit none
    logical(lp) removeit
    integer i,k,n,je(:),t1,t2,j,kr,m(:,:)

     
    removeit=my_true
    t1=0;    t2=0;
    je(k)=je(k)-1
    do i=1,n,2
       if(coast(i)) cycle 
     j=(i+1)/2
     t1=t1+abs(je(i)-je(i+1)+m(j,kr))
     t2=t2+abs(je(i)-je(i+1)-m(j,kr))
    enddo
      if(t1==0.or.t2==0) removeit=my_false
    je(k)=je(k)+1

    end subroutine check_resonance

    subroutine check_resonance_spin(k,n,je,kr,ms,m,removeit) 
!#internal: normal
!# This routine identifies terms in a spin-orbital orbital vector field that
!# are left per user's request.
!# This is used if a spin resonance family is to be left in the map.
!# See Sec.6.4 of Springer book.
    implicit none
    logical(lp) removeit
    integer i,k,n,je(:),t1,t2,j,kr,ms(:),m(:,:)

     
    removeit=my_true
    t1=0;    t2=0;
 
    do i=1,n,2
       if(coast(i)) cycle 
     j=(i+1)/2
     t2=t2+abs(je(i)-je(i+1)-m(j,kr))
     t1=t1+abs(je(i)-je(i+1)+m(j,kr))
    enddo
!        if(k==1) then
!         t1=t1+iabs(-spin_def_tune-ms(kr))
!         t2=t2+iabs(-spin_def_tune+ms(kr))
!        elseif(k==3) then
!         t1=t1+iabs(spin_def_tune-ms(kr))
!         t2=t2+iabs(spin_def_tune+ms(kr))
!        else
!         t1=t1+iabs(ms(kr))
!         t2=t2+iabs(ms(kr))
!        endif
       if(k==1) then
        if(ms(kr)>0) then
!         t2=t2+iabs(spin_def_tune)
           if(t2==0) removeit=my_false
        elseif(ms(kr)<0) then
!         t1=t1+iabs(spin_def_tune)
           if(t1==0) removeit=my_false
        endif
        elseif(k==3) then
        if(ms(kr)>0) then
!         t1=t1+iabs(spin_def_tune)
           if(t1==0) removeit=my_false
        elseif(ms(kr)<0) then
!         t2=t2+iabs(spin_def_tune)
           if(t2==0) removeit=my_false
        endif
        else
         t1=t1+iabs(ms(kr))
         t2=t2+iabs(ms(kr))
          if(t1==0.or.t2==0) removeit=my_false
        endif



    end subroutine check_resonance_spin

   SUBROUTINE C_AVERAGE(F,A,F_FLOQUET) 
!#internal: normal
!# This routine averages a function F using the canonical trnsformation A.
!# This is explained in Sec.2.2.1 of my Springer book.
!# The result expressed in phasor's basis is in F_FLOQUET
    IMPLICIT NONE
    TYPE(c_damap) A
    TYPE(c_TAYLOR) F,F_FLOQUET,FQ
    complex(dp) v
    INTEGER I
    logical(lp) removeit
    integer, allocatable :: je(:)

    allocate(je(nv))  
         
     call alloc(fq)
    
    F_FLOQUET=(F*A)*TO_PHASOR(-1) 

    i=1
    do while(.true.) 
     call  c_cycle(F_FLOQUET,i,v ,je); if(i==0) exit;   
     call check_kernel(0,A%n,je,removeit) 
     if(.not.removeit) then
       fq=fq+(v.cmono.je)
     endif
    enddo
     
    f_floquet=fq
    call kill(fq)

    deallocate(je)

   END SUBROUTINE C_AVERAGE 


  function c_expflo_fac(h,x)   !,eps,nrmax)
    implicit none
    ! DOES EXP( \VEC{H} ) X = Y
    integer i,localmaster
    type(c_damap) c_expflo_fac
    type(c_damap),optional, intent(in):: x
    type(c_factored_lie), intent(in) :: h
    IF(.NOT.C_STABLE_DA) then
     c_expflo_fac%v%i=0
     RETURN
     endif

    localmaster=c_master
    if(present(x)) then
     c_expflo_fac%n=x%n
    else
     c_expflo_fac%n=nd2
    endif
    call c_assmap(c_expflo_fac)

    if(present(x)) then
     c_expflo_fac=x
    else
     c_expflo_fac=1
    endif

    if(h%dir==1) then
      do i=h%n,1,-1
       c_expflo_fac=texp(h%f(i),c_expflo_fac)
      enddo   
    else
      do i=1,h%n
       c_expflo_fac=texp(h%f(i),c_expflo_fac)
      enddo   
    endif
        if(present(x)) c_expflo_fac%x0=x%x0

    c_master=localmaster

   end function c_expflo_fac


  function c_expflo_fac_inv(h,x)   !,eps,nrmax)
    implicit none
    ! DOES EXP( \VEC{H} ) X = Y
    integer i,localmaster
    type(c_damap) c_expflo_fac_inv
    type(c_damap),optional, intent(in):: x
    type(c_factored_lie), intent(in) :: h
    IF(.NOT.C_STABLE_DA) then
     c_expflo_fac_inv%v%i=0
     RETURN
     endif

    localmaster=c_master
    if(present(x)) then
     c_expflo_fac_inv%n=x%n
    else
     c_expflo_fac_inv%n=nd2
    endif
    call c_assmap(c_expflo_fac_inv)

    if(present(x)) then
     c_expflo_fac_inv=x
    else
     c_expflo_fac_inv=1
    endif

    if(h%dir==-1) then
      do i=h%n,1,-1
       c_expflo_fac_inv=texp(-h%f(i),c_expflo_fac_inv)
      enddo   
    else
      do i=1,h%n
       c_expflo_fac_inv=texp(-h%f(i),c_expflo_fac_inv)
      enddo   
    endif

        if(present(x)) c_expflo_fac_inv%x0=x%x0
    c_master=localmaster

  end function c_expflo_fac_inv


      function c_add_map(s1,s2)   !,eps,nrmax)
    implicit none
 
    integer i,localmaster
    type(c_damap) c_add_map
    type(c_damap), intent(in):: s1
    type(c_damap), intent(in) :: s2
    type(c_taylor) tt
      localmaster=c_master

 
     c_add_map%n=s1%n
    call alloc(tt)
    call c_assmap(c_add_map)
    
    do i=1,s1%n
      tt=s1%v(i)+ s2%v(i) 
      c_add_map%v(i)=tt
    enddo
    c_add_map%s=s1%s+ s2%s 
    c_add_map%e_ij=s1%e_ij+ s2%e_ij 

    c_master=localmaster
 call kill(tt)
  end function c_add_map

      function c_sub_map(s1,s2)   !,eps,nrmax)
    implicit none
 
    integer i,localmaster
    type(c_damap) c_sub_map
    type(c_damap), intent(in):: s1
    type(c_damap), intent(in) :: s2
    type(c_taylor) tt

      localmaster=c_master

 
     c_sub_map%n=s1%n
    call alloc(tt)
    call c_assmap(c_sub_map)
    
    do i=1,s1%n
! etienne was completely wrong fixed 1/5/2016
      tt=s1%v(i)- s2%v(i) 
      c_sub_map%v(i)=tt
    enddo
    c_sub_map%s=s1%s- s2%s 
    c_sub_map%e_ij=s1%e_ij- s2%e_ij 

    c_master=localmaster
 call kill(tt)
  end function c_sub_map

 


      function c_1_vf_q(s1,c)   !,eps,nrmax)
    implicit none
  !  extracts mt-1 as a vector field including spin
    integer i,localmaster,ci,k
    type(c_vector_field) c_1_vf_q
    type(c_damap), intent(in):: s1
    type(c_taylor) tt
    real(dp) norm
    integer, optional :: c
      localmaster=c_master
     ci=1
     if(present(c)) ci=c
     c_1_vf_q%n=s1%n
    call alloc(tt)
    call c_ass_vector_field(c_1_vf_q)
    
    do i=1,s1%n
      tt=s1%v(i)-ci*(1.d0.cmono.i)
      c_1_vf_q%v(i)=tt
    enddo
    
    call c_full_norm_quaternion(s1%q,i,norm)
    if(i/=0) then
     do k=1,3
     c_1_vf_q%q%x(k)=s1%q%x(k)
    enddo
     c_1_vf_q%q%x(0)=0.0_dp  !s1%q%x(0)-1.0_dp

    endif
    
    c_master=localmaster
 call kill(tt)
  end function c_1_vf_q

      function c_1_map(s1,c)   !,eps,nrmax)
    implicit none
 
    integer i,j,localmaster,ci
    type(c_damap) c_1_map
    type(c_damap), intent(in):: s1
    type(c_taylor) tt
    real(dp) norm
    integer, optional :: c
      localmaster=c_master
     ci=-1
     if(present(c)) ci=c
     c_1_map%n=s1%n
    call alloc(tt)
    call c_assmap(c_1_map)
     c_1_map=s1
    do i=1,s1%n
      tt=s1%v(i)+ci*(1.d0.cmono.i)
      c_1_map%v(i)=tt
    enddo
    
    call c_full_norm_spin(s1%s,i,norm)
    if(i/=0) then
      do j=1,3
           c_1_map%s%s(j,j)=c_1_map%s%s(j,j)+ci 
      enddo
    endif
    
    c_master=localmaster
 call kill(tt)
  end function c_1_map


      function c_add_vf(s1,s2)   !,eps,nrmax)
    implicit none
 
    integer i,localmaster
    type(c_vector_field) c_add_vf
    type(c_vector_field), intent(in):: s1
    type(c_vector_field), intent(in) :: s2
    type(c_taylor) tt
      localmaster=c_master

 
     c_add_vf%n=s1%n
    call alloc(tt)
    call c_ass_vector_field(c_add_vf)
    
    do i=1,s1%n
      tt=s1%v(i)+ s2%v(i) 
      c_add_vf%v(i)=tt
    enddo
 
    c_add_vf%q=s1%q+s2%q 

     
    c_add_vf%nrmax=max(s1%nrmax,s2%nrmax)
    c_add_vf%eps=min(s1%eps,s2%eps)
    c_master=localmaster
 call kill(tt)
  end function c_add_vf

      function c_sub_vf(s1,s2)   !,eps,nrmax)
    implicit none
 
    integer i,localmaster
    type(c_vector_field) c_sub_vf
    type(c_vector_field), intent(in):: s1
    type(c_vector_field), intent(in) :: s2
    type(c_taylor) tt
      localmaster=c_master

 
     c_sub_vf%n=s1%n
    call alloc(tt)
    call c_ass_vector_field(c_sub_vf)
    
    do i=1,s1%n
      tt=s1%v(i)- s2%v(i) 
      c_sub_vf%v(i)=tt
    enddo
 
    c_sub_vf%q=s1%q-s2%q 

     
    c_sub_vf%nrmax=max(s1%nrmax,s2%nrmax)
    c_sub_vf%eps=min(s1%eps,s2%eps)
    c_master=localmaster
 call kill(tt)
  end function c_sub_vf


    FUNCTION real_mul_map( r,S1 )
    implicit none
    TYPE (c_damap) real_mul_map
    real(dp),intent(in):: r
    TYPE (c_damap), INTENT (IN) :: S1
    integer localmaster,i,j

    IF(.NOT.C_STABLE_DA) then
     real_mul_map%v%i=0
     RETURN
     endif
    localmaster=c_master

    !    call check(s1)
    real_mul_map%n=s1%n
    call c_assmap(real_mul_map)
    
    do i=1,s1%n
     real_mul_map%v(i)=r*s1%v(i)
    enddo


    do i=1,3
    do j=1,3
     real_mul_map%s%s(i,j)=r*s1%s%s(i,j)
    enddo
    enddo

    c_master=localmaster

  END FUNCTION real_mul_map

    FUNCTION real_mul_vec( r,S1 )
    implicit none
    TYPE (c_vector_field) real_mul_vec
    real(dp),intent(in):: r
    TYPE (c_vector_field), INTENT (IN) :: S1
    integer localmaster,i

    IF(.NOT.C_STABLE_DA) then
     real_mul_vec%v%i=0
     RETURN
     endif
    localmaster=c_master

    !    call check(s1)
    real_mul_vec%n=s1%n
    call c_ass_vector_field(real_mul_vec)
    
    do i=1,s1%n
     real_mul_vec%v(i)=r*s1%v(i)
    enddo

 

    do i=0,3
     real_mul_vec%q%x(i)=r*s1%q%x(i)
    enddo


    real_mul_vec%nrmax=s1%nrmax
    real_mul_vec%eps=s1%eps

    c_master=localmaster

  END FUNCTION real_mul_vec
  
    FUNCTION complex_mul_vec( r,S1 )
    implicit none
    TYPE (c_vector_field) complex_mul_vec
    complex(dp),intent(in):: r
    TYPE (c_vector_field), INTENT (IN) :: S1
    integer localmaster,i

    IF(.NOT.C_STABLE_DA) then
     complex_mul_vec%v%i=0
     RETURN
     endif
    localmaster=c_master

    !    call check(s1)
    complex_mul_vec%n=s1%n
    call c_ass_vector_field(complex_mul_vec)
    
    do i=1,s1%n
     complex_mul_vec%v(i)=r*s1%v(i)
    enddo

 

    do i=0,3
     complex_mul_vec%q%x(i)=r*s1%q%x(i)
    enddo

    complex_mul_vec%nrmax=s1%nrmax
    complex_mul_vec%eps=s1%eps

    c_master=localmaster

  END FUNCTION complex_mul_vec 

    FUNCTION map_mul_vec_q( r,S1 )
    implicit none
    TYPE (c_vector_field) map_mul_vec_q
    type(c_damap),intent(in):: r
    TYPE (c_vector_field), INTENT (IN) :: S1
    integer localmaster,i,k
    TYPE (c_vector_field) f_orb
    TYPE (c_quaternion) alpha_inv
    IF(.NOT.C_STABLE_DA) then
     map_mul_vec_q%v%i=0
     RETURN
     endif
    localmaster=c_master
    map_mul_vec_q%n=s1%n
    call c_ass_vector_field(map_mul_vec_q)
     f_orb%n=s1%n
     call alloc(f_orb)
     call alloc(alpha_inv)
     f_orb=s1


     alpha_inv=r%q**(-1)

     f_orb=map_mul_vec( r,f_orb )
     f_orb%q=0.0_dp
      map_mul_vec_q=f_orb

     do i=0,3
      f_orb%q%x(i)=f_orb*alpha_inv%x(i)
    enddo
    map_mul_vec_q%q=(alpha_inv*(s1%q*r)+f_orb%q)*r%q

    map_mul_vec_q%nrmax=s1%nrmax
    map_mul_vec_q%eps=s1%eps
    call kill(alpha_inv)
    call kill(f_orb)
    c_master=localmaster
  END FUNCTION map_mul_vec_q


    FUNCTION map_mul_vec( r,S1 )
    implicit none
    TYPE (c_vector_field) map_mul_vec
    type(c_damap),intent(in):: r
    type(c_damap) ri,r0
    TYPE (c_vector_field), INTENT (IN) :: S1
    integer localmaster,i,k

    IF(.NOT.C_STABLE_DA) then
     map_mul_vec%v%i=0
     RETURN
     endif
    localmaster=c_master
    call alloc(ri,r0)
     r0=r
     ri=r0**(-1)
 
    !    call check(s1)
    map_mul_vec%n=s1%n
    call c_ass_vector_field(map_mul_vec)

    map_mul_vec=0
   
    do i=1,s1%n
    do k=1,s1%n
     map_mul_vec%v(k)=s1%v(i)*(ri%v(k).d.i)+map_mul_vec%v(k)
    enddo
    enddo
    do k=1,s1%n
     map_mul_vec%v(k)=map_mul_vec%v(k)*r0
    enddo

    map_mul_vec%nrmax=s1%nrmax
    map_mul_vec%eps=s1%eps

    map_mul_vec%q=s1%q
    c_master=localmaster
    call kill(ri,r0)
  END FUNCTION map_mul_vec

    function exp_ad(h,x)  !  exp(Lie bracket)
    implicit none
    ! DOES EXP( \VEC{H} ) X = Y

    integer i,j,localmaster
    type(c_vector_field) exp_ad,ft
    type(c_vector_field), intent(in):: x
    type(c_vector_field), intent(in) :: h
    real(dp) prec,xnorm,r,xnorma 
    IF(.NOT.C_STABLE_DA) then
     exp_ad%v%i=0
     RETURN
     endif 

prec=1.d-8
    localmaster=c_master

 
     exp_ad%n=x%n
     ft%n=x%n
     call alloc(ft)
    call c_ass_vector_field(exp_ad)

     exp_ad=x
   !  c=1.0_dp
     ft=x
     
         xnorm=0.0_dp
          do j=1,ft%n
             r=full_abs(ft%v(j))
             xnorm=xnorm+r
          enddo
      
     do i=1,x%nrmax

         ft=(1.0_dp/i)*(h.lb.ft)
         exp_ad=exp_ad+ ft
         
          xnorma=0.0_dp
          do j=1,ft%n
             r=full_abs(ft%v(j))
             xnorma=xnorma+r
          enddo        
      
          if(xnorma>=xnorm.and.i > 20) then

             exit 
          endif    
          xnorm=xnorma      
     enddo
    
 


     call kill(ft)
     
    c_master=localmaster
 
  end function exp_ad


  function c_expflo_map(h,x)   !,eps,nrmax)
    implicit none
    ! DOES EXP( \VEC{H} ) X = Y

    integer i,j,localmaster
    type(c_damap) c_expflo_map
    type(c_damap), optional, intent(in):: x
    type(c_vector_field), intent(in) :: h

    IF(.NOT.C_STABLE_DA) then
     c_expflo_map%v%i=0
     RETURN
     endif 

    localmaster=c_master

    if(present(x)) then
     c_expflo_map%n=x%n
    else
     c_expflo_map%n=nd2
    endif
 
    call c_assmap(c_expflo_map)

    if(present(x)) then
     c_expflo_map=x
    else
     c_expflo_map=1
    endif


    do i=1,c_expflo_map%n
     c_expflo_map%v(i)=texp(h,c_expflo_map%v(i))
    enddo
 
   c_expflo_map%q=texp(h,c_expflo_map%q)
        if(present(x)) c_expflo_map%x0=x%x0    
c_master=localmaster
 
  end function c_expflo_map


  function c_expflo(h,x)   !,eps,nrmax)
    implicit none
    ! DOES EXP( \VEC{H} ) X = Y
    logical(lp) more
    integer i ,localmaster
    type(c_taylor) c_expflo
    type(c_taylor), intent(in):: x
    type(c_taylor) b1,b2,b3,b4
    type(c_vector_field), intent(in) :: h
    real(dp) coe,r,rbefore
    IF(.NOT.C_STABLE_DA) then
     c_expflo%i=0
     RETURN
     endif 

    localmaster=c_master

    call ass(c_expflo)

    call alloc(b1)
    call alloc(b2)
    call alloc(b3)
    call alloc(b4)
    b4=x
    b1=x

    more=.true.
    rbefore=1e30_dp
    do i=1,h%nrmax
       coe=1.0_dp/REAL(i,kind=DP)
 
        b2=coe*b1
 
!       call dacmu(b1,coe,b2)
        b1=h*b2
 
!       call daflo(h,b2,b1)
        b3=b1+b4
!       call daadd(b4,b1,b3)
       r=full_abs(b1)
 !      call daabs(b1,r)
       if(more) then
          if(r.gt.h%eps) then
             rbefore=r
             goto 100
          else
             rbefore=r
             more=.false.
          endif
       else
          if(r.ge.rbefore) then
             c_expflo=b3
            c_master=localmaster
!             call dacop(b3,y)
             call kill(b4)
             call kill(b3)
             call kill(b2)
             call kill(b1)
             return
          endif
          rbefore=r
       endif
100    continue
       b4=b3
!       call dacop(b3,b4)
       
    enddo
    if(lielib_print(2)==1) then
       write(6,'(a6,1x,G21.14,1x,a25)') ' NORM ',h%eps,' NEVER REACHED IN EXPFLO '
    endif
    c_expflo=b3
!    call dacop(b3,y)
    call kill(b4)
    call kill(b3)
    call kill(b2)
    call kill(b1)

    c_master=localmaster
    
  end function c_expflo
  
   subroutine c_flofacg(xy0,h,epso,n)
!#internal: manipulation
!# Use preferably function c_logf or its interface
    implicit none
    ! GENERAL ONE EXPONENT FACTORIZATION
    logical(lp) more
    integer i,k,nrmax,k1,n1,nn,n2
!    integer,dimension(:)::xy,h
!    integer,dimension(ndim2)::x,v,w,t,z
     type(c_damap), intent(inout) :: xy0
     type(c_vector_field), intent(inout) :: h
     type(c_vector_field) t,z,z3,e2,e3,e4,e5
     type(c_damap)  x,v,xy
     real(dp), optional :: epso
     integer, optional :: n
    integer,dimension(lnv)::jj
    real(dp) eps,epsone,r,xn,xnbefore,xnorm,xnorm1
    
    n1=1
    if(present(n)) n1=n
    
    if(.not.c_stable_da) return
    jj=0
    jj(1)=1
    !
    xy%n=xy0%n
    call alloc(xy);
    x%n=xy%n;v%n=xy%n; z%n=xy%n; t%n=xy%n; 
    call alloc(x);call alloc(v); call alloc(z);call alloc(t);
    z3%n=xy%n;call alloc(z3);   
    
        e2%n=xy%n;call alloc(e2);        
        e3%n=xy%n;call alloc(e3);   
        e4%n=xy%n;call alloc(e4);  
        e5%n=xy%n;call alloc(e5);   


    v=1
    do n2=2,2     ! put two for extra convergence
    do nn=1,n1   
      !  write(6,*) n2,nn
        if(n2==1) then
          do k1=1,xy0%n
           xy%v(k1)=v%v(k1)+ nn*((xy0%v(k1).sub.1)-v%v(k1))/n1
          enddo
        else
          do k1=1,xy0%n
           xy%v(k1)=(xy0%v(k1).sub.1) + nn*(xy0%v(k1)-(xy0%v(k1).sub.1))/n1
          enddo
        endif
          
    xnorm1=0.0_dp
    do i=1,xy%n
       r=full_abs(xy%v(i))
       xnorm1=xnorm1+r
    enddo
    if(present(epso)) then
     epsone= epso
    else
     epsone= xnorm1/epso_factor
    endif
    xnbefore=1e36_dp
    more=.false.
!    eps=1e-5_dp
 
    nrmax=1000
    xn=1e4_dp
    if(lielib_print(3)==1) write(6,*) "epsone,xnorm1 ",epsone,xnorm1


    if(epsone>0.0_dp) then  !epsone>zero
       do k=1,nrmax
       if(.not.(check_stable.and.C_STABLE_DA)) then
        write(6,*) " unstable in map logarithm "
        check_stable=.false.
        call kill(x)
        call kill(xy)
        call kill(e2)
        call kill(e3)
        call kill(e4)
        call kill(e5)
        call kill(z3)
        call kill(v)
        call kill(t)
        call kill(z)
        return
         return
       endif
          t=-h

          x=texp(t,xy)
          do k1=1,xy%n
           t%v(k1)=x%v(k1)-v%v(k1)
          enddo
!if(k==1) write(6,*)k, full_abs(t%v(1)),full_abs(t%v(2))


          if(xn.lt.epsone) then  !!! Quasi quadratic after some convergence obtained
             if(lielib_print(3)==1) then
                write(6,'(a28,i4,g21.14,1x,g21.14)') " Norm of   CBH  iteration # ",k,xn,epsone 
             endif
  
             do k1=1,xy%n
              e2%v(k1)=-0.5_dp*(t*t%v(k1))
             enddo
             do k1=1,xy%n
              e3%v(k1)=-0.5_dp*(e2*t%v(k1))-(1.0_dp/6.0_dp)*(t*e2%v(k1))
             enddo
         
          do k1=1,xy%n
           t%v(k1)=t%v(k1)+e2%v(k1)+e3%v(k1) ! t is now tau
          enddo
             
                z=h.lb.t   ! <h,tau>    
                e2=h.lb.z  ! <h,<h,tau>>
                e3=t.lb.z  ! <tau,<h,tau>>
             do k1=1,xy%n
               t%v(k1)=t%v(k1)+0.5_dp*z%v(k1) +(1.0_dp/12.0_dp)*(e2%v(k1)-e3%v(k1)) 
             enddo
                 e4=t.lb.e2  ! (tau,<h,<h,tau>>> (0)
              do k1=1,xy%n
               t%v(k1)=t%v(k1)-(1.0_dp/24.0_dp)*e4%v(k1) 
             enddo
            if(extra_terms_log) then     
                e4=h.lb.e2  ! <h,<h,<h,tau>>> (1)
                e5=h.lb.e4  
              do k1=1,xy%n
               t%v(k1)=t%v(k1)-(1.0_dp/720.0_dp)*e5%v(k1)
             enddo
              e5=t.lb.e4      !  (2)
              do k1=1,xy%n
               t%v(k1)=t%v(k1)+(1.0_dp/360.0_dp)*e5%v(k1)
             enddo
               e4=t.lb.e2  ! <t,<h,<h,tau>>>
               e5=h.lb.e4  ! <h,<t,<h,<h,tau>>>  (3)
              do k1=1,xy%n
               t%v(k1)=t%v(k1)-(1.0_dp/120.0_dp)*e5%v(k1)
             enddo
                e4=h.lb.e3  ! <h,<t,<h,tau>>>
               e5=t.lb.e4  ! <t,<h,<t,<h,tau>>>   (4)
              do k1=1,xy%n
               t%v(k1)=t%v(k1)+(1.0_dp/120.0_dp)*e5%v(k1)
             enddo
               e4=t.lb.e3  ! <t,<t,<h,tau>>>
               e5=h.lb.e4  ! <h,<t,<t,<h,tau>>>   (5)
              do k1=1,xy%n
               t%v(k1)=t%v(k1)-(1.0_dp/360.0_dp)*e5%v(k1)
             enddo
               e5=t.lb.e4  ! <t,<t,<t,<h,tau>>>   (6)
              do k1=1,xy%n
               t%v(k1)=t%v(k1)+(1.0_dp/720.0_dp)*e5%v(k1)
             enddo
           endif
          endif    !!!! end of quasi quadratic 

                   do k1=1,xy%n
                    h%v(k1)=t%v(k1)+h%v(k1) 
                   enddo  

          xnorm=0.0_dp
          do i=1,xy%n
             r=full_abs(t%v(i))
             xnorm=xnorm+r
          enddo
          xn=xnorm/xnorm1
          if(xn.ge.epsone.and.(lielib_print(3)==1)) then
             write(6,'(a28,i4,g21.14,1x,g21.14)') " Norm of linear iteration # ",k,xn,epsone
          endif
 !         if(xn.lt.eps.or.more) then
         if(xn.lt.epsone.or.more) then

             more=.true.
             if(xn.ge.xnbefore) goto 1000
             xnbefore=xn
          endif
       enddo

1000   continue
       ! write(6,*) "k",k   
 else  !epsone>zero

       do k=1,nint(abs(epsone))-1
         if(lielib_print(3)==1) then
              xnorm=0.0_dp
              do i=1,xy%n
                 r=full_abs(t%v(i))
                 xnorm=xnorm+r
              enddo
              xn=xnorm/xnorm1
              write(6,'(a20,g21.14)') " norm of correction ",xn
         endif
             t=-h
        
            x=texp(t,xy)

                   do k1=1,xy%n
                    t%v(k1)=x%v(k1)-v%v(k1) 
                   enddo  


                   do k1=1,xy%n
                    h%v(k1)=t%v(k1)+h%v(k1) 
                   enddo  

       enddo
    endif
 !   if(lielib_print(3)==1) WRITE(6,*) " K ", K,epsone
    if(lielib_print(3)==1) then
       write(6,'(a11,i4)') " iteration " , k-1
    endif


    enddo  ! nn

    enddo  ! n2

    call kill(x)
    call kill(xy)
    call kill(e2)
    call kill(e3)
    call kill(e4)
    call kill(e5)
    call kill(z3)
    call kill(v)
    call kill(t)
    call kill(z)
    return
  end subroutine c_flofacg

 

!!!!!!!!!!!!!!!!!!!!! spin

  subroutine c_find_n0(s0,n0,linear) 
!#general: manipulation & normal
!# This routine finds the invariant vector n0
!# if the rotation s0 of type c_spinmatrix.
!# If linear is true, only the constant part of 
!# s0 is used. Thus s0_constant =exp(theta n0.L) where n0 is a 
!# vector of length one. This is the n0 of the normal form.
!# If linear is false, then s0=exp(theta n0.L) but this is NOT
!# the invariant ISF n. See Sec.6.3.2. where I explained why this n0
!# cannot be the ISF vector n.
    implicit none
    type(c_spinmatrix),intent(inout) :: s0
    type(c_spinor), intent(inout) :: n0
    type(c_taylor)  norm0
    type(c_taylor)  det,detm
    type(c_spinmatrix) ss
    logical(lp), optional :: linear
    logical(lp) linear0
    integer i,is,j

    linear0=my_false
    if(present(linear)) linear0=linear 
    call alloc(det)
    call alloc(detm)
    call alloc(norm0)
    call alloc(ss)
    
!    ss=s0

    
   do i=1,3
       do j=1,3
       if(linear0) then  
         ss%s(i,j)=(s0%s(i,j).sub.'0')
       else
         ss%s(i,j)=s0%s(i,j)
       endif
       enddo
    enddo



    do i=1,3
       ss%s(i,i)=ss%s(i,i)-1.0_dp
    enddo

    det=(ss%s(2,2)*ss%s(3,3)-ss%s(2,3)*ss%s(3,2))
 
    is=1
    detm=(ss%s(1,1)*ss%s(3,3)-ss%s(1,3)*ss%s(3,1))
 
    if(abs(detm)>=abs(det)) then
       det=detm
       is=2
    endif
 
    detm=ss%s(1,1)*ss%s(2,2)-ss%s(1,2)*ss%s(2,1)
    if(abs(detm)>=abs(det)) then
       det=detm
       is=3
    endif
 

    n0%v(is)=1.0_dp
    if(is==1) then
       n0%v(2)=(-ss%s(3,3)*ss%s(2,1)+ss%s(2,3)*ss%s(3,1))/det
       n0%v(3)=(-ss%s(2,2)*ss%s(3,1)+ss%s(2,1)*ss%s(3,2))/det
    elseif(is==2) then
       n0%v(1)=(-ss%s(3,3)*ss%s(1,2)+ss%s(3,2)*ss%s(1,3))/det
       n0%v(3)=(-ss%s(1,1)*ss%s(3,2)+ss%s(1,2)*ss%s(3,1))/det
    else
       n0%v(1)=(-ss%s(2,2)*ss%s(1,3)+ss%s(2,3)*ss%s(1,2))/det
       n0%v(2)=(-ss%s(1,1)*ss%s(2,3)+ss%s(1,3)*ss%s(2,1))/det
    endif

     norm0=sqrt(n0%v(1)**2+n0%v(2)**2+n0%v(3)**2)

 

    do i=1,3
       n0%v(i)=n0%v(i)/norm0
    enddo


    call kill(det)
    call kill(detm)
    call kill(norm0)
    call kill(ss)

  end subroutine c_find_n0

  subroutine c_n0_to_nr(n0,nr)  
!#general: manipulation & normal
!# Puts an vector n0 into the eigenoperators of L_y
!# nt%v(2)=n0%v(2)
!# nt%v(1)=n0%v(1)-i_*n0%v(3) ! coefficient of  1/2(L_x + i L_z) 
!# nt%v(3)=n0%v(1)+i_*n0%v(3) ! coefficient of  1/2(L_x - i L_z)
!# The inverse routine is c_nr_to_n0

    implicit none
    TYPE(c_spinor), INTENT(INout) :: n0,nr
    TYPE(c_spinor) nt

    call alloc(nt)

     nt%v(2)=n0%v(2)
     nt%v(1)=n0%v(1)-i_*n0%v(3) ! coefficient of  1/2(L_x + i L_z) 
     nt%v(3)=n0%v(1)+i_*n0%v(3) ! coefficient of  1/2(L_x - i L_z)
     nr=nt   
    call kill(nt)

  end subroutine c_n0_to_nr

  subroutine c_nr_to_n0(nr,n0) 
!#general: manipulation & normal
!# Puts an eigenoperators vector nr into the cartesian basis
!# (L_x,L_y,L_z).
!# nt%v(2)=nr%v(2)
!# nt%v(1)=(nr%v(1)+nr%v(3))/2.0_dp    ! coefficient of L_x 
!# nt%v(3)=i_*(nr%v(1)-nr%v(3))/2.0_dp ! coefficient of L_z
!# The inverse routine is c_n0_to_nr.
    implicit none
    TYPE(c_spinor), INTENT(INout) :: n0,nr
    TYPE(c_spinor) nt

    call alloc(nt)

     nt%v(2)=nr%v(2)
     nt%v(1)=(nr%v(1)+nr%v(3))/2.0_dp    ! coefficient of L_x
     nt%v(3)=i_*(nr%v(1)-nr%v(3))/2.0_dp ! coefficient of L_z
     n0=nt 
    call kill(nt)

  end subroutine c_nr_to_n0

  subroutine c_find_om_da(S,om,n) 
!#restricted : normal
!# This routine finds om such that the c_spinmatrix S= exp(om.L).
!# However it does it in "n" steps if n is specified.
!# Otherwise it does it in no steps, no is the order of the complex
!# TPSA package. 
!# This routine works only if the matrix S does not have constant parts.
!# It is used in the normal form algorithm for the spin.
    implicit none
    TYPE(c_spinmatrix), INTENT(INout) :: S
    TYPE(c_spinor), INTENT(INout) :: om
    integer, optional :: n
    TYPE(c_spinmatrix) h
    integer i,nr
    TYPE(c_spinmatrix) dh,dhn,di
    complex(dp) c,cl
    real(dp) depsb,depsa
    logical pre
    call alloc(h)
    call alloc(dh)
    call alloc(dhn)
    call alloc(di)
    !  this only works with a da-map
    nr=no
      pre=present(n)
    if(present(n)) then 
      nr=n
    endif
    dh=S
    h=0
    dhn=1
    dh%s(1,1)=dh%s(1,1)-1.0_dp
    dh%s(2,2)=dh%s(2,2)-1.0_dp
    dh%s(3,3)=dh%s(3,3)-1.0_dp
    depsb=1.e38_dp
    c=1.0_dp
    do i=1,nr
       dhn=dhn*dh
       cl=c/i
       di=cl*dhn
       h=h+di
       c=-c
        call c_norm_spinmatrix(di,depsa)
       if(lielib_print(3)==1) write(6,*) i,depsa
       if(pre.and.i>nr/2) then
        if(depsa>=depsb) exit
           depsb=depsa
       endif
    enddo

    if(pre.and.i>nr-1) then
      write(6,*) " did not converged in c_find_om_da ",i
      write(6,*) " Norms ",depsb,depsa
      stop
     endif

    om%v(1)=h%s(3,2)
    om%v(2)=h%s(1,3)
    om%v(3)=h%s(2,1) 

    call kill(h)
    call kill(dh)
    call kill(dhn)
    call kill(di)
  end subroutine c_find_om_da

  subroutine c_find_as(n22,a) 
!#general : normal & manipulation
!# Find the c_spinmatrix "a" such that
!# e_y = (0,1,0)= a**(-1)*n0
!# because a*exp(theta n0.L)*a**(-1)= exp(theta (a*n0).L). See Sec.6.5.2.

    implicit none
    type(c_spinor), intent(inout) ::  n22 
    type(c_spinmatrix), intent(inout) ::  a 
    type(c_spinor)  n1 ,n3,n2 
    type(c_taylor)   s,n
    real(dp) x
    integer i,is

    call alloc(n1)
    call alloc(n2)
    call alloc(n3)
    call alloc(s,n)

    n2=n22

    x=n2.dot.n2

    x=sqrt(x)

    do i=1,3
     n2%v(i)=n2%v(i)/x
    enddo

    ! here we find smallest value of n2
    is=2
    if(abs(n2%v(1))< abs(n2%v(2))) is=1

    if(is==1) then
       if(abs(n2%v(3))<abs(n2%v(1))) is=3
    else
       if(abs(n2%v(3))<abs(n2%v(2))) is=3
    endif

    !  put n1 in along that value
    do i=1,3
       n1%v(i)=0.0_dp
    enddo
    n1%v(is)=1.0_dp

    s=n2%v(is)*n1%v(is)

    n=0.0_dp
    do i=1,3
       n1%v(i)=n1%v(i)-s*n2%v(i)
       n=n1%v(i)**2+n
    enddo
    do i=1,3
       n1%v(i)=n1%v(i)/sqrt(n)
    enddo

    n3%v(1)=n1%v(2)*n2%v(3)-n1%v(3)*n2%v(2)
    n3%v(2)=n1%v(3)*n2%v(1)-n1%v(1)*n2%v(3)
    n3%v(3)=n1%v(1)*n2%v(2)-n1%v(2)*n2%v(1)

    n=0.0_dp
    do i=1,3
       n=n3%v(i)**2+n
    enddo
    do i=1,3
       n3%v(i)=n3%v(i)/sqrt(n)
    enddo

 !   if(spin_normal_position==2) then
 if(abs(n1%v(1))>abs(n3%v(1))) then
       do i=1,3
          a%s(i,1)=n1%v(i)
          a%s(i,2)=n2%v(i)
          a%s(i,3)=n3%v(i)
       enddo
else
x=n3%v(1).sub.'0'

if(x<0) then
       do i=1,3
          a%s(i,1)=-n3%v(i)
          a%s(i,2)=n2%v(i)
          a%s(i,3)=n1%v(i)
       enddo
else
       do i=1,3
          a%s(i,1)=n3%v(i)
          a%s(i,2)=n2%v(i)
          a%s(i,3)=-n1%v(i)
       enddo
endif

endif



    call kill(n1)
    call kill(n2)
    call kill(n3)
    call kill(s,n)
  end subroutine c_find_as

  subroutine c_inv_as(m,mi)  
!#internal : manipulation
!# Trivially inverts an O(3) rotation by taking the transpose.
!# mi=m**(-1)=m^t.
    implicit none
    integer i,j
    TYPE (c_spinmatrix), intent (inout):: m,mi
    TYPE (c_spinmatrix) n

   call alloc(n)

    do i=1,3
       do j=1,3
          n%s(j,i)=m%s(i,j)
       enddo
    enddo

    do i=1,3
       do j=1,3
          mi%s(i,j)=n%s(i,j)
       enddo
    enddo

   call kill(n)

  end subroutine c_inv_as



  subroutine c_find_spin_angle(S,tune,radian)  
!#restricted : normal
!# Find the tune of a spin rotation S
!# which is a rotation around the y-axis.
    implicit none

    TYPE (c_spinmatrix), intent (inout):: S
    TYPE (c_taylor) tune    
    logical(lp),optional :: radian
    logical(lp) rad
     
     rad=my_true
      if(present(radian)) rad=radian
     tune=S%s(1,1)+i_*S%s(1,3)

     tune=-i_*log(tune) 

     if(.not.rad) tune=tune/twopi

  end subroutine c_find_spin_angle

  function c_log_spinmatrix(s,exact,n) ! spin routine
    implicit none
    TYPE(c_spinor) c_log_spinmatrix
    TYPE(c_spinmatrix), INTENT(INout) :: s
    type(c_taylor) tune 
    type(c_damap) as
    logical(lp), optional :: exact
    integer, optional :: n
    logical(lp) exa,useq
    integer localmaster
    real(dp) d
    integer k,n1
    IF(.NOT.C_STABLE_DA) then
     c_log_spinmatrix%v%i=0
     RETURN
     endif
useq=use_quaternion
use_quaternion=.false.
    exa=my_false
    n1=no
   if(present(exact)) exa=exact    
    localmaster=c_master

    call c_ass_spinor(c_log_spinmatrix)

    if(present(n)) then
      call c_find_om_da(s,c_log_spinmatrix,n)
    else    

     if(.not.exa) then
        call c_norm_spin(s,k,EPS=d)
        if(d<1.d-2) then
         n1=no+20000
         k=1
        endif
     else
      k=0
     endif
! write(6,*) " k, d ",k,d
!pause 
!etienne
    if(k==1) then
      call c_find_om_da(s,c_log_spinmatrix,n1)
    else
        call alloc(as)
        call alloc(tune)
        
        as=1

        call c_find_n0(s,c_log_spinmatrix)

        call c_find_as(c_log_spinmatrix,AS%s)

        AS%s= AS%s**(-1)*s*as%s


        call c_find_spin_angle(as%s,tune)

 
        c_log_spinmatrix=tune*c_log_spinmatrix


          call kill(as)
          call kill(tune)
     endif

  endif
         c_master=localmaster
use_quaternion=useq
  end function c_log_spinmatrix

!  now useless, use quaternion if needed
! function c_spinor_spinmatrix(h_axis,ds) ! spin routine
!    implicit none
!    TYPE(c_spinmatrix) c_spinor_spinmatrix
!    TYPE(c_spinmatrix), INTENT(IN) :: DS
!    TYPE(c_spinor), INTENT(IN) :: h_axis

 !   integer localmaster
 !   TYPE(c_spinmatrix) dh 
 !   real(dp) eps,norm1,norm2
 !
 !   logical check
 !   IF(.NOT.C_STABLE_DA) then
 !    c_spinor_spinmatrix%s%i=0
 !    RETURN
 !    endif
 !
 !    localmaster=c_master
 !
 !     call c_ass_spinmatrix(c_spinor_spinmatrix)
 !     call alloc(dh)
 !
 !   dh%s(3,1)=-h_axis%v(2)
 !   dh%s(2,1)=h_axis%v(3)
 !   dh%s(1,3)=h_axis%v(2)
 !   dh%s(3,2)=h_axis%v(1)
 !   dh%s(1,2)=-h_axis%v(3)
 !   dh%s(2,3)=-h_axis%v(1)
    
!    c_spinor_spinmatrix=dh*ds
!    Lie operator order
!     c_spinor_spinmatrix=ds*dh 




!    call kill(dh)
!
!     c_master=localmaster
!
!   end function c_spinor_spinmatrix

 

function c_vector_field_quaternion(h,ds) ! spin routine
    implicit none
    TYPE(c_quaternion) c_vector_field_quaternion
    TYPE(c_quaternion), INTENT(IN) :: DS
    TYPE(c_vector_field), INTENT(IN) :: h
    integer  nmax
    integer i,j,localmaster
    IF(.NOT.C_STABLE_DA) then
     c_vector_field_quaternion%x(1)%i=0
     RETURN
     endif

     localmaster=c_master

      call c_ass_quaternion(c_vector_field_quaternion)
      
      do i=0,3
        c_vector_field_quaternion%x(i)=h*ds%x(i)
      enddo
  ! order reversed for compositional map considerations
       c_vector_field_quaternion=c_vector_field_quaternion+ds*h%q
 
     c_master=localmaster

   end function c_vector_field_quaternion



  function c_exp_spinmatrix(h_axis,ds) ! spin routine
    implicit none
    TYPE(c_spinmatrix) c_exp_spinmatrix
    TYPE(c_spinmatrix),optional, INTENT(INout) :: DS
    TYPE(c_spinor), INTENT(IN) :: h_axis
    integer  nmax
    integer i,localmaster
    TYPE(c_spinmatrix) dh,dhn,dr,dst
    real(dp) eps,norm1,norm2
    complex(dp) c
    logical check
    IF(.NOT.C_STABLE_DA) then
     c_exp_spinmatrix%s%i=0
     RETURN
     endif

    localmaster=c_master

    call c_ass_spinmatrix(c_exp_spinmatrix)




    check=.true.
    eps=1.d-5
    nmax=1000

    call alloc(dh)
    call alloc(dhn)
    call alloc(dr)
 

    !  this  works with a  tpsa-map

 
     c_exp_spinmatrix=1
  
    dh=0
    dh%s(2,1)=h_axis%v(3)
    dh%s(1,3)=h_axis%v(2)
    dh%s(3,2)=h_axis%v(1)
    dh%s(1,2)=-h_axis%v(3)
    dh%s(3,1)=-h_axis%v(2)
    dh%s(2,3)=-h_axis%v(1)

    dhn=1
    c=1.0_dp
    norm1=mybig
    do i=1,nmax
       dhn=dhn*dh
       c=c/i

       dr=c_exp_spinmatrix
       c_exp_spinmatrix=c_exp_spinmatrix+c*dhn 

       dr=c_exp_spinmatrix+(-1.0_dp,0.0_dp)*dr

       call c_full_norm_spinmatrix(dr,norm2)
 
       if(check) then
          if(norm2<eps.and.i>10) then
             check=.false.
          endif
       else
          if(norm2>=norm1) exit
       endif
       norm1=norm2
    enddo

    if(i>nmax-10) then
       write(6,*) "no convergence in c_exp_spinmatrix, enter 0 to stop "
       read(5,*) norm1
       if(norm1==0)  stop 1066
    endif
    if(present(ds)) c_exp_spinmatrix=c_exp_spinmatrix*ds

    call kill(dh)
    call kill(dhn)
    call kill(dr)
 
    c_master=localmaster
  end   function c_exp_spinmatrix

  function c_exp_quaternion(h_axis,ds) ! spin routine
    implicit none
    TYPE(c_quaternion) c_exp_quaternion
    TYPE(c_quaternion),optional, INTENT(INout) :: DS
    TYPE(c_quaternion), INTENT(IN) :: h_axis
    integer  nmax
    integer i,localmaster,k
    TYPE(c_quaternion) dh,dhn,dr,dst
    real(dp) eps,norm1,norm2
    complex(dp) c
    logical check
    IF(.NOT.C_STABLE_DA) then
     c_exp_quaternion%x(1)%i=0
     RETURN
     endif

    localmaster=c_master

    call c_ass_quaternion(c_exp_quaternion)

    check=.true.
    eps=1.d-5
    nmax=1000

    call alloc(dh)
    call alloc(dhn)
    call alloc(dr)
 
     c_exp_quaternion=1.0_dp
  
    dh=h_axis
 

    dhn=1.0_dp
    c=1.0_dp
    norm1=mybig
    do i=1,nmax
       dhn=dhn*dh
       c=1.0_dp/i
       dhn=c*dhn

       dr=c_exp_quaternion

       c_exp_quaternion=c_exp_quaternion+dhn 

       dr=c_exp_quaternion+(-1.0_dp,0.0_dp)*dr

       call c_full_norm_quaternion(dr,k,norm2)


       if(check) then
          if(norm2<eps.and.i>10) then
             check=.false.
          endif
       else
          if(norm2>=norm1) exit
       endif
       norm1=norm2
    enddo

    if(i>nmax-10) then
       write(6,*) "no convergence in c_exp_quaternion, enter 0 to stop "
       read(5,*) norm1
       if(norm1==0)  stop 1066
    endif
    if(present(ds)) c_exp_quaternion=c_exp_quaternion*ds

    call kill(dh)
    call kill(dhn)
    call kill(dr)
 
    c_master=localmaster
  end   function c_exp_quaternion



 function c_exp_vectorfield_on_quaternion(h,ds) ! spin routine
    implicit none
    TYPE(c_quaternion) c_exp_vectorfield_on_quaternion
    TYPE(c_quaternion), INTENT(INout) :: DS
    TYPE(c_vector_field), INTENT(IN) :: h
    integer  nmax,k
    integer i,localmaster
    TYPE(c_quaternion) dh,dr
    real(dp) eps,norm1,norm2
    complex(dp) c
    logical check
    IF(.NOT.C_STABLE_DA) then
     c_exp_vectorfield_on_quaternion%x(1)%i=0
     RETURN
     endif

    localmaster=c_master

    call c_ass_quaternion(c_exp_vectorfield_on_quaternion)




    check=.true.
    eps=1.d-10
    nmax=1000

    call alloc(dh)
 
    call alloc(dr)
 

    !  this  works with a  tpsa-map

 
     c_exp_vectorfield_on_quaternion=ds
  

   ! if(present(ds)) then 
     dh=ds
   ! else
   !  dh=1
   !  dh=h.lb.dh
   ! endif

 
    c=1.0_dp
    norm1=mybig
    do i=1,nmax

       dh=h*dh
       c=c/i

       dr=c_exp_vectorfield_on_quaternion
       c_exp_vectorfield_on_quaternion=c_exp_vectorfield_on_quaternion+c*dh 

       dr=c_exp_vectorfield_on_quaternion+(-1.0_dp,0.0_dp)*dr

       call c_full_norm_quaternion(dr,k,norm2)
 
       if(check) then
          if(norm2<eps.and.i>10) then
             check=.false.
          endif
       else
         if(norm2>=norm1) exit
       endif
       norm1=norm2
    enddo

    if(i>nmax-10) then
       write(6,*) "no convergence in c_exp_vectorfield_on_quaternion, enter 0 to stop "
       read(5,*) norm1
       if(norm1==0)  stop 1066
    endif


    call kill(dh)
    call kill(dr)
 
    c_master=localmaster
  end   function c_exp_vectorfield_on_quaternion


  subroutine c_full_norm_damap(m,norm) ! spin routine
    implicit none
    TYPE(c_damap), INTENT(IN) :: m
    real(dp) norm
    integer i,j

    norm=0.0_dp
    do i=1,m%n
          norm=norm+full_abs( m%v(i) )
    enddo
    if(use_quaternion) then
       do i=0,3
                norm=norm+full_abs( m%q%x(i) )
       enddo
    else
    do i=1,3
       do j=1,3
          norm=norm+full_abs( m%s%s(i,j) )
       enddo
    enddo
    endif
  end subroutine c_full_norm_damap


  subroutine c_full_norm_spin_map(m,k,norm) ! spin routine
    implicit none
    TYPE(c_damap), INTENT(IN) :: m
    real(dp) norm
    integer i,k

    norm=0.0_dp
    if(use_quaternion) then
     call c_full_norm_quaternion(m%q,k,norm)
    else
     call c_full_norm_spin(m%s,k,norm)
    endif

  end subroutine c_full_norm_spin_map

  subroutine c_full_norm_spinmatrix(s,norm) ! spin routine
    implicit none
    TYPE(c_spinmatrix), INTENT(IN) :: S
    real(dp) norm
    integer i,j

    norm=0.0_dp

    do i=1,3
       do j=1,3
          norm=norm+full_abs( s%s(i,j) )
       enddo
    enddo

  end subroutine c_full_norm_spinmatrix

  subroutine c_full_norm_quaternion(q,k,norm) ! spin routine
    implicit none
    TYPE(c_quaternion), INTENT(IN) :: q
    real(dp) norm,normn,nr
    integer i,k
    k=-1
    norm=0.0_dp

    do i=1,3
 
          norm=norm+full_abs( q%x(i) )
        
    enddo
       normn=norm
       norm=norm+full_abs( q%x(0) )

    if(norm==0) k=0
     nr=q%x(0) 
    if(norm==1.and.normn==0.and.nr==1) k=1

  end subroutine c_full_norm_quaternion

  subroutine c_norm_spinmatrix(s,norm) ! spin routine
    implicit none
    TYPE(c_spinmatrix), INTENT(IN) :: S
    real(dp) norm
    integer i,j

    norm=0.0_dp

    do i=1,3
       do j=1,3
          norm=norm+abs( s%s(i,j) )
       enddo
    enddo

  end subroutine c_norm_spinmatrix

  subroutine c_full_norm_vector_field(s,norm) ! spin routine
    implicit none
    TYPE(c_vector_field), INTENT(IN) :: S
    real(dp) norm,norms
    integer i,k

    norm=0.0_dp

    do i=1,s%n
          norm=norm+full_abs( s%v(i) )
    enddo
    norms=0.0_dp
    if(use_quaternion) then
     call c_full_norm_quaternion(s%q,k,norms) 
    endif
    norm=norm+norms
  end subroutine c_full_norm_vector_field

  subroutine c_full_norm_spinor(s,norm) ! spin routine
    implicit none
    TYPE(c_spinor), INTENT(IN) :: S
    real(dp) norm
    integer i

    norm=0.0_dp

    do i=1,3
          norm=norm+full_abs( s%v(i) )
    enddo

  end subroutine c_full_norm_spinor

  subroutine c_full_norm_fourier(s,norm) ! spin routine
    implicit none
    TYPE(c_vector_field_fourier), INTENT(IN) :: S
    real(dp) norm,nor1
    integer i

    norm=0.0_dp

    do i=-n_fourier,n_fourier
         call c_full_norm_vector_field(s%f(i),nor1)
          norm=norm+nor1
    enddo

  end subroutine c_full_norm_fourier

!!!! radiation 

 subroutine c_check_rad(e_ij,rad_in)
    implicit none
    logical(lp), INTENT(INout) :: rad_in
    complex(dp) e_ij(6,6)
    integer i,j
    real(dp) norm

    rad_in=my_true
    norm=0.0_dp
    do i=1,6
       do j=1,6
          norm=norm+abs(e_ij(i,j))
       enddo
    enddo

    if(norm==0.0_dp) then
        rad_in=.false.
    endif

  end  subroutine c_check_rad

!!!!!!!!!!!!!!!!!!!!!!! eig6


 subroutine c_eig6(fm,reval,aieval,revec,aievec)
!*
    implicit none
    !**************************************************************************

    !  Diagonalization routines of NERI

    !ccccccccccccccccc
    !
    !  this routine finds the eigenvalues and eigenvectors
    !  of the full matrix fm.
    !  the eigenvectors are normalized so that the real and
    !  imaginary part of vectors 1, 3, and 5 have +1 antisymmetric
    !  product:
    !      revec1 J aivec1 = 1 ; revec3 J aivec3 = 1 ;
    !      revec5 J aivec5 = one
    !  the eigenvectors 2 ,4, and 6 have the opposite normalization.
    !  written by F. Neri, Feb 26 1986.
    !
    integer jet,nn,i,i1,ilo,ihi,mdim,info
    real(dp),dimension(ndim2t)::reval,aieval,ort
    real(dp),dimension(ndim2t,ndim2t)::revec,aievec,fm,aa,vv

    if(.not.c_%stable_da) return


    !  copy matrix to temporary storage (the matrix aa is destroyed)
    do i=1,nd2harm  !nd2t !-ndc2t
       do i1=1,nd2harm  !,nd2t !-ndc2t
          aa(i1,i) = fm(i1,i)
       enddo
    enddo
    
    ilo = 1
    ihi = nd2harm  !nd2t !-ndc2t
    mdim = ndim2t
    nn = nd2harm   !nd2t !-ndc2t
    !  compute eigenvalues and eigenvectors using double
    !  precision Eispack routines:
    call ety(mdim,nn,ilo,ihi,aa,ort)
    call etyt(mdim,nn,ilo,ihi,aa,ort,vv)
    call ety2(mdim,nn,ilo,ihi,aa,reval,aieval,vv,info)
    
    if ( info .ne. 0 ) then
       write(6,*) '  ERROR IN C_EIG6'
       stop 999
    endif
    !      call neigv(vv,pbkt)
    
    do i=1,ndharm    !ndt !-ndct
       do jet=1,nd2harm   !nd2t !-ndc2t
          revec(jet,2*i-1)=vv(jet,2*i-1)
          revec(jet,2*i)=vv(jet,2*i-1)
          aievec(jet,2*i-1)=vv(jet,2*i)
          aievec(jet,2*i)=-vv(jet,2*i)
       enddo
    enddo
    do i=1,nd2harm    !nd2t !-ndc2t
       if(abs(reval(i)**2+aieval(i)**2 -1.0_dp).gt.1e-10_dp) then
           if(lielib_print(4)==1) then
             write(6,*) ' EIG6: Eigenvalues off the unit circle!'
             write(6,*) sqrt(reval(i)**2+aieval(i)**2)
          endif
       endif
    enddo
    return
  end subroutine c_eig6

  subroutine ety(nm,n,low,igh,a,ort)
!*
    implicit none
    !
    !     this subroutine is a translation of the algol procedure orthes,
    !     num. math. 12, 349-368(1968) by martin and wilkinson.
    !     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
    !
    !     given a real general matrix, this subroutine
    !     reduces a submatrix situated in rows and columns
    !     low through igh to upper hessenberg form by
    !     orthogonal similarity transformations.
    !
    !     on input-
    !
    !        nm must be set to the row dimension of two-dimensional
    !          array parameters as declared in the calling program
    !          dimension statement,
    !
    !        n is the order of the matrix,
    !
    !        low and igh are integers determined by the balancing
    !          subroutine  balanc.  if  balanc  has not been used,
    !          set low=1, igh=n,
    !
    !        a contains the input matrix.
    !
    !     on output-
    !
    !        a contains the hessenberg matrix.  information about
    !          the orthogonal transformations used in the reduction
    !          is stored in the remaining triangle under the
    !          hessenberg matrix,
    !
    !        ort contains further information about the transformations.
    !          only elements low through igh are used.
    !
    !     fortran routine by b. s. garbow
    !     modified by filippo neri.
    !
    !
    integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
    real(dp),dimension(:,:)::a
    real(dp),dimension(:)::ort
     
    real(dp) f,g,h,scale
    if(.not.c_stable_da) return

    la = igh - 1
    kp1 = low + 1
    if (la .lt. kp1) go to 200
    !
    do m = kp1, la
       h = 0.0_dp
       ort(m) = 0.0_dp
       scale = 0.0_dp
       !     ********** scale column (algol tol then not needed) **********
       do i = m, igh
          scale = scale + abs(a(i,m-1))
       enddo
       !
       if (scale .eq. 0.0_dp) go to 180
       mp = m + igh
       !     ********** for i=igh step -1 until m do -- **********
       do ii = m, igh
          i = mp - ii
          ort(i) = a(i,m-1) / scale
          h = h + ort(i) * ort(i)
       enddo
       !
       g = -sign(SQRT(h),ort(m))
       h = h - ort(m) * g
       ort(m) = ort(m) - g
       !     ********** form (i-(u*ut)/h) * a **********
       do j = m, n
          f = 0.0_dp
          !     ********** for i=igh step -1 until m do -- **********
          do ii = m, igh
             i = mp - ii
             f = f + ort(i) * a(i,j)
          enddo
          !
          f = f / h
          !
          do i = m, igh
             a(i,j) = a(i,j) - f * ort(i)
          enddo
          !
       enddo
       !     ********** form (i-(u*ut)/h)*a*(i-(u*ut)/h) **********
       do i = 1, igh
          f = 0.0_dp
          !     ********** for j=igh step -1 until m do -- **********
          do jj = m, igh
             j = mp - jj
             f = f + ort(j) * a(i,j)
          enddo
          !
          f = f / h
          !
          do j = m, igh
             a(i,j) = a(i,j) - f * ort(j)
          enddo
          !
       enddo
       !
       ort(m) = scale * ort(m)
       a(m,m-1) = scale * g
180    continue
    enddo
    !
200 return
    !     ********** last card of ety **********
  end subroutine ety

  subroutine etyt(nm,n,low,igh,a,ort,z)
!*
    implicit none
    !
    !     this subroutine is a translation of the algol procedure ortrans,
    !     num. math. 16, 181-204(1970) by peters and wilkinson.
    !     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
    !
    !     this subroutine accumulates the orthogonal similarity
    !     transformations used in the reduction of a real general
    !     matrix to upper hessenberg form by  ety.
    !
    !     on input-
    !
    !        nm must be set to the row dimension of two-dimensional
    !          array parameters as declared in the calling program
    !          dimension statement,
    !
    !        n is the order of the matrix,
    !
    !        low and igh are integers determined by the balancing
    !          subroutine  balanc.  if  balanc  has not been used,
    !          set low=1, igh=n,
    !
    !        a contains information about the orthogonal trans-
    !          formations used in the reduction by  orthes
    !          in its strict lower triangle,
    !
    !          ort contains further information about the trans-
    !          formations used in the reduction by  ety.
    !          only elements low through igh are used.
    !
    !     on output-
    !
    !        z contains the transformation matrix produced in the
    !          reduction by  ety,
    !
    !        ort has been altered.
    !
    !     fortran routine by b. s. garbow.
    !     modified by f. neri.
    !
    !
    integer i,j,n,kl,mm,mp,nm,igh,low,mp1
    real(dp) g
    real(dp),dimension(:)::ort
    real(dp),dimension(:,:)::a
    real(dp),dimension(:,:)::z
    if(.not.c_stable_da) return

    !     ********** initialize z to identity matrix **********
    do i = 1, n
       !
       do j = 1, n
          z(i,j) = 0.0_dp
       enddo
       !
       z(i,i) = 1.0_dp
    enddo
    !
    kl = igh - low - 1
    if (kl .lt. 1) go to 200
    !     ********** for mp=igh-1 step -1 until low+1 do -- **********
    do mm = 1, kl
       mp = igh - mm
       if (a(mp,mp-1) .eq. 0.0_dp) go to 140
       mp1 = mp + 1
       !
       do i = mp1, igh
          ort(i) = a(i,mp-1)
       enddo
       !
       do j = mp, igh
          g = 0.0_dp
          !
          do i = mp, igh
             g = g + ort(i) * z(i,j)
          enddo
          !     ********** divisor below is negative of h formed in orthes.
          !                double division avoids possible underflow **********
          g = (g / ort(mp)) / a(mp,mp-1)
          !
          do i = mp, igh
             z(i,j) = z(i,j) + g * ort(i)
          enddo
          !
       enddo
       !
140    continue
    enddo
    !
200 return
    !     ********** last card of etyt **********
  end subroutine etyt

  subroutine ety2(nm,n,low,igh,h,wr,wi,z,ierr)
!*
    implicit none
    !
    !
    !
    !     this subroutine is a translation of the algol procedure hqr2,
    !     num. math. 16, 181-204(1970) by peters and wilkinson.
    !     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
    !
    !     this subroutine finds the eigenvalues and eigenvectors
    !     of a real upper hessenberg matrix by the qr method.  the
    !     eigenvectors of a real general matrix can also be found
    !     if  elmhes  and  eltran  or  orthes  and  ortran  have
    !     been used to reduce this general matrix to hessenberg form
    !     and to accumulate the similarity transformations.
    !
    !     on input-
    !
    !        nm must be set to the row dimension of two-dimensional
    !          array parameters as declared in the calling program
    !          dimension statement,
    !
    !        n is the order of the matrix,
    !
    !        low and igh are integers determined by the balancing
    !          subroutine  balanc.  if  balanc  has not been used,
    !          set low=1, igh=n,
    !
    !        h contains the upper hessenberg matrix,
    !
    !        z contains the transformation matrix produced by  eltran
    !          after the reduction by  elmhes, or by  ortran  after the
    !          reduction by  orthes, if performed.  if the eigenvectors
    !          of the hessenberg matrix are desired, z must contain the
    !          identity matrix.
    !
    !     on output-
    !
    !        h has been destroyed,
    !
    !        wr and wi contain the real and imaginary parts,
    !          respectively, of the eigenvalues.  the eigenvalues
    !          are unordered except that complex conjugate pairs
    !          of values appear consecutively with the eigenvalue
    !          having the positive imaginary part first.  if an
    !          error exit is made, the eigenvalues should be correct
    !          for indices ierr+1,...,n,
    !
    !        z contains the real and imaginary parts of the eigenvectors.
    !          if the i-th eigenvalue is real, the i-th column of z
    !          contains its eigenvector.  if the i-th eigenvalue is complex
    !          with positive imaginary part, the i-th and (i+1)-th
    !          columns of z contain the real and imaginary parts of its
    !          eigenvector.  the eigenvectors are unnormalized.  if an
    !          error exit is made, none of the eigenvectors has been found,
    !
    !        ierr is set to
    !          zero       for normal return,
    !          j          if the j-th eigenvalue has not been
    !                     determined after 200 iterations.
    !
    !     arithmetic is real(dp). complex division
    !     is simulated by routin etdiv.
    !
    !     fortran routine by b. s. garbow.
    !     modified by f. neri.
    !
    !
    logical(lp) notlas
    integer i,j,k,l,m,n,en,ii,jj,ll,mm,na,nm,nn,igh,its,low,mp2,enm2,ierr
    real(dp) p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,z3r,z3i
    real(dp),dimension(:)::wr,wi
    real(dp),dimension(:,:)::h,z
    if(.not.c_stable_da) return

    !     ********** machep is a machine dependent parameter specifying
    !                the relative precision of floating point arithmetic.
    !
    !                **********
    !     machep = r1mach(4)
    !
    ierr = 0
    norm = 0.0_dp
    k = 1
    !     ********** store roots isolated by balanc
    !                and compute matrix norm **********
    do i = 1, n
       !
       do j = k, n
          norm = norm + abs(h(i,j))
       enddo
       !
       k = i
       if (i .ge. low .and. i .le. igh) go to 50
       wr(i) = h(i,i)
       wi(i) = 0.0_dp
50     continue
    enddo
    !
    en = igh
    t = 0.0_dp
    !     ********** search for next eigenvalues **********
60  if (en .lt. low) go to 340
    its = 0
    na = en - 1
    enm2 = na - 1
    !     ********** look for single small sub-diagonal element
    !                for l=en step -1 until low do -- **********
70  do ll = low, en
       l = en + low - ll
       if (l .eq. low) go to 100
       s = abs(h(l-1,l-1)) + abs(h(l,l))
       if (s .eq. 0.0_dp) s = norm
       if (abs(h(l,l-1)) .le. machep * s) go to 100
    enddo
    !     ********** form shift **********
100 x = h(en,en)
    if (l .eq. en) go to 270
    y = h(na,na)
    w = h(en,na) * h(na,en)
    if (l .eq. na) go to 280
    if (its .eq. 200) go to 1000
    if (its .ne. 10 .and. its .ne. 20) go to 130
    !     ********** form exceptional shift **********
    t = t + x
    !
    do i = low, en
       h(i,i) = h(i,i) - x
    enddo
    !
    s = abs(h(en,na)) + abs(h(na,enm2))
    x = 0.75_dp * s
    y = x
    w = -0.4375_dp * s * s
130 its = its + 1
    !     ********** look for two consecutive small
    !                sub-diagonal elements.
    !                for m=en-2 step -1 until l do -- **********
    do mm = l, enm2
       m = enm2 + l - mm
       zz = h(m,m)
       r = x - zz
       s = y - zz
       p = (r * s - w) / h(m+1,m) + h(m,m+1)
       q = h(m+1,m+1) - zz - r - s
       r = h(m+2,m+1)
       s = abs(p) + abs(q) + abs(r)
       p = p / s
       q = q / s
       r = r / s
       if (m .eq. l) go to 150
       if (abs(h(m,m-1)) * (abs(q) + abs(r)) .le. machep * abs(p) * (abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1)))) go to 150
    enddo
    !
150 mp2 = m + 2
    !
    do i = mp2, en
       h(i,i-2) = 0.0_dp
       if (i .eq. mp2) go to 160
       h(i,i-3) = 0.0_dp
160    continue
    enddo
    !     ********** double qr step involving rows l to en and
    !                columns m to en **********
    do k = m, na
       notlas = k .ne. na
       if (k .eq. m) go to 170
       p = h(k,k-1)
       q = h(k+1,k-1)
       r = 0.0_dp
       if (notlas) r = h(k+2,k-1)
       x = abs(p) + abs(q) + abs(r)
       if (x .eq. 0.0_dp) go to 260
       p = p / x
       q = q / x
       r = r / x
170    s = sign(SQRT(p*p+q*q+r*r),p)
       if (k .eq. m) go to 180
       h(k,k-1) = -s * x
       go to 190
180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
190    p = p + s
       x = p / s
       y = q / s
       zz = r / s
       q = q / p
       r = r / p
       !     ********** row modification **********
       do j = k, n
          p = h(k,j) + q * h(k+1,j)
          if (.not. notlas) go to 200
          p = p + r * h(k+2,j)
          h(k+2,j) = h(k+2,j) - p * zz
200       h(k+1,j) = h(k+1,j) - p * y
          h(k,j) = h(k,j) - p * x
       enddo
       !
       j = min0(en,k+3)
       !     ********** column modification **********
       do i = 1, j
          p = x * h(i,k) + y * h(i,k+1)
          if (.not. notlas) go to 220
          p = p + zz * h(i,k+2)
          h(i,k+2) = h(i,k+2) - p * r
220       h(i,k+1) = h(i,k+1) - p * q
          h(i,k) = h(i,k) - p
       enddo
       !     ********** accumulate transformations **********
       do i = low, igh
          p = x * z(i,k) + y * z(i,k+1)
          if (.not. notlas) go to 240
          p = p + zz * z(i,k+2)
          z(i,k+2) = z(i,k+2) - p * r
240       z(i,k+1) = z(i,k+1) - p * q
          z(i,k) = z(i,k) - p
       enddo
       !
260    continue
    enddo
    !
    go to 70
    !     ********** one root found **********
270 h(en,en) = x + t
    wr(en) = h(en,en)
    wi(en) = 0.0_dp
    en = na
    go to 60
    !     ********** two roots found **********
280 p = (y - x) / 2.0_dp
    q = p * p + w
    zz = SQRT(abs(q))
    h(en,en) = x + t
    x = h(en,en)
    h(na,na) = y + t
    if (q .lt. 0.0_dp) go to 320
    !     ********** real pair **********
    zz = p + sign(zz,p)
    wr(na) = x + zz
    wr(en) = wr(na)
    if (zz .ne. 0.0_dp) wr(en) = x - w / zz
    wi(na) = 0.0_dp
    wi(en) = 0.0_dp
    x = h(en,na)
    s = abs(x) + abs(zz)
    p = x / s
    q = zz / s
    r = SQRT(p*p+q*q)
    p = p / r
    q = q / r
    !     ********** row modification **********
    do j = na, n
       zz = h(na,j)
       h(na,j) = q * zz + p * h(en,j)
       h(en,j) = q * h(en,j) - p * zz
    enddo
    !     ********** column modification **********
    do i = 1, en
       zz = h(i,na)
       h(i,na) = q * zz + p * h(i,en)
       h(i,en) = q * h(i,en) - p * zz
    enddo
    !     ********** accumulate transformations **********
    do i = low, igh
       zz = z(i,na)
       z(i,na) = q * zz + p * z(i,en)
       z(i,en) = q * z(i,en) - p * zz
    enddo
    !
    go to 330
    !     ********** complex pair **********
320 wr(na) = x + p
    wr(en) = x + p
    wi(na) = zz
    wi(en) = -zz
330 en = enm2
    go to 60
    !     ********** all roots found.  backsubstitute to find
    !                vectors of upper triangular form **********
340 if (norm .eq. 0.0_dp) go to 1001
    !     ********** for en=n step -1 until 1 do -- **********
    do nn = 1, n
       en = n + 1 - nn
       p = wr(en)
       q = wi(en)
       na = en - 1
       if (q.lt.0) goto 710
       if (q.eq.0) goto 600
       if (q.gt.0) goto 800
       !     ********** real vector **********
600    m = en
       h(en,en) = 1.0_dp
       if (na .eq. 0) go to 800
       !     ********** for i=en-1 step -1 until 1 do -- **********
       do ii = 1, na
          i = en - ii
          w = h(i,i) - p
          r = h(i,en)
          if (m .gt. na) go to 620
          !
          do j = m, na
             r = r + h(i,j) * h(j,en)
          enddo
          !
620       if (wi(i) .ge. 0.0_dp) go to 630
          zz = w
          s = r
          go to 700
630       m = i
          if (wi(i) .ne. 0.0_dp) go to 640
          t = w
          if (w .eq. 0.0_dp) t = machep * norm
          h(i,en) = -r / t
          go to 700
          !     ********** solve real equations **********
640       x = h(i,i+1)
          y = h(i+1,i)
          q = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i)
          t = (x * s - zz * r) / q
          h(i,en) = t
          if (abs(x) .le. abs(zz)) go to 650
          h(i+1,en) = (-r - w * t) / x
          go to 700
650       h(i+1,en) = (-s - y * t) / zz
700       continue
       enddo
       !     ********** end real vector **********
       go to 800
       !     ********** complex vector **********
710    m = na
       !     ********** last vector component chosen imaginary so that
       !                eigenvector matrix is triangular **********
       if (abs(h(en,na)) .le. abs(h(na,en))) go to 720
       h(na,na) = q / h(en,na)
       h(na,en) = -(h(en,en) - p) / h(en,na)
       go to 730
       ! 720    z3 = cmplx(zero,-h(na,en)) / cmplx(h(na,na)-p,q)
       !        h(na,na) = real(z3,kind=dp)
       !        h(na,en) = aimag(z3)
720    call etdiv(z3r,z3i,0.0_dp,-h(na,en),h(na,na)-p,q)
       h(na,na) = z3r
       h(na,en) = z3i
730    h(en,na) = 0.0_dp
       h(en,en) = 1.0_dp
       enm2 = na - 1
       if (enm2 .eq. 0) go to 800
       !     ********** for i=en-2 step -1 until 1 do -- **********
       do ii = 1, enm2
          i = na - ii
          w = h(i,i) - p
          ra = 0.0_dp
          sa = h(i,en)
          !
          do j = m, na
             ra = ra + h(i,j) * h(j,na)
             sa = sa + h(i,j) * h(j,en)
          enddo
          !
          if (wi(i) .ge. 0.0_dp) go to 770
          zz = w
          r = ra
          s = sa
          go to 790
770       m = i
          if (wi(i) .ne. 0.0_dp) go to 780
          !           z3 = cmplx(-ra,-sa) / cmplx(w,q)
          !           h(i,na) = real(z3,kind=dp)
          !           h(i,en) = aimag(z3)
          call etdiv(z3r,z3i,-ra,-sa,w,q)
          h(i,na) = z3r
          h(i,en) = z3i
          go to 790
          !     ********** solve complex equations **********
780       x = h(i,i+1)
          y = h(i+1,i)
          vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q
          vi = (wr(i) - p) * 2.0_dp * q
          if (vr .eq. 0.0_dp .and. vi .eq. 0.0_dp) vr = machep * norm  * (abs(w) + abs(q) + abs(x) + abs(y) + abs(zz))
          !           z3 = cmplx(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra) / cmplx(vr,vi)
          !           h(i,na) = real(z3,kind=dp)
          !           h(i,en) = aimag(z3)
          call etdiv(z3r,z3i,x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi)
          h(i,na) = z3r
          h(i,en) = z3i
          if (abs(x) .le. abs(zz) + abs(q)) go to 785
          h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
          h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
          go to 790
          ! 785       z3 = cmplx(-r-y*h(i,na),-s-y*h(i,en)) / cmplx(zz,q)
          !           h(i+1,na) = real(z3,kind=dp)
          !           h(i+1,en) = aimag(z3)
785       call etdiv(z3r,z3i,-r-y*h(i,na),-s-y*h(i,en),zz,q)
          h(i+1,na) = z3r
          h(i+1,en) = z3i
790       continue
       enddo
       !     ********** end complex vector **********
800    continue
    enddo
    !     ********** end back substitution.
    !                vectors of isolated roots **********
    do i = 1, n
       if (i .ge. low .and. i .le. igh) go to 840
       !
       do j = i, n
          z(i,j) = h(i,j)
       enddo
       !
840    continue
    enddo
    !     ********** multiply by transformation matrix to give
    !                vectors of original full matrix.
    !                for j=n step -1 until low do -- **********
    do jj = low, n
       j = n + low - jj
       m = min0(j,igh)
       !
       do i = low, igh
          zz = 0.0_dp
          !
          do k = low, m
             zz = zz + z(i,k) * h(k,j)
          enddo
          !
          z(i,j) = zz

       enddo
    enddo
    !
    go to 1001
    !     ********** set error -- no convergence to an
    !                eigenvalue after 200 iterations **********
1000 ierr = en
1001 return
    !     ********** last card of ety2 **********
  end subroutine ety2

  subroutine etdiv(a,b,c,d,e,f)
!*
    implicit none
    !   computes the complex division
    !     a + ib = (c + id)/(e + if)
    !  very slow, but tries to be as accurate as
    !  possible by changing the order of the
    !  operations, so to avoid under(over)flow
    !  problems.
    !  Written by F. Neri Feb. 12 1986
    !
    integer flip
    real(dp) a,b,c,d,e,f,s,t,cc,dd,ee,ff,temp
    if(.not.c_stable_da) return

    flip = 0
    cc = c
    dd = d
    ee = e
    ff = f
    if( abs(f).ge.abs(e) ) then
       ee = f
       ff = e
       cc = d
       dd = c
       flip = 1
    endif
    s = 1.0_dp/ee
    t = 1.0_dp/(ee+ ff*(ff*s))
    if ( abs(ff) .ge. abs(s) ) then
       temp = ff
       ff = s
       s = temp
    endif
    if( abs(dd) .ge. abs(s) ) then
       a = t*(cc + s*(dd*ff))
    else if ( abs(dd) .ge. abs(ff) ) then
       a = t*(cc + dd*(s*ff))
    else
       a = t*(cc + ff*(s*dd))
    endif
    if ( abs(cc) .ge. abs(s)) then
       b = t*(dd - s*(cc*ff))
    else if ( abs(cc) .ge. abs(ff)) then
       b = t*(dd - cc*(s*ff))
    else
       b = t*(dd - ff*(s*cc))
    endif
    if (flip.ne.0 ) then
       b = -b
    endif
    return
  end subroutine etdiv


subroutine ohmi_factor(a_t,z,r,ok,mf)
implicit none
type(c_damap), intent(inout) :: a_t, z,r
type(c_damap) at,h
type(damap) ma
integer,optional :: mf
integer mf0,i
logical ok
real(dp) norm

mf0=0
if(present(mf)) mf0=mf

call alloc(ma)
call alloc(at,h)

at=a_t
z=1
r=1

 ma=at
call checksymp(ma,norm)

write(6,*) " norm 1",norm

call get_6d_disp(at,h)
call get_6d_ohmi(at,h,z,mf,ok)

at=z**(-1)*a_t
 ma=at
call checksymp(ma,norm)

write(6,*) " norm 2",norm

write(6,*) " teng also",norm
read(5,*) i
if(i==1.and.ok) call get_4d_disp0(at,r,ok)


call kill(at,h)
call kill(ma)
end subroutine ohmi_factor

subroutine get_4d_disp0(a_t,r,ok)  !,h1,sig)
implicit none 
type(c_damap), intent(inout) ::a_t,r
type(c_taylor) h1(2,2),sig(2,2),mu,tc
type(c_taylor) a12(2,2),a22(2,2),a22d(2,2)
integer i,j,nd2n
integer, allocatable :: je(:)
logical ok

call alloc_nn(h1)
call alloc_nn(sig)
call alloc_nn(a12)
call alloc_nn(a22)
call alloc_nn(a22d)
call alloc(mu,tc)


nd2n=6
allocate(je(nd2n))

je=0
do i=1,2
do j=1,2
 je(2+j)=1
 a12(i,j)=a_t%v(i).par.je
 je(2+j)=0
enddo
enddo
do i=1,2
do j=1,2
 je(2+j)=1
 a22(i,j)=a_t%v(i+2).par.je
 je(2+j)=0
enddo
enddo

call   dagger_22(a22,a22d)
call matmul_nn(a12,a22d,h1)
call matmul_nn(a22,a22d,sig)
!call print(a_t%v(1),6)
!call print(a_t%v(2),6)
!call print(a12,6)
! pause 756
!call print(a_t%v(3),6)
!call print(a_t%v(4),6)
!call print(a22,6)
!pause 980
!call print(h1,6)
!pause 981
!call print(sig,6)
if(real(sig(1,1).sub.'0')<=0.0_dp) then
 ok=.false.
 deallocate(je)
 call kill_nn(h1)
 call kill_nn(sig)
 call kill_nn(a12)
 call kill_nn(a22)
 call kill_nn(a22d)
 call kill(mu,tc)
endif

mu=sqrt(sig(1,1))
tc=1.0_dp/mu
call  matmulr_nn(h1,h1,tc)   !mu D1
call   dagger_22(h1,sig)
tc=-1.0_dp
call  matmulr_nn(sig,sig,tc) 
!
r=0
! 1 1 block
do i=1,4
 r%v(i)=mu*(1.0_dp.cmono.i)
enddo
! 3 3
r%v(5)=1.0_dp.cmono.5
r%v(6)=1.0_dp.cmono.6

! 1 2
do i=1,2
do j=1,2
 r%v(i)=r%v(i) + h1(i,j)*(1.0_dp.cmono.(j+2))
enddo
enddo
! 2 1

do i=1,2
do j=1,2
 r%v(i+2)=r%v(i+2) + sig(i,j)*(1.0_dp.cmono.(j))
enddo
enddo




deallocate(je)
call kill_nn(h1)
call kill_nn(sig)
call kill_nn(a12)
call kill_nn(a22)
call kill_nn(a22d)
call kill(mu,tc)

 end subroutine get_4d_disp0 


subroutine get_6d_disp(a_t,h)
implicit none 
type(c_damap), intent(inout) :: a_t, h
type(c_taylor) disp(6),disp_ave0(6)
integer kp,i,n,j,k,k0
integer, allocatable :: je(:)
type(c_damap) r0
complex(dp) w

h=0
call alloc(r0)
call alloc(disp)
call alloc(disp_ave0)

allocate(je(nv))
je=0
do i=1,6
 disp(i)=a_t%v(i)*c_phasor()
enddo


do kp=1,6
       j=1

        do while(.true.) 
          call  c_cycle(disp(kp),j,w ,je); if(j==0) exit;
          k0=0;
          do n=1,2

!!!   keep only the terms in the third plane

         k0=k0+abs(je(2*n))+abs(je(2*n-1))

          enddo
         if(k0==0) disp_ave0(kp)=disp_ave0(kp)+ (w.cmono.je)

       enddo

enddo


!!!! These are all the dispersion function a la Ripken in my Nishikawa paper
do i=1,6
 h%v(i)=disp_ave0(i)
enddo


!!!!  Here I set the initial transverse conditions to be zero
h=h*ci_phasor()*a_t**(-1)
 ! " Dispersion and zeta in terms of initial delta "
deallocate(je)





call kill(r0)
call kill(disp)
call kill(disp_ave0)

 end subroutine get_6d_disp 


subroutine get_6d_ohmi(a_t,h,z,mf,ok)
implicit none 
type(c_damap), intent(inout) :: h,z,a_t
type(c_taylor) h1(2,2),h2(2,2),sig(2,2),sigma,rho,sigmai,det1,det2,lam,tc
type(c_taylor) h1d(2,2),h2d(2,2),t(2,2)
type(c_vector_field) vf,vfs
integer, allocatable :: je(:)
integer i,j,mf,kll
logical ok
type(c_damap) a_cs,q,at
real(dp) norm

ok=.true.
 
vf%n=0;vfs%n=0;

call  alloc_nn(h1)
call  alloc_nn(h2)
call  alloc_nn(h1d)
call  alloc_nn(h2d)
call  alloc_nn(sig)
call  alloc_nn(t)
 call  alloc(det2,lam,tc)
 call  alloc(det1)
call alloc(vf);call alloc(vfs);
call alloc(sigma,rho,sigmai)   
call alloc(a_cs,at,q); 
 
allocate(je(6))
je=0
do i=1,2
do j=1,2
 je(4+j)=1
 h1(i,j)=h%v(i).par.je
 je(4+j)=0
enddo
enddo
do i=1,2
do j=1,2
 je(4+j)=1
 h2(i,j)=h%v(i+2).par.je
 je(4+j)=0
enddo
enddo
do i=1,2
do j=1,2
 je(4+j)=1
 sig(i,j)=h%v(i+4).par.je
 je(4+j)=0
enddo
enddo
 



sigma=sig(1,1)
if(real(sigma.sub.'0')<=0.0_dp) then
 ok=.false.
call  kill_nn(h1)
call  kill_nn(h2)
call  kill_nn(h1d)
call  kill_nn(h2d)
call  kill_nn(sig)
call  kill_nn(t)
 call  kill(det2,lam,tc)
 call  kill(det1)
call kill(vf);call kill(vfs);
call kill(a_cs,at,q); 
call kill(sigma,rho,sigmai)  


 return
endif
rho=sqrt(sigma)
sigmai=1.0_dp/sigma

call  matmulr_nn(h1,h1,sigmai)
call  matmulr_nn(h2,h2,sigmai)
call  dagger_22(h1,h1d)
call  dagger_22(h2,h2d)

det1=h1(1,1)*h1(2,2)-h1(1,2)*h1(2,1)
det2=h2(1,1)*h2(2,2)-h2(1,2)*h2(2,1)
lam=rho**2/(1.0_dp+rho)

z=0
! 1 1 block
t(1,1)=1.0_dp;t(2,2)=1.0_dp;t(1,2)=0.0_dp;t(2,1)=0.0_dp
tc=1.0_dp-lam*det1
call  matmulr_nn(t,t,tc)
do i=1,2
do j=1,2
 z%v(i)=z%v(i) + t(i,j)*(1.0_dp.cmono.j)
enddo
enddo
! 2 2
t(1,1)=1.0_dp;t(2,2)=1.0_dp;t(1,2)=0.0_dp;t(2,1)=0.0_dp
tc=1.0_dp-lam*det2
call  matmulr_nn(t,t,tc)
do i=1,2
do j=1,2
 z%v(i+2)=z%v(i+2) + t(i,j)*(1.0_dp.cmono.(j+2))
enddo
enddo
! 3 3
t(1,1)=1.0_dp;t(2,2)=1.0_dp;t(1,2)=0.0_dp;t(2,1)=0.0_dp
tc=rho
call  matmulr_nn(t,t,tc)
do i=1,2
do j=1,2
 z%v(i+4)=z%v(i+4) + t(i,j)*(1.0_dp.cmono.(j+4))
enddo
enddo
! 1 2
call matmul_nn(h1,h2d,t)
tc=-lam
call  matmulr_nn(t,t,tc)

do i=1,2
do j=1,2
 z%v(i)=z%v(i) + t(i,j)*(1.0_dp.cmono.(j+2))
enddo
enddo
! 2 1
call matmul_nn(h2,h1d,t)
tc=-lam
call  matmulr_nn(t,t,tc)

do i=1,2
do j=1,2
 z%v(i+2)=z%v(i+2) + t(i,j)*(1.0_dp.cmono.(j))
enddo
enddo

! 1 3
 
tc=rho
call  matmulr_nn(h1,t,tc)

do i=1,2
do j=1,2
 z%v(i)=z%v(i) + t(i,j)*(1.0_dp.cmono.(j+4))
enddo
enddo

! 2 3
 
tc=rho
call  matmulr_nn(h2,t,tc)

do i=1,2
do j=1,2
 z%v(i+2)=z%v(i+2) + t(i,j)*(1.0_dp.cmono.(j+4))
enddo
enddo

! 3 1
 
tc=-rho
call  matmulr_nn(h1d,t,tc)

do i=1,2
do j=1,2
 z%v(i+4)=z%v(i+4) + t(i,j)*(1.0_dp.cmono.(j))
enddo
enddo

! 3 2
 
tc=-rho
call  matmulr_nn(h2d,t,tc)

do i=1,2
do j=1,2
 z%v(i+4)=z%v(i+4) + t(i,j)*(1.0_dp.cmono.(j+2))
enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
a_cs=z
q=z

at=a_t
 
do kll=1,no+2

! kll=1+kll

at=a_cs**(-1)*at

h=0
 call get_6d_disp(at,h)
h%v(5)=0
h%v(6)=0


vf=0
do i=1,4
vf%v(i)=h%v(i)
enddo

tc=getpb_from_transverse(vf,vfs)


call c_full_norm_damap(h,norm)

if(mf/=0) write(mf,*) "norm in Ohmi ",norm

 



a_cs=exp(vfs)
 
q=q*a_cs
 
  if(mf/=0) then
  write(mf,*) " Dispersion and zeta in terms of initial delta ",kll
 

  call print(h,mf)

endif
!write(6,*) " more "
!!read(5,*) kkk

 

enddo
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

z=q

deallocate(je)
call  kill_nn(h1)
call  kill_nn(h2)
call  kill_nn(h1d)
call  kill_nn(h2d)
call  kill_nn(sig)
call  kill_nn(t)
 call  kill(det2,lam,tc)
 call  kill(det1)
call kill(vf);call kill(vfs);
call kill(a_cs,at,q); 
call kill(sigma,rho,sigmai) 
end subroutine get_6d_ohmi


subroutine teng_edwards_a1(a1,R_TE,CS_TE,COSLIKE,t_e)
!#general: normal
!# This is the famous Teng-Edwards factorisation 
!# a1= R_TE * CS_TE
!# If the flag t_e is false, then the factorisation failed.
!# This is NOT a global representation of a 4x4 symplectic matrix
!# and it fails for large coupling. 
!# The parameter coslike is true if the R_TE(1,1) entry of the matrix
!# is between -1 and 1. If it is "hyperbolic", then coslike is false.
!# I am not a big fan of this factorisation, but it works for small couplings.


    implicit none
    type(c_damap) , intent(inout) :: a1,R_TE,CS_TE

    type(c_damap) g,s1,s1i

 
    integer i,j

    type(c_taylor) m(4,4),unpb
    type(c_taylor) at(2,2),bt(2,2),ct(2,2),dt(2,2),ati(2,2),bti(2,2),alpha,det
    logical(lp) t_e,COSLIKE
    real(dp) alpha0,epsone
    type(c_vector_field) h

    if(nd2t/=4) then
     write(6,*) " The number of harmonic planes in the orbital part must be 4 "
     t_e=my_false
    return
    endif
    
    call alloc_nn(m)
    call alloc_nn(at)
    call alloc_nn(bt)
    call alloc_nn(ct)
    call alloc_nn(dt)
    call alloc_nn(ati)
    call alloc_nn(bti)
    call alloc(alpha,det)


       t_e=my_true
       call  copy_damap_matrix(a1,m)
       call  copy_matrix_matrix(m(1:2,1:2),at)
       call  copy_matrix_matrix(m(1:2,3:4),ct)
       call  copy_matrix_matrix(m(3:4,1:2),dt)
       call  copy_matrix_matrix(m(3:4,3:4),bt)
        
       call invert_22(at,ati)
       call invert_22(bt,bti)
       if(.not.c_%STABLE_DA) then
        t_e=my_false
       endif 


       call matmul_nn(dt,ati,ati,sc=-1.0_dp)
       call matmul_nn(ati,ct,ct)
       call matmul_nn(ct,bti,ct)
       if(.not.c_%STABLE_DA) then
        t_e=my_false
        goto 888
       endif 


       alpha=ct(1,1)
       alpha0=alpha

       if(alpha0<=-1.0_dp) then
        t_e=my_false
        goto 888
       endif
        
       det=sqrt(1.0_dp/(1.0_dp+alpha))




       if(alpha0>=0.0_dp) then
          COSLIKE=my_true
       else
          ! det=sqrt(one/(one-alpha))
          COSLIKE=my_false
       endif

 
        CS_TE=0
       do i=1,2
          do j=1,2
             CS_TE%v(i)=at(i,j)*(1.0_dp.cmono.j)/det+CS_TE%v(i)
             CS_TE%v(i+2)=bt(i,j)*(1.0_dp.cmono.(j+2))/det+CS_TE%v(i+2)
          enddo
       enddo

       !  The rotation matrix is created but it may not have the correct path length
       !dependence
       if(c_%ndpt/=0.and.t_e) then   
          g%n=nv
          call alloc(g)
          call alloc(h)
          call alloc(unpb)
          call alloc(s1)
          call alloc(s1i)
          epsone=-no

          do i=1,nv
             g%v(i)=1.0_dp.cmono.i
          enddo
          g%v(ndpt)=0.0_dp

          do i=1,nd2t
             s1%v(i) = CS_TE%v(i)*g
          enddo


           s1%v(ndptb)=1.0_dp.cmono.ndptb
           s1%v(ndpt)=1.0_dp.cmono.ndpt


           CS_TE%v(ndptb)=1.0_dp.cmono.ndptb
           CS_TE%v(ndpt) =1.0_dp.cmono.ndpt
!!!!!
          s1i=s1**(-1)
          s1i=s1i*CS_TE    ! s1i is completely nonlinear.
          call c_flofacg(s1i,h,epsone)

          alpha0=1.0_dp
          if(mod(ndpt,2)==0) alpha0=-1.0_dp
          do i=1,nd2t
             h%v(i)=alpha0*(h%v(i).d.ndpt)
          enddo
          call c_int_partial(h,unpb,2)
          !  un%pb=un%vector  ! this is the longitudinal part

             s1i%v(ndptb)=s1i%v(ndptb)+unpb



   
          CS_TE=s1*s1i


          call kill(s1i)
          call kill(s1)
          call kill(unpb)
          call kill(g)

          do i=7,nd2
           cs_te%v(i)=1.0_dp.cmono.i
          enddo

       else


         do i=5,nd2
          cs_te%v(i)=1.0_dp.cmono.i
         enddo

       endif
         888 continue
       if(.not.t_e) then       
        c_%STABLE_DA=my_true
        cs_te=0
        R_TE=0
        write(6,*) " Teng-Edwards is crap : Too much coupling! "
       else       
        R_TE=a1*cs_TE**(-1)
       endif  



    call kill_nn(m)
    call kill_nn(at)
    call kill_nn(bt)
    call kill_nn(ct)
    call kill_nn(dt)
    call kill_nn(ati)
    call kill_nn(bti)
    call kill(alpha,det)
end subroutine teng_edwards_a1

  subroutine c_int_partial(v,h,nd0)
    implicit none
    ! IF SCA=-one
    !     \VEC{V}.GRAD   = J GRAD H . GRAD = :H:
    !
    ! IF SCA=one
    !     \VEC{V}.GRAD  = GRAD H . GRAD
    integer i,nd0
    type(c_vector_field) v
    type(c_taylor) h
    type(c_taylor) b4,b3,b2,b1
    type(c_damap) x
 
    if(.not.c_%stable_da) return


    nd_used=nd0
    call alloc(x)
    call alloc(b4,b3,b2,b1)

    x=1

    do i=1,nd_used
       call cfu(v%v(2*i-1),dlie,b3)
       call cfu(v%v(2*i),dlie,b1)
       b2=b1*x%v(2*i-1)
       b1=b3*x%v(2*i)
       b3=b2-b1
       b2=b3+b4
       b4=b2
    enddo
    h=b4


    call kill(b4,b3,b2,b1)
    call kill(x)


  end subroutine c_int_partial

  complex(dp) function dlie(j)
    implicit none
    integer i
    !      INTEGER J(NTT)
    integer,dimension(:)::j
    dlie=0.0_dp  
    if(.not.c_%stable_da) return


    do i=1,nd_used
       dlie=REAL(j(2*i-1)+j(2*i),kind=DP)+dlie
    enddo
    dlie=dlie+1.0_dp
    dlie=1.0_dp/dlie
    return
  end function dlie

  subroutine copy_damap_matrix(mi,a)
    implicit none
    type(c_taylor), intent(inout) :: a(4,4)
    type(c_damap), intent(in) :: mi
    type(c_damap) m

    integer i,j,nt

    integer, allocatable :: jl(:)

    call alloc(m)

    m=mi

    nt=4

    allocate(jl(nt))
    jl=0
    do i=1,nt
       do j=1,nt
          jl(j)=1
          a(i,j)=m%v(i).par.jl
          jl(j)=0
       enddo
    enddo

    call kill(m)
    deallocate(jl)

  end subroutine copy_damap_matrix

 subroutine invert_22(a,ai)
!*
    implicit none
    type(c_taylor) a(2,2),ai(2,2),t(2,2)
    type(c_taylor) det
    call alloc_nn(t)
    call alloc(det)

    det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
    t(1,1)=a(2,2)/det
    t(2,2)=a(1,1)/det
    t(1,2)=-a(1,2)/det
    t(2,1)=-a(2,1)/det

    call copy_matrix_matrix(t,ai)

    call kill_nn(t)
    call kill(det)

  end subroutine invert_22

 subroutine dagger_22(a,ai)
!*
    implicit none
    type(c_taylor) a(2,2),ai(2,2),t(2,2),at(2,2)

    call alloc_nn(t)
    call alloc_nn(at)

    t(1,2)=1.0_dp
    t(2,1)=-1.0_dp

    at(1,1)=a(1,1)
    at(2,2)=a(2,2)
    at(1,2)=a(2,1)
    at(2,1)=a(1,2)

    call matmul_nn(t,at,at)
    call matmul_nn(at,t,at,-1.0_dp)

    call copy_matrix_matrix(at,ai)

    call kill_nn(t)
    call kill_nn(at)

  end subroutine dagger_22

  subroutine matmulr_33(m,mo,sc0)
!*
!  mo=sc m.n
    implicit none
    type(c_taylor) m(:,:),mo(:,:)
    real(dp) sc0
    integer i,j,k

    do i=1,size(m,dim=1)
       do j=1,size(m,dim=2)
             mo(i,j)=m(i,j)*sc0
       enddo
    enddo
  end subroutine matmulr_33

  subroutine matmult_33(m,mo,sc0)
!*
!  mo=sc m.n
    implicit none
    type(c_taylor) m(:,:),mo(:,:)
    type(c_taylor) sc0
    integer i,j,k

    do i=1,size(m,dim=1)
       do j=1,size(m,dim=2)
             mo(i,j)=m(i,j)*sc0
       enddo
    enddo
  end subroutine matmult_33

  subroutine matmul_33(m,n,mo,sc)
!*
!  mo=sc m.n
    implicit none
    type(c_taylor) m(:,:),n(:,:),mo(:,:)
    type(c_taylor), allocatable :: a(:,:)
    real(dp), optional :: sc
    real(dp) sc0
    integer i,j,k
    sc0=1.0_dp
    allocate(a(size(m,dim=1),size(n,dim=2)))

    call alloc_nn(a)

    do i=1,size(m,dim=1)
       do j=1,size(m,dim=2)
          do k=1,size(n,dim=2)
             a(i,k)=m(i,j)*n(j,k)+a(i,k)
          enddo
       enddo
    enddo

    if(present(sc)) sc0=sc
    do i=1,size(mo,dim=1)
       do j=1,size(mo,dim=2)
          mo(i,j)=sc0*a(i,j)
       enddo
    enddo
    call kill_nn(a)
    deallocate(a)
  end subroutine matmul_33

  subroutine alloc_33t(a)
!*
    implicit none
    type(c_taylor) a(:,:)
    integer i,j

    do i=1,size(a,dim=1)
       do j=1,size(a,dim=2)
          call alloc(a(i,j))
       enddo
    enddo

  end subroutine alloc_33t

  subroutine print_33t(a,mf,prec)
!*
    implicit none
    type(c_taylor) a(:,:)
    integer i,j,mf
    real(dp), optional :: prec

    do i=1,size(a,dim=1)
       do j=1,size(a,dim=2)
          write(mf,*) i,j
          call print(a(i,j),mf,prec)
       enddo
    enddo

  end subroutine print_33t


  subroutine kill_33t(a)
!*
    implicit none
    type(c_taylor) a(:,:)
    integer i,j

    do i=1,size(a,dim=1)
       do j=1,size(a,dim=2)
          call kill(a(i,j))
       enddo
    enddo

  end subroutine kill_33t

  subroutine copy_matrix_matrix(ma,a)
!*
    implicit none
    type(c_taylor), intent(inout) :: a(:,:)
    type(C_taylor), intent(in) :: ma(:,:)
    integer i,j

    do i=1,size(a,dim=1)
       do j=1,size(a,dim=2)
          a(i,j)=ma(i,j)
       enddo
    enddo

  end subroutine copy_matrix_matrix




subroutine extract_linear_from_normalised(m,a1,phi1,f1,f2,integer_part)
    implicit none
    type(c_damap) , intent(inout) :: m,a1,phi1
    type(c_vector_field) , intent(inout) :: f1,f2
    real(dp), intent(IN)::    integer_part(:)
    type(c_damap) b1,a
    complex(dp) v
 
    integer i,j,k,kr
    integer, allocatable :: je(:)

        

    call alloc(a)
    call alloc(b1) 

     a=to_phasor()*m*from_phasor() 

    allocate(je(nv))
     je=0
 
    
 
!! extract the linear part as a function parameters (delta included)
     
      do i=1,nd2
       j=1
        do while(.true.) 


          call  c_cycle(a%v(i),j,v ,je); if(j==0) exit;   
          kr=0
          do k=1,nd2
           if(k==ndpt) cycle
           kr=je(k)+kr
          enddo
          if(i==ndptb) then
           if(kr<=2) then
            b1%v(i)=b1%v(i)+(v.cmono.je)
           endif
          else
           if(kr==1) then
            b1%v(i)=b1%v(i)+(v.cmono.je)
           endif
          endif
       enddo
     enddo


     if(ndpt/=0)  then 
      b1%v(ndpt)=1.0_dp.cmono.ndpt
      f1%v(ndptb)=b1%v(ndptb)-(1.0_dp.cmono.ndptb)
     endif

     a1=a*b1**(-1)

     phi1=b1

     do i=1,nd2t
      f1%v(i)=(log(b1%v(i).d.i).cut.c_%no) *(1.0_dp.cmono.i)
     enddo
     do i=1,nd2t/2
      f1%v(2*i-1)=-(i_*twopi*integer_part(i).cmono.(2*i-1))+ f1%v(2*i-1)
      f1%v(2*i)=(i_*twopi*integer_part(i).cmono.(2*i))+ f1%v(2*i)
     enddo

      f2=log(a1)

     f2=transform_vector_field_by_map(f2,to_phasor())     
     f1=transform_vector_field_by_map(f1,to_phasor())


    deallocate(je)
     call kill(b1) 


end subroutine extract_linear_from_normalised

subroutine extract_a0(a,a0)
!#internal: manipulation
!# This routines extracts a0: the full fixed point map.
    implicit none
    type(c_damap) , intent(inout) :: a,a0
 
    type(c_damap) a0t,at,ai,x,v
    type(c_taylor)  t1
    integer i
    integer np_pos


    ! a0 is for the fixed point
    
    at%n=nv
    call alloc(at)
    call alloc(ai)
    call alloc(a0t)
    call alloc(t1)
    call alloc(x)
    call alloc(v)    
     ai=a
    
!! extract the dispersion
     if(ndpt/=0) at%v(ndpt)=1.0_dp.cmono.ndpt
     np_pos=nv-np+1
     do i=np_pos,nv
      at%v(i)=1.0_dp.cmono.i      
     enddo     
     
     a0t=ai*at

 

!!! Force the symplectic condition in the time variable
    if(ndpt/=0) then 
         t1=0
         do i=1,ndt
            x%v(2*i)  =(-1)**(2*i-1)*(a0t%v(2*i-1))
            x%v(2*i-1)=(-1)**(2*i  )*(a0t%v(2*i))
            v%v(2*i)=   x%v(2*i).d.ndpt
            v%v(2*i-1)= x%v(2*i-1).d.ndpt
         enddo
         do i=1,ndt
             t1=-(1.0_dp.cmono.(2*i-1))*v%v(2*i-1)+t1  ! first order
             t1=-(1.0_dp.cmono.(2*i))*v%v(2*i)+t1      ! first order
             t1=-0.5_dp*(x%v(2*i-1)*v%v(2*i)-x%v(2*i)*v%v(2*i-1))+t1  ! second order
         enddo
         t1=(-1)**ndpt*t1
         a0t%v(ndptb)=(1.0_dp.cmono.ndptb)+t1 !!! effect on  time added to identity map in the time-energy plane
         a0t%v(ndpt)=1.0_dp.cmono.ndpt
  
    endif


     do i=1,nd2
      if(i/=ndpt.and.i/=ndptb) a0t%v(i)=a0t%v(i)+(1.0_dp.cmono.i)      
     enddo  
    if(ndpt/=0) then
     if(.not.time_lie_choice) then
         t1=0.0_dp
       do i=1,nd2-ndc2t
         t1=t1+(1.0_dp.cmono.i)*(a0t%v(ndptb).d.i)
       enddo
        a0t%v(ndptb)=(1.0_dp.cmono.ndptb)+t1 
     endif
    endif

 
    a0=a0t
    a=a0**(-1)*a

 


    call kill(at)
    call kill(ai)
    call kill(t1) 
    call kill(a0t)  
    call kill(x)
    call kill(v)  
end subroutine extract_a0

subroutine extract_only_a0(a,a0)
!#internal: manipulation
!# This routines extracts a0: the full fixed point map.
    implicit none
    type(c_damap) , intent(inout) :: a,a0
 
    type(c_damap) a0t,at,ai,x,v
    type(c_taylor)  t1
    integer i
    integer np_pos


    ! a0 is for the fixed point
    
    at%n=nv
    call alloc(at)
    call alloc(ai)
    call alloc(a0t)
    call alloc(t1)
    call alloc(x)
    call alloc(v)    
     ai=a
    
!! extract the dispersion
     if(ndpt/=0) at%v(ndpt)=1.0_dp.cmono.ndpt
     np_pos=nv-np+1
     do i=np_pos,nv
      at%v(i)=1.0_dp.cmono.i      
     enddo     
     
     a0t=ai*at

 

!!! Force the symplectic condition in the time variable
    if(ndpt/=0) then 
         t1=0
         do i=1,ndt
            x%v(2*i)  =(-1)**(2*i-1)*(a0t%v(2*i-1))
            x%v(2*i-1)=(-1)**(2*i  )*(a0t%v(2*i))
            v%v(2*i)=   x%v(2*i).d.ndpt
            v%v(2*i-1)= x%v(2*i-1).d.ndpt
         enddo
         do i=1,ndt
             t1=-(1.0_dp.cmono.(2*i-1))*v%v(2*i-1)+t1  ! first order
             t1=-(1.0_dp.cmono.(2*i))*v%v(2*i)+t1      ! first order
             t1=-0.5_dp*(x%v(2*i-1)*v%v(2*i)-x%v(2*i)*v%v(2*i-1))+t1  ! second order
         enddo
         t1=(-1)**ndpt*t1
         a0t%v(ndptb)=(1.0_dp.cmono.ndptb)+t1 !!! effect on  time added to identity map in the time-energy plane
         a0t%v(ndpt)=1.0_dp.cmono.ndpt
  
    endif


     do i=1,nd2
      if(i/=ndpt.and.i/=ndptb) a0t%v(i)=a0t%v(i)+(1.0_dp.cmono.i)      
     enddo  


 
    a0=a0t
    a=a0**(-1)*a

 !   if(ndpt/=0) then
 !        call print(a%v(ndpt),6)
 !        at=0
 !        at%v(ndpt)=(1.0_dp.cmono.ndpt)
 !        ai=a*at
 !        at=1
 !        at%v(ndpt)=ai%v(ndpt)
 !        call print(at%v(ndpt),6)
 !        a=at**(-1)*a
 !        a0=a0*at
 !   endif
 


    call kill(at)
    call kill(ai)
    call kill(t1) 
    call kill(a0t)  
    call kill(x)
    call kill(v)  
end subroutine extract_only_a0

 subroutine extract_a1(a,a1,phi1)
!#internal: manipulation
!# This routines extracts a1: the full linear canonical transformation.
!# Teng-Edward A_12=0 is imposed  or, for fun, Anti-Teng-Edwards A_21=0
!# based on the flag courant_snyder_teng_edwards (defaulted to true).
!# See Sec.7.4 of my Springer book. Phi1 is the phase advanced induced
!# by the choice of canonical transformation. Sec 2.4 gives a review in 1-d-f.
    implicit none
    type(c_damap) , intent(inout) :: a,a1
    type(c_damap),optional , intent(inout) :: phi1
    type(c_damap) b1 
    type(c_taylor)  cphi,sphi,t
    complex(dp) v
 
    integer i,j,k,kr
    integer, allocatable :: je(:)

 
     call alloc(b1) 
    call alloc(cphi,sphi,t)
 
    allocate(je(nv))
     je=0
 
    
 
!! extract the linear part as a function parameters (delta included)
     
      do i=1,nd2
       j=1
        do while(.true.) 


          call  c_cycle(a%v(i),j,v ,je); if(j==0) exit;   
          kr=0
          do k=1,nd2
           if(k==ndpt) cycle
           kr=je(k)+kr
          enddo
          if(kr==1) then
           b1%v(i)=b1%v(i)+(v.cmono.je)
          endif
       enddo
     enddo
     if(ndpt/=0)  then 
       b1%v(ndpt)=b1%v(ndpt)+(1.0_dp.cmono.ndpt)
 
       j=1
        do while(.true.) 


          call  c_cycle(a%v(ndptb),j,v ,je); if(j==0) exit;   
          kr=0
          do k=1,nd2t
           kr=je(k)+kr
          enddo
          if(kr==2) then
           b1%v(ndptb)=b1%v(ndptb)+(v.cmono.je)
          endif
       enddo
            
     endif
    
 
     a=b1**(-1)*a

     a1=b1

!  imposed Teng-Edward A_12=0 or, for fun, Anti-Teng-Edwards A_21=0

    if(present(phi1)) then

     phi1=1
      do i=1,nd2/2
      if((i<=ndt).or.(i>nd-rf)) then
       if(courant_snyder_teng_edwards) then
        t=sqrt((b1%v(2*i-1).d.(2*i))**2+(b1%v(2*i-1).d.(2*i-1))**2)
        cphi=(b1%v(2*i-1).d.(2*i-1))/t
        sphi=(b1%v(2*i-1).d.(2*i))/t
       else
        t=sqrt((b1%v(2*i).d.(2*i-1))**2+(b1%v(2*i).d.(2*i))**2)
        cphi=(b1%v(2*i).d.(2*i))/t
        sphi=-(b1%v(2*i).d.(2*i-1))/t
       endif
        phi1%v(2*i-1) =cphi*(1.0_dp.cmono.(2*i-1))+sphi*(1.0_dp.cmono.(2*i))
        phi1%v(2*i)   =cphi*(1.0_dp.cmono.(2*i))-sphi*(1.0_dp.cmono.(2*i-1))

!!! Corrects the time variable with the d/dp_t term
       if(ndpt/=0) then         
        t= cphi + i_*sphi
         t=-i_*log(t)
         t=-(t.d.ndpt)/(2.0_dp,0.0_dp)
         phi1%v(ndptb)=phi1%v(ndptb)+(-1)**(ndptb)*t*((1.0_dp.cmono.(2*i-1))**2+(1.0_dp.cmono.(2*i))**2)
       endif

       endif
      enddo

    

            a1=a1*phi1**(-1)
            a=phi1*a*phi1**(-1) 

 
  endif
    
 


    deallocate(je)
     call kill(b1) 
    call kill(cphi,sphi,t)
end subroutine extract_a1

 subroutine extract_only_a1(a,a1)
!#internal: manipulation
!# This routines extracts a1

    implicit none
    type(c_damap) , intent(inout) :: a,a1
    type(c_damap) b1 
    type(c_taylor)  t
    complex(dp) v
 
    integer i,j,k,kr
    integer, allocatable :: je(:)

 
     call alloc(b1) 
    call alloc(t)
 
    allocate(je(nv))
     je=0
 
    
 
!! extract the linear part as a function parameters (delta included)
     
      do i=1,nd2
       j=1
        do while(.true.) 


          call  c_cycle(a%v(i),j,v ,je); if(j==0) exit;   
          kr=0
          do k=1,nd2

           kr=je(k)+kr
          enddo
          if(kr==1) then
           b1%v(i)=b1%v(i)+(v.cmono.je)
          endif
       enddo
     enddo

    
 
     a=b1**(-1)*a

     a1=b1

!  imposed Teng-Edward A_12=0 or, for fun, Anti-Teng-Edwards A_21=0

 


    deallocate(je)
     call kill(b1) 
    call kill(t)
end subroutine extract_only_a1



subroutine extract_a2(a,phi2)
    implicit none
    type(c_damap) , intent(inout) :: a
    type(c_damap),optional , intent(inout) :: phi2
    type(c_damap) b1 
    type(c_vector_field) h,hr,hf

    real(dp) eps
    complex(dp) v
 
    integer j,k,kr
    integer, allocatable :: je(:)
    logical(lp) removeit
 
     if(.not.present(phi2)) return
 
     call alloc(b1) 
     call alloc(h) 
     call alloc(hr) 
     call alloc(hf) 



    allocate(je(nv))
     je=0


 
    b1=c_simil(from_phasor(-1),a,1)

 
    eps=-no-1  !1.0e-7_dp

 !! Here I assume that in a single Lie exponent representation there are NO kernel terms
 !do kr=2,no
  do kr=1,no

 hr=0
    call c_flofacg(b1,h,eps)


  do k=1,nd2
  j=1
        do while(.true.) 

          call  c_cycle(h%v(k),j,v ,je); if(j==0) exit;
          call check_kernel(k,nd2,je,removeit)
           if(.not.removeit) then
             hr%v(k)=hr%v(k)+(v.cmono.je)
           endif

       enddo
    hf%v(k)=hf%v(k)+hr%v(k)
 enddo
   b1=texp(-hr,b1) 

 enddo  

   a=c_simil(from_phasor(-1),b1,-1)
   phi2=texp(hf) 
   phi2=c_simil(from_phasor(-1),phi2,-1)
   


    deallocate(je)
  


     call kill(b1) 
     call kill(h) 
     call kill(hr) 
     call kill(hf) 

end subroutine extract_a2


 ! spin routine
 subroutine factor_ely_rest(as_xyz,as_y,r_y,n_y,n_r) !r_y in normal units, tune in phasors   as_xyz=  r_y o as_xyz 
    implicit none
    type(c_damap), intent(inout) :: as_xyz,as_y
    type(c_damap) , intent(inout) :: r_y
    type(c_spinor), intent(inout) ::  n_y
    type(c_spinor), intent(inout) ::  n_r
    type(c_damap) temp,as_nl,rot_y
    type(c_spinor) n_expo,n_tune,tune0
    type(c_taylor) tr
    type(c_taylor) t
    integer i
    integer  nmax
    real(dp) eps,norm1,norm2,d
    logical check
!!!  original as_xyz = as_xyz*r_y = a_y*a_nl*r_y  on exit
    check=.true.
    eps=1.d-6
    nmax=1000
 
    call alloc(n_expo)
    call alloc(n_tune)
    call alloc(tune0)
    call alloc(temp,as_nl,rot_y)
    call alloc(tr)
    call alloc(t)
    tune0=0
    as_y=as_xyz
       as_y= as_y 
    norm1=mybig
    rot_y=1
    temp=1
    check=.true.
    norm1=mybig

    !   call c_full_norm_spin(as_xyz%s,i,EPS=d)
      call c_full_norm_spin_map(as_xyz,i,d)
 
      eps=d*eps
!      write(6,*) eps


!pause 1256

    do i=1,nmax
   !    call find_exponent_only(as_y,n_expo)

 

       n_expo=log(as_y%s,exact=my_false)




       !     call dalog_spinor_8(as_y,n_expo)
       norm2=0.0_dp
  !     call c_n0_to_nr(n_expo,n_tune)  ! not necessary
 
       n_tune=n_expo
 
       n_tune%v(1)=0.0_dp
       n_tune%v(3)=0.0_dp
 
 
    !   call cfu(n_tune%v(2),c_phase_shift,n_tune%v(2))
  
 
        tune0%v(2)=tune0%v(2)+n_tune%v(2) 
        norm2=FULL_ABS(n_tune%v(2)) 
 
       temp%s=exp(n_tune)
 
 
       rot_y=rot_y*temp
    !   call inv_as(temp%s%s,temp%s%s)
       as_y=temp**(-1)*as_y
   !    write(6,*) norm1, norm2,i
       if(check) then
          if(norm2<eps.and.i>10) then
             check=.false.
          endif
       else
          if(norm2>=norm1) exit   !.and.i>40) exit
       endif
       norm1=norm2
    enddo

    if(i>nmax-10) then
       write(6,*) "no convergence in remove_y_rot "
       !stop 1067
    endif
    
  !  as_xyz= as_y
 
      r_y= rot_y
 
       n_r=log(as_y%s,exact=my_false)
 
       n_y=0
  
       n_y%v(2)= tune0%v(2)   ! in phasors
 

    call kill(n_expo)
    call kill(temp,as_nl,rot_y)
    call kill(tune0)
    call kill(tr)
    call kill(n_tune)
    call kill(t)

  end subroutine factor_ely_rest

  ! spin routine
 subroutine c_remove_y_rot(as_xyz,r_y,tune) !r_y in normal units, tune in phasors   as_xyz= as_xyz o r_y
    implicit none
    type(c_damap), intent(inout) :: as_xyz
    type(c_damap), optional , intent(inout) ::r_y
    type(c_taylor), optional ,intent(inout) ::  tune
    type(c_damap) temp,as_y,as_nl,rot_y
    type(c_spinor) n_expo,n_tune,tune0
    type(c_taylor) tr
    type(taylor) si,co
     real(dp) si0,co0
    type(c_taylor) t
    type(c_quaternion) qnr
    type(q_linear) q,qr
    integer i
    integer  nmax
    real(dp) eps,norm1,norm2,d,dt,aq
    logical check
!!!  original as_xyz = as_xyz*r_y = a_y*a_nl*r_y  on exit
    check=.true.
    eps=1.d-9
    nmax=1000
 
    call alloc(n_expo)
    call alloc(n_tune)
    call alloc(tune0)
    call alloc(temp,as_y,as_nl,rot_y)
    call alloc(tr)
    call alloc(si,co)
    call alloc(t)
    call alloc(qnr)
    tune0=0
    as_y=as_xyz
       as_y=from_phasor(-1)*as_y*from_phasor()
    norm1=mybig
    rot_y=1
    temp=1
    check=.true.
    norm1=mybig
dt=0
!       call c_full_norm_spin(as_xyz%s,i,EPS=d)
      call c_full_norm_spin_map(as_xyz,i,d)

      eps=d*eps
!      write(6,*) eps

!pause 1256

    do i=1,nmax
   !    call find_exponent_only(as_y,n_expo)

 
      if(use_quaternion) then
 
      
      if(i==1.and.qphase) then
          q=1
          q=as_y%q
 !  make sure isf not below y plane
       aq=q%q(0,0)**2-(q%q(1,0)**2+q%q(2,0)**2+q%q(3,0)**2)
if(aq<0) then
         temp%q=1   !  = i 
         as_y%q=as_y%q*temp%q
         q=as_y%q
endif
 
            aq=-atan2(real(q%q(2,0)),real(q%q(0,0)))
            temp%q=1.0_dp
            temp%q%x(0)= cos(aq)
            temp%q%x(2)= -sin(aq)

            n_tune%v(1)=0.0_dp
            n_tune%v(3)=0.0_dp
            n_tune%v(2)=-aq*2.0_dp
 
     else

            temp%q=1.0_dp
            si0=as_y%q%x(2)
            co0=as_y%q%x(0)
!call print(as_y%q,16)
            n_expo=as_y%q 
             call c_n0_to_nr(n_expo,n_tune)  ! not necessary
  
            n_tune=n_expo
  
            n_tune%v(1)=0.0_dp
            n_tune%v(3)=0.0_dp
  
  
            call cfu(n_tune%v(2),c_phase_shift,n_tune%v(2))
             qnr=n_tune
           temp%q=exp(qnr)
            n_tune%v(2)=n_tune%v(2)*2.0_dp
  endif
 
       else
            n_expo=log(as_y%s,exact=my_false)
            !     call dalog_spinor_8(as_y,n_expo)
            norm2=0.0_dp
             call c_n0_to_nr(n_expo,n_tune)  ! not necessary
  
            n_tune=n_expo
  
            n_tune%v(1)=0.0_dp
            n_tune%v(3)=0.0_dp
  
  
            call cfu(n_tune%v(2),c_phase_shift,n_tune%v(2))
              temp%s=exp(n_tune)
      endif


            dt=dt+n_tune%v(2)
            tune0%v(2)=tune0%v(2)+n_tune%v(2) 


        norm2=FULL_ABS(n_tune%v(2)) 
 
  !      call c_nr_to_n0(n_tune,n_expo) 
       !  enddo




 
 
       rot_y=temp*rot_y
    !   call inv_as(temp%s%s,temp%s%s)
       as_y=as_y*temp**(-1)
   !    write(6,*) norm1, norm2,i
if(use_quaternion) then
       if(check) then
          if(norm2<eps) then
             check=.false.
          endif
       else
          if(norm2>=norm1) goto 123   !.and.i>40) exit
       endif
else
       if(check) then
          if(norm2<eps) then
             check=.false.
          endif
       else

          if(norm2>=norm1) exit   !.and.i>40) exit
       endif
endif
       norm1=norm2

    enddo

    if(i>nmax-10) then
       write(6,*) "no convergence in remove_y_rot ",norm2,norm1,eps,i,norm2>=norm1
       !stop 1067
    endif
    
123    as_xyz=from_phasor()*as_y*from_phasor(-1)

 
    if(present(r_y)) then
      r_y=from_phasor()*rot_y*from_phasor(-1) 
    endif
       if(present(tune)) then
        tune=0.0_dp
        tune=spin_def_tune*tune0%v(2)   ! in phasors
       endif

 

    call kill(n_expo)
    call kill(temp,as_y,as_nl,rot_y)
    call kill(tune0)
    call kill(si,co)
    call kill(tr)
    call kill(n_tune)
    call kill(t)
    call kill(qnr)


  end subroutine c_remove_y_rot

  complex(dp) function c_phase_shift(j)
    implicit none

    !      INTEGER J(NTT)
    integer,dimension(:)::j

    logical(lp) removeit

    c_phase_shift=0.d0
    if(.not.c_%stable_da) return


    call check_kernel(0,nd2,j,removeit)

    if(.not.removeit) c_phase_shift=1.0_dp

    return
  end function c_phase_shift
  !
  complex(dp) function c_zeroth(j)
    implicit none
    integer i,k
    !       keeps in regular harmonic planes zeroth order
    integer,dimension(:)::j

 
    c_zeroth=0.d0
    if(.not.c_%stable_da) return

    k=0
    do i=1,nd2t
    k=k+abs(j(i))
    enddo

    if(k==0) c_zeroth=1.0_dp

    return
  end function c_zeroth
!!!!!!!!!!!!!!!!!!!!! tracking unitary !!!!!!!!!!!!!!!!!!!!!!!!

 ! spin routine
subroutine produce_orthogonal(n,ray,s)
implicit none
type(c_spinor),intent(inout) :: n
type(c_spinmatrix),intent(inout) ::  s
type(c_ray),intent(inout):: ray
type(c_spinor) nt

call alloc(nt)

nt=n.o.ray

 s=exp(nt)

call kill(nt)
end subroutine produce_orthogonal 




 ! spin routine


 ! spin routine
 subroutine orthogonalise_ray(ray)
 implicit none
 type(c_ray),intent(inout):: ray
 real(dp)n
 integer i
   
      n=0.0_dp
    do i=1,3
      n=n+abs(ray%s1(i))**2    
    enddo
      ray%s1=ray%s1/sqrt(n)


   
    n=0.0_dp

    DO I=1,3
       n=n+ray%s1(i)*ray%s2(i)
    ENDDO

    DO I=1,3
    ray%s2(i)=ray%s2(i)-n*ray%s1(i)
    ENDDO
     
    n=0.0_dp
    do i=1,3
      n=n+abs(ray%s2(i))**2    
    enddo

      ray%s2=ray%s2/sqrt(n)

        ray%s3(1)=ray%s1(2)*ray%s2(3)-ray%s1(3)*ray%s2(2)
        ray%s3(2)=ray%s1(3)*ray%s2(1)-ray%s1(1)*ray%s2(3)
        ray%s3(3)=ray%s1(1)*ray%s2(2)-ray%s1(2)*ray%s2(1)

      n=0.0_dp
    do i=1,3
      n=n+abs(ray%s3(i))**2    
    enddo
      ray%s3=ray%s3/sqrt(n)

  end subroutine orthogonalise_ray 

  SUBROUTINE  c_IdentityEQUALVECfourier(S2,S1)
!*
    implicit none
    type (c_vector_field_fourier),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1
    integer i
    IF(.NOT.C_STABLE_DA) RETURN
     if(s1/=0) then 
       write(6,*) "c_IdentityEQUALVECfourier"
      stop
     endif
     do i=-s2%n,s2%n
      s2%f(i)=0
     enddo

  END SUBROUTINE c_IdentityEQUALVECfourier

SUBROUTINE  equal_c_vector_field_fourier(S1,s2)
!*
    implicit none
    type (c_vector_field_fourier),INTENT(INOUT) :: S1
    type (c_vector_field_fourier),INTENT(IN) :: S2
    INTEGER i

    do i=-n_fourier,n_fourier
     s1%f(i)=s2%f(i)
    enddo
    s1%n=s2%n

end SUBROUTINE  equal_c_vector_field_fourier

  SUBROUTINE  alloc_c_vector_field_fourier(S1,N)
!*
    implicit none
    type (c_vector_field_fourier),INTENT(INOUT) :: S1
    integer,optional :: n
    INTEGER i
   ! logical existed

   ! existed=.false.
  !  if(s1%n/=0) existed=.true.

   ! if(present(n)) then
   !  s1%n=n
   ! else
     s1%n=n_fourier
   ! endif
     
     if(associated(s1%f) ) then
      deallocate(s1%f)
      nullify(s1%f)
     endif
     allocate(s1%f(-s1%n:s1%n))

    do i=-s1%n,s1%n
      s1%f(i)%n=0
     if(present(n)) s1%f(i)%n=n
     call alloc(s1%f(i))
    enddo


  END SUBROUTINE alloc_c_vector_field_fourier

  SUBROUTINE  kill_c_vector_field_fourier(S1)
!*
    implicit none
    type (c_vector_field_fourier),INTENT(INOUT) :: S1
    INTEGER i
 

   ! if(present(n)) then
   !  s1%n=n
   ! else
   !  s1%n=0
   ! endif
     



    do i=-s1%n,s1%n
     call kill(s1%f(i))
    enddo
           s1%n=0
     deallocate(s1%f)

  END SUBROUTINE kill_c_vector_field_fourier

  subroutine transform_vector_field_fourier_by_map(s1,s2,m)
    implicit none
    TYPE (c_vector_field_fourier), INTENT (INout) :: S1 ,s2
    TYPE (c_damap), INTENT (IN) :: m
 
    integer i


     IF(.NOT.C_STABLE_DA) then
     RETURN
     endif
!!! The tunes are stored for the nonlinear normal form recursive algorithm


     do i=-n_fourier,n_fourier
      s2%f(i)=transform_vector_field_by_map(s1%f(i),m)
     enddo
   
end  subroutine transform_vector_field_fourier_by_map

! etienne

subroutine exp_vector_field_fourier(F,H,K,nlin)
    implicit none
    TYPE (c_vector_field_fourier), INTENT (INout) :: F,H,K
    TYPE (c_vector_field_fourier) s3,t,dF_dt,dhs
    integer, optional :: nlin
    INTEGER I ,nl
    real(dp) norm(3)
    complex(dp) fac,fac1
    nl=0
    if(present(nlin)) nl=nlin
    call alloc(s3) 
    call alloc(t)
    call alloc(dhs)
    call alloc(dF_dt)
    
    call ddt_vector_field_fourier(F,dF_dt) ! (1)

    s3=H ; t=H ; dhs=dF_dt;
    fac=1.0_dp;    fac1=1.0_dp;

    if(nl/=0) write(6,*) "Priting Iterations for convergence check "
    do i=1,no+nl   ! extra terms for cheap convergence
     fac=1.0_dp/i
     
      call bra_vector_field_fourier(F,t,t) ! (2a)
      call mulc_vector_field_fourier(t,fac,t)  !  t=fac*t
      call add_vector_field_fourier(s3,t,s3) ! (2b) 
      call bra_vector_field_fourier(F,dF_dt,dF_dt) ! (3a)
      fac1=1.0_dp/(i+1)
      call mulc_vector_field_fourier(dF_dt,fac1,dF_dt)
      call add_vector_field_fourier(dhs,dF_dt,dhs) ! (3b)
       if(i>no.and.i>no+nl-10) then
         call c_full_norm_fourier(s3,norm(1))
         call c_full_norm_fourier(t,norm(2))
         call c_full_norm_fourier(dF_dt,norm(3))
         write(6,'(i4,1x,3(g23.16,1x))')i, norm
       endif
    enddo
      fac=-1.d0
      call add_vector_field_fourier(s3,dhs,s3,fac) ! (3c)
     K=s3
     call kill(s3)
     call kill(t)
     call kill(dhs)
     call kill(dF_dt)

end subroutine exp_vector_field_fourier


subroutine ddt_vector_field_fourier(s1,ds1)
    implicit none
    TYPE (c_vector_field_fourier), INTENT (INout) :: S1,ds1
    INTEGER I
 
    ds1%f(0)=0
    do i=1,n_fourier
     ds1%f(i)= (i*i_)*s1%f(i)
     ds1%f(-i)= (-i*i_)*s1%f(-i)
    enddo

end subroutine ddt_vector_field_fourier

subroutine print_vector_field_fourier(s1,mf)
    implicit none
    TYPE (c_vector_field_fourier), INTENT (INout) :: S1 
    INTEGER I,mf
 
     write(mf,*) 0,"th mode"
     call print(s1%f(0),mf,dospin=.false.)
    do i=1,n_fourier
     write(mf,*) i,"th mode"
     call print(s1%f(i),mf,dospin=.false.)
     call print(s1%f(-i),mf,dospin=.false.)
    enddo

end subroutine print_vector_field_fourier

subroutine print_poisson_bracket_fourier(s1,mf)
    implicit none
    TYPE (c_vector_field_fourier), INTENT (INout) :: S1 
    type(c_taylor) h
    INTEGER I,mf
    call alloc(h)
     write(mf,*) 0,"th mode"
!     call print(s1%f(0),mf)
      h=cgetpb(s1%f(0))
      call print(h,mf)
    do i=1,n_fourier
     write(mf,*) i,"th mode"
      h=cgetpb(s1%f(i))
      call print(h,mf)
      h=cgetpb(s1%f(-i))
      call print(h,mf)
  !   call print(s1%f(i),mf)
  !   call print(s1%f(-i),mf)
    enddo
    call kill(h)
end subroutine print_poisson_bracket_fourier

subroutine bra_vector_field_fourier(s1,s2,r)
    implicit none
    TYPE (c_vector_field_fourier), INTENT (INout) :: S1,s2,r
    TYPE (c_vector_field_fourier) s3
    INTEGER I,J
    call alloc(s3)

     DO I=-n_fourier,n_fourier
     DO J=-n_fourier,n_fourier
      IF(ABS(I+J)>n_fourier) CYCLE
      S3%F(I+J)=S3%F(I+J) + (s1%f(i).lb.s2%f(j))
     ENDDO
     ENDDO

     r=s3
     call kill(s3)

end subroutine bra_vector_field_fourier

subroutine add_vector_field_fourier(s1,s2,r,fac)
    implicit none
    TYPE (c_vector_field_fourier), INTENT (INout) :: S1,s2,r
    TYPE (c_vector_field_fourier) s3
    complex(dp), optional :: fac
    INTEGER I
    call alloc(s3)

     if(present(fac) ) then
     DO I=-n_fourier,n_fourier

      S3%F(I)=s1%f(i)+fac*s2%f(i)

     ENDDO
     else
     DO I=-n_fourier,n_fourier

      S3%F(I)=s1%f(i)+s2%f(i)

     ENDDO
     endif
     r=s3
     call kill(s3)

end subroutine add_vector_field_fourier

subroutine mulc_vector_field_fourier(s1,fac,r)
    implicit none
    TYPE (c_vector_field_fourier), INTENT (INout) :: S1,r
    TYPE (c_vector_field_fourier) s3
    complex(dp) fac
    INTEGER I
    call alloc(s3)

     DO I=-n_fourier,n_fourier

     S3%F(I)=fac*s1%f(i)

     ENDDO

     r=s3
     call kill(s3)

end subroutine mulc_vector_field_fourier

  SUBROUTINE  c_clean_vector_field_fourier(S1,S2,prec,r)
    implicit none
    type (c_vector_field_fourier),INTENT(INOUT)::S2
    type (c_vector_field_fourier), intent(INOUT):: s1
    real(dp) prec
    type(c_ray), optional :: r
    integer i

    do i=-n_fourier,n_fourier
       call c_clean_vector_field(s1%f(i),s2%f(i),prec,r)
    enddo

  END SUBROUTINE c_clean_vector_field_fourier

  SUBROUTINE  c_clean_taylors(S1,S2,prec,r)
    implicit none
    type (c_taylor),INTENT(INOUT)::S2(:)
    type (c_taylor), intent(INOUT):: s1(:)
    real(dp) prec
    type(c_ray), optional :: r
    integer i

    do i=lbound(s1,1),ubound(s1,1)
       call clean(s1(i),s2(i),prec,r)
    enddo

  END SUBROUTINE c_clean_taylors

  SUBROUTINE  c_evaluate_vector_field_fourier(S1,theta,S2)
    implicit none
    type (c_vector_field_fourier),INTENT(INOUT)::S1
    type (c_vector_field), intent(INOUT):: s2
    real(dp) theta
    integer i

    S2=0
    do i=-n_fourier,n_fourier
       S2=S2+(exp(i_*i*theta))*S1%f(i)
    enddo

  END SUBROUTINE c_evaluate_vector_field_fourier

  subroutine normalise_vector_field_fourier(H,F,K,F1)
    implicit none
    TYPE (c_vector_field_fourier), INTENT (INout) :: H,F,K
    TYPE (c_vector_field_fourier),optional, INTENT (INout) :: F1
    TYPE (c_taylor) temp
 
    integer ki,n,m,j,l,o,nl,i1
    complex(dp), allocatable :: eg(:)
    integer, allocatable :: je(:)
    complex(dp) v,lam
    logical(lp) removeit
    type(c_vector_field_fourier) ht,H1
     IF(.NOT.C_STABLE_DA) then
     RETURN
     endif
     call alloc(temp)
     call alloc(ht)
     call alloc(h1)
     F=0
     H1=H
     ht=H
     n=H%f(0)%n
     i1=2
     if(present(F1)) i1=1
       allocate(eg(n),je(nv)); 
       eg=0.0_dp
       je=0
      do ki=1,n
       if(coast(ki)) then
        eg(ki)=0
       else
        je=0
        je(ki)=1
        eg(ki)=H%f(0)%v(ki).sub.je   ! (1)
       endif  
      enddo

    i1=2 ;if(present(F1)) i1=1 ;nl=0; nl=n_extra;

    do o=i1,no

        ht=H1
       IF(O>1) call exp_vector_field_fourier(F,Ht,Ht)   ! (2)

    do m=-n_fourier,n_fourier
    
      do ki=1,n

       j=1
        do while(.true.) 
          temp=ht%f(m)%v(ki).sub.o
          call  c_cycle(temp,j,v ,je); if(j==0) exit;

     if(m/=0) then
        removeit=.true.
     else
          call check_kernel(ki,n,je,removeit)
     endif
           if(removeit) then
             lam=-i_*m                  ! (3a)
             je(ki)=je(ki)-1
               do l=1,n 
               if(coast(l)) cycle 
                 lam=lam-eg(l)*je(l)    ! (3b)
               enddo
             je(ki)=je(ki)+1
              F%f(m)%v(ki)=F%f(m)%v(ki)-(v.cmono.je)/lam  ! (4)
            endif

        enddo  ! over monomials
       enddo  ! over vector index
      enddo ! over fourier mode
      IF(o==1) THEN   
        call exp_vector_field_fourier(F,Ht,H1,nlin=nl)  ! (5)
        F1=F; F=0;NL=0;
      ENDIF
    enddo  ! over order o
     ht=H1
       call exp_vector_field_fourier(F,Ht,Ht)   ! (6)
     K=ht

   deallocate(eg,je)
   call kill(temp)
   call kill(ht)
   call kill(h1)
 
end  subroutine normalise_vector_field_fourier

 subroutine normalise_vector_field_fourier_factored(H)
    implicit none
    TYPE (c_vector_field_fourier), INTENT (INout) :: H
    TYPE (c_taylor) temp
 
    integer k,n,m,j,l,o,nl
    complex(dp), allocatable :: eg(:)
    integer, allocatable :: je(:)
    complex(dp) v,lam
    logical(lp) removeit
    type(c_vector_field_fourier) Ht,F
     IF(.NOT.C_STABLE_DA) then
     RETURN
     endif

     call alloc(temp)
     call alloc(Ht)
     call alloc(F)

     f=0
     Ht=H
     n=H%f(0)%n

       allocate(eg(n),je(nv)); 
       eg=0.0_dp
       je=0
      do k=1,n
       if(coast(k)) then
        eg(k)=0
       else
        je=0
        je(k)=1
        eg(k)=H%f(0)%v(k).sub.je   ! (1)
       endif  
      enddo

    do o=1,no
    nl=0
    if(o==1) nl=n_extra
!     ht=s1
!     call exp_vector_field_fourier(s2,ht,ht) ! (2)
     F=0
    do m=-n_fourier,n_fourier
                                   ! (2)
      do k=1,n
       j=1
        do while(.true.) 
          temp=Ht%f(m)%v(k).sub.o
          call  c_cycle(temp,j,v ,je); if(j==0) exit;

             if(m/=0) then
                removeit=.true.
             else
                  call check_kernel(k,n,je,removeit)
             endif


           if(removeit) then


             lam=-i_*m                     ! (3a)
             je(k)=je(k)-1
               do l=1,n 
               if(coast(l)) cycle 
                 lam=lam-eg(l)*je(l)       ! (3b)
               enddo
             je(k)=je(k)+1
              F%f(m)%v(k)=F%f(m)%v(k)-(v.cmono.je)/lam  ! (4)
            endif

        enddo  ! over the monomials
       enddo  ! over the components of the vector
      enddo ! over the fourier modes

       call exp_vector_field_fourier(F,Ht,Ht,nlin=nl)   ! (4)

    enddo  ! over order o

   ! ht=s1
   ! call exp_vector_field_fourier(s2,ht,ht)                 ! (5)
   ! s3=ht

   H=Ht                               

   deallocate(eg,je)
   call kill(temp)
   call kill(Ht)
   call kill(F)
 
end  subroutine normalise_vector_field_fourier_factored

subroutine symplectify_for_sethna(m,ms,a1,a2,eps_and_norm)
implicit none
TYPE(c_damap),intent(inout):: m,ms
TYPE(c_damap),optional,intent(inout):: a1,a2
real(dp),optional:: eps_and_norm
TYPE(c_damap) mt,l,b1,b2
type(damap) mm
type(c_vector_field) f,fs
complex(dp) v
type(c_taylor) t,dt
real(dp),allocatable::  mat(:,:), matt(:,:),S(:,:),id(:,:)
integer i,j,k,n(11),nv,nd2,al,ii,a
integer, allocatable :: je(:)
real(dp) dm,norm,normb,norma

call c_get_indices(n,0)
nv=n(4)
nd2=n(3)

call alloc(mt,l,b1,b2);call alloc(f);call alloc(fs);
call alloc(mm);call alloc(t,dt)

b1=1
b2=1
if(present(a1)) b1=a1
if(present(a2)) then
b2=a2
else
 b2=b1
endif
m=b2**(-1)*m*b1
allocate(mat(m%n,m%n),matt(m%n,m%n),S(m%n,m%n),id(m%n,m%n))
mat=0
matt=0
id=0
S=0
do i=1,nd2/2
 S(2*I-1,2*I)=1 ; S(2*I,2*I-1)=-1;
 Id(2*I-1,2*I-1)=1 ; id(2*I,2*I)=1;
enddo


if(nv-nd2==0) then
mt=m
 call c_factor_map(mt,l,f,dir=-1)  
mat=l
else
mat=m
l=mat
mt=l**(-1)*m
f=log(mt)
endif



matt=transpose(mat)

mm=l


! constructing Furman's contracting matrix from my review sec.3.8.2


call checksymp(mm,normb)

do i=1,1000
mat=matmul((1.5e0_dp*id+0.5e0_dp*matmul(mat,matmul(s,matmul(matt,S)))),mat)
matt=transpose(mat)
mm=matt
call checksymp(mm,norm)
if(norm>=normb) exit
 normb=norm
enddo
if(i>=999) then
 Write(6,*) "procedure did not converge in symplectify_for_sethna"
stop 
endif

l=mat
mm=l
fs=0

! Integrating a symplectic operator using the hypercube's diagonal

allocate(je(nv))
je=0
do i=1,f%n

       j=1

        do while(.true.) 

          call  c_cycle(f%v(i),j,v ,je); if(j==0) exit;
         dm=1
         do ii=1,nd2
          dm=dm+je(ii)
         enddo
        t=v.cmono.je
        do a=1,nd2
         dt=t.d.a
        do al=1,nd2
        do k=1,nd2
          fs%v(al)=fs%v(al)+s(a,al)*s(k,i)*(id(k,a)*t+(1.0_dp.cmono.k)*dt)/dm
        enddo ! k
        enddo ! al
        enddo ! a
        enddo

enddo

mt=exp(fs,l)


if(present(eps_and_norm)) then
 m=m-(m.sub.0)
 
 l=m-mt
 norma=0.0_dp
 k=0

do i=1,l%n

       j=1

        do while(.true.) 

          call  c_cycle(l%v(i),j,v ,je); if(j==0) exit;
       
          normb=abs(m%v(i).sub.je)
          norm=abs(v)
  !        if(norm>eps_and_norm) then
 !         norma=norm/normb+norma
           norma=norm+norma
      !    write(16,*) je
      !    write(16,*) norma,norm,normb
          k=k+1
   !       endif


        enddo

enddo
!write(6,*) k
eps_and_norm=norma/k

endif



ms=mt


ms=b2*ms*b1**(-1)
m=b2*m*b1**(-1)


deallocate(mat,matt,S,id)
deallocate(je)
call kill(mt,l,b1,b2); call kill(f);call kill(fs);
call kill(mm);call kill(t,dt)

end subroutine symplectify_for_sethna

subroutine nth_root(m,mr,n)
implicit none
TYPE(c_damap) m,mr
TYPE(c_damap) dm
type(c_vector_field) f
integer, parameter :: nn=36
real(dp) mt(nn,nn),ma(6,6),bb(nn),im(nn,nn),mm(nn,nn)
integer i,n,nt,k,a,b,j
real(dp) xn,norm
integer ind(6,6),k1(nn),k2(nn)

if(n==1) then
mr=m
return
endif
call alloc(dm)
call alloc(f)
k=0
do i=1,c_%nd2
do j=1,c_%nd2
 k=k+1
 ind(i,j)=k
k1(k)=i
k2(k)=j
enddo
enddo


f=c_logf_spin(m)

xn=1.0_dp/n
f=xn*f

dm=exp(f)



nt=nn  

ma=dm


do i=1,c_%nd2
do j=1,c_%nd2
 bb(ind(i,j))=m%e_ij(i,j)
do a=1,c_%nd2
do b=1,c_%nd2
 mt(ind(i,j),ind(a,b))=ma(i,a)*ma(j,b)
enddo
enddo
enddo
enddo

im=0;mm=0;
do i=1,nt
 im(i,i)=1
enddo


do i=0,n-1
mm=im+mm
im=matmul(im,mt)
enddo

 
  call matinv(mm,im,nt,nt,i)
if(i/=0) then
 write(6,*) "Could not inverse in nth_root",i
 stop
endif
bb=matmul(im,bb)
 
do i=1,nt
 dm%e_ij(k1(i),k2(i))=bb(i)
enddo

mr=dm

call kill(f)
call kill(dm)

end subroutine nth_root

!!!!!   stuff for arrays of nodes

 subroutine alloc_node_array_tpsa(a)
 implicit none
 type(node_array), allocatable :: a(:)
 integer i

 do i=1,size(a)
   call alloc(a(i)%f)
   call alloc(a(i)%m)
 enddo

 end  subroutine alloc_node_array_tpsa

 subroutine kill_node_array_tpsa(a)
 implicit none
 type(node_array), allocatable :: a(:)
 integer i

 do i=1,size(a)
   call kill(a(i)%f)
   call kill(a(i)%m)
 enddo

 end  subroutine kill_node_array_tpsa

 subroutine kill_node_array(a)
 implicit none
 type(node_array), allocatable :: a(:)
 integer i

 do i=1,size(a)
   deallocate(a(i)%pos,a(i)%f)
   deallocate(a(i)%v,a(i)%vmax,a(i)%err,a(i)%s)
 enddo

 end  subroutine kill_node_array

 subroutine alloc_node_array(a,n,m)
 implicit none
 type(node_array), allocatable :: a(:)
 integer i,n,m

 allocate(a(n))

 do i=1,n
   allocate(a(i)%pos,a(i)%f,a(i)%m)
   allocate(a(i)%v,a(i)%vmax,a(i)%err,a(i)%s(m))
   a(i)%s=0.0_dp
   a(i)%v=0.0_dp;a(i)%vmax=1.d38;
   a(i)%pos=0
 enddo

 end  subroutine alloc_node_array

!!!!!!!!!!!!!!!!!!   Yu Li Hua  factorization   !!!!!!!!!!!!!!!!!! 

subroutine get_c_yu_w(n,yu,a0,a1,a2,ugiven)
implicit none
type(c_normal_form), intent(inout) ::  n
type(c_yu_w), intent(inout) :: yu
type(c_damap),optional ::  ugiven,a0,a1,a2
type(c_damap) u,b0,b1,bn,ui
type(c_vector_field) f
type(c_taylor) t,p
integer i,j,k
complex(dp), allocatable :: mu(:)
integer, allocatable :: js(:) 

!!! 

allocate(mu(nd2t),js(nd2t))



!nd2t # harmonic planes
if(yu%n/=0) then
 call kill(yu)
endif
call alloc(yu)

 !   

call alloc(u,b0,b1,bn,ui)
call alloc(f);call alloc(t,p);

if(present(ugiven) ) then
 call c_canonise(n%a_t,u,a0=b0,a1=b1,a2=bn)
ui=ci_phasor()*bn*c_phasor()
u=ui**(-1)
f=n%ker
f=u*f
 u=ugiven
else
 call c_canonise(n%a_t,u,a0=b0,a1=b1,a2=bn)
ui=ci_phasor()*bn*c_phasor()
u=ui**(-1)
f=n%ker
f=u*f
endif




do i=1,nd2t/2
js=0
js(2*i-1)=1
 mu(2*i-1)=f%v(2*i-1).sub.js
js=0
js(2*i)=1
 mu(2*i)=f%v(2*i).sub.js
enddo


js=0

do i=1,size(yu%w,1)   ! should be nd2t # harmonic planes 
  yu%w(i,0)=u%v(i)
  t=yu%w(i,0)
  p=1.d0
 do j=1,yu%n
   t=f*t
   p=p*(1.d0+(1.d0.cmono.1))
   yu%w(i,j)=t 
  do k=0,j-1
   js(1)=j-k
  yu%w(i,j)=yu%w(i,j)-(p.sub.js)*mu(i)**(j-k)*yu%w(i,k)
 enddo
enddo
enddo

if(.not.present(ugiven)) then

if(present(a0) ) a0=b0
if(present(a1) ) a1=b1
if(present(a2) ) a2=bn

endif


call kill(u,b0,b1,bn,ui)
call kill(f);call kill(t,p);
deallocate(mu,js)
end subroutine get_c_yu_w

subroutine transform_c_yu_w(yu,a)
implicit none
type(c_yu_w), intent(inout) :: yu
type(c_damap) a 
integer i,j

do i=1,size(yu%w,1)
 do j=0,yu%n
  yu%w(i,j)=yu%w(i,j)*a
 enddo
enddo

end subroutine transform_c_yu_w

!!!!!!!!!!!!!!!!!!   End of Yu Li Hua  factorization   !!!!!!!!!!!!!!!!!! 
subroutine c_fast_canonise(u,u_c,phase,damping,q_c,q_ptc,q_rot,spin_tune ,dospin)
implicit none
type(c_damap), intent(inout) ::  u,u_c
real(dp), optional, intent(inout) :: phase(:),damping(:)
real(dp), optional, intent(inout) :: spin_tune(2)
type(q_linear), optional :: q_c,q_ptc,q_rot
real(dp) b(6,6),b0(6,6),ri(6,6),ang,damp(3),t,cphi,sphi,s(6,6),aq,daq
type(q_linear) q ,qr,qc,qrot
complex(dp) cri(6,6)
integer i
logical dos
logical, optional :: dospin

dos=.false.
if(present(dospin)) dos=dospin
s=0
b0=0
do i=1,nd
b0(2*i-1,2*i-1)=1
b0(2*i,2*i)=1
s(2*i-1,2*i)=1 
s(2*i,2*i-1)=-1 
enddo
if(present(q_rot) ) then 
qrot=0  ! actually makes identity
endif
b=0

ri=0
b=u
 
if(ndpt/=0)  call extract_a0_mat(b,b0)

s=matmul(matmul(b,s),transpose(b))

do i=1,ndt
    damp(i)=sqrt(abs(s(2*i-1,2*i)))
enddo

!det=FindDet(b(1:nd2,1:nd2), nd2)**(1.0_dp/nd2)
!write(6,*) damp

 
      do i=1,ndt
      if((i<=ndt)) then  !.or.(i>nd-rf)) then
 !    damp(i)=sqrt(b(2*i-1,2*i-1)*b(2*i,2*i)-b(2*i-1,2*i)*b(2*i,2*i-1))
       if(courant_snyder_teng_edwards) then
        t=sqrt(b(2*i-1,2*i-1)**2+b(2*i-1,2*i)**2)
        cphi=b(2*i-1,2*i-1)/t
        sphi=b(2*i-1,2*i)/t
       else
        t=sqrt(b(2*i,2*i-1)**2+b(2*i,2*i)**2)
        cphi=b(2*i,2*i)/t
        sphi=-b(2*i,2*i-1)/t
       endif
       ri(2*i-1,2*i-1)=cphi /damp(i)
       ri(2*i,2*i)=cphi/damp(i)
       ri(2*i-1,2*i)=-sphi /damp(i)
       ri(2*i,2*i-1)=sphi /damp(i)
  
    endif
 
 if(present(phase)) then
     ang=-atan2(sphi,cphi)
  phase(i)=phase(i)-ang/twopi
 endif
  if(present(damping)) then
  damping(i)=damping(i)-log(damp(i))
 endif
if(present(q_rot) ) then 
 qrot%mat(2*i-1,2*i-1)=ri(2*i-1,2*i-1)*damp(i)**2
 qrot%mat(2*i,2*i)=ri(2*i,2*i)*damp(i)**2
 qrot%mat(2*i-1,2*i)=-ri(2*i-1,2*i)*damp(i)**2
 qrot%mat(2*i,2*i-1)=-ri(2*i,2*i-1)*damp(i)**2
endif
      enddo



      if(ndpt/=0) then
        ri(5,5)=1
        ri(6,6)=1
        ri(ndptb,ndpt)=- b(ndptb,ndpt)
        if(mod(ndpt,2)==0) then
         i=ndpt/2
        else
         i=ndptb/2
        endif
if(present(phase))       phase(i)=phase(i)+b(ndptb,ndpt)
if(present(q_rot) ) then 
 qrot%mat(ndptb,ndpt)=-ri(ndptb,ndpt)
 qrot%mat(5,5)=ri(5,5)
 qrot%mat(6,6)=ri(6,6)
endif
      endif

       s=matmul(b,ri)
       s=matmul(b0,s)
    u_c=s
if(use_quaternion.and.dos) then
q=1
 q=u%q
!  make sure isf not below y plane
       aq=q%q(0,0)**2-(q%q(1,0)**2+q%q(2,0)**2+q%q(3,0)**2)
if(aq<0) then
          qr=1   !  = i 
          qr=0.0_dp
          qr%q(1,0)=1.0_dp
          q=q*qr
endif
 qr=1
 qr=0.0_dp

 aq=-atan2(real(q%q(2,0)),real(q%q(0,0)))

 qr%q(0,0)= cos(aq)
 qr%q(2,0)= sin(aq)

 daq=0
 if(ndpt/=0) then  

  daq=(q%q(0,ndpt)*qr%q(2,0)+q%q(2,ndpt)*qr%q(0,0))/(q%q(2,0)*qr%q(2,0)-q%q(0,0)*qr%q(0,0))
  qr%q(0,ndpt)=  -daq*qr%q(2,0)
  qr%q(2,ndpt)= daq*qr%q(0,0)
 endif
 qc=q*qr
 if(present(spin_tune)) then
  spin_tune(1)=spin_tune(1)+aq/pi   ! changed 2018.11.01
 endif
cri=ri
qc=qc*cri
 u_c%q=qc 
if(present(q_c) ) then 
q_c=1 
cri=inv_symplectic66(s)
q_c=qc*cri
endif

endif
 
 if(present(spin_tune)) then
  spin_tune(2)=spin_tune(2)+daq/pi   ! changed 2018.11.01
 endif

if(present(q_ptc) ) then 
q_ptc=u_c 
endif
if(present(q_rot) ) then 
 qrot%q=qr%q 
 qrot%q(1:3,0:6)=-qrot%q(1:3,0:6)
endif

if(present(q_rot) ) then 
 q_rot=qrot*q_rot
endif
 
end subroutine c_fast_canonise

subroutine extract_a0_mat(a,a0)
!#internal: manipulation
!# This routines extracts a0: the full fixed point map.
    implicit none
    real(dp), intent(inout) :: a(6,6),a0(6,6)
    type(c_damap) aa,aa0
    call alloc(aa,aa0)
    aa=a
    aa0=a0
    call extract_a0(aa,aa0)
    a=aa
    a0=aa0
      call kill(aa,aa0)
end subroutine extract_a0_mat




!Function to find the determinant of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Description: The subroutine is based on two key points:
!1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!row operations (column operations would work as well) are used
!to convert the matrix into upper traingular form
!2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
!
REAL(dp) FUNCTION FindDet(mat, n)
    IMPLICIT NONE
    REAL(dp), DIMENSION(n,n) :: matrix,mat
    INTEGER, INTENT(IN) :: n
    REAL(dp) :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    matrix=mat
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
    
    !Calculate determinant by finding product of diagonal elements
    FindDet = l
    DO i = 1, n
        FindDet = FindDet * matrix(i,i)
    END DO
    
END FUNCTION FindDet

  END MODULE  c_tpsa
