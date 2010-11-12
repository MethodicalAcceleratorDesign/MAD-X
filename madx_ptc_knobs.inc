!_____________________________________________________________



type tablepush_poly
   character(20) :: tabname  !table name  to put the coefficients
   logical       :: pushtab  ! if true it should be pushed to table
   character(20) :: colname  !column name to put the coefficients
   integer       :: index ! column in 2D array of taylors of the result
   integer       :: element  !defines the element of 6D function
   character(10) :: monomial !defiens the monomial f.g. 100000 is coeff of x, and 020300 is coeff of px^2*py^3
end type tablepush_poly
!describes result

type resultdescr
   integer       :: index ! column in 2D array of taylors of the result
   character(10) :: name ! name of the parameter defined by the user
end type resultdescr
!_____________________________________________________________

type mapelement
   integer       :: element  !defines the element of 6D function
   character(10) :: monomial !defiens the monomial f.g. 100000 is coeff of x, and 020300 is coeff of px^2*py^3
end type mapelement
!_____________________________________________________________

type mapelresult
   type (resultdescr)  :: descr
   type (mapelement)   :: r
end type mapelresult
!_____________________________________________________________


!defines a file
type file
   character(20) :: name  !file name
   integer       :: fd !file descriptor
end type file
!_____________________________________________________________

type dumpmap_file
   type(file)    :: file
   type(mapelement)   :: r
end type dumpmap_file


integer, parameter                ::   beta11=1
integer, parameter                ::   beta12=2
integer, parameter                ::   beta13=3

integer, parameter                ::   beta21=4
integer, parameter                ::   beta22=5
integer, parameter                ::   beta23=6

integer, parameter                ::   beta31=7
integer, parameter                ::   beta32=8
integer, parameter                ::   beta33=9

integer, parameter                ::   alfa11=10
integer, parameter                ::   alfa12=11
integer, parameter                ::   alfa13=12

integer, parameter                ::   alfa21=13
integer, parameter                ::   alfa22=14
integer, parameter                ::   alfa23=15

integer, parameter                ::   alfa31=16
integer, parameter                ::   alfa32=17
integer, parameter                ::   alfa33=18

integer, parameter                ::   gama11=19
integer, parameter                ::   gama12=20
integer, parameter                ::   gama13=21

integer, parameter                ::   gama21=22
integer, parameter                ::   gama22=23
integer, parameter                ::   gama23=24

integer, parameter                ::   gama31=25
integer, parameter                ::   gama32=26
integer, parameter                ::   gama33=27

! additional parameters for the Twiss derivatives w.r.t delta_p:

integer, parameter  :: beta11p = 28
integer, parameter  ::  beta12p = 29
integer, parameter  :: beta13p = 30

integer, parameter  :: beta21p = 31
integer, parameter  :: beta22p = 32
integer, parameter  :: beta23p = 33

integer, parameter  :: beta31p = 34
integer, parameter  :: beta32p = 35
integer, parameter  :: beta33p = 36

integer, parameter  :: alfa11p = 37
integer, parameter  :: alfa12p = 38
integer, parameter  :: alfa13p = 39

integer, parameter  :: alfa21p = 40
integer, parameter  :: alfa22p = 41
integer, parameter  :: alfa23p = 42

integer, parameter  :: alfa31p = 43
integer, parameter  :: alfa32p = 44
integer, parameter  :: alfa33p = 45

integer, parameter  :: gama11p = 46
integer, parameter  :: gama12p = 47
integer, parameter  :: gama13p = 48

integer, parameter  :: gama21p = 49
integer, parameter  :: gama22p = 50
integer, parameter  :: gama23p = 51

integer, parameter  :: gama31p = 52
integer, parameter  :: gama32p = 53
integer, parameter  :: gama33p = 54

!---

integer, parameter                ::   disp1=55
integer, parameter                ::   disp2=56
integer, parameter                ::   disp3=57
integer, parameter                ::   disp4=58

integer, parameter                ::   disp1p = 59
integer, parameter                ::   disp2p = 60
integer, parameter                ::   disp3p = 61
integer, parameter                ::   disp4p = 62

! second derivatives of dispersion
integer, parameter		  ::   disp1p2 = 63
integer, parameter		  ::   disp2p2 = 64
integer, parameter 		  ::   disp3p2 = 65
integer, parameter                ::   disp4p2 = 66
! third order derivatives of dispersion
integer, parameter		  ::   disp1p3 = 67
integer, parameter		  ::   disp2p3 = 68
integer, parameter                ::   disp3p3 = 69
integer, parameter                ::   disp4p3 = 70


integer, parameter                ::   mu1=71 ! 63+8 dispersion derivatives
integer, parameter                ::   mu2=72
integer, parameter                ::   mu3=73 ! 65+8 dispersion derivatives

integer, parameter                ::   kn_x  =74 !kn prefix is added to limit crashes with the other modules name space
integer, parameter                ::   kn_px =75
integer, parameter                ::   kn_y  =76
integer, parameter                ::   kn_py =77
integer, parameter                ::   kn_t  =78
integer, parameter                ::   kn_dp =79 ! btw nothing to do with ntwisses



integer, parameter                ::   ntwisses = 79 ! 73 till mu3 + 6 phase-space variables
integer, parameter                ::   ninparams = 27 	!possible initial conditions as parameters
							! (unused)



! increase string size from 6 to 7 to accomodate Twiss derivatives w.r.t deltap
CHARACTER(7), parameter          ::   twissnames(ntwisses) = (/ &
     'beta11 ','beta12 ','beta13 ',& ! space-padding to reach 7
     'beta21 ','beta22 ','beta23 ',&
     'beta31 ','beta32 ','beta33 ',&
     'alfa11 ','alfa12 ','alfa13 ',&
     'alfa21 ','alfa22 ','alfa23 ',&
     'alfa31 ','alfa32 ','alfa33 ',&
     'gama11 ','gama12 ','gama13 ',&
     'gama21 ','gama22 ','gama23 ',&
     'gama31 ','gama32 ','gama33 ',&
     'beta11p','beta12p','beta13p',&
     'beta21p','beta22p','beta23p',&
     'beta31p','beta32p','beta33p',&
     'alfa11p','alfa12p','alfa13p',&
     'alfa21p','alfa22p','alfa23p',&
     'alfa31p','alfa32p','alfa33p',&
     'gama11p','gama12p','gama13p',&
     'gama21p','gama22p','gama23p',&    
     'gama31p','gama32p','gama33p',&
     'disp1  ','disp2  ','disp3  ','disp4  ',&
     'disp1p ','disp2p ','disp3p ','disp4p ',&
     'disp1p2','disp2p2','disp3p2','disp4p2',&
     'disp1p3','disp2p3','disp3p3','disp4p3',&
     'mu1    ','mu2    ','mu3    ',&
     'x      ','px     ',&
     'y      ','py     ',&
     'pt     ','t      '/)
