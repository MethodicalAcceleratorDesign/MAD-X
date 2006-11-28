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

integer, parameter                ::   kn_x  =28 !kn prefix is added to limit crashes with the other modules name space
integer, parameter                ::   kn_px =29
integer, parameter                ::   kn_y  =30
integer, parameter                ::   kn_py =31
integer, parameter                ::   kn_t  =32
integer, parameter                ::   kn_dp =33



integer, parameter                ::   ntwisses = 33

CHARACTER*(6), parameter          ::   twissnames(ntwisses) = (/ &
     'beta11','beta12','beta13',&
     'beta21','beta22','beta23',&
     'beta31','beta32','beta33',&
     'alfa11','alfa12','alfa13',&
     'alfa21','alfa22','alfa23',&
     'alfa31','alfa32','alfa33',&
     'gama11','gama12','gama13',&
     'gama21','gama22','gama23',&
     'gama31','gama32','gama33',&
     'x     ','px    ',&
     'y     ','py    ',&
     'pt    ','t     '          /)
