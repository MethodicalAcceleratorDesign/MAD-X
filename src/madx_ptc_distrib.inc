
type momentdef
   integer        :: iarray(6) !defiens the moment as integer array f.g. 200000 is <x*x>
   integer        :: index     ! column in 2D array of taylors of the result, zero if not parameteric
   character(20)  :: table     !column name, if not provided it is "iarray(1)-iarray(2)-iarray(3)"
   character(17)  :: column    !column name, if not provided it is "iarray(1)-iarray(2)-iarray(3)"
end type momentdef


integer, parameter                ::   distr_gauss=1
integer, parameter                ::   distr_flat5=2
integer, parameter                ::   distr_flat56=3
integer, parameter                ::   distr_user=4
