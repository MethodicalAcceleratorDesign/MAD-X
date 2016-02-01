module madx_ptc_normal_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  ! madx_ptc_normal module
  ! Frank Schmidt (CERN)
  ! updates by Piotr Skowronski (CERN)
  use madx_ptc_module
  use c_TPSA
!  use madx_ptc_twiss
  
  !_________________________________________________________________
  implicit none
  save

  public :: ptc_normal

  
  private display_table_results
  
  character(1000), private  :: whymsg
  
contains
  !____________________________________________________________________________________________


  SUBROUTINE ptc_normal()
    USE ptc_results
    USE madx_ptc_intstate_module
    implicit none
    logical(lp) closed_orbit,normal,maptable
    integer no,mynd2,npara,nda,icase,flag_index,why(9)
    integer i,ii,j1,k,l,starti
    integer n_rows,row,n_haml,n_gnfu,nres,mynres,n1,n2
    integer,external :: select_ptc_idx, minimum_acceptable_order, &
         string_from_table_row, double_from_table_row, double_to_table_row
    real(dp) x(6),deltap0,deltap,dt
    integer :: indexa(4)
    integer :: row_haml(101)
    integer :: index1(1000,2)
    real(kind(1d0)) get_value,val_ptc,tmp
    character(len = 5) name_var
    type(real_8) y(6)
    type(real_8) :: theAscript(6) ! used here to compute dispersion's derivatives
    type(c_damap)  :: c_Map
    type(c_vector_field) vf_kernel, vf
    type(c_taylor)  :: g_io
    type(c_normal_form), target :: theNormalForm, theNormalForm_t2 ! normal from type 2 for hamiltonian terms
    
    !--------------------------------------------------------;----------------------

    if(universe.le.0.or.EXCEPTION.ne.0) then
       call fort_warn('return from ptc_normal: ',' no universe created')
       return
    endif
    if(index_mad.le.0.or.EXCEPTION.ne.0) then
       call fort_warn('return from ptc_normal: ',' no layout created')
       return
    endif
    
    ! do the check early so it leaves immediately and does not leave allocated objects
    no = get_value('ptc_normal ','no ')
    
    if (no < minimum_acceptable_order()) then
      print*, "minimum acceptable order: ", minimum_acceptable_order()
      call seterrorflag(11,"ptc_normal ","Order of calculation is not sufficient to calculate required parameters")
      call fort_warn('ptc_normal: ',&
           'Order of calculation (parameter no) is not sufficient to calculate required parameters')
      return
    endif
    
    call cleartables() !defined in madx_ptc_knobs
    
    nda=0

    icase = get_value('ptc_normal ','icase ')
    deltap0 = get_value('ptc_normal ','deltap ')
    
    !
    deltap = zero
    call my_state(icase,deltap,deltap0)
    CALL UPDATE_STATES

    x(:)=zero
    if(mytime) then
       call Convert_dp_to_dt (deltap, dt)
    else
       dt=deltap
    endif
    if(icase.eq.5) x(5)=dt
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Closed orbit search
    !!    if requested
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    closed_orbit = get_value('ptc_normal ','closed_orbit ') .ne. 0
    if(closed_orbit) then
      ! pass starting point for closed orbit search
       x(1)=get_value('ptc_normal ','x ')
       x(2)=get_value('ptc_normal ','px ')
       x(3)=get_value('ptc_normal ','y ')
       x(4)=get_value('ptc_normal ','py ')
       x(6)=get_value('ptc_normal ','t ')
       x(5)=x(5)+get_value('ptc_normal ','pt ')

       call find_orbit(my_ring,x,1,default,c_1d_7)
       CALL write_closed_orbit(icase,x)
    endif


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! INIT PTC
    !!    and allocate map
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !new complex PTC
    call init_all(default,no,nda,BERZ,mynd2,npara)
    c_verbose=.false.
    
    i_piotr(:) = 0
    i_piotr(1)=1; i_piotr(2)=3; i_piotr(3)=5;
    
    c_normal_auto=.true.;

    call alloc(y)
    y=npara
    Y=X

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Track the map
    !!    and allocate map
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    c_%watch_user=.true.
    call track(my_ring,y,1,default)
    if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
       write(whymsg,*) 'DA got unstable in tracking: PTC msg: ',messagelost
       call fort_warn('ptc_normal: ',whymsg)
       call seterrorflag(10,"ptc_normal ",whymsg)
       return
    endif
    call PRODUCE_APERTURE_FLAG(flag_index)
    if(flag_index/=0) then
       call ANALYSE_APERTURE_FLAG(flag_index,why)
       write(whymsg,*) "out of aperture ",why
       call fort_warn('ptc_normal: ',whymsg)
       call seterrorflag(10,"ptc_normal ",whymsg)
       
       CALL kill(y)
       c_%watch_user=.false.
       return
    endif
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Write the map 
    !!    to fort.18
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    c_%watch_user=.false.
    !if (getdebug()>1)
    print77=.false.
    read77 =.false.

    call daprint(y,18)


    maptable = get_value('ptc_normal ','maptable ') .ne. 0
    if(maptable) then
       call makemaptable(y,no)
    endif


    normal = get_value('ptc_normal ','normal ') .ne. 0
    if(normal) then

       !------ Find the number of occurences of the attribute 'haml'

       n_rows = select_ptc_idx()
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !! 
       !! HERE WE DO NORMAL FORM TYPE 1
       !!     Q (+chroma+anharmon) and DX
       !! 
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       call alloc(theNormalForm) 
       call alloc(c_Map)
   
       c_Map = y;

       if (getdebug() > 0) print*,"Normal Form Type 1"
       
       call  c_normal(c_Map,theNormalForm)       ! (4) !! HERE WE DO NORMAL FORM TYPE 1
       
       if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
          write(whymsg,*) 'DA got unstable in Normal Form: PTC msg: ',messagelost
          call fort_warn('ptc_normal: ',whymsg)
          call seterrorflag(10,"ptc_normal ",whymsg);
          return
       endif

       if (getdebug() > 1) then
          write(19,'(/a/)') 'NORMAL: Dispersion, First and Higher Orders'
          call daprint(theNormalForm%A1,19)
       endif
       

       call alloc(vf_kernel)
       vf_kernel=0
       call flatten_c_factored_lie(theNormalForm%ker,vf_kernel)

       if (getdebug() > 1) then
          write(19,'(/a/)') 'NORMAL:  Tunes, Chromaticities and Anharmonicities'

          call daprint(theNormalForm%A_t,19)! WARNING, this does not work with 6D dispersion, need to use Chao-Sands prescription
          call daprint(vf_kernel,19) ! tunes


       endif


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
       !! 
       !! HERE WE DO NORMAL FORM TYPE 2
       !!       HAML and GNFU
       !! It has resonanses left in the normal form
       !! which are defined with normal%M

       call alloc(theNormalForm_t2) 
       
       nres = 0       
       n_haml = 0
       n_gnfu = 0
       
       ! read HAML's and GNUF's that user asks for
       if (n_rows > 0) then
          do row = 1,n_rows
             name_var=" "
             k = string_from_table_row('normal_results ', 'name ', row, name_var)


             if (name_var(:4) .eq. 'gnfu') then
                if (getdebug() > 2) print*,"ptc_normal: adding gnfu"
                n_gnfu = n_gnfu + 1
             endif
             if (name_var(:4) .eq. 'haml') then
                n_haml = n_haml + 1
                row_haml(n_haml) = row
             endif
          enddo

          if (n_haml > 0) then
             
             !call alloc(pbrh)
             
             do j1 =1,n_haml
                row = row_haml(j1)
                k = double_from_table_row("normal_results ", "value ", row, doublenum)
                mynres = int(doublenum)
                row = row_haml(j1) - 3*mynres + 2
                starti = 1
                
                if (j1 .eq. 1) then
                   k = double_from_table_row("normal_results ", "order1 ", row, doublenum) !ordrer in X
                   indexa(1) = int(doublenum)
                   k = double_from_table_row("normal_results ", "order2 ", row, doublenum) !ordrer in PX
                   indexa(2) = int(doublenum)
                   k = double_from_table_row("normal_results ", "order3 ", row, doublenum) !ordrer in Y
                   indexa(3) = int(doublenum)
                   k = double_from_table_row("normal_results ", "order4 ", row, doublenum) !ordrer in PY
                   indexa(4) = int(doublenum)

                   index1(1,1) = indexa(1) - indexa(2)
                   index1(1,2) = indexa(3) - indexa(4)
                   
                   
                   !n_t2%m(1,1)= index1(1,1)
                   !n_t2%m(2,1)= index1(1,2)
                   !if(c_%nd2.eq.6) n_t2%m(3,1)= indexa(3)
                   call fort_warn('ptc_normal: ',' Hamiltonian terms not yet implemented with the new complex PTC')
                   
                   
                   nres = 1
                   starti = 2
                endif
                !============ nres is the number of resonances to be set

                if (mynres .ge. starti) then
                   do i = starti,mynres
                      ii = row + 3*(i-1)
                      k = double_from_table_row("normal_results ", "order1 ", ii, doublenum)
                      indexa(1) = int(doublenum)
                      k = double_from_table_row("normal_results ", "order2 ", ii, doublenum)
                      indexa(2) = int(doublenum)
                      k = double_from_table_row("normal_results ", "order3 ", ii, doublenum)
                      indexa(3) = int(doublenum)
                      k = double_from_table_row("normal_results ", "order4 ", ii, doublenum)
                      indexa(4) = int(doublenum)
                      n1 = indexa(1) - indexa(2)
                      n2 = indexa(3) - indexa(4)
                      do l = 1,nres
                         if (n1 .eq. index1(l,1) .and. n2 .eq. index1(l,2)) goto 100
                      enddo
                      nres = nres + 1
                      index1(nres,1) = n1
                      index1(nres,2) = n2
                      !n_t2%m(1,nres)= n1
                      !n_t2%m(2,nres)= n2
                      !if(c_%nd2.eq.6) n_t2%m(3,nres)= indexa(3)
	  call fort_warn('ptc_normal: ',' Hamiltonian terms not yet implemented with the new complex PTC')
	  
100                   continue
                   enddo
                endif
             enddo
             
             !n_t2%nres = nres
             call fort_warn('ptc_normal: ',' Hamiltonian terms not yet implemented with the new complex PTC')
          endif

       endif
       !------------------------------------------------------------------------
       
       if (n_haml > 0) then
           
           if (getdebug() > 0) print*,"Normal Form Type 2 (for Hamiltonian Terms)"
           
           call fort_warn('ptc_normal: ',' Hamiltonian terms not yet implemented with the new complex PTC')

           !call alloc(theNormalForm_t2)
           !call  c_normal(c_Map,theNormalForm)


           if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
              write(whymsg,*) 'DA got unstable in Normal Form: PTC msg: ',messagelost
              call fort_warn('ptc_normal: ',whymsg)
              call seterrorflag(10,"ptc_normal ",whymsg)
              return
           endif


       endif
       
       if (n_gnfu > 0) then

          !!!!!!!!!!!!!!!!!!!!!!
          !Generating functions

          !n%g is the vecotor field for the transformation 
          !from resonance basis (action angle coordinate system, (x+ipx),(x-ipx)) back to cartesion X,Y
          !the ndim polynomials need to be flattened to get RDT's
          call alloc(vf);
          call flatten_c_factored_lie(theNormalForm%G,vf)
          
          call alloc(g_io);
          
          call equal_c_tayls(g_io,-cgetpb(vf))
          
          !g_io =-cgetpb(vf)
          !call putGnormaltable(g_io)

       endif
       
       !------ get values and store them in the table 'normal_results' ---------
       if (n_rows > 0) then
          do row = 1,n_rows
             name_var=" "
             k = string_from_table_row("normal_results ", "name ", row, name_var)

             if (name_var(:3) .eq. 'ham') then
               val_ptc = double_from_normal_t2(name_var, row, icase)  !!! HERE IS THE RETRIVAL
             else
               val_ptc = double_from_normal_t1(name_var, row, icase)  !!! HERE IS THE RETRIVAL
             endif
             
             if (name_var(:4) .ne. 'haml' .and. name_var(:4) .ne. 'gnfu') then
               k = double_to_table_row("normal_results ", "value ", row, val_ptc)
             endif  
             
          enddo
       endif
       
       call alloc(theAscript)
       theAscript = X+theNormalForm%A_t;
       
       call putusertable(1,'$end$ ',dt ,dt,y,theAscript)
       call kill(theAscript)
       


       if (n_gnfu > 0) then
          call kill(g_io)
          call kill(vf)
       endif   

       !if (n_haml > 0) then 
         !clean after hamiltonian when implemented
       !endif  
 

       call kill(theNormalForm)
       !call kill(n_t2)

    endif
    

    
    CALL kill(y)

    close(18)
    
    if (getdebug() > 1) then
      close(19)
    endif  

! f90flush is not portable, and useless...
!    call f90flush(18,my_false)
!    call f90flush(19,my_false)
  
  contains
  !________________________________________________________________________________

     FUNCTION double_from_normal_t1(name_var,row,icase)
       USE ptc_results
       implicit none
       logical(lp) name_l
       integer,intent(IN) ::  row,icase
       real(dp) double_from_normal_t1, d_val, d_val1, d_val2
       complex(dp) :: c_val
       integer ii,i1,i2,jj
       integer j,k,ind(6)
       integer double_from_table_row
       character(len = 4)  name_var
       character(len = 2)  name_var1
       character(len = 3)  name_var2

       name_l = .false.
       double_from_normal_t1 = zero

       name_var1 = name_var
       SELECT CASE (name_var1)
       CASE ("dx")
          ind(:)=0
          k = double_from_table_row("normal_results ", "order1 ", row, doublenum)
          ind(5) = int(doublenum)
          if (ind(5) == 0) ind(5) = 1
          ind(6) = 0
          d_val = theNormalForm%A1%V(1).sub.ind
       CASE ('dy')
          ind(:)=0
          k = double_from_table_row("normal_results ", "order1 ", row, doublenum)
          ind(5) = int(doublenum)
          if (ind(5) == 0) ind(5) = 1
          ind(6) = 0
          d_val = theNormalForm%A1%V(3).sub.ind
       CASE ('q1')
          ind(:)=0
          ind(1)=1 !kernel has extra exponent at the plane variable, q1=v(1).sub.'100000'
          c_val=vf_kernel%v(1).sub.ind 
          d_val = -aimag(c_val)/(2.*pi)

       CASE ('q2')
          ind(:)=0
          ind(3)=1 !kernel has extra exponent at the plane variable, q1=v(1).sub.'100000'
          c_val=vf_kernel%v(3).sub.ind 
          d_val = -aimag(c_val)/(2.*pi)
       CASE DEFAULT
          name_l = .true.
       END SELECT
       if (name_l) then
          name_l = .false.
          name_var2 = name_var
          SELECT CASE (name_var2)
          CASE ('dpx')
             ind(:)=0
             k = double_from_table_row("normal_results ", "order1 ", row, doublenum)
             ind(5) = int(doublenum)
             if (ind(5) == 0) ind(5) = 1
             ind(6) = 0
             d_val = theNormalForm%A1%V(2).sub.ind
          CASE ('dpy')
             ind(:)=0
             k = double_from_table_row("normal_results ", "order1 ", row, doublenum)
             ind(5) = int(doublenum)
             if (ind(5) == 0) ind(5) = 1
             ind(6) = 0
             d_val = theNormalForm%A1%V(4).sub.ind
          CASE ('dq1')
             k = double_from_table_row("normal_results ", "order1 ", row, doublenum)
             ind(:)=0
             ind(1)=1
             ind(5) = int(doublenum)
             if (ind(5) == 0) ind(5) = 1
             ind(6) = 0
             c_val=vf_kernel%v(1).sub.ind 
             d_val = -aimag(c_val)/(2.*pi)
          CASE ('dq2')
             k = double_from_table_row("normal_results ", "order1 ", row, doublenum)
             ind(:)=0
             ind(3)=1
             ind(5) = int(doublenum)
             if (ind(5) == 0) ind(5) = 1
             ind(6) = 0
             c_val=vf_kernel%v(3).sub.ind 
             d_val = -aimag(c_val)/(2.*pi)
          CASE DEFAULT
             name_l = .true.
          END SELECT
       endif
       if (name_l) then
          SELECT CASE (name_var)
          CASE ('anhx')
             k = double_from_table_row("normal_results ", "order1 ", row, doublenum)
             do j = 1,2
                ind(j) =  int(doublenum)
             enddo
             k = double_from_table_row("normal_results ", "order2 ", row, doublenum)
             do j = 3,4
                ind(j) = int(doublenum)
             enddo
             k = double_from_table_row("normal_results ", "order3 ", row, doublenum)
             ind(5) = int(doublenum)
             ind(6) = 0
             ind(1) = ind(1) + 1
             c_val=vf_kernel%v(1).sub.ind 
             d_val = -aimag(c_val)/(2.*pi)
             
             
          CASE ('anhy')
             k = double_from_table_row("normal_results ", "order1 ", row, doublenum)
             do j = 1,2
                ind(j) = int(doublenum)
             enddo
             k = double_from_table_row("normal_results ", "order2 ", row, doublenum)
             do j = 3,4
                ind(j) = int(doublenum)
             enddo
             k = double_from_table_row("normal_results ", "order3 ", row, doublenum)
             ind(5) = int(doublenum)
             ind(6) = 0
             ind(3) = ind(3) + 1
             c_val=vf_kernel%v(3).sub.ind 
             d_val = -aimag(c_val)/(2.*pi)
          CASE ('eign')
             ii=(icase/2)*2
             k = double_from_table_row("normal_results ", "order1 ", row, doublenum)
             i1 = int(doublenum)
             if(i1.gt.ii) call aafail('return from double_from_normal_t1: ',' wrong # of eigenvectors')
             jj=0
             if(i1.eq.5) jj=6
             if(i1.eq.6) i1=5
             if(jj.eq.6) i1=jj
             k = double_from_table_row("normal_results ", "order2 ", row, doublenum)
             i2 = int(doublenum)
             if(i2.gt.ii) call aafail('return from double_from_normal_t1: ',' eigenvectors too many components')
             jj=0
             if(i2.eq.5) jj=6
             if(i2.eq.6) i2=5
             if(jj.eq.6) i2=jj
             ind(:)=0
             ind(i2)=1
             double_from_normal_t1 = theNormalForm%A_t%V(i1).sub.ind
             if(mytime.and.i2.eq.6) double_from_normal_t1 = -double_from_normal_t1
             RETURN
           CASE ('gnfc')
              k = double_from_table_row("normal_results ", "order1 ", row, doublenum)
              ind(1) = int(doublenum)
              k = double_from_table_row("normal_results ", "order2 ", row, doublenum)
              ind(2) = int(doublenum)
              k = double_from_table_row("normal_results ", "order3 ", row, doublenum)
              ind(3) = int(doublenum)
              k = double_from_table_row("normal_results ", "order4 ", row, doublenum)
              ind(4) = int(doublenum)
              ind(5) = 0
              ind(6) = 0
              c_val = g_io.sub.ind
              d_val = -real(c_val)
              double_from_normal_t1 = d_val
              RETURN
           CASE ('gnfs')
              k = double_from_table_row("normal_results ", "order1 ", row, doublenum)
              ind(1) = int(doublenum)
              k = double_from_table_row("normal_results ", "order2 ", row, doublenum)
              ind(2) = int(doublenum)
              k = double_from_table_row("normal_results ", "order3 ", row, doublenum)
              ind(3) = int(doublenum)
              k = double_from_table_row("normal_results ", "order4 ", row, doublenum)
              ind(4) = int(doublenum)
              ind(5) = 0
              ind(6) = 0
              c_val = g_io.sub.ind
              d_val = -imag(c_val)
              double_from_normal_t1 = d_val
              RETURN
           CASE ('gnfa')
              k = double_from_table_row("normal_results ", "order1 ", row, doublenum)
              ind(1) = int(doublenum)
              k = double_from_table_row("normal_results ", "order2 ", row, doublenum)
              ind(2) = int(doublenum)
              k = double_from_table_row("normal_results ", "order3 ", row, doublenum)
              ind(3) = int(doublenum)
              k = double_from_table_row("normal_results ", "order4 ", row, doublenum)
              ind(4) = int(doublenum)
              ind(5) = 0
              ind(6) = 0
              c_val = g_io.sub.ind
              d_val1 = imag(c_val)
              d_val2 = real(c_val)
              double_from_normal_t1 = SQRT(d_val1**2 + d_val2**2)
              RETURN
           CASE ('gnfu')
              double_from_normal_t1 = zero ! should never arrive here
              RETURN
          CASE DEFAULT
             print *,"--Error in the table normal_results-- Unknown input: ",name_var
          END SELECT
       endif
       double_from_normal_t1 = d_val*(factorial(ind(1))*factorial(ind(3))*factorial(ind(5)))
     END FUNCTION double_from_normal_t1
   !________________________________________________
   !Extraction of normal type 2 variables: hamiltonian and generating functions

     FUNCTION double_from_normal_t2(name_var,row,icase)
       USE ptc_results
       implicit none
       logical(lp) name_l
       integer,intent(IN) ::  row,icase
       real(dp) double_from_normal_t2, d_val, d_val1, d_val2
       integer ii,i1,i2,jj
       integer j,k,ind(6)
       integer double_from_table_row
       character(len = 4)  name_var
       character(len = 2)  name_var1
       character(len = 3)  name_var2

       double_from_normal_t2 = zero

       name_var1 = name_var

      ! if (LEN_TRIM(name_var))

       SELECT CASE (name_var)
       CASE ('hamc')
          k = double_from_table_row("normal_results ", "order1 ", row, doublenum)
          ind(1) = int(doublenum)
          k = double_from_table_row("normal_results ", "order2 ", row, doublenum)
          ind(2) = int(doublenum)
          k = double_from_table_row("normal_results ", "order3 ", row, doublenum)
          ind(3) = int(doublenum)
          k = double_from_table_row("normal_results ", "order4 ", row, doublenum)
          ind(4) = int(doublenum)
          ind(5) = 0
          ind(6) = 0
          !d_val = pbrh%cos%h.sub.ind
          d_val = zero
          call fort_warn('double_from_normal_t2: ',' Hamiltonian terms not yet implemented with the new complex PTC')
          double_from_normal_t2 = d_val
          RETURN
       CASE ('hams')
          k = double_from_table_row("normal_results ", "order1 ", row, doublenum)
          ind(1) = int(doublenum)
          k = double_from_table_row("normal_results ", "order2 ", row, doublenum)
          ind(2) = int(doublenum)
          k = double_from_table_row("normal_results ", "order3 ", row, doublenum)
          ind(3) = int(doublenum)
          k = double_from_table_row("normal_results ", "order4 ", row, doublenum)
          ind(4) = int(doublenum)
          ind(5) = 0
          ind(6) = 0
          !d_val = pbrh%sin%h.sub.ind
          d_val = zero
          call fort_warn('double_from_normal_t2: ',' Hamiltonian terms not yet implemented with the new complex PTC')
          double_from_normal_t2 = d_val
          RETURN
       CASE ('hama')
          k = double_from_table_row("normal_results ", "order1 ", row, doublenum)
          ind(1) = int(doublenum)
          k = double_from_table_row("normal_results ", "order2 ", row, doublenum)
          ind(2) = int(doublenum)
          k = double_from_table_row("normal_results ", "order3 ", row, doublenum)
          ind(3) = int(doublenum)
          k = double_from_table_row("normal_results ", "order4 ", row, doublenum)
          ind(4) = int(doublenum)
          ind(5) = 0
          ind(6) = 0
          !d_val1 = pbrh%cos%h.sub.ind
          !d_val2 = pbrh%sin%h.sub.ind
          !double_from_normal_t2 = SQRT(d_val1**2 + d_val2**2)
          double_from_normal_t2 = zero
          call fort_warn('double_from_normal_t2: ',' Hamiltonian terms not yet implemented with the new complex PTC')
          
          RETURN
       CASE ('haml')
          double_from_normal_t2 = zero
          RETURN
       CASE DEFAULT
          print *,"--Error in the table normal_results-- Unknown input: ",name_var
       END SELECT

       double_from_normal_t2 = d_val*(factorial(ind(1))*factorial(ind(3))*factorial(ind(5)))
     END FUNCTION double_from_normal_t2

  END subroutine ptc_normal
!_________________________________________________________________

  SUBROUTINE display_table_results()
    implicit none
    integer, external :: select_ptc_idx, string_from_table_row, double_from_table_row
    character(len = 4) name_var
    integer row,k
    integer :: ord(3)

    print *,"Variable name  Order 1  order 2  order 3        Value      "
    do row = 1 , select_ptc_idx()
       name_var=" "
       k = string_from_table_row("normal_results ", "name ", row, name_var)
       k = double_from_table_row("normal_results ", "order1 ", row, doublenum)
       ord(1) = int(doublenum)
       k = double_from_table_row("normal_results ", "order2 ", row, doublenum)
       ord(2) = int(doublenum)
       k = double_from_table_row("normal_results ", "order3 ", row, doublenum)
       ord(3) = int(doublenum)
       k = double_from_table_row("normal_results ", "value ", row, doublenum)
       WRITE(*,100) name_var,ord(1),ord(2),ord(3),doublenum
    enddo
100 FORMAT(3X,A4,14X,I1,8X,I1,8X,I1,5X,f25.18)
  END SUBROUTINE display_table_results
END MODULE madx_ptc_normal_module
