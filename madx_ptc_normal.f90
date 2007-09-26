module madx_ptc_normal_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  ! madx_ptc_normal module
  ! Frank Schmidt (CERN)
  !
  use madx_ptc_module
  !_________________________________________________________________
  implicit none

  public

  private double_from_ptc_normal,display_table_results

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
         string_from_table, double_from_table, result_from_normal
    real(dp) x(6),deltap0,deltap,dt
    integer :: indexa(4)
    integer :: row_haml(101)
    integer :: index1(1000,2)
    real(kind(1d0)) get_value,val_ptc
    character(len = 5) name_var
    !------------------------------------------------------------------------------

    if(universe.le.0) then
       call fort_warn('return from ptc_normal: ',' no universe created')
       return
    endif
    if(index_mad.le.0) then
       call fort_warn('return from ptc_normal: ',' no layout created')
       return
    endif

    nda=0

    icase = get_value('ptc_normal ','icase ')
    deltap0 = get_value('ptc_normal ','deltap ')

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
    closed_orbit = get_value('ptc_normal ','closed_orbit ') .ne. 0
    if(closed_orbit) then
       call find_orbit(my_ring,x,1,default,c_1d_7)
       CALL write_closed_orbit(icase,x)
    endif

    no = get_value('ptc_normal ','no ')

    call init(default,no,nda,BERZ,mynd2,npara)


    call alloc(y)
    y=npara
    Y=X

    c_%watch_user=.true.
    call track(my_ring,y,1,default)
    if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
       call fort_warn('ptc_normal: ','DA got unstable')
       call seterrorflag(10,"ptc_normal ","DA got unstable ");
       return
    endif
    call PRODUCE_APERTURE_FLAG(flag_index)
    if(flag_index/=0) then
       call ANALYSE_APERTURE_FLAG(flag_index,why)
       Write(6,*) "ptc_normal unstable (map production)-programs continues "
       Write(6,*) why ! See produce aperture flag routine in sd_frame
       CALL kill(y)
       c_%watch_user=.false.
       return
    endif
    c_%watch_user=.false.
    !if (getdebug()>1)
    print77=.false.
    read77 =.false.

    call daprint(y,18)

    maptable = get_value('ptc_normal ','maptable ') .ne. 0
    if(maptable) then
       call makemaptable(y)
    endif


    normal = get_value('ptc_normal ','normal ') .ne. 0
    if(normal) then
       call alloc(n)

       !------ Find the number of occurences of the attribute 'haml'

       n_rows = select_ptc_idx()
       n_haml = 0
       n_gnfu = 0

       if (n_rows > 0) then
          do row = 1,n_rows
             name_var=" "
             k = string_from_table('normal_results ', 'name ', row, name_var)


             if (name_var(:4) .eq. 'gnfu') n_gnfu = n_gnfu + 1
             if (name_var(:4) .eq. 'haml') then
                n_haml = n_haml + 1
                row_haml(n_haml) = row
             endif
          enddo


          if (n_gnfu > 0) call alloc(pbrg)

          if (n_haml > 0) then
             call alloc(pbrh)
             do j1 =1,n_haml
                row = row_haml(j1)
                k = double_from_table("normal_results ", "value ", row, doublenum)
                mynres = int(doublenum)
                row = row_haml(j1) - 3*mynres + 2
                starti = 1
                if (j1 .eq. 1) then
                   k = double_from_table("normal_results ", "order1 ", row, doublenum)
                   indexa(1) = int(doublenum)
                   k = double_from_table("normal_results ", "order2 ", row, doublenum)
                   indexa(2) = int(doublenum)
                   k = double_from_table("normal_results ", "order3 ", row, doublenum)
                   indexa(3) = int(doublenum)
                   k = double_from_table("normal_results ", "order4 ", row, doublenum)
                   indexa(4) = int(doublenum)
                   index1(1,1) = indexa(1) - indexa(2)
                   index1(1,2) = indexa(3) - indexa(4)
                   n%m(1,1)= index1(1,1)
                   n%m(2,1)= index1(1,2)
                   if(ndim2.eq.6) n%m(3,1)= indexa(3)
                   nres = 1
                   starti = 2
                endif
                !============ nres is the number of resonances to be set

                if (mynres .ge. starti) then
                   do i = starti,mynres
                      ii = row + 3*(i-1)
                      k = double_from_table("normal_results ", "order1 ", ii, doublenum)
                      indexa(1) = int(doublenum)
                      k = double_from_table("normal_results ", "order2 ", ii, doublenum)
                      indexa(2) = int(doublenum)
                      k = double_from_table("normal_results ", "order3 ", ii, doublenum)
                      indexa(3) = int(doublenum)
                      k = double_from_table("normal_results ", "order4 ", ii, doublenum)
                      indexa(4) = int(doublenum)
                      n1 = indexa(1) - indexa(2)
                      n2 = indexa(3) - indexa(4)
                      do l = 1,nres
                         if (n1 .eq. index1(l,1) .and. n2 .eq. index1(l,2)) goto 100
                      enddo
                      nres = nres + 1
                      index1(nres,1) = n1
                      index1(nres,2) = n2
                      n%m(1,nres)= n1
                      n%m(2,nres)= n2
                      if(ndim2.eq.6) n%m(3,nres)= indexa(3)
100                   continue
                   enddo
                endif
             enddo
             n%nres = nres
          endif

       endif
       !------------------------------------------------------------------------


       n=y
       if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
          call fort_warn('ptc_normal: ','Fatal Error: DA in NormalForm got unstable')
          stop
       endif
       if (n_gnfu > 0) pbrg = n%a%pb
       if (n_haml > 0) pbrh = n%normal%pb
       if (getdebug() > 1) then
          write(19,'(/a/)') 'Dispersion, First and Higher Orders'
          call daprint(n%A1,19)
       endif

       !------ get values and store them in the table 'normal_results' ---------


       n_rows = select_ptc_idx()
       if (no < minimum_acceptable_order()) then
          print*, "minimum acceptable order: ", minimum_acceptable_order()
          call seterrorflag(11,"ptc_normal ","Order of calculation is not sufficient to calculate required parameters")
          call fort_warn('ptc_normal: ',&
               'Order of calculation (parameter no) is not sufficient to calculate required parameters')
          return
       endif

       if (n_rows > 0) then
          do row = 1,n_rows
             name_var=" "
             k = string_from_table("normal_results ", "name ", row, name_var)
             val_ptc = double_from_ptc_normal(name_var,row,icase)
             if (name_var(:4) .ne. 'haml'.and.name_var(:4) .ne. 'gnfu')    &
                  call double_to_table_row("normal_results ", "value ", row, val_ptc)
          enddo
       endif

       !------------------------------------------------------------------------

       if (getdebug() > 1) then
          write(19,'(/a/)') 'Tunes, Chromaticities and Anharmonicities'
          !  call daprint(n%A_t,19)
          !  call daprint(n%A,19)
          call daprint(n%dhdj,19)
          !  call daprint(pbrh,19)
       endif

       if (n_gnfu > 0) call kill(pbrg)
       if (n_haml > 0) call kill(pbrh)
       call kill(n)
    endif
    CALL kill(y)
    call f90flush(18,my_false)
    call f90flush(19,my_false)

  END subroutine ptc_normal
  !________________________________________________________________________________

  FUNCTION double_from_ptc_normal(name_var,row,icase)
    USE ptc_results
    implicit none
    logical(lp) name_l
    integer,intent(IN) ::  row,icase
    real(dp) double_from_ptc_normal, d_val, d_val1, d_val2
    integer ii,i1,i2,jj
    integer j,k,ind(6)
    integer double_from_table
    character(len = 4)  name_var
    character(len = 2)  name_var1
    character(len = 3)  name_var2

    name_l = .false.
    double_from_ptc_normal = zero

    name_var1 = name_var
    SELECT CASE (name_var1)
    CASE ("dx")
       ind(:)=0
       k = double_from_table("normal_results ", "order1 ", row, doublenum)
       ind(5) = int(doublenum)
       if (ind(5) == 0) ind(5) = 1
       ind(6) = 0
       d_val = n%A1%V(1).sub.ind
    CASE ('dy')
       ind(:)=0
       k = double_from_table("normal_results ", "order1 ", row, doublenum)
       ind(5) = int(doublenum)
       if (ind(5) == 0) ind(5) = 1
       ind(6) = 0
       d_val = n%A1%V(3).sub.ind
    CASE ('q1')
       ind(:)=0
       d_val = n%dhdj%V(3).sub.ind
    CASE ('q2')
       ind(:)=0
       d_val = n%dhdj%V(4).sub.ind
    CASE DEFAULT
       name_l = .true.
    END SELECT
    if (name_l) then
       name_l = .false.
       name_var2 = name_var
       SELECT CASE (name_var2)
       CASE ('dpx')
          ind(:)=0
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(5) = int(doublenum)
          if (ind(5) == 0) ind(5) = 1
          ind(6) = 0
          d_val = n%A1%V(2).sub.ind
       CASE ('dpy')
          ind(:)=0
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(5) = int(doublenum)
          if (ind(5) == 0) ind(5) = 1
          ind(6) = 0
          d_val = n%A1%V(4).sub.ind
       CASE ('dq1')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(:)=0
          ind(5) = int(doublenum)
          if (ind(5) == 0) ind(5) = 1
          ind(6) = 0
          d_val = n%dhdj%V(3).sub.ind
       CASE ('dq2')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(:)=0
          ind(5) = int(doublenum)
          if (ind(5) == 0) ind(5) = 1
          ind(6) = 0
          d_val = n%dhdj%V(4).sub.ind
       CASE DEFAULT
          name_l = .true.
       END SELECT
    endif
    if (name_l) then
       SELECT CASE (name_var)
       CASE ('anhx')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          do j = 1,2
             ind(j) =  int(doublenum)
          enddo
          k = double_from_table("normal_results ", "order2 ", row, doublenum)
          do j = 3,4
             ind(j) = int(doublenum)
          enddo
          k = double_from_table("normal_results ", "order3 ", row, doublenum)
          ind(5) = int(doublenum)
          ind(6) = 0
          d_val = n%dhdj%V(3).sub.ind
       CASE ('anhy')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          do j = 1,2
             ind(j) = int(doublenum)
          enddo
          k = double_from_table("normal_results ", "order2 ", row, doublenum)
          do j = 3,4
             ind(j) = int(doublenum)
          enddo
          k = double_from_table("normal_results ", "order3 ", row, doublenum)
          ind(5) = int(doublenum)
          ind(6) = 0
          d_val = n%dhdj%V(4).sub.ind
       CASE ('hamc')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(1) = int(doublenum)
          k = double_from_table("normal_results ", "order2 ", row, doublenum)
          ind(2) = int(doublenum)
          k = double_from_table("normal_results ", "order3 ", row, doublenum)
          ind(3) = int(doublenum)
          k = double_from_table("normal_results ", "order4 ", row, doublenum)
          ind(4) = int(doublenum)
          ind(5) = 0
          ind(6) = 0
          d_val = pbrh%cos%h.sub.ind
          double_from_ptc_normal = d_val
          RETURN
       CASE ('hams')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(1) = int(doublenum)
          k = double_from_table("normal_results ", "order2 ", row, doublenum)
          ind(2) = int(doublenum)
          k = double_from_table("normal_results ", "order3 ", row, doublenum)
          ind(3) = int(doublenum)
          k = double_from_table("normal_results ", "order4 ", row, doublenum)
          ind(4) = int(doublenum)
          ind(5) = 0
          ind(6) = 0
          d_val = pbrh%sin%h.sub.ind
          double_from_ptc_normal = d_val
          RETURN
       CASE ('hama')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(1) = int(doublenum)
          k = double_from_table("normal_results ", "order2 ", row, doublenum)
          ind(2) = int(doublenum)
          k = double_from_table("normal_results ", "order3 ", row, doublenum)
          ind(3) = int(doublenum)
          k = double_from_table("normal_results ", "order4 ", row, doublenum)
          ind(4) = int(doublenum)
          ind(5) = 0
          ind(6) = 0
          d_val1 = pbrh%cos%h.sub.ind
          d_val2 = pbrh%sin%h.sub.ind
          double_from_ptc_normal = SQRT(d_val1**2 + d_val2**2)
          RETURN
       CASE ('haml')
          double_from_ptc_normal = zero
          RETURN
       CASE ('gnfc')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(1) = int(doublenum)
          k = double_from_table("normal_results ", "order2 ", row, doublenum)
          ind(2) = int(doublenum)
          k = double_from_table("normal_results ", "order3 ", row, doublenum)
          ind(3) = int(doublenum)
          k = double_from_table("normal_results ", "order4 ", row, doublenum)
          ind(4) = int(doublenum)
          ind(5) = 0
          ind(6) = 0
          d_val = pbrg%cos%h.sub.ind
          double_from_ptc_normal = d_val
          RETURN
       CASE ('gnfs')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(1) = int(doublenum)
          k = double_from_table("normal_results ", "order2 ", row, doublenum)
          ind(2) = int(doublenum)
          k = double_from_table("normal_results ", "order3 ", row, doublenum)
          ind(3) = int(doublenum)
          k = double_from_table("normal_results ", "order4 ", row, doublenum)
          ind(4) = int(doublenum)
          ind(5) = 0
          ind(6) = 0
          d_val = pbrg%sin%h.sub.ind
          double_from_ptc_normal = d_val
          RETURN
       CASE ('gnfa')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(1) = int(doublenum)
          k = double_from_table("normal_results ", "order2 ", row, doublenum)
          ind(2) = int(doublenum)
          k = double_from_table("normal_results ", "order3 ", row, doublenum)
          ind(3) = int(doublenum)
          k = double_from_table("normal_results ", "order4 ", row, doublenum)
          ind(4) = int(doublenum)
          ind(5) = 0
          ind(6) = 0
          d_val1 = pbrg%cos%h.sub.ind
          d_val2 = pbrg%sin%h.sub.ind
          double_from_ptc_normal = SQRT(d_val1**2 + d_val2**2)
          RETURN
       CASE ('gnfu')
          double_from_ptc_normal = zero
          RETURN
       CASE ('eign')
          ii=(icase/2)*2
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          i1 = int(doublenum)
          if(i1.gt.ii) call aafail('return from double_from_ptc_normal: ',' wrong # of eigenvectors')
          jj=0
          if(i1.eq.5) jj=6
          if(i1.eq.6) i1=5
          if(jj.eq.6) i1=jj
          k = double_from_table("normal_results ", "order2 ", row, doublenum)
          i2 = int(doublenum)
          if(i2.gt.ii) call aafail('return from double_from_ptc_normal: ',' eigenvectors too many components')
          jj=0
          if(i2.eq.5) jj=6
          if(i2.eq.6) i2=5
          if(jj.eq.6) i2=jj
          ind(:)=0
          ind(i2)=1
          double_from_ptc_normal = n%A_t%V(i1).sub.ind
          if(mytime.and.i2.eq.6) double_from_ptc_normal = -double_from_ptc_normal
          RETURN
       CASE DEFAULT
          print *,"--Error in the table normal_results-- Unknown input: ",name_var
       END SELECT
    endif
    double_from_ptc_normal = d_val*(factorial(ind(1))*factorial(ind(3))*factorial(ind(5)))
  END FUNCTION double_from_ptc_normal

  SUBROUTINE display_table_results()
    implicit none
    integer, external :: select_ptc_idx, string_from_table, double_from_table
    character(len = 4) name_var
    integer row,k
    integer :: ord(3)

    print *,"Variable name  Order 1  order 2  order 3        Value      "
    do row = 1 , select_ptc_idx()
       name_var=" "
       k = string_from_table("normal_results ", "name ", row, name_var)
       k = double_from_table("normal_results ", "order1 ", row, doublenum)
       ord(1) = int(doublenum)
       k = double_from_table("normal_results ", "order2 ", row, doublenum)
       ord(2) = int(doublenum)
       k = double_from_table("normal_results ", "order3 ", row, doublenum)
       ord(3) = int(doublenum)
       k = double_from_table("normal_results ", "value ", row, doublenum)
       WRITE(*,100) name_var,ord(1),ord(2),ord(3),doublenum
    enddo
100 FORMAT(3X,A4,14X,I1,8X,I1,8X,I1,5X,f25.18)
  END SUBROUTINE display_table_results
END MODULE madx_ptc_normal_module
