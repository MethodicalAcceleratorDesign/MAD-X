!/*
! * Copyright(C) 2008 by Lingyun Yang
! * lingyun(.dot.]yang@gmail.com
! * http://www.lingyunyang.com
! *
! * Please get permission from Lingyun Yang before you redistribute this file.
! *
! * Version: $Id: c_tpsa_interface.F90,v 1.13 2009-06-24 09:06:43 frs Exp $
! */


module dabnew
  !$$$$  use da_arrays
  use dabnew_b !$$$$
  implicit none
  public
  ! integer,private,parameter:: lsw=1

  ! integer,private,parameter::nmax=400,lsw=1
  ! real(dp),private,parameter::tiny=c_1d_20
  character(120),private :: line

#ifdef _WIN32_DLL
  !DEC$ ATTRIBUTES DLLIMPORT :: ad_tra, ad_shift, ad_print, ad_save_block, ad_read_block
  !DEC$ ATTRIBUTES DLLIMPORT :: ad_fill_ran, ad_nvar, ad_length, ad_derivative
  !DEC$ ATTRIBUTES DLLIMPORT :: ad_subst, ad_cos, ad_sin, ad_log, ad_exp, ad_sqrt, ad_abs
  !DEC$ ATTRIBUTES DLLIMPORT :: ad_div_c, ad_c_div, ad_mult_const, ad_add_const
  !DEC$ ATTRIBUTES DLLIMPORT :: ad_div, ad_mult, ad_sub, ad_reset, ad_pok, ad_pek
  !DEC$ ATTRIBUTES DLLIMPORT :: ad_var, ad_truncate, ad_const, ad_count, ad_free, ad_add
  !DEC$ ATTRIBUTES DLLIMPORT :: ad_copy, ad_clean, ad_alloc, ad_reserve, ad_init, ad_elem

  DLL_IMPORT ad_tra, ad_shift, ad_print, ad_save_block, ad_read_block
  DLL_IMPORT ad_fill_ran, ad_nvar, ad_length, ad_derivative
  DLL_IMPORT ad_subst, ad_cos, ad_sin, ad_log, ad_exp, ad_sqrt, ad_abs
  DLL_IMPORT ad_div_c, ad_c_div, ad_mult_const, ad_add_const
  DLL_IMPORT ad_div, ad_mult, ad_sub, ad_reset, ad_pok, ad_pek
  DLL_IMPORT ad_var, ad_truncate, ad_const, ad_count, ad_free, ad_add
  DLL_IMPORT ad_copy, ad_clean, ad_alloc, ad_reserve, ad_init, ad_elem ,ad_resetvars
#endif

  private trx_cpp


contains

  !done
  subroutine daini(nd,nv,k)
    implicit none
    !    integer, intent(in) :: nd,nv,nd2,k
    integer nd,nv,k,last_nv
    if(lingyun_yang) then !%%%%
       if(last_tpsa==1.and.lielib_print(10)==0) then
          call ad_nvar(last_nv)
          call ad_resetvars(last_nv)
       endif
       call danum0(nd,nv)
       call ad_init(nv, nd)
       call ad_reserve(lda)
       !     write(6,*) " LDA Lingyun ", lda
       last_tpsa=1
    else !%%%%
       call daini_b(nd,nv,k)
       !           last_tpsa=2   done in daini_b
    endif !%%%%
  end subroutine daini

  !done
  subroutine daall0(i)
    implicit none
    !    integer, intent(in) :: i
    integer  i
    if(lingyun_yang) then !%%%%
       call ad_alloc(i)
    else !%%%%
       call daall0_b(i)
    endif !%%%%
  end subroutine daall0

  !done
  subroutine daall1(i,ccc,n0,nv)
    implicit none
    !    integer, intent(in) :: i
    integer i,n0, nv
    character(10) ccc

    if(lingyun_yang) then !%%%%
       call ad_alloc(i)
    else !%%%%
       call daall1_b(i,ccc,n0,nv)
    endif !%%%%

  end subroutine daall1

  subroutine daclean(ina,value)
    implicit none
    integer ina
    !
    real(dp) value
    !
    !
    if(lingyun_yang) then !%%%%
       call ad_clean(ina, value)
    else !%%%%
       call daclean_b(ina, value)
    endif !%%%%
    return
  end subroutine daclean

  !done
  subroutine daadd(i,j,k)
    implicit none
    !    integer, intent(in) :: i,j,k
    integer  i,j,k
    integer itmp
    !write(*,*) k
    ! if j!= k and i != k, then optimize
    if(lingyun_yang) then !%%%%
       call ad_alloc(itmp)
       call ad_copy(i, itmp)
       call ad_add(itmp, j)
       call ad_copy(itmp, k)
       call ad_free(itmp)
    else !%%%%
       call daadd_b(i,j,k)
    endif !%%%%
  end subroutine daadd

  !done
  subroutine dacon(i,r)
    implicit none
    !    integer, intent(in) :: i
    !    real(dp), intent(in) :: r
    integer  i
    real(dp)  r

    if(lingyun_yang) then !%%%%
       call ad_const(i, r)
    else !%%%%
       call dacon_b(i,r)
    endif !%%%%
  end subroutine dacon

  !done
  subroutine dadal(idal,j)
    implicit none
    integer  j
    !    integer, intent(in) :: j
    integer,dimension(:)::idal
    integer k
    if(lingyun_yang) then !%%%%
       do k=1,j
          call ad_free(idal(k))
       enddo
    else !%%%%
       call dadal_b(idal,j)
    endif !%%%%
  end subroutine dadal

  !done
  subroutine dadal1(idal)
    implicit none
    integer  idal
    !    integer, intent(inout) :: idal
    if(lingyun_yang) then !%%%%
       call ad_free(idal)
    else !%%%%
       call dadal1_b(idal)
    endif !%%%%
  end subroutine dadal1

  !done
  subroutine count_da(idal)
    implicit none
    !    integer, intent(inout) :: idal
    integer  idal
    !call ad_all(i)
    if(lingyun_yang) then !%%%%
       call ad_count(idal)
    else !%%%%
       call count_da_b(idal)
    endif !%%%%
  end subroutine count_da

  !done
  subroutine davar(ina,ckon,i)
    implicit none
    !    integer, intent(in) :: ina,i
    integer ina,i
    real(dp), intent(in) :: ckon
    if(lingyun_yang) then !%%%%
       call ad_var(ina, ckon, i-1)
    else !%%%%
       call davar_b(ina,ckon,i)
    endif !%%%%

  end subroutine davar

  ! done
  subroutine danot(not)
    implicit none
    !    integer, intent(in) :: not
    integer  not
    if(lingyun_yang) then !%%%%
       print *, 'ERROR: This is not used in new TPSA routines.'
       STOP
    else !%%%%
       call danot_b(not)
    endif !%%%%
  end subroutine danot

  !done
  subroutine datrunc(ina, imd, inb)
    implicit none
    !    integer, intent(in) :: ina,imd,inb
    integer  ina,imd,inb
    if(lingyun_yang) then !%%%%
       call ad_copy(ina, inb)
       ! truncate imd and above.
       call ad_truncate(inb, imd)
    else !%%%%
       call datrunc_b(ina, imd, inb)
    endif !%%%%
  end subroutine datrunc

  ! done
  subroutine daeps(deps)
    implicit none
    real(dp) deps
    !    real(dp), intent(inout) :: deps

    if(lingyun_yang) then !%%%%
       print *, 'ERROR: We use machine dependent eps instead.'
       if (deps<zero) then
          deps = eps
       else
          STOP
       endif
    else !%%%%
       call daeps_b(deps)
    endif !%%%%

  end subroutine daeps

  ! done
  subroutine dapek(ina,jv,cjj)
    implicit none
    !    integer,intent(in) :: ina
    integer  ina
    real(dp),intent(inout) :: cjj
    integer,dimension(:)::jv     ! 2002.12.4
    if(lingyun_yang) then !%%%%
       call ad_pek(ina, jv, size(jv), cjj)
    else !%%%%
       call dapek_b(ina,jv,cjj)
    endif !%%%%
  end subroutine dapek

  !done
  subroutine dapok(ina,jv,cjj)
    implicit none
    integer  ina
    !    integer,intent(in) :: ina
    real(dp),intent(in) :: cjj
    integer,dimension(:)::jv     ! 2002.12.4
    if(lingyun_yang) then !%%%%
       call ad_pok(ina, jv, size(jv), cjj)
    else !%%%%
       call dapok_b(ina,jv,cjj)
    endif !%%%%
  end subroutine dapok

  !done
  subroutine daclr(ina)
    implicit none
    !    integer,intent(in) :: ina
    integer  ina
    if(lingyun_yang) then !%%%%
       call ad_reset(ina)
    else !%%%%
       call daclr_b(ina)
    endif !%%%%
  end subroutine daclr

  !done
  subroutine dacop(ina,inb)
    implicit none
    !    integer,intent(in) :: ina,inb
    integer  ina,inb
    if(lingyun_yang) then !%%%%
       call ad_copy(ina, inb)
    else !%%%%
       call dacop_b(ina,inb)
    endif !%%%%
  end subroutine dacop

  !done
  subroutine dasub(ina,inb,inc)
    implicit none
    !    integer,intent(in) :: ina,inb,inc
    integer  ina,inb,inc
    integer i1,i2
    if(lingyun_yang) then !%%%%
       call ad_alloc(i1)
       call ad_alloc(i2)
       call ad_copy(ina, i1)
       call ad_copy(inb, i2)
       call ad_sub(i1, i2)
       call ad_copy(i1, inc)
       call ad_free(i1)
       call ad_free(i2)
    else !%%%%
       call dasub_b(ina,inb,inc)
    endif !%%%%
  end subroutine dasub

  !done
  subroutine damul(ina,inb,inc)
    implicit none
    !    integer,intent(in) :: ina,inb,inc
    integer  ina,inb,inc
    integer i1, i2
    if(lingyun_yang) then !%%%%
       call ad_alloc(i1)
       call ad_alloc(i2)
       call ad_copy(ina, i1)
       call ad_copy(inb, i2)
       call ad_mult(i1, i2, inc)
       call ad_free(i1)
       call ad_free(i2)
    else !%%%%
       call damul_b(ina,inb,inc)
    endif !%%%%
  end subroutine damul

  !done
  subroutine dadiv(ina,inb,inc)
    implicit none
    !    integer,intent(in) :: ina,inb,inc
    integer  ina,inb,inc
    integer i1, i2

    if(lingyun_yang) then !%%%%
       call ad_alloc(i1)
       call ad_alloc(i2)
       call ad_copy(ina, i1)
       call ad_copy(inb, i2)
       call ad_div(i1, i2, inc)
       call ad_free(i1)
       call ad_free(i2)
    else !%%%%
       call dadiv_b(ina,inb,inc)
    endif !%%%%

  end subroutine dadiv

  !done
  subroutine dacad(ina,ckon,inc)
    implicit none
    integer ina,inc
    !    integer,intent(in) :: ina,inc
    real(dp), intent(in):: ckon
    if(lingyun_yang) then !%%%%
       call ad_copy(ina, inc)
       call ad_add_const(inc, ckon)
    else !%%%%
       call dacad_b(ina,ckon,inc)
    endif !%%%%
  end subroutine dacad

  !done
  subroutine dacsu(ina,ckon,inc)
    implicit none
    !    integer,intent(in) :: ina,inc
    integer  ina,inc
    real(dp)  ckon
    call dacad(ina, -ckon, inc)
  end subroutine dacsu

  !done
  subroutine dasuc(ina,ckon,inc)
    implicit none
    integer  ina,inc
    real(dp)  ckon
    if(lingyun_yang) then !%%%%
       call ad_copy(ina, inc)
       call ad_mult_const(inc, -1.d0)
       call dacad(inc, ckon, inc)
    else !%%%%
       call dasuc_b(ina,ckon,inc)
    endif !%%%%
  end subroutine dasuc

  !done
  subroutine dacmu(ina,ckon,inc)
    implicit none
    integer  ina,inc
    integer it1
    real(dp)  ckon
    if(lingyun_yang) then !%%%%
       call ad_alloc(it1)
       call ad_copy(ina, it1)
       call ad_mult_const(it1, ckon)
       call ad_copy(it1, inc)
       call ad_free(it1)
    else !%%%%
       call dacmu_b(ina,ckon,inc)
    endif !%%%%
  end subroutine dacmu

  !done
  subroutine dadic(ina,ckon,inc)
    implicit none
    integer ina,inc
    real(dp)  ckon
    integer i1

    if(lingyun_yang) then !%%%%
       call ad_alloc(i1)
       call ad_c_div(ina, ckon, i1)
       call ad_copy(i1, inc)
       call ad_free(i1)
    else !%%%%
       call dadic_b(ina,ckon,inc)
    endif !%%%%
  end subroutine dadic

  !done
  subroutine dacdi(ina,ckon,inc)
    implicit none
    integer  ina,inc
    real(dp)  ckon
    integer i1

    if(lingyun_yang) then !%%%%
       call ad_alloc(i1)
       call ad_copy(ina, i1)
       call ad_div_c(i1, ckon)
       call ad_copy(i1, inc)
       call ad_free(i1)
    else !%%%%
       call dacdi_b(ina,ckon,inc)
    endif !%%%%
  end subroutine dacdi

  !done
  subroutine daabs(ina,r)
    implicit none
    integer  ina
    real(dp)  r
    if(lingyun_yang) then !%%%%
       call ad_abs(ina, r)
    else !%%%%
       call daabs_b(ina,r)
    endif !%%%%

  end subroutine daabs

  subroutine dalin(ina,afac,inb,bfac,inc)
    implicit none
    integer  ina,inb,inc
    integer it1,it2
    real(dp)  afac,bfac
    !write(*,*) "Not implemented"
    !call exit

    if(lingyun_yang) then !%%%%
       call ad_alloc(it1)
       call ad_alloc(it2)
       call ad_copy(ina,it1)
       call ad_copy(inb,it2)
       call ad_mult_const(it1,afac)
       call ad_mult_const(it2,bfac)
       call ad_copy(it1, inc)
       call ad_add(inc, it2)
       call ad_free(it1)
       call ad_free(it2)
    else !%%%%
       call dalin_b(ina,afac,inb,bfac,inc)
    endif !%%%%
  end subroutine dalin

  subroutine dafun(cf,ina,inc)
    implicit none
    integer  ina,inc
    character(4)  cf
    INTEGER IPAUSE,MYPAUSES
    if(lingyun_yang) then !%%%%
       select case (cf)
       case('INV ')
          !call exit
          call ad_c_div(ina, 1.d0, inc)
       case('SQRT')
          call ad_sqrt(ina, inc)
       case('EXP ')
          call ad_exp(ina, inc)
       case('LOG ')
          call ad_log(ina, inc)
       case('SIN ')
          call ad_sin(ina, inc)
       case('COS ')
          call ad_cos(ina, inc)
       case('SINH')
          write(*,*) "BUG in sinh"
          STOP
       case('COSH')
          write(*,*) "BUG in cosh"
          STOP
       case default
          write(line,'(a28,1x,a4)')  'ERROR, UNSOPPORTED FUNCTION ',cf
          ipause=mypauses(35,line)
       end select
       !
    else !%%%%
       call dafun_b(cf,ina,inc)
    endif !%%%%
  end subroutine dafun

  subroutine dacct(ma,ia,mb,ib,mc,ic)
    implicit none
    integer,dimension(:)::ma,mb,mc
    integer  ia,ib,ic
    integer, allocatable :: c(:)
    integer i

    if(lingyun_yang) then !%%%%

       if((size(ma)<ia).or.(size(mb)<ib).or.(size(mc)<ic))then
          write(6,*) "Error caught in interface dacct 1"
          stop
       endif
       if(ia/=ic)then
          write(6,*) "Error caught in interface dacct 2"
          stop
       endif

       !    if(ib/=nv)then
       !     write(6,*) "Error caught in interface dacct 3"
       !    stop
       !    endif

       allocate(c(ic))

       do i=1,ic
          call ad_alloc(c(i))
       enddo

       do i=1,ia
          !write(*,*) ib
          call ad_subst(ma(i), mb, ib, c(i))
       enddo

       do i=1,ic
          call ad_copy(c(i),mc(i))
          call ad_free(c(i))
       enddo

       deallocate(c)

    else !%%%%
       call dacct_b(ma,ia,mb,ib,mc,ic)
    endif !%%%%


  end subroutine dacct

  subroutine mtree(mb,ib,mc,ic)
    implicit none
    integer,dimension(:)::mb,mc
    integer  ib,ic
    integer i
    if(lingyun_yang) then !%%%%
       if(ib/=ic) then
          write(6,*) "In mtree not compatible to C++ TPSA "
          stop 666
       endif
       do i=1,ib
          call dacop(mb(i),mc(i))
       enddo
    else !%%%%
       call  mtree_b(mb,ib,mc,ic)
    endif !%%%%

  end subroutine mtree

  subroutine ppushstore(mc,nd2,coef,ml,mv)
    implicit none
    !
    integer nd2
    integer,dimension(:), intent(in)::mc
    integer,dimension(:), intent(out)::ml,mv
    real(dp),dimension(:),intent(out)::coef


    if(lingyun_yang) then !%%%%
       mv=0
       ml=0
       coef=0
       write(*,*) "ppushstore should be called using the LBNL Version of Berz TPSA"
       STOP 666
    else !%%%%
       call   ppushstore_b(mc,nd2,coef,ml,mv)
    endif !%%%%
  end subroutine ppushstore

  subroutine ppushGETN(mc,ND2,ntot)
    implicit none
    !
    integer ntot,ND2
    integer,dimension(:), intent(inout)::mc
    !
    if(lingyun_yang) then !%%%%
       write(*,*) "ppushGETN should be called using the LBNL Version of Berz TPSA"
       STOP 666
    else !%%%%
       call   ppushGETN_b(mc,ND2,ntot)
    endif !%%%%
  end subroutine ppushGETN

  subroutine ppush1(mc,xi,xf)
    implicit none
    integer mc
    real(dp) xf
    real(dp),dimension(:)::xi

    if(lingyun_yang) then !%%%%
       call trx_cpp(mc,xi,xf)
    else !%%%%
       call   ppush1_b(mc,xi,xf)
    endif !%%%%
  end subroutine ppush1

  !  used for slow pushing
  subroutine trx_cpp(h,xi,xf)  !,rh,y)
    implicit none
    real(dp) xf
    real(dp),dimension(:)::xi
    integer,dimension(1)::rh,hs
    integer i,h,rh1
    integer,dimension(lnv)::y
    integer, allocatable :: jv(:)

    if(.not.c_%stable_da) return
    hs(1)=h
    allocate(jv(c_%nv))
    jv=0
    call ad_alloc(rh1)
    rh(1)=rh1


    do i=1,c_%nv
       call ad_alloc(y(i))
       call ad_const(y(i), xi(i))
    enddo

    call dacct(hs,1,y,c_%nv,rh,1)
    rh1=rh(1)

    call ad_pek(rh1, jv, c_%nv, xf)


    do i=1,c_%nv
       call ad_free(y(i))
    enddo
    call ad_free(rh1)
    deallocate(jv)

    return
  end subroutine trx_cpp

  subroutine dainv(ma,ia,mc,ib)
    implicit none
    integer,dimension(:)::ma,mc
    integer,allocatable::jj(:),ms(:),ml(:),mb(:)
    real(dp),allocatable ::aa(:,:), ai(:,:)
    real(dp),allocatable::ac(:)
    integer i,j,k,ia,ib,ier
    real(dp) amsjj

    if(lingyun_yang) then !%%%%
       allocate(jj(ib),ms(ib),ml(ib))
       allocate(aa(ib,ib), ai(ib,ib))
       allocate(ac(ib))
       allocate(mb(ib))
       jj=0
       do i=1,ib
          call ad_alloc(mb(i))
          call ad_alloc(ms(i))
          call ad_alloc(ml(i))
       enddo

       do i=1,ia
          call dapek(ma(i), jj, ac(i))
          call dapok(ma(i), jj, zero)
       enddo

       do i=1,ia
          call dacon(mb(i),zero)
          do j=1,ib
             jj(j)=1
             call dapek(ma(i), jj, aa(i,j))
             call dapok(ma(i), jj, zero)
             jj(j)=0
          enddo
          call dacmu(ma(i), -one, ma(i))
       enddo


       call matinv(aa,ai,ia,ia,ier)
       if (ier.eq.132) then
          write(*,*) "Can not inverse matrix:"
          write(*,*) aa
          STOP
       endif
       do i=1,ib
          do j=1,ib
             do k=1,ib
                jj(k)=0
             enddo
             jj(j)=1
             call dapok(mb(i),jj,ai(i,j))
             call dapok(ml(i),jj,ai(i,j))
          enddo
       enddo
       !
       !
       do i=2,c_%no
          call dacct(ma,ia,mb,ib,ms,ia)
          do j=1,ib
             do k=1,ib
                jj(k)=0
             enddo
             jj(j)=1
             call dapek(ms(j),jj,amsjj)
             call dapok(ms(j),jj,amsjj+one)
          enddo
          call dacct(ml,ia,ms,ia,mb,ib)
       enddo
       do i=1,ia
          call dacmu(ma(i),-one,ma(i))
          do j=1,ia
             do k=1,ib
                jj(k)=0
             enddo
             call dapok(ma(i),jj,ac(i))
             jj(j)=1
             call dapok(ma(i),jj,aa(i,j))
          enddo
       enddo

       do i=1,ib
          call ad_copy(mb(i), mc(i))
          call ad_free(mb(i))
          call ad_free(ml(i))
          call ad_free(ms(i))
       enddo

       deallocate(mb)
       deallocate(jj,ms,ml)
       deallocate(aa, ai)
       deallocate(ac)
    else !%%%%
       call   dainv_b(ma,ia,mc,ib)
    endif !%%%%

  end subroutine dainv

  subroutine dapin(ma,ia,mb,ib,jx)
    implicit none
    integer ia,ib,i,nv
    integer,dimension(:)::ma,mb,jx
    !
    integer,allocatable :: me(:), mn(:),mi(:),jj(:)

    if(lingyun_yang) then !%%%%
       allocate(me(ib),mn(ib),mi(ib),jj(ib))

       ! me = const = zero, Identity matrix, linear part
       do i=1,ib
          jj(i) = 0
       enddo

       do i=1,ib
          call ad_alloc(me(i))
          call ad_alloc(mn(i))
          call dapok(me(i), jj, zero)
          jj(i)=1
          call dapok(me(i),jj, one)
          jj(i) = 0
       enddo
       ! if jx(i) == 0) mn(i)=me(i)
       ! if(jx(i)/=0) mn(i)=ma(i)
       if (jx(i).eq.0) then
          call ad_copy(me(i), mn(i))
       else
          call ad_copy(ma(i), mn(i))
       endif
       ! dainv(mn,nv,mi,nv)
       call dainv(mn,nv,mi,nv)
       ! if jx(i) == 0) me(i)=ma(i)
       !  dacct(me,nv,mi,nv,mb)
       do i=1,ib
          if (jx(i).eq.0) then
             call ad_copy(ma(i),me(i))
          endif
       enddo
       call dacct(me,ib,mi,ib,mb,ib)

       deallocate(me, mn, mi,jj)
    else !%%%%
       call   dapin_b(ma,ia,mb,ib,jx)
    endif !%%%%

  end subroutine dapin

  subroutine dader(idif,ina,inc)
    implicit none
    integer idif,ina,inc,itmp
    if(lingyun_yang) then !%%%%
       call ad_alloc(itmp)
       call ad_derivative(ina,idif-1,itmp)
       call ad_copy(itmp, inc)
       call ad_free(itmp)
    else !%%%%
       call dader_b(idif,ina,inc)
    endif !%%%%
  end subroutine dader

  subroutine dacfuR(ina,fun,inc)
    implicit none
    integer ina,inc
    integer, allocatable :: j1d(:)
    real(dp), allocatable :: v(:)
    integer ns,nvar,i

    interface
       function fun(abc)
         use precision_constants
         implicit none
         complex(dp) fun
         integer,dimension(:)::abc
       end function fun
    end interface

    if(lingyun_yang) then !%%%%
       call ad_length(ina,ns)

       call ad_nvar(nvar)
       allocate(v(ns),j1d(ns*nvar))

       call ad_read_block(ina, v, j1d, ns)

       do i=1,ns
          if(v(i)/=zero) then
             v(i)=real(fun(j1d((i-1)*nvar+1: i*nvar)),kind=dp)*v(i)
          endif
       enddo

       call ad_save_block(inc, v, j1d, ns)
       deallocate(v,j1d)
    else !%%%%
       call   dacfuR_b(ina,fun,inc)
    endif !%%%%



  end subroutine dacfuR

  subroutine dacfuI(ina,fun,inc)
    implicit none
    integer ina,inc
    integer, allocatable :: j1d(:)
    real(dp), allocatable :: v(:)
    integer ns,nvar,i

    interface
       function fun(abc)
         use precision_constants
         implicit none
         complex(dp) fun
         integer,dimension(:)::abc
       end function fun
    end interface

    if(lingyun_yang) then !%%%%
       call ad_length(ina,ns)

       call ad_nvar(nvar)
       allocate(v(ns),j1d(ns*nvar))

       call ad_read_block(ina, v, j1d, ns)

       do i=1,ns
          if(v(i)/=zero) then
             v(i)=aimag(fun(j1d((i-1)*nvar+1: i*nvar)))*v(i)
          endif
       enddo

       call ad_save_block(inc, v, j1d, ns)
       deallocate(v,j1d)
    else !%%%%
       call dacfuI_b(ina,fun,inc)
    endif !%%%%

  end subroutine dacfuI

  ! done
  subroutine dacfu(ina,fun,inc)
    implicit none
    integer ina,inc
    integer, allocatable :: j1d(:)
    real(dp), allocatable :: v(:)
    integer ns,nvar,i
    !real(dp),external::fun
    interface
       function fun(abc)
         use precision_constants
         implicit none
         real(dp) fun
         integer,dimension(:)::abc
       end function fun
    end interface
    if(lingyun_yang) then !%%%%
       !pause 1
       call ad_length(ina,ns)
       ! write(6,*) ns

       !pause 2
       call ad_nvar(nvar)
       ! write(6,*) nvar
       !pause 3
       allocate(v(ns),j1d(ns*nvar))
       !call ad_print(ina)

       call ad_read_block(ina, v, j1d, ns)

       do i=1,ns
          if(v(i)/=zero) then
             v(i)=fun(j1d((i-1)*nvar+1: i*nvar))*v(i)
          endif
       enddo

       call ad_save_block(inc, v, j1d, ns)
       deallocate(v,j1d)
    else !%%%%
       call dacfu_b(ina,fun,inc)
    endif !%%%%

  end subroutine dacfu

  ! done
  !  subroutine GET_C_J(ina,I,C,J)
  !    implicit none
  !    INTEGER I,nv,ns,k,ina
  !    integer, dimension(lnv)::j
  !    integer, allocatable :: j1d(:)
  !    real(dp), allocatable :: v(:)
  !    real(dp) C
  !    if(lingyun_yang) then !%%%%
  !       call ad_length(ina, ns)
  !       call ad_nvar(nv)
  !       allocate(j1d(nv*ns))
  !       call ad_read_block(ina, v, j1d, ns)
  !       C = v(i)
  !       do k=1,nv
  !          j(k) = j1d((i-1)*nv+k)
  !       enddo
  !       deallocate(j1d)
  !    else !%%%%
  !       call GET_C_J_b(ina,I,C,J)
  !    endif !%%%%
  !  end  subroutine GET_C_J

  subroutine dapri(ina,iunit)
    implicit none
    !    INTEGER ina,iunit
    integer i,ina,ipresent,illa,iunit,iii,ioa,cn
    real(dp) value
    integer, allocatable :: j(:)
    character(10) c10,k10

    if(lingyun_yang) then !%%%%
       ipresent=1
       call dacycle(ina,ipresent,value,illa)

       allocate(j(c_%nv))
       cn=0
       do i=1,illa
          call dacycle(ina,i,value,illa,j)

          if(abs(value)<=eps) then
             cn=cn+1
          endif
       enddo

       if(cn==illa) then
          illa=0
       endif


       j=0
       write(iunit,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)') " Lingyun  ",', NO =',c_%no,', NV =',c_%nv,', INA =',ina,&
            '*********************************************'
       if(illa.ne.0) write(iunit,'(A)') '    I  COEFFICIENT          ORDER   EXPONENTS'
       if(illa.eq.0) write(iunit,'(A)') '   ALL COMPONENTS ZERO '
       c10='      NO ='
       k10='      NV ='
       i=0


       !    write(iunit,'(A10,I6,A10,I6)') c10,c_%no,k10,c_%nv
       cn=0
       do i=1,illa
          call dacycle(ina,i,value,illa,j)

          ioa=0
          do iii=1,c_%nv
             ioa= j(iii)+ioa
          enddo
          if(abs(value)>eps) then
             cn=cn+1
             write(iunit,'(I6,2X,G21.14,I5,4X,18(2i2,1X))') cn,value,ioa,(j(iii),iii=1,c_%nv)
             !ETIENNE
             write(iunit,*) value
          endif
       enddo
       j=0
       write(iunit,'(A)') '                                      '
       deallocate(j)
    else !%%%%
       call dapri_b(ina,iunit)
    endif !%%%%


  end subroutine dapri

  subroutine dapri77(ina,iunit)
    implicit none
    !    INTEGER ina,iunit
    integer i,ina,ipresent,illa,iunit,iii,ioa,cn
    real(dp) value
    integer, allocatable :: j(:)
    character(10) c10,k10

    if(lingyun_yang) then !%%%%
       ipresent=1
       call dacycle(ina,ipresent,value,illa)

       allocate(j(c_%nv))
       cn=0
       do i=1,illa
          call dacycle(ina,i,value,illa,j)

          if(abs(value)<=eps) then
             cn=cn+1
          endif
       enddo

       if(cn==illa) then
          illa=0
       endif


       j=0
       write(iunit,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)') " Lingyun  ",', NO =',c_%no,', NV =',c_%nv,', INA =',ina,&
            '*********************************************'
       if(illa.ne.0) write(iunit,'(A)') '    I  COEFFICIENT          ORDER   EXPONENTS'
       if(illa.eq.0) write(iunit,'(A)') '   ALL COMPONENTS ZERO '
       c10='      NO ='
       k10='      NV ='
       i=0


       write(iunit,'(A10,I6,A10,I6)') c10,c_%no,k10,c_%nv
       cn=0
       do i=1,illa
          call dacycle(ina,i,value,illa,j)

          ioa=0
          do iii=1,c_%nv
             ioa= j(iii)+ioa
          enddo
          if(abs(value)>eps) then
             cn=cn+1
             write(iunit,501) ioa,value,(j(iii),iii=1,c_%nv)
          endif
       enddo
501    format(' ', i3,1x,g23.16,1x,100(1x,i2))
502    format(' ', i5,1x,g23.16,1x,100(1x,i2))
       j=0
       write(iunit,502) -cn,zero,(j(iii),iii=1,c_%nv)
       deallocate(j)
    else !%%%%
       call dapri77_b(ina,iunit)
    endif !%%%%

  end subroutine dapri77

  subroutine dashift(ina,inc,ishift)
    implicit none
    integer ina,inc,ishift,itmp
    real(dp) eps
    if(lingyun_yang) then !%%%%
       call ad_alloc(itmp)
       eps = 1e-20
       call ad_shift(ina, ishift, itmp, eps)
       call ad_copy(itmp, inc)
       call ad_free(itmp)
    else !%%%%
       call dashift_b(ina,inc,ishift)
    endif !%%%%
  end  subroutine dashift

  subroutine darea(ina,iunit)
    implicit none
    integer ina,iunit
    integer ii,iin,i,nno,io,iwarin,iwarnv,iwarno,io1
    character(10) c10
    integer,dimension(lnv)::j
    real(dp) c

    if(lingyun_yang) then !%%%%
       read(iunit,'(A10)') c10
       read(iunit,'(18X,I4)') nno
       read(iunit,'(A10)') c10
       read(iunit,'(A10)') c10
       read(iunit,'(A10)') c10

       iin = 0
       !
10     continue
       iin = iin + 1
       read(iunit,'(I6,2X,G21.14,I5,4X,18(2i2,1X))') ii,c,io,(j(i),i=1,c_%nv)
       if(ii.eq.0) goto 20
       !ETIENNE
       read(iunit,*) c
       !ETIENNE
       if(ii.ne.iin) then
          iwarin = 1
       endif
       io1 = 0
       do i=1,c_%nv
          io1 = io1 + j(i)
       enddo
       !
       if(io1.ne.io) then
          iwarnv = 1
          goto 10
       endif
       if(io.gt.c_%no) then
          !        IF(IWARNO.EQ.0) PRINT*,'WARNING IN DAREA, FILE ',
          !    *              'CONTAINS HIGHER ORDERS THAN VECTOR '
          iwarno = 1
          goto 10
       endif
       !
       call ad_pok(ina, j, c_%nv, c)
       goto 10
       !
20     continue

       !
       return
    else !%%%%
       call darea_b(ina,iunit)
    endif !%%%%


  end subroutine darea



  subroutine darea77(ina,iunit)
    implicit none
    integer ina,iunit
    integer nojoh,nvjoh,ii,iche,k,i
    character(10) c10,k10
    integer,dimension(lnv)::j
    real(dp) c

    if(lingyun_yang) then !%%%%
       read(iunit,'(A10)') c10
       read(iunit,'(A10)') c10
       read(iunit,'(A10)') c10
       read(iunit,'(A10)') c10
       read(iunit,'(A10)') c10
       read(iunit,'(A10,I6,A10,I6)') c10,nojoh,k10,nvjoh

10     continue
       read(iunit,*) ii,c,(j(k),k=1,nvjoh)
       if(ii.lt.0) goto 20

       do i=c_%nv+1,nvjoh
          if(j(i).ne.0) goto 10
       enddo
       iche=0
       do i=1,c_%nv
          iche=iche+j(i)
       enddo
       if(iche.gt.c_%no) goto 10

       call ad_pok(ina, j, c_%nv, c)

       goto 10

20     continue

    else !%%%%
       call darea77_b(ina,iunit)
    endif !%%%%

  end subroutine darea77

  ! done
  !  subroutine dainf(inc,inoc,invc,ipoc,ilmc,illc)
  !    implicit none
  !    integer inc,inoc,invc,ipoc,ilmc,illc
  !    if(lingyun_yang) then !%%%%
  !      write(6,*) " FILL_UNI of Sagan only in Berz"
  !      STOP 666
  !    else !%%%%
  !     call dainf_b(inc,inoc,invc,ipoc,ilmc,illc)
  !    endif !%%%%
  !  end subroutine dainf

  subroutine datra(idif,ina,inc)
    implicit none
    integer idif,ina,inc,itmp
    if(lingyun_yang) then !%%%%
       call ad_alloc(itmp)
       call ad_tra(ina,idif-1,itmp)
       call ad_copy(itmp, inc)
       call ad_free(itmp)
    else !%%%%
       call datra_b(idif,ina,inc)
    endif !%%%%
  end subroutine datra

  ! done
  subroutine daran(ina,cm,xran)
    implicit none
    integer ina
    real(dp) cm,xran

    if(lingyun_yang) then !%%%%
       call ad_fill_ran(ina,cm,xran)
    else !%%%%
       call daran_b(ina,cm,xran)
    endif !%%%%

  end subroutine daran

  subroutine dacycle(ina,ipresent,value,illa,j)
    implicit none
    integer illa,ina,ipresent
    integer,optional,dimension(:)::j
    real(dp) value

    if(lingyun_yang) then !%%%%
       call ad_length(ina, illa)
       if(.not.present(j)) return

       call ad_elem(ina, ipresent, j, value)
    else !%%%%
       call dacycle_b(ina,ipresent,value,illa,j)
    endif !%%%%
  end subroutine dacycle

  subroutine dawritefile(ina)
    implicit none
    integer ina,n, nvar,i
    integer, allocatable :: j(:)
    real(dp),allocatable :: v(:)
    ! get real length of ina, without zeros in the tail
    call ad_length(ina, n)
    ! how many
    call ad_nvar(nvar)
    allocate(j(n*nvar),v(n))
    call ad_read_block(ina, v, j, n)
    do i=1,n
       write(6,'(1x,e15.8,5(1x,i4))') v(i),j((i-1)*nvar+1:i*nvar)
    enddo

    ! for compare
    call ad_print(ina)

    deallocate(j,v)
  end subroutine dawritefile

  subroutine dareadfile(ina)
    implicit none
    integer ina
  end subroutine dareadfile

end module dabnew !$$$$
