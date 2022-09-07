subroutine soddin(ierr)
! Copied from Frank Schmidt hrr March 2022
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !  The name SODD stands for "Second Order Detuning and Distortion"
  !           ----             -      -     -            -
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !
  !  Programm SODD consists of 3 parts:
  !
  !  1.) Subroutine Detune calculates the detuning function terms in
  !      first and second order in the strength of the multipoles.
  !  2.) Subroutine Distort1 calculates the distortion function and
  !      Hamiltonian terms in first order in the strength of the
  !      multipoles.
  !  3.) Subroutine Distort2 calculates the distortion function and
  !      Hamiltonian terms in second order in the strength of the
  !      multipoles.
  !
  !
  !  Theory:  Analytical Calculation of Smear and Tune Shift  SSC-232
  !  Authors: J. Bengtsson and J. Irwin
  !  Date:    February 1990
  !
  !  Program Author: Frank Schmidt, CERN SL-Group
  !  Date:           November 1998 - January 1999
  !  Maximum Order:  parameter mmul=11 (22-Pole)
  !  Comments to:    Frank.Schmidt@cern.ch
  !
  !---------------------------------------------------------------------
  !
  !  Requirements:
  !  -------------
  !
  !   The program expects a file "fc.34" which contains 8 columns:
  !   - Position of Multipole [m]
  !   - Name of Multipole
  !   - Multipole Type (>0 erect Multipole; <0 skew Multipole);
  !      eg "3" is an erect sextupole, "-5" is a skew decapole
  !   - Single particle strength `a la SIXTRACK
  !   - Horizontal Betafunction
  !   - Vertical Betafunction
  !   - Horizontal Phase-advance
  !   - Vertical Phase-advance
  !
  !   The last line in this file is special: All entries are given at
  !   the end of the accelerator with the name "END" and the artificial
  !   multipole type "100".
  !
  !   This file can be automatically produced with SIXTRACK when the
  !   linear optics parameter are calculated with the "LINE" block
  !   (see web locaction: http://wwwslap.cern.ch/~frs/). It
  !   is also produced directly from MAD with the DOOM program by
  !   invoking the MAD-SIXTRACK convertor using the `special' flag (see
  !   web locaction: http://wwwslap.cern.ch/~hansg/doom/doom_six.html
  !   for details).
  !
  !   The program needs roughly 50Mbytes in particular for the second
  !   order distortion function due to pairs of multipoles up to 11th
  !   order.
  !
  !---------------------------------------------------------------------
  !
  !  Program Manual:
  !  ---------------
  !
  !  a.) Program selection:
  !      Specify the program:
  !      - iprog=1,3,5,7 => Run Program Detune
  !      - iprog=2,3,6,7 => Run Program Distort1
  !      - iprog=4,5,6,7 => Run Program Distort2
  !  b.) Order selection:
  !      Specify low (n1) and high (n2) limit of orders to be studied
  !      erect and skew elements are denoted with positive and negative
  !      values respectively; eg.:
  !      n1=-4, n2=10 => the orders (-4,-3,-2 and 3-10) will be treated
  !  c.) Position Range:
  !      Specify a position range (etl1,etl2) in [m], the multipoles
  !      located between this range will be used for the calculation
  !  d.) Print_out Switch:
  !      Output files meant for spreadsheets, i.e. no comments. Human
  !      readable output is found in fort.6
  !      iu_on = 0   => No printout
  !      iu_on = 1   => Output at the end of the position range
  !      iu_on = 2   => At each Multipole in the position range
  !      In Detail:
  !
  !      Program Detune
  !
  !                  Multipole
  !                  Strength
  !      Unit iu_on  Order          Contents in each column
  !      ----------------------------------------------------------------
  !      70     1      1            multipole order
  !                          (hor.,ver.) plane => (1,2)
  !                          hor. or ver. detuning
  !                          order of horizontal Emittance
  !                          order of vertical Emittance
  !      ----------------------------------------------------------------
  !      71     2      1            multipole order
  !                          (hor.,ver.) plane => (1,2)
  !                          hor. or ver. detuning
  !                          order of horizontal Emittance
  !                          order of vertical Emittance
  !      ----------------------------------------------------------------
  !      72     1      2            first multipole order
  !                          second multipole order
  !                          horizontal detuning
  !                          order of horizontal Emittance
  !                          order of vertical Emittance
  !      ----------------------------------------------------------------
  !      73     1      2      first multipole order
  !                          second multipole order
  !                          vertical detuning
  !                          order of horizontal Emittance
  !                          order of vertical Emittance
  !      ---------------------------------------------------------------
  !
  !      Program Distort1
  !
  !                  Multipole
  !                  Strength
  !      Unit iu_on  Order          Contents in each column
  !      ----------------------------------------------------------------
  !      74     1      1            multipole order
  !                          cosine part of distortion
  !                          sine part of distortion
  !                          amplitude of distortion
  !                          j
  !                          k
  !                          l
  !                          m
  !      ----------------------------------------------------------------
  !      75     1      1            multipole order
  !                          cosine part of Hamiltonian
  !                          sine part of Hamiltonian
  !                          amplitude of Hamiltonian
  !                          j
  !                          k
  !                          l
  !                          m
  !      ----------------------------------------------------------------
  !      76     2      1            multipole order
  !                          appearance number in position range
  !                          number of resonance
  !                          position
  !                          cosine part of Hamiltonian
  !                          sine part of Hamiltonian
  !                          amplitude of Hamiltonian
  !                          j
  !                          k
  !                          l
  !                          m
  !      ----------------------------------------------------------------
  !      77     2      1            multipole order
  !                          appearance number in position range
  !                          number of resonance
  !                          position
  !                          cosine part of distortion
  !                          sine part of distortion
  !                          amplitude of distortion
  !                          j
  !                          k
  !                          l
  !                          m
  !      ----------------------------------------------------------------
  !
  !      Program Distort2
  !
  !                  Multipole
  !                  Strength
  !      Unit iu_on  Order          Contents in each column
  !      ----------------------------------------------------------------
  !      78     1      2            first multipole order
  !                                 second multipole order
  !                                 cosine part of distortion
  !                          sine part of distortion
  !                          amplitude of distortion
  !                          j
  !                          k
  !                          l
  !                          m
  !      ----------------------------------------------------------------
  !      79     1      2            first multipole order
  !                                 second multipole order
  !                                 cosine part of Hamiltonian
  !                          sine part of Hamiltonian
  !                          amplitude of Hamiltonian
  !                          j
  !                          k
  !                          l
  !                          m
  !      ----------------------------------------------------------------
  !
  !---------------------------------------------------------------------
  use sodd
  implicit none

  integer ierr

  integer ier,iprog,n3
  double precision etl20

  !--- szcompar is the size of the arrays
  !--- returned by the routine comm_para
  integer szcompar
  parameter (szcompar = 100)

  !--- szchara is the size of the character strings char_a
  !--- returned by the routine comm_para
  integer szchara
  parameter (szchara = 400)

  integer detun,disto1,disto2,prend,prblup
  integer nint, ndble, k, int_arr(szcompar), char_l(szcompar)
  double precision d_arr(szcompar)
  character * (szchara) char_a

  !---------------------------------------------------------------------

  ierr = 0
  j70 = 0
  j71 = 0
  j72 = 0
  j73 = 0
  j74 = 0
  j75 = 0
  j76 = 0
  j77 = 0
  j78 = 0
  j79 = 0
  table_size_70 = 2
  table_size_71 = 2
  table_size_72 = 2
  table_size_73 = 2
  table_size_74 = 2
  table_size_75 = 2
  table_size_76 = 2
  table_size_77 = 2
  table_size_78 = 2
  table_size_79 = 2

  !--- get detun

  detun = 0
  call comm_para('detune ', nint, ndble, k, int_arr, d_arr,         &
       &char_a, char_l)
  if (nint .gt. 0) detun = int_arr(1)

  !--- get disto1

  disto1 = 0
  call comm_para('distort1 ', nint, ndble, k, int_arr, d_arr,       &
       &char_a, char_l)
  if (nint .gt. 0) disto1 = int_arr(1)

  !--- get disto2

  disto2 = 0
  call comm_para('distort2 ', nint, ndble, k, int_arr, d_arr,       &
       &char_a, char_l)
  if (nint .gt. 0) disto2 = int_arr(1)

  iprog = detun + disto1*2 + disto2*4
  if(iprog.le.0.or.iprog.gt.7) call prror(0,14,iprog)
  if(iprog.eq.1.or.iprog.eq.3.or.iprog.eq.5.or.iprog.eq.7)          &
       &write(6,*) ' Program <DETUNE> will be executed '
  if(iprog.eq.2.or.iprog.eq.3.or.iprog.eq.6.or.iprog.eq.7)          &
       &write(6,*) ' Program <DISTORT1> will be executed '
  if(iprog.eq.4.or.iprog.eq.5.or.iprog.eq.6.or.iprog.eq.7)          &
       &write(6,*) ' Program <DISTORT2> will be executed '


  char_a = ' '
  call comm_para('multipole_order_range ', nint, ndble, k, int_arr, &
       &d_arr,char_a, char_l)
  n1 = int_arr(1)
  n2 = int_arr(2)
  if(n1.gt.n2) then
     n3=n1
     n1=n2
     n2=n3
  endif
  if(abs(n1).gt.mmul.or.abs(n2).gt.mmul) call prror(0,15,mmul)

  !--- get etl1 & etl2 (start position & stop position)

  etl1 = 0.
  etl2 = 0.
  call comm_para('start_stop ', nint, ndble, k, int_arr, d_arr,     &
       &char_a, char_l)
  if (ndble .gt. 0) then
     etl1 = d_arr(1)
     etl2 = d_arr(2)
  endif

  print *,"Start and End Position in [m]: ", etl1,etl2

  if(etl1.lt.0) etl1=-etl1
  if(etl2.lt.0) etl2=-etl2
  if(etl2.lt.etl1) then
     etl20=etl2
     etl1=etl2
     etl2=etl20
  endif

  !--- get no print

  iu_on = 0
  call comm_para('noprint ', nint, ndble, k, int_arr, d_arr,        &
       &char_a, char_l)

  if (int_arr(1) .eq. 0) then

     !--- get prend

     prend = 0
     call comm_para('print_at_end ', nint, ndble, k, int_arr, d_arr, &
          &char_a, char_l)
     if (nint .gt. 0) prend = int_arr(1)

     !--- get prblup

     prblup = 0
     call comm_para('print_all ', nint, ndble, k, int_arr, d_arr,    &
          &char_a, char_l)
     if (nint .gt. 0) prblup = int_arr(1)

     iu_on = prend + prblup*2
  endif

  if(iu_on.lt.0.or.iu_on.gt.3) iu_on=0
  if(iu_on.eq.0) write(6,*) ' No print_out requested '
  if(iu_on.eq.1.or.iu_on.eq.3) then
     write(6,*) ' Print_out at end of structure '
     if(detun.eq.1) then
        open(70,file='detune_1_end',form='formatted',status='unknown',&
             &iostat=ier)
        if(ier.ne.0) call prror(0,16,70)
        write(70,10070)
        open(72,file='detune_2_hor',form='formatted',status='unknown',&
             &iostat=ier)
        if(ier.ne.0) call prror(0,16,72)
        write(72,10072)
        open(73,file='detune_2_ver',form='formatted',                 &
             &status='unknown',iostat=ier)
        if(ier.ne.0) call prror(0,16,73)
        write(73,10073)
     endif
     if(disto1.eq.1) then
        open(74,file='distort_1_f_end',form='formatted',              &
             &status='unknown',iostat=ier)
        if(ier.ne.0) call prror(0,16,74)
        write(74,10074)
        open(75,file='distort_1_h_end',form='formatted',              &
             &status='unknown',iostat=ier)
        if(ier.ne.0) call prror(0,16,75)
        write(75,10074)
     endif
     if(disto2.eq.1) then
        open(78,file='distort_2_f_end',form='formatted',              &
             &status='unknown',iostat=ier)
        if(ier.ne.0) call prror(0,16,78)
        write(78,10078)
        open(79,file='distort_2_h_end',form='formatted',              &
             &status='unknown',iostat=ier)
        if(ier.ne.0) call prror(0,16,79)
        write(79,10078)
     endif
  endif
  if(iu_on.eq.2.or.iu_on.eq.3) then
     write(6,*) ' Print_out at each multipole '
     if(detun.eq.1) then
        open(71,file='detune_1_all',form='formatted',                 &
             &status='unknown',iostat=ier)
        if(ier.ne.0) call prror(0,16,71)
        write(71,10070)
     endif
     if(disto1.eq.1) then
        open(76,file='distort_1_f_all',form='formatted',              &
             &status='unknown',iostat=ier)
        if(ier.ne.0) call prror(0,16,76)
        write(76,10076)
        open(77,file='distort_1_h_all',form='formatted',              &
             &status='unknown',iostat=ier)
        if(ier.ne.0) call prror(0,16,77)
        write(77,10076)
     endif
  endif
  open(34,file='fc.34',form='formatted',status='unknown',           &
       &recl=8192,iostat=ier)
  if(ier.ne.0) call prror(0,16,34)
  !-----------------------------------------------------------------------
  !  Initialisation
  !-----------------------------------------------------------------------
  call init
  call timest(tlim)
  call timex(tim1)
  !-----------------------------------------------------------------------
  !  Read Data
  !-----------------------------------------------------------------------
  call readdat
  !-----------------------------------------------------------------------

  call timex(tim2)
  write(6,10000) tim2-tim1
  if(iprog.eq.1.or.iprog.eq.3.or.iprog.eq.5.or.iprog.eq.7)          &
       &call detune
  call timex(tim2)
  if(iprog.eq.2.or.iprog.eq.3.or.iprog.eq.6.or.iprog.eq.7)          &
       &call distort1
  call timex(tim2)
  if(iprog.eq.4.or.iprog.eq.5.or.iprog.eq.6.or.iprog.eq.7)          &
       &call distort2
  call timex(tim6)
  write(6,10010) tim6-tim1

  continue
  if(detun.eq.1) then
     close(70,status='keep')
     close(71,status='keep')
     close(72,status='keep')
     close(73,status='keep')
  endif
  if(disto1.eq.1) then
     close(74,status='keep')
     close(75,status='keep')
     close(76,status='keep')
     close(77,status='keep')
  endif
  if(disto2.eq.1) then
     close(78,status='keep')
     close(79,status='keep')
  endif

10000 format(/80('-')//'Reading Data took ',f10.3,' second(s)',         &
       &' of Execution Time'//80('-')//)
10010 format(/80('-')//'Program <SODD> used a total of ',f10.3,         &
       &' second(s) of Execution Time'//80('-')//)
10070 format('MulOrd ',' Plane',4x,'Hor/Vert detuning',3x,'Hinv Vinv')
10072 format('MulOrd1 ','MulOrd2 ',3x,'Horizontal detuning',            &
       &2x,'Hinv Vinv')
10073 format('MulOrd1 ','MulOrd2 ',3x,'Vertical detuning  ',            &
       &2x,'Hinv Vinv')
10074 format('MulOrd',11x,'Cosine',17x,'Sine',19x,'Amplitude',          &
       &8x,'j',3x,'k',3x,'l',3x,'m',1x)
10076 format('MulOrd',3x,'Loc',2x,'Res',6x,'Position[m]',14x,'Cosine',  &
       &17x,'Sine',19x,'Amplitude',                                       &
       &8x,'j',3x,'k',3x,'l',3x,'m',1x)
10078 format('MulOrd1',' MulOrd2',10x,'Cosine',17x,'Sine',19x,          &
       &'Amplitude',8x,'j',3x,'k',3x,'l',3x,'m',1x)
  return
end subroutine soddin

subroutine detune
  !-----------------------------------------------------------------------
  !--Program Detune
  !-----------------------------------------------------------------------
  use sodd
  implicit none
  integer i,ia,iaa,iah,iamin,ib,ic,id,ih1,ii,imu,imuty,j,jj,jstep,k,&
       &k1,k2,k3,l,l0,l1,l2,l3,m
  double precision betdx,betdxi,fac,facd0,facd1,phix,phiy,          &
       &sigde,sinar,tij
  dimension tij(0:mmul2,-mmul2:mmul2)
  dimension betdxi(mmultx),betdx(mmultx)
  dimension imu(2)
  !---------------------------------------------------------------------
  iamin=2
  jstep=2
  !-----------------------------------------------------------------------
  !  Read Data
  !-----------------------------------------------------------------------
  do i=0,mmul2
     do j=-mmul2,mmul2
        tij(i,j)=zero
     enddo
  enddo
  !---------------------------------------------------------------------
  !--
  !--Detuning First Order in Multipole Strength
  !--
  !---------------------------------------------------------------------
  do ia=2,mmul,jstep
     imu(1)=ia
     imu(2)=0
     iah=ia/2
     if(ityc(ia).gt.0) then
        write(6,10000) 'Structure has: ',ityc(ia),' single ',         &
             &comment(ia)
        ih1=ia/2-1
        do i=1,ityc(ia)
           fac=-two/dble(ia)*bstr(ia,i)
           sigde=-one
           do j=1,iah
              jj=j-1
              sigde=-sigde
              det1(1,j)=sigde*fac*cosav(jj)*cosav(iah-jj)*dble(iah-jj)* &
                   &fact(ia,1+jj*2)*beta(1,ia,i)**(iah-jj)*beta(2,ia,i)**jj
              det(1,ia,iah-j,jj)=det(1,ia,iah-j,jj)+det1(1,j)
              if(j.eq.iah) then
                 sigde=-sigde
                 det1(2,iah)=sigde*fac*dble(iah)*cosav(iah)*             &
                      &beta(2,ia,i)**iah
                 det(2,ia,0,iah-1)=det(2,ia,0,iah-1)+det1(2,iah)
              endif
           enddo
           if(iu_on.gt.0) call detwri(0,imu,ih1)
        enddo
        !---------------------------------------------------------------------
        !--
        !--Printing First Order Terms
        !--
        !---------------------------------------------------------------------
        call detwri(1,imu,ih1)
     endif
  enddo
  call timex(tim3)
  write(6,10030) tim3-tim2
  !---------------------------------------------------------------------
  !--
  !--Calculation of Detuning Poisson Brackets
  !--
  !-----Needed for Second Order Terms
  !--
  !---------------------------------------------------------------------
  do imuty=1,2
     call timex(tim3)
     !---------------------------------------------------------------------
     !--
     !--Detuning Second Order in Multipole Strength
     !--
     !---------------------------------------------------------------------
     do ia=iamin,mmul
        do ib=ia,mmul,jstep
           if(imuty.eq.1) ic=ia
           if(imuty.eq.1) id=ib
           if(imuty.eq.2) ic=-ia
           if(imuty.eq.2) id=-ib
           if(ityc(ic).gt.0.and.ityc(id).gt.0) then
              imu(1)=ic
              imu(2)=id
              do k=1,mmultx
                 ihv(k)=0
                 iplane(k,1)=0
                 iplane(k,2)=0
                 betexp(k,1)=zero
                 betexp(k,2)=zero
                 betexp(k,3)=zero
                 betexp(k,4)=zero
                 fac2(k)=zero
                 do m=1,mmult
                    itij(k,m,1)=0
                    itij(k,m,2)=0
                    itij(k,m,3)=0
                    fac4(k,m)=zero
                 enddo
              enddo
              call caldet2(ia,ib,imuty)
              do ii=1,ityc(ic)
                 do k=1,icc
                    betdxi(k)=(beta(1,ic,ii)**betexp(k,1))*               &
                         &(beta(2,ic,ii)**betexp(k,2))
                 enddo
                 do jj=1,ityc(id)
                    if(ii.eq.1.and.jj.eq.1) then
                       if(ia.ne.ib) then
                          write(6,10050) 'Pairs of: ',ityc(ic),comment(ic), &
                               &' and: ',ityc(id),comment(id)
                       else
                          write(6,10060) 'Quadratic Effect of: ',           &
                               &ityc(ic),comment(ic)
                       endif
                    endif
                    !--Phase Differences
                    phix=abs(phi(1,id,jj)-phi(1,ic,ii))-qx
                    phiy=abs(phi(2,id,jj)-phi(2,ic,ii))-qy
                    !--Phase Factor
                    if(imuty.eq.1) iaa=0
                    if(imuty.eq.2) iaa=1
                    do i=ia+iaa,0,-2
                       do j=i-ia,ia-i,2
                          if(i.ne.0.or.j.gt.0) then
                             sinar=sin(dble(i)*qx+dble(j)*qy)
                             if(abs(sinar).gt.pieni) then
                                tij(i,j)=cos(dble(i)*phix+dble(j)*phiy)/sinar
                             else
                                write(6,*) '  Warning: Results for detuning ',&
                                     &'probably wrong; sinar= ',sinar
                                tij(i,j)=zero
                             endif
                          endif
                       enddo
                    enddo
                    !--Detuning Calculation
                    facd0=bstr(ic,ii)*bstr(id,jj)
                    do k=1,icc
                       k1=iplane(k,1)
                       k2=iplane(k,2)
                       k3=ihv(k)
                       if(k3.eq.1.or.(k3.eq.2.and.k1.eq.0)) then
                          betdx(k)=betdxi(k)*                               &
                               &(beta(1,id,jj)**betexp(k,3))*                                     &
                               &(beta(2,id,jj)**betexp(k,4))
                          facd1=zero
                          l0=icd(k)
                          do l=1,l0
                             l1=itij(k,l,1)
                             l2=itij(k,l,2)
                             l3=itij(k,l,3)
                             facd1=facd1+fac4(k,l)*(tij(l2,l3)+              &
                                  &dble(l1)*tij(l2,-l3))
                          enddo
                          det(k3,ic,k1,k2)=det(k3,ic,k1,k2)+facd0*          &
                               &fac2(k)*facd1*betdx(k)
                       endif
                    enddo
                    do i=ia+iaa,0,-2
                       do j=i-ia,ia-i,2
                          if(i.ne.0.or.j.gt.0) then
                             tij(i,j)=zero
                          endif
                       enddo
                    enddo
                 enddo
              enddo
              !---------------------------------------------------------------------
              !--
              !--Printing Second Order Terms
              !--
              !---------------------------------------------------------------------
              ih1=(ia+ib)/2-2
              call detwri(2,imu,ih1)
           endif
        enddo
     enddo
     call timex(tim4)
     write(6,10090) tim4-tim3
  enddo
  call timex(tim5)
  write(6,10100) tim5-tim2
  !---------------------------------------------------------------------
10000 format(//80('-')//5x,a,i6,2a)
10030 format(/'First Order Detuning Calculation took ',                 &
       &f10.3,' second(s)',' of Execution Time'//80('-')//)
10050 format(//80('-')//10x,a,i6,2a,i6,a)
10060 format(//80('-')//10x,a,i6,a)
10090 format(/80('-')/'Second Order Detuning Calculation took ',        &
       &f10.3,' second(s)',' of Execution Time'//80('-')//)
10100 format(/80('-')//'Program <DETUNE> used ',f10.3,                  &
       &' second(s) of Execution Time'//80('-')//)
  return
end subroutine detune
subroutine distort1
  !-----------------------------------------------------------------------
  !--Program Distort1 (first order)
  !-----------------------------------------------------------------------
  use sodd
  implicit none
  integer i,ia,iadd,ib,ic,id,id1,id11,id110,id12,id120,idiff1,      &
       &idiff2,ie,if,ii,iia,iia2,ivar,j
  double precision distc,dists,dstc,dsts,facc,facs,fact1,factb1,    &
       &factb2,fpre,hamc,hams,hmc,hms,signii,sinar,tc,ts,amp
  integer data_size
  character*16 table_name
  integer int_to_write(11)
  double precision double_to_write(11)

  dimension distc(-mmul:mmul,mh),dists(-mmul:mmul,mh)
  dimension dstc(mh),dsts(mh),ivar(-mmul:mmul),fact1(mh)
  dimension hamc(-mmul:mmul,mh),hams(-mmul:mmul,mh)
  !---------------------------------------------------------------------
  do i=1,mh
     fact1(i)=zero
     do j=-mmul,mmul
        ivar(j)=0
        ifacta(1,j,i)=0
        ifacta(2,j,i)=0
        ifact4(1,j,i)=0
        ifact4(2,j,i)=0
        ifact4(3,j,i)=0
        ifact4(4,j,i)=0
        factb(1,j,i)=zero
        factb(2,j,i)=zero
        fact0(j,i)=zero
        distc(j,i)=zero
        dists(j,i)=zero
        hamc(j,i)=zero
        hams(j,i)=zero
     enddo
  enddo
  iadd=2
  idiff1=1
  idiff2=2
  ivar(1)=1
  ivar(-1)=1
  !--Number of Resonances
  do i=2,mmul
     ivar(i)=iadd
     ivar(-i)=iadd
     if(mod(i,2).eq.0) then
        idiff1=idiff1+idiff2
        idiff2=idiff2+1
     endif
     iadd=iadd+idiff1
  enddo
  !--Calculation of the various Factors for Multipoles -mmul -> +mmul
  do ia=-mmul,mmul
     if(ityc(ia).gt.0) then
        iia=iabs(ia)
        fpre=one/dble((iia*2**(iia+1)))
        ib=0
        signii=-one
        if(ia.lt.0) then
           iia2=iia-1
        else
           iia2=iia
        endif
        do ic=iia2,0,-2
           signii=signii*(-one)
           ib=ib+1
           id=iia-ic
           id1=id+1
           ifacta(1,ia,ib)=ic
           ifacta(2,ia,ib)=id
           ifact4(1,ia,ib)=ic
           ifact4(2,ia,ib)=0
           ifact4(3,ia,ib)=id
           ifact4(4,ia,ib)=0
           factb1=abs(ifacta(1,ia,ib))/two
           factb2=abs(ifacta(2,ia,ib))/two
           factb(1,ia,ib)=factb1
           factb(2,ia,ib)=factb2
           fact0(ia,ib)=signii*fpre*fact(iia,id1)
           if(id.gt.0.and.ic.gt.0) then
              ib=ib+1
              ifacta(1,ia,ib)=ic
              ifacta(2,ia,ib)=-id
              ifact4(1,ia,ib)=ic
              ifact4(2,ia,ib)=0
              ifact4(3,ia,ib)=0
              ifact4(4,ia,ib)=id
              factb(1,ia,ib)=factb1
              factb(2,ia,ib)=factb2
              fact0(ia,ib)=signii*fpre*fact(iia,id1)
           endif
           id11=1
           id110=0
           do ie=ic,0,-2
              if(id.ne.0.or.ie.ne.0) then
                 if(ie.ne.ic) then
                    id110=id110+1
                    ib=ib+1
                    id11=id11+1
                    ifacta(1,ia,ib)=ie
                    ifacta(2,ia,ib)=id
                    ifact4(1,ia,ib)=ie+id110
                    ifact4(2,ia,ib)=id110
                    ifact4(3,ia,ib)=id
                    ifact4(4,ia,ib)=0
                    factb(1,ia,ib)=factb1
                    factb(2,ia,ib)=factb2
                    fact0(ia,ib)=signii*fpre*fact(iia,id1)*               &
                         &fact(ic,id11)
                    if(id.gt.0.and.ie.gt.0) then
                       ib=ib+1
                       ifacta(1,ia,ib)=ie
                       ifacta(2,ia,ib)=-id
                       ifact4(1,ia,ib)=ie+id110
                       ifact4(2,ia,ib)=id110
                       ifact4(3,ia,ib)=0
                       ifact4(4,ia,ib)=id
                       factb(1,ia,ib)=factb1
                       factb(2,ia,ib)=factb2
                       fact0(ia,ib)=signii*fpre*fact(iia,id1)*             &
                            &fact(ic,id11)
                    endif
                 endif
                 id12=1
                 id120=0
                 do if=id-2,0,-2
                    if(if.ne.0.or.ie.ne.0) then
                       id120=id120+1
                       ib=ib+1
                       id12=id12+1
                       ifacta(1,ia,ib)=ie
                       ifacta(2,ia,ib)=if
                       ifact4(1,ia,ib)=ie+id110
                       ifact4(2,ia,ib)=id110
                       ifact4(3,ia,ib)=if+id120
                       ifact4(4,ia,ib)=id120
                       factb(1,ia,ib)=factb1
                       factb(2,ia,ib)=factb2
                       fact0(ia,ib)=signii*fpre*fact(iia,id1)*             &
                            &fact(ic,id11)*fact(id,id12)
                       if(if.gt.0.and.ie.gt.0) then
                          ib=ib+1
                          ifacta(1,ia,ib)=ie
                          ifacta(2,ia,ib)=-if
                          ifact4(1,ia,ib)=ie+id110
                          ifact4(2,ia,ib)=id110
                          ifact4(3,ia,ib)=id120
                          ifact4(4,ia,ib)=if+id120
                          factb(1,ia,ib)=factb1
                          factb(2,ia,ib)=factb2
                          fact0(ia,ib)=signii*fpre*fact(iia,id1)*           &
                               &fact(ic,id11)*fact(id,id12)
                       endif
                    endif
                 enddo
                 id120=0
              endif
           enddo
           id110=0
        enddo
        call sortres(2,ia,ib)
     endif
  enddo
  !--Calculation of I Order Distortion Function
  do iia=2,mmul
     do ii=1,2
        if(ii.eq.1) ia=iia
        if(ii.eq.2) ia=-iia
        if(ityc(ia).gt.0) then
           write(6,10100) 'First Order Contribution of: ',             &
                &ityc(ia),comment(ia)
           do i=1,ityc(ia)
              do j=1,ivar(ia)
                 !--Factor of Machine Parameters
                 fact1(j)=bstr(ia,i)*fact0(ia,j)*                        &
                      &beta(1,ia,i)**factb(1,ia,j)*beta(2,ia,i)**factb(2,ia,j)
                 !--Phase Factor
                 tc=sin(ifacta(1,ia,j)*(phi(1,ia,i)-qx)+                 &
                      &ifacta(2,ia,j)*(phi(2,ia,i)-qy))
                 ts=cos(ifacta(1,ia,j)*(phi(1,ia,i)-qx)+                 &
                      &ifacta(2,ia,j)*(phi(2,ia,i)-qy))
                 !--Generating Function Calculation
                 sinar=sin(ifacta(1,ia,j)*qx+ifacta(2,ia,j)*qy)
                 if(abs(sinar).gt.pieni) then
                    dstc(j)=fact1(j)*tc/sinar
                    dsts(j)=fact1(j)*ts/sinar
                 else
                    write(6,*) '  Warning: Results for distortion ',      &
                         &'function probably wrong; sinar= ',sinar
                    dstc(j)=zero
                    dsts(j)=zero
                 endif
                 distc(ia,j)=distc(ia,j)+dstc(j)
                 dists(ia,j)=dists(ia,j)+dsts(j)
              enddo
              !--Write Generating Function Calculation of each Element
              if(iu_on.eq.2.or.iu_on.eq.3) then
                 table_name ='distort_1_f_all '
                 data_size = ivar(ia)
                 call fit_table(table_name,data_size,table_size_76)
                 do j=1,data_size
                    amp = dsqrt(dstc(j)*dstc(j)+dsts(j)*dsts(j))
                    write(76,10000) ia,i,j,etl(ia,i),dstc(j),dsts(j),amp, &
                         &ifact4(1,ia,j),ifact4(2,ia,j),                                    &
                         &ifact4(3,ia,j),ifact4(4,ia,j)
                    data_size = ivar(ia)
                    int_to_write(1) = ia
                    int_to_write(2) = i
                    int_to_write(3) = j
                    double_to_write(4) = etl(ia,i)
                    double_to_write(5) = dstc(j)
                    double_to_write(6) = dsts(j)
                    double_to_write(7) = amp
                    int_to_write(8) = ifact4(1,ia,j)
                    int_to_write(9) = ifact4(2,ia,j)
                    int_to_write(10) = ifact4(3,ia,j)
                    int_to_write(11) = ifact4(4,ia,j)
                    call write_table(table_name,11,int_to_write,          &
                         &double_to_write)
                    j76 = j76 + 1
                    call augment_count(table_name)
                 enddo
              endif
              !--Calculate and write Hamiltonian of each Element
              if(iu_on.eq.2.or.iu_on.eq.3) then
                 table_name ='distort_1_h_all '
                 data_size = ivar(ia)
                 call fit_table(table_name,data_size,table_size_77)
                 do j=1,data_size
                    facc=one-cos(two*(ifacta(1,ia,j)*qx+                  &
                         &ifacta(2,ia,j)*qy))
                    facs=sin(two*(ifacta(1,ia,j)*qx+ifacta(2,ia,j)*qy))
                    hmc=facc*dstc(j)+facs*dsts(j)
                    hms=facc*dsts(j)-facs*dstc(j)
                    double_to_write(5)=hmc
                    double_to_write(6)=hms
                    amp = dsqrt(hmc*hmc+hms*hms)
                    write(77,10000) ia,i,j,etl(ia,i),double_to_write(5),  &
                         &double_to_write(6),amp,ifact4(1,ia,j),ifact4(2,ia,j),             &
                         &ifact4(3,ia,j),ifact4(4,ia,j)
                    int_to_write(1) = ia
                    int_to_write(2) = i
                    int_to_write(3) = j
                    double_to_write(4) = etl(ia,i)
                    double_to_write(7) = amp
                    int_to_write(8) = ifact4(1,ia,j)
                    int_to_write(9) = ifact4(2,ia,j)
                    int_to_write(10) = ifact4(3,ia,j)
                    int_to_write(11) = ifact4(4,ia,j)
                    call write_table(table_name,11,int_to_write,          &
                         &double_to_write)
                    j77 = j77 + 1
                    call augment_count(table_name)
                 enddo
              endif
           enddo
           !--Write total Generating Function
           if(iu_on.eq.1.or.iu_on.eq.3) then
              write(6,10020)
              table_name ='distort_1_f_end '
              data_size = ivar(ia)
              call fit_table(table_name,data_size,table_size_74)
              do j=1,data_size
                 amp = sqrt(distc(ia,j)*distc(ia,j)+dists(ia,j)*         &
                      &dists(ia,j))
                 write(6,10040) distc(ia,j),dists(ia,j),amp,             &
                      &ifact4(1,ia,j),ifact4(2,ia,j),ifact4(3,ia,j),                     &
                      &ifact4(4,ia,j)
                 write(74,10010) ia,distc(ia,j),dists(ia,j),amp,         &
                      &ifact4(1,ia,j),ifact4(2,ia,j),ifact4(3,ia,j),                     &
                      &ifact4(4,ia,j)
                 int_to_write(1) = ia
                 double_to_write(2) = distc(ia,j)
                 double_to_write(3) = dists(ia,j)
                 double_to_write(4) = amp
                 int_to_write(5) = ifact4(1,ia,j)
                 int_to_write(6) = ifact4(2,ia,j)
                 int_to_write(7) = ifact4(3,ia,j)
                 int_to_write(8) = ifact4(4,ia,j)
                 call write_table(table_name,8,int_to_write,             &
                      &double_to_write)
                 j74 = j74 + 1
                 call augment_count(table_name)
              enddo
              !--Calculate and write total Hamiltonian
              write(6,10030)
              table_name ='distort_1_h_end '
              data_size = ivar(ia)
              call fit_table(table_name,data_size,table_size_75)
              do j=1,data_size
                 facc=one-cos(two*(ifacta(1,ia,j)*qx+ifacta(2,ia,j)*qy))
                 facs=sin(two*(ifacta(1,ia,j)*qx+ifacta(2,ia,j)*qy))
                 hamc(ia,j)=facc*distc(ia,j)+facs*dists(ia,j)
                 hams(ia,j)=facc*dists(ia,j)-facs*distc(ia,j)
                 double_to_write(2)=hamc(ia,j)
                 double_to_write(3)=hams(ia,j)
                 amp = sqrt(hamc(ia,j)*hamc(ia,j)+hams(ia,j)*hams(ia,j))
                 write(6,10040) double_to_write(2),double_to_write(3),   &
                      &amp,ifact4(1,ia,j),ifact4(2,ia,j),                                &
                      &ifact4(3,ia,j),ifact4(4,ia,j)
                 write(75,10010) ia,double_to_write(2),                  &
                      &double_to_write(3),amp,ifact4(1,ia,j),ifact4(2,ia,j),             &
                      &ifact4(3,ia,j),ifact4(4,ia,j)
                 int_to_write(1) = ia
                 double_to_write(4) = amp
                 int_to_write(5) = ifact4(1,ia,j)
                 int_to_write(6) = ifact4(2,ia,j)
                 int_to_write(7) = ifact4(3,ia,j)
                 int_to_write(8) = ifact4(4,ia,j)
                 call write_table(table_name,8,int_to_write,             &
                      &double_to_write)
                 j75 = j75 + 1
                 call augment_count(table_name)
              enddo
           endif
        endif
     enddo
  enddo
  call timex(tim5)
  write(6,10110) tim5-tim2
  !---------------------------------------------------------------------
10000 format(i4,i7,i5,4(1pe23.15),4i4)
10010 format(i4,4x,3(1pe23.15),6i4)
10020 format(/80('-')/10x,'Distortionfunction Terms (first order)'/     &
       &80('-')/4x,'Cosine-Term',8x,'Sine-Term',8x,'Amplitude',7x,'J',    &
       &5x,'K',5x,'L',5x,'M'/80('-'))
10030 format(/80('-')/10x,'Hamiltonian Terms (first order)'/80('-')/    &
       &4x,'Cosine-Term',8x,'Sine-Term',8x,'Amplitude',7x,'J',            &
       &5x,'K',5x,'L',5x,'M'/80('-'))
10040 format(3(1pe17.9),4i6)
10100 format(//80('-')//10x,a,i8,a)
10110 format(/80('-')//'Program <DISTORT1> used ',f10.3,                &
       &' second(s) of Execution Time'//80('-')//)
  return
end subroutine distort1
subroutine distort2
  !-----------------------------------------------------------------------
  !--Program Distort2 (second order)
  !-----------------------------------------------------------------------
  use sodd
  implicit none
  integer i,ia,ia0,ib,ib0,ic,ic0,ic00,ica,ica2,id,id0,id00,ide,     &
       &ii,imax,imuty,itij0,jj,k,k1,k2,l,l0,l3,l4,ll3,ll4,m,mde,num
  double precision arg,argij,betdx,betdxi,dl12q,dl34q,efact,ephix1, &
       &ephiy1,esin34,etrgij,facc,facd0,facd1,facd10,facp,facs,ffff,phix1,&
       &phix2,phix20,phiy1,phiy2,phiy20,sin34,spij,trgij
  integer data_size
  character*16 table_name
  dimension efact(mmultx,mmult),ll3(mmultx,mmult),ll4(mmultx,mmult)
  dimension betdxi(mmultx),betdx(mmultx)
  dimension facd1(mmultx,mmult,2)
  dimension itij0(mmultx,mmult,4),facd10(mmultx,mmult,2)
  dimension mde(mmul)
  integer int_to_write(11)
  double precision double_to_write(11)
  !---------------------------------------------------------------------
  !--Set the rest of the Variables to zero
  do i=1,mmult
     do ii=1,mmultx
        facd1(ii,i,1)=zero
        facd1(ii,i,2)=zero
     enddo
  enddo
  do i=1,mmul
     mde(i)=0
  enddo
  ii=0
  do i=4,mmul,2
     ii=ii+1
     mde(ii)=i
  enddo
  !---------------------------------------------------------------------
  !--
  !--Calculation of Driving-term Poisson Brackets
  !--
  !-----Needed for Second Order Terms
  !--
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !--
  !--Driving & Hamiltonian Terms Second Order in Multipole Strength
  !--
  !---------------------------------------------------------------------
  do ic00=n1,n2
     ic=ic00
     do id00=-iabs(ic),iabs(ic)
        id=id00
        ia=iabs(ic)
        ib=iabs(id)
        !--Difference between erect and skew
        !--imuty=1 ==> "erect";      imuty=2 ==> "skew";
        !--imuty=3 ==> "erect-skew"; imuty=4 ==> "skew-erect"
        if(ic.gt.0.and.id.gt.0) imuty=1
        if(ic.lt.0.and.id.lt.0) imuty=2
        if(ic.gt.0.and.id.lt.0) imuty=3
        if(ic.lt.0.and.id.gt.0) imuty=4
        call timex(tim3)
        if(ityc(ic).gt.0.and.ityc(id).gt.0) then
           icc=0
           imax=ia+ib-2
           do k=1,mmultx
              betexp(k,1)=zero
              betexp(k,2)=zero
              betexp(k,3)=zero
              betexp(k,4)=zero
              icd(k)=0
              fac2(k)=zero
              do m=1,mmult
                 isig(k,m)=0
                 itij(k,m,1)=0
                 itij(k,m,2)=0
                 itij(k,m,3)=0
                 itij(k,m,4)=0
                 itij(k,m,5)=0
                 itij(k,m,6)=0
                 itij0(k,m,1)=0
                 itij0(k,m,2)=0
                 itij0(k,m,3)=0
                 itij0(k,m,4)=0
                 fac4(k,m)=zero
              enddo
           enddo
           !--
           !--The ide(n0=0,yes=1) switch checks if the average(detuning)terms
           !--create contributions to the distortion function in higher order
           !--A full loop is needed
           !--
           ide=0
           if(imuty.eq.1.and.ia.ne.ib) then
              do k=1,mmul
                 if(ib.eq.mde(k)) ide=1
              enddo
           endif
           if(ide.eq.1) then
              !--
              !--Extra (ide=1) Loop starts here
              !--
              ia0=ia
              ib0=ib
              ic0=ic
              id0=id
              ia=ib0
              ib=ia0
              ic=id0
              id=ic0
              call caldt2(ia,ib,imuty)
              do k=1,mmultx
                 do m=1,mmult
                    itij0(k,m,1)=itij(k,m,3)
                    itij0(k,m,2)=itij(k,m,4)
                    itij0(k,m,3)=itij(k,m,5)
                    itij0(k,m,4)=itij(k,m,6)
                 enddo
              enddo
              do ii=1,ityc(ic)
                 !--Beta term of first set of multipoles
                 do k=1,icc
                    betdxi(k)=(beta(1,ic,ii)**betexp(k,1))*               &
                         &(beta(2,ic,ii)**betexp(k,2))
                 enddo
                 do jj=1,ityc(id)
                    facd0=two*bstr(ic,ii)*bstr(id,jj)
                    !--Phase Differences (phix2,phiy2) of secondary resonances change
                    !--sign for reversed order of the position of multipoles
                    spij=one
                    if(phi(1,id,jj).lt.phi(1,ic,ii)) spij=-one
                    !--Drivingterm Calculation
                    do k=1,icc
                       !--Beta term of second set of multipoles multiplied to first
                       betdx(k)=betdxi(k)*                                 &
                            &(beta(1,id,jj)**betexp(k,3))*                                     &
                            &(beta(2,id,jj)**betexp(k,4))
                       l0=icd(k)
                       !--Cosine and sine part of driving term
                       do l=1,l0
                          if(isig(k,l).eq.1) then
                             dl12q=itij(k,l,1)*qx+itij(k,l,2)*qy
                             if(abs(dl12q).lt.pieni) dl12q=pieni
                             l3=itij(k,l,3)-itij(k,l,4)
                             l4=itij(k,l,5)-itij(k,l,6)
                             dl34q=l3*qx+l4*qy
                             if(abs(dl34q).lt.pieni) dl34q=pieni
                             sin34=sin(dl34q)
                             phix1=-qx
                             phiy1=-qy
                             trgij=sin(dl12q)
                             if(phi(1,id,jj).eq.phi(1,ic,ii)) then
                                phix2=phi(1,id,jj)+qx
                                phiy2=phi(2,id,jj)+qy
                             else
                                phix2=phi(1,id,jj)-spij*qx
                                phiy2=phi(2,id,jj)-spij*qy
                             endif
                             arg=l3*phix1+l4*phiy1+itij(k,l,1)*              &
                                  &phix2+itij(k,l,2)*phiy2
                             facp=fac4(k,l)*facd0*fac2(k)*                   &
                                  &betdx(k)/trgij/sin34
                             facd10(k,l,1)=facd10(k,l,1)+facp*sin(arg)
                             facd10(k,l,2)=facd10(k,l,2)+facp*cos(arg)
                          endif
                       enddo
                    enddo
                 enddo
              enddo
              ia=ia0
              ib=ib0
              ic=ic0
              id=id0
           endif
           !---------------------------------------------------------------------
           !--
           !--Main loop
           !--
           !---------------------------------------------------------------------
           icc=0
           do k=1,mmultx
              betexp(k,1)=zero
              betexp(k,2)=zero
              betexp(k,3)=zero
              betexp(k,4)=zero
              icd(k)=0
              fac2(k)=zero
              do m=1,mmult
                 isig(k,m)=0
                 itij(k,m,1)=0
                 itij(k,m,2)=0
                 itij(k,m,3)=0
                 itij(k,m,4)=0
                 itij(k,m,5)=0
                 itij(k,m,6)=0
                 fac4(k,m)=zero
              enddo
           enddo
           call caldt2(ia,ib,imuty)
           do k=1,icc
              l0=icd(k)
              do l=1,l0
                 dl12q=itij(k,l,1)*qx+itij(k,l,2)*qy
                 if(abs(dl12q).lt.pieni) dl12q=pieni
                 etrgij=sin(dl12q)
                 ll3(k,l)=itij(k,l,3)-itij(k,l,4)
                 ll4(k,l)=itij(k,l,5)-itij(k,l,6)
                 dl34q=ll3(k,l)*qx+ll4(k,l)*qy
                 if(abs(dl34q).lt.pieni) dl34q=pieni
                 esin34=sin(dl34q)
                 efact(k,l)=fac4(k,l)/etrgij/esin34
              enddo
           enddo
           do ii=1,ityc(ic)
              !--Beta term of first set of multipoles
              do k=1,icc
                 betdxi(k)=(beta(1,ic,ii)**betexp(k,1))*                 &
                      &(beta(2,ic,ii)**betexp(k,2))
              enddo
              do jj=1,ityc(id)
                 facd0=bstr(ic,ii)*bstr(id,jj)
                 if(ic.ne.id.and.ic.ne.-id) facd0=facd0*two
                 if(ii.eq.1.and.jj.eq.1) then
                    if(ic.ne.id) then
                       write(6,10090) 'Pairs of: ',ityc(ic),comment(ic),   &
                            &' and: ',ityc(id),comment(id)
                    else
                       write(6,10100) 'Quadratic Effect of: ',             &
                            &ityc(ic),comment(ic)
                    endif
                 endif
                 !--Phase Differences (phix2,phiy2) of secondary resonances change
                 !--sign for reversed order of the position of multipoles
                 spij=one
                 if(phi(1,id,jj).lt.phi(1,ic,ii)) spij=-one
                 phix20=phi(1,id,jj)-phi(1,ic,ii)-spij*qx
                 phiy20=phi(2,id,jj)-phi(2,ic,ii)-spij*qy
                 ephix1=phi(1,ic,ii)-qx
                 ephiy1=phi(2,ic,ii)-qy
                 !--Drivingterm Calculation
                 do k=1,icc
                    !--Beta term of second set of multipoles multiplied to first
                    betdx(k)=betdxi(k)*                                   &
                         &(beta(1,id,jj)**betexp(k,3))*                                     &
                         &(beta(2,id,jj)**betexp(k,4))
                    ffff=facd0*fac2(k)*betdx(k)
                    l0=icd(k)
                    num=0
                    do l=1,l0
                       num=num+isig(k,l)
                    enddo
                    if (num.eq.0) then
                       !--Cosine and sine part of driving term
                       do l=1,l0
                          arg=ll3(k,l)*ephix1+ll4(k,l)*ephiy1+              &
                               &itij(k,l,1)*phix20+itij(k,l,2)*phiy20
                          facp=ffff*efact(k,l)
                          facd1(k,l,1)=facd1(k,l,1)+facp*sin(arg)
                          facd1(k,l,2)=facd1(k,l,2)+facp*cos(arg)
                       enddo
                    else
                       do l=1,l0
                          if(isig(k,l).eq.0) then
                             arg=ll3(k,l)*ephix1+ll4(k,l)*ephiy1+            &
                                  &itij(k,l,1)*phix20+itij(k,l,2)*phiy20
                             facp=ffff*efact(k,l)
                             facd1(k,l,1)=facd1(k,l,1)+facp*sin(arg)
                             facd1(k,l,2)=facd1(k,l,2)+facp*cos(arg)
                          else
                             dl12q=itij(k,l,1)*qx+itij(k,l,2)*qy
                             if(abs(dl12q).lt.pieni) dl12q=pieni
                             l3=itij(k,l,3)-itij(k,l,4)
                             l4=itij(k,l,5)-itij(k,l,6)
                             dl34q=l3*qx+l4*qy
                             if(abs(dl34q).lt.pieni) dl34q=pieni
                             sin34=sin(dl34q)
                             phix1=-qx
                             phiy1=-qy
                             trgij=sin(dl12q)
                             phix2=phi(1,id,jj)
                             phiy2=phi(2,id,jj)
                             if(ic.eq.id.and.phi(1,id,jj).eq.phi(1,ic,ii))   &
                                  &then
                                trgij=trgij/cos(dl12q)
                             else
                                phix2=phix2-spij*qx
                                phiy2=phiy2-spij*qy
                             endif
                             arg=l3*phix1+l4*phiy1+itij(k,l,1)*              &
                                  &phix2+itij(k,l,2)*phiy2
                             facp=fac4(k,l)*facd0*fac2(k)*                   &
                                  &betdx(k)/trgij/sin34
                             !--Factor 2 comes from trigonometric functions
                             if(ic.eq.id) facp=two*facp
                             facd1(k,l,1)=facd1(k,l,1)+facp*sin(arg)
                             facd1(k,l,2)=facd1(k,l,2)+facp*cos(arg)
                          endif
                       enddo
                    endif
                 enddo
              enddo
           enddo
           !--Collect a total of ica non-zero terms
           ica=0
           do k1=1,mmultx
              do k2=1,mmult
                 if(abs(facd1(k1,k2,1)).gt.pieni.or.                     &
                      &abs(facd1(k1,k2,2)).gt.pieni) then
                    ica=ica+1
                    if(ica.gt.mmultf) call prror(3,10,mmultf)
                    facd2(1,ica)=facd1(k1,k2,1)
                    facd2(2,ica)=facd1(k1,k2,2)
                    facd1(k1,k2,1)=zero
                    facd1(k1,k2,2)=zero
                    do ii=1,4
                       ifacd2(ii,ica)=itij(k1,k2,ii+2)
                    enddo
                 endif
              enddo
           enddo
           !--Collect also the additional detuning induced terms
           if(ide.eq.1) then
              do k1=1,mmultx
                 do k2=1,mmult
                    if(abs(facd10(k1,k2,1)).gt.pieni.or.                  &
                         &abs(facd10(k1,k2,2)).gt.pieni) then
                       ica=ica+1
                       if(ica.gt.mmultf) call prror(3,10,mmultf)
                       facd2(1,ica)=facd10(k1,k2,1)
                       facd2(2,ica)=facd10(k1,k2,2)
                       facd10(k1,k2,1)=zero
                       facd10(k1,k2,2)=zero
                       do ii=1,4
                          ifacd2(ii,ica)=itij0(k1,k2,ii)
                       enddo
                    endif
                 enddo
              enddo
           endif
           !--Add up all contributions related to a specific resonance
           do k1=1,ica
              do k2=k1+1,ica
                 if(ifacd2(1,k2).eq.ifacd2(1,k1).and.                    &
                      &ifacd2(2,k2).eq.ifacd2(2,k1).and.                                 &
                      &ifacd2(3,k2).eq.ifacd2(3,k1).and.                                 &
                      &ifacd2(4,k2).eq.ifacd2(4,k1)) then
                    facd2(1,k1)=facd2(1,k1)+facd2(1,k2)
                    facd2(2,k1)=facd2(2,k1)+facd2(2,k2)
                    facd2(1,k2)=zero
                    facd2(2,k2)=zero
                    do ii=1,4
                       ifacd2(ii,k2)=0
                    enddo
                 endif
              enddo
           enddo
           ica2=0
           !--Collect the final total of ica2 non-zero resonance terms
           do k1=1,ica
              if(abs(facd2(1,k1)).gt.pieni.or.                          &
                   &abs(facd2(2,k1)).gt.pieni) then
                 ica2=ica2+1
                 if(ica2.ne.k1) then
                    facd2(1,ica2)=facd2(1,k1)
                    facd2(2,ica2)=facd2(2,k1)
                    facd2(1,k1)=zero
                    facd2(2,k1)=zero
                    do ii=1,4
                       ifacd2(ii,ica2)=ifacd2(ii,k1)
                       ifacd2(ii,k1)=0
                    enddo
                 endif
              endif
           enddo
           !--Order the resonance terms
           call sortres(3,imax,ica2)
           write(6,10060)
           table_name ='distort_2_f_end '
           data_size = ica2
           call fit_table(table_name,data_size,table_size_78)
           do k1=1,data_size
              argij=two*((ifacd2(1,k1)-ifacd2(2,k1))*qx+                &
                   &(ifacd2(3,k1)-ifacd2(4,k1))*qy)
              facc=one-cos(argij)
              facs=sin(argij)
              ham(1,k1)=facc*facd2(1,k1)+facs*facd2(2,k1)
              ham(2,k1)=facc*facd2(2,k1)-facs*facd2(1,k1)
              double_to_write(5)=dsqrt(facd2(1,k1)*facd2(1,k1)+         &
                   &facd2(2,k1)*facd2(2,k1))
              write(6,10080) facd2(1,k1),facd2(2,k1),double_to_write(5),&
                   &ifacd2(1,k1),ifacd2(2,k1),ifacd2(3,k1),ifacd2(4,k1)
              if(iu_on.eq.1.or.iu_on.eq.3)                              &
                   &write(78,10110) ic,id,facd2(1,k1),facd2(2,k1),                    &
                   &double_to_write(5),ifacd2(1,k1),ifacd2(2,k1),                     &
                   &ifacd2(3,k1),ifacd2(4,k1)
              int_to_write(1) = ic
              int_to_write(2) = id
              double_to_write(3) = facd2(1,k1)
              double_to_write(4) = facd2(2,k1)
              int_to_write(6) = ifacd2(1,k1)
              int_to_write(7) = ifacd2(2,k1)
              int_to_write(8) = ifacd2(3,k1)
              int_to_write(9) = ifacd2(4,k1)
              call write_table(table_name,9,int_to_write,               &
                   &double_to_write)
              j78 = j78 + 1
              call augment_count(table_name)
              facd2(1,k1)=zero
              facd2(2,k1)=zero
           enddo
           write(6,10070)
           table_name ='distort_2_h_end '
           data_size = ica2
           call fit_table(table_name,data_size,table_size_79)
           do k1=1,data_size
              double_to_write(5)=dsqrt(ham(1,k1)*ham(1,k1)+             &
                   &ham(2,k1)*ham(2,k1))
              write(6,10080) ham(1,k1),ham(2,k1),double_to_write(5),    &
                   &ifacd2(1,k1),ifacd2(2,k1),ifacd2(3,k1),ifacd2(4,k1)
              if(iu_on.eq.1.or.iu_on.eq.3)                              &
                   &write(79,10110) ic,id,ham(1,k1),ham(2,k1),                        &
                   &double_to_write(5),ifacd2(1,k1),ifacd2(2,k1),                     &
                   &ifacd2(3,k1),ifacd2(4,k1)
              int_to_write(1) = ic
              int_to_write(2) = id
              double_to_write(3) = ham(1,k1)
              double_to_write(4) = ham(2,k1)
              int_to_write(6) = ifacd2(1,k1)
              int_to_write(7) = ifacd2(2,k1)
              int_to_write(8) = ifacd2(3,k1)
              int_to_write(9) = ifacd2(4,k1)
              call write_table(table_name,9,int_to_write,               &
                   &double_to_write)
              j79 = j79 + 1
              call augment_count(table_name)
              do ii=1,4
                 ifacd2(ii,k1)=0
                 ham(1,k1)=zero
                 ham(2,k1)=zero
              enddo
           enddo
           call timex(tim4)
           write(6,10030) tim4-tim3
           write(6,10020) tim4-tim2
        endif
     enddo
  enddo
  !---------------------------------------------------------------------
  call timex(tim5)
  write(6,10040) tim5-tim2
  return
  !---------------------------------------------------------------------
10020 format(/80('-')//'Intermediate Time in  <DISTORT2>: ',            &
       &f10.3,' second(s)',' of Execution Time'                           &
       &//80('-')//)
10030 format(/80('-')//'Second Order Drivingterm Calculation took ',    &
       &f10.3,' second(s)',' of Execution Time'//80('-')//)
10040 format(/80('-')//'Program <DISTORT2> used ',f10.3,                &
       &' second(s) of Execution Time'//80('-')//)
10060 format(/80('-')/10x,'Distortionfunction Terms (second order)'/    &
       &80('-')/4x,'Cosine-Term',8x,'Sine-Term',8x,'Amplitude',7x,        &
       &'J',5x,'K',5x,'L',5x,'M'/80('-'))
10070 format(/80('-')/10x,'Hamiltonian Terms (second order)'/80('-')/   &
       &4x,'Cosine-Term',8x,'Sine-Term',8x,'Amplitude',7x,                &
       &'J',5x,'K',5x,'L',5x,'M'/80('-'))
10080 format(3(1pe17.9),4i6)
10090 format(//80('-')//10x,a,i8,2a,i8,a)
10100 format(//80('-')//10x,a,i8,a)
10110 format(i4,4x,i4,4x,3(1pe23.15),4i4)
end subroutine distort2
subroutine caldet2(ia,ib,imuty)
  !---------------------------------------------------------------------
  !--
  !--Calculation of the poisson bracket to derive the detuning function
  !--Second order in the strength of the multipoles
  !--
  !---------------------------------------------------------------------
  use sodd
  implicit none
  integer ia,ib,ic,id,ido,if2,ihoff,ihw,ii,ii1,ii11,iiend,iii,iii1, &
       &ijnu,ijnu1,imuty,ipl,iti,ivoff,jj,jj1,jj11,jjend,jjj,jjj1,k,k1,kk,&
       &l,l1,m,mend,mexp,mij,mm,nend,nexp,nij,nn,n
  double precision betl,betxp,fac,fac1,fac3,fcc2,fcc4,signii, &
       &signjj
  !---------------------------------------------------------------------
  dimension betl(mmultx,4),betxp(2,mmultx,4)
  dimension ihw(2,mmult),ipl(2,mmult,2),iti(2,mmultx,mmult,6)
  dimension fcc2(2,mmultx),fcc4(2,mmultx,mmult)
  do k=1,mmultx
     fcc2(1,k)=zero
     fcc2(2,k)=zero
     do m=1,mmult
        do nn=1,2
           ihw(nn,k)=0
           ipl(nn,k,1)=0
           ipl(nn,k,2)=0
           iti(nn,k,m,1)=0
           iti(nn,k,m,2)=0
           iti(nn,k,m,3)=0
           betl(k,1)=zero
           betl(k,2)=zero
           betl(k,3)=zero
           betl(k,4)=zero
           betxp(nn,k,1)=zero
           betxp(nn,k,2)=zero
           betxp(nn,k,3)=zero
           betxp(nn,k,4)=zero
           fcc4(nn,k,m)=zero
        enddo
     enddo
  enddo
  fac=-two/dble(4*2**((ia+ib)/2)*ia*ib)
  if(ia.ne.ib) fac=fac*two
  iiend=(ia+2)/2
  jjend=(ib+2)/2
  signii=-one
  if(imuty.eq.2.and.2*iiend.eq.ia+2) iiend=iiend-1
  do ii=1,iiend
     ii1=(ii-1)*2
     if(imuty.eq.2) ii1=ii1+1
     signii=signii*(-one)
     signjj=-one
     if(imuty.eq.2.and.2*jjend.eq.ib+2) jjend=jjend-1
     do jj=1,jjend
        if2=0
        jj1=(jj-1)*2
        if(imuty.eq.2) jj1=jj1+1
        signjj=signjj*(-one)
        fac1=fac*signii*signjj*fact(ia,1+ii1)*fact(ib,1+jj1)
        iii=ia-ii1
        jjj=ib-jj1
        ijnu=(ii-1)*jjend+jj
        betl(ijnu,1)=dble(iii)/two
        betl(ijnu,2)=dble(ii1)/two
        betl(ijnu,3)=dble(jjj)/two
        betl(ijnu,4)=dble(jj1)/two
        ihoff=0
        ivoff=0
        if(iii.eq.0.and.jjj.ne.0) then
           fac1=fac1*fact(jjj,jjj/2+1)
           ihoff=1
        endif
        if(iii.ne.0.and.jjj.eq.0) then
           fac1=fac1*fact(iii,iii/2+1)
           ihoff=1
        endif
        if(ii1.eq.0.and.jj1.ne.0) then
           fac1=fac1*fact(jj1,jj1/2+1)
           ivoff=1
        endif
        if(ii1.ne.0.and.jj1.eq.0) then
           fac1=fac1*fact(ii1,ii1/2+1)
           ivoff=1
        endif
        k=(iii+jjj)/2
        k1=k-1
        l=(ii1+jj1)/2
        l1=l-1
        mij=min(iii,jjj)
        nij=min(ii1,jj1)
        mend=(mij+2)/2
        nend=(nij+2)/2
        iii1=0
        jjj1=0
        ii11=0
        jj11=0
        if(iii.gt.mij) iii1=(iii+2)/2-mend
        if(jjj.gt.mij) jjj1=(jjj+2)/2-mend
        if(ii1.gt.nij) ii11=(ii1+2)/2-nend
        if(jj1.gt.nij) jj11=(jj1+2)/2-nend
        !---------------------------------------------------------------------
        !--
        !--Deriving the Derivatives
        !--
        !---------------------------------------------------------------------
        if(k.gt.0) then
           if(l.gt.0.and.ivoff.eq.0.and.ihoff.eq.0) then
              !--
              !--Qx Term first Derivative by Iy for the Detuning by Ix
              !--
              ido=1
              !--Finding Double Appearance
              do ijnu1=1,ijnu-1
                 if(                                                     &
                      &betl(ijnu1,1).eq.betl(ijnu,3).and.                                &
                      &betl(ijnu1,2).eq.betl(ijnu,4).and.                                &
                      &betl(ijnu1,3).eq.betl(ijnu,1).and.                                &
                      &betl(ijnu1,4).eq.betl(ijnu,2)) then
                    ido=0
                    if2=if2+1
                    if(if2.gt.2) call prror(1,13,if2)
                    fcc2(if2,ijnu1)=two*fcc2(if2,ijnu1)
                 endif
              enddo
              if(ido.eq.1) then
                 if2=if2+1
                 if(if2.gt.2) call prror(1,13,if2)
                 ihw(if2,ijnu)=1
                 ipl(if2,ijnu,1)=k1
                 ipl(if2,ijnu,2)=l1
                 betxp(if2,ijnu,1)=betl(ijnu,1)
                 betxp(if2,ijnu,2)=betl(ijnu,2)
                 betxp(if2,ijnu,3)=betl(ijnu,3)
                 betxp(if2,ijnu,4)=betl(ijnu,4)
                 fcc2(if2,ijnu)=fac1*dble(k*ii1*jj1)/dble(2**(k1+l1))
                 do m=1,mend
                    mexp=mij-(m-1)*2
                    fac3=fact(iii,iii1+m)*fact(jjj,jjj1+m)
                    do n=1,nend
                       nexp=nij-(n-1)*2
                       nn=(m-1)*nend+n
                       fcc4(if2,ijnu,nn)=                                  &
                            &fac3*fact(ii1,ii11+n)*fact(jj1,jj11+n)*                           &
                            &dble(nexp)*(one/dble(ii1)+one/dble(jj1))/two
                       if(mexp.ne.0) then
                          iti(if2,ijnu,nn,1)=-1
                          iti(if2,ijnu,nn,2)=mexp
                          iti(if2,ijnu,nn,3)=nexp
                       else
                          iti(if2,ijnu,nn,1)=0
                          iti(if2,ijnu,nn,2)=0
                          iti(if2,ijnu,nn,3)=nexp
                       endif
                    enddo
                 enddo
              endif
           endif
           if(k1.gt.0.and.ihoff.eq.0) then
              !--
              !--Qx Term first Derivative by Ix for the Detuning by Ix
              !--
              ido=1
              !--Finding Double Appearance
              do ijnu1=1,ijnu-1
                 if(                                                     &
                      &betl(ijnu1,1).eq.betl(ijnu,3).and.                                &
                      &betl(ijnu1,2).eq.betl(ijnu,4).and.                                &
                      &betl(ijnu1,3).eq.betl(ijnu,1).and.                                &
                      &betl(ijnu1,4).eq.betl(ijnu,2)) then
                    ido=0
                    if2=if2+1
                    if(if2.gt.2) call prror(1,13,if2)
                    fcc2(if2,ijnu1)=two*fcc2(if2,ijnu1)
                 endif
              enddo
              if(ido.eq.1) then
                 if2=if2+1
                 if(if2.gt.2) call prror(1,13,if2)
                 ihw(if2,ijnu)=1
                 ipl(if2,ijnu,1)=k1-1
                 ipl(if2,ijnu,2)=l
                 betxp(if2,ijnu,1)=betl(ijnu,1)
                 betxp(if2,ijnu,2)=betl(ijnu,2)
                 betxp(if2,ijnu,3)=betl(ijnu,3)
                 betxp(if2,ijnu,4)=betl(ijnu,4)
                 fcc2(if2,ijnu)=fac1*dble(k1*iii*jjj)/                   &
                      &dble(2**(k1+l1))
                 !--Terms Iy**n where n>0
                 if(ii1.ne.0.and.jj1.ne.0) then
                    do n=1,nend
                       nexp=nij-(n-1)*2
                       fac3=fact(ii1,ii11+n)*fact(jj1,jj11+n)
                       do m=1,mend
                          mexp=mij-(m-1)*2
                          mm=(n-1)*mend+m
                          fcc4(if2,ijnu,mm)=                                &
                               &fac3*fact(iii,iii1+m)*fact(jjj,jjj1+m)*                           &
                               &dble(mexp)*(one/dble(iii)+one/dble(jjj))/two
                          if(abs(fcc4(if2,ijnu,mm)).gt.pieni) then
                             if(nexp.ne.0) then
                                iti(if2,ijnu,mm,1)=1
                                iti(if2,ijnu,mm,2)=mexp
                                iti(if2,ijnu,mm,3)=nexp
                             else if(mexp.ne.0) then
                                iti(if2,ijnu,mm,1)=0
                                iti(if2,ijnu,mm,2)=mexp
                                iti(if2,ijnu,mm,3)=0
                             endif
                          endif
                       enddo
                    enddo
                    !--Terms Iy**0=1
                 else
                    do m=1,mend
                       mexp=mij-(m-1)*2
                       fcc4(if2,ijnu,m)=                                   &
                            &fact(iii,iii1+m)*fact(jjj,jjj1+m)*                                &
                            &dble(mexp)*(one/dble(iii)+one/dble(jjj))/two
                       if(abs(fcc4(if2,ijnu,m)).gt.pieni) then
                          iti(if2,ijnu,m,1)=0
                          iti(if2,ijnu,m,2)=mexp
                          iti(if2,ijnu,m,3)=0
                       endif
                    enddo
                 endif
              endif
           else if(ihoff.eq.1.and.ivoff.eq.0) then
              !--
              !--Horizontal Term only one x, first Derivative by Iy and then Ix
              !--
              ido=1
              !--Finding Double Appearance
              do ijnu1=1,ijnu-1
                 if(                                                     &
                      &betl(ijnu1,1).eq.betl(ijnu,3).and.                                &
                      &betl(ijnu1,2).eq.betl(ijnu,4).and.                                &
                      &betl(ijnu1,3).eq.betl(ijnu,1).and.                                &
                      &betl(ijnu1,4).eq.betl(ijnu,2)) then
                    ido=0
                    if2=if2+1
                    if(if2.gt.2) call prror(1,13,if2)
                    fcc2(if2,ijnu1)=two*fcc2(if2,ijnu1)
                 endif
              enddo
              if(ido.eq.1) then
                 if2=if2+1
                 if(if2.gt.2) call prror(1,13,if2)
                 ihw(if2,ijnu)=1
                 ipl(if2,ijnu,1)=k1
                 ipl(if2,ijnu,2)=l1
                 betxp(if2,ijnu,1)=betl(ijnu,1)
                 betxp(if2,ijnu,2)=betl(ijnu,2)
                 betxp(if2,ijnu,3)=betl(ijnu,3)
                 betxp(if2,ijnu,4)=betl(ijnu,4)
                 fcc2(if2,ijnu)=fac1*dble(k*ii1*jj1)/dble(2**(k1+l1))
                 do n=1,nend
                    nexp=nij-(n-1)*2
                    fcc4(if2,ijnu,n)=                                     &
                         &fact(ii1,ii11+n)*fact(jj1,jj11+n)*                                &
                         &dble(nexp)*(one/dble(ii1)+one/dble(jj1))/two
                    iti(if2,ijnu,n,1)=0
                    iti(if2,ijnu,n,2)=0
                    iti(if2,ijnu,n,3)=nexp
                 enddo
              endif
           else if(iii.eq.1.and.jjj.eq.1) then
              !--
              !--Vertical Term for uneven Multipoles
              !--
              if2=if2+1
              if(if2.gt.2) call prror(1,13,if2)
              ihw(if2,ijnu)=2
              ipl(if2,ijnu,1)=0
              ipl(if2,ijnu,2)=l-1
              betxp(if2,ijnu,1)=betl(ijnu,1)
              betxp(if2,ijnu,2)=betl(ijnu,2)
              betxp(if2,ijnu,3)=betl(ijnu,3)
              betxp(if2,ijnu,4)=betl(ijnu,4)
              fcc2(if2,ijnu)=fac1*dble(l*iii*jjj)/dble(2**(l-1))
              do n=1,nend
                 fac3=fact(ii1,ii11+n)*fact(jj1,jj11+n)
                 nexp=nij-(n-1)*2
                 fcc4(if2,ijnu,n)=fac3*                                  &
                      &(one/dble(iii)+one/dble(jjj))/two
                 if(nexp.ne.0) then
                    iti(if2,ijnu,n,1)=1
                    iti(if2,ijnu,n,2)=1
                    iti(if2,ijnu,n,3)=nexp
                 else
                    iti(if2,ijnu,n,1)=0
                    iti(if2,ijnu,n,2)=1
                    iti(if2,ijnu,n,3)=nexp
                 endif
              enddo
           endif
           !--
           !--Vertical Term for even Multipoles
           !--
        else if(ivoff.eq.0) then
           if2=if2+1
           if(if2.gt.2) call prror(1,13,if2)
           ihw(if2,ijnu)=2
           ipl(if2,ijnu,1)=0
           ipl(if2,ijnu,2)=l1-1
           betxp(if2,ijnu,1)=betl(ijnu,1)
           betxp(if2,ijnu,2)=betl(ijnu,2)
           betxp(if2,ijnu,3)=betl(ijnu,3)
           betxp(if2,ijnu,4)=betl(ijnu,4)
           fcc2(if2,ijnu)=fac1*dble(l1*ii1*jj1)/dble(2**(l-2))
           do n=1,nend
              nexp=nij-(n-1)*2
              fcc4(if2,ijnu,n)=                                         &
                   &fact(ii1,ii11+n)*fact(jj1,jj11+n)*                                &
                   &dble(nexp)*(one/dble(ii1)+one/dble(jj1))/two
              iti(if2,ijnu,n,1)=0
              iti(if2,ijnu,n,2)=0
              iti(if2,ijnu,n,3)=nexp
           enddo
        endif
     enddo
  enddo
  !--Write out Parameters properly
  ic=0
  do nn=1,2
     do kk=1,mmultx
        if(abs(fcc2(nn,kk)).gt.pieni) then
           id=0
           ic=ic+1
           if(ic.gt.mmultx) call prror(1,6,mmultx)
           ihv(ic)=ihw(nn,kk)
           iplane(ic,1)=ipl(nn,kk,1)
           iplane(ic,2)=ipl(nn,kk,2)
           betexp(ic,1)=betxp(nn,kk,1)
           betexp(ic,2)=betxp(nn,kk,2)
           betexp(ic,3)=betxp(nn,kk,3)
           betexp(ic,4)=betxp(nn,kk,4)
           fac2(ic)=fcc2(nn,kk)
           do mm=1,mmult
              if(abs(fcc4(nn,kk,mm)).gt.pieni) then
                 id=id+1
                 if(id.gt.mmult) call prror(1,7,mmult)
                 itij(ic,id,1)=iti(nn,kk,mm,1)
                 itij(ic,id,2)=iti(nn,kk,mm,2)
                 itij(ic,id,3)=iti(nn,kk,mm,3)
                 fac4(ic,id)=fcc4(nn,kk,mm)
              endif
           enddo
           icd(ic)=id
        endif
     enddo
     icc=ic
  enddo
  return
end subroutine caldet2
subroutine caldt2(ia,ib,imuty)
  !---------------------------------------------------------------------
  !--
  !--Calculation of the poisson bracket to derive the distortion function
  !--Second order in the strength of the multipoles
  !--
  !---------------------------------------------------------------------
  use sodd
  implicit none
  integer ia,ib,ic,ica,icol1,icol2,icol3,icol4,icol5,icol6,id,ido1, &
       &ido2,ido3,ido4,ihoff,ihoff1,ihoff2,ihvall,ii,ii1,iiend,iii,ijnu,  &
       &imuty,isiga,iti,ivoff,ivoff1,ivoff2,jj,jj1,jjend,jjj,k,k1,k2,k3,  &
       &k4,k41,k412s,k412s1,k412s2,k41s,k42,k42s,k43,k434s,k434s1,k434s2, &
       &k44,k5,kk,l1,m,mm,nn,nnanf,nnend,nsig
  double precision betl,diab,fac11,fac12,fac21,fac22,faca,facc,     &
       &facc1,facc2,fcc2,fcc4,signii,signjj
  dimension iti(2,mmultx,mmult,6),isiga(2,mmultx,mmult)
  dimension fcc2(2,mmultx),fcc4(2,mmultx,mmult)
  dimension betl(mmultx,4)
  !---------------------------------------------------------------------
  do k=1,mmultx
     fcc2(1,k)=zero
     fcc2(2,k)=zero
     do m=1,mmult
        do nn=1,2
           iti(nn,k,m,1)=0
           iti(nn,k,m,2)=0
           iti(nn,k,m,3)=0
           iti(nn,k,m,4)=0
           iti(nn,k,m,5)=0
           iti(nn,k,m,6)=0
           isiga(nn,k,m)=0
           betl(k,1)=zero
           betl(k,2)=zero
           betl(k,3)=zero
           betl(k,4)=zero
           fcc4(nn,k,m)=zero
        enddo
     enddo
  enddo
  diab=dble(ia+ib)/two
  facc=one/four/(two**diab)/dble(ia*ib)
  !--Prefactor I (General factors)
  iiend=(ia+2)/2
  jjend=(ib+2)/2
  !--(iiend,jjend)=>Total number of terms in (x+iy)**n
  !--(i,j) is on left or right side of the Poisson bracket respectively
  signii=-one
  if((imuty.eq.4.or.imuty.eq.2).and.2*iiend.eq.ia+2) iiend=iiend-1
  do ii=1,iiend
     ii1=(ii-1)*2
     if(imuty.eq.4.or.imuty.eq.2) ii1=ii1+1
     signii=signii*(-one)
     signjj=-one
     if((imuty.eq.3.or.imuty.eq.2).and.2*jjend.eq.ib+2) jjend=jjend-1
     do jj=1,jjend
        jj1=(jj-1)*2
        !--(ii1,jj1)=>Increment of powers (0->{(iiend,jjend)-1}*2)
        !--           goes in steps of 2
        !--These are the exponents of the second variable y
        if(imuty.eq.3.or.imuty.eq.2) jj1=jj1+1
        signjj=signjj*(-one)
        facc1=facc*signii*signjj*fact(ia,1+ii1)*fact(ib,1+jj1)
        !--Prefactor II (signs, binominal coefficients)
        iii=ia-ii1
        jjj=ib-jj1
        !--(iii,jjj)=>Exponent of the first variable x
        ijnu=(ii-1)*jjend+jj
        if(ijnu.gt.mmultx) call prror(3,6,mmultx)
        !--inju=># of current case treated of a total of (iiend*jjend)
        betl(ijnu,1)=dble(iii)/two
        betl(ijnu,2)=dble(ii1)/two
        betl(ijnu,3)=dble(jjj)/two
        betl(ijnu,4)=dble(jj1)/two
        ihoff=0
        ihoff1=0
        ihoff2=0
        ivoff=0
        ivoff1=0
        ivoff2=0
        if(iii.eq.0.and.jjj.eq.0) ihoff=1
        if(iii.eq.0.and.jjj.ne.0.and.ihoff.eq.0) ihoff1=1
        if(iii.ne.0.and.jjj.eq.0.and.ihoff.eq.0) ihoff2=1
        if(ii1.eq.0.and.jj1.eq.0) ivoff=1
        if(ii1.eq.0.and.jj1.ne.0.and.ivoff.eq.0) ivoff1=1
        if(ii1.ne.0.and.jj1.eq.0.and.ivoff.eq.0) ivoff2=1
        ihvall=ihoff+ihoff1+ihoff2+ivoff+ivoff1+ivoff2
        !--(ihoff,ivoff)=0/1 (x,y) on (at least on one side/on no sides) of
        !--                          Poisson bracket
        !--(ihoff1,ivoff1)=0/1 (x,y) on (both sides/on left side) of
        !--                          Poisson bracket
        !--(ihoff2,ivoff2)=0/1 (x,y) on (both sides/on right side) of
        !--                          Poisson bracket
        facc2=-facc1/(two**(diab-one))/four
        !--Prefactor III (Invariant->emittance) and -1/4d0 from evaluation of
        !--             trigonometric functions
        !---------------------------------------------------------------------
        !--
        !--Deriving the Derivatives
        !--
        !---------------------------------------------------------------------
        nnanf=1
        nnend=2
        if(ihoff.eq.1.or.ihoff1.eq.1.or.ihoff2.eq.1) then
           nnanf=2
           nnend=2
        endif
        if(ivoff.eq.1.or.ivoff1.eq.1.or.ivoff2.eq.1) nnend=1
        if(ihvall.eq.0) then
           nnanf=1
           nnend=2
        endif
        !--Differentiation nn=1 => horizontal nn=2 => vertical
        !-- ==> exchange planes (ido1,ido2 =>  ido3,ido4)
        !-- ==> exchange terms of secondary resonances (icol1 => icol2)
        !-- ==> exchange primary resonance terms (icol3,icol4 => icol5,icol6)
        do nn=nnanf,nnend
           if(nn.eq.1) then
              ido1=iii
              ido2=jjj
              ido3=ii1
              ido4=jj1
              icol1=1
              icol2=2
              icol3=3
              icol4=4
              icol5=5
              icol6=6
              !--nsig see below
              nsig=3
           endif
           if(nn.eq.2) then
              ido1=ii1
              ido2=jj1
              ido3=iii
              ido4=jjj
              icol1=2
              icol2=1
              icol3=5
              icol4=6
              icol5=3
              icol6=4
              !--nsig is set to 4 for the vertical differentiation as the
              !--horizontal stays always positive
              nsig=4
           endif
           fcc2(nn,ijnu)=facc2*dble(ido1*ido2)
           ica=0
           !--(k1,k2) always differentiation nn=1 horizontal nn=2 vertical
           do k1=0,ido1,2
              k41=ido1-k1
              !--k41 see k44
              fac11=dble(k41)/dble(ido1)*fact(ido1,1+k1/2)
              fac12=fact(ido1,1+k1/2)
              do k2=0,ido2,2
                 k42=ido2-k2
                 !--k42 see k44
                 fac21=fact(ido2,1+k2/2)
                 fac22=dble(k42)/dble(ido2)*fact(ido2,1+k2/2)
                 !--(l1) sums over the 3-4 possible sign combinations
                 do l1=1,nsig
                    faca=sign(l1,1)*fac11*fac21-                          &
                         &sign(l1,2)*fac12*fac22
                    if(k41.eq.0) faca=-sign(l1,2)*fac12*fac22
                    if(k42.eq.0) faca=sign(l1,1)*fac11*fac21
                    !--Condition 1 & 5 clear
                    !--2 Condition avoids duplication of l1 by l2
                    !--3 Condition avoids duplication of l1 by l3
                    !--4 Condition avoids duplication of l1 with any other
                    if((abs(faca).gt.pieni).and.                          &
                         &(l1.ne.2.or.k41.ne.0).and.                                        &
                         &(l1.ne.3.or.k42.ne.0).and.                                        &
                         &(l1.ne.4.or.(k41.ne.0.and.k42.ne.0)).and.                         &
                         &((k41+k42).ge.1)) then
                       !--(k3,k4) plane with no differentiation nn=1 vertical nn=2 horizontal
                       !--((ido3+3),(ido4+3)) allow for constant term in the
                       !--undifferentiated part
                       do k3=0,2*ido3,2
                          do k4=0,2*ido4,2
                             !--k43,k44 control the secondary resonance
                             !--i.e. the same that appear in first order
                             k43=ido3-k3
                             k44=ido4-k4
                             !--(k41s,k42s) takes care of the sign of k41 and k42
                             k41s=sign(l1,1)*k41
                             k42s=sign(l1,2)*k42
                             !--(k412s) sum of the 2 terms with differentiation
                             k412s=k41s+k42s
                             k412s1=(ido1+ido2+k412s)/2-1
                             k412s2=k412s1-k412s
                             !--(k434s) sum of the 2 terms without differentiation
                             k434s=k43+k44
                             k434s1=(ido3+ido4+k434s)/2
                             k434s2=k434s1-k434s
                             !--1. (0,0) Secondary Resonance has to be excluded
                             !--2. Select resonances in the desired nomenclature
                             !--at least the sum must be larger than 0
                             !--and the resonance should be written in descending order starting
                             !--with the horizontal one
                             if(                                             &
                                  &(k42s.ne.0.or.k44.ne.0).and.                                      &
                                  &(k412s.gt.0.or.k434s.gt.0).and.                                   &
                                  &((nn.eq.1.and.k412s1.ge.k412s2).or.                               &
                                  &(nn.eq.2.and.k434s1.ge.k434s2))) then
                                !--Program check: finds duplicated cases
                                do k5=1,ica
                                   if(                                         &
                                        &iti(nn,ijnu,k5,icol1).eq.k42s.and.                                &
                                        &iti(nn,ijnu,k5,icol2).eq.k44.and.                                 &
                                        &iti(nn,ijnu,k5,icol3).eq.k412s1.and.                              &
                                        &iti(nn,ijnu,k5,icol4).eq.k412s2.and.                              &
                                        &iti(nn,ijnu,k5,icol5).eq.k434s1.and.                              &
                                        &iti(nn,ijnu,k5,icol6).eq.k434s2                                   &
                                        &) call prror(3,9,2)
                                enddo
                                ica=ica+1
                                if(ica.gt.mmult) call prror(3,5,mmult)
                                fcc4(nn,ijnu,ica)=faca*                       &
                                     &fact(ido3,1+k3/2)*fact(ido4,1+k4/2)
                                iti(nn,ijnu,ica,icol1)=k42s
                                iti(nn,ijnu,ica,icol2)=k44
                                iti(nn,ijnu,ica,icol3)=k412s1
                                iti(nn,ijnu,ica,icol4)=k412s2
                                iti(nn,ijnu,ica,icol5)=k434s1
                                iti(nn,ijnu,ica,icol6)=k434s2
                                if(k41.eq.0.and.k43.eq.0) isiga(nn,ijnu,ica)=1
                             endif
                          enddo
                       enddo
                    endif
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  !--Write out Parameters properly
  ic=0
  do nn=1,2
     do kk=1,mmultx
        if(abs(fcc2(nn,kk)).gt.pieni) then
           id=0
           ic=ic+1
           if(ic.gt.mmultx) call prror(3,6,mmultx)
           betexp(ic,1)=betl(kk,1)
           betexp(ic,2)=betl(kk,2)
           betexp(ic,3)=betl(kk,3)
           betexp(ic,4)=betl(kk,4)
           fac2(ic)=fcc2(nn,kk)
           do mm=1,mmult
              if(abs(fcc4(nn,kk,mm)).gt.pieni) then
                 id=id+1
                 if(id.gt.mmult) call prror(3,7,mmult)
                 if(isiga(nn,kk,mm).eq.1) then
                    isig(ic,id)=1
                    isiga(nn,kk,mm)=0
                 endif
                 itij(ic,id,1)=iti(nn,kk,mm,1)
                 itij(ic,id,2)=iti(nn,kk,mm,2)
                 itij(ic,id,3)=iti(nn,kk,mm,3)
                 itij(ic,id,4)=iti(nn,kk,mm,4)
                 itij(ic,id,5)=iti(nn,kk,mm,5)
                 itij(ic,id,6)=iti(nn,kk,mm,6)
                 fac4(ic,id)=fcc4(nn,kk,mm)
              endif
           enddo
           icd(ic)=id
        endif
     enddo
     icc=ic
  enddo
  return
end subroutine caldt2
subroutine init
  !---------------------------------------------------------------------
  !--Initilisation
  !---------------------------------------------------------------------
  use sodd
  implicit none
  integer i,ii,j,k
  !---------------------------------------------------------------------
  tlim=1e7
  zero=dble(0d0)
  one=1d0
  two=2d0
  four=4d0
  pi=atan(one)*four
  pi2=pi*two
  pihi=two/pi
  comment(2)=' Erect Quadrupoles'
  comment(-2)=' Skew Quadrupoles '
  comment(3)=' Erect Sextupoles '
  comment(-3)=' Skew Sextupoles '
  comment(4)=' Erect Octupoles '
  comment(-4)=' Skew Octupoles '
  comment(5)=' Erect Decapoles '
  comment(-5)=' Skew Decapoles '
  comment(6)=' Erect Dodecapoles'
  comment(-6)=' Skew Dodecapoles '
  comment(7)=' Erect 14th-Pole '
  comment(-7)=' Skew 14th-Pole '
  comment(8)=' Erect 16th-Pole '
  comment(-8)=' Skew 16th-Pole '
  comment(9)=' Erect 18th-Pole '
  comment(-9)=' Skew 18th-Pole '
  comment(10)=' Erect 20th-Pole '
  comment(-10)=' Skew 20th-Pole '
  comment(11)=' Erect 22th-Pole '
  comment(-11)=' Skew 22th-Pole '
  !--
  !--Binominals
  !--
  cosav(0)=one
  do i=2,mmul,2
     ii=i/2
     cosav(ii)=cosav(ii-1)*dble(i-1)/dble(i)
  enddo
  do j=0,mmul
     do i=0,mmul2
        fact(j,i)=one
     enddo
  enddo
  do i=2,mmul
     do j=2,mmul
        if(j.ge.i) then
           fact(j,i)=fact(j-1,i)+fact(j-1,i-1)
        endif
     enddo
  enddo
  do i=-mmul,mmul
     do j=0,mmul
        do k=0,mmul
           det(1,i,j,k)=zero
           det(2,i,j,k)=zero
        enddo
     enddo
  enddo
  do i=0,mmul
     det1(1,i)=zero
     det1(2,i)=zero
  enddo
  !--
  !--sign array
  !--
  sign(1,1)=one
  sign(1,2)=one
  sign(2,1)=-one
  sign(2,2)=one
  sign(3,1)=one
  sign(3,2)=-one
  sign(4,1)=-one
  sign(4,2)=-one
  return
end subroutine init
subroutine readdat
  !---------------------------------------------------------------------
  !--
  !--Reading Lattice Information from File # 34
  !--
  !---------------------------------------------------------------------
  use sodd
  implicit none
  integer i,ich1,ich2,ii,isize,itype,j,jj,ndum,num
  double precision betx,bety,et,phix0,phiy0,str
  character*16 name
  character*200 ch
  character*205 ch1
  double precision factor
  !---------------------------------------------------------------------

  do i=1,nblz
     do j=-mmul,mmul
        ityc(j)=0
        bstr(j,i)=zero
        strn(j,i)=" "
        beta(1,j,i)=zero
        beta(2,j,i)=zero
        sbeta(1,j,i)=zero
        sbeta(2,j,i)=zero
        phi(1,j,i)=zero
        phi(2,j,i)=zero
        etl(j,i)=zero
     enddo
  enddo
  do ii=1,10000000
     read(34,'(A150)',end=900) ch
     ich1=0
     ich2=0
     do jj=1,200
        if(ich1.eq.0.and.ch(jj:jj).ne.' ') ich1=1
        if(ich1.eq.1.and.ch(jj:jj).eq.' ') ich1=2
        if(ich1.eq.2.and.ch(jj:jj).ne.' ') then
           ich1=3
           ich2=jj
        endif
        if(ich1.eq.3.and.ch(jj:jj).eq.' ') then
           ch1(1:200)=ch(1:ich2-1)//''''//ch(ich2:jj-1)//''''//        &
                &ch(jj:200)//' / '
           goto 210
        endif
     enddo
210  continue
     read(ch1,*) et,name,itype,str,betx,bety,phix0,phiy0
     if(itype.eq.100) then
        qx=phix0*pi
        qy=phiy0*pi
        goto 5
     endif
     if(et.ge.etl1.and.et.le.etl2) then
        if(abs(str).ge.pieni.and.itype.ge.n1.and.itype.le.n2) then
           ityc(itype)=ityc(itype)+1
           isize=ityc(itype)
           if(isize.gt.nblz) call prror(3,3,nblz)
           etl(itype,isize)=et
           strn(itype,isize)=name
           factor = 10**(3*(iabs(itype)-2))
           bstr(itype,isize)=str*factor
           beta(1,itype,isize)=betx
           beta(2,itype,isize)=bety
           sbeta(1,itype,isize)=dsqrt(betx)
           sbeta(2,itype,isize)=dsqrt(bety)
           phi(1,itype,isize)=phix0*two*pi
           phi(2,itype,isize)=phiy0*two*pi
        endif
     endif
  enddo
  call prror(3,4,num)
5 continue
  write(6,*) ' Qx= ',qx/pi,' Qy= ',qy/pi
  close(34)
  return
900 close(34)
  call prror(3,1,ndum)
end subroutine readdat
subroutine detwri(icase,imu,ih1)
  !---------------------------------------------------------------------
  !--Organize Writing for Detune1
  !---------------------------------------------------------------------
  use sodd
  implicit none
  integer i,ic,icase,id,ih1,imu,j,k
  character*16 table_name
  integer int_to_write(11)
  double precision double_to_write(11)
  dimension imu(2)
  integer data_size
  !---------------------------------------------------------------------
  ic=imu(1)
  id=imu(2)
  if(icase.eq.0) then

     !--- Zero order multiple strength order

     !---  handle unit 71 (detune_1_all)

     if(iu_on.eq.2.or.iu_on.eq.3) then
        table_name = 'detune_1_all '
        data_size = 2*ih1+2+j71
        call fit_table(table_name,data_size,table_size_71)
        do i=ih1,0,-1
           write(71,10040) ic,1,det(1,ic,i,ih1-i)/pi2,i,ih1-i
           int_to_write(1) = ic
           int_to_write(2) = 1
           double_to_write(3) = det(1,ic,i,ih1-i)/pi2
           int_to_write(4) = i
           int_to_write(5) = ih1-i
           call write_table(table_name,5,int_to_write,double_to_write)
           j71 = j71 + 1
           call augment_count(table_name)
        enddo
        do i=ih1-1,0,-1
           write(71,10040) ic,2,                                       &
                &det(1,ic,i,ih1-i)*dble(ih1-i)/dble(i+1)/pi2,i+1,ih1-i-1
           int_to_write(2) = 2
           double_to_write(3) = det(1,ic,i,ih1-i)*dble(ih1-i)/dble(i+1)&
                &/pi2
           int_to_write(4) = i+1
           int_to_write(5) = ih1-i-1
           call write_table(table_name,5,int_to_write,double_to_write)
           j71 = j71 + 1
           call augment_count(table_name)
        enddo
        write(71,10040) ic,2,det(2,ic,0,ih1)/pi2,0,ih1
        double_to_write(3) = det(2,ic,0,ih1)/pi2
        int_to_write(4) = 0
        int_to_write(5) = ih1
        call write_table(table_name,5,int_to_write,double_to_write)
        j71 = j71 + 1
        call augment_count(table_name)
     endif

     !--- First order multiple strength order

     !---  handle unit 70 (detune_1_end)

  else if(icase.eq.1) then
     write(6,10000)
     table_name = 'detune_1_end '
     data_size = 2*ih1+2+j70
     call fit_table(table_name,data_size,table_size_70)
     do i=ih1,0,-1
        write(6,10010) det(1,ic,i,ih1-i)/pi2,i,ih1-i
        if(iu_on.eq.1.or.iu_on.eq.3) then
           write(70,10040) ic,1,det(1,ic,i,ih1-i)/pi2,i,ih1-i
           int_to_write(1) = ic
           int_to_write(2) = 1
           double_to_write(3) = det(1,ic,i,ih1-i)/pi2
           int_to_write(4) = i
           int_to_write(5) = ih1-i
           call write_table(table_name,5,int_to_write,double_to_write)
           j70 = j70 + 1
           call augment_count(table_name)
        endif
     enddo
     write(6,10020)
     do i=ih1-1,0,-1
        write(6,10010) det(1,ic,i,ih1-i)*dble(ih1-i)/dble(i+1)/pi2,   &
             &i+1,ih1-i-1
        if(iu_on.eq.1.or.iu_on.eq.3) then
           write(70,10040) ic,2,det(1,ic,i,ih1-i)*dble(ih1-i)/dble(i+1)&
                &/pi2,i+1,ih1-i-1
           int_to_write(2) = 2
           double_to_write(3) = det(1,ic,i,ih1-i)*dble(ih1-i)/dble(i+1)&
                &/pi2
           int_to_write(4) = i+1
           int_to_write(5) = ih1-i-1
           call write_table(table_name,5,int_to_write,double_to_write)
           j70 = j70 + 1
           call augment_count(table_name)
        endif
     enddo
     write(6,10010) det(2,ic,0,ih1)/pi2,0,ih1
     if(iu_on.eq.1.or.iu_on.eq.3) then
        write(70,10040) ic,2,det(2,ic,0,ih1)/pi2,0,ih1
        double_to_write(3) = det(2,ic,0,ih1)/pi2
        int_to_write(4) = 0
        int_to_write(5) = ih1
        call write_table(table_name,5,int_to_write,double_to_write)
        j70 = j70 + 1
        call augment_count(table_name)
     endif
  else if(icase.eq.2) then

     !--- Second order multiple strength order

     !---  handle unit 72 (detune_2_hor)

     if(iu_on.eq.1.or.iu_on.eq.3) then
        write(6,10000)
        table_name = 'detune_2_hor '
        data_size = ih1+1+j72
        call fit_table(table_name,data_size,table_size_72)
     endif
     do i=ih1,0,-1
        write(6,10010) det(1,ic,i,ih1-i)/pi2,i,ih1-i
        if(iu_on.eq.1.or.iu_on.eq.3) then
           write(72,10050) ic,id,det(1,ic,i,ih1-i)/pi2,i,ih1-i
           int_to_write(1) = ic
           int_to_write(2) = id
           double_to_write(3) = det(1,ic,i,ih1-i)/pi2
           int_to_write(4) = i
           int_to_write(5) = ih1-i
           call write_table(table_name,5,int_to_write,double_to_write)
           j72 = j72 + 1
           call augment_count(table_name)
        endif
     enddo

     !---  handle unit 73 (detune_2_ver)

     if(iu_on.eq.1.or.iu_on.eq.3) then
        write(6,10020)
        table_name = 'detune_2_ver '
        data_size = ih1+2+j73
        call fit_table(table_name,data_size,table_size_73)
     endif
     do i=ih1-1,0,-1
        write(6,10010) det(1,ic,i,ih1-i)*dble(ih1-i)/                 &
             &dble(i+1)/pi2,i+1,ih1-i-1
        if(iu_on.eq.1.or.iu_on.eq.3) then
           write(73,10050) ic,id,                                      &
                &det(1,ic,i,ih1-i)*dble(ih1-i)/dble(i+1)/pi2,i+1,ih1-i-1
           int_to_write(1) = ic
           int_to_write(2) = id
           double_to_write(3) = det(1,ic,i,ih1-i)*dble(ih1-i)/dble(i+1)&
                &/pi2
           int_to_write(4) = i+1
           int_to_write(5) = ih1-i-1
           call write_table(table_name,5,int_to_write,double_to_write)
           j73 = j73 + 1
           call augment_count(table_name)
        endif
     enddo
     write(6,10010) det(2,ic,0,ih1)/pi2,0,ih1
     if(iu_on.eq.1.or.iu_on.eq.3) then
        write(73,10050) ic,id,                                        &
             &det(2,ic,0,ih1)/pi2,0,ih1
        double_to_write(3) = det(2,ic,0,ih1)/pi2
        int_to_write(4) = 0
        int_to_write(5) = ih1
        call write_table(table_name,5,int_to_write,double_to_write)
        j73 = j73 + 1
        call augment_count(table_name)
     endif
     do i=-mmul,mmul
        do j=0,mmul
           do k=0,mmul
              det(1,i,j,k)=zero
              det(2,i,j,k)=zero
           enddo
        enddo
     enddo
     do i=0,mmul
        det1(1,i)=zero
        det1(2,i)=zero
     enddo
  endif
  return
  !---------------------------------------------------------------------
10000 format(/80('-')/8x,'Horizontal Coefficient',4x,'E-x',3x,'E-y'/    &
       &80('-'))
10010 format(10x,1pe20.12,2i6)
10020 format(/80('-')/8x,'  Vertical Coefficient',4x,'E-x',3x,'E-y'/    &
       &80('-'))
10040 format(i4,3x,i4,2x,1pe23.15,i4,1x,i4)
10050 format(i4,4x,i4,3x,1pe23.15,1x,i4,1x,i4)

end subroutine detwri
subroutine sortres(icase,im1,ic)
  !---------------------------------------------------------------------
  !--
  !--Order the resonance terms
  !--
  !---------------------------------------------------------------------
  use sodd
  implicit none
  integer ic,icase,idum1,idum2,idum3,idum4,idum5,idum6,im,im1,is,   &
       &is0,k0,k1,k2,k3,k4,l2,l4,m
  double precision dum1,dum2,dum3,dum4
  !---------------------------------------------------------------------
  im=iabs(im1)
  is0=0
  is=0
  do k0=im,0,-1
     do k1=k0,0,-1
        do k2=k1,0,-1
           l2=k1-k2
           do k3=k0-k1,k0-k1,-1
              do k4=k3,0,-1
                 l4=k3-k4
                 do m=1,ic
                    if(icase.eq.2.and.                                    &
                         &k2.eq.ifact4(1,im1,m).and.l2.eq.ifact4(2,im1,m)                   &
                         &.and.k4.eq.ifact4(3,im1,m).and.                                   &
                         &l4.eq.ifact4(4,im1,m)) then
                       is0=is0+1
                       dum1=fact0(im1,is0)
                       dum2=factb(1,im1,is0)
                       dum3=factb(2,im1,is0)
                       idum1=ifacta(1,im1,is0)
                       idum2=ifacta(2,im1,is0)
                       idum3=ifact4(1,im1,is0)
                       idum4=ifact4(2,im1,is0)
                       idum5=ifact4(3,im1,is0)
                       idum6=ifact4(4,im1,is0)
                       fact0(im1,is0)=fact0(im1,m)
                       factb(1,im1,is0)=factb(1,im1,m)
                       factb(2,im1,is0)=factb(2,im1,m)
                       ifacta(1,im1,is0)=ifacta(1,im1,m)
                       ifacta(2,im1,is0)=ifacta(2,im1,m)
                       ifact4(1,im1,is0)=ifact4(1,im1,m)
                       ifact4(2,im1,is0)=ifact4(2,im1,m)
                       ifact4(3,im1,is0)=ifact4(3,im1,m)
                       ifact4(4,im1,is0)=ifact4(4,im1,m)
                       fact0(im1,m)=dum1
                       factb(1,im1,m)=dum2
                       factb(2,im1,m)=dum3
                       ifacta(1,im1,m)=idum1
                       ifacta(2,im1,m)=idum2
                       ifact4(1,im1,m)=idum3
                       ifact4(2,im1,m)=idum4
                       ifact4(3,im1,m)=idum5
                       ifact4(4,im1,m)=idum6
                    endif
                    if(icase.eq.3.and.                                    &
                         &k2.eq.ifacd2(1,m).and.l2.eq.ifacd2(2,m)                           &
                         &.and.k4.eq.ifacd2(3,m).and.l4.eq.ifacd2(4,m)) then
                       is=is+1
                       dum1=facd2(1,is)
                       dum2=facd2(2,is)
                       dum3=ham(1,is)
                       dum4=ham(2,is)
                       idum1=ifacd2(1,is)
                       idum2=ifacd2(2,is)
                       idum3=ifacd2(3,is)
                       idum4=ifacd2(4,is)
                       facd2(1,is)=facd2(1,m)
                       facd2(2,is)=facd2(2,m)
                       ham(1,is)=ham(1,m)
                       ham(2,is)=ham(2,m)
                       ifacd2(1,is)=ifacd2(1,m)
                       ifacd2(2,is)=ifacd2(2,m)
                       ifacd2(3,is)=ifacd2(3,m)
                       ifacd2(4,is)=ifacd2(4,m)
                       facd2(1,m)=dum1
                       facd2(2,m)=dum2
                       ham(1,m)=dum3
                       ham(2,m)=dum4
                       ifacd2(1,m)=idum1
                       ifacd2(2,m)=idum2
                       ifacd2(3,m)=idum3
                       ifacd2(4,m)=idum4
                    endif
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  return
end subroutine sortres
subroutine prror(ipro,ier,num)
  implicit none
  integer ier,ipro,num
  !-----------------------------------------------------------------------
  !--Error Output
  !-----------------------------------------------------------------------

  if(ipro.eq.0)  write(6,*) ' Error occurred in Reading User Input '
  if(ipro.eq.1)  write(6,*) ' Error occurred in Program Detune, ',  &
       &' which calculates the first and second Order Detuning.'
  if(ipro.eq.2)  write(6,*) ' Error occurred in Program Distort1, ',&
       &' which calculates the Distortion Function in first Order.'
  if(ipro.eq.3)  write(6,*) ' Error occurred in Program Distort2, ',&
       &' which calculates the Distortion Function in second Order.'
  goto(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160),ier
10 write(6,*) ' File 34 empty or corrupted - Program stops'
  stop 1
20 write(6,*) ' File 34 cannot be read - Program stops'
  stop 2
30 write(6,*) ' Number of Elements nblz: ',num,' too small!'
  stop 3
40 write(6,*) ' Input File fc.34 corrupted '
  stop 4
50 write(6,*) ' Too many Cases: in Derivation mmult:',num
  stop 5
60 write(6,*) ' Too many Cases: ic terms, increase mmult: ',num
  stop 6
70 write(6,*) ' Too many Cases: id terms, increase mmultx: ',num
  stop 7
80 write(6,*) ' Number of Resonance Terms exceeded, increase ',      &
       &'mmultx: ',num
  stop 8
90 write(6,*) ' Program Error unphysical Terms in Case: ',num
  stop 9
100 write(6,*) ' Maximum Number of Resonances: ',num,' too large ',   &
       &'==> increase mmultf'
  stop 10
110 write(6,*) ' Inconsistent Number: ',num,' of Resonances for ',    &
       &'the mixed Case'
  stop 11
120 write(6,*) ' Final Ordering of Resonances is not possible as ',   &
       &'their Number: ',num,' is too large.'
  stop 12
130 write(6,*) ' There are only 2 derivatives possible - program',    &
       &' error'
  stop 13
140 write(6,*) ' Program Specifier iprog: ',num,' must lie between',  &
       &' 1 and 7 '
  stop 14
150 write(6,*) ' Analysis is restricted to order: ',num
  stop 15
160 write(6,*) ' Unit: ',num,' could not be opened '
  stop 16
end subroutine prror

subroutine fit_table(table_name,data_size,tab_size)
  implicit none
  integer data_size,tab_size,ll
  character table_name*(*)
  character*16 tab_name

  tab_name = table_name
  ll = len(tab_name)
  tab_name(ll:ll) = ' '
10 if(data_size .gt. tab_size) then
     call double_table(tab_name)
     tab_size = tab_size*2
     goto 10
  endif

  return
end subroutine fit_table

subroutine ertab(k,table_name,column,row)
  implicit none
  integer k,row
  character table_name*(*),column*(*)
  if (k .eq.0) return
  if (k .eq. -1) print *,"table ",table_name,"  does not exist"
  if (k .eq. -2) print *," in table ",table_name,"column ",column,  &
       &" does not exist"
  if (k .eq. -3) print *,"in table ",table_name,"row",row,          &
       &" does not exist"

  return
end subroutine ertab

subroutine write_table(table_name,table_type,int_to_write,        &
     &double_to_write)
  implicit none
  integer table_type,int_to_write(11)
  integer k,table_type_index(11)
  integer j70,j71,j72,j73,j74,j75,j76,j77,j78,j79
  integer tab_types_5(5),tab_types_8(8),tab_types_9(9),             &
       &tab_types_11(11)
  double precision double_to_write(11)
  character table_name*(*)
  character*16  name_5(5),name_8(8),name_9(9),name_11(11)
  common/indeces/j70,j71,j72,j73,j74,j75,j76,j77,j78,j79
  data table_type_index/0,0,0,0,1,0,0,2,3,0,4/
  data tab_types_5/1,1,2,1,1/
  data tab_types_8/1,2,2,2,1,1,1,1/
  data tab_types_9/1,1,2,2,2,1,1,1,1/
  data tab_types_11/1,1,1,2,2,2,2,1,1,1,1/
  data name_5/'multipoleorder ','plane ','detuning ','h_inv_order ',&
       &'v_inv_order '/
  data name_8/'multipoleorder ','cosine ','sine ','amplitude ',     &
       &'j ','k ','l ','m '/
  data name_9/'multipoleorder1 ','multipoleorder2 ','cosine ',      &
       &'sine ','amplitude ','j ','k ','l ','m '/
  data name_11/'multipoleorder ','location ','resonance ',          &
       &'position[m] ','cosine ','sine ','amplitude ','j ','k ','l ','m '/

  goto(50,80,90,110) table_type_index(table_type)
50 do k = 1,5
     if(tab_types_5(k) .eq. 2) then
        call double_to_table_curr(table_name,name_5(k),double_to_write(k))
     else
        call double_to_table_curr(table_name,name_5(k),                    &
             &dble(int_to_write(k)))
     endif
  enddo
  go to 1000

80 do k = 1,8
     if(tab_types_8(k) .eq. 2) then
        call double_to_table_curr(table_name,name_8(k),double_to_write(k))
     else
        call double_to_table_curr(table_name,name_8(k),                    &
             &dble(int_to_write(k)))
     endif
  enddo
  go to 1000

90 do k = 1,9
     if(tab_types_9(k) .eq. 2) then
        call double_to_table_curr(table_name,name_9(k),double_to_write(k))
     else
        call double_to_table_curr(table_name,name_9(k),                    &
             &dble(int_to_write(k)))
     endif
  enddo
  goto 1000

110 do k = 1,11
     if(tab_types_11(k) .eq. 2) then
        call double_to_table_curr(table_name,name_11(k),double_to_write(k))
     else
        call double_to_table_curr(table_name,name_11(k),                   &
             &dble(int_to_write(k)))
     endif
  enddo

1000 return
end subroutine write_table

      subroutine timex(r1)
        implicit none
        real r1,timestart,timenow
        common /mytimes/timestart
        save
        call timest(0.0)
        call cpu_time(timenow)
        r1=timenow-timestart
        return
      end subroutine timex

      subroutine timest(r1)
        implicit none
        real r1,timestart
        logical start
        common /mytimes/timestart
        data start /.false./
        save
        if (.not.start) then
           start=.true.
           call cpu_time(timestart)
        endif
        return
      end subroutine timest
