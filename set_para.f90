module set_para

contains

  SUBROUTINE set_PARAMETERS(R,nt,iorder,IFAM,inda,scale)
    !Strength of Multipole of order iorder as parameter
    USE S_TRACKING

    IMPLICIT NONE
    integer ipause, mypause
    logical(lp) ok
    INTEGER iorder,i,j,jj,k,lstr,IFAM,tot,nt,inda,min1
    INTEGER,parameter::ipara=100
    real(dp) scale(ipara),value
    character(20) str
    CHARACTER(3) STR1
    character(10),dimension(10)::multname
    type(layout) r
    type(fibre), POINTER :: current
    INTEGER,ALLOCATABLE,dimension(:)::DAFAM
    INTEGER,ALLOCATABLE,dimension(:,:)::FAM
    real(dp),ALLOCATABLE,dimension(:)::SFAM
    multname=(/"Dipole    ","Quadrupole","Sextupole ","Octupole  ","Decapole  ",&
         "Dodecapole","14-Pole   ","16-Pole   ","18-Pole   ","20-Pole   "/)

    ALLOCATE(FAM(IFAM,0:R%N),DAFAM(IFAM),SFAM(IFAM))

    min1=0
    if(iorder.lt.0) then
       min1=1
       iorder=-iorder
    endif
    DO I=1,IFAM
       OK=.TRUE.
       DO WHILE(OK)
          TOT=0
          if(min1.eq.0) WRITE(6,*) " Identify ",multname(iorder)
          if(min1.eq.1) WRITE(6,*) " Identify ","SKEW-"//multname(iorder)
          READ(5,*) STR
          STR=TRIM(ADJUSTL(STR))
          LSTR=LEN_TRIM (STR)
          current=>r%start
          DO J=1,R%N
             IF(current%MAG%NAME==STR.and.current%MAG%P%NMUL==iorder) THEN
                print*,current%MAG%P%NMUL,current%MAG%KIND
                TOT=TOT+1
                FAM(I,TOT)=J
             ENDIF
             current=>current%next
          ENDDO
          WRITE(6,*) TOT," Is that OK? YES or NO?"
          READ (5,*) STR1
          STR1=TRIM(ADJUSTL(STR1))
          IF(STR1(1:1)=='Y'.OR.STR1(1:1)=='y') THEN
             OK=.FALSE.
             inda=inda+1
             if(inda.gt.100) then
                write(6,*) " Problem: Only ",ipara," Parameters allowed"
                ipause=mypause(2002)
             endif
             DAFAM(I)=inda
             WRITE(6,*) " Give Scaling Factor, '0' uses Default"
             read(5,*) value
             if(value==0) then
                WRITE(6,*) " Take Default Scaling Value : ",scale(inda)
                SFAM(I)=scale(inda)
             else
                SFAM(I)=value
             endif
          ENDIF
       ENDDO

       FAM(I,0)=TOT
       current=r%start
       DO JJ=1,FAM(I,0)
          J=FAM(I,JJ)
          ! ALLOCATION GYMNASTIC IF Multipole NOT YET ALLOCATED
          IF(current%MAGP%P%NMUL<iorder) THEN
             CALL KILL(current%MAGP%BN,current%MAGP%P%NMUL)
             CALL KILL(current%MAGP%AN,current%MAGP%P%NMUL)
             current%MAGP%P%NMUL=iorder
             DEALLOCATE(current%MAGP%BN)
             DEALLOCATE(current%MAGP%AN)
             CALL ALLOC(current%MAGP%BN,iorder)
             CALL ALLOC(current%MAGP%AN,iorder)
             ALLOCATE(current%MAGP%BN(iorder),current%MAGP%AN(iorder))
             DO K=1,current%MAG%P%NMUL
                current%MAGP%BN(K)=current%MAG%BN(K)
                current%MAGP%AN(K)=current%MAG%AN(K)
             ENDDO
             DEALLOCATE(current%MAG%BN)
             DEALLOCATE(current%MAG%AN)
             ALLOCATE(current%MAG%BN(iorder),current%MAG%AN(iorder))
             call equal(current%MAG,current%MAGP)
          ENDIF
          if(min1.eq.0) then
             current%MAGP%BN(iorder)%I=NT+I
             current%MAGP%BN(iorder)%KIND=3
          else
             current%MAGP%AN(iorder)%I=NT+I
             current%MAGP%AN(iorder)%KIND=3
          endif
          current=>current%next
       ENDDO
    ENDDO

    current=r%start
    DO I=1,IFAM
       DO JJ=1,1
          J=FAM(I,JJ)
          if(min1.eq.0) WRITE(6,*)  current%MAG%NAME,' ', current%MAG%BN(iorder)
          if(min1.eq.1) WRITE(6,*)  current%MAG%NAME,' ', current%MAG%AN(iorder)
          current=>current%next
       ENDDO
    ENDDO

    DEALLOCATE(FAM,STAT=I)
    !    WRITE(6,*) I
    DEALLOCATE(DAFAM,STAT=I)
    !    WRITE(6,*) I
    DEALLOCATE(SFAM,STAT=I)
    !    WRITE(6,*) I

  end subroutine set_PARAMETERS

end module set_para
