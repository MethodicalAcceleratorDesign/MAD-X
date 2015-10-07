subroutine res_index(skew,mynorder,myn1,myn2,indexa,mynres)
  use resindexfi
  implicit none
  logical skew
  integer mynorder,myn1,myn2,indexa(mnres,4),mynres
  integer i,j,n1,n2,n3,n4
  !---------------------------------------------------------------------
  
  !print*, skew,mynorder,myn1,myn2,mynres

  do i=1,mnres
     do j=1,4
        indexa(i,j)=0
     enddo
  enddo
  if(mynorder.le.0.or.mynorder.gt.mymorder) then
     print*," Routine res_index cannot fill index array!"
     print*," mynorder variable out of range (1 -- ",mymorder,")"
     return
  endif
  n1=mynorder
  n2=0
  n3=0
  n4=0
  mynres=0
  call myindex(skew,mynres,myn1,myn2,n1,n2,n3,n4,indexa)
  if(mynres.gt.mnres) goto 999
1 n1=n1-1
  if(n1.lt.0) return
  n2=mynorder-n1
  if(n2.gt.n1) then
     n2=n1
  elseif(n2.eq.n1.and.n1+n2.eq.mynorder) then
     n2=n2-1
  else
     call myindex(skew,mynres,myn1,myn2,n1,n2,n3,n4,indexa)
     if(mynres.gt.mnres) goto 999
     n2=n2-1
  endif
2 continue
  n3=mynorder-n1-n2
  n4=0
  call myindex(skew,mynres,myn1,myn2,n1,n2,n3,n4,indexa)
  if(mynres.gt.mnres) goto 999
  if(n3.gt.0) then
3    n3=n3-1
     n4=n4+1
     if((n1.eq.0.or.n1.eq.n2).and.(n3.le.n4)) then
        if(n2.eq.0) then
           n3=0
           n4=0
           goto 1
        else
           n2=n2-1
           goto 2
        endif
     endif
     call myindex(skew,mynres,myn1,myn2,n1,n2,n3,n4,indexa)
     if(mynres.gt.mnres) goto 999
     if(n3.gt.0) then
        goto 3
     else if(n2.eq.0) then
        n3=0
        n4=0
        goto 1
     else
        n2=n2-1
        goto 2
     endif
  endif
999 continue
  return
end subroutine res_index

subroutine myindex(skew,mynres,myn1,myn2,n1,n2,n3,n4,indexa)
  use resindexfi
  implicit none
  logical odd
  logical skew
  integer mynres,myn1,myn2,n1,n2,n3,n4,indexa(mnres,4)
  integer no,nd1,nd2
  odd=.false.
  no=n1+n2+n3+n4
  if(mod(no,2).ne.0) odd=.true.
  nd1=n1-n2
  nd2=n3-n4
  if(.not.skew) then
     if(odd.and.(nd1.eq.0.or.mod(nd1,2).eq.0)) return
     if(.not.odd.and.mod(nd1,2).ne.0) return
  else
     if(odd.and.mod(nd1,2).ne.0) return
     if(.not.odd.and.(nd1.eq.0.or.mod(nd1,2).eq.0)) return
  endif
  if((myn1.eq.0.and.myn2.eq.0).or.((myn1.eq.nd1).and.(myn2.eq.nd2))) then
     mynres=mynres+1
     if(mynres.gt.mnres) then
        print*," Maximum number: ",mnres," of resonance too small"
        return
     endif
     indexa(mynres,1)=n1
     indexa(mynres,2)=n2
     indexa(mynres,3)=n3
     indexa(mynres,4)=n4
  endif
  return
end subroutine myindex
