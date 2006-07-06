module madx_ptc_eplacement_module
  !This module enables the user to PLACE an element in arbitrary position 
  use madx_keywords
  use madx_ptc_intstate_module, only : getdebug
  use madx_ptc_module, only: my_ring;
  implicit none
  save
  private

  !============================================================================================
  !  PUBLIC INTERFACE
  public                                      :: place_element
  public                                      :: printlayout_rootm

  !============================================================================================
  !  PRIVATE
  !    data structures
  integer, parameter                          :: black  = 2
  integer, parameter                          :: red    = 2
  integer, parameter                          :: green  = 3
  integer, parameter                          :: blue   = 4
  integer, parameter                          :: yellow = 5
  integer, parameter                          :: magenta= 6
  integer, parameter                          :: cyan = 7
  integer, parameter                          :: darkgreen = 8
  integer, parameter                          :: dgrey = 14
  integer, parameter                          :: grey = 16
  integer, parameter                          :: lgrey = 18
  !    routines
  private                                     :: rot
  private                                     :: printfframes
  !============================================================================================

contains
  !____________________________________________________________________________________________

  subroutine place_element(elementidx,refframe)
    implicit none
    integer              :: elementidx
    integer              :: refframe
    integer              :: j, i
    integer              :: restart_sequ,advance_node
    real(kind(1d0))      :: get_value
    type(fibre), pointer :: p
    type(fibre), pointer :: mvf !moved fibre
    real(dp)             :: ent(3,3),a(3)
    logical              :: onlyposition, onlyorientation, autoplace
    
    ent(:,:) = 0
    ent(1,1) = 1
    ent(2,2) = 1
    ent(3,3) = 1
    
    a(:) = 0
    nullify(mvf)
    
    j=restart_sequ()
    j=0
    p=>my_ring%start
    
    elementidx = elementidx + 1
    if (getdebug()>2) then
      print*, "I am in placeelement: Element index is ", elementidx
      print*, "refframe is", refframe
    endif  
    
    if ( (elementidx .lt. 1) .and. (elementidx .gt. my_ring%n) ) then
      call fort_warn("place_element","element out of range of the current layout")
      return
    endif
    
    do 

      j=j+1
      
      if (elementidx == j) then
        exit
      endif
      
      if(advance_node().eq.0) then
        return
      endif

      p=>p%next
    
    enddo

    if (getdebug() > 1 ) then
       print*,"Found element no. ", elementidx," named ", p%mag%name, &
             &" of kind ", p%mag%kind, mytype(p%mag%kind)
    endif
    
    onlyorientation = get_value('ptc_eplacement ','onlyorientation ') .ne. 0
    if (onlyorientation .eqv. .false.) then
       a(1) = get_value('ptc_eplacement ','x ')
       a(2) = get_value('ptc_eplacement ','y ')
       a(3) = get_value('ptc_eplacement ','z ')
       write(6,'(a20, 3(f12.10,1x))') 'Read position ',a
       
       if (refframe ==  1) then
         !reference frame of the passed parameters is the current magnet position
         if (getdebug() > 2) then
           print*,"Reference frame: the current magnet position"
         endif

         a = rot(P%chart%f%ent,a)
   
         write(6,'(a65, 3f8.4)') 'Position tranformed from the current magnet ref frame to global',a
         
         a(1) = a(1) +  P%chart%f%a(1)
         a(2) = a(2) +  P%chart%f%a(2)
         a(3) = a(3) +  P%chart%f%a(3)
         
         
       endif

       if (refframe ==  2) then
         !reference frame of the passed parameters is the end face of the preceding magnet
         if (getdebug() > 2) then
           print*,"Reference frame: the end face of the preceding magnet"
         endif

         a = rot(P%previous%chart%f%exi,a)
         
         a(1) = a(1) +  P%previous%chart%f%b(1)
         a(2) = a(2) +  P%previous%chart%f%b(2)
         a(3) = a(3) +  P%previous%chart%f%b(3)
         
       endif
       
    else
       a =  P%chart%f%a !copy current position
    endif
        
    onlyposition = get_value('ptc_eplacement ','onlyposition ') .ne. 0
    if (onlyposition .eqv. .false.) then
       call readrotmatrix() !reads ent matrix

       if (refframe ==  1) then
         ent = rotm(P%chart%f%ent,ent)
       endif
       
       if (refframe ==  2) then
         ent = rotm(P%previous%chart%f%exi,ent)
       endif
       
    else
       ent = P%chart%f%exi
    endif
    
    
    call survey(p,ENT,a)
    

    autoplace = get_value('ptc_eplacement ','autoplacedownstream ') .ne. 0

    if (autoplace) then
      call survey(my_ring,j,my_ring%n)
    else
      CALL FIND_PATCH(P%PREVIOUS,P,NEXT=my_true,ENERGY_PATCH=MY_FALSE)
      CALL FIND_PATCH(P,P%NEXT,NEXT=MY_FALSE,ENERGY_PATCH=MY_FALSE)
    endif  


    mvf=>p


    do i=j,my_ring%n-1
       CALL FIND_PATCH(P,P%next,NEXT=MY_TRUE,ENERGY_PATCH=MY_FALSE)
       P=>P%NEXT
    ENDDO
!    CALL FIND_PATCH(P,P%next,NEXT=my_false,ENERGY_PATCH=MY_FALSE)
    
!    print*, "Exiting placeelement"
  contains 
    subroutine readrotmatrix
      implicit none
      real(dp)             :: phi, theta
      real(dp)             :: sinphi, cosphi, sintheta, costheta
       phi   = get_value('ptc_eplacement ','phi ')
       theta = get_value('ptc_eplacement ','theta ')

       write(6,'(a20, 2f8.4)') 'Read rotations ',phi, theta
       sinphi = sin(phi)
       cosphi = cos(phi)
       sintheta = sin(theta)
       costheta = cos(theta)

       ent(1,1) = cosphi 
       ent(1,2) = 0
       ent(1,3) = sinphi

       ent(2,1) =  sintheta*sinphi
       ent(2,2) =  costheta
       ent(2,3) = -sintheta*cosphi

       ent(3,1) = -costheta*sinphi
       ent(3,2) =  sintheta
       ent(3,3) =  costheta*cosphi

       write(6,*) " Rotation MATRIX  "
       write(6,'(3f8.4)') ent(1,:)
       write(6,'(3f8.4)') ent(2,:)
       write(6,'(3f8.4)') ent(3,:)
    end subroutine readrotmatrix
  end subroutine place_element


  !____________________________________________________________________________________________
  !____________________________________________________________________________________________
  !____________________________________________________________________________________________

  subroutine printlayout_rootm(filenameIA)
    implicit none
    include 'twissa.fi'
    integer   filenameIA(*)
    integer       :: i  !iterator + tmp
    integer       :: mf !macro file descriptor
    integer       :: xmin, xmax, ymin, ymax
    character(48) :: filename
    character(48) :: fctname
    TYPE(LAYOUT),pointer :: r
    type(fibre), pointer :: p
    real(dp)     :: z, a(3)
    
    
    r=>my_ring
    
    filename = charconv(filenameIA)
    
!    print*, "I am in printlayout_rootm: Macro name is is ", filename

    call kanalnummer(mf)
    open(unit=mf,file=filename)
    
    i=index(filename,'.C')
    fctname = filename(1:i-1)
    
    if (getdebug() > 2 ) then
      print*, ".C found at ",i," function name is ", fctname
    endif  

    write(mf,*) "void ", fctname,'()'
    write(mf,*) "{"
    write(mf,*) 'TBRIK* l;'
    write(mf,*) 'TShape* s;'
    write(mf,*) 'TNode* mn;'
    write(mf,*) 'TNode* n;'
    write(mf,*) 'TPoints3DABC* pts;'
    write(mf,*) 'Double_t rotmatrix[9];'
    write(mf,*) 'TRotMatrix* m;'
    
    write(mf,*) 'TCanvas* c = new TCanvas("c","PTC Layout",10,10,800,600);'

    write(mf,*) 's = new TBRIK("START","START","void",0.01,0.01,0.01);'
    write(mf,*) 's->SetLineColor(2);'
    write(mf,*) 'mn = new TNode("NODE1","NODE1","START");'
    write(mf,*) 'mn->cd();'


    
    
    xmin = 0
    xmax = 0
    ymin = 0
    ymax = 0
    p=>r%start
    do i=1,r%n
      
      if ( P%mag%p%f%a(1) .lt. xmin) xmin = P%mag%p%f%a(1)
      if ( P%mag%p%f%b(1) .lt. xmin) xmin = P%mag%p%f%b(1)
      if ( P%mag%p%f%a(3) .gt. ymax) ymax = P%mag%p%f%a(3)
      if ( P%mag%p%f%b(3) .gt. ymax) ymax = P%mag%p%f%b(3)
      
      P=>P%NEXT
    enddo
    xmin = xmin -1
    ymin = ymin -1

    xmax = xmax +1
    ymax = ymax +1
    
    write(mf,*) 'c->Range(', xmin, ',', ymin, ',', xmax, ',', ymax,');'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    p=>r%start
    do i=1,r%n
      write(mf,*) 'cout<<"',p%mag%name,'"<<endl;'
      if (getdebug() > 2) then
        print*, i,p%mag%name
        print*, 'Edges: ', p%mag%P%EDGE(1), p%mag%P%EDGE(2)
      endif  
      
      if (p%mag%l == 0) goto 100;
      
      select case(p%mag%kind)
       
       case(kind1)
         call drawtube(p,mf,0.05_dp,lgrey)

       case(kind11)
         call drawtube(p,mf,0.05_dp,magenta)

       case(kind10)
         call drawsbend(p,mf)

       case(kind20)
         call drawboxm(p,mf,magenta)

       case(kind21)
         call drawtube(p,mf,0.25_dp,yellow)
         
       case(kind16)
         if (getdebug() > 3) then 
           print*, "KIND16: likemad is ", p%mag%k16%likemad
           print*, "KIND16: bn(0) ", p%mag%bn(0), " bn(1)", p%mag%bn(1), " bn(2)", p%mag%bn(2)
           print*, "KIND16: an(0) ", p%mag%an(0), " an(1)", p%mag%an(1), " an(2)", p%mag%an(2)
         endif
         if (p%mag%bn(1) /= zero ) then
           a(1)= P%mag%p%f%ent(3,1)*p%mag%l/two + P%mag%p%f%a(1)
           a(2)= P%mag%p%f%ent(3,2)*p%mag%l/two + P%mag%p%f%a(2)
           a(3)= P%mag%p%f%ent(3,3)*p%mag%l/two + P%mag%p%f%a(3)
           call drawbox(p,mf,P%mag%p%f%ent,a,blue)
         else    
           if (p%mag%bn(2) .gt. zero ) then
             call drawboxm(p,mf,red)  !QUAD foc
           else
             call drawboxm(p,mf,green)!QUAD defoc
           endif 
         endif
         
         
       case default
         print*, "################# "
         print*, "################# "
         print*, "Unrecognized kind ", p%mag%kind, mytype(p%mag%kind)
         print*, "################# "
         print*, "################# "
      end select
      
      
      
100   continue      
      P=>P%NEXT
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

    write(mf,*) 'mn->Draw();'
    write(mf,*) 'c->Update();'
    write(mf,*) 'c->GetViewer3D("ogl");'

    write(mf,*) "}"
    close(mf)

  contains 

    subroutine  drawsbend(p,mf)
      implicit none
      type(fibre), pointer :: p
      integer      :: mf !macro file descriptor
      integer      :: color = blue
      character(9) :: fname
      character(8) :: nodname
      character(7) :: mtxname
      real(dp)     :: x,y,z,r, xx,yy, phi
      !To be finished - need more information about SBEND frames in PTC
      
      
      if (P%mag%p%f%a(2) .ne. P%mag%p%f%b(2) ) then
       print*, "Not able yet to drow horizonthally skewed SBEND. DRAWING AS RBEND"
       call drawboxm(p,mf,color)
       return
      endif

      x = P%mag%p%f%b(1) - P%mag%p%f%a(1)
      y = P%mag%p%f%b(3) - P%mag%p%f%a(3)

      z = (x*x + y*y)/four !square of the half of distance between ent and exi

      x = P%mag%p%f%a(1) + x/two !x coord of the half of distance between ent and exi
      y = P%mag%p%f%a(3) + y/two !z coord of the half of distance between ent and exi

      !!!!!!!!!!!!!!!!!!!
      
      x = x - P%mag%p%f%o(1)
      y = y - P%mag%p%f%o(3)
      
      y = sqrt(x*x + y*y)

      if (y == 0) then
       print*, "All three reference frames are inline. DRAWING AS RBEND"
       call drawboxm(p,mf,color)
       return
      endif

      r = (z + y*y)/(two*y)
      
      phi = two*arcsin(z/r)
!      print*, "R is ", r," phi is ", phi," z is ", z," y ",y
      
      write(fname,'(a5,i4.4)') 'SBEND',i
      write(mf,*) 's = new TTUBS("',fname, '","',  fname,'","void",',&
                      &  r-0.25_dp,',',r+0.25_dp,',0.25,0,',phi,');'
      write(mf,*) 's->SetLineColor(',color,');'
    
      write(mtxname,'(a3,i4.4)') 'mtx',i
      call setmatrix(P%mag%p%f%mid,mtxname,mf)
      
      write(nodname,'(a4,i4.4)') 'NODE',i
      write(mf,*) 'n = new TNode("',nodname,'","',  nodname,'",s,', &
                                 & x,',',y, ',',z,',m);'
    
    end subroutine  drawsbend

    subroutine  drawbox(p,mf,m,a,color)
      implicit none
      type(fibre), pointer :: p
      integer      :: mf !macro file descriptor
      real(dp)     :: m(3,3)
      real(dp)     :: a(3)
      integer      :: color
      character(9) :: fname
      character(8) :: nodname
      character(7) :: mtxname
      real(dp)     :: x,y,z

      write(fname,'(a5,i4.4)') 'RECTA',i
      write(mf,*) 's = new TBRIK("',fname, '","',  fname,'","void",0.5,0.5,',p%mag%l/2_dp,');'
      write(mf,*) 's->SetLineColor(',color,');'

      x = a(1)
      y = a(2)
      z = a(3)
      

      write(mtxname,'(a3,i4.4)') 'mtx',i
      call setmatrix(m,mtxname,mf)
      
      write(nodname,'(a4,i4.4)') 'NODE',i
      write(mf,*) 'n = new TNode("',nodname,'","',  nodname,'",s,', &
                                 & x,',',y, ',',z,',m);'
      
    end subroutine  drawbox  

    subroutine  drawboxm(p,mf,color)
      implicit none
      type(fibre), pointer :: p
      integer      :: mf !macro file descriptor
      integer      :: color

      call drawbox(p,mf, P%mag%p%f%mid, P%mag%p%f%o, color)
      
    end subroutine  drawboxm
    
    
    subroutine  drawtube(p,mf,r,color)
      implicit none
      type(fibre), pointer :: p
      integer      :: mf !macro file descriptor
      real(dp)     :: r
      integer      :: color
      character(9) :: fname
      character(8) :: nodname
      character(7) :: mtxname
      real(dp)     :: x,y,z
  
      write(fname,'(a5,i4.4)') 'DRIFT',i
      write(mf,*) 's = new TTUBE("',fname, '","',  fname,'","void",',r,',',p%mag%l/2_dp,');'
      write(mf,*) 's->SetLineColor(',color,');'
      
      x = P%mag%p%f%o(1)
      y = P%mag%p%f%o(2)
      z = P%mag%p%f%o(3)
      

      write(mtxname,'(a3,i4.4)') 'mtx',i
      call setmatrix(P%mag%p%f%mid,mtxname,mf)
      
      write(nodname,'(a4,i4.4)') 'NODE',i
      write(mf,*) 'n = new TNode("',nodname,'","',  nodname,'",s,', &
                                 & x,',',y, ',',z,',m);'
      
      
      
    end subroutine  drawtube

    subroutine setmatrix(m,name,mf)
      implicit none
      real(dp)     :: m(3,3)
      character(6) :: name
      integer      :: mf !macro file descriptor
      
      
      write(mf,*) 'rotmatrix[0] = ', m(1,1),';'
      write(mf,*) 'rotmatrix[1] = ', m(1,2),';'
      write(mf,*) 'rotmatrix[2] = ', m(1,3),';'

      write(mf,*) 'rotmatrix[3] = ', m(2,1),';'
      write(mf,*) 'rotmatrix[4] = ', m(2,2),';'
      write(mf,*) 'rotmatrix[5] = ', m(2,3),';'

      write(mf,*) 'rotmatrix[6] = ', m(3,1),';'
      write(mf,*) 'rotmatrix[7] = ', m(3,2),';'
      write(mf,*) 'rotmatrix[8] = ', m(3,3),';'
      write(mf,*) 'm = new TRotMatrix("',name,'","',name,'",rotmatrix);'
      
    end subroutine setmatrix
    

  end subroutine printlayout_rootm

  subroutine printfframes(p)
    implicit none
    type(fibre), pointer :: p
    write(6,*) " Actual magnet positioning  "
    write(6,*) "  "
    write(6,*) " Entrance origin A(3) "
    write(6,*) P%mag%p%f%a
    write(6,*) " Entrance frame (i,j,k) basis in the ent(3,3) array "
    write(6,*) P%mag%p%f%ent(1,:)
    write(6,*) P%mag%p%f%ent(2,:)
    write(6,*) P%mag%p%f%ent(3,:)
    write(6,*) " Middle origin O(3) "
    write(6,*) P%mag%p%f%o
    write(6,*) " Middle frame (i,j,k) basis in the ent(3,3) array "
    write(6,*) P%mag%p%f%mid(1,:)
    write(6,*) P%mag%p%f%mid(2,:)
    write(6,*) P%mag%p%f%mid(3,:)
    write(6,*) " Exit origin B(3) "
    write(6,*) P%mag%p%f%B
    write(6,*) " Exit frame (i,j,k) basis in the ent(3,3) array "
    write(6,*) P%mag%p%f%exi(1,:)
    write(6,*) P%mag%p%f%exi(2,:)
    write(6,*) P%mag%p%f%exi(3,:)
  end subroutine printfframes

  function rot(m,a)
    implicit none
    real(dp)             :: rot(3)
    real(dp)             :: a(3),m(3,3)
    
    rot(1) = m(1,1)*a(1) + m(2,1)*a(2) + m(3,1)*a(3)
    rot(2) = m(1,2)*a(1) + m(2,2)*a(2) + m(3,2)*a(3)
    rot(3) = m(1,3)*a(1) + m(2,3)*a(2) + m(3,3)*a(3)

    if (getdebug() > 3) then
      write(6,*) " Vector Rotation  "
      write(6,'(a26, f8.4)')             " ", a(1)
      write(6,'(a26, f8.4)')             " ", a(2)
      write(6,'(a26, f8.4)')             " ", a(3)
      write(6,*) "     x "
      write(6,'(3f8.4,a2,f8.4)') m(:,1), " ", rot(1)
      write(6,'(3f8.4,a2,f8.4)') m(:,2), " ", rot(2)
      write(6,'(3f8.4,a2,f8.4)') m(:,3), " ", rot(3)
    endif
        
  end function   
  !____________________________________________________________________________________________

  function rotm(a,m)
    implicit none
    real(dp)             :: rotm(3,3)
    real(dp)             :: a(3,3),m(3,3)
    
    rotm(1,1) = a(1,1)*m(1,1) + a(1,2)*m(2,1) + a(1,3)*m(3,1)
    rotm(2,1) = a(2,1)*m(1,1) + a(2,2)*m(2,1) + a(2,3)*m(3,1)
    rotm(3,1) = a(3,1)*m(1,1) + a(3,2)*m(2,1) + a(3,3)*m(3,1)

    rotm(1,2) = a(1,1)*m(1,2) + a(1,2)*m(2,2) + a(1,3)*m(3,2)
    rotm(2,2) = a(2,1)*m(1,2) + a(2,2)*m(2,2) + a(2,3)*m(3,2)
    rotm(3,2) = a(3,1)*m(1,2) + a(3,2)*m(2,2) + a(3,3)*m(3,2)

    rotm(1,3) = a(1,1)*m(1,3) + a(1,2)*m(2,3) + a(1,3)*m(3,3)
    rotm(2,3) = a(2,1)*m(1,3) + a(2,2)*m(2,3) + a(2,3)*m(3,3)
    rotm(3,3) = a(3,1)*m(1,3) + a(3,2)*m(2,3) + a(3,3)*m(3,3)

    
    if (getdebug() > 3) then
      write(6,*) " 2 Rotation MATRIX product "
      write(6,'(a26, 3f8.4)')             " ", a(:,1)
      write(6,'(a26, 3f8.4)')             " ", a(:,2)
      write(6,'(a26, 3f8.4)')             " ", a(:,3)
      write(6,*) "     x "
      write(6,'(3f8.4,a2,3f8.4)') m(:,1), " ", rotm(:,1)
      write(6,'(3f8.4,a2,3f8.4)') m(:,2), " ", rotm(:,2)
      write(6,'(3f8.4,a2,3f8.4)') m(:,3), " ", rotm(:,3)
    endif
    
  end function 
  
  
end module madx_ptc_eplacement_module
