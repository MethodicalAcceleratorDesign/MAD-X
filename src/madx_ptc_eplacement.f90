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
  integer, parameter                          :: black  = 1
  integer, parameter                          :: red    = 2
  integer, parameter                          :: green  = 3
  integer, parameter                          :: blue   = 4
  integer, parameter                          :: yellow = 5
  integer, parameter                          :: magenta= 6
  integer, parameter                          :: cyan = 7
  integer, parameter                          :: darkgreen = 8
  integer, parameter                          :: violet = 9
  integer, parameter                          :: color_n_sext = 46
  integer, parameter                          :: color_s_sext = 30
  integer, parameter                          :: dgrey = 14
  integer, parameter                          :: grey = 16
  integer, parameter                          :: lgrey = 18
  integer, parameter                          :: color_of_ghost = 19

  !    routines
  private                                     :: rot
  private                                     :: printfframes
  
  
  external :: fort_warn, read_ptc_command77, gino_ptc_command77, gettrack
  
  
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
    real(dp),target      :: a(3), ange(3)
    real(dp)             :: phi, theta, psi
    real(dp),target      :: idm(3,3),ent(3,3)
    real(dp),pointer     :: base(:,:)
    logical              :: onlyposition, onlyorientation, autoplace, surveyall

    idm(:,:) = zero
    idm(1,1) = one
    idm(2,2) = one
    idm(3,3) = one

    a(:) = 0

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
       if (getdebug() > 2 ) then
          write(6,'(a, 3(f13.10,1x))') 'ptc_eplacement: Read position ',a
       endif
    endif

    select case(refframe)
    case(0)
       if (getdebug() > 2) then
          print*,"ptc_eplacement: Reference frame: Global Coordinate System"
       endif
       a = a - p%chart%f%a
       base => idm
    case(1)
       !reference frame of the passed parameters is the current magnet position
       if (getdebug() > 2) then
          print*,"ptc_eplacement: Reference frame: the current magnet position"
       endif

       ent = p%chart%f%ent  !we need to make copy because it may be overwritten at the time of rotation
       base => ent
    case(2)
       !reference frame of the passed parameters is the end face of the preceding magnet
       if (getdebug() > 2) then
          print*,"ptc_eplacement: Reference frame: the end face of the preceding magnet"
       endif
       base => P%previous%chart%f%exi
    case default
       refframe = 0
       call fort_warn("ptc_eplacement","Such reference frame is not supported. Using global")
       a = a - p%chart%f%a
       base => idm
    end select




    onlyposition = get_value('ptc_eplacement ','onlyposition ') .ne. 0
    if (onlyposition .eqv. .false.) then

       phi   = get_value('ptc_eplacement ','phi ')
       theta = get_value('ptc_eplacement ','theta ')
       psi   = get_value('ptc_eplacement ','psi ')
       if (getdebug() > 2 ) then
          write(6,'(a20, 2f8.4)') 'Read rotations ',phi, theta
       endif


       if (refframe /=  1) then
          !brings it to the ref system
          CALL COMPUTE_ENTRANCE_ANGLE(p%chart%f%ent,base,ANGE)

          if (getdebug() > 2 ) then
             write(6,'(a, 3(f13.10,1x))') 'ptc_eplacement: R0 Computed entr angles ',ange
          endif

          CALL ROTATE_Fibre(p,p%chart%f%a,ange)

       endif

       ange(1) = theta
       ange(2) = phi
       ange(3) = psi

       CALL ROTATE_Fibre(p,p%chart%f%a,ange,1,base)

    endif

    if (onlyorientation .eqv. .false.) then
       if (refframe ==  2) then ! face translation of the moved one to the end of the previous one
          CALL TRANSLATE_Fibre(p,p%previous%chart%f%b - p%chart%f%a)
       endif
       CALL TRANSLATE_Fibre(p,a,1,base)
    endif


    p%patch=-1
    p%patch= 0
    CALL FIND_PATCH(P%PREVIOUS,P,NEXT=my_true,ENERGY_PATCH=MY_FALSE)

    autoplace = get_value('ptc_eplacement ','autoplacedownstream ') .ne. 0
    if (autoplace) then
       if (getdebug() > 2) then
          print*,"ptc_eplacement: autoplacedownstream=true: running survey"
       endif
       call survey(my_ring,j+1,my_ring%n)
    else
       if (getdebug() > 2) then
          print*,"ptc_eplacement: autoplacedownstream=false: finding patch"
       endif
       CALL FIND_PATCH(P,P%NEXT,NEXT=MY_FALSE,ENERGY_PATCH=MY_FALSE)
    endif


    j=j+1
    p=>p%next

    do i=j,my_ring%n-1

       p%patch =-1
       p%patch = 0
       CALL FIND_PATCH(P,P%next,NEXT=MY_TRUE,ENERGY_PATCH=MY_FALSE)
       P=>P%NEXT
    ENDDO

    surveyall = get_value('ptc_eplacement ','surveyall ') .ne. 0
    if (surveyall) then
       if (getdebug() > 2) then
          print*,"ptc_eplacement: surveyall=true"
       endif
       call survey(my_ring)
    endif

    !    CALL FIND_PATCH(P,P%next,NEXT=my_false,ENERGY_PATCH=MY_FALSE)

    !    print*, "Exiting placeelement"
  end subroutine place_element


  !____________________________________________________________________________________________
  !____________________________________________________________________________________________
  !____________________________________________________________________________________________

  subroutine printlayout_rootm(filenameIA)
    use twissafi
    implicit none
    integer   filenameIA(*)
    integer       :: i  !iterator + tmp
    integer       :: mf !macro file descriptor
    integer       :: xmin, xmax, ymin, ymax
    character(48) :: filename
    character(48) :: fctname
    TYPE(LAYOUT),pointer :: r
    type(fibre), pointer :: p
    real(dp)     :: a(3)
    integer      :: nmul
    character(48) charconv

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

    write(mf,*) '#ifndef __MAKECINT__'
    write(mf,*) ' #ifndef __CINT__'
    write(mf,*) ''
    write(mf,*) ' #include "TROOT.h"'
    write(mf,*) ' #include "TCanvas.h"'
    write(mf,*) ' #include "Riostream.h"'
    write(mf,*) ''
    write(mf,*) ' #include "TBRIK.h"'
    write(mf,*) ' #include "TShape.h"'
    write(mf,*) ' #include "TNode.h"'
    write(mf,*) ' #include "TCanvas.h"'
    write(mf,*) ' #include "TGLViewer.h"'
    write(mf,*) ' #include "TPoints3DABC.h"'
    write(mf,*) ' #include "TTUBE.h"'
    write(mf,*) ' #include "TRotMatrix.h"'
    write(mf,*) ''
    write(mf,*) ' #endif'
    write(mf,*) '#endif'
    write(mf,*) ''
    write(mf,*) ''
    write(mf,*) "void ", fctname,'()'
    write(mf,*) "{"
    write(mf,*) 'TShape* s;'
    write(mf,*) 'TNode* mn;'
    write(mf,*) 'TNode* n;'
    write(mf,*) 'Double_t rotmatrix[9];'
    write(mf,*) 'TRotMatrix* m;'
    write(mf,*) ''
    write(mf,*) 'gSystem->Load("libGed");'
    write(mf,*) 'gSystem->Load("libRGL");'
    write(mf,*) ''
    write(mf,*) 'TCanvas* c = new TCanvas("c","PTC Layout",10,10,800,600);'
    write(mf,*) ''
    write(mf,*) 's = new TBRIK("START","START","void",0.01,0.01,0.01);'
    write(mf,*) 's->SetLineColor(2);'
    write(mf,*) 'mn = new TNode("NODE1","NODE1","START");'
    write(mf,*) 'mn->cd();'
    write(mf,*) ''




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

    write(mf,*) ''
    write(mf,*) 'c->Range(', xmin, ',', ymin, ',', xmax, ',', ymax,');'
    write(mf,*) ''
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    p=>r%start
    do i=1,r%n
       write(mf,*)
       write(mf,*) '//cout<<',i,'<<" ',p%mag%name,'"<<endl;'

       !print*, i,p%mag%name
       if (getdebug() > 2) then
          print*, i,p%mag%name,' of kind ',p%mag%kind
          print*, 'Edges: ', p%mag%P%EDGE(1), p%mag%P%EDGE(2)
       endif


       if (p%mag%l == zero) then
          !print*, i,p%mag%name(1:3), p%mag%kind
          if(p%mag%name(1:3) == 'M__') then
            call drawtube(p,mf,1.0_dp,lgrey)
          endif

          goto 100;
       endif

       select case(p%mag%kind)

       case(kind1)
          call drawtube(p,mf,0.05_dp,lgrey)

       case(kind11) !monitor
          call drawtube(p,mf,0.05_dp,magenta)

       case(kind12) !H monitor
          call drawtube(p,mf,0.05_dp,magenta)

       case(kind13) !V monitor
          call drawtube(p,mf,0.05_dp,magenta)

       case(kind14) !instrument
          call drawtube(p,mf,0.05_dp,magenta)

       case(kind10) !SBEND
          call drawsbend(p,mf)

       case(32) !SBEND (Multipole) with model 1 ( Kick Drift Kick) and NON-exact
          !print*,"This is type 32"
          call drawsbend(p,mf)

       case(kind20) !stright exact bend
          call drawboxm(p,mf,magenta)

       case(kind15)!septum (electric)
          call drawboxm(p,mf,cyan)

       case(kind21) !TW CAV
          call drawtube(p,mf,0.25_dp,yellow)

       case(kind4) !TW CAV
          call drawtube(p,mf,0.25_dp,yellow) !TW CAV

       case(kind5) !solenoid
          call drawtube(p,mf,0.25_dp,darkgreen)

       case(kind16,kind7)
          if ( p%mag%kind==kind16 .and. getdebug() > 3) then
             print*, "KIND16: likemad is ", p%mag%k16%likemad
             print*, "KIND16: bn(0) ", p%mag%bn(0), " bn(1)", p%mag%bn(1), " bn(2)", p%mag%bn(2)
             print*, "KIND16: an(0) ", p%mag%an(0), " an(1)", p%mag%an(1), " an(2)", p%mag%an(2)
          endif

          nmul = p%mag%p%nmul
          if (p%mag%bn(1) /= zero ) then
             !BEND
             a(1)= P%mag%p%f%ent(3,1)*p%mag%l/two + P%mag%p%f%a(1)
             a(2)= P%mag%p%f%ent(3,2)*p%mag%l/two + P%mag%p%f%a(2)
             a(3)= P%mag%p%f%ent(3,3)*p%mag%l/two + P%mag%p%f%a(3)
             call drawbox(p,mf,P%mag%p%f%ent,a,blue)
          elseif (nmul < 2)  then
             !not powered element or very high order multipole
             call drawboxm(p,mf,color_of_ghost)
             cycle

          elseif (   (p%mag%bn(2) /= zero) ) then
             !QUAD
             if (p%mag%bn(2) .gt. zero ) then
                call drawboxm(p,mf,red)  !QUAD foc
             else
                call drawboxm(p,mf,green)!QUAD defoc
             endif

          elseif (nmul < 3)  then
             !not powered element or very high order multipole
             call drawboxm(p,mf,color_of_ghost)
             goto 100

          elseif ( (p%mag%bn(3) /= zero) .or. (p%mag%an(3) /= zero)  ) then
             !SEXTUPOLE
             call drawboxm(p,mf,color_n_sext)

          elseif (nmul < 4)  then
             !not powered element or very high order multipole
             call drawboxm(p,mf,color_of_ghost)
             goto 100
          elseif ( (p%mag%bn(4) /= zero) .or. (p%mag%an(4) /= zero) ) then

             !OCTUPOLE
             print*,"OCTUPOLE ",p%mag%bn(4),p%mag%an(4)
             call drawboxm(p,mf,color_s_sext)

          else
             !not powered element or very high order multipole
             call drawboxm(p,mf,color_of_ghost)
          endif


       case default
          print*, "################# "
          print*, "################# "
          print*, "Unrecognized kind ", p%mag%kind, mytype(p%mag%kind)
          print*, "################# "
          print*, "################# "
       end select



100    continue
       P=>P%NEXT
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(mf,*) ''
    write(mf,*) 'mn->Draw("ogl");'
    write(mf,*) 'TGLViewer * v = (TGLViewer *)c->GetViewer3D();'

    write(mf,*) ''
    write(mf,*) "}"
    close(mf)

  contains

    subroutine  drawsbend(p,mf)
      implicit none
      type(fibre), pointer :: p
      integer      :: mf !macro file descriptor
      integer      :: color = blue
      character(10) :: fname
      character(9) :: nodname
      character(8) :: mtxname
      real(dp)     :: x,y,z,r, phi
      !To be finished - need more information about SBEND frames in PTC


      if (P%mag%p%f%a(2) .ne. P%mag%p%f%b(2) ) then
         print*, "Not able yet to draw horizonthally skewed SBEND. DRAWING AS RBEND"
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
         print*,i,p%mag%name, "All three reference frames are inline. DRAWING AS RBEND"
         call drawboxm(p,mf,color)
         return
      endif

      r = (z + y*y)/(two*y)
      if (r > 100000) then
         print*,i,p%mag%name, "SBEND curvature is almost null. DRAWING AS RBEND"
         call drawboxm(p,mf,color)
         return
      endif

      phi = two*arcsin(z/r)
      print*, "R is ", r," phi is ", phi," z is ", z," y ",y

      write(fname,'(a5,i5.5)') 'SBEND',i
      write(mf,*) 's = new TTUBS("',fname, '","',  fname,'","void",',&
           &  r-0.25_dp,',',r+0.25_dp,',0.25,0,',phi,');'
      write(mf,*) 's->SetLineColor(',color,');'

      write(mtxname,'(a3,i5.5)') 'mtx',i
      !      print*, fname, " ",mtxname

      call setmatrix(P%mag%p%f%mid,mtxname,mf)

      write(nodname,'(a4,i5.5)') 'NODE',i
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
      character(10) :: fname
      character(9) :: nodname
      character(8) :: mtxname
      real(dp)     :: x,y,z

      write(fname,'(a5,i5.5)') 'RECTA',i
      write(mf,*) 's = new TBRIK("',fname, '","',  fname,'","void",0.5,0.5,',p%mag%l/two,');'
      write(mf,*) 's->SetLineColor(',color,');'

      x = a(1)
      y = a(2)
      z = a(3)


      write(mtxname,'(a3,i5.5)') 'mtx',i
      !      print*, fname," ",mtxname

      call setmatrix(m,mtxname,mf)

      write(nodname,'(a4,i5.5)') 'NODE',i
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
      character(10) :: fname
      character(9) :: nodname
      character(8) :: mtxname
      real(dp)     :: x,y,z

      write(fname,'(a5,i5.5)') 'DRIFT',i
      write(mf,*) 's = new TTUBE("',fname, '","',  fname,'","void",',r,',',p%mag%l/two,');'
      write(mf,*) 's->SetLineColor(',color,');'

      x = P%mag%p%f%o(1)
      y = P%mag%p%f%o(2)
      z = P%mag%p%f%o(3)


      write(mtxname,'(a3,i5.5)') 'mtx',i
      !      print*, fname," ", mtxname

      call setmatrix(P%mag%p%f%mid,mtxname,mf)

      write(nodname,'(a4,i5.5)') 'NODE',i
      write(mf,*) 'n = new TNode("',nodname,'","',  nodname,'",s,', &
           & x,',',y, ',',z,',m);'



    end subroutine  drawtube

    subroutine setmatrix(m,name,mf)
      implicit none
      real(dp)     :: m(3,3)
      character(8) :: name
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

  end function rot
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

  end function rotm
  !____________________________________________________________________________________________
! nowhere used
!   subroutine rotmatrixfromeuler(phi, theta, psi, ent)
!     implicit none
!     real(dp)             :: phi, theta, psi
!     real(dp)             :: ent(3,3)
!     real(dp)             :: sinphi, cosphi, sintheta, costheta, sinpsi, cospsi
! 
!     psi = zero ! we do not support rotations around the magnet axis yet
! 
! 
!     sinphi = sin(phi)
!     cosphi = cos(phi)
!     sintheta = sin(theta)
!     costheta = cos(theta)
!     sinpsi = sin(psi)
!     cospsi = cos(psi)
! 
!     ent(1,1) = cosphi*cospsi    ! cosphi
!     ent(1,2) = cospsi*sinpsi    ! 0
!     ent(1,3) = sinphi           ! sinphi
! 
!     ent(2,1) =  sintheta*sinphi*cospsi  - sinpsi*costheta   !-?         ! sintheta*sinphi
!     ent(2,2) =  costheta*cospsi         + sinpsi*sintheta*sinphi ! -?   ! costheta 
!     ent(2,3) = -sintheta*cosphi                                         !-sintheta*cosphi
! 
!     ent(3,1) = -costheta*sinphi*cospsi  + sinpsi*costheta    !-?        !-costheta*sinphi  
!     ent(3,2) =  sintheta*cospsi         + sinpsi*costheta*sinphi  !-?   ! sintheta
!     ent(3,3) =  costheta*cosphi                                         ! costheta*cosphi
! 
!     write(6,*) " Rotation MATRIX  "
!     write(6,'(3f8.4)') ent(1,:)
!     write(6,'(3f8.4)') ent(2,:)
!     write(6,'(3f8.4)') ent(3,:)
!   end subroutine rotmatrixfromeuler
! 

end module madx_ptc_eplacement_module
!____________________________________________________________________________________________
!____________________________________________________________________________________________
!____________________________________________________________________________________________
!____________________________________________________________________________________________


function w_ptc_getnelements()
  use madx_ptc_module, only: my_ring;
  implicit none
  integer ::  w_ptc_getnelements
  integer :: n

  n = 0
  if (associated(my_ring)) then
     n = my_ring%n
  endif

  w_ptc_getnelements = n

end function w_ptc_getnelements
!____________________________________________________________________________________________

subroutine w_ptc_getelname(n)
  use madx_ptc_module, only: my_ring;
  implicit none
  integer n

  if (associated(my_ring)) then
     n = my_ring%n
  endif

end subroutine w_ptc_getelname
