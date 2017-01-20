module madx_ptc_setcavs_module
  USE madx_keywords
  USE madx_ptc_intstate_module
  implicit none
  public

  ! flag for debugging ranges from 0 (no debug printout) to 10 (the most detailed)
  logical(lp), public                     :: cavsareset   = .false.
  ! flag that indicates if cavities were already set for the current setup

  !********************************************************************************************
  !********************************************************************************************
  !********************************************************************************************

contains
  subroutine setcavities(my_ring, maxaccel)
    implicit                none
    type(layout),target  :: my_ring
    type(internal_state) :: localis ! internal state to be use in this routine = intstate+totalpath+time
    logical(lp)          :: maxaccel
    real(dp)             :: charge    ! charge of an accelerated particle
    !      use madx_keywords
    integer              :: i,j!,currentelement=1       !iterators
    integer              :: mf1,mf2
    type(fibre), pointer :: p     ! skowron: temporary variable: current fibre
    type(work)           :: startfen
    type(work)           :: nfen      ! New Fibre ENergy
    integer, pointer     :: poscav(:) !array keeping indexes of cavities
    real(dp),allocatable :: phasecav(:) !array keeping phases of cavities
    real(dp)             :: patchprecision=c_1d_8
    logical(lp)          :: patchenergy=.true.
    logical(lp)          :: patchnext=.true.
    real(dp)             :: prevbeta0  !just a temporary real variable
    real (dp)            :: x(1:6)   ! track vector -
    ! here we always use closed orbit track, that is all its relative coordinates are 0
    real(dp)             :: position=zero !synchronous particle position
    real(kind(1d0))     :: get_value
    !------------------------------------------------------
100 format (a20, f10.4, a10, f10.4, a10, f10.4)
110 format (8f10.4, l2, i3)
120 format (i16, f16.3, f16.4, e16.3, e16.3, f16.3)
130 format (a12, i5, a6, a30, a6, f10.4, a20, f10.4)
    !------------------------------------------------------


    !Below we enforce that x(6) is cT, and it is time of flight from the start
    !we use time T=x(6)/ctime to find the time of arrival to a cavity so we can adjust its phase optimally
    
    !This is not needed, and even to contrary, it creates a bug: 
    !if the tracking is done without totalpath, cavities should be tuned to such x(6)=0 gives max acceleration

    localis = getintstate()
    
!    localis = localis - nocavity0 + totalpath0
    localis = localis - delta0 - only_4d0 - nocavity0 + totalpath0
    if (getdebug() > 1) then
       print *, "I am in setcavities "
       call print(localis,6)
    endif

    patchnext=.true.

    charge = get_value('beam ', "charge ")

    if (getdebug() > 1) then
       call kanalnummer(mf1)
       open(unit=mf1,file='sychrpart.txt')
       call kanalnummer(mf2)
       open(unit=mf2,file='twcavsettings.txt')
       write(mf2,'(6a16)') "!ElNo     ","Ref.Momentum","Phase","Frequency [Hz]","Voltage","DeltaE"
    endif

    nfen = 0
    startfen = 0
    x(:)=zero

    call locate_all_twcav(my_ring,poscav)
    if ( getdebug() > 2 ) write(6,*) "There are ", size(poscav), " Cavities in the line."
    if ( size(poscav) == 0) then
       if (getdebug() > 1) then
          close(mf1);close(mf2);
       endif
       return
    endif

    allocate(phasecav(size(poscav)))

    !  Here is tracking element by element
    p=>my_ring%start

    startfen=p  !setting up start energy for record
    nfen=p      ! current fibre energy

    if ( getdebug() > 1 ) then
        print *, 'c_%feed_p0c = ', c_%feed_p0c
    endif

    if ( getdebug() > 2 ) write (*,*) 'START TRACKING TILL THE FIRST CAVITY'

    i = 1

    do j=1,size(poscav)

       if ( getdebug() > 2 ) then
          write (*,*) 'Current cavity no is j=',j
          write (*,*) 'Setting beam momentum AND tracking ', nfen%p0c,' till the cavity (',poscav(j),')'
       endif


       do i=i,poscav(j)-1 !from the current element i to the current cavity j

          p = nfen   ! set current reference energy
          call track(my_ring,x,i,i+1,localis)

          if ( .not. c_%stable_da) then
             call fort_warn('setcavities: ','DA got unstable')
             call seterrorflag(10,"setcavities ","DA got unstable");

             deallocate(poscav);
             deallocate(phasecav);
             if (getdebug() > 1) then
                close(mf1);close(mf2);
             endif
             return
          endif

          if ( getdebug()>1 ) then
             write (6,*) ' i=',i,' name=',p%mag%name, &
                  ' beta0 ', nfen%beta0, &
                  ' newpos ', x(6), &
                  ' Current energy ',nfen%energy
             write(6,'(6f8.4)') x
          endif

          if (getdebug() > 1) then
             write (mf1,*) ' '
             write (mf1,130) 'i=',i,' name=',p%mag%name,' p0c=',p%mag%p%p0c, ' Current energy ',nfen%energy
             write (mf1,'(6f8.4)') x
          endif

          p=>p%next
       enddo

       if (associated(p%next)) then
          if ( (p%next%mag%kind==kind21) .or. (p%next%mag%kind==kind4) ) then
             write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
             write(6,*) "!!!                                     !!!"
             write(6,*) "!!!    CONESCUTIVE        CAVITIES      !!!"
             write(6,*) "!!!                                     !!!"
             write(6,*) "!!!     Currently it is forbiden        !!!"
             write(6,*) "!!! to place one cavity after another.  !!!"
             write(6,*) "!!! Plese insert a marker between them. !!!"
             write(6,*) "!!!                                     !!!"
             write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
             call aafail('setcavities:','Found consecutive cavities; insert a marker. Program stops')
          endif
       endif

       !AT THIS POINT WE HAVE ARRIVED TO A CAVITY
       ! p point to this cavity

       ! set the reference energy in this cavity

       p = nfen

       if ( getdebug() > 1 ) then
          write (6,130) 'i=',i,' name=',p%mag%name,' p0c=',p%mag%p%p0c, ' Current energy ',nfen%energy
          write (6,'(6f8.4)') x
       endif

       !TUNE CAVITY
       call setcavity(p,x,phasecav(j),charge,maxaccel)

       if (getdebug() > 1) then
          write(mf2,120) poscav(j), p%mag%p%p0c, p%mag%phas*c_360/twopi, p%mag%freq, p%mag%volt, p%mag%delta_e
       endif

       !TRACK CAVITY
       call track(my_ring,x,poscav(j),poscav(j)+1,localis)

       if ( .not. c_%stable_da) then
          call fort_warn('setcavities: ','DA got unstable')
          call seterrorflag(10,"setcavities ","DA got unstable");

          deallocate(poscav);
          deallocate(phasecav);
          if (getdebug() > 1) then
             close(mf1);close(mf2);
          endif
          return
       endif

       if (getdebug() > 1) then
          write (mf1,*) ' '
          write (mf1,130) 'poscav(j)=',poscav(j),' name=',p%mag%name,' p0c=',p%mag%p%p0c, ' Current energy ',nfen%energy
          write (mf1,'(6f8.4)') x
       endif

       if ( getdebug() > 2 ) then
          write(6,'(a, 6f12.8)') ' Track parameters after cavity ',x
          write(*,100) 'Old Fibre: energy=',nfen%energy,' momentum=',nfen%p0c,' kinetic=',nfen%kinetic
       endif

       ! GET NEW ENERGY AFTER THE CAVITY
       prevbeta0 = nfen%beta0

       nfen= x(5)*nfen%p0c

       if ( getdebug() > 2 ) then
          write(6,100) 'New Fibre: energy=',nfen%energy,' momentum=',nfen%p0c,' kinetic=',nfen%kinetic
          write(6,110) nfen
          write(6,'(a10, f8.4)') 'Relative E increase', (nfen%energy-startfen%energy)/startfen%kinetic
       endif

       if (getdebug()>1) write (6,*) 'beta0 ', prevbeta0, ' oldpos ', position, ' newpos ', x(6), ' Current energy ',nfen%energy
       position = x(6)

       !PATCH THE NEXT ELEMENT ON ENTRANCE
       p%next = nfen

       if ( getdebug() > 2 ) write (*,*) 'Finding patch for j=',j,' ',p%mag%name
       call find_patch(p,next=patchnext,ENERGY_PATCH=patchenergy,PREC=patchprecision)

       i=poscav(j)+1
       p=>p%next

       !from this point on we do not need to calculate TOF cause there is no further cavs to set
    enddo

    if ( getdebug() > 2 ) then
       write (*,*) 'Loop over cavities done'
       write (*,*) 'Current element is ', p%mag%name
       write (*,*) 'Doing loop from the first element after the last cavity to the END'
    endif

    do i=i,my_ring%n !setting beam energies to the end of line
       p = nfen
       call track(my_ring,x,i,i+1,localis)

       if ( .not. c_%stable_da) then
          call fort_warn('setcavities: ','DA got unstable')
          call seterrorflag(10,"setcavities ","DA got unstable");

          deallocate(poscav);
          deallocate(phasecav);
          if (getdebug() > 1) then
             close(mf1);close(mf2);
          endif
          return
       endif

       if (getdebug() > 1) then
          write (mf1,*) ' '
          write (mf1,130) 'i=',i,' name=',p%mag%name,' p0c=',p%mag%p%p0c, ' Current energy ',nfen%energy
          write (mf1,'(6f8.4)') x
       endif

       if ( getdebug() > 1 ) then
          write(6,*) ' i=',i,' name=',p%mag%name, &
               ' beta0 ', nfen%beta0, &
               ' newpos ', x(6), &
               ' Current energy ',nfen%energy
          write(6,'(6f8.4)') x
       endif

       p=>p%next
    enddo


    if (getdebug() > 1) then
       write (mf1,*) ' '
       write (mf1,*) 'END'
       write (mf1,'(6f8.4)') x

       write(6,*) 'PARAMETERS AT THE END OF LINE:'
       write(6,'(a, 6f8.4)') ' Track parameters ',x
       write(*,100) 'START energy=',startfen%energy,' momentum=',startfen%p0c,' kinetic=',startfen%kinetic
       write(6,100) 'END energy=',nfen%energy,' momentum=',nfen%p0c,' kinetic=',nfen%kinetic
       write(6,110) nfen
       write(6,'(a10, f8.4)') 'Relative E increase', (nfen%energy-startfen%energy)/startfen%kinetic
    endif

    ! FINISHED HERE

    patchnext=.false.
    p=>my_ring%start

    do i=1,my_ring%n
       if ( associated(p%mag) .eqv. .false.) then
          if (getdebug() > 1 ) then
              print *, 'Fibre no. ',i,' has no mag assigned to it'
          endif
          cycle
       endif
       if ( getdebug() > 2 ) then
          write(6,*) 'Name: ', p%mag%name, ' Kind: ', p%mag%kind
       endif

       if( (p%mag%kind == kind21) .or. (p%mag%kind == kind4)) then

          if(p%next%patch%energy==1) then
             p%patch%energy=2
             p%next%patch%energy=0
          endif

          if ( getdebug() > 1 ) then
             write (6,*) 'Cavity ',i,' name ',p%mag%name,' phase ', p%mag%phas,' Volt ',p%mag%volt, &
                  & ' length ', p%mag%l
             write(6,*) 'DELTAE ', p%mag%DELTA_E
          endif


       endif
       p=>p%next
    enddo

    cavsareset = .true. !module field indicating that cavities were set appriopriately
    deallocate(poscav);
    deallocate(phasecav);

    if (getdebug() > 1) then
       close(mf1);close(mf2);
    endif

    !****************************************************************************************
    !*********  E N D   O F   PTC_TRACKCAVS  ************************************************
    !****************************************************************************************
    !________________________________________________________________________________________

  contains  ! what follows are internal subroutines of ptc_trackline
    !____________________________________________________________________________________________

    subroutine setcavity(f, x, phase_rel, charge, ene)
      implicit none
      type(fibre),intent(inout):: f         ! fiber -> here must be a cavity, i.e. kind21 (tw) or kind4 (rf)
      real(dp)                 :: x(6)      ! reference particle coordinates (closed orbit for a circular machine)
      real(dp)                 :: phase_rel ! final relative phase
      real(dp)                 :: phase_old ! for printing old phase
      real(dp)                 :: charge    ! charge of an particle
      logical(lp)              :: ene       ! switches if cavity should always maximally accelerate
      ! the reference track; lag is calculated
      !      logical(lp)              :: givendene = .false. ! makes cavity always accelerate about a given value;
      !      integer(4)               :: tmp
      ! volt is calculated; lag and freq is preserved
      real(dp)                 :: de_mev ! delta energy
      real(dp)                 :: arrivtime !time of arrival
!      type(internal_state)     :: globalis ! internal state to be use in the tracking
      
      
      
      arrivtime = x(6)/clight
      
      if (getdebug()>2) then
          print *, 'arrivtime = ', arrivtime*1e9, ' ns'
      endif

      if( (f%mag%kind/=kind21) .and. (f%mag%kind/=kind4) ) then
         write(6,*) " fatal error: not a Cavity "
         call aafail('setcavity', 'fatal error: not a Cavity. Program stops')
      endif

      if ( f%mag%kind==kind21) then
         if(f%mag%cav21%psi/=zero) then
            write(6,*) " warning: backwards wave present ",f%mag%cav21%psi
            f%mag%cav21%psi=zero   ! removing backward waves
            f%magp%cav21%psi=zero   ! removing backward waves
         endif
      endif

!      globalis = getintstate()
!      print*, "Total path is ", globalis%totalpath
!      chargesign = charge/abs(charge)
      
      phase_old = f%mag%phas
      
      if(ene) then

         if ( getdebug() > 2 ) then
            print*,"MAX ACCEL MODE"
            de_mev=f%mag%volt*f%mag%l
            write(*,*) '   Max Energy to gain: ', de_mev, ' MeV, x(6)', x(6)
         endif
         ! nunmber of radians from the lauch of the synchr particle to its arrival to the cavity
         f%mag%phas = pi/two - twopi*f%mag%freq*arrivtime - f%mag%lag ! here we tune to be on the crest and then we add the lag
         f%magp%phas= f%mag%phas
         phase_rel=f%mag%phas

      else
         if ( getdebug() > 2 ) then
           print*,"REGULAR MODE, not max accel"
         endif  
         
         if( f%mag%kind == kind21  ) then
         
           f%mag%phas = - f%mag%lag
           f%magp%phas= f%mag%phas
           phase_rel=f%mag%phas
         
         else !(f%mag%kind==kind4) 
          
           f%mag%phas = - f%mag%lag
           f%magp%phas= f%mag%phas
           phase_rel=f%mag%phas
         
         endif

      endif


      !    write (*,*) 'energy (t/f)? :',ene, 'charge: ', charge
      if ( getdebug() > 1 ) then
         write(6,*) 'Cavity settings:'
         write(6,*)                    '    Name   ', f%mag%name
         write(6,'(a12,f12.5,a10,l1)') '    Charge ', charge,' max ene? : ',ene
         write(6,'(a12,f12.5,a10)')    '    Volt ',   f%mag%volt,' MV '
         write(6,'(a12,f12.5,a10)')    '    DELTAE ', f%mag%delta_e, ' GeV '
         write(6,'(a12,f12.5,a10)')    '    Length ', f%mag%l,' m'
         write(6,'(a12,f12.3,a10)')    '    Phase ',  f%mag%phas, ' rad '
         write(6,'(a12,f12.3,a10)')    '    Old Ph ', phase_old, ' rad '
         write(6,'(a12,f12.0,a10)')    '    Freq ',   f%mag%freq, ' Hz '
         write(6,'(a12,f12.5,a10,f12.4,a10)') '    Lag ',    f%mag%lag/twopi*c_360,' deg ', f%mag%lag,' rad '
         write(6,'(a12,f12.5,a10)')    '    P0c ',    f%mag%p%p0c, 'GeV/c'
      endif

    end subroutine setcavity
    !____________________________________________________________________________________________

    subroutine locate_all_twcav(r,pos)
      implicit none
      type(layout), target, intent(inout) :: r
      type(fibre), pointer:: p
      integer, pointer ::  pos(:)
      integer i,ic
      ic=0
      p=>r%start
      do i=1,r%n
         if( (p%mag%kind==kind21) .or. (p%mag%kind==kind4) ) then
            if(p%mag%freq/=zero) then
               ic=ic+1
            endif
         endif
         p=>p%next
      enddo
      allocate(pos(ic))
      pos=0
      ic=0
      p=>r%start
      do i=1,r%n
         if( (p%mag%kind==kind21) .or. (p%mag%kind==kind4) ) then
            if(p%mag%freq/=zero) then
               ic=ic+1
               pos(ic)=i
            endif
         endif
         p=>p%next
      enddo

    end subroutine locate_all_twcav

  end subroutine setcavities

end module madx_ptc_setcavs_module


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         if (f%mag%delta_e > 0) then
!            kf= -f%mag%delta_e/f%mag%l/f%dir/charge
!            if(givendene) then
!               !push about given energy, calculates the needed voltage keeping lag and freq
!               f%mag%volt=kf
!               f%magp%volt=kf
!            else
!               !uses simply the peak voltage and phase
!               f%mag%volt=sign(one,kf*f%mag%volt) * f%mag%volt
!               f%magp%volt=f%mag%volt
!            endif
!            f%mag%phas=pi/two-twopi*f%mag%freq*arrivtime-c_%phase0
!            f%magp%phas=f%mag%phas
!         endif
!         phase_rel=f%mag%phas+twopi*f%mag%freq*arrivtime

!       tmp = f%mag%phas/twopi
!       print *, tmp, f%mag%phas/twopi
!       f%mag%phas = f%mag%phas - twopi*dble(tmp)
!       f%magp%phas=f%mag%phas
!       print *, f%mag%phas
