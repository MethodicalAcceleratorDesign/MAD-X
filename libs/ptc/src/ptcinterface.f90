!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----These subroutines are called from C to manage PTC
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine ptc_init(p_in_file)

        USE madx_ptc_module
!        use accel_ptc
        IMPLICIT NONE

        character p_in_file*128 
        TYPE(LAYOUT), POINTER :: o_ring
        integer i,j
        real(dp) x(6)
        logical exists
        character*20 file0
  !      type(mad_universe),target:: orbit_universe
  !      M_U=>orbit_universe
        
        write (*,*) "=============PTC INIT=====START==========="
!    allocate(m_u)
!    call set_up_universe(m_u)
!    allocate(m_t)
!    call set_up_universe(m_t)
  call ptc_ini_no_append()

!old        allocate(m_u)
!old        call set_up_universe(m_u)
!        call APPEND_EMPTY_LAYOUT(m_u)
!        call READ_INTO_VIRGIN_LAYOUT(m_u%start,p_in_file,lmax=lmax)
           

          N_CAV4_F=3
          CALL  READ_AND_APPEND_VIRGIN_general(m_u,p_in_file,lmax0=lmax)
!          call READ_INTO_VIRGIN_LAYOUT(m_u%start,p_in_file,lmax=lmax)


           file0="pre_orbit_set.txt"
          INQUIRE (FILE = file0, EXIST = exists)

          if(exists) then
            call  read_ptc_command77(file0)
          endif
        o_ring=>m_U%end

        call MAKE_NODE_LAYOUT(o_ring)
        if(lmax==zero) then
         write(6,*) " Error lmax = 0 "
!         pause
         stop 777
        endif


        if(lmax>0) then
           CALL ORBIT_MAKE_NODE_LAYOUT(o_ring,my_true)
        else
           lmax=-lmax
           CALL ORBIT_MAKE_NODE_LAYOUT(o_ring,my_false)
        endif

        write(6,*) " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" 
        write(6,*) " ORBIT_USE_ORBIT_UNITS ",my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS
        write(6,*) " ORBIT_N_NODE ",my_ORBIT_LATTICE%ORBIT_N_NODE
        write(6,*) " ORBIT_WARNING ",my_ORBIT_LATTICE%ORBIT_WARNING
        write(6,*) " ORBIT_OMEGA ",my_ORBIT_LATTICE%ORBIT_OMEGA
        write(6,*) " ORBIT_HARMONIC ",my_ORBIT_LATTICE%ORBIT_HARMONIC
        write(6,*) " ORBIT_P0C ",my_ORBIT_LATTICE%ORBIT_P0C
        write(6,*) " ORBIT_KINETIC ",my_ORBIT_LATTICE%ORBIT_KINETIC
        write(6,*) " ORBIT_MASS_IN_AMU ",my_ORBIT_LATTICE%ORBIT_MASS_IN_AMU
        write(6,*) " ORBIT_L ",my_ORBIT_LATTICE%ORBIT_L
        write(6,*) " ORBIT_GAMMAT ",my_ORBIT_LATTICE%ORBIT_GAMMAT
        write(6,*) " ORBIT_CHARGE ",my_ORBIT_LATTICE%ORBIT_CHARGE
        write(6,*) " ORBIT_LMAX AND LMAX ",my_ORBIT_LATTICE%ORBIT_LMAX,LMAX
        write(6,*) " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" 
!        MY_STATE=>my_orbit_lattice%state
!        MY_ORBIT_STATE=>MY_STATE
CAVITY_TOTALPATH=0

      return
      end  

!===================================================
!reads additional ptc commands from file and execute them inside ptc
!===================================================
      subroutine ptc_script(p_in_file)

        USE orbit_ptc
        USE pointer_lattice
        IMPLICIT NONE

        character p_in_file*128 

        CALL read_ptc_command(p_in_file)

     end  subroutine ptc_script

!===================================================
!get initial twiss at entrance of the ring
!===================================================
      subroutine ptc_get_twiss_init(bx,by,ax,ay,ex,epx)

        USE orbit_ptc
        IMPLICIT NONE
        REAL(DP) bx,by,ax,ay,ex,epx

        INTEGER i 
        REAL(DP) length 

        call GET_N_NODE(i)
        call GET_info(i,length,bx,by,ax,ay,ex,epx)

      return
      end


!===================================================
!get number of PTC ORBIT nodes, harmonic number, 
!     length of the ring, and gamma transition
!===================================================
      subroutine ptc_get_ini_params(nNodes,nHarm,lRing,gammaT)

        USE orbit_ptc
        IMPLICIT NONE
        REAL(DP) lRing,gammaT
        INTEGER nNodes,nHarm

        call GET_HARMONIC(nHarm)
        call GET_CIRCUMFERENCE(lRing)
        call GET_N_NODE(nNodes)
        call GET_GAMMAT(gammaT)

      return
      end

!===================================================
!get synchronous particle params mass, charge, and energy
!===================================================
      subroutine ptc_get_syncpart(mass,charge,kin_energy)

        USE orbit_ptc
        IMPLICIT NONE
        REAL(DP) mass, kin_energy
!        real(dp) charge
        integer charge

        call GET_MASS_AMU(mass)
        
        call GET_kinetic(kin_energy)
        call GET_CHARGE(charge)

      return
      end

!===================================================
!get twiss and length for a node with index
!===================================================
      subroutine ptc_get_twiss_for_node(node_index, node_length, bx,by,ax,ay,ex,epx)

        USE orbit_ptc
        IMPLICIT NONE
        REAL(DP) node_length, bx,by,ax,ay,ex,epx
        INTEGER node_index

        INTEGER i
        i = node_index + 1
        call GET_info(i,node_length,bx,by,ax,ay,ex,epx)

      return
      end


!===================================================
!track 6D coordinates through the PTC-ORBIT node
!===================================================
!       subroutine ptc_track_particle(node_index, x,xp,y,yp,phi,dE)
!
!         USE orbit_ptc
!         IMPLICIT NONE
!         REAL(DP) x,xp,y,yp,phi,dE
!         INTEGER node_index
!         INTEGER i
!        
!         i = node_index + 1
!        
!         call PUT_RAY(x,xp,y,yp,phi,dE)
!         call TRACK_ONE_NODE(i)
!         call GET_RAY(x,xp,y,yp,phi,dE)
!
!       return
!       end subroutine orbit_ptc_track_particle
!
!===========================================================
! This subroutine should be called before particle tracking.
!  It specifies the type of the task that will be performed
!  in ORBIT before particle tracking for particular node.
!  i_task = 0 - do not do anything
!  i_task = 1 - energy of the sync. particle changed
!===========================================================
   SUBROUTINE ptc_get_task_type(i_node,i_task)
     
    USE orbit_ptc
    IMPLICIT NONE
    INTEGER  i_node,i_task
    INTEGER  i_node1
    
    i_node1 = i_node + 1
    call GET_task(i_node1,i_task)

   END SUBROUTINE  ptc_get_task_type

!===================================================
!It returns the lowest frequency of the RF cavities
!This will be used to transform phi to z[m] 
!===================================================
   SUBROUTINE ptc_get_omega(x)
     
     USE orbit_ptc
     IMPLICIT NONE
     REAL(DP) x
     call GET_omega(x)

   END SUBROUTINE ptc_get_omega

!===================================================
!It reads the acceleration table into the ptc code
!===================================================
   subroutine ptc_read_accel_table(p_in_file)

     IMPLICIT NONE
     character p_in_file*128 
     write(6,*) " This is just an ordinary Script now! "
     
     call read_ptc_command77(p_in_file)
       
   END SUBROUTINE ptc_read_accel_table

SUBROUTINE ptc_synchronous_set(i_node)
     IMPLICIT NONE
     integer i_node
 !  write(6,*) " Not needed anymore "
end SUBROUTINE ptc_synchronous_set

SUBROUTINE ptc_synchronous_after(i_node)
     use pointer_lattice
     IMPLICIT NONE
     integer i_node
   
   if(i_node==-100) then
      write(6,*) " ********************************************************* "
      write(6,*) "  "
      write(6,*) " Orbit State is Printed "
      if(accelerate) then
       write(6,*) "Acceleration is turned ON "
      else
       write(6,*) "Acceleration is turned OFF "
      endif
      write(6,*) "  "
      call print(my_ORBIT_LATTICE%state,6)
      write(6,*) " ********************************************************* "
   else
   ! write(6,*) " ptc_synchronous_after not needed anymore "
   endif
end SUBROUTINE ptc_synchronous_after

