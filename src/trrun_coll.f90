subroutine trcoll(flag, apx, apy, turn, sum, part_id, last_turn,  &
     last_pos, last_orbit, z, ntrk,al_errors,offx,offy)

  use twiss0fi
  use name_lenfi
  implicit none

  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   test for collimator aperture limits.                               *
  ! input:                                                               *
  !   flag      (integer)   aperture type flag:                          *
  !                         1: elliptic, 2: rectangular                  *
  !   apx       (double)    x aperture or half axis                      *
  !   apy       (double)    y aperture or half axis                      *
  !   turn      (integer)   current turn number.                         *
  !   sum       (double)    accumulated length.                          *
  ! input/output:                                                        *
  !   part_id   (int array) particle identification list                 *
  !   last_turn (int array) storage for number of last turn              *
  !   last_pos  (dp. array) storage for last position (= sum)            *
  !   last_orbit(dp. array) storage for last orbit                       *
  !   z(6,*)    (double)    track coordinates: (x, px, y, py, t, pt).    *
  !   ntrk      (integer) number of surviving tracks.                    *
  !----------------------------------------------------------------------*
  integer flag,turn,part_id(*),last_turn(*),ntrk,i,n,nn
  double precision apx,apy,sum,last_pos(*),last_orbit(6,*),z(6,*),  &
       one,al_errors(align_max),offx,offy
  parameter(one=1d0)
  character(name_len) aptype

  n = 1
10 continue
  do i = n, ntrk

     !---- Is particle outside aperture?
     if (flag .eq. 1.and.((z(1,i)-al_errors(11)- offx)/apx)**2             &
          +((z(3,i)-al_errors(12)- offy) / apy)**2 .gt. one) then
        go to 99
     else if(flag .eq. 2                                             &
          .and. (abs(z(1,i)-al_errors(11)- offx) .gt. apx                         &
          .or. abs(z(3,i)-al_errors(12)- offy) .gt. apy)) then
        go to 99
        !***  Introduction of marguerite : two ellipses
     else if(flag .eq. 3.and. ((z(1,i)-al_errors(11)- offx) / apx)**2      &
          + ((z(3,i)-al_errors(12)- offy) / apy)**2 .gt. one .and.                &
          ((z(1,i)-al_errors(11)- offx) / apy)**2 +                               &
          ((z(3,i)-al_errors(12)- offy) / apx)**2 .gt. one) then
        go to 99
     endif
     go to 98
99   n = i
     nn=name_len
     call node_string('apertype ',aptype,nn)
     call trkill(n, turn, sum, ntrk, part_id,                        &
          last_turn, last_pos, last_orbit, z, aptype)
     goto 10
98   continue
  enddo
end subroutine trcoll

subroutine trcoll1(flag, apx, apy, turn, sum, part_id, last_turn,  &
     last_pos, last_orbit, z, ntrk,al_errors, apr,offx,offy)

  use twiss0fi
  use name_lenfi
  implicit none

  !----------------------------------------------------------------------*
  ! Similar with trcoll, for racetrack type aperture
  !-------------Racetrack type , added by Yipeng SUN 21-10-2008---
  !----------------------------------------------------------------------*
  integer flag,turn,part_id(*),last_turn(*),ntrk,i,n,nn
  double precision apx,apy,sum,last_pos(*),last_orbit(6,*),z(6,*),  &
       one,al_errors(align_max),apr,offx,offy
  parameter(one=1d0)
  character(name_len) aptype

  n = 1
10 continue
  do i = n, ntrk

     !---- Is particle outside aperture?
     if (flag .eq. 4                                                 &
          .and. (abs(z(1,i)-al_errors(11)- offx)) .gt. (apr+apx)                  &
          .or. abs(z(3,i)-al_errors(12)- offy) .gt. (apy+apr) .or.                &
          ((((abs(z(1,i)-al_errors(11)- offx)-apx)**2+                            &
          (abs(z(3,i)-al_errors(12)- offy)-apy)**2) .gt. apr**2)                  &
          .and. (abs(z(1,i)-al_errors(11)- offx)) .gt. apx                        &
          .and. abs(z(3,i)-al_errors(12)- offy) .gt. apy)) then
        go to 99
     endif
     go to 98
99   n = i
     nn=name_len
     call node_string('apertype ',aptype,nn)
     call trkill(n, turn, sum, ntrk, part_id,                        &
          last_turn, last_pos, last_orbit, z, aptype)
     goto 10
98   continue
  enddo
end subroutine trcoll1
subroutine trkill(n, turn, sum, jmax, part_id,                    &
     last_turn, last_pos, last_orbit, z, aptype)

  use name_lenfi
  implicit none

  !hbu--- kill particle:  print, modify part_id list
  logical recloss
  integer i,j,n,turn,part_id(*),jmax,last_turn(*),get_option
  double precision sum, z(6,*), last_pos(*), last_orbit(6,*),       &
       torb(6) !, theta, node_value, st, ct, tmp
  character(name_len) aptype
  !hbu
  character(name_len) el_name

  recloss = get_option('recloss ') .ne. 0

  !!--- As elements might have a tilt we have to transform back
  !!--- into the original coordinate system!
  !      theta = node_value('tilt ')
  !      if (theta .ne. 0.d0)  then
  !          st = sin(theta)
  !          ct = cos(theta)
  !!--- rotate trajectory (at exit)
  !            tmp = z(1,n)
  !            z(1,n) = ct * tmp - st * z(3,n)
  !            z(3,n) = ct * z(3,n) + st * tmp
  !            tmp = z(2,n)
  !            z(2,n) = ct * tmp - st * z(4,n)
  !            z(4,n) = ct * z(4,n) + st * tmp
  !      endif

  last_turn(part_id(n)) = turn
  last_pos(part_id(n)) = sum
  do j = 1, 6
     torb(j) = z(j,n)
     last_orbit(j,part_id(n)) = z(j,n)
  enddo

  !hbu
  call element_name(el_name,len(el_name))
  !hbu
  write(6,'(''particle #'',i6,'' lost turn '',i6,''  at pos. s ='', f10.2,'' element='',a,'' aperture ='',a)') &
       part_id(n),turn,sum,el_name,aptype
  print *,"   X=",z(1,n),"  Y=",z(3,n),"  T=",z(5,n)

  if(recloss) then
     call tt_ploss(part_id(n),turn,sum,torb,el_name)
  endif

  do i = n+1, jmax
     part_id(i-1) = part_id(i)
     do j = 1, 6
        z(j,i-1) = z(j,i)
     enddo
  enddo
  jmax = jmax - 1

end subroutine trkill
subroutine tt_ploss(npart,turn,spos,orbit,el_name)

  use name_lenfi
  implicit none

  !hbu added spos
  !----------------------------------------------------------------------*
  !--- purpose: enter lost particle coordinates in table                 *
  !    input:                                                            *
  !    npart  (int)           particle number                            *
  !    turn   (int)           turn number                                *
  !    spos    (double)       s-coordinate when loss happens             *
  !    orbit  (double array)  particle orbit                             *
  !    orbit0 (double array)  reference orbit                            *
  !----------------------------------------------------------------------*
  integer npart,turn,j
  double precision orbit(6),tmp,tt,tn
  double precision energy,get_value
  character(36) table
  character(name_len) el_name
  !hbu
  double precision spos
  !hbu
  character(4) vec_names(7)
  !hbu
  data vec_names / 'x', 'px', 'y', 'py', 't', 'pt', 's' /
  data table / 'trackloss' /

  tn = npart
  tt = turn

  energy = get_value('probe ','energy ')


  ! the number of the current particle
  call double_to_table_curr(table, 'number ', tn)
  ! the number of the current turn
  call double_to_table_curr(table, 'turn ', tt)
  !hbu spos

  call double_to_table_curr(table,vec_names(7),spos)

  do j = 1, 6
     tmp = orbit(j)
     call double_to_table_curr(table, vec_names(j), tmp)
  enddo
  call double_to_table_curr(table, 'e ', energy)
  call string_to_table_curr(table, 'element ', el_name)

  call augment_count(table)
end subroutine tt_ploss
