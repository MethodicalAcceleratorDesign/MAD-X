subroutine tthacdip(track,ktrack,turn)

  implicit none

  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Track a set of trajectories through a thin h ac dipole (zero l.)   *
  ! Input/output:                                                        *
  !   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
  !   KTRACK    (integer) number of surviving tracks.                    *
  !----------------------------------------------------------------------*
  ! Added by Yipeng SUN on 11 Nov 2009                                   *
  !----------------------------------------------------------------------*
  integer itrack,ktrack,turn,turn1,turn2,turn3,turn4
  double precision bi2gi2,dl,el,omega,dtbyds,phirf,pt,rff,betas,    &
       gammas,rfl,rfv,track(6,*),clight,twopi,vrf,deltap,deltas,pc,pc0,  &
       get_variable,node_value,get_value,one,two,half,ten3m,ten6p,px
  !      double precision px,py,ttt,beti,el1
  parameter(one=1d0,two=2d0,half=5d-1,ten3m=1d-3,ten6p=1d6)

  !---- Initialize
  clight=get_variable('clight ')
  twopi=get_variable('twopi ')

  !---- Fetch data.
  rfv = node_value('volt ')
  rff = node_value('freq ')
  rfl = node_value('lag ')
  pc0 = get_value('beam ','pc ')

  turn1 = node_value('ramp1 ')
  turn2 = node_value('ramp2 ')
  turn3 = node_value('ramp3 ')
  turn4 = node_value('ramp4 ')
  !---- Set up.
  omega = rff * twopi
  vrf   = 300 * rfv * ten3m / pc0
  phirf = rfl * twopi

  if (turn .lt. turn1)  then
     vrf = 0
  else if (turn .ge. turn1 .and. turn .lt. turn2) then
     vrf = (turn-turn1) * vrf / (turn2-turn1)
  else if (turn .ge. turn2 .and. turn .lt. turn3) then
     vrf = vrf
  else if (turn .ge. turn3 .and. turn .lt. turn4) then
     vrf = (turn4-turn) * vrf / (turn4-turn3)
  else
     vrf = 0
  endif
  !      if (turn .le. 10 .or. turn .gt. 1990) then
  !       print*," turn: ",turn, " vrf: ", vrf
  !      endif
  do itrack = 1, ktrack
     px  = track(2,itrack)                                           &
          + vrf * sin(phirf + omega * turn)

     !---- track(2,jtrk) = track(2,jtrk)
     !        pt = track(6,itrack)
     !        pt = pt + vrf * sin(phirf - omega * track(5,itrack))

     track(2,itrack) = px
  enddo
end subroutine tthacdip
subroutine ttvacdip(track,ktrack,turn)

  implicit none

  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Track a set of trajectories through a thin v ac dipole (zero l.)   *
  ! Input/output:                                                        *
  !   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
  !   KTRACK    (integer) number of surviving tracks.                    *
  !----------------------------------------------------------------------*
  ! Added by Yipeng SUN on 11 Nov 2009                                   *
  !----------------------------------------------------------------------*
  integer itrack,ktrack,turn,turn1,turn2,turn3,turn4
  double precision bi2gi2,dl,el,omega,dtbyds,phirf,pt,rff,betas,    &
       gammas,rfl,rfv,track(6,*),clight,twopi,vrf,deltap,deltas,pc,pc0,  &
       get_variable,node_value,get_value,one,two,half,ten3m,ten6p,py
  !      double precision px,py,ttt,beti,el1
  parameter(one=1d0,two=2d0,half=5d-1,ten3m=1d-3,ten6p=1d6)

  !---- Initialize
  clight=get_variable('clight ')
  twopi=get_variable('twopi ')

  !---- Fetch data.
  rfv = node_value('volt ')
  rff = node_value('freq ')
  rfl = node_value('lag ')
  pc0 = get_value('beam ','pc ')

  turn1 = node_value('ramp1 ')
  turn2 = node_value('ramp2 ')
  turn3 = node_value('ramp3 ')
  turn4 = node_value('ramp4 ')
  !---- Set up.
  omega = rff * twopi
  vrf   = 300 * rfv * ten3m / pc0
  phirf = rfl * twopi


  if (turn .lt. turn1)  then
     vrf = 0
  else if (turn .ge. turn1 .and. turn .lt. turn2) then
     vrf = (turn-turn1) * vrf / (turn2-turn1)
  else if (turn .ge. turn2 .and. turn .lt. turn3) then
     vrf = vrf
  else if (turn .ge. turn3 .and. turn .lt. turn4) then
     vrf = (turn4-turn) * vrf / (turn4-turn3)
  else
     vrf = 0
  endif
  !      if (turn .le. 10 .or. turn .gt. 1990) then
  !       print*," turn: ",turn, " vrf: ", vrf
  !      endif
  do itrack = 1, ktrack
     py  = track(4,itrack)                                           &
          + vrf * sin(phirf + omega * turn)

     !---- track(2,jtrk) = track(2,jtrk)
     !        pt = track(6,itrack)
     !        pt = pt + vrf * sin(phirf - omega * track(5,itrack))

     track(4,itrack) = py
  enddo
end subroutine ttvacdip
