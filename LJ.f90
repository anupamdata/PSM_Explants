program repulsion
  use mod_param
  use mod_force 
  use mod_phase 
  use mtmod

  implicit none
  integer*4::ip,jp,inp
  integer*4::iloop,nloop
  real*8::rip,xip,yip,rin,xin,yin
  real*8::rnum,delta0,maxf,tot_t
!! 
  
  call read_input
  mot0 = Df*gm;  ! Single cell motility
  R0 = dble(Npr0)*a;  ! Radius of sphere.
  R1 = dble(Npr1)*a;  ! Radius of sphere (Cell+Matrix).
  Ly = dble(Npy)*a; !PBC in Y, with system length Ly
  Lx = dble(Npx)*a;
  Lx0 = 20.0*a;
  rc = 3.0d0*a;   ! Cut-off for LJ, here it is being used to determine the density in the neighbourhood of cells.
  area_c = pi*rc*rc
  if(two_phase) Nrs = int(pi*(R1**2-R0**2)*rhos)
  if(TB_1D) Nrs = Npy

!!
  write(cnpar,'(g8.0)') 3+Npy*ndata
  formp = '('//trim(adjustl(cnpar))//'es16.3E2)'
  call open_track
!!

! Allocate the arrays
  allocate(xx(Npt),yy(Npt),vv(Npt),Pr(Npt),xix(Npt),xiy(Npt))
  allocate(dtime(Npt)) !Time of last division
  allocate(xx0(Npt),yy0(Npt),mot(Npt),ttp(Npt))
  allocate(xrs(Nrs),yrs(Nrs))
  allocate(xrs0(Nrs),yrs0(Nrs))
  allocate(fxrs(Nrs),fyrs(Nrs))
  allocate(fx(Npt),fy(Npt))
  allocate(Aa(Npt),rho(Npt))
  allocate(ran(3))
  if (Clock) then 
    allocate(theta(Npt),omg(Npt),eps_theta(Npt),fth(Npt))
  endif
  call sgrnd(iseed)
  xx = -2.0d0; yy = 0.0d0;
  dtime = 0.0d0;
  gmm = gm;
  ! y dependent friction factor 
  mot = 0;  ! Inititalize the motility as 0
  Aa = a; !Aamin;    ! Inititalize the attraction potential Aa0

  if(nrun==1) then
  if(two_phase) call init_matrix
  if(TB_1D) call init_TB    
  if(Clock) call init_phase_TB
  else
! Read from old trajectory file.
     call read_traj
  endif
!
!! Time marching 
  rhotemp = rho0;
  delta0 = delta;
  do it = 1,Nt;
    tt = dble(it)*delta;
    tau_current = tau*dexp(-tt/taud)
    tau_currentw = tau*dexp(-tt/taudw)
!
!! Calculate the force acting on the particle. f(rij) = - [\partial U/\partial r]
    fx = 0.0d0;  fy = 0.0d0;
    fxrs = 0.0d0;  fyrs = 0.0d0;
!!
! First for the substrate
!
  if(two_phase) call cal_force_matrix
  if(TB_1D) call cal_force_TB
!!
! Second for the cells
  if(Clock) then 
    fth = 0.0d0;
    call cal_force_phase_cell
  else
    call cal_force_cell
  endif

!! Adaptive time step
   delta = delta0
   maxf = (max(maxval(abs(fx)),maxval(abs(fy)))*(delta/gm))/a;
   if(maxf.gt.0.1d0) then
     nloop = (ceiling(maxf/0.1d0))
     delta = delta0/dble(nloop)
     do iloop = 1,nloop
    fx = 0.0d0;  fy = 0.0d0;
    fxrs = 0.0d0;  fyrs = 0.0d0;
       if(two_phase) call cal_force_matrix
       if(TB_1D) call cal_force_TB
       if(Clock) then 
         fth = 0.0d0;
         call cal_force_phase_cell
       else
         call cal_force_cell
       endif
       if(two_phase) call evolve_matrix
       if(TB_1D) call evolve_TB
       if(Clock) then 
         call evolve_phase_cell
       else
         call evolve_cell
       endif
     enddo
   write(99,*)delta,delta0
     goto 75
   endif

!! Evolving the particles 
!
!! Evolve the matrix
    if(two_phase)  call evolve_matrix
!! Evolve the TB
    if(TB_1D)  call evolve_TB
!
!! Evolve the rest of the cells 
   ! Deterministic process
    if(Clock) then 
      call evolve_phase_cell
    else
      call evolve_cell
    endif
   ! Wiener process
75    delta = delta0
    call evolve_cell_Weiner

!!-- Cell proliferation
    if(stress_dependent) call cell_proliferation
    if((uniform_prol).and.((mod(it,int(0.1d0*tauG/delta)).eq.0))) call cell_proliferation_uniform
    if(TB_1D)  then
!      call cell_proliferation_TB
     if(mod(it,250).eq.0) print*,'00',sum(xx(1:Npsm))/Npsm,sum(yy(1:Npsm))/Npsm
    endif
!!
!!
  if(mod(it,1000).eq.0)call write_track
  !call write_track
!!
  enddo

  call close_track

end program repulsion



