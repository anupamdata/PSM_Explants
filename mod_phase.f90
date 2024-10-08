!! Module to define tha parameters
!! ---------------------------------------------------------
!!
module mod_phase
  use mod_param
  use mtmod
  implicit none
  real*8::thLJ,thij,sigma_theta
  real*8::omginf

  save
  !
contains
!! ---------------SUBROUINES ---------------
!!
!! Intialize Cell phase
!- For TB Case
subroutine init_phase_TB
!!
  implicit none
  integer::ip
!!
  sigma_theta = Df*rhom*tau/a; ! Length scale of PSM
  omginf = omgL/(1-exp(-Df*rhom*tau/(a*sigma_theta)));
  eps_theta = eps_theta0;

! Assign the phase and intrinsic frequency to the  PSM 
  do ip = 1,Npsm;
    omg(ip) = omginf*(1.0d0-exp(-(xx(ip)-(Lx0-Df*rhom*tau/a))/sigma_theta))
    theta(ip) = (grnd()-0.5) + omgL*(a*xx(ip)/Df*rhom);          !omgL x/v_TB
  enddo;

endsubroutine init_phase_TB

!!
! Cell phase interaction
!!
subroutine cal_force_phase_cell
  implicit none
  integer::ip,jp
  real*8::aij,aij1,aij2,tfact,avg_thij
!!  
    do ip = 1,np;
      ! Only repulsion from matrix
      do jp = 1,Nrs;
          xij = xrs(jp) - xx(ip);
        if(two_phase) then
          yij = yrs(jp) - yy(ip);
          rij = dsqrt(xij**2+yij**2);
        if(rij<a) then
          fx(ip) = fx(ip) + kas*(rij-a)*xij/rij;
          fy(ip) = fy(ip) + kas*(rij-a)*yij/rij;
        endif
        elseif(TB_1D) then
          yij = yrs(jp) - yy(ip);
          if(TB_Wall) then;else;    yij=yij-dnint(yij/Ly)*Ly;    endif
          rij = dsqrt(xij**2+yij**2);
          !fLJ = -4.d0*epst*((r12*Ar/(rij)**(r12+2)-a6*Aamin/(rij)**(a6+2))) 
          fLJ = 0.0d0
          !if(rij<rlj) then
          if(rij<sigma_min) then
            aij  = sigma_min;
            fLJ = -(2.0d0*alphar*epst*aij*aij/rij**4)*((aij/rij)**2-1.0d0)**2 !Only repulsion 
          endif
          fx(ip) = fx(ip) + fLJ*xij;
          fy(ip) = fy(ip) + fLJ*yij;
        endif
      enddo
!
      do jp = 1,np;
! To avoid self counting
        if(jp==ip) then
        else
          yij = yy(jp)-yy(ip);
          if(TB_1D)  then  
          if(TB_Wall) then;else;    yij=yij-dnint(yij/Ly)*Ly;    endif
          endif          
          xij = xx(jp)-xx(ip);
          rij = dsqrt(xij**2+yij**2);
          rij = max(rij,rmin)
! Cell-cell interaction
          thLJ = 0.0d0
          fLJ = 0.0d0
          if(rij<rlj) then
            thij = theta(jp)-theta(ip);
            avg_thij = (theta(jp)+theta(ip))*0.5d0;
            !aij1 = a+(1.0d0-dcos(thij))*0.30d0*(rlj-a)+(1.0d0-dcos(avg_thij))*0.70d0*(rlj-a) 
            aij1 = a+((1.0d0-dcos(thij))*0.30d0+(1.0d0-dcos(avg_thij))*0.70d0)*(rlj-a) 
            aij2 = rlj
            tfact = 0.0d0;!dexp(-(ttp(ip))/tau_current)
            aij = aij2*tfact+aij1*(1-tfact)
            fLJ = -(2.0d0*alphar*eps0/rij**4)*((rlj/rij)**2-1.0d0)*(3.0d0*(aij*rlj/rij)**2-2.0d0*rlj*rlj-aij**2)
            thLJ = eps_theta(ip)*(1.0d0/rij**2-1.0d0/rlj**2)*dsin(thij);
          endif
          fx(ip) = fx(ip) + fLJ*xij;
          fy(ip) = fy(ip) + fLJ*yij;
          fth(ip) = fth(ip) + thLJ;
        endif  
      enddo
    enddo
end subroutine cal_force_phase_cell
!!


!! Evolve the cell: Deterministic process
subroutine evolve_phase_cell
  implicit none
  integer::ip
   ! Deterministic process
    do ip = 1,np;
      omg(ip) = omginf*(1.0d0-exp(-(xx(ip)-(Lx0-Df*rhom*tau_currentw/a))/sigma_theta))
      !omg(ip) = omginf*(1.0d0-exp(-(xx(ip)-(Lx0-Df*rhom*tau/a))/sigma_theta))
      if((Lx0-xx(ip)).gt.Df*rhom*tau/a) omg(ip) = 0.0d0;
      xx(ip) = xx(ip) + delta*(fx(ip))/gm ! This term can be evolved using second or third order schemes as well
      yy(ip) = yy(ip) + delta*(fy(ip))/gm ! This term can be evolved using second or third order schemes as well
      theta(ip) = theta(ip) + delta*(omg(ip)+fth(ip)) ! This term can be evolved using second or third order schemes as well
    enddo
!!
end subroutine evolve_phase_cell
!!

!! Evolve the cell: Wiener process
subroutine evolve_phase_Weiner
  implicit none
  integer::ip
  real*8::rnumx,rnumy
   ! Wiener process
    do ip = 1,np;
    if (FGF_grad) mot(ip) = mot0*dexp(-(ttp(ip))/tau_current)
      ttp(ip) = ttp(ip) + delta
!64    rnumx = gasdev(); rnumy = gasdev();
!      !if(max(abs(rnumx),abs(rnumy))>4.0d0) 
!      if(max(abs(rnumx),abs(rnumy))>8.0d0) goto 64
      rnumx = (grnd()-0.5d0);   rnumy = (grnd()-0.5d0) !White noise 
      xix(ip) = dsqrt(2.d0*mot(ip)*gm)*rnumx
      xiy(ip) = dsqrt(2.d0*mot(ip)*gm)*rnumy
      xx(ip) = xx(ip) + dsqrt(delta)*(xix(ip))/gm ! Noise term; evolved at the end with a sqrt(delta) term
      yy(ip) = yy(ip) + dsqrt(delta)*(xiy(ip))/gm ! Noise term; evolved at the end with a sqrt(delta) term
    enddo
    if(TB_1D) then
      do ip = 1,np;
        if(TB_Wall) then;
          if(yy(ip).lt.0.5d0*a) yy(ip) = a-yy(ip)  !Reflecting wall condition at y = 0
          if(yy(ip).gt.(Npy-0.5d0)*a) yy(ip) = 2.0d0*Npy-a-yy(ip)  !Reflecting wall condition at y = Npy*a
        else
          yy(ip) = yy(ip) - floor(yy(ip)/Ly)*Ly ! Periodic Boundary condition
        endif
        if(xx(ip).lt.0.5d0*a) xx(ip) = a-xx(ip)  !Reflecting wall condition at x = 0
      enddo
    endif



end subroutine evolve_phase_Weiner

!!
end module mod_phase
