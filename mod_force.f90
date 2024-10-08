!! Module to define tha parameters
!! ---------------------------------------------------------
!!
module mod_force
  use mod_param
  use mtmod
  implicit none
  save
  !
contains
!! ---------------SUBROUINES ---------------
!!
!! Intialize Cells
!- For Matrix Case
subroutine init_matrix
  implicit none
  integer::ip,ii,jj,nrs0
  real*8::rrnd,trnd
! Initialize the position of the particle in the chain (boundary of the tissue)
  dth = (2.0*pi)/dble(Nrs)
  is = 0
  do ii = 1,nint(R1-R0);  ! Matrix 
    nrs0 = int(2.0*pi*(R0+ii-1)*rhos)
    dth = (2.0*pi)/dble(nrs0)
    do jj = 1,nrs0;  ! Matrix 
      angle = dble(jj)*dth+(ii-1)*dth*0.5
      is = is + 1
      xrs(is) = dble(R0+ii-1)*dcos(angle)
      yrs(is) = dble(R0+ii-1)*dsin(angle)
    enddo
  enddo
  Nrs = is;

!!- Random packing
!  do is = 1,Nrs;  ! Matrix 
!62  rrnd = (R1-R0)*grnd()+R0;  trnd = (2.0*pi)*grnd()     ! For compact packing
!    xrs(is) = rrnd*cos(trnd);
!    yrs(is) = rrnd*sin(trnd);
!    dr0 = 2.0*a;
!    do ii = 1,is-1;
!      yij = yrs(ii)-yrs(is); 
!      xij = xrs(ii)-xrs(is);
!      rij = dsqrt(xij**2+yij**2);
!      if(rij<dr0) dr0 = rij
!    enddo
!    if(dr0<a) goto 62
!  enddo; 
  xrs0 = xrs;
  yrs0 = yrs;
  call write_track
! Populate the rest of PSM 
  ! When whole spheroid is populated
   Npsm = int(pi*(R0-a)**2/a*rho0) 
  ! When only a bilayer shell/ring is populated
  ! Npsm = ceiling(2.0*pi*(R0-3.0*a)*2.0/a*rho0)
  do np = 1,Npsm;
63  rrnd = (R0-a)*grnd();  trnd = (2.0*pi)*grnd()     ! For compact packing
!63  rrnd = (R0-4.0*a)+2.0d0*grnd();  trnd = (2.0*pi)*grnd()
    xx(np) = rrnd*cos(trnd);
    yy(np) = rrnd*sin(trnd);
    dr0 = 2.0*a;
    do ii = 1,np-1;
      yij = yy(np)-yy(ii); 
      xij = xx(np)-xx(ii);
      rij = dsqrt(xij**2+yij**2);
      if(rij<dr0) dr0 = rij
    enddo
    if(dr0<a) goto 63
  enddo; 
  np = Npsm;
  do ip = 1,np; 
  dtime(ip) = dble(ip-np)/dble(np)*tauG ;
  enddo
  call write_track
! The motility of cells decreases as the cells move away from tailbud
   mot = mot0; ttp = 0.0d0;
   if (FGF_grad) then
   do ip = 1,np;
     ttp(ip) = (R0-xx(ip))/Df
     mot(ip) = mot0*dexp(-ttp(ip)/tau)
   enddo
   endif

endsubroutine init_matrix

!!
!- For TB Case
subroutine init_TB
  implicit none
  integer::ip,ii
  real*8::rnumx,rnumy,Ly0

! Initialize the position of the first Npy particle in the chain (Tail-bud)
    Ny0 = Npy*Aratio/(Aratio+1)/2;
    Nx0 = Npy/(Aratio+1)/2;
    Lx0 = dble(Nx0/2);
    Ly0 = dble(Ny0/2)
  do ip = 1,Ny0;    
!!
    xrs(ip) = dble(Nx0/2);
    yrs(ip) = dble(ip-Ny0/2);

    xrs(Nx0+Ny0+ip) = -dble(Nx0/2);
    yrs(Nx0+Ny0+ip) = -dble(ip-Ny0/2);

!!
  enddo
    do ip = 1,Nx0;    
!!
    xrs(Ny0+ip) = -dble(ip-Nx0/2);
    yrs(Ny0+ip) = dble(Ny0/2);

    xrs(2*Ny0+Nx0+ip) = dble(ip-Nx0/2);
    yrs(2*Ny0+Nx0+ip) = -dble(Ny0/2);
!!
  enddo

  do ip = 1,Npy;
    write(17,*)xrs(ip),yrs(ip)
  enddo
  Nrs = Npy;
  call write_track
! Populate the rest of PSM 
  Npsm = int((Nx0-1)*(Ny0-1)*rho0)
  print*,Npsm
  do np = 1,Npsm;
63  rnumx = grnd();  rnumy = grnd()
    xx(np) = dble(Nx0-1)*(rnumx-0.5d0);
    yy(np) = dble(Ny0-1)*(rnumy-0.5d0);
    dr0 = 2.0*a;
    !print*,yy(np),xx(np),np,dr0 
    do ii = 1,np-1;
      yij = yy(np)-yy(ii); yij=yij-dnint(yij/Ly)*Ly
      xij = xx(np)-xx(ii);
      rij = dsqrt(xij**2+yij**2);
      if(rij<dr0) dr0 = rij
    enddo
    !print*,yy(np),xx(np),np,dr0 
    if(dr0<a) goto 63
  enddo;
  np = Npsm;
  call write_track
  do ip = 1,Npsm;
    write(18,*)xx(ip),yy(ip)
  enddo
! The motility of cells decreases as the cells move away from tailbud
   mot = mot0; ttp = 0.0d0;
   if (FGF_grad) then
   do ip = 1,np;
     ttp(ip) = (Ly0-yy(ip))/Df
     mot(ip) = mot0*dexp(-ttp(ip)/tau)
   enddo
   endif

endsubroutine init_TB

!!
! Force acting on the cell aggregate 
subroutine cal_force_matrix
  implicit none
  integer::ip,jp
    do is = 1,Nrs;
! Soft repulsion
      do jp = 1,Nrs;
  ! To avoid self counting
        if(jp==is) then
        else
          yij = yrs(jp)-yrs(is); 
          xij = xrs(jp)-xrs(is);
          rij = dsqrt(xij**2+yij**2);
          rij = max(rij,rmin)
       !   if(rij<a) then
       !     fxrs(is) = fxrs(is) + kss*(rij-a)*xij/rij;
       !     fyrs(is) = fyrs(is) + kss*(rij-a)*yij/rij;
       !   endif
          !fLJ = -4.d0*epsm*((rn*Ar/(rij)**(rn+2)-an*Aa0/(rij)**(an+2))) 
! Matrix-matrix LJ interaction
          fLJ = 0.0d0
          if(rij<rlj) then
          fLJ = -(2.0d0*alphar*epsm/rij**4)*((rlj/rij)**2-1.0d0)*(3.0d0*(a*rlj/rij)**2-2.0d0*rlj*rlj-a*a)
          endif
! Matrix-matrix LJ interaction
          fxrs(is) = fxrs(is) + fLJ*xij;
          fyrs(is) = fyrs(is) + fLJ*yij;
        endif  
      enddo
! Link rearrangement
       yij = yrs0(is)-yrs(is); 
       xij = xrs0(is)-xrs(is);
       rij = dsqrt(xij**2+yij**2);
       if(rij>ls_cut) then
         xrs0(is) = xrs(is);
         yrs0(is) = yrs(is);
       endif
! Elastic Contribution
       fxrs(is) = fxrs(is) + Es*xij;
       fyrs(is) = fyrs(is) + Es*yij;
! Inter species repulsion
      do jp = 1,np;
        yij = yy(jp)-yrs(is); 
        xij = xx(jp)-xrs(is);
        rij = dsqrt(xij**2+yij**2);
        if(rij<a) then
          fxrs(is) = fxrs(is) + kas*(rij-a)*xij/rij;
          fyrs(is) = fyrs(is) + kas*(rij-a)*yij/rij;
        endif
      enddo
    enddo
end subroutine cal_force_matrix
!!

! Force acting on the TB 
subroutine cal_force_TB
  implicit none
  integer::ip,jp
  real*8::rin,rip,xin,xip,yin,yip,aij
    do is = 1,Nrs;
! Spring force calculation    
      if(is==1) then
        yip = yrs(is)-yrs(is+1);       yin = yrs(is)-yrs(Npy);
        xip = xrs(is)-xrs(is+1);       xin = xrs(is)-xrs(Npy);
      else if(is==Npy) then
        yip = yrs(is)-yrs(1);       yin = yrs(is)-yrs(is-1);
        xip = xrs(is)-xrs(1);          xin = xrs(is)-xrs(is-1);
      else
        yip = yrs(is)-yrs(is+1);       yin = yrs(is)-yrs(is-1);
        xip = xrs(is)-xrs(is+1);       xin = xrs(is)-xrs(is-1);
      end if
    if(TB_Wall) then 
      if(is==1) then
        yin = a/dsqrt(2.0d0);  xin = a/dsqrt(2.0d0);
      else if(is==Npy) then
        yip = -a/dsqrt(2.0d0);  xip = -a/dsqrt(2.0d0);
      endif
    endif
      !if(yin.lt.0.0d0) yin = yin+Ly
      !if(yip.gt.0.0d0) yip = yip-Ly
      !if(yin.gt.5.0d0) yin = 5.d0*a;
      !if(yip.lt.5.0d0) yip = -5.d0*a;
      rip = dsqrt(xip**2+yip**2);  rin = dsqrt(xin**2+yin**2);
      fxrs(is) = fxrs(is) - kc*((rip-a)*xip/rip + (rin-a)*xin/rin);
      fyrs(is) = fyrs(is) - kc*((rip-a)*yip/rip + (rin-a)*yin/rin);
! Inter species repulsion
      do jp = 1,np;
        yij = yy(jp)-yrs(is); 
        if(TB_Wall) then;else;    yij=yij-dnint(yij/Ly)*Ly;    endif
        xij = xx(jp)-xrs(is);
        rij = dsqrt(xij**2+yij**2);
! Tailbud to rest; force calculation
        !fLJ = -4.d0*epst*((r12*Ar/(rij)**(r12+2)-a6*Aamin/(rij)**(a6+2)))
        fLJ = 0.0d0
        if(rij<sigma_min) then
          aij = sigma_min;
          fLJ = -(2.0d0*alphar*epst*aij*aij/rij**4)*((aij/rij)**2-1.0d0)**2 !Only repulsion 
        endif
        fxrs(is) = fxrs(is) + fLJ*xij;
        fyrs(is) = fyrs(is) + fLJ*yij;
        !if(rij<a) then
        !  fxrs(is) = fxrs(is) + kas*(rij-a)*xij/rij;
        !  fyrs(is) = fyrs(is) + kas*(rij-a)*yij/rij;
        !endif
      enddo
    enddo
end subroutine cal_force_TB
!!

! Force acting on the cell aggregate 
subroutine cal_force_cell
  implicit none
  integer::ip,jp
  real*8::aij
      rho = 0.0d0
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
          if(rij<sigma_min) then
            aij  = sigma_min;
            fLJ = -(2.0d0*alphar*epst*aij*aij/rij**4)*((aij/rij)**2-1.0d0)**2 !Only repulsion 
          endif
          fx(ip) = fx(ip) + fLJ*xij;
          fy(ip) = fy(ip) + fLJ*yij;
        endif
      enddo
!!  
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
          fLJ = 0.0d0
          if(rij<rlj) then
            !! Modified by Anupam to bring in repulsion
            aij = a;
            !fLJ = -(2.0d0*alphar*eps0/rij**4)*((rlj/rij)**2-1.0d0)*(3.0d0*(aij*rlj/rij)**2-2.0d0*rlj*rlj-aij**2)
            aij  = sigma_min;
            fLJ = -(2.0d0*alphar*eps0*aij*aij/rij**4)*((aij/rij)**2-1.0d0)**2 !Only repulsion 
          endif
          fx(ip) = fx(ip) + fLJ*xij;
          fy(ip) = fy(ip) + fLJ*yij;
        endif  
      enddo
    enddo
end subroutine cal_force_cell
!!

!! Evolve the matrix 
subroutine evolve_matrix
  implicit none
  integer::ip
    do ip = 1,Nrs;
      xrs(ip) = xrs(ip) + delta*(fxrs(ip))/mu ! This term can be evolved using second or third order schemes as well
      yrs(ip) = yrs(ip) + delta*(fyrs(ip))/mu ! This term can be evolved using second or third order schemes as well
    enddo
end subroutine evolve_matrix
!!

!! Evolve the TB 
subroutine evolve_TB
  implicit none
  integer::ip
  if(TB_Wall) then;
    do ip = 1,Nrs;
      xrs(ip) = xrs(ip) + delta*(fxrs(ip))/mu ! This term can be evolved using second or third order schemes as well
      yrs(ip) = yrs(ip) + delta*(fyrs(ip))/mu ! This term can be evolved using second or third order schemes as well
    enddo
    yrs(1) = 0.5d0*a;    yrs(Nrs) = (dble(Npy)-0.5d0)*a;
  else;
    do ip = 1,Nrs;
      xrs(ip) = xrs(ip) + delta*(fxrs(ip))/mu ! This term can be evolved using second or third order schemes as well
      yrs(ip) = yrs(ip) + delta*(fyrs(ip))/mu ! This term can be evolved using second or third order schemes as well
      !yrs(ip) = yrs(ip) - floor(yrs(ip)/Ly)*Ly ! Periodic Boundary condition
    enddo
  endif
end subroutine evolve_TB
!!

!! Evolve the cell: Deterministic process
subroutine evolve_cell
  implicit none
  integer::ip
   ! Deterministic process
    do ip = 1,np;
      xx(ip) = xx(ip) + delta*(fx(ip))/gm ! This term can be evolved using second or third order schemes as well
      yy(ip) = yy(ip) + delta*(fy(ip))/gm ! This term can be evolved using second or third order schemes as well
    enddo
!!
end subroutine evolve_cell
!!

!! Evolve the cell: Wiener process
subroutine evolve_cell_Weiner
  implicit none
  integer::ip
  real*8::rnumx,rnumy
   ! Wiener process
    do ip = 1,np;
    !if (FGF_grad) mot(ip) = mot0*dexp(-(ttp(ip))/(tau_current))
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
!    if(TB_1D) then
!      do ip = 1,np;
!        if(TB_Wall) then;
!          if(yy(ip).lt.0.5d0*a) yy(ip) = a-yy(ip)  !Reflecting wall condition at y = 0
!          if(yy(ip).gt.(Npy-0.5d0)*a) yy(ip) = 2.0d0*Npy-a-yy(ip)  !Reflecting wall condition at y = Npy*a
!        else
!          yy(ip) = yy(ip) - floor(yy(ip)/Ly)*Ly ! Periodic Boundary condition
!        endif
!        if(xx(ip).lt.0.5d0*a) xx(ip) = a-xx(ip)  !Reflecting wall condition at x = 0
!      enddo
!    endif



end subroutine evolve_cell_Weiner
!!

!!-- Cell proliferation
subroutine cell_proliferation
  implicit none
  integer::ip,ii,iir,np0
  real*8::Unp,dPr
  real*8::dstep
    !rnum = grnd() !White noise for position of particle in y direction
! Check the density if cell generation is needed.
! Cell generation; add one cell, if there is enough space. 
    do ip = 1,np;
        np0 = np
        iir = floor((np0-1)*grnd())+1 !White noise for position of particle in y direction
      if((tt-dtime(iir))>=tauG) then
        angle = (2.0*pi)*grnd()
        xx(np+1) = (xx(iir))+a*cos(angle) 
        yy(np+1) = (yy(iir))+a*sin(angle) 
        np = np+1;
        ttp(np) = 0.d0;
        dr0 = 2.0*a;
        Unp = 0.0d0;
        do ii = 1,np-1;
          yij = yy(np)-yy(ii); 
          xij = xx(np)-xx(ii);
          rij = dsqrt(xij**2+yij**2);
          if(rij<dr0) dr0 = rij
          if(rij<rc) then
            Unp = Unp + 4.d0*eps0*((Ar/(rij)**(rn)-0.5d0*(Aa(ip)+Aa(ii))/(rij)**(an)))
          endif
        enddo
!  Loop for the wall
        do ii = 1,nrs;
          yij = yy(np)-yrs(ii); 
          xij = xx(np)-xrs(ii);
          rij = dsqrt(xij**2+yij**2);
          if(rij<dr0) dr0 = rij
          if(rij<a) then
            Unp = Unp + 0.5d0*kas*(rij-a)**2;
          endif
        enddo
        dPr = dexp(-Unp/mot0)
        dstep = grnd()
        
        !if(dr0<0.9*a) then
        if(dstep>dPr) then
          np = np-1;  
        else
          print*,dstep,dPr,dr0,np,it
          dtime(iir) = tt
          dtime(np) = tt
        endif
        rhotemp = dble(np)/(pi*R0**2)
      endif
    enddo !Cell proliferation
end subroutine cell_proliferation


!!-- Cell proliferation
subroutine cell_proliferation_uniform
  implicit none
  integer::ip,ii,iir,np0
  real*8::Unp,dPr
  real*8::dstep
    !rnum = grnd() !White noise for position of particle in y direction
! Check the density if cell generation is needed.
! Cell generation; add one cell, if there is enough space. 
        np0 = np
    do ip = 1,np;
        ii = floor((np0-1)*grnd())+1 !White noise for picking the random cell
 print*,ii,dtime(ii)
      if((tt-dtime(ii))>=tauG) then
        angle = (2.0*pi)*grnd()
        xx(np+1) = (xx(ii))+a*cos(angle)
        yy(np+1) = (yy(ii))+a*sin(angle)
        np = np+1;
        ttp(np) = 0.d0;
        print*,dstep,dPr,dr0,np,it
        dtime(iir) = tt
        dtime(np) = tt
      endif
    if(np.gt.np0) goto34
    enddo !Cell proliferation fix
34  print*,np0,np
end subroutine cell_proliferation_uniform


!!-- Cell proliferation for TB
subroutine cell_proliferation_TB
  implicit none
  integer::ip,ii,inp,np0
  real*8::rnum
    !rnum = grnd() !White noise for position of particle in y direction
! Check the density if cell generation is needed.
! Cell generation; add one cell, if there is enough space. 
      rhotemp = dble(np)/(Lx0*dble(Npy))
      np0 = np;
    do ip = 1,Nrs/3 ;
    if(rhotemp>rhom) then
      goto 65
    else
      rnum = grnd() !White noise for position of particle in y direction
      np = np+1;
      yy(np) = dble(Npy)*a*rnum;
      inp = ceiling(yy(np))
      xx(np) = xrs(inp)-a;
      ttp(np) = 0.d0;
      dr0 = 2.0*a;
      do ii = 1,np-1;
        yij = yy(np)-yy(ii); yij=yij-dnint(yij/Ly)*Ly
        xij = xx(np)-xx(ii);
        rij = dsqrt(xij**2+yij**2);
        if(rij<dr0) dr0 = rij
      enddo
      if(dr0<a) then
         np = np-1;
      else
         if(Clock) theta(np) = theta(np0-1)+0.1*(grnd()-0.5d0);
      endif
      rhotemp = dble(np)/(Lx0*dble(Npy))
      ! print*,'rhotemp',rhotemp
    endif
    enddo
65    print*,tt,np,it,dr0,rhotemp
end subroutine cell_proliferation_TB

!!
end module mod_force
