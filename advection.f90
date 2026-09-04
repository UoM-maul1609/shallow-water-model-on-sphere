	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>advection routines
    module advection
    use numerics_type
    private
    public :: lax_wendroff_sphere, lax_wendroff_ll, dissipation, smagorinsky, &
              smagorinsky_spherical_stress, spherical_stress_divergence
    contains
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>advects a scalar field on the sphere 
	!>@param[in] ip: number of east-west points
	!>@param[in] jp: ditto for north-south
	!>@param[in] o_halo: halos required for advection scheme
	!>@param[in] dt:  timestep
	!>@param[in] dx,dy:  dx,dy
	!>@param[in] g: gravity
	!>@param[inout] u,v,h: prognostic variables
	!>@param[in] hs: surface height
	!>@param[in] re: radius of planet
	!>@param[in] theta: latitude
	!>@param[in] f_cor: Coriolis parameter
	!>solves the 1-d advection equation:
	!>\f$ \frac{\partial \psi}{\partial t} + \frac{\partial u \psi}{\partial x} = 0 \f$
    subroutine lax_wendroff_sphere(ip,jp,o_halo,dt,dx,dy,g,u,v,h,hs,re,theta,f_cor)

		use numerics_type
		implicit none
		integer(i4b), intent(in) :: ip,jp,o_halo
		real(wp), intent(in) :: dt, g, re
		real(wp), intent(in), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: &
																		hs, f_cor, dx, dy
		real(wp), dimension(1-o_halo:jp+o_halo), intent(in) :: theta
		real(wp), intent(inout), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: &
																				h, u, v
																				
		! local variables:
		real(wp), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: & 
					dy1, v1, h1, vh, uh, vh1, Ux, Uy, Vx, Vy, Vy2
		real(wp), dimension(1:ip,1:jp) :: &
									    uh_new, vh_new, h_new
										
		real(wp), dimension(0:ip,1:jp) :: h_mid_xt, uh_mid_xt, Ux_mid_xt, Vx_mid_xt
		real(wp), dimension(0:ip,1:jp) :: vh_mid_xt
		real(wp), dimension(1:ip,0:jp) :: h_mid_yt, uh_mid_yt, vh_mid_yt, c_mid_yt, &
		 									Uy_mid_yt, Vy_mid_yt, Vy_mid_yt2
		integer(i4b) :: j, i
			
		
		do j=1-o_halo,jp+o_halo
			dy1(:,j)=dy(:,j)*cos(theta(j))
			v1(:,j)=v(:,j)*cos(theta(j))
			h1(:,j)=h(:,j)*cos(theta(j))
		enddo
		do j=0,jp
			c_mid_yt(:,j)=cos(0.5_wp*(theta(j+1)+theta(j)))
		enddo
		
		uh=u *h
		vh=v*h
		vh1=v1*h
		
		! continuity equation (calculate mid-point values at 0.5*dt):
		h_mid_xt = 0.5_wp*(h(1:ip+1,1:jp)+h(0:ip,1:jp)) &
		  -(0.5_wp*dt/(0.5_wp*(dx(1:ip+1,1:jp)+dx(0:ip,1:jp)))) &
		  *(uh(1:ip+1,1:jp)-uh(0:ip,1:jp))
		  
		h_mid_yt = 0.5_wp*(h(1:ip,1:jp+1)+h(1:ip,0:jp)) &
		  -(0.5_wp*dt/(0.5_wp*(dy1(1:ip,1:jp+1)+dy1(1:ip,0:jp)))) &
		  *(vh1(1:ip,1:jp+1)-vh1(1:ip,0:jp))

		! v-phi, or u momentum equation (calculate mid-point values at 0.5*dt):
		Ux = uh*u + g*h**2*(0.5_wp)
		Uy = uh*v1
		uh_mid_xt(0:ip,:) = 0.5_wp*(uh(1:ip+1,1:jp)+uh(0:ip,1:jp)) &
		  -(0.5_wp*dt/(0.5_wp*(dx(1:ip+1,1:jp)+dx(0:ip,1:jp))))* &
		  	(Ux(1:ip+1,1:jp)-Ux(0:ip,1:jp)) &
		  +0.125_wp*dt*(f_cor(1:ip+1,1:jp)+f_cor(0:ip,1:jp))* &
		  	(vh(1:ip+1,1:jp)+vh(0:ip,1:jp))
		
		uh_mid_yt = 0.5_wp*(uh(1:ip,1:jp+1)+uh(1:ip,0:jp)) &
		  -(0.5_wp*dt/(0.5_wp*(dy1(1:ip,1:jp+1)+dy1(1:ip,0:jp))))* &
		  	(Uy(1:ip,1:jp+1)-Uy(1:ip,0:jp)) &
		  +0.125_wp*dt*(f_cor(1:ip,1:jp+1)+f_cor(1:ip,0:jp))* &
		  	(vh(1:ip,1:jp+1)+vh(1:ip,0:jp))



		! v-theta, or v momentum equation (calculate mid-point values at 0.5*dt):
		Vx = uh*v
		Vy = vh1*v
		Vy2 = 0.5_wp*g*h**2
		vh_mid_xt(0:ip,1:jp) = 0.5_wp*(vh(1:ip+1,1:jp)+vh(0:ip,1:jp)) &
		  -(0.5_wp*dt/(0.5_wp*(dx(1:ip+1,1:jp)+dx(0:ip,1:jp))))*(Vx(1:ip+1,1:jp)-Vx(0:ip,1:jp)) &
		  -0.125_wp*dt*(f_cor(1:ip+1,1:jp)+f_cor(0:ip,1:jp))*(uh(1:ip+1,1:jp)+uh(0:ip,1:jp))

		vh_mid_yt(1:ip,0:jp) = 0.5_wp*(vh(1:ip,1:jp+1)+vh(1:ip,0:jp)) &
		  -(0.5_wp*dt/(0.5_wp*(dy1(1:ip,1:jp+1)+dy1(1:ip,0:jp))))*(Vy(1:ip,1:jp+1)-Vy(1:ip,0:jp)) &
		  -(0.5_wp*dt/(dy(1:ip,0:jp)))*(Vy2(1:ip,1:jp+1)-Vy2(1:ip,0:jp)) &
		  -0.125_wp*dt*(f_cor(1:ip,1:jp+1)+f_cor(1:ip,0:jp))*(uh(1:ip,1:jp+1)+uh(1:ip,0:jp))

! 		calculate mid-point value of cos (theta)
! 		c_mid_yt=cos(0.5.*(THETA(:,2:end)+THETA(:,1:end-1)));
! 
! 
		! Now use the mid-point values to predict the values at the next timestep
		! continuity:
		h_new = h(1:ip,1:jp) &
		  - (dt/(0.5_wp*(dx(1:ip,1:jp)+dx(0:ip-1,1:jp))))*(uh_mid_xt(1:ip,1:jp)-uh_mid_xt(0:ip-1,1:jp)) &
		  - (dt/(0.5_wp*(dy1(1:ip,1:jp)+dy1(1:ip,0:jp-1)))) * &
		  (vh_mid_yt(1:ip,1:jp)*c_mid_yt(1:ip,1:jp)-vh_mid_yt(1:ip,0:jp-1)*c_mid_yt(1:ip,0:jp-1))


		! u-momentum equation:
		Ux_mid_xt = uh_mid_xt*uh_mid_xt/h_mid_xt + 0.5_wp*g*h_mid_xt**2
		Uy_mid_yt = uh_mid_yt*vh_mid_yt/h_mid_yt*c_mid_yt
		uh_new = uh(1:ip,1:jp) &
		  - (dt/(0.5_wp*(dx(1:ip,1:jp)+dx(0:ip-1,1:jp))))*  (Ux_mid_xt(1:ip,1:jp)-Ux_mid_xt(0:ip-1,1:jp)) &
		  - (dt/(0.5_wp*(dy1(1:ip,1:jp)+dy1(1:ip,0:jp-1))))*(Uy_mid_yt(1:ip,1:jp)-Uy_mid_yt(1:ip,0:jp-1))


		! v-momentum equation:
		Vx_mid_xt = uh_mid_xt*vh_mid_xt/h_mid_xt
		Vy_mid_yt = vh_mid_yt*vh_mid_yt/h_mid_yt*c_mid_yt
		Vy_mid_yt2 = 0.5_wp*g*h_mid_yt**2
		vh_new = vh(1:ip,1:jp) &
		  - (dt/(0.5_wp*(dx(1:ip,1:jp)+dx(0:ip-1,1:jp))))*(Vx_mid_xt(1:ip,1:jp)-Vx_mid_xt(0:ip-1,1:jp)) &
		  - (dt/(0.5_wp*(dy1(1:ip,1:jp)+dy1(1:ip,0:jp-1))))*(Vy_mid_yt(1:ip,1:jp)-Vy_mid_yt(1:ip,0:jp-1)) &
		  - (dt/(dy(1:ip,0:jp-1) ))* &
		  (Vy_mid_yt2(1:ip,1:jp)-Vy_mid_yt2(1:ip,0:jp-1))


		! add on Coriolis and contribution of orography to pressure gradient:
		uh_new=uh_new  +dt*.5_wp*(f_cor(1:ip,1:jp)*v(1:ip,1:jp) - &
			g*(hs(2:ip+1,1:jp)-hs(0:ip-1,1:jp))/(2._wp*dx(1:ip,1:jp)))* &
			(h(1:ip,1:jp)+h_new)

		vh_new=vh_new  -dt*.5_wp*(f_cor(1:ip,1:jp)*u(1:ip,1:jp) + &
			g*(hs(1:ip,2:jp+1)-hs(1:ip,0:jp-1))/(2._wp*dy1(1:ip,1:jp)))* &
			(h(1:ip,1:jp)+h_new)

! 		uh_new=uh_new  +dt*.5_wp*(f_cor(1:ip,1:jp)*(vh_mid_xt(1:ip,1:jp)+vh_mid_xt(0:ip-1,1:jp)) )
! 
! 		vh_new=vh_new  -dt*.5_wp*(f_cor(1:ip,1:jp)*(uh_mid_xt(1:ip,1:jp)+uh_mid_xt(0:ip-1,1:jp)) )


		! re-calculate u and v.
		u(1:ip,1:jp) = uh_new/(h_new)
		v(1:ip,1:jp) = vh_new/h_new
		h(1:ip,1:jp) = h_new
																				
! 		do j=1,jp
! 			do i=1,ip
! 				
! 			
! 			enddo
! 		enddo																				


	end subroutine lax_wendroff_sphere
	
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>advects a scalar field on a ll grid 
	!>@param[in] ip: number of east-west points
	!>@param[in] jp: ditto for north-south
	!>@param[in] o_halo: halos required for advection scheme
	!>@param[in] dt:  timestep
	!>@param[in] g: gravity
	!>@param[inout] u,v,h: prognostic variables
	!>@param[in] hs: surface height
	!>@param[in] re: radius of planet
	!>@param[in] theta: latitude
	!>@param[in] thetan: latitude - staggered
	!>@param[in] dtheta: latitude step
	!>@param[in] dthetan: latitude step - staggered
	!>@param[in] phi: phi
	!>@param[in] phin: phin
	!>@param[in] dphi: dphi
	!>@param[in] dphin: dphin
	!>@param[in] f_cor: Coriolis parameter
	!>@param[in] recqdq - for efficiency
	!>@param[in] recqdp - for efficiency
	!>@param[in] recqdp_s - for efficiency
	!>@param[in] recqdq_s - for efficiency
	!>@param[in] redq_s - for efficiency
	!>@param[in] redq - for efficiency
	!>@param[in] cq - for efficiency
	!>@param[in] cq_s - for efficiency
	!>@param[in] momentum_metric_terms - 0 legacy/off; 1 spherical momentum curvature terms
	!>solves the 1-d advection equation:
	!>\f$ \frac{\partial \psi}{\partial t} + \frac{\partial u \psi}{\partial x} = 0 \f$
    subroutine lax_wendroff_ll(ip,jp,o_halo,dt,g,u,v,h,hs,re,&
    		theta,thetan,dtheta,dthetan, phi, phin, dphi, dphin, f_cor, &
    		recqdq, recqdp, recqdp_s, recqdq_s, redq_s, redq, cq, cq_s, coriolis_scheme, momentum_metric_terms)

		use numerics_type
		implicit none
		integer(i4b), intent(in) :: ip,jp,o_halo, coriolis_scheme, momentum_metric_terms
		real(wp), intent(in) :: dt, g, re
		real(wp), intent(in), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: &
																		hs, f_cor, &
    				recqdq, recqdp, recqdp_s, recqdq_s, redq_s, redq, &
    				cq, cq_s
		real(wp), dimension(1-o_halo:jp+o_halo), intent(in) :: &
											theta,thetan, dtheta, dthetan
		real(wp), dimension(1-o_halo:ip+o_halo), intent(in) :: &											
											phi, phin, dphi, dphin
		real(wp), intent(inout), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: &
																				h, u, v
																				
		! local variables:
		real(wp), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: & 
					dy1, v1, h1, vh, uh, vh1, Ux, Uy, Vx, Vy, Vy2, &
                    metric_u, metric_v
		real(wp), dimension(1:ip,1:jp) :: &
									    uh_new, vh_new, h_new, &
									    cor_alpha, cor_denom, cor_rhs_u, cor_rhs_v
										
		real(wp), dimension(0:ip,1:jp) :: h_mid_xt, uh_mid_xt, Ux_mid_xt, Vx_mid_xt
		real(wp), dimension(0:ip,1:jp) :: vh_mid_xt
		real(wp), dimension(1:ip,0:jp) :: h_mid_yt, uh_mid_yt, vh_mid_yt,  &
		 									Uy_mid_yt, Vy_mid_yt, Vy_mid_yt2
		integer(i4b) :: j, i
			
		
		v1=v*cq ! cq=cos(theta)
		h1=h*cq
		
		uh=u *h
		vh=v*h
		vh1=v1*h

        ! Spherical curvature/metric source terms for the conservative
        ! eastward/northward momenta. Keep them disabled by default for
        ! exact backwards compatibility with the original model.
        metric_u = 0._wp
        metric_v = 0._wp
        select case (momentum_metric_terms)
        case (0)
            ! Legacy behaviour: no momentum curvature terms.
        case (1)
            do j=1-o_halo,jp+o_halo
                metric_u(:,j) =  uh(:,j)*v(:,j)*tan(theta(j))/re
                metric_v(:,j) = -uh(:,j)*u(:,j)*tan(theta(j))/re
            enddo
        case default
            write(*,*) 'ERROR: unknown momentum_metric_terms = ', momentum_metric_terms
            write(*,*) '       valid values are 0 (legacy/off) and 1 (spherical/on)'
            stop 1
        end select
		
		! continuity equation (calculate mid-point values at 0.5*dt):
		h_mid_xt = 0.5_wp*(h(1:ip+1,1:jp)+h(0:ip,1:jp)) &
		  -(0.5_wp*dt/(recqdp_s(0:ip,1:jp))) &
		  *(uh(1:ip+1,1:jp)-uh(0:ip,1:jp))
		  
		h_mid_yt = 0.5_wp*(h(1:ip,1:jp+1)+h(1:ip,0:jp)) &
		  -(0.5_wp*dt/(recqdq_s(1:ip,0:jp))) &
		  *(vh1(1:ip,1:jp+1)-vh1(1:ip,0:jp))

		! v-phi, or u momentum equation (calculate mid-point values at 0.5*dt):
		Ux = uh*u + g*h**2*(0.5_wp)
		Uy = uh*v1
		uh_mid_xt(0:ip,:) = 0.5_wp*(uh(1:ip+1,1:jp)+uh(0:ip,1:jp)) &
		  -(0.5_wp*dt/(recqdp_s(0:ip,1:jp)))* &
		  	(Ux(1:ip+1,1:jp)-Ux(0:ip,1:jp)) &
		  +0.125_wp*dt*(f_cor(1:ip+1,1:jp)+f_cor(0:ip,1:jp))* &
		  	(vh(1:ip+1,1:jp)+vh(0:ip,1:jp))
		
		uh_mid_yt = 0.5_wp*(uh(1:ip,1:jp+1)+uh(1:ip,0:jp)) &
		  -(0.5_wp*dt/(recqdq_s(1:ip,0:jp)))* &
		  	(Uy(1:ip,1:jp+1)-Uy(1:ip,0:jp)) &
		  +0.125_wp*dt*(f_cor(1:ip,1:jp+1)+f_cor(1:ip,0:jp))* &
		  	(vh(1:ip,1:jp+1)+vh(1:ip,0:jp))



		! v-theta, or v momentum equation (calculate mid-point values at 0.5*dt):
		Vx = uh*v
		Vy = vh1*v
		Vy2 = 0.5_wp*g*h**2
		vh_mid_xt(0:ip,1:jp) = 0.5_wp*(vh(1:ip+1,1:jp)+vh(0:ip,1:jp)) &
		  -(0.5_wp*dt/(recqdp_s(0:ip,1:jp)))*(Vx(1:ip+1,1:jp)-Vx(0:ip,1:jp)) &
		  -0.125_wp*dt*(f_cor(1:ip+1,1:jp)+f_cor(0:ip,1:jp))*(uh(1:ip+1,1:jp)+uh(0:ip,1:jp))

		vh_mid_yt(1:ip,0:jp) = 0.5_wp*(vh(1:ip,1:jp+1)+vh(1:ip,0:jp)) &
		  -(0.5_wp*dt/(recqdq_s(1:ip,0:jp)))*(Vy(1:ip,1:jp+1)-Vy(1:ip,0:jp)) &
		  -(0.5_wp*dt/(redq_s(1:ip,0:jp)))*(Vy2(1:ip,1:jp+1)-Vy2(1:ip,0:jp)) &
		  -0.125_wp*dt*(f_cor(1:ip,1:jp+1)+f_cor(1:ip,0:jp))*(uh(1:ip,1:jp+1)+uh(1:ip,0:jp))

        ! Add half-step spherical momentum metric sources at the faces.
        ! The arithmetic mean of the adjacent cell-centred source is used,
        ! multiplied by dt/2 for the Lax-Wendroff predictor.
        if (momentum_metric_terms == 1) then
            uh_mid_xt(0:ip,1:jp) = uh_mid_xt(0:ip,1:jp) + 0.25_wp*dt* &
                (metric_u(1:ip+1,1:jp)+metric_u(0:ip,1:jp))
            uh_mid_yt(1:ip,0:jp) = uh_mid_yt(1:ip,0:jp) + 0.25_wp*dt* &
                (metric_u(1:ip,1:jp+1)+metric_u(1:ip,0:jp))
            vh_mid_xt(0:ip,1:jp) = vh_mid_xt(0:ip,1:jp) + 0.25_wp*dt* &
                (metric_v(1:ip+1,1:jp)+metric_v(0:ip,1:jp))
            vh_mid_yt(1:ip,0:jp) = vh_mid_yt(1:ip,0:jp) + 0.25_wp*dt* &
                (metric_v(1:ip,1:jp+1)+metric_v(1:ip,0:jp))
        endif

! 		calculate mid-point value of cos (theta)
! 		c_mid_yt=cos(0.5.*(THETA(:,2:end)+THETA(:,1:end-1)));
! 
! 
		! Now use the mid-point values to predict the values at the next timestep
		! continuity:
		h_new = h(1:ip,1:jp) &
		  - (dt/(recqdp(0:ip-1,1:jp)))*(uh_mid_xt(1:ip,1:jp)-uh_mid_xt(0:ip-1,1:jp)) &
		  - (dt/(recqdq(1:ip,0:jp-1))) * &
		  (vh_mid_yt(1:ip,1:jp)*cq_s(1:ip,1:jp)-vh_mid_yt(1:ip,0:jp-1)*cq_s(1:ip,0:jp-1))


		! u-momentum equation:
		Ux_mid_xt = uh_mid_xt*uh_mid_xt/h_mid_xt + 0.5_wp*g*h_mid_xt**2
		Uy_mid_yt = uh_mid_yt*vh_mid_yt/h_mid_yt*cq_s(1:ip,0:jp)
		uh_new = uh(1:ip,1:jp) &
		  - (dt/(recqdp(0:ip-1,1:jp)))*  (Ux_mid_xt(1:ip,1:jp)-Ux_mid_xt(0:ip-1,1:jp)) &
		  - (dt/(recqdq(1:ip,0:jp-1)))*(Uy_mid_yt(1:ip,1:jp)-Uy_mid_yt(1:ip,0:jp-1))


		! v-momentum equation:
		Vx_mid_xt = uh_mid_xt*vh_mid_xt/h_mid_xt
		Vy_mid_yt = vh_mid_yt*vh_mid_yt/h_mid_yt*cq_s(1:ip,0:jp)
		Vy_mid_yt2 = 0.5_wp*g*h_mid_yt**2
		vh_new = vh(1:ip,1:jp) &
		  - (dt/(recqdp(0:ip-1,1:jp)))*(Vx_mid_xt(1:ip,1:jp)-Vx_mid_xt(0:ip-1,1:jp)) &
		  - (dt/(recqdq(1:ip,0:jp-1)))*(Vy_mid_yt(1:ip,1:jp)-Vy_mid_yt(1:ip,0:jp-1)) &
		  - (dt/(redq(1:ip,0:jp-1) ))* &
		  (Vy_mid_yt2(1:ip,1:jp)-Vy_mid_yt2(1:ip,0:jp-1))


        ! Full-step spherical momentum metric source.  This mirrors the
        ! original source treatment: old-time velocity with mean layer depth
        ! 0.5*(h^n+h^{n+1}).  With coriolis_scheme=1 this increment is part of
        ! the non-Coriolis RHS subsequently coupled to the CN Coriolis solve.
        if (momentum_metric_terms == 1) then
            do j=1,jp
                uh_new(:,j) = uh_new(:,j) + dt*0.5_wp* &
                    (h(1:ip,j)+h_new(:,j))*u(1:ip,j)*v(1:ip,j)*tan(theta(j))/re
                vh_new(:,j) = vh_new(:,j) - dt*0.5_wp* &
                    (h(1:ip,j)+h_new(:,j))*u(1:ip,j)*u(1:ip,j)*tan(theta(j))/re
            enddo
        endif

		! Coriolis/source corrector.  Keep the original scheme as the default
		! for backwards compatibility; scheme 1 uses a Crank-Nicolson
		! corrector for the conservative momenta.
		select case (coriolis_scheme)
		case (0)
			! Legacy/original source update.  This is intentionally kept in
			! the same form as the original model.
			uh_new=uh_new  +dt*.5_wp*(f_cor(1:ip,1:jp)*v(1:ip,1:jp) - &
				g*(hs(2:ip+1,1:jp)-hs(0:ip-1,1:jp))/(recqdp(1:ip,1:jp)+recqdp(0:ip-1,1:jp)))* &
				(h(1:ip,1:jp)+h_new)

			vh_new=vh_new  -dt*.5_wp*(f_cor(1:ip,1:jp)*u(1:ip,1:jp) + &
				g*(hs(1:ip,2:jp+1)-hs(1:ip,0:jp-1))/(redq(1:ip,1:jp)+redq(1:ip,0:jp-1)))* &
				(h(1:ip,1:jp)+h_new)

		case (1)
			! Add the orographic pressure-gradient source first.
			uh_new = uh_new - dt*0.5_wp*g* &
				(hs(2:ip+1,1:jp)-hs(0:ip-1,1:jp)) / &
				(recqdp(1:ip,1:jp)+recqdp(0:ip-1,1:jp)) * &
				(h(1:ip,1:jp)+h_new)

			vh_new = vh_new - dt*0.5_wp*g* &
				(hs(1:ip,2:jp+1)-hs(1:ip,0:jp-1)) / &
				(redq(1:ip,1:jp)+redq(1:ip,0:jp-1)) * &
				(h(1:ip,1:jp)+h_new)

			! Crank-Nicolson Coriolis corrector for the conservative momenta:
			!   uh^{n+1} = uh* + alpha (vh^n + vh^{n+1})
			!   vh^{n+1} = vh* - alpha (uh^n + uh^{n+1})
			! with alpha = f*dt/2.
			cor_alpha = 0.5_wp*dt*f_cor(1:ip,1:jp)
			cor_denom = 1._wp + cor_alpha**2

			cor_rhs_u = uh_new + cor_alpha*vh(1:ip,1:jp)
			cor_rhs_v = vh_new - cor_alpha*uh(1:ip,1:jp)

			uh_new = (cor_rhs_u + cor_alpha*cor_rhs_v) / cor_denom
			vh_new = (cor_rhs_v - cor_alpha*cor_rhs_u) / cor_denom

		case default
			write(*,*) 'ERROR: unknown coriolis_scheme = ', coriolis_scheme
			write(*,*) '       valid values are 0 (legacy) and 1 (Crank-Nicolson)'
			stop 1
		end select


		! re-calculate u and v.
		u(1:ip,1:jp) = uh_new/(h_new)
		v(1:ip,1:jp) = vh_new/h_new
		h(1:ip,1:jp) = h_new


	end subroutine lax_wendroff_ll

	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates del2 of prognostic variable
	!>@param[in] ip: number of east-west points
	!>@param[in] jp: ditto for north-south
	!>@param[in] o_halo: halos required for advection scheme
	!>@param[in] dt:  timestep
	!>@param[in] f: prognostic variable
	!>@param[inout] delsq: delsq of f
	!>@param[in] re: radius of planet
	!>@param[in] theta: latitude
	!>@param[in] thetan: latitude - staggered
	!>@param[in] dtheta: latitude step
	!>@param[in] dthetan: latitude step - staggered
	!>@param[in] phi: phi
	!>@param[in] phin: phin
	!>@param[in] dphi: dphi
	!>@param[in] dphin: dphin
	!>@param[in] recq: for efficiency
	!>@param[in] cq_s: for efficiency
	!>@param[in] dp1: for efficiency
	!>@param[in] dq: for efficiency
	!>calculates del**2:
	!>\f$ visterm = \frac{1}{re^2\cos\theta}
	!>  \frac{\partial }{\partial \theta} 
	!>\left(\cos\theta\frac{\partial f}{\partial\theta} \right) + 
	!> \frac{1}{re^2\cos^2\theta}\frac{\partial^2 f}{\partial \phi ^2}\f$
    subroutine dissipation(ip,jp,o_halo,dt,f,delsq,re,&
    		theta,thetan,dtheta,dthetan, phi, phin, dphi, dphin, &
    		recq, cq_s, dp1, dq)

		use numerics_type
		implicit none
		integer(i4b), intent(in) :: ip,jp,o_halo
		real(wp), intent(in) :: dt, re
		real(wp), dimension(1-o_halo:jp+o_halo), intent(in) :: &
											theta,thetan, dtheta, dthetan
		real(wp), dimension(1-o_halo:ip+o_halo), intent(in) :: &											
											phi, phin, dphi, dphin
		real(wp), intent(in), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: &
																				f, &
															recq, cq_s, dp1, dq
		real(wp), intent(inout), dimension(1:ip,1:jp) :: delsq
		! local variables:
		integer(i4b) :: j, i
			
		
		
		! calculate del^2 using 2nd order difference 
		! (central difference of forward and backward):
		delsq(1:ip,1:jp)  =1._wp/(re*recq(1:ip,1:jp))* &
			( cq_s(1:ip,1:jp)*(f(1:ip,2:jp+1)-f(1:ip,1:jp))/dq(1:ip,1:jp) - &
			  cq_s(1:ip,0:jp-1)*(f(1:ip,1:jp)-f(1:ip,0:jp-1))/dq(1:ip,0:jp-1) ) / &
			  dq(1:ip,1:jp)
			  
		delsq(1:ip,1:jp)  = delsq(1:ip,1:jp) + &
			1._wp/(recq(1:ip,1:jp)**2._wp)* &
			( (f(2:ip+1,1:jp)-f(1:ip,1:jp))/dp1(1:ip,1:jp) - &
			  (f(1:ip,1:jp)-f(0:ip-1,1:jp))/dp1(0:ip-1,1:jp) ) / &
			  dp1(1:ip,1:jp)

	end subroutine dissipation

	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates smagorinsky-lilly viscosity
	!>@param[in] ip: number of east-west points
	!>@param[in] jp: ditto for north-south
	!>@param[in] o_halo: halos required for advection scheme
	!>@param[in] cvis:  coefficient for viscosity
	!>@param[in] u,v: u and v winds
	!>@param[inout] vis: viscosity
	!>@param[in] re: radius of planet
	!>@param[in] recq: for efficiency
	!>@param[in] dp1: for efficiency
	!>@param[in] dq: for efficiency
	!>calculates smagorinsky-lilly viscosity:
	!>\f$ visco = C_s^2\Delta x\Delta y|S|\f$
    subroutine smagorinsky(ip,jp,o_halo,cvis,u,v,vis,re,&
    		recq, dp1, dq)

		use numerics_type
		implicit none
		integer(i4b), intent(in) :: ip,jp,o_halo
		real(wp), intent(in) :: cvis, re
		real(wp), intent(in), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: &
																		u,v, &
															recq, dp1, dq
		real(wp), intent(inout), dimension(1:ip,1:jp) :: vis
		! local variables:
		integer(i4b) :: j, i
			
		
		
		! calculate viscosity using centred differences:
		vis(1:ip,1:jp) = cvis**2._wp*re*recq(1:ip,1:jp)*dp1(1:ip,1:jp)*dq(1:ip,1:jp)* &
		sqrt( ( (u(2:ip+1,1:jp)-u(0:ip-1,1:jp))/ &
			(recq(1:ip,1:jp)*(dp1(0:ip-1,1:jp)+dp1(1:ip,1:jp))) )**2._wp + &
			( (v(1:ip,2:jp+1)-v(1:ip,0:jp-1))/ &
			(re*(dq(1:ip,1:jp)+dq(1:ip,0:jp-1))) )**2._wp + &
			0.5_wp*(  &
			(u(1:ip,2:jp+1)-u(1:ip,0:jp-1))/ &
			(re*(dq(1:ip,1:jp)+dq(1:ip,0:jp-1)))+ &
			(v(2:ip+1,1:jp)-v(0:ip-1,1:jp))/ &
			(recq(1:ip,1:jp)*(dp1(1:ip,1:jp)+dp1(0:ip-1,1:jp))) &
			)**2._wp )

	end subroutine smagorinsky


	!>@author
	!>Paul J. Connolly / spherical SGS update
	!>@brief
	!>Construct a thickness-weighted spherical Smagorinsky stress tensor.
	!>
	!>The physical east/north strain components are
	!>  S_ee = (1/(a cos(theta))) du/dlambda - v tan(theta)/a
	!>  S_nn = (1/a) dv/dtheta
	!>  S_en = 0.5[(1/a) du/dtheta + (1/(a cos(theta))) dv/dlambda
	!>                 + u tan(theta)/a].
	!>Only the deviatoric part is used in the SGS stress, so resolved
	!>horizontal divergence is not treated as a bulk-viscous deformation.
	!>The stress is tau_ij = 2 h nu_t S'_ij, with
	!>nu_t = cvis^2 Delta_x Delta_y |S'| and |S'|=sqrt(2 S':S').
    subroutine smagorinsky_spherical_stress(ip,jp,o_halo,cvis,h,u,v,&
            tau_uu,tau_uv,tau_vv,vis,re,theta,recq,dp1,dq)

        use numerics_type
        implicit none
        integer(i4b), intent(in) :: ip,jp,o_halo
        real(wp), intent(in) :: cvis, re
        real(wp), dimension(1-o_halo:jp+o_halo), intent(in) :: theta
        real(wp), intent(in), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: &
                h,u,v,recq,dp1,dq
        real(wp), intent(inout), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: &
                tau_uu,tau_uv,tau_vv
        real(wp), intent(inout), dimension(1:ip,1:jp) :: vis

        real(wp), dimension(1:ip) :: dudx,dvdy,dudy,dvdx, &
                see,snn,sen,sdev_ee,sdev_nn,strain_mag
        real(wp) :: metric_curv
        integer(i4b) :: j

        ! Initialise halos too: MPI exchange fills internal-neighbour halos;
        ! physical latitude halos are set explicitly in the driver.
        tau_uu = 0._wp
        tau_uv = 0._wp
        tau_vv = 0._wp
        vis = 0._wp

        do j=1,jp
            dudx = (u(2:ip+1,j)-u(0:ip-1,j)) / &
                    (recq(1:ip,j)*(dp1(1:ip,j)+dp1(0:ip-1,j)))
            dvdx = (v(2:ip+1,j)-v(0:ip-1,j)) / &
                    (recq(1:ip,j)*(dp1(1:ip,j)+dp1(0:ip-1,j)))
            dudy = (u(1:ip,j+1)-u(1:ip,j-1)) / &
                    (re*(dq(1:ip,j)+dq(1:ip,j-1)))
            dvdy = (v(1:ip,j+1)-v(1:ip,j-1)) / &
                    (re*(dq(1:ip,j)+dq(1:ip,j-1)))

            metric_curv = tan(theta(j))/re
            see = dudx - v(1:ip,j)*metric_curv
            snn = dvdy
            sen = 0.5_wp*(dudy + dvdx + u(1:ip,j)*metric_curv)

            ! In two horizontal dimensions, removing half the trace gives
            ! S'_ee = (S_ee-S_nn)/2 and S'_nn=-S'_ee.
            sdev_ee = 0.5_wp*(see-snn)
            sdev_nn = -sdev_ee
            strain_mag = sqrt(2._wp*sdev_ee**2 + 2._wp*sdev_nn**2 + &
                              4._wp*sen**2)

            ! Delta^2 = Delta_x Delta_y on the local spherical grid.
            vis(1:ip,j) = cvis**2 * re*recq(1:ip,j)*dp1(1:ip,j)*dq(1:ip,j) * &
                          strain_mag

            tau_uu(1:ip,j) = 2._wp*h(1:ip,j)*vis(1:ip,j)*sdev_ee
            tau_vv(1:ip,j) = 2._wp*h(1:ip,j)*vis(1:ip,j)*sdev_nn
            tau_uv(1:ip,j) = 2._wp*h(1:ip,j)*vis(1:ip,j)*sen
        enddo

    end subroutine smagorinsky_spherical_stress


	!>@brief
	!>Spherical divergence of a symmetric horizontal SGS stress tensor.
	!>Returns tendencies of conservative momenta hu and hv.  Written in a
	!>face-flux form for the derivative pieces, with the required spherical
	!>tensor metric term retained in the northward component.
    subroutine spherical_stress_divergence(ip,jp,o_halo,tau_uu,tau_uv,tau_vv,&
            tend_u,tend_v,re,theta,thetan,recq,cq,cq_s,dp1,dq)

        use numerics_type
        implicit none
        integer(i4b), intent(in) :: ip,jp,o_halo
        real(wp), intent(in) :: re
        real(wp), dimension(1-o_halo:jp+o_halo), intent(in) :: theta,thetan
        real(wp), intent(in), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: &
                tau_uu,tau_uv,tau_vv,recq,cq,cq_s,dp1,dq
        real(wp), intent(out), dimension(1:ip,1:jp) :: tend_u,tend_v

        real(wp), dimension(1:ip) :: uu_e,uu_w,uv_e,uv_w, &
                uv_n,uv_s,vv_n,vv_s
        integer(i4b) :: j

        do j=1,jp
            ! Arithmetic face stresses give a centred, conservative flux
            ! difference while requiring only the model's existing one-cell halo.
            uu_e = 0.5_wp*(tau_uu(1:ip,j)+tau_uu(2:ip+1,j))
            uu_w = 0.5_wp*(tau_uu(0:ip-1,j)+tau_uu(1:ip,j))
            uv_e = 0.5_wp*(tau_uv(1:ip,j)+tau_uv(2:ip+1,j))
            uv_w = 0.5_wp*(tau_uv(0:ip-1,j)+tau_uv(1:ip,j))
            uv_n = 0.5_wp*(tau_uv(1:ip,j)+tau_uv(1:ip,j+1))
            uv_s = 0.5_wp*(tau_uv(1:ip,j-1)+tau_uv(1:ip,j))
            vv_n = 0.5_wp*(tau_vv(1:ip,j)+tau_vv(1:ip,j+1))
            vv_s = 0.5_wp*(tau_vv(1:ip,j-1)+tau_vv(1:ip,j))

            ! (div tau)_east = 1/(a cos phi) d tau_ee/dlambda
            !                + 1/(a cos^2 phi) d(tau_en cos^2 phi)/dphi
            tend_u(1:ip,j) = (uu_e-uu_w)/(recq(1:ip,j)*dp1(1:ip,j)) + &
                (uv_n*cq_s(1:ip,j)**2 - uv_s*cq_s(1:ip,j-1)**2) / &
                (re*cq(1:ip,j)**2*dq(1:ip,j))

            ! (div tau)_north = 1/(a cos phi) d tau_en/dlambda
            !                 + 1/(a cos phi) d(tau_nn cos phi)/dphi
            !                 + tau_ee tan(phi)/a
            tend_v(1:ip,j) = (uv_e-uv_w)/(recq(1:ip,j)*dp1(1:ip,j)) + &
                (vv_n*cq_s(1:ip,j) - vv_s*cq_s(1:ip,j-1)) / &
                (re*cq(1:ip,j)*dq(1:ip,j)) + &
                tau_uu(1:ip,j)*tan(theta(j))/re
        enddo

    end subroutine spherical_stress_divergence

	end module advection
