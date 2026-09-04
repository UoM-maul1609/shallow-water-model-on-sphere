	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>drivers for the shallow water model
    module drivers
    use numerics_type
    !use variables
    private
    public :: model_driver
    contains
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calls IO and runs one time-step of model
	!>@param[in] ip: number of east-west levels on global grid
	!>@param[in] ipp: number of east-west levels on this PE
	!>@param[in] jp: ditto for south-north
	!>@param[in] jpp: ditto for south-north
	!>@param[in] ntim: number of time-levels
	!>@param[in] f: rotation rate
	!>@param[in] re: radius of planet
	!>@param[in] g: gravity
	!>@param[in] rho: density of fluid
	!>@param[in] dphi: step in longitude
	!>@param[in] dtheta: step in latitude
	!>@param[in] dphin: step in longitude - staggered
	!>@param[in] dthetan: step in latitude
	!>@param[in] f_cor: coriolis parameter
	!>@param[inout] h: depth of fluid
	!>@param[in] hs: height of surface above reference
	!>@param[inout] u,v: winds
	!>@param[inout] height: height of fluid
	!>@param[in] dt, dx,dy, x, y: grids
	!>@param[in] phi, theta, phin, thetan: grids and staggered grids
	!>@param[in] recqdp - for efficiency
	!>@param[in] recqdp_s - for efficiency
	!>@param[in] recqdq_s - for efficiency
	!>@param[in] redq_s - for efficiency
	!>@param[in] redq - for efficiency
	!>@param[in] recq - for efficiency
	!>@param[in] cq_s - for efficiency
	!>@param[in] cq - for efficiency
	!>@param[in] dp1 - for efficiency
	!>@param[in] dq - for efficiency
	!>@param[in] recqdq - for efficiency
	!>@param[in] u_nudge: wind to nudge to
	!>@param[inout] new_file: flag for if this is a new file
	!>@param[in] outputfile: netcdf output
	!>@param[in] output_interval: interval for output (s)
	!>@param[in] nudge: logical if we want to nudge
	!>@param[in] nudge_tau: time-scale (s) for nudging
	!>@param[in] subgrid_model - 1 is constant viscosity, 2 is smagorinsky approach
	!>@param[in] viscous_dissipation: add dissipation term
	!>@param[in] dissipate_h: add dissipation term to h-field
	!>@param[in] vis: viscosity
	!>@param[in] cvis: smagorinsky parameter
	!>@param[in] vis_eq: viscosity in equatorial region (needed for stability)
	!>@param[in] lat_eq: latitude north and south over which to apply vis_eq
	!>@param[in] dims,id, world_process, ring_comm: mpi variables
    subroutine model_driver(ip,ipp, jp,jpp, ntim, f, &
				re, g, rho, dphi, dtheta, dphin, dthetan, &
				f_cor,h,hs, u, v, &
				height, dt, dx, dy, x, y, &
				phi, theta, phin, thetan, &
				recqdp, recqdp_s, recqdq_s, redq_s, redq, &
    			recq, cq_s, cq, dp1, dq,recqdq, &
			    u_nudge,o_halo, &
				ipstart, jpstart, coords, &
				new_file,outputfile, output_interval, nudge, nudge_tau, &
				subgrid_model, viscous_dissipation, dissipate_h,vis, cvis, &
				vis_eq, lat_eq, coriolis_scheme, momentum_metric_terms, smagorinsky_scheme, &
                lat_boundary_scheme, sponge_south_width, sponge_north_width, &
                sponge_south_timescale, sponge_north_timescale, slat, nlat, &
				dims,id, world_process, rank, ring_comm)
		use numerics_type
		use mpi_module
		use advection

		implicit none
		logical, intent(inout) :: new_file
		logical, intent(in) :: nudge, viscous_dissipation, dissipate_h
		integer(i4b), intent(in) :: ip,ipp, jp,jpp, ntim, o_halo, ipstart, jpstart, &
									subgrid_model, coriolis_scheme, momentum_metric_terms, smagorinsky_scheme, &
                                    lat_boundary_scheme
		integer(i4b), intent(in) :: id, world_process, ring_comm, rank
		integer(i4b), dimension(2), intent(in) :: coords, dims
		character (len=*), intent(in) :: outputfile
		real(wp), intent(in) :: f, re, g, rho, dt, output_interval
		real(wp), dimension(1-o_halo:ipp+o_halo), intent(in) :: phi, phin, dphi, dphin
		real(wp), dimension(1-o_halo:jpp+o_halo), intent(in) :: theta, thetan, u_nudge, &
																dtheta, dthetan
		real(wp), dimension(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), &
					intent(in) :: f_cor, &
    				recqdp, recqdp_s, recqdq_s, redq_s, redq, &
    				recq, cq_s, cq, dp1, dq, recqdq
		real(wp), dimension(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), &
					intent(inout) :: h, hs, u, v, height
		real(wp), dimension(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), &
					intent(in) :: dx, dy, x, y
		real(wp), intent(in) :: vis, nudge_tau, cvis, lat_eq, vis_eq, &
                            sponge_south_width, sponge_north_width, &
                            sponge_south_timescale, sponge_north_timescale, slat, nlat
					
		! locals:		
		integer(i4b) :: n, cur=1, j, error, rank2
        logical :: south_wall, north_wall
		real(wp) :: time, time_last_output, output_time
		real(wp), dimension(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo) :: &
				u_old, v_old, h_old, u_sgs, v_sgs, h_sgs, &
                tau_uu, tau_uv, tau_vv, u_bc_ref, eta_bc_ref
		real(wp), dimension(1:ipp,1:jpp) :: delsq, vort, visco, &
                sgs_mom_u, sgs_mom_v, mom_u_tmp, mom_v_tmp, &
                h_sponge_ref, u_sponge_ref, v_sponge_ref
		

		time_last_output=-output_interval
		output_time=output_interval
		rank2=dims(1)*dims(2)

        ! Latitude boundary setup.  Scheme 0 deliberately preserves the
        ! legacy frozen physical ghost cells.  Scheme 1 uses reflecting
        ! free-slip ghost cells and exact zero normal mass flux at the
        ! Lax-Wendroff latitude faces.
        select case (lat_boundary_scheme)
        case (0)
            south_wall=.false.
            north_wall=.false.
        case (1)
            south_wall=(coords(2) == 0)
            north_wall=(coords(2) == dims(2)-1)
        case default
            write(*,*) 'ERROR: unknown lat_boundary_scheme = ',lat_boundary_scheme
            write(*,*) '       valid values are 0 (legacy) and 1 (free-slip)'
            stop 1
        end select

        if (lat_boundary_scheme == 1) then
            if (sponge_south_width < 0._wp .or. sponge_north_width < 0._wp .or. &
                sponge_south_timescale < 0._wp .or. sponge_north_timescale < 0._wp) then
                write(*,*) 'ERROR: sponge widths/timescales must be >= 0'
                stop 1
            endif

            ! Save the initial model state as the sponge reference.
            h_sponge_ref=h(1:ipp,1:jpp)
            u_sponge_ref=u(1:ipp,1:jpp)
            v_sponge_ref=v(1:ipp,1:jpp)

            ! Ensure the first diagnostic/advection call sees free-slip halos.
            call exchange_halos(ring_comm, id, ipp, jpp, o_halo, h)
            call exchange_halos(ring_comm, id, ipp, jpp, o_halo, u)
            call exchange_halos(ring_comm, id, ipp, jpp, o_halo, v)
            call exchange_halos(ring_comm, id, ipp, jpp, o_halo, hs)
            ! Construct the balanced wall once at t=0, then freeze that
            ! background continuation as the reference state. Subsequent
            ! timesteps reflect perturbations about this reference rather than
            ! recomputing an instantaneous gradient-wind slope from evolving u.
            call apply_balanced_free_slip_lat_halos(ipp,jpp,o_halo,h,hs,u,v,theta,f_cor, &
                re,g,momentum_metric_terms,coords,dims)
            u_bc_ref = u
            eta_bc_ref = h + hs
        endif

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! time-loop                                                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		do n=1,ntim	
		
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! write netcdf variables                                                     !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			time=real(n-1,wp)*dt
			if (time-time_last_output >= output_interval) then
				if (id==world_process) &
					print *,'output no ',cur,' at time (hrs) ', &
						time/3600._wp,n,' steps of ',ntim
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! calculate diagnostics for output: vorticity, etc                       !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				call diagnostics(ipp,jpp,o_halo,dt,u,v, vort,re,&
						theta,thetan,dtheta,dthetan, phi, phin, dphi, dphin, &
						recq, cq_s, cq, dp1, dq)
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				call output(new_file,outputfile,cur,ip,ipp,ipstart,jp,jpp,jpstart, &
							o_halo, &
							time,phi,theta, &
							u_nudge, f_cor, height, h, u, v, vort, &
							id, world_process, rank2, ring_comm)
				time_last_output=time
				cur=cur+1
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			
			
			

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! advance solution 1 time-step                                               !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			h_old=h
			u_old=u
			v_old=v
			call lax_wendroff_ll(ipp,jpp,o_halo,dt,g,u,v,h,hs,re,&
	    		theta,thetan,dtheta,dthetan, phi, phin, dphi, dphin, f_cor, &
    			recqdq, recqdp, recqdp_s, recqdq_s, redq_s, redq, cq, cq_s, &
                coriolis_scheme, momentum_metric_terms, south_wall, north_wall)	    		
! 			call lax_wendroff_sphere(ipp,jpp,o_halo,dt,dx,dy,g,u,v,h,hs,re,theta,f_cor)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! nudge                                                                      !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if (nudge) then
				do j=1,jpp
					! mid-point rule:
					u(1:ipp,j)=u(1:ipp,j)+ &
						(u_nudge(j)- &
						(0.5_wp*(u(1:ipp,j)+u_old(1:ipp,j)))/real(1,wp) ) &
						/nudge_tau * dt
! 					v(1:ipp,j)=v(1:ipp,j)+&
! 						(0._wp- &
! 						sum(0.5_wp*(v(1:ipp,j)+v_old(1:ipp,j)))/real(ipp,wp) ) &
! 						/nudge_tau *dt

! 					Derived by integrating du/dt=(u_nudge-u)/tau
! 					u(1:ipp,j)=u_nudge(j)- &
! 						(u_nudge(j)- &
! 						(u(1:ipp,j))/real(1,wp) ) * &
! 						exp(-dt/nudge_tau )
! 					v(1:ipp,j)=0._wp- &
! 						(0._wp- &
! 						(v(1:ipp,j))/real(1,wp) ) * &
! 						exp(-dt/nudge_tau )
					
				enddo
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! calculate dissipation: mid-point rule                                      !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if (viscous_dissipation) then
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! halo exchanges                                                         !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				call exchange_halos(ring_comm, id, ipp, jpp, o_halo, u)
				call exchange_halos(ring_comm, id, ipp, jpp, o_halo, v)
                if (lat_boundary_scheme == 1) then
                    call apply_reference_vector_free_slip_lat_halos(ipp,jpp,o_halo,u,v,u_bc_ref,coords,dims)
                endif
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				select case(subgrid_model)
				case (1)
					! Legacy constant scalar viscosity.
					call dissipation(ipp,jpp,o_halo,dt,0.5_wp*(u_old+u), delsq,re,&
						theta,thetan,dtheta,dthetan, phi, phin, dphi, dphin, &
						recq, cq_s, dp1, dq)
					u(1:ipp,1:jpp)=u(1:ipp,1:jpp)+dt*delsq*vis

					call dissipation(ipp,jpp,o_halo,dt,0.5_wp*(v_old+v), delsq,re,&
						theta,thetan,dtheta,dthetan, phi, phin, dphi, dphin, &
						recq, cq_s, dp1, dq)
					v(1:ipp,1:jpp)=v(1:ipp,1:jpp)+dt*delsq*vis

				case (2)
					select case (smagorinsky_scheme)
					case (0)
						! Legacy Smagorinsky: scalar Laplacian of u and v.
						call dissipation(ipp,jpp,o_halo,dt,0.5_wp*(u_old+u), delsq,re,&
							theta,thetan,dtheta,dthetan, phi, phin, dphi, dphin, &
							recq, cq_s, dp1, dq)
						call smagorinsky(ipp,jpp,o_halo,cvis,0.5_wp*(u_old+u),&
							0.5_wp*(v_old+v),visco,re,recq, dp1, dq)
						u(1:ipp,1:jpp)=u(1:ipp,1:jpp)+dt*delsq*visco

						call dissipation(ipp,jpp,o_halo,dt,0.5_wp*(v_old+v), delsq,re,&
							theta,thetan,dtheta,dthetan, phi, phin, dphi, dphin, &
							recq, cq_s, dp1, dq)
						v(1:ipp,1:jpp)=v(1:ipp,1:jpp)+dt*delsq*visco

					case (1)
						! Spherical, thickness-weighted, conservative SGS stress.
						! Evaluate h, u, v at the old/new midpoint, matching the
						! time-centred treatment used by the legacy viscosity.
						call exchange_halos(ring_comm, id, ipp, jpp, o_halo, h)
                        if (lat_boundary_scheme == 1) then
                            call apply_reference_height_lat_halos(ipp,jpp,o_halo,h,hs,eta_bc_ref,coords,dims)
                        endif
						u_sgs = 0.5_wp*(u_old+u)
						v_sgs = 0.5_wp*(v_old+v)
						h_sgs = 0.5_wp*(h_old+h)

						call smagorinsky_spherical_stress(ipp,jpp,o_halo,cvis,h_sgs,u_sgs,v_sgs,&
							tau_uu,tau_uv,tau_vv,visco,re,theta,recq,dp1,dq)

						! Exchange the SGS stresses.  At the physical latitude edges
						! use a zero-normal-SGS-stress boundary condition: making the
						! ghost stress the negative of the adjacent cell makes the
						! face-centred stress exactly zero.
						call exchange_halos(ring_comm, id, ipp, jpp, o_halo, tau_uu)
						call exchange_halos(ring_comm, id, ipp, jpp, o_halo, tau_uv)
						call exchange_halos(ring_comm, id, ipp, jpp, o_halo, tau_vv)
						if (coords(2) == 0) then
							tau_uu(1:ipp,0) = -tau_uu(1:ipp,1)
							tau_uv(1:ipp,0) = -tau_uv(1:ipp,1)
							tau_vv(1:ipp,0) = -tau_vv(1:ipp,1)
						endif
						if (coords(2) == dims(2)-1) then
							tau_uu(1:ipp,jpp+1) = -tau_uu(1:ipp,jpp)
							tau_uv(1:ipp,jpp+1) = -tau_uv(1:ipp,jpp)
							tau_vv(1:ipp,jpp+1) = -tau_vv(1:ipp,jpp)
						endif

						call spherical_stress_divergence(ipp,jpp,o_halo,tau_uu,tau_uv,tau_vv,&
							sgs_mom_u,sgs_mom_v,re,theta,thetan,recq,cq,cq_s,dp1,dq)

						! Apply the SGS term to conservative momenta.  h itself is
						! not diffused by this closure.
						mom_u_tmp = h(1:ipp,1:jpp)*u(1:ipp,1:jpp) + dt*sgs_mom_u
						mom_v_tmp = h(1:ipp,1:jpp)*v(1:ipp,1:jpp) + dt*sgs_mom_v
						u(1:ipp,1:jpp) = mom_u_tmp/h(1:ipp,1:jpp)
						v(1:ipp,1:jpp) = mom_v_tmp/h(1:ipp,1:jpp)

					case default
						write(*,*) 'ERROR: unknown smagorinsky_scheme = ', smagorinsky_scheme
						write(*,*) '       valid values are 0 (legacy) and 1 (spherical stress)'
						stop 1
					end select

				case default
					print *,'error subgrid ',subgrid_model
				end select

				! Existing extra equatorial v viscosity.  In the new SGS path
				! delsq has not otherwise been evaluated for v, so do it here.
				if ((subgrid_model == 2) .and. (smagorinsky_scheme == 1) .and. &
					(abs(vis_eq) > tiny(1._wp))) then
					call exchange_halos(ring_comm, id, ipp, jpp, o_halo, v)
					call dissipation(ipp,jpp,o_halo,dt,0.5_wp*(v_old+v), delsq,re,&
						theta,thetan,dtheta,dthetan, phi, phin, dphi, dphin, &
						recq, cq_s, dp1, dq)
				endif

				do j=1,jpp
					if((theta(j) >-lat_eq*pi/180._wp) .and. &
						 (theta(j) < lat_eq*pi/180._wp)) then

						v(1:ipp,j)=v(1:ipp,j)+dt*delsq(1:ipp,j)*vis_eq* &
							cos(theta(j)*90._wp/lat_eq)
					endif
				enddo

				! Preserve the existing optional h diffusion for the constant-
				! viscosity model only.  The conservative Smagorinsky closure
				! acts on momentum, not layer thickness.
				if (dissipate_h .and. (subgrid_model == 1)) then
					call exchange_halos(ring_comm, id, ipp, jpp, o_halo, h)
                    if (lat_boundary_scheme == 1) then
                        call apply_even_lat_halos(ipp,jpp,o_halo,h,coords,dims)
                    endif
					call dissipation(ipp,jpp,o_halo,dt,0.5_wp*(h_old+h), delsq,re,&
						theta,thetan,dtheta,dthetan, phi, phin, dphi, dphin, &
						recq, cq_s, dp1, dq)
					h(1:ipp,1:jpp)=h(1:ipp,1:jpp)+dt*delsq*vis
				endif

			endif	    	
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



            ! Optional latitude sponge. A zero width or zero timescale disables
            ! that side. The damping ramp is quadratic in distance into the zone.
            if (lat_boundary_scheme == 1) then
                call apply_latitude_sponge(ipp,jpp,o_halo,dt,h,u,v, &
                    h_sponge_ref,u_sponge_ref,v_sponge_ref,theta,slat,nlat, &
                    sponge_south_width,sponge_north_width, &
                    sponge_south_timescale,sponge_north_timescale)
            endif

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! halo exchanges                                                             !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			call exchange_halos(ring_comm, id, ipp, jpp, o_halo, h)
			call exchange_halos(ring_comm, id, ipp, jpp, o_halo, u)
			call exchange_halos(ring_comm, id, ipp, jpp, o_halo, v)
            if (lat_boundary_scheme == 1) then
                call apply_reference_free_slip_lat_halos(ipp,jpp,o_halo,h,hs,u,v, &
                    u_bc_ref,eta_bc_ref,coords,dims)
            endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
! 			if(coords(2)==(dims(2)-1)) then
! 				v(:,jpp+1:jpp+o_halo)=0._wp
! 				do j=jpp+1,jpp+o_halo
! 					u(:,j)=u_nudge(j)
! 				enddo
! 				h(:,jpp+o_halo)=h(:,jpp)
! 			endif
! 			if(coords(2)==0) then
! 				v(:,1-o_halo)=0._wp
! 				do j=1-o_halo,0
! 					u(:,j)=u_nudge(j)
! 				enddo
! 				h(:,1-o_halo)=h(:,1)
! 			endif
			

		enddo
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








		
	end subroutine model_driver
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	

    ! Balanced impermeable free-slip latitude ghosts.  Tangential velocity
    ! is reflected evenly, normal velocity oddly, and the total free-surface
    ! height eta=h+hs is extrapolated with the local geostrophic/gradient-wind
    ! slope.  This avoids imposing the unphysical d(h+hs)/dtheta=0 condition
    ! at a rotating free-slip wall.
    subroutine apply_balanced_free_slip_lat_halos(ip,jp,o_halo,h,hs,u,v,theta,f_cor, &
                                                   re,g,momentum_metric_terms,coords,dims)
        use numerics_type
        implicit none
        integer(i4b), intent(in) :: ip,jp,o_halo,momentum_metric_terms
        integer(i4b), dimension(2), intent(in) :: coords,dims
        real(wp), intent(in) :: re,g
        real(wp), dimension(1-o_halo:jp+o_halo), intent(in) :: theta
        real(wp), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo), intent(in) :: f_cor
        real(wp), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo), intent(inout) :: h,hs,u,v
        integer(i4b) :: k, ji, jg
        real(wp) :: theta_face, dtheta_pair
        real(wp), dimension(1-o_halo:ip+o_halo) :: f_face, u_face, deta_dtheta, eta_i

        if (momentum_metric_terms /= 0 .and. momentum_metric_terms /= 1) then
            write(*,*) 'ERROR: unknown momentum_metric_terms in latitude BC = ',momentum_metric_terms
            stop 1
        endif

        if (coords(2) == 0) then
            do k=1,o_halo
                ji=k
                jg=1-k
                ! Free-slip velocity reflection.
                u(:,jg)=u(:,ji)
                v(:,jg)=-v(:,ji)
                ! Continue static topography evenly across the artificial wall.
                hs(:,jg)=hs(:,ji)

                theta_face=0.5_wp*(theta(jg)+theta(ji))
                dtheta_pair=theta(jg)-theta(ji)
                f_face=0.5_wp*(f_cor(:,jg)+f_cor(:,ji))
                u_face=u(:,ji)
                deta_dtheta=-(re*f_face*u_face)/g
                if (momentum_metric_terms == 1) then
                    deta_dtheta=deta_dtheta-(u_face*u_face*tan(theta_face))/g
                endif
                eta_i=h(:,ji)+hs(:,ji)
                h(:,jg)=eta_i+deta_dtheta*dtheta_pair-hs(:,jg)
            enddo
        endif

        if (coords(2) == dims(2)-1) then
            do k=1,o_halo
                ji=jp+1-k
                jg=jp+k
                u(:,jg)=u(:,ji)
                v(:,jg)=-v(:,ji)
                hs(:,jg)=hs(:,ji)

                theta_face=0.5_wp*(theta(jg)+theta(ji))
                dtheta_pair=theta(jg)-theta(ji)
                f_face=0.5_wp*(f_cor(:,jg)+f_cor(:,ji))
                u_face=u(:,ji)
                deta_dtheta=-(re*f_face*u_face)/g
                if (momentum_metric_terms == 1) then
                    deta_dtheta=deta_dtheta-(u_face*u_face*tan(theta_face))/g
                endif
                eta_i=h(:,ji)+hs(:,ji)
                h(:,jg)=eta_i+deta_dtheta*dtheta_pair-hs(:,jg)
            enddo
        endif
    end subroutine apply_balanced_free_slip_lat_halos

    ! Reference-state free-slip latitude ghosts.  The initial balanced
    ! background continuation is kept fixed, while perturbations of tangential
    ! velocity and free-surface height are reflected evenly.  Normal velocity
    ! is reflected oddly.  This avoids imposing instantaneous gradient-wind
    ! balance on evolving disturbances at a rigid wall.
    subroutine apply_reference_free_slip_lat_halos(ip,jp,o_halo,h,hs,u,v, &
                                                    u_ref,eta_ref,coords,dims)
        use numerics_type
        implicit none
        integer(i4b), intent(in) :: ip,jp,o_halo
        integer(i4b), dimension(2), intent(in) :: coords,dims
        real(wp), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo), intent(in) :: &
            u_ref,eta_ref
        real(wp), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo), intent(inout) :: &
            h,hs,u,v
        integer(i4b) :: k,ji,jg

        if (coords(2) == 0) then
            do k=1,o_halo
                ji=k
                jg=1-k
                u(:,jg)=u_ref(:,jg) + (u(:,ji)-u_ref(:,ji))
                v(:,jg)=-v(:,ji)
                h(:,jg)=eta_ref(:,jg) + ((h(:,ji)+hs(:,ji))-eta_ref(:,ji)) - hs(:,jg)
            enddo
        endif
        if (coords(2) == dims(2)-1) then
            do k=1,o_halo
                ji=jp+1-k
                jg=jp+k
                u(:,jg)=u_ref(:,jg) + (u(:,ji)-u_ref(:,ji))
                v(:,jg)=-v(:,ji)
                h(:,jg)=eta_ref(:,jg) + ((h(:,ji)+hs(:,ji))-eta_ref(:,ji)) - hs(:,jg)
            enddo
        endif
    end subroutine apply_reference_free_slip_lat_halos

    subroutine apply_reference_vector_free_slip_lat_halos(ip,jp,o_halo,u,v,u_ref,coords,dims)
        use numerics_type
        implicit none
        integer(i4b), intent(in) :: ip,jp,o_halo
        integer(i4b), dimension(2), intent(in) :: coords,dims
        real(wp), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo), intent(in) :: u_ref
        real(wp), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo), intent(inout) :: u,v
        integer(i4b) :: k,ji,jg
        if (coords(2) == 0) then
            do k=1,o_halo
                ji=k; jg=1-k
                u(:,jg)=u_ref(:,jg) + (u(:,ji)-u_ref(:,ji))
                v(:,jg)=-v(:,ji)
            enddo
        endif
        if (coords(2) == dims(2)-1) then
            do k=1,o_halo
                ji=jp+1-k; jg=jp+k
                u(:,jg)=u_ref(:,jg) + (u(:,ji)-u_ref(:,ji))
                v(:,jg)=-v(:,ji)
            enddo
        endif
    end subroutine apply_reference_vector_free_slip_lat_halos

    subroutine apply_reference_height_lat_halos(ip,jp,o_halo,h,hs,eta_ref,coords,dims)
        use numerics_type
        implicit none
        integer(i4b), intent(in) :: ip,jp,o_halo
        integer(i4b), dimension(2), intent(in) :: coords,dims
        real(wp), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo), intent(in) :: hs,eta_ref
        real(wp), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo), intent(inout) :: h
        integer(i4b) :: k,ji,jg
        if (coords(2) == 0) then
            do k=1,o_halo
                ji=k; jg=1-k
                h(:,jg)=eta_ref(:,jg) + ((h(:,ji)+hs(:,ji))-eta_ref(:,ji)) - hs(:,jg)
            enddo
        endif
        if (coords(2) == dims(2)-1) then
            do k=1,o_halo
                ji=jp+1-k; jg=jp+k
                h(:,jg)=eta_ref(:,jg) + ((h(:,ji)+hs(:,ji))-eta_ref(:,ji)) - hs(:,jg)
            enddo
        endif
    end subroutine apply_reference_height_lat_halos

    subroutine apply_vector_free_slip_lat_halos(ip,jp,o_halo,u,v,coords,dims)
        use numerics_type
        implicit none
        integer(i4b), intent(in) :: ip,jp,o_halo
        integer(i4b), dimension(2), intent(in) :: coords,dims
        real(wp), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo), intent(inout) :: u,v
        integer(i4b) :: k
        if (coords(2) == 0) then
            do k=1,o_halo
                u(:,1-k)=u(:,k)
                v(:,1-k)=-v(:,k)
            enddo
        endif
        if (coords(2) == dims(2)-1) then
            do k=1,o_halo
                u(:,jp+k)=u(:,jp+1-k)
                v(:,jp+k)=-v(:,jp+1-k)
            enddo
        endif
    end subroutine apply_vector_free_slip_lat_halos

    subroutine apply_even_lat_halos(ip,jp,o_halo,a,coords,dims)
        use numerics_type
        implicit none
        integer(i4b), intent(in) :: ip,jp,o_halo
        integer(i4b), dimension(2), intent(in) :: coords,dims
        real(wp), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo), intent(inout) :: a
        integer(i4b) :: k
        if (coords(2) == 0) then
            do k=1,o_halo
                a(:,1-k)=a(:,k)
            enddo
        endif
        if (coords(2) == dims(2)-1) then
            do k=1,o_halo
                a(:,jp+k)=a(:,jp+1-k)
            enddo
        endif
    end subroutine apply_even_lat_halos

    ! Exact exponential relaxation toward the initial state. Widths are degrees
    ! latitude and taus are the e-folding seconds at the physical boundary.
    subroutine apply_latitude_sponge(ip,jp,o_halo,dt,h,u,v,h_ref,u_ref,v_ref,theta, &
            slat,nlat,south_width,north_width,south_tau,north_tau)
        use numerics_type
        implicit none
        integer(i4b), intent(in) :: ip,jp,o_halo
        real(wp), intent(in) :: dt,slat,nlat,south_width,north_width,south_tau,north_tau
        real(wp), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo), intent(inout) :: h,u,v
        real(wp), dimension(1:ip,1:jp), intent(in) :: h_ref,u_ref,v_ref
        real(wp), dimension(1-o_halo:jp+o_halo), intent(in) :: theta
        integer(i4b) :: j
        real(wp) :: lat_deg,xi,rate,damp

        do j=1,jp
            lat_deg=theta(j)*180._wp/pi
            rate=0._wp
            if (south_width > 0._wp .and. south_tau > 0._wp) then
                if (lat_deg < slat+south_width) then
                    xi=max(0._wp,min(1._wp,(slat+south_width-lat_deg)/south_width))
                    rate=rate+xi*xi/south_tau
                endif
            endif
            if (north_width > 0._wp .and. north_tau > 0._wp) then
                if (lat_deg > nlat-north_width) then
                    xi=max(0._wp,min(1._wp,(lat_deg-(nlat-north_width))/north_width))
                    rate=rate+xi*xi/north_tau
                endif
            endif
            if (rate > 0._wp) then
                damp=exp(-dt*rate)
                h(1:ip,j)=h_ref(:,j)+damp*(h(1:ip,j)-h_ref(:,j))
                u(1:ip,j)=u_ref(:,j)+damp*(u(1:ip,j)-u_ref(:,j))
                v(1:ip,j)=v_ref(:,j)+damp*(v(1:ip,j)-v_ref(:,j))
            endif
        enddo
    end subroutine apply_latitude_sponge

	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the vorticity field, centred difference 
	!>@param[in] ip: number of east-west points
	!>@param[in] jp: ditto for north-south
	!>@param[in] o_halo: halos required for advection scheme
	!>@param[in] dt:  timestep
	!>@param[in] u,v: wind fields
	!>@param[inout] vort: vorticity field
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
	!>@param[in] cq: for efficiency
	!>@param[in] dp1: for efficiency
	!>@param[in] dq: for efficiency
	!>solves the 1-d advection equation:
	!>\f$ \zeta _r = \frac{1}{\cos\theta}
	!> \left( \frac{\partial u\cos\theta}{\partial \theta} -
	!> \frac{\partial v}{\partial \phi}\right)\f$
    subroutine diagnostics(ip,jp,o_halo,dt,u,v,vort,re,&
    		theta,thetan,dtheta,dthetan, phi, phin, dphi, dphin, &
    		recq, cq_s, cq, dp1, dq)

		use numerics_type
		implicit none
		integer(i4b), intent(in) :: ip,jp,o_halo
		real(wp), intent(in) :: dt, re
		real(wp), dimension(1-o_halo:jp+o_halo), intent(in) :: &
											theta,thetan, dtheta, dthetan
		real(wp), dimension(1-o_halo:ip+o_halo), intent(in) :: &											
											phi, phin, dphi, dphin
		real(wp), intent(in), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: &
																				u,v, &
										recq, cq_s, cq, dp1, dq
		real(wp), intent(inout), dimension(1:ip,1:jp) :: vort
		! local variables:
		integer(i4b) :: j, i
			
		
				
		! calculate relative vorticity 
		! (central difference ):
		vort(1:ip,1:jp)  =-1._wp/(recq(1:ip,1:jp))* &
			( (u(1:ip,2:jp+1)*cq(1:ip,2:jp+1)- &
			   u(1:ip,0:jp-1)*cq(1:ip,0:jp-1))/dq(1:ip,1:jp) - &
			  (v(2:ip+1,1:jp)-v(0:ip-1,1:jp))/dp1(1:ip,1:jp) ) 
			  

	end subroutine diagnostics
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	
	


	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>outputs variables to NetCDF file using MPI
	!>@param[inout] new_file: flag if this is a new file
	!>@param[in] outputfile: outputfilename
	!>@param[in] n: time-level
	!>@param[in] ip: number of east-west global grid
	!>@param[in] ipp: number of east-west levels on this PE
	!>@param[in] ipstart: start of i index on global grid
	!>@param[in] jp: ditto for south-north
	!>@param[in] jpp: ditto for south-north
	!>@param[in] jpstart: start of j index on global grid
	!>@param[in] o_halo: halo
	!>@param[in] time: time (s)
	!>@param[in] phi: longitude
	!>@param[in] theta: latitude
	!>@param[in] u_nudge: winds to nudge to
	!>@param[in] f_cor: Coriolis parameter
	!>@param[in] height: height of fluid
	!>@param[in] u: u-wind
	!>@param[in] v: v-wind
	!>@param[in] vort: vorticity
	!>@param[in] id: id
	!>@param[in] world_process: world_process
	!>@param[in] rank: rank
	!>@param[in] ring_comm: ring_comm
	subroutine output(new_file,outputfile,n,ip,ipp,ipstart,jp,jpp,jpstart,o_halo, &
					time,phi,theta, &
					u_nudge, f_cor, height, h, u, v, vort, &
				    id, world_process, rank, ring_comm)
	
		use netcdf
		use mpi
		use variables, only : MPI_INTEGER9

		implicit none
		logical, intent(inout) :: new_file
		character (len=*), intent(in) :: outputfile
		integer(i4b), intent(in) :: n, ip, ipp, ipstart, jp, jpp, jpstart, o_halo
		real(wp), intent(in) :: time
		real(wp), dimension(1-o_halo:ipp+o_halo), intent(in) :: phi
		real(wp), dimension(1-o_halo:jpp+o_halo), intent(in) :: theta, u_nudge
		
		real(wp), dimension(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), &
					intent(in) :: f_cor, height, h, u, v
		real(wp), dimension(1:ipp,1:jpp), intent(in) :: vort
		
		integer(i4b), intent(in) :: id ,world_process, rank, ring_comm
	
		integer(i4b) :: ncid, x_dimid, nx_dimid, ny_dimid, error, varid,a_dimid, id_go
		integer(i4b) :: i, tag1
		logical :: var


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! perform a blocking recv to wait for message from main process, 				 !
		! before carrying on															 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(id .ne. world_process) then
			tag1=id
			call MPI_Recv(var,1, MPI_LOGICAL, world_process, &
				tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE,error)
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	
		if((id==world_process) .and. new_file) then
			! open the file
		
			call check( nf90_create(outputfile, NF90_CLOBBER, ncid) )

			! define dimensions (netcdf hands back a handle)
			call check( nf90_def_dim(ncid, "times", NF90_UNLIMITED, x_dimid) )
			call check( nf90_def_dim(ncid, "ip", ip, nx_dimid) )
			call check( nf90_def_dim(ncid, "jp", jp, ny_dimid) )


			! close the file, freeing up any internal netCDF resources
			! associated with the file, and flush any buffers
			call check( nf90_close(ncid) )
		
			! now define some variables, units, etc
			call check( nf90_open(outputfile, NF90_WRITE, ncid) )
			
			
			! define mode
			call check( nf90_redef(ncid) )

			! define variable: time
			call check( nf90_def_var(ncid, "time", NF90_REAL, &
						(/x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "time", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "seconds") )
						
						
			! define variable: phi
			call check( nf90_def_var(ncid, "phi", NF90_REAL, &
						(/nx_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "phi", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "radians") )

			! define variable: theta
			call check( nf90_def_var(ncid, "theta", NF90_REAL, &
						(/ny_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "theta", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "radians") )

			! define variable: u_nudge
			call check( nf90_def_var(ncid, "u_nudge", NF90_REAL, &
						(/ny_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "u_nudge", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "m/s") )

			! define variable: f_cor
			call check( nf90_def_var(ncid, "f_cor", NF90_REAL, &
						(/nx_dimid, ny_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "f_cor", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "s**-1") )

			! define variable: height
			call check( nf90_def_var(ncid, "height", NF90_REAL, &
						(/nx_dimid, ny_dimid, x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "height", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "m") )

			! define variable: h
			call check( nf90_def_var(ncid, "h", NF90_REAL, &
						(/nx_dimid, ny_dimid, x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "h", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "m") )

			! define variable: u
			call check( nf90_def_var(ncid, "u", NF90_REAL, &
						(/nx_dimid, ny_dimid, x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "u", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "m/s") )

			! define variable: v
			call check( nf90_def_var(ncid, "v", NF90_REAL, &
						(/nx_dimid, ny_dimid, x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "v", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "m/s") )

			! define variable: vort
			call check( nf90_def_var(ncid, "vort", NF90_REAL, &
						(/nx_dimid, ny_dimid, x_dimid/), varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(ncid, "vort", a_dimid) )
			! units
			call check( nf90_put_att(ncid, a_dimid, &
					   "units", "s**-1") )



			! exit define mode
			call check( nf90_enddef(ncid) )
			
			
			call check( nf90_close(ncid) )

			new_file=.false.
		endif
	
! 	



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! now send messages from the main process to all other processes                 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(id == world_process) then
			do i=1,rank-1
				tag1=i
				call MPI_Send(var, 1, MPI_LOGICAL, i, &
						tag1, MPI_COMM_WORLD, error)
			enddo
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	
	

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! perform a blocking recv to wait for message from main process,                 !
		! before carrying on                               								 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if(id .ne. world_process) then
			tag1=id
			call MPI_Recv(id_go,1, MPI_INTEGER9, id-1, &
				tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE,error)
		else
			id_go=world_process ! lets us go for first run
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		
		
		
		

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! ****WRITE****																	 !			
		! now we can write to file - each PE writes its own segment						 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call check( nf90_open(outputfile, NF90_WRITE, ncid) )
		
		if(n == 1) then
			! write variable: phi
			call check( nf90_inq_varid(ncid, "phi", varid ) )
			call check( nf90_put_var(ncid, varid, phi(1:ipp), &
						start = (/1+ipstart/)))	

			! write variable: theta
			call check( nf90_inq_varid(ncid, "theta", varid ) )
			call check( nf90_put_var(ncid, varid, theta(1:jpp), &
						start = (/1+jpstart/)))	
			! write variable: u_nudge
			call check( nf90_inq_varid(ncid, "u_nudge", varid ) )
			call check( nf90_put_var(ncid, varid, u_nudge(1:jpp), &
						start = (/1+jpstart/)))	
			! write variable: f_cor
			call check( nf90_inq_varid(ncid, "f_cor", varid ) )
			call check( nf90_put_var(ncid, varid, f_cor(1:ipp,1:jpp), &
						start = (/1+ipstart,1+jpstart/)))	
			! write variable: height
			call check( nf90_inq_varid(ncid, "height", varid ) )
			call check( nf90_put_var(ncid, varid, height(1:ipp,1:jpp), &
						start = (/1+ipstart,1+jpstart,1/)))	
		endif

		if(id==world_process) then
			! write variable: time
			call check( nf90_inq_varid(ncid, "time", varid ) )
			call check( nf90_put_var(ncid, varid, time, &
						start = (/n/)))	
	    endif
	    
	    
		! write variable: h
		call check( nf90_inq_varid(ncid, "h", varid ) )
		call check( nf90_put_var(ncid, varid, h(1:ipp,1:jpp), &
					start = (/1+ipstart,1+jpstart,n/)))	

		! write variable: u
		call check( nf90_inq_varid(ncid, "u", varid ) )
		call check( nf90_put_var(ncid, varid, u(1:ipp,1:jpp), &
					start = (/1+ipstart,1+jpstart,n/)))	
		! write variable: v
		call check( nf90_inq_varid(ncid, "v", varid ) )
		call check( nf90_put_var(ncid, varid, v(1:ipp,1:jpp), &
					start = (/1+ipstart,1+jpstart,n/)))	
		! write variable: vort
		call check( nf90_inq_varid(ncid, "vort", varid ) )
		call check( nf90_put_var(ncid, varid, vort(1:ipp,1:jpp), &
					start = (/1+ipstart,1+jpstart,n/)))	
		call check( nf90_close(ncid) )
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! perform a send, to essentially allow next PE to resume and start the write     !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if((id == id_go).and.((id+1).lt.rank)) then
			tag1=id+1
			call MPI_Send(id+1, 1, MPI_INTEGER9, id+1, &
						tag1, MPI_COMM_WORLD, error)
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		if (rank > 1 ) then
			! send to world_process to complete ring
			tag1=2010
			if( ((id+1).eq.rank) ) then
				call MPI_Send(var, 1, MPI_LOGICAL, world_process, &
							tag1, MPI_COMM_WORLD, error)			
			endif


			! receive at world_process to complete ring
			if((id==world_process) ) then
				call MPI_Recv(var,1, MPI_LOGICAL,rank-1, &
					tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE,error)
			endif
		endif	
	


	end subroutine output
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! HELPER ROUTINE                                                       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine check(status)
		use netcdf
		use numerics_type
		integer(i4b), intent ( in) :: status

		if(status /= nf90_noerr) then
			print *, trim(nf90_strerror(status))
			stop "Stopped"
		end if
	end subroutine check
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	end module drivers
