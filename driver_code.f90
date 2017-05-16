	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>drivers for the shallow water model
    module drivers
    use nrtype
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
				vis_eq, lat_eq, &
				dims,id, world_process, rank, ring_comm)
		use nrtype
		use mpi_module
		use advection

		implicit none
		logical, intent(inout) :: new_file
		logical, intent(in) :: nudge, viscous_dissipation, dissipate_h
		integer(i4b), intent(in) :: ip,ipp, jp,jpp, ntim, o_halo, ipstart, jpstart, &
									subgrid_model
		integer(i4b), intent(in) :: id, world_process, ring_comm, rank
		integer(i4b), dimension(2), intent(in) :: coords, dims
		character (len=*), intent(in) :: outputfile
		real(sp), intent(in) :: f, re, g, rho, dt, output_interval
		real(sp), dimension(1-o_halo:ipp+o_halo), intent(in) :: phi, phin, dphi, dphin
		real(sp), dimension(1-o_halo:jpp+o_halo), intent(in) :: theta, thetan, u_nudge, &
																dtheta, dthetan
		real(sp), dimension(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), &
					intent(in) :: f_cor, hs, &
    				recqdp, recqdp_s, recqdq_s, redq_s, redq, &
    				recq, cq_s, cq, dp1, dq, recqdq
		real(sp), dimension(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), &
					intent(inout) :: h, u, v, height
		real(sp), dimension(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), &
					intent(in) :: dx, dy, x, y
		real(sp), intent(in) :: vis, nudge_tau, cvis, lat_eq, vis_eq
					
		! locals:		
		integer(i4b) :: n, cur=1, j, error
		real(sp) :: time, time_last_output, output_time
		real(sp), dimension(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo) :: &
				u_old, v_old, h_old
		real(sp), dimension(1:ipp,1:jpp) :: delsq, vort, visco
		

		time_last_output=-output_interval
		output_time=output_interval

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! time-loop                                                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		do n=1,ntim	
		
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! write netcdf variables                                                     !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			time=real(n-1,sp)*dt
			if (time-time_last_output >= output_interval) then
				if (id==world_process) &
					print *,'output no ',cur,' at time (hrs) ', &
						time/3600._sp,n,' steps of ',ntim
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
							id, world_process, rank, ring_comm)
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
    			recqdq, recqdp, recqdp_s, recqdq_s, redq_s, redq, cq, cq_s)	    		
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
						(0.5_sp*(u(1:ipp,j)+u_old(1:ipp,j)))/real(1,sp) ) &
						/nudge_tau * dt
! 					v(1:ipp,j)=v(1:ipp,j)+&
! 						(0._sp- &
! 						sum(0.5_sp*(v(1:ipp,j)+v_old(1:ipp,j)))/real(ipp,sp) ) &
! 						/nudge_tau *dt

! 					Derived by integrating du/dt=(u_nudge-u)/tau
! 					u(1:ipp,j)=u_nudge(j)- &
! 						(u_nudge(j)- &
! 						(u(1:ipp,j))/real(1,sp) ) * &
! 						exp(-dt/nudge_tau )
! 					v(1:ipp,j)=0._sp- &
! 						(0._sp- &
! 						(v(1:ipp,j))/real(1,sp) ) * &
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
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				
				
				
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! dissipate u                                                            !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				call dissipation(ipp,jpp,o_halo,dt,0.5_sp*(u_old+u), delsq,re,&
					theta,thetan,dtheta,dthetan, phi, phin, dphi, dphin, &
					recq, cq_s, dp1, dq)
										
				select case(subgrid_model)
				case (1)
					u(1:ipp,1:jpp)=u(1:ipp,1:jpp)+dt*delsq*vis			
				case (2)
					call smagorinsky(ipp,jpp,o_halo,cvis,0.5_sp*(u_old+u),&
									0.5_sp*(v_old+v),visco,re,recq, dp1, dq)
					u(1:ipp,1:jpp)=u(1:ipp,1:jpp)+dt*delsq*visco			
				case default
					print *,'error subgrid ',subgrid_model
				end select
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




				


				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! dissipate v                                                            !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				call dissipation(ipp,jpp,o_halo,dt,0.5_sp*(v_old+v), delsq,re,&
					theta,thetan,dtheta,dthetan, phi, phin, dphi, dphin, &
					recq, cq_s, dp1, dq)
								
				select case(subgrid_model)
				case (1)
					v(1:ipp,1:jpp)=v(1:ipp,1:jpp)+dt*delsq*vis
				case (2)
					v(1:ipp,1:jpp)=v(1:ipp,1:jpp)+dt*delsq*visco			
				case default
					print *,'error subgrid ',subgrid_model
				end select
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


				do j=1,jpp
					if((theta(j) >-lat_eq*pi/180._sp) .and. &
						 (theta(j) < lat_eq*pi/180._sp)) then

						v(1:ipp,j)=v(1:ipp,j)+dt*delsq(1:ipp,j)*vis_eq* &
							cos(theta(j)*90._sp/lat_eq)
					endif
				enddo

				
					
							
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! dissipate h                                                            !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				if (dissipate_h .and. (subgrid_model == 1)) then
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					! halo exchanges                                                     !
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					call exchange_halos(ring_comm, id, ipp, jpp, o_halo, h)
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

					call dissipation(ipp,jpp,o_halo,dt,0.5_sp*(h_old+h), delsq,re,&
						theta,thetan,dtheta,dthetan, phi, phin, dphi, dphin, &
						recq, cq_s, dp1, dq)
						
						
					h(1:ipp,1:jpp)=h(1:ipp,1:jpp)+dt*delsq*vis
				endif
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



			endif	    		
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! halo exchanges                                                             !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			call exchange_halos(ring_comm, id, ipp, jpp, o_halo, h)
			call exchange_halos(ring_comm, id, ipp, jpp, o_halo, u)
			call exchange_halos(ring_comm, id, ipp, jpp, o_halo, v)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
! 			if(coords(2)==(dims(2)-1)) then
! 				v(:,jpp+1:jpp+o_halo)=0._sp
! 				do j=jpp+1,jpp+o_halo
! 					u(:,j)=u_nudge(j)
! 				enddo
! 				h(:,jpp+o_halo)=h(:,jpp)
! 			endif
! 			if(coords(2)==0) then
! 				v(:,1-o_halo)=0._sp
! 				do j=1-o_halo,0
! 					u(:,j)=u_nudge(j)
! 				enddo
! 				h(:,1-o_halo)=h(:,1)
! 			endif
			

		enddo
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








		
	end subroutine model_driver
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	
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

		use nrtype
		implicit none
		integer(i4b), intent(in) :: ip,jp,o_halo
		real(sp), intent(in) :: dt, re
		real(sp), dimension(1-o_halo:jp+o_halo), intent(in) :: &
											theta,thetan, dtheta, dthetan
		real(sp), dimension(1-o_halo:ip+o_halo), intent(in) :: &											
											phi, phin, dphi, dphin
		real(sp), intent(in), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: &
																				u,v, &
										recq, cq_s, cq, dp1, dq
		real(sp), intent(inout), dimension(1:ip,1:jp) :: vort
		! local variables:
		integer(i4b) :: j, i
			
		
				
		! calculate relative vorticity 
		! (central difference ):
		vort(1:ip,1:jp)  =-1._sp/(recq(1:ip,1:jp))* &
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
		real(sp), intent(in) :: time
		real(sp), dimension(1-o_halo:ipp+o_halo), intent(in) :: phi
		real(sp), dimension(1-o_halo:jpp+o_halo), intent(in) :: theta, u_nudge
		
		real(sp), dimension(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), &
					intent(in) :: f_cor, height, h, u, v
		real(sp), dimension(1:ipp,1:jpp), intent(in) :: vort
		
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


	


	end subroutine output
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! HELPER ROUTINE                                                       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine check(status)
		use netcdf
		use nrtype
		integer(i4b), intent ( in) :: status

		if(status /= nf90_noerr) then
			print *, trim(nf90_strerror(status))
			stop "Stopped"
		end if
	end subroutine check
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	end module drivers
