	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>initialisation for the shallow water model
    module initialisation
    use numerics_type

    private
    public :: allocate_and_set
    contains
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Allocate and set arrays                                                            !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>allocate arrays on each PE, and initialise them
	!>@param[inout] ip for this grid
	!>@param[inout] jp for this grid
	!>@param[inout] ntim - number of time-steps
	!>@param[inout] f - rotation rate of planet
	!>@param[inout] re - radius of planet
	!>@param[inout] g - gravity of planet
	!>@param[inout] rho - density
	!>@param[inout] dphi - step in longitude
	!>@param[inout] dtheta - step in latitude
	!>@param[inout] dphin - step in longitude - staggered
	!>@param[inout] dthetan - step in latitude
	!>@param[inout] f_cor - coriolis parameter
	!>@param[inout] h - depth of fluid
	!>@param[inout] hs - height of surface above a reference
	!>@param[inout] u - east-west wind
	!>@param[inout] v - north-south wind
	!>@param[inout] height - total height of fluid
	!>@param[inout] dt - time-step
	!>@param[inout] dx - east-west step
	!>@param[inout] dy - north-south step
	!>@param[inout] x - east-west distance
	!>@param[inout] y - north-south distance
	!>@param[inout] phi - longitude
	!>@param[inout] theta - latitude
	!>@param[inout] phin - staggered longitude
	!>@param[inout] thetan - staggerd latitude
	!>@param[inout] recqdp - for efficiency
	!>@param[inout] recqdp_s - for efficiency
	!>@param[inout] recqdq_s - for efficiency
	!>@param[inout] redq_s - for efficiency
	!>@param[inout] redq - for efficiency
	!>@param[inout] recq - for efficiency
	!>@param[inout] cq_s - for efficiency
	!>@param[inout] cq - for efficiency
	!>@param[inout] dp1 - for efficiency
	!>@param[inout] dq - for efficiency
	!>@param[inout] recqdq - for efficiency
	!>@param[inout] u_nudge - wind to nudge to
	!>@param[in] o_halo - number of halos required
	!>@param[inout] ipstart - start of i indexing for this PE
	!>@param[inout] jpstart - start of j indexing for this PE
	!>@param[inout] coords - coordinates of cartesian topology
	!>@param[in] inputfile - netcdf file of saturn winds
	!>@param[in] add_random_height_noise - add noise to get going
	!>@param[in] height_noise_scheme - 0 legacy grid-cell noise; 1 correlated physical-scale noise
	!>@param[in] height_noise_amplitude - RMS height perturbation in metres for scheme 1
	!>@param[in] height_noise_corr_length - Gaussian correlation sigma in metres for scheme 1
	!>@param[in] initially_geostrophic - diagnose balanced winds from height after perturbation
	!>@param[in] momentum_metric_terms - 0 legacy/geostrophic; 1 spherical curvature/gradient-wind
	!>@param[in] initial_winds - flag: saturn, or jet?
	!>@param[in] ideal jet parameters: u_jet, theta_jet, h_jet
	!>@param[in] ip - global ip (all pes)
	!>@param[in] jp - global jp (all pes)
	!>@param[in] wind_factor - factor to multiply mean wind by
	!>@param[in] wind_shift - amount to shift mean wind by
	!>@param[in] wind_reduce - amount to reduce mean wind by
	!>@param[in] runtime - length of run
	!>@param[in] dt_nm - time-step
	!>@param[in] grav - gravity read in from namelist
	!>@param[in] rho_nm - density of fluid read in from namelist
	!>@param[in] re_nm - radius of planet read in from namelist
	!>@param[in] rotation_period_hours - rotation period of planet
	!>@param[in] scale_height - scale height of atmosphere
	!>@param[inout] slat - most southerly latitude
	!>@param[inout] nlat - most northerly latitude
	!>@param[in] slat_thresh - most southerly latitude
	!>@param[in] nlat_thresh - most northerly latitude
	!>@param[in] dims - dimensions of cartesian topology
	!>@param[in] id - id of this PE
	!>@param[in] comm2d - communicator for cartesian topology
	subroutine allocate_and_set(ipp,jpp,ntim, f, &
				re, g, rho, dphi, dtheta, dphin,dthetan, &
				f_cor,h,hs, u, v, &
				height, dt, dx, dy, x, y, &
				phi, theta, phin, thetan, &
				recqdp, recqdp_s, recqdq_s, redq_s,redq, &
				recq, cq_s, cq, dp1, dq, &
				recqdq, &
				u_nudge, o_halo, ipstart, jpstart, coords, &
				inputfile, add_random_height_noise, &
                height_noise_scheme, height_noise_amplitude, height_noise_corr_length, &
				initially_geostrophic, momentum_metric_terms, initial_winds, &
				u_jet, theta_jet, h_jet, &
				ip, jp, &
				wind_factor, wind_shift, wind_reduce, runtime, &
				dt_nm, grav, rho_nm, re_nm, &
				rotation_period_hours, scale_height, slat, nlat, &
				slat_thresh, nlat_thresh, dims, id, comm2d)
				
		use numerics_type
		use mpi
		use netcdf
		use numerics, only : find_pos, poly_int
		use random, only : random_normal
		use mpi_module
		
		implicit none
		! grid variables to be set
		integer(i4b), intent(inout) :: ipp, jpp, ntim, o_halo
    	real(wp), intent(inout) :: f, re, g, rho, dt
    	real(wp), intent(inout), allocatable, dimension(:,:) :: &
    									f_cor, h, hs, u, v, height, &
    											dx, dy, x, y, &
    				recqdp, recqdp_s, recqdq_s, redq_s, redq, &
    				recq, cq_s, cq, dp1, dq, recqdq
    	real(wp), intent(inout), allocatable, dimension(:) :: &
    									phi, theta, phin, thetan,u_nudge, &
    									dphi, dtheta, dphin, dthetan

		! namelist variables used to set the grid
		character (len=*), intent(in) :: inputfile
		logical, intent(in) :: add_random_height_noise, initially_geostrophic
		integer(i4b), intent(in) :: initial_winds, ip, jp, momentum_metric_terms, &
                                  height_noise_scheme
		real(wp), intent(in) :: wind_factor, wind_shift, wind_reduce, runtime, &
							dt_nm, grav, rho_nm, re_nm, &
                            height_noise_amplitude, height_noise_corr_length, &
							rotation_period_hours, scale_height, &
							slat_thresh, nlat_thresh, &
							u_jet, theta_jet, h_jet
		real(wp), intent(inout) :: slat, nlat
		integer(i4b), dimension(2), intent(in) :: dims
		integer(i4b), dimension(2), intent(inout) :: coords
		integer(i4b), intent(inout) :: ipstart, jpstart
		integer(i4b), intent(in) :: id, comm2d
		! locals
		integer(i4b) :: iloc, error, AllocateStatus, ncid, varid1,varid2, dimid, nlats, &
						i, j
		real(wp), dimension(:), allocatable :: latitude, wind
		real(wp), dimension(:,:), allocatable :: u_base
        real(wp), dimension(:,:), allocatable :: noise_raw, noise_tmp, noise_corr
		real(wp) :: var, dummy, delta_omega, slat_thresh2, nlat_thresh2, &
                    pgrad_y, pgrad_x, pgrad_y_base, kcurv, balance_freq, f_eff, &
                    lat_model, lat_sample, &
                    noise_mean, noise_rms, noise_sum, noise_sumsq, noise_weight, &
                    sigma_i, sigma_j, dx_noise, dy_noise, wgt, lat_global, &
                    band_south, band_north
		! for random number:
		real(wp) :: r
		real(wp), dimension(10,10) :: rs
		integer(i4b) :: k, nbottom, ntop, tag1, ii, jj, iii, jjj, &
                                  radius_i, radius_j, n_noise
		integer(i4b), allocatable, dimension(:) :: seed
		
		
		

		if(id>=dims(1)*dims(2)) return

		! start of initialisation code
		! 1. basic equalities:
		re=re_nm
		g=grav
		rho=rho_nm
		dt=dt_nm

		! 2. scalar formulae:
		ntim=ceiling(runtime/dt)
		f=2._wp*PI / (rotation_period_hours*3600._wp)
		
		! 3. arrays:




		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! find the number of grid points on each PE                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call MPI_CART_COORDS(comm2d, id, 2, coords, error)
! 		print *,'Coords of ',id,' are ',coords

		! number of grid points in all but last:
		ipp = floor(real(ip,wp)/real(dims(1),wp)) 
		ipstart = ipp*(coords(1))     
		if(coords(1) == (dims(1)-1)) then
			ipp=ip-(dims(1)-1)*ipp ! number of grid points in last
		endif
		! number of grid points in all but last:
		jpp = floor(real(jp,wp)/real(dims(2),wp))      
		jpstart = jpp*(coords(2))     
		if(coords(2) == (dims(2)-1)) then
			jpp=jp-(dims(2)-1)*jpp ! number of grid points in last
		endif
! 		print *,ip,jp,ipp,jpp,ipstart, jpstart,coords
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! allocate arrays                                                                !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		allocate( f_cor(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( h(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( hs(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( u(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( u_base(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( v(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( height(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( x(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( y(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dx(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dy(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		
		allocate( recqdp(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( recqdp_s(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( recqdq_s(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( redq_s(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( redq(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( recq(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( cq_s(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( cq(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dp1(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dq(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( recqdq(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	

		allocate( phi(1-o_halo:ipp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dphi(1-o_halo:ipp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( phin(1-o_halo:ipp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dphin(1-o_halo:ipp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( theta(1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dtheta(1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( thetan(1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dthetan(1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( u_nudge(1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! read netcdf file                                                               !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
		! the file.
		call check( nf90_open(inputfile, NF90_NOWRITE, ncid) )

		call check( nf90_inq_dimid(ncid, "nlats", dimid) )
		call check( nf90_inquire_dimension(ncid, dimid, len = nlats) )
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! allocate arrays                                                                !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		allocate( latitude(1:nlats), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( wind(1:nlats), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		

		! Get the varid of the latitude variable, based on its name.
		call check( nf90_inq_varid(ncid, "latitude", varid1) )
		! Get the varid of the wind variable, based on its name.
		call check( nf90_inq_varid(ncid, "wind", varid2) )
		
		call check( nf90_get_var(ncid, varid1, latitude, start = [1] ) )
		call check( nf90_get_var(ncid, varid2, wind, start = [1] ) )
		
		! Close the file, freeing all resources.
		call check( nf90_close(ncid) )
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		
		
		



		! interpolate to grid, do finite diffs, pass halos, etc
		dphi=(2._wp*PI) / real(ip,wp) ! lon: ip distinct periodic cells
		dphin=(2._wp*PI) / real(ip,wp) ! lon: ip distinct periodic cells

		! set up longitude array:
		phi=dphi(1)*(/(i,i=ipstart+1-o_halo-1,ipstart+ipp+o_halo-1)/)
		phin=phi+dphi/2._wp

		nlat=min(nlat,nlat_thresh)
		slat=max(slat,slat_thresh)
			
		! latitude array:	
		dtheta=(nlat-slat) &
				/ real(jp-1,wp) * PI/180._wp  ! lat
				
		! deal with singularity at the poles:
! 		if((coords(2) == (dims(2)-1))  .and. (nlat .gt. 0._wp)) then
! 			dtheta(jpp) = 2._wp*(90._wp-nlat)*PI/180._wp
! 		endif
! 		if((coords(2) == 0)  .and. (slat .lt. 0._wp)) then
! 			dtheta(0) = 2._wp*(90._wp+slat)*PI/180._wp
! 		endif
				


		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! set up latitude array need mpi to add them up:                                 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call MPI_CART_SHIFT( comm2d, 1, 1, nbottom, ntop, error)
		if(nbottom .ne. -1) then
			tag1=001
			call MPI_Recv(theta(1-o_halo), 1, MPI_REAL8, nbottom, &
					tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE,error)
		endif

		if(coords(2) == 0) then
			theta(0) = slat*PI/180._wp-dtheta(0)
		endif

		do i=1,jpp+o_halo
			theta(i)=theta(i-1)+dtheta(i-1)
		enddo
		
		if(ntop .ne. -1) then
			tag1=001
			call MPI_Send(theta(jpp), 1, MPI_REAL8, ntop, &
				tag1, MPI_COMM_WORLD,error)
		endif		
! 		theta=dtheta(10)*(/(i,i=jpstart+1-o_halo-1,jpstart+jpp+o_halo-1)/) + slat*PI/180._wp
		
		

		
		thetan=theta+dtheta/2._wp
		do i=1-o_halo,jpp
			dthetan(i)=thetan(i+1)-thetan(i)
		enddo
		
		! if this is the top then set dtheta:

		if(ntop == -1 ) then
			dthetan(jpp+1)=2._wp*(90._wp*PI/180._wp -dthetan(jpp))
		endif 
		if(ntop /= -1 ) then
			call MPI_Recv(dthetan(jpp+o_halo), 1, MPI_REAL8, ntop, &
					tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE,error)
		endif 
		if(nbottom /= -1 ) then
			call MPI_Send(dthetan(1), 1, MPI_REAL8, nbottom, &
					tag1, MPI_COMM_WORLD,error)
		endif 	
! 		dthetan=(nlat-slat) / real(jp-1,wp) * PI/180._wp  ! lat
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		select case (initial_winds)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! saturn winds:                                                              !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			case (1) 
			! Saturn mean-wind profile.  The three wind controls act on the
			! observed/reference profile itself, before the rotation-frame correction:
			!
			!   U_new(phi) = wind_factor * U_ref(phi - wind_shift) - wind_reduce
			!
			! wind_shift is in degrees: positive shifts profile features northward,
			! negative shifts them southward.  wind_reduce is a literal subtraction
			! in m/s (so sufficiently weak eastward winds can become westward).
			delta_omega=2._wp*PI/(3600._wp)*(1._wp/10.656_wp-1._wp/rotation_period_hours)
			do i=1-o_halo,jpp+o_halo
				lat_model=asin(sin(theta(i)))*180._wp/PI
				lat_sample=lat_model-wind_shift

				iloc=find_pos(latitude(1:nlats),lat_sample)
				iloc=min(nlats-1,iloc)
				iloc=max(1,iloc)

				! Linear interpolation of the reference wind at the shifted latitude.
				! Preserve the legacy endpoint clamping outside the supplied profile.
				call poly_int(latitude(iloc:iloc+1), wind(iloc:iloc+1), &
					min(max(lat_sample,latitude(nlats)),latitude(1)), var,dummy)
				
				if((theta(i)*180._wp/pi>90._wp) .or. &
					(theta(i)*180._wp/pi<-90._wp) ) var=-var

				! Apply requested Saturn-wind transformations in profile space.
				var=wind_factor*var
				var=var-wind_reduce

				! Convert from the reference System III period used by the wind data
				! to the model's chosen rotation period.  Do not scale/reduce this term.
				var=var+delta_omega*re*cos(theta(i))
			
				u_nudge(i)=var
			enddo
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! ideal jet:                                                                 !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			case (2) 
			do i=1-o_halo,jpp+o_halo
				u_nudge(i)=u_jet* exp(-0.5_wp*(theta(i)-theta_jet*pi/180._wp)**2._wp &
				 				/ (h_jet*pi/180._wp)**2._wp )
			enddo
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			case default
				print *,'error initial_winds',initial_winds
				stop
		end select
		
		! calculate Coriolis param:
		do i=1-o_halo,ipp+o_halo
			f_cor(i,:)=2._wp*f*sin(theta)			
		enddo
! 		if(coords(2)==0) f_cor(:,1-o_halo)=-f_cor(:,1-o_halo)
! 		if(coords(2)==dims(2)-1) f_cor(:,jpp+o_halo)=-f_cor(:,jpp+o_halo)
			
		! calculate x, y, dx, dy:
		do i=1-o_halo,ipp+o_halo
			x(i,:)=phi(i)*cos(asin(sin(theta)))*re			
		enddo
		do i=1-o_halo,ipp+o_halo
			y(i,:)=(theta)*re			
		enddo
		
		do i=1-o_halo,ipp+o_halo
			dx(i,:)=dphi(i)*cos(asin(sin(theta)))*re			
		enddo
		do i=1-o_halo,ipp+o_halo
			dy(i,:)=(dtheta)*re			
		enddo
		
		
		! set u and v:
		v(:,:)=0._wp
		do j=1-o_halo,jpp+o_halo
			u(:,j)=u_nudge(j)
		enddo
		u_base(:,:)=u(:,:)
		
		! set surface to zero.
		hs(:,:)=0._wp




		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! calculate the height field from winds	- MPI needed to span sub-domains		 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
		call MPI_CART_COORDS(comm2d, id, 2, coords, error)
		call MPI_CART_SHIFT( comm2d, 1, 1, nbottom, ntop, error)	
		
		! if the y coordinate is not the most northerly
		! then receive from domain north of here:
		if ( (coords(2)+1) /= dims(2) ) then
			tag1=2010
			call MPI_Recv(height(:,jpp+1:jpp+o_halo),& ! the data packet to receive into
				ipp+2*o_halo, & ! size of the data packet
				MPI_REAL8, &
				ntop, & ! receive from above (north)
				tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE,error)
		
		endif		
		! if most northerly grid point (note set the jpp+1 point)
		if( (coords(2)+1) == dims(2) ) then
			height(:,jpp+1)=scale_height
		endif

		
        select case (momentum_metric_terms)
        case (0)
            ! Legacy geostrophic balance used by the original model.
            do j=jpp,0,-1
                height(:,j)=height(:,j+1)+ &
                        0.25_wp*(f_cor(:,j+1)+f_cor(:,j))* &
                        (u(:,j+1)+u(:,j))* &
                        re/g*(dtheta(j))
            enddo
        case (1)
            ! Gradient-wind balance for a zonal flow on a sphere:
            !   g/re dh/dtheta = -f*u - u^2*tan(theta)/re.
            ! Integrate southward from the prescribed northern height,
            ! using midpoint values between adjacent latitude rows.
            do j=jpp,0,-1
                height(:,j)=height(:,j+1)+ &
                        ( 0.25_wp*(f_cor(:,j+1)+f_cor(:,j))* &
                          (u(:,j+1)+u(:,j))*re + &
                          0.25_wp*(u(:,j+1)+u(:,j))**2*tan(thetan(j)) ) * &
                        dtheta(j)/g
            enddo
        case default
            write(*,*) 'ERROR: unknown momentum_metric_terms = ', momentum_metric_terms
            write(*,*) '       valid values are 0 (legacy/off) and 1 (spherical/on)'
            stop 1
        end select
		! if the y coordinate is not most southerly
		if ( (coords(2)+1) /= 1 ) then
			tag1=2010
			call MPI_Send(height(:,1:1),& ! the data packet to send
				ipp+2*o_halo, & ! size of the data packet
				MPI_REAL8, &
				nbottom, & ! send to below (south)
				tag1, MPI_COMM_WORLD,error)
		
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		
		
		
		
		
		
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! calculate and add noise														 !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
		if (add_random_height_noise) then
            call random_seed(size=k)
            allocate(seed(1:k))
            seed(:)=2
            call random_seed(put=seed)

            select case (height_noise_scheme)
            case (0)
                ! Legacy behaviour: independent grid-cell Gaussian height noise,
                ! with the original nominal ~1000 m amplitude expression.
                select case (initial_winds)
                case (1)
                    do j=1,jp
                        do i=1,ip
                            r=random_normal()
                            if((i > ipstart) .and. (i <=ipstart+ipp) &
                                .and. (j > jpstart) .and. (j <= jpstart+jpp) ) then
                                if ((theta(j-jpstart)*180._wp/PI) > 75._wp &
                                    .and. (theta(j-jpstart)*180._wp/PI) < 80._wp) then
                                    height(i-ipstart,j-jpstart) = &
                                        height(i-ipstart,j-jpstart) + &
                                        r*1000.e0_wp*0.6e5_wp/height(i-ipstart,j-jpstart)
                                endif
                            endif
                        enddo
                    enddo
                case (2)
                    do j=1,jp
                        do i=1,ip
                            r=random_normal()
                            if((i > ipstart) .and. (i <=ipstart+ipp) &
                                .and. (j > jpstart) .and. (j <= jpstart+jpp) ) then
                                if ((theta(j-jpstart)*180._wp/PI) > (theta_jet-h_jet*3._wp) &
                                    .and. (theta(j-jpstart)*180._wp/PI) < (theta_jet+h_jet*3._wp)) then
                                    height(i-ipstart,j-jpstart) = &
                                        height(i-ipstart,j-jpstart) + &
                                        r*1000.e0_wp*0.6e5_wp/height(i-ipstart,j-jpstart)
                                endif
                            endif
                        enddo
                    enddo
                case default
                    print *,'error initial_winds',initial_winds
                    stop
                end select

            case (1)
                ! Resolution-independent correlated perturbation.  Generate the
                ! same global Gaussian field on every MPI rank, smooth it with a
                ! Gaussian kernel whose standard deviation is specified in metres,
                ! then normalise the perturbation over the active jet band so that
                ! height_noise_amplitude is its RMS height perturbation in metres.
                if (height_noise_amplitude < 0._wp) then
                    write(*,*) 'ERROR: height_noise_amplitude must be >= 0'
                    stop 1
                endif
                if (height_noise_corr_length < 0._wp) then
                    write(*,*) 'ERROR: height_noise_corr_length must be >= 0'
                    stop 1
                endif

                select case (initial_winds)
                case (1)
                    band_south=75._wp
                    band_north=80._wp
                case (2)
                    band_south=theta_jet-3._wp*h_jet
                    band_north=theta_jet+3._wp*h_jet
                case default
                    print *,'error initial_winds',initial_winds
                    stop
                end select

                allocate(noise_raw(1:ip,1:jp), noise_tmp(1:ip,1:jp), &
                         noise_corr(1:ip,1:jp), STAT=AllocateStatus)
                if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

                do j=1,jp
                    do i=1,ip
                        noise_raw(i,j)=random_normal()
                    enddo
                enddo

                if (height_noise_corr_length > 0._wp) then
                    ! Meridional Gaussian smoothing.
                    dy_noise=re*abs((nlat-slat)*PI/180._wp/real(jp-1,wp))
                    sigma_j=height_noise_corr_length/max(dy_noise,tiny(1._wp))
                    radius_j=min(jp-1,max(1,ceiling(3._wp*sigma_j)))
                    do j=1,jp
                        do i=1,ip
                            noise_sum=0._wp
                            noise_weight=0._wp
                            do jj=-radius_j,radius_j
                                jjj=j+jj
                                if (jjj < 1 .or. jjj > jp) cycle
                                wgt=exp(-0.5_wp*(real(jj,wp)/sigma_j)**2)
                                noise_sum=noise_sum+wgt*noise_raw(i,jjj)
                                noise_weight=noise_weight+wgt
                            enddo
                            noise_tmp(i,j)=noise_sum/noise_weight
                        enddo
                    enddo

                    ! Zonal Gaussian smoothing.  The number of grid points in the
                    ! kernel varies with latitude so the physical correlation scale
                    ! remains approximately constant on the lat-lon grid.
                    do j=1,jp
                        lat_global=(slat+(nlat-slat)*real(j-1,wp)/real(jp-1,wp))*PI/180._wp
                        dx_noise=re*abs(cos(lat_global))*(2._wp*PI/real(ip,wp))
                        sigma_i=height_noise_corr_length/max(dx_noise,tiny(1._wp))
                        radius_i=min(ip/2,max(1,ceiling(3._wp*sigma_i)))
                        do i=1,ip
                            noise_sum=0._wp
                            noise_weight=0._wp
                            do ii=-radius_i,radius_i
                                iii=modulo(i-1+ii,ip)+1
                                wgt=exp(-0.5_wp*(real(ii,wp)/sigma_i)**2)
                                noise_sum=noise_sum+wgt*noise_tmp(iii,j)
                                noise_weight=noise_weight+wgt
                            enddo
                            noise_corr(i,j)=noise_sum/noise_weight
                        enddo
                    enddo
                else
                    noise_corr=noise_raw
                endif

                ! Remove the band mean and scale to the requested RMS amplitude.
                noise_sum=0._wp
                noise_sumsq=0._wp
                n_noise=0
                do j=1,jp
                    lat_global=slat+(nlat-slat)*real(j-1,wp)/real(jp-1,wp)
                    if (lat_global > band_south .and. lat_global < band_north) then
                        do i=1,ip
                            noise_sum=noise_sum+noise_corr(i,j)
                            noise_sumsq=noise_sumsq+noise_corr(i,j)**2
                            n_noise=n_noise+1
                        enddo
                    endif
                enddo
                if (n_noise <= 0) then
                    write(*,*) 'ERROR: no grid points in height-noise latitude band'
                    stop 1
                endif
                noise_mean=noise_sum/real(n_noise,wp)
                noise_rms=sqrt(max(0._wp,noise_sumsq/real(n_noise,wp)-noise_mean**2))
                if (noise_rms <= tiny(1._wp) .and. height_noise_amplitude > 0._wp) then
                    write(*,*) 'ERROR: correlated height-noise RMS is zero'
                    stop 1
                endif

                do j=1,jp
                    lat_global=slat+(nlat-slat)*real(j-1,wp)/real(jp-1,wp)
                    if (lat_global > band_south .and. lat_global < band_north) then
                        if (j > jpstart .and. j <= jpstart+jpp) then
                            do i=1,ip
                                if (i > ipstart .and. i <= ipstart+ipp) then
                                    if (height_noise_amplitude > 0._wp) then
                                        height(i-ipstart,j-jpstart)=height(i-ipstart,j-jpstart) + &
                                            height_noise_amplitude* &
                                            (noise_corr(i,j)-noise_mean)/noise_rms
                                    endif
                                endif
                            enddo
                        endif
                    endif
                enddo

                deallocate(noise_raw,noise_tmp,noise_corr)

            case default
                write(*,*) 'ERROR: unknown height_noise_scheme = ',height_noise_scheme
                stop 1
            end select

            deallocate(seed)
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		



		! set halos in Coriolis array
		call exchange_halos(comm2d, id, ipp, jpp, o_halo, f_cor)
		! set halos in height array
		call exchange_halos(comm2d, id, ipp, jpp, o_halo, height)


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
		! calculate u and v from height field (with noise) - no message passing needed
		! if halos in height are set correctly
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
		if(initially_geostrophic) then
            select case (momentum_metric_terms)
            case (0)
                ! Legacy geostrophic diagnosis after the random height
                ! perturbation has been added.
                do j=1,jpp
                    u(:,j)=-g*(height(:,j+1)-height(:,j-1)) / &
                            (re*(dtheta(j)+dtheta(j-1))*f_cor(:,j))
                enddo
                do i=1,ipp
                    v(i,:)=g*(height(i+1,:)-height(i-1,:)) / &
                            (re*(dphi(i)+dphi(i-1))*cos(theta(:))*f_cor(i,:))
                enddo

            case (1)
                ! The zonally symmetric base jet and height field are already in
                ! exact gradient-wind balance.  After adding an arbitrary 2-D
                ! height perturbation, an exact local nonlinear gradient-wind
                ! inversion is not guaranteed to have a real solution.  Diagnose
                ! the perturbation winds from the balance linearised about the
                ! known base jet U:
                !
                !   (f + 2*k*U) u' + P'_y = 0,   k=tan(theta)/re
                !   (f +   k*U) v'       = P'_x
                !
                ! where P'_y is the total meridional pressure-gradient
                ! acceleration minus that of the gradient-wind-balanced base state.
                do j=1,jpp
                    kcurv=tan(theta(j))/re
                    do i=1-o_halo,ipp+o_halo
                        pgrad_y=g*(height(i,j+1)-height(i,j-1)) / &
                                (re*(dtheta(j)+dtheta(j-1)))
                        pgrad_y_base=-f_cor(i,j)*u_base(i,j) - &
                                     kcurv*u_base(i,j)**2
                        balance_freq=f_cor(i,j)+2._wp*kcurv*u_base(i,j)
                        if (abs(balance_freq) <= tiny(1._wp)) then
                            write(*,*) 'ERROR: zero linear gradient-wind frequency'
                            write(*,*) ' i,j,theta(deg),frequency = ', &
                                i,j,theta(j)*180._wp/pi,balance_freq
                            stop 1
                        endif
                        u(i,j)=u_base(i,j) - &
                               (pgrad_y-pgrad_y_base)/balance_freq
                    enddo
                enddo

                ! The base height has no zonal gradient, so P_x is entirely the
                ! perturbation pressure gradient.  Linearised zonal momentum
                ! balance gives (f + k*U) v' = P'_x.
                do j=1,jpp
                    kcurv=tan(theta(j))/re
                    do i=1,ipp
                        pgrad_x=g*(height(i+1,j)-height(i-1,j)) / &
                                (re*(dphi(i)+dphi(i-1))*cos(theta(j)))
                        f_eff=f_cor(i,j)+kcurv*u_base(i,j)
                        if (abs(f_eff) <= tiny(1._wp)) then
                            write(*,*) 'ERROR: zero effective Coriolis in linear gradient-wind initialisation'
                            write(*,*) ' i,j,theta(deg),f_eff = ', &
                                i,j,theta(j)*180._wp/pi,f_eff
                            stop 1
                        endif
                        v(i,j)=pgrad_x/f_eff
                    enddo
                enddo
            end select
		endif	
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		


		! set halos in u array
		call exchange_halos(comm2d, id, ipp, jpp, o_halo, u)
		! set halos in v array
		call exchange_halos(comm2d, id, ipp, jpp, o_halo, v)


		! set h array
		h(:,:)=height(:,:)-hs(:,:)
		
		
		! set some variables for efficiency
		do j=1-o_halo,jpp+o_halo
			recqdp(:,j)=re*cos(theta(j))*dphi(:)
			recqdp_s(:,j)=re*cos(theta(j))*dphin(:)
			recqdq_s(:,j)=re*cos(thetan(j))*dthetan(j)
			redq_s(:,j)=re*dthetan(j)
			redq(:,j)=re*dtheta(j)

			recq(:,j)=re*cos(theta(j))
			recqdq(:,j)=re*cos(theta(j))*dtheta(j)
			cq_s(:,j)=cos(thetan(j))
			cq(:,j)=cos(theta(j))
			dp1(:,j)=dphi(:)
			dq(:,j)=dtheta(j)
		enddo
		
		

		deallocate(wind)
		deallocate(latitude)
		
	end subroutine allocate_and_set
	

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

	end module initialisation
	
