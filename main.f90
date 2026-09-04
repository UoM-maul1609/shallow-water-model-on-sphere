	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2017
	!>@brief
	!>Shallow Water Model on Sphere (SWMOS): 
	!>Solves the shallow water equations on a rotating sphere for a lat / lon grid:
	!> <br> <b>Continuity:</b> <br>
	!>\f$	\frac{\partial h}{\partial t}+\frac{1}{R_e\cos\theta}
    !>           \frac{\partial }{\partial \theta} \left(hv_\theta\cos\theta\right)+
    !>           \frac{1}{R_e\cos\theta}
    !> \frac{\partial}{\partial \phi} \left(hv_\phi\right)=0 
    !>\f$
    !>
	!> <br> <b>horizontal momentum 1:</b> <br>
	!>\f$ \frac{\partial v_\phi h}{\partial t} + 
	!> \frac{1}{R_e\cos\theta}\frac{\partial }{\partial \theta}
	!> \left(v_\phi v_\theta h\cos\theta \right) + 
    !>     \frac{1}{R_e\cos\theta}\frac{\partial }{\partial \phi}
    !> \left(v_\phi ^2 h+\frac{gh^2}{2}\right) =
    !>      -\frac{gh}{R_e\cos\theta}
    !> \frac{\partial }{\partial \phi} \left(H\right)+fhv_\theta
    !>\f$
    !>
	!> <br><br>
	!> <br> <b>horizontal momentum 2:</b> <br>
	!>\f$ \frac{\partial v_\theta h}{\partial t} + 
	!> \frac{1}{R_e\cos\theta}\frac{\partial }{\partial \theta}
	!> \left(v_\theta ^2 h+\frac{gh^2\cos\theta}{2}\right) + 
    !>      \frac{1}{R_e\cos\theta}\frac{\partial }{\partial \phi}
    !> \left(v_\phi v_\theta h\right) =
    !>       -\frac{gh}{R_e}\frac{\partial }{\partial \theta} \left(H \right)-fhv_\phi
    !>\f$
    !>
	!> <br><br>
	!> compile using the Makefile (note requires netcdf) and then run using: <br>
	!> mpiexec -n 4 ./main.exe namelist.in
	!> <br><br>
	!> (namelist used for initialisation).
	!> <br><br>



	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>main programme reads in information, allocates arrays, then calls the model driver

    program main
        use variables
        use mpi
        use mpi_module
        use initialisation
        use drivers
        use numerics_type, only : wp
        
        implicit none
        character (len=200) :: nmlfile = ' '
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! namelist for run variables                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        namelist /run_vars/ nm1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! MPI initialisation                                                   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call MPI_Init ( mp1%error )
		call MPI_Comm_rank ( MPI_COMM_WORLD, mp1%id, mp1%error )
		call MPI_Comm_size ( MPI_COMM_WORLD, mp1%rank, mp1%error )
		mp1%wtime = MPI_Wtime ( )	
		print *,'MPI running with ID: ',mp1%id,' and rank ',mp1%rank
		call mpi_define(MPI_INTEGER9)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! read in namelists													   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call getarg(1,nmlfile)
        open(8,file=nmlfile,status='old', recl=80, delim='apostrophe')
        read(8,nml=run_vars)
        close(8)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		grid1%o_halo=1


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Block until processors have synced	     						   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call block_ring(MPI_COMM_WORLD,mp1%id,world_process,mp1%rank)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Set-up the Cartesian topology										   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Let MPI choose a factorisation that uses every process.  Keep
		! reorder=.true.; all Cartesian neighbour communication uses ring_comm.
		mp1%dims=[0_i4b,0_i4b]
		call MPI_Dims_create(mp1%rank,mp1%ndim,mp1%dims,mp1%error)
		mp1%dx=mp1%dims(1)
		mp1%dy=mp1%dims(2)

		if(mp1%id == world_process) then
			print *,'Cartesian topology: ',mp1%dx, mp1%dy
		endif

		mp1%periods=[.true.,.false.]
		! cart topo:
		call MPI_CART_CREATE( MPI_COMM_WORLD, mp1%ndim, mp1%dims, &
							mp1%periods, mp1%reorder,mp1%ring_comm, mp1%error )
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Allocate and initialise arrays									   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call allocate_and_set(grid1%ip,grid1%jp,grid1%ntim, grid1%f, &
				grid1%re, grid1%g, grid1%rho, grid1%dphi, grid1%dtheta, &
				grid1%dphin, grid1%dthetan, &
				grid1%f_cor,grid1%h,grid1%hs, grid1%u, grid1%v, &
				grid1%height, grid1%dt,grid1%dx, grid1%dy, grid1%x, grid1%y, &
				grid1%phi, grid1%theta, grid1%phin, grid1%thetan, &
   				grid1%recqdp, grid1%recqdp_s, grid1%recqdq_s, grid1%redq_s, grid1%redq, &
    				grid1%recq, grid1%cq_s, grid1%cq, grid1%dp1, grid1%dq, &
    				grid1%recqdq, &
 				grid1%u_nudge,grid1%o_halo, &
				grid1%ipstart, grid1%jpstart, grid1%coords, &
				nm1%inputfile, nm1%add_random_height_noise, &
                nm1%height_noise_scheme, nm1%height_noise_amplitude, &
                nm1%height_noise_corr_length, &
				nm1%initially_geostrophic, nm1%momentum_metric_terms, nm1%initial_winds, &
				nm1%u_jet, nm1%theta_jet, nm1%h_jet, &
				nm1%ip, nm1%jp, &
				nm1%wind_factor, nm1%wind_shift, nm1%wind_reduce, nm1%runtime, &
				nm1%dt, nm1%grav, nm1%rho, nm1%re, &
				nm1%rotation_period_hours, nm1%scale_height, nm1%slat, &
				nm1%nlat, nm1%slat_thresh, nm1%nlat_thresh, &
				mp1%dims, mp1%id, mp1%ring_comm)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Block until processors have synced	     						   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call block_ring(MPI_COMM_WORLD,mp1%id,world_process,mp1%rank)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Driver code: time-loop, advance solution, output	   				   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 		call model_driver(nm1%ip,grid1%ip,nm1%jp,grid1%jp,grid1%ntim, grid1%f, &
				grid1%re, grid1%g, grid1%rho, grid1%dphi, grid1%dtheta, &
				grid1%dphin, grid1%dthetan, &
				grid1%f_cor,grid1%h,grid1%hs, grid1%u, grid1%v, &
				grid1%height, grid1%dt, grid1%dx, grid1%dy, grid1%x, grid1%y, &
				grid1%phi, grid1%theta, grid1%phin, grid1%thetan, &
   				grid1%recqdp, grid1%recqdp_s, grid1%recqdq_s, grid1%redq_s, grid1%redq, &
    			grid1%recq, grid1%cq_s, grid1%cq, grid1%dp1, grid1%dq, &
    			grid1%recqdq, &
				grid1%u_nudge,grid1%o_halo, &
				grid1%ipstart, grid1%jpstart, grid1%coords, &
				io1%new_file, nm1%outputfile, nm1%output_interval, &
				nm1%nudge,nm1%nudge_timescale,nm1%nudge_scheme, &
				nm1%subgrid_model, nm1%viscous_dissipation, &
				nm1%dissipate_h,nm1%vis,nm1%cvis, &
				nm1%vis_eq,nm1%lat_eq, nm1%coriolis_scheme, nm1%momentum_metric_terms, &
                nm1%smagorinsky_scheme, nm1%lat_boundary_scheme, &
                nm1%sponge_south_width, nm1%sponge_north_width, &
                nm1%sponge_south_timescale, nm1%sponge_north_timescale, &
                nm1%slat, nm1%nlat, &
				mp1%dims,mp1%id, world_process, mp1%rank, mp1%ring_comm)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Terminate MPI											    		   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call MPI_Finalize ( mp1%error )
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end program main



