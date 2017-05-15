	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>variables for the shallow water model on sphere
    module variables
    use nrtype
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>variables and types for the thermal cloud model

    implicit none
    
		!>@brief
		!>main model prognostic variables
        type grid
            ! variables for grid
            integer(i4b) :: ip, jp, ntim, o_halo, ipstart, jpstart
            integer(i4b), dimension(2) :: coords
            real(sp) :: f, re, g, rho, dt
            real(sp), dimension(:,:), allocatable :: f_cor, h, hs, u,v, height, &
            										 dx, dy, x, y, &
    				recqdp, recqdp_s, recqdq_s, redq_s, redq, &
    				recq, cq_s, cq, dp1, dq
            										 
            real(sp), dimension(:), allocatable :: phi, theta, phin, thetan, u_nudge, &
            	dphi, dtheta, dphin, dthetan
        end type grid



    											
				
	

		!>@brief
		!>variables for namelist input
        type namelist_input
            character (len=200) :: inputfile='input'
            character (len=200) :: outputfile='output'
            logical :: add_random_height_noise, &
            			initially_geostrophic, &
            			viscous_dissipation, &
            			dissipate_h, nudge, restart
            integer(i4b) :: initial_winds, ip, jp
            real(sp) :: wind_factor, wind_shift, wind_reduce, vis, &
            			runtime, dt, output_interval, &
            			grav, rho, Re, rotation_period_hours, scale_height, &
            			slat, nlat, slat_thresh, nlat_thresh, nudge_timescale, &
            			u_jet, theta_jet, h_jet
        end type namelist_input



		!>@brief
		!>variables for mpi
        type mpi_vars
        	integer(i4b) :: rank, id, error, dx, dy
        	real(sp) :: wtime
        	logical, dimension(2) :: periods
        	logical :: reorder=.true.
        	integer(i4b) :: ndim=2
        	integer(i4b), dimension(2) :: dims
        	integer(i4b) :: ring_comm
        end type mpi_vars



        
		!>@brief
		!>variables for NetCDF file output
        type io
            ! variables for io
            integer(i4b) :: ncid, varid, x_dimid, y_dimid, z_dimid, &
                            dimids(2), a_dimid, xx_dimid, yy_dimid, &
                            zz_dimid, i_dimid, j_dimid, k_dimid, nq_dimid, nprec_dimid
            integer(i4b) :: icur=1
            logical :: new_file=.true.
        end type io




		! declare a namelist type
		type(namelist_input) :: nm1
		! declare an mpi type
		type(mpi_vars) :: mp1
        ! declare a grid type
        type(grid) :: grid1
        ! declare an io type
        type(io) :: io1


		integer(i4b) :: world_process=0  
		integer(i4b) :: MPI_INTEGER9      
        
		
    end module variables




