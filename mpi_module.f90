	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>mpi routines for shallow water model
    module mpi_module
    use nrtype
    use mpi
    
    private
    public :: mpi_define, block_ring, exchange_halos
    
	contains
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! define some types                                                                  !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>define some types to be used in the model
	!>@param[inout] MPI_INTEGER9_DEF - type to be defined as integer kind=9
	subroutine mpi_define(MPI_INTEGER9_DEF)
		implicit none
		integer(i4b), intent(inout) :: MPI_INTEGER9_DEF
		
		integer(i4b) :: error
		
		call MPI_TYPE_CREATE_F90_INTEGER (9, MPI_INTEGER9_DEF, error)
	end subroutine mpi_define
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! exchange halos for a variable using Cartesian topology                             !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>define some types to be used in the model
	!>@param[in] comm2, id, ipp, jpp, o_halo
	!>@param[inout] array: the array to exchange_halos on
	subroutine exchange_halos(comm2d, id, ipp, jpp, o_halo, array)
		implicit none
		integer(i4b), intent(in) :: comm2d, id, ipp, jpp, o_halo
		real(sp), intent(inout), &
			 dimension(1-o_halo:o_halo+ipp,1-o_halo:o_halo+jpp) :: array
		
		integer(i4b) :: error, nbrleft, nbrright, nbrbottom, nbrtop, tag1, &
						request
		integer(i4b), dimension(MPI_STATUS_SIZE) :: status
		
		! Find the processors neighbours
		call MPI_CART_SHIFT( comm2d, 0, 1, nbrleft, nbrright, error)	
		call MPI_CART_SHIFT( comm2d, 1, 1, nbrbottom, nbrtop, error)
		
! 		print *,id,nbrleft, nbrright, nbrbottom, nbrtop

		
		! now receive data from top, bottom, left, and right
		if (nbrleft /= id) then
			tag1=010
			! send from left (specify destination):
			call MPI_Issend(array(ipp+1-o_halo:ipp,1:jpp), jpp, MPI_REAL8, nbrright, &
				tag1, MPI_COMM_WORLD, request,error)
			! receive from left (specify source):
			call MPI_Recv(array(1-o_halo:0,1:jpp), jpp, MPI_REAL8, nbrleft, &
				tag1, MPI_COMM_WORLD, status,error)
			call MPI_Wait(request, status, error)
		else
			array(1-o_halo:0,1:jpp)=array(ipp+1-o_halo:ipp,1:jpp)
		endif
		
		
		
		if (nbrright /= id) then
			tag1=010
			! send from right (specify destination):
			call MPI_Issend(array(1:o_halo,1:jpp), jpp, MPI_REAL8, nbrleft, &
				tag1, MPI_COMM_WORLD, request,error)
			! receive from right (specify source):
			call MPI_Recv(array(ipp+1:ipp+o_halo,1:jpp), jpp, MPI_REAL8, nbrright, &
				tag1, MPI_COMM_WORLD, status,error)
			call MPI_Wait(request, status, error)
		else
			array(ipp+1:ipp+o_halo,1:jpp)=array(1:o_halo,1:jpp)
		endif
		
		
		if ((nbrtop /= id)) then
			tag1=110
			! send from top (specify destination):
			call MPI_Issend(array(1:ipp,jpp+1-o_halo:jpp), ipp, MPI_REAL8, nbrtop, &
				tag1, MPI_COMM_WORLD, request,error)
			! receive from top (specify source):
			call MPI_Recv(array(1:ipp,1-o_halo:0), ipp, MPI_REAL8, nbrbottom, &
				tag1, MPI_COMM_WORLD, status,error)
			call MPI_Wait(request, status, error)
		else
!			array(1:ipp,1-o_halo:0)=array(1:ipp,jpp-o_halo:jpp)
		endif

		if ((nbrbottom /= id)) then
			tag1=110
			! send from bottom (specify destination):
			call MPI_Issend(array(1:ipp,1:o_halo), ipp, MPI_REAL8, nbrbottom, &
				tag1, MPI_COMM_WORLD, request,error)
			! receive from bottom (specify source):
			call MPI_Recv(array(1:ipp,jpp+1:jpp+o_halo), ipp, MPI_REAL8, nbrtop, &
				tag1, MPI_COMM_WORLD, status,error)
			call MPI_Wait(request, status, error)
		else
!			array(1:ipp,jpp+1:jpp+o_halo)=array(1:ipp,1:o_halo)
		endif
		
		
					
	end subroutine exchange_halos
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	
	
	
	
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Block via ring                                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>block by sending / receiving a message around the ring
	!>@param[in] ring_comm - comm of the cart topology
	!>@param[in] id - id of this process
	!>@param[in] world_process - id of world process
	!>@param[in] rank - rank of mpi job
	subroutine block_ring(ring_comm,id,world_process,rank)
		implicit none
		integer(i4b), intent(in) :: ring_comm, id, world_process, rank
		integer(i4b) :: error, tag1=2010
		character (len=3) :: mesg='Yo!'
		
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! essentially blocks until all processors catch up					   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call MPI_Barrier(ring_comm, error)
! 		if (id .ne. world_process ) then
! 			! processors except 0 are waiting to recv from previous pe:
! 			call MPI_Recv(mesg, len(mesg), MPI_CHARACTER, id-1, &
! 				tag1, ring_comm, MPI_STATUS_IGNORE,error)
! 		endif
! 		if ( (world_process+1) .ne. rank ) then ! so we don't send a message to ourselves!
! 			! processor 0 will send here first (as not waiting)
! 			call MPI_Send(mesg, len(mesg), MPI_CHARACTER, mod(id+1,rank), &
! 					tag1, ring_comm, error)
! 			! lastly receive message from last process
! 			if(id == world_process) then
! 				! processor 0 waiting to recv from last pe in ring
! 				call MPI_Recv(mesg, len(mesg), MPI_CHARACTER, rank-1, &
! 					tag1, ring_comm, MPI_STATUS_IGNORE,error)
! 			endif
! 		endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end subroutine block_ring
	
	
	
	
	end module mpi_module
	
	