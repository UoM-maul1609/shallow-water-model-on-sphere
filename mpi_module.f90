	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>mpi routines for shallow water model
    module mpi_module
    use numerics_type
    use mpi
    implicit none
    
#if VAR_TYPE==0
		integer(i4b), parameter :: MPIREAL=MPI_REAL4
#endif
#if VAR_TYPE==1
		integer(i4b), parameter :: MPIREAL=MPI_REAL8
#endif
    private
    public :: mpi_define, block_ring, exchange_halos, MPIREAL
    
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
		real(wp), intent(inout), &
			 dimension(1-o_halo:o_halo+ipp,1-o_halo:o_halo+jpp) :: array

		integer(i4b) :: error, nbrleft, nbrright, nbrbottom, nbrtop
		integer(i4b) :: i,j,k,nxbuf,nybuf
		real(wp), allocatable :: sendbuf(:), recvbuf(:)

		! Neighbour ranks returned by MPI_CART_SHIFT are ranks in comm2d, so all
		! communication below must also use comm2d.  This is essential when the
		! Cartesian communicator was created with reorder=.true.
		call MPI_CART_SHIFT(comm2d, 0, 1, nbrleft, nbrright, error)
		call MPI_CART_SHIFT(comm2d, 1, 1, nbrbottom, nbrtop, error)

		! Longitude halo columns are strided Fortran array sections.  Pack them
		! explicitly rather than passing a non-contiguous column to MPI as though
		! it were contiguous memory.
		nybuf=jpp*o_halo
		allocate(sendbuf(nybuf),recvbuf(nybuf))

		k=0
		do j=1,jpp
			do i=ipp-o_halo+1,ipp
				k=k+1; sendbuf(k)=array(i,j)
			enddo
		enddo
		call MPI_Sendrecv(sendbuf,nybuf,MPIREAL,nbrright,10, &
		                  recvbuf,nybuf,MPIREAL,nbrleft,10,comm2d,MPI_STATUS_IGNORE,error)
		k=0
		do j=1,jpp
			do i=1-o_halo,0
				k=k+1; array(i,j)=recvbuf(k)
			enddo
		enddo

		k=0
		do j=1,jpp
			do i=1,o_halo
				k=k+1; sendbuf(k)=array(i,j)
			enddo
		enddo
		call MPI_Sendrecv(sendbuf,nybuf,MPIREAL,nbrleft,11, &
		                  recvbuf,nybuf,MPIREAL,nbrright,11,comm2d,MPI_STATUS_IGNORE,error)
		k=0
		do j=1,jpp
			do i=ipp+1,ipp+o_halo
				k=k+1; array(i,j)=recvbuf(k)
			enddo
		enddo
		deallocate(sendbuf,recvbuf)

		! Latitude edge strips are also packed explicitly because 1:ipp does not
		! span the full declared first dimension when longitude halos are present.
		nxbuf=ipp*o_halo
		allocate(sendbuf(nxbuf),recvbuf(nxbuf))

		k=0
		do j=jpp-o_halo+1,jpp
			do i=1,ipp
				k=k+1; sendbuf(k)=array(i,j)
			enddo
		enddo
		call MPI_Sendrecv(sendbuf,nxbuf,MPIREAL,nbrtop,20, &
		                  recvbuf,nxbuf,MPIREAL,nbrbottom,20,comm2d,MPI_STATUS_IGNORE,error)
		if (nbrbottom /= MPI_PROC_NULL) then
			k=0
			do j=1-o_halo,0
				do i=1,ipp
					k=k+1; array(i,j)=recvbuf(k)
				enddo
			enddo
		endif

		k=0
		do j=1,o_halo
			do i=1,ipp
				k=k+1; sendbuf(k)=array(i,j)
			enddo
		enddo
		call MPI_Sendrecv(sendbuf,nxbuf,MPIREAL,nbrbottom,21, &
		                  recvbuf,nxbuf,MPIREAL,nbrtop,21,comm2d,MPI_STATUS_IGNORE,error)
		if (nbrtop /= MPI_PROC_NULL) then
			k=0
			do j=jpp+1,jpp+o_halo
				do i=1,ipp
					k=k+1; array(i,j)=recvbuf(k)
				enddo
			enddo
		endif
		deallocate(sendbuf,recvbuf)

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
	
	