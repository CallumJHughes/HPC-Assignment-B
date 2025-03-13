program main
	use mpi
	implicit none
	integer, parameter  :: dp = selected_real_kind(15,300)

	integer :: N, i, k, sysHeight, sysWidth, t, steps
	real(kind=dp), dimension(:,:), allocatable :: system, systemChunk
	real(kind=dp), dimension(:), allocatable :: haloLeft, haloRight
	real(kind=dp) :: neighbourUp, neighbourDown, neighbourRight, neighbourLeft, phi
	real(kind=dp) :: lOuter, lInner, rOuter, rInner, avgBefore, avgAfter, systemSumBefore, systemSumAfter
	integer :: num_procs, num_elements_per_proc, rank, ierror, boundaryRequestSend, boundaryRequestRecv
	logical :: completed, allCompleted

	!****************!
	! Initialise MPI !
	!****************!
	call MPI_Init(ierror)
	if (ierror/=0) print *, "ERROR in initialising MPI"

	call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
	call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierror)

	!*****************************!
	! Assign values to parameters !
	!*****************************!
	steps = 0
	phi = 1000 ! Potential of inner, conductive cylinder

	! Define dimensions of physical system
	rInner = 5.0
	rOuter = 10.0
	lInner = 5.0
	lOuter = 20.0

	sysHeight = int(lOuter) ! Height of grid of data points
	sysWidth = int(rOuter) ! Width of grid of data points
	num_elements_per_proc = (sysHeight*sysWidth) / num_procs ! Number of elements in each data slice

	completed = .FALSE.

	!*******************************************!
	! Define initial system in master processor !
	!*******************************************!
	if (rank == 0) then
		print *, "Master processor (0) speaking"
		print *, "Number of processors: ", num_procs
		print *, "Number of elements per processor: ", num_elements_per_proc
		allocate(system(sysHeight,sysWidth)) ! Only allocate overall system in master processor

		! Initialise system elements to zero		
		system = 0 

		! Assign elements  holding information on inner metal conductor to phi
		do i = 1, int(lInner)
			do k = 1, int(rInner)
				system(i,k) = phi
			end do
		end do

		do i = 1, sysHeight
			print *, system(i,:)
		end do

		print *, "=============================================================="
	end if

	!************************************************************************!
	! Define buffers holding a subset of the system array for each processor !
	!************************************************************************!
	allocate(systemChunk(sysHeight,sysWidth/num_procs))
	print *, "Buffers allocate of size", size(systemChunk), "for processor", rank
	print *, "=============================================================="

	!*****************************!
	! Send data to each processor !
	!*****************************!
	call MPI_Scatter(system,num_elements_per_proc,MPI_DOUBLE_PRECISION,systemChunk,num_elements_per_proc,MPI_DOUBLE_PRECISION & 
		& ,0,MPI_COMM_WORLD,ierror)
	if (ierror/=0) print *, "ERROR in scattering data"
	if (ierror==0) print *, "Data scattered to processor", rank

	!***************************************!
	! Allocate both halos in each processor !
	!***************************************!
	allocate(haloLeft(sysHeight))
	allocate(haloRight(sysHeight))

	!******************************************************************************!
	! Initialise avgBefore and avgAfter to random values so 'while' loop can start !
	!******************************************************************************!
	!avgBefore = 1
	!avgAfter = 10

	!***************************************************!
	! Run simulation until equilibrium has been reached !
	!***************************************************!
	time_loop: do t = 1, 1
		steps = steps + 1

		!if ((avgAfter-avgBefore) .LT. 1.0E-12_dp) completed = .TRUE.
		call MPI_Allreduce(sum(systemChunk),systemSumBefore,1,MPI_DOUBLE_PRECISION, &
			& MPI_SUM,MPI_COMM_WORLD,ierror)
		if (ierror/=0) stop "ERROR in reducing"

		!**************************************************!
		! Calculate average of system before being updated !
		!**************************************************!
		avgBefore = sum(systemChunk)/size(systemChunk)
		!print *, "Average before",rank,":", avgBefore

		!**********************************************!
		! Sending halo data to neighbouring processors !
		!**********************************************!
		if (rank < num_procs-1) then
			! Send boundary to the right
			call MPI_Isend(systemChunk(:,sysWidth/num_procs),sysHeight,MPI_DOUBLE_PRECISION, &
				& rank+1,1,MPI_COMM_WORLD,boundaryRequestSend,ierror)
			if (ierror/=0) stop "ERROR in sending right boundary data in processor"

			! Receive boundary data from the right
			call MPI_Irecv(haloRight,sysHeight,MPI_DOUBLE_PRECISION,rank+1,2,MPI_COMM_WORLD,boundaryRequestRecv,ierror)
			if (ierror/=0) stop "ERROR in receiving right boundary data in processor"
		end if

		if (rank > 0) then
			! Send boundary data to the left
			call MPI_Isend(systemChunk(:,1),sysHeight,MPI_DOUBLE_PRECISION, &
				& rank-1,2,MPI_COMM_WORLD,boundaryRequestSend,ierror)
			if (ierror/=0) stop "ERROR in sending left boundary data in processor"

			! Receive boundary data from left
			call MPI_Irecv(haloLeft,sysHeight,MPI_DOUBLE_PRECISION,rank-1,1,MPI_COMM_WORLD,boundaryRequestRecv,ierror)
			if (ierror/=0) stop "ERROR in receiving left boundary data"
		end if

		!************************************************!
		! Check it is safe to update halo-dependent data !
		!************************************************!
		if (rank < num_procs-1) then
			! Check data was sent to the right
			call MPI_Wait(boundaryRequestSend,MPI_STATUS_IGNORE,ierror)
			if (ierror/=0) stop "ERROR with MPI_Wait sending right"

			! Check data was received from right
			call MPI_Wait(boundaryRequestRecv,MPI_STATUS_IGNORE,ierror)
			if (ierror/=0) stop "ERROR with MPI_Wait reveiving right"
		endif
		if (rank > 0) then
			call MPI_wait(boundaryRequestRecv,MPI_STATUS_IGNORE,ierror)
			if (ierror/=0) stop "ERROR with MPI_Wait receiving left"
			call MPI_wait(boundaryRequestSend,MPI_STATUS_IGNORE,ierror)
			if (ierror/=0) stop "ERROR with MPI_Wait sending left"
		endif

		!**************************************!
		! Assigns halo data to edge processors !
		!**************************************!
		if (rank == 0) haloLeft = 0 ! Assigns halo left data in left processor to 0 (Left-hand side boundary of system)

		if (rank == num_procs-1) haloRight = 0 ! Assigns halo right data in right processor to 0 (Right-hand side boundary of system)

		print *, "Halo right for",rank,":",haloRight(:)
		print *, "Halo left for",rank,":",haloLeft(:)
		

		!******************************!
		! Update halo-independent data !
		!******************************!
		do i = 1, sysHeight
			do k = 1, (sysWidth/num_procs)
				if (systemChunk(i,k) == phi) then
					cycle
				else if (k == 1 .AND. rank == 0) then
					call FindNeighbours
					systemChunk(i,k) = UpdatedPotential(CalcPotR0(neighbourUp,neighbourDown,neighbourRight),systemChunk(i,k))
				else
					call FindNeighbours
					if (rank == 0 .AND. k == sysWidth/num_procs) then
						print *, "NR for",i,"rank",rank,":",neighbourRight
					else if (rank == 1 .AND. k == 1) then
						print *, "NL for",i,"rank",rank,":",neighbourLeft
					end if
					systemChunk(i,k) = UpdatedPotential(CalcPotU(neighbourUp,neighbourDown,neighbourRight,neighbourLeft),systemChunk(i,k))
				end if
			end do
		end do				

		

		!***************************!
		! Upate halo-dependent data !
		!***************************!
		!do i = 1, sysHeight
		!	! Update left-hand side boundary data
		!	k = 1
		!	if (systemChunk(i,k) == phi) then
		!		cycle
		!	else
		!		call FindNeighbours
		!		if (rank == 0) then
		!			systemChunk(i,k) = UpdatedPotential(CalcPotR0(neighbourUp,neighbourDown,neighbourRight),systemChunk(i,k))
		!		else if (rank /= 0) then
		!			systemChunk(i,k) = UpdatedPotential(CalcPotU(neighbourUp,neighbourDown,neighbourRight,neighbourLeft),systemChunk(i,k))
		!		end if
		!	end if
!
		!	! Update right-hand side boundary data
		!	k = sysWidth/num_procs
		!	if (systemChunk(i,k) == phi) then
		!		cycle
		!	else
		!		call FindNeighbours
		!		systemChunk(i,k) = UpdatedPotential(CalcPotU(neighbourUp,neighbourDown,neighbourRight,neighbourLeft),systemChunk(i,k))
		!	end if
		!end do

		!*************************************************!
		! Calculate average of system after being updated !
		!*************************************************!
		avgAfter = sum(systemChunk)/size(systemChunk)
		!print *, "Average after",rank,":", avgAfter

		! -14
		!if ((avgAfter-avgBefore) .LT. 1.0E-12_dp) completed = .TRUE.
		!call MPI_Allreduce(completed,allCompleted,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierror)
		!if (ierror/=0) stop "ERROR in reducing"
		!if (allCompleted) exit time_loop

		call MPI_Allreduce(sum(systemChunk),systemSumAfter,1,MPI_DOUBLE_PRECISION, &
			& MPI_SUM,MPI_COMM_WORLD,ierror)
		if (ierror/=0) stop "ERROR in reducing"

		if (systemSumAfter-systemSumBefore .LT. 1.0E-4) exit time_loop

	end do time_loop

	print *, "Steps:",steps

	if (allCompleted) then
		print *, "Processor",rank,"done"
	else
		print *, "Processor",rank,"not done - need more iterations"
	end if

	!call MPI_Barrier(MPI_COMM_WORLD,ierror)
	!if (ierror/=0) stop "ERROR in barrier"
	!print *, "Processor",rank,"done"

	!*******************************!
	! Gather data to root processor !
	!*******************************!
	call MPI_Gather(systemChunk,num_elements_per_proc,MPI_DOUBLE_PRECISION, &
		& system,num_elements_per_proc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
	if (ierror/=0) stop "ERROR in gathering data"

	if (rank == 0) then
		do i = 1, sysHeight
			print *, system(i,:)
		end do
	end if

	!***************!
	! Finalizes MPI !
	!***************!
	call MPI_Finalize(ierror)
	if (ierror/=0) print *, "ERROR in finalizing MPI"

contains

	!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!! FUNCTIONS !!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!

	function CalcPotU(neighbourUp,neighbourDown,neighbourRight,neighbourLeft)
	!!! Function to calculate the potential using the potential of the current point's neigbours !!!
		real(kind=dp) :: CalcPotU
		real(kind=dp), intent(in) :: neighbourUp, neighbourDown, neighbourRight, neighbourLeft

		CalcPotU = 0.25*(neighbourUp+neighbourDown+neighbourRight+neighbourLeft) + (1.0/(8.0*k))*(neighbourRight-neighbourLeft)
	end function

	function CalcPotR0(neighbourUp,neighbourDown,neighbourRight)
	!!! Function to calculate the potential when k = 0 (i.e. r = 0) !!!
		real(kind=dp) :: CalcPotR0
		real(kind=dp), intent(in) :: neighbourUp, neighbourDown, neighbourRight

		CalcPotR0 = (2.0/3.0)*(neighbourRight) + (1/6)*(neighbourUp+neighbourDown)
	end function

	function UpdatedPotential(U, phiOld)
	!!! Function to update the new value for the potential at that point !!!
		real(kind=dp) :: UpdatedPotential
		real(kind=dp), intent(in) :: U, phiOld
		real(kind=dp), parameter :: w = 1 ! Defines the rate of convergence to the solution

		UpdatedPotential = phiOld + w*(U - phiOld)
	end function

	!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!! SUBROUTINES !!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!

	subroutine FindNeighbours
	!!! Subroutine to assign values to neighbourUp, neighbourDown, neighbourRight, and neighbourLeft !!!
	!!! Also, checks for boundary conditions and assigns values accordingly !!!
		if (i .EQ. sysHeight) then
			neighbourUp = 0
		else
			neighbourUp = systemChunk(i+1,k)
		end if

		if (i .EQ. 1) then
			neighbourDown = 0
		else
			neighbourDown = systemChunk(i-1,k)
		end if

		if (k .EQ. sysWidth/num_procs) then
			neighbourRight = haloRight(i)
		else
			neighbourRight = systemChunk(i,k+1)
		end if

		if (k .EQ. 1) then
			neighbourLeft = haloLeft(i)
		else
			neighbourLeft = systemChunk(i,k-1)
		end if
	end subroutine

end program main
