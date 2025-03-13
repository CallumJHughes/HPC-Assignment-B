program main
	use mpi
	implicit none
	integer, parameter  :: dp = selected_real_kind(15,300)

	!***************!
	! MPI Variables !
	!***************!
	integer 	  :: num_procs, num_elements_per_proc, rank, ierror, boundaryRequestSend, boundaryRequestRecv
	real(kind=dp) :: systemSumBefore, systemSumAfter, start, finish

	!*********************!
	! Algorithm Variables !
	!*********************!
	real(kind=dp), dimension(:,:), allocatable :: system, systemChunk
	real(kind=dp), dimension(:), allocatable   :: haloLeft, haloRight
	real(kind=dp) 							   :: neighbourUp, neighbourDown, neighbourRight, neighbourLeft
	real(kind=dp) 							   :: lOuter, lInner, rOuter, rInner, phi
	integer 								   :: i, k, t, steps
	integer									   :: N, sysHeight, sysWidth, scaleFactor
	logical 								   :: allCompleted

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
	scaleFactor = 1
	phi = 1000 ! Potential of inner, conductive cylinder

	! Define dimensions of physical system
	rInner = 1.0
	rOuter = 10.0
	lInner = 5.0
	lOuter = 20.0

	sysHeight = int(lOuter)*scaleFactor ! Height of grid of data points
	sysWidth = int(rOuter)*scaleFactor ! Width of grid of data points
	num_elements_per_proc = (sysHeight*sysWidth) / num_procs ! Number of elements in each data slice

	allCompleted = .FALSE. ! Check to see if max limit of steps has been reached

	!*******************************************!
	! Define initial system in master processor !
	!*******************************************!
	if (rank == 0) then
		print *, "=============================================================="
		print *, "Master processor (0) speaking"
		print *, "Number of processors: ", num_procs
		print *, "Number of elements per processor: ", num_elements_per_proc
		allocate(system(sysHeight,sysWidth)) ! Only allocate overall system in master processor

		system = 0 ! Initialise system elements to zero 

		! Assign elements holding information on inner metal conductor to phi
		do i = 1, int(lInner)*scaleFactor
			do k = 1, int(rInner)*scaleFactor
				system(i,k) = phi
			end do
		end do

		print *, "=============================================================="
	end if

	!************************************************************************!
	! Define buffers holding a subset of the system array for each processor !
	!************************************************************************!
	allocate(systemChunk(sysHeight,1+(rank*(sysWidth/num_procs)):(sysWidth/num_procs)*(rank+1)))
	print *, "Buffers allocate of size", size(systemChunk), "for processor", rank

	!*****************************!
	! Send data to each processor !
	!*****************************!
	call MPI_Scatter(system,num_elements_per_proc,MPI_DOUBLE_PRECISION,systemChunk,num_elements_per_proc,MPI_DOUBLE_PRECISION & 
		& ,0,MPI_COMM_WORLD,ierror)
	if (ierror/=0) print *, "ERROR in scattering data"
	if (ierror==0) print *, "Data scattered to processor", rank
	print *, "=============================================================="

	!***************************************!
	! Allocate both halos in each processor !
	!***************************************!
	allocate(haloLeft(sysHeight))
	allocate(haloRight(sysHeight))

	start = MPI_Wtime()

	!***************************************************!
	! Run simulation until equilibrium has been reached !
	!***************************************************!
	time_loop: do t = 1, 1000000
		steps = steps + 1

		!**************************************************!
		! Calculate average of system before being updated !
		!**************************************************!
		call MPI_Allreduce(sum(systemChunk),systemSumBefore,1,MPI_DOUBLE_PRECISION, &
			& MPI_SUM,MPI_COMM_WORLD,ierror)
		if (ierror/=0) stop "ERROR in reducing"

		!**********************************************!
		! Sending halo data to neighbouring processors !
		!**********************************************!
		if (rank < num_procs-1) then
			! Send boundary to the right
			call MPI_Isend(systemChunk(:,(sysWidth/num_procs)*(rank+1)),sysHeight,MPI_DOUBLE_PRECISION, &
				& rank+1,1,MPI_COMM_WORLD,boundaryRequestSend,ierror)
			if (ierror/=0) stop "ERROR in sending right boundary data in processor"

			! Receive boundary data from the right
			call MPI_Irecv(haloRight,sysHeight,MPI_DOUBLE_PRECISION,rank+1,2,MPI_COMM_WORLD,boundaryRequestRecv,ierror)
			if (ierror/=0) stop "ERROR in receiving right boundary data in processor"
		end if

		if (rank > 0) then
			! Send boundary data to the left
			call MPI_Isend(systemChunk(:,1+(rank*(sysWidth/num_procs))),sysHeight,MPI_DOUBLE_PRECISION, &
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

		!**************************************************!
		! Upate all data in system chunk for each processor!
		!**************************************************!
		do k = 1 + (rank*(sysWidth/num_procs)), ((sysWidth/num_procs)*(rank+1))
			do i = 1, sysHeight
				if (systemChunk(i,k) == phi) then
					cycle
				else if (k == 1 .AND. rank == 0) then
					call FindNeighbours
					systemChunk(i,k) = UpdatedPotential(CalcPotR0(neighbourUp,neighbourDown,neighbourRight),systemChunk(i,k))
				else
					call FindNeighbours
					systemChunk(i,k) = UpdatedPotential(CalcPotU(neighbourUp,neighbourDown,neighbourRight,neighbourLeft),systemChunk(i,k))
				end if
			end do
		end do

		!*************************************************!
		! Calculate average of system after being updated !
		!*************************************************!
		call MPI_Allreduce(sum(systemChunk),systemSumAfter,1,MPI_DOUBLE_PRECISION, &
			& MPI_SUM,MPI_COMM_WORLD,ierror)
		if (ierror/=0) stop "ERROR in reducing"

		!**********************!
		! Test for convergence !
		!**********************!
		if (systemSumAfter-systemSumBefore .LT. 1.0E-3) then
			allCompleted = .TRUE.
			exit time_loop
		end if

		!*******************************************************!
		! Ensure each processor has reached end of current step !
		!*******************************************************!
		call MPI_Barrier(MPI_COMM_WORLD,ierror)
		if (ierror/=0) stop "ERROR in Barrier in loop"

	end do time_loop

	finish = MPI_Wtime() ! Finds finishing time of the algorithm

	print *, "Processor", rank, "took", finish-start, "seconds and", steps, "steps"

	!******************************************************************!
	! Ensures all processors have reached this stage before continuing !
	!******************************************************************!
	call MPI_Barrier(MPI_COMM_WORLD,ierror)
	if (ierror/=0) stop "ERROR in barrier"

	!********************************************************************************!
	! Checks whether algorithm was finished due to convergence or max limit of steps !
	!********************************************************************************!
	if (allCompleted) then
		print *, "Processor",rank,"done"
	else
		print *, "Processor",rank,"not done - need more iterations"
	end if

	!*******************************!
	! Gather data to root processor !
	!*******************************!
	call MPI_Gather(systemChunk,num_elements_per_proc,MPI_DOUBLE_PRECISION, &
		& system,num_elements_per_proc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
	if (ierror/=0) stop "ERROR in gathering data"

	!********************************************************!
	! Prints results to standard output using root processor !
	!********************************************************!
	if (rank == 0) then
		print *, "Point A:", system(int(8.5*scaleFactor),1) 				   ! Point(7.5,0)
		print *, "Point B:", system(int(6*scaleFactor),int(7*scaleFactor))	   ! Point(5,6)
		print *, "Point C:", system(int(13.5*scaleFactor),int(6*scaleFactor))  ! Point(12.5,5)
		do i = 1, sysHeight
			print *, system(i,1)
		end do
		print *, "=============================================================="
		do k = 1, sysWidth
			print *, system(1,i)
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

		CalcPotU = 0.25*(neighbourUp+neighbourDown+neighbourRight+neighbourLeft) + (1.0/(8.0*real(k)))*(neighbourRight-neighbourLeft)
	end function

	function CalcPotR0(neighbourUp,neighbourDown,neighbourRight)
	!!! Function to calculate the potential when k = 0 (i.e. r = 0) !!!
		real(kind=dp) :: CalcPotR0
		real(kind=dp), intent(in) :: neighbourUp, neighbourDown, neighbourRight

		CalcPotR0 = (2.0/3.0)*(neighbourRight) + (1.0/6.0)*(neighbourUp+neighbourDown)
	end function

	function UpdatedPotential(U, phiOld)
	!!! Function to update the new value for the potential at that point !!!
		real(kind=dp) :: UpdatedPotential
		real(kind=dp), intent(in) :: U, phiOld
		real(kind=dp), parameter :: w = 1.5 ! Defines the rate of convergence to the solution

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

		if (k .EQ. (sysWidth/num_procs)*(rank+1)) then
			neighbourRight = haloRight(i)
		else
			neighbourRight = systemChunk(i,k+1)
		end if

		if (k .EQ. (1+(rank*(sysWidth/num_procs)))) then
			neighbourLeft = haloLeft(i)
		else
			neighbourLeft = systemChunk(i,k-1)
		end if
	end subroutine

end program main
