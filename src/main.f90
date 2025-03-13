program main
!use mpi
implicit none
	integer, parameter  :: dp = selected_real_kind(15,300)

	integer :: N, i, k, sysHeight, sysWidth, t, step
	real(kind=dp), dimension(:,:), allocatable :: system
	real(kind=dp) :: neighbourUp, neighbourDown, neighbourRight, neighbourLeft, phi
	real(kind=dp) :: lOuter, lInner, rOuter, rInner, avgBefore, avgAfter

	step = 0

	! Define parameters for dimensions of conductors
	rInner = 5.0
	rOuter = 10.0
	lInner = 5.0
	lOuter = 20.0

	! Potential of inner, conductive cylinder
	phi = 1000

	! Define height and width of grid point, then allocate (maybe should multiply by 100 later on)
	sysHeight = int(lOuter)*20
	sysWidth = int(rOuter)*20
	allocate(system(sysHeight,sysWidth))

	! Initialise system elements all to 0
	system = 0

	! Assign elements  holding information on inner metal conductor to phi
	do i = 1, int(lInner)
		do k = 1, int(rInner)
			system(i,k) = phi
		end do
	end do

	! Initialise avgBefore and avgAfter to random values so 'while' loop can start
	avgBefore = 1
	avgAfter = 10

	! Update system until equilibrium has been reached
	do while ((avgAfter-avgBefore) .GT. 1.0E-14)
		step = step + 1
		avgBefore = sum(system)/size(system)
		do i = 1, sysHeight
			do k = 1, sysWidth
				if (i .LE. int(lInner) .AND. k .LE. int(rInner)) then
					cycle
				else
					call FindNeighbours
					if (k .EQ. 1) then
						system(i,k) = UpdatedPotential(CalcPotR0(neighbourUp,neighbourDown,neighbourRight),system(i,k))
					else
						system(i,k) = UpdatedPotential(CalcPotU(neighbourUp,neighbourDown,neighbourRight,neighbourLeft),system(i,k))
					end if
				end if
			end do
		end do
	avgAfter = sum(system)/size(system)
	end do

	! Prints final system
	!do i = 1, sysHeight
	!	print *, system(i,:)
	!end do

	!print *, system(0,:)
	print *, "Top right point:", system(sysHeight,sysWidth)
	print *, "Average:", sum(system)/size(system)
	print *, "steps:", step

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
			neighbourUp = system(i+1,k)
		end if

		if (i .EQ. 1) then
			neighbourDown = 0
		else
			neighbourDown = system(i-1,k)
		end if

		if (k .EQ. sysWidth) then
			neighbourRight = 0
		else
			neighbourRight = system(i,k+1)
		end if

		if (k .EQ. 1) then
			neighbourLeft = 0
		else
			neighbourLeft = system(i,k-1)
		end if
	end subroutine

end program main
