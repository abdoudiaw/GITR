
module boundaries
    use share 
    use utilities
    implicit none

    contains

    subroutine boundaryConditions()
        use share
        implicit none
        ! create file to store particles that leave the simulation and append to it
        open(unit=13, file="distributionZmin.out", status="unknown", action="write", position="append")
        do j = 1, npart
            if (x(j) > xmax .or. x(j) < xmin .or. y(j) > ymax .or. y(j) < ymin .or. z(j) > zmax .or. z(j) < zmin) then
            isActive(j) = .false.
            write ( *, *) " particle ", j, " is out of bounds"
            write(13,*) x(j), y(j), z(j), vx(j), vy(j), vz(j)
            end if
        end do
    end subroutine boundaryConditions
end module boundaries
