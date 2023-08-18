module parametersIO
        use share
        use utilities, only: ran3
        implicit none
        private
        public :: read_input, WriteIonData, WriteInitialConditions, Initialization, Finalization
        ! parameters
        !
        contains
        subroutine read_input()
        implicit none
        !    
        ! Read the initial conditions
        integer :: iunit, ierr
        open(newunit = iunit, file = 'input.in', status = 'old', action = 'read', iostat = ierr)
        if (ierr == 0) then
            read(iunit,*) npart
            read(iunit,*) dens_back
            read(iunit,*) zbar_back
            read(iunit,*) pmass_back
            read(iunit,*) tempi_back
            read(iunit,*) tempe_back
            read(iunit,*) pmass_imp
            read(iunit,*) zbar_imp
            read(iunit,*) tempi_imp
            read(iunit,*) dt
            read(iunit,*) steps
            read(iunit,*) snaps
            read(iunit,*) sheath_factor
            read(iunit,*) Bmagnitude
            read(iunit,*) alpha
            read(iunit, *) xmin
            read(iunit, *) xmax
            read(iunit, *) ymin
            read(iunit, *) ymax
            read(iunit, *) zmin
            read(iunit, *) zmax
            close(iunit)
            write(*,*) 'npart', npart
            write(*,*) 'dens_back', dens_back
            write(*,*) 'zbar_back', zbar_back
            write(*,*) 'pmass_back', pmass_back
            write(*,*) 'tempi_back', tempi_back
            write(*,*) 'tempe_back', tempe_back
            write(*,*) 'pmass_imp', pmass_imp
            write(*,*) 'zbar_imp', zbar_imp
            write(*,*) 'tempi_imp', tempi_imp
            write(*,*) 'dt', dt
            write(*,*) 'steps', steps
            write(*,*) 'snaps', snaps
            write(*,*) 'sheath_factor', sheath_factor
            write(*,*) 'Bmagnitude', Bmagnitude
            write(*,*) 'alpha', alpha
            write(*,*) 'xmin', xmin
            write(*,*) 'xmax', xmax
            write(*,*) 'ymin', ymin
            write(*,*) 'ymax', ymax
            write(*,*) 'zmin', zmin
            write(*,*) 'zmax', zmax

        else
            write(*, *) '*** Error reading the input file *** '
            write(*, *) 'Aborting...'
            stop
        endif
        close(iunit)

        q_over_m  = zbar_imp / pmass_imp
        vti = sqrt(k * tempi_imp * 11600. /( pmass_imp  * m_p))
        omega_ci = q_over_m * Bmagnitude * (e / m_p)
        larmor_radius = vti / omega_ci
        Debye = sqrt(eps0 * k * tempe_back * 11600. / (dens_back * e * e))
        larmor_over_debye_radius = larmor_radius / Debye

        write(*,*) 'larmor_radius', larmor_radius
        write(*,*) 'Debye', Debye
        write(*,*) 'larmor_over_debye_radius', larmor_over_debye_radius
        
    end subroutine read_input

    subroutine WriteInitialConditions()

        implicit none   
        integer ::  iseed
        iseed = -123456
        sigma = 1.0
        ! move this into own module
        do i = 1, npart
            r1 = ran3(iseed)
            r2 = ran3(iseed)
            local_vx(i) = sigma * sqrt(-2.0_8 * log(r1)) * cos(2.0_8 * pi * r2)
            r1 = ran3(iseed)
            r2 = ran3(iseed)
            local_vy(i) = sigma * sqrt(-2.0_8 * log(r1)) * cos(2.0_8 * pi * r2)
            r1 = ran3(iseed)
            r2 = ran3(iseed)
            local_vz(i) = sigma * sqrt(-2.0_8 * log(r1)) * cos(2.0_8 * pi * r2)

            ! Define direction vectors
            xdir(1) = cos(alpha * pi / 180.0_8)
            xdir(2) = 0.0_8
            xdir(3) = -sin(alpha * pi / 180.0_8)

            ydir(1) = 0.0_8
            ydir(2) = 1.0_8
            ydir(3) = 0.0_8

            zdir(1) = sin(alpha * pi / 180.0_8)
            zdir(2) = 0.0_8
            zdir(3) = cos(alpha * pi / 180.0_8)

            ! Calculate velocity components
            vx0(i) = local_vx(i)
            vy0(i) = local_vy(i)
            vz0 (i)= local_vz(i)

            vx(i) = xdir(1) * vx0(i) + ydir(1) * vy0(i) + zdir(1) * vz0(i)
            vy(i) = xdir(2) * vx0(i) + ydir(2) * vy0(i) + zdir(2) * vz0(i)
            vz(i) = xdir(3) * vx0(i) + ydir(3) * vy0(i) + zdir(3) * vz0(i)

            x(i) = 0.0_8
            y(i) = 0.0_8
            z(i) = sheath_factor * larmor_radius 
        end do

        ! Write the initial conditions to file
        open(unit=12, status='unknown', file='particleSource.out')
        ! Open the file for writing
        ! Write the header
        write(12, '(a)') 'ITEM: TIMESTEP'
        write(12, '(i0)') 0
        write(12, '(a)') 'ITEM: NUMBER OF ATOMS'
        write(12, '(i0)') npart
        write(12, '(a)') 'ITEM: BOX BOUNDS pp pp pp'
        write(12, '(2f8.3)') xmin, xmax
        write(12, '(2f8.3)') ymin, ymax
        write(12, '(2f8.3)') zmin, zmax
        write(12, '(a)') 'ITEM: ATOMS id type x y z vx vy vz'
        ! create a unique id for each particle
        do j = 1, npart
            write(12,*) j,1,x(j), y(j), z(j), vx(j), vy(j), vz(j)
        end do
        ! close(12)
    end subroutine WriteInitialConditions


    subroutine WriteIonData()
        use share
        implicit none        
        !
        open(unit=14, status='unknown', file='iqp.out')
        ! Open the file for writing
        ! Write the header
        write(14, '(a)') 'ITEM: TIMESTEP'
        write(14, '(i0)') loop - 1 ! loop starts at 1, so subtract 1
        write(14, '(a)') 'ITEM: NUMBER OF ATOMS'
        write(14, '(i0)') npart
        write(14, '(a)') 'ITEM: BOX BOUNDS pp pp pp'
        write(14, '(2f8.3)') xmin, xmax
        write(14, '(2f8.3)') ymin, ymax
        write(14, '(2f8.3)') zmin, zmax
        write(14, '(a)') 'ITEM: ATOMS id type x y z vx vy vz'
        ! create a unique id for each particle
        do j = 1, npart
            write(14,*) j,1,x(j), y(j), z(j), vx(j), vy(j), vz(j)
        end do
        ! close(14)
        
        end subroutine WriteIonData


          subroutine Initialization()
            use share
            implicit none
            ! Read in the input file and initial conditions
            write(*,*) 'Reading in input file...'
            call read_input()
            allocate(isActive(npart))
            isActive = .true.
            write(*,*) 'Reading in particle starting positions and velocities...'
            ! write done if successful
            call WriteInitialConditions()
            write(*,*) 'Done reading input file...'
          end subroutine Initialization

          subroutine Finalization()
            use share
            implicit none
            write(*,*) 'Done simulating....'
            close(14)
            close(12)
            close(13)
            deallocate(isActive)
          end subroutine Finalization



end module parametersIO


