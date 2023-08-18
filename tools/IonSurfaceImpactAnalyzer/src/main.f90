program IonSurfaceImpactAnalysis
   use share
   use parametersIO, only: read_input, WriteIonData, WriteInitialConditions, Initialization, Finalization 
   use pusher, only: borisPusher
   use boundaries, only: boundaryConditions
   use utilities, only: CompactArrays
   implicit none

   ! Initialization: read in the input file and initial particles conditions
   call Initialization()
   !
   do loop = 1, steps
      write (*,*) 'step', loop
      call borisPusher()
      call boundaryConditions()
      call CompactArrays()
   !  
      call WriteIonData()
   !
       if (npart == 0) then
          write(*,*) 'No more particles.  Stopping.'
          exit
       end if
   end do ! loop

   call Finalization()

 end program IonSurfaceImpactAnalysis
!
! End of file
 
