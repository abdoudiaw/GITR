# Link previously defined CMake compilation "targets" together as needed

# link source targets
target_link_libraries( surfaceReactions )
target_link_libraries( elementaryProcesses  )
target_link_libraries( interpolater thrust )
target_link_libraries( flags libconfig thrust  )
target_link_libraries( utils libconfig thrust interpolater netcdf  )
target_link_libraries( utils pusher thrust )
target_link_libraries( geometryCheck sheathModel pusher )

# Improvement: Conditionally link based on whether the GITR_USE_<component> clause is enabled
target_link_libraries( GITR 
                      interpolater
                      elementaryProcesses
                      netcdf 
                      libconfig   
                      utils
                      sheathModel
                      pusher
                      surfaceReactions
                      geometryCheck
                      cli11
                      configInterface 
                      flags
                      )
if( GITR_USE_CUDA )
  target_link_libraries( GITR 
                         CUDA::cudart )
  foreach( target IN LISTS gpu_targets gpu_test_targets )
    target_link_libraries( ${target} thrust CUDA::cudart )
  endforeach()

endif()
