# Only OpenMP or Cuda can be specified
if( GITR_USE_CUDA AND GITR_USE_OPENMP )
  message( FATAL_ERROR "Both GITR_USE_CUDA and GITR_USE_OPENMP are set. Please select one" )
endif()

# Handle external dependencies
include( ExternalProject )

set( prefix "${CMAKE_BINARY_DIR}/external" )

if( APPLE )
  set( suffix ".dylib" )
else()
  set( suffix ".so" )
endif()

set( CMAKE_BUILD_WITH_INSTALL_RPATH True )

# ensure shared dependency libs are discoverable at load-time
set( CMAKE_INSTALL_RPATH
     "${CMAKE_BINARY_DIR}"
      "${prefix}/indicators_install/lib"
     "${prefix}/libconfig_install/lib"
     "${prefix}/netcdf-c-install/lib"
     "${prefix}/netcdf-cxx4-install/lib" )

set( dependencies "" )

include( CMake/CLI11.cmake )

if( GITR_USE_CUDA )
  include( CMake/cuda.cmake ) 
  list ( APPEND dependencies cuda::cudart )
endif()
include( CMake/indicators.cmake )
include( CMake/thrust.cmake ) 
include( CMake/libconfig.cmake ) 
include( CMake/netcdf.cmake ) 