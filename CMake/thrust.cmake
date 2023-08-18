
# set the thrust backend accelerator
if( GITR_USE_CUDA )

  add_compile_definitions( THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA )
else()

  add_compile_definitions( THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP )

endif()

# If CUDA is specified, use the thrust that ships with the CUDA package
if( GITR_USE_CUDA )
  set( THRUST_INCLUDE_DIR "${CUDAToolkit_INCLUDE_DIRS}" CACHE PATH "" FORCE )

else()

  set( THRUST_INCLUDE_DIR "${prefix}/thrust" CACHE PATH "" FORCE )

  if( NOT EXISTS ${THRUST_INCLUDE_DIR} )
    
    message( "thrust will be downloaded..." )
    set( thrust_url "https://github.com/NVIDIA/thrust.git" )
    set( download_command
         git clone --depth 1 --branch 1.15.0
         ${thrust_url} 
         ${prefix}/thrust )

    ExternalProject_Add( thrust_download
                         PREFIX ${prefix} 
                         DOWNLOAD_COMMAND ${download_command}
                         CONFIGURE_COMMAND ""
                         BUILD_COMMAND ""
                         INSTALL_COMMAND "" )

    endif()
endif()

add_library( thrust INTERFACE )
if( TARGET thrust_download )
  add_dependencies( thrust thrust_download )
endif()

include_directories( ${THRUST_INCLUDE_DIR} )
target_include_directories( thrust INTERFACE ${THRUST_INCLUDE_DIR} )

list( APPEND dependencies thrust )

