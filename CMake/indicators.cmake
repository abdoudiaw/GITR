# Add indicators

set( INDICATORS_INCLUDE_DIR 
     "${prefix}/indicators_install/include" 
     CACHE PATH "" FORCE )

set( INDICATORS_LIBRARY 
     "${prefix}/indicators_install/lib/indicators${suffix}" 
     CACHE FILEPATH "" FORCE )

if( NOT EXISTS ${INDICATORS_INCLUDE_DIR} OR
    NOT EXISTS ${INDICATORS_LIBRARY} )

  message( "indicators will be downloaded..." )

  set( indicators_url "https://github.com/p-ranav/indicators.git" )

  if( EXISTS ${prefix}/indicators )
    set( download_command "" )
  else()
    set( download_command 
         git clone --depth 1 
         ${indicators_url}
         ${prefix}/indicators )
  endif()

  set( configure_command
       ${CMAKE_COMMAND}
       -S ${prefix}/indicators
       -B ${prefix}/indicators_build
       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
       -DCMAKE_INSTALL_PREFIX=${prefix}/indicators_install )

    ExternalProject_Add( indicators_download
                         PREFIX ${prefix}
                         DOWNLOAD_COMMAND ${download_command}
                         CONFIGURE_COMMAND ${configure_command}
                         BUILD_BYPRODUCTS ${INDICATORS_LIBRARY}
                         BUILD_COMMAND ${CMAKE_COMMAND} --build ${prefix}/indicators_build -- -j
                         INSTALL_COMMAND ${CMAKE_COMMAND} --install ${prefix}/indicators_build ) 
endif()

add_library( indicators INTERFACE )

if( TARGET indicators_download )
  add_dependencies( indicators indicators_download )
endif()

include_directories( ${INDICATORS_INCLUDE_DIR} )

target_include_directories( indicators INTERFACE 
                            ${INDICATORS_INCLUDE_DIR} )

target_link_libraries( indicators INTERFACE
                       ${INDICATORS_LIBRARY} )

list( APPEND dependencies indicators )
