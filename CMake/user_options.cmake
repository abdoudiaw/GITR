# String description for each option
set( description "(no description added - see define_options.cmake)" )
set( GITR_USE_CUDA 0 CACHE STRING "${description}" FORCE )
set( GITR_USE_DOUBLE 1 CACHE STRING "${description}" FORCE )
# set (GITR_USE_INDICATION 0 CACHE STRING "${description}" FORCE )

add_compile_definitions( 
        USE_CUDA=${GITR_USE_CUDA}
        USE_DOUBLE=${GITR_USE_DOUBLE}
        FIELD_ALIGNED_VALUES=0
         )
