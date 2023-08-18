
# Serial modules
add_executable(GITR src/main.cpp)
target_include_directories(GITR PUBLIC include)
if(GITR_USE_CUDA)
    set_source_files_properties(src/main.cpp PROPERTIES LANGUAGE CUDA)
    set_target_properties(GITR PROPERTIES LINKER_LANGUAGE CUDA)
    set_property(TARGET GITR PROPERTY CUDA_SEPARABLE_COMPILATION ON)
    target_compile_options(GITR PRIVATE --expt-relaxed-constexpr)
endif()

# CPU-only targets
set(non_gpu_targets
    utils
    flags
    configInterface
    )

# Conditionally compile as GPU targets
set(gpu_targets
    surfaceReactions
    interpolater
    elementaryProcesses
    sheathModel
    pusher
    geometryCheck
    )

# Include all targets for non-GPU build
if(NOT GITR_USE_CUDA)
    set(non_gpu_targets ${non_gpu_targets} ${gpu_targets})
endif()

# Compile non_gpu_targets
foreach(component IN LISTS non_gpu_targets)
    add_library(${component} src/${component}.cpp)
    target_include_directories(${component} PUBLIC include)
endforeach()

# Compile gpu_targets
if(GITR_USE_CUDA)
    foreach(component IN LISTS gpu_targets)
        add_library(${component} src/${component}.cpp)
        set_source_files_properties(src/${component}.cpp PROPERTIES LANGUAGE CUDA)
        set_target_properties(${component} PROPERTIES COMPILE_FLAGS "-dc")
        target_include_directories(${component} PUBLIC include)
        target_compile_options(${component} PRIVATE --expt-relaxed-constexpr)
    endforeach()
endif()
