# CMake project file for IonSurfaceImpactAnalyzer

cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(IonSurfaceImpactAnalyzer VERSION 28.0)

ENABLE_LANGUAGE( Fortran )


add_executable(lmp share.f90
                        utilities.f90
                        parametersIO.f90
                        pusher.f90
                        boundaries.f90
                        main.f90)


target_include_directories(lmp PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )
