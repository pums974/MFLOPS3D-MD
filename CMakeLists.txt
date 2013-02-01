# -------------------------------------------------------------------------------------------------------------------------
# General options
# -------------------------------------------------------------------------------------------------------------------------

# Top-Level CmakeLists.txt
PROJECT( mflops3d-md Fortran)

#enable_language (Fortran)

# Check version of cmake
cmake_minimum_required(VERSION 2.6)

# to see command lines
#SET(CMAKE_VERBOSE_MAKEFILE ON)

# to add Color on CMake Process
#SET(CMAKE_COLOR_MAKEFILE ON)

# 
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/includes)

# Type of build 
message( STATUS "Build type: " ${CMAKE_BUILD_TYPE})

# -------------------------------------------------------------------------------------------------------------------------
# Compiler options
# -------------------------------------------------------------------------------------------------------------------------

# ifort compiler (LINUX)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
message( STATUS "Fortran compiler: " ${CMAKE_Fortran_COMPILER} )
set(NETCDFLIB "-I/opt/netcdf-4.1.2.ifort/include /opt/netcdf-4.1.2.ifort/lib/libnetcdff.a /opt/netcdf-4.1.2.ifort/lib/libnetcdf.a")
set(NETCDFLIB "${NETCDFLIB}  -L/usr/lib64/ -lhdf5_hl -lhdf5 -lz -lcurl")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -pc64 -xW -axp -cm -w -vec_report0 -assume byterecl -mcmodel=medium -heap-arrays ")
set(CMAKE_Fortran_FLAGS_DEBUG " -O0 -pc64 -g -assume byterecl -mcmodel=medium  -heap-arrays -check ")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ")
set(LAPACKLIB "-Wl,--start-group /opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/mkl/lib/intel64/libmkl_core.a /opt/intel/mkl/lib/intel64/libmkl_sequential.a  -Wl,--end-group")
set(DOUBLE "-r8")
include_directories(/opt/netcdf-4.1.2.ifort/include)
endif ()

# gfortran compiler (LINUX)
if (CMAKE_Fortran_COMPILER MATCHES "gfortran")
message( STATUS "Fortran compiler: " ${CMAKE_Fortran_COMPILER} )
set(NETCDFLIB "-I/opt/netcdf/netcdf4/include -L/opt/netcdf/netcdf4/lib -lnetcdff -L/opt/pnetcdf/pnetcdf-1.2.0/lib -L/opt/netcdf/netcdf4/hdf5.1.8.8/lib -lnetcdf -lpnetcdf -lhdf5_hl -lhdf5 -lm -lz -lcurl")
set(NETCDFLIB "${NETCDFLIB}")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -ffast-math -funroll-loops -std=legacy -fopenmp")
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -Warray-bounds -fbounds-check -std=legacy -Wall -Wno-unused ")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ")
set(DOUBLE "-fdefault-real-8")
#set(LAPACKLIB "-Wl,--start-group /usr/lib64/libblas.a /usr/lib64/liblapack.a  -Wl,--end-group")
#set(LAPACKLIB "-Wl,--start-group /opt/acml/acml-4.4.0.gfortran/gfortran64/lib/libacml.a  -Wl,--end-group")
set(LAPACKLIB "-Wl,--start-group /opt/acml/acml-5.0.0.gfortran/gfortran64_mp/lib/libacml_mp.a  -Wl,--end-group")
include_directories(/usr/lib64/gfortran/modules/)
endif ()

# gfortran compiler (LINUX)
if (CMAKE_Fortran_COMPILER MATCHES "mpif90")
#if (CMAKE_Fortran_COMPILER MATCHES "tau_f90.sh")
message( STATUS "Fortran compiler: " ${CMAKE_Fortran_COMPILER} )
set(NETCDFLIB "-I/opt/netcdf/netcdf4/include -L/opt/netcdf/netcdf4/lib -lnetcdff -L/opt/pnetcdf/pnetcdf-1.2.0/lib -L/opt/netcdf/netcdf4/hdf5.1.8.8/lib -lnetcdf -lpnetcdf -lhdf5_hl -lhdf5 -lm -lz -lcurl")
set(NETCDFLIB "${NETCDFLIB}")
#set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -ffast-math -funroll-loops -fstack-arrays -mcmodel=medium -std=legacy -fopenmp -cpp")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -funroll-loops -fstack-arrays -mcmodel=medium -fopenmp -cpp")
#set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -ffast-math -funroll-loops  -mcmodel=medium -fopenmp -cpp")
#set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -ffast-math -funroll-loops  -mcmodel=medium -cpp")
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -Warray-bounds -fbounds-check -mcmodel=medium -std=legacy -Wall -Wno-unused  -cpp")
#set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -Warray-bounds -fbounds-check -std=legacy  -mcmodel=medium -Wall -Wno-unused -fopenmp -cpp")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ")
set(DOUBLE "-fdefault-real-8")
#set(LAPACKLIB "-Wl,--start-group /usr/lib64/libblas.a /usr/lib64/liblapack.a  -Wl,--end-group")
#set(LAPACKLIB "-Wl,--start-group /opt/acml/acml-5.0.0.gfortran/gfortran64/lib/libacml.a  -Wl,--end-group")
set(LAPACKLIB "-Wl,--start-group /opt/acml/acml-5.0.0.gfortran/gfortran64_mp/lib/libacml_mp.a  -Wl,--end-group")
#set(PETSCLIB "-L/opt/petsc/petsc-3.3-p2.openmpi-1.4.5.gcc-4.7.0.debug/lib/ -lpetsc -lHYPRE  -lsuperlu_dist_3.0 -ldmumps -lmumps_common -lparmetis -lmetis -lpord -lscalapack -lblacs -lstdc++ -lml")
set(PETSCLIB "-L/opt/petsc/petsc-3.3-p2.release/lib/ -lpetsc -lHYPRE -lsuperlu_dist_3.1 -ldmumps -lmumps_common -lparmetis -lmetis -lpord -lscalapack -lblacs -lstdc++ -lml")
set(X11LIB "-lX11")
#include_directories(/usr/lib64/gfortran/modules/)
include_directories(/opt/netcdf/netcdf4/include/ /opt/petsc/petsc-3.3-p2.release/include)
endif ()

if (${CMAKE_BUILD_TYPE} MATCHES "Release")
message( STATUS "Fortran compiler flags:" ${CMAKE_Fortran_FLAGS_RELEASE})
elseif (${CMAKE_BUILD_TYPE} MATCHES "Debug")
message( STATUS "Fortran compiler flags:" ${CMAKE_Fortran_FLAGS_DEBUG})
endif()

message( STATUS "Netcdf library flags:" ${NETCDFLIB})

# -------------------------------------------------------------------------------------------------------------------------
# Libraries
# -------------------------------------------------------------------------------------------------------------------------

# class color_print
SET(COLOR color)
ADD_LIBRARY(${COLOR} STATIC src/utils/color_print.f90 )

# class precision
SET(PREC prec)
ADD_LIBRARY(${PREC} STATIC src/precision/precision.f90 )

# class command line
SET(CMDL cmdl)
ADD_LIBRARY(${CMDL} STATIC src/utils/command_line.f90 ) 
target_link_libraries(${CMDL} ${PREC})

# class io
SET(IO io)
ADD_LIBRARY( ${IO} STATIC src/class_io/class_io.f90) 
target_link_libraries(${IO} ${PREC} ${NETCDFLIB} ${MD})

# class dierckx
#SET(DIERCKX dierckx)
#ADD_LIBRARY( ${DIERCKX} STATIC src/class_filter/dierckx.f) 
#target_link_libraries(${DIERCKX})
#set_target_properties(${DIERCKX} PROPERTIES COMPILE_FLAGS ${DOUBLE})

# class derivatives
SET(DER der)
ADD_LIBRARY( ${DER} STATIC src/class_derivatives/minpack.f
  src/class_derivatives/class_derivatives_coefficient.f90 src/class_derivatives/class_derivatives.f90 )
target_link_libraries(${DER} ${PREC} ${LAPACKLIB})

# class solver_1d
SET(SOLVER_1D solver_1d)
ADD_LIBRARY( ${SOLVER_1D} STATIC src/class_solver/class_solver_coefficient_1d.f90 src/class_solver/class_solver_1d.f90 )
target_link_libraries(${SOLVER_1D} ${PREC} ${LAPACKLIB})

# class solver_3d
SET(SOLVER_3D solver_3d)
ADD_LIBRARY( ${SOLVER_3D} STATIC src/class_solver/class_solver_coefficient_3d.f90 src/class_solver/class_solver_3d.f90 )
target_link_libraries(${SOLVER_3D} ${PREC} ${LAPACKLIB})

# class filter
#SET(FILTER filter)
#ADD_LIBRARY( ${FILTER} STATIC src/class_filter/class_filter.f90)
#target_link_libraries(${FILTER} ${DIERCKX} ${PREC})

# class md
SET(MD md)
ADD_LIBRARY( ${MD} STATIC src/class_md/class_md.f90  )
target_link_libraries(${MD} ${PREC} ${PETSCLIB} ${X11LIB} ${LAPACKLIB})

# class libfield 
SET(FIELD field)
ADD_LIBRARY( ${FIELD} STATIC src/class_field/class_field.f90) 
target_link_libraries(${FIELD} ${PREC} ${IO} ${DER} ${FILTER} ${SOLVER_1D} ${SOLVER_3D} ${MD})

# class mesh
SET(MESH mesh)
ADD_LIBRARY( ${MESH} STATIC src/class_mesh/class_mesh.f90) 
target_link_libraries(${MESH} ${FILTER} ${PREC} ${DER} ${SOLVER_1D} ${SOLVER_3D})

# class navier3d 
SET(NAVIER3D navier3d)
ADD_LIBRARY( ${NAVIER3D} STATIC src/class_navier/class_navier_3D.f90) 
target_link_libraries(${NAVIER3D} ${FIELD} ${MESH} ${COLOR} ${CMDL})

# -------------------------------------------------------------------------------------------------------------------------
# Test executables
# -------------------------------------------------------------------------------------------------------------------------

# test solver_1d
SET(TESTSOLVER_1D testsolver_1d)
ADD_EXECUTABLE( ${TESTSOLVER_1D} src/main/test_solver_1d.f90 src/main/parameters.f90)
target_link_libraries(${TESTSOLVER_1D} ${FIELD} ${MESH} ${LAPACKLIB} ${COLOR})

# test solver_3d
SET(TESTSOLVER_3D testsolver_3d)
ADD_EXECUTABLE( ${TESTSOLVER_3D} src/main/test_solver_3d.f90 src/main/parameters.f90)
target_link_libraries(${TESTSOLVER_3D} ${FIELD} ${MESH} ${LAPACKLIB} ${COLOR} ${CMDL})

# test md
SET(TESTMD testmd)
ADD_EXECUTABLE( ${TESTMD} src/main/test_md.f90 src/main/parameters.f90)
target_link_libraries(${TESTMD} ${MD} ${FIELD} ${MESH} ${LAPACKLIB} ${COLOR} ${CMDL})

# test derivatives
SET(TESTDERO testdero)
ADD_EXECUTABLE( ${TESTDERO} src/main/test_derivatives_omp.f90 src/main/parameters.f90)
target_link_libraries(${TESTDERO} ${FIELD} ${MESH} ${LAPACKLIB} ${COLOR} ${CMDL})

# test derivatives
SET(TESTDER testder)
ADD_EXECUTABLE( ${TESTDER} src/main/test_derivatives.f90 src/main/parameters.f90)
target_link_libraries(${TESTDER} ${FIELD} ${MESH} ${LAPACKLIB} ${COLOR})

# test io
SET(TESTIO testio)
ADD_EXECUTABLE( ${TESTIO} src/main/test_io.f90 src/main/parameters.f90)
target_link_libraries(${TESTIO} ${FIELD} ${MESH} ${COLOR})

# test filter
#SET(TESTFILTER testfilter)
#ADD_EXECUTABLE( ${TESTFILTER} src/main/test_filter.f90)
#target_link_libraries(${TESTFILTER} ${FIELD} ${MESH} ${COLOR})


# navier3d
SET(TESTNAV testnav)
ADD_EXECUTABLE( ${TESTNAV} src/main/navier3d.f90 src/main/parameters.f90)
target_link_libraries(${TESTNAV} ${NAVIER3D} ${MD} ${FIELD} ${MESH} ${COLOR} ${CMDL})


