# 
# Experimental CMake configuration script for Elmer
# 
# Authors:
#  Sami Ilvonen, CSC - IT Center for Science, Ltd.
#  Phil Weir, NUMA Engineering Services, Ltd.
#
# First public beta release 11th of October, 2013
#
# Note that this version is highly experimental and includes
# only a very limited subset of features that are supported
# by the original autotools scripts.
#

PROJECT(Elmer Fortran C CXX)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/Modules")

ENABLE_TESTING()

SET(ELMER_FEM_MAJOR_VERSION 7)
SET(ELMER_FEM_MINOR_VERSION 0)
SET(ELMER_FEM_VERSION
  ${ELMER_FEM_MAJOR_VERSION}.${ELMER_FEM_MINOR_VERSION})

SET(VERSION ${ELMER_FEM_VERSION})
SET(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
SET(MUMPS_INCLUDE "/usr/include")
FIND_PACKAGE(MPI REQUIRED)
ADD_DEFINITIONS(${MPI_Fortran_COMPILE_FLAGS})
INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${MUMPS_INCLUDE})
LINK_DIRECTORIES(${MPI_Fortran_LIBRARIES})

#FIND_PACKAGE(CUDA)
#SET(HAVE_CUDA ${CUDA_FOUND})
#IF(HAVE_CUDA)
#   ADD_SUBDIRECTORY(cuda_hutiter)
#ENDIF(HAVE_CUDA)

FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)

#FIND_PACKAGE(Trilinos)
#SET(HAVE_TRILINOS ${Trilinos_FOUND})
FIND_PACKAGE(Mumps)
SET(HAVE_MUMPS ${Mumps_FOUND})
MESSAGE(${HAVE_MUMPS})


# Check if Fortran compiler supports contiguous keyword
INCLUDE(testContiguous)
IF(CMAKE_Fortran_COMPILER_SUPPORTS_CONTIGUOUS)
  ADD_DEFINITIONS(-DCONTIG=,CONTIGUOUS)
ELSE()
  ADD_DEFINITIONS(-DCONTIG=)
ENDIF()

SET(FC_STDCALLBULL " ")
SET(FC_FUNC FC_GLOBAL)
SET(FC_FUNC_ FC_GLOBAL_)

INCLUDE(FortranCInterface)
FortranCInterface_HEADER(FCMangle.h MACRO_NAMESPACE "FC_")
FortranCInterface_VERIFY(CXX)

SET(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/bin/fortran CACHE PATH "Build directory for Fortran modules."
)
SET(CMAKE_Fortran_MODULE_DESTINATION
    ${CMAKE_INSTALL_PREFIX}/include/elmer/f90 CACHE PATH "Destination for Fortran modules."
)

EXEC_PROGRAM(
    "git"
    ${CMAKE_CURRENT_SOURCE_DIR}
    ARGS "describe" "--tags" "--always"
    OUTPUT_VARIABLE GIT_REVISION)
FILE(WRITE ${PROJECT_BINARY_DIR}/gitrev "${GIT_REVISION}")
INSTALL(FILES ${PROJECT_BINARY_DIR}/gitrev DESTINATION share/elmersolver)

ADD_SUBDIRECTORY(umfpack)
ADD_SUBDIRECTORY(hutiter)
ADD_SUBDIRECTORY(eio)
ADD_SUBDIRECTORY(matc)
ADD_SUBDIRECTORY(fem)
ADD_SUBDIRECTORY(mathlibs)
ADD_SUBDIRECTORY(elmergrid)

MESSAGE(STATUS "------------------------------------------------")
MESSAGE(STATUS "  BLAS library:   " ${BLAS_LIBRARIES})
MESSAGE(STATUS "  LAPACK library: " ${LAPACK_LIBRARIES})
MESSAGE(STATUS "------------------------------------------------")
MESSAGE(STATUS "  MPI Fortran:    " ${MPI_Fortran_FOUND})
MESSAGE(STATUS "  MPI Fortran compiler:    " ${MPI_Fortran_COMPILER} 
  ${MPI_Fortran_COMPILER_FLAGS})
MESSAGE(STATUS "  MPI Fortran include dir: " ${MPI_Fortran_INCLUDE_PATH})
MESSAGE(STATUS "  MPI Fortran libraries:   " ${MPI_Fortran_LIBRARIES})
MESSAGE(STATUS "  MPI Fortran link flags:  " ${MPI_Fortran_LINK_FLAGS})
MESSAGE(STATUS "------------------------------------------------")
MESSAGE(STATUS "  MPI C:             " ${MPI_C_FOUND})
MESSAGE(STATUS "  MPI C compiler:    " ${MPI_C_COMPILER} 
  ${MPI_C_COMPILER_FLAGS})
MESSAGE(STATUS "  MPI C include dir: " ${MPI_Fortran_INCLUDE_PATH})
MESSAGE(STATUS "  MPI C libraries:   " ${MPI_Fortran_LIBRARIES})
MESSAGE(STATUS "  MPI C flags:       " ${MPI_Fortran_LINK_FLAGS})
MESSAGE(STATUS "------------------------------------------------")
