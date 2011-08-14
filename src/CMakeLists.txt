cmake_minimum_required(VERSION 2.8)
# Note, the optional FORTRAN specification is only requied due to LAPACK below.  If we were to use
# CLAPACK (i.e. f2c) exclusivly we wouldn't need this.  I didn't want to write a FindF2C.cmake file
# so we'll just keep with the fortran libraries.
project (fit_rv C Fortran)

set(CMAKE_BUILD_TYPE Release)

#set(MULTINEST_ROOT_DIR /home/bkloppen/lib/MultiNest_v2.12)
find_package(MULTINEST REQUIRED)
include_directories(${MULTINEST_INCLUDE_DIRS})
set(LIBS ${LIBS} ${MULTINEST_LIBRARIES})

if(${CMAKE_Fortran_COMPILER_ID} EQUAL "Intel")
    set(USING_INTEL True)
endif()


set(MPI_COMPILER mpif90)
find_package(MPI)
if(MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
    set(LIBS ${LIBS} ${MPI_LIBRARY})
    set(CMAKE_Fortran_COMPILER ${MPI_COMPILER})
    set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} "-nofor_main")
endif(MPI_FOUND)

# We need lapack, but the cluster doesn't have it installed in the
# right location, workaround
#set(LIB ${LIB} /usr/lib64 -llapack)

find_package(LAPACK REQUIRED)
include_directories(${LAPACK_INCLUDE_DIRS})
set(LIBS ${LIBS} ${LAPACK_liBRARIES})

# Requires the math library
set(LIBS ${LIBS} m)

#file(GLOB_RECURSE SOURCE . *.c)
file(GLOB SOURCE . *.c)

add_executable (fit_rv ${SOURCE})

SET_TARGET_PROPERTIES(fit_rv PROPERTIES LINKER_LANGUAGE Fortran)

target_link_libraries(fit_rv ${LIBS})
