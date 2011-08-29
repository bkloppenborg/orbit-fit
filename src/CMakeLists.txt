cmake_minimum_required(VERSION 2.8)
# Note, the optional FORTRAN specification is only requied due to LAPACK below.  If we were to use
# CLAPACK (i.e. f2c) exclusivly we wouldn't need this.  I didn't want to write a FindF2C.cmake file
# so we'll just keep with the fortran libraries.
project (orbit_fit CXX Fortran)

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)

set(MPI_COMPILER mpif90)
find_package(MPI)
if(MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
    set(CMAKE_Fortran_COMPILER ${MPI_COMPILER})
    set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} "-DMPI")
endif(MPI_FOUND)

# We need lapack, but the cluster doesn't have it installed in the
# right location, workaround
#list(APPEND LIBS /usr/lib64 -llapack)

find_package(LAPACK REQUIRED)
include_directories(${LAPACK_INCLUDE_DIRS})

# Requires the math library
set(LIBS ${LIBS} m)

# Explicitly include the orbit and textio libraries/include directories
include_directories(.)
include_directories(${CMAKE_SOURCE_DIR}/lib/textio)
include_directories(${CMAKE_SOURCE_DIR}/lib/orbital_motion)
include_directories(${CMAKE_SOURCE_DIR}/lib/multinest)

# Now specify the source and output for the individual binaries:
file(GLOB GEN_DATA_SOURCE gen_data/*.cpp gen_data/*.c)
add_executable(gendata common.cpp ${GEN_DATA_SOURCE})
SET_TARGET_PROPERTIES(gendata PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(gendata orbit textio)

file(GLOB FIT_RV_SOURCE fit_rv/*.cpp)
add_executable(fitrv common.cpp multinest_inf.cpp ${FIT_RV_SOURCE})
SET_TARGET_PROPERTIES(fitrv PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(fitrv orbit textio multinest ${MPI_LIBRARY} ${LAPACK_LIBRARIES})

file(GLOB FIT_AST_SOURCE fit_astrometry/*.cpp)
add_executable(fitast common.cpp multinest_inf.cpp ${FIT_AST_SOURCE})
SET_TARGET_PROPERTIES(fitast PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(fitast orbit textio multinest ${MPI_LIBRARY} ${LAPACK_LIBRARIES})

