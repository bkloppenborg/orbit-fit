# - Try to find MULTINEST.
# Variables used by this module:
#  MULTINEST_ROOT_DIR     - MULTINEST root directory
# Variables defined by this module:
#  MULTINEST_FOUND        - system has MULTINEST
#  MULTINEST_INCLUDE_DIR  - the MULTINEST include directory (cached)
#  MULTINEST_INCLUDE_DIRS - the MULTINEST include directories
#                         (identical to MULTINEST_INCLUDE_DIR)
#  MULTINEST_LIBRARY      - the MULTINEST library (cached)
#  MULTINEST_LIBRARIES    - the MULTINEST libraries
#                         (identical to MULTINEST_LIBRARY)

# Copyright (C) 2009
# ASTRON (Netherlands Institute for Radio Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
#
# This file is part of the LOFAR software suite.
# The LOFAR software suite is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The LOFAR software suite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
#
# $Id$

#TODO: we need to make a special rule for the intel compiler to not call find_package
# for blas/lapack
if(NOT MULTINEST_FOUND)

  find_path(MULTINEST_INCLUDE_DIR multinest.h
    PATHS ${MULTINEST_ROOT_DIR} PATH_SUFFIXES include)
  find_library(MULTINEST_LIBRARY NAMES multinest libmultinest
    PATHS ${MULTINEST_ROOT_DIR} PATH_SUFFIXES lib)

  if(${CMAKE_Fortran_COMPILER_ID} EQUAL "Intel")
    set(USING_INTEL True)
  endif()

  if(NOT ${USING_INTEL})
    find_library(LAPACK_LIBRARY lapack)
    mark_as_advanced(MULTINEST_INCLUDE_DIR MULTINEST_LIBRARY LAPACK_LIBRARY)
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(MULTINEST DEFAULT_MSG
    MULTINEST_LIBRARY LAPACK_LIBRARY MULTINEST_INCLUDE_DIR)

    set(MULTINEST_INCLUDE_DIRS ${MULTINEST_INCLUDE_DIR})
    set(MULTINEST_LIBRARIES "-L"${MULTINEST_ROOT_DIR} ${MULTINEST_LIBRARY} ${LAPACK_LIBRARY})

  else()
    mark_as_advanced(MULTINEST_INCLUDE_DIR MULTINEST_LIBRARY)
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(MULTINEST DEFAULT_MSG
    MULTINEST_LIBRARY MULTINEST_INCLUDE_DIR)

    set(MULTINEST_INCLUDE_DIRS ${MULTINEST_INCLUDE_DIR})
    set(MULTINEST_LIBRARIES ${MULTINEST_LIBRARY} "-llapack -lblas")


  endif()

  #include(FindPackageHandleStandardArgs)
  #find_package_handle_standard_args(MULTINEST DEFAULT_MSG
  #  MULTINEST_LIBRARY LAPACK_LIBRARY MULTINEST_INCLUDE_DIR)

  #set(MULTINEST_INCLUDE_DIRS ${MULTINEST_INCLUDE_DIR})
  #set(MULTINEST_LIBRARIES ${MULTINEST_LIBRARY} ${LAPACK_LIBRARY})

endif(NOT MULTINEST_FOUND)
