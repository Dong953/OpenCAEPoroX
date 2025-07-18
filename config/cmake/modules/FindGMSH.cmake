# Once done this will define
#  GMSH_FOUND        - System has GMSH
#  GMSH_DIR          - The GMSH directory
#  GMSH_INCLUDE_DIRS - The GMSH include directories
#  GMSH_LIBRARIES    - The libraries needed to use GMSH
#
#  Jindong Wang
#  07/17/2025

set(GMSH_DIR "${GMSH_DIR}")

# Check for header file
find_path(GMSH_INCLUDE_DIRS gmsh.h
	HINTS ${GMSH_DIR}/include $ENV{GMSH_DIR}/include
	DOC "Directory where the GMSH header is located")
mark_as_advanced(GMSH_INCLUDE_DIRS)

# Check for FASP library
find_library(GMSH_LIBRARIES gmsh 
	HINTS ${GMSH_DIR}/lib $ENV{GMSH_DIR}/lib ${PROJECT_SOURCE_DIR}/gmsh/lib
	DOC "The GMSH library")
mark_as_advanced(GMSH_LIBRARIES)

# Collect libraries
set(GMSH_LIBRARIES ${GMSH_LIBRARIES})
# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMSH
	"GMSH could not be found. Check GMSH_DIR."
	GMSH_LIBRARIES GMSH_INCLUDE_DIRS)


