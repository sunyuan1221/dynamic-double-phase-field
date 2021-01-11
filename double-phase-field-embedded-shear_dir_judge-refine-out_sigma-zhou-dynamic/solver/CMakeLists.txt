# ==============================================================================
# GEOCENTRIC SOLVER CMake script
# ==============================================================================
MESSAGE(STATUS "Configuring Geocentric")

# ------------------------------------------------------------------------------
# Modify The Following Sections To Describe Your Specific Setup
# ------------------------------------------------------------------------------

# Set target
SET(TARGET "solver")

# Add source files (may need to change according to directory organization)
FILE(GLOB_RECURSE TARGET_SRC "src/*.cc" "../shared/src/*.cc")
FILE(GLOB_RECURSE TARGET_INC "src/*.h"  "../shared/include/*.h")
INCLUDE_DIRECTORIES(include ../shared/include)

# ------------------------------------------------------------------------------
# Typically No Need To Modify Anything Beyond This Point
# ------------------------------------------------------------------------------

# Require cmake 2.8.8+
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

# When installing, include rpaths to avoid have user set LD_LIBRARY_PATH.
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Look for deal.II
FIND_PACKAGE(deal.II 8.4 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\nCould not find a suitably recent version of deal.II.\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake "
    "or set an environment variable \"DEAL_II_DIR\" that contains a path to a "
    "sufficiently recent version of deal.II."
    )
ENDIF()

MESSAGE(STATUS "Found deal.II version ${DEAL_II_PACKAGE_VERSION} at '${deal.II_DIR}'")

# Check trilinos and p4est are available
SET(_DEALII_GOOD ON)
IF(NOT DEAL_II_WITH_P4EST)
    MESSAGE(SEND_ERROR
      "\n-- deal.II was built without support for p4est!\n"
      )
    SET(_DEALII_GOOD OFF)
ENDIF()
IF(NOT DEAL_II_WITH_TRILINOS)
    MESSAGE(SEND_ERROR
      "\n-- deal.II was built without support for Trilinos!\n"
      )
    SET(_DEALII_GOOD OFF)
ENDIF()
IF (NOT _DEALII_GOOD)
  MESSAGE(FATAL_ERROR
    "\nGeocentric requires a deal.II installation built with support for Trilinos and p4est but one or both of these appears to be missing!\n"
    )
ENDIF()

# Invoke convenient autopilot from deal.ii
DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()

# Installation rule(s)
INSTALL_TARGETS(/bin ${TARGET})


