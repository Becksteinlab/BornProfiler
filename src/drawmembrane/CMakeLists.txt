# CMakeLists.txt: Written by Oliver Beckstein <orbeckst@gmail.com> (2010)
# Placed in the Public Domain (2010)
#
# CMake project file for draw_membrane2a.c [Oliver Beckstein (2010)]
# draw_membrane2a.c  is based on draw_membrane2.c [Michael Grabe (2008)]
# from http://www.poissonboltzmann.org/apbs/examples/potentials-of-mean-force/the-polar-solvation-potential-of-mean-force-for-a-helix-in-a-dielectric-slab-membrane
# and some code from draw_membrane4.c (from apbsmem-1.04.zip) [Michael Grabe (2010)]

cmake_minimum_required (VERSION 2.4)
project (draw_membrane)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation in user's HOME directory
if (NOT HAVE_CMAKE_INSTALL_PREFIX)
   set(CMAKE_INSTALL_PREFIX $ENV{HOME} CACHE PATH
       "Choose the installation directory." FORCE)
   set(HAVE_CMAKE_INSTALL_PREFIX 1 CACHE INTERNAL
       "hack so that users can set PREFIX but start with our default" FORCE)			      
endif (NOT HAVE_CMAKE_INSTALL_PREFIX)

# find zlib in the system
include(${CMAKE_ROOT}/Modules/FindZLIB.cmake)

if (NOT ZLIB_FOUND)
   if (CMAKE_VERBOSE_MAKEFILE)
      message("**********************************************")
      message("  Necessary zlib library not found")
      message("  Install zlib through your package manager")
      message("  and look for zlib-devel or similar")
      message("**********************************************")
   endif (CMAKE_VERBOSE_MAKEFILE)
endif (NOT ZLIB_FOUND)

include_directories(${ZLIB_INCLUDE_DIRS})
list(APPEND LIBS ${ZLIB_LIBRARIES})

# new draw_membrane2a needs both libz and the maths lib
# (math should come at the end)
find_library(M_LIB m)
list(APPEND LIBS ${M_LIB})
mark_as_advanced(M_LIB)

# build executables
add_executable(draw_membrane2a draw_membrane2a.c)

# linking
# (need libm and libz for new 2a)
target_link_libraries (draw_membrane2a ${LIBS})

# install executables
install (TARGETS draw_membrane2a  RUNTIME DESTINATION bin)

