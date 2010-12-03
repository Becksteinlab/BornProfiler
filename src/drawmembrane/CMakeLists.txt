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
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      )
endif (NOT CMAKE_BUILD_TYPE)

# default installation in user's HOME directory
# (although ccmake still insists on /usr/local)
set (CMAKE_INSTALL_PREFIX $ENV{HOME} CACHE STRING
      "Choose the installation directory.")

# build executables
add_executable (draw_membrane2a draw_membrane2a.c)

# linking
# (need libm and libz for new 2a)
target_link_libraries (draw_membrane2a m z)

# install executables
install (TARGETS draw_membrane2a  RUNTIME DESTINATION bin)

