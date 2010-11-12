# CMakeLists.txt: Written by Oliver Beckstein <orbeckst@gmail.com> (2010)
# Placed in the Public Domain (2010)
#
# CMake project file for 
# (1) draw_membrane2.c (c) 2008 Michael Grabe
# Downloaded from http://www.poissonboltzmann.org/apbs/examples/potentials-of-mean-force/the-polar-solvation-potential-of-mean-force-for-a-helix-in-a-dielectric-slab-membrane
# (2) draw_membrane4.c (c) 2010 Michael Grabe (from apbsmem-1.04.zip)
# (3) draw_membrane2a.c is draw_membrane2.c with modifications by
#     Oliver Beckstein (2010)

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
add_executable (draw_membrane2 draw_membrane2.c)
add_executable (draw_membrane2a draw_membrane2a.c)
add_executable (draw_membrane4 draw_membrane4.c)

# linking
# (need libm and libz for new 2a)
target_link_libraries (draw_membrane2 m)
target_link_libraries (draw_membrane2a m z)
target_link_libraries (draw_membrane4 m)

# install executables
#install (TARGETS draw_membrane2 draw_membrane2a draw_membrane4 RUNTIME DESTINATION bin)
install (TARGETS draw_membrane2a  RUNTIME DESTINATION bin)

