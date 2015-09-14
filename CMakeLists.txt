cmake_minimum_required(VERSION 2.6)

project(PUROR_LIB)


###############################################################
#Bootstrap search for libraries 
# (We need to find cmake modules in Gadgetron)
###############################################################
find_path(GADGETRON_CMAKE_MODULES FindGadgetron.cmake HINTS
$ENV{GADGETRON_HOME}/cmake
/usr/local/gadgetron)

if (NOT GADGETRON_CMAKE_MODULES)
  MESSAGE(FATAL_ERROR "GADGETRON_CMAKE_MODULES cannot be found. 
   Try to set GADGETRON_HOME environment variable.")
endif(NOT GADGETRON_CMAKE_MODULES)

set(CMAKE_MODULE_PATH ${GADGETRON_CMAKE_MODULES})
###############################################################

find_package(Ismrmrd REQUIRED)
find_package(Gadgetron REQUIRED)
find_package(Boost REQUIRED)
find_package(ACE REQUIRED)
#find_package(MPI REQUIRED)

#include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_INSTALL_PREFIX ${GADGETRON_HOME})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -fopenmp -O3") #-O3
set(CMAKE_CXX_COMPILER "/opt/rh/devtoolset-3/root/usr/bin/c++")
#set(CMAKE_CXX_COMPILER "/usr/lib64/openmpi/bin/mpicxx")
#set(CMAKE_CXX_COMPILER "/usr/lib64/openmpi/bin/mpicc")
set(CMAKE_C_COMPILER "/opt/rh/devtoolset-3/root/usr/bin/cc")

#set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
#set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})



include_directories(include " ~/home/tim/ismrmrd/")
INCLUDE_DIRECTORIES(${ACE_INCLUDE_DIR} 
     ${Boost_INCLUDE_DIR}
     ${GADGETRON_INCLUDE_DIR}
     ${ISMRMRD_INCLUDE_DIR}
     ${ISMRMRD_SCHEMA_DIR}
     ${ISMRMRD_XSD_INCLUDE_DIR}
     ${ISMRMRD_REAL_INCLUDE}
     )

set(SOURCES "${SOURCES} ${ISMRMRD_REAL_INCLUDE}/ismrmrd.h")


LINK_DIRECTORIES(${GADGETRON_LIB_DIR} )

ADD_LIBRARY(gadgetron_PUROR SHARED AccGadget.cpp ContrastGadget.cpp Mask.cpp UnwrapGadget.cpp GetMaskGadget.cpp PUROR.cpp Unwrap3DGadget.cpp HPFGadget.cpp ExtractMagGadget.cpp)

TARGET_LINK_LIBRARIES(gadgetron_PUROR #${MPI_LIBRARIES}
                      optimized ${ACE_LIBRARIES} 
                      debug ${ACE_DEBUG_LIBRARY}
			    gadgetron_toolbox_cpufft)


INSTALL (FILES PUROR.h Mask.h GetMaskGadget.h ContrastGadget.h UnwrapGadget.h AccGadget.h Unwrap3DGadget.h HPFGadget.h ExtractMagGadget.h DESTINATION include)

INSTALL(TARGETS gadgetron_PUROR DESTINATION lib)


