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
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
###############################################################

find_package(ISMRMRD REQUIRED)
find_package(Gadgetron REQUIRED)
find_package(Boost REQUIRED)
find_package(ACE REQUIRED)

set(CMAKE_INSTALL_PREFIX ${GADGETRON_HOME})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -fopenmp -O3 -Wunused-variable -Wunused-but-set-variable -Wunused-value -Wfloat-equal -Wuseless-cast") #-O3

INCLUDE_DIRECTORIES(${ACE_INCLUDE_DIR} 
     ${Boost_INCLUDE_DIR}
     ${GADGETRON_INCLUDE_DIR}
     ${ISMRMRD_INCLUDE_DIR}
     ${ISMRMRD_SCHEMA_DIR}
     ${ISMRMRD_XSD_INCLUDE_DIR}
     ${ISMRMRD_REAL_INCLUDE}
     )

LINK_DIRECTORIES(${GADGETRON_LIB_DIR} )

ADD_LIBRARY(gadgetron_PUROR SHARED	IEVChannelSumGadget.cpp
					ContrastGadget.cpp
					Mask.cpp
					UnwrapGadget.cpp
					GetMaskGadget.cpp
					PUROR.cpp
					ReorderGadget.cpp
					AddMag.cpp
					AddMetaData.cpp
					HPFGadget.cpp
					PhaseMaskGadget.cpp
					SWIGadget.cpp
					DecImageIndex.cpp
					FlipGadget.cpp
					)

TARGET_LINK_LIBRARIES(gadgetron_PUROR #${MPI_LIBRARIES}
                      optimized ${ACE_LIBRARIES} 
                      debug ${ACE_DEBUG_LIBRARY}
			    gadgetron_toolbox_cpufft)


INSTALL (FILES 				PUROR.h
					gadgetron_PUROR_export.h
					Mask.h
					GetMaskGadget.h
					ContrastGadget.h
					UnwrapGadget.h
					IEVChannelSumGadget.h
					ReorderGadget.h
					HPFGadget.h
					PhaseMaskGadget.h
					SWIGadget.h
					AddMag.h
					AddMetaData.h
					DecImageIndex.h
					FlipGadget.h
					DESTINATION include)

INSTALL(TARGETS gadgetron_PUROR DESTINATION lib)


