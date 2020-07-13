cmake_minimum_required (VERSION 2.6)
project (PhyloNaut)
OPTION(BUILD_STATIC "build static (otherwise build shared)" ON)
file(GLOB SOURCES *.cpp TripartitionScorer/*.cpp CladeExtractor/*.cpp Analysis/*.cpp)

if(WIN32)
add_definitions(/std:c++14)
INCLUDE_DIRECTORIES("C:/Users/pranjal/src/phylokit/install/x64-Release/include")
else()
add_definitions(-std=c++14 -g -O3 -Wall)
include_directories("../../install/include")
link_directories("../../install/lib")
endif()

option(ENABLE_PROFILING
  "Enable google profiler"
  OFF)




if(ENABLE_PROFILING)
  add_definitions(-DENABLE_PROFILING)

  SET(Boost_USE_STATIC_LIBS OFF)
  SET(Boost_USE_MULTITHREADED ON)
  SET(Boost_USE_STATIC_RUNTIME OFF)

else()
endif(ENABLE_PROFILING)

IF (BUILD_STATIC)
  ADD_LIBRARY(phylonaut_static STATIC ${SOURCES})
  SET_TARGET_PROPERTIES(phylonaut_static PROPERTIES OUTPUT_NAME phylonaut)  	
  install(TARGETS phylonaut_static
		DESTINATION lib)
  TARGET_LINK_LIBRARIES(phylonaut_static phylokit)
  ENDIF()
INCLUDE_DIRECTORIES("boost/" )

if(WIN32)
LINK_DIRECTORIES("C:/Users/pranjal/src/phylokit/install/x64-Release/lib/")

endif()

if(NOT(BUILD_STATIC))
ADD_LIBRARY(phylonaut_shared SHARED ${SOURCES} )

SET_TARGET_PROPERTIES(phylonaut_shared PROPERTIES OUTPUT_NAME phylonaut) 
TARGET_LINK_LIBRARIES(phylonaut_shared phylokit)
install(TARGETS phylonaut_shared
		DESTINATION lib)
		endif()


file(GLOB HEADERS *.hpp TripartitionScorer/*.hpp Analysis/*.hpp CladeExtractor/*.hpp)

install(FILES ${HEADERS} DESTINATION include/phylonaut)
