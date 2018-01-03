cmake_minimum_required (VERSION 2.6)
project (PhyloNaut)

file(GLOB SOURCES *.cpp TripartitionScorer/*.cpp CladeExtractor/*.cpp Analysis/*.cpp)

add_definitions(-std=c++14 -g -O3 -Wall)

option(ENABLE_PROFILING
  "Enable google profiler"
  OFF)




if(ENABLE_PROFILING)
  add_definitions(-DENABLE_PROFILING)

  SET(Boost_USE_STATIC_LIBS OFF)
  SET(Boost_USE_MULTITHREADED ON)
  SET(Boost_USE_STATIC_RUNTIME OFF)

else()
  ADD_LIBRARY(phylonaut_static STATIC ${SOURCES})
  SET_TARGET_PROPERTIES(phylonaut_static PROPERTIES OUTPUT_NAME phylonaut)  	
  install(TARGETS phylonaut_static
		DESTINATION lib)
  TARGET_LINK_LIBRARIES(phylonaut_static phylokit)


endif(ENABLE_PROFILING)


ADD_LIBRARY(phylonaut_shared SHARED ${SOURCES} )
SET_TARGET_PROPERTIES(phylonaut_shared PROPERTIES OUTPUT_NAME phylonaut) 
TARGET_LINK_LIBRARIES(phylonaut_shared phylokit)
install(TARGETS phylonaut_shared
		DESTINATION lib)


file(GLOB HEADERS *.hpp TripartitionScorer/*.hpp Analysis/*.hpp CladeExtractor/*.hpp)

install(FILES ${HEADERS} DESTINATION include/phylonaut)
