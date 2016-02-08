cmake_minimum_required(VERSION 2.8.6)

project (general)

set (PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

set(PROJECT_general_SRCS
  acploc.cpp
  bias.cpp
  comparscen.cpp
  conf.cpp
  data.cpp
  estimparam.cpp
  general.cpp
  header.cpp
  history.cpp
  matrices.cpp
  mesutils.cpp
  particleset.cpp
  particuleC.cpp
  randforest.cpp
  randomgenerator.cpp
  reftable.cpp
  simfile.cpp
  sumstat.cpp
  modchec.cpp
)

set(PROJECT_general_HEADERS
  acploc.hpp
  bias.hpp
  comparscen.hpp
  conf.hpp
  data.hpp
  estimparam.hpp
  global.hpp
  header.hpp
  history.hpp
  matrices.hpp
  mesutils.hpp
  modchec.hpp
  particleset.hpp
  particuleC.hpp
  randforest.hpp
  randomgenerator.hpp
  reftable.hpp
  simfile.hpp
)
source_group("Include Files" FILES ${PROJECT_general_HEADERS})

set (PROJECT_initRNGproject_SRCS
  initRNGproject.cpp
  mesutils.cpp
  matrices.cpp
  randomgenerator.cpp)

add_executable(general ${PROJECT_general_SRCS} ${PROJECT_general_HEADERS})
add_executable(initRNGproject ${PROJECT_initRNGproject_SRCS})

set_directory_properties(PROPERTIES COMPILE_OPTIONS "-g;-O3")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

target_link_libraries(general dcmt)
target_link_libraries(initRNGproject dcmt)    

if (MSVC)
  target_link_libraries(general wingetopt)
  target_link_libraries(initRNGproject wingetopt)
endif()


