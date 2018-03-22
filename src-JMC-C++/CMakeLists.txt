cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project (general)

set(PROJECT_general_SRCS
  acploc.cpp
  bias.cpp
  comparscen.cpp
  conf.cpp
  data.cpp
  estimparam.cpp
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
  statdefs.cpp
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
  statdefs.hpp
  )


source_group("Include Files" FILES ${PROJECT_general_HEADERS})
add_executable(general general.cpp ${PROJECT_general_SRCS} ${PROJECT_general_HEADERS})
target_link_libraries(general dcmt)
target_compile_features(general PUBLIC cxx_std_17)

find_package(OpenMP)
if (OpenMP_FOUND)
#   add_library(OpenMP::OpenMP_BLOUBH SHARED IMPORTED)
   target_link_libraries(general OpenMP::OpenMP_CXX) 
#  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

install(TARGETS general RUNTIME DESTINATION .)
