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

set (PROJECT_initRNGproject_SRCS
  initRNGproject.cpp
  mesutils.cpp
  matrices.cpp
  randomgenerator.cpp)

add_compile_options(-fopenmp -g -O3)
add_executable(general ${PROJECT_general_SRCS})
target_link_libraries(general -ldcmt -lgomp)

add_executable(initRNGproject ${PROJECT_initRNGproject_SRCS})
target_link_libraries(initRNGproject -ldcmt -lgomp)
