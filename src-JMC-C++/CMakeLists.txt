cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project (general)

set(PROJECT_general_SRCS
  acploc.cpp
  data.cpp
  header.cpp
  history.cpp
  matrices.cpp
  mesutils.cpp
  particleset.cpp
  particuleC.cpp
  randomgenerator.cpp
  reftable.cpp
  simfile.cpp
  statdefs.cpp
  sumstat.cpp
  modchec.cpp
  estimparam.cpp
  )

set(PROJECT_general_HEADERS
  acploc.hpp
  data.hpp
  header.hpp
  history.hpp
  matrices.hpp
  mesutils.hpp
  particleset.hpp
  particuleC.hpp
  randomgenerator.hpp
  reftable.hpp
  simfile.hpp
  statdefs.hpp
  modchec.hpp
  estimparam.hpp
  )

source_group("Include Files" FILES ${PROJECT_general_HEADERS})
add_executable(general general.cpp ${PROJECT_general_SRCS} ${PROJECT_general_HEADERS})
target_link_libraries(general dcmt)

if (MSVC)
  target_link_libraries(general wingetopt)
endif()

find_package(OpenMP)

if(${OpenMP_Found})
  target_link_libraries(general OpenMP::OpenMP_CXX OpenMP::OpenMP_C)
endif()