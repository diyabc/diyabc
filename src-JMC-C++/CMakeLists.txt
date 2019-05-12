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
  randomgenerator.hpp
  reftable.hpp
  simfile.hpp
  statdefs.hpp
  )


IF(MSVC)
  add_compile_options(/permissive-)
  # https://stackoverflow.com/questions/14172856/compile-with-mt-instead-of-md-using-cmake
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
#else()
#  link_libraries(c++ c++abi)
endif()

source_group("Include Files" FILES ${PROJECT_general_HEADERS})
add_executable(general general.cpp ${PROJECT_general_SRCS} ${PROJECT_general_HEADERS})
target_link_libraries(general dcmt)
target_compile_features(general PUBLIC cxx_std_17)

if (MSVC)
  target_link_libraries(general wingetopt)
endif()

find_package(Threads)
target_link_libraries(general Threads::Threads)

find_package(OpenMP)
if (OpenMP_FOUND)
   target_link_libraries(general OpenMP::OpenMP_CXX) 
endif()

#install(TARGETS general RUNTIME DESTINATION .)
