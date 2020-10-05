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

if (MAKE_STATIC_EXE)
  if (APPLE AND MAKE_STATIC_EXE)
    target_link_libraries(general /usr/local/lib/libomp.a)
  endif()

  if (UNIX AND NOT APPLE) 
    target_link_libraries(general iomp5 pthread dl c++ c++abi)
  endif() 
endif()

if (MSVC OR (NOT MAKE_STATIC_EXE))
  find_package(OpenMP REQUIRED)
  target_link_libraries(general OpenMP::OpenMP_CXX)
endif()