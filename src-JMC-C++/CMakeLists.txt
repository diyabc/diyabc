cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project (general)
enable_testing()

set (PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

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
  )
source_group("Include Files" FILES ${PROJECT_general_HEADERS})

set (PROJECT_initRNGproject_SRCS
  initRNGproject.cpp
  mesutils.cpp
  matrices.cpp
  randomgenerator.cpp)

add_executable(general general.cpp ${PROJECT_general_SRCS} ${PROJECT_general_HEADERS})
add_executable(general-test general-test.cpp ${PROJECT_general_SRCS} ${PROJECT_general_HEADERS})
add_executable(initRNGproject ${PROJECT_initRNGproject_SRCS})

find_package(Boost COMPONENTS filesystem REQUIRED)

add_test(NAME f3reich WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test/f3reich/pool COMMAND general-test)

if (NOT (CMAKE_CXX_COMPILER_ID MATCHES "Intel"))
  target_compile_features(general PRIVATE cxx_right_angle_brackets)
  target_compile_features(general-test PRIVATE cxx_right_angle_brackets)
  target_compile_features(initRNGproject PRIVATE cxx_right_angle_brackets)
endif()

set_directory_properties(PROPERTIES COMPILE_OPTIONS "-g;-O3")

target_link_libraries(general dcmt)
target_link_libraries(general-test dcmt ${Boost_LIBRARIES})
target_link_libraries(initRNGproject dcmt)    

if (MSVC)
  target_link_libraries(general wingetopt)
  target_link_libraries(general-test wingetopt)
  target_link_libraries(initRNGproject wingetopt)
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
  if ((CMAKE_CXX_COMPILER_ID MATCHES "Clang") AND (CMAKE_SYSTEM_NAME MATCHES "Windows"))
    set (OpenMP_C_FLAGS "${OpenMP_C_FLAGS} -fnoopenmp-use-tls")
    set (OpenMP_CXX_FLAGS "${OpenMP_CXX_FLAGS} -fnoopenmp-use-tls")
    target_link_libraries(general -lomp)
    target_link_libraries(general-test -lomp)
    target_link_libraries(initRNGproject -lomp) 
  endif()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
