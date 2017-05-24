cmake_minimum_required(VERSION 2.8.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project (dcmt)
set(DCMT_SOURCE_DIR dcmt)
include_directories("${DCMT_SOURCE_DIR}/include")

set(PROJECT_dcmt_SRCS
  ${DCMT_SOURCE_DIR}/lib/init.c
  ${DCMT_SOURCE_DIR}/lib/check32.c
  ${DCMT_SOURCE_DIR}/lib/prescr.c
  ${DCMT_SOURCE_DIR}/lib/mt19937.c
  ${DCMT_SOURCE_DIR}/lib/eqdeg.c
  ${DCMT_SOURCE_DIR}/lib/seive.c
  ${DCMT_SOURCE_DIR}/lib/genmtrand.c
  )

add_library(dcmt ${PROJECT_dcmt_SRCS})

project (diyabc)
enable_testing()

IF(MSVC)
  add_subdirectory (wingetopt)
  link_directories (${dcmt_BINARY_DIR} ${wingetopt_BINARY_DIR})
  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Za")
ELSE() 
  link_directories (${dcmt_BINARY_DIR})
ENDIF()

project (dcmt)


add_subdirectory (src-JMC-C++)






