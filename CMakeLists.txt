cmake_minimum_required(VERSION 3.12.0)

cmake_policy(SET CMP0048 NEW)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project (dcmt)

IF(MSVC AND MAKE_STATIC_EXE)
  add_compile_options(
    $<$<CONFIG:>:/MT> #---------|
    $<$<CONFIG:Debug>:/MTd> #---|-- Statically link the runtime libraries
    $<$<CONFIG:Release>:/MT> #--|
  )
endif()  


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

IF(MSVC AND MAKE_STATIC_EXE)
  add_compile_options(
    $<$<CONFIG:>:/MT> #---------|
    $<$<CONFIG:Debug>:/MTd> #---|-- Statically link the runtime libraries
    $<$<CONFIG:Release>:/MT> #--|
  )
endif()  


IF(MSVC)
  add_subdirectory (wingetopt)
  link_directories (${dcmt_BINARY_DIR} ${wingetopt_BINARY_DIR})
ELSE() 
  link_directories (${dcmt_BINARY_DIR})
ENDIF()

add_subdirectory (src-JMC-C++)

include(CTest)
enable_testing()

add_test(NAME functionaltest COMMAND "Rscript" "--vanilla tests.R $<TARGET_FILE:general>" WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/tests")