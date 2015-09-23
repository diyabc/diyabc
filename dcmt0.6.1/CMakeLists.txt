cmake_minimum_required(VERSION 2.8.6)

project (dcmt)

include_directories("${PROJECT_SOURCE_DIR}/include")

set(PROJECT_dcmt_SRCS
  ${PROJECT_SOURCE_DIR}/lib/init.c
  ${PROJECT_SOURCE_DIR}/lib/check32.c
  ${PROJECT_SOURCE_DIR}/lib/prescr.c
  ${PROJECT_SOURCE_DIR}/lib/mt19937.c
  ${PROJECT_SOURCE_DIR}/lib/eqdeg.c
  ${PROJECT_SOURCE_DIR}/lib/seive.c
  ${PROJECT_SOURCE_DIR}/lib/genmtrand.c
  )

add_library(dcmt ${PROJECT_dcmt_SRCS})




