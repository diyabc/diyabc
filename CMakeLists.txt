cmake_minimum_required(VERSION 2.8.6)

project (diyabc)


add_subdirectory (dcmt0.6.1)

IF(MSVC)
  add_subdirectory (wingetopt)
  link_directories (${dcmt_BINARY_DIR} ${wingetopt_BINARY_DIR})
  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Za")
ELSE() 
  link_directories (${dcmt_BINARY_DIR})
ENDIF()

add_subdirectory (src-JMC-C++)







