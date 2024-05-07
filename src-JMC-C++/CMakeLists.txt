cmake_minimum_required(VERSION 3.12.0)

project (general VERSION 1.0.4)

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
  bias.cpp
  conf.cpp
  comparscen.cpp
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
  bias.hpp
  conf.hpp
  comparscen.hpp
  )


source_group("Include Files" FILES ${PROJECT_general_HEADERS})
add_executable(general general.cpp ${PROJECT_general_SRCS} ${PROJECT_general_HEADERS})
target_link_libraries(general dcmt)

if (MSVC)
  target_link_libraries(general wingetopt)
endif()

if (MAKE_STATIC_EXE)
  if (MAKE_STATIC_EXE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif()
  if (UNIX AND NOT APPLE) 
    target_link_libraries(general iomp5 pthread dl c++ c++abi)
  endif() 
endif()

if (MSVC OR APPLE OR (NOT MAKE_STATIC_EXE))
  if (APPLE)
    set(OpenMP_ROOT "/usr/local/opt/libomp")
  endif()
  find_package(OpenMP REQUIRED)
  target_link_libraries(general OpenMP::OpenMP_CXX)
  include_directories(${OpenMP_CXX_INCLUDE_DIRS})
endif()

cmake_policy(SET CMP0087 NEW)

# install(CODE [[
#   file(TIMESTAMP $<TARGET_FILE:general> build_date "%Y%m%d-%H%M")
#   file(INSTALL 
#        DESTINATION ${CMAKE_BINARY_DIR}/install
#        FILES $<TARGET_FILE:general>
#        RENAME diyabc-${CMAKE_SYSTEM_NAME}-${build_date})
#   message("Build date : ${build_date}")
# ]])

install(TARGETS general RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/install)

if(MSVC)
  install(CODE [[
    message("Target to install : $<TARGET_FILE:general>")
    file(GET_RUNTIME_DEPENDENCIES
        RESOLVED_DEPENDENCIES_VAR _r_deps
        UNRESOLVED_DEPENDENCIES_VAR _u_deps
        PRE_INCLUDE_REGEXES "^.*[vV][cC]([oO][mM][pP]|[rR][uU][nN][tT][iI][mM][eE])140.*$"
        PRE_EXCLUDE_REGEXES "^.*$"
        EXECUTABLES  $<TARGET_FILE:general>)
    foreach(_file ${_r_deps})
      file(INSTALL
        DESTINATION ${CMAKE_BINARY_DIR}/install
        FILES "${_file}"
      )
    endforeach()
    list(LENGTH _u_deps _u_length)
    if("${_u_length}" GREATER 0)
      message(WARNING "Unresolved dependencies detected!")
    endif()
        ]]
  )
  install(FILES ${resolved_deps}
          DESTINATION ${CMAKE_BINARY_DIR}/install)

endif()
