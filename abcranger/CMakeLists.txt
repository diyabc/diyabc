project (abcranger)
set(RANGER_SOURCE_DIR ../ranger/src)
include_directories("${RANGER_SOURCE_DIR}")
set(PROJECT_abcranger_SRCS
  ${RANGER_SOURCE_DIR}/Data.cpp
  ${RANGER_SOURCE_DIR}/DataChar.cpp
  ${RANGER_SOURCE_DIR}/DataDouble.cpp
  ${RANGER_SOURCE_DIR}/DataFloat.cpp
  ${RANGER_SOURCE_DIR}/Forest.cpp
  ${RANGER_SOURCE_DIR}/ForestClassification.cpp
  ${RANGER_SOURCE_DIR}/ForestRegression.cpp
  ${RANGER_SOURCE_DIR}/Tree.cpp
  ${RANGER_SOURCE_DIR}/TreeClassification.cpp
  ${RANGER_SOURCE_DIR}/TreeRegression.cpp
  ${RANGER_SOURCE_DIR}/utility.cpp
  )
  
  add_executable(abcranger mainrangertest.cpp abcranger.cpp ${PROJECT_abcranger_SRCS})
  target_compile_features(general PRIVATE cxx_right_angle_brackets cxx_lambdas)