cmake_minimum_required(VERSION 3.15)

project(
    RadialModel
    VERSION 0.0.1
    LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)    

# add in NumericConcepts
include(FetchContent)
FetchContent_Declare(
  NumericConcepts
  GIT_REPOSITORY https://github.com/da380/NumericConcepts.git
  GIT_TAG main
)
FetchContent_MakeAvailable(NumericConcepts)


# add in NumericConcepts
include(FetchContent)
FetchContent_Declare(
  Dimensions
  GIT_REPOSITORY https://github.com/da380/Dimensions.git
  GIT_TAG main
)
FetchContent_MakeAvailable(Dimensions)

# add in Eigen3
FetchContent_Declare(
  Eigen3
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG master
)
FetchContent_MakeAvailable(Eigen3)

# add in Interpolation library. 
FetchContent_Declare(
  Interpolation
  GIT_REPOSITORY git@github.com:da380/Interpolation.git
  GIT_TAG newBuildInterface
)
FetchContent_MakeAvailable(Interpolation)


add_library(${PROJECT_NAME} INTERFACE)
target_link_libraries(${PROJECT_NAME} INTERFACE NumericConcepts Dimensions Interpolation)


if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    add_subdirectory(examples)
endif()

include(GNUInstallDirs)

target_include_directories(
    ${PROJECT_NAME}
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
              $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)








