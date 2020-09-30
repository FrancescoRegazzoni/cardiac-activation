cmake_minimum_required( VERSION 3.10 )

project( cardiac-activation )

if( NOT CMAKE_BUILD_TYPE )
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package( Boost 1.41 REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )

find_package( Eigen3 3.3 REQUIRED NO_MODULE )

add_executable( run_model main.cpp sarcomere.cpp model_RDQ18.cpp model_RDQ20_MF.cpp model_RDQ20_SE.cpp )

target_link_libraries( run_model Eigen3::Eigen )