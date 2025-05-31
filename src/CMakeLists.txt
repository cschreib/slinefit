cmake_minimum_required(VERSION 2.6...4.0)
project(slinefit C CXX)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release (-DCMAKE_BUILD_TYPE=...)")
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(NO_REFLECTION ON)
set(NO_FFTW ON)
set(NO_GSL ON)
set(NO_LAPACK ON)
set(NO_WCSLIB ON)

find_package(vif REQUIRED)

# Include dependencies
message(STATUS ${VIF_INCLUDE_DIRS})
include_directories(${VIF_INCLUDE_DIRS})

# Build slinefit
add_executable(slinefit slinefit.cpp)
target_link_libraries(slinefit ${VIF_LIBRARIES})
install(TARGETS slinefit DESTINATION bin)
