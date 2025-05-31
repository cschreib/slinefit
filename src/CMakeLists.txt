cmake_minimum_required(VERSION 3.8...4.0)
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

# Build slinefit
add_executable(slinefit slinefit.cpp)
target_include_directories(slinefit PRIVATE ${VIF_INCLUDE_DIRS})
target_link_libraries(slinefit PRIVATE ${VIF_LIBRARIES})
target_compile_features(slinefit PRIVATE cxx_std_11)
install(TARGETS slinefit DESTINATION bin)
