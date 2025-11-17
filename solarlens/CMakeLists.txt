cmake_minimum_required(VERSION 3.16)
project(SolarLens VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Space-grade compile flags
add_compile_options(
    -Wall -Wextra -Werror -Wpedantic
    -fno-exceptions -fno-rtti
    -O2 -march=native
    -D_GNU_SOURCE
)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Source files
set(SOURCES
    src/main.cpp
    src/gravitational_lens.cpp
    src/exoplanet_detector.cpp
    src/formation_control.cpp
    src/navigation.cpp
    src/deep_space_relay.cpp
    src/power_thermal.cpp
)

# Executable
add_executable(solarlens ${SOURCES})

# Math library
target_link_libraries(solarlens m)

# Installation
install(TARGETS solarlens DESTINATION bin)

# Print build information
message(STATUS "SolarLens Mission Control System")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
