
cmake_minimum_required(VERSION 3.10)
project(SignalProcessing)

set(CMAKE_CXX_STANDARD 17)

# Find Python (for matplotlibcpp)
find_package(Python3 COMPONENTS Development REQUIRED)

# Include directories
include_directories(${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS} .)

# Sources
set(SOURCES examples.cpp sig_utils.cpp transfer.cpp)

# Executable
add_executable(signal_processing ${SOURCES})

# Link FFTW and Python
target_link_libraries(signal_processing fftw3 ${Python3_LIBRARIES})

