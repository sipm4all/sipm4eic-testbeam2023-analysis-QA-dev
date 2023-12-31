# CMakeLists.txt

cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(hough CUDA CXX)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Set CUDA architecture (change according to your GPU architecture)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} --gpu-architecture=sm_75)

# Find CUDA toolkit
find_package(CUDA REQUIRED)
find_package(Boost COMPONENTS program_options system REQUIRED)
find_package(ROOT COMPONENTS RIO REQUIRED)

# Add CUDA include directories
include_directories(${CUDA_INCLUDE_DIRS} ${ROOT_INCLUDE_DIR})

# Add your source files (main.cpp and cuda.cu)
set(SOURCE_FILES
    main.cc
    cuda.cu
)

# Add executable and specify source files
cuda_add_executable(hough ${SOURCE_FILES})
set_target_properties(hough PROPERTIES CUDA_ARCHITECTURES "75")
target_link_libraries(hough ${CUDA_LIBRARIES} ${Boost_LIBRARIES} ${ROOT_LIBS})
install(TARGETS hough RUNTIME DESTINATION bin)
