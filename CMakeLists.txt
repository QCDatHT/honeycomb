# ============================
# BASIC DEFINITIONS
# ============================

# Define minimum version of cmake
cmake_minimum_required(VERSION 3.10)

# Define project name and its language
project(Honeycomb C)

# Define code version
set(Honeycomb_VERSION 1.0.1)

# Define c++ standard and issue all the warning demanded by this standard
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
add_compile_options(-Wall -Wextra -Werror -fPIC -O3)
add_compile_options(
   -march=native
   -fno-math-errno
   -fno-trapping-math
   -fno-rounding-math
   -fexcess-precision=fast
   -fcx-limited-range
   -fomit-frame-pointer
   -funroll-loops
)

# -fopenmp -DUSE_OMP

# Check the processor architecture
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|i[3-6]86|amd64")
   # x86-specific compiler flags
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mfpmath=sse -msse2")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
   # ARM-specific compiler flags
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mcpu=native -mtune=native")
endif()

if(NOT DEFINED CMAKE_MACOSX_RPATH)
   set(CMAKE_MACOSX_RPATH 0)
endif()

if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.27.0")
   cmake_policy(SET CMP0148 OLD)
endif()

# ============================
# MAIN
# ============================

# Define include directory
include_directories(${PROJECT_SOURCE_DIR}/inc)

# Generate list of source files
file(GLOB_RECURSE source_files src/*)

# Define target library
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
add_library(Honeycomb SHARED ${source_files})

# Define libraries to be linked
target_link_libraries(Honeycomb)

# ============================
# Enable testing and build test codes
add_subdirectory(tests)

# ============================
# CONFIGURATION
# ============================

# Configure version to file
configure_file("${PROJECT_SOURCE_DIR}/inc/Honeycomb/version.h.in" "${PROJECT_SOURCE_DIR}/inc/Honeycomb/version.h")

# Configure configuration script
set(prefix ${CMAKE_INSTALL_PREFIX})
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")
set(cflags "-I${prefix}/include")
set(ldflags "-L${prefix}/lib -lHoneycomb")

configure_file("${PROJECT_SOURCE_DIR}/bin/Honeycomb-config.in" "${PROJECT_BINARY_DIR}/bin/Honeycomb-config")

# ============================
# INSTALLATION
# ============================

# Installation of code
install(TARGETS Honeycomb RUNTIME DESTINATION bin LIBRARY DESTINATION lib)
install(DIRECTORY inc/honeycomb DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/bin/Honeycomb-config DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
