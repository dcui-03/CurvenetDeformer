cmake_minimum_required(VERSION 3.20)
project(ProfileMover LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(PROFILEMOVER_WITH_POLYSCOPE "Build the Polyscope frontend executable" ON)

if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64)$" AND NOT DEFINED CMAKE_OSX_ARCHITECTURES)
  set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "Build architectures for macOS" FORCE)
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  execute_process(
    COMMAND uname -r
    OUTPUT_VARIABLE PROFILEMOVER_UNAME_R
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(PROFILEMOVER_UNAME_R MATCHES "WSL|Microsoft")
    set(PROFILEMOVER_IS_WSL TRUE)
    message(STATUS "Configuring ProfileMover for WSL. GUI execution requires a working X11/Wayland bridge.")
  endif()
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(nlohmann_json CONFIG QUIET)
if(NOT TARGET nlohmann_json::nlohmann_json)
  find_path(NLOHMANN_JSON_INCLUDE_DIR nlohmann/json.hpp)
  if(NOT NLOHMANN_JSON_INCLUDE_DIR)
    message(FATAL_ERROR
      "Could not find nlohmann/json.hpp. Install nlohmann_json or provide its include path.")
  endif()
  add_library(nlohmann_json::nlohmann_json INTERFACE IMPORTED)
  set_target_properties(nlohmann_json::nlohmann_json PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${NLOHMANN_JSON_INCLUDE_DIR}"
  )
endif()

include(cmake/libigl.cmake)

set(PROFILEMOVER_CORE_SOURCES
  src/IO/io.cpp
  src/curvenet/components/control.cpp
  src/curvenet/components/spline.cpp
  src/curvenet/curvenet.cpp
  src/dcurvenet/components/dvert.cpp
  src/dcurvenet/components/dsegment.cpp
  src/dcurvenet/components/dspline.cpp
  src/dcurvenet/dcurvenet.cpp
  src/dcurvenet/pdcurvenet.cpp
  src/mesh/cutmesh.cpp
  src/mesh/mesh.cpp
  src/mesh/polyHE/polyHE.cpp
  src/profilemover/profilemover.cpp
  src/utils/decUtils.cpp
  src/utils/jsonUtils.cpp
  src/utils/utils.cpp
)

add_library(profilemover_core ${PROFILEMOVER_CORE_SOURCES})
target_include_directories(profilemover_core PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
)
target_link_libraries(profilemover_core PUBLIC
  Eigen3::Eigen
  nlohmann_json::nlohmann_json
  igl::core
)

if(PROFILEMOVER_WITH_POLYSCOPE)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope/CMakeLists.txt")
    add_subdirectory(deps/polyscope)
  else()
    message(FATAL_ERROR
      "Polyscope dependency not found at deps/polyscope.\n"
      "Clone it with:\n"
      "  git clone --recursive https://github.com/nmwsharp/polyscope.git deps/polyscope")
  endif()

  add_executable(profile_mover src/main.cpp)
  target_link_libraries(profile_mover PRIVATE
    profilemover_core
    polyscope
  )
endif()
