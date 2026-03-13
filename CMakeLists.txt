cmake_minimum_required(VERSION 3.20)
project(ProfileFormer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(PROFILEFORMER_WITH_POLYSCOPE "Build Polyscope frontend executable" ON)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(nlohmann_json CONFIG QUIET)
if(NOT TARGET nlohmann_json::nlohmann_json)
  find_path(NLOHMANN_JSON_INCLUDE_DIR nlohmann/json.hpp)
  if(NOT NLOHMANN_JSON_INCLUDE_DIR)
    message(FATAL_ERROR "Could not find nlohmann/json.hpp. Install nlohmann_json (apt/conda) or set include path.")
  endif()
  add_library(nlohmann_json::nlohmann_json INTERFACE IMPORTED)
  set_target_properties(nlohmann_json::nlohmann_json PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${NLOHMANN_JSON_INCLUDE_DIR}")
endif()

include(cmake/libigl.cmake)

add_library(profileformer_core
  src/curvenet/components/control.cpp
  src/curvenet/components/spline.cpp
  src/curvenet/curvenet.cpp
  src/dcurvenet/components/dvert.cpp
  src/dcurvenet/components/dsegment.cpp
  src/dcurvenet/components/dspline.cpp
  src/dcurvenet/dcurvenet.cpp
  src/dcurvenet/pdcurvenet.cpp
  src/profileformer/profileformer.cpp
  src/utils/jsonUtils.cpp
  src/utils/utils.cpp
)

target_include_directories(profileformer_core PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

target_link_libraries(profileformer_core PUBLIC
  Eigen3::Eigen
  nlohmann_json::nlohmann_json
  igl::core
)

if(PROFILEFORMER_WITH_POLYSCOPE)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope/CMakeLists.txt")
    add_subdirectory(deps/polyscope)
  else()
    message(FATAL_ERROR
      "Polyscope dependency not found at deps/polyscope.\n"
      "For WSL/Linux, clone it:\n"
      "  git clone --recursive https://github.com/nmwsharp/polyscope.git deps/polyscope")
  endif()

  add_executable(profile_former src/main.cpp)
  target_link_libraries(profile_former PRIVATE
    profileformer_core
    polyscope
  )
endif()
list(REMOVE_ITEM PROJECT_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")

add_library(pm_core ${PROJECT_SOURCES})
target_include_directories(pm_core PUBLIC
  #"${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope/deps/args"
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

# Link igl (and the glfw module) to your project
target_link_libraries(pm_core PUBLIC 
  igl::glfw
  ## Other modules you could link to
  # igl::embree
  # igl::imgui
  # igl::opengl
  # igl::stb
  # igl::predicates
  # igl::xml
  # igl_copyleft::cgal
  # igl_copyleft::comiso
  # igl_copyleft::core
  # igl_copyleft::cork
  #igl_copyleft::tetgen
  # igl_restricted::matlab
  # igl_restricted::mosek
  #igl_restricted::triangle
  polyscope
  #geometry-central
)

# Add your project files
add_executable(profile_mover src/main.cpp)

# 
target_link_libraries(profile_mover PRIVATE pm_core)
