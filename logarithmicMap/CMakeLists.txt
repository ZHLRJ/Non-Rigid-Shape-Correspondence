cmake_minimum_required(VERSION 3.10.0)

project(logarithmicMap)

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Print the build type
# set(CMAKE_BUILD_TYPE "Debug")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
endif()
message(STATUS "cmake build type: ${CMAKE_BUILD_TYPE}")

### Configure the compiler
# This is a basic, decent setup that should do something sane on most compilers

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  SET(BASE_CXX_FLAGS "-std=c++11 -Wall -Wextra")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces -Wno-unused-private-field")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Setting clang-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=3 -fcolor-diagnostics")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g3 -fsanitize=address -fno-limit-debug-info")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("Setting gcc-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g3")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
  endif()
  
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  
  # https://gitlab.kitware.com/cmake/cmake/-/issues/22484
  include(CheckCXXCompilerFlag)
  unset(COMPILER_SUPPORTS_MARCH_NATIVE CACHE)
  CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
  if(COMPILER_SUPPORTS_MARCH_NATIVE)
    SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native -DNDEBUG")
  else()
    SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -mcpu=apple-m1 -DNDEBUG") # Apple M1
  endif()
else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()

# Edit files to build here.
set(SRCS 
  src/main.cpp
  src/heat-method.cpp
  ../../utils/src/colormap.cpp
  ../../core/src/geometry.cpp
)


add_subdirectory(../../extern/geometry-central deps/geometry-central)
add_subdirectory(../../extern/polyscope deps/polyscope)
add_subdirectory(../../extern/googletest deps/googletest)

# Change the executable name by changing "main" in the lines below.
add_executable(p10-main "${SRCS}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main)
target_include_directories(p10-main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(p10-main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# add the args.hxx project which we use for command line args
target_include_directories(p10-main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../extern/polyscope/deps/args")
target_link_libraries(p10-main geometry-central polyscope)

# Logarithmic Map  checkerMap
set(checkerMap_SRCS
        src/checkerMap.cpp
        ../../core/src/geometry.cpp
        )
add_executable(checkerMap "${checkerMap_SRCS}")
target_include_directories(checkerMap PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_link_libraries(checkerMap geometry-central polyscope)


set(KarcherMean_SRCS
        src/KarcherMean.cpp
        ../../core/src/geometry.cpp
        )
add_executable(KarcherMean "${KarcherMean_SRCS}")
target_include_directories(KarcherMean PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_link_libraries(KarcherMean geometry-central polyscope)


add_executable(MMP  src/MMP.cpp)
target_include_directories(MMP PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_link_libraries(MMP geometry-central polyscope)
# geodesic path
set(drawGeoPath
      src/geodesicpath.cpp
      src/heat-method.cpp
      ../../utils/src/colormap.cpp
      ../../core/src/geometry.cpp
    )

# Change the executable name by changing "main" in the lines below.
add_executable(drawGeodesicPath ${drawGeoPath})
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main)
target_include_directories(drawGeodesicPath PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(drawGeodesicPath PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# add the args.hxx project which we use for command line args
target_include_directories(drawGeodesicPath PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../extern/polyscope/deps/args")
target_link_libraries(drawGeodesicPath geometry-central polyscope)