cmake_minimum_required(VERSION 3.26)

include(FetchContent)

project(phototools)

set(CMAKE_CXX_STANDARD 17)

# Add custom CMake module directory to path.
list( APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/" )

add_subdirectory(src/)

# Enable versioning support for vcpkg.
#
# see: https://devblogs.microsoft.com/cppblog/take-control-of-your-vcpkg-dependencies-with-versioning-support/
set(VCPKG_FEATURE_FLAGS "versions")

# used for loading and processing RAW images.
find_package(libraw CONFIG REQUIRED)

# used for loading and encoding JPEG images.
# recommend using libjpeg-turbo for speed.
find_package(JPEG REQUIRED)
find_package(libjpeg-turbo CONFIG REQUIRED)

# PNG decoding and encoding
find_package(PNG REQUIRED)

# Colorspace engine.
# see: https://github.com/mm2/Little-CMS
find_package(lcms2 CONFIG REQUIRED)

add_executable(phototools main.cpp
        image.cpp
        image.h
        raw_image.cpp
        raw_image.h
        src/engine/core/math/matrix.cpp
        src/engine/core/math/matrix.h
        src/engine/core/color.cpp
        src/engine/core/color.h
        src/engine/core/image.cpp
        src/engine/core/image.h)

# ====================
# LibRaw

# Fix for error:
#   Imported target "libraw::raw" includes non-existent path
#
#    "/Users/kat/CLionProjects/phototools/cmake-build-debug-llvm-clang/vcpkg_installed/arm64-osx/libraw"
#
#   in its INTERFACE_INCLUDE_DIRECTORIES.  Possible reasons include:
#
#   * The path was deleted, renamed, or moved to another location.
#
#   * An install or uninstall procedure did not complete successfully.
#
#   * The installation package was faulty and references files it does not
#   provide.
#
# TODO: fix in upstream vcpkg port?
get_target_property(LIBRAW_WRONG_INTERFACE_INCLUDE_DIRECTORIES libraw::raw INTERFACE_INCLUDE_DIRECTORIES)
get_filename_component(LIBRAW_WRONG_INTERFACE_INCLUDE_DIRECTORIES_PARENT
        "${LIBRAW_WRONG_INTERFACE_INCLUDE_DIRECTORIES}" DIRECTORY) # get parent directory
get_filename_component(LIBRAW_WRONG_INTERFACE_INCLUDE_DIRECTORIES_FILE
        "${LIBRAW_WRONG_INTERFACE_INCLUDE_DIRECTORIES}" NAME) # get parent directory
message(STATUS "Fixed libraw include directory ${LIBRAW_WRONG_INTERFACE_INCLUDE_DIRECTORIES_PARENT}/include/${LIBRAW_WRONG_INTERFACE_INCLUDE_DIRECTORIES_FILE}")
set_target_properties(libraw::raw PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${LIBRAW_WRONG_INTERFACE_INCLUDE_DIRECTORIES_PARENT}/include/${LIBRAW_WRONG_INTERFACE_INCLUDE_DIRECTORIES_FILE}"
) # correctly inject include/ into the last paht component

#target_link_libraries(phototools_engine_libraw PRIVATE libraw::raw_r) # thread safe vs libraw::raw
target_link_libraries(phototools PRIVATE libraw::raw) # thread safe vs libraw::raw

# ====================
# Rawspeed


# Allow rawspeed CMake list to fetch googletest.
set(ALLOW_DOWNLOADING_GOOGLETEST ON)
set(ALLOW_DOWNLOADING_GOOGLEBENCHMARK ON)

# NOTE: may need LLVM Clang on Apple to work. Apple GCC and Clang don't seem to provide it.
set(WITH_OPENMP ON)

# TODO: should be ON for safety, but throws a ton of unsafe errors currently.
set(RAWSPEED_ENABLE_WERROR OFF)

FetchContent_Declare(rawspeed
        GIT_REPOSITORY https://github.com/darktable-org/rawspeed.git
        GIT_TAG        a56970ceb8e1ad76987d1c549bb07c861d65eea3 # latest "stable" branch commit
        #0452b3f628a7e6c1b3dc9c8ba8128abcc3f05731 # v3.6
)

FetchContent_MakeAvailable(rawspeed)



if(WITH_OPENMP AND OPENMP_FOUND)
    target_link_libraries(phototools INTERFACE ${OpenMP_CXX_FLAGS})
endif()


target_link_libraries(phototools PUBLIC
        rawspeed
)


# Copy the camera metadata to the destination directory.
add_custom_command(
        TARGET phototools POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${rawspeed_SOURCE_DIR}/data/cameras.xml
        ${CMAKE_BINARY_DIR}/cameras.xml) # NOTE: not CURRENT_BINARY_DIR so it goes in the correct location.

# =====
# libjpeg-turbo
target_include_directories(phototools PRIVATE JPEG::JPEG)
target_link_libraries(phototools PRIVATE $<IF:$<TARGET_EXISTS:libjpeg-turbo::turbojpeg>,libjpeg-turbo::turbojpeg,libjpeg-turbo::turbojpeg-static>)

# =====
# libpng
target_link_libraries(phototools PRIVATE PNG::PNG)

# =====
# littleCMS (colorspace management engine)
target_link_libraries(phototools PRIVATE lcms2::lcms2)

# ====================
# Other

#target_link_libraries(phototools PUBLIC
#        #phototools_engine
#        #phototools_engine_rawspeed
#        #phototools_engine_libraw
#)

# Copy the camera metadata to the destination directory.
add_custom_command(
        TARGET phototools POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/DSCF1374.RAF
        ${CMAKE_CURRENT_BINARY_DIR}/DSCF1374.RAF)
