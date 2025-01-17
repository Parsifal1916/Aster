cmake_minimum_required(VERSION 3.15)
project(Aster VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

if(APPLE)
    list(APPEND CMAKE_PREFIX_PATH 
        "/usr/local/opt/glfw"
        "/usr/local/opt/glew"
        "/opt/homebrew/opt/glfw"
        "/opt/homebrew/opt/glew"
    )
    
    include_directories(
        "/usr/local/include"
        "/opt/homebrew/include"
    )
    
    link_directories(
        "/usr/local/lib"
        "/opt/homebrew/lib"
    )
elseif(WIN32)
    if(DEFINED ENV{VCPKG_ROOT})
        set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
            CACHE STRING "Vcpkg toolchain file")
    endif()
endif()

# Trova OpenGL
find_package(OpenGL REQUIRED)

# Trova GLEW
if(WIN32)
    find_package(GLEW REQUIRED)
else()
    find_package(GLEW REQUIRED)
    if(NOT GLEW_FOUND)
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(GLEW REQUIRED glew)
    endif()
endif()

find_package(glfw3 REQUIRED)
if(NOT glfw3_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GLFW3 REQUIRED glfw3)
endif()

file(GLOB_RECURSE ASTER_SOURCES src/*.cpp)
file(GLOB_RECURSE ASTER_HEADERS include/Aster/*.h include/Aster/*.hpp)

add_library(Aster ${ASTER_SOURCES} ${ASTER_HEADERS})

target_include_directories(Aster 
    PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        ${OPENGL_INCLUDE_DIR}
        ${GLEW_INCLUDE_DIRS}
        ${GLFW3_INCLUDE_DIRS}
)

if(WIN32)
    target_link_libraries(Aster 
        PUBLIC
            OpenGL::GL
            GLEW::GLEW
            glfw
    )
elseif(APPLE)
    target_link_libraries(Aster 
        PUBLIC
            "-framework OpenGL"
            ${GLEW_LIBRARIES}
            glfw
    )
else()
    target_link_libraries(Aster 
        PUBLIC
            OpenGL::GL
            GLEW::GLEW
            glfw
            ${CMAKE_DL_LIBS}
    )
endif()

set_target_properties(Aster PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    PUBLIC_HEADER "${ASTER_HEADERS}"
    POSITION_INDEPENDENT_CODE ON
    DEBUG_POSTFIX "d"
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    "AsterConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/AsterConfig.cmake.in"
    "${PROJECT_BINARY_DIR}/AsterConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Aster
)

# Installazione
install(FILES include/Aster.hpp 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Aster
)

install(DIRECTORY include/Aster
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)

install(TARGETS Aster
    EXPORT AsterTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Aster
)

install(EXPORT AsterTargets
    FILE AsterTargets.cmake
    NAMESPACE Aster::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Aster
)

install(FILES
    "${PROJECT_BINARY_DIR}/AsterConfig.cmake"
    "${PROJECT_BINARY_DIR}/AsterConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Aster
)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    find_package(Catch2 REQUIRED)

    add_executable(AsterTests tests/main_tests.cpp)
    
    target_link_libraries(AsterTests PRIVATE
        Catch2::Catch2
        Catch2::Catch2WithMain
        Aster
    )

    if(NOT MSVC)
        target_compile_options(AsterTests PRIVATE 
            -O3 
            -mtune=native 
            -march=native 
            -funroll-loops 
            -ffast-math 
            -Wno-narrowing 
            -Wno-aggressive-loop-optimizations
        )
    endif()

    include(CTest)
    include(Catch)
    catch_discover_tests(AsterTests)
endif()
