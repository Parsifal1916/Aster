cmake_minimum_required(VERSION 3.15)
project(Aster VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)

find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(glfw3 REQUIRED)

file(GLOB_RECURSE ASTER_HEADERS include/Aster/*.h include/Aster/*.hpp)

add_library(Aster INTERFACE)

target_include_directories(Aster INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(Aster INTERFACE
    OpenGL::GL
    GLEW::GLEW
    glfw
)

# Configurazione per l'installazione
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

# Installazione degli header
install(FILES include/Aster.hpp 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/Aster
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Target di installazione
install(TARGETS Aster
    EXPORT AsterTargets
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
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

# Parte dei test...
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    find_package(Catch2 REQUIRED)

    file(GLOB_RECURSE TEST_HEADERS tests/*.h)
    add_executable(AsterTests tests/main_tests.cpp)

    target_link_libraries(AsterTests PRIVATE
        Catch2::Catch2
        Catch2::Catch2WithMain
        OpenGL::GL
        GLEW::GLEW
        glfw
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
