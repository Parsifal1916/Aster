cmake_minimum_required(VERSION 3.15)
project(Aster VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)

include(GNUInstallDirs)
include(CMakeFindDependencyMacro)

find_dependency(OpenGL REQUIRED)
find_dependency(GLEW REQUIRED)
find_dependency(glfw3 REQUIRED)
find_dependency(OpenCL REQUIRED)

file(GLOB_RECURSE ASTER_SOURCES src/*.cpp)
file(GLOB_RECURSE ASTER_HEADERS include/Aster/.h include/Aster/*.hpp)

add_library(Aster ${ASTER_SOURCES} ${ASTER_HEADERS})

target_include_directories(Aster PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(Aster PRIVATE
    OpenGL::GL
    GLEW::GLEW
    glfw
    OpenCL::OpenCL
    ${CMAKE_DL_LIBS}
)

set_target_properties(Aster PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    PUBLIC_HEADER "${ASTER_HEADERS}"
    POSITION_INDEPENDENT_CODE ON

    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}"
)

if(NOT MSVC)
target_compile_options(Aster PRIVATE 
    -Wno-aggressive-loop-optimizations
    -O3
    -lOpenCL
)
endif()

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

install(FILES include/Aster.hpp 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/Aster
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
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

set(LIBRARY_PATH "${CMAKE_INSTALL_PREFIX}/include/Aster/thirdparty")
target_compile_definitions(Aster PRIVATE LIBRARY_PATH="${LIBRARY_PATH}")

file(COPY ${CMAKE_SOURCE_DIR}/include/Aster/thirdparty/roboto.ttf
     DESTINATION ${CMAKE_BINARY_DIR}/include/Aster/thirdparty)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    find_package(Catch2 REQUIRED)

    file(GLOB_RECURSE TEST_HEADERS tests/*.h)
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
