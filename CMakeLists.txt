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

install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Aster/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Aster
)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    # Find Catch2 for testing
    find_package(Catch2 REQUIRED)

    file(GLOB_RECURSE ASTER_TEST_SOURCES tests/*.cpp)
    add_executable(AsterTests ${ASTER_TEST_SOURCES})

    target_link_libraries(AsterTests PRIVATE
        Catch2::Catch2
        Catch2::Catch2WithMain
        OpenGL::GL
        GLEW::GLEW
        glfw
        Aster
    )

    target_include_directories(AsterTests PRIVATE ${CMAKE_SOURCE_DIR}/include)

    if(MSVC)
        target_compile_options(AsterTests PRIVATE 
            /O2 
            /arch:AVX2
        )
    else()
        target_compile_options(AsterTests PRIVATE 
            -O3 
            -mtune=native 
            -march=native 
            -funroll-loops 
            -ffast-math 
            -fcx-limited-range 
            -Wno-narrowing 
            -Wno-aggressive-loop-optimizations
        )
    endif()

    include(CTest)
    include(Catch)
    catch_discover_tests(AsterTests)
endif()