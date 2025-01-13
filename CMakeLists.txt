cmake_minimum_required(VERSION 3.15)
project(Aster VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)

file(GLOB_RECURSE ASTER_HEADERS include/Aster/*.h include/Aster/*.hpp)

add_library(Aster INTERFACE)

target_include_directories(Aster INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    find_package(SFML  REQUIRED COMPONENTS system window graphics audio network)
    find_package(Catch2 REQUIRED)

    file(GLOB_RECURSE ASTER_TEST_SOURCES tests/*.cpp)
    add_executable(AsterTests ${ASTER_TEST_SOURCES})

    target_link_libraries(AsterTests PRIVATE
        Catch2::Catch2
        Catch2::Catch2WithMain
        sfml-system
        sfml-window
        sfml-graphics
        Aster
    )

    target_include_directories(AsterTests PRIVATE ${CMAKE_SOURCE_DIR}/include)

    target_compile_options(AsterTests PRIVATE -O3 -mtune=native -march=native -funroll-loops -ffast-math -fcx-limited-range -Wno-narrowing -Wno-aggressive-loop-optimizations)

    include(CTest)
    include(Catch)
    catch_discover_tests(AsterTests)
endif()
