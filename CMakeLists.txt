cmake_minimum_required(VERSION 3.15)
project(Aster VERSION 1.7.4 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

if(APPLE)
    set(CMAKE_MACOSX_RPATH ON)
endif()

option(BUILD_SHARED_LIBS "Build shared libraries" ON)

include(GNUInstallDirs)
include(CMakeFindDependencyMacro)

find_dependency(OpenGL REQUIRED)
find_dependency(GLEW REQUIRED)
find_dependency(glfw3 REQUIRED)
find_dependency(OpenCL REQUIRED)
find_dependency(TBB REQUIRED)
find_dependency(Threads REQUIRED)


file(GLOB_RECURSE ASTER_SOURCES src/*.cpp)
file(GLOB_RECURSE ASTER_HEADERS include/Aster/.h include/Aster/*.hpp)

add_library(Aster SHARED ${ASTER_SOURCES} ${ASTER_HEADERS})

target_include_directories(Aster PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(Aster PUBLIC
    OpenGL::GL
    GLEW::GLEW
    glfw
    OpenCL::OpenCL
    Threads::Threads
    TBB::tbb
    ${CMAKE_DL_LIBS}
)

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

set_target_properties(Aster PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    PUBLIC_HEADER "${ASTER_HEADERS}"
    POSITION_INDEPENDENT_CODE ON
    LIBRARY_OUTPUT_DIRECTORY "${LIBRARY_OUTPUT_PATH}"
    ARCHIVE_OUTPUT_DIRECTORY "${LIBRARY_OUTPUT_PATH}"
    RUNTIME_OUTPUT_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}"
)

if(NOT MSVC)
    target_compile_options(Aster PUBLIC
        -Wno-aggressive-loop-optimizations
        -O3
        -lTBB
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

if(APPLE)
    target_compile_options(Aster PRIVATE -Wno-nullability-completeness)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    find_package(Catch2 REQUIRED)

    file(GLOB_RECURSE TEST_HEADERS tests/*.h)
    add_executable(AsterTests tests/main_tests.cpp)
    
    set_target_properties(AsterTests PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}"
    )

    if(APPLE)
        set_target_properties(AsterTests PROPERTIES
            BUILD_WITH_INSTALL_RPATH TRUE
            INSTALL_RPATH "@loader_path/../lib"
        )
    else()
        set_target_properties(AsterTests PROPERTIES
            BUILD_WITH_INSTALL_RPATH TRUE
            INSTALL_RPATH "$ORIGIN/../lib"
        )
    endif()
    
    target_link_libraries(AsterTests PRIVATE
        Catch2::Catch2
        Catch2::Catch2WithMain
        Aster
        OpenCL::OpenCL
    )

    if(APPLE)
        target_compile_options(AsterTests PRIVATE -Wno-nullability-completeness)
        
        add_custom_command(TARGET AsterTests POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/lib
            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:Aster> ${CMAKE_BINARY_DIR}/lib/
            COMMENT "Copying Aster library to lib directory"
        )
    endif()

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
    
    if(NOT APPLE) 
        add_custom_command(TARGET AsterTests POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/lib
            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:Aster> ${CMAKE_BINARY_DIR}/lib/
            COMMENT "Copying Aster library to lib directory"
        )
    endif()
    
    catch_discover_tests(AsterTests)
endif()