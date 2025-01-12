cmake_minimum_required(VERSION 3.15)
project(Aster VERSION 1.0 LANGUAGES CXX)

# Impostazioni del compilatore
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)


# Raccolta dei file header
file(GLOB_RECURSE ASTER_HEADERS include/Aster/*.h include/Aster/*.hpp)

# Creazione di un target INTERFACE per la libreria header-only
add_library(Aster INTERFACE)

# Specifica le directory di inclusione
target_include_directories(Aster INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Imposta le opzioni di compilazione per chi utilizza questa libreria
target_compile_options(Aster INTERFACE -O3)

# Installa gli header
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
