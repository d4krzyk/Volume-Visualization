cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(VolumeRenderer)

# Ustawienia C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Znalezienie bibliotek VTK
find_package(VTK COMPONENTS 
  CommonColor
  CommonCore
  CommonDataModel
  FiltersCore
  FiltersPoints
  InteractionStyle
  RenderingContextOpenGL2
  RenderingVolume
  FiltersGeometry
  RenderingVolumeOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
  RenderingVolume
  RenderingVolumeOpenGL2
  GUISupportQt
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "Vol: Unable to find the VTK build folder.")
endif()

# Znalezienie Qt (w zależności od wersji)
find_package(Qt5 REQUIRED COMPONENTS Core Gui Widgets OpenGL Concurrent)

if (NOT Qt5_FOUND)
  message(FATAL_ERROR "Qt5 not found!")
endif()

# Konfiguracja automatycznego generowania plików MOC/UIC/RCC
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# Zapobieganie błędom związanym z długimi ścieżkami w Windows
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")

# Dodanie pliku wykonywalnego
add_executable(VolumeRenderer MACOSX_BUNDLE Vol.cxx)

# Powiązanie bibliotek Qt5 i VTK
target_link_libraries(VolumeRenderer PRIVATE 
  Qt5::Core
  Qt5::Gui
  Qt5::Widgets
  Qt5::OpenGL
  Qt5::Concurrent
  ${VTK_LIBRARIES}
)

# vtk_module_autoinit jest wymagane
vtk_module_autoinit(
  TARGETS VolumeRenderer
  MODULES ${VTK_LIBRARIES}
)
