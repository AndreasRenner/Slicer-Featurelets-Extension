cmake_minimum_required(VERSION 2.8.9)

project(Slicer-Featurelets)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://www.example.com/Slicer/Extensions/Slicer-Featurelets")
set(EXTENSION_CATEGORY "Registration")
set(EXTENSION_CONTRIBUTORS "Andreas Renner (AnyWare Corp.)")
set(EXTENSION_DESCRIPTION "This Extension uses regularly distributed image subvolumes to achiev multimodal deformable registration.")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/Slicer-Featurelets.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/Slicer-Featurelets/Screenshots/1.png")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(Featurelets-Registration)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
