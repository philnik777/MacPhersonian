cmake_minimum_required(VERSION 3.20)

project(MacPhersonians)

add_executable(lower_cones creating_all_oriented_matroids/lower_cones.cpp
                           creating_all_oriented_matroids/OMs.cpp)

add_executable(translate_finschi_representatives creating_all_oriented_matroids/translate_Finschi_representatives.cpp
                                                 creating_all_oriented_matroids/OMs.cpp)

target_compile_options(lower_cones PRIVATE
                       -Wno-vla-cxx-extension # FIXME: Remove all VLAs
                      )

target_compile_options(translate_finschi_representatives PRIVATE
                       -Wno-vla-cxx-extension # FIXME: Remove all VLAs
                      )

set_target_properties(lower_cones PROPERTIES
                      CXX_STANDARD 23
                      CXX_STANDARD_REQUIRED On
                      CXX_EXTENSIONS Off)

set_target_properties(translate_finschi_representatives PROPERTIES
                      CXX_STANDARD 23
                      CXX_STANDARD_REQUIRED On
                      CXX_EXTENSIONS Off)
