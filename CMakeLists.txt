cmake_minimum_required(VERSION 3.20)

project(MacPhersonians)

add_executable(lower_cones creating_all_oriented_matroids/lower_cones.cpp
                           creating_all_oriented_matroids/OMs.cpp)

add_executable(translate_finschi_representatives creating_all_oriented_matroids/translate_Finschi_representatives.cpp
                                                 creating_all_oriented_matroids/OMs.cpp)
