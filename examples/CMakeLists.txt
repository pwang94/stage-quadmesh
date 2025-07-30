cmake_minimum_required(VERSION 3.12...3.30)

project(ultimaille-examples LANGUAGES C CXX)
# or...
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_executable(chart chart.cpp)
add_executable(tache_to_poly tache_to_poly.cpp)
add_executable(framefield framefield.cpp)
add_executable(tache_singu tache_singu.cpp)
add_executable(recentre recentre.cpp)
add_executable(annot_poly annot_poly.cpp)
add_executable(suppr_non_disque suppr_non_disque.cpp)



target_link_libraries(chart ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
target_link_libraries(tache_to_poly ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
target_link_libraries(framefield ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
target_link_libraries(tache_singu ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
target_link_libraries(recentre ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
target_link_libraries(annot_poly ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
target_link_libraries(suppr_non_disque ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)



# Copy asset files to build directory
configure_file(assets/joint_sub_sub.geogram assets/joint_sub_sub.geogram COPYONLY)

