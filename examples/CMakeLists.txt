cmake_minimum_required(VERSION 3.12...3.30)

project(ultimaille-examples LANGUAGES C CXX)
# or...
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# add_executable(naive_subdivision naive_subdivision.cpp)
# add_executable(taches taches.cpp)
# add_executable(tachequad tachequad.cpp)
# add_executable(midpoint midpoint.cpp)
# add_executable(tachepoly tachepoly.cpp)
# add_executable(tachepoly2 tachepoly2.cpp)
add_executable(midpointpoly midpointpoly.cpp)
add_executable(chart chart.cpp)
# add_executable(axe_median axe_median.cpp)
add_executable(tache_to_poly tache_to_quad.cpp)
# add_executable(poly_to_poly2 poly_to_poly2.cpp)
# add_executable(poly_aligne poly_aligne.cpp)
add_executable(framefield framefield.cpp)
add_executable(tache_singu tache_singu.cpp)
add_executable(recentre recentre.cpp)
# add_executable(poly_to_quad poly_to_quad.cpp)
add_executable(annot_poly annot_poly.cpp)
add_executable(suppr_non_disque suppr_non_disque.cpp)

#set(CMAKE_PREFIX_PATH "C:\\Program Files (x86)")
#add_library(osqp::osqpstatic STATIC IMPORTED)
#set_target_properties(osqp::osqpstatic PROPERTIES
#    IMPORTED_LOCATION "C:/Program Files (x86)/osqp/lib/osqp.lib"
#    INTERFACE_INCLUDE_DIRECTORIES "C:/Program Files (x86)/osqp/include"
#)
# find_package(osqp REQUIRED)


# target_link_libraries(taches ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
# target_link_libraries(tachequad ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
# target_link_libraries(midpoint ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
# target_link_libraries(tachepoly ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
# target_link_libraries(tachepoly2 ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
# target_link_libraries(midpointpoly ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
target_link_libraries(chart ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
# target_link_libraries(axe_median ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
target_link_libraries(tache_to_poly ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
# target_link_libraries(poly_to_poly ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
# target_link_libraries(poly_to_poly2 ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
# target_link_libraries(poly_aligne ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
target_link_libraries(framefield ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
target_link_libraries(tache_singu ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
target_link_libraries(recentre ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
# target_link_libraries(poly_to_quad ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
target_link_libraries(annot_poly ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
target_link_libraries(suppr_non_disque ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)

# target_include_directories(framefield PRIVATE ${osqp_INCLUDE_DIRS})


# Copy asset files to build directory
configure_file(assets/catorus.geogram assets/catorus.geogram COPYONLY)
configure_file(assets/tet_13225.geogram assets/tet_13225.geogram COPYONLY)
configure_file(assets/quarterdisk_tri.mesh assets/quarterdisk_tri.mesh COPYONLY)
#sphere
configure_file(assets/sphere.geogram assets/sphere.geogram COPYONLY)
configure_file(assets/spherep.geogram assets/spherep.geogram COPYONLY)
configure_file(assets/joint_with_tets.obj assets/joint_with_tets.obj COPYONLY)

configure_file(assets/cube_tri.geogram assets/cube_tri.geogram COPYONLY)
configure_file(assets/joint_sub.geogram assets/joint_sub.geogram COPYONLY)
configure_file(assets/joint_sub_sub.geogram assets/joint_sub_sub.geogram COPYONLY)
configure_file(assets/B46sub.geogram assets/B46sub.geogram COPYONLY)
configure_file(assets/B72sub.geogram assets/B72sub.geogram COPYONLY)
configure_file(assets/charts_to_polygon.geogram assets/charts_to_polygon.geogram COPYONLY)
configure_file(assets/charts.geogram assets/charts.geogram COPYONLY)
configure_file(assets/sortie.geogram assets/sortie.geogram COPYONLY)
configure_file(assets/cubesub.geogram assets/cubesub.geogram COPYONLY)
#configure_file(../build/examples/Debug/a_transformer.geogram assets/a_transformer.geogram COPYONLY)

