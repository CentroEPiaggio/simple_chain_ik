find_package(orocos_kdl REQUIRED)
message("orocos_kdl_LIBRARIES = ${orocos_kdl_LIBRARIES} - version = ${orocos_kdl_VERSION}")
# if(${orocos_kdl_VERSION} VERSION_GREATER 1.3)
#     message(FATAL_ERROR "The version of ${orocos_kdl_LIBRARIES} is >= 1.4: you should now remove the internal kdl library needed before (the whole my_kdl subdirectory) and remove dependency from it in the solvers!")
# endif()

add_library(my_kdl_lib
    my_chainiksolvervel_wdls.cpp
)
target_link_libraries(my_kdl_lib
    ${orocos_kdl_LIBRARIES}
)
