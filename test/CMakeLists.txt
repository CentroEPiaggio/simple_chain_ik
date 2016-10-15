## Declare a C++ executable
add_executable(simple_chain_ik_test simple_chain_ik_test.cpp)

## Add cmake target dependencies of the executable
## same as for the library above
add_dependencies(simple_chain_ik_test ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Specify libraries to link a library or executable target against
target_link_libraries(simple_chain_ik_test
    simple_chain_ik_solver
    ${catkin_LIBRARIES}
)

add_executable(chain_and_solvers_test chain_and_solvers_test_using_vito.cpp)
add_dependencies(chain_and_solvers_test ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(chain_and_solvers_test
    chain_and_solvers
    ${catkin_LIBRARIES}
)

add_executable(eigen_indexig_test eigen_indexig_test.cpp)
target_link_libraries(eigen_indexig_test ${catkin_LIBRARIES})
