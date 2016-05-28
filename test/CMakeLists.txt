## Declare a C++ executable
add_executable(qb_legs_ik_test qb_legs_ik_test.cpp)

## Add cmake target dependencies of the executable
## same as for the library above
add_dependencies(qb_legs_ik_test ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Specify libraries to link a library or executable target against
target_link_libraries(qb_legs_ik_test
    qb_legs
    ${catkin_LIBRARIES}
)
