# Test scenarios:
#   Two tools, many parameters for each
#       - Raremetalworker
#       - Raremetal


set(RMW_BIN ${CMAKE_BINARY_DIR}/raremetalworker)
set(RM_BIN ${CMAKE_BINARY_DIR}/raremetal)


#### RAREMETAL TUTORIAL
# Make test assets (inputs) available in build folder
file(COPY datasets/raremetal_tutorial
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


## Prepare the basic directory structure (in the build folder) for tests to run
# Each test has a base folder with `inputs`, `expected`, and `outputs`
# For a first-stage analysis, there may be no special inputs. Outputs will be wiped on every test run.
set(BASE_TEST_FOLDER ${CMAKE_CURRENT_BINARY_DIR}/test_tut_rmw)
file(MAKE_DIRECTORY ${BASE_TEST_FOLDER})
file(MAKE_DIRECTORY ${BASE_TEST_FOLDER}/inputs)
file(MAKE_DIRECTORY ${BASE_TEST_FOLDER}/expected)

file(COPY raremetalworker/tutorial/expected
     DESTINATION ${BASE_TEST_FOLDER})

add_test(NAME test_tut_rmw
         COMMAND ${CMAKE_COMMAND}
            -DRMW_BIN=${RMW_BIN}
            -DBASE_TEST_FOLDER=${BASE_TEST_FOLDER}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/raremetalworker/script.cmake)


set(BASE_TEST_FOLDER ${CMAKE_CURRENT_BINARY_DIR}/test_tut_rm)
file(MAKE_DIRECTORY ${BASE_TEST_FOLDER})
file(MAKE_DIRECTORY ${BASE_TEST_FOLDER}/inputs)
file(MAKE_DIRECTORY ${BASE_TEST_FOLDER}/expected)

file(COPY raremetal/tutorial/expected
        DESTINATION ${BASE_TEST_FOLDER})
file(COPY raremetal/tutorial/expected
        DESTINATION ${BASE_TEST_FOLDER})
file(COPY raremetal/tutorial/inputs
        DESTINATION ${BASE_TEST_FOLDER})

add_test(NAME test_tut_rm
        COMMAND ${CMAKE_COMMAND}
        -DRM_BIN=${RM_BIN}
        -DBASE_TEST_FOLDER=${BASE_TEST_FOLDER}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/raremetal/script.cmake)
