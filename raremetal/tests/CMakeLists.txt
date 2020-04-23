add_executable(tests-raremetal
  tests.cpp
  ../src/GroupFromAnnotation.cpp
  ../src/GroupFromAnnotation.h
  ../src/Meta.cpp
  ../src/Meta.h
  ../src/MetaUtility.cpp
  ../src/MetaUtility.h
  ../src/My_mvt.h
  ../src/QuadProg.cpp
  ../src/QuadProg.h
  ../src/SummaryFileReader.cpp
  ../src/SummaryFileReader.h
  ../src/WritePDF.cpp
  ../src/WritePDF.h
)
target_include_directories(tests-raremetal PRIVATE ../../otherLib/eigen-3.2.0)
target_link_libraries(tests-raremetal libraremetal)
target_compile_definitions(tests-raremetal PRIVATE VERSION="${PROJECT_VERSION}")