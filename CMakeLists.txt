cmake_minimum_required(VERSION 3.2)
project(raremetal VERSION 4.15.1 LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 11)

find_package(Threads)
find_library(LIB_STATGEN StatGen)
find_library(LIB_SAVVY savvy)
find_library(LIB_ZSTD zstd)
find_library(LIB_HTS hts)
find_library(LIB_Z z)
find_library(LIB_RMATH Rmath)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

add_library(libraremetal
            libRareMetal/mvtnorm/randomF77.c
            libRareMetal/pdf/PDF.cpp
            libRareMetal/pdf/PDF.h
            libRareMetal/pdf/PDFchartaxis.cpp
            libRareMetal/pdf/PDFchartaxis.h
            libRareMetal/pdf/PDFchartbar.cpp
            libRareMetal/pdf/PDFchartbar.h
            libRareMetal/pdf/PDFchartbasics.cpp
            libRareMetal/pdf/PDFchartbasics.h
            libRareMetal/pdf/PDFchartlegend.cpp
            libRareMetal/pdf/PDFchartlegend.h
            libRareMetal/pdf/PDFchartline.cpp
            libRareMetal/pdf/PDFchartline.h
            libRareMetal/pdf/PDFchartmarker.cpp
            libRareMetal/pdf/PDFchartmarker.h
            libRareMetal/pdf/PDFchartobject.cpp
            libRareMetal/pdf/PDFchartobject.h
            libRareMetal/pdf/PDFfont.cpp
            libRareMetal/pdf/PDFfont.h
            libRareMetal/pdf/PDFgrid.cpp
            libRareMetal/pdf/PDFgrid.h
            libRareMetal/pdf/PDFgridcell.cpp
            libRareMetal/pdf/PDFgridcell.h
            libRareMetal/pdf/PDFhistogram.cpp
            libRareMetal/pdf/PDFhistogram.h
            libRareMetal/pdf/PDFinfo.cpp
            libRareMetal/pdf/PDFinfo.h
            libRareMetal/pdf/PDFlinechart.cpp
            libRareMetal/pdf/PDFlinechart.h
            libRareMetal/pdf/PDFlinechartwithpolygon.cpp
            libRareMetal/pdf/PDFlinechartwithpolygon.h
            libRareMetal/pdf/PDFmanhattan.cpp
            libRareMetal/pdf/PDFmanhattan.h
            libRareMetal/pdf/PDFpage.cpp
            libRareMetal/pdf/PDFpage.h
            libRareMetal/pdf/PDFpageobject.cpp
            libRareMetal/pdf/PDFpageobject.h
            libRareMetal/src/Kinship.cpp
            libRareMetal/src/Kinship.h
            libRareMetal/src/KinshipX.cpp
            libRareMetal/src/KinshipX.h
            libRareMetal/src/MathCholesky.cpp
            libRareMetal/src/MathCholesky.h
            libRareMetal/src/MathDeriv.cpp
            libRareMetal/src/MathDeriv.h
            libRareMetal/src/MathFloatVector.cpp
            libRareMetal/src/MathFloatVector.h
            libRareMetal/src/MathGenMin.cpp
            libRareMetal/src/MathGenMin.h
            libRareMetal/src/MathGold.cpp
            libRareMetal/src/MathGold.h
            libRareMetal/src/MathNormal.cpp
            libRareMetal/src/MathNormal.h
            libRareMetal/src/MathSVD.cpp
            libRareMetal/src/MathSVD.h
            libRareMetal/src/Rmath.h
            libRareMetal/src/TraitTransformations.cpp
            libRareMetal/src/TraitTransformations.h
            libRareMetal/mvtnorm/mvt.f)

set_target_properties(libraremetal PROPERTIES OUTPUT_NAME raremetal)
target_include_directories(libraremetal PUBLIC libRareMetal/pdf libRareMetal/src)
target_link_libraries(libraremetal ${LIB_STATGEN} ${LIB_HTS} ${CMAKE_THREAD_LIBS_INIT} ${LIB_Z} ${LIB_RMATH})

add_executable(raremetal
#               raremetal/src/Calculate_mvt_pvalue.h
#               raremetal/src/group.cpp
               raremetal/src/GroupFromAnnotation.cpp
               raremetal/src/GroupFromAnnotation.h
               raremetal/src/Main.cpp
               raremetal/src/Meta.cpp
               raremetal/src/Meta.h
               raremetal/src/MetaUtility.cpp
               raremetal/src/MetaUtility.h
#               raremetal/src/MixChidist.h
#               raremetal/src/My_mvt.h
#               raremetal/src/qfc.h
               raremetal/src/QuadProg.cpp
               raremetal/src/QuadProg.h
#               raremetal/src/Rmath.h
               raremetal/src/SummaryFileReader.cpp
               raremetal/src/SummaryFileReader.h
#               raremetal/src/WriteLog.h
               raremetal/src/WritePDF.cpp
               raremetal/src/WritePDF.h)

target_include_directories(raremetal PRIVATE otherLib/eigen-3.2.0)
target_link_libraries(raremetal libraremetal)
target_compile_definitions(raremetal PRIVATE VERSION="${PROJECT_VERSION}")

add_executable(raremetalworker
               raremetalworker/src/AutoFit.cpp
               raremetalworker/src/AutoFit.h
               raremetalworker/src/AutoFit2.cpp
               raremetalworker/src/AutoFit2.h
               #raremetalworker/src/CheckRef.cpp
               #raremetalworker/src/CheckRef.h
               raremetalworker/src/DataQC.cpp
               raremetalworker/src/DataQC.h
               raremetalworker/src/FastFit.cpp
               raremetalworker/src/FastFit.h
               raremetalworker/src/GroupFromAnnotation.cpp
               raremetalworker/src/GroupFromAnnotation.h
               raremetalworker/src/KinshipEmp.cpp
               raremetalworker/src/KinshipEmp.h
               raremetalworker/src/logistic.cpp
               raremetalworker/src/logistic.h
               raremetalworker/src/Main.cpp
               raremetalworker/src/OutputKin.cpp
               raremetalworker/src/OutputKin.h
               raremetalworker/src/PreMeta.cpp
               raremetalworker/src/PreMeta.h
               #raremetalworker/src/Rmath.h
               #raremetalworker/src/snpHWE.h
               raremetalworker/src/TransformResiduals.cpp
               raremetalworker/src/TransformResiduals.h
               #raremetalworker/src/WriteLog.h
               raremetalworker/src/WritePDF.cpp
               raremetalworker/src/WritePDF.h)

target_include_directories(raremetalworker PRIVATE otherLib/eigen-3.2.0)
target_link_libraries(raremetalworker libraremetal ${LIB_SAVVY} ${LIB_ZSTD})
target_compile_definitions(raremetalworker PRIVATE VERSION="${PROJECT_VERSION}")

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
    add_subdirectory(raremetal/tests)
endif()