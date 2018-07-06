cmake_minimum_required(VERSION 3.5)
cmake_policy(SET CMP0023 NEW)

project(AD3)

if(UNIX)
	link_libraries(pthread)
endif()

include_directories(".")

set(HEADER_FILES
	examples/cpp/parsing/FactorTree.h
	ad3/Factor.h
	ad3/FactorDense.h
	ad3/FactorGraph.h
	ad3/GenericFactor.h
	ad3/MultiVariable.h
	ad3/Utils.h
	examples/cpp/dense/FactorSequence.h
	examples/cpp/parsing/FactorTree.h
	examples/cpp/parsing/FactorGrandparentHeadAutomaton.h
	examples/cpp/parsing/FactorHeadAutomaton.h
	examples/cpp/parsing/FactorTree.h
)

add_library(ad3 STATIC
		ad3/Factor.cpp
		ad3/FactorGraph.cpp
		ad3/GenericFactor.cpp
		ad3/Utils.cpp
		examples/cpp/parsing/FactorTree.cpp
		${HEADER_FILES}
)
target_include_directories(ad3 PUBLIC Eigen)

set_property(TARGET ad3 PROPERTY POSITION_INDEPENDENT_CODE TRUE)

link_libraries(ad3)


add_executable(ad3_multi
	ad3_multi.cpp
	examples/cpp/parsing/FactorTree.cpp
	examples/cpp/parsing/FactorTree.h
)
target_include_directories(ad3_multi
	PRIVATE
		examples/cpp/dense
		examples/cpp/parsing
		examples/cpp/logic
		examples/cpp/summarization
)

if(NOT SKBUILD)
    add_executable(simple_grid examples/cpp/dense/simple_grid.cpp examples/cpp/dense/FactorSequence.h)
    add_executable(simple_parser
        examples/cpp/parsing/simple_parser
        examples/cpp/parsing/FactorGrandparentHeadAutomaton.h
        examples/cpp/parsing/FactorHeadAutomaton.h
        examples/cpp/parsing/FactorTree.h
    )
    add_executable(simple_coref examples/cpp/logic/simple_coref.cpp)
endif()

if(NOT SKBUILD)
    install(TARGETS ad3
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)
    install(FILES ad3/Factor.h DESTINATION include/ad3)
    install(FILES ad3/FactorDense.h DESTINATION include/ad3)
    install(FILES ad3/FactorGraph.h DESTINATION include/ad3)
    install(FILES ad3/GenericFactor.h DESTINATION include/ad3)
    install(FILES ad3/MultiVariable.h DESTINATION include/ad3)
    install(FILES ad3/Utils.h DESTINATION include/ad3)
    install(FILES examples/cpp/parsing/FactorTree.h DESTINATION include/ad3/examples/cpp/parsing)
    install(FILES examples/cpp/dense/FactorSequence.h DESTINATION include/ad3/examples/cpp/dense)
    install(FILES examples/cpp/parsing/FactorTree.h DESTINATION include/ad3/examples/cpp/parsing)
    install(FILES examples/cpp/parsing/FactorGrandparentHeadAutomaton.h DESTINATION include/ad3/examples/cpp/parsing)
    install(FILES examples/cpp/parsing/FactorHeadAutomaton.h DESTINATION include/ad3/examples/cpp/parsing)
    install(FILES examples/cpp/parsing/FactorTree.h DESTINATION include/ad3/examples/cpp/parsing)
else()
    if(WIN32)
        set(CONDA_SUBDIR "Library/")
    else()
        set(CONDA_SUBDIR "")
    endif()
    install(TARGETS ad3
            LIBRARY DESTINATION "${CONDA_SUBDIR}lib"
            ARCHIVE DESTINATION "${CONDA_SUBDIR}lib")
    install(FILES ad3/Factor.h DESTINATION "${CONDA_SUBDIR}include/ad3")
    install(FILES ad3/FactorDense.h DESTINATION "${CONDA_SUBDIR}include/ad3")
    install(FILES ad3/FactorGraph.h DESTINATION "${CONDA_SUBDIR}include/ad3")
    install(FILES ad3/GenericFactor.h DESTINATION "${CONDA_SUBDIR}include/ad3")
    install(FILES ad3/MultiVariable.h DESTINATION "${CONDA_SUBDIR}include/ad3")
    install(FILES ad3/Utils.h DESTINATION "${CONDA_SUBDIR}include/ad3")
    install(FILES examples/cpp/parsing/FactorTree.h DESTINATION "${CONDA_SUBDIR}include/ad3/examples/cpp/parsing")
    install(FILES examples/cpp/dense/FactorSequence.h DESTINATION "${CONDA_SUBDIR}include/ad3/examples/cpp/dense")
    install(FILES examples/cpp/parsing/FactorTree.h DESTINATION "${CONDA_SUBDIR}include/ad3/examples/cpp/parsing")
    install(FILES examples/cpp/parsing/FactorGrandparentHeadAutomaton.h
            DESTINATION "${CONDA_SUBDIR}include/ad3/examples/cpp/parsing")
    install(FILES examples/cpp/parsing/FactorHeadAutomaton.h
            DESTINATION "${CONDA_SUBDIR}include/ad3/examples/cpp/parsing")
    install(FILES examples/cpp/parsing/FactorTree.h DESTINATION "${CONDA_SUBDIR}include/ad3/examples/cpp/parsing")

    find_package(PythonInterp)
    find_package(PythonLibs)
    find_package(Cython REQUIRED)
    find_package(NumPy REQUIRED)
    find_package(PythonExtensions REQUIRED)

    add_subdirectory(python)
endif()
