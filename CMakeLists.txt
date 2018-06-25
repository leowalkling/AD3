cmake_minimum_required(VERSION 3.0)

project(AD3)

if (UNIX)
	link_libraries(pthread)
endif()

include_directories(".")

add_library(ad3 STATIC
		ad3/Factor.cpp
		ad3/FactorGraph.cpp
		ad3/GenericFactor.cpp
		ad3/Utils.cpp
		examples/cpp/parsing/FactorTree.cpp
		# Headers
		ad3/Factor.h
		ad3/FactorDense.h
		ad3/FactorGraph.h
		ad3/GenericFactor.h
		ad3/MultiVariable.h
		ad3/Utils.h
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


add_executable(simple_grid examples/cpp/dense/simple_grid.cpp)
add_executable(simple_parser examples/cpp/parsing/simple_parser)
add_executable(simple_coref examples/cpp/logic/simple_coref.cpp)

install(
    TARGETS ad3 ad3_multi simple_grid simple_parser simple_coref
    ARCHIVE DESTINATION lib
	RUNTIME DESTINATION .
	LIBRARY DESTINATION lib
)
install(TARGETS ad3
        LIBRARY DESTINATION ad3
		ARCHIVE DESTINATION ad3)
