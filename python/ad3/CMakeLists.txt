include_directories(${CMAKE_CURRENT_LIST_DIR})
add_cython_target(python_base base.pyx CXX)
add_cython_target(python_extensions extensions.pyx CXX)
add_cython_target(python_factor_graph factor_graph.pyx CXX)

add_library(python_base MODULE ${python_base})
python_extension_module(base
                        LINKED_MODULES_VAR linked_modules_list
                        FORWARD_DECL_MODULES_VAR fdecl_module_list)

add_library(python_extensions MODULE ${python_extensions})
python_extension_module(extensions
                        LINKED_MODULES_VAR linked_modules_list
                        FORWARD_DECL_MODULES_VAR fdecl_module_list)

add_library(python_factor_graph MODULE ${python_factor_graph})
python_extension_module(factor_graph
                        LINKED_MODULES_VAR linked_modules_list
                        FORWARD_DECL_MODULES_VAR fdecl_module_list)

install(TARGETS python_base python_extensions python_factor_graph
        ARCHIVE DESTINATION python/ad3
        LIBRARY DESTINATION python/ad3)

python_modules_header(modules
                      FORWARD_DECL_MODULES_LIST ${fdecl_module_list})
