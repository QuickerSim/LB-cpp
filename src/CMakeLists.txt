list(APPEND lb_files
    main.cpp
    core.hpp
    lbblock.cpp
    lbblock.hpp
    quadtree.cpp
    quadtree.hpp
    mpiblock.cpp
    mpiblock.hpp
    treeid.cpp
    treeid.hpp
) 

add_executable(lb ${lb_files})

if(CMAKE_BUILD_TYPE MATCHES Profiler)
    target_link_libraries(lb easy_profiler)
endif(CMAKE_BUILD_TYPE MATCHES Profiler)
