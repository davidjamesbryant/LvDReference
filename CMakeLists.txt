cmake_minimum_required(VERSION 3.20)
project(LvDReference LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ------------------------------------------------------------
# Shared include directories
# ------------------------------------------------------------
set(LVD_INCLUDE_DIRS
        src/Phylib
        src/DecompositionTree
        src/StandardLikelihood
        src/Eigen
)

# ------------------------------------------------------------
# Shared sources built as a static library
# ------------------------------------------------------------
add_library(lvd_core STATIC
        src/DecompositionTree/decompositionTree.cpp
        src/DecompositionTree/decompositionLikelihood.cpp
        src/StandardLikelihood/standardLikelihood.cpp
        src/Phylib/sequences/phylip_seq.cpp
)

target_include_directories(lvd_core PUBLIC ${LVD_INCLUDE_DIRS})

target_compile_options(lvd_core PUBLIC
        $<$<CONFIG:Debug>:-g>
        $<$<CONFIG:Release>:-O3 -Wno-deprecated-declarations>
)

# Uncomment if you want native optimisation locally (not portable):
# target_compile_options(lvd_core PUBLIC $<$<CONFIG:Release>:-march=native>)

# ------------------------------------------------------------
# Executables
# ------------------------------------------------------------
foreach(TARGET RunSimulation TreeHeights LvDLikelihoods TidyTree)
    add_executable(${TARGET} src/${TARGET}.cpp)
    target_link_libraries(${TARGET} PRIVATE lvd_core)
endforeach()
