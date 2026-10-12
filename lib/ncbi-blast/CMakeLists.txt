# Vendored NCBI BLAST Karlin mathematics (C, upstream-faithful).
# Original filenames kept so `diff -u` against
# github.com/ncbi/ncbi-cxx-toolkit-public stays readable.
add_library(ncbi_blast OBJECT
        blast_stat.c
        ncbi_math.c
        )
target_include_directories(ncbi_blast PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# Upstream C, keep quiet; never inherit C++ -fno-exceptions/-std flags.
set_target_properties(ncbi_blast PROPERTIES COMPILE_FLAGS "-w")
