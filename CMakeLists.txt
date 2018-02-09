cmake_minimum_required( VERSION 3.9.1 )

PROJECT( new_app )
ENABLE_TESTING()
find_package( Boost 1.53 REQUIRED COMPONENTS program_options )
if (MSVC)
  # prevent auto linking
  add_definitions( -DBOOST_ALL_NO_LIB=1 )
endif()

# global flags
set( CMAKE_CXX_STANDARD 14 )
if (CMAKE_COMPILER_IS_GNUCC OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
  add_definitions(-Wall -Wextra -Werror)  # be conservative about checks
  if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # some helpful additional checks suggested by http://blog.llvm.org/2013/09/clang-warnings.html
    add_definitions( -Wheader-guard -Wlogical-not-parentheses -Wduplicate-enum -Wdangling-field
      -Wstring-plus-int -Wstring-conversion -Wsizeof-array-argument -Wsometimes-uninitialized
      -Woverloaded-shift-op-parentheses -Wloop-analysis )
    add_definitions(-stdlib=libc++)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++ -lc++abi")
  endif()
elseif (MSVC11)
  # VC11 workaround: emplace_back requires variadic template emulation, for which VC11 has
  # a default limit of 5 arguments due to some memory explosion issue.  Raising this to 6
  # allows our calls to emplace_back (to name one example) to succeed
  add_definitions(-D_VARIADIC_MAX=6)
endif()

# for packages without a standard CMake find script:
set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

# packages relied on by multiple directories
find_package(TCL 8.6)

add_subdirectory( parser )
add_subdirectory( gui )
add_subdirectory( db )
add_subdirectory( apps )
if( TCL_FOUND )
  add_subdirectory( tclint )
endif()
add_subdirectory( analysis )
