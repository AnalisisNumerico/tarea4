# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.9

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jeanpaul/Code/c++/tarea4

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jeanpaul/Code/c++/tarea4/cmake-build-debug

# Include any dependencies generated for this target.
include test/CMakeFiles/tester.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/tester.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/tester.dir/flags.make

test/CMakeFiles/tester.dir/testAllocator.cpp.o: test/CMakeFiles/tester.dir/flags.make
test/CMakeFiles/tester.dir/testAllocator.cpp.o: ../test/testAllocator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jeanpaul/Code/c++/tarea4/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/tester.dir/testAllocator.cpp.o"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tester.dir/testAllocator.cpp.o -c /home/jeanpaul/Code/c++/tarea4/test/testAllocator.cpp

test/CMakeFiles/tester.dir/testAllocator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tester.dir/testAllocator.cpp.i"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jeanpaul/Code/c++/tarea4/test/testAllocator.cpp > CMakeFiles/tester.dir/testAllocator.cpp.i

test/CMakeFiles/tester.dir/testAllocator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tester.dir/testAllocator.cpp.s"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jeanpaul/Code/c++/tarea4/test/testAllocator.cpp -o CMakeFiles/tester.dir/testAllocator.cpp.s

test/CMakeFiles/tester.dir/testAllocator.cpp.o.requires:

.PHONY : test/CMakeFiles/tester.dir/testAllocator.cpp.o.requires

test/CMakeFiles/tester.dir/testAllocator.cpp.o.provides: test/CMakeFiles/tester.dir/testAllocator.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/tester.dir/build.make test/CMakeFiles/tester.dir/testAllocator.cpp.o.provides.build
.PHONY : test/CMakeFiles/tester.dir/testAllocator.cpp.o.provides

test/CMakeFiles/tester.dir/testAllocator.cpp.o.provides.build: test/CMakeFiles/tester.dir/testAllocator.cpp.o


test/CMakeFiles/tester.dir/testLU.cpp.o: test/CMakeFiles/tester.dir/flags.make
test/CMakeFiles/tester.dir/testLU.cpp.o: ../test/testLU.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jeanpaul/Code/c++/tarea4/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object test/CMakeFiles/tester.dir/testLU.cpp.o"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tester.dir/testLU.cpp.o -c /home/jeanpaul/Code/c++/tarea4/test/testLU.cpp

test/CMakeFiles/tester.dir/testLU.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tester.dir/testLU.cpp.i"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jeanpaul/Code/c++/tarea4/test/testLU.cpp > CMakeFiles/tester.dir/testLU.cpp.i

test/CMakeFiles/tester.dir/testLU.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tester.dir/testLU.cpp.s"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jeanpaul/Code/c++/tarea4/test/testLU.cpp -o CMakeFiles/tester.dir/testLU.cpp.s

test/CMakeFiles/tester.dir/testLU.cpp.o.requires:

.PHONY : test/CMakeFiles/tester.dir/testLU.cpp.o.requires

test/CMakeFiles/tester.dir/testLU.cpp.o.provides: test/CMakeFiles/tester.dir/testLU.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/tester.dir/build.make test/CMakeFiles/tester.dir/testLU.cpp.o.provides.build
.PHONY : test/CMakeFiles/tester.dir/testLU.cpp.o.provides

test/CMakeFiles/tester.dir/testLU.cpp.o.provides.build: test/CMakeFiles/tester.dir/testLU.cpp.o


test/CMakeFiles/tester.dir/testMain.cpp.o: test/CMakeFiles/tester.dir/flags.make
test/CMakeFiles/tester.dir/testMain.cpp.o: ../test/testMain.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jeanpaul/Code/c++/tarea4/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object test/CMakeFiles/tester.dir/testMain.cpp.o"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tester.dir/testMain.cpp.o -c /home/jeanpaul/Code/c++/tarea4/test/testMain.cpp

test/CMakeFiles/tester.dir/testMain.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tester.dir/testMain.cpp.i"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jeanpaul/Code/c++/tarea4/test/testMain.cpp > CMakeFiles/tester.dir/testMain.cpp.i

test/CMakeFiles/tester.dir/testMain.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tester.dir/testMain.cpp.s"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jeanpaul/Code/c++/tarea4/test/testMain.cpp -o CMakeFiles/tester.dir/testMain.cpp.s

test/CMakeFiles/tester.dir/testMain.cpp.o.requires:

.PHONY : test/CMakeFiles/tester.dir/testMain.cpp.o.requires

test/CMakeFiles/tester.dir/testMain.cpp.o.provides: test/CMakeFiles/tester.dir/testMain.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/tester.dir/build.make test/CMakeFiles/tester.dir/testMain.cpp.o.provides.build
.PHONY : test/CMakeFiles/tester.dir/testMain.cpp.o.provides

test/CMakeFiles/tester.dir/testMain.cpp.o.provides.build: test/CMakeFiles/tester.dir/testMain.cpp.o


test/CMakeFiles/tester.dir/testMatrix.cpp.o: test/CMakeFiles/tester.dir/flags.make
test/CMakeFiles/tester.dir/testMatrix.cpp.o: ../test/testMatrix.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jeanpaul/Code/c++/tarea4/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object test/CMakeFiles/tester.dir/testMatrix.cpp.o"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tester.dir/testMatrix.cpp.o -c /home/jeanpaul/Code/c++/tarea4/test/testMatrix.cpp

test/CMakeFiles/tester.dir/testMatrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tester.dir/testMatrix.cpp.i"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jeanpaul/Code/c++/tarea4/test/testMatrix.cpp > CMakeFiles/tester.dir/testMatrix.cpp.i

test/CMakeFiles/tester.dir/testMatrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tester.dir/testMatrix.cpp.s"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jeanpaul/Code/c++/tarea4/test/testMatrix.cpp -o CMakeFiles/tester.dir/testMatrix.cpp.s

test/CMakeFiles/tester.dir/testMatrix.cpp.o.requires:

.PHONY : test/CMakeFiles/tester.dir/testMatrix.cpp.o.requires

test/CMakeFiles/tester.dir/testMatrix.cpp.o.provides: test/CMakeFiles/tester.dir/testMatrix.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/tester.dir/build.make test/CMakeFiles/tester.dir/testMatrix.cpp.o.provides.build
.PHONY : test/CMakeFiles/tester.dir/testMatrix.cpp.o.provides

test/CMakeFiles/tester.dir/testMatrix.cpp.o.provides.build: test/CMakeFiles/tester.dir/testMatrix.cpp.o


test/CMakeFiles/tester.dir/testQR.cpp.o: test/CMakeFiles/tester.dir/flags.make
test/CMakeFiles/tester.dir/testQR.cpp.o: ../test/testQR.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jeanpaul/Code/c++/tarea4/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object test/CMakeFiles/tester.dir/testQR.cpp.o"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tester.dir/testQR.cpp.o -c /home/jeanpaul/Code/c++/tarea4/test/testQR.cpp

test/CMakeFiles/tester.dir/testQR.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tester.dir/testQR.cpp.i"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jeanpaul/Code/c++/tarea4/test/testQR.cpp > CMakeFiles/tester.dir/testQR.cpp.i

test/CMakeFiles/tester.dir/testQR.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tester.dir/testQR.cpp.s"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jeanpaul/Code/c++/tarea4/test/testQR.cpp -o CMakeFiles/tester.dir/testQR.cpp.s

test/CMakeFiles/tester.dir/testQR.cpp.o.requires:

.PHONY : test/CMakeFiles/tester.dir/testQR.cpp.o.requires

test/CMakeFiles/tester.dir/testQR.cpp.o.provides: test/CMakeFiles/tester.dir/testQR.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/tester.dir/build.make test/CMakeFiles/tester.dir/testQR.cpp.o.provides.build
.PHONY : test/CMakeFiles/tester.dir/testQR.cpp.o.provides

test/CMakeFiles/tester.dir/testQR.cpp.o.provides.build: test/CMakeFiles/tester.dir/testQR.cpp.o


test/CMakeFiles/tester.dir/testSolver.cpp.o: test/CMakeFiles/tester.dir/flags.make
test/CMakeFiles/tester.dir/testSolver.cpp.o: ../test/testSolver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jeanpaul/Code/c++/tarea4/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object test/CMakeFiles/tester.dir/testSolver.cpp.o"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tester.dir/testSolver.cpp.o -c /home/jeanpaul/Code/c++/tarea4/test/testSolver.cpp

test/CMakeFiles/tester.dir/testSolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tester.dir/testSolver.cpp.i"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jeanpaul/Code/c++/tarea4/test/testSolver.cpp > CMakeFiles/tester.dir/testSolver.cpp.i

test/CMakeFiles/tester.dir/testSolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tester.dir/testSolver.cpp.s"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jeanpaul/Code/c++/tarea4/test/testSolver.cpp -o CMakeFiles/tester.dir/testSolver.cpp.s

test/CMakeFiles/tester.dir/testSolver.cpp.o.requires:

.PHONY : test/CMakeFiles/tester.dir/testSolver.cpp.o.requires

test/CMakeFiles/tester.dir/testSolver.cpp.o.provides: test/CMakeFiles/tester.dir/testSolver.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/tester.dir/build.make test/CMakeFiles/tester.dir/testSolver.cpp.o.provides.build
.PHONY : test/CMakeFiles/tester.dir/testSolver.cpp.o.provides

test/CMakeFiles/tester.dir/testSolver.cpp.o.provides.build: test/CMakeFiles/tester.dir/testSolver.cpp.o


# Object files for target tester
tester_OBJECTS = \
"CMakeFiles/tester.dir/testAllocator.cpp.o" \
"CMakeFiles/tester.dir/testLU.cpp.o" \
"CMakeFiles/tester.dir/testMain.cpp.o" \
"CMakeFiles/tester.dir/testMatrix.cpp.o" \
"CMakeFiles/tester.dir/testQR.cpp.o" \
"CMakeFiles/tester.dir/testSolver.cpp.o"

# External object files for target tester
tester_EXTERNAL_OBJECTS =

test/tester: test/CMakeFiles/tester.dir/testAllocator.cpp.o
test/tester: test/CMakeFiles/tester.dir/testLU.cpp.o
test/tester: test/CMakeFiles/tester.dir/testMain.cpp.o
test/tester: test/CMakeFiles/tester.dir/testMatrix.cpp.o
test/tester: test/CMakeFiles/tester.dir/testQR.cpp.o
test/tester: test/CMakeFiles/tester.dir/testSolver.cpp.o
test/tester: test/CMakeFiles/tester.dir/build.make
test/tester: src/libanpi.a
test/tester: /usr/lib/x86_64-linux-gnu/libboost_filesystem.so
test/tester: /usr/lib/x86_64-linux-gnu/libboost_system.so
test/tester: /usr/lib/x86_64-linux-gnu/libboost_unit_test_framework.so
test/tester: test/CMakeFiles/tester.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jeanpaul/Code/c++/tarea4/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX executable tester"
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tester.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/tester.dir/build: test/tester

.PHONY : test/CMakeFiles/tester.dir/build

test/CMakeFiles/tester.dir/requires: test/CMakeFiles/tester.dir/testAllocator.cpp.o.requires
test/CMakeFiles/tester.dir/requires: test/CMakeFiles/tester.dir/testLU.cpp.o.requires
test/CMakeFiles/tester.dir/requires: test/CMakeFiles/tester.dir/testMain.cpp.o.requires
test/CMakeFiles/tester.dir/requires: test/CMakeFiles/tester.dir/testMatrix.cpp.o.requires
test/CMakeFiles/tester.dir/requires: test/CMakeFiles/tester.dir/testQR.cpp.o.requires
test/CMakeFiles/tester.dir/requires: test/CMakeFiles/tester.dir/testSolver.cpp.o.requires

.PHONY : test/CMakeFiles/tester.dir/requires

test/CMakeFiles/tester.dir/clean:
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test && $(CMAKE_COMMAND) -P CMakeFiles/tester.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/tester.dir/clean

test/CMakeFiles/tester.dir/depend:
	cd /home/jeanpaul/Code/c++/tarea4/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jeanpaul/Code/c++/tarea4 /home/jeanpaul/Code/c++/tarea4/test /home/jeanpaul/Code/c++/tarea4/cmake-build-debug /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test /home/jeanpaul/Code/c++/tarea4/cmake-build-debug/test/CMakeFiles/tester.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/tester.dir/depend

