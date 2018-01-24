# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.9

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_COMMAND = /home/tzii/Programs/clion-2017.3.1/bin/cmake/bin/cmake

# The command to remove a file.
RM = /home/tzii/Programs/clion-2017.3.1/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/tzii/Projects/bstat

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/tzii/Projects/bstat

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/home/tzii/Programs/clion-2017.3.1/bin/cmake/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/home/tzii/Programs/clion-2017.3.1/bin/cmake/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/tzii/Projects/bstat/CMakeFiles /home/tzii/Projects/bstat/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/tzii/Projects/bstat/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named bstat

# Build rule for target.
bstat: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 bstat
.PHONY : bstat

# fast build rule for target.
bstat/fast:
	$(MAKE) -f CMakeFiles/bstat.dir/build.make CMakeFiles/bstat.dir/build
.PHONY : bstat/fast

main.o: main.c.o

.PHONY : main.o

# target to build an object file
main.c.o:
	$(MAKE) -f CMakeFiles/bstat.dir/build.make CMakeFiles/bstat.dir/main.c.o
.PHONY : main.c.o

main.i: main.c.i

.PHONY : main.i

# target to preprocess a source file
main.c.i:
	$(MAKE) -f CMakeFiles/bstat.dir/build.make CMakeFiles/bstat.dir/main.c.i
.PHONY : main.c.i

main.s: main.c.s

.PHONY : main.s

# target to generate assembly for a file
main.c.s:
	$(MAKE) -f CMakeFiles/bstat.dir/build.make CMakeFiles/bstat.dir/main.c.s
.PHONY : main.c.s

src/hmm.o: src/hmm.c.o

.PHONY : src/hmm.o

# target to build an object file
src/hmm.c.o:
	$(MAKE) -f CMakeFiles/bstat.dir/build.make CMakeFiles/bstat.dir/src/hmm.c.o
.PHONY : src/hmm.c.o

src/hmm.i: src/hmm.c.i

.PHONY : src/hmm.i

# target to preprocess a source file
src/hmm.c.i:
	$(MAKE) -f CMakeFiles/bstat.dir/build.make CMakeFiles/bstat.dir/src/hmm.c.i
.PHONY : src/hmm.c.i

src/hmm.s: src/hmm.c.s

.PHONY : src/hmm.s

# target to generate assembly for a file
src/hmm.c.s:
	$(MAKE) -f CMakeFiles/bstat.dir/build.make CMakeFiles/bstat.dir/src/hmm.c.s
.PHONY : src/hmm.c.s

src/util.o: src/util.c.o

.PHONY : src/util.o

# target to build an object file
src/util.c.o:
	$(MAKE) -f CMakeFiles/bstat.dir/build.make CMakeFiles/bstat.dir/src/util.c.o
.PHONY : src/util.c.o

src/util.i: src/util.c.i

.PHONY : src/util.i

# target to preprocess a source file
src/util.c.i:
	$(MAKE) -f CMakeFiles/bstat.dir/build.make CMakeFiles/bstat.dir/src/util.c.i
.PHONY : src/util.c.i

src/util.s: src/util.c.s

.PHONY : src/util.s

# target to generate assembly for a file
src/util.c.s:
	$(MAKE) -f CMakeFiles/bstat.dir/build.make CMakeFiles/bstat.dir/src/util.c.s
.PHONY : src/util.c.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... bstat"
	@echo "... main.o"
	@echo "... main.i"
	@echo "... main.s"
	@echo "... src/hmm.o"
	@echo "... src/hmm.i"
	@echo "... src/hmm.s"
	@echo "... src/util.o"
	@echo "... src/util.i"
	@echo "... src/util.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system
