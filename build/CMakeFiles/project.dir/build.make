# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /snap/cmake/858/bin/cmake

# The command to remove a file.
RM = /snap/cmake/858/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/antoine/Documents/Viscous/CPU

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/antoine/Documents/Viscous/CPU/build

# Include any dependencies generated for this target.
include CMakeFiles/project.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/project.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/project.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/project.dir/flags.make

CMakeFiles/project.dir/src/project.c.o: CMakeFiles/project.dir/flags.make
CMakeFiles/project.dir/src/project.c.o: ../src/project.c
CMakeFiles/project.dir/src/project.c.o: CMakeFiles/project.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/antoine/Documents/Viscous/CPU/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/project.dir/src/project.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/project.dir/src/project.c.o -MF CMakeFiles/project.dir/src/project.c.o.d -o CMakeFiles/project.dir/src/project.c.o -c /home/antoine/Documents/Viscous/CPU/src/project.c

CMakeFiles/project.dir/src/project.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/project.dir/src/project.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/antoine/Documents/Viscous/CPU/src/project.c > CMakeFiles/project.dir/src/project.c.i

CMakeFiles/project.dir/src/project.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/project.dir/src/project.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/antoine/Documents/Viscous/CPU/src/project.c -o CMakeFiles/project.dir/src/project.c.s

CMakeFiles/project.dir/src/shaders.c.o: CMakeFiles/project.dir/flags.make
CMakeFiles/project.dir/src/shaders.c.o: ../src/shaders.c
CMakeFiles/project.dir/src/shaders.c.o: CMakeFiles/project.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/antoine/Documents/Viscous/CPU/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/project.dir/src/shaders.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/project.dir/src/shaders.c.o -MF CMakeFiles/project.dir/src/shaders.c.o.d -o CMakeFiles/project.dir/src/shaders.c.o -c /home/antoine/Documents/Viscous/CPU/src/shaders.c

CMakeFiles/project.dir/src/shaders.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/project.dir/src/shaders.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/antoine/Documents/Viscous/CPU/src/shaders.c > CMakeFiles/project.dir/src/shaders.c.i

CMakeFiles/project.dir/src/shaders.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/project.dir/src/shaders.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/antoine/Documents/Viscous/CPU/src/shaders.c -o CMakeFiles/project.dir/src/shaders.c.s

CMakeFiles/project.dir/src/viscous.c.o: CMakeFiles/project.dir/flags.make
CMakeFiles/project.dir/src/viscous.c.o: ../src/viscous.c
CMakeFiles/project.dir/src/viscous.c.o: CMakeFiles/project.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/antoine/Documents/Viscous/CPU/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/project.dir/src/viscous.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/project.dir/src/viscous.c.o -MF CMakeFiles/project.dir/src/viscous.c.o.d -o CMakeFiles/project.dir/src/viscous.c.o -c /home/antoine/Documents/Viscous/CPU/src/viscous.c

CMakeFiles/project.dir/src/viscous.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/project.dir/src/viscous.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/antoine/Documents/Viscous/CPU/src/viscous.c > CMakeFiles/project.dir/src/viscous.c.i

CMakeFiles/project.dir/src/viscous.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/project.dir/src/viscous.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/antoine/Documents/Viscous/CPU/src/viscous.c -o CMakeFiles/project.dir/src/viscous.c.s

CMakeFiles/project.dir/src/window.c.o: CMakeFiles/project.dir/flags.make
CMakeFiles/project.dir/src/window.c.o: ../src/window.c
CMakeFiles/project.dir/src/window.c.o: CMakeFiles/project.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/antoine/Documents/Viscous/CPU/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/project.dir/src/window.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/project.dir/src/window.c.o -MF CMakeFiles/project.dir/src/window.c.o.d -o CMakeFiles/project.dir/src/window.c.o -c /home/antoine/Documents/Viscous/CPU/src/window.c

CMakeFiles/project.dir/src/window.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/project.dir/src/window.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/antoine/Documents/Viscous/CPU/src/window.c > CMakeFiles/project.dir/src/window.c.i

CMakeFiles/project.dir/src/window.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/project.dir/src/window.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/antoine/Documents/Viscous/CPU/src/window.c -o CMakeFiles/project.dir/src/window.c.s

# Object files for target project
project_OBJECTS = \
"CMakeFiles/project.dir/src/project.c.o" \
"CMakeFiles/project.dir/src/shaders.c.o" \
"CMakeFiles/project.dir/src/viscous.c.o" \
"CMakeFiles/project.dir/src/window.c.o"

# External object files for target project
project_EXTERNAL_OBJECTS =

project: CMakeFiles/project.dir/src/project.c.o
project: CMakeFiles/project.dir/src/shaders.c.o
project: CMakeFiles/project.dir/src/viscous.c.o
project: CMakeFiles/project.dir/src/window.c.o
project: CMakeFiles/project.dir/build.make
project: /usr/lib/x86_64-linux-gnu/libpython3.8.so
project: /usr/lib/x86_64-linux-gnu/libGLEW.so
project: /usr/lib/x86_64-linux-gnu/libglfw.so
project: /usr/lib/x86_64-linux-gnu/libOpenGL.so
project: /usr/lib/x86_64-linux-gnu/libGLX.so
project: /usr/lib/x86_64-linux-gnu/libGLU.so
project: CMakeFiles/project.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/antoine/Documents/Viscous/CPU/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking C executable project"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/project.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/project.dir/build: project
.PHONY : CMakeFiles/project.dir/build

CMakeFiles/project.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/project.dir/cmake_clean.cmake
.PHONY : CMakeFiles/project.dir/clean

CMakeFiles/project.dir/depend:
	cd /home/antoine/Documents/Viscous/CPU/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/antoine/Documents/Viscous/CPU /home/antoine/Documents/Viscous/CPU /home/antoine/Documents/Viscous/CPU/build /home/antoine/Documents/Viscous/CPU/build /home/antoine/Documents/Viscous/CPU/build/CMakeFiles/project.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/project.dir/depend

