# .-----------------------------------------------------------------------.
# | SFcollapse1D                                                          |
# | Gravitational collapse of scalar fields in spherical symmetry         |
# |                                                                       |
# | Copyright (c) 2020, Leonardo Werneck                                  |
# |                                                                       |
# | This program is free software: you can redistribute it and/or modify  |
# | it under the terms of the GNU General Public License as published by  |
# | the Free Software Foundation, either version 3 of the License, or     |
# | (at your option) any later version.                                   |
# |                                                                       |
# | This program is distributed in the hope that it will be useful,       |
# | but WITHOUT ANY WARRANTY; without even the implied warranty of        |
# | MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
# | GNU General Public License for more details.                          |
# |                                                                       |
# | You should have received a copy of the GNU General Public License     |
# | along with this program.  If not, see <https://www.gnu.org/licenses/>.|
# .-----------------------------------------------------------------------.
#
# Set the src directory
SRC_DIR = ./src

# Set the doc directory
DOC_DIR = ./doc

# Set the obj directory
OBJ_DIR = ./obj

# Set the C++ compiler
CXX = g++

# Set the C++ compiler flags
CXXFLAGS = -Wall -O2 -march=native -fopenmp

# Set the objects - When a new source file is added to the code,
# an object with the same base name must be added to the list below
OBJECTS = SFcollapse1D.o      \
	  grid.o              \
          gridfunction.o      \
          initial_condition.o \
          lapse_rescaling.o   \
          newton_raphson.o    \
          evolution.o         \
          utilities.o

# Set the objects path - When a new source file is added to the code,
# an object with the same base name must be added to the list below
OBJ_PATHS = $(OBJ_DIR)/SFcollapse1D.o      \
	    $(OBJ_DIR)/grid.o              \
            $(OBJ_DIR)/gridfunction.o      \
            $(OBJ_DIR)/initial_condition.o \
            $(OBJ_DIR)/lapse_rescaling.o   \
            $(OBJ_DIR)/newton_raphson.o    \
            $(OBJ_DIR)/evolution.o         \
            $(OBJ_DIR)/utilities.o

all: SFcollapse1D

SFcollapse1D: $(OBJECTS) out_dir
	$(CXX) $(CXXFLAGS) $(OBJ_PATHS) -o SFcollapse1D

SFcollapse1D.o: $(SRC_DIR)/SFcollapse1D.cpp $(SRC_DIR)/logo.hpp $(SRC_DIR)/macros.hpp $(SRC_DIR)/grid.hpp $(SRC_DIR)/gridfunction.hpp obj_dir
	$(CXX) $(CXXFLAGS) -o $(OBJ_DIR)/SFcollapse1D.o -c $(SRC_DIR)/SFcollapse1D.cpp

grid.o: $(SRC_DIR)/grid.cpp $(SRC_DIR)/grid.hpp $(SRC_DIR)/macros.hpp obj_dir
	$(CXX) $(CXXFLAGS) -o $(OBJ_DIR)/grid.o -c $(SRC_DIR)/grid.cpp

gridfunction.o: $(SRC_DIR)/gridfunction.cpp $(SRC_DIR)/gridfunction.hpp $(SRC_DIR)/grid.hpp obj_dir
	$(CXX) $(CXXFLAGS) -o $(OBJ_DIR)/gridfunction.o -c $(SRC_DIR)/gridfunction.cpp

initial_condition.o: $(SRC_DIR)/initial_condition.cpp $(SRC_DIR)/macros.hpp $(SRC_DIR)/grid.hpp $(SRC_DIR)/gridfunction.hpp obj_dir
	$(CXX) $(CXXFLAGS) -o $(OBJ_DIR)/initial_condition.o -c $(SRC_DIR)/initial_condition.cpp

lapse_rescaling.o: $(SRC_DIR)/lapse_rescaling.cpp $(SRC_DIR)/macros.hpp $(SRC_DIR)/grid.hpp $(SRC_DIR)/gridfunction.hpp obj_dir
	$(CXX) $(CXXFLAGS) -o $(OBJ_DIR)/lapse_rescaling.o -c $(SRC_DIR)/lapse_rescaling.cpp

newton_raphson.o: $(SRC_DIR)/newton_raphson.cpp $(SRC_DIR)/macros.hpp $(SRC_DIR)/grid.hpp $(SRC_DIR)/gridfunction.hpp obj_dir
	$(CXX) $(CXXFLAGS) -o $(OBJ_DIR)/newton_raphson.o -c $(SRC_DIR)/newton_raphson.cpp

evolution.o: $(SRC_DIR)/evolution.cpp $(SRC_DIR)/macros.hpp $(SRC_DIR)/grid.hpp $(SRC_DIR)/gridfunction.hpp obj_dir
	$(CXX) $(CXXFLAGS) -o $(OBJ_DIR)/evolution.o -c $(SRC_DIR)/evolution.cpp

utilities.o: $(SRC_DIR)/utilities.cpp $(SRC_DIR)/utilities.hpp $(SRC_DIR)/macros.hpp $(SRC_DIR)/grid.hpp obj_dir
	$(CXX) $(CXXFLAGS) -o $(OBJ_DIR)/utilities.o -c $(SRC_DIR)/utilities.cpp

obj_dir:
	mkdir -p obj

out_dir:
	mkdir -p out

clean:
	rm -rf $(OBJ_PATHS) SFcollapse1D

realclean:
	rm -rf $(OBJ_PATHS) SFcollapse1D *.txt out/*.dat animations/*.gif out/ obj/
