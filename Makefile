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
OBJECTS = SFcollapse1D.o \
	  grid.o         \
          gridfunction.o \
          evolution.o    \
          utilities.o

# Set the objects path - When a new source file is added to the code,
# an object with the same base name must be added to the list below
OBJ_PATHS = $(OBJ_DIR)/SFcollapse1D.o \
	    $(OBJ_DIR)/grid.o         \
            $(OBJ_DIR)/gridfunction.o \
            $(OBJ_DIR)/evolution.o    \
            $(OBJ_DIR)/utilities.o

# Set header files for the SFcollapse1D.cpp source file
SFcollapse1D_headers = $(SRC_DIR)/logo.hpp         \
                       $(SRC_DIR)/grid.hpp         \
                       $(SRC_DIR)/macros.hpp       \
                       $(SRC_DIR)/utilities.hpp    \
                       $(SRC_DIR)/evolution.hpp    \
                       $(SRC_DIR)/gridfunction.hpp

# Set header files for the grid.cpp source file
grid_headers = $(SRC_DIR)/macros.hpp \
               $(SRC_DIR)/grid.hpp

# Set header files for the gridfunction.cpp source file
gridfunction_headers = $(SRC_DIR)/grid.hpp         \
                       $(SRC_DIR)/utilities.hpp    \
                       $(SRC_DIR)/gridfunction.hpp

# Set header files for the evolution.cpp source file
evolution_headers = $(SRC_DIR)/grid.hpp         \
                    $(SRC_DIR)/macros.hpp       \
                    $(SRC_DIR)/utilities.hpp    \
                    $(SRC_DIR)/evolution.hpp    \
                    $(SRC_DIR)/gridfunction.hpp

# Set header files for the utilities.cpp source file
utilities_headers = $(SRC_DIR)/grid.hpp         \
                    $(SRC_DIR)/macros.hpp       \
                    $(SRC_DIR)/utilities.hpp    \
                    $(SRC_DIR)/gridfunction.hpp

all: SFcollapse1D

SFcollapse1D: $(OBJECTS) out_dir
	$(CXX) $(CXXFLAGS) $(OBJ_PATHS) -o SFcollapse1D

SFcollapse1D.o: $(SRC_DIR)/SFcollapse1D.cpp $(SFcollapse1D_headers) obj_dir
	$(CXX) $(CXXFLAGS) -o $(OBJ_DIR)/SFcollapse1D.o -c $(SRC_DIR)/SFcollapse1D.cpp

grid.o: $(SRC_DIR)/grid.cpp $(grid_headers) obj_dir
	$(CXX) $(CXXFLAGS) -o $(OBJ_DIR)/grid.o -c $(SRC_DIR)/grid.cpp

gridfunction.o: $(SRC_DIR)/gridfunction.cpp $(gridfunction_headers) obj_dir
	$(CXX) $(CXXFLAGS) -o $(OBJ_DIR)/gridfunction.o -c $(SRC_DIR)/gridfunction.cpp

evolution.o: $(SRC_DIR)/evolution.cpp $(evolution_headers) obj_dir
	$(CXX) $(CXXFLAGS) -o $(OBJ_DIR)/evolution.o -c $(SRC_DIR)/evolution.cpp

utilities.o: $(SRC_DIR)/utilities.cpp $(utilities_headers) obj_dir
	$(CXX) $(CXXFLAGS) -o $(OBJ_DIR)/utilities.o -c $(SRC_DIR)/utilities.cpp

obj_dir:
	mkdir -p obj

out_dir:
	mkdir -p out

clean:
	rm -rf $(OBJ_PATHS) SFcollapse1D

realclean:
	rm -rf $(OBJ_PATHS) SFcollapse1D *.txt out/*.dat animations/*.gif doc/*.pdf out/ obj/
