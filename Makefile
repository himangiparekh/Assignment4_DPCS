# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17

# Directories
SRC_DIR = .
OBJ_DIR = obj
BIN_DIR = bin

# Source files
SOURCES = $(SRC_DIR)/polyhedron.cpp $(SRC_DIR)/main.cpp \
          $(SRC_DIR)/input.cpp $(SRC_DIR)/transformations.cpp \
          $(SRC_DIR)/projections.cpp

# Object files
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SOURCES))

# Headers (for dependency tracking)
HEADERS = $(SRC_DIR)/polyhedron.h $(SRC_DIR)/input.h \
          $(SRC_DIR)/transformations.h $(SRC_DIR)/projections.h

# Final executable
EXECUTABLE = $(BIN_DIR)/polyhedron_program

# Default target - build the system
all: $(EXECUTABLE)

# Rule to link the final executable
$(EXECUTABLE): $(OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(EXECUTABLE)
	@echo "Build complete: $(EXECUTABLE)"

# Rule to compile source files into object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean compiled files
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Clean and rebuild
rebuild: clean all

# Phony targets (not actual files)
.PHONY: all clean rebuild
