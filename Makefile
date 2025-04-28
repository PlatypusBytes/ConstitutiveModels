# Detect OS
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
    SO_EXT = .so
    RM = rm -rf
    MKDIR = mkdir -p
else ifeq ($(OS),Windows_NT)  # Git Bash, MSYS, or native make on Windows
    SO_EXT = .dll
    RM = del /Q /S
    MKDIR = mkdir
endif

CC = gcc
CFLAGS = -fPIC -Wall -Wextra -O2 -g
LDFLAGS = -static-libgcc -static-libstdc++ -shared

# Output directories
BUILD_DIR = build_C
LIB_DIR = $(BUILD_DIR)/lib

# Source directories
C_MODELS_DIR = c_models
ELASTIC_LAWS_DIR = $(C_MODELS_DIR)/elastic_laws
LINEAR_ELASTIC_TENS_DIR = $(C_MODELS_DIR)/linear_elastic_with_vertical_tens_cutoff
MATSUOKA_NAKAI_DIR = $(C_MODELS_DIR)/matsuoka_nakai
YIELD_SURFACES_DIR = $(C_MODELS_DIR)/yield_surfaces

# Include directories
INCLUDES = -I$(C_MODELS_DIR) -I$(ELASTIC_LAWS_DIR) -I$(YIELD_SURFACES_DIR)

# Common object files
COMMON_OBJS = \
	$(BUILD_DIR)/globals.o \
	$(BUILD_DIR)/utils.o \
	$(BUILD_DIR)/stress_utils.o \
	$(BUILD_DIR)/hookes_law.o

# Target shared libraries
LIBS = \
	$(LIB_DIR)/linear_elastic_with_vertical_tens_cutoff$(SO_EXT) \
	$(LIB_DIR)/matsuoka_nakai$(SO_EXT)

# Default target
all: directories $(LIBS)

# Create necessary directories
directories:
	$(MKDIR) $(BUILD_DIR)
	$(MKDIR) $(LIB_DIR)

# Common object files compilation
$(BUILD_DIR)/globals.o: $(C_MODELS_DIR)/globals.c $(C_MODELS_DIR)/globals.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/utils.o: $(C_MODELS_DIR)/utils.c $(C_MODELS_DIR)/utils.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/stress_utils.o: $(C_MODELS_DIR)/stress_utils.c $(C_MODELS_DIR)/stress_utils.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/hookes_law.o: $(ELASTIC_LAWS_DIR)/hookes_law.c $(ELASTIC_LAWS_DIR)/hookes_law.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Matsuoka-Nakai yield surface
$(BUILD_DIR)/matsuoka_nakai_surface.o: $(YIELD_SURFACES_DIR)/matsuoka_nakai_surface.c $(YIELD_SURFACES_DIR)/matsuoka_nakai_surface.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Linear elastic with vertical tension cutoff model
$(BUILD_DIR)/linear_elastic_with_vertical_tens_cutoff.o: $(LINEAR_ELASTIC_TENS_DIR)/linear_elastic_with_vertical_tens_cutoff.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Matsuoka-Nakai model
$(BUILD_DIR)/matsuoka_nakai.o: $(MATSUOKA_NAKAI_DIR)/matsuoka_nakai.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Shared libraries
$(LIB_DIR)/linear_elastic_with_vertical_tens_cutoff$(SO_EXT): $(COMMON_OBJS) $(BUILD_DIR)/linear_elastic_with_vertical_tens_cutoff.o
	$(CC) $(LDFLAGS) -o $@ $^

$(LIB_DIR)/matsuoka_nakai$(SO_EXT): $(COMMON_OBJS) $(BUILD_DIR)/matsuoka_nakai_surface.o $(BUILD_DIR)/matsuoka_nakai.o
	$(CC) $(LDFLAGS) -o $@ $^

# Clean target
clean:
	$(RM) $(BUILD_DIR)

# Help target
help:
	@echo "Available targets:"
	@echo "  all       - Build all shared libraries (default)"
	@echo "  clean     - Remove all build artifacts"
	@echo "  help      - Show this help message"

.PHONY: all directories clean help