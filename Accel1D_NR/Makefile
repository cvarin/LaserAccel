################################################################################
BIN  = Accel1D.exe
BUILD_DIR = ./build
DATA_DIR  = ./output
SRC_DIR   = ./src ./src/nr

################### Fetch names ################################################
VPATH = $(SRC_DIR)
SRC   = $(foreach DIR, $(SRC_DIR), $(wildcard $(DIR)/*.c))
HEADS = $(foreach DIR, $(SRC_DIR), $(wildcard $(DIR)/*.h))
INCLUDES = $(addprefix -I./, $(SRC_DIR) . )
xnames = $(notdir $(SRC))
NAMES = $(subst .c,,$(xnames))
OBJS  = $(addsuffix .o, $(NAMES))
LINKS = $(addprefix $(BUILD_DIR)/,$(OBJS))
DATA  = $(wildcard $(DATA_DIR)/*)

################### Commands ###################################################
CPP   = g++
FLAGS = -fPIC -O3 -Wall
# définir la macro NDEBUG pour supprimer les "assert()".
RM    = rm -fr

################### Compiling and linking ######################################
all: $(BIN)

$(BIN): $(LINKS) $(HEADS)
	@echo "LINKING: $(CPP) $(OBJS) -o $(BIN)"
	@$(CPP) $(LINKS) -o $(BIN)

$(BUILD_DIR)/%.o: %.c $(HEADS)
	@echo "COMPILING: $(CPP) $(FLAGS) -c $< -o $@"
	@$(CPP) $(FLAGS) $(INCLUDES) -c $< -o $@

################### Run ########################################################
.PHONY: r run
r: run
run: $(BIN)
	./$(BIN)

################### Clean ######################################################
.PHONY: clean c
c: clean
clean:
	${RM} $(LINKS) $(BIN) $(DATA)

.PHONY: clean_data cd
cd: clean_data
clean_data:
	$(RM) $(DATA)

################### END of FILE ################################################

