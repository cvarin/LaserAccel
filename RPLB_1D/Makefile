################################################################################
BIN  = prog.exe
BUILD_DIR = ./build
DATA_DIR  = ./output
SRC_DIR   = ./src
################### Fetch names ################################################
SRC   = $(wildcard $(SRC_DIR)/*.c)
xnames = $(subst $(SRC_DIR)/,,$(SRC))
NAMES = $(subst .c,,$(xnames))
HEADS = $(wildcard $(SRC_DIR)/*.h)
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
	$(CPP) $(LINKS) -o $(BIN) $(LIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(HEADS)
	$(CPP) $(FLAGS) -c $< -o $@

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

