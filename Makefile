MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile)
MPI_LINK_FLAGS = $(shell mpicc --showme:link)

TARGET_DIR 	= ./bin
TARGET		= mandelbrot
INCDIR		= ./src/include
SRCDIR		= ./src
OBJDIR		= ./src/obj
CC			= gcc
CFLAGS		= -Wall -Ofast $(MPI_COMPILE_FLAGS) -fopenmp-simd -march=native 
INCLUDES	= -I$(INCDIR)
LINKER		= gcc
LIBS		= -lm -lpng -fopenmp -lavformat -lavcodec -lavutil -lpthread -lswscale
LFLAGS		= -Wall -I. $(LIBS) $(MPI_LINK_FLAGS)

SRCS := $(wildcard $(SRCDIR)/*.c)
DEPS := $(patsubst $(SRCDIR)/%.c, $(INCDIR)/%.h, $(SRCS))
OBJS := $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SRCS))

dir_guard=@mkdir -p $(@D)

all: $(TARGET_DIR)/$(TARGET)
	@echo "Compilation of "$<" finished!"

$(TARGET_DIR)/$(TARGET): $(OBJS)
	$(dir_guard)
	@$(LINKER) $(OBJS) $(LFLAGS) -o $@
	@echo "Linking completed"

$(OBJS): $(OBJDIR)/%.o : $(SRCDIR)/%.c $(DEPS)
	$(dir_guard)
	@$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@
	@echo "Compiled "$<" successfully"

.PHONY: clean
clean:
	@rm -f $(OBJS)
	@rm -rf $(OBJDIR)
	@echo "Object files removed"
	@rm -f $(TARGET)
	@rm -rf $(TARGET_DIR)
	@echo "Executable removed"
	@echo "Cleaning completed!"

print-%  : ; @echo $* = $($*)
