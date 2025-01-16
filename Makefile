CC := gcc
CFLAGS := -Wall -Wextra -Werror -std=gnu11 -MMD -MP -march=native
LFLAGS :=
OFLAGS := -Ofast
EXEC := whirlpool_gd
BUILDDIR := ./build
SRCDIR := ./src
SRC_PAT := $(SRCDIR)/%.c
OBJ_PAT := $(BUILDDIR)/%.o

SRCS := $(wildcard $(SRCDIR)/*.c)
OBJS := $(SRCS:$(SRC_PAT)=$(OBJ_PAT))
DEPS := $(OBJS:.o=.d)

.PHONY: run_and_plot plot run build clean assembly

run_and_plot: run plot

plot:
	@python plot_c_result.py

run: build
	@$(EXEC).exe

build: $(OBJS)
	$(CC) $(LFLAGS) -o $(EXEC) $(OBJS)

debug: OFLAGS = -g
debug: build

$(OBJ_PAT): $(SRC_PAT)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(OFLAGS) -o $@ -c $<

clean:
	rm -f *.exe $(BUILDDIR)/*
	
-include $(DEPS)