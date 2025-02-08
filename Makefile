CC := gcc
CFLAGS := -Wall -Wextra -Werror -std=gnu11 -MMD -MP -march=native
LFLAGS :=
OFLAGS := -Ofast
EXEC := whirlpool_gd
BUILDDIR := ./build
SRCDIR := ./src
SRC_PAT := $(SRCDIR)/%.c
OBJ_PAT := $(BUILDDIR)/%.o

SRCS := $(shell find $(SOURCEDIR) -name '*.c')
OBJS := $(subst $(SRCDIR),$(BUILDDIR),$(SRCS:.c=.o))
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

build/%.o: src/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(OFLAGS) -o $@ -c $<

clean:
	rm -rf ./*.exe $(BUILDDIR)/*
	
-include $(DEPS)