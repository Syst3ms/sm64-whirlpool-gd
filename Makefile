CC = gcc
CFLAGS = -Wall -Wextra -Werror -MMD -MP
LFLAGS =
OFLAGS = -O2
SRCS = main.c util.c math_funcs.c
EXEC = whirlpool_gd
BUILDDIR = ./build

OBJS = $(SRCS:%.c=$(BUILDDIR)/%.o)
DEPS = $(OBJS:.o=.d)

.PHONY: run_and_plot plot run build clean

run_and_plot: run plot

plot:
	@python plot_c_result.py
	
run: build
	@./$(EXEC).exe

build: $(OBJS)
	$(CC) $(LFLAGS) $(OFLAGS) -o $(EXEC) $(OBJS)

debug: OFLAGS = -g
debug: build

$(BUILDDIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(OFLAGS) -o $@ -c $<

clean:
	rm -f *.exe $(BUILDDIR)/*
	
-include $(DEPS)