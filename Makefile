CC = gcc
CFLAGS = -Wall -Wextra -Werror -O2 -MMD -MP
LFLAGS =
SRCS = main.c util.c
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
	@$(CC) $(LFLAGS) -o $(EXEC) $(OBJS)

$(BUILDDIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -f *.exe $(BUILDDIR)
	
-include $(DEPS)