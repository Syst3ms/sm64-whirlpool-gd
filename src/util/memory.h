#pragma once
#include <stddef.h>
#include <mm_malloc.h>

#include "../types.h"

void * _mm_realloc(void *aligned_ptr, size_t size, size_t align);

struct history init_history(size_t initial_size);
void store_into_history(struct data *d, struct history *mem);
void free_history(struct history mem);