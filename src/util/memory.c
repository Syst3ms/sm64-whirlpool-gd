#include <stdio.h>

#include "./memory.h"

void * _mm_realloc(void *aligned_ptr, size_t size, size_t align)
{
  void *malloc_ptr = ((void**) aligned_ptr)[-1];

  malloc_ptr = realloc(malloc_ptr, size + align);
  if (!malloc_ptr)
    return ((void *) 0);

  /* Align  We have at least sizeof (void *) space below malloc'd ptr. */
  aligned_ptr = (void *) (((size_t) malloc_ptr + align)
			  & ~((size_t) (align) - 1));

  /* Store the original pointer just before p.  */	
  ((void **) aligned_ptr) [-1] = malloc_ptr;

  return aligned_ptr;
}

struct history init_history(size_t initial_size) {
    struct history mem = {
        initial_size,
        0,
        malloc(initial_size * POINTS * sizeof(v2d))
    };

    if (mem.pts == NULL) {
        puts("Could not allocate arrays for history!");
        exit(1);
    }
    
    return mem;
}

void store_into_history(struct data *d, struct history *mem) {
    size_t offset = POINTS * mem->next;
    for (size_t i = 0; i < POINTS; i++) {
        mem->pts[offset + i] = d->points[i].pos;
    }

    if (++mem->next >= mem->size) {
        mem->size += 1000;
        printf("Memory exceeded, reallocing to %d\n", mem->size);
        mem->pts = realloc(mem->pts, mem->size * POINTS * sizeof(v2d));
        if (mem->pts == NULL) {
            puts("Could not realloc!");
            exit(1);
        }
    }
}

void free_history(struct history mem) {
    free(mem.pts);
}