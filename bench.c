#include <stdio.h>

#include "rtree.h"

int main() {
  point_t pts[8] = {
    (point_t){ .x = -4, .y = 1 },
    (point_t){ .x = -3, .y = 0 },
    (point_t){ .x = -2, .y = 1 },
    (point_t){ .x = -1, .y = 0 },
    (point_t){ .x = 1, .y = -1 },
    (point_t){ .x = 2, .y = 1 },
    (point_t){ .x = 3, .y = -1 },
    (point_t){ .x = 4, .y = 1 },
  };

  tree_t tree = tree_new(8, pts);

  for (size_t i = 0; i < tree.npts; i++) {
    printf("tree.pts[%zu] = (%lf, %lf)\n", i, tree.pts[i].x, tree.pts[i].y);
  }

  return 0;
}
