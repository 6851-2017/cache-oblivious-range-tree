#include <stdio.h>

#include "rtree.h"

int main() {
  point_t pts[4] = {
    (point_t){ .x = -2, .y = 1 },
    (point_t){ .x = -1, .y = 0 },
    (point_t){ .x = 1, .y = -1 },
    (point_t){ .x = 2, .y = 1 }
  };

  tree_t tree = tree_new(4, pts);

  for (size_t i = 0; i < 4; i++) {
    printf("tree.pts[%zu] = (%lf, %lf)\n", i, tree.pts[i].x, tree.pts[i].y);
  }

  return 0;
}
