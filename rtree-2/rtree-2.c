/**
 * rtree-2.c
 *
 * Static cache-oblivious range tree for 2-sided orthogonal range searching.
 */

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include <stdio.h>

#include "rtree.h"

#define max(a, b) (a > b ? a : b)
#define min(a, b) (a < b ? a : b)
#define swap(a, b) {a ^= b; b ^= a; a ^= b;}

#define ALPHA 2
// S_CUTOFF determines the max size of the last S
#define S_CUTOFF 2

/**
 * get_sparse_y
 *
 * args:
 *   n (size_t): number of points
 *   ptr (point_t*): list of points, sorted by x then y
 *
 * returns the maximum y for which (*, <=y) is sparse in ptr
 */
double get_sparse_y(size_t n, point_t *ptr) {
  const size_t sparse_n = (n - 1) / ALPHA;
  point_t s_pts[sparse_n];  // Sparse points
  double max_y = -INFINITY;
  double max_max_y = INFINITY;
  for (point_t *p_itr = ptr, *s_itr = s_pts; p_itr < ptr + n; p_itr++) {
    // We already pruned points with y values <= to this one
    if (s_itr->y >= max_max_y) {
      continue;
    }

    // Just take the first sparse_n points...
    if (s_itr < s_pts + sparse_n) {
      *s_itr = *p_itr;
      max_y = max(max_y, p_itr->y);
      s_itr++;
      continue;
    }

    // Now check if the current point will move max_y down
    if (p_itr->y > max_y) {
      max_max_y = min(max_max_y, p_itr->y);
      continue;
    }

    // Compute the new max y
    double prev_max_y = max_y;
    double min_x = INFINITY;  // Find min x corresponding to max y
    max_y = -INFINITY;
    for (point_t *itr = s_pts; itr < s_itr; itr++) {
      if (itr->y < prev_max_y) {
        if (itr->y >= max_y) {
          max_y = itr->y;
          min_x = min(min_x, itr->x);
        }
      }
    }
    max_max_y = min(max_max_y, prev_max_y);

    // Drop elements based on new max y
    point_t *next_s_itr = s_pts;
    for (point_t *itr = s_pts; itr < s_itr; itr++) {
      if (itr->x <= min_x && itr->y <= max_y) {
        *next_s_itr = *itr;
        next_s_itr++;
      }
    }
    s_itr = next_s_itr;
  }
  return max_y;
}

/**
 * get_sparse_x
 *
 * args:
 *   n (size_t): number of points
 *   ptr (point_t*): list of points, sorted by x then y
 *
 * returns the maximum x for which (<=x, *) is sparse in ptr
 */
double get_sparse_x(size_t n, point_t *ptr) {
  double y = get_sparse_y(n, ptr);
  double x = ptr[0].x - 1;
  double next_x = ptr[0].x;
  size_t count = 0;
  size_t next_count = 1;

  // Iterate through by x
  for (point_t *itr = ptr + 1; itr < ptr + n && next_count * ALPHA < n; itr++) {
    if (itr->x > next_x) {
      x = next_x;
      count = next_count;
    }

    next_count += (itr->y <= y);
    next_x = itr->x;
  }

  if (next_count * ALPHA < n) {
    x = next_x;
  }

  return x;
}

/**
 * get_p
 *
 * args:
 *   n (size_t): number of points in s
 *   s (point_t*): the parent S_{i-1} from which to generate the P_{i-1}, sorted
 *                 by x then y
 *   p (point_t*): the p buffer, should have at least n / ALPHA slots, sorted by
 *                 x then y
 *
 * returns the number of elements in p
 */
size_t get_p(size_t n, point_t *s, point_t *p) {
  double max_x = get_sparse_x(n, s);
  size_t i;
  for (i = 0; i < n && s[i].x <= max_x; i++) {
    p[i] = s[i];
  }
  return i;
}

/**
 * get_s
 *
 * args:
 *   n (size_t): number of points in s0
 *   s0 (point_t*): the parent S_{i-1} from which to generate the S_i, sorted by
 *                  x then y
 *   s1 (point_t*): the S_i buffer, which should have at least n slots, sorted
 *                  by x then y
 *
 * returns the number of elements in s1
 */
size_t get_s(size_t n, point_t *s0, point_t *s1) {
  double max_x = get_sparse_x(n, s0);
  double max_y = get_sparse_y(n, s0);
  point_t *itr1 = s1;
  for (point_t *itr0 = s0; itr0 < s0 + n; itr0++) {
    if (itr0->x > max_x || itr0->y <= max_y) {
      *itr1 = *itr0;
      itr1++;
    }
  }

  return (size_t)(itr1 - s1);
}

tree_t tree_new(size_t n, point_t *data) {
  // Bit trick to round up to next power of 2
  size_t bst_size = n;
  bst_size--;
  bst_size |= bst_size >> 1;
  bst_size |= bst_size >> 2;
  bst_size |= bst_size >> 4;
  bst_size |= bst_size >> 8;
  bst_size |= bst_size >> 16;
  bst_size |= bst_size >> 32;
  bst_size++;

  // ALPHA / (ALPHA - 1) + 1 to account for rounding errors
  size_t pts_size = (ALPHA / (ALPHA - 1) + 1) * n;

  tree_t tree = (tree_t) {
    .pts = (point_t*)malloc(pts_size * (sizeof *tree.pts)),
    .bst = (bst_t*)malloc(bst_size * (sizeof *tree.bst)),
    .npts = pts_size,
    .nbst = bst_size,
  };

  // Prepare for vEB order construction
  const size_t MAX_S_N = n;
  point_t s_buf0[MAX_S_N];
  point_t s_buf1[MAX_S_N];
  uintptr_t s0 = (uintptr_t)s_buf0;  // Alias so we can swap buffers
  uintptr_t s1 = (uintptr_t)s_buf1;  // Alias so we can swap buffers
  memcpy((point_t*)s0, data, sizeof s_buf0);

  // Construct in vEB order for cache-obliviousness
  size_t n_rem = n;  // # elements left to put in array
  point_t *pptr = tree.pts;  // Pointer to next location in array to fill
  volatile int counter = 0;
  while (n_rem > S_CUTOFF) {
    ++counter;
    point_t *p0 = pptr;  // p = P_{i-1}
    pptr += get_p(n_rem, (point_t*)s0, p0);
    n_rem = get_s(n_rem, (point_t*)s0, (point_t*)s1);
    swap(s0, s1);
  }

  // Tail case of the last few elements, i.e. the final S_i
  if (n_rem) {
    memcpy(pptr, (point_t*)s0, n_rem * (sizeof *pptr));
  }

  return tree;
}

result_t tree_query(tree_t *tree, point_t *p) {
}
