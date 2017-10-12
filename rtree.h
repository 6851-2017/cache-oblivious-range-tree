/**
 * rtree.h
 *
 * Header for orthogonal 2D range searching with static trees
 */

#pragma once

#include <stdlib.h>

typedef struct {
  double x, y;
} point_t;

typedef struct {
  point_t *pts;
  size_t num;
} result_t;

typedef struct bst_s {
  struct bst_s *left;
  struct bst_s *right;
  struct bst_s *parent;
  point_t *min;
  point_t *max;
} bst_t;

typedef struct {
  point_t *pts;
  bst_t *bst;
  size_t num;
} tree_t;


/**
 * tree_new
 *
 * args:
 *   n (size_t): number of points in data
 *   data (point_t*): the n points to statically store, sorted by x then y
 *
 * returns a tree_t, the range tree representing data
 */
tree_t tree_new(size_t n, point_t *data);

/**
 * tree_query2
 *
 * 2-sided tree range query
 *
 * args:
 *   tree (tree_t*): the tree_t to search
 *   p (point_t*): the point representing the 2-sided query
 *
 * returns a result_t, which contains the set of points in tree with both
 *   coordinates <= the coordinates of p
 */
result_t tree_query(tree_t *tree, point_t *p);
