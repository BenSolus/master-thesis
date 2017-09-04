/* ------------------------------------------------------------
 * This is the file "main.c" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      main.c
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#include "basic.h"
#include "singquad2d.h"

#include <iostream>

const uint px[][3] = { { 0, 1, 2 }, // x = ( a, b, c ), y = ( a, b, c ) or
                                    // x = ( a, b, c ), y = ( d, e, f )
                       /* Identical vertex cases */
                       /* Both contain a */
                       { 0, 1, 2 }, // x = ( a, b, c ), y = ( a, d, e )
                       { 0, 1, 2 }, // x = ( a, b, c ), y = ( d, a, e )
                       { 0, 1, 2 }, // x = ( a, b, c ), y = ( d, e, a )
                       /* Both contain b */
                       { 1, 0, 2 }, // x = ( a, b, c ), y = ( b, d, e )
                       { 1, 0, 2 }, // x = ( a, b, c ), y = ( d, b, e )
                       { 1, 0, 2 }, // x = ( a, b, c ), y = ( d, e, b )
                       /* Both contain c */
                       { 2, 0, 1 }, // x = ( a, b, c ), y = ( c, d, e )
                       { 2, 0, 1 }, // x = ( a, b, c ), y = ( d, c, e )
                       { 2, 0, 1 }, // x = ( a, b, c ), y = ( d, e, c )
                       /* Identical edge cases */
                       /* Both contain a and b */
                       { 0, 1, 2 }, // x = ( a, b, c ), y = ( a, b, d )
                       { 0, 1, 2 }, // x = ( a, b, c ), y = ( b, a, d )
                       { 0, 1, 2 }, // x = ( a, b, c ), y = ( b, d, a )
                       { 0, 1, 2 }, // x = ( a, b, c ), y = ( a, d, b )
                       { 0, 1, 2 }, // x = ( a, b, c ), y = ( d, a, b )
                       { 0, 1, 2 }, // x = ( a, b, c ), y = ( d, b, a )
                       /* Both contain b and c */
                       { 1, 2, 0 }, // x = ( a, b, c ), y = ( b, c, d )
                       { 1, 2, 0 }, // x = ( a, b, c ), y = ( c, b, d )
                       { 1, 2, 0 }, // x = ( a, b, c ), y = ( c, d, b )
                       { 1, 2, 0 }, // x = ( a, b, c ), y = ( b, d, c )
                       { 1, 2, 0 }, // x = ( a, b, c ), y = ( d, b, c )
                       { 1, 2, 0 }, // x = ( a, b, c ), y = ( d, c, b )
                       /* Both contain a and c */
                       { 0, 2, 1 }, // x = ( a, b, c ), y = ( a, c, d )
                       { 0, 2, 1 }, // x = ( a, b, c ), y = ( c, a, d )
                       { 0, 2, 1 }, // x = ( a, b, c ), y = ( c, d, a )
                       { 0, 2, 1 }, // x = ( a, b, c ), y = ( a, d, c )
                       { 0, 2, 1 }, // x = ( a, b, c ), y = ( d, a, c )
                       { 0, 2, 1 }, // x = ( a, b, c ), y = ( d, c, a )
};

const uint py[][3] = { { 0, 1, 2 }, // x = ( a, b, c ), y = ( d, e, f )
                                    // x = ( a, b, c ), y = ( d, e, f )
                       /* Identical vertex cases */
                       /* Both contain a */
                       { 0, 1, 2 }, // x = ( a, b, c ), y = ( a, b, d )
                       { 1, 0, 2 }, // x = ( a, b, c ), y = ( d, a, e )
                       { 2, 0, 1 }, // x = ( a, b, c ), y = ( d, a, e )
                       /* Both contain b */
                       { 0, 1, 2 }, // x = ( a, b, c ), y = ( b, d, e )
                       { 1, 0, 2 }, // x = ( a, b, c ), y = ( b, d, e )
                       { 2, 0, 1 }, // x = ( a, b, c ), y = ( d, e, b )
                       /* Both contain c */
                       { 0, 1, 2 }, // x = ( a, b, c ), y = ( d, e, b )
                       { 1, 0, 2 }, // x = ( a, b, c ), y = ( d, c, e )
                       { 2, 0, 1 }, // x = ( a, b, c ), y = ( d, e, c )
                       /* Identical edge cases */
                       /* Both contain a and b */
                       { 0, 1, 2 }, // x = ( a, b, c ), y = ( a, b, d )
                       { 1, 0, 2 }, // x = ( a, b, c ), y = ( b, a, d )
                       { 2, 0, 1 }, // x = ( a, b, c ), y = ( b, d, a )
                       { 0, 2, 1 }, // x = ( a, b, c ), y = ( a, d, b )
                       { 1, 2, 0 }, // x = ( a, b, c ), y = ( d, a, b )
                       { 2, 1, 0 }, // x = ( a, b, c ), y = ( d, b, a )
                       /* Both contain b and c */
                       { 0, 1, 2 }, // x = ( a, b, c ), y = ( b, c, d )
                       { 1, 0, 2 }, // x = ( a, b, c ), y = ( c, b, d )
                       { 2, 0, 1 }, // x = ( a, b, c ), y = ( c, d, b )
                       { 0, 2, 1 }, // x = ( a, b, c ), y = ( b, d, c )
                       { 1, 2, 0 }, // x = ( a, b, c ), y = ( d, b, c )
                       { 2, 1, 0 }, // x = ( a, b, c ), y = ( d, c, b )
                       /* Both contain a and c */
                       { 0, 1, 2 }, // x = ( a, b, c ), y = ( a, c, d )
                       { 1, 0, 2 }, // x = ( a, b, c ), y = ( c, a, d )
                       { 2, 0, 1 }, // x = ( a, b, c ), y = ( c, d, a )
                       { 0, 2, 1 }, // x = ( a, b, c ), y = ( a, d, c )
                       { 1, 2, 0 }, // x = ( a, b, c ), y = ( d, a, c )
                       { 2, 1, 0 }, // x = ( a, b, c ), y = ( d, c, a )
                     };

int
main(int argc, char** argv)
{
  (void) argc;
  (void) argv;

  pcsingquad2d sq = build_singquad2d(NULL, 2, 4);

  const uint x[3] = { 0, 1, 2 };
  const uint y[3] = { 3, 2, 1 };

  const uint p = (uint) (x[0] == y[0]) + (x[0] == y[1]) + (x[0] == y[2]) +
                        (x[1] == y[0]) + (x[1] == y[1]) + (x[1] == y[2]) +
                        (x[2] == y[0]) + (x[2] == y[1]) + (x[2] == y[2]);

  uint nq;
  real base;

  uint xp[3] = { 0, 0, 0 };
  uint yp[3] = { 0, 0, 0 };

  real *xq, *yq, *w;

  select_quadrature_singquad2d(sq, x, y, xp, yp, &xq, &yq, &w, &nq, &base);

  std::cout << "Number of identical vertices: " << p << "\n";

  /* Offsets: 0 for p = 0 or p = 3 and 1 + (p - 1) * 9 for p = 1 or p = 2 */
  uint i = ((p == 1) || (p == 2)) * (1 + (p - 1) * 9);

  /* i == 0: dist or id case */

  /* i == 1 ... 9: identical vertex case */
  /* i == 1: (p == 1) && (x[0] == y[0]) */
  i += ((p == 1) && (x[0] == y[1]));     // i == 2
  i += ((p == 1) && (x[0] == y[2])) * 2; // i == 3
  i += ((p == 1) && (x[1] == y[0])) * 3; // i == 4
  i += ((p == 1) && (x[1] == y[1])) * 4; // i == 5
  i += ((p == 1) && (x[1] == y[2])) * 5; // i == 6
  i += ((p == 1) && (x[2] == y[0])) * 6; // i == 7
  i += ((p == 1) && (x[2] == y[1])) * 7; // i == 8
  i += ((p == 1) && (x[2] == y[2])) * 8; // i == 9

  /* i == 10 ... 18: identical edge case */
  /* Both contain a and b */
  /* i == 10: (p == 2) && (x[0] == y[0]) && (x[1] == y[1]) */
  i += ((p == 2) && (x[0] == y[1]) && (x[1] == y[0]));      // i == 11
  i += ((p == 2) && (x[0] == y[2]) && (x[1] == y[0])) *  2; // i == 11
  i += ((p == 2) && (x[0] == y[0]) && (x[1] == y[2])) *  3; // i == 13
  i += ((p == 2) && (x[0] == y[1]) && (x[1] == y[2])) *  4; // i == 14
  i += ((p == 2) && (x[0] == y[2]) && (x[1] == y[1])) *  5; // i == 15

  i += ((p == 2) && (x[1] == y[0]) && (x[2] == y[1])) *  6; // i == 16
  i += ((p == 2) && (x[1] == y[1]) && (x[2] == y[0])) *  7; // i == 17
  i += ((p == 2) && (x[1] == y[2]) && (x[2] == y[0])) *  8; // i == 18
  i += ((p == 2) && (x[1] == y[0]) && (x[2] == y[2])) *  9; // i == 19
  i += ((p == 2) && (x[1] == y[1]) && (x[2] == y[2])) * 10; // i == 20
  i += ((p == 2) && (x[1] == y[2]) && (x[2] == y[1])) * 11; // i == 21

  i += ((p == 2) && (x[0] == y[0]) && (x[2] == y[1])) * 12; // i == 22
  i += ((p == 2) && (x[0] == y[1]) && (x[2] == y[0])) * 13; // i == 23
  i += ((p == 2) && (x[0] == y[2]) && (x[2] == y[0])) * 14; // i == 24
  i += ((p == 2) && (x[0] == y[0]) && (x[2] == y[2])) * 15; // i == 25
  i += ((p == 2) && (x[0] == y[1]) && (x[2] == y[2])) * 16; // i == 26
  i += ((p == 2) && (x[0] == y[1]) && (x[2] == y[2])) * 16; // i == 26

  std::cout << "Permutation number " << i << "\n";
  std::cout << "Permutation indices of x: ";

  for(uint i(0); i < 3; ++i)
    std::cout << xp[i] << " ";

  std::cout << "\nPermutation indices of y: ";

  for(uint i(0); i < 3; ++i)
    std::cout << yp[i] << " ";

  std::cout << "\n";
}