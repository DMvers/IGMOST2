#include "project.h"
#include "parameter.h"

extern Par par;

/* This file contains a number of additional functions that are used during the simulation. */

double max(double a, double b)
{
  /* returns whichever is the highest number of a and b */
  double c = a > b ? a : b;
  return c;
}

double min(double a, double b)
{
  /* returns whichever is the lowest number of a and b */
  double c = a < b ? a : b;
  return c;
}

int idu = -1;
double Uniform()
{
  /* Function that returns a random number between 0 and 1.
   * Implements Knuth's substrative method, see Numerical Recipes.
   */

  static int inext, inextp;
  static long ma[56];
  static int iff = 0;
  long mj, mk;
  int g, ii, k;
  if (idu < 0 || iff == 0)
  {
    iff = 1;
    mj = abs(par.randseed - (idu < 0 ? -idu : idu));
    mj %= MBIG;
    ma[55] = mj;
    mk = 1;
    for (g = 1; g <= 54; g++)
    {
      ii = (21 * g) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if (mk < MZ)
        mk += MBIG;
      mj = ma[ii];
    }
    for (k = 1; k <= 4; k++)
      for (g = 1; g <= 55; g++)
      {
        ma[g] -= ma[1 + (g + 30) % 55];
        if (ma[g] < MZ)
          ma[g] += MBIG;
      }
    inext = 0;
    inextp = 31;
    idu = 1;
  }
  if (++inext == 56)
    inext = 1;
  if (++inextp == 56)
    inextp = 1;
  mj = ma[inext] - ma[inextp];
  if (mj < MZ)
    mj += MBIG;
  ma[inext] = mj;
  return mj * FAC;
}

int RandNum(int n)
{
  /* Function that returns a random integer between 1 and n, including 1 and n.
   * Uses Uniform() internally.
   */
  int temp;
  if (n <= 0)
    return (0);
  else
  {
    temp = n * Uniform() + 1;
    return (temp);
  }
}

void Shufflevector(int *vector, int size)
{
  /* Function that randomly shuffles a vector or an array.
   * Uses RandNum(n) and thus Uniform() internally.
   */
  int xa, xb, pos = 0;
  int tussen;
  for (xa = 0; xa < size; xa++)
  {
    xb = pos + RandNum(size) - 1;
    tussen = vector[xa];
    vector[xa] = vector[xb];
    vector[xb] = tussen;
    pos++;
    size--;
  }
}
