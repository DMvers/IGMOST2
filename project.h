#ifndef FBA_SPATIAL_
#define FBA_SPATIAL_

#include <gsl/gsl_randist.h>
#include <glpk.h>
#include <vector>

#define NBS 9
#define LYS 1
#define MBIG 1000000000
#define MZ 0
#define FAC (1.0 / MBIG)

double max(double a, double b); // returns whichever is the highest number of a and b
double min(double a, double b); // returns whichever is the lowest number of a and b

double Uniform();                          // Function that returns a random number between 0 and 1.
int RandNum(int n);                        // Function that returns a random integer between 1 and n.
void Shufflevector(int *vector, int size); // Function that randomly shuffles a vector or an array.

#endif
