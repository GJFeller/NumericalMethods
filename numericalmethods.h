#ifndef NUMERICALMETHODS_H
#define NUMERICALMETHODS_H

#include <stdio.h>
#include <stdlib.h>

class NumericalMethods
{
public:

    enum PivotType{ NORMAL, PARTIAL };

    NumericalMethods();

    void quadraticGaussJordan(int n, float* A, float* b, float* x, PivotType pivotType);

    void LUFatoration(int n, float* A, float *L, float *U);

private:
    void swapLines(int n, float* Ab, int lin1, int lin2);
    void print_mat(int r, int c, float* m);

    void normalPivot(int n, float* Ab);
    void partialPivot(int n, float* Ab);



};

#endif // NUMERICALMETHODS_H
