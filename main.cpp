#include <QApplication>
#include <QDebug>
#include <stdio.h>
#include <stdlib.h>
#include "numericalmethods.h"

int main() {

    NumericalMethods test;
    float A[4][4] = {
        {3, 2, 0, 1},
        {9, 8, -3, 4},
        {-6, 4, -8, 0},
        {3, -8, 3, -4}
    };

    /*A[0][0] = 3;
    A[0][1] = 2;
    A[0][2] = 0;
    A[0][3] = 1;
    A[1][0] = 9;
    A[1][1] = 8;
    A[1][2] = -3;
    A[1][3] = 4;
    A[2][0] = -6;
    A[2][1] = 4;
    A[2][2] = -8;
    A[2][3] = 0;
    A[3][0] = 3;
    A[3][1] = -8;
    A[3][2] = 3;
    A[3][3] = -4;*/

    /*printf("A:\n");
    for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                    printf("\t%f", A[i][j]);
            }
            printf("\n");
    }*/
    float* Amatrix = (float *) malloc(sizeof (float) * 4 * 4);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            Amatrix[i*4+j] = A[i][j];
        }
    }
    /*printf("Amatrix:\n");
    for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                    printf("\t%f", Amatrix[i*4+j]);
            }
            printf("\n");
    }*/


    float b[4] = {3, 6, -16, 18};
    float* x = (float *) malloc(sizeof (float) * 4);

    //test.quadraticGaussJordan(4, Amatrix, b, x, NumericalMethods::PivotType::PARTIAL);

    float *L = nullptr;
    float *U = nullptr;
    test.LUFatoration(4, Amatrix, L, U);

    printf("Solution:\n");
    for(int i = 0; i < 4; i++) {
        printf("\t%f", x[i]);
    }

    return 0;
}
