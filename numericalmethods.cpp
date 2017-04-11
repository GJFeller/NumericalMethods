#include "numericalmethods.h"

NumericalMethods::NumericalMethods()
{
    //ctor
}

void NumericalMethods::quadraticGaussJordan(int n, float* A, float* b, float* x, PivotType pivotType)
{
    float *matrix =  (float *) malloc(sizeof (float) * n * (n+1));

    //Ab copy to matrix
    for(int i = 0; i < n; i++){
        for(int j = 0; j <= n; j++) {
            if(j != n) {
                matrix[i*(n+1)+j] = A[i*n+j];
            } else {
                matrix[i*(n+1)+j] = b[i];
            }
        }
    }
    printf("Ab matrix:\n");
    print_mat(n, n+1, matrix);

    switch(pivotType) {
        case PARTIAL:
            partialPivot(n, matrix);
            break;
        default:
            normalPivot(n, matrix);
    }

    printf("Triangular matrix:\n");
    print_mat(n, n+1, matrix);
    // Substitution
    for(int i = n-1; i >= 0; i--){
        for(int j = i+1; j < n; j++){
            matrix[i*(n+1)+n] = matrix[i*(n+1)+n] - x[j]*matrix[i*(n+1)+j];
        }
        x[i] = matrix[i*(n+1)+n]/matrix[i*(n+1)+i];
    }
    /*printf("Solution:\n");
    for(int i = 0; i < 4; i++) {
        printf("\t%f", x[i]);
    }*/
    return;
}

void NumericalMethods::print_mat(int r, int c, float* m)
{
        for (int i = 0; i < r; i++) {
                for (int j = 0; j < c; j++) {
                        printf("\t%f", m[i*c+j]);
                }
                printf("\n");
        }
        return;
}

void NumericalMethods::normalPivot(int n, float *Ab){
    // Triangulation
    for(int i = 0; i < n-1; i++) {
        for(int j = i+1; j < n; j++) {
            float mult = -Ab[j*(n+1)+i]/Ab[i*(n+1)+i];
            for(int k = i; k <= n; k++) {
                Ab[j*(n+1)+k] = mult * Ab[i*(n+1)+k] + Ab[j*(n+1)+k];
            }
        }
    }
}

void NumericalMethods::partialPivot(int n, float *Ab){
    // Triangulation
    for(int i = 0; i < n-1; i++) {
        //Find the maximum value for pivot
        float max = Ab[i*(n+1)+i];
        int maxLine = i;
        for(int lin = i+1; lin < n; lin++) {
            if(Ab[lin*(n+1)+i] > max) {
                max = Ab[lin*(n+1)+i];
                maxLine = lin;
            }
        }
        if(maxLine != i) {
            swapLines(n+1, Ab, i, maxLine);
        }
        for(int j = i+1; j < n; j++) {
            float mult = -Ab[j*(n+1)+i]/Ab[i*(n+1)+i];
            for(int k = i; k <= n; k++) {
                Ab[j*(n+1)+k] = mult * Ab[i*(n+1)+k] + Ab[j*(n+1)+k];
            }
        }
    }
}

void NumericalMethods::swapLines(int n, float *Ab, int lin1, int lin2) {
    for(int j = 0; j < n; j++) {
        float copy = Ab[lin1*(n)+j];
        Ab[lin1*(n)+j] = Ab[lin2*(n)+j];
        Ab[lin2*(n)+j] = copy;
    }
}

void NumericalMethods::LUFatoration(int n, float *A, float *L, float *U){
    if(L == nullptr || U == nullptr) {
        L = (float *) malloc(sizeof (float) * n * n);
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                if(i == j) {
                    L[i*n+j] = 1;
                } else {
                    L[i*n+j] = 0;
                }
            }
        }
        U = (float *) malloc(sizeof (float) * n * n);
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                U[i*n+j] = A[i*n+j];
            }
        }
    }
    /*
    // Triangulation
    for(int i = 0; i < n-1; i++) {
        //Find the maximum value for pivot
        float max = U[i*(n)+i];
        int maxLine = i;
        for(int lin = i+1; lin < n; lin++) {
            if(U[lin*(n)+i] > max) {
                max = U[lin*(n)+i];
                maxLine = lin;
            }
        }
        if(maxLine != i) {
            swapLines(n, U, i, maxLine);
        }
        for(int j = i+1; j < n; j++) {
            float mult = -U[j*(n)+i]/U[i*(n)+i];
            L[j*n+i] = -mult;
            for(int k = i; k <= n; k++) {
                U[j*(n)+k] = mult * U[i*(n)+k] + U[j*(n)+k];
            }
        }
    }*/
    // Triangulation
    for(int i = 0; i < n-1; i++) {
        for(int j = i+1; j < n; j++) {
            float mult = -U[j*(n)+i]/U[i*(n)+i];
            L[j*n+i] = -mult;
            for(int k = i; k <= n; k++) {
                U[j*(n)+k] = mult * U[i*(n)+k] + U[j*(n)+k];
            }
        }
    }

    printf("L:\n");
    print_mat(n, n, L);
    printf("U:\n");
    print_mat(n, n, U);

}
