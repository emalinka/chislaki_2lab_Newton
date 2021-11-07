//
//  main.cpp
//  lab2newton
//
//  Created by Elizaveta on 7.11.21.
//

#include <iostream>
#include <math.h>
#include <iomanip>


void PrintMatrix(double** matrix_A, int size){ //вывод матрицы
    for (int i = 0; i < size; i++){
        for (int j = 0; j <= size; j++)
            std::cout << std::setw(15) << matrix_A[i][j];
        std::cout << '\n';
    }
}

void FillVectorF(double* vector_F, double* vector_X){
    vector_F[0] = 1.5 * vector_X[0] * vector_X[0] * vector_X[0] - vector_X[1] * vector_X[1] - 1;
    vector_F[1] = vector_X[0]*vector_X[1] * vector_X[1] * vector_X[1] - vector_X[1] - 4;
}

void FillMatrixJacobi_I(double** matrix_Jacobi,  double* vector_X){
    matrix_Jacobi[0][0] = 4.5 * vector_X[0] * vector_X[0];
    matrix_Jacobi[0][1] = -2 * vector_X[1];
    matrix_Jacobi[1][0] = vector_X[1] * vector_X[1] * vector_X[1];
    matrix_Jacobi[1][1] = 3 * vector_X[0] * vector_X[1] * vector_X[1] - 1;
}

void FillMatrixJacobi_II(double** matrix_Jacobi, double* vector_F, double* vector_X, int size, double increment_M){
    double F1, F2;
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++){
            vector_X[j] += increment_M;
            FillVectorF(vector_F, vector_X);
            F1 = vector_F[i];
            vector_X[j] -= increment_M;
            FillVectorF(vector_F, vector_X);
            F2 = vector_F[i];
            matrix_Jacobi[i][j] = (F1 - F2) / increment_M;
            matrix_Jacobi[i][size] = -vector_F[i];
        }
}


void GaussMethod(double** matrix_Jacobi, double* vector_X, double* vector_X2, int size){
    for (int i = 0; i < size; i++){
        double max = fabs(matrix_Jacobi[i][i]);
        int index = i;
        for (int j = i; j < size; j++)
            if (fabs(matrix_Jacobi[j][i]) > max){
                max = fabs(matrix_Jacobi[j][i]);
                index = j;
            }

        if (index != i){
            double* temp = matrix_Jacobi[i];
            matrix_Jacobi[i] = matrix_Jacobi[index];
            matrix_Jacobi[index] = temp;
        }

        double mainelement = matrix_Jacobi[i][i];
        for (int j = 0; j < size + 1; j++)
            matrix_Jacobi[i][j] /= mainelement;
        

        for (int j = i + 1; j < size; j++){
            double temp = matrix_Jacobi[j][i];
            for (int k = i; k < size + 1; k++)
                matrix_Jacobi[j][k] -= matrix_Jacobi[i][k] * temp;
        }
        //cout << "\nTransformed matrix:\n";
        //PrintMatrix(matrix_Jacobi, size);
    }

    for (int i = size - 1; i > 0; i--){
        double temp = matrix_Jacobi[i][size];
        for (int j = i - 1; j >= 0; j--)
            matrix_Jacobi[j][size] -= matrix_Jacobi[j][i] * temp;
        //cout << "\n___Reversed Gauss method___\n";
        //PrintMatrix(matrix_Jacobi, size);
    }

    for (int i = 0; i < size; ++i)
        vector_X2[i] = matrix_Jacobi[i][size];

    for (int i = 0; i < size; i++)
        vector_X[i] += matrix_Jacobi[i][size];
}

void NewtonMethod(double** matrix_Jacobi, double* vector_F, double* vector_X, double* vector_X2, int size, int iter, double increment_M, double eps){
    
    std::cout << std::setw(10) << "X1" << std::setw(10) << "X2" << std::setw(17) << "Delta_1" << std::setw(17)
    << "Delta_2" << std::setw(6) << "Iter" << '\n';
       
    double delta1, delta2,max;
    int c = 0;
    FillMatrixJacobi_II(matrix_Jacobi, vector_F, vector_X, size, increment_M);

//    PrintMatrix(matrix_Jacobi, size);
    std::cout << '\n';
    while (true){
        FillVectorF(vector_F, vector_X);
        
        FillMatrixJacobi_I(matrix_Jacobi, vector_X);

        
        
        GaussMethod(matrix_Jacobi, vector_X, vector_X2, size);
        max = 0;

        FillVectorF(vector_F, vector_X);
        for (int i = 0; i < size; i++)
            if (fabs(vector_F[i]) > max)
                max = fabs(vector_F[i]);

        delta1 = max;
        max = 0;

        for (int i = 0; i < size; i++){
            if (fabs(vector_X[i]) < 1 && fabs(vector_X2[i]) > max)
                max = fabs(vector_X2[i]);
            if (fabs(vector_X[i]) >= 1 && fabs(vector_X2[i] / vector_X[i]) > max)
                max = fabs(vector_X2[i] / vector_X[i]);
        }

        delta2 = max;
        std::cout << std::setw(10) << vector_X[0] << std::setw(10) << vector_X[1] << std::setw(17) << delta1
            << std::setw(17) << delta2 << std::setw(6) << ++c << '\n';
        if (delta1 <= eps && delta2 <= eps || c >= iter)
            break;
    }
}


int main() {
    std::cout << "Newton's method \n";
    double eps = 1e-9;
    double increment_M = 0.01;
    
    std::cout << "Введите количество итераций: ";
    int iter;
    std::cin >> iter;
    int size = 2;
    double* vector_F = new double[size];
    
    double* vector_X = new double[size];
    double* vector_X2 = new double[size];
    std::cout << "\nВведите 2 начальных приближения: ";
    for (int i = 0; i < size; i++)
        std::cin >> vector_X[i];
    
    double** matrix_Jacobi = new double* [size];
    for (int i = 0; i < size; i++)
        matrix_Jacobi[i] = new double[size + 1];

    NewtonMethod(matrix_Jacobi, vector_F, vector_X, vector_X2, size, iter, increment_M, eps);

    delete[] vector_F;
    delete[] vector_X;
    delete[] vector_X2;
    for (int i = 0; i < size; i++)
        delete[] matrix_Jacobi[i];
    delete[] matrix_Jacobi;

    return 0;
}
