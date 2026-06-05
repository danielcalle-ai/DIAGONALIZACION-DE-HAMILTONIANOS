#include "diagonalizacion.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdexcept>

using namespace Eigen;
using namespace std;

ResultadosDiagonalizacion diagonalizarHamiltoniano(const vector<vector<double>>& H, const vector<double>& pesos) {
    int N = H.size();
    ResultadosDiagonalizacion res;
    
    if (N == 0 || H[0].size() != N) {
        throw invalid_argument("La matriz del Hamiltoniano debe ser cuadrada.");
    }

    MatrixXd H_eigen(N, N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            H_eigen(i, j) = H[i][j];
        }
    }
    
    SelfAdjointEigenSolver<MatrixXd> solver(H_eigen);
    if (solver.info() != Success) {
        throw runtime_error("Error critico: La diagonalizacion no convergio.");
    }
    
    VectorXd En = solver.eigenvalues();
    MatrixXd Cn = solver.eigenvectors();
    
    res.autovalores.resize(N);
    res.autovectores.assign(N, vector<double>(N));
    res.funcionesOnda.assign(N, vector<double>(N));
    
    for (int i = 0; i < N; i++) {
        res.autovalores[i] = En(i);
    }
    
    for (int n = 0; n < N; n++) {
        for (int i = 0; i < N; i++) {
            res.autovectores[i][n] = Cn(i, n);
            if (pesos[i] > 1e-15) {
                res.funcionesOnda[n][i] = Cn(i, n) / sqrt(pesos[i]);
            } else {
                res.funcionesOnda[n][i] = 0.0;
            }
        }
    }
    
    normalizarFuncionesOnda(res.funcionesOnda, pesos);
    
    return res;
}

void normalizarFuncionesOnda(vector<vector<double>>& funcionesOnda, const vector<double>& pesos) {
    int N_estados = funcionesOnda.size();
    int N_puntos = pesos.size();
    
    for (int n = 0; n < N_estados; n++) {
        double norma = 0.0;
        for (int i = 0; i < N_puntos; i++) {
            norma += funcionesOnda[n][i] * funcionesOnda[n][i] * pesos[i];
        }
        norma = sqrt(norma);
        
        if (norma > 1e-12) {
            for (int i = 0; i < N_puntos; i++) {
                funcionesOnda[n][i] /= norma;
            }
        }
    }
}

void imprimirResultados(const ResultadosDiagonalizacion& res, int numEstados) {
    int N = res.autovalores.size();
    int mostrar = min(numEstados, N);
    
    cout << "\n=== AUTOVALORES DE ENERGIA CALCULADOS ===" << endl;
    for (int n = 0; n < mostrar; n++) {
        cout << "Estado E[" << n << "] = " << res.autovalores[n] << endl;
    }
}
