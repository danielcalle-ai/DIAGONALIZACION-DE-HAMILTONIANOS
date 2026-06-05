#ifndef DIAGONALIZACION_H
#define DIAGONALIZACION_H

#include <vector>
#include <Eigen/Dense>

struct ResultadosDiagonalizacion {
    std::vector<double> autovalores;                 // E[n]
    std::vector<std::vector<double>> autovectores;   // C[i][n]
    std::vector<std::vector<double>> funcionesOnda;  // psi[n][i]
};

ResultadosDiagonalizacion diagonalizarHamiltoniano(
    const std::vector<std::vector<double>>& H,
    const std::vector<double>& pesos           
);

void normalizarFuncionesOnda(
    std::vector<std::vector<double>>& funcionesOnda,
    const std::vector<double>& pesos
);

void imprimirResultados(const ResultadosDiagonalizacion& res, int numEstados = 5);

#endif
