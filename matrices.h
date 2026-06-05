#ifndef MATRICES_H
#define MATRICES_H

#include "nodos_pesos.h"
#include <vector>
#include <string>
#include <functional>

struct ParametrosFisicos {
    double hbar;   // ħ reducido
    double m;      // masa efectiva
    double h;      // parámetro de escala (no usado, por compatibilidad)
};

struct ParametrosTI {
    double r_a;           // radio del TI [m]
    double r_b;           // radio externo del QD [m]
    double B;             // campo magnético externo [T]
    double theta_bar;     // parámetro magnetoeléctrico θ̄
    double epsilon1;      // permitividad del TI
    double epsilon2;      // permitividad de GaAs
    double epsilon0;      // permitividad del vacío
    double mu0;           // permeabilidad del vacío
    double V0;            // profundidad del potencial de confinamiento [J]
    double q;             // carga del electrón [C]
    double c;             // velocidad de la luz [m/s]
};

using FuncionPotencial2D = std::function<double(double r, double u)>;

std::vector<double> escalamientoYmapeo_r(const std::vector<double>& nodos_originales, double r_a, double r_b);
std::vector<double> escalamientoYmapeo_u(const std::vector<double>& nodos_originales);

std::vector<std::vector<double>> construirMatrizCinetica(
    int N_r, int N_u,
    const std::vector<double>& nodos_r,
    const std::vector<double>& nodos_u,
    const ParametrosFisicos& params,
    double r_a, double r_b); // Añadidos radios físicos para corregir la derivada radial

std::vector<std::vector<double>> construirMatrizPotencial(
    int N_r, int N_u,
    const std::vector<double>& nodos_r,
    const std::vector<double>& nodos_u,
    FuncionPotencial2D V);

std::vector<std::vector<double>> construirHamiltoniano(
    int N_r, int N_u,
    const std::vector<double>& nodos_r_originales,
    const std::vector<double>& nodos_u_originales,
    const ParametrosFisicos& params,
    FuncionPotencial2D V,
    double r_a, double r_b); // Añadidos para el pasaje de la escala real del core-shell

#endif
