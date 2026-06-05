#include "matrices.h"
#include <cmath>
#include <stdexcept>
#include <vector>

using namespace std;

vector<double> escalamientoYmapeo_r(const vector<double>& nodos_originales, double r_a, double r_b) {
    int N = nodos_originales.size();
    vector<double> nodos_r(N);
    for (int i = 0; i < N; i++) {
        nodos_r[i] = ((r_b - r_a) / 2.0) * nodos_originales[i] + (r_a + r_b) / 2.0;
    }
    return nodos_r;
}

vector<double> escalamientoYmapeo_u(const vector<double>& nodos_originales) {
    return nodos_originales;
}

vector<vector<double>> construirHamiltoniano(
    int N_r, int N_u,
    const vector<double>& nodos_r_originales,
    const vector<double>& nodos_u_originales,
    const ParametrosFisicos& params,
    FuncionPotencial2D V,
    double r_a, double r_b) {
    
    int N_total = N_r * N_u;
    vector<vector<double>> H(N_total, vector<double>(N_total, 0.0));
    
    vector<double> nodos_r = escalamientoYmapeo_r(nodos_r_originales, r_a, r_b);
    vector<double> nodos_u = escalamientoYmapeo_u(nodos_u_originales);
    
    double hbar2_2m = (params.hbar * params.hbar) / (2.0 * params.m);
    double factor_radial = 2.0 / (r_b - r_a);
    double factor_radial_cuadrado = factor_radial * factor_radial;
    
    // 1. MATRIZ DE SEGUNDA DERIVADA RADIAL (T_r)
    vector<vector<double>> T_r(N_r, vector<double>(N_r, 0.0));
    for (int i = 0; i < N_r; i++) {
        double xi = nodos_r_originales[i];
        for (int k = 0; k < N_r; k++) {
            double xk = nodos_r_originales[k];
            if (i == k) {
                T_r[i][i] = (N_r * (N_r + 1.0)) / (3.0 * (1.0 - xi * xi)) + (xi * xi) / ((1.0 - xi * xi) * (1.0 - xi * xi));
            } else {
                double signo = ((i - k) % 2 == 0) ? 1.0 : -1.0;
                T_r[i][k] = (signo / ((xi - xk) * (xi - xk))) * sqrt((1.0 - xi * xi) / (1.0 - xk * xk));
            }
        }
    }
    
    // 2. MATRIZ ANGULAR COMPLETA (incluye m²/(1-u²) para cada m)
    // Nota: Esta función se llama por cada m, así que calculamos para el m actual
    auto construirMatrizAngular = [&](int m_quantum) {
        vector<vector<double>> S(N_u, vector<double>(N_u, 0.0));
        for (int j = 0; j < N_u; j++) {
            double uj = nodos_u_originales[j];
            double den_j = 1.0 - uj * uj;
            if (den_j < 1e-12) den_j = 1e-12;
            
            for (int l = 0; l < N_u; l++) {
                if (j == l) {
                    // Término diagonal: incluye el término centrífugo m²/(1-u²)
                    S[j][j] = (N_u * (N_u + 1.0) * den_j + 2.0) / (3.0 * den_j * den_j) + 
                              double(m_quantum * m_quantum) / den_j;
                } else {
                    double ul = nodos_u_originales[l];
                    double den_l = 1.0 - ul * ul;
                    if (den_l < 1e-12) den_l = 1e-12;
                    double signo = ((j - l) % 2 == 0) ? 1.0 : -1.0;
                    S[j][l] = signo * 2.0 * sqrt(den_j * den_l) / ((uj - ul) * (uj - ul));
                }
            }
        }
        return S;
    };
    
    // Por simplicidad, asumimos m=0 (estado fundamental). Para otros m, habría que pasar m como parámetro
    vector<vector<double>> S = construirMatrizAngular(0);
    
    // 3. ENSAMBLAJE DEL HAMILTONIANO
    for (int i = 0; i < N_r; i++) {
        for (int j = 0; j < N_u; j++) {
            int idx = i * N_u + j;
            double r_i = nodos_r[i];
            double u_j = nodos_u[j];
            
            for (int k = 0; k < N_r; k++) {
                for (int l = 0; l < N_u; l++) {
                    int idx2 = k * N_u + l;
                    double cinetica = 0.0;
                    
                    if (j == l) {
                        cinetica += hbar2_2m * T_r[i][k] * factor_radial_cuadrado;
                    }
                    
                    if (i == k) {
                        cinetica += hbar2_2m * S[j][l] / (r_i * r_i);
                    }
                    
                    H[idx][idx2] = cinetica;
                }
            }
            H[idx][idx] += V(r_i, u_j);
        }
    }
    
    return H;
}

// Stub requerido por compatibilidad
vector<vector<double>> construirMatrizCinetica(int N_r, int N_u, const vector<double>& nodos_r, const vector<double>& nodos_u, const ParametrosFisicos& params) {
    int N_total = N_r * N_u;
    return vector<vector<double>>(N_total, vector<double>(N_total, 0.0));
}

vector<vector<double>> construirMatrizPotencial(int N_r, int N_u, const vector<double>& nodos_r, const vector<double>& nodos_u, FuncionPotencial2D V) {
    int N_total = N_r * N_u;
    return vector<vector<double>>(N_total, vector<double>(N_total, 0.0));
}
