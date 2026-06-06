#include "matrices.h"
#include <cmath>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <iomanip>

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
    
    cout << "\n[DEBUG] Parámetros de construcción:" << endl;
    cout << "  ℏ²/2m = " << hbar2_2m << endl;
    cout << "  Factor radial = " << factor_radial << endl;
    cout << "  Factor radial² = " << factor_radial_cuadrado << endl;
    cout << "  Dominio: [" << r_a << ", " << r_b << "]" << endl;
    
    // 1. MATRIZ DE SEGUNDA DERIVADA RADIAL (T_r)
    vector<vector<double>> T_r(N_r, vector<double>(N_r, 0.0));
    for (int i = 0; i < N_r; i++) {
        double xi = nodos_r_originales[i];
        double denom_i = (1.0 - xi * xi);
        if (fabs(denom_i) < 1e-14) denom_i = 1e-14;
        
        for (int k = 0; k < N_r; k++) {
            double xk = nodos_r_originales[k];
            double denom_k = (1.0 - xk * xk);
            if (fabs(denom_k) < 1e-14) denom_k = 1e-14;
            
            if (i == k) {
                T_r[i][i] = (N_r * (N_r + 1.0)) / (3.0 * denom_i) + (xi * xi) / (denom_i * denom_i);
            } else {
                double signo = ((i - k) % 2 == 0) ? 1.0 : -1.0;
                T_r[i][k] = (signo / ((xi - xk) * (xi - xk))) * sqrt(denom_i / denom_k);
            }
        }
    }
    
    // 2. MATRIZ ANGULAR COMPLETA
    auto construirMatrizAngular = [&](int m_quantum) {
        vector<vector<double>> S(N_u, vector<double>(N_u, 0.0));
        for (int j = 0; j < N_u; j++) {
            double uj = nodos_u_originales[j];
            double den_j = 1.0 - uj * uj;
            if (fabs(den_j) < 1e-14) den_j = 1e-14;
            
            for (int l = 0; l < N_u; l++) {
                if (j == l) {
                    S[j][j] = (N_u * (N_u + 1.0) * den_j + 2.0) / (3.0 * den_j * den_j) + 
                              double(m_quantum * m_quantum) / den_j;
                } else {
                    double ul = nodos_u_originales[l];
                    double den_l = 1.0 - ul * ul;
                    if (fabs(den_l) < 1e-14) den_l = 1e-14;
                    double signo = ((j - l) % 2 == 0) ? 1.0 : -1.0;
                    S[j][l] = signo * 2.0 * sqrt(den_j * den_l) / ((uj - ul) * (uj - ul));
                }
            }
        }
        return S;
    };
    
    vector<vector<double>> S = construirMatrizAngular(0);
    
    cout << "[DEBUG] Matrices de derivadas calculadas" << endl;
    
    // 3. ENSAMBLAJE CORRECTO DEL HAMILTONIANO
    cout << "[DEBUG] Ensamblando Hamiltoniano..." << endl;
    
    for (int i = 0; i < N_r; i++) {
        for (int j = 0; j < N_u; j++) {
            int idx = i * N_u + j;
            double r_i = nodos_r[i];
            double u_j = nodos_u[j];
            
            // PASO 1: Inicializar con potencial en diagonal
            double V_ij = V(r_i, u_j);
            H[idx][idx] = V_ij;
            
            // PASO 2: Sumar términos cinéticos
            for (int k = 0; k < N_r; k++) {
                for (int l = 0; l < N_u; l++) {
                    int idx2 = k * N_u + l;
                    double cinetica = 0.0;
                    
                    // Término radial (solo si j == l)
                    if (j == l) {
                        cinetica += hbar2_2m * T_r[i][k] * factor_radial_cuadrado;
                    }
                    
                    // Término angular (solo si i == k, y dividido por r²)
                    if (i == k) {
                        if (r_i > 1e-10) {
                            cinetica += hbar2_2m * S[j][l] / (r_i * r_i);
                        } else {
                            cinetica += 0.0;  // En r=0, el término angular se anula
                        }
                    }
                    
                    // Asignar correctamente
                    if (idx == idx2) {
                        // Diagonal: SUMA cinética al potencial
                        H[idx][idx] += cinetica;
                    } else {
                        // Fuera de diagonal: solo cinética
                        H[idx][idx2] = cinetica;
                    }
                }
            }
        }
    }
    
    cout << "[DEBUG] Hamiltoniano ensamblado" << endl;
    
    // ========== VERIFICACIÓN DE HERMITICIDAD ==========
    cout << "\n[DEBUG] VERIFICACIÓN DE HERMITICIDAD" << endl;
    cout << "========================================" << endl;
    
    double max_asimetria = 0.0;
    double max_elemento = 0.0;
    int pos_max_i = 0, pos_max_j = 0;
    
    for (int i = 0; i < N_total; i++) {
        for (int j = i + 1; j < N_total; j++) {
            double h_ij = H[i][j];
            double h_ji = H[j][i];
            double diff = fabs(h_ij - h_ji);
            double max_elem_pair = max(fabs(h_ij), fabs(h_ji));
            
            if (diff > max_asimetria) {
                max_asimetria = diff;
                pos_max_i = i;
                pos_max_j = j;
            }
            
            if (max_elem_pair > max_elemento) {
                max_elemento = max_elem_pair;
            }
        }
    }
    
    cout << "Máximo elemento en H: " << scientific << setprecision(6) << max_elemento << endl;
    cout << "Máxima asimetría: " << scientific << setprecision(6) << max_asimetria << endl;
    cout << "Posición: H[" << pos_max_i << "][" << pos_max_j << "]" << endl;
    
    if (max_asimetria > 0) {
        double rel_error = max_asimetria / (max_elemento + 1e-15);
        cout << "Error relativo: " << scientific << setprecision(3) << rel_error << endl;
        
        if (rel_error < 1e-10) {
            cout << "✓ Matriz es hermitiana (error < 1e-10)" << endl;
        } else if (rel_error < 1e-8) {
            cout << "✓ Matriz es aproximadamente hermitiana (error < 1e-8)" << endl;
        } else if (rel_error < 1e-6) {
            cout << "⚠ Matriz tiene pequeños errores de hermiticidad (error < 1e-6)" << endl;
        } else {
            cout << "✗ ADVERTENCIA: Matriz NO es hermitiana (error > 1e-6)" << endl;
            cout << "  H[" << pos_max_i << "][" << pos_max_j << "] = " << H[pos_max_i][pos_max_j] << endl;
            cout << "  H[" << pos_max_j << "][" << pos_max_i << "] = " << H[pos_max_j][pos_max_i] << endl;
        }
    }
    
    cout << "========================================\n" << endl;
    
    // ========== ESTADÍSTICAS DE LA MATRIZ ==========
    cout << "[DEBUG] ESTADÍSTICAS DE LA MATRIZ" << endl;
    cout << "================================" << endl;
    
    double sum_diag = 0.0;
    double min_diag = 1e10;
    double max_diag = -1e10;
    int nnz = 0;  // Número de elementos no cero
    
    for (int i = 0; i < N_total; i++) {
        sum_diag += H[i][i];
        min_diag = min(min_diag, H[i][i]);
        max_diag = max(max_diag, H[i][i]);
        
        for (int j = 0; j < N_total; j++) {
            if (fabs(H[i][j]) > 1e-15) nnz++;
        }
    }
    
    cout << "Traza (sum diagonal): " << scientific << setprecision(6) << sum_diag << endl;
    cout << "Min diagonal: " << min_diag << endl;
    cout << "Max diagonal: " << max_diag << endl;
    cout << "Elementos no-cero: " << nnz << " / " << N_total*N_total 
         << " (" << fixed << setprecision(2) << (100.0*nnz)/(N_total*N_total) << "%)" << endl;
    cout << "================================\n" << endl;
    
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
