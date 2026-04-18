#include "nodos_pesos.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <cmath>
#include <stdexcept>

using namespace std;

// Implementación de mallaLegendre
void mallaLegendre::construir(int N) {
    nodos.resize(N);
    pesos.resize(N);
    
    gsl_integration_glfixed_table* t = gsl_integration_glfixed_table_alloc(N);
    
    for(int i = 0; i < N; i++) {
        double xi, wi;
        gsl_integration_glfixed_point(-1.0, 1.0, i, &xi, &wi, t);
        nodos[i] = xi;
        pesos[i] = wi;
    }
    
    gsl_integration_glfixed_table_free(t);
}

const vector<double>& mallaLegendre::getNodos() const { return nodos; }
const vector<double>& mallaLegendre::getPesos() const { return pesos; }

// Laguerre
void mallaLaguerre::construir(int N) {
    nodos.resize(N);
    pesos.resize(N);
    
    double alpha = 0.0;//laguerre estandar

    //construir matriz tridiagonal de jacobi
    
    gsl_matrix* J = gsl_matrix_alloc(N, N);
    gsl_matrix_set_zero(J);
    
    for (int i = 0; i < N; i++) {
      //diagonal
        gsl_matrix_set(J, i, i, 2.0 * i + 1.0 + alpha);
        
        if (i < N-1) {
            double beta = sqrt((i + 1.0) * (alpha + i + 1.0));
            gsl_matrix_set(J, i, i+1, beta);
            gsl_matrix_set(J, i+1, i, beta);
        }
    }
    //calcular autovalores y auto vectores
    gsl_vector* eval = gsl_vector_alloc(N);
    gsl_matrix* evec = gsl_matrix_alloc(N, N);
    
    gsl_eigen_symmv_workspace* work = gsl_eigen_symmv_alloc(N);
    gsl_eigen_symmv(J, eval, evec, work);
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
    //los autovalores son los nodos
    for (int i = 0; i < N; i++) {
      nodos[i] = gsl_vector_get(eval, i); //Pesos: (primer componente del autovector)^2 * factor de normalización
        double v0 = gsl_matrix_get(evec, i, 0);
        pesos[i] = v0 * v0 * tgamma(alpha + 1.0);
    }
    
    gsl_matrix_free(J);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_eigen_symmv_free(work);
}

const vector<double>& mallaLaguerre::getNodos() const { return nodos; }
const vector<double>& mallaLaguerre::getPesos() const { return pesos; }

// Implementación de mallaHermite
void mallaHermite::construir(int N) {
    nodos.resize(N);
    pesos.resize(N);
    
    //crear matris tridiagonal
    gsl_matrix* J = gsl_matrix_alloc(N, N);
    gsl_matrix_set_zero(J);
    
    for (int i = 0; i < N; i++) {
        gsl_matrix_set(J, i, i, 0.0);
        
        if (i < N-1) {
            double beta = sqrt((i+1) / 2.0);
            gsl_matrix_set(J, i, i+1, beta);
            gsl_matrix_set(J, i+1, i, beta);
        }
    }
    
    gsl_vector* eval = gsl_vector_alloc(N);
    gsl_matrix* evec = gsl_matrix_alloc(N, N);
    
    gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(N);
    gsl_eigen_symmv(J, eval, evec, w);
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
    //los autovalores son nodos
    for (int i = 0; i < N; i++) {
        nodos[i] = gsl_vector_get(eval, i);
        double v0 = gsl_matrix_get(evec, i, 0);
        pesos[i] = v0 * v0 * sqrt(M_PI);
    }
    
    gsl_matrix_free(J);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_eigen_symmv_free(w);
}

const vector<double>& mallaHermite::getNodos() const { return nodos; }
const vector<double>& mallaHermite::getPesos() const { return pesos; }

// Función factory
unique_ptr<malla> crearmalla(const string& tipo_malla) {
    if (tipo_malla == "Legendre")
        return make_unique<mallaLegendre>();
    else if (tipo_malla == "Laguerre")
        return make_unique<mallaLaguerre>();
    else if (tipo_malla == "Hermite")
        return make_unique<mallaHermite>();
    else
        throw invalid_argument("No esta esa malla");
}
