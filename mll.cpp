#include "nodos_pesos.h"
#include <gsl/gsl_integration.h>
#include <cmath>
#include <stdexcept>
// te llegan commits?
using namespace std;

void mallaLegendre::construir(int N) {
    nodos.resize(N); pesos.resize(N);
    gsl_integration_glfixed_table* t = gsl_integration_glfixed_table_alloc(N);
    for(int i = 0; i < N; i++) gsl_integration_glfixed_point(-1.0, 1.0, i, &nodos[i], &pesos[i], t);
    gsl_integration_glfixed_table_free(t);
}

void mallaLaguerre::construir(int N) {
    nodos.resize(N); pesos.resize(N);
    gsl_integration_fixed_workspace* w = gsl_integration_fixed_alloc(gsl_integration_fixed_laguerre, N, 0.0, 0.0, 0.0, 0.0);
    for(int i = 0; i < N; i++) { nodos[i] = gsl_integration_fixed_nodes(w)[i]; pesos[i] = gsl_integration_fixed_weights(w)[i]; }
    gsl_integration_fixed_free(w);
}

void mallaHermite::construir(int N) {
    nodos.resize(N); pesos.resize(N);
    gsl_integration_fixed_workspace* w = gsl_integration_fixed_alloc(gsl_integration_fixed_hermite, N, 0.0, 0.0, 0.0, 0.0);
    for(int i = 0; i < N; i++) { nodos[i] = gsl_integration_fixed_nodes(w)[i]; pesos[i] = gsl_integration_fixed_weights(w)[i]; }
    gsl_integration_fixed_free(w);
}

const vector<double>& mallaLegendre::getNodos() const { return nodos; }
const vector<double>& mallaLegendre::getPesos() const { return pesos; }

const vector<double>& mallaLaguerre::getNodos() const { return nodos; }
const vector<double>& mallaLaguerre::getPesos() const { return pesos; }

const vector<double>& mallaHermite::getNodos() const { return nodos; }
const vector<double>& mallaHermite::getPesos() const { return pesos; }

unique_ptr<malla> crearmalla(const string& tipo_malla) {
    if (tipo_malla == "Legendre") return make_unique<mallaLegendre>();
    if (tipo_malla == "Laguerre") return make_unique<mallaLaguerre>();
    if (tipo_malla == "Hermite") return make_unique<mallaHermite>();
    throw invalid_argument("Malla no soportada.");
}
