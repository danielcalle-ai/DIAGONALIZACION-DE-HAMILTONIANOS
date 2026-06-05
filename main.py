import numpy as np
import matplotlib.pyplot as plt
import pipeline_cuantico_cpp

alpha = 1/137.0

# Parametros del paper
r_a_nm = 5.0
r_b_nm = 10.0
V0_meV = 0.0
epsilon1 = 4.0
epsilon2 = 12.53
epsilon0 = 8.854e-12
mu0 = 1.2566e-6
c = 3.0e8
m_efectiva = 0.067
carga_q = 1.602e-19
B_Tesla = 3.0

# Conversion a unidades atomicas
nm_a_bohr = 0.0529177
meV_a_ua = 1.0 / 27211.4

r_a = r_a_nm / nm_a_bohr
r_b = r_b_nm / nm_a_bohr
V0 = V0_meV * meV_a_ua

# ========== CONTROLES ==========
TI_encendido = False  # True = con aislante, False = sin aislante
B_encendido = False    # True = con campo, False = sin campo
# ===============================

theta_param = (11.0 * alpha) if TI_encendido else 0.0
B_valor = B_Tesla if B_encendido else 0.0

print(f"\n--- CONFIGURACION ---")
print(f"Aislante TI: {'ON' if TI_encendido else 'OFF'}")
print(f"Campo B: {B_valor} T")
print(f"V0: {V0_meV} meV = {V0:.6f} u.a.")
print(f"r_a = {r_a:.2f} u.a., r_b = {r_b:.2f} u.a.")

# Ejecutar
resultados = pipeline_cuantico_cpp.ejecutar_pipeline_cuantico(
    50, 50, r_a, r_b, V0, B_valor, theta_param,
    epsilon1, epsilon2, epsilon0, mu0, c, m_efectiva, carga_q, 0.0
)

energias = np.array(resultados.energias)
funciones_onda = resultados.funciones_onda
nodos_r = np.array(resultados.nodos_r)
nodos_u = np.array(resultados.nodos_u)

print(f"\n--- ENERGIAS [u.a.] ---")
for i in range(min(5, len(energias))):
    print(f"E{i} = {energias[i]:.8f}")

# Reconstruccion de la funcion de onda
Nrad = len(nodos_r)
Nang = len(nodos_u)
psi_vector = np.array(funciones_onda[0])
psi_2d = np.zeros((Nang, Nrad))

for i in range(Nrad):
    for j in range(Nang):
        idx = i * Nang + j
        if idx < len(psi_vector):
            psi_2d[j, i] = psi_vector[idx]

nodos_r_nm = nodos_r * nm_a_bohr
thetas = np.arccos(nodos_u)

idx_r = np.argsort(nodos_r_nm)
idx_u = np.argsort(thetas)

R_mesh, Theta_mesh = np.meshgrid(nodos_r_nm[idx_r], thetas[idx_u])
psi_ordenada = psi_2d[idx_u, :][:, idx_r]
psi_vis = np.abs(psi_ordenada)
if np.max(psi_vis) > 0:
    psi_vis /= np.max(psi_vis)

# Grafico polar
fig, ax = plt.subplots(figsize=(9, 8), subplot_kw=dict(projection='polar'))
ax.contourf(Theta_mesh, R_mesh, psi_vis, levels=100, cmap='viridis')
ax.set_ylim(0, r_b_nm)
ax.set_title(f"Estado fundamental | E0 = {energias[0]:.6f} u.a. | TI={'ON' if TI_encendido else 'OFF'}, B={B_valor}T")
plt.show()
