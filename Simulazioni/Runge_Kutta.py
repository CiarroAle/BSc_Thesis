import numpy as np
import matplotlib.pyplot as plt

# Definizione dei parametri
beta1 = 1 / 15
beta2 = 1 / 2
gamma1 = 1 / 12
gamma2 = 1 / 10
N = 10000  # Popolazione totale

# Tempo di simulazione
T = 120  # giorni
dt = 1.0  # passo temporale
n_steps = int(T / dt)

# Inizializzazione delle variabili
s = np.zeros(n_steps)
i1 = np.zeros(n_steps)
i2 = np.zeros(n_steps)
r1 = np.zeros(n_steps)
r2 = np.zeros(n_steps)
t = np.linspace(0, T, n_steps)

# Condizioni iniziali
s[0] = N - 1010
i1[0] = 1000
i2[0] = 10  # Introdurre una piccola popolazione iniziale infetta da Omicron
r1[0] = 0
r2[0] = 0


# Funzione per calcolare i tassi di variazione
def derivatives(s, i1, i2, r1, r2):
    ds_dt = - (beta1 / N) * s * i1 - (beta2 / N) * s * i2
    di1_dt = (beta1 / N) * s * i1 - gamma1 * i1
    di2_dt = (beta2 / N) * s * i2 - gamma2 * i2
    dr1_dt = gamma1 * i1
    dr2_dt = gamma2 * i2
    return ds_dt, di1_dt, di2_dt, dr1_dt, dr2_dt


# Metodo di Runge-Kutta a 4 stadi (RK4)
for k in range(n_steps - 1):
    ds_dt1, di1_dt1, di2_dt1, dr1_dt1, dr2_dt1 = derivatives(s[k], i1[k], i2[k], r1[k], r2[k])

    s_temp = s[k] + 0.5 * ds_dt1 * dt
    i1_temp = i1[k] + 0.5 * di1_dt1 * dt
    i2_temp = i2[k] + 0.5 * di2_dt1 * dt
    r1_temp = r1[k] + 0.5 * dr1_dt1 * dt
    r2_temp = r2[k] + 0.5 * dr2_dt1 * dt

    ds_dt2, di1_dt2, di2_dt2, dr1_dt2, dr2_dt2 = derivatives(s_temp, i1_temp, i2_temp, r1_temp, r2_temp)

    s_temp = s[k] + 0.5 * ds_dt2 * dt
    i1_temp = i1[k] + 0.5 * di1_dt2 * dt
    i2_temp = i2[k] + 0.5 * di2_dt2 * dt
    r1_temp = r1[k] + 0.5 * dr1_dt2 * dt
    r2_temp = r2[k] + 0.5 * dr2_dt2 * dt

    ds_dt3, di1_dt3, di2_dt3, dr1_dt3, dr2_dt3 = derivatives(s_temp, i1_temp, i2_temp, r1_temp, r2_temp)

    s_temp = s[k] + ds_dt3 * dt
    i1_temp = i1[k] + di1_dt3 * dt
    i2_temp = i2[k] + di2_dt3 * dt
    r1_temp = r1[k] + dr1_dt3 * dt
    r2_temp = r2[k] + dr2_dt3 * dt

    ds_dt4, di1_dt4, di2_dt4, dr1_dt4, dr2_dt4 = derivatives(s_temp, i1_temp, i2_temp, r1_temp, r2_temp)

    s[k + 1] = s[k] + (dt / 6) * (ds_dt1 + 2 * ds_dt2 + 2 * ds_dt3 + ds_dt4)
    i1[k + 1] = i1[k] + (dt / 6) * (di1_dt1 + 2 * di1_dt2 + 2 * di1_dt3 + di1_dt4)
    i2[k + 1] = i2[k] + (dt / 6) * (di2_dt1 + 2 * di2_dt2 + 2 * di2_dt3 + di2_dt4)
    r1[k + 1] = r1[k] + (dt / 6) * (dr1_dt1 + 2 * dr1_dt2 + 2 * dr1_dt3 + dr1_dt4)
    r2[k + 1] = r2[k] + (dt / 6) * (dr2_dt1 + 2 * dr2_dt2 + 2 * dr2_dt3 + dr2_dt4)

# Plot dei risultati
plt.figure(figsize=(12, 8))
plt.plot(t, s, label='Suscettibili', color='blue')
plt.plot(t, i1, label='Infetti Delta', color='red')
plt.plot(t, i2, label='Infetti Omicron', color='green')
plt.plot(t, r1, label='Rimossi Delta', color='purple')
plt.plot(t, r2, label='Rimossi Omicron', color='orange')

plt.xlabel('Tempo (giorni)')
plt.ylabel('Numero di individui')
plt.title('Modello SIR con Varianti Delta e Omicron - Metodo RK4')
plt.legend()
plt.grid(True)
plt.show()
