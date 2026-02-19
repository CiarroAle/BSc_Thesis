import numpy as np
import matplotlib.pyplot as plt

# Definizione dei parametri
beta1 = 1/12
beta2 = 1/2
gamma1 = 1/12
gamma2 = 1/10
N = 10000 # Popolazione totale

# Tempo di simulazione
T = 120  # giorni
dt = 0.1  # passo temporale
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

# Metodo di Eulero esplicito per risolvere il sistema continuo
for k in range(n_steps - 1):
    ds_dt = - (beta1 / N) * s[k] * i1[k] - (beta2 / N) * s[k] * i2[k]
    di1_dt = (beta1 / N) * s[k] * i1[k] - gamma1 * i1[k]
    di2_dt = (beta2 / N) * s[k] * i2[k] - gamma2 * i2[k]
    dr1_dt = gamma1 * i1[k]
    dr2_dt = gamma2 * i2[k]

    s[k + 1] = s[k] + ds_dt * dt
    i1[k + 1] = i1[k] + di1_dt * dt
    i2[k + 1] = i2[k] + di2_dt * dt
    r1[k + 1] = r1[k] + dr1_dt * dt
    r2[k + 1] = r2[k] + dr2_dt * dt

# Plot dei risultati
plt.figure(figsize=(12, 8))
plt.plot(t, s, label='Suscettibili', color='blue')
plt.plot(t, i1, label='Infetti Delta', color='red')
plt.plot(t, i2, label='Infetti Omicron', color='green')
plt.plot(t, r1, label='Rimossi Delta', color='purple')
plt.plot(t, r2, label='Rimossi Omicron', color='orange')

plt.xlabel('Tempo (giorni)')
plt.ylabel('Numero di individui')
plt.title('Modello SIR con Varianti Delta e Omicron')
plt.legend()
plt.grid(True)
plt.show()
