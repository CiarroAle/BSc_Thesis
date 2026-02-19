import numpy as np
import matplotlib.pyplot as plt

# Definizione dei parametri
beta1 = 1/15
beta2 = 1/2
gamma1 = 1/12
gamma2 = 1/10
N = 10000  # Popolazione totale

# Condizioni iniziali
S = N - 1010
I1 = 1000
I2 = 10
R1 = 0
R2 = 0

# Tempo di simulazione
T_max = 40
t = 0

# Inizializzazione delle liste per tracciare i risultati
time_points = [t]
S_values = [S]
I1_values = [I1]
I2_values = [I2]
R1_values = [R1]
R2_values = [R2]

# Inizializzazione delle liste per tracciare le nuove infezioni giornaliere
time_points_daily = np.arange(0, T_max+1)
h = len(time_points_daily)
new_delta_infections = np.zeros(h)
new_omicron_infections = np.zeros(h)
new_delta_removals = np.zeros(h)
new_omicron_removals = np.zeros(h)

while t < T_max-1 and (I1 > 0 or I2 > 0):
    # Calcolo dei tassi di evento
    rate_infection1 = beta1 * S * I1 / N
    rate_infection2 = beta2 * S * I2 / N
    rate_recovery1 = gamma1 * I1
    rate_recovery2 = gamma2 * I2
    total_rate = rate_infection1 + rate_infection2 + rate_recovery1 + rate_recovery2

    if total_rate == 0:
        break

    # Determinazione del tempo fino al prossimo evento
    delta_t = np.random.exponential(1 / total_rate)
    t += delta_t
    day = int(np.floor(t))

    # Determinazione del tipo di evento
    r = np.random.uniform(0, total_rate)
    if r < rate_infection1: # infezione con delta
        S -= 1
        I1 += 1
        new_delta_infections[day] += 1

    elif r < rate_infection1 + rate_infection2: # infezione con omicron
        S -= 1
        I2 += 1
        new_omicron_infections[day] += 1

    elif r < rate_infection1 + rate_infection2 + rate_recovery1: # rimozione da delta
        I1 -= 1
        R1 += 1
        new_delta_removals[day] += 1

    else: # rimozione da omicron
        I2 -= 1
        R2 += 1
        new_omicron_removals[day] += 1


    # Registrazione dei risultati
    time_points.append(t)
    S_values.append(S)
    I1_values.append(I1)
    I2_values.append(I2)
    R1_values.append(R1)
    R2_values.append(R2)


# operazioni componente per componente
new_delta_infections_effective = new_delta_infections - new_delta_removals
new_omicron_infections_effective = new_omicron_infections - new_omicron_removals




# Plot dei risultati
plt.figure(figsize=(12, 8))
plt.plot(time_points, S_values, label='Suscettibili', color='blue')
plt.plot(time_points, I1_values, label='Infetti Delta', color='red')
plt.plot(time_points, I2_values, label='Infetti Omicron', color='green')
plt.plot(time_points, R1_values, label='Rimossi Delta', color='purple')
plt.plot(time_points, R2_values, label='Rimossi Delta', color='yellow')

plt.xlabel('Tempo (giorni)')
plt.ylabel('Numero di individui')
plt.title('Algoritmo di Gillespie per Modello SIR con Varianti Delta e Omicron')
plt.legend()
plt.grid(True)
plt.show()

# Calcolo della probabilità di infezione con Omicron rispetto a Delta
P_omicron = np.array(I2_values) / (np.array(I1_values) + np.array(I2_values))

plt.figure(figsize=(12, 8))
plt.plot(time_points, P_omicron, label='Probabilità di infezione con Omicron', color='green')
plt.xlabel('Tempo (giorni)')
plt.ylabel('Probabilità')
plt.title('Probabilità di infezione con Omicron nel tempo')
plt.legend()
plt.grid(True)
plt.show()

# calcolo della probabilità di infezione con Omicron rispetto a Delta AL GIORNO
P_omicron_daily = np.where(new_delta_infections + new_omicron_infections > 0,
                           new_omicron_infections / (new_delta_infections + new_omicron_infections),
                           0)

plt.figure(figsize=(12, 8))
plt.plot(time_points_daily, P_omicron_daily, label='Probabilità di infezione con Omicron', color='green')
plt.xlabel('Tempo (giorni)')
plt.ylabel('Probabilità')
plt.title('Probabilità giornaliera di infezione con Omicron rispetto a Delta al giorno')
plt.legend()
plt.grid(True)
plt.show()




