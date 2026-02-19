# Simulazioni SIR Variant COVID-19

Questo progetto analizza la diffusione di due varianti del virus COVID-19 tramite il modello epidemiologico SIR. Il lavoro confronta modelli deterministici e stocastici per studiare il passaggio di dominanza tra varianti

## Contenuto della cartella

I codici Python implementano gli algoritmi descritti della tesi:
- **Metodo di Eulero**: Simulazione deterministica con passo *dt = 0.1*.
- **Rounge-Kutta (RK4)**: Modello deterministico ottimizzato per efficienza e precisione.
- **Algoritmo di Gillespie**: Simulazione stocastica che modella le probabilità di infezione casuali.

## Risultati principali
- La variante 2 diventa predominante se il suo numero di riproduzione $R_0$ è molto maggiore della variante 1, anche se quest'ultima è ancora in crescita.
- I modelli stocastici evidenziano come la probabilità di infettarsi con la variante 2 aumenti significativamente nel tempo.

## Requisiti

- Python 3.x
- Librerie: `numpy`, `matplotlib`
