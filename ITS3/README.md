# ITS3 - Formato dati

I dati acquisiti con la sorgente di ⁵⁵Fe (limitatamente a quelli utilizzati nell'analisi) sono contenuti nel file `data.root`.

Il file contiene 16 `TTree`, ciascuno chiamato `pixelN`, con `N = 0, 1, ..., 15`.

## Pixel centrali (oscilloscopio)

I `TTree` relativi ai pixel centrali, letti con l'oscilloscopio, hanno le seguenti branch:

| Branch          | Tipo                         |
|-----------------|------------------------------|
| `baseline`      | `double`                     |
| `baseline_rms`  | `double`                     |
| `baseline_err`  | `double`                     |
| `underline`     | `double`                     |
| `underline_rms` | `double`                     |
| `underline_err` | `double`                     |
| `signal`        | `double`                     |
| `CFD_times`     | struct di 10 `double`        |

Le branch contengono un valore `double` per ogni evento, ad eccezione di `CFD_times`, che è una struct di 10 valori `double` (`t0`, `t10`, ..., `t90`).

Per semplicità, il nome della branch che permette di accedere al singolo valore di tempo è:

```
CFD_times.tN    (con N = 0, 10, ..., 90)
```

## Pixel esterni (DAQ board)

I `TTree` relativi ai pixel esterni, letti con la DAQ board, hanno le seguenti branch di valori `double`:

| Branch     | Tipo     |
|------------|----------|
| `baseline` | `double` |
| `rms`      | `double` |
| `p2p`      | `double` |