# Come leggere i file

Ho sonno e non mi stava venendo in mente nulla di più furbo

Il file CV_data.root ha un tree per ogni pixel testato, le branch sono V, C e i relativi errori, ogni leaf è un punto

Il file IV_data.root ha un tree per ogni pixel testato, le branch sono T, V, I e gli errori su I e V, ogni leaf è un punto di temperatura e un vector per il resto contenente tutti i dati a quella temperatura. Il vector è ordinato quindi non dovremmo avere problemi