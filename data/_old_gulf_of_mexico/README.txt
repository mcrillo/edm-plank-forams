
GOM_original.xlsx
Original data, sent to me by Lukas Jonkers on the 18.12.2017 


GOM_temp.csv
sampling_start : day the bottle (sample) was opened
sampling_end : day the bottle (sample) was closed
resolution_days : number of days sampled (= days between sampling_end and sampling_start)
time_step : time step for EDM use (sequencial 1 to 203, resolution 7 days)
G_menardii - G_falconensis: species, values are in (number of shells)*(meter^-2)*(day-1)
temp_K : mean of the daily sea surface temperature for the sampling interval (sampling_end - sampling_start). Data from http://monitor.cicsnc.org/obs4MIPs/data/OISST/Daily/
temp_oC : temperature in Celsius (= temp_K - 273.15)

Dados de abundancias (fluxos de conchas) das especies:

"NA_gap": são os casos em que a garrafinha foi perdida, um gap 'verdadeiro' pois não temos nenhum dado para este time_step 
Solução: a maioria dos NA_gap são apenas um time-step (interpolação linear provavelmente ok). Porem temos um caso de NA_gap para 6 time-steps seguidos (linhas 98 ate 103) - como resolver? talvez usar "Expectation–maximization algorithm".

"NA_resolution": são os casos nos quais a resolução amostrada foi maior que a estabelecida para a analise EDM. Nós estabelecemos uma resolução de 7 dias (cada time-step equivale a 7 dias). Algumas amostras são valores médios para 14 dias. Teremos então que dividir esta amostra entre duas linhas de 7 dias. Eu ja dupliquei cada amostra de 14 dias e a linha abaixo de cada uma eu preenchi com "NA_resolution" (note que as datas de amostragem se repetem 2 a 2). 
Solução: uma solução simples eh assumir que o fluxo de conchas eh constante para o intervalo de 14 dias, e então colocar valores idênticos nas duas linhas que dividem o intervalo de 14 dias em dois de 7 dias. Porem, olhando os dados, da pra ver que as abundância variam bastante de uma semana para outra.  
Outra solução eh talvez estimar para cada espécie uma distribuição de abundância no tempo (contar quantos time-steps tem abundância x, y, z) e então com essa informação sortear (?) como a abundância de 14 dias se distribuira em dois intervalos de 7 dias.



plots folder:
plot das series temporais de todas as espécies juntas
plots das series temporais para cada espécie:
- fluxo 
- log do fluxo
- first difference do fluxo