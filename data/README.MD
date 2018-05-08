
Excel files explanation:


GOM_original.xlsx
Original data, sent to me by Lukas Jonkers on the 18.12.2017 


GOM_original_sst.csv
Original da with extra columns about the timing of each sampling period (see below), plus Sea Surface Temperature data (sst). The problem of ’GOM_original_sst.csv' is that the time-steps have different lengths (see column "resolution": sometimes 7 days and sometimes 14 days).
New columns:
open: day sediment trap sample opened
close: day sediment trap sample closed
middle: day in between open and close
doy: Day Of the Year of the middle date (i.e., 2/2/2018 = 33)
cum_days: cumulative days of the whole sampling period
resolution: number of days in between open and close
next_open: number of days in between two samplings events (i.e., two rows) = in between the last close and the next open)
set: Sea Surface Temperature data - mean sst for the days in between open and close. Data from http://monitor.cicsnc.org/obs4MIPs/data/OISST/Daily/

GOM_NAs.csv
GOM_original_sst.csv spreadsheet, but now all time-steps are equivalent to 7 days (added column: time-step). 
Two types of NAs: 
(1) ‘NA_gap' são os NAs verdadeiros (i.e., período de tempo que não temos amostras). 
(2) ’NA_resolution' são amostras que foram coletadas num periodo de 14 dias (em vez de 7 dias como a maioria da series). Então, para colocar cada linha da matriz numa mesma resolução de 7 dias (mesmo tamanho de time-steps), eu dividi as amostras de 14 dias em duas de 7, e coloquei NA_resolution em uma delas. Porem, nós temos dados coletados neste NA_resolution: a linha de cima (que foi coletada no mesmo período - mesmo 'open' e 'close'). Então em vez de deixar o NA no NA_resolution, o melhor seria simplesmente copiar a linha de cima. O Brenno disse que talvez isso seja problemático pro CCM, pois teremos algumas linhas na tabelas que serão idênticas. 

GOM_NA_only.csv 
Adicionada pelo Rafa, com todos os NA_gap e NA_resolution da GOM_NAs.csv como NA.

GOM_NA_gap.csv
Adicionada pela Marina: Na_gap = NA, mas NA_resolution preenchidos com a linha de cima (linhas idênticas) 



———————

Mais explicação sobre os NA_gap e NA_resolution (da planilha GOM_NAs.csv)

"NA_gap": são os casos em que a garrafinha foi perdida, um gap 'verdadeiro' pois não temos nenhum dado para este time_step 
Solução: a maioria dos NA_gap são apenas um time-step (interpolação linear provavelmente ok). Porem temos um caso de NA_gap para 6 time-steps seguidos (linhas 98 ate 103) - como resolver? talvez usar "Expectation–maximization algorithm".

"NA_resolution": são os casos nos quais a resolução amostrada foi maior que a estabelecida para a analise EDM. Nós estabelecemos uma resolução de 7 dias (cada time-step equivale a 7 dias). Algumas amostras são valores médios para 14 dias. Teremos então que dividir esta amostra entre duas linhas de 7 dias. Eu ja dupliquei cada amostra de 14 dias e a linha abaixo de cada uma eu preenchi com "NA_resolution" (note que as datas de amostragem se repetem 2 a 2). 
Solução: uma solução simples eh assumir que o fluxo de conchas eh constante para o intervalo de 14 dias, e então colocar valores idênticos nas duas linhas que dividem o intervalo de 14 dias em dois de 7 dias. Porem, olhando os dados, da pra ver que as abundância variam bastante de uma semana para outra.  
Outra solução eh talvez estimar para cada espécie uma distribuição de abundância no tempo (contar quantos time-steps tem abundância x, y, z) e então com essa informação sortear (?) como a abundância de 14 dias se distribuira em dois intervalos de 7 dias.
