%Cria um dataset para a condição de LoS com várias BSs
%Autor: Paulo Francisco
clc;clear;close all
Fc=28e9;
c=300;
b=[
    30, 10, 3;
    50, 10, 3;
    30, 30, 3;
    50, 30, 3;
    10, 10, 3;
    70, 10, 3;
    10, 30, 3;
    70, 30, 3;
    ]';
m=[40, 20, 1]';

% mostra(b,m)

los=1;
Nt=64;              % number of TX antennas
Nr=Nt;              % number of RX antennas
N=15;               % number of subcarriers
B=100;              % total BW in MHz
Ns=20;              % number of sybols sent
arranjo=1;          %1-URA, 2-UCA

toas=buscaToAs (b,m)/c;
aods=buscaAoDsLos(b,m);

rng(1);
K=8;
sn=-20:5:20;
nruns=100;
for k=1:K
    toa=toas(k);
    aod=aods(k,:);
    aoa=[0.2618 1.6581];
    for i=1:size(sn,2)
        snr=sn(i);
        for j=1:nruns
            [H,y,w,yp] = channelModeling(toa, aod, aoa, Nr, Nt, B, N,Fc,Ns,snr,arranjo);
            [ToAs_BfLoS(k,j,i), AoDs_BfLoS(k,:,j,i), ~,~]=DCS_SOMP_Adaptativo(Nt, Nr, Ns, N, B, c, y, w, 1,Fc,arranjo,6,0);
        end
    end
end
save("dataSet_apendice_9Ant.mat","ToAs_BfLoS","AoDs_BfLoS","sn");