%Cria um dataset contendo estimativas de parâmetros dos caminhos NLoS
%@Autor: Paulo Francisco
clc;clear;close all
Fc=28e9;
c=300;
b=[-8, 0, 5]';
m=[7, 10, 1]';
s1=[-10, 4, 3];
s2=[10, 8, 3];
s3=[-10, 8, 4];
s4=[10, 4, 4];
s=[]';
los=1;
Nt=64;              % number of TX antennas
Nr=Nt;              % number of RX antennas
N=10;               % number of subcarriers
B=100;              % total BW in MHz
Ns=20;              % number of sybols sent
arranjo=1;          %1-URA, 2-UCA

%% Determinação de parâmetros reais
[ToAs, AoDs, AoAs, L]=getTrueParameters(b,m,s,los,c);

AoAs=aod_aoa_Swap(AoAs);
[ToAs, AoDs, AoAs]=orderParameters(ToAs, AoDs, AoAs);

alpha=[ToAs,AoDs,AoAs];
empilha=0; %Para mostrar resultados de cada iteração

rng(1);
sn=-20:5:20;
Nruns=1000;
tic
for i=1:size(sn,2)
    snr=sn(i);    
    for j=1:Nruns % QUANTIDADE DE RUNS        
        [H,y,w,yp] = channelModeling(ToAs, AoDs, AoAs, Nr, Nt, B, N,Fc,Ns,snr,arranjo);        
        %Estimação de parâmetros com Beamforming Adaptativo
        [ToAs_BfLoS(:,j,i), ADos_BfLoS(:,:,j,i), AoAs_BfLoS(:,:,j,i),wo]=DCS_SOMP_beam(Nt, Nr, Ns, N, B, c, y, w, L,Fc,arranjo,6,empilha);        
    end
end
tocsave("dataSet_4Caminhos.mat","ToAs_BfLoS","ADos_BfLoS","AoAs_BfLoS","sn");
