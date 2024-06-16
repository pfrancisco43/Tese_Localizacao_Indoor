%Executa uma execução do DCS-SOMP e do DCS-SOMP, mostrando o resultado de
%cada iteração
%Autor: Paulo Francisco
clc;clear;close all
Fc=28e9;
c=300; % em us
b=[-8, 0, 5]';
m=[7, 10, 1]';
s=[]';
los=1;
Nt=64;              % number of TX antennas
Nr=Nt;              % number of RX antennas
N=10;               % number of subcarriers
B=100;              % total BW in MHz
Ns=20;              % number of sybols sent
arranjo=1;          %1-URA, 2-UCA

%% Determinação de parâmetros reias
[ToAs, AoDs, AoAs, L]=getTrueParameters(b,m,s,los,c);
AoAs=aod_aoa_Swap(AoAs);

%% Modelagem do canal
snr=10;
rng(1);
[H,y,w,yp] = channelModeling(ToAs, AoDs, AoAs, Nr, Nt, B, N,Fc,Ns,snr,arranjo);

%% Estimação de parâmetros com DCS-SOMP Simples
[ToAs_Dc, AoDs_Dc, AoAs_Dc,~]=DCS_SOMP_Simples(Nt, Nr, Ns, N, B, c, y, w, L,Fc,arranjo);
AoAs_Dc=aod_aoa_Swap(AoAs_Dc);

%% Estimação de parâmetros com DCS-SOMP Adaptativo
empilha=1;
it=10;
[ToAs_Ad, AoDs_Ad, AoAs_Ad,~]=DCS_SOMP_Adaptativo(Nt, Nr, Ns, N, B, c, y, w, L,Fc,arranjo,it,empilha);
AoAs_Ad=aod_aoa_Swap(AoAs_Ad);

AoAs=aod_aoa_Swap(AoAs);

AoDs=rad2deg(AoDs);
AoAs=rad2deg(AoAs);

AoDs_Dc=rad2deg(AoDs_Dc);
AoAs_Dc=rad2deg(AoAs_Dc);
AoDs_Ad=rad2deg(AoDs_Ad);
AoAs_Ad=rad2deg(AoAs_Ad);
%% plots - ToA
figure
t=length(ToAs_Ad);
semilogy([1:t],ones(1,t)*ToAs,'--k','linewidth',2)
hold on
plot([1:t],ones(1,t)*ToAs_Dc,'--b','linewidth',1)
plot([1:t],ToAs_Ad,'-^r','linewidth',1)
xlabel('Iteração')
ylabel('ToA ($\mu$s)','Interpreter','latex')
grid on
legend('Valor real','DCS-SOMP','DCS-SOMP Adaptativo')
 yticks([min(ToAs_Ad) max(ToAs_Ad)])

%% plots - AoD-Az
figure
semilogy([1:t],ones(1,t)*AoDs(1,1),'--k','linewidth',2)
hold on
semilogy([1:t],ones(1,t)*AoDs_Dc(1,1),'--b','linewidth',1)
semilogy([1:t],AoDs_Ad(:,1),'-^r','linewidth',1)
xlabel('Iteração')
ylabel('AoD$^{az}$ ($^\circ$)','Interpreter','latex')
grid on
legend('Valor real','DCS-SOMP','DCS-SOMP Adaptativo')

%% plots - AoD-El
figure
semilogy([1:t],ones(1,t)*AoDs(1,2),'--k','linewidth',2)
hold on
semilogy([1:t],ones(1,t)*AoDs_Dc(1,2),'--b','linewidth',1)
semilogy([1:t],AoDs_Ad(:,2),'-^r','linewidth',1)
xlabel('Iteração')
ylabel('AoD$^{el}$ ($^\circ$)','Interpreter','latex')
grid on
legend('Valor real','DCS-SOMP','DCS-SOMP Adaptativo')

%% plots - AoA-Az
figure
semilogy([1:t],ones(1,t)*AoAs(1,1),'--k','linewidth',2)
hold on
semilogy([1:t],ones(1,t)*AoAs_Dc(1,1),'--b','linewidth',1)
semilogy([1:t],AoAs_Ad(:,1),'-^r','linewidth',1)
xlabel('Iteração')
ylabel('AoA$^{az}$ ($^\circ$)','Interpreter','latex')
grid on
legend('Valor real','DCS-SOMP','DCS-SOMP Adaptativo')

%% plots - AoA-El
figure
semilogy([1:t],ones(1,t)*AoAs(1,2),'--k','linewidth',2)
hold on
semilogy([1:t],ones(1,t)*AoAs_Dc(1,2),'--b','linewidth',1)
semilogy([1:t],AoAs_Ad(:,2),'-^r','linewidth',1)
xlabel('Iteração')
ylabel('AoA$^{el}$ ($^\circ$)','Interpreter','latex')
grid on
legend('Valor real','DCS-SOMP','DCS-SOMP Adaptativo')

