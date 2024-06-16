clc;clear;close all
Fc=28e9;
c=300;
it=4;
b=[-8, 0, 5]';
m=[7, 10, 1]';
s1=[-10, 4, 3];
s2=[10, 8, 3];
s3=[-10, 8, 4];
s4=[10, 4, 4];
s=[s1;s2;s3;s4]';
los=0;
Nt=64;            % number of TX antennas
Nr=Nt;             % number of RX antennas
N=10;               % number of subcarriers
B=100;             % total BW in MHz
Ns=20;             % number of sybols sent
arranjo=1;         %1-URA, 2-UCA

%% Determinação de parâmetros reias
[ToAs, AoDs, AoAs, L]=getTrueParameters(b,m,s,los,c);
AoAs=aod_aoa_Swap(AoAs);

%Selecionando o caminho 3
ToAs=ToAs(3);
AoDs=AoDs(3,:);
AoAs=AoAs(3,:);
L=1;

alpha=[ToAs,AoDs,AoAs];

sn=-20:1:0;
rng(1);
Nruns=1000;
empilha=0;
for i=1:size(sn,2)
    snr=sn(i);
    for j=1:Nruns % QUANTIDADE DE RUNS
        [H,y,w,yp] = channelModeling(ToAs, AoDs, AoAs, Nr, Nt, B, N,Fc,Ns,snr,arranjo);

        % Estimação de parâmetros com DCS-SOMP Simples        
        [ToAs_Dc, AoDs_Dc, AoAs_Dc,~]=DCS_SOMP_Simples(Nt, Nr, Ns, N, B, c, y, w, L,Fc,arranjo);
        e_t1(j)=norm(ToAs_Dc-ToAs);
        e_aodaz1(j)=norm(AoDs_Dc(1,1)-AoDs(1,1));
        e_aodel1(j)=norm(AoDs_Dc(1,2)-AoDs(1,2));
        e_aoaaz1(j)=norm(AoAs_Dc(1,1)-AoAs(1,1));
        e_aoael1(j)=norm(AoAs_Dc(1,2)-AoAs(1,2));

        %Estimação de parâmetros com DCS-SOMP Adaptativo
        [ToAs_Ad, AoDs_Ad, AoAs_Ad,~]=DCS_SOMP_Adaptativo(Nt, Nr, Ns, N, B, c, y, w, L,Fc,arranjo,it,empilha);
        e_t2(j)=norm(ToAs_Ad(end)-ToAs);
        e_aodaz2(j)=norm(AoDs_Ad(end,1)-AoDs(1,1));
        e_aodel2(j)=norm(AoDs_Ad(end,2)-AoDs(1,2));
        e_aoaaz2(j)=norm(AoAs_Ad(end,1)-AoAs(1,1));
        e_aoael2(j)=norm(AoAs_Ad(end,2)-AoAs(1,2));

        %Estimação de parâmetros com DCS-SOMP e beamforming
        [ToAs_Bf, AoDs_Bf, AoAs_Bf,wo]=DCS_SOMP_beam(Nt, Nr, Ns, N, B, c, y, w, L,Fc,arranjo,snr,alpha,yp,H,empilha);
        e_t3(j)=norm(ToAs_Bf(end)-ToAs);
        e_aodaz3(j)=norm(AoDs_Bf(end,1)-AoDs(1,1));
        e_aodel3(j)=norm(AoDs_Bf(end,2)-AoDs(1,2));
        e_aoaaz3(j)=norm(AoAs_Bf(end,1)-AoAs(1,1));
        e_aoael3(j)=norm(AoAs_Bf(end,2)-AoAs(1,2));
    end
    et_t1(i)=sqrt(mean(e_t1.^2));
    et_aodaz1(i)=sqrt(mean(e_aodaz1.^2));
    et_aodel1(i)=sqrt(mean(e_aodel1.^2));
    et_aoaaz1(i)=sqrt(mean(e_aoaaz1.^2));
    et_aoael1(i)=sqrt(mean(e_aoael1.^2));

    et_t2(i)=sqrt(mean(e_t2.^2));
    et_aodaz2(i)=sqrt(mean(e_aodaz2.^2));
    et_aodel2(i)=sqrt(mean(e_aodel2.^2));
    et_aoaaz2(i)=sqrt(mean(e_aoaaz2.^2));
    et_aoael2(i)=sqrt(mean(e_aoael2.^2));

    et_t3(i)=sqrt(mean(e_t3.^2));
    et_aodaz3(i)=sqrt(mean(e_aodaz3.^2));
    et_aodel3(i)=sqrt(mean(e_aodel3.^2));
    et_aoaaz3(i)=sqrt(mean(e_aoaaz3.^2));
    et_aoael3(i)=sqrt(mean(e_aoael3.^2));
end

%% plots - ToA
figure
hold on
semilogy(sn,et_t1,'b*-')
semilogy(sn,et_t2,'r^--')
semilogy(sn,et_t3,'mo--')
legend('DCS-SOMP','DCS-SOMP Adaptativo','Beamforming Adaptativo');
xlabel('SNR')
ylabel('ToA - RMSE ($\mu$s)','Interpreter','latex')
grid on;

%% plots - AoD az
figure
hold on
semilogy(sn,et_aodaz1*180/pi,'b*-')
semilogy(sn,et_aodaz2*180/pi,'r^--')
semilogy(sn,et_aodaz3*180/pi,'mo--')
legend('DCS-SOMP','DCS-SOMP Adaptativo','Beamforming Adaptativo');
xlabel('SNR')
ylabel('AoD$^{az}$ - RMSE ($^\circ$)','Interpreter','latex')
grid on;

%% plots - AoD el
figure
hold on
semilogy(sn,et_aodel1*180/pi,'b*-')
semilogy(sn,et_aodel2*180/pi,'r^--')
semilogy(sn,et_aodel3*180/pi,'mo--')
legend('DCS-SOMP','DCS-SOMP Adaptativo','Beamforming Adaptativo');
xlabel('SNR')
ylabel('AoD$^{el}$ - RMSE ($^\circ$)','Interpreter','latex')
grid on;

%% plots - AoA az
figure
hold on
semilogy(sn,et_aoaaz1*180/pi,'b*-')
semilogy(sn,et_aoaaz2*180/pi,'r^--')
semilogy(sn,et_aoaaz3*180/pi,'mo--')
legend('DCS-SOMP','DCS-SOMP Adaptativo','Beamforming Adaptativo');
xlabel('SNR')
ylabel('AoA$^{az}$ - RMSE ($^\circ$)','Interpreter','latex')
grid on;

%% plots - AoA el
figure
hold on
semilogy(sn,et_aoael1*180/pi,'b*-')
semilogy(sn,et_aoael2*180/pi,'r^--')
semilogy(sn,et_aoael3*180/pi,'mo--')
legend('DCS-SOMP','DCS-SOMP Adaptativo','Beamforming Adaptativo');
xlabel('SNR')
ylabel('AoA$^{el}$ - RMSE ($^\circ$)','Interpreter','latex')
grid on;
%save("simulacao1.mat","et_t1","et_t2","et_t3","et_aodaz1","et_aodaz2","et_aodaz3","et_aoaaz1","et_aoaaz2","et_aoaaz3","et_aodel1","et_aodel2","et_aodel3","et_aoael1","et_aoael2","et_aoael3")

