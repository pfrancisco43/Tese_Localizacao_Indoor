%Compara os métodos sob condição de LoS e várias BSs
%RMSE em relação a variação da SNR
clc;close all; clear;
load dataSet_apendice
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

L=5;
cs=[1 2 3 5 7]';
bs=b(:,cs);

b=b(:,cs);
nruns=100;
for snr=1:size(sn,2)  
    a=reshape(ToAs_BfLoS(cs,:,snr),[],1);
    s_toas=sqrt(var(a)); %Variância na estimação de ToA    
    for j=1:nruns
        toas=ToAs_BfLoS(cs,j,snr)*c;       
        aods=AoDs_BfLoS(cs,:,j,snr);       
        
        [eA(j),eP(j)]=estimadorToAAoD(b,toas,aods,m);
        [eEGC(j)]=estimadorEGC_TOA (m,b,toas);
        [eTaylor(j)]=estimadorTaylorToA(m,b,toas,s_toas);
        [eIMTL(j)]= estimadorIMTL(m,b,toas);
        
    end
    meA(snr)=sqrt(mean(eA.^2));
    meP(snr)=sqrt(mean(eP.^2));
    meEGC(snr)=sqrt(mean(eEGC.^2));
    meTaylor(snr)=sqrt(mean(eTaylor.^2));
    meIMTL(snr)=sqrt(mean(eIMTL.^2));    
end
vetor=sn;
semilogy(vetor,meA,'^-r')
hold on;
semilogy(vetor,meP,'s-m')
semilogy(vetor,meEGC,'o-b')
semilogy(vetor,meTaylor,'*-k')
semilogy(vetor,meIMTL,'x-','Color','#77AC30')

%ylim([-1 15])
grid on;
xlabel('SNR (dB)')
ylabel('RMSE (m)')
legend('Proposta - ToA + AoD e Média Simples', 'Proposta - ToA + AoD e Média Ponderada','Proposta - EGC','Método de Taylor','Método IMTL');