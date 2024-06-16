%Compara os métodos sob condição de LoS e várias BSs
%RMSE em relação a variação da Qdt de BSs
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
cs=[1 2 3 4 5 6 7 8]';
bs=b(:,cs);

nruns=100;
snr=5;
bo=b;
for i=1:8
    cs=1:i;
    b=bo(:,cs);
    a=reshape(ToAs_BfLoS(cs,:,snr),[],1);
    s_toas=sqrt(var(a)); %Variância na estimação de ToA
    for j=1:nruns
        toas=ToAs_BfLoS(cs,j,snr)*c;
        aods=AoDs_BfLoS(cs,:,j,snr);        

        [eA(j),eP(j)]=estimadorToAAoD(b,toas,aods,m);
        if i>2
            [eEGC(j)]=estimadorEGC_TOA (m,b,toas);
            [eTaylor(j)]=estimadorTaylorToA(m,b,toas,s_toas);
            [eIMTL(j)]= estimadorIMTL(m,b,toas);
        end

    end
    meA(i)=sqrt(mean(eA.^2));
    meP(i)=sqrt(mean(eP.^2));
    if i>2
        meEGC(i-2)=sqrt(mean(eEGC.^2));
        meTaylor(i-2)=sqrt(mean(eTaylor.^2));
        meIMTL(i-2)=sqrt(mean(eIMTL.^2));
    end
end
vetor=3:8;
semilogy([1:8],meA,'^-r')
hold on;
semilogy(1:8,meP,'s-m')
semilogy(vetor,meEGC,'o-b')
semilogy(vetor,meTaylor,'*-k')
semilogy(vetor,meIMTL,'x-','Color','#77AC30')

%ylim([-1 15])
grid on;
xlabel('Quantidade de Caminhos')
ylabel('RMSE (m)')
legend('Proposta - ToA + AoD e Média Simples', 'Proposta - ToA + AoD e Média Ponderada','Proposta - EGC','Método de Taylor','Método IMTL');