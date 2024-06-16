%Compara os métodos sob condição de LoS e várias BSs
%CDF
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
cs=[1 2 3 4 5 6 7]';
bs=b(:,cs);

b=b(:,cs);
nruns=100;
kk=1;
for snr=[2 4:6]
    a=reshape(ToAs_BfLoS(cs,:,snr),[],1);
    s_toas=sqrt(var(a)); %Variância na estimação de ToA    
    for j=1:nruns
        toas=ToAs_BfLoS(cs,j,snr)*c;       
        aods=AoDs_BfLoS(cs,:,j,snr);        
        
        [eA(kk),eP(kk)]=estimadorToAAoD(b,toas,aods,m);
        [eEGC(kk)]=estimadorEGC_TOA (m,b,toas);
        [eTaylor(kk)]=estimadorTaylorToA(m,b,toas,s_toas);
        [eIMTL(kk)]= estimadorIMTL(m,b,toas); 
        kk=kk+1;
    end
end
sqrt(mean(eA).^2)
sqrt(mean(eP).^2)
sqrt(mean(eEGC).^2)
sqrt(mean(eTaylor).^2)
sqrt(mean(eIMTL).^2)

eA(nruns*0.95)
eP(nruns*0.95)
eEGC(nruns*0.95)
eTaylor(nruns*0.95)
eIMTL(nruns*0.95)

eA=sort(eA);
eP=sort(eP);
eEGC=sort(eEGC);
eTaylor=sort(eTaylor);
eIMTL=sort(eIMTL);

nruns=kk-1;
y=linspace(0,1,nruns);%para o gráfico apenas
semilogx(eA,y,'-r','LineWidth',2)
hold on;
semilogx(eP,y,'-m','LineWidth',2)
semilogx(eEGC,y,'-b','LineWidth',2)
semilogx(eTaylor,y,'-k','LineWidth',2)
semilogx(eIMTL,y,'-','LineWidth',2,'Color','#77AC30')

semilogx([min(eA),max(eEGC)],[0.95,0.95],'--','Color','#4DBEEE','LineWidth',2)
grid on;
xlabel('Erro (m)')
ylabel('CDF')
legend('Proposta - ToA + AoD e Média Simples', 'Proposta - ToA + AoD e Média Ponderada','Proposta - EGC','Método de Taylor','Método IMTL','Probabilidade em 95%');
