%Gera a CDF e gera a CDF para um cenário com a probabilidade de LoS
%determinada segundo a 3GPP 38901
clc;clear;close all
%busca das bases
load dataSet_CaminhoLoS %Base los beam
load dataSet_4Caminhos_Nlos %base com 4 caminhos nlos beam

%%%%%%só pra saber o que foi usado no dataSet
b=[-8, 0, 5]';
m=[7, 10, 1]';
s1=[-10, 4, 3];
s2=[10, 8, 3];
s3=[-10, 8, 4];
s4=[10, 4, 4];
s=[s1;s2;s3;s4]';
c=300;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nruns=1000;
tic
prob=DeterminaProbLos(b,m); %probabilidade de LoS
%prob=0.6;
ql=0;
for i=1:Nruns
    los=rand < prob; %Será LoS?
    snr=randi([1 9],1,1); %Aleatório para snr
    j=randi([1 100],1,1); %Aleatório qual medição do dataset

    a=reshape(ToAs_Bf(:,:,snr),[],1);
    s_toas=sqrt(var(a)); %Variância na estimação de ToA
    a=reshape(AoDs_Bf(:,:,:,snr),[],1);
    s_angs=sqrt(var(a)); %Variância na estimação dos angulos

    if los
        toaL=ToAs_BfLoS(:,j,snr)*c + sqrt(0.01^2)*randn;
        aodL=ADos_BfLoS(:,:,j,snr)+ sqrt(0.01^2)*randn;
        [eT(i),eV(i),eH(i)]=estimadorToAAoD(b,toaL,aodL,m);
        ql=ql+1;
    else
        toasN=ToAs_Bf(:,j,snr)*c+ sqrt(0.01^2)*randn;
        aodsN=AoDs_Bf(:,:,j,snr)+ sqrt(0.01^2)*randn;
        aoasN=AoAs_Bf(:,:,j,snr)+ sqrt(0.01^2)*randn;
        aoasN=aod_aoa_Swap(aoasN);
        [e1,eH(i),eV(i),ss,rr,kk,posi,posf,esi,esf,eT(i),cr,ee1,ee2]=estimadorIntersecLinhasNovo(toasN,aoasN,aodsN,b,m,s,s_angs,s_toas);

    end
end
toc
eT=sort(eT);
eH=sort(eH);
eV=sort(eV);

sqrt(mean(eT.^2));
sqrt(mean(eH.^2));
sqrt(mean(eV.^2));
ql
eT(Nruns*0.95)
eH(Nruns*0.95)
eV(Nruns*0.95)

y=linspace(0,1,Nruns);%para o gráfico apenas
semilogx(eT,y,'-k','LineWidth',2)
hold on
semilogx(eH,y,'-b','LineWidth',2)
semilogx(eV,y,'-m','LineWidth',2)
semilogx([min(eT),max(eT)],[0.95,0.95],'--','Color','red','LineWidth',2)
grid on;
xlabel('Erro (m)')
ylabel('CDF')
legend('Erro Conjunto','Erro Horizontal','Erro Vertical','Acurácia (95%)');
