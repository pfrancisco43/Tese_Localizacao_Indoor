%Gera gráficos comparativos de RMSE para as propostas em comparação com
%métodos da literatura - Variação da quantidade de caminhols
%autor: Paulo Francisco
clc;clear;close all
load dataSet_CaminhoLoS
load dataSet_4Caminhos_Nlos
%%%%%%só pra saber o que foi usado no dataSet
%Devem estar na mesma ordem da estimação
b=[-8, 0, 5]';
m=[7, 10, 1]';
s1=[-10, 4, 3];
s2=[10, 8, 3];
s3=[-10, 8, 4];
s4=[10, 4, 4];
S=[s1;s2;s3;s4]';
c=300;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nruns=size(ToAs_BfLoS,2);

cam=2:4;
snr=1;
qnt=4;

for i=2:qnt
    s=S(:,1:i);
    a=reshape(ToAs_Bf(1:i,:,snr),[],1);
    s_toas=sqrt(var(a)); %Variância na estimação de ToA
    a=reshape(AoDs_Bf(1:i,:,:,snr),[],1);
    s_angs=sqrt(var(a)); %Variância na estimação dos angulos
    for j=1:Nruns % QUANTIDADE DE RUNS
        toaL=ToAs_BfLoS(1,j,snr)*c;
        aodL=ADos_BfLoS(1,:,j,snr);
        aoaL=AoAs_BfLoS(1,:,j,snr);

        toasN=ToAs_Bf(1:i,j,snr)*c;
        aodsN=AoDs_Bf(1:i,:,j,snr);
        aoasN=AoAs_Bf(1:i,:,j,snr);
        aoasN=aod_aoa_Swap(aoasN);

        [eL(j),~]=estimadorToAAoD(b,toaL,aodL,m);

        [e1(j),eh,ev,ss,rr,kk,posi,posf,esi,esf,e2(j),cr,ee1(j),ee2(j)]=estimadorIntersecLinhasNovo(toasN,aoasN,aodsN,b,m,s,s_angs,s_toas);
    
        [es(j)]=Shikur.Shikur(toasN, aodsN, aoasN, b, m, s, s_toas, s_angs);
        [~,~,ew(j),ewe(j)]=AlgoritmoHenk.main(b,s,m,s_toas,s_angs);
        et(j)=tobias.tobias(toasN,aoasN,aodsN,b,s,m);
    end
    etL(i-1)=sqrt(mean(eL.^2));
    et1(i-1)=sqrt(mean(e1.^2));
    et2(i-1)=sqrt(mean(e2.^2));
    ets(i-1)=sqrt(mean(es.^2));
    etw(i-1)=sqrt(mean(ew.^2));
    ett(i-1)=sqrt(mean(et.^2));

    etee1(i-1)=sqrt(mean(ee1.^2));
    etee2(i-1)=sqrt(mean(ee2.^2));
    etewe(i-1)=sqrt(mean(ewe.^2));
end
figure
hold on
semilogy(cam,etL,'bX--');
semilogy(cam,et1,'bX-');
semilogy(cam,et2,'d-','Color','#77AC30');
semilogy(cam,ets,'r^-');
semilogy(cam,etw,'ms-');
semilogy(cam,ett,'ko-');
xticks([2 3 4])
legend('Método LoS','Proposta 1','Proposta 2','(Shikur; Weber, 2014)','(Wymeersch, 2018)','(Wei; Palleit; Weber,2011)');
xlabel('Quantidade de Caminhos NLoS')
ylabel('Localização do UE - RMSE (m)');
grid on;

figure
hold on
semilogy(cam,etee1,'bX-');
semilogy(cam,etee2,'d-','Color','#77AC30');
semilogy(cam,etewe,'ms-');
xticks([2 3 4])
legend('Proposta 1','Proposta 2','(Wymeersch, 2018)');
xlabel('Quantidade de Caminhos NLoS')
ylabel('Localização dos SCs - RMSE (m)');
grid on;