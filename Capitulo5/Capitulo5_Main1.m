%Gera gráficos comparativos de RMSE para as propostas em comparação com
%métodos da literatura - Variação da SNR
%autor: Paulo Francisco
clc;clear;close all
%busca das bases
load dataSet_CaminhoLoS %Base los com Beamforming Adaptativo
load dataSet_4Caminhos_Nlos %base com 4 caminhos nlos com Beamforming Adaptativo
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

Nruns=100;
for i=1:size(sn,2) 
    a=reshape(ToAs_Bf(:,:,i),[],1);
    s_toas=sqrt(var(a)); %Variância na estimação de ToA
    a=reshape(AoDs_Bf(:,:,:,i),[],1);
    s_angs=sqrt(var(a)); %Variância na estimação dos angulos
    for j=1:Nruns % QUANTIDADE DE RUNS
        toaL=ToAs_BfLoS(:,j,i)*c;
        aodL=ADos_BfLoS(:,:,j,i);

        toasN=ToAs_Bf(:,j,i)*c;
        aodsN=AoDs_Bf(:,:,j,i);
        aoasN=AoAs_Bf(:,:,j,i);
        aoasN=aod_aoa_Swap(aoasN);%Pegar o oposto do AoA, na base está invertido

        [eL(j),~]=estimadorToAAoD(b,toaL,aodL,m);

        [e1(j),eh,ev,ss,rr,kk,posi,posf,esi,esf,e2(j),cr,ee1(j),ee2(j)]=estimadorIntersecLinhasNovo(toasN,aoasN,aodsN,b,m,s,s_angs,s_toas);
    
        [es(j)]=Shikur.Shikur(toasN, aodsN, aoasN, b, m, s, s_toas, s_angs);
        [~,~,ew(j),ewe(j)]=AlgoritmoHenk.main(b,s,m,s_toas,s_angs);
        et(j)=tobias.tobias(toasN,aoasN,aodsN,b,s,m);

    end
    etL(i)=sqrt(mean(eL.^2));
    et1(i)=sqrt(mean(e1.^2));
    et2(i)=sqrt(mean(e2.^2));
    ets(i)=sqrt(mean(es.^2));
    etw(i)=sqrt(mean(ew.^2));
    ett(i)=sqrt(mean(et.^2));

    etee1(i)=sqrt(mean(ee1.^2));
    etee2(i)=sqrt(mean(ee2.^2));
    etewe(i)=sqrt(mean(ewe.^2));
end

figure
hold on
semilogy(sn,etL,'bX--');
semilogy(sn,et1,'bX-');
semilogy(sn,et2,'d-','Color','#77AC30');
semilogy(sn,ets,'r^-');
semilogy(sn,etw,'ms-');
semilogy(sn,ett,'ko-');
legend({'Método LoS','Proposta 1','Proposta 2','(Shikur; Weber, 2014)','(Wymeersch, 2018)','(Wei; Palleit; Weber,2011)'},'NumColumns',2);
xlabel('SNR (dB)')
ylabel('Localização do UE - RMSE (m)');
grid on;

figure
hold on
semilogy(sn,etee1,'bX-');
semilogy(sn,etee2,'d-','Color','#77AC30');
semilogy(sn,etewe,'ms-');
legend('Proposta 1','Proposta 2','(Wymeersch, 2018)');
xlabel('SNR (dB)')
ylabel('Localização dos SCs - RMSE (m)');
grid on;