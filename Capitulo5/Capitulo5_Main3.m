%Gera gráficos comparativos para as propostas em comparação com
%métodos da literatura - CDF e BoxPlot
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
s=S;
c=300;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nruns=100;
for i=1:size(sn,2)
    a=reshape(ToAs_Bf(:,:,i),[],1);
    s_toas=sqrt(var(a)); %Variância na estimação de ToA
    a=reshape(AoDs_Bf(:,:,:,i),[],1);
    s_angs=sqrt(var(a)); %Variância na estimação dos angulos
    for j=1:Nruns % QUANTIDADE DE RUNS
        toaL=ToAs_BfLoS(1,j,i)*c;
        aodL=ADos_BfLoS(1,:,j,i);

        toasN=ToAs_Bf(:,j,i)*c;
        aodsN=AoDs_Bf(:,:,j,i);
        aoasN=AoAs_Bf(:,:,j,i);
        aoasN=aod_aoa_Swap(aoasN);

        [eL(i,j),~]=estimadorToAAoD(b,toaL,aodL,m);

        [e1(i,j),eh,ev,ss,rr,kk,posi,posf,esi,esf,e2(i,j),cr,ee1(j),ee2(j)]=estimadorIntersecLinhasNovo(toasN,aoasN,aodsN,b,m,s,s_angs,s_toas);

        [es(i,j)]=Shikur.Shikur(toasN, aodsN, aoasN, b, m, s, s_toas, s_angs);
        [~,~,ew(i,j),ewe(j)]=AlgoritmoHenk.main(b,s,m,s_toas,s_angs);
        et(i,j)=tobias.tobias(toasN,aoasN,aodsN,b,s,m);
    end    
end
eL=sort(reshape(eL,1,[]));
e1=sort(reshape(e1,1,[]));
e2=sort(reshape(e2,1,[]));
es=sort(reshape(es,1,[]));
ew=sort(reshape(ew,1,[]));
et=sort(reshape(et,1,[]));

y=linspace(0,1,Nruns*9);
figure
hold on
semilogx(eL,y,'-b','LineWidth',2);
semilogx(e1,y,'--b','LineWidth',2);
semilogx(e2,y,'-b','LineWidth',2,'Color','#77AC30');
semilogx(es,y,'-r','LineWidth',2);
semilogx(ew,y,'-m','LineWidth',2);
semilogx(et,y,'-k','LineWidth',2);

semilogx([min(eL),max(et)],[0.95,0.95],'--','Color','#EDB120','LineWidth',3)

legend('Método LoS','Proposta 1','Proposta 2','(Shikur; Weber, 2014)','(Wymeersch, 2018)','(Wei; Palleit; Weber,2011)');
grid on;
xlabel('Erro (m)')
ylabel('CDF')

sqrt(mean(eL).^2)
sqrt(mean(e1).^2)
sqrt(mean(e2).^2)
sqrt(mean(es).^2)
sqrt(mean(ew).^2)
sqrt(mean(et).^2)

%95%
eL(855)
e1(855)
e2(855)
es(855)
ew(855)
et(855)


figure
boxplot([eL',e1' e2' es' ew' et'],{'Método LoS','Proposta 1','Proposta 2','(Shikur; Weber, 2014)','(Wymeersch, 2018)','(Wei; Palleit; Weber,2011)'})
ylabel('Erro (m)')
grid on