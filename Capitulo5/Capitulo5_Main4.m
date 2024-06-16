%Gera gráfico ilustrativo (alvo) da acurácia, mostrando:
%Exatidão e precisão
%Autor: Paulo Francisco
clc;clear;close all
load dataSet_CaminhoLoS
%load dataSet_CaminhoLoS_Adaptativo
load dataSet_4Caminhos_Nlos
% load dataSet_4Caminhos_NLoSBeamBeta085 %base com beam beta 085, 200 amostras

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

snr=9; %Escolhe uma das SNRs

a=reshape(ToAs_Bf(:,:,snr),[],1);
s_toas=sqrt(var(a)); %Variância na estimação de ToA
a=reshape(AoDs_Bf(:,:,:,snr),[],1);
s_angs=sqrt(var(a)); %Variância na estimação dos angulos
for j=1:Nruns % QUANTIDADE DE RUNS
    toaL=ToAs_BfLoS(1,j,snr)*c;
    aodL=ADos_BfLoS(1,:,j,snr);
     
    toasN=ToAs_Bf(:,j,snr)*c;
    aodsN=AoDs_Bf(:,:,j,snr);
    aoasN=AoAs_Bf(:,:,j,snr);
    aoasN=aod_aoa_Swap(aoasN);

    [~,~,~,eL(:,j)]=estimadorToAAoD(b,toaL,aodL,m);

    [e1(j),eh,ev,ss,rr,kk,posi(:,j),posf(:,j),esi,esf,e2(j),cr,ee1(j),ee2(j)]=estimadorIntersecLinhasNovo(toasN,aoasN,aodsN,b,m,s,s_angs,s_toas);

    [~,~,es(:,j)]=Shikur.Shikur(toasN, aodsN, aoasN, b, m, s, s_toas, s_angs);
    [~,~,~,~,ew(:,j)]=AlgoritmoHenk.main(b,s,m,s_toas,s_angs);
    [~,et(:,j)]=tobias.tobias(toasN,aoasN,aodsN,b,s,m);
end

p = nsidedpoly(1000, 'Center', [m(1) m(2)], 'Radius', 16);
plot(p, 'FaceColor', '#A9A9A9', 'EdgeColor', '#A9A9A9')
axis equal
hold on
p = nsidedpoly(1000, 'Center', [m(1) m(2)], 'Radius', 13);
plot(p, 'FaceColor', 'w', 'EdgeColor', 'w')
p = nsidedpoly(1000, 'Center', [m(1) m(2)], 'Radius', 10);
plot(p, 'FaceColor', '#A9A9A9', 'EdgeColor', '#A9A9A9')
p = nsidedpoly(1000, 'Center', [m(1) m(2)], 'Radius', 7);
plot(p, 'FaceColor', 'w', 'EdgeColor', 'w')
p = nsidedpoly(1000, 'Center', [m(1) m(2)], 'Radius', 4);
plot(p, 'FaceColor', '#A9A9A9', 'EdgeColor', '#A9A9A9')
p = nsidedpoly(1000, 'Center', [m(1) m(2)], 'Radius', 2);
plot(p, 'FaceColor', 'w', 'EdgeColor', 'w')
t=6;
p1=plot3(eL(1,:),eL(2,:),eL(3,:),'Xb','markersize',t);
p2=plot3(posi(1,:),posi(2,:),posi(3,:),'>','markersize',t,'Color','#8B4513');
p3=plot3(posf(1,:),posf(2,:),posf(3,:),'d','markersize',t,'Color','#77AC30');

p4=plot3(es(1,:),es(2,:),es(3,:),'^r','markersize',t);
p5=plot3(ew(1,:),ew(2,:),ew(3,:),'ms','markersize',t);
p6=plot3(et(1,:),et(2,:),et(3,:),'ko','markersize',t);
%plot3 (m(1),m(2),m(3),'.r','markersize',30);
xlim([-10 23])
ylim([-6 26])
xlabel('x')
ylabel('y')
zlabel('z')


s=determinaPrecisao(eL);
s1=sprintf('Método LoS                        [%.2f]',s);
s1=strrep(s1,'.',',');
s=determinaPrecisao(posi);
s2=sprintf('Proposta 1                          [%.2f]',s);
s2=strrep(s2,'.',',');
s=determinaPrecisao(posf);
s3=sprintf('Proposta 2                          [%.2f]',s);
s3=strrep(s3,'.',',');

s=determinaPrecisao(es);
s4=sprintf('(Shikur; Weber, 2014)        [%.2f]',s);
s4=strrep(s4,'.',',');
s=determinaPrecisao(ew);
s5=sprintf('(Wymeersch, 2018)            [%.2f]',s);
s5=strrep(s5,'.',',');
s=determinaPrecisao(et);
s6=sprintf('(Wei; Palleit; Weber,2011) [%.2f]',s);
s6=strrep(s6,'.',',');

leg=legend([p1,p2,p3,p4,p5,p6],s1,s2,s3,s4,s5,s6);
title(leg,'Método      [Precisão]')

function p=determinaPrecisao(a)
    % Calcula a variância ao longo de cada componente
    std_devs = std(a, 0, 2);  % 0 para o fator de normalização N-1, 2 para operar ao longo das colunas
    % Calcula a média das variâncias das três componentes
    p = mean(std_devs);
end