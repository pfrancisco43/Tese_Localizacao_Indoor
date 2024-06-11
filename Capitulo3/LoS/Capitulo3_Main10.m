%Compara técnicas híbridas: distância e ângulos

clc;clear;close all;
sigma.AoA_az = 0.05;    %[rad]  desvio padrão para AoA_az
sigma.AoA_el = 0.01;    %[rad]  desvio padrão para AoA_el
sigma.AoD_az = 0.01;    %[rad]  desvio padrão para AoD_az
sigma.AoD_el = 0.01;    %[rad]  desvio padrão para AoD_el
sigma.ToA = 0.1;        %[m]    desvio padrão para ToA
sigma.TDoA = 0.1;       %[m]    desvio padrão para TDoA

q=[1];
eqmtaylors=[];
crlbs=[];
qits=[];
qiis=[];
nsim=1000;
for i=1:length(q)
    tic
    ns=7;  %Qtd de BSs
    sig1=sigma.ToA^2;
    sig2=sigma.AoA_az^2;
    nd=3; %dimensões    
    EM=[70; 0; 0]; %Posição do móvel 3d
    for qs=1:nsim
        ES=geraBSs(ns,nd); %Posições das BSs
        toas=buscaToAs(ES,EM);
        ruidoToA=sqrt(sig1)*randn(ns,1);
        toas=toas+ruidoToA;
        
        tdoas=buscaTDoAs(ES,EM);
        ruidoTDoA=sqrt(sig1)*randn(ns-1,1);
        tdoas=tdoas+ruidoTDoA;
        
        aoas=buscaAoAs(ES,EM);
        ruidoAoA=sqrt(sig2)*randn(ns,2);
        aoas=aoas+ruidoAoA;
        
        [pos,~,qi]=estimadorTaylorToA(EM,ES,toas,sig1);
        eqmtoa(qs)=norm(EM-pos);
        
        [pos,~,qi]=estimadorTaylorToAAoA(EM,ES,toas, aoas,sig1,sig2);
        eqmtaoa(qs)=norm(EM-pos);
        
        [pos,~,qi]=estimadorTaylorTDoA(EM,ES,tdoas,sig1);
        eqmtdoa(qs)=norm(EM-pos);
        
        [pos,~,qi]=estimadorTaylorTDoAAoA(EM,ES,tdoas, aoas,sig1,sig2);
        eqmtdaoa(qs)=norm(EM-pos);
        
    end
    toc
    eqmtoas=sort(eqmtoa);
    eqmtaoas=sort(eqmtaoa);
    eqmtdoas=sort(eqmtdoa);
    eqmtdaoas=sort(eqmtdaoa);
end
yy=linspace(0,1,nsim);
semilogx(eqmtoas,yy,'b-','LineWidth',1);
hold on
semilogx(eqmtaoas,yy,'m--','LineWidth',1);
semilogx(eqmtdoas,yy,'k-','LineWidth',1);
semilogx(eqmtdaoas,yy,'r--','LineWidth',1);


ylabel('CDF')
xlabel('Erro (m)')
grid on
legend('ToA','ToA + AoA','TDoA','TDoA + AoA');

