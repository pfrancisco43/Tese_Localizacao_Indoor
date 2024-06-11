%Compara técnicas híbridas: distância e ângulos

clc;clear;close all;
sigma.AoA_az = 0.05;    %[rad]  desvio padrão para AoA_az
sigma.AoA_el = 0.01;    %[rad]  desvio padrão para AoA_el
sigma.AoD_az = 0.01;    %[rad]  desvio padrão para AoD_az
sigma.AoD_el = 0.01;    %[rad]  desvio padrão para AoD_el
sigma.ToA = 0.1;        %[m]    desvio padrão para ToA
sigma.TDoA = 0.1;       %[m]    desvio padrão para TDoA

q=[0.01 0.05 0.1 0.15 0.2];
eqmtaylors=[];
crlbs=[];
qits=[];
qiis=[];
eqmtoas=[];
eqmtaoas=[];
eqmtdoas=[];
eqmtdaoas=[];
for i=1:length(q)
    ns=4;% Qtd de BSs
    sig1=sigma.ToA;
    sig2=q(i);
    nd=3; %dimensões    
    EM=[70; 0; 0]; %Posição do móvel 3d
    ES=geraBSs(ns,nd); %Posições das BSs
    for qs=1:1000        
        toas=buscaToAs(ES,EM);
        ruidoToA=sqrt(sig1^2)*randn(ns,1);
        toas=toas+ruidoToA;
        
        tdoas=buscaTDoAs(ES,EM);
        ruidoTDoA=sqrt(sig1^2)*randn(ns-1,1);
        tdoas=tdoas+ruidoTDoA;
        
        aoas=buscaAoAs(ES,EM);
        ruidoAoA=sqrt(sig2^2)*randn(ns,2);
        aoas=aoas+ruidoAoA;
        
        [pos]=estimadorTaylorToA(EM,ES,toas,sig1);
        eqmtoa(qs)=norm(EM-pos);
        
        [pos]=estimadorTaylorToAAoA(EM,ES,toas, aoas,sig1,sig2);
        eqmtaoa(qs)=norm(EM-pos);
        
        [pos]=estimadorTaylorTDoA(EM,ES,tdoas,sig1);
        eqmtdoa(qs)=norm(EM-pos);
        
        pos=estimadorTaylorTDoAAoA(EM,ES,tdoas, aoas,sig1,sig2);
        eqmtdaoa(qs)=norm(EM-pos);
        
    end
    eqmtoas=[eqmtoas sqrt(mean(eqmtoa.^2))];
    eqmtaoas=[eqmtaoas sqrt(mean(eqmtaoa.^2))];
    eqmtdoas=[eqmtdoas sqrt(mean(eqmtdoa.^2))];
    eqmtdaoas=[eqmtdaoas sqrt(mean(eqmtdaoa.^2))];
end
semilogy(q,eqmtoas,'bO-','LineWidth',1);
hold on
semilogy(q,eqmtaoas,'md--','LineWidth',1);
semilogy(q,eqmtdoas,'ks-','LineWidth',1);
semilogy(q,eqmtdaoas,'r*--','LineWidth',1);


ylabel('RMSE (m)')
xlabel('\sigma_{AoA} (rad)')
grid on
legend('ToA','ToA + AoA','TDoA','TDoA + AoA');

