%Compara IMTL com série de Taylor

clc;clear;close all;
sigma.AoA_az = 0.01;    %[radian]  standard deviation of DOA in azimuth
sigma.AoA_el = 0.01;    %[radian]  standard deviation of DOA in elevation
sigma.AoD_az = 0.01;    %[radian]  standard deviation of DOD in azimuth
sigma.AoD_el = 0.01;    %[radian]  standard deviation of DOD in elevation
sigma.ToA = 0.1;        %[m]       standard deviation of TOA in elevation
sigma.TDoA = 0.1;        %[m]       standard deviation of TOA in elevation

q=[0.1 0.5 1 1.5 2];
eqmtaylors=[];
crlbs=[];
qits=[];
qiis=[];
eqmimtls=[];
for i=1:length(q)
    ns=6;  %Qtd de BS
    sig=q(i);
    nd=3; %dimensões    
    EM=[70; 0; 0]; %Posição do móvel 3d
    for qs=1:1000
        ES=geraBSs(ns,nd); %Gera posições das BSs
        toas=buscaToAs(ES,EM);
        ruidoToA=sqrt(sig^2)*randn(ns,1);
        toas=toas+ruidoToA;
        [pos,crlb(qs),qit(qs)]=estimadorTaylorToA(EM,ES,toas,sig);
        eqmtaylor(qs)=norm(EM-pos);
        
        [pos,qii(qs)]=estimadorIMTL(ES,toas);
        eqmimtl(qs)=norm(EM-pos);
    end
    eqmtaylors=[eqmtaylors sqrt(mean(eqmtaylor.^2))];
    crlbs=[crlbs sqrt(mean(crlb.^2))];
    qits=[qits mean(qit)];
    qiis=[qiis mean(qii)];
    eqmimtls=[eqmimtls sqrt(mean(eqmimtl.^2))];
end
semilogy(q,eqmtaylors,'bO-','LineWidth',1);
hold on
semilogy(q,crlbs,'r--','LineWidth',1);
semilogy(q,eqmimtls,'md-','LineWidth',1);


ylabel('RMSE (m)')
xlabel('\sigma_{ToA}')
grid on
legend('Método de Taylor','CRLB','IMTL');

