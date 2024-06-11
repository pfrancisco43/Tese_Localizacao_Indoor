%Compara ToA com TDoA

clc;clear;close all;
sigma.ToA = 0.1;  %[m] desvio padrão para ToA
sigma.TDoA = 0.1; %[m] desvio padrão para TDoA

q=[5,6,7,8,9,10];
eqmtaylortoas=[];
eqmtaylortdoas=[];
eqmtayloralls=[];
crlbs=[];
qits=[];
qitds=[];
qialls=[];
for i=1:length(q)
    ns=q(i);  %Qtd de BSs
    sig=sigma.ToA;
    nd=3; %dimensões   
    EM=[70; 0; 0]; %Posição do móvel 3d    
    for qs=1:1000
        ES=geraBSs(ns,nd); %Gera posições das BSs
        toas=buscaToAs(ES,EM);
        ruidoToA=sqrt(sig^2)*randn(ns,1);
        toas=toas+ruidoToA;
        tdoas=buscaTDoAs(ES,EM);
        ruidoTDoA=sqrt(sig^2)*randn(ns-1,1);
        tdoas=tdoas+ruidoTDoA;
       
        [pos,crlb(qs),qit(qs)]=estimadorTaylorToA(EM,ES,toas,sig);
        eqmtaylortoa(qs)=norm(EM-pos);
      
        
        [pos,~,qitd(qs)]=estimadorTaylorTDoA(EM,ES,tdoas,sig);
        eqmtaylortdoa(qs)=norm(EM-pos);
        
        [pos,~,qiall(qs)]=estimadorTaylorToATDoA(EM,ES,toas,tdoas,sig);
        eqmtaylorall(qs)=norm(EM-pos);
    end    
    eqmtaylortoas=[eqmtaylortoas sqrt(mean(eqmtaylortoa.^2))];
    eqmtaylortdoas=[eqmtaylortdoas sqrt(mean(eqmtaylortdoa.^2))];
    eqmtayloralls=[eqmtayloralls sqrt(mean(eqmtaylorall.^2))];
    crlbs=[crlbs sqrt(mean(crlb.^2))];
    qits=[qits mean(qit)];
    qitds=[qitds mean(qitd)];
    qialls=[qialls mean(qiall)];
end
semilogy(q,eqmtaylortoas,'bO-','LineWidth',1);
hold on
semilogy(q,eqmtaylortdoas,'mX-','LineWidth',1);
semilogy(q,eqmtayloralls,'kd-','LineWidth',1);
semilogy(q,crlbs,'r--','LineWidth',1);

ylabel('RMSE (m)')
xlabel('Quantidade de BSs')
grid on
legend('ToA','TDoA','ToA+TDoA','CRLB');

