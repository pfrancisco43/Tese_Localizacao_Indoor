%Compara IMTL com série de Taylor

clc;clear;close all;
EM=[70; 0; 0];       %Posição do móvel 3d
sToa = 0.1;          %desvio padrão para ToA
nsim=1000;

q=[5,6,7,8,9,10];

for i=1:length(q)
    ns=q(i); %Qtd de BS    
    nd=3; %dimensões    
    for qs=1:nsim  
        EB=geraBSs(ns,nd); %Gera posições das BSs
        toas=buscaToAs(EB,EM);
        ruidoToA=sqrt(sToa^2)*randn(ns,1);
        toas=toas+ruidoToA;
        [pos,crlb(qs),qit(qs)]=estimadorTaylorToA(EM,EB,toas,sToa);
        et(qs)=norm(EM-pos);
         
        [pos,qii(qs)]=estimadorIMTL(EB,toas);
        ei(qs)=norm(EM-pos);
    end
    rmset(i)=sqrt(mean(et.^2));
    crlbs(i)=sqrt(mean(crlb.^2));
    mqit(i)=mean(qit);
    mqii(i)=mean(qii);
    rmsei(i)=sqrt(mean(ei.^2));
end
semilogy(q,rmset,'bO-','LineWidth',1);
hold on
semilogy(q,crlbs,'r--','LineWidth',1);
semilogy(q,rmsei,'md-','LineWidth',1);

ylabel('RMSE (m)')
xlabel('Quantidade de BSs')
grid on
legend('Método de Taylor','CRLB','IMTL');

