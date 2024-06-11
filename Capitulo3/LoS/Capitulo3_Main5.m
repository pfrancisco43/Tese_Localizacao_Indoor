%Técnicas que usam RSS
clc;clear;close all;
sigma.ToA = 0.1;  %[m] desvio padrão para ToA
sigma.RSS = 1;    %[dB] desvio padrão para RSS
p0=30;
beta=3;

q=[5,6,7,8,9,10];
eqmtoas=[];
eqmrsss=[];
eqmdiss=[];
eqmwcls=[];
qits=[];
qirs=[];
qids=[];
EM=[10; 10; 0]; %Posição do móvel 3d
for i=1:length(q)
    ns=q(i); %Qtd de BSs
    sig1=sigma.ToA;
    sig2=sigma.RSS;
    nd=3; %dimensões    
    for qs=1:1000
        ES=geraBSs(ns,nd); %Posições das BSs
        toas=buscaToAs(ES,EM);
        ruidoToA=sqrt(sig1^2)*randn(ns,1);
        toas=toas+ruidoToA;
        rss=buscaRSS(ES,EM,p0,beta);
        ruidoRSS=sqrt(sig2^2)*randn(ns,1);
        rss=rss+ruidoRSS;
        dis=rss2dis(rss)';
        [pos,crlb(qs),qit(qs)]=estimadorTaylorToA(EM,ES,toas,sig1);
        eqmtoa(qs)=norm(EM-pos);
        
        [pos,~,qir(qs)]=estimadorTaylorRSS(EM,ES,rss,sig2,p0,beta);
        eqmrss(qs)=norm(EM-pos);
        
        [pos,qid(qs)]=estimadorTaylorDis(ES,dis,EM,sig1);
        eqmdis(qs)=norm(EM-pos);
        
        [pos]=estimadorWCL(ES,rss');
        eqmwcl(qs)=norm(EM-pos');
    end
    eqmtoas=[eqmtoas sqrt(mean(eqmtoa.^2))];
    eqmrsss=[eqmrsss sqrt(mean(eqmrss.^2))];
    eqmdiss=[eqmdiss sqrt(mean(eqmdis.^2))];
    eqmwcls=[eqmwcls sqrt(mean(eqmwcl.^2))];
    qits=[qits mean(qit)];
    qirs=[qirs mean(qir)];
    qids=[qids mean(qid)];
end
semilogy(q,eqmtoas,'bO-','LineWidth',1);
hold on
semilogy(q,eqmrsss,'rs-','LineWidth',1);
semilogy(q,eqmdiss,'md-','LineWidth',1);
semilogy(q,eqmwcls,'k*-','LineWidth',1);

ylabel('RMSE (m)')
xlabel('Quantidade de BSs')
grid on
legend('ToA','RSS-PL','Distância-PL','WCL');

