%CDF - Compara ToA com TDoA
clc;clear;close all;
sigma.ToA = 0.1;  %[m] desvio padrão para ToA
sigma.TDoA = 0.1; %[m] desvio padrão para TDoA

q=[7];
eqmtaylortoas=[];
eqmtaylortdoas=[];
eqmtayloralls=[];
crlbs=[];
qits=[];
qitds=[];
qialls=[];
for i=1:length(q)
    ns=q(i); %Qtd de BS
    sig=sigma.ToA;
    nd=3; %dimensões   
    EM=[70; 0; 0]; %Posição do móvel 3d
    
    nr=1000;
    for qs=1:nr
        ES=geraBSs(ns,nd); %Posições das BSs
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
end
yy=linspace(0,1,nr);
semilogx(sort(eqmtaylortoa),yy,'b-','LineWidth',1);
hold on
semilogx(sort(eqmtaylortdoa),yy,'m-','LineWidth',1);
semilogx(sort(eqmtaylorall),yy,'k-','LineWidth',1);

ylabel('CDF')
xlabel('Erro (m)')
grid on
legend('ToA','TDoA','ToA+TDoA');

