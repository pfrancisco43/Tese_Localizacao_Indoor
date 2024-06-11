clc;close all; clear;
EB=[0;0;10];
EM=[20;0;0];

sigmaT=1;
sigmaA=0.05;
L=5; %Quantidade de caminhos

nsim=1000;
vetor=linspace(0.01,0.08,10);
tam=size(vetor,2);
for i=1:tam
    for si=1:nsim
        %%%VARIANDO O TOA%%%
        sigmaA=vetor(i);
        
        ES=geraEspalhadores(L,3);
        
        toas=buscaToAsNlos (EB,ES,EM);
        ruidoToA=sqrt(sigmaT^2)*randn(L,1);
        toas=toas+ruidoToA;
        
        aods=buscaAoDs(EB,ES);
        ruidoAoD=sqrt(sigmaA^2)*randn(L,2);
        aods=aods+ruidoAoD;
        
        aoas=buscaAoAs(ES,EM);
        ruidoAoA=sqrt(sigmaA^2)*randn(L,2);
        aoas=aoas+ruidoAoA;
        
        et(si)=tobias.tobias(toas,aoas,aods,EB,ES,EM);        
        [~,~,eh(si)]=AlgoritmoHenk.main(EB,ES,EM,sigmaT,sigmaA);        
        [es(si), crlb(si)]=Shikur.Shikur(toas, aods, aoas, EB, EM, ES, sigmaT, sigmaA);                
    end    
    met(i)=mean(et);   
    meh(i)=mean(eh);
    mes(i)=mean(es);
    mcrlb(i)=mean(crlb);
   
end

semilogy(vetor,met,'^-r')
hold on;
semilogy(vetor,meh,'s-m')
semilogy(vetor,mes,'o-k')
semilogy(vetor,mcrlb,'--','Color','#EDB120','LineWidth',2)

grid on;
xlabel('\sigma_{AoA/AoD} (rad)')
ylabel('RMSE (m)')
legend('(Wei; Palleit; Weber,2011)','(Wymeersch, 2018)','(Shikur; Weber, 2014)','CRLB')


