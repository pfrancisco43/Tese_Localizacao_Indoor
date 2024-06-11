clc;close all; clear;
EB=[0;0;10];
EM=[20;0;0];

sigmaT=2;
sigmaA=0.05;
L=4; %Quantidade de caminhos

nsim=1000;
vetor=[1];
tam=size(vetor,2);
%tic
for i=1:tam
    for si=1:nsim
                       
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
    met(i)=mean(et)   
    meh(i)=mean(eh)
    mes(i)=mean(es)   
end
%toc
et=sort(et);
eh=sort(eh);
es=sort(es);

y=linspace(0,1,nsim);
semilogx(et,y,'-r','LineWidth',2)
hold on;
semilogx(eh,y,'-m','LineWidth',2)
semilogx(es,y,'-k','LineWidth',2)

semilogx([min(es),max(es)],[0.95,0.95],'--','Color','#EDB120','LineWidth',2)

grid on;
xlabel('Erro (m)')
ylabel('CDF')
legend('(Wei; Palleit; Weber,2011)','(Wymeersch, 2018)','(Shikur; Weber, 2014)')
