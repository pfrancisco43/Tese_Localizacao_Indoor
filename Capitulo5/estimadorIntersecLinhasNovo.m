%Esta função prepara os pontos para a estimação do ponto de intersecção
%Autor: Paulo Francisco
function [e1,eh,ev,s,r,k,posi,posf,esi,esf,e2,cr,EstError1,EstError2]=estimadorIntersecLinhasNovo(toasN,aoasN,aodsN,EB,EM,ES,s_angs,s_toas)
    
    [s,r]=determinaPts(aodsN,aoasN,toasN,EB);   
    
    posi=lineIntersect3D_2(s,r);
    s=s';
    r=r';
    posi=posi';
    [esi,k]=determinaEspalhadores (EB,posi,aoasN,aodsN,toasN);
    [e2,cr,posf,esf]=estimadorAjusteFinalMelhorado(EB,EM,ES,toasN, aodsN, aoasN,s_toas,s_angs,posi,esi);
    
    eh=norm(EM(1:2)-posf(1:2));
    ev=norm(EM(3)-posf(3));

    e1=norm(posi-EM);%Somente intesecção de Linhas

    EstError1=determinaErro(esi,ES);
    EstError2=determinaErro(esf,ES);
end

function eq=determinaErro(y,ES)
    for i=1:size(ES,2)
        e(i)=norm(y(:,i)-ES(:,i));
    end
    eq=mean(e);
end

function [y r]=determinaEspalhadores(EB,EM,aoas,aods,toas)
    ns=length(toas);
    pi=[EB';EM'];
    for i=1:ns
        s(:,i)=determinaPontoS(aods(i,1),aods(i,2),toas(i),EB);
        r(:,i)=determinaPonto2(aoas(i,1),aoas(i,2),toas(i),EM);
        pf=[s(:,i)';r(:,i)'];
        [y(:,i), dis]=lineIntersect3D(pi,pf);
    end   
end

%função que determina os pontos que passam o espalhador e EM
function [s,r]=determinaPts(aods,aoas,toas,EB)
    for i=1:size(aods,1)
        s(i,:)=determinaPontoS(aods(i,1),aods(i,2),toas(i),EB);
        r(i,:)=determinaPontoR(aoas(i,1),aoas(i,2),toas(i),EB); %Caso NLos
    end
end

function p=determinaPontoS(aa,ae,toa,EB)
    m=[sin(ae)*cos(aa);
        sin(ae)*sin(aa);
        cos(ae)];
    p=toa*m+EB;
end

function p=determinaPontoR(aa,ae,toa,EB)    
    m=[sin(pi-ae)*cos(aa-pi);
        sin(pi-ae)*sin(aa-pi);
        cos(pi-ae)];
    p=toa*m+EB;
end

function p=determinaPonto2(aa,ae,toa,EM)
    m=[sin(pi-ae)*cos(aa-pi);
        sin(pi-ae)*sin(aa-pi);
        cos(pi-ae)];
    p=EM-toa*m;
end