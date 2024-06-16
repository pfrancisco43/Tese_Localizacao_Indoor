%Determina a posição da MS com base no ToA e AoD
%Autor: Paulo Francisco
function [eA,eH,eV,pos]=estimadorToAAoD(ES,toas,aods,EM)
    ns=size(ES,2);
    x=zeros(3,ns);  
    
    for i=1:ns
        x(1,i)=toas(i)*sin(aods(i,2))*cos(aods(i,1))+ES(1,i);
        x(2,i)=toas(i)*sin(aods(i,2))*sin(aods(i,1))+ES(2,i);
        x(3,i)=toas(i)*cos(aods(i,2))+ES(3,i);       
    end

    pos=mean(x,2);
    eA=norm(pos-EM);
    eH=norm(pos(1:2)-EM(1:2));
    eV=norm(pos(3)-EM(3));
    
    toas=1./toas';
    pos(1)=(sum(x(1,:).*toas))/sum(toas);
    pos(2)=(sum(x(2,:).*toas))/sum(toas);
    pos(3)=(sum(x(3,:).*toas))/sum(toas);   
end