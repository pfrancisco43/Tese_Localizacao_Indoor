%replicação do artigo AOD/AOA/TOA-based 3D Positioning in NLOS
%Multipath Environments
%Autores: WEI, X.; PALLEIT, N.; WEBER, T.

function [eq, pos]=tobias(toas,aoas,aods,EB,ES,EM)    
    k=0;
    ns=size(ES,2);
    for i=1:ns
        k=k+1;
        M(k,1)=(cos(aods(i,2)) + cos(aoas(i,2)));
        M(k,2)=0;
        M(k,3)=-(sin(aods(i,2)) * cos(aods(i,1)) + sin(aoas(i,2)) * cos(aoas(i,1)));
        b(k,1)=toas(i) * (sin(aods(i,2)) * cos(aods(i,1)) * cos(aoas(i,2)) - sin(aoas(i,2)) * cos(aoas(i,1)) * cos(aods(i,2)));
        
        k=k+1;
        M(k,1)=0;
        M(k,2)=(cos(aods(i,2)) + cos(aoas(i,2)));
        M(k,3)= -(sin(aods(i,2)) * sin(aods(i,1)) + sin(aoas(i,2)) * sin(aoas(i,1))); 
        b(k,1)=toas(i) * (sin(aods(i,2)) * sin(aods(i,1)) * cos(aoas(i,2)) - sin(aoas(i,2)) * sin(aoas(i,1)) * cos(aods(i,2)));
    end
    
    pos=(pinv(M'*M)*M')*b+EB;
    eq=norm(pos-EM);
end