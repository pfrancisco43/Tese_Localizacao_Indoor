%Determina o ToA
%Autor: Paulo Franicsco
function toas=buscaToAsNlos (EB,ES,EM) 
    ns=size(ES,2);
    toas=zeros(ns,1);
    for i=1:ns
        d1=norm(EB-ES(:,i));   
        d2=norm(EM-ES(:,i)); 
        toas(i,1)=d1+d2;      
    end    
end
    
        