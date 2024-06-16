function tdoas=buscaTDoAShikur (EB,ES,EM)    
    %No cenário los para todos os espalhadores
    db=0;
    ns=size(ES,2);
    for i=1:ns
        d1=norm(EB-ES(:,i));   
        d2=norm(EM-ES(:,i)); 
        toa=d1+d2;  
        tdoas(i,1)=toa-db;
    end    
end
    
        