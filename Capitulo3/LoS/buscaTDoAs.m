function tdoas=buscaTDoAs (BS,MS)    
    %No cenário los para todas as BSs
    EB=BS(:,1);
    BS=BS(:,2:end);
    db=norm(MS-EB);
    ns=size(BS,2);
    for i=1:ns
        ds=norm(MS-BS(:,i));  
        tdoas(i,1)=ds-db;
    end    
end
    
        