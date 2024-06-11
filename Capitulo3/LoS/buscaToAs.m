function toas=buscaToAs (BS,MS)    
    %No cenário los para todoas as BSs
    ns=size(BS,2);
    for i=1:ns
        toas(i,1)=norm(MS-BS(:,i));       
    end    
end
    
        