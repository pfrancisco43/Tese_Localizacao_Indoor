%No cenário los para todos as BSs
function rsss=buscaRSS (BS,x,p0,beta) 
    ns=size(BS,2);
    for i=1:ns
        rsss(i,1)=p0-10*beta*log10(norm(x-BS(:,i)));
    end    
end
    
        