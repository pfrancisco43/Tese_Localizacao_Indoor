%No cenário los para todas as BSs
function aoas=buscaAoAs (BS,MS)
    ns=size(BS,2);
    for i=1:ns
        aoas(i,:)=[atan((MS(2)-BS(2,i))/(MS(1)-BS(1,i))) acos((MS(3)-BS(3,i))/norm(BS(:,i)-MS))];
    end    
end
    
        