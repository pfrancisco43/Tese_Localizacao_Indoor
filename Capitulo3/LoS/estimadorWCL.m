%Centroid ponderado pelo rss
function pos=estimadorWCL(BS,p)
    %p deve ser um vetor linha
    nd=size(BS,1);
    for i=1:nd
        a=p.*BS(i,:);
        pos(i)=sum(a)/sum(p);
    end    
end