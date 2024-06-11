%Transforma RSS em distância
function d=rss2dis(rss)
    p0=30;
    beta=3;
    for i =1:length(rss)
        d(i)=10^((p0-rss(i))/(10*beta));       
    end
end