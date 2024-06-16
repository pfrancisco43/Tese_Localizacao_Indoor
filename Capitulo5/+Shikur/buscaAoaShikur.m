function aoas=buscaAoaShikur(ES,EM)
    qs=size(ES,2);
    M=EM;
    for i=1:qs
        S=ES(:,i);
        az=(pi/2) * (1 - sign(S(1) - M(1))) + atan((S(2) - M(2)) / (S(1) - M(1)));
        el=(pi/2) - atan((S(3) - M(3)) / sqrt((S(1) - M(1))^2 + (S(2) - M(2))^2));
        
        aoas(i,:)=[az el];
    end