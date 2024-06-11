function aods=buscaAodShikur(ES,EB)
    qs=size(ES,2);
    B=EB;
    for i=1:qs
        S=ES(:,i);
        az=(pi/2) * (1 - sign(S(1) - B(1))) + atan((S(2) - B(2)) / (S(1) - B(1)));
        el=(pi/2) - atan((S(3)-B(3)) / sqrt((S(1) - B(1))^2 + (S(2) - B(2))^2));
        
        aods(i,:)=[az el];
    end