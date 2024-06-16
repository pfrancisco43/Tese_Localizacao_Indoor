%Determina a probabilidade em conformidade com a TR38901 da 3GPP
%Autor: Paulo Francisco
function PrLos=DeterminaProbLos(b,m)
    hut=1;
    hbs=b(3);
    d2=norm(b(1:2)-m(1:2));  %####
    %Para cenário inF-SL
    dclutter=10; %em metros
    r=30/100; % em %
    hc=4; %Altura dos obstáculos
    Ksub=-(dclutter/(log(1-r)))*((hbs-hut)/(hc-hut));
    PrLos=exp(-(d2/Ksub));
end
