function ToAE=estimateToA(h,Ts,c,N)
    %% Estima o ToA com base no período e quantidade de subportadoras
    %Autor: Paulo Francisco
    d=-mean(diff(phase(h)))*(N*Ts)*c/(2*pi);
    if (d<0)
        d=d+N*Ts*c;
    end
    ToAE=d/c;
end