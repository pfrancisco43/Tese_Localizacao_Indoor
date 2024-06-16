%Chama o método DCS-SOMP Simples (Não considera iterações, it=0)
%Autor: Paulo Francisco
function [ToAs_E, AoDs_E, AoAs_E,its]=DCS_SOMP_Simples(Nt, Nr, Ns, N, B, c, y, x, L,fc,arranjo)
    L_az=10;
    L_el=L_az;
    it=0;
    lim=1e-6;
    [ToAs_E, AoDs_E, AoAs_E,its]=parameterEstimation (Nt, Nr, Ns, N, B, c, y, x, L, L_az, L_el, it,lim,fc,arranjo);

    ToAs_E=ToAs_E(end,:)';
    
    rAoD=zeros(L, 2);
    rAoA=zeros(L, 2);
    for k = 1:L
            rAoD(k, :) = AoDs_E(end, :, k);
            rAoA(k, :) = AoAs_E(end, :, k);
    end
    AoDs_E=rAoD;
    AoAs_E=rAoA;

    %Order the parameters
    [ToAs_E, AoDs_E, AoAs_E]=orderParameters(ToAs_E, AoDs_E, AoAs_E);
end