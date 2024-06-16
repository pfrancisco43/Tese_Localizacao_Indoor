%Ordena os parâmetros usando o ToA como chave de busca
%Autor: Paulo Francisco
function [ToAs_E, AoDs_E, AoAs_E]=orderParameters(toas, aods, aoas)
    [ToAs_E,i] = sort(toas);
    AoDs_E = aods(i,:);
    AoAs_E = aoas(i,:);
end