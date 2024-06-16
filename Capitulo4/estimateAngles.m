function [aodEst,aoaEst]=estimateAngles(h,rangeAz_aod,rangeEl_aod,rangeAz_aoa,rangeEl_aoa)
    %% Estima os ângulos considerando a posição da correlação máxima
    %Autor: Paulo Francisco
    NrA=numel(rangeAz_aod);
    NrE=numel(rangeEl_aod);
    id_1=ceil(h/(NrA*NrE));
    if mod(id_1,NrE)== 0
        id_A=id_1/NrE;
        id_E=mod(id_1,NrE)+NrE;
    else
        id_A=fix(id_1/NrE+1);
        id_E=mod(id_1,NrE);
    end
    aod_az=rangeAz_aod(id_A);
    aod_el=rangeEl_aod(id_E);
    aodEst=[aod_az aod_el];

    id_1=h-(id_1-1)*(NrA*NrE);
    if mod(id_1,NrE)== 0
        id_A=id_1/NrE;
        id_E=mod(id_1,NrE)+NrE;
    else
        id_A=fix(id_1/NrE+1);
        id_E=mod(id_1,NrE);
    end
    aoa_az=rangeAz_aoa(id_A);
    aoa_el=rangeEl_aoa(id_E);
    aoaEst=[aoa_az aoa_el];
end