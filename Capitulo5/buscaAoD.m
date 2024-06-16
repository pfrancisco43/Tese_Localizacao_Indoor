%Ângulo que sai de uma EB para várias ESs
%Autor: Paulo Franicsco
function aods=buscaAoD (Eb,Es)
    %No cenário los para todos os espalhadores
    ns=size(Es,2);
    aods=zeros(ns,2);
    for i=1:ns
        aod_az=azimute(Eb,Es(:,i));
        aod_el=(pi/2) - atan((Es(3,i)-Eb(3))/sqrt((Es(1,i)-Eb(1))^2 + (Es(2,i)-Eb(2))^2));%αl
        aods(i,:)=[aod_az aod_el];
    end
end

function a=azimute(EB, ES)
    if ES(1) < EB(1)
        a=atan((ES(2)-EB(2))/(ES(1)-EB(1)))+pi;
    else
        a=atan((ES(2)-EB(2))/(ES(1)-EB(1)));
    end
end

