%é o ângulo de chegada na EM vindo do Espalhador ou EB
%Autor: Paulo Franicsco
function aods=buscaAoA(Es,Em)   
    ns=size(Es,2);
    aods=zeros(ns,2);
    for i=1:ns
        aod_az=azimute(Es(:,i), Em);
        aod_el=(pi/2) - atan((Es(3,i)-Em(3))/sqrt((Es(1,i)-Em(1))^2 + (Es(2,i)-Em(2))^2));%αl
        aods(i,:)=[aod_az aod_el];
    end
end

function a=azimute(ES, EM)
    if ES(1) < EM(1)
        a=atan((ES(2)-EM(2))/(ES(1)-EM(1)))+pi;
    else
        a=atan((ES(2)-EM(2))/(ES(1)-EM(1)));       
    end    
end

