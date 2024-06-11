%código que replica os resultados encontrados em:
%A Simple and Efficient Iterative Method for Toa Localization
%Para localização em ambiente LOS com sensores conhecidos
function [pos, qi] = estimadorIMTL(ES,toas)
    x=[1;1;1]; %chute inicial    
    K=100; %Limite de iterações
    eps=10e-4; %erro
    M=size(ES,2);
    e=zeros(M,K);
    for k=1:K  %iterações        
        S=0;
        for i=1:M %até a qnt de espalhadores
            a=norm(x-ES(:,i));
            b=toas(i)/a;
            c=(x-ES(:,i));
            d=b*c+ES(:,i);
            S=S+d;
            e(i,k)=(toas(i)-norm(x-ES(:,i)))/toas(i);
        end
        x=(1/M)*S; %novo x
        if k==1
            ea=0; %erro anterior
        else
            ea=e(:,k-1); %erro anterior
        end
        if norm(ea - e(:,k)) <= eps
            break;  
        end
    end
    qi=k;
    pos=x;
end

