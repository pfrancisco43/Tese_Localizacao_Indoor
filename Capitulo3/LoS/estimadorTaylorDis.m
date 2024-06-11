function [pos,qi]=estimadorTaylorDis(SP,dis,EM,sigma)
    %Estima a localização com base na distncia usando Taylor
        
    ns=size(SP,2);
    x=[1;1;1];
    delt = eye(size(SP,2));
    Q = (sigma^2)*delt;
    C=3;%Quantidade de colunas (parametros a estimar)
    L=ns;%quantidade de medições;
    G=zeros(L,C);
    
    Fx=geraFx(x,SP);%Busca medições com base no primeiro chute, não acrescenta ruido
    tol=1e-4;
    for kk=1:100%quantidade de iterações
        qi=kk;
        for i = 1:L
            %             for j = 1:C
            %                 f=Ja(i,j);
            %                 ex=subs(f,[mx my mz sx sy sz],[x',SP(1,i),SP(2,i),SP(3,i)]);
            %                 G(i,j)=eval(ex);
            %             end
            G(i,1)=(2*x(1) - 2*SP(1,i))/(2*((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2)^(1/2));
            G(i,2)=(2*x(2) - 2*SP(2,i))/(2*((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2)^(1/2));
            G(i,3)=(2*x(3) - 2*SP(3,i))/(2*((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2)^(1/2));
        end
        h=dis-Fx;
        delta=(pinv(G' * inv(Q) *G) * G' * inv(Q))*h;        
        xAnt=x;
        x=x+delta;
        
        erro=norm(xAnt-x);
        if(erro<=tol)            
            break
        end
        Fx=geraFx(x,SP);%busca medidas com base no valor atualizado
    end
    erro=norm(EM-x);
    if(erro>1000) %Verifica não convergência e adiciona um erro médio do método
        x=EM+15;
    end
    pos=x;
end

function [Fx]=geraFx(x,ES)
    ns=size(ES,2);
    Fx=zeros(ns,1);
    for i=1:ns
        Fx(i)=norm(x-ES(:,i));
    end
end