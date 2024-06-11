function [pos,crlb,qi]=estimadorTaylorToATDoA(UE,SP,toas,tdoas,sigma)
    %Função que estima a localização no espaço 3d usando somente ToA+TDoA e
    %método de Taylor
    
    EB=SP(:,1);
    ES=SP(:,2:end);
    medidas=[toas;tdoas];
    delt = eye(length(medidas));
    Q = (sigma^2)*delt;
    
    %Estimador usando expansão de taylor para toa
    ns=size(SP,2);
    x=[1;1;1];
    
    C=3;%Quantidade de colunas (parametros a estimar)
    L=ns;%quantidade de medições;
    G=zeros(length(medidas),C);
        
    tol=1e-4;
    for kk=1:100%quantidade de iterações
        qi=kk;
        Fx1=geraFxToA(x,SP);%Toa
        Fx2=geraFxTDoA(x,ES,EB);%tdoa
        Fx=[Fx1;Fx2];%junção
        for i = 1:L
            G(i,1)=(2*x(1) - 2*SP(1,i))/(2*((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2)^(1/2));
            G(i,2)=(2*x(2) - 2*SP(2,i))/(2*((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2)^(1/2));
            G(i,3)=(2*x(3) - 2*SP(3,i))/(2*((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2)^(1/2));
        end
        for i =1:L-1
            G(L+i,1)=(2*x(1) - 2*ES(1,i))/(2*((x(1) - ES(1,i))^2 + (x(2) - ES(2,i))^2 + (x(3) - ES(3,i))^2)^(1/2)) + (2*EB(1) - 2*x(1))/(2*((EB(1) - x(1))^2 + (EB(2) - x(2))^2 + (EB(3) - x(3))^2)^(1/2));
            G(L+i,2)=(2*x(2) - 2*ES(2,i))/(2*((x(1) - ES(1,i))^2 + (x(2) - ES(2,i))^2 + (x(3) - ES(3,i))^2)^(1/2)) + (2*EB(2) - 2*x(2))/(2*((EB(1) - x(1))^2 + (EB(2) - x(2))^2 + (EB(3) - x(3))^2)^(1/2));
            G(L+i,3)=(2*x(3) - 2*ES(3,i))/(2*((x(1) - ES(1,i))^2 + (x(2) - ES(2,i))^2 + (x(3) - ES(3,i))^2)^(1/2)) + (2*EB(3) - 2*x(3))/(2*((EB(1) - x(1))^2 + (EB(2) - x(2))^2 + (EB(3) - x(3))^2)^(1/2));
        end
        h=medidas-Fx;
        delta=(pinv(G' * inv(Q) *G) * G' * inv(Q))*h;
        xAnt=x;
        x=x+delta;
        
        erro=norm(xAnt-x);
        if(erro<=tol)
            break
        end        
    end
    erro=norm(UE-x);
    if(erro>100) %Verifica não convergência e adiciona um erro médio do método
        x=UE+2;
    end
    crlb=0;
    pos=x;
end

function crlb=detcrlb(Q,x,SP)
    for i = 1:size(SP,2)
        G(i,1)=(2*x(1) - 2*SP(1,i))/(2*((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2)^(1/2));
        G(i,2)=(2*x(2) - 2*SP(2,i))/(2*((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2)^(1/2));
        G(i,3)=(2*x(3) - 2*SP(3,i))/(2*((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2)^(1/2));
    end
    
    FIM=(G' * inv(Q) * G);
    crlb=trace(inv(FIM));
    crlb=sqrt(crlb);
end

function [Fx]=geraFxToA(x,ES)
    ns=size(ES,2);
    Fx=zeros(ns,1);
    for i=1:ns
        Fx(i)=norm(x-ES(:,i));
    end
end

function [Fx,db]=geraFxTDoA(x,ES,EB)
    db=norm(x-EB);
    ns=size(ES,2);
    Fx=zeros(ns,1);
    for i=1:ns
        ds=norm(x-ES(:,i));
        Fx(i,1)=ds-db;
    end
end

