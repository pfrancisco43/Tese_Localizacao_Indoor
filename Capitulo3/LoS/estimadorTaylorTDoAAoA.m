function [pos,crlb,qi]=estimadorTaylorTDoAAoA(UE,SP,tdoas, aoas,sigmaTdoa,sigmaAoa)
    %Função que estima a localização no espaço 3d usando somente TDoA+AoA e
    %o método de Taylor
    
    ns=size(SP,2);
    medidas=[tdoas;aoas(:,1);aoas(:,2)];
    sigmaTdoa=repmat(sigmaTdoa,ns-1,1);
    sigmaAoa=repmat(sigmaAoa,ns*2,1);
    sigma=[sigmaTdoa;sigmaAoa];
    delt = eye(ns*3-1);
    Q = (sigma.^2).*delt;    
    
    x=[1;1;1];    
    tol=1e-4;    
    for kk=1:1000%quantidade de iterações
        qi=kk;
        Fx=geraFx(x,SP);%busca medidas com base no valor atualizado
        G=[buscaGTdoA(x,SP);buscaGAoAAz(x,SP);buscaGAoAEl(x,SP)];
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

function G=buscaGTdoA(x,ES)
    EB=ES(:,1);
    ES=ES(:,2:end);    
    for i=1:size(ES,2)
        G(i,1)=(2*x(1) - 2*ES(1,i))/(2*((x(1) - ES(1,i))^2 + (x(2) - ES(2,i))^2 + (x(3) - ES(3,i))^2)^(1/2)) + (2*EB(1) - 2*x(1))/(2*((EB(1) - x(1))^2 + (EB(2) - x(2))^2 + (EB(3) - x(3))^2)^(1/2));
        G(i,2)=(2*x(2) - 2*ES(2,i))/(2*((x(1) - ES(1,i))^2 + (x(2) - ES(2,i))^2 + (x(3) - ES(3,i))^2)^(1/2)) + (2*EB(2) - 2*x(2))/(2*((EB(1) - x(1))^2 + (EB(2) - x(2))^2 + (EB(3) - x(3))^2)^(1/2));
        G(i,3)=(2*x(3) - 2*ES(3,i))/(2*((x(1) - ES(1,i))^2 + (x(2) - ES(2,i))^2 + (x(3) - ES(3,i))^2)^(1/2)) + (2*EB(3) - 2*x(3))/(2*((EB(1) - x(1))^2 + (EB(2) - x(2))^2 + (EB(3) - x(3))^2)^(1/2));
    end
end

function G=buscaGAoAAz(x,SP)
    for i=1:size(SP,2)
        G(i,1)=-(x(2) - SP(2,i))/((x(1) - SP(1,i))^2*((x(2) - SP(2,i))^2/(x(1) - SP(1,i))^2 + 1));
        G(i,2)=1/((x(1) - SP(1,i))*((x(2) - SP(2,i))^2/(x(1) - SP(1,i))^2 + 1));
        G(i,3)=0;
    end
end

function G=buscaGAoAEl(x,SP)
    for i=1:size(SP,2)
        G(i,1)=((x(3) - SP(3,i))*(2*x(1) - 2*SP(1,i)))/(2*(1 - (x(3) - SP(3,i))^2/((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2))^(1/2)*((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2)^(3/2));
        G(i,2)=((x(3) - SP(3,i))*(2*x(2) - 2*SP(2,i)))/(2*(1 - (x(3) - SP(3,i))^2/((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2))^(1/2)*((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2)^(3/2));
        G(i,3)=-(1/((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2)^(1/2) - ((x(3) - SP(3,i))*(2*x(3) - 2*SP(3,i)))/(2*((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2)^(3/2)))/(1 - (x(3) - SP(3,i))^2/((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2))^(1/2);
    end
end

function [Fx]=geraFx(x,ES)
    EB=ES(:,1);
    SP=ES(:,2:end);
    ns=size(SP,2);
    db=norm(x-EB);
    Fx1=zeros(ns,1);
    for i=1:ns
        ds=norm(x-SP(:,i));
        Fx1(i)=ds-db;
    end
    
    ns=size(ES,2);
    Fx2=zeros(ns,1);
    Fx3=zeros(ns,1);
    for i=1:ns
        Fx2(i)=atan((x(2)-ES(2,i))/(x(1)-ES(1,i)));
        Fx3(i)=acos((x(3)-ES(3,i))/norm(ES(:,i)-x));
    end
    Fx=[Fx1;Fx2;Fx3];
end

