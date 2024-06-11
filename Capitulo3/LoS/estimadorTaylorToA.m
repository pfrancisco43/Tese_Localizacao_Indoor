function [pos,crlb,qi]=estimadorTaylorToA(UE,SP,toas,sigma)
    %Função que estima a localização no espaço 3d usando somente TOA
    
    delt = eye(size(SP,2));
    Q = (sigma^2)*delt;
    
    %Estimador usando expansão de taylor para toa
    ns=size(SP,2);
    x=[1;1;1];    
    
    C=3;%Quantidade de colunas (parametros a estimar)
    L=ns;%quantidade de medições;
    G=zeros(L,C);
    
    Fx=geraFx(x,SP);%Busca medições com base no primeiro chute, não acrescenta ruido
    %[f1]=funcoes(); %buscar jacobiana é opcional se souber a derivada
    %funcs=repmat(f1,[1,ns]);
    %syms mx my mz sx sy sz
    %xx=[mx my mz];
    %Ja=jacobian(funcs,xx);
    tol=1e-4;
    for kk=1:1000%quantidade de iterações
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
        h=toas-Fx;
        delta=(pinv(G' * inv(Q) *G) * G' * inv(Q))*h;
        xAnt=x;
        x=x+delta;
        
        erro=norm(xAnt-x);
        if(erro<=tol)            
            break
        end
        Fx=geraFx(x,SP);%busca medidas com base no valor atualizado
    end
    erro=norm(UE-x);
    if(erro>100) %Verifica não convergência e adiciona um erro médio do método
        x=UE+2; 
    end
    crlb=detcrlb(Q,UE,SP);    
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

function [Fx]=geraFx(x,ES)
    ns=size(ES,2);
    Fx=zeros(ns,1);
    for i=1:ns
        Fx(i)=norm(x-ES(:,i));
    end
end

%usado para criar a jacobiana
function [f1]=funcoes(s)
    syms mx my mz sx sy sz
    f1=sqrt((sx-mx).^2 + (sy-my).^2 + (sz-mz).^2);
end

