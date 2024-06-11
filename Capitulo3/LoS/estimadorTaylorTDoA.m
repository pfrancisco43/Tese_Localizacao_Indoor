function [pos,crlb,qi]=estimadorTaylorTDoA(UE,ES,tdoas,sigma)
    %Função que estima a localização no espaço 3d usando somente TDoA com o
    %método de Taylor
    
    EB=ES(:,1);
    ES=ES(:,2:end);
    ns=size(ES,2);
    delt = eye(ns);
    Q = (sigma^2)*delt;    
    
    x=[1;1;1];
    
    C=3;%Quantidade de colunas (parametros a estimar)
    L=ns;%quantidade de medições;
    G=zeros(L,C);
    
         Fx=geraFx(x,ES,EB);%Busca medições com base no primeiro chute, não acrescenta ruido
    %     [f1]=funcoes(); %buscar jacobiana é opcional se souber a derivada
    %     funcs=repmat(f1,[1,ns]);
    %     syms mx my mz sx sy sz bx by bz
    %     xx=[mx my mz];
    %     Ja=jacobian(funcs,xx);
    tol=1e-4;
    for kk=1:1000%quantidade de iterações
        qi=kk;
        for i = 1:L
            %             for j = 1:C
            %                 f=Ja(i,j);
            %                 ex=subs(f,[mx my mz sx sy sz],[x',SP(1,i),SP(2,i),SP(3,i)]);
            %                 G(i,j)=eval(ex);
            %             end
            
            G(i,1)=(2*x(1) - 2*ES(1,i))/(2*((x(1) - ES(1,i))^2 + (x(2) - ES(2,i))^2 + (x(3) - ES(3,i))^2)^(1/2)) + (2*EB(1) - 2*x(1))/(2*((EB(1) - x(1))^2 + (EB(2) - x(2))^2 + (EB(3) - x(3))^2)^(1/2));
            G(i,2)=(2*x(2) - 2*ES(2,i))/(2*((x(1) - ES(1,i))^2 + (x(2) - ES(2,i))^2 + (x(3) - ES(3,i))^2)^(1/2)) + (2*EB(2) - 2*x(2))/(2*((EB(1) - x(1))^2 + (EB(2) - x(2))^2 + (EB(3) - x(3))^2)^(1/2));
            G(i,3)=(2*x(3) - 2*ES(3,i))/(2*((x(1) - ES(1,i))^2 + (x(2) - ES(2,i))^2 + (x(3) - ES(3,i))^2)^(1/2)) + (2*EB(3) - 2*x(3))/(2*((EB(1) - x(1))^2 + (EB(2) - x(2))^2 + (EB(3) - x(3))^2)^(1/2));
            
        end
        h=tdoas-Fx;
        delta=(pinv(G' * inv(Q) *G) * G' * inv(Q))*h;
        xAnt=x;
        x=x+delta;
        
        erro=norm(xAnt-x);
        if(erro<=tol)
            break
        end
        Fx=geraFx(x,ES,EB);%busca medidas com base no valor atualizado
    end
    erro=norm(UE-x);
    if(erro>100) %Verifica não convergência e adiciona um erro médio do método
        x=UE+20; 
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

function [Fx,db]=geraFx(x,ES,EB)
    db=norm(x-EB);
    ns=size(ES,2);
    Fx=zeros(ns,1);
    for i=1:ns
        ds=norm(x-ES(:,i));
        Fx(i,1)=ds-db;
    end
end

%usado para criar a jacobiana
function [f1]=funcoes()
    syms mx my mz sx sy sz bx by bz
    f1=sqrt((sx-mx).^2 + (sy-my).^2 + (sz-mz).^2)-sqrt((bx-mx).^2 + (by-my).^2 + (bz-mz).^2);
end

