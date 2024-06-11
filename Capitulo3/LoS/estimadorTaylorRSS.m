function [pos,crlb,qi]=estimadorTaylorRSS(EM,ES,rsss,sigma,p00,beta0)
    %Método de Taylor com RSS
    
    delt = eye(size(ES,2));
    Q = (sigma^2)*delt;
    
    ns=size(ES,2);
    x=[1;1;1];
    
    C=3;%Quantidade de colunas (parametros a estimar)
    L=ns;%quantidade de medições;
    G=zeros(L,C);
    
%     [f1]=funcoes(); %buscar jacobiana é opcional se souber a derivada
%     funcs=repmat(f1,[1,ns]);
%     syms mx my mz sx sy sz p0 bet
%     xx=[mx my mz];
%     Ja=jacobian(funcs,xx);
    tol=1e-4;
    for kk=1:1000%quantidade de iterações
        Fx=geraFx(x,ES,p00,beta0);%Busca medições com base nos valores estimados, não acrescenta ruido
        qi=kk;
        for i = 1:L
%             for j = 1:C
%                 f=Ja(i,j);
%                 ex=subs(f,[mx my mz sx sy sz p0 bet],[x',ES(1,i),ES(2,i),ES(3,i), p00, beta0]);
%                 G(i,j)=eval(ex);
%             end
            G(i,1)=-(5*beta0*(2*x(1) - 2*ES(1,i)))/(log(10)*((x(1) - ES(1,i))^2 + (x(2) - ES(2,i))^2 + (x(3) - ES(3,i))^2));
            G(i,2)=-(5*beta0*(2*x(2) - 2*ES(2,i)))/(log(10)*((x(1) - ES(1,i))^2 + (x(2) - ES(2,i))^2 + (x(3) - ES(3,i))^2));
            G(i,3)=-(5*beta0*(2*x(3) - 2*ES(3,i)))/(log(10)*((x(1) - ES(1,i))^2 + (x(2) - ES(2,i))^2 + (x(3) - ES(3,i))^2));
        end
        h=rsss-Fx;
        delta=(pinv(G' * inv(Q) *G) * G' * inv(Q))*h;
        xAnt=x;
        x=x+delta;
        
        erro=norm(xAnt-x);
        if(erro<=tol)            
            break
        end       
    end
    erro=norm(EM-x);
    if(erro>1000) %Verifica não convergência e adiciona um erro médio do método
        x=EM+15;
        eqm=norm(x-EM);       
    end
    crlb=0;    
    pos=x;    
end

function crlb=detcrlb(Q,x,ES)
    for i = 1:size(ES,2)
        G(i,1)=(2*x(1) - 2*ES(1,i))/(2*((x(1) - ES(1,i))^2 + (x(2) - ES(2,i))^2 + (x(3) - ES(3,i))^2)^(1/2));
        G(i,2)=(2*x(2) - 2*ES(2,i))/(2*((x(1) - ES(1,i))^2 + (x(2) - ES(2,i))^2 + (x(3) - ES(3,i))^2)^(1/2));
        G(i,3)=(2*x(3) - 2*ES(3,i))/(2*((x(1) - ES(1,i))^2 + (x(2) - ES(2,i))^2 + (x(3) - ES(3,i))^2)^(1/2));
    end
    FIM=(G' * inv(Q) * G);
    crlb=trace(inv(FIM));
    crlb=sqrt(crlb);
end

function [Fx]=geraFx(x,ES,p0,bet)
    ns=size(ES,2);
    Fx=zeros(ns,1);
    for i=1:ns
        Fx(i)=p0-10*bet*log10(norm(x-ES(:,i)));
    end
end

%usado para criar a jacobiana
function [f1]=funcoes()
    syms mx my mz sx sy sz p0 bet
    f1=p0-10*bet*log10(sqrt((sx-mx).^2 + (sy-my).^2 + (sz-mz).^2));    
end

