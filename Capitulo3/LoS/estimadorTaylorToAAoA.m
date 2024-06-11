function [pos,crlb,qi]=estimadorTaylorToAAoA(UE,SP,toas, aoas,sigmaToa,sigmaAoa)
    %Função que estima a localização no espaço 3d usando somente ToA + AoA
    %com método de Taylor
    ns=size(SP,2);
    medidas=[toas;aoas(:,1);aoas(:,2)];
    sigmaToa=repmat(sigmaToa,ns,1);
    sigmaAoa=repmat(sigmaAoa,ns*2,1);
    sigma=[sigmaToa;sigmaAoa];
    delt = eye(ns*3);
    Q = (sigma.^2).*delt;    
    
    x=[1;1;1];
    
    C=3;%Quantidade de colunas (parametros a estimar)
    L=ns*3;%quantidade de medições;
    G=zeros(L,C);
    
   
%     [f1]=funcoes(); %buscar jacobiana é opcional se souber a derivada
%     funcs1=repmat(f1(1),[1,ns]);
%     funcs2=repmat(f1(2),[1,ns]);
%     funcs3=repmat(f1(3),[1,ns]);
%     funcs=[funcs1,funcs2,funcs3];
%     syms mx my mz sx sy sz
%     xx=[mx my mz];
%     Ja=jacobian(funcs,xx);
     tol=1e-4;
%     si=1;
    for kk=1:1000%quantidade de iterações
        qi=kk;
        Fx=geraFx(x,SP);%busca medidas com base no valor atualizado
%         for i = 1:L
%             for j = 1:C
%                 f=Ja(i,j);
%                 ex=subs(f,[mx my mz sx sy sz],[x',SP(1,si),SP(2,si),SP(3,si)]);
%                 G(i,j)=eval(ex);
%             end
%             si=si+1;
%             if si==ns+1
%                 si=1;
%             end
%         end
        G=[buscaGToA(x,SP);buscaGAoAAz(x,SP);buscaGAoAEl(x,SP)];
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

function G=buscaGToA(x,SP)
    for i=1:size(SP,2)
        G(i,1)=(2*x(1) - 2*SP(1,i))/(2*((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2)^(1/2));
        G(i,2)=(2*x(2) - 2*SP(2,i))/(2*((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2)^(1/2));
        G(i,3)=(2*x(3) - 2*SP(3,i))/(2*((x(1) - SP(1,i))^2 + (x(2) - SP(2,i))^2 + (x(3) - SP(3,i))^2)^(1/2));
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

function crlb=detcrlb(Q,x,SP)
    for i = 1:size(SP,2)
       
    end
    FIM=(G' * inv(Q) * G);
    crlb=trace(inv(FIM));
    crlb=sqrt(crlb);
end

function [Fx]=geraFx(x,ES)
    ns=size(ES,2);
    Fx1=zeros(ns,1);
    Fx2=zeros(ns,1);
    Fx3=zeros(ns,1);
    for i=1:ns
        Fx1(i)=norm(x-ES(:,i));
%         Fx2(i)=atan2(ES(2,i)-x(2),ES(1,i)-x(1));
%         Fx3(i)=asin((ES(3,i)-x(3))/Fx1(i));
        Fx2(i)=atan((x(2)-ES(2,i))/(x(1)-ES(1,i)));
        Fx3(i)=acos((x(3)-ES(3,i))/norm(ES(:,i)-x));
    end
    Fx=[Fx1;Fx2;Fx3];
end

%usado para criar a jacobiana
function [f]=funcoes()
    syms mx my mz sx sy sz
    f1=sqrt((sx-mx).^2 + (sy-my).^2 + (sz-mz).^2);
    %f2=atan2(sy-my,sx-mx);%aoa az
    f2=atan((my-sy)/(mx-sx));
    %f3=asin((sz-mz)/sqrt((sx-mx).^2 + (sy-my).^2 + (sz-mz).^2));%aoa el
    f3=acos((mz-sz)/sqrt((sx-mx).^2 + (sy-my).^2 + (sz-mz).^2));
    f=[f1,f2,f3];
end

