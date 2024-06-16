%Esta função minimiza o erro de estimação baseado no método de Gauss-Newton
%ou Taylor de primeira ordem
%Autor: Paulo Francisco
function [erro,crlb,posEM,posES]=estimadorAjusteFinalMelhorado(EB,EM,ES,toas, aods, aoas,sigmaToa,sigmaAngs,Xi,Si)   
    ns=size(ES,2);
    sigmaAoa=sigmaAngs;
    sigmaAod=sigmaAngs;
    medidas=[toas;aods(:,1);aods(:,2);aoas(:,1);aoas(:,2)];
    sigmaToa=repmat(sigmaToa,ns,1);
    sigmaAoa=repmat(sigmaAoa,ns*2,1);
    sigmaAod=repmat(sigmaAod,ns*2,1);
    sigma=[sigmaToa;sigmaAoa;sigmaAod];
    
    delt = eye(ns*5);
    Q = (sigma.^2).*delt;    
    
    x=Xi; %o primeiro é o da EM
    for i =1: ns
        x=[x;Si(:,i)];
    end
    xo=x;
    
    C=(ns*3)+3; %Quantidade de colunas (número de espalhadores *3 + 3 da EM)
    L=ns*5;%quantidade de medições;
    G=zeros(L,C);
    
    % [f1]=funcoes(); %buscar jacobiana é opcional se souber a derivada
    % funcs1=repmat(f1(1),[1,ns]);
    % funcs2=repmat(f1(2),[1,ns]);
    % funcs3=repmat(f1(3),[1,ns]);
    % funcs4=repmat(f1(4),[1,ns]);
    % funcs5=repmat(f1(5),[1,ns]);
    % funcs=[funcs1,funcs2,funcs3,funcs4,funcs5];
    % syms mx my mz sx sy sz bx by bz;
    % %xx=[mx my mz sx sy sz];
    % xx=[mx my mz];
    % Ja=jacobian(funcs,xx);
    
    tol=1e-4; 
 
    for kk=1:1000%quantidade de iterações        
        Fx=geraFx(x,EB);%busca medidas com base no valor atualizado
        
        G=[buscaGToA(x,EB);buscaGAoDAz(x,EB);buscaGAoDEl(x,EB);buscaGAoAAz(x);buscaGAoAEl(x)];
        
        h=medidas-Fx;
        %delta=(pinv(G' *G) * G')*h; %Gauss
        delta=(pinv(G' * inv(Q) *G) * G' * inv(Q))*h; %Usando a covariância - MLE
        xAnt=x;
        x=x+delta;
        
        erro=norm(xAnt-x);
        if(erro<=tol || erro > 1000)
            break
        end
        
    end
 
    if(erro>10)
        x=xo;       
    end
    
    posEM=x(1:3);
    posES=x(4:end);
    ns=size(posES,1)/3;%qnt espalhadores
    posES=reshape(posES,[3,ns]);
   
    erro=norm(posEM-EM);
    
    crlb=0;
end

function G=buscaGToA(X,B)
    x=X(1:3);
    ES=X(4:end);    
    ns=size(ES,1)/3;%qnt espalhadores
    C=(ns+1)*3;%qnt de colunas da jacobiana
    ES=reshape(ES,[3,ns]);

    mx=x(1);
    my=x(2);
    mz=x(3);
    bx=B(1);
    by=B(2);
    bz=B(3);
    
    G=zeros(ns,C);
    for i=1:ns
        sx=ES(1,i);
        sy=ES(2,i);
        sz=ES(3,i);
        
        G(i,1)=(2*mx - 2*sx)/(2*((mx - sx)^2 + (my - sy)^2 + (mz - sz)^2)^(1/2));               
        G(i,2)=(2*my - 2*sy)/(2*((mx - sx)^2 + (my - sy)^2 + (mz - sz)^2)^(1/2));               
        G(i,3)=(2*mz - 2*sz)/(2*((mx - sx)^2 + (my - sy)^2 + (mz - sz)^2)^(1/2));               
        
        %Aqui é zerado onde não tem nada a ver
        G(i,i*3+1)= - (2*mx - 2*sx)/(2*((mx - sx)^2 + (my - sy)^2 + (mz - sz)^2)^(1/2)) - (2*bx - 2*sx)/(2*((bx - sx)^2 + (by - sy)^2 + (bz - sz)^2)^(1/2));
        G(i,i*3+2)= - (2*my - 2*sy)/(2*((mx - sx)^2 + (my - sy)^2 + (mz - sz)^2)^(1/2)) - (2*by - 2*sy)/(2*((bx - sx)^2 + (by - sy)^2 + (bz - sz)^2)^(1/2));
        G(i,i*3+3)= - (2*mz - 2*sz)/(2*((mx - sx)^2 + (my - sy)^2 + (mz - sz)^2)^(1/2)) - (2*bz - 2*sz)/(2*((bx - sx)^2 + (by - sy)^2 + (bz - sz)^2)^(1/2));
    end
end

function G=buscaGAoAAz(X)
    x=X(1:3);
    ES=X(4:end);
    ns=size(ES,1)/3;%qnt espalhadores
    C=(ns+1)*3;%qnt de colunas da jacobiana
    ES=reshape(ES,[3,ns]);
    
    mx=x(1);
    my=x(2);
    mz=x(3);

    G=zeros(ns,C);
    for i=1:ns
        sx=ES(1,i);
        sy=ES(2,i);
        sz=ES(3,i);
        
        G(i,1)=pi*dirac(mx - sx) - (my - sy)/((mx - sx)^2*((my - sy)^2/(mx - sx)^2 + 1));               
        G(i,2)=1/((mx - sx)*((my - sy)^2/(mx - sx)^2 + 1));
        G(i,3)=0;
        %Aqui é zerado onde não tem nada a ver        
        G(i,i*3+1)= (my - sy)/((mx - sx)^2*((my - sy)^2/(mx - sx)^2 + 1)) - pi*dirac(mx - sx);
        G(i,i*3+2)= -1/((mx - sx)*((my - sy)^2/(mx - sx)^2 + 1));
        G(i,i*3+3)= 0;
    end
end

function G=buscaGAoAEl(X)
    x=X(1:3);
    ES=X(4:end);
    ns=size(ES,1)/3;%qnt espalhadores
    C=(ns+1)*3;%qnt de colunas da jacobiana
    ES=reshape(ES,[3,ns]);
    
    mx=x(1);
    my=x(2);
    mz=x(3);

    G=zeros(ns,C);
    for i=1:size(ES,2)
        sx=ES(1,i);
        sy=ES(2,i);
        sz=ES(3,i);
        
        G(i,1)=-((mz - sz)*(2*mx - 2*sx))/(2*((mx - sx)^2 + (my - sy)^2)^(3/2)*((mz - sz)^2/((mx - sx)^2 + (my - sy)^2) + 1));
        G(i,2)=-((mz - sz)*(2*my - 2*sy))/(2*((mx - sx)^2 + (my - sy)^2)^(3/2)*((mz - sz)^2/((mx - sx)^2 + (my - sy)^2) + 1));
        G(i,3)=1/(((mx - sx)^2 + (my - sy)^2)^(1/2)*((mz - sz)^2/((mx - sx)^2 + (my - sy)^2) + 1));
        %Aqui é zerado onde não tem nada a ver
        G(i,i*3+1)= ((mz - sz)*(2*mx - 2*sx))/(2*((mx - sx)^2 + (my - sy)^2)^(3/2)*((mz - sz)^2/((mx - sx)^2 + (my - sy)^2) + 1));
        G(i,i*3+2)= ((mz - sz)*(2*my - 2*sy))/(2*((mx - sx)^2 + (my - sy)^2)^(3/2)*((mz - sz)^2/((mx - sx)^2 + (my - sy)^2) + 1));
        G(i,i*3+3)= -1/(((mx - sx)^2 + (my - sy)^2)^(1/2)*((mz - sz)^2/((mx - sx)^2 + (my - sy)^2) + 1));
    end
end

function G=buscaGAoDAz(X,B)
    x=X(1:3);
    ES=X(4:end);
    ns=size(ES,1)/3;%qnt espalhadores
    C=(ns+1)*3;%qnt de colunas da jacobiana
    ES=reshape(ES,[3,ns]);
    
    mx=x(1);
    my=x(2);
    mz=x(3);
    bx=B(1);
    by=B(2);
    bz=B(3);
    
    G=zeros(ns,C);
    for i=1:ns
        sx=ES(1,i);
        sy=ES(2,i);
        sz=ES(3,i);
        
        G(i,1)=0;
        G(i,2)=0;
        G(i,3)=0;
        %Aqui é zerado onde não tem nada a ver
        if dirac(bx - sx) ~= 0
            hhh=4;
        end
        G(i,i*3+1)= (by - sy)/((bx - sx)^2*((by - sy)^2/(bx - sx)^2 + 1)) - pi*dirac(bx - sx);
        G(i,i*3+2)= -1/((bx - sx)*((by - sy)^2/(bx - sx)^2 + 1));
        G(i,i*3+3)= 0;
    end
end

function G=buscaGAoDEl(X,B)
    x=X(1:3);
    ES=X(4:end);
    ns=size(ES,1)/3;%qnt espalhadores
    C=(ns+1)*3;%qnt de colunas da jacobiana
    ES=reshape(ES,[3,ns]);
    
    mx=x(1);
    my=x(2);
    mz=x(3);
    bx=B(1);
    by=B(2);
    bz=B(3);
    
    G=zeros(ns,C);
    for i=1:ns
        sx=ES(1,i);
        sy=ES(2,i);
        sz=ES(3,i);
        
        G(i,1)=0;
        G(i,2)=0;
        G(i,3)=0;
        %Aqui é zerado onde não tem nada a ver
        G(i,i*3+1)= ((bz - sz)*(2*bx - 2*sx))/(2*((bx - sx)^2 + (by - sy)^2)^(3/2)*((bz - sz)^2/((bx - sx)^2 + (by - sy)^2) + 1));
        G(i,i*3+2)= ((bz - sz)*(2*by - 2*sy))/(2*((bx - sx)^2 + (by - sy)^2)^(3/2)*((bz - sz)^2/((bx - sx)^2 + (by - sy)^2) + 1));
        G(i,i*3+3)= -1/(((bx - sx)^2 + (by - sy)^2)^(1/2)*((bz - sz)^2/((bx - sx)^2 + (by - sy)^2) + 1));
    end
end

function [Fx]=geraFx(X,EB)
    x=X(1:3);
    ES=X(4:end);
    ns=size(ES,1)/3;
    ES=reshape(ES,[3,ns]);
    Fx1=buscaToAsNlos (EB,ES,x);
    
    aods=buscaAoD(EB,ES);
    Fx2=aods(:,1);
    Fx3=aods(:,2);
    
    aoas=buscaAoA(ES,x);
    Fx4=aoas(:,1);
    Fx5=aoas(:,2);
    
    Fx=[Fx1;Fx2;Fx3;Fx4;Fx5];
end

%usado para criar a jacobiana
function [f]=funcoes()
    syms mx my mz sx sy sz bx by bz
    f1=sqrt((sx-bx).^2 + (sy-by).^2 + (sz-bz).^2) + sqrt((sx-mx).^2 + (sy-my).^2 + (sz-mz).^2);
    
    f2=(pi/2) * (1 - sign(sx - bx)) + atan((sy - by) / (sx - bx));
    f3=(pi/2) - atan((sz - bz) / sqrt((sx - bx)^2 + (sy - by)^2));
    
    f4=(pi/2) * (1 - sign(sx - mx)) + atan((sy - my) / (sx - mx));
    f5=(pi/2) - atan((sz - mz) / sqrt((sx - mx)^2 + (sy - my)^2));
    
    f=[f1,f2,f3,f4,f5];
end


function crlb=detcrlb(Q,EM,ES,EB)
    ns=size(ES,2);
    x=[EM;reshape(ES,[ns*3,1])]; %o primeiro é o da EM
    
    G=[buscaGToA(x,EB);buscaGAoDAz(x,EB);buscaGAoDEl(x,EB);buscaGAoAAz(x);buscaGAoAEl(x)];
    FIM=(G' * inv(Q) * G);

    FIM_inv = inv(FIM); % Inverte a matriz de informação de Fisher

    crlbs = diag(FIM_inv); % Extrai a diagonal da matriz inversa, que são os CRLBs individuais
    crlbs = sqrt(crlbs); % Tira a raiz quadrada para obter a estimativa do erro padrão

    cr=trace(inv(FIM));
    crlb=sqrt(cr);
   
end
