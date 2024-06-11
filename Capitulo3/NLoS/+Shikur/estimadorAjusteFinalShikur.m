function [erro,crlb]=estimadorAjusteFinalShikur(EB,EM,ES,toas, aods, aoas,sigmaToa,sigmaAngs,Xi)
    
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
    
    
    x=Xi;%veio do método geométrico
    
    C=3;%Quantidade de colunas (parametros a estimar)
    L=ns*5;%quantidade de medições;
    G=zeros(L,C);
%     
%     
%         [f1]=funcoes(); %buscar jacobiana é opcional se souber a derivada
%         funcs1=repmat(f1(1),[1,ns]);
%         funcs2=repmat(f1(2),[1,ns]);
%         funcs3=repmat(f1(3),[1,ns]);
%         funcs4=repmat(f1(4),[1,ns]);
%         funcs5=repmat(f1(5),[1,ns]);
%         funcs=[funcs1,funcs2,funcs3,funcs4,funcs5];
%         syms mx my mz sx sy sz bx by bz;
%         xx=[mx my mz];
%         Ja=jacobian(funcs,xx);
    tol=1e-4;
    %     si=1;
    for kk=1:1000%quantidade de iterações
        qi=kk;
        Fx=geraFx(x,ES,EB);%busca medidas com base no valor atualizado
        
        G=[buscaGToA(x,ES);buscaGAoDAz(x,ES);buscaGAoDEl(x,ES);buscaGAoAAz(x,ES);buscaGAoAEl(x,ES)];
        
        h=medidas-Fx;
        delta=(pinv(G' * inv(Q) *G) * G' * inv(Q))*h;
        xAnt=x;
        x=x+delta;
        
        erro=norm(xAnt-x);
        if(erro<=tol || erro > 1000)
            break
        end
        
    end
    
    if(erro>1000)
        pos=Xi;
    else
        pos=x;
    end    
    erro=norm(pos-EM);
    
    crlb=detcrlb(Q,EM,ES);
end

function G=buscaGToA(x,ES)
    mx=x(1);
    my=x(2);
    mz=x(3);    
    for i=1:size(ES,2)
        sx=ES(1,i);
        sy=ES(2,i);
        sz=ES(3,i);
        
        G(i,1)=(2*mx - 2*sx)/(2*((mx - sx)^2 + (my - sy)^2 + (mz - sz)^2)^(1/2));
        G(i,2)=(2*my - 2*sy)/(2*((mx - sx)^2 + (my - sy)^2 + (mz - sz)^2)^(1/2));
        G(i,3)=(2*mz - 2*sz)/(2*((mx - sx)^2 + (my - sy)^2 + (mz - sz)^2)^(1/2));
    end
end

function G=buscaGAoAAz(x,ES)
    mx=x(1);
    my=x(2);
    mz=x(3);
    
    for i=1:size(ES,2)
        sx=ES(1,i);
        sy=ES(2,i);
        sz=ES(3,i);
        
        G(i,1)=pi*dirac(mx - sx) - (my - sy)/((mx - sx)^2*((my - sy)^2/(mx - sx)^2 + 1));
        G(i,2)=1/((mx - sx)*((my - sy)^2/(mx - sx)^2 + 1));
        G(i,3)=0;
    end
end

function G=buscaGAoAEl(x,ES)
    mx=x(1);
    my=x(2);
    mz=x(3);
    
    for i=1:size(ES,2)
        sx=ES(1,i);
        sy=ES(2,i);
        sz=ES(3,i);
        
        G(i,1)=-((mz - sz)*(2*mx - 2*sx))/(2*((mx - sx)^2 + (my - sy)^2)^(3/2)*((mz - sz)^2/((mx - sx)^2 + (my - sy)^2) + 1));
        G(i,2)=-((mz - sz)*(2*my - 2*sy))/(2*((mx - sx)^2 + (my - sy)^2)^(3/2)*((mz - sz)^2/((mx - sx)^2 + (my - sy)^2) + 1));
        G(i,3)=1/(((mx - sx)^2 + (my - sy)^2)^(1/2)*((mz - sz)^2/((mx - sx)^2 + (my - sy)^2) + 1));
    end
end

function G=buscaGAoDAz(x,ES)
    for i=1:size(ES,2)
        G(i,1)=0;
        G(i,2)=0;
        G(i,3)=0;
    end
end

function G=buscaGAoDEl(x,ES)
    for i=1:size(ES,2)
        G(i,1)=0;
        G(i,2)=0;
        G(i,3)=0;
    end
end

function [Fx]=geraFx(x,ES,EB)
    
    Fx1=Shikur.buscaTDoAShikur (EB,ES,x);
    
    aods=Shikur.buscaAodShikur(ES,EB);
    Fx2=aods(:,1);
    Fx3=aods(:,2);
    
    aoas=Shikur.buscaAoaShikur(ES,x);
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


function crlb=detcrlb(Q,x,ES)
    G=[buscaGToA(x,ES);buscaGAoDAz(x,ES);buscaGAoDEl(x,ES);buscaGAoAAz(x,ES);buscaGAoAEl(x,ES)];
    FIM=(G' * inv(Q) * G);
    cr=trace(inv(FIM));
    crlb=sqrt(cr);
    if(cr<crlb)
        crlb=cr;
    end
end
