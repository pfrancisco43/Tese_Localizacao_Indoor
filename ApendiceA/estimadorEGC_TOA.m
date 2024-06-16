%Estima localização em um cenário LOS 
%Usa equação geral da circunferência e não precisa de chute inicial
%Autor: Paulo Francisco
function [e,pos]=estimadorEGC_TOA (EM,ES,toas)
    ns=size(ES,2);    
    if(ns==4)
        ns=ns-1;
    end
    dim=size(ES,1);
    func='interEq';
    if(dim==3)
        func='interEq3D';
    end
    
    h=nchoosek((1:ns),3);%Cria matriz de possibilidades
    m=size(h,1);
    A=[]; b=[];
    %busca as equações com base em 3 espalhadores e vai montando o conjunto
    %de equações
    for i=1:m
        [AA,bb]=feval(func,ES(:,h(i,1)),ES(:,h(i,2)),ES(:,h(i,3)),toas(h(i,1)),toas(h(i,2)),toas(h(i,3)));
        A=[A;AA];
        b=[b;bb];
    end
    pos=pinv(A'*A)*A'*b;%Resolve o sistema de equaçoes sobre-determinado
    %pos=[0,0,0]';
    %OPCIONAL - Minimização do erro
    % Configuração das opções de otimização
    options = optimoptions('lsqnonlin', 'Display', 'off', 'Algorithm', 'trust-region-reflective','OptimalityTolerance',0.19);
    % % Função de erro para minimização
    error_function = @(x) arrayfun(@(i) norm(x - ES(:, i)) - toas(i), 1:size(ES, 2));
    % % Chamada ao otimizador
    posf= lsqnonlin(error_function, pos, [], [], options);

    posf(3)=rand;
    e=norm(posf-EM);
end

%Usando 3 BSs e 3 equações, em 2D
function [A,b]=interEq(s1,s2,s3,r1,r2,r3)
    A1=[-2*s1(1)+2*s2(1) , -2*s1(2)+2*s2(2)];
    A2=[-2*s1(1)+2*s3(1) , -2*s1(2)+2*s3(2)];
    b1=r1^2-r2^2-s1(1)^2+s2(1)^2-s1(2)^2+s2(2)^2;
    b2=r1^2-r3^2-s1(1)^2+s3(1)^2-s1(2)^2+s3(2)^2;
    A=[A1;A2];
    b=[b1;b2];
end

%Usando 3 BSs e 2 equações
function [A,b]=interEq3D(s1,s2,s3,r1,r2,r3)
    A1=[-2*s1(1)+2*s2(1) , -2*s1(2)+2*s2(2), -2*s1(3)+2*s2(3)];
    A2=[-2*s1(1)+2*s3(1) , -2*s1(2)+2*s3(2), -2*s1(3)+2*s3(3)];
    b1=r1^2-r2^2 -s1(1)^2+s2(1)^2 -s1(2)^2+s2(2)^2 -s1(3)^2+s2(3)^2;
    b2=r1^2-r3^2 -s1(1)^2+s3(1)^2 -s1(2)^2+s3(2)^2 -s1(3)^2+s3(3)^2;
    A=[A1;A2];
    b=[b1;b2];
end

%Usando 3 BSs e 3 equações
function [A,b]=interEq3D_2(s1,s2,s3,r1,r2,r3)
    A1=[s2(1)-s1(1) , s2(2)-s1(2), s2(3)-s1(3)];
    A2=[s3(1)-s1(1) , s3(2)-s1(2), s3(3)-s1(3)];
    A3=[s3(1)-s2(1) , s3(2)-s2(2), s3(3)-s2(3)];
    b1=(r1^2-r2^2 -s1(1)^2+s2(1)^2 -s1(2)^2+s2(2)^2 -s1(3)^2+s2(3)^2)/2;
    b2=(r1^2-r3^2 -s1(1)^2+s3(1)^2 -s1(2)^2+s3(2)^2 -s1(3)^2+s3(3)^2)/2;
    b3=(r2^2-r3^2 -s2(1)^2+s3(1)^2 -s2(2)^2+s3(2)^2 -s2(3)^2+s3(3)^2)/2;
    A=[A1;A2;A3];
    b=[b1;b2;b3];
end

