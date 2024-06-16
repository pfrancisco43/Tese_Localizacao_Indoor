function [indice,h]=DCSSOMPRef(Y,A,L)
    %Versão simplificada do DCS-SOMP que executa apenas uma vez
    %L=1 em todos os casos

    K=size(Y,2);        % Canais
    N=size(Y,1);        % Observações
    M=size(A,2);        
    if (L<=0)
        L=N;
    end
    % Inicialização
    R=Y;
    psi=zeros(N,L,K);
    indices=zeros(1,L);
    columns=zeros(N,L,K);
    betamatrix=zeros(L,K);

    for counter=1:L
        % Encontra a correlação máxima
        cost=zeros(1,M);
        for m=1:M
            for k=1:K
                v=abs(A(:,m,k)'*R(:,k))/norm(A(:,m,k));
                cost(m)=cost(m)+v;
            end
        end

        [maxval,maxi]=max(cost);

        indice=maxi;

        for k=1:K
            % Ortogonalização
            columns(:,counter,k)=A(:,maxi,k);
            omega=A(:,maxi,k);
            psi(:,counter,k)=omega;
                     
            beta=(psi(:,counter,k)'*R(:,k)/(norm(psi(:,counter,k)))^2);
            betamatrix(counter,k)=beta;
        end
    end

    % Retorno para o cálculo do ToA
    h=zeros(L,K);
    for k=1:K
        [Q,Rqr]=qr(columns(:,:,k),0);
        h(:,k)=inv(Rqr)*Q'*psi(:,:,k)*betamatrix(:,k);
    end