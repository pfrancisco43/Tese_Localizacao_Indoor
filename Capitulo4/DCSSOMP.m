function [aodEst,aoaEst,toaEst,its]=DCSSOMP(Y,A,L,rangeAz_aod,rangeEl_aod,rangeAz_aoa,rangeEl_aoa,x,Nt, Nr, Ts, c,it,lim,fc,arranjo)
    % Este método recebe a matriz de busca e o sinal vetorizado
    % O métode determina a covariância máxima e o resídua
    % O método executa enquanto tiver caminhos válidos
    % O método pode ser executado na sua versão simples (sem considerar iterações)
    % O Método pode ser executado na sua versão adaptativa (considerando
    % mais iterações), nesse caso é executada a função fining

    K=size(Y,2);        % Canais
    N=size(Y,1);        % Observações por canal
    M=size(A,2);        % Tamanho da esparcidade
    if (L<=0)
        L=N;
    end
    % Inicialização
    R=Y; % Resíduo
    psi=zeros(N,L,K);
    indices=zeros(1,L);
    columns=zeros(N,L,K);
    betamatrix=zeros(L,K);

    aodEst=zeros(it+1,2,L);
    aoaEst=zeros(it+1,2,L);
    toaEst=zeros(it+1,L);

    
    ratio_prev = inf;
    limiarQueda = 0.2; % Limiar para a queda da relação
    %% Determina a correlação máxima entre resíduo e matriz de detecção
    for l=1:L        
        cost=zeros(1,M);
        for m=1:M
            for k=1:K
                v=abs(A(:,m,k)'*R(:,k))/norm(A(:,m,k));
                cost(m)=cost(m)+v;
            end
        end

        [maxval,maxi]=max(cost);
        mean_cost = mean(cost);
        ratio = maxval / mean_cost;
        %% Determina se ainda há caminhos válidos
        % Verifica a condição de parada baseada na relação
        % if ratio_prev~= inf
        %     PercQueda=(ratio_prev - ratio)/ratio_prev;
        %     if PercQueda > limiarQueda
        %         break; % Parar se a relação cai abaixo de um certo percentual da razão anterior
        %     end
        % end
        % ratio_prev = ratio;  % Atualiza a razão anterior
        %plot(cost)

        %% Estima os ângulos com base na correlação máxima (2D-AoD, 2D-AoA)
        [aodEi,aoaEi]=estimateAngles(maxi,rangeAz_aod,rangeEl_aod,rangeAz_aoa,rangeEl_aoa);
        aodEst(1,:,l)=aodEi;
        aoaEst(1,:,l)=aoaEi;

        %% Verifica se executará o método simples ou adaptativo
        if it > 0
            Omega=zeros(size(A,1),size(A,2),k);
            [aodEstr,aoaEstr,toaEstr,its(l)]=fining(aodEi,aoaEi,rangeAz_aod,rangeEl_aod,rangeAz_aoa,rangeEl_aoa,x,Omega, Nt, Nr,R,Ts,c,it,lim,fc,arranjo);
            aodEst(:,:,l)=aodEstr;
            aoaEst(:,:,l)=aoaEstr;
            toaEst(2:end,l)=toaEstr';
        else
            its=0;
        end

        %% Cálculo do resíduo do sinal
        indices(l)=maxi;
        for k=1:K
            % 3. orthogonalize
            columns(:,l,k)=A(:,maxi,k);
            omega=A(:,maxi,k);
            psi(:,l,k)=omega;
            for counter2=1:l-1
                psi(:,l,k)=psi(:,l,k)-(psi(:,counter2,k)'*omega)*psi(:,counter2,k)/(norm(psi(:,counter2,k)))^2;
            end
            
            % Atualiza o resíduo
            lambda=0.1;
            beta=(psi(:,l,k)'*R(:,k)/(norm(psi(:,l,k)))^2+lambda);
            betamatrix(l,k)=beta;
            R(:,k)=R(:,k)-(beta)*psi(:,l,k);
        end
    end

    %% Determinação do ToA
    h=zeros(L,K);
    for k=1:K
        [Q,Rqr]=qr(columns(:,:,k),0);
        h(:,k)=inv(Rqr)*Q'*psi(:,:,k)*betamatrix(:,k);
    end
    for l=1:L
        toaEst(1,l)=estimateToA(h(l,:),Ts,c,K);
    end
end

%% Função para o DCS-SOMP Adaptativo
function [aodEstr,aoaEstr,toaEstr,its]=fining(aodEst,aoaEst,rangeAz_aod,rangeEl_aod,rangeAz_aoa,rangeEl_aoa,x,Omega, Nt, Nr,Y,Ts,c,it,lim,fc,arranjo)
    K=size(Omega,3);
    I=it;
    aodEstr(1,:)=aodEst;
    aoaEstr(1,:)=aoaEst;
    regis=false;
    its=I;
    for i=1:I
        %% Atualiza os candidatos
        difA=rangeAz_aod(2)-rangeAz_aod(1);

        ia1=aodEst(1)-difA;
        fa1=aodEst(1)+difA;
        if aodEst(1)==rangeAz_aod(1)
            ia1=aodEst(1);
        end
        if aodEst(1)==rangeAz_aod(end)
            fa1=aodEst(1);
        end

        ia2=aoaEst(1)-difA;
        fa2=aoaEst(1)+difA;
        if aoaEst(1)==rangeAz_aoa(1)
            ia2=aoaEst(1);
        end
        if aoaEst(1)==rangeAz_aoa(end)
            fa2=aoaEst(1);
        end

        difE=rangeEl_aod(2)-rangeEl_aod(1);
        ie1=aodEst(2)-difE;
        fe1=aodEst(2)+difE;
        if aodEst(2)==rangeEl_aod(1)
            ie1=aodEst(2);
        end
        if aodEst(2)==rangeEl_aod(end)
            fe1=aodEst(2);
        end

        ie2=aoaEst(2)-difE;
        fe2=aoaEst(2)+difE;
        if aoaEst(2)==rangeEl_aoa(1)
            ie2=aoaEst(2);
        end
        if aoaEst(2)==rangeEl_aoa(end)
            fe2=aoaEst(2);
        end


        if i==1
            L_az=6;
            L_el=6;
        else
            L_az=numel(rangeAz_aod);
            L_el=numel(rangeEl_aod);
        end

        %% Atualiza a matriz de detecção
        [Omega, rangeAz_aod, rangeEl_aod, rangeAz_aoa, rangeEl_aoa]=createDictionary(ia1,fa1,ie1,fe1,ia2,fa2,ie2,fe2,L_az,L_el,x,Omega, Nt, Nr,fc,arranjo);
        [ir,hr]=DCSSOMPRef(Y,Omega,1); %Chama o refinamento adaptativo, L=1

        %############### Estima os ângulos e o ToA novamente
        [aodEst,aoaEst]=estimateAngles(ir,rangeAz_aod,rangeEl_aod,rangeAz_aoa,rangeEl_aoa);
        toaEst=estimateToA(hr,Ts,c,K);
        %###############
        if abs(aodEst(1)-aodEstr(i,1)) < lim && regis==false
            its=i;
            regis=true;
        end
        aodEstr(i+1,:)=aodEst;
        aoaEstr(i+1,:)=aoaEst;
        toaEstr(i)=toaEst;
    end
end