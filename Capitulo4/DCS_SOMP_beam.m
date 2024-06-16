function [ToAs_E, AoDs_E, AoAs_E,wo,yo]=DCS_SOMP_beam(Nt, Nr, Ns, N, B, c, y, w, L,Fc,arranjo,snr,alpha,yp,H,empilha)
    [ToAs_E, AoDs_E, AoAs_E,~]=DCS_SOMP_Adaptativo(Nt, Nr, Ns, N, B, c, y, w, L,Fc,arranjo,6,empilha);

    % Parametros atuais
    toas=ToAs_E;
    aods=AoDs_E;    
    aoas=AoAs_E;

    % toas=alpha(1);
    % aods=alpha(2:3);
    % aoas=alpha(4:5);

    %Otimiza w e y
    [wo,yo,ToAs_E,AoDs_E,AoAs_E]=lms (w, y, toas, aods, aoas, Nr, Nt, B, N, Fc, Ns, snr, arranjo,c,L,H,yp);

end

function [wo,yo,toas,aods,aoas] = lms(w, y, ToA, AoD, AoA, Nr, Nt, B, N, fc, Ns, snr, arranjo,c,L,Hp,yp)
    % Atualização dos Pesos do Feixe 

    % wo=exp(1j*rand(Nt,Ns,N)*2*pi);
    toas=ToA;
    aods=AoD;
    aoas=AoA;
    tol=1;
    
    aant=inf(1,L*5);
    
    yo=y;
    wo=w;

    mu = 0.0001;  % Taxa de aprendizagem
    %wo=zeros(Nt,Ns,N);
    beta=0.85; %Conhecimento do canal - Parâmetro ajustável
    nIt=1;%Parâmetro ajustável (quantidade de iterações)
    for i=1:nIt
        [Ha,~] = getH_y(wo,toas, aods, aoas, Nr, Nt, B, N,fc,Ns,snr,arranjo);
        Hf = beta * Hp + (1 - beta) * Ha; 
        aatu=[toas',reshape(aods,1,L*2),reshape(aoas,1,L*2)];
        erro=norm(aant-aatu,'fro');
        if erro<tol
            break;
        end
        aant=aatu;

        % Estrutura de laços para ajuste adaptativo
        ym=zeros(Nt,Ns,N);
        for n = 1:N  % Percorre cada subportadora
            for k = 1:Ns  % Percorre cada símbolo
                for j = 1:nIt  % Iterações para ajustar w para este símbolo e esta subportadora
                    % Cálculo do erro
                    e = yo(:, k, n) - Hf(:,:,n) * wo(:, k, n);

                    % Atualização de w usando LMS
                    wo(:, k, n) = wo(:, k, n) + mu * Hf(:,:,n)' * e;
                    ym(:, k, n)=Hf(:,:,n) * wo(:, k, n);

                    % Condição de parada (opcional)
                    v=norm(e);
                    if v < 1e-1                        
                        break;
                    end
                end
            end
        end

        yo=ym;
        [ToAs_E, AoDs_E, AoAs_E,~]=DCS_SOMP_Adaptativo(Nt, Nr, Ns, N, B, c, yo, wo, L,fc,arranjo,6,0);
        toas=ToAs_E;
        aods=AoDs_E;
        aoas=AoAs_E;
    end    
end