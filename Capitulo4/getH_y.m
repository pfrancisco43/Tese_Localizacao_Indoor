function [H,y] = getH_y(w,ToA, AoD, AoA, Nr, Nt, B, N,fc,Ns,snr,arranjo)
    % Parâmetros adicionais
    Ts = 1 / B; % Período de amostragem
    L = size(ToA, 1); % Quantidade de caminhos
        
    % Inicialização da matriz de resposta do canal
    H = zeros(Nr, Nt, N);    
    % Loop sobre cada subportadora
    h=10;   
    for n = 1:N
        % Loop sobre cada caminho
        for l = 1:L
            % Calcular vetor resposta para o caminho atual
            toa=ToA(l);
            aodaz=AoD(l, 1);
            aodel=AoD(l, 2);
            aoaaz=AoA(l, 1);
            aoael=AoA(l, 2);
            Ar = sqrt(Nr) * getResponse(fc, Nr, aoaaz, aoael,arranjo);
            At = sqrt(Nt) * getResponse(fc, Nt, aodaz, aodel,arranjo);
            % Calcular matriz de canal para o caminho atual e subportadora atual
            H(:,:,n) = H(:,:,n) + Ar * h* exp(-1j * 2 * pi * toa * (n-1) / (N * Ts)) * At';
        end
    end

    y=zeros(Nr,Ns,N);
           
    SNR_linear = 10^(snr/10);
    sigma=1/SNR_linear;
    for k=1:Ns
        for n=1:N 
            noise=sigma/sqrt(2)*(randn(Nr,1)+1j*randn(Nr,1));
            y(:,k,n)=H(:,:,n)*w(:,k,n);     
        end
    end
end