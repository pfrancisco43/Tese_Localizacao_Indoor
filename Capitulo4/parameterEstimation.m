% Author: Paulo Francisco
% Objective: Estimation of the parameters
% Syntax:
%       [ToAs_E, AoDs_E, AoAs_E]=parameterEstimation (Nt, Nr, Ns, N, B, c, y, x, L, L_az, L_el, H, fin)
% Inputs:
%       Nt, Nr, Ns: Antennas Tx, Rx and number of symbols
%       N, B, c: Numbert of subcariers, BW and propagation speed
%       y, x, L: signal transmited, original symbols and number of paths
%       L_az, L_el, H, fin: Length of Grid sensing, Channel model and option of fine tunning
%
% Outputs:
%       ToAs_E, AoDs_E, AoAs_E - Parameters estimated
%
function [ToAs_E, AoDs_E, AoAs_E,its]=parameterEstimation (Nt, Nr, Ns, N, B, c, y, x, L, L_az, L_el, it, lim,fc,arranjo)
    Ts=1/B;  % Sampling period in us


    %Grid for the general scenario:
    % rangeAz_aod=linspace(10,170,L_az)*pi/180;
    % rangeEl_aod=linspace(91,170,L_el)*pi/180;
    % rangeAz_aoa=[linspace(10,170,L_az)]*pi/180;
    % rangeEl_aoa=linspace(91,170,L_el)*pi/180;
    
    %%Grid - In the specific simulation scenario (limit the grid for increase speed):
    rangeAz_aod=linspace(12,117,L_az)*pi/180;
    rangeEl_aod=linspace(93,115,L_el)*pi/180;
    rangeAz_aoa=linspace(6,147,L_az)*pi/180;
    rangeEl_aoa=linspace(96,120,L_el)*pi/180;

    %Grid - In the specific simulation scenario (limit the grid for increase speed):
    % rangeAz_aod=linspace(10,120,L_az)*pi/180;
    % rangeEl_aod=linspace(20,130,L_el)*pi/180;
    % rangeAz_aoa=linspace(100,210,L_az)*pi/180;
    % rangeEl_aoa=linspace(10,120,L_el)*pi/180;

    %% Create dictionary
    k=1;
    NrA=length(rangeAz_aoa);
    NrE=length(rangeEl_aoa);
    Ur=zeros(Nr,NrA*NrE);
    for i=1:NrA
        azi=rangeAz_aoa(i);
        for j=1:NrE
            ele=rangeEl_aoa(j);
            Ur(:,k)=sqrt(Nr) * getResponse(fc, Nr, azi, ele,arranjo);
            k=k+1;
        end
    end
    
    k=1;
    Ut=zeros(Nt,NrA*NrE);
    for i=1:NrA
        azi=rangeAz_aod(i);
        for j=1:NrE
            ele=rangeEl_aod(j);
            Ut(:,k)=sqrt(Nt) * getResponse(fc, Nt, azi, ele,arranjo);
            k=k+1;
        end
    end

    %% Vectorize and generation of the basis
    yb=zeros(Nr*Ns,N);
    Omega=zeros(Nr*Ns,(NrA*NrE)^2,N);     
    for n=1:N
        yb(:,n)=reshape(y(:,:,n),Nr*Ns,1);
        Omega(:,:,n)=kronMult((Ut'*x(:,:,n)).',Ur);
        %Omega(:,:,n)=kronMult(Ut,Ur);
        %Omega(:,:,n)=kronMult((Ut'*zeros(Nr,Ns)).',Ur);
    end

    %% run DCS-SOMP
    [AoDs_E,AoAs_E,ToAs_E,its]=DCSSOMP(yb,Omega,L,rangeAz_aod,rangeEl_aod,rangeAz_aoa,rangeEl_aoa,x,Nt, Nr, Ts, c,it, lim,fc,arranjo);
end

function X = kronMult(A,B)

    %KRON Kronecker product.
    %   kron(A,B) returns the Kronecker product of two matrices A and B, of
    %   dimensions I-by-J and K-by-L respectively. The result is an I*K-by-J*L
    %   block matrix in which the (i,j)-th block is defined as A(i,j)*B.
    %   Version: 06/02/2011
    %   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
    [I, J] = size(A);
    [K, L] = size(B);
    if ~issparse(A) && ~issparse(B)

        % Both matrices are dense.
        A = reshape(A,[1 I 1 J]);
        B = reshape(B,[K 1 L 1]);
        X = reshape(bsxfun(@times,A,B),[I*K J*L]);

    else

        % One of the matrices is sparse.
        [ia,ja,sa] = find(A);
        [ib,jb,sb] = find(B);
        ix = bsxfun(@plus,K*(ia(:)-1).',ib(:));
        jx = bsxfun(@plus,L*(ja(:)-1).',jb(:));

        % The @and operator is slightly faster for logicals.
        if islogical(sa) && islogical(sb)
            X = sparse(ix,jx,bsxfun(@and,sb(:),sa(:).'),I*K,J*L);
        else
            X = sparse(ix,jx,double(sb(:))*double(sa(:).'),I*K,J*L);
        end
    end
end