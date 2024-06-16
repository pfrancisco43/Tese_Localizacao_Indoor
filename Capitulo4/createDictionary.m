%Cria o dicionário para cada iteração
%Autor: Paulo Francisco
function [Omega, rangeAz_aod, rangeEl_aod, rangeAz_aoa, rangeEl_aoa]=createDictionary(ia1,fa1,ie1,fe1,ia2,fa2,ie2,fe2,L_az,L_el,x,Omega,Nt, Nr,fc,arranjo)
    
    rangeAz_aod=linspace(ia1,fa1,L_az);
    rangeEl_aod=linspace(ie1,fe1,L_el);
    
    rangeAz_aoa=linspace(ia2,fa2,L_az);   
    rangeEl_aoa=linspace(ie2,fe2,L_el);

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
    %

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
    aa=size(Omega,1);    
    N=size(Omega,3);
    Omega=zeros(aa,(NrA*NrE)^2,N);
    for n=1:N        
        Omega(:,:,n)=kronMult((Ut'*x(:,:,n)).',Ur);
    end
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