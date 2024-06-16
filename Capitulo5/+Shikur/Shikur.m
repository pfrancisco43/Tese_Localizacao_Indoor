%replicação do artigo DOA/AOD/AOA Localization in NLOS Environments.
%Autores: SHIKUR, B. Y.; WEBER, T. 

function [e, crlb, Xi]=Shikur(tdoas, aods, aoas, EB, EM, ES, sigmaToa, sigmaAngs)
    xb = EB(1);
    yb = EB(2);
    zb = EB(3);
    %valores obtidos pelo algoritmo SAGE (valores do artigo)
    ad = aods(:,1);%ψl
    ed = aods(:,2);%αl
    aa = aoas(:,1);% Φl
    ea = aoas(:,2);%βl
    dcl = tdoas;%δl
    
    %variáveis da formulação geométrica
    L=size(ES,2);
    al1=zeros(1,L);
    al2=zeros(1,L);
    al3=zeros(1,L);
    al4=zeros(1,L);
    al5=zeros(1,L);
    
       
    %calclando valores para cada caminho
    for i=1:L
        al1(i) = sin(ed(i)) * cos(ad(i)) + sin(ea(i)) * cos(aa(i));
        al2(i) = sin(ed(i)) * sin(ad(i)) + sin(ea(i)) * sin(aa(i));
        al3(i) = cos(ed(i)) + cos(ea(i));
        %al4(i) = sin(ea(i)) * sin(ed(i)) * sin(ad(i)) - aa(i);%Erro no
        %artigo, Desenvolvi novamente e ficou:
        al4(i)=sin(ea(i))*cos(aa(i))*al2(i) - sin(ea(i))*sin(aa(i))*al1(i);
        al5(i) = sin(ea(i)) * cos(aa(i)) * cos(ed(i)) - cos(ea(i)) * sin(ed(i)) * cos(ad(i));        
    end
    
    d1=0;
    %confere(EB,EM,ES,aods,aoas,tdoas,al1,al2,al3);
    %confere2(EM,EB,al1,al2,al3,al4,al5,dcl,d1);
    
    
    A=zeros(4,2*L);%a quantidade de expressões é 2x a quantidade de caminhos
    b=zeros(1,L*2);
    k=1;
    for i=1:2:L*2
        a1=[-al2(k) al1(k) 0 -al4(k)]';%são duas equações para cada caminho
        a2=[-al3(k) 0 al1(k) -al5(k)]';
        A(:,i)=a1;
        A(:,i+1)=a2;
        
        b1=-xb * al2(k) + yb * al1(k) + dcl(k) * al4(k);%são 2 valores de b para cada caminho
        b2=-xb * al3(k) + zb * al1(k) + dcl(k) * al5(k);
        b(1,i:i+1)=[b1 b2];
        k=k+1;
    end
    
    A=A';
    b=b';
    r=A\b;
    %confere2(r(1:3),EB,al1,al2,al3,al4,al5,dcl,r(4))
    
    %##################################################
    Xi=r(1:3);
    [~,crlb]=Shikur.estimadorAjusteFinalShikur(EB,EM,ES,tdoas, aods, aoas,sigmaToa,sigmaAngs,Xi);
   e=norm(Xi-EM);
    
    
end  
    %só pra conferir as captações
    %gera os valores pra verificar os parâmetros
function confere(EB,EM,ES,aods,aoas,tdoas,al1,al2,al3)
    ad = aods(:,1);%ψl
    ed = aods(:,2);%αl
    aa = aoas(:,1);% Φl
    ea = aoas(:,2);%βl
    dcl = tdoas;%δl
    s=ES;
    m=EM;
    b=EB;
    ns=size(ES,2);
    d1=0;
    for i=1:ns
        dbl=sqrt((s(1,i)-b(1))^2 + (s(2,i)-b(2))^2 + (s(3,i)-b(3))^2);
        dml=sqrt((s(1,i)-m(1))^2 + (s(2,i)-m(2))^2 + (s(3,i)-m(3))^2);
        dl=dbl+dml;
        S(1) = b(1) + dbl* sin(ed(i)) * cos(ad(i));
        S(2) = b(2) + dbl* sin(ed(i)) * sin(ad(i));
        S(3) = b(3) + dbl* cos(ed(i));
        M(1)=s(1,i)-(dl-dbl)*sin(ea(i))*cos(aa(i));
        MM(1)=b(1)+dbl*al1(i)-(d1+dcl(i))*sin(ea(i))*cos(aa(i));
        
        M(2)=s(2,i)-(dl-dbl)*sin(ea(i))*sin(aa(i));
        MM(2)=b(2)+dbl*al2(i)-(d1+dcl(i))*sin(ea(i))*sin(aa(i));
        
        M(3)=s(3,i)-(dl-dbl)*cos(ea(i));
        MM(3)=b(3)+dbl*al3(i)-(d1+dcl(i))*cos(ea(i));
        DBL=(m(1)-b(1)+(d1+dcl(i))*sin(ea(i))*cos(aa(i)))/al1(i);
    end
end

function confere2(EM,EB,al1,al2,al3,al4,al5,dcl,d1)
    m=EM;
    b=EB;
    ns=size(al1,2);    
    for i=1:ns
        aux1=-m(1)*al2(i) + m(2)*al1(i);
        aux2=-b(1)*al2(i) + b(2)*al1(i) + (d1+dcl(i))*al4(i);
        
        aux3=-m(1)*al3(i) + m(3)*al1(i);
        aux4=-b(1)*al3(i) + b(3)*al1(i) + (d1+dcl(i))*al5(i);
    end
end

    