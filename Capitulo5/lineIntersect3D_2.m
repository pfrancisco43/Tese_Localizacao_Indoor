%Detemina a menor distância entre várias retas no espaço
%Autor: Paulo Francisco
function [P_intersect] = lineIntersect3D_2(PA,PB)   
    Si = PB - PA; %N lines described as vectors
    L=size(PA,1);
    for i=1:L
        ss=sqrt(sum(Si(i,:).^2));
        for j=1:3
            Si(i,j)=Si(i,j)/ss;
        end
    end
    for i=1:3     
        ss=0;
        for j=1:3
            if i==j
                H(i,j)=sum(Si(:,j).^2-1);
                ss=ss+PA(:,j).*(Si(:,j).^2-1);
            else
                H(i,j)=sum(Si(:,i).*Si(:,j));
                ss=ss+PA(:,j).*(Si(:,i).*Si(:,j));
            end            
        end      
        c(i,1)=sum(ss);
    end    
    P_intersect = (H\c)';
end