%Gera posições aleatórias para as BSs
function [SP]=geraBSs(NS,nd)
    %NS=NS-1;
    x=rand(1,NS)*100-50;    %eixo x
    %x=[0 x];
    y=rand(1,NS)*100-50;    %eixo y
    %y=[0 y];
    z=rand(1,NS)*40-20;     %eixo z
    %z=[40 z];
    
    SP=[x; y];    
    if(nd==3)
         SP=[x; y; z];
    end    
end
    
