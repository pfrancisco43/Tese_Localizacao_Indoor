%Inverte ângulos AoD <-> AoA
%Autor: Paulo Francisco
function out=aod_aoa_Swap(in)
    out(:,1)=in(:,1)+pi;
    for i=1:size(out,1)
        if out(i,1)>1.5*pi
            out(i,1)=out(i,1)-2*pi;
        end
    end
    out(:,2)=pi-in(:,2);
end

