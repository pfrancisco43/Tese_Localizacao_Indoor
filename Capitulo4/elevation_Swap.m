%Inverte o angulo de elevação
function out=elevation_Swap(in)
    out=in;
    out(:,2)=pi-in(:,2);
end

