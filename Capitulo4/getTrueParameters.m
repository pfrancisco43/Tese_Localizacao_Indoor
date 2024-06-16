% Author: Paulo Francisco
% Objective: Calcule true parameters
% Syntax:
%       [ToAs, AoDs, AoAs]=getParameters(b,m,s,los)
% Inputs:
%       b, m, s - BS, MS, and SCs coordinates
%       los - if exists los
%       c - propagation speed
% Outputs:
%       ToAs, AoDs, AoAs - Parameters calculated

function [ToAs, AoDs, AoAs, L]=getTrueParameters(b,m,s,los,c)
    K=size(s,2);
    ToAs=zeros(K,1);
    AoDs=zeros(K,2);
    AoAs=zeros(K,2);
    i=1;
    if los==1
        ToAs=zeros(K+1,1);
        AoDs=zeros(K+1,2);
        AoAs=zeros(K+1,2);

        ToAs(i,1)=getToA (b,m);
        AoDs(i,:)=getAoD (b,m);
        AoAs(i,:)=getAoA (b,m);
        i=i+1;
    end
    if K>=1
        ToAs(i:end,1)=getToA (b,m,s);
        AoDs(i:end,:)=getAoD (b,s);
        AoAs(i:end,:)=getAoA (s,m);
    end
    ToAs=ToAs./c;
    L=numel(ToAs); 
end

function toa=getToA (b,m,s)
    if nargin == 2
        toa=norm(b-m);
    else
        toa=zeros(size(s,2),1);
        for i=1:size(s,2)
            toa(i,1)=norm(b-s(:,i))+norm(s(:,i)-m);
        end
    end
end

function aoas=getAoA(s,m)
    ns=size(s,2);
    aoas=zeros(ns,2);
    for i=1:ns
        aod_az=(pi/2) * (1 - sign(s(1,i) - m(1))) + atan((s(2,i) - m(2)) / (s(1,i) - m(1)));
        aod_el=(pi/2) - atan((s(3,i)-m(3))/sqrt((s(1,i)-m(1))^2 + (s(2,i)-m(2))^2));
        aoas(i,:)=[aod_az aod_el];
    end
end

function aods=getAoD (b,s)
    ns=size(s,2);
    aods=zeros(ns,2);
    for i=1:ns
        aod_az=(pi/2) * (1 - sign(s(1,i) - b(1))) + atan((s(2,i) - b(2)) / (s(1,i) - b(1)));
        aod_el=(pi/2) - atan((s(3,i)-b(3))/sqrt((s(1,i)-b(1))^2 + (s(2,i)-b(2))^2));
        aods(i,:)=[aod_az aod_el];
    end
end