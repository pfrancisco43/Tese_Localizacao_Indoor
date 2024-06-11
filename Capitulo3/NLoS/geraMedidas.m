function [m,mr]=geraMedidas(BS,UE,SP,covMedidas,B,alfa)

[~,QS]=size(SP);
m=zeros(QS,5);
mr=zeros(QS,5);
for i=1:QS
    %ToA
    SPi=SP(:,i);
    TOA=norm(SPi-UE)+norm(BS-SPi)+B;
    %AoD
    d_SP_BS=norm(BS-SPi);
    AOD=[atan2(SPi(2),SPi(1)) asin((SPi(3)-BS(3))./d_SP_BS)];
    %AoA
    AOA=[atan2(SPi(2)-UE(2),SPi(1)-UE(1))-alfa asin((SPi(3)-UE(3))/(norm(SPi-UE)))];    
    m(i,:)=[TOA AOD AOA];
    n=sqrtm(covMedidas)*randn(5,1);
    n=n';
    mr(i,:)=m(i,:)+n;
end