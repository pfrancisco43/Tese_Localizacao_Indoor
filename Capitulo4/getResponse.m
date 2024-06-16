%Resposta, ou matriz de direção, com base no arranjo de antenas
%Autor: Paulo Francisco
function resp = getResponse(fc, N, azimute, elevacao,arranjo)
    % Velocidade da luz (m/s)
    c = 3e8;
    % Comprimento de onda (m)
    lambda = c / fc;

    %%Tipos de arranjo
    switch(arranjo)
        case 1
            resp = getURA(N, azimute, elevacao, lambda);
        case 2
            resp=getUCA(azimute,elevacao,N,fc);
        case 3
            resp=getURA2(N, azimute, elevacao, lambda);
        case 4
            resp = getSVURA(fc, lambda, N,azimute,elevacao);
    end
end

%ToolBox Matlab
function resp=getSVURA(fc, lam, N,az,el)    
    array = phased.URA('Size',[sqrt(N) sqrt(N)],'ElementSpacing',0.5*lam);
    steervec = phased.SteeringVector('SensorArray',array);    
    
    resp = steervec(fc,[rad2deg(az); rad2deg(el)]);
    
    resp=resp/sqrt(N);
end

%Artigo: Gain and Phase Autocalibration for Uniform Rectangular Arrays
function resp = getURA(N, az, el, lambda)
    M=sqrt(N);
    L=M;
    m=1:M;
    n=1:L;
    k=(2*pi)/lambda;
    d=lambda/2;
    
    Ux=sin(el)*cos(az);
    Uy=sin(el)*sin(az);

    % Ux=cos(el)*cos(az);
    % Uy=sin(el)*cos(az);

    ax=exp(1i*k*d*(m-1)*Ux);
    ay=exp(1i*k*d*(n-1)*Uy);

    resp=kron(ax,ay)';
    resp=(1/sqrt(N))*resp;

end

%Antenna Array Topologies for mmWave Massive MIMO Systems: Spectral Efficiency Analysis
function resp=getUCA(az,el,N,fc)    
    n=1:N-1;
    c=physconst('LightSpeed');
    lamb=c/fc;
    k=(2*pi)/lamb;
    dc=lamb/2;

    r=(N-1)*dc/(2*pi);

    % ang=sin(n);
    ang=((2*pi)/(N))*n;

    resp=exp(1i*k*r*sin(el)*cos(az-ang));
    resp=[1,resp]./N;
    resp=resp'./sqrt(N);
end

%Versão 2 do URA
function resp = getURA2(N, az, el, lambda)
    N_H = sqrt(N);
    N_V = sqrt(N);
    X_R = zeros(1,N);
    [Y_R,Z_R] = meshgrid(0:N_H-1,0:N_V-1);
    RPos = [X_R;Y_R(:).';Z_R(:).'];

    sv = [(sin(el).*cos(az)).';...
        (sin(el).*sin(az)).';...
        (cos(el)).'];

    k=(2*pi)/lambda;
    %k=pi;
    d=lambda/2;
    resp = (1/sqrt(N))*exp(1i*k*d*RPos.'*sv);
    % m=0:N_H-1;
    % n=0:N_V-1;
    % resp = exp(1i*k*d*(m*sin(el)*sin(az)+n*cos(el)))

end