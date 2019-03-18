%% Assignment 3 
%% Part 1
% For Part 1 , we have simulating the electrons in an electric field. This
% field will create a constant force on each electron and thus accelerate.


clc
clear

global C



C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665; %metres (32.1740 ft) per s²

%              ELECTRIC FIELD
% Voltage of 0.l V across x dimension
% Using V = E * d ,
%   We get E = 0.1 / 200 = 500 uV/m
eField = 0.1 / (200);

%                      FORCE  & ACCELERATION
% E = F / q , we know the q = C.q_0
force = eField * C.q_0;
eAcceleration = force / (C.m_0 * 0.26);





nSim = 150;
noe = 1000;
r2 = randi(360,noe,1);
xbound = 200;
ybound = 100;
x = randi(200,noe,1);
y = randi(100,noe,1);
vth = sqrt((C.kb * 300)/(C.m_0 * 0.26));
vx = vth * cos(r2) ;
vy = vth * sin(r2);

nPlot = 10;
colourArray= rand(nPlot,1);

MFP = vth * 0.2 * 10^-12;

pScat = 1 - exp((-35 * 10^-16)/(0.2 * 10^-12));

tMatrix = zeros (noe);


for t = 1:nSim
    vxc = vx; % create copy of vx
    vyc = vy; % create copy of vy
    [n,m] = size(vx);
    [n1,m1] = size(vy);
    
    %%randomly permutation of positions in vx and vy%%%
    idx = randperm(n);
    randomvx = vx;
    randomvx(idx,1)= vx (:,1) ;
    
    idy = randperm(n1);
    randomvy = vy;
    randomvy(idy,1) = vy(:,1);
    
    
    %Modelling scattering%%%%%%
    rScatter= rand(noe,1);
    
    % this gives 1s and 0s. 1 means it scatters
    tempScatter = rScatter < pScat;
    randomvx = tempScatter .* randomvx; % not scattered are 0s
    randomvy = tempScatter .* randomvy ; % not scattered are 0s
    
    %not scattered
    notScatter = rScatter >= pScat;
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    vx = vx .* notScatter; % the scattered vx are now 0
    vy = vy .* notScatter; % scattered vy = 0
    
    vx = vx + randomvx;
    vy = vy + randomvy;
    
    
    
    %%%%%%%%%%%%%%
    xc = x; % x copy
    yc = y; % y copy
    
    
    %Reflecting for y bounds%
    temp = y >= ybound ;
    temp1 = y < ybound ;
    
    
    temp = temp * -1;
    
    tempHigher = temp + temp1;
    
    
    temp2 = y <= 0;
    temp3 = y > 0;
    
    temp2 = temp2 * -1;
    tempLower = temp2 + temp3;
    
    vy = vy .* tempHigher;
    vy = vy .* tempLower;
    
    %%%%%%%%%%%%%%%%%%%
    
    % when x > 200%%%%%
    tempx1 = x <= 200;
    
    x = x .* tempx1;
    %%%%%%%%%%%%%%%%%%
    
    %%When x goes less than zero , come from 200 %%%%%
    
    tempx2 = x < -0.1;
    
    
    tempx2 = tempx2 * 200;
    tempxFinal = x + tempx2;
    
    x = tempxFinal;
    
    %%%%%%%%%%%%%%%%%%%
    dx = vx * (1/2000000);
    dy = vy * (1/2000000);
    
    x = x + dx;
    y = y + dy;
    vsq = (vy).^2 + (vx).^2 ;
    average = mean(vsq);
    
    for q = 1:1:nPlot
        plotx(q) = x(q);
        ploty(q) = y(q);
        
    end
    tMatrix = ((vsq * 0.26 *  C.m_0)/C.kb);

    
    
    
    figure(1)
    scatter(plotx,ploty,3,colourArray);
    axis([0 200 0 100]);
    xlabel("x");
    ylabel("y");
    title ("Electron simulation");
    vx =  vx + (eAcceleration * (1/20000));
    pause(0.01);
    hold on
    
    figure(2)
    I = average *noe * eField * C.q_0;
    scatter(t, I,'b.');
    axis tight;
    title("Current density");
    hold on;
    
end





[X,Y] = meshgrid (x , y);
f1 = scatteredInterpolant(x,y,tMatrix);
Z = f1(X,Y);
figure (3);
mesh(X,Y,Z);
title('Temperature plot');
xlabel('X Position');
ylabel('Y Position');
zlabel('Temparature(K)');

figure(4)
hist3([x y] , [50 50])
title("Electron Density Map")
