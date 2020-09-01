%% INITIALIZE MATLAB
clear all;
clc;
close all;
format long;

%% DASHBOARD

% Initial parameters; defined by paper
Pe = 10579.5; % psia
P_AF = 8856.2;
Pw = 7998.9;
r_AF = 5.9; % ft
re = 626.1; % ft
Tr = 250; % degrees F
rw = 0.583; % ft
phi = .3;
k = .25; % darcy
PI = 3.1; % bbl/psi
Pb = 4000; % psi

% not defined by the paper; based on matching the initial pressure profile
mu = 4.5; % reservoir fluid viscosity, cp; 1.3
h = PI/7.08/k*mu*log(re/rw); % reservoir thickness, ft; 55

% TPR Parameters
gamma_o = 0.85;
rho_o = gamma_o*62.43; % lb/ft3
depth = 8000; % ft
P_wh = 300; % psi
D_tub = 3/12; % ft
e = 0.05;

Pwf = 4000:1000:9000;
for y = 1:length(Pwf);
   q_TPR(y) = calcFlowrate(depth, rho_o, P_wh, mu, D_tub, Pwf(y), e); 
end

% Increments
t_End = 8; % hrs
dt = .25; % hrs
dr = 0.05; % ft
M = round((r_AF - rw)/dr) + 1; % for damage radius
M1 = round((re - rw)/dr) + 1; % for drainage radius
N = t_End/dt +1;
t = 0:dt:t_End; % hrs

% flow-rate
q = 5000:1000:15500; % bbl/day

for m = 4

% Tuning parameters
alpha = 1/3;
beta = 13e-5;
gamma = beta/alpha;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Calculating Hydraulic Radius %% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dg = 6.5; % micron 
rH = dg/6*(phi/(1-phi));

% rH = 0.0314*(k*1000/phi)^.5;
dH = 4*rH; % micron
dAP = alpha*dH; % micron

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Calculating No-damage Pressure Profile %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


P = zeros(M1,N);
P(1,1) = 7998.9;
for i = 1:M1
    r(i) = dr*i;
    A_init(i) = 2*pi*r(i)*phi*h; % sqft
end

for i = 1:M1-1
    P(i+1,1) = P(1,1) + q(m)*mu/(7.08*k*h)*log(r(i+1)/rw);
end

% % plotting initial pressure profile vs. radius
figure
set(gcf, 'Position', [100, 100, 1100, 450])
subplot(1,2,1);
plot(r, P(:,1), '*r' , 'LineWidth', 2);
grid on;
title('Pressure Profile at t = 0');
xlabel('Radial Distance, ft');
ylabel('Pressure, psia');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Particle Size Distribution (acc. the paper) (Step 5) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSD data
AsphalteneFraction = [0.21 0.165 .125 .09 .07 .05 .04 .03 .015 .01 .03 .06 .046 .025 .015 0.009 0.007 0.002 0.001 0.0005];
ParticleDiameter = [0.00001 0.02 0.05 0.08 0.11 0.16 0.2 0.21 0.24 0.25 0.3 0.4 0.56 0.77 1 1.28 1.5 1.75 1.9 2];
% ref: Fig. 4 of the paper

% % PSD plot
subplot(1,2,2);
plot(ParticleDiameter, AsphalteneFraction, 'g', 'LineWidth', 2);
title('Particle Size Distribution');
xlabel('Particle Diameter, micron');
ylabel('Asphaltene Fraction');
grid on

% PSD integral
for o = 1:length(ParticleDiameter)
    if ParticleDiameter(o) >= dAP
        PD(o) = ParticleDiameter(o);
    end
end
PD1 = nonzeros(PD)';
nn = length(PD) - length(PD1) + 1;
f_trap = trapz(PD1, AsphalteneFraction(nn:length(AsphalteneFraction)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Calculating S (Step 4) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reservoir Fluid Molecular Weight
component_MoleFraction = [0.013887 0.533045 0.047187 0.039397 0.024263 0.036136 0.058189 0.062541 0.071256 0.021870 0.057985 0.007795 0.001417 0.019344 0.005688];
component_MolecularWeight = [42.117 16 30 44 74.839 119.657 133.514 170.318 257.588 265.684 330 453.624 466.771 380.000 475.000];
component_Density = [2 1 2 493 626 670 725 750 1200 1200 1400 1400 1400 561 1200]; % kg/m3, acc. Wikipedia

rho_RF = sum(component_Density.*component_MoleFraction); 
MW_RF = sum(component_MoleFraction.*component_MolecularWeight); % reservoir fluid mean molecular volume, gr/mol
MW_Asphaltene = 475; % gr/mole
% ref: Table 4-Reservoir Oil Characterization, SPE 37252
% also: Table-2 of the paper

wt_fraction = [.0001 .001 .001119 .0011115 .001101 0.0010845 .0001];
mole_fraction = wt_fraction.*MW_RF/MW_Asphaltene;
Pressure = [1000 3000 3500 4500 5000 5500 11000];
% ref: slide 23 of lecture 2

% Fit: 
[xData, yData] = prepareCurveData( Pressure, mole_fraction );

% Set up fittype and options.
ft = fittype( 'poly2' );

% Fit model to data.
S = fit( xData, yData, ft ); % weight fraction of asphaltene to reservoir fluid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Calculating M_RF using q (Step 4) and finally Calculating A_AP (Step 6) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculating molar volumes of asphaltene and reservoir fluid
vA = MW_Asphaltene/1.2; % cc/mole, asphaltene density equal to 1.2 g/cc
v_RF = MW_RF*1000/rho_RF; % cc/mole
A_AP = zeros(M,N);
% calculating the main parameters
    for j = 1:N
        M_RF(j) = q(m)*(dt/24)*j/v_RF*0.158987294928*1e6; % mole of reservoir fluid
    %     Main Loop
        for i = M:-1:1
            A_AP(i,j) = S(P(i,j))*f_trap*M_RF(j)*vA*gamma*6/dH*10.7639; % sqft    
            DOD(i,j) = 1/(1 - (A_AP(i,j)/A_init(i))); % (Step 7)
            k_dam(i,j) = k/DOD(i,j);                % (Step 8)
            phi_dam(i,j) = phi/DOD(i,j);      
            dP_dam(i,j) = q(m)*mu/1.127/k/A_init(i)*DOD(i,j)*dr;
            if j<N
                P(i,j+1) = P(i,j) - dP_dam(i,j);
                P_wf(m,j+1) = P(find(r==.75), j+1);
                dPs(:,j+1) = P(:,j) - P(:,j+1);
                skin(:,j+1) = dPs(:,j+1)*7.08*k*h/q(m)/mu;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Plotting the Results %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pressure profile vs. r
figure
set(gcf, 'Position', [100, 100, 700, 600])
subplot(2,2,1); 
plot(r(find(r==.5):M), P(find(r==.5):M,1:5:N) ,'LineWidth', 2)
xlabel('Radial Distance, ft');
ylabel('Pressure, psia');
legend('t=0', 't=1.25 hr', 't=2.5 hr', 't=3.75 hr', 't=5 hr', 't=6.25 hr', 't=7.5 hr' ); 
title('Post-Damage Pressure Profile');
grid on

% Permeability vs. r
subplot(2,2,2);
% figure
plot(r(find(r==0.5):M), k_dam(find(r==.5):M,4:5:N), 'LineWidth', 2)
xlabel('Radial Distance, ft');
ylabel('Permeability, Darcy');
legend( 't=.75 hr', 't=2 hr', 't=3.25 hr', 't=4.5 hr', 't=5.75 hr','t=7 hr' ); 
title('Post-Damage Permeability Profile');
grid on

% Porosity vs. r 
subplot(2,2,3);
% figure
plot(r(find(r==.5):M), phi_dam(find(r==.5):M,4:5:N), 'LineWidth', 2)
xlabel('Radial Distance, ft');
ylabel('Porosity, fraction');
legend( 't=.75 hr', 't=2 hr', 't=3.25 hr', 't=4.5 hr', 't=5.75 hr','t=7 hr' ); 
title('Post-Damage Porosity Profile');
grid on

% Skin vs. Time
% figure
subplot(2,2,4);
plot(t, skin(find(r==.5),:), '*' ,'LineWidth', 2)
xlabel('Time, hr');
ylabel('Skin');
title('Skin Factor vs. Time');
grid on
