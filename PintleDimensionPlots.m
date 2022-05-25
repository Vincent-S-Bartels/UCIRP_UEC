%% Pintle Dimension Plots
clear; clc; close all;
%Given Parameters (Change these)
%% Engine Constraints:
OFratio=3.05;%enter O/F wt ratio 2.50:.05:3.50
n=12;
radius_throat =.65;    % inches --later converted to metric
mach_throat = 1;    % mach
pressure_sea = 101325;  % Pascal
pressure_design=7.950e4;
R = 8.3144598;  % Joules/(mol*Kelvin)
gravity = 9.80665;  % m/s^2

%THERMO PARAMETERS BELOW ONLY VALID FOR 500PSIA CHAMBER AND 50=Pc/Pe
pressure_chamber = 3.792e6; % Pascal 550 psi

T=readtable('CEA_Proccessed/CEAParameters(550).xlsx');
temperature_chamber = T{n,14};  % Kelvin
pressure_exit = T{n,13}*10^5;   % Pascal
molarmass_chamber = T{n,8}; % g/mol
molarmass_throat = T{n,9};  % g/mol
molarmass_exit = T{n,10};   % g/mol
gamma_chamber = T{n,5};    % cp/cv
gamma_throat = T{n,6}; % cp/cv
gamma_exit = T{n,7};   % cp/cv
rho = T{n,2}; %
R_chamber = (R/molarmass_chamber)*1000; %Joules/(Kg*K)

% Throat Calculations
area_throat = pi*((radius_throat*(.0254))^2);
massflow_throat=((area_throat*pressure_chamber*mach_throat)/...
    (sqrt(temperature_chamber*R_chamber)))*((sqrt(gamma_throat*...
    (1+((gamma_throat-1)/2)*mach_throat^2)^...
    (-(gamma_throat+1)/(gamma_throat-1))))); % mass flow rate through an orifice at various mach numbers (choked is mach_throat = 1)
%% Preasure Iterate;
pressureLossCoef = 0.05; % Lose due to piping etc.
pressureDropCoef = 0.1:0.005:0.20;   %percent drop in pressure across injector
preInjectorPressure = (1./(1-pressureDropCoef)).*pressure_chamber; %PSI on propellant side of injector
tankPressure = (1./(1-pressureLossCoef)).*preInjectorPressure; %PSI required in tanks using given values and chamber pressure

delPInject = abs((preInjectorPressure - pressure_chamber));
Cd = 0.5;

%% Values from CEA
mdotTotal=massflow_throat;
mdotO2=mdotTotal/(1+(1/OFratio));
mdotCH4=mdotTotal-mdotO2;

rho_O2=1141;
rho_CH4=422.6;

T_CH4=111;
T_O2=90;
R_O2=259.84;

T_b_CH4=111.63;
T_cr_CH4=190.55;
P_cr_CH4=4.595e6;
T_r_CH4=T_CH4/T_cr_CH4;

mu_CH4=1225e-7;     %https://www.sciencedirect.com/science/article/pii/0031891473902577
mu_O2=6.93e-6;      %https://www.engineeringtoolbox.com/oxygen-O2-dynamic-kinematic-viscosity-temperature-pressure-d_2081.html
sig_CH4=12.909e-3;  %http://www.ddbst.com/en/EED/PCP/SFT_C1051.php

Qdot02 = mdotO2/rho_O2;
QdotCH4 = mdotCH4/rho_CH4;
QdotTotal = Qdot02+QdotCH4;
mtoinch=39.3701;


%D_pt = (0.25:0.05:0.5)*0.0254;
D_pt = .375*0.0254;
D_pr=.25*0.0254;
%theta_pt = 45;
theta_pt = 30:5:60;
v_CH4 = zeros(length(delPInject),1);
v_O2 = zeros(length(delPInject), 1);
t_ann = zeros(length(delPInject), 1);
for j = 1:1:length(theta_pt)
    for i = 1:1:length(delPInject)
        A_inletO2(i,j) = mdotO2/(Cd*sqrt(2*delPInject(i)*rho_O2));
        A_inletCH4(i,j)=mdotCH4/(Cd*sqrt(2*delPInject(i)*rho_CH4));
        v_CH4(i) = mdotCH4/(A_inletCH4(i,j)*rho_CH4);
        v_O2(i) = mdotO2/(A_inletO2(i,j)*rho_O2);

        D_cg=2*sqrt((A_inletCH4(i,j)/pi)+(D_pr/2)^2);
        t_post=(D_pt-D_cg)/2;      %Thickness of post between D_cg and D_post
        D_post=D_cg+2*t_post;
        D_outer(i)=2*sqrt((A_inletO2(i,j)./pi)+(D_post/2).^2);
        t_ann(i)=(D_outer(i)-D_post)/2;
        t_cg=(D_cg-D_pr)/2;

        R_pt=D_pt/2;
        R_pr=D_pr/2;
        R_cg=D_cg/2;
        R_post=D_post/2;
        R_outer=D_outer(i)/2;
        L_open(i,j) = (R_post-sqrt(R_post^2-A_inletCH4(i,j)*(sind(theta_pt(j))/pi)))/(sind(theta_pt(j))*cosd(theta_pt(j)))*mtoinch;
    end
end
%% Plots
close all;
subplot(2,2,1)
hold on
for j = 1:1:length(theta_pt)
    plot(pressureDropCoef, L_open(:,j));

end
xlabel("Pressure drop across injector (PSI)")
ylabel('Opening distance in inches')
legend({'30','35', '40','45','50','55','60'}, 'Location','best')
title('Opening distance vs Pressure with varying pintle angles')
hold off

subplot(2,2,2)
hold on
plot(pressureDropCoef, v_O2)
plot(pressureDropCoef, v_CH4)
legend({'O_2', 'CH_4'}, 'Location','best')
xlabel("Pressure drop across injector (PSI)")
ylabel('Velocity of Propellants (m/s)')
hold off

subplot(2,2,3)
hold on
plot(pressureDropCoef, t_ann*mtoinch);
xlabel("Pressure drop across injector (PSI)")
ylabel('Annular Distance in Inches')
hold off


subplot(2,2,4)
hold on
plot(pressureDropCoef, D_outer*mtoinch);
hold off

