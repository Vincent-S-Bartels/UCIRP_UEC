%% Project Bulldog Engine Calculator
    %Created by Owen Trimble
    %Improved by Zackey Sahebzada
    %Revised by Brian Devine

clear; clc; close all;
%Given Parameters (Change these)

OFratio=input('Enter O/F wt. ratio (2.50:.05:3.50):');    %enter O/F wt ratio 2.50:.05:3.50
OFr=2.5:.05:3.5;
n=1;

while OFr(n)~=OFratio
  n=n+1;
  if n>numel(OFr)
    OFratio=input("Please Input Valid OF Ratio (2.50:.05:3.50): ");
    n=1;
  end
end

radius_throat =input('Enter throat RADIUS in inches:');    % inches --later converted to metric
mach_throat = 1;    % mach
pressure_sea = 101325;  % Pascal
pressure_design=7.950e4;
R = 8.3144598;  % Joules/(mol*Kelvin)
gravity = 9.80665;  % m/s^2

%THERMO PARAMETERS BELOW ONLY VALID FOR 500PSIA CHAMBER AND 50=Pc/Pe
pressure_chamber = 3.792e6; % Pascal 550 psi

T=readtable('CEAParameters(550).xlsx');

temperature_chamber = T{n,14};  % Kelvin
pressure_exit = T{n,13}*10^5;   % Pascal
molarmass_chamber = T{n,8}; % g/mol
molarmass_throat = T{n,9};  % g/mol
molarmass_exit = T{n,10};   % g/mol
gamma_chamber = T{n,5};    % cp/cv
gamma_throat = T{n,6}; % cp/cv
gamma_exit = T{n,7};   % cp/cv
rho = T{n,2}; %

%% Calculate the Gas Constants & Heat Capacity
R_chamber = (R/molarmass_chamber)*1000; %Joules/(Kg*K)
R_throat = (R/molarmass_throat)*1000; %Joules/(Kg*K)
R_exit = (R/molarmass_exit)*1000; %Joules/(Kg*K)

cp_chamber = (gamma_chamber*R_chamber)/(gamma_chamber-1); %Joules/(Kg*K)
cp_throat = (gamma_throat*R_throat)/(gamma_throat-1);%Joules/(Kg*K)
cp_exit = (gamma_exit*R_exit)/(gamma_exit-1);%Joules/(Kg*K)
cv_chamber = (R_chamber)/(gamma_chamber-1); %Joules/(Kg*K)
cv_throat = (R_throat)/(gamma_throat-1); %Joules/(Kg*K)
cv_exit = (R_exit)/(gamma_exit-1);
%% Throat Calculations
area_throat = pi*((radius_throat*(.0254))^2);
massflow_throat=((area_throat*pressure_chamber*mach_throat)/...
    (sqrt(temperature_chamber*R_chamber)))*((sqrt(gamma_throat*...
    (1+((gamma_throat-1)/2)*mach_throat^2)^...
    (-(gamma_throat+1)/(gamma_throat-1))))); % mass flow rate through an orifice at various mach numbers (choked is mach_throat = 1)

%% Calculate Exit Conditions
mach_exit = sqrt((2/(gamma_exit-1))*((pressure_chamber/pressure_exit)^...
    ((gamma_exit-1)/gamma_exit)-1));
%Solving for Mach Exit using StagnationPressure/Pressure Exit formula

temperature_exit = temperature_chamber/(1+((gamma_exit-1)/2)*mach_exit^2);

epsilon = (1/mach_exit)*((2/(gamma_exit+1))*(1+(((gamma_exit-1)*mach_exit^2)/2)))^...
    ((gamma_exit+1)/(2*(gamma_exit-1)));
%as a function of gamma_exit mach_exit, epsilon=A_exit/A_Throat by definition.

area_exit = epsilon*area_throat;
radius_exit = sqrt(area_exit/pi);
velocity_exit = mach_exit*sqrt(gamma_exit*R_exit*temperature_exit);
%epsilon = area_exit/area_throat;

%% Calculate Thrust

thrust_design = massflow_throat*velocity_exit+...
               (pressure_exit-pressure_design)*area_exit; % Assuming that design has perfect expansion at pressure altitude, otherwise change pressure to match
thrust_sea = massflow_throat*velocity_exit+...
            (pressure_exit-pressure_sea)*area_exit;

%% Performance Parameters
c_star = pressure_chamber*area_throat/massflow_throat;
c_star_verify = sqrt(R_chamber*temperature_chamber*(((gamma_chamber+1)/2)^((gamma_chamber+1)/(gamma_chamber-1)))*(1/gamma_chamber)); %combustion chamber performance parameter
thrust_coefficient = thrust_design/(pressure_chamber*area_throat);
thrust_coefficient_verify = sqrt(((2*gamma_exit^2)/(gamma_exit-1))*((2/(gamma_exit+1))^((gamma_exit+1)/(gamma_exit-1)))*(1-((pressure_exit/pressure_chamber)^((gamma_exit-1)/gamma_exit))))+((pressure_exit-pressure_sea)/pressure_chamber)*(area_exit/area_throat); % Assuming that design has perfect expansion at pressure altitude, otherwise change pressure to match
ISP = thrust_design/(massflow_throat*gravity);
ISP_verify = velocity_exit/gravity;

%% Injector Calc

%Created by Vincent S. Bartels

%% Chosen Design Constraints:

pressureLossCoef = 0.05; % Lose due to piping etc.
pressureDropCoef = 0.15;   %percent drop in pressure across injector
preInjectorPressure = (1/(1-pressureDropCoef))*pressure_chamber; %PSI on propellant side of injector
tankPressure = (1/(1-pressureLossCoef))*preInjectorPressure; %PSI required in tanks using given values and chamber pressure

delPInject = abs((preInjectorPressure - pressure_chamber));
Cd = 0.8;

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

%% Can 1/4" tubing handle the N2 Volumetric Flow Rate?

mdot_throat = (((0.152^2/4*pi)*600*1)/...
(sqrt(580*661.98)))*((sqrt(1.4*...
(1+((1.4-1)/2)*1^2)^...
(-(1.4+1)/(1.4-1)))));
volumetric = mdot_throat/(0.0725*35.3147);
chillin = volumetric>QdotTotal;

if chillin==1
  fprintf('We are Chillin\n\n');
else
  fprintf('We are NOT Chillin\n\n');
end

%% Injector Sizing

A_inletO2=mdotO2/(Cd*sqrt(2*delPInject*rho_O2));
A_inletCH4=mdotCH4/(Cd*sqrt(2*delPInject*rho_CH4));

v_CH4=mdotCH4/(A_inletCH4*rho_CH4);
v_O2=mdotO2/(A_inletO2*rho_O2);

D_pt=8e-3;        %Pintle Tip Diameter
D_pr=.375*D_pt;   %Diameter of Pintle Rod
theta_pt=input('Enter Pintle Angle (Degs): ');      %Pintle Angle (Deg)

D_cg=2*sqrt((A_inletCH4/pi)+(D_pr/2)^2);
t_post=(D_pt-D_cg)/2;      %Thickness of post between D_cg and D_post
D_post=D_cg+2*t_post;
D_outer=2*sqrt((A_inletO2/pi)+(D_post/2)^2);
t_ann=(D_outer-D_post)/2;
t_cg=(D_cg-D_pr)/2;

R_pt=D_pt/2;
R_pr=D_pr/2;
R_cg=D_cg/2;
R_post=D_post/2;
R_outer=D_outer/2;

%L_open=((v_O2*mdotCH4*rho_O2)/(v_CH4*mdotO2*rho_CH4))*((D_pt+t_ann)/D_pt)*t_ann;
L_open=(R_post-sqrt(R_post^2-A_inletCH4*(sind(theta_pt)/pi)))/(sind(theta_pt)*cosd(theta_pt));
mtoinch=39.3701;

fprintf('Pintle Tip Radius (in): %f\n',R_pt*mtoinch);
fprintf('Pintle Rod Radius (in): %f\n',R_pr*mtoinch);
fprintf('Pintle Center Gap Radius (in): %f\n',R_cg*mtoinch);
fprintf('Pintle Post Radius (in): %f\n',R_post*mtoinch);
fprintf('Pintle Outer Radius (in): %f\n',R_outer*mtoinch);
fprintf('Opening Distance of the Pintle (in): %f\n\n',L_open*mtoinch);

%% Injector Performance Calculations

K=(v_O2*t_ann)/(v_CH4*L_open);
TMR=(mdotCH4*v_CH4*cosd(theta_pt))/(mdotO2*v_O2+mdotCH4*v_CH4*sind(theta_pt));             %Total Momentum Ratio
alph=atand(TMR);  %Spray Angle of Propellants
We=(rho_O2*L_open*(v_O2-v_CH4)^2)/sig_CH4;                                                 %Weber Number

%% Chamber Geometry

xi=(90-theta_pt)/90;
s=1.15+1.35*xi;
p=1.3+.9*xi;
q=3.455-.225*xi;
SMD=L_open*xi^-1*exp(4-(q*We^.1));
Large_gamma=sqrt(((gamma_chamber+1)/gamma_chamber)^((gamma_chamber+1)/(gamma_chamber-1)));
u_droplet=(mdotCH4*v_CH4+mdotO2*v_O2)/mdotTotal;
