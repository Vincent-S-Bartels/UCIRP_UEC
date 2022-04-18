%% Engine Calculator
%Created by Owen Trimble
%Improved by Zackey Sahebzada
%Revised by Brian Devine

clear; clc; close all;
%% Given Parameters (Change these)

OFratio=input('Enter O/F wt. ratio (2.50:.05:3.50):');    %enter O/F wt ratio 2.50:.05:3.50
radius_throat =input('Enter throat RADIUS in inches:');    % inches --later converted to metric 
mach_throat = 1;    % mach
pressure_sea = 101325;  % Pascal
R = 8.3144598;  % Joules/(mol*Kelvin)
gravity = 9.80665;  % m/s^2

%THERMO PARAMETERS BELOW ONLY VALID FOR 500PSIA CHAMBER AND 50=Pc/Pe
pressure_chamber = 3.447e6; % Pascal

T=readtable('CEA_Proccessed/CEAParameters(500).xlsx')

switch OFratio
    case 2.50		
		temperature_chamber = T{1,14};  % Kelvin
		pressure_exit = T{1,13}*10^5;   % Pascal
		molarmass_chamber = T{1,8}; % g/mol
		molarmass_throat = T{1,9};  % g/mol
		molarmass_exit = T{1,10};   % g/mol	
        gamma_chamber = T{1,5};    % cp/cv
        gamma_throat = T{1,6}; % cp/cv
        gamma_exit = T{1,7};   % cp/cv
        rho = T{1,2}; %
    case 2.55		
		temperature_chamber = T{2,14};  % Kelvin
		pressure_exit = T{2,13}*10^5;   % Pascal
		molarmass_chamber = T{2,8}; % g/mol
		molarmass_throat = T{2,9};  % g/mol
		molarmass_exit = T{2,10};   % g/mol
        gamma_chamber = T{2,5};    % cp/cv
        gamma_throat = T{2,6}; % cp/cv
        gamma_exit = T{2,7};   % cp/cv
    case 2.60		
		temperature_chamber = T{3,14};  % Kelvin
		pressure_exit = T{3,13}*10^5;   % Pascal
		molarmass_chamber = T{3,8}; % g/mol
		molarmass_throat = T{3,9};  % g/mol
		molarmass_exit = T{3,10};   % g/mol
        gamma_chamber = T{3,5};    % cp/cv
        gamma_throat = T{3,6}; % cp/cv
        gamma_exit = T{3,7};   % cp/cv
    case 2.65		
		temperature_chamber = T{4,14};  % Kelvin
		pressure_exit = T{4,13}*10^5;   % Pascal
		molarmass_chamber = T{4,8}; % g/mol
		molarmass_throat = T{4,9};  % g/mol
		molarmass_exit = T{4,10};   % g/mol
        gamma_chamber = T{4,5};    % cp/cv
        gamma_throat = T{4,6}; % cp/cv
        gamma_exit = T{4,7};   % cp/cv
    case 2.70		
		temperature_chamber = T{5,14};  % Kelvin
		pressure_exit = T{5,13}*10^5;   % Pascal
		molarmass_chamber = T{5,8}; % g/mol
		molarmass_throat = T{5,9};  % g/mol
		molarmass_exit = T{5,10};   % g/mol
        gamma_chamber = T{5,5};    % cp/cv
        gamma_throat = T{5,6}; % cp/cv
        gamma_exit = T{5,7};   % cp/cv
    case 2.75		
		temperature_chamber = T{6,14};  % Kelvin
		pressure_exit = T{6,13}*10^5;   % Pascal
		molarmass_chamber = T{6,8}; % g/mol
		molarmass_throat = T{6,9};  % g/mol
		molarmass_exit = T{6,10};   % g/mol 
        gamma_chamber = T{6,5};    % cp/cv
        gamma_throat = T{6,6}; % cp/cv
        gamma_exit = T{6,7};   % cp/cv
    case 2.80		
		temperature_chamber = T{7,14};  % Kelvin
		pressure_exit = T{7,13}*10^5;   % Pascal
		molarmass_chamber = T{7,8}; % g/mol
		molarmass_throat = T{7,9};  % g/mol
		molarmass_exit = T{7,10};   % g/mol
        gamma_chamber = T{7,5};    % cp/cv
        gamma_throat = T{7,6}; % cp/cv
        gamma_exit = T{7,7};   % cp/cv
    case 2.85		
		temperature_chamber = T{8,14};  % Kelvin
		pressure_exit = T{8,13}*10^5;   % Pascal
		molarmass_chamber = T{8,8}; % g/mol
		molarmass_throat = T{8,9};  % g/mol
		molarmass_exit = T{8,10};   % g/mol 
        gamma_chamber = T{8,5};    % cp/cv
        gamma_throat = T{8,6}; % cp/cv
        gamma_exit = T{8,7};   % cp/cv
    case 2.90		
		temperature_chamber = T{9,14};  % Kelvin
		pressure_exit = T{9,13}*10^5;   % Pascal
		molarmass_chamber = T{9,8}; % g/mol
		molarmass_throat = T{9,9};  % g/mol
		molarmass_exit = T{9,10};   % g/mol 
        gamma_chamber = T{9,5};    % cp/cv
        gamma_throat = T{9,6}; % cp/cv
        gamma_exit = T{9,7};   % cp/cv
    case 2.95		
		temperature_chamber = T{10,14};  % Kelvin
		pressure_exit = T{10,13}*10^5;   % Pascal
		molarmass_chamber = T{10,8}; % g/mol
		molarmass_throat = T{10,9};  % g/mol
		molarmass_exit = T{10,10};   % g/mol
        gamma_chamber = T{10,5};    % cp/cv
        gamma_throat = T{10,6}; % cp/cv
        gamma_exit = T{10,7};   % cp/cv
    case 3.00		
		temperature_chamber = T{11,14};  % Kelvin
		pressure_exit = T{11,13}*10^5;   % Pascal
		molarmass_chamber = T{11,8}; % g/mol
		molarmass_throat = T{11,9};  % g/mol
		molarmass_exit = T{11,10};   % g/mol 
        gamma_chamber = T{11,5};    % cp/cv
        gamma_throat = T{11,6}; % cp/cv
        gamma_exit = T{11,7};   % cp/cv
     case 3.05		
		temperature_chamber = T{12,14};  % Kelvin
		pressure_exit = T{12,13}*10^5;   % Pascal
		molarmass_chamber = T{12,8}; % g/mol
		molarmass_throat = T{12,9};  % g/mol
		molarmass_exit = T{12,10};   % g/mol 
        gamma_chamber = T{12,5};    % cp/cv
        gamma_throat = T{12,6}; % cp/cv
        gamma_exit = T{12,7};   % cp/cv
    case 3.10		
		temperature_chamber = T{13,14};  % Kelvin
		pressure_exit = T{13,13}*10^5;   % Pascal
		molarmass_chamber = T{13,8}; % g/mol
		molarmass_throat = T{13,9};  % g/mol
		molarmass_exit = T{13,10};   % g/mol
        gamma_chamber = T{13,5};    % cp/cv
        gamma_throat = T{13,6}; % cp/cv
        gamma_exit = T{13,7};   % cp/cv
    case 3.15		
		temperature_chamber = T{14,14};  % Kelvin
		pressure_exit = T{14,13}*10^5;   % Pascal
		molarmass_chamber = T{14,8}; % g/mol
		molarmass_throat = T{14,9};  % g/mol
		molarmass_exit = T{14,10};   % g/mol
        gamma_chamber = T{14,5};    % cp/cv
        gamma_throat = T{14,6}; % cp/cv
        gamma_exit = T{14,7};   % cp/cv
    case 3.20		
		temperature_chamber = T{15,14};  % Kelvin
		pressure_exit = T{15,13}*10^5;   % Pascal
		molarmass_chamber = T{15,8}; % g/mol
		molarmass_throat = T{15,9};  % g/mol
		molarmass_exit = T{15,10};   % g/mol
        gamma_chamber = T{15,5};    % cp/cv
        gamma_throat = T{15,6}; % cp/cv
        gamma_exit = T{15,7};   % cp/cv
    case 3.25		
		temperature_chamber = T{16,14};  % Kelvin
		pressure_exit = T{16,13}*10^5;   % Pascal
		molarmass_chamber = T{16,8}; % g/mol
		molarmass_throat = T{16,9};  % g/mol
		molarmass_exit = T{16,10};   % g/mol
        gamma_chamber = T{16,5};    % cp/cv
        gamma_throat = T{16,6}; % cp/cv
        gamma_exit = T{16,7};   % cp/cv
    case 3.30		
		temperature_chamber = T{17,14};  % Kelvin
		pressure_exit = T{17,13}*10^5;   % Pascal
		molarmass_chamber = T{17,8}; % g/mol
		molarmass_throat = T{17,9};  % g/mol
		molarmass_exit = T{17,10};   % g/mol
        gamma_chamber = T{17,5};    % cp/cv
        gamma_throat = T{17,6}; % cp/cv
        gamma_exit = T{17,7};   % cp/cv
    case 3.35		
		temperature_chamber = T{18,14};  % Kelvin
		pressure_exit = T{18,13}*10^5;   % Pascal
		molarmass_chamber = T{18,8}; % g/mol
		molarmass_throat = T{18,9};  % g/mol
		molarmass_exit = T{18,10};   % g/mol
        gamma_chamber = T{18,5};    % cp/cv
        gamma_throat = T{18,6}; % cp/cv
        gamma_exit = T{18,7};   % cp/cv
    case 3.40		
		temperature_chamber = T{19,14};  % Kelvin
		pressure_exit = T{19,13}*10^5;   % Pascal
		molarmass_chamber = T{19,8}; % g/mol
		molarmass_throat = T{19,9};  % g/mol
		molarmass_exit = T{19,10};   % g/mol
        gamma_chamber = T{19,5};    % cp/cv
        gamma_throat = T{19,6}; % cp/cv
        gamma_exit = T{19,7};   % cp/cv
    case 3.45		
		temperature_chamber = T{20,14};  % Kelvin
		pressure_exit = T{20,13}*10^5;   % Pascal
		molarmass_chamber = T{20,8}; % g/mol
		molarmass_throat = T{20,9};  % g/mol
		molarmass_exit = T{20,10};   % g/mol
        gamma_chamber = T{20,5};    % cp/cv
        gamma_throat = T{20,6}; % cp/cv
        gamma_exit = T{20,7};   % cp/cv
    case 3.50		
		temperature_chamber = T{21,14}; % Kelvin
		pressure_exit = T{21,13}*10^5;  % Pascal
		molarmass_chamber = T{21,8};    % g/mol
		molarmass_throat = T{21,9}; % g/mol
		molarmass_exit = T{21,10};  % g/mol
        gamma_chamber = T{21,5};    % cp/cv
        gamma_throat = T{21,6}; % cp/cv
        gamma_exit = T{21,7};   % cp/cv
end			
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
%% Calculate Throat Area
area_throat = pi*((radius_throat*(.0254))^2);
%% Calculate Mass Flow Rate
massflow_throat = ((area_throat*pressure_chamber*mach_throat)/...
    (sqrt(temperature_chamber*R_chamber)))*((sqrt(gamma_throat*...
    (1+((gamma_throat-1)/2)*mach_throat^2)^...
    (-(gamma_throat+1)/(gamma_throat-1))))) % mass flow rate through an orifice at various mach numbers (choked is mach_throat = 1)

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
    (pressure_exit-pressure_exit)*area_exit % Assuming that design has perfect expansion at pressure altitude, otherwise change pressure to match
thrust_sea = massflow_throat*velocity_exit+...
    (pressure_exit-pressure_sea)*area_exit;
%% Performance Parameters
c_star = pressure_chamber*area_throat/massflow_throat; % combustion chamber performance parameter
c_star_verify = sqrt(R_chamber*temperature_chamber/gamma_chamber*(((gamma_chamber+1)/2)...
    ^((gamma_throat+1)/(gamma_throat-1))));
thrust_coefficient = thrust_design/(pressure_chamber*area_throat);
thrust_coefficient_verify = sqrt(((2*gamma_exit^2)/(gamma_exit-1))*((2/(gamma_exit+1))^((gamma_exit+1)/(gamma_exit-1)))*(1-((pressure_exit/pressure_chamber)^((gamma_exit-1)/gamma_exit))))+((pressure_exit-pressure_exit)/pressure_chamber)*(area_exit/area_throat); % Assuming that design has perfect expansion at pressure altitude, otherwise change pressure to match
ISP = thrust_design/(massflow_throat*gravity)
ISP_verify = velocity_exit/gravity;



%% Injector Calc Addendum 
%Created by Vincent S. Bartels & Brian Devine
% Chosen Design Constraints:

rhoO2 = 1141; %kg/m3
rhoCH4 = 422.6;
Cd = 0.8;

% chosen delP (roughly 13% drop accross injector) tank pressure ~620 psi
% chamber pressure of 550 psi
delP = 70*6894.76 % psi to pascal, 


%% mass flow of prop/oxid
mdotTotal = massflow_throat;
mdot02 = mdotTotal / (1 + (1/OFratio));
mdotCH4 = mdot02 / OFratio;
Qdot02 = mdot02/rhoO2;
QdotCH4 = mdotCH4/rhoCH4;
QdotTotal = Qdot02 + QdotCH4;

%% Inlets:

areaO2 = mdot02/(Cd * sqrt(2*delP * rhoO2)); % m^2
vO2 = Cd*sqrt(2*delP/rhoO2); % m/s

areaCH4 = mdotCH4/(Cd*sqrt(2*delP*rhoCH4));
vCH4 = Cd*sqrt(2 * delP/rhoCH4);


%% Pintle specific addendum
%list of variables from (Middle East)
c_star;
Dcp = 0; % central gap diameter
Dip = 0; % inner body diameter
Dp = 0;  % pintle support thickness 
Dt = 0;  % pintle dip diamter
K1 = 0;  % dimensionless parameter
K2 = 0;  % dimensionless parameter
K3 = 0;  % dimensionless parameter
Laod = 0; % actual opening distance
Lopen = 0; % pintle opening (y axis from pintle slant edge to post
Ltod = 0; %tracing opening distance
Rcp = 0; %center radius post
tag = 0; %thickness of annual gap
tt = 0; %tip thickness 
thetaPintle = 45; %pintle angle in degress
thetaShadow = 0; %shadow angle
K = (rhoO2 * vO2^2) / (rhoCH4*vCH4^2) * (dO2/dCH4);



