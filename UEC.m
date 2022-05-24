%% UEC
% The goal of this code is to describe the optimal performance of a potential
% engine with only the Oxidizer Fuel Ratio (OF by weight) and throat radius
clear all; clc; close all;

%OFratio=input('Enter O/F wt. ratio (2.50:.05:3.50):');    %enter O/F wt ratio 2.50:.05:3.50
OFratio = 3.05;
OFr=2.5:.05:3.5;
n=1;

while OFr(n)~=OFratio
    n=n+1;
    if n>numel(OFr)
        OFratio=input("Please Input Valid OF Ratio (2.50:.05:3.50): ");
        n=1;
    end
end

radius_throat = 0.5; %inches
%radius_throat =input('Enter throat RADIUS in inches:');    % inches --later converted to metric
mach_throat = 1;    % mach
pressure_sea = 101325;  % Pascal
pressure_design=7.950e4;
R = 8.3144598;  % Joules/(mol*Kelvin)
gravity = 9.80665;  % m/s^2

if input('Enter Chamber Pressure (500 or 550 PSI)\nEnter: ') == 500
    T=readtable('CEA_Proccessed/CEAParameters(500).xlsx');
    pressure_chamber = 500 *  6894.75729; % psi to pascal

else
    T=readtable('CEA_Proccessed/CEAParameters(550).xlsx');
    pressure_chamber = 550 *  6894.75729; % psi to pascal

end


temperature_chamber = T{n,14};  % Kelvin
pressure_exit = T{n,13}*10^5;   % Pascal
molarmass_chamber = T{n,8}; % g/mol
molarmass_throat = T{n,9};  % g/mol
molarmass_exit = T{n,10};   % g/mol
gamma_chamber = T{n,5};    % cp/cv
gamma_throat = T{n,6}; % cp/cv
gamma_exit = T{n,7};   % cp/cv
rho = T{n,2}; %