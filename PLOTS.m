clear; clc; close all;

%% PLOTS
disp('Give plots of ISP and Thrust for given radius and Chamber psi: ');

%% Constants
mach_throat = 1;    % mach
pressure_sea = 101325;  % Pascal
R = 8.3144598;  % Joules/(mol*Kelvin)
gravity = 9.80665;  % m/s^2
%% Input
chamberPressureInput = input('Enter chamber pressure [500 or 550]: ');
if chamberPressureInput == 500
    T=readtable('CEA_Proccessed/CEAParameters(500).xlsx');
    pressure_chamber = 500 *  6894.75729; % psi to pascal

elseif chamberPressureInput == 550
    T=readtable('CEA_Proccessed/CEAParameters(550).xlsx');
    pressure_chamber = 550 *  6894.75729; % psi to pascal

end

pressure_design=pressure_chamber/50; %P_TOT/P_exit = 50;
%% Iterative variables
throatRadiusIterate =0.45:0.01:0.55; % IN INCHES
OFratioIterate = 2.5:0.05:3.5; % OF RATIO

ISPStoreSealevel = zeros(length(OFratioIterate),length(throatRadiusIterate)); % VARIABLE TO STORE ISP
ISPStoreDesign = ISPStoreSealevel;

thrustStoreSea = ISPStoreSealevel; % VARIABLE TO STORE THRUST
thrustStoreDesign = ISPStoreDesign;

cStarStore = ISPStoreSealevel;  % VARIABLE TO STORE C*
massFlowStore = ISPStoreSealevel;   % Variable to store mDot
exitRadiusStore = ISPStoreSealevel;

%% The loop
for j = 1:1:length(throatRadiusIterate)
    radius_throat = throatRadiusIterate(j);

    for i = 1:1:length(OFratioIterate)
        n = i;
        %% Table Reading:

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
            (sqrt(temperature_chamber*R_throat)))*((sqrt(gamma_throat*...
            (1+((gamma_throat-1)/2)*mach_throat^2)^...
            (-(gamma_throat+1)/(gamma_throat-1))))); % mass flow rate through an orifice at various mach numbers (choked is mach_throat = 1)
        massFlowStore(i,j) = massflow_throat;
        %% Calculate Exit Conditions
        mach_exit = sqrt((2/(gamma_exit-1))*((pressure_chamber/pressure_exit)^...
            ((gamma_exit-1)/gamma_exit)-1));
        %Solving for Mach Exit using StagnationPressure/Pressure Exit formula

        %temperature_exit = temperature_chamber/(1+((gamma_exit-1)/2)*mach_exit^2);
        temperature_exit = T{n, 16};
        epsilon = (1/mach_exit)*((2/(gamma_exit+1))*(1+(((gamma_exit-1)*mach_exit^2)/2)))^...
            ((gamma_exit+1)/(2*(gamma_exit-1)));
        %as a function of gamma_exit mach_exit, epsilon=A_exit/A_Throat by definition.

        area_exit = epsilon*area_throat;
        radius_exit = sqrt(area_exit/pi);
        exitRadiusStore(i,j) = radius_exit;
        velocity_exit = mach_exit*sqrt(gamma_exit*R_exit*temperature_exit);
        %epsilon = area_exit/area_throat;
        %% Calculate Thrust

        thrust_design = massflow_throat*velocity_exit+...
            (pressure_exit-pressure_design)*area_exit; % Assuming that design has perfect expansion at pressure altitude, otherwise change pressure to match
        thrust_sea = massflow_throat*velocity_exit+...
            (pressure_exit-pressure_sea)*area_exit;
        thrustStoreSea(i,j) = thrust_sea;
        thrustStoreDesign(i,j) = thrust_design;
        %% Performance Parameters
        c_star = pressure_chamber*area_throat/massflow_throat;
        c_star_verify = sqrt(R_chamber*temperature_chamber*(((gamma_chamber+1)/2)^((gamma_chamber+1)/(gamma_chamber-1)))*(1/gamma_chamber)); %combustion chamber performance parameter
        cStarStore(i,j) = c_star;
        thrust_coefficient = thrust_design/(pressure_chamber*area_throat);
        thrust_coefficient_verify = sqrt(((2*gamma_exit^2)/(gamma_exit-1))*((2/(gamma_exit+1))^((gamma_exit+1)/(gamma_exit-1)))*(1-((pressure_exit/pressure_chamber)^((gamma_exit-1)/gamma_exit))))+((pressure_exit-pressure_sea)/pressure_chamber)*(area_exit/area_throat); % Assuming that design has perfect expansion at pressure altitude, otherwise change pressure to match
        ISP_SeaLevel = thrust_sea/(massflow_throat*gravity);
        ISP_Design = velocity_exit/gravity;
        ISPStoreSealevel(i,j) = ISP_SeaLevel;
        ISPStoreDesign(i,j) = ISP_Design;
    end
end

%% Plots
close all;
% f = figure(1);
% f.Position =  [300 300 1000 563];
% hold on
% plot((OFratioIterate), ISPStoreDesign(:,length(throatRadiusIterate)));
% plot((OFratioIterate), ISPStoreSealevel(:,length(throatRadiusIterate)));
% xlabel('OF Ratio');
% ylabel('ISP at Design Altitude (ft)')
% title('ISP at Design Altitude vs OF Ratio [Various Throat Radius]')
% yyaxis right
% plot((OFratioIterate), thrustStoreDesign(:, length(throatRadiusIterate)));
% plot((OFratioIterate), thrustStoreSea(:, length(throatRadiusIterate)));
% legend({'Design Alt.', 'Sea Level', 'Design Alt.', 'Sea Level'}, 'Location','best');
% hold off
mtoinch=39.3701;
f = figure(1);
f.Position = [20 20, 1500 845];
%% Thrust Plot
subplot(2,2,1)
hold on
title('Sea Level Thrust vs OF Ratio [Various R_t]');
xlabel('OF Ratio')
ylabel('Thrust in Newtons')
for j = 1:1:(length(throatRadiusIterate))
    plot(OFratioIterate, thrustStoreDesign(:,j));
end
legend({'0.45', '0.46', '0.47', '0.48', '0.49', '0.50','0.51','0.52','0.53','0.54','0.55'}, ...
    'Location','best', 'NumColumns',2)
hold off
%% C* and ISP Plot
subplot(2,2,2)
hold on
title('Characteristic Velocity and ISP')
plot(OFratioIterate, cStarStore(:,j))
xlabel('OF Ratio')
ylabel('m/s')
yyaxis right
plot(OFratioIterate, ISPStoreSealevel(:, j));
plot(OFratioIterate, ISPStoreDesign(:, j));
ylabel('ISP (s)');
legend({'C*', 'ISP Sea', 'ISP Design'}, 'Location','best')
hold off
%% Mass Flow Rate Plot
subplot(2,2,3)
hold on
title('Mass Flow Rate vs OF Ratio [Various R_t]');
for j =1:1:length(throatRadiusIterate)
    plot(OFratioIterate , massFlowStore(:,j))
end
xlabel('OF Ratio')
ylabel('kg/s')
legend({'0.45', '0.46', '0.47', '0.48', '0.49', '0.50','0.51','0.52','0.53','0.54','0.55'}, ...
    'Location','best', 'NumColumns', 2)
hold off

subplot(2,2,4)
hold on
title('Exit Area vs OF [Various R_t]')
for j =1:1:length(throatRadiusIterate)
    plot(OFratioIterate, exitRadiusStore(:,j)*mtoinch)
end
xlabel('OF Ratio')
ylabel('Radius in Inches')
legend({'0.45', '0.46', '0.47', '0.48', '0.49', '0.50','0.51','0.52','0.53','0.54','0.55'}, ...
    'Location','best', 'NumColumns', 2)
hold off

for j = 1:1:2
    f = figure(j+1 );
    f.Position =  [20 20 1500 845];

    subplot(2,2,1);
    hold on
    titleText = sprintf('Thrust vs OF for (R_t = %0.2f)', throatRadiusIterate(j));
    title(titleText)
    plot(OFratioIterate, thrustStoreSea(:, j));
    plot(OFratioIterate, thrustStoreDesign(:,j));
    xlabel('OF Ratio');
    ylabel('Thrust in Newtons')
    yyaxis right
    plot((OFratioIterate), ISPStoreSealevel(:,1));
    plot((OFratioIterate), ISPStoreDesign(:, 1));
    ylabel('ISP in Seconds')
    legend({'Sea Level', 'Design Alt.','Sea Level', 'Design Alt.'}, ...
        'Location','best');

    hold off

    subplot(2,2,2);
    hold on
    titleText = sprintf('Mass Flow Rate vs OF Ratio (R_t = %0.2f)', throatRadiusIterate(j));
    title(titleText);
    plot(OFratioIterate,  massFlowStore(:, j));
    xlabel('Of Ratio')
    ylabel('Mass flow in kg/s')
    ylim([(massFlowStore(1, j) -0.05) (massFlowStore(length(OFratioIterate), j) +0.05)])
    legend({'Mass Flow'}, 'Location','best');
    hold off

    subplot(2,2,3);
    hold on
    titleText = sprintf('Exit Radius vs OF Ratio (R_t = %0.2f)', throatRadiusIterate(j));
    title(titleText)
    plot(OFratioIterate, exitRadiusStore(:, j) * mtoinch)
    xlabel('Of Ratio')
    ylabel('Exit Radius in INCHES')
    legend({'R_e'}, 'Location','best');
    hold off

    subplot(2,2,4);
    hold on
    titleText = sprintf('C* vs OF Ratio(R_t = %0.2f)', throatRadiusIterate(j));
    title(titleText)
    plot(OFratioIterate, cStarStore(:, j))
    xlabel('Of Ratio')
    ylabel('C* in m/s')
    legend({'C*'}, 'Location','best');
    hold off

end

