%% CEA Processing Script: Steps to follow
% 1. Input as a txt file the CEA tabulation 
% 2. Open and change the read file name to the new txt file in
% 'ingestData.py'
% 3. Run ingestData.py with chosen input and output file names
% 4. Modify the names within this matlab script to the appropriate file
% names and locations
% 5. Run the script
% 6. OPTIONAL: Pray.
%
% Created by Vincent S. Bartels and Zackey Sahebzada
format longG
clc;
clear;
%% Dimensions of input
x = readmatrix('CEA_Raw/CEA550.csv');
[i, j] = size(x);

%% Dimensions of output
k = i/3 - 1;
l = j*3 + 1;
y = zeros(k,l);
%% need to enter the proper OF Range
OF = linspace(2.5, 3.55, k+1)';
for m =1: length(OF)
    y(m,1) = OF(m);
end
xOffset = 0;

for a = 1:(i/3)
    for m = 1:3
        rho(m) = x(m + xOffset*3,1);
        gam(m) = x(m + xOffset*3,2);
        molMass(m) = x(m + xOffset*3,3);
        press(m) = x(m + xOffset*3,4);
        temp(m) = x(m + xOffset*3,5);
        ae(m) = x(m + xOffset*3, 6);
        mach(m) = x(m + xOffset*3, 7);
        sonV(m) = x(m + xOffset*3, 8);
    end

    for m = 1:3
        y(a,m+1) = rho(m);
        y(a,m+4) = gam(m);
        y(a,m+7) = molMass(m);
        y(a,m+10) = press(m);
        y(a, m+13) = temp(m);
        y(a, m+16) = ae(m);
        y(a, m+19) = mach(m);
        y(a, m+22) = sonV(m);
    end
    xOffset = xOffset + 1;
end


for m = 1:3:i
    chamber((m+2)/3,:) = x(m,:);
    throat((m+2)/3,:) = x(m+1,:);
    exit((m+2)/3,:) = x(m+2,:);
end

T = array2table(y);
% Default heading for the columns will be A1, A2 and so on. 
% You can assign the specific headings to your table in the following manner
T.Properties.VariableNames(1:l) = {'O/F', 'rho chamb', 'rho throat', 'rho exit' ...
    ,'gamma chamb', 'gamma throat', 'gamma exit', 'molecular weight M chamber' ...
    ,'molecular weight M throat', 'molecular weight M exit' ...
    ,'pressure chamber', 'pressure throat', 'pressure exit', ...
    'temp. chamber', 'temp. throat', 'temp. exit', 'area ratio chamber', ...
    'area ratio throat', 'area ratio exit', 'mach chamber', 'mach throat' ...
    ,'mach exit', 'sonV chamb', 'sonV throat', 'sonV exit'};
writetable(T,'CEA_Proccessed/CEA550Formatted.xlsx')

