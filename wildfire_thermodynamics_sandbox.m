% Wildfire thermodynamics sandbox
clear all; close all
% -------- WILDFIRE params -----------

h_f         = 1.87e7; % J/kg heat of combustion
conv_frac   = 0.8; % Fraction of heat release delivered to convection
FMC         = 0.95; % Just fix a high-but-dry-ish value for now
l_v         = 2.257e6;  % latent heat of vaporization (J kg^-1)
cp_w        = 4190;  % specific heat of liquid water (J/kg/K)
dT_w        = 75;   % Assumed amount of heating in fuel water

%  --> one approach focusing on variations in spread rate and consumption
% r_burn = [100 500 1000 5000]; % Radius of the burned area perimeter. 
% frac_r_fireline = 1/3; % Assumption of how much of the perimeter is active at any one time
% 
n_fires = 1001;
w_vec = linspace(1,20,n_fires)'; % kg/m2 - forest fuel density
ros_vec = logspace(-1,0.6,n_fires)'; % m/s - fire spread rates, 0.1 - 4 m/s

% --> another approach focusing on fire growth
n_steps = 1001;

r_burn = logspace(1,3.7,n_steps)'; % m
frac_r_fireline = 1/3;

w = [1 5]; % kg/m2
ros = [0.1]; % m/s 
ros_hi = 1;
w_hi = 10;

% --------- VOLCANO PARAMS -----------

T_m = 850;  % T (C)
c_p = 1250; % J/kg/K

% MER's, kg/s
Q_m = [1e5 % Eyja low
       1e6; % Eyja peak
       1e7 % Taal 2020
       1e8]; % Grimsvotn 2011

% Assumed plume source params
rho_0 = 3; % kg/m3
u_0 = [40 80 150 300]';
   
figDir = '~/Kahuna/phd-docs/group-stuff/my-talks/2024-02-28/figures';
%% Calcs

% Volcanoes are easy...

H_v = Q_m .* c_p .* T_m;  
r_vent = (Q_m ./ (rho_0.*u_0 .* pi)).^(1/2);

% Wildfires are also kinda easy here...
h_f_adjusted = h_f - (l_v + dT_w.*cp_w);

L   = r_burn .* (2*pi) .* frac_r_fireline; % Fireline length
I_b = (1 - FMC./(1+FMC)) .* h_f .* w .* ros; % Byram's fire intensity

H_f = conv_frac .* I_b .* L; % Total sensible heat flux - worth accounting for sensible/latent heat loss to FMC as well?

H_f_extreme = conv_frac .* (1 - FMC./(1+FMC)) .* h_f .* w_hi .* ros_hi .* L;

I_b_vec = (1 - FMC./(1+FMC)) .* h_f .* [w w_hi] .* ros_vec;
Q_c = [w w_hi] .* ros_vec .* L; % "Mass combustion rate" kg/s

%% Plots

figure
pv = scatter(r_vent,H_v,'s','filled');
hold on
pa = plot(r_burn, H_f);
pb = plot(r_burn, H_f_extreme);
set(gca,'YScale','log','XScale','log')
xlabel('Plume source radius (m)')
ylabel('Sensible heat flux (W)')
grid on
legend([pv; pa; pb],{'Volcanoes','w = 1 kg/m^2, ROS = 0.1 m/s','w = 5 kg/m^2, ROS = 0.1 m/s','w = 10 kg/m^2, ROS = 1 m/s'})
printpdf('Sensible_heat_flux',figDir,[14 10])

dy = 0.1; 
figure
% tightSubplot(2,1,1,[],dy)
plot(nan,nan)
pa = plot(ros_vec,I_b_vec);
set(gca,'YScale','log','XScale','log')
ylabel('Byram''s Fire Intensity (W/m)')
xlabel('Rate of spread (m/s)')
grid on
hold on
pb = plot(ros_vec([1 end]), 4e6*[1 1], '--k');
legend([pa; pb],{'w = 1 kg/m^2','w = 5 kg/m^2','w = 10 kg/m^2','Suppression Limit'},...
    'location','northwest')
printpdf('I_b',figDir,[14 10])

% tightSubplot(2,1,2,[],dy)
% plot(ros_vec.*L, Q_c)
% ylabel('Mass combustion rate (kg/s)')
% xlabel('Area burn rate (m^2/s)')
% set(gca,'YScale','log','XScale','log')
% grid on
