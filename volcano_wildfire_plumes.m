%% ==================================================================
%        Volcanic/wildfire plume model comparisons
% ===================================================================
clear all; close all 

% Set up code path
codeDir = '~/code/research-projects/';
addpath(genpath(fullfile(codeDir,'hydroVolc/')))

% ========== Volcano setup ============

% -------------
R_d        = 287;   % gas constant of dry air (J/kg/K) 
R_v        = 461;   % gas constant of volcanic gas (water) (J/kg/K) 
epsilon     = R_d/R_v;
T0          = 273.15;

% ------ General plume input -----

% Might sweep over these a bit later, will want to invert for input...
Q = 1e5; % Start with a basic low-ish value
% T = [850 850 850 225] + T0;
% u0 = [55 75 55 10];
% r0 = [60 100 60 100];
% n0 = [0.08 0.05 0.08 0.99];
% lambda = [1e-2 1e-2 1e-5 1e-2];
% ll = {'Small volcanic','Large Volcanic','Small V., low cond.','Small wildfire'};

% 3 egs
% T = [850 850 225] + T0;
% u0 = [55 75 10];
% r0 = [60 100 500];
% n0 = [0.08 0.05 0.99];
% lambda = [1e-2 1e-2 1e-2];
% ll = {'Small volcanic','Moderate Volcanic','Large wildfire'};

% Simple eg.s
T = [850 150] + T0;
u0 = [75 10];
r0 = [100 300];
n0 = [0.05 0.99];
lambda = [1e-2 1e-2];
rho_m  = [2400 1200];
ll = {'Moderate Volcanic','Wildfire'};

% Key mass and heat flux
pI.T0   = 850 + 273.15;
pI.r_0  = 100;
pI.u0   = 35; 
pI.n_0  = 0.05;

% pI.atmo = fullfile(codeDir,'hydroVolc/atmoFiles/atmprofile.mat');
% pI.atmo = fullfile(codeDir,'hydroVolc/atmoFiles/atm_ERAreanalysis_Shinmoedake2011_absWind.mat');
pI.atmo = fullfile(codeDir,'hydroVolc/atmoFiles/atm_ERAreanalysis_Tungarahua2014_01_absWind.mat');
pI.D = 3;


% -------------- Wildfire setup ----------------




figDir = '~/Kahuna/phd-docs/research/wildfire/plumeFafo';
printfigs = false;
%%

q0 = zeros(length(u0));
for qi = 1:length(u0)
    pI.T0  = T(qi);
    pI.u_0 = u0(qi);
    pI.r_0 = r0(qi);
    pI.n_0 = n0(qi);
    pI.rho_m = rho_m(qi);
    pI.lambda = lambda(qi);
    
    dat(qi).pI = getPlumeSource(pI);
    
    dat(qi).pO = hmodel(dat(qi).pI,true);

    q0(qi) = dat(qi).pO.m(1);
    
    
%     plotAtmoProfile(pI.atmo)
    
    
    % First approx on u_0
%     rho_B0  = (volcSource.n_0/rho_v0 + (1-n_0)/rho_s)^(-1);
end

    % Temporary altitude range
%     z = dat(end).pO.z; 
    z = (0:100:2e4)';
    load(pI.atmo)
    atmo = interpAtmoArray(dat(1).pI.atmo,z);
    aP   = getAtmoProps(atmo,z);
    htropo  = findTPheight(atmo(:,2)/1e3,atmo(:,3)) * 1e3;
    Ptropo  = 10.^interp1(atmo(:,2),log10(atmo(:,1)),htropo,'pchip','extrap');

    % Plot plume(s)
co = [0.8500    0.3250    0.0980
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330];
    
plotPlumePairMS(dat,co(1:length(u0),:))
if printfigs; printpdf('Volc_fire_plumeOutput',figDir,[18 12]); end
skewt(aP.P/100,aP.T-273.15,aP.Tdew) % , varargin)
tp = plot(xlim,Ptropo/100*[1 1],'--k','LineWidth',1.3);
tp(2) = plot(nan,nan,'-k','LineWidth',1.3);
tp(3) = plot(nan,nan,'LineWidth',1.3,'Color',[0 0.4470 0.7410]);

for qi = 1:length(u0)
    Pq = 10.^interp1(atmo(:,2),log10(atmo(:,1)),dat(qi).pO.z,'pchip','extrap');
    Pb = 10.^interp1(atmo(:,2),log10(atmo(:,1)),dat(qi).pO.hb,'pchip','extrap');
    
    % Estimate dew point - NEEDS FIXING
    moistureProps = satProps(Pq,dat(qi).pO.theta,1.0); % assume to get initial props
    w_a = dat(qi).pO.m_v./dat(qi).pO.m_d;
    rh = w_a./moistureProps.w_s;
    % rh(rh>1) = 1;
    moistureProps = satProps(Pq,dat(qi).pO.theta,1.0);
    
%     qz = rh .* epsilon .* moistureProps.es./(Pq - moistureProps.es); % = w_a?
%     chi = log(P/100 .* qz ./ (6.112.*(epsilon + qz)));
%     props.T_dew = 243.5 .* props.chi ./ (17.67 - props.chi);

%     plot(Pq,dat(qi).pO.theta
    skewt(Pq/100,dat(qi).pO.theta-T0,dat(qi).pO.theta-T0,'DrawLines',0,'Dwptopts',{'LineWidth',1.7,'Color',co(qi,:)})
    plot(xlim,Pb/100*[1 1],'--','LineWidth',1.3,'Color',co(qi,:))
    
    % Temp dummies for legend
    pp(qi) = plot(nan,nan,'LineWidth',1.7,'Color',co(qi,:));
end

legend([tp pp],[{'Tropopause','Atmospheric T','Dew Pt. T'} ll])
set(gca,'Fontsize',14)
%% 
T = 400;
P0 = 101325;


es_0   = 100.*6.112.*exp(17.67.*(T-273.15)./(T +243.5-273.15)); % saturation vapor pressure for water
dp = max([P0-es_0  2*es_0],[],2); % Avoid a div0 when abs. pressure becomes very small
w_s    = 1./epsilon .* es_0 ./ dp; % w_a = r_h*w_s, mass mixing ratio of water vapor to dry air at saturation

%% Quick calcs for water vapour
T = 273.16:1:373.15;
P = 101325*ones(size(T));
rh = [0.2:0.2:1.0];
T0 = 273.15;


nr = 3;
figure
for ai = 1:nr
    ax(ai) = subplot(nr,1,ai);
end
hold(ax,'on')
co = get(gca,'ColorOrder');
for ri = 1:length(rh)

    V(ri) = satProps(P,T,rh(ri).*ones(size(T)));
    
    es_exp = exp(34.494 - (4924.99) ./ (V(ri).T+237.1-T0)) ./ (V(ri).T + 105-T0).^1.57;
    
    plot(ax(1),V(ri).T - T0, V(ri).es,'Color',co(ri,:));    
    plot(ax(1),V(ri).T - T0, es_exp,'--','Color',co(ri,:));    
    plot(ax(2),V(ri).T - T0, V(ri).w_s*1000);
    plot(ax(3),V(ri).T - T0, V(ri).T_dew);
    
end
ylabel(ax(1),'Sat. Vapour Pressure (Pa)')
ylabel(ax(2),'Vapour mixing ratio (g/kg)')
ylabel(ax(3),'Dew Point T (C)')
xlabel(ax(nr),'T (C)')
grid(ax,'on')
set(ax(1),'YScale','log')