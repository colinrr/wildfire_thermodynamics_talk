function props = satProps(P,T,rh,iapws)
% TEMPORARY FUNCTION to get some reasonable properties
% P = Pressure (Pa)
% T = temp (K)
% rh = rel hum (fraction)
% iapws = use iapws - slower, but much more accurate

R_d        = 287;   % gas constant of dry air (J/kg/K) 
R_v        = 461;   % gas constant of volcanic gas (water) (J/kg/K) 
T0         = 273.15;
epsilon   = R_d/R_v;

props.P = P;
props.T = T;
props.RH = rh;

% props.Tsat_fofP = saturationTemperature(props.P);
props.es   = 100.*6.112.*exp(17.67.*(T-T0)./(T + 243.5 - T0)); % saturation vapor pressure for water
dp     = max([P-props.es  2*props.es],[],2); % Avoid a div0 when abs. pressure becomes very small
props.w_s    = 1./epsilon .* props.es ./ dp; % mass mixing ratio of water vapor to dry air at saturation
props.w_a   = rh.*props.w_s; % mass mixing ratio of water vapor to dry air

qz = rh .* epsilon .* props.es./(P - props.es); % = w_a?

props.chi = log(P/100 .* qz ./ (6.112.*(epsilon + qz)));

props.T_dew = 243.5 .* props.chi ./ (17.67 - props.chi);

props.rho_B  = P./(R_v * T) .* (1 + props.w_a) ./ (props.w_a + epsilon); % Bulk density


end