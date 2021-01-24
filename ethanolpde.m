function derivs = ethanolpde(t, curr_states)
% derivs = derivates for the differents odes
% params = vector of parameters (gamma, c, C_p-gas, C_p-cat,epsilon_bed, v_z
% U, R_r, T_wall, n_eff, rho_cat, rho_eff,dH)
%args = vector of args [T_i-1, step_size, C_CO2_i-1, C_CO_i-1, C_H2O_i-1, C_H2_i-1, C_MeOH_i-1] 
% Currstates = [T, C]

[m, ~] = size(curr_states);
curr_states = reshape(curr_states,[m(1)/6,6]);
[m, n] = size(curr_states);
derivs = zeros(m,n);
step_size = 1 / m;

eps_bed = 0.37;
eps=0.528;
gamma=eps_bed + (1 - eps_bed)*eps;
v_z=1; 
%-------m/s---------%
eta_eff=0.1;
rho_cat=3960;
%----- kg/m^3----------%
%-----since the unit of r is mol/g/s, here rho_cat is supposed to have a unit of g/m^3---------%
%------the above statement is derived from equation (23) to make sure the
%units of both sides of = is equal.

c = sum(curr_states(2:end));
C_pgas=2500;
%-----J/kg/k--------%
C_pcat=1000;
%----J/kg/K---------%

U=50;
%-----W/m^2/K-------%
R_r=6.3e-2;
%----- m ---------%
T_wall=200+273.15;
%-------200~centigrade,noted by jingyu----------%
%-------473.15~K, noted by jingyu----------%

deltaH_rxnD=90.6;
deltaH_rxnW=-41.1;
deltaH_rxnR=49.5;
%----the unit of delta is kJ/mol--------%


% For the energy balance eq.:
%const1 = gamma*c*C_pgas + rho_cat*C_pcat;
%const2 = eps_bed*v_z*c*C_pgas;
%const3 = 2*U/R_r;
%const4 = eta_eff*rho_cat;

T = curr_states(:, 1);
C_CO2 = curr_states(:, 2);
C_CO = curr_states(:, 3);
C_H2O = curr_states(:, 4);
C_MeOH = curr_states(:, 5);
C_H2 = curr_states(:, 6);

R = 8.314e-3;

%----same value with the unit of either kPa*m^3/mol/K or KJ/mol/K-------%

K_eqD = @(T) 1.718e14.*exp(-95418e-3./(R*T));
K_eqW = @(T) 9.543e-3.*exp(39876e-3./(R*T));
K_eqR = @(T) 1.849e10.*exp(-56087e-3./(R*T));

%-------the unit of K_eq is kPa^x--------------%
%-------here the unit of R is KJ/mol/K-------% 
%-------the unit of T is always K---------%

P_MeOH = @(c_MeOH,T) c_MeOH.*(R*T);
P_CO = @(c_CO, T) c_CO.*(R*T);
P_CO2 = @(c_CO2, T) c_CO2.*(R*T);
P_H2 = @(c_H2, T) c_H2.*(R*T);
P_H2O = @(c_H2O, T) c_H2O.*(R*T);
%-----the unit of P is kPa---------%
%-----here the unit of R is kPa*m^3/mol/K-------%
%-----the unit of c is mol/m^3--------%

E_qD= @(c_CO, c_H2, c_MeOH, T) 1 - (P_CO(c_CO, T).*P_H2(c_H2, T).^2)./(K_eqD(T).*P_MeOH(c_MeOH,T));
E_qW= @(c_CO2, c_H2, c_CO, c_H2O, T) 1 - (P_CO2(c_CO2, T).*P_H2(c_H2, T))./(K_eqW(T).*P_CO(c_CO, T).*P_H2O(c_H2O, T));
E_qR= @(c_CO2, c_H2, c_MeOH, c_H2O, T) 1 - (P_CO2(c_H2, T).*P_H2(c_H2, T).^3)./(K_eqR(T).*P_MeOH(c_MeOH, T).*P_H2O(c_H2O, T));
%------E has no unit --------% 

k_D=@(T) 1.12*exp(-76/(R*T))';
k_W=@(T) 0.0023*exp(-50/(R*T))';
k_R=@(T) 6.75*exp(-81/(R*T))';
%-------the unit of k is same with Pre-exponential factor (mol/g/s kPa^x)--------%
%------here the unit of R is KJ/mol/K-------% 

r_D = @(c_CO, c_H2, c_MeOH, T) k_D(T).*P_MeOH(c_MeOH,T).*E_qD(c_CO, c_H2, c_MeOH, T);
r_W = @(c_CO2, c_H2, c_CO, c_H2O, T) k_W(T).*P_CO(c_CO, T).*E_qW(c_CO2, c_H2, c_CO, c_H2O, T);
r_R = @(c_CO2, c_H2, c_MeOH, c_H2O, T) k_R(T).*P_MeOH(c_MeOH,T).*E_qR(c_CO2, c_H2, c_MeOH, c_H2O, T);
%-----the unit of r is the unit of k products kPa (mol/g/s *kPa^(x+1))------%
%------based on analyze of shijie, here x is supposed to be -1 ---------%

r_CO2 = @(c_CO2, c_H2, c_CO, c_H2O, c_MeOH, T) r_W(c_CO2, c_H2, c_CO, c_H2O, T)+r_R(c_CO2, c_H2, c_MeOH, c_H2O, T);
r_CO = @(c_CO2, c_H2, c_CO, c_H2O, c_MeOH, T) r_D(c_CO, c_H2, c_MeOH, T)+r_W(c_CO2, c_H2, c_CO, c_H2O, T);
r_H2O = @(c_CO2, c_H2, c_CO, c_H2O, c_MeOH, T) -r_CO2(c_CO2, c_H2, c_CO, c_H2O, c_MeOH, T);
r_MeOH = @(c_CO2, c_H2, c_CO, c_H2O, c_MeOH, T) -(r_CO2(c_CO2, c_H2, c_CO, c_H2O, c_MeOH, T)+r_CO(c_CO2, c_H2, c_MeOH, c_CO, c_H2O, T));
r_H2= @(c_CO2, c_H2, c_CO, c_H2O, c_MeOH, T) 3*r_CO2(c_CO2, c_H2, c_CO, c_H2O, c_MeOH, T)+2*r_CO(c_CO2, c_H2, c_MeOH, c_CO, c_H2O, T);


sum_ = deltaH_rxnD * r_D(C_CO, C_H2, C_MeOH, T) + ...
deltaH_rxnW * r_W(C_CO2, C_H2, C_CO, C_H2O, T) + deltaH_rxnR * r_R(C_CO2, C_H2, C_MeOH, C_H2O, T);

derivs(2:m-1, 1) = -((eps_bed*v_z*c*C_pgas)/(gamma*c*C_pgas + rho_cat*C_pcat))*((T(3:end) - T(2:end-1))./step_size) + ...
    (1/(gamma*c*C_pgas + rho_cat*C_pcat))*((2*U/R_r)*(T_wall*ones(m-2, 1) - T(2:end-1)) - eta_eff*rho_cat*sum_(2:end-1));
%derivs(2:m-1, 1) = -(const2/const1) * ((T(3:end) - T(2:end-1))./step_size) + ...
%(const3/const1) * (T_wall*ones(m-2, 1) - T(2:end-1)) - ((const4/const1) * sum_(2:end-1));

% Don't include edges
C_CO2_ = C_CO2(2:end-1);
C_H2_ = C_H2(2:end-1);
C_CO_ = C_CO(2:end-1);
C_H2O_ = C_H2O(2:end-1);
C_MeOH_ = C_MeOH(2:end-1);


derivs(2:m-1, 2) = (eps_bed * v_z / gamma) * (curr_states(3:end, 2) ... 
- curr_states(2:end-1, 2))/step_size - (eta_eff*rho_cat/gamma) * r_CO(C_CO2_, C_H2_, C_CO_, C_H2O_, C_MeOH_, T(2:end-1));

derivs(2:m-1, 3) = (eps_bed * v_z / gamma) * (curr_states(3:end, 3) ... 
- curr_states(2:end-1, 3))/step_size - (eta_eff*rho_cat/gamma) * r_H2O(C_CO2_, C_H2_, C_CO_, C_H2O_, C_MeOH_, T(2:end-1));

derivs(2:m-1, 4) = (eps_bed * v_z / gamma) * (curr_states(3:end, 4) ... 
- curr_states(2:end-1, 4))/step_size - (eta_eff*rho_cat/gamma) * r_H2(C_CO2_, C_H2_, C_CO_, C_H2O_, C_MeOH_, T(2:end-1));

derivs(2:m-1, 5) = (eps_bed * v_z / gamma) * (curr_states(3:end, 5) ... 
- curr_states(2:end-1, 5))/step_size - (eta_eff*rho_cat/gamma) * r_MeOH(C_CO2_, C_H2_, C_CO_, C_H2O_, C_MeOH_, T(2:end-1));

derivs(2:m-1, 6) = (eps_bed * v_z / gamma) * (curr_states(3:end, 6) ... 
- curr_states(2:end-1, 6))/step_size - (eta_eff*rho_cat/gamma) * r_CO2(C_CO2_, C_H2_, C_CO_, C_H2O_, C_MeOH_, T(2:end-1));
derivs = reshape(derivs, [m(1)*6, 1]);