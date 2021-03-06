% Equation to solve
function [c_pde,f,s] = pdefun(x,t,u,dudx)
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

R = 8.314e-3;

%----same value with the unit of either kPa*m^3/mol/K or KJ/mol/K-------%
K_eqD = @(T) 1.718e14.*exp(-95418e-3./(R*T));
K_eqW = @(T) 9.543e-3.*exp(39876e-3./(R*T));
K_eqR = @(T) 1.849e10.*exp(-56087e-3./(R*T));

%-------the unit of K_eq is kPa^x--------------%
%-------here the unit of R is KJ/mol/K-------% 
%-------the unit of T is always K---------%
P_H2O = @(c_H2O, T) c_H2O.*(R*T);
P_MeOH = @(c_MeOH,T) c_MeOH.*(R*T);
P_CO = @(c_CO, T) c_CO.*(R*T);
P_CO2 = @(c_CO2, T) c_CO2.*(R*T);
P_H2 = @(c_H2, T) c_H2.*(R*T);

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
r_H2 = @(c_CO2, c_H2, c_CO, c_H2O, c_MeOH, T) 3*r_CO2(c_CO2, c_H2, c_CO, c_H2O, c_MeOH, T)+2*r_CO(c_CO2, c_H2, c_MeOH, c_CO, c_H2O, T);
r_H2O = @(c_CO2, c_H2, c_CO, c_H2O, c_MeOH, T) -r_CO2(c_CO2, c_H2, c_CO, c_H2O, c_MeOH, T);
r_MeOH = @(c_CO2, c_H2, c_CO, c_H2O, c_MeOH, T) -(r_CO2(c_CO2, c_H2, c_CO, c_H2O, c_MeOH, T)+r_CO(c_CO2, c_H2, c_MeOH, c_CO, c_H2O, T));

% sum_ is the sum in the last term of the energy balance:
sum_ = @(c_CO2, c_H2, c_CO, c_H2O, c_MeOH, T) deltaH_rxnD * r_D(c_CO, c_H2, c_MeOH, T) + ...
deltaH_rxnW * r_W(c_CO2, c_H2, c_CO, c_H2O, T) + deltaH_rxnR * r_R(c_CO2, c_H2, c_MeOH, c_H2O, T);

c = sum(u(2:end));

% Eqs. 23 & 24 rewritten on the form c*dudt = df/dx + s*F(T,u), which is
% the form needed to use the pdepe-func.
c_pde = ones(6,1);

f = [-u(1)*((eps_bed*v_z*c*C_pgas)/(gamma*c*C_pgas + rho_cat*C_pcat)); 
     -u(2:end)*(eps_bed * v_z / gamma)]; 

s = [(1/(gamma*c*C_pgas + rho_cat*C_pcat))*((2*U/R_r)*(T_wall - u(1)) - eta_eff*rho_cat*sum_(u(2), u(3), u(4), u(5), u(6), u(1)));
    (eta_eff*rho_cat/gamma)*r_CO2(u(2), u(3), u(4), u(5), u(6), u(1));
    (eta_eff*rho_cat/gamma)*r_H2(u(2), u(3), u(4), u(5), u(6), u(1));
    (eta_eff*rho_cat/gamma)*r_CO(u(2), u(3), u(4), u(5), u(6), u(1));
    (eta_eff*rho_cat/gamma)*r_H2O(u(2), u(3), u(4), u(5), u(6), u(1));
    (eta_eff*rho_cat/gamma)*r_MeOH(u(2), u(3), u(4), u(5), u(6), u(1))];

% Not sure these work as I want them to work... I don't think they will
% let the concentrations increase again whenever they have become 0 once.
%f(u<0)=0;
%s(u<0)=0;

end
