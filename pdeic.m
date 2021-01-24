% The order is: u0 = [Temperature, C_CO2, C_H2, C_CO, C_H2O, C_MeOH]

function u0 = pdeic(x) 
u0=[500;1;1;1;2.5e4;2.5e4];
end