% Boundary Conditions

function [pl, ql, pr, qr] = pdebc(xl,ul,xr,ur,t)
pl = [ul(1)-500; 0; 0; 0; ul(5)-2.5e4; ul(6)-2.5e4];
ql = [0; 1; 1; 1; 0; 0]; 
pr = [0; 0; 0; 0; 0; 0];
qr = [1; 1; 1; 1; 1; 1];
end
