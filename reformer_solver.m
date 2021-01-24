% https://se.mathworks.com/help/matlab/math/solve-system-of-pdes.html
x = linspace(0,0.45,45);
t = linspace(0,1,200);

m = 0;
sol = pdepe(m,@pdefun,@pdeic,@pdebc,x,t);

plot(t, sol(:,20,1))