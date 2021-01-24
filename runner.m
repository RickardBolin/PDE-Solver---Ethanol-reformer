m = 10;
% Initial conditions
curr_states = zeros(m, 6);
curr_states(:, 1) = ones(m,1) * 500;
curr_states(:, 2) = ones(m,1) * 10;
curr_states(:, 3) = ones(m,1) * 10;
curr_states(:, 4) = ones(m,1) * 10;
curr_states(:, 5) = ones(m,1) * 2.6e4;
curr_states(:, 6) = ones(m,1) * 2.6e4;

% Boundary conditions
curr_states(1, 1) = 473.15;
curr_states(1, 2) = 1;
curr_states(1, 3) = 1;
curr_states(1, 4) = 1;
curr_states(1, 5) = 2.5e4;
curr_states(1, 6) = 2.5e4;

t0 = [0, 0.01];

[tsol, ysol] = ode23s(@(t, y) ethanolpde(t, y), t0, curr_states);
plot(linspace(0,1,10), ysol(end, 1:10))