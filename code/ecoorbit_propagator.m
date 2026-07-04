% EcoOrbit — Physics Layer: Numerical Orbital Decay Propagator
% Forward-Euler integration of da/dt = -rho * v * a * Cd * A / m  (Eq. 4)
% Exports trajectory.csv consumed by the Python intelligence layer.
%
% Reference scenario: 250 kg LEO small satellite, 500 km circular orbit.
% Run:  >> ecoorbit_propagator

% ---- Constants (single source of truth; mirror of ecoorbit.py) ----------
mu    = 3.986e14;      % Earth gravitational parameter, m^3/s^2
Re    = 6378e3;        % Earth equatorial radius, m
rho0  = 6.97e-13;      % density at 500 km, kg/m^3
H     = 50e3;          % scale height, m
h_ref = 500e3;         % reference altitude for rho0, m
h_re  = 120e3;         % effective re-entry boundary, m
dt    = 3600;          % integration step, s (1 hour)

% ---- Spacecraft parameters ----------------------------------------------
m  = 250;              % mass, kg
Cd = 2.2;              % drag coefficient
A  = 3.2;              % reference cross-section, m^2  (beta = 35.5 kg/m^2)
h0 = 500e3;            % initial altitude, m

% ---- Algorithm 1: propagation -------------------------------------------
a = Re + h0;  v = sqrt(mu / a);  t = 0;
T = [];                % trajectory log: [t_days, h_km, v_kms, Fd_N]

while (a - Re) > h_re
    h   = a - Re;
    rho = rho0 * exp(-(h - h_ref) / H);        % Eq. (2)
    Fd  = 0.5 * Cd * A * rho * v^2;            % Eq. (1)
    da  = -rho * v * a * Cd * A / m * dt;      % Eq. (4)
    T(end+1, :) = [t/86400, h/1e3, v/1e3, Fd]; %#ok<SAGROW>
    a = a + da;
    v = sqrt(mu / a);
    t = t + dt;
end

fprintf('Numerical RUL: %.1f days (%.0f steps)\n', T(end,1), size(T,1));

% ---- King-Hele analytical cross-check ------------------------------------
beta = m / (Cd * A);
t_kh = beta * H / (rho0 * sqrt(mu/(Re+h0)) * (Re+h0)) ...
       * (1 - exp(-(h0 - h_re)/H)) / 86400;
fprintf('King-Hele analytical: %.1f days (delta = %+.1f%%)\n', ...
        t_kh, 100*(t_kh - T(end,1))/T(end,1));

% ---- Figures: altitude, velocity, drag ------------------------------------
figure('Name', 'EcoOrbit — MATLAB Physics Layer', 'Color', 'w');

subplot(2, 2, 1);
plot(T(:,1), T(:,2), 'Color', [0.06 0.46 0.43], 'LineWidth', 2); hold on;
plot([0 T(end,1)], [120 120], '--r');
text(10, 150, 'Re-entry (120 km)', 'Color', 'r');
xlabel('Time (days)'); ylabel('Altitude (km)');
title('Altitude vs Time'); grid on;

subplot(2, 2, 2);
plot(T(:,1), T(:,3), 'Color', [0.06 0.46 0.43], 'LineWidth', 2);
xlabel('Time (days)'); ylabel('Velocity (km/s)');
title('Orbital Velocity vs Time'); grid on;

subplot(2, 2, 3);
semilogy(T(:,1), T(:,4), 'Color', [0.06 0.46 0.43], 'LineWidth', 2);
xlabel('Time (days)'); ylabel('Drag force (N, log)');
title('Atmospheric Drag vs Time'); grid on;

subplot(2, 2, 4);
plot(T(:,1), T(:,2), 'Color', [0.06 0.46 0.43], 'LineWidth', 2); hold on;
plot([0 t_kh], [T(1,2) 120], '--', 'Color', [0.79 0.54 0.02], 'LineWidth', 1.5);
legend('Numerical (Euler)', sprintf('King-Hele endpoint: %.0f d', t_kh), ...
       'Location', 'southwest');
xlabel('Time (days)'); ylabel('Altitude (km)');
title(sprintf('Validation: \\Delta = %+.1f%%', 100*(t_kh - T(end,1))/T(end,1)));
grid on;

% Save a copy of the figure next to the CSV (works in MATLAB and Octave)
print('matlab_decay_plots.png', '-dpng', '-r120');

% ---- Export for Python intelligence layer --------------------------------
fid = fopen('trajectory.csv', 'w');
fprintf(fid, 't_days,altitude_km,velocity_km_s,drag_N\n');
fprintf(fid, '%.3f,%.3f,%.5f,%.6e\n', T');
fclose(fid);
disp('Exported trajectory.csv');
