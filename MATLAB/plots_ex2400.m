close all; clc; clearvars
% sim1 = sim('excGBS2', 15);
% sim1 = sim('ex2400',15);
 sim1 = sim('ex2400_elev'); %ELEVATOR
% sim1 = sim('ex2400_ail'); %AILERON
% sim1 = sim('ex2400_rud'); %RUDDER

%% Figure style settings
FONT_SIZE   = 18;
LABEL_SIZE  = 19;
LEGEND_SIZE = 16;
LINE_WIDTH  = 1.6;
AXES_WIDTH  = 1.1;
FIG_WIDTH_CM  = 16;
FIG_HEIGHT_PER_TILE_CM = 3.4;

set(0, 'DefaultAxesFontSize', FONT_SIZE);
set(0, 'DefaultTextFontSize', LABEL_SIZE);
set(0, 'DefaultLegendFontSize', LEGEND_SIZE);
set(0, 'DefaultLineLineWidth', LINE_WIDTH);
set(0, 'DefaultAxesLineWidth', AXES_WIDTH);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

GBS3_OUT = readmatrix('GBS_3Out.csv');
% MIRAR LAS COLUMNAS DE LAS DEFLEXIONES
% Columnas reales en GBS_3Out.csv
col.time = 1;
col.vt_fps = 22;
col.udot = 27;
col.vdot = 28;
col.wdot = 29;
col.p = 11;
col.q = 12;
col.r = 13;
col.pdot = 14;
col.qdot = 15;
col.rdot = 16;
col.alt_asl_ft = 62;
col.alt_agl_ft = 63;
col.phi_deg = 64;
col.theta_deg = 65;
col.psi_deg = 66;
col.elevator_cmd = 145;
col.pitch_trim = 146;
col.elevator_ctrl = 147;
col.aileron_ctrl = 150;
col.rudder_ctrl = 155;
col.ail_cmd_norm = 2;
col.left_ail_deg = 6;
col.right_ail_deg = 7;

t = GBS3_OUT(:, col.time);

%% CRUISE / POSICION
h_asl = GBS3_OUT(:, col.alt_asl_ft);
h_agl = GBS3_OUT(:, col.alt_agl_ft);
v = GBS3_OUT(:, col.vt_fps);
wrap180 = @(x) mod(x + 180, 360) - 180;
phi = wrap180(GBS3_OUT(:, col.phi_deg));
theta = wrap180(GBS3_OUT(:, col.theta_deg));
psi = wrap180(GBS3_OUT(:, col.psi_deg));

fig = new_shared_x_figure(3, FIG_WIDTH_CM, FIG_HEIGHT_PER_TILE_CM);

ax = nexttile;
plot(t, h_asl); grid on;
ylabel('$h$ (ft)');
xlim([0.0, 20]);
format_tile(ax, false);

ax = nexttile;
plot(t, v * 0.3048); grid on;
ylabel('$v$ (m/s)');
xlim([0.0, 20]);
format_tile(ax, false);

ax = nexttile;
plot(t, phi); grid on;
hold on
plot(t, theta);
plot(t, psi);
hold off
ylabel('(deg)');
legend('$\phi$', '$\theta$', '$\psi$', 'Location', 'best');
xlim([0.0, 20]);
format_tile(ax, true);

%% ROLL
ax1 = GBS3_OUT(:, col.udot);
p1 = deg2rad(GBS3_OUT(:, col.p));
d_a1 = rad2deg(GBS3_OUT(:, col.aileron_ctrl));
pdot1 = deg2rad(GBS3_OUT(:, col.pdot));

% delta_a_prop_deg = -0.50668; THIS IS TO ADD THE PROPELLER'S EFFECT
% d_a1 = d_a1 + delta_a_prop_deg;

fig = new_shared_x_figure(4, FIG_WIDTH_CM, FIG_HEIGHT_PER_TILE_CM);

ax = nexttile;
plot(t, d_a1); grid on;
ylabel('$\delta_a$ (deg)');
%ylim([-1, 1]);
xlim([0.0, 20]);
format_tile(ax, false);

ax = nexttile;
plot(t, p1); grid on;
ylabel('$p$ (rad/s)');
xlim([0.11, 20]);
format_tile(ax, false);

ax = nexttile;
plot(t, pdot1); grid on;
ylabel('$\dot{p}$ (rad/s$^2$)');
xlim([0.11, 20]);
format_tile(ax, false);

ax = nexttile;
plot(t, ax1 * 0.3048); grid on;
ylabel('$a_x$ (m/s$^2$)');
xlim([0.11, 20]);
format_tile(ax, true);

%% PITCH
az1 = GBS3_OUT(:, col.wdot);
q1 = deg2rad(GBS3_OUT(:, col.q));
d_e1 = rad2deg(GBS3_OUT(:, col.elevator_ctrl));
qdot1 = deg2rad(GBS3_OUT(:, col.qdot));

fig = new_shared_x_figure(4, FIG_WIDTH_CM, FIG_HEIGHT_PER_TILE_CM);

ax = nexttile;
plot(t, d_e1); grid on;
ylabel('$\delta_e$ (deg)');
%ylim([-1, 1]);
xlim([0.0, 20]);
format_tile(ax, false);

ax = nexttile;
plot(t, q1); grid on;
ylabel('$q$ (rad/s)');
xlim([0.11, 20]);
format_tile(ax, false);

ax = nexttile;
plot(t, qdot1); grid on;
ylabel('$\dot{q}$ (rad/s$^2$)');
xlim([0.11, 20]);
format_tile(ax, false);

ax = nexttile;
plot(t, az1 * 0.3048); grid on;
ylabel('$a_z$ (m/s$^2$)');
xlim([0.11, 20]);
format_tile(ax, true);

%% YAW
ay1 = GBS3_OUT(:, col.vdot);
r1 = deg2rad(GBS3_OUT(:, col.r));
d_r1 = rad2deg(GBS3_OUT(:, col.rudder_ctrl));
rdot1 = deg2rad(GBS3_OUT(:, col.rdot));

fig = new_shared_x_figure(4, FIG_WIDTH_CM, FIG_HEIGHT_PER_TILE_CM);

ax = nexttile;
plot(t, d_r1); grid on;
ylabel('$\delta_r$ (deg)');
% ylim([-1, 1]);
xlim([0.0, 20]);
format_tile(ax, false);

ax = nexttile;
plot(t, r1); grid on;
ylabel('$r$ (rad/s)');
xlim([0.11, 20]);
format_tile(ax, false);

ax = nexttile;
plot(t, rdot1); grid on;
ylabel('$\dot{r}$ (rad/s$^2$)');
xlim([0.11, 20]);
format_tile(ax, false);

ax = nexttile;
plot(t, ay1 * 0.3048); grid on;
ylabel('$a_y$ (m/s$^2$)');
xlim([0.11, 20]);
format_tile(ax, true);

%% Functions (plot's style)
function fig = new_shared_x_figure(n_tiles, width_cm, height_per_tile_cm)
% Makes the figure wide and short instead of long and narrow
fig = figure();
set(fig, 'Units', 'centimeters', ...
    'Position', [2, 2, width_cm, height_per_tile_cm * n_tiles]);
tiledlayout(n_tiles, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
end

function format_tile(ax, is_last)
% Applies the shared x-axis
if is_last
    xlabel(ax, '$t$ (s)');
else
    set(ax, 'XTickLabel', []);
end
end
