close all; clc; clear all

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

%% RUDDER FRF
% Validate the rudder response using the JSBSim frequency response

sim1 = sim('ex2400_rud');

%% 1. Read JSBSim output
file_name = 'GBS_3Out.csv';
OUT = readmatrix(file_name);

time_col = 1;
rudder_col = 9;      % Rudder Position (deg)
r_col = 13;          % R (deg/s)
beta_col = 80;       % Beta (deg)
V_col = 22;          % V_Total (ft/s)
rho_col = 48;        % Rho (slugs/ft^3)
theta_col = 65;      % Theta (deg)
alpha_col = 79;      % Alpha (deg)

t = OUT(:, time_col);
dt = median(diff(t));
fs = 1 / dt;

dr = deg2rad(OUT(:, rudder_col));
r_meas = deg2rad(OUT(:, r_col));
beta_meas = deg2rad(OUT(:, beta_col));

%% 2. Remove trim offset
% The chirp XML excites from 2 to 22 s.
t_trim_ini = 1.0;
t_trim_fin = 2.0;
t_chirp_ini = 2.0;
t_chirp_fin = 22.0;

idx_trim = t >= t_trim_ini & t < t_trim_fin;
idx_chirp = t >= t_chirp_ini & t <= t_chirp_fin;

dr_trim = mean(dr(idx_trim));
r_trim = mean(r_meas(idx_trim));
beta_trim = mean(beta_meas(idx_trim));

dr1 = dr(idx_chirp) - dr_trim;
r1 = r_meas(idx_chirp) - r_trim;
beta1 = beta_meas(idx_chirp) - beta_trim;

fprintf('EX2400_TRIM rudder lateral validation\n')
fprintf('FRF interval: %.2f s to %.2f s, N = %d samples\n', t_chirp_ini, t_chirp_fin, nnz(idx_chirp))
fprintf('Rudder RMS: %.4f rad = %.3f deg\n\n', rms(dr1), rad2deg(rms(dr1)))

%% 3. JSBSim frequency response
frf_band = [0.1 5.0];
lateral_fit_band = [0.2 3.0];
f = logspace(log10(frf_band(1)), log10(frf_band(2)), 500);

win_length = min(round(5.0 * fs), nnz(idx_chirp));
win = hamming(win_length);
noverlap = round(0.80 * win_length);

Hr_jsb = tfestimate(dr1, r1, win, noverlap, f, fs);
Hbeta_jsb = tfestimate(dr1, beta1, win, noverlap, f, fs);

%% 4. Aircraft data from ex2400_trim.xml
g = 32.174;                         % ft/s^2
S = 5.68183775;                     % ft^2
b = 7.87401575;                     % ft
m = 5.5116 / g;                     % slugs

Ixx = 0.2529;                       % slug*ft^2
Izz = 0.47654398;                   % slug*ft^2
Ixz = 0.01629;                      % slug*ft^2

V0 = mean(OUT(idx_trim, V_col));
rho = mean(OUT(idx_trim, rho_col));
alpha0 = deg2rad(mean(OUT(idx_trim, alpha_col)));
theta0 = deg2rad(mean(OUT(idx_trim, theta_col)));

qbar0 = 0.5 * rho * V0^2;
qbarS = qbar0 * S;
qbarS_V = qbarS / V0;

u0 = V0 * cos(alpha0);
w0 = V0 * sin(alpha0);

% JSBSim applies aerodynamic forces at AERORP
empty_wt = 3.7479;
tank_wt = 1.7637;
cg_empty_in = [12.7953; 0.0000; -1.9685];
tank_in = [11.1000; 0.0000; -1.9685];
aerorp_in = [17.2875; 0.0000; -1.9685];
cg_total_in = (empty_wt*cg_empty_in + tank_wt*tank_in) / (empty_wt + tank_wt);

rx = (cg_total_in(1) - aerorp_in(1)) / 12.0;
rz = (cg_total_in(3) - aerorp_in(3)) / 12.0;

%% 5. Lateral-directional derivatives from ex2400_trim.xml
% Values in JSBSim body-axis convention.
CYb = -0.5316963;
CYp = 0.0;
CYr = 0.0;

Clb = -0.5440491;
Clp = -0.7119247;
Clr = 0.0582024;

Cnb = 0.3087214;
Cnp = -0.0982054;
Cnr = -0.0798044;

CYdr = -0.2169216;
Cldr = -0.0470489;
Cndr = 0.0834980;

%% 6. Analytical lateral model
Ixx_s = Ixx * cos(alpha0)^2 + Izz * sin(alpha0)^2 + Ixz * sin(2 * alpha0);
Izz_s = Ixx * sin(alpha0)^2 + Izz * cos(alpha0)^2 - Ixz * sin(2 * alpha0);
Ixz_s = Ixz * (cos(alpha0)^2 - sin(alpha0)^2) - 0.5 * (Ixx - Izz) * sin(2 * alpha0);

Yv = qbarS_V * CYb;
Yp = 0.5 * qbarS_V * b * CYp;
Yr = 0.5 * qbarS_V * b * CYr;

Lv = qbarS_V * b * Clb;
Lp = 0.5 * qbarS_V * b^2 * Clp;
Lr = 0.5 * qbarS_V * b^2 * Clr;

Nv = qbarS_V * b * Cnb;
Np = 0.5 * qbarS_V * b^2 * Cnp;
Nr = 0.5 * qbarS_V * b^2 * Cnr;

Lv = Lv - rz * Yv;
Lp = Lp - rz * Yp;
Lr = Lr - rz * Yr;

Nv = Nv + rx * Yv;
Np = Np + rx * Yp;
Nr = Nr + rx * Yr;

Ix_ = (Ixx_s * Izz_s - Ixz_s^2) / Izz_s;
Iz_ = (Ixx_s * Izz_s - Ixz_s^2) / Ixx_s;
Izx_ = Ixz_s / (Ixx_s * Izz_s - Ixz_s^2);

A_lat = zeros(4,4);
A_lat(1,1) = Yv / m;
A_lat(1,2) = Yp / m + w0;
A_lat(1,3) = Yr / m - u0;
A_lat(1,4) = g * cos(theta0);

A_lat(2,1) = Lv / Ix_ + Izx_ * Nv;
A_lat(2,2) = Lp / Ix_ + Izx_ * Np;
A_lat(2,3) = Lr / Ix_ + Izx_ * Nr;

A_lat(3,1) = Izx_ * Lv + Nv / Iz_;
A_lat(3,2) = Izx_ * Lp + Np / Iz_;
A_lat(3,3) = Izx_ * Lr + Nr / Iz_;

A_lat(4,2) = 1;
A_lat(4,3) = tan(theta0);

Ydr = qbar0 * S * CYdr;
Ldr = qbar0 * S * b * Cldr;
Ndr = qbar0 * S * b * Cndr;

Ldr = Ldr - rz * Ydr;
Ndr = Ndr + rx * Ydr;

B_lat = zeros(4,1);
B_lat(1) = Ydr / m;
B_lat(2) = Ldr / Ix_ + Izx_ * Ndr;
B_lat(3) = Izx_ * Ldr + Ndr / Iz_;

%% 7. Analytical transfer functions
Cr = [0 0 1 0];
Cbeta = [1/V0 0 0 0];

Hr_base = zeros(numel(f), 1);
Hbeta_base = zeros(numel(f), 1);

I4 = eye(4);
for k = 1:numel(f)
    s = 1i * 2*pi*f(k);
    x_over_dr = (s*I4 - A_lat) \ B_lat;

    Hr_base(k) = Cr * x_over_dr;
    Hbeta_base(k) = Cbeta * x_over_dr;
end

%% 8. Identify FRF-equivalent derivatives
% Each derivative is multiplied by a scaling factor. For every candidate, A and B matrices are rebuilt and compared with
% the two JSBSim frequency responses and the closest match with them is selected
derivative_names = {'CYb'; 'Clb'; 'Clp'; 'Clr'; 'Cnb'; ...
                    'Cnp'; 'Cnr'; 'CYdr'; 'Cldr'; 'Cndr'};
base_derivatives = [CYb; Clb; Clp; Clr; Cnb; ...
                    Cnp; Cnr; CYdr; Cldr; Cndr];

ratio_id = ones(10,1);
ratio_values = 0.90:0.005:1.10;
regularization = 0.20;

idx_fit = f >= lateral_fit_band(1) & f <= lateral_fit_band(2);
Hr_jsb_fit = Hr_jsb(idx_fit);
Hbeta_jsb_fit = Hbeta_jsb(idx_fit);
Hr_jsb_fit = Hr_jsb_fit(:);
Hbeta_jsb_fit = Hbeta_jsb_fit(:);

for sweep = 1:4
    for ratio_number = 1:10
        best_ratio = ratio_id(ratio_number);
        best_error = inf;

        for ratio_candidate = ratio_values
            ratio_test = ratio_id;
            ratio_test(ratio_number) = ratio_candidate;

            CYb_s = CYb * ratio_test(1);
            Clb_s = Clb * ratio_test(2);
            Clp_s = Clp * ratio_test(3);
            Clr_s = Clr * ratio_test(4);
            Cnb_s = Cnb * ratio_test(5);
            Cnp_s = Cnp * ratio_test(6);
            Cnr_s = Cnr * ratio_test(7);
            CYdr_s = CYdr * ratio_test(8);
            Cldr_s = Cldr * ratio_test(9);
            Cndr_s = Cndr * ratio_test(10);

            Yv_s = qbarS_V * CYb_s;
            Lv_s = qbarS_V * b * Clb_s;
            Lp_s = 0.5 * qbarS_V * b^2 * Clp_s;
            Lr_s = 0.5 * qbarS_V * b^2 * Clr_s;
            Nv_s = qbarS_V * b * Cnb_s;
            Np_s = 0.5 * qbarS_V * b^2 * Cnp_s;
            Nr_s = 0.5 * qbarS_V * b^2 * Cnr_s;

            Lv_s = Lv_s - rz * Yv_s;
            Nv_s = Nv_s + rx * Yv_s;

            A_test = zeros(4,4);
            A_test(1,1) = Yv_s / m;
            A_test(1,2) = w0;
            A_test(1,3) = -u0;
            A_test(1,4) = g * cos(theta0);
            A_test(2,1) = Lv_s / Ix_ + Izx_ * Nv_s;
            A_test(2,2) = Lp_s / Ix_ + Izx_ * Np_s;
            A_test(2,3) = Lr_s / Ix_ + Izx_ * Nr_s;
            A_test(3,1) = Izx_ * Lv_s + Nv_s / Iz_;
            A_test(3,2) = Izx_ * Lp_s + Np_s / Iz_;
            A_test(3,3) = Izx_ * Lr_s + Nr_s / Iz_;
            A_test(4,2) = 1;
            A_test(4,3) = tan(theta0);

            Ydr_s = qbar0 * S * CYdr_s;
            Ldr_s = qbar0 * S * b * Cldr_s - rz * Ydr_s;
            Ndr_s = qbar0 * S * b * Cndr_s + rx * Ydr_s;

            B_test = zeros(4,1);
            B_test(1) = Ydr_s / m;
            B_test(2) = Ldr_s / Ix_ + Izx_ * Ndr_s;
            B_test(3) = Izx_ * Ldr_s + Ndr_s / Iz_;

            Hr_test = zeros(numel(f),1);
            Hbeta_test = zeros(numel(f),1);

            for k = 1:numel(f)
                s = 1i * 2*pi*f(k);
                response = (s*I4 - A_test) \ B_test;
                Hr_test(k) = Cr * response;
                Hbeta_test(k) = Cbeta * response;
            end

            if any(~isfinite(Hr_test)) || any(~isfinite(Hbeta_test))
                continue
            end

            Hr_test_fit = Hr_test(idx_fit);
            Hbeta_test_fit = Hbeta_test(idx_fit);
            err_r_test = norm(Hr_test_fit(:) - Hr_jsb_fit) / ...
                max(norm(Hr_jsb_fit), eps);
            err_beta_test = norm(Hbeta_test_fit(:) - Hbeta_jsb_fit) / ...
                max(norm(Hbeta_jsb_fit), eps);

            total_error = err_r_test^2 + err_beta_test^2 + ...
                regularization * sum((ratio_test - 1).^2);

            if total_error < best_error
                best_error = total_error;
                best_ratio = ratio_candidate;
            end
        end

        ratio_id(ratio_number) = best_ratio;
    end
end

frf_equivalent_derivatives = base_derivatives .* ratio_id;

% Build the final FRF-equivalent model using the identified factors.
CYb_s = frf_equivalent_derivatives(1);
Clb_s = frf_equivalent_derivatives(2);
Clp_s = frf_equivalent_derivatives(3);
Clr_s = frf_equivalent_derivatives(4);
Cnb_s = frf_equivalent_derivatives(5);
Cnp_s = frf_equivalent_derivatives(6);
Cnr_s = frf_equivalent_derivatives(7);
CYdr_s = frf_equivalent_derivatives(8);
Cldr_s = frf_equivalent_derivatives(9);
Cndr_s = frf_equivalent_derivatives(10);

Yv_s = qbarS_V * CYb_s;
Lv_s = qbarS_V * b * Clb_s;
Lp_s = 0.5 * qbarS_V * b^2 * Clp_s;
Lr_s = 0.5 * qbarS_V * b^2 * Clr_s;
Nv_s = qbarS_V * b * Cnb_s;
Np_s = 0.5 * qbarS_V * b^2 * Cnp_s;
Nr_s = 0.5 * qbarS_V * b^2 * Cnr_s;

Lv_s = Lv_s - rz * Yv_s;
Nv_s = Nv_s + rx * Yv_s;

A_frf = zeros(4,4);
A_frf(1,1) = Yv_s / m;
A_frf(1,2) = w0;
A_frf(1,3) = -u0;
A_frf(1,4) = g * cos(theta0);
A_frf(2,1) = Lv_s / Ix_ + Izx_ * Nv_s;
A_frf(2,2) = Lp_s / Ix_ + Izx_ * Np_s;
A_frf(2,3) = Lr_s / Ix_ + Izx_ * Nr_s;
A_frf(3,1) = Izx_ * Lv_s + Nv_s / Iz_;
A_frf(3,2) = Izx_ * Lp_s + Np_s / Iz_;
A_frf(3,3) = Izx_ * Lr_s + Nr_s / Iz_;
A_frf(4,2) = 1;
A_frf(4,3) = tan(theta0);

Ydr_s = qbar0 * S * CYdr_s;
Ldr_s = qbar0 * S * b * Cldr_s - rz * Ydr_s;
Ndr_s = qbar0 * S * b * Cndr_s + rx * Ydr_s;

B_frf = zeros(4,1);
B_frf(1) = Ydr_s / m;
B_frf(2) = Ldr_s / Ix_ + Izx_ * Ndr_s;
B_frf(3) = Izx_ * Ldr_s + Ndr_s / Iz_;

Hr_frf = zeros(numel(f),1);
Hbeta_frf = zeros(numel(f),1);

for k = 1:numel(f)
    s = 1i * 2*pi*f(k);
    response = (s*I4 - A_frf) \ B_frf;
    Hr_frf(k) = Cr * response;
    Hbeta_frf(k) = Cbeta * response;
end

%% 9. Tables
derivative_table = table( ...
    derivative_names, ...
    base_derivatives, ...
    frf_equivalent_derivatives, ...
    abs((frf_equivalent_derivatives - base_derivatives) ./ base_derivatives) * 100, ...
    'VariableNames', {'Derivative','Base_VSPAERO_JSBSim', ...
    'FRF_Equivalent','Difference_pct'});

Hr_base_fit = Hr_base(idx_fit);
Hbeta_base_fit = Hbeta_base(idx_fit);
Hr_base_fit = Hr_base_fit(:);
Hbeta_base_fit = Hbeta_base_fit(:);

err_r = norm(Hr_base_fit - Hr_jsb_fit) / max(norm(Hr_jsb_fit), eps);
err_beta = norm(Hbeta_base_fit - Hbeta_jsb_fit) / max(norm(Hbeta_jsb_fit), eps);
Hr_frf_fit = Hr_frf(idx_fit);
Hbeta_frf_fit = Hbeta_frf(idx_fit);
err_r_frf = norm(Hr_frf_fit(:) - Hr_jsb_fit) / max(norm(Hr_jsb_fit), eps);
err_beta_frf = norm(Hbeta_frf_fit(:) - Hbeta_jsb_fit) / max(norm(Hbeta_jsb_fit), eps);

error_table = table( ...
    {'r/dr'; 'beta/dr'}, ...
    [err_r; err_beta], ...
    [err_r_frf; err_beta_frf], ...
    'VariableNames', {'FRF','Base_Model_Error','FRF_Equivalent_Error'});

disp('Base derivatives and FRF-equivalent derivatives:')
disp(derivative_table)

disp('Base and FRF-equivalent model errors:')
disp(error_table)

%% 10. Plots

FONT_SIZE   = 18;
LABEL_SIZE  = 19;
LEGEND_SIZE = 16;
LINE_WIDTH  = 1.6;
AXES_WIDTH  = 1.1;

set(0, 'DefaultAxesFontSize', FONT_SIZE);
set(0, 'DefaultTextFontSize', LABEL_SIZE);
set(0, 'DefaultLegendFontSize', LEGEND_SIZE);
set(0, 'DefaultAxesLineWidth', AXES_WIDTH);

figure('Units', 'centimeters', 'Position', [2, 2, 16, 12]);
tiledlayout(2,1,'TileSpacing','tight','Padding','compact')

ax1 = nexttile;
semilogx(f, 20*log10(abs(Hr_jsb)), 'k', ...
         f, 20*log10(abs(Hr_base)), '--', 'LineWidth', LINE_WIDTH)
grid on
set(ax1, 'TickLabelInterpreter', 'latex', 'XTickLabel', [])
ylabel('$|r/\delta_r|$ (dB)', 'Interpreter', 'latex')
% title('Rudder FRF validation', 'Interpreter', 'latex')
legend('JSBSim FRF', 'Analytical model', ...
    'Interpreter', 'latex', 'Location', 'best')

ax2 = nexttile;
semilogx(f, 20*log10(abs(Hbeta_jsb)), 'k', ...
         f, 20*log10(abs(Hbeta_base)), '--', 'LineWidth', LINE_WIDTH)
grid on
set(ax2, 'TickLabelInterpreter', 'latex')
ylabel('$|\beta/\delta_r|$ (dB)', 'Interpreter', 'latex')
xlabel('Frequency, $f$ (Hz)', 'Interpreter', 'latex')
legend('JSBSim FRF', 'Analytical model', ...
    'Interpreter', 'latex', 'Location', 'best')

linkaxes([ax1, ax2], 'x')
