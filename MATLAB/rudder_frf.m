close all; clc; clear all

%% HAMADA_AIL_EX2400_TRIM
% Validate the aileron response using the JSBSim frequency response

sim1 = sim('ex2400_ail');

%% 1. Read JSBSim output
file_name = 'GBS_3Out.csv';
OUT = readmatrix(file_name);

time_col = 1;
aileron_col = 6;     % Left Aileron Position (deg)
p_col = 11;          % P (deg/s)
r_col = 13;          % R (deg/s)
V_col = 22;          % V_Total (ft/s)
rho_col = 48;        % Rho (slugs/ft^3)
theta_col = 65;      % Theta (deg)
alpha_col = 79;      % Alpha (deg)

t = OUT(:, time_col);
dt = median(diff(t));
fs = 1/dt;

da = deg2rad(OUT(:, aileron_col));
p_meas = deg2rad(OUT(:, p_col));
r_meas = deg2rad(OUT(:, r_col));

%% 2. Remove trim offset
% The chirp XML excites from 2 to 22 s.
t_trim_ini = 1.0;
t_trim_fin = 2.0;
t_chirp_ini = 2.0;
t_chirp_fin = 22.0;

idx_trim = t >= t_trim_ini & t < t_trim_fin;
idx_chirp = t >= t_chirp_ini & t <= t_chirp_fin;

da_trim = mean(da(idx_trim));
p_trim = mean(p_meas(idx_trim));
r_trim = mean(r_meas(idx_trim));

da1 = da(idx_chirp) - da_trim;
p1 = p_meas(idx_chirp) - p_trim;
r1 = r_meas(idx_chirp) - r_trim;

fprintf('EX2400_TRIM aileron lateral validation\n')
fprintf('FRF interval: %.2f s to %.2f s, N = %d samples\n', t_chirp_ini, t_chirp_fin, nnz(idx_chirp))
fprintf('Aileron RMS: %.4f rad = %.3f deg\n\n', rms(da1), rad2deg(rms(da1)))

%% 3. JSBSim frequency response
frf_band = [0.1 5.0];
lateral_fit_band = [0.2 3.0];
f = logspace(log10(frf_band(1)), log10(frf_band(2)), 500);

win_length = min(round(5.0 * fs), nnz(idx_chirp));
win = hamming(win_length);
noverlap = round(0.80 * win_length);

Hp_jsb = tfestimate(da1, p1, win, noverlap, f, fs);

%% 4. Aircraft data from ex2400_trim.xml
g = 32.174;                         % ft/s^2
S = 5.68183775;                     % ft^2
b = 7.87401575;                     % ft
m = 5.5116/g;                     % slugs

Ixx = 0.2529;                       % slug*ft^2
Izz = 0.47654398;                   % slug*ft^2
Ixz = 0.01629;                      % slug*ft^2

V0 = mean(OUT(idx_trim, V_col));
rho = mean(OUT(idx_trim, rho_col));
alpha0 = deg2rad(mean(OUT(idx_trim, alpha_col)));
theta0 = deg2rad(mean(OUT(idx_trim, theta_col)));

phidot1 = p1 + r1 * tan(theta0);
Hphidot_jsb = tfestimate(da1, phidot1, win, noverlap, f, fs);
Hphi_jsb = Hphidot_jsb(:) ./ (1i * 2*pi*f(:));

qbar0 = 0.5 * rho * V0^2;
qbarS = qbar0 * S;
qbarS_V = qbarS/V0;

u0 = V0 * cos(alpha0);
w0 = V0 * sin(alpha0);

% JSBSim applies aerodynamic forces at AERORP
empty_wt = 3.7479;
tank_wt = 1.7637;
cg_empty_in = [12.7953; 0.0000; -1.9685];
tank_in = [11.1000; 0.0000; -1.9685];
aerorp_in = [17.2875; 0.0000; -1.9685];
cg_total_in = (empty_wt*cg_empty_in + tank_wt*tank_in)/(empty_wt + tank_wt);

rx = (cg_total_in(1) - aerorp_in(1))/12.0;
rz = (cg_total_in(3) - aerorp_in(3))/12.0;

%% 5. Lateral-directional derivatives from ex2400_trim.xml
%Values in JSBSim body-axis convention.

CYb = -0.5316963;
CYp = 0.0;
CYr = 0.0;

Clb = -0.5440491;
Clp = -0.7119247;
Clr = 0.0582024;

Cnb = 0.3087214;
Cnp = -0.0982054;
Cnr = -0.0798044;

CYda = -0.0174547;
ClDa = -0.3650137;
Cnda = -0.0011450;

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

Ix_ = (Ixx_s * Izz_s - Ixz_s^2)/Izz_s;
Iz_ = (Ixx_s * Izz_s - Ixz_s^2)/Ixx_s;
Izx_ = Ixz_s/(Ixx_s * Izz_s - Ixz_s^2);

A_lat = zeros(4,4);
A_lat(1,1) = Yv/m;
A_lat(1,2) = Yp/m + w0;
A_lat(1,3) = Yr/m - u0;
A_lat(1,4) = g * cos(theta0);

A_lat(2,1) = Lv/Ix_ + Izx_ * Nv;
A_lat(2,2) = Lp/Ix_ + Izx_ * Np;
A_lat(2,3) = Lr/Ix_ + Izx_ * Nr;

A_lat(3,1) = Izx_ * Lv + Nv/Iz_;
A_lat(3,2) = Izx_ * Lp + Np/Iz_;
A_lat(3,3) = Izx_ * Lr + Nr/Iz_;

A_lat(4,2) = 1;
A_lat(4,3) = tan(theta0);

Yda = qbar0 * S * CYda;
Lda = qbar0 * S * b * ClDa;
Nda = qbar0 * S * b * Cnda;

Lda = Lda - rz * Yda;
Nda = Nda + rx * Yda;

B_lat = zeros(4,1);
B_lat(1) = Yda/m;
B_lat(2) = Lda/Ix_ + Izx_ * Nda;
B_lat(3) = Izx_ * Lda + Nda/Iz_;

%% 7. Analytical transfer functions
Cp = [0 1 0 0];
Cphi = [0 0 0 1];

Hp_base = zeros(numel(f), 1);
Hphi_base = zeros(numel(f), 1);

I4 = eye(4);
for k = 1:numel(f)
    s = 1i * 2*pi*f(k);
    x_over_da = (s*I4 - A_lat) \ B_lat;

    Hp_base(k) = Cp * x_over_da;
    Hphi_base(k) = Cphi * x_over_da;
end

%% 8. Identify FRF-equivalent derivatives
% Each derivative is multiplied by a scaling factor. For every candidate, A and B matrices are rebuilt and compared with
% the two JSBSim frequency responses and the closest match with them is selected

derivative_names = {'CYb'; 'Clb'; 'Clp'; 'Clr'; 'Cnb'; ...
                    'Cnp'; 'Cnr'; 'CYda'; 'ClDa'; 'Cnda'};
base_derivatives = [CYb; Clb; Clp; Clr; Cnb; ...
                    Cnp; Cnr; CYda; ClDa; Cnda];

ratio_id = ones(10,1);
ratio_values = 0.90:0.005:1.10;
regularization = 0.20;

idx_fit = f >= lateral_fit_band(1) & f <= lateral_fit_band(2);
Hp_jsb_fit = Hp_jsb(idx_fit);
Hphi_jsb_fit = Hphi_jsb(idx_fit);
Hp_jsb_fit = Hp_jsb_fit(:);
Hphi_jsb_fit = Hphi_jsb_fit(:);

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
            CYda_s = CYda * ratio_test(8);
            ClDa_s = ClDa * ratio_test(9);
            Cnda_s = Cnda * ratio_test(10);

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
            A_test(1,1) = Yv_s/m;
            A_test(1,2) = w0;
            A_test(1,3) = -u0;
            A_test(1,4) = g * cos(theta0);
            A_test(2,1) = Lv_s/Ix_ + Izx_ * Nv_s;
            A_test(2,2) = Lp_s/Ix_ + Izx_ * Np_s;
            A_test(2,3) = Lr_s/Ix_ + Izx_ * Nr_s;
            A_test(3,1) = Izx_ * Lv_s + Nv_s/Iz_;
            A_test(3,2) = Izx_ * Lp_s + Np_s/Iz_;
            A_test(3,3) = Izx_ * Lr_s + Nr_s/Iz_;
            A_test(4,2) = 1;
            A_test(4,3) = tan(theta0);

            Yda_s = qbar0 * S * CYda_s;
            Lda_s = qbar0 * S * b * ClDa_s - rz * Yda_s;
            Nda_s = qbar0 * S * b * Cnda_s + rx * Yda_s;

            B_test = zeros(4,1);
            B_test(1) = Yda_s/m;
            B_test(2) = Lda_s/Ix_ + Izx_ * Nda_s;
            B_test(3) = Izx_ * Lda_s + Nda_s/Iz_;

            Hp_test = zeros(numel(f),1);
            Hphi_test = zeros(numel(f),1);

            for k = 1:numel(f)
                s = 1i * 2*pi*f(k);
                response = (s*I4 - A_test) \ B_test;
                Hp_test(k) = Cp * response;
                Hphi_test(k) = Cphi * response;
            end

            if any(~isfinite(Hp_test)) || any(~isfinite(Hphi_test))
                continue
            end

            Hp_test_fit = Hp_test(idx_fit);
            Hphi_test_fit = Hphi_test(idx_fit);
            err_p_test = norm(Hp_test_fit(:) - Hp_jsb_fit)/...
                max(norm(Hp_jsb_fit), eps);
            err_phi_test = norm(Hphi_test_fit(:) - Hphi_jsb_fit)/...
                max(norm(Hphi_jsb_fit), eps);

            total_error = err_p_test^2 + err_phi_test^2 + ...
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

% Build the final FRF-equivalent model using the identified factors
CYb_s = frf_equivalent_derivatives(1);
Clb_s = frf_equivalent_derivatives(2);
Clp_s = frf_equivalent_derivatives(3);
Clr_s = frf_equivalent_derivatives(4);
Cnb_s = frf_equivalent_derivatives(5);
Cnp_s = frf_equivalent_derivatives(6);
Cnr_s = frf_equivalent_derivatives(7);
CYda_s = frf_equivalent_derivatives(8);
ClDa_s = frf_equivalent_derivatives(9);
Cnda_s = frf_equivalent_derivatives(10);

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
A_frf(1,1) = Yv_s/m;
A_frf(1,2) = w0;
A_frf(1,3) = -u0;
A_frf(1,4) = g * cos(theta0);
A_frf(2,1) = Lv_s/Ix_ + Izx_ * Nv_s;
A_frf(2,2) = Lp_s/Ix_ + Izx_ * Np_s;
A_frf(2,3) = Lr_s/Ix_ + Izx_ * Nr_s;
A_frf(3,1) = Izx_ * Lv_s + Nv_s/Iz_;
A_frf(3,2) = Izx_ * Lp_s + Np_s/Iz_;
A_frf(3,3) = Izx_ * Lr_s + Nr_s/Iz_;
A_frf(4,2) = 1;
A_frf(4,3) = tan(theta0);

Yda_s = qbar0 * S * CYda_s;
Lda_s = qbar0 * S * b * ClDa_s - rz * Yda_s;
Nda_s = qbar0 * S * b * Cnda_s + rx * Yda_s;

B_frf = zeros(4,1);
B_frf(1) = Yda_s/m;
B_frf(2) = Lda_s/Ix_ + Izx_ * Nda_s;
B_frf(3) = Izx_ * Lda_s + Nda_s/Iz_;

Hp_frf = zeros(numel(f),1);
Hphi_frf = zeros(numel(f),1);

for k = 1:numel(f)
    s = 1i * 2*pi*f(k);
    response = (s*I4 - A_frf) \ B_frf;
    Hp_frf(k) = Cp * response;
    Hphi_frf(k) = Cphi * response;
end

%% 9. Tables
derivative_table = table( ...
    derivative_names, ...
    base_derivatives, ...
    frf_equivalent_derivatives, ...
    abs((frf_equivalent_derivatives - base_derivatives) ./ base_derivatives) * 100, ...
    'VariableNames', {'Derivative','Base_VSPAERO_JSBSim', ...
    'FRF_Equivalent','Difference_pct'});

Hp_base_fit = Hp_base(idx_fit);
Hphi_base_fit = Hphi_base(idx_fit);
Hp_base_fit = Hp_base_fit(:);
Hphi_base_fit = Hphi_base_fit(:);

err_p = norm(Hp_base_fit - Hp_jsb_fit)/max(norm(Hp_jsb_fit), eps);
err_phi = norm(Hphi_base_fit - Hphi_jsb_fit)/max(norm(Hphi_jsb_fit), eps);
Hp_frf_fit = Hp_frf(idx_fit);
Hphi_frf_fit = Hphi_frf(idx_fit);
err_p_frf = norm(Hp_frf_fit(:) - Hp_jsb_fit)/max(norm(Hp_jsb_fit), eps);
err_phi_frf = norm(Hphi_frf_fit(:) - Hphi_jsb_fit)/max(norm(Hphi_jsb_fit), eps);

error_table = table( ...
    {'p/da'; 'phi/da'}, ...
    [err_p; err_phi], ...
    [err_p_frf; err_phi_frf], ...
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
semilogx(f, 20*log10(abs(Hp_jsb)), 'k', ...
         f, 20*log10(abs(Hp_base)), '--', 'LineWidth', LINE_WIDTH)
grid on
set(ax1, 'TickLabelInterpreter', 'latex', 'XTickLabel', [])
ylabel('$|p/\delta_a|$ (dB)', 'Interpreter', 'latex')
%title('Aileron FRF validation', 'Interpreter', 'latex')
legend('JSBSim FRF', 'Analytical model', ...
    'Interpreter', 'latex', 'Location', 'best')

ax2 = nexttile;
semilogx(f, 20*log10(abs(Hphi_jsb)), 'k', ...
         f, 20*log10(abs(Hphi_base)), '--', 'LineWidth', LINE_WIDTH)
grid on
set(ax2, 'TickLabelInterpreter', 'latex')
ylabel('$|\phi/\delta_a|$ (dB)', 'Interpreter', 'latex')
xlabel('Frequency, $f$ (Hz)', 'Interpreter', 'latex')
legend('JSBSim FRF', 'Analytical model', ...
    'Interpreter', 'latex', 'Location', 'best')

linkaxes([ax1, ax2], 'x')
