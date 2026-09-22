close all, clear all, clc
%% FRF ELEVATOR
% Validate the elevator response using the JSBSim frequency response.
sim1 = sim('ex2400_elev');

%% 1. Read JSBSim output
file_name = 'GBS_3Out.csv';
OUT = readmatrix(file_name);

fid = fopen(file_name, 'r');
header_line = fgetl(fid);
fclose(fid);
headers = strsplit(header_line, ',');

time_col = 1;
elevator_col = 147;
q_col = 12;
az_col = 35;

t = OUT(:, time_col);
dt = median(diff(t));
fs = 1/dt;

de = OUT(:, elevator_col);              % elevator deflection, rad
q_meas = deg2rad(OUT(:, q_col));        % pitch rate, rad/s
az_meas = -OUT(:, az_col)/32.174;     % vertical acceleration, g

%% 2. Remove trim offset
% The frequency response must describe perturbations around trim
t_trim_ini = 1.0;
t_trim_fin = 2.0;
t_chirp_ini = 2.0;
t_chirp_fin = 22.0;

idx_trim = t >= t_trim_ini & t < t_trim_fin;
idx_chirp = t >= t_chirp_ini & t <= t_chirp_fin;

de1 = de(idx_chirp) - mean(de(idx_trim));
q1 = q_meas(idx_chirp) - mean(q_meas(idx_trim));
az1 = az_meas(idx_chirp) - mean(az_meas(idx_trim));

fprintf('EX2400_TRIM short-period validation\n')
fprintf('Elevator RMS: %.4f rad = %.3f deg\n\n', rms(de1), rad2deg(rms(de1)))

%% 3. JSBSim frequency response
frf_band = [0.1 5.0];               % Hz, chirp excitation range
short_period_band = [1.5 4.0];      % Hz, focused short-period validation zone
f = logspace(log10(frf_band(1)), log10(frf_band(2)), 500);

win_length = min(round(5.0 * fs), nnz(idx_chirp));
win = hamming(win_length);
noverlap = round(0.80 * win_length);

Hq_jsb = tfestimate(de1, q1, win, noverlap, f, fs);
Haz_jsb = tfestimate(de1, az1, win, noverlap, f, fs);

%% 4. Aircraft data from ex2400_trim.xml
g = 32.174;                         % ft/s^2
S = 5.68183775;                     % ft^2
c = 0.73774278;                     % ft
m = 5.5116/g;                     % slugs
Iyy = 0.19724398;                      % slug*ft^2
V0 = 29.1577 * 1.687809857;         % ft/s
alpha0 = deg2rad(0.972105);    % rad
theta0 = deg2rad(0.972105);    % rad
de_trim = 0.000057;              % rad

% Density at 100 m, ISA atmosphere.
h_m = 100;
[~, ~, ~, rho_si] = atmosisa(h_m);
rho = rho_si * 0.00194032033;       % kg/m^3 to slug/ft^3

% Use the trimmed condition from the JSBSim output when available.
% This keeps the analytical model centered on the same point as the FRF.
V_col = find(strcmp(headers, 'V_{Total} (ft/s)'), 1);
alpha_col = find(strcmp(headers, 'Alpha (deg)'), 1);
theta_col = find(strcmp(headers, 'Theta (deg)'), 1);
rho_col = find(strcmp(headers, 'Rho (slugs/ft^3)'), 1);
thrust_col = find(contains(headers, 'Thrust (engine 0 in lbs)'), 1);

if ~isempty(V_col)
    V0 = mean(OUT(idx_trim, V_col));
end
if ~isempty(alpha_col)
    alpha0 = deg2rad(mean(OUT(idx_trim, alpha_col)));
end
if ~isempty(theta_col)
    theta0 = deg2rad(mean(OUT(idx_trim, theta_col)));
end
if ~isempty(rho_col)
    rho = mean(OUT(idx_trim, rho_col));
end
de_trim = mean(de(idx_trim));
thrust_trim = 0.0;
if ~isempty(thrust_col)
    thrust_trim = mean(OUT(idx_trim, thrust_col));
end

% Longitudinal derivatives in JSBSim convention.
CLalpha = 5.5414465;
CLq     = 13.9286835;
CLde    = 0.6111759;
CLadot  = 0.9815960;
CL0     = 0.2457510;

CDo     = 0.0154292;
K       = 0.0253324;
CDde    = 0.0067909;

Cm0     = 0.2548651;
Cmalpha = -2.8119101;
Cmq     = -26.7082377;
Cmde    = -2.1799963;
Cmadot  = 0.4032911;

% Arm from AERORP to CG. JSBSim applies lift/drag at AERORP.
empty_wt = 3.7479;
tank_wt = 1.7637;
cg_empty_in = [12.7953; 0; -1.9685];
tank_in = [11.1000; 0; -1.9685];
aerorp_in = [17.2875; 0; -1.9685];
prop_in = [24.8031; 0; 4.6063];

cg_total_in = (empty_wt*cg_empty_in + tank_wt*tank_in)/(empty_wt + tank_wt);
rx = (cg_total_in(1) - aerorp_in(1))/12.0;
rz = (cg_total_in(3) - aerorp_in(3))/12.0;
rz_prop = (cg_total_in(3) - prop_in(3))/12.0;

%% 5. Identify FRF-equivalent derivatives
% This is a comparison between the VSPAERO derivatives and the equivalent derivatives from the FRF.
derivative_names = {'CLalpha'; 'CLq'; 'Cmalpha'; 'Cmq'; 'CLde'; 'Cmde'};
base_derivatives = [CLalpha; CLq; Cmalpha; Cmq; CLde; Cmde];

ratio_id = ones(6, 1);
ratio_values = 0.90:0.005:1.10;
idx_fit = f >= short_period_band(1) & f <= short_period_band(2);
regularization = 0.20;
best_total_error = inf;

for sweep = 1:4
    for ratio_number = 1:6
        best_ratio_for_this_derivative = ratio_id(ratio_number);

        for ratio_candidate = ratio_values
            ratio_test = ratio_id;
            ratio_test(ratio_number) = ratio_candidate;

            CLalpha_s = CLalpha * ratio_test(1);
            CLq_s = CLq * ratio_test(2);
            Cmalpha_s = Cmalpha * ratio_test(3);
            Cmq_s = Cmq * ratio_test(4);
            CLde_s = CLde * ratio_test(5);
            Cmde_s = Cmde * ratio_test(6);

            qbar0 = 0.5 * rho * V0^2;
            qbarS = qbar0 * S;

            u0 = V0 * cos(alpha0);
            w0 = V0 * sin(alpha0);

            ca = cos(alpha0);
            sa = sin(alpha0);
            ci2vel = c/(2 * V0);

            dalpha_du = -w0/V0^2;
            dalpha_dw =  u0/V0^2;
            dqbar_du = rho * u0;
            dqbar_dw = rho * w0;

            dalphadot_dudot = -w0/V0^2;
            dalphadot_dwdot =  u0/V0^2;

            CL_trim = CL0 + CLalpha_s * alpha0 + CLde_s * de_trim;
            CD_trim = CDo + K * CL_trim^2 + CDde * de_trim;
            Cm_trim = Cm0 + Cmalpha_s * alpha0 + Cmde_s * de_trim;

            FX_trim = -CD_trim * ca + CL_trim * sa;
            FZ_trim = -CD_trim * sa - CL_trim * ca;

            dCL_du = CLalpha_s * dalpha_du;
            dCD_du = 2 * K * CL_trim * dCL_du;
            dCm_du = Cmalpha_s * dalpha_du;

            dCL_dw = CLalpha_s * dalpha_dw;
            dCD_dw = 2 * K * CL_trim * dCL_dw;
            dCm_dw = Cmalpha_s * dalpha_dw;

            dCL_dq = CLq_s * ci2vel;
            dCD_dq = 2 * K * CL_trim * dCL_dq;
            dCm_dq = Cmq_s * ci2vel;

            dCL_de = CLde_s;
            dCD_de = 2 * K * CL_trim * dCL_de + CDde;
            dCm_de = Cmde_s;

            dCL_dudot = CLadot * ci2vel * dalphadot_dudot;
            dCD_dudot = 2 * K * CL_trim * dCL_dudot;
            dCm_dudot = Cmadot * ci2vel * dalphadot_dudot;

            dCL_dwdot = CLadot * ci2vel * dalphadot_dwdot;
            dCD_dwdot = 2 * K * CL_trim * dCL_dwdot;
            dCm_dwdot = Cmadot * ci2vel * dalphadot_dwdot;

            dFX_du = -ca*dCD_du + sa*dCL_du + (CD_trim*sa + CL_trim*ca)*dalpha_du;
            dFZ_du = -sa*dCD_du - ca*dCL_du + (-CD_trim*ca + CL_trim*sa)*dalpha_du;

            dFX_dw = -ca*dCD_dw + sa*dCL_dw + (CD_trim*sa + CL_trim*ca)*dalpha_dw;
            dFZ_dw = -sa*dCD_dw - ca*dCL_dw + (-CD_trim*ca + CL_trim*sa)*dalpha_dw;

            dFX_dq = -ca*dCD_dq + sa*dCL_dq;
            dFZ_dq = -sa*dCD_dq - ca*dCL_dq;

            dFX_de = -ca*dCD_de + sa*dCL_de;
            dFZ_de = -sa*dCD_de - ca*dCL_de;

            dFX_dudot = -ca*dCD_dudot + sa*dCL_dudot;
            dFZ_dudot = -sa*dCD_dudot - ca*dCL_dudot;

            dFX_dwdot = -ca*dCD_dwdot + sa*dCL_dwdot;
            dFZ_dwdot = -sa*dCD_dwdot - ca*dCL_dwdot;

            Xu = S * (dqbar_du * FX_trim + qbar0 * dFX_du);
            Xw = S * (dqbar_dw * FX_trim + qbar0 * dFX_dw);
            Xq = qbarS * dFX_dq;
            Xde = qbarS * dFX_de;
            Xudot = qbarS * dFX_dudot;
            Xwdot = qbarS * dFX_dwdot;

            Zu = S * (dqbar_du * FZ_trim + qbar0 * dFZ_du);
            Zw = S * (dqbar_dw * FZ_trim + qbar0 * dFZ_dw);
            Zq = qbarS * dFZ_dq;
            Zde = qbarS * dFZ_de;
            Zudot = qbarS * dFZ_dudot;
            Zwdot = qbarS * dFZ_dwdot;

            Mu = S*c*(dqbar_du*Cm_trim + qbar0*dCm_du) + rz*Xu - rx*Zu;
            Mw = S*c*(dqbar_dw*Cm_trim + qbar0*dCm_dw) + rz*Xw - rx*Zw;
            Mq = qbarS*c*dCm_dq + rz*Xq - rx*Zq;
            Mde = qbarS*c*dCm_de + rz*Xde - rx*Zde;
            Mudot = qbarS*c*dCm_dudot + rz*Xudot - rx*Zudot;
            Mwdot = qbarS*c*dCm_dwdot + rz*Xwdot - rx*Zwdot;

            % Propulsive linearization around trim. With constant-power, thrust decreases when airspeed
            % increases: dT/dV ~= -T0/V0.
            dTdV = -thrust_trim/V0;
            dTdu = dTdV * u0/V0;
            dTdw = dTdV * w0/V0;

            Xu = Xu + dTdu;
            Xw = Xw + dTdw;
            Mu = Mu + rz_prop * dTdu;
            Mw = Mw + rz_prop * dTdw;

            lhs = [m - Xudot, -Xwdot,    0;
                   -Zudot,    m - Zwdot, 0;
                   -Mudot,    -Mwdot,    Iyy];

            rhs_A = [Xu, Xw, Xq - m*w0, -m*g*cos(theta0);
                     Zu, Zw, Zq + m*u0, -m*g*sin(theta0);
                     Mu, Mw, Mq,         0];

            rhs_B = [Xde; Zde; Mde];

            A = zeros(4,4);
            B = zeros(4,1);
            A(1:3,:) = lhs \ rhs_A;
            B(1:3) = lhs \ rhs_B;
            A(4,3) = 1;

            Cq = [0 0 1 0];
            Dq = 0;
            Caz = -(A(2,:) - [0 0 u0 0])/g;
            Daz = -B(2)/g;

            Hq_model = zeros(numel(f), 1);
            Haz_model = zeros(numel(f), 1);
            I4 = eye(4);

            for k = 1:numel(f)
                s = 1i * 2*pi*f(k);
                Hq_model(k) = Cq * ((s*I4 - A) \ B) + Dq;
                Haz_model(k) = Caz * ((s*I4 - A) \ B) + Daz;
            end

            if any(~isfinite(Hq_model)) || any(~isfinite(Haz_model))
                continue
            end

            Hq_model_fit = Hq_model(idx_fit);
            Hq_jsb_fit = Hq_jsb(idx_fit);
            Haz_model_fit = Haz_model(idx_fit);
            Haz_jsb_fit = Haz_jsb(idx_fit);

            err_q = sqrt(sum(abs(Hq_model_fit(:) - Hq_jsb_fit(:)).^2)/...
                max(sum(abs(Hq_jsb_fit(:)).^2), eps));
            err_az = sqrt(sum(abs(Haz_model_fit(:) - Haz_jsb_fit(:)).^2)/...
                max(sum(abs(Haz_jsb_fit(:)).^2), eps));

            total_error = err_q^2 + err_az^2 + regularization * sum((ratio_test - 1).^2);

            if total_error < best_total_error
                best_total_error = total_error;
                best_ratio_for_this_derivative = ratio_candidate;
            end
        end

        ratio_id(ratio_number) = best_ratio_for_this_derivative;
    end
end

frf_equivalent_derivatives = base_derivatives .* ratio_id;

%% 6. Build base model for validation plots
CLalpha_s = CLalpha;
CLq_s = CLq;
Cmalpha_s = Cmalpha;
Cmq_s = Cmq;
CLde_s = CLde;
Cmde_s = Cmde;

qbar0 = 0.5 * rho * V0^2;
qbarS = qbar0 * S;

u0 = V0 * cos(alpha0);
w0 = V0 * sin(alpha0);

ca = cos(alpha0);
sa = sin(alpha0);
ci2vel = c/(2 * V0);

dalpha_du = -w0/V0^2;
dalpha_dw =  u0/V0^2;
dqbar_du = rho * u0;
dqbar_dw = rho * w0;

dalphadot_dudot = -w0/V0^2;
dalphadot_dwdot =  u0/V0^2;

CL_trim = CL0 + CLalpha_s * alpha0 + CLde_s * de_trim;
CD_trim = CDo + K * CL_trim^2 + CDde * de_trim;
Cm_trim = Cm0 + Cmalpha_s * alpha0 + Cmde_s * de_trim;

FX_trim = -CD_trim * ca + CL_trim * sa;
FZ_trim = -CD_trim * sa - CL_trim * ca;

dCL_du = CLalpha_s * dalpha_du;
dCD_du = 2 * K * CL_trim * dCL_du;
dCm_du = Cmalpha_s * dalpha_du;

dCL_dw = CLalpha_s * dalpha_dw;
dCD_dw = 2 * K * CL_trim * dCL_dw;
dCm_dw = Cmalpha_s * dalpha_dw;

dCL_dq = CLq_s * ci2vel;
dCD_dq = 2 * K * CL_trim * dCL_dq;
dCm_dq = Cmq_s * ci2vel;

dCL_de = CLde_s;
dCD_de = 2 * K * CL_trim * dCL_de + CDde;
dCm_de = Cmde_s;

dCL_dudot = CLadot * ci2vel * dalphadot_dudot;
dCD_dudot = 2 * K * CL_trim * dCL_dudot;
dCm_dudot = Cmadot * ci2vel * dalphadot_dudot;

dCL_dwdot = CLadot * ci2vel * dalphadot_dwdot;
dCD_dwdot = 2 * K * CL_trim * dCL_dwdot;
dCm_dwdot = Cmadot * ci2vel * dalphadot_dwdot;

dFX_du = -ca*dCD_du + sa*dCL_du + (CD_trim*sa + CL_trim*ca)*dalpha_du;
dFZ_du = -sa*dCD_du - ca*dCL_du + (-CD_trim*ca + CL_trim*sa)*dalpha_du;

dFX_dw = -ca*dCD_dw + sa*dCL_dw + (CD_trim*sa + CL_trim*ca)*dalpha_dw;
dFZ_dw = -sa*dCD_dw - ca*dCL_dw + (-CD_trim*ca + CL_trim*sa)*dalpha_dw;

dFX_dq = -ca*dCD_dq + sa*dCL_dq;
dFZ_dq = -sa*dCD_dq - ca*dCL_dq;

dFX_de = -ca*dCD_de + sa*dCL_de;
dFZ_de = -sa*dCD_de - ca*dCL_de;

dFX_dudot = -ca*dCD_dudot + sa*dCL_dudot;
dFZ_dudot = -sa*dCD_dudot - ca*dCL_dudot;

dFX_dwdot = -ca*dCD_dwdot + sa*dCL_dwdot;
dFZ_dwdot = -sa*dCD_dwdot - ca*dCL_dwdot;

Xu = S * (dqbar_du * FX_trim + qbar0 * dFX_du);
Xw = S * (dqbar_dw * FX_trim + qbar0 * dFX_dw);
Xq = qbarS * dFX_dq;
Xde = qbarS * dFX_de;
Xudot = qbarS * dFX_dudot;
Xwdot = qbarS * dFX_dwdot;

Zu = S * (dqbar_du * FZ_trim + qbar0 * dFZ_du);
Zw = S * (dqbar_dw * FZ_trim + qbar0 * dFZ_dw);
Zq = qbarS * dFZ_dq;
Zde = qbarS * dFZ_de;
Zudot = qbarS * dFZ_dudot;
Zwdot = qbarS * dFZ_dwdot;

Mu = S*c*(dqbar_du*Cm_trim + qbar0*dCm_du) + rz*Xu - rx*Zu;
Mw = S*c*(dqbar_dw*Cm_trim + qbar0*dCm_dw) + rz*Xw - rx*Zw;
Mq = qbarS*c*dCm_dq + rz*Xq - rx*Zq;
Mde = qbarS*c*dCm_de + rz*Xde - rx*Zde;
Mudot = qbarS*c*dCm_dudot + rz*Xudot - rx*Zudot;
Mwdot = qbarS*c*dCm_dwdot + rz*Xwdot - rx*Zwdot;

% Propulsive linearization around trim. With the constant-power
% assumption, thrust decreases when airspeed increases: dT/dV ~= -T0/V0.
dTdV = -thrust_trim/V0;
dTdu = dTdV * u0/V0;
dTdw = dTdV * w0/V0;

Xu = Xu + dTdu;
Xw = Xw + dTdw;
Mu = Mu + rz_prop * dTdu;
Mw = Mw + rz_prop * dTdw;

lhs = [m - Xudot, -Xwdot,    0;
       -Zudot,    m - Zwdot, 0;
       -Mudot,    -Mwdot,    Iyy];

rhs_A = [Xu, Xw, Xq - m*w0, -m*g*cos(theta0);
         Zu, Zw, Zq + m*u0, -m*g*sin(theta0);
         Mu, Mw, Mq,         0];

rhs_B = [Xde; Zde; Mde];

A_base = zeros(4,4);
B_base = zeros(4,1);
A_base(1:3,:) = lhs \ rhs_A;
B_base(1:3) = lhs \ rhs_B;
A_base(4,3) = 1;

Cq = [0 0 1 0];
Dq = 0;
Caz = -(A_base(2,:) - [0 0 u0 0])/g;
Daz = -B_base(2)/g;

Hq_base = zeros(numel(f), 1);
Haz_base = zeros(numel(f), 1);
I4 = eye(4);
for k = 1:numel(f)
    s = 1i * 2*pi*f(k);
    Hq_base(k) = Cq * ((s*I4 - A_base) \ B_base) + Dq;
    Haz_base(k) = Caz * ((s*I4 - A_base) \ B_base) + Daz;
end

%% 7. Print comparison tables
Hq_base_fit = Hq_base(idx_fit);
Hq_jsb_fit = Hq_jsb(idx_fit);
Haz_base_fit = Haz_base(idx_fit);
Haz_jsb_fit = Haz_jsb(idx_fit);

err_q_base = sqrt(sum(abs(Hq_base_fit(:) - Hq_jsb_fit(:)).^2)/...
    max(sum(abs(Hq_jsb_fit(:)).^2), eps));
err_az_base = sqrt(sum(abs(Haz_base_fit(:) - Haz_jsb_fit(:)).^2)/...
    max(sum(abs(Haz_jsb_fit(:)).^2), eps));

eig_base = eig(A_base);
eig_base_complex = eig_base(imag(eig_base) > 0);
[~, idx_base_sp] = max(imag(eig_base_complex));
sp_base = eig_base_complex(idx_base_sp);

fd_base_hz = imag(sp_base)/(2*pi);
fn_base_hz = abs(sp_base)/(2*pi);
zeta_base = -real(sp_base)/abs(sp_base);

f_fit = f(idx_fit);
[~, idx_q_peak] = max(abs(Hq_jsb_fit));
[~, idx_az_peak] = max(abs(Haz_jsb_fit));

fq_frf_hz = f_fit(idx_q_peak);
faz_frf_hz = f_fit(idx_az_peak);
q_peak_db = 20*log10(abs(Hq_jsb_fit(idx_q_peak)));
az_peak_db = 20*log10(abs(Haz_jsb_fit(idx_az_peak)));

derivative_table = table( ...
    derivative_names, ...
    base_derivatives, ...
    frf_equivalent_derivatives, ...
    abs((frf_equivalent_derivatives - base_derivatives) ./ base_derivatives) * 100, ...
    'VariableNames', {'Derivative','Base_VSPAERO_JSBSim','FRF_Equivalent','Difference_pct'});

fprintf('Short-period validation band: %.1f-%.1f Hz\n\n', short_period_band(1), short_period_band(2))

disp('Base derivatives and FRF-equivalent derivatives:')
disp(derivative_table)

%% 8. Plots

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

Hq_jsb_plot = smoothdata(20*log10(abs(Hq_jsb)), 'movmedian', 7);
Haz_jsb_plot = smoothdata(20*log10(abs(Haz_jsb)), 'movmedian', 7);
Hq_base_plot = 20*log10(abs(Hq_base));
Haz_base_plot = 20*log10(abs(Haz_base));

ax1 = nexttile;
semilogx(f, Hq_jsb_plot, 'k', ...
         f, Hq_base_plot, '--', 'LineWidth', LINE_WIDTH)
grid on
set(ax1, 'TickLabelInterpreter', 'latex', 'XTickLabel', [])
ylabel('$|q/\delta_e|$ (dB)', 'Interpreter', 'latex')
% title('Elevator FRF validation', 'Interpreter', 'latex')
legend('JSBSim FRF', 'Analytical model', ...
    'Interpreter', 'latex', 'Location', 'best')

ax2 = nexttile;
semilogx(f, Haz_jsb_plot, 'k', ...
         f, Haz_base_plot, '--', 'LineWidth', LINE_WIDTH)
grid on
set(ax2, 'TickLabelInterpreter', 'latex')
ylabel('$|a_z/\delta_e|$ (dB)', 'Interpreter', 'latex')
xlabel('Frequency, $f$ (Hz)', 'Interpreter', 'latex')
legend('JSBSim FRF', 'Analytical model', ...
    'Interpreter', 'latex', 'Location', 'best')

linkaxes([ax1, ax2], 'x')

exportgraphics(gcf, 'elev0_corrected.png', 'Resolution', 300)
