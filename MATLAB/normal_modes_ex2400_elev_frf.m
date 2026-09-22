close all; clc; clear all

%% ELEVATOR FRF-EQUIVALENT LONGITUDINAL MODES

%% Aircraft and trim data
g = 32.174;                    % ft/s^2

% Geometry
S = 5.68183775;                % ft^2
b = 7.87401575;                % ft
c = 0.73774278;                % ft

% Mass and inertia in body axes
m = 5.5116 / g;                % slugs, total mass: 1.7 kg empty + 0.8 kg battery
Ixx = 0.2529;                  % slug*ft^2
Iyy = 0.19724398;              % slug*ft^2 
Izz = 0.47654398;              % slug*ft^2
Ixz = 0.01629;                 % slug*ft^2

% Trim condition
h0_ft = 328.084;               % ft
V0 = 29.1577 * 1.687809857;    % ft/s
alpha0 = deg2rad(0.972105);    % rad, trim from ex2400_trim
theta0 = deg2rad(0.972105);    % rad, trim from ex2400_trim

% Trimmed controls
de0 = 0.000057;                 % rad, trim from ex2400_trim
da0 = -0.000444;                % rad, trim from ex2400_trim
dr0 = -0.000001;                % rad, trim from ex2400_trim

% Atmosphere
h0_m = h0_ft * 0.3048;
[~, ~, ~, rho_si] = atmosisa(h0_m);
rho = rho_si * 0.00194032033;  % kg/m^3 to slug/ft^3

%% Aerodynamic coefficients from ex2400_trim.xml
% Original (VSPAERO) values, for reference:
%   CLalpha=5.5414465  CLq=13.9286835  Cmalpha=-2.8119101
%   Cmq=-26.7082377    CLde=0.6111759  Cmde=-2.1799963
CLalpha = 5.2921;
CLq     = 13.1630;
CLadot  = 0.9815960;
CL0_base = 0.2457510;

CDo     = 0.0154292;
K       = 0.0253324;
CDde    = 0.0067909;

Cm0     = 0.2548651;
Cmalpha = -2.8400;
Cmq     = -24.0370;
Cmde    = -2.2345;
Cmadot  = 0.4032911;
CLde    = 0.63562;

CYb     = -0.5316963;

% Lateral derivatives converted from VSPAERO to JSBSim body axes
Clb     = -0.5440491;
Clp     = -0.7119247;
Clr     =  0.0582024;

Cnb     =  0.3087214;
Cnp     = -0.0982054;
Cnr     = -0.0798044;

%% JSBSim AERORP to CG arm
% JSBSim applies aerodynamic forces at AERORP and transfers the moments to the CG
empty_wt = 3.7479;
tank_wt = 1.7637;
cg_empty_in = [12.7953; 0.0000; -1.9685];
tank_in = [11.1000; 0.0000; -1.9685];
aerorp_in = [17.2875; 0.0000; -1.9685];
cg_total_in = (empty_wt*cg_empty_in + tank_wt*tank_in) / (empty_wt + tank_wt);

rx = (cg_total_in(1) - aerorp_in(1)) / 12.0;
rz = (cg_total_in(3) - aerorp_in(3)) / 12.0;

%% Inertia in stability axes
Ixx_s = Ixx * cos(alpha0)^2 + Izz * sin(alpha0)^2 + Ixz * sin(2 * alpha0);
Izz_s = Ixx * sin(alpha0)^2 + Izz * cos(alpha0)^2 - Ixz * sin(2 * alpha0);
Ixz_s = Ixz * (cos(alpha0)^2 - sin(alpha0)^2) - 0.5 * (Ixx - Izz) * sin(2 * alpha0);
Iyy_s = Iyy;

%% Longitudinal stability derivatives
qbar0 = 0.5 * rho * V0^2;
qbarS = qbar0 * S;
qbarS_V = qbarS / V0;

u0 = V0 * cos(alpha0);
w0 = V0 * sin(alpha0);

ca = cos(alpha0);
sa = sin(alpha0);
ci2vel = c / (2 * V0);

dalpha_du = -w0 / V0^2;
dalpha_dw =  u0 / V0^2;
dqbar_du = rho * u0;
dqbar_dw = rho * w0;

dalphadot_dudot = -w0 / V0^2;
dalphadot_dwdot =  u0 / V0^2;

CL0 = CL0_base + CLalpha * alpha0 + CLde * de0;
CD0 = CDo + K * CL0^2 + CDde * de0;
Cm_trim = Cm0 + Cmalpha * alpha0 + Cmde * de0;

FX0 = -CD0 * ca + CL0 * sa;
FZ0 = -CD0 * sa - CL0 * ca;

%% Derivatives of CL, CD and Cm
% With respect to u
dCL_du = CLalpha * dalpha_du;
dCD_du = 2 * K * CL0 * dCL_du;
dCm_du = Cmalpha * dalpha_du;

% With respect to w
dCL_dw = CLalpha * dalpha_dw;
dCD_dw = 2 * K * CL0 * dCL_dw;
dCm_dw = Cmalpha * dalpha_dw;

% With respect to q
dCL_dq = CLq * ci2vel;
dCD_dq = 2 * K * CL0 * dCL_dq;
dCm_dq = Cmq * ci2vel;

% With respect to elevator
dCL_de = CLde;
dCD_de = 2 * K * CL0 * dCL_de + CDde;
dCm_de = Cmde;

% With respect to udot
dCL_dudot = CLadot * ci2vel * dalphadot_dudot;
dCD_dudot = 2 * K * CL0 * dCL_dudot;
dCm_dudot = Cmadot * ci2vel * dalphadot_dudot;

% With respect to wdot
dCL_dwdot = CLadot * ci2vel * dalphadot_dwdot;
dCD_dwdot = 2 * K * CL0 * dCL_dwdot;
dCm_dwdot = Cmadot * ci2vel * dalphadot_dwdot;

%% Derivatives of X and Z
dFX_du = -ca*dCD_du + sa*dCL_du + (CD0*sa + CL0*ca)*dalpha_du;
dFZ_du = -sa*dCD_du - ca*dCL_du + (-CD0*ca + CL0*sa)*dalpha_du;

dFX_dw = -ca*dCD_dw + sa*dCL_dw + (CD0*sa + CL0*ca)*dalpha_dw;
dFZ_dw = -sa*dCD_dw - ca*dCL_dw + (-CD0*ca + CL0*sa)*dalpha_dw;

dFX_dq = -ca*dCD_dq + sa*dCL_dq;
dFZ_dq = -sa*dCD_dq - ca*dCL_dq;

dFX_de = -ca*dCD_de + sa*dCL_de;
dFZ_de = -sa*dCD_de - ca*dCL_de;

dFX_dudot = -ca*dCD_dudot + sa*dCL_dudot;
dFZ_dudot = -sa*dCD_dudot - ca*dCL_dudot;

dFX_dwdot = -ca*dCD_dwdot + sa*dCL_dwdot;
dFZ_dwdot = -sa*dCD_dwdot - ca*dCL_dwdot;

Xu = S * (dqbar_du * FX0 + qbar0 * dFX_du);
Xw = S * (dqbar_dw * FX0 + qbar0 * dFX_dw);
Xq = qbarS * dFX_dq;
Xde = qbarS * dFX_de;
Xudot = qbarS * dFX_dudot;
Xwdot = qbarS * dFX_dwdot;

Zu = S * (dqbar_du * FZ0 + qbar0 * dFZ_du);
Zw = S * (dqbar_dw * FZ0 + qbar0 * dFZ_dw);
Zq = qbarS * dFZ_dq;
Zde = qbarS * dFZ_de;
Zudot = qbarS * dFZ_dudot;
Zwdot = qbarS * dFZ_dwdot;

%% Derivatives of M
Mu = S*c*(dqbar_du*Cm_trim + qbar0*dCm_du) + rz*Xu - rx*Zu;
Mw = S*c*(dqbar_dw*Cm_trim + qbar0*dCm_dw) + rz*Xw - rx*Zw;
Mq = qbarS*c*dCm_dq + rz*Xq - rx*Zq;
Mde = qbarS*c*dCm_de + rz*Xde - rx*Zde;
Mudot = qbarS*c*dCm_dudot + rz*Xudot - rx*Zudot;
Mwdot = qbarS*c*dCm_dwdot + rz*Xwdot - rx*Zwdot;

%% Longitudinal matrix
% States: x_long = [u, w, q, theta]'
lhs_long = [m - Xudot, -Xwdot,       0;
            -Zudot,    m - Zwdot,    0;
            -Mudot,    -Mwdot,       Iyy_s];

rhs_long = [Xu, Xw, Xq - m*w0, -m*g*cos(theta0);
            Zu, Zw, Zq + m*u0, -m*g*sin(theta0);
            Mu, Mw, Mq,         0];

A_long = zeros(4,4);
A_long(1:3,:) = lhs_long \ rhs_long;
A_long(4,1) = 0;
A_long(4,2) = 0;
A_long(4,3) = 1;
A_long(4,4) = 0;

%% Lateral-directional stability derivatives

CYp = 0;
CYr = 0;

Yv = qbarS_V * CYb;
Yp = 0.5 * qbarS_V * b * CYp;
Yr = 0.5 * qbarS_V * b * CYr;

Lv = qbarS_V * b * Clb;
Lp = 0.5 * qbarS_V * b^2 * Clp;
Lr = 0.5 * qbarS_V * b^2 * Clr;

Nv = qbarS_V * b * Cnb;
Np = 0.5 * qbarS_V * b^2 * Cnp;
Nr = 0.5 * qbarS_V * b^2 * Cnr;

% Moment transfer due to side force applied at AERORP
Lv = Lv - rz * Yv;
Lp = Lp - rz * Yp;
Lr = Lr - rz * Yr;

Nv = Nv + rx * Yv;
Np = Np + rx * Yp;
Nr = Nr + rx * Yr;

Ix_ = (Ixx_s * Izz_s - Ixz_s^2) / Izz_s;
Iz_ = (Ixx_s * Izz_s - Ixz_s^2) / Ixx_s;
Izx_ = Ixz_s / (Ixx_s * Izz_s - Ixz_s^2);

%% Lateral-directional matrix
% States: x_lat = [v, p, r, phi]'
A_lateral = zeros(4,4);

A_lateral(1,1) = Yv / m;
A_lateral(1,2) = Yp / m + w0;
A_lateral(1,3) = Yr / m - u0;
A_lateral(1,4) = g * cos(theta0);

A_lateral(2,1) = Lv / Ix_ + Izx_ * Nv;
A_lateral(2,2) = Lp / Ix_ + Izx_ * Np;
A_lateral(2,3) = Lr / Ix_ + Izx_ * Nr;
A_lateral(2,4) = 0;

A_lateral(3,1) = Izx_ * Lv + Nv / Iz_;
A_lateral(3,2) = Izx_ * Lp + Np / Iz_;
A_lateral(3,3) = Izx_ * Lr + Nr / Iz_;
A_lateral(3,4) = 0;

A_lateral(4,1) = 0;
A_lateral(4,2) = 1;
A_lateral(4,3) = tan(theta0);
A_lateral(4,4) = 0;

%% Eigenvalues and eigenvectors
[evec_long, eval_long] = eig(A_long);
[evec_lat, eval_lat] = eig(A_lateral);

lambda_lon = diag(eval_long);
lambda_lat = diag(eval_lat);

real_lon = real(lambda_lon);
imag_lon = imag(lambda_lon);
wn_lon = sqrt(real_lon.^2 + imag_lon.^2);
zeta_lon = -real_lon ./ max(wn_lon, eps);

period_lon = nan(size(lambda_lon));
osc_lon = abs(imag_lon) > 1e-9;
period_lon(osc_lon) = 2*pi ./ abs(imag_lon(osc_lon));

real_lat = real(lambda_lat);
imag_lat = imag(lambda_lat);
wn_lat = sqrt(real_lat.^2 + imag_lat.^2);
zeta_lat = -real_lat ./ max(wn_lat, eps);

period_lat = nan(size(lambda_lat));
osc_lat = abs(imag_lat) > 1e-9;
period_lat(osc_lat) = 2*pi ./ abs(imag_lat(osc_lat));

%% Results
disp('A_longitudinal (FRF-equivalent)')
disp(A_long)

disp('A_lateral')
disp(A_lateral)

disp('Longitudinal eigenvalues (FRF-equivalent):')
disp(lambda_lon)

longitudinal_data = table(lambda_lon, real_lon, imag_lon, wn_lon, zeta_lon, period_lon, ...
    'VariableNames', {'Eigenvalue','Real','Imag','NaturalFreq','DampingRatio','Period'});

disp('Longitudinal data (FRF-equivalent)')
disp(longitudinal_data)

disp('Lateral-directional eigenvalue')
disp(lambda_lat)

lateral_data = table(lambda_lat, real_lat, imag_lat, wn_lat, zeta_lat, period_lat, ...
    'VariableNames', {'Eigenvalue','Real','Imag','NaturalFreq','DampingRatio','Period'});

disp('Lateral-directional data')
disp(lateral_data)

%% Pole plots
figure();
plot(real(lambda_lon), imag(lambda_lon), 'x', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
xlabel('Real part (1/s)');
ylabel('Imag part (rad/s)');
title('Longitudinal Poles (FRF-equivalent)');
