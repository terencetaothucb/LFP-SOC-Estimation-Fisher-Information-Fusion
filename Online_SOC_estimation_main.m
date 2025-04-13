%% Summary:

%  Inputs: measured terminal voltage and current

%  Author: Junzhe Shi

%  Date created       : 3/18/2023

%  Date last modified : 02/09/2025

%  Modified by        : Junzhe Shi

clear all
close all

labelsize = 20;
legendsize = 20;
titlesize = 20;
linewidth = 3;
set(0,'DefaultFigureWindowStyle','docked')

%% Load configurations
% Load secod order filter
load('Configuration Set/second_order_matrix.mat');

% load ocv soc model
load('Configuration Set/Real_SOC_OCV_Curve.mat');
load('Configuration Set/OCV_H_to_SOC_map.mat')
load('Configuration Set/SOC_H_to_dSOCdOCV_map.mat');

%% Load test data
load('Test Data Set/Long_flat_zone_data_25_degree.mat');
% load('Test Data Set/Long_flat_zone_data_10_degree.mat');
% load('Test Data Set/With_constant_current_25_degree.mat');

True_data = battery_tests{1};

%% Error set up
SOC_guess_intial = 1; % [%]
I_bias = 0; % [A]
Add_Voltage_Quantization_error = false; 

%% Check data field
initalSOC = min(1,True_data.SOC(1));

if(~isfield(True_data,'Current_measured'))
    True_data.Current_measured = True_data.Current;
    True_data.VT_measured = True_data.VT;
end

%% System parameter set up
Cp = 1.2; % cpacity [Ah]
dt = 1; % sampling time [s]
N = length(True_data.SOC); % total time time length

%% Quantization
V_supply = 5; % Supply voltage in volts
n_bits = 10; % Number of bits in ADC
n_levels = 2^n_bits; % Calculate number of levels
quantization_step = V_supply / (n_levels - 1); 

%% initial SOC estimation parameters
% parameter indentification
linear_regression_window_size = 300; % window size
linear_regression_input = zeros(6,linear_regression_window_size);
linear_regression_output = zeros(1,linear_regression_window_size);
VT_x_next = zeros(size(A0, 1), 1);
I_x_next= zeros(size(A0, 1), 1);
VT_dot_x_next= zeros(size(A1, 1), 1);
I_dot_x_next= zeros(size(A1, 1), 1);
VT_2dot_x_next= zeros(size(A2, 1), 1);
I_2dot_x_next= zeros(size(A2, 1), 1);
OCV_est = zeros(N,1);

% Hysteresis temrs
CC = 550; % Hysteresis fatting parameter
tau1 = 200 * 0.05; % current filter time constant
a1 = exp(-abs(1/tau1)); % current filter discretized factor
I1_est = zeros(N,1); % filter current
I_prev = 0; % pervious current
H_est = zeros(N,1);
if (True_data.VT(1) > 3.29)
    H_est(1) = 1;  % start from charging
else
    H_est(1) = -1; % start from discharging
end
% Fisher information matrix 
sigm_V = 0.01; % voltage measurement variance
sigm_V_s = sigm_V * sigm_V;
S = zeros(linear_regression_window_size, 6);
S_prev = zeros(linear_regression_window_size, 6);
F = 1/sigm_V_s * S' * S;
F_inv = inv(F);
largeNumber = 1e10;
F_inv(isinf(F_inv)) = largeNumber;% Replace Inf with the large number
cov_SOCm_save = zeros(N, 1);

% SOC
SOC = SOC_guess_intial;
SOC_est = zeros(N,1); 
SOC_m= zeros(N,1);
P_SOC = 2; % inital SOC covariance
sigm_vv_SOC = (1 / Cp / 3600)^2 * 0.0002;% processing noise, Y = aX + b -> Var(Y) = a^2 * Var(X) 
sigm_ww_OCV = zeros(N,1);
P_SOC_save = zeros(N,1);

%% Run Online estimation
Current = True_data.Current';
V_measure = True_data.VT';

tic
for i = 1: N
    
    % add current bias
    I = Current(i) + I_bias;

    % add voltage quantization error
    VT = V_measure(i);

    if(Add_Voltage_Quantization_error)
        VT = round(VT / quantization_step) * quantization_step;
    end

    % apply the discretized transfter function to to get processed voltage
    % and current data
    VT_x = VT_x_next;
    I_x = I_x_next;
    VT_dot_x = VT_dot_x_next;
    I_dot_x = I_dot_x_next;
    VT_2dot_x = VT_2dot_x_next;
    I_2dot_x = I_2dot_x_next;

    VT_x_next = A0 * VT_x + B0 * VT;
    VT_y = C0 * VT_x + D0 * VT;
    I_x_next = A0 * I_x + B0 * I;
    I_y = C0 * I_x + D0 * I;

    VT_dot_x_next = A1 * VT_dot_x + B1 * VT;
    VT_dot_y = C1 * VT_dot_x + D1 * VT;
    I_dot_x_next = A1 * I_dot_x + B1 * I;
    I_dot_y = C1 * I_dot_x + D1 * I;

    VT_2dot_x_next = A2 * VT_2dot_x + B2 * VT;
    VT_2dot_y = C2 * VT_2dot_x + D2 * VT;
    I_2dot_x_next = A2 * I_2dot_x + B2 * I;
    I_2dot_y = C2 * I_2dot_x + D2 * I;

    % run linear regression
    x = [1, -I_2dot_y, -I_dot_y, -I_y, -VT_2dot_y, -VT_dot_y]';
    rg_index = mod(i,linear_regression_window_size) + 1;
    linear_regression_input(:,rg_index) = x;
    linear_regression_output(1,rg_index) = VT_y;
    thetaest = linear_regression_input' \ linear_regression_output';
    
    % update the hysteresis term
    if i > 1
        I1_est(i) = a1 * I1_est(i - 1) + (1 - a1) * I_prev; % I > 0 is discharge, 
        % A low pass filter is used for the current sign here to
        % eliminate the impact of measurement noise on the charging and discharging directions 
        H_est(i) = exp(-abs(I/CC)) * H_est(i - 1) + (1 - exp(-abs(I/CC))) * sign(-I1_est(i));
    end

    % get estimated OCV for parameter estimation
    OCV_est(i) = thetaest(1);
    OCV_est(i) = max(OCV_range(1),min(OCV_range(end),OCV_est(i))); % keep the OCV in the constraint

    % get SOC_{OCV-H} from the OCV-H-SOC inversion
    SOC_m(i) = interp2(OCV_range, H_range,OCV_H_to_SOC_map', OCV_est(i), H_est(i))/100;
    SOC_m(i) = max(0,min(1,SOC_m(i))); % keep the SOC in the constraint

    % update the sensitivity matrix
    S(1,:) = x';
    S(2:end,:) = S_prev(1:end - 1,:);
    S_prev = S;

    % calculate the fisher information matrix
    F = 1/sigm_V_s * S' * S;
    F_inv = inv(F);
    F_inv(isinf(F_inv)) = largeNumber; % replace Inf with the large number to avoid singular issues

    % get the cov of OCV_{OCV-H}
    sigm_ww_OCV(i) = F_inv(1);

    % get the dsoc/docv value with the given soc and hysteresis term
    dSOCdOCV = interp2(SOC_range, H_range, SOC_H_to_dSOCdOCV_map', SOC * 100, H_est(i))/100;

    % get the cov of SOC_{OCV-H}
    sigm_ww_SOC = dSOCdOCV^2 * sigm_ww_OCV(i);

    % estimate SOC using KF
    SOC = SOC - I / Cp / 3600; % Coulomb counting
    H_SOC = 1; % H matrix in the KF

    % run KF after the 5 steps to ensure the linear regression work
    if i >= length(thetaest) * 5
        % measurement covariance update
        P_SOC = P_SOC + sigm_vv_SOC;
        P_SOC_save(i) = P_SOC;

        % Kalman gain calculation
        K_SOC = P_SOC * H_SOC' * inv(H_SOC * P_SOC * H_SOC' + sigm_ww_SOC);

        % state update
        SOC = SOC + K_SOC * (SOC_m(i) - H_SOC*SOC);

        % covariance update 
        P_SOC = (eye(size(P_SOC)) - K_SOC*H_SOC) * P_SOC * (eye(size(P_SOC)) - K_SOC*H_SOC)' + K_SOC * sigm_ww_SOC * K_SOC';
    else
        K_SOC = P_SOC * H_SOC' * inv(H_SOC * P_SOC * H_SOC' + sigm_ww_SOC);
    end

    SOC = max(0,min(1,SOC)); % keep the SOC in the constraint

    % save data
    I_prev = I;
    SOC_est(i) = SOC;
    cov_SOCm_save(i) = sigm_ww_SOC;
end

fprintf('Total solving time: %i\n',toc);
Total_solving_time = toc;

%% Plot results
range = 1:i;
figure;
subplot(2,1,1)
plot(range/ 3600,SOC_m(range) .* 100,'-','LineWidth',2)
hold on
plot([0, range]'/ 3600,[SOC_guess_intial; SOC_est(range) ].* 100,'LineWidth',2)
hold on
plot(range/ 3600,True_data.SOC(range) .* 100,'LineWidth',2)
ylabel('SOC [%]','FontWeight','bold','FontSize',labelsize)
legend('SOC_{OCV-h}','SOC_{Est}','SOC_{True}','FontWeight','bold','FontSize',legendsize)
ax = gca;
ax.FontSize = legendsize;
grid on
subplot(2,1,2)
plot([0, range]' ./ 3600,abs([initalSOC - SOC_guess_intial;(True_data.SOC(range) - SOC_est(range))] .* 100),'LineWidth',2)
xlabel('Time [h]','FontWeight','bold','FontSize',labelsize)
ylabel('Abs Error [%]','FontWeight','bold','FontSize',labelsize)
legend('Error','FontWeight','bold','FontSize',legendsize)
grid on
ax = gca;
ax.FontSize = legendsize;


figure;
subplot(2,1,1)
semilogy([range]' ./ 3600,cov_SOCm_save(range),'LineWidth',linewidth)
ylabel('Cov (log scale)','FontWeight','bold','FontSize',labelsize)
legend('Cov_{SOC_{OCV-H}}','FontWeight','bold','FontSize',legendsize)
grid on
ax = gca;
ax.FontSize = legendsize;
subplot(2,1,2)
plot([range]' ./ 3600,Current(range),'--','LineWidth',2)
xlabel('Time [h]','FontWeight','bold','FontSize',labelsize)
ylabel('Crate','FontWeight','bold','FontSize',labelsize)
legend('Current','FontWeight','bold','FontSize',legendsize)
ax = gca;
ax.FontSize = legendsize;
grid on

RMSE = rmse(True_data.SOC(range) .* 100, SOC_est(range) .* 100)

