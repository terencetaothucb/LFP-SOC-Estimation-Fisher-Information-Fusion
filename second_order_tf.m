clear all
% Define the parameters lambda_0, lambda_1, and the sampling time T

lambda_0 = 0.4; 
lambda_1 = 9;
Ts = 1; % Sampling period

s = tf('s');
Gs2 = (lambda_0*s*s)/(s*s + lambda_1 *s + lambda_0);
Gs1 = (lambda_0*s)/(s*s + lambda_1 *s + lambda_0);
Gs0 = (lambda_0)/(s*s + lambda_1 *s + lambda_0);

% Discretize using zero-order hold (ZOH)
Gz2 = c2d(Gs2, Ts, 'zoh');
Gz1 = c2d(Gs1, Ts, 'zoh');
Gz0 = c2d(Gs0, Ts, 'zoh');
% Get the state-space representation of the discrete transfer function
[A2, B2, C2, D2] = ssdata(Gz2);
[A1, B1, C1, D1] = ssdata(Gz1);
[A0, B0, C0, D0] = ssdata(Gz0);



% Print the state-space matrices which represent the difference equations
save('Configuration Set/second_order_matrix.mat', 'A0', 'A1', 'A2', 'B0', 'B1', 'B2', 'C0', 'C1', 'C2', 'D0', 'D1', 'D2');


