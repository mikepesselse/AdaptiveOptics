%% Filtering assignment
% Mike Pesselse [4300564]
% Bart de Jong  [4367146]

clearvars
close all

%% Load data
load systemMatrices.mat
load turbulenceData.mat

%% Method 1
% Exercise 1.6
C = zeros(49,49);
var_e1 = zeros(1, 20);
mean_e = zeros(49, 20);
var_e1_nc = zeros(1, 20);

for j = 1:20
    phik = phiSim{j};
    
    for i = 1:length(phik)
        C = C + phik(:,i)*phik(:,i)';
    end
    
    C = C/length(phik);    
    
    [var_e1(j)] = AOloopRW(G, H, C, sigmae, phik);
    [var_e1_nc(j)] = AOloop_nocontrol(phik, sigmae, H, G);
end

figure; bar(1:20, [var_e1_nc; var_e1]);
title('Exercise 1: Mean variance for every dataset')b
legend('Without control', 'With control')


%% Method 2
% Exercise 2.6

% Initialise vectors
C_phi_0 = zeros(49, 49, 20);
C_phi_1 = zeros(49, 49, 20);
A       = zeros(49, 49, 20);
Cw      = zeros(49, 49, 20);
K       = zeros(49, 72, 20);

for j = 1:20
    
    phik = phiSim{j};
    
    % Calculate C_phi_0
    phicov = zeros(49, 49);
    for i=1:5000
        phicov = phicov + phiSim{j}(:,i)*phiSim{j}(:,i)';
    end
    C_phi_0(:, :, j) = phicov/5000;
    
    %Calculate C_phi_1
    phicov = zeros(49, 49);
    for i=1:4999
        phicov = phicov + phiSim{j}(:,i+1)*phiSim{j}(:,i)';
    end
    C_phi_1(:, :, j) = phicov/4999;
    
    [A(:,:,j), Cw(:,:,j), K(:,:,j)] = computeKalmanAR(C_phi_0(:, :, j), C_phi_1(:, :, j), G, sigmae);
    [var_e2(j)] = AOloopAR(G, H, C_phi_0(:, :, j), sigmae, A(:, :, j), Cw(:, :, j), K(:, :, j), phik);
end

figure; bar(1:20, [var_e1_nc; var_e1; var_e2]);
title('Exercise 2: Mean variance for every dataset')
legend('Without control', 'With control', 'With Kalman control')


%% Method 3
% Exercise 3.5
n = 100;
p = 6;
s = 5;

As = zeros(n, n, 20);
Cs = zeros(72, n, 20);
Ks = zeros(n, 72, 20);

for j=1:20
    
    s_o = G*phiIdent{j} + sigmae*randn(length(G), 5000);
    [As(:,:,j), Cs(:,:,j), Ks(:,:,j)] = SubId(s_o, s, n);
    
    [var_e3(j)]        = AOloopSID(G, H, As(:, :, j), Cs(:, :, j), Ks(:, :, j), sigmae, phiSim{j});
end

figure; bar(1:20, [var_e1_nc; var_e1; var_e2; var_e3]);
title('Exercise 2: Mean variance for every dataset')
legend('Without control', 'With control', 'With Kalman control', 'With SID')
%% Functions
function [var_e]        = AOloopRW(G, H, C, sigma, phiSim)

phi     = phiSim - ones(49, 1)*mean(phiSim, 1);
[m, p]  = size(phi);

P = C*G'/(G*C*G'+sigma^2*eye(length(G)));   % [49x72]

epsilon = zeros(m, p);
u       = zeros(m, p);
s       = zeros(length(G),p);

for k = 2:p
    epsilon(:,k) = phi(:,k)-H*u(:,k-1);
    s(:,k) = G*epsilon(:,k) + sigma*randn(length(G), 1);
    delta_u = H \ P * s(:,k);
    u(:,k) = delta_u + u(:,k-1);
end

mean_e1 = mean(epsilon, 1);     % Calculate mean for each timestep
var_e = var(epsilon - ones(m, 1)*mean_e1, 0, 2);
var_e = mean(var_e, 1);

end

function [A, Cw, K]     = computeKalmanAR(C_phi_0, C_phi_1, G, sigma)

% If idare is not working, use the following code:
% A = C_phi_1 / C_phi_0;
% Cw = C_phi_0 - A * C_phi_0 * A';
% [~,~,K] = dare(A',G', Cw, sigma^2*eye(72));
% K = K';

A  = C_phi_1/C_phi_0;
Cw = C_phi_0 - A*C_phi_0*A';
%Cw = (Cw+Cw')/2;
[~, K, ~] = idare(A', G', Cw, sigma^2*eye(length(G)), 0, []);
K = K';
end

function [var_e]        = AOloopAR(G, H, C_phi_0, sigma, A, Cw, K, phiSim)

%Define dimension
phi     = phiSim;
[m, p]  = size(phi);

epsilon = zeros(m, p);
u       = zeros(m, p);
s       = zeros(length(G),p);

epsilonreal = zeros(m, p);

% Closed loop

epsilonreal(:,1) = phi(:,1);
s(:, 1) = G*epsilonreal(:,1) + sigma*randn(length(G), 1);
s(:, 2) = G*epsilonreal(:,2) + sigma*randn(length(G), 1);
u(:, 1) = H'*H \ H' * K*s(:,1);

for k = 2:p-1
    u(:,k) = H'*H \ H' * ((A-K*G) * epsilon(:,k) + A*H*u(:,k-1) + K*s(:,k));
    epsilon(:,k+1) = (A-K*G) * epsilon(:, k) - H*u(:,k) + A*H*u(:,k-1) + K*s(:,k);
    epsilonreal(:,k+1) = phi(:, k+1) - H*u(:, k);
    s(:, k+1) = G*epsilonreal(:,k+1) + sigma*randn(length(G), 1);
end

mean_e1_real = mean(epsilonreal, 1);     % Calculate mean for each timestep
var_e = var(epsilonreal - ones(m, 1)*mean_e1_real, 0, 2);
var_e = mean(var_e, 1);

end

function [As, Cs, Ks]   = SubId(s_o, s, n)

%%%%%%%%%%%%
Ncol = size(s_o, 2);
sid = s_o(:, 1:Ncol);
[sid_r, sid_c] = size(sid);
N = Ncol -s +2;

Y_0 = [];
Y_s = [];
for k=1:s
    Y_vec = sid(:, k:k+N-2);
    Y_0 = [Y_0; Y_vec];
end

for k=s+1:2*s
    Y_vec = sid(:, k:N+k-s-2);
    Y_s = [Y_s; Y_vec];
end

%plot of the singular values of S method 1
[~,SY0,~] = svd(Y_0);
figure; semilogy(diag(SY0), 'xr');

Zn = Y_0(:, 1:length(Y_s));    % Resolve concatenation errors
RQmat = [Zn; Y_s];

[~,R] = qr(RQmat');
R = R';

R11 = R(1:s*sid_r, 1:s*sid_r);
R21 = R(s*sid_r+1 : 2*s*sid_r, 1:s*sid_r);

SVDmat    = R21/R11*Zn;
[~, S, V] = svd(SVDmat);

%plot of the singular values of S method 2
[S_r,~] = size(S);
figure('Position', [100, 100, 1500, 600])
scatter(1:S_r,diag(S),'*')
title('Singular Values, s = ?')
xlim([0, 200])
grid on

XsN = sqrtm(S(1:n, 1:n))*V(:, 1:n)';
F = XsN(:, 1:end-1);
xA = XsN(:, 2:end);
yC = Y_s(1:sid_r, 1:end-1);
yC = sid(:, s+1:end-s);

As = xA * F' /(F*F');    % Frobenius norm = y*f'*(f*f')^-1
Cs = yC * F' /(F*F');    % Frobenius norm

Ws = XsN(:, 2:end)         - As*XsN(:, 1:end-1);
Vs = Y_s(1:sid_r, 1:end-1) - Cs*XsN(:, 1:end-1);

Q_hat = Ws*Ws'/(Ncol-s);
R_hat = Vs*Vs'/(Ncol-s);
S_hat = Ws*Vs'/(Ncol-s);

[~, K, ~] = idare(As', Cs', Q_hat, R_hat, S_hat, []);
Ks = K';

end

function [var_e]        = AOloopSID(G, H, As, Cs, Ks, sigmae, phiSim)

[U, S, V] = svd(G);
rankS = rank(S);
Gamma = V(:, 1:rankS)/S(1:rankS, 1:rankS)*U(:, 1:rankS)';

phi = phiSim;
[m, ~]  = size(phi);

epsilon(:,1) = phi(:,1);
s(:, 1) = G*epsilon(:,1) + sigmae*randn(length(G), 1);
u(:, 1) = H'*H \ H' * Gamma * Cs *Ks*s(:,1);
x = zeros(length(As), size(phi, 2)-2);

for k = 2:size(phi, 2)-1
    epsilon(:, k) = phi(:, k) - H*u(:, k-1);
    s(:, k) = G*epsilon(:,k) + sigmae*randn(length(G), 1);
    x(:, k+1) = (As - Ks*Cs)*x(:, k) + Ks * (s(:, k) + G*H*u(:, k-1));
    u(:, k) = H'*H \ H' * Gamma*Cs*x(:,k+1);
    %u(:, k) = (H^-1) * Gamma*Cs*x(:,k+1);
end

mean_e1_real = mean(epsilon, 1); % Calculate mean for each timestep
var_e = var(epsilon - ones(m, 1)*mean_e1_real, 0, 2);
var_e = mean(var_e);

end