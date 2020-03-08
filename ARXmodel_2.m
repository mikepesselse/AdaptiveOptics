%% MATLAB Exericse II: Filtering and Identification
% For the course SC42025, Filtering and Identification,
% By Mike Pesselse [4300564],
% On 18-12-2019.

%% Initialisation of data and parameters
clc; clearvars; close all;
load('Homework.mat')

%% 1) ARX Model
% See below for the ARX function 

%% 2) subspace identification Model
% See below for the subspace identification function 

%% 3)Identification of a 6th order model
% In question 3 a 6th order model is identified using four couples of
% inputs an outputs of an unknown dynamical system: (u1, y1), (u2, y2),
% (u3, y3) and (u4, y4)
%% 3.1) Which sets are useful for identification?
% The 2 input sets that are useful for identification should be 
% persistently exciting. According to definition 10.1 on page 358 of the
% textbook, the Hankel matrix of the input set should have a rank equal to 
% the order of the model. So in this case we should check which of the
% four input sets have a rank of 6. This seems to be input set 1 and input
% set 4. 

% Define order of the model 
n = 6; 

% Create Hankel matrices for the four input singal
u1_H = hankel(u1(1:n), u1(n:end));
u2_H = hankel(u2(1:n), u2(n:end));
u3_H = hankel(u3(1:n), u3(n:end));
u4_H = hankel(u4(1:n), u4(n:end));

% Determine the rank of the Hankel input matrices 
rank_u1_H = rank(u1_H);
rank_u2_H = rank(u2_H);
rank_u3_H = rank(u3_H);
rank_u4_H = rank(u4_H);

% Display rank of matrices
disp(['rank of Hankel input matrix 1 is ', num2str(rank_u1_H)])
disp(['rank of Hankel input matrix 2 is ', num2str(rank_u2_H)])
disp(['rank of Hankel input matrix 3 is ', num2str(rank_u3_H)])
disp(['rank of Hankel input matrix 4 is ', num2str(rank_u4_H)])

%% 3.2) How to define s in mysubid?
% s should be chosen such that it is larger than n, According to page 295
% of the textbook. (9.2). Take for example s=2n.

%% 3.3) Why is S an useful output in mysubid? 
% S are the singular values of Y0,s,N Π⊥U0,s,N
% The singular values (S) give an indication on how to choose the order of 
% the model. The number of non-zero singular values is the order of
% the system.


%% 3.4) Identifty 6th order model of the system
% Using both created mysubid and myarx functions we identify the 6th order 
% model.

% Define s for mysubid
s = 30; % 

% Use mysubid function for subspace identification of the useful sets
% Set 1
[At_1, Bt_1, Ct_1, Dt_1, x0_1, S1] = mysubid(y1, u1, s, n); 
sys1 = ss(At_1, Bt_1, Ct_1, Dt_1, []);

% Set 4
[At_4, Bt_4, Ct_4, Dt_4, x0_4, S4] = mysubid(y4, u4, s, n); 
sys4 = ss(At_4, Bt_4, Ct_4, Dt_4, []);

% use myarx function for identification of ARX model of the useful sets
% Set 1
[aest_1, best_1] = myarx(y1, u1, n); 
tf_arx_1 = tf(best_1', aest_1', []);

% Set 4
[aest_4, best_4] = myarx(y4, u4, n); 
tf_arx_4 = tf(best_4', aest_4', []);

%% 3.5) Bode plots of the 6th order transfer functions
% When comparing the bode plots of the 6th order transfer function
% identifications of question 4, we can see that the mysubid model for set
% 4 tracks the original system the best. Overall mysubid is preforming
% better than myarx. 

% Also, increasing s seems to improve the mysubid models. 
bode_options = bodeoptions;
bode_options.Grid = 'on';
bode(tfse, sys1, sys4, tf_arx_1, tf_arx_4, bode_options)
legend('System','SubId 1', 'SubId 4', 'ARX 1', 'ARX 4')
title('3.5) Bode plots of the 6th order transfer functions')

%% 4)Identification of a 6th order model using vector u0 and y0
% In question 4 the 6th order model is identified using the vectors u0 and
% y0 of the file homework.mat 
%% 4.1) Estimate order of y0
% The singular values obtained by choosing two different values of n and s 
% show that the order of the system is 2, due to the large gap between 
% the other singular values. 

% Singular Values for s = 20 
[At0, Bt0, Ct0, Dt0, x0t0, S0] = mysubid(y0, u0, 20, 2);
figure; 
semilogy(diag(S0),'ob');
title('4.1) Estimate order of y0')
legend('Singular Values')
 
%% 4.2) Suppose a bias error of 0.2
% In this case it can be seen from the plotted singular values that there
% are relevant differences with what was found at 4.1. We can see there are
% three singular values and then a big gap to the others. This means that
% y00 has order 3. 

% Define y00
y00 = y0 + 0.2;

% Singular Values for s = 20
[At00, Bt00, Ct00, Dt00, x0t00, S00] = mysubid(y00, u0, 20, 3);
figure; 
semilogy(diag(S00),'or');
title('4.2) Estimate order of y00')
legend('Singular Values')

%% 4.3) Is the system found at 4.2 reachable?
% The system found at 4.2 is not reachable since we have an uncontrollable
% mode and The poles can be seen below. 

% Poles of At00
Poles_At00 = eig(At00);
Poles_At00

% Putting the ss in controllability form shows that we have an
% uncontrollable at 1, thus not reachable
[Abar,Bbar,Cbar,T,k]= ctrbf(At00, Bt00, Ct00);
Abar

%% 1) ARX Model

function [aest, best] = myarx(y, u, n)
%MYARX: identification of a SISO discrete-time ARX model
%   The function will return the numerator (best) and denominator (aest) 
%   coefficients of the discrete- time SISO ARX model from a given input 
%   (column) vector u and its corresponding output vector y. 
%   The user should also specify the desired order (n) for the model.
%   Mike Pesselse - 4300564

N = length(y); 

% Pre determine size of phi
% phi = zeros(2*n,N-n);

for i = n+1:N                   % timestep k 
    j = (i-1:-1:i-n);           % inputs for phi
    phi(:,i-n) = [y(j);u(j)];   % create phi
end

% Determine Y for useful data points y
Y_N = y(n+1:N);

% Optimal Theta
Theta = inv((1/N)*phi*(phi'))*(1/N)*phi*Y_N;

% Create coefficients of the denominator aest and numerator best 
aest = [1; Theta(1:n)];
best = [0; Theta(n+1:end)];

end      

%% 2) subspace identification Model

function [At, Bt, Ct, Dt, x0t, S] = mysubid(y, u, s, n)
%MYSUBID:  subspace identification of a SISO state-space model
%   The function will return the state-space matrices (At, Bt, Ct and Dt) 
%   and the initial state x0t from a given input (column) vector u and
%   its corresponding output vector y. The user should also specify the 
%   desired order (n) for the model and the parameter s 
%   (the number of rows of the block-Hankel matrices Y0,s,N and U0,s,N 
%   used in subspace identification). An additional output (S) is a vector 
%   containing the singular values of Y0,s,N Π⊥U0,s,N .
%   Mike Pesselse - 4300564


% Build Y and U matrices
N = length(y)-s+1;
Yh = hankel(y(1:s), y(s:N+s-1));
Uh = hankel(u(1:s), u(s:N+s-1));

% Calculate projection matrix
pie = eye(N) - (Uh'*(inv(Uh*(Uh')))*Uh);

% Calculate SVD so we get U and S matrix
[U,S,~] = svd(Yh*pie);

% Select correct terms from U
Un = U(:,1:n);

% At and Ct output
Ct = Un(1,:);
At = Un(1:s-1,:)\Un(2:s,:);

% Bt, Dt and X0t output
[Bt,Dt,x0t] = subidhelp(y,u,At,Ct);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I first build mysubid using QR factorization but using page 306 of the
% textbook. However this does not seems to work or produce the same result
% as the method from above. I still dont understand why, nobody was able to
% explain me. 

% % Build Y and U matrices
% [r,c] = size(y);
% Y = y(1:r-s+1,:)'; 
% U = u(1:r-s+1,:)';
% 
% for i=2:s
%     Y = [Y;y(i:r-s+i)']; 
%     U = [U;u(i:r-s+i)'];
% end
% 
% % QR factorization on U and Y
% [~,R] = qr([U;Y]);
% 
% % Select R_22 out of R
% R_T = R';
% R_22 = R_T(s+1:end, s+1:2*s);
% 
% % SVD of R_22 as is done on page 306 of the textbook (9.2.4.3)
% [U_svd, S_svd, ~] = svd(R_22);
% 
% At = U_svd(1:c*(s-1), 1:n) \ U_svd((c+1):(s*c), 1:n);
% Ct = U_svd(1:c,1:n);
% S = diag(S_svd);
% [Bt,Dt,x0t] = subidhelp(y,u,At,Ct);

end
