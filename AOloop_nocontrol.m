function [ var_eps ] = AOloop_nocontrol(phik,sigmae,H,G)
% Example of online AO simulation for open_loop measurements
% IN
% phik  : incoming turbulence wavefront
% sigmae: measurement noise parameter
% H     : influence matrix mapping the wavefront on the mirror
% G     : measurement matrix 
% OUT
% var_eps : mean variance of the residual wavefront


n = size(H,1);      % dimension lifted wavefront
ns = size(G,1);     % dimension lifted sensor slopes
T = length(phik);   % number of temporal phase points

epsk = zeros(n,T);  % residual wavefront
eps_piston_removed = zeros(n,T); % residual wavefront with mean removed
sk = zeros(ns,T);   % slopes measurements
var_eps = zeros(T,1);

for k = 1:T-1
    epsk(:,k+1) = phik(:,k+1);
    eps_piston_removed(:,k+1) = epsk(:,k+1)-mean(epsk(:,k+1)); 
    sk(:,k+1) = G*epsk(:,k+1) + sigmae*randn(ns,1);
    var_eps(k+1) = var(eps_piston_removed(:,k+1));
end
var_eps = mean(var_eps);

end



