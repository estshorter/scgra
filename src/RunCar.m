rehash; clear all; close all;

% ‰Šú„’è‰ğ‚Ì—pˆÓ
N = 50;

x_est = zeros(2,N+1);
u_est = zeros(2,N+1);

pi_est = 2*sqrt(2);
umax = 0.5;
tmp = 0;
t = 0:1/N:1;
for ts=1:N/2
    t = (ts-1)/N;
    u_est(1,ts) = umax;
    u_est(2,ts) = asin((u_est(1,ts))/1);%sqrt(1 - u_est(1,ts)^2);
    x_est(2,ts) = pi_est * umax * t;
    x_est(1,ts) = pi_est^2 * 0.5 * umax * t^2;
end

for ts=N/2+1:N+1
    t = (ts-1)/N;
    u_est(1,ts) = - umax;
    u_est(2,ts) = asin((u_est(1,ts))/1);%sqrt(1 - u_est(1,ts)^2);
    x_est(2,ts) = pi_est * umax * (1-t);
    x_est(1,ts) = umax * pi_est^2 * (t- 0.5*t^2) - umax/4*pi_est^2;
end
% plot(x_est(1,:))

% load('x_car.mat')
% load('u_car.mat')
% load('pi_car.mat')
% x_est = x;
% u_est = u;
% pi_est = pi;

% 
self = Car();
self.Init();
self.SetInitialGuess(x_est,u_est,pi_est);
self.Calc()