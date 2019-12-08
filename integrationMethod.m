%Values in simulation:
%   x1_init = 0.01
%   R1: 1.5
%   Alpha12: 1.1
%   K1 = 1.2

%   x2_init = 0.02
%   R2: 1.6
%   Alpha21: 1.4
%   K2 = 1.3

clc
clear

load X1
load X2

% hold on
% title("Analyzed Timeseries")
% scatter(X1)
% scatter(X2)

S1 = X1.data;
S2 = X2.data;
m = length(S1)-1;
n = 2;

% Solving for Species 1

d1 = zeros(m, 1);
xdash1 = zeros(m, n);

for i = (1:m)
    d1(i) = (S1(i+1) - S1(i));
    xdash1(i,1) = (S1(i+1) + S1(i))/2;
    xdash1(i,2) = ((S1(i+1))^2 + (S1(i))^2)/2;
    xdash1(i,3) = (S1(i+1)*S2(i+1) + S1(i)*S2(i))/2;
end

a1 = inv(transpose(xdash1)*xdash1)*transpose(xdash1)*d1

% Solving for Species 2

d2 = zeros(m, 1);
xdash2 = zeros(m, n);

for i = (1:m)
    d2(i) = (S2(i+1) - S2(i));
    xdash2(i,1) = (S2(i+1) + S2(i))/2;
    xdash2(i,2) = ((S2(i+1))^2 + (S2(i))^2)/2;
    xdash2(i,3) = (S1(i+1)*S2(i+1) + S1(i)*S2(i))/2;
end

a2 = inv(transpose(xdash2)*xdash2)*transpose(xdash2)*d2

%Resimulate the results [NOT Working]
tspan = [0 10];
x_init = [0.01 0.02];

[TEMP1, TEMP2] = odefun(dx1,dx2,a1,a2)
[dx1, dx2] = ode45(@(dx1, dx2) odefun(dx1,dx2,a1,a2), tspan, x_init);

 

function [dx1, dx2] = odefun(x1, x2, a1, a2)    
    dx1 = x1*(a1(1) + a1(2)*x1 +a1(3)*x2);
    dx2 = x2*(a2(1) + a2(2)*x2 +a2(3)*x1);
end


