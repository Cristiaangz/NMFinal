%Values in simulation:
%   R1: 1.5
%   Alpha12: 1.1
%   K1 = 1.2

%   R2: 1.6
%   Alpha21: 1.4
%   K2 = 1.3

clc
clear

load X1
load X2

% hold on
% title("Analyzed Timeseries")
% plot(X1)
% plot(X2)

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
 
