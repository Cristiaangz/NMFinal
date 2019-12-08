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
n = 3;

% Solving for timeseries X1

d = zeros(m, 1);
xdash = zeros(m, n);

for i = (1:m)
    d(i) = (S1(i+1) - S1(i));
    xdash(i,1) = (S1(i+1) + S1(i))/2;
    xdash(i,2) = ((S1(i+1))^2 + (S1(i))^2)/2;
    xdash(i,3) = (S1(i+1)*S2(i+1) + S1(i)*S2(i))/2;
end

% There is a problem when we try to get the A matrix this arises from the
% fact that the rank below is 2 instead of the expected 3 which makes it
% non-invertible.

rank(transpose(xdash)*xdash)
