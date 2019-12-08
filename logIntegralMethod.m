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
xdash1 = ones(m, n);

for i = (1:m)
    d1(i) = (log(S1(i+1)) - log(S1(i)));
    xdash1(i,2) = (S1(i+1) + S1(i))/2;
    xdash1(i,3) = (S2(i+1) + S2(i))/2;
end


a1 = inv(transpose(xdash1)*xdash1)*transpose(xdash1)*eye(m,1)

% Solving for Species 2

d2 = zeros(m, 1);
xdash2 = ones(m, n);

for i = (1:m)
    d2(i) = (log(S2(i+1)) - log(S2(i)));
    xdash2(i,2) = (S2(i+1) + S2(i))/2;
    xdash2(i,3) = (S1(i+1) + S1(i))/2;
end


a2 = inv(transpose(xdash2)*xdash2)*transpose(xdash2)*eye(m,1)
