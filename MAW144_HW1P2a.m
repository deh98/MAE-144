clear; close all; clc;

%% MAE144 HW1 Problem2
a = RR_poly([-1 1 -3 3 -6 6],1);
b = RR_poly([-2 2 -5 5],1);
f = RR_poly([-1,-1,-3,-3,-6,-6,-20,-20,-20,-20,-20,-20],1);
[x,y] = RR_diophantine(a,b,f)
test = trim(a*x + b*y)
residual1 = norm(f - test)