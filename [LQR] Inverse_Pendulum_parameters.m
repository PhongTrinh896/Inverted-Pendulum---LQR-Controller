m = 0.1;
M = 0.136;
c = 0.63;
b = 0.00007892;
l = 0.2;
J = 0.0007176;
g = 9.81;
kb = 0.031;
kt = 0.031;
Rm = 12.5;
r = 0.006;
pi = 3.141592654;
k_swing = 12;
alp = J*(M+m) + M*m*l^2;
% Ma tran he so trang thai
A = [ 0         0       1       0;
      0         0       0       1;
      0         ((m^2)*l^2*g)/alp       -(J+m*l^2)*(c+kt*kb/(Rm*r^2))/alp       -b*l*m/alp;
      0         (M+m)*m*g*l/alp         -l*m*(c + kt*kb/(Rm*r^2))/alp           -(M+m)*b/alp];
B = [0;         0;          (J+m*l^2)*kt/(alp*Rm*r);        l*m*kt/(alp*Rm*r)];
% Ma tran he so Q, R
Q = [5000   0   0   0;
    0       1500    0   0;
    0       0       100   0;
    0       0       0   1];
R = 1;
K = lqr(A,B,Q,R)
