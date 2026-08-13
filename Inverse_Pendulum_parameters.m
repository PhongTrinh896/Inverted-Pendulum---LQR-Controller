% ======== He so cho mo hinh con lac/ Inverted pendulum parameters =============
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
% ============== Ma tran he so trang thai tuyen tinh hoa / Linearized state matrices =============
A = [ 0         0       1       0;
      0         0       0       1;
      0         ((m^2)*l^2*g)/alp       -(J+m*l^2)*(c+kt*kb/(Rm*r^2))/alp       -b*l*m/alp;
      0         (M+m)*m*g*l/alp         -l*m*(c + kt*kb/(Rm*r^2))/alp           -(M+m)*b/alp];
B = [0;         0;          (J+m*l^2)*kt/(alp*Rm*r);        l*m*kt/(alp*Rm*r)];
C = [1 0 0 0; 0 1 0 0];
% ========= Thiet ke bo dieu khien LQR / LQR controller design ============
Q = [1000   0   0   0;
    0       1500    0   0;
    0       0       100   0;
    0       0       0   1];
R = 1;
P = care(A, B, Q, R); 
K = lqr(A,B,Q,R);

% ========= Thiet ke bo loc Kalman / Kalman filter design ============
G = eye(4);
Qn = 0.0001*G;
Rn = [0.0001  0 ; 0  0.0001];
Kk = lqe(A, G, C, Qn, Rn);