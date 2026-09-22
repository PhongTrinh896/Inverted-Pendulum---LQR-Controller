% ======== He so cho mo hinh con lac / Inverted pendulum parameters (DISCRETE) ========
m = 0.04467; % p.A.L
M = 0.1;
c = 0.63;
b = 0.00007892;
l = 0.1; % Chieu dai nua con lac
% Duong kinh thanh: d = 6 mm
% KLR: p = 7900

J = 0.000596;
g = 9.81;
kb = 0.031;
kt = 0.031;
Rm = 12.5;
r = 0.006;
pi = 3.141592654;
k_swing = 42;
alp = J*(M+m) + M*m*l^2;

% ============== Ma tran he so trang thai tuyen tinh hoa (LIEN TUC) / Continuous-time state matrices =============
A = [ 0         0       1       0;
      0         0       0       1;
      0         ((m^2)*l^2*g)/alp       -(J+m*l^2)*(c+kt*kb/(Rm*r^2))/alp       -b*l*m/alp;
      0         (M+m)*m*g*l/alp         -l*m*(c + kt*kb/(Rm*r^2))/alp           -(M+m)*b/alp];
B = [0;         0;          (J+m*l^2)*kt/(alp*Rm*r);        l*m*kt/(alp*Rm*r)];
C = [1 0 0 0; 0 1 0 0];
D = [0; 0];

% =====================================================================
% ============ RỜI RẠC HÓA / DISCRETIZATION (ZOH) =====================
% =====================================================================
% Ts = chu ky lay mau tren STM32 (giay). Chinh lai cho phu hop voi
% toc do doc encoder / cap nhat PWM thuc te cua ban.
Ts = 0.005;   % 5 ms -> 200 Hz (gia tri de xuat, co the tinh chinh)

sys_c = ss(A, B, C, D);
sys_d = c2d(sys_c, Ts, 'zoh');   % Zero-Order Hold discretization

Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;

fprintf('--- He roi rac (Ts = %.4f s) ---\n', Ts);
disp('Ad ='); disp(Ad);
disp('Bd ='); disp(Bd);

% ========= Thiet ke bo dieu khien LQR RỜI RẠC / Discrete LQR ============
Q = [5000   0   0   0;
    0       10    0   0;
    0       0       1   0;
    0       0       0   1];
R = 1;

Kd = dlqr(Ad, Bd, Q, R);   % Ma tran hoi tiep trang thai roi rac
fprintf('--- Discrete LQR gain Kd ---\n');
disp(Kd);

% (Tuy chon) kiem tra on dinh vong kin roi rac
eig_cl = eig(Ad - Bd*Kd);
fprintf('--- Cuc vong kin (phai nam trong duong tron don vi) ---\n');
disp(eig_cl);
disp('Do lon cuc dai:'); disp(max(abs(eig_cl)));

% ========= Thiet ke bo loc Kalman RỜI RẠC / Discrete Kalman filter ============
G  = eye(4);
Qn = 0.0001*G;          % Nhieu qua trinh (process noise covariance)
Rn = [0.0001 0; 0 0.0001]; % Nhieu do (measurement noise covariance)

% dlqe tra ve gain Kalman roi rac tuong ung voi lqe ban dau
Kkd = dlqe(Ad, G, Cd, Qn, Rn);
fprintf('--- Discrete Kalman gain Kkd ---\n');
disp(Kkd);

