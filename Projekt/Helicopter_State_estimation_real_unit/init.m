l_cm = 0.015;
m_heli = 0.479;
j_eq_p = 0.0172;
j_eq_y = 0.0210;
g = 9.81;
K_pp = 0.0556;
K_yy = 0.21084;
K_py = 0.005;
K_yp = 0.15;
B_p = 0.01;
B_y = 0.08;

teta_op = -10*pi/180;
psi_op = pi/2;
wteta_op = 0;
wpsi_op = 0;

V_mp_op =K_yy*m_heli*g*cos(teta_op)*l_cm/(K_yy*K_pp - K_yp*K_py);
v_my_op = K_yp*V_mp_op/K_yy;

A_c = [0, 0, 1, 0;
    0, 0, 0, 1;
    m_heli*g*l_cm*sin(teta_op)/(j_eq_p+m_heli*l_cm^2), 0, -B_p/(j_eq_p+m_heli*l_cm^2), 0;
    0, 0, 0, -B_y/(j_eq_y + m_heli*(cos(teta_op)^2)*l_cm^2)];

B_c = [0, 0;
    0, 0;
    K_pp/(j_eq_p + m_heli*l_cm^2), -K_py/(j_eq_p + m_heli*l_cm^2);
    K_yp/(j_eq_y + m_heli*(cos(teta_op)^2)*l_cm^2), -K_yy/(j_eq_y + m_heli*(cos(teta_op)^2)*l_cm^2)];

C_c = [1,0,0,0;
       0,1,0,0];

D_c = zeros(2);

%desired_poles = [-2, -3, -4, -5].*(1/5);

%L = place(A_c', C_c', desired_poles);
%L = L';

G = eye(4);
H = zeros(2,4);

plant = ss(A_c, [B_c,G], C_c, 0);

W = diag([7e1, 7e1, 8e1, 6e1]);
V = diag([0.1, 0.1]);

[kalmf, L, P, M] = kalman(plant, W, V);




