clc; clear;
format short g

% ------------------------------------------------------------------------
% Inertia Coordinate Frame Convention
%
%        x (up)
%        ^
%        |
%        |
%        O ----> y (left)
%
% • Origin O is at the stance foot contact point.
% • +x axis points UPWARD (vertical direction).
% • +y axis points to the LEFT (horizontal direction).
% • Angles are measured CLOCKWISE and taken as POSITIVE.

%% Symbols
syms q1 q2 q3 p_h p_v dp_h dp_v dq1 dq2 dq3 real
syms r m Mh Mt l g real 

q = [q1; q2; q3;];
dq = [dq1; dq2; dq3;];


%% Helper transforms (2D homogeneous, 3x3)
Rot = @(th) [cos(th) -sin(th) 0;
             sin(th)  cos(th) 0;
             0        0       1];

TransX = @(a) [1 0 a;
               0 1 0;
               0 0 1];

%% --- Forward kinematics (SYMBOLIC) ---

% Stance leg: origin -> link1 start rotation, then translate along x by r
T_o_l1s = Rot(q1);
T_l1s_l1e = TransX(r);
T_o_l1e = T_o_l1s * T_l1s_l1e;              % hip end (top of stance leg)

% Torso: from hip end rotate by (pi - q3) then translate by l
T_l1e_torso_s = Rot(pi - q3);
T_torso_s_torso_e = TransX(l);
T_o_torso_e = T_o_l1e * T_l1e_torso_s * T_torso_s_torso_e;

% m1: midpoint of stance leg
T_l1s_m1 = TransX(r/2);
T_o_m1 = T_o_l1s * T_l1s_m1;

% Swing leg: from hip end rotate by (pi - q2) then translate by r/2 to m2
T_l1e_sw_s = Rot(pi - q2);
T_sw_s_m2 = TransX(r/2);
T_sw_s_m2_e = TransX(r);
T_o_m2 = T_o_l1e * T_l1e_sw_s * T_sw_s_m2;
T_o_m2_e = T_o_l1e * T_l1e_sw_s * T_sw_s_m2_e;

%% --- Extract positions (2x1) ---
pHip   = T_o_l1e(1:2, 3);
pTorso = T_o_torso_e(1:2, 3);
pM1    = T_o_m1(1:2, 3);
pM2    = T_o_m2(1:2, 3);
pM2End    = T_o_m2_e(1:2, 3);

%% --- Jacobians (2x3) ---
J_Hip   = simplify(jacobian(pHip,   q));
J_Torso = simplify(jacobian(pTorso, q));
J_M1    = simplify(jacobian(pM1,    q));
J_M2    = simplify(jacobian(pM2,    q));
J_M2End    = simplify(jacobian(pM2End,    q));

%% Display symbolic results
disp("pHip = ");   disp(pHip);
disp("pTorso = "); disp(pTorso);
disp("pM1 = ");    disp(pM1);
disp("pM2 = ");    disp(pM2);

disp("J_Hip = ");   disp(J_Hip);
disp("J_Torso = "); disp(J_Torso);
disp("J_M1 = ");    disp(J_M1);
disp("J_M2 = ");    disp(J_M2);

%% --- Optional: evaluate at your test values ---
r_val  = 1.0;
l_val  = 0.5;
q1_val = pi/6;
q2_val = 2*pi/3;
q3_val = pi/4;

subs_list = [r, l, q1, q2, q3];
vals_list = [r_val, l_val, q1_val, q2_val, q3_val];

pHip_num   = double(subs(pHip,   subs_list, vals_list));
pTorso_num = double(subs(pTorso, subs_list, vals_list));
pM1_num    = double(subs(pM1,    subs_list, vals_list));
pM2_num    = double(subs(pM2,    subs_list, vals_list));

J_Hip_num   = double(subs(J_Hip,   subs_list, vals_list));
J_Torso_num = double(subs(J_Torso, subs_list, vals_list)); 
J_M1_num    = double(subs(J_M1,    subs_list, vals_list));
J_M2_num    = double(subs(J_M2,    subs_list, vals_list));

disp("Numeric pHip, pTorso, pM1, pM2:");
disp(pHip_num); disp(pTorso_num); disp(pM1_num); disp(pM2_num);

disp("Numeric Jacobians:");
disp(J_Hip_num); disp(J_Torso_num); disp(J_M1_num); disp(J_M2_num);

dq1_val = 0.5;
dq2_val = 0.3;
dq3_val = 0.2;

dpHip = J_Hip_num * [dq1_val; dq2_val; dq3_val];
dpTorso = J_Torso_num * [dq1_val; dq2_val; dq3_val];
dpM1 = J_M1_num * [dq1_val; dq2_val; dq3_val];
dpM2 = J_M2_num * [dq1_val; dq2_val; dq3_val];

disp("Velocities in cartesian space:");
disp(dpHip); disp(dpTorso); disp(dpM1); disp(dpM2);

% masses 
m_val = 1.0; 
Mh_val = 2.0; 
Mt_val = 3.0;
M_tot = Mh_val + Mt_val + 2*m_val;
pCoM = (Mh_val*pHip + Mt_val*pTorso + m_val*pM1 + m_val*pM2)/M_tot;
J_PCoM    = simplify(jacobian(pCoM,    q));
disp("J_PCoM = ");    disp(J_PCoM);
J_PCoM_num    = double(subs(J_PCoM,    subs_list, vals_list));
pCoM_num    = double(subs(pCoM,    subs_list, vals_list));
dpCoM = J_PCoM_num * [dq1_val; dq2_val; dq3_val];
disp("CoM pose:");
disp(pCoM_num)
disp("CoM velocity:");
disp(dpCoM)

% kinetic energy equations

D = ((m*J_M1'*J_M1)+(m*J_M2'*J_M2)+(Mt*J_Torso'*J_Torso)+(J_Hip'*Mh*J_Hip));
disp(D)
KE = 0.5*(dq')*(D)*dq;

PE = g *((Mt*pTorso(1))+(Mh*pHip(1))+(m*pM1(1))+(m*pM2(1)));

coriolis_matrix(D, q, dq);

%TEST CODE BELOW FOR THE 3 FOR LOOPS


function C = coriolis_matrix(D, q, qdot)
% CORIOLIS_MATRIX Computes the Coriolis and centrifugal matrix C(q, qdot)
% using Christoffel symbols of the first kind.
%
% Inputs:
%   D     - n x n symbolic inertia/mass matrix (function of q)
%   q     - n x 1 symbolic vector of generalized coordinates
%   qdot  - n x 1 symbolic vector of generalized velocities
%
% Output:
%   C     - n x n symbolic Coriolis matrix

n = length(q);  % number of generalized coordinates
C = sym(zeros(n, n));

C_logic = sym(zeros(n,n));
tic
for k = 1:n
    for j = 1:n
        c_kj = sym(0);
        for i = 1:n
            % Christoffel symbol: (1/2)*(dd_kj/dq_i + dd_ki/dq_j - dd_ij/dq_k)
            christoffel = (1/2) * ( diff(D(k,j), q(i)) ...
                                  + diff(D(k,i), q(j)) ...
                                  - diff(D(i,j), q(k)) );
            c_kj = c_kj + christoffel * qdot(i);
        end
        C(k,j) = simplify(c_kj);
    end
end
elapsed_time = toc;
% disp("Brute force elapsed time")
% disp(elapsed_time)
% disp("C_bruteforce")
% disp(C)

global dict
dict = containers.Map('KeyType', 'char', 'ValueType', 'double');
% tic
for k = 1:n
    for j = 1:n
        c_kj = sym(0);
        for i = 1:n
            christoffel = sym(0);
            if i ~= j
                key = [k j i];
                key_str = mat2str(key);
                if isKey(dict, key_str)
                    christoffel = dict(key);
                
                else
                    % Christoffel symbol: (1/2)*(dd_kj/dq_i + dd_ki/dq_j - dd_ij/dq_k)
                    christoffel = (1/2) * ( diff(D(k,j), q(i)) ...
                                  + diff(D(k,i), q(j)) ...
                                  - diff(D(i,j), q(k)) );
                end
            else
                christoffel = (1/2) * ( diff(D(k,j), q(i)) ...
                                  + diff(D(k,i), q(j)) ...
                                  - diff(D(i,j), q(k)) );
            end
            c_kj = c_kj + christoffel * qdot(i);
        end
        C_logic(k,j) = simplify(c_kj);
    end
end
% elapsed_time = toc;
% disp("Logic elapsed time")
% disp(elapsed_time)
% disp("C_logic: ")
% disp(C_logic)

end

%%%%%%%%%%%%%%%%%% What is control input matrix? %%%%%%%%%%%%%%%%%%%%%%%%%%
B = [0, 0; ...
     1, 0; ...
     0, 1];

% Write 3 link model to file
% Inputs:
%       q
%       dq
%       params
%
% Outputs: 
%       D: Inertia matrix
%       C: Coriolis matrix
%       G: Gravity matrix
%       B: 
%
write_symbolic_term_to_mfile(q,dq,params,D,C,G,B)

%-------------------------------------------------------------------------%
%%%% Impact map

% Using same psotion vectors as above, but taking partial with respect to qe
% instead

% Extended configuration variables
p_e = [p_h; p_v];

qe = [q; p_h; p_v];
dqe = [dq; dp_h; dp_v];

% Extended position
pMh_e = pMh + p_e;%bruh what p_e needs to
pMt_e = pMt + p_e;
pm1_e = pm1 + p_e;
pm2_e = pm2 + p_e;
P2e = P2 + p_e;

% Extended velocities
J_mh_e   = simplify(jacobian(pMh_e,   qe));
J_mt_e   = simplify(jacobian(pMt_e, qe));
J_m1_e   = simplify(jacobian(pm1_e,    qe));
J_m2_e   = simplify(jacobian(pm2_e,    qe));

vMh_e = J_mh_e * dqe;
vMt_e = J_mt_e * dqe;
vm1_e = J_m1_e * dqe;
vm2_e = J_m2_e * dqe;

K_Mhe = 0.5*(dqe')*(Mh*J_mh_e'*J_mh_e)*dqe;
K_Mte = 0.5*(dqe')*(Mt*J_mt_e'*J_mt_e)*dqe;
K_m1e = 0.5*(dqe')*(Mh*J_m1_e'*J_mh_e)*dqe;
K_m2e = 0.5*(dqe')*(Mh*J_m2_e'*J_mh_e)*dqe;

Ke = K_m1e + K_Mhe + K_Mte + K_m2e;

% Extended inertia matrix
De = (Mh*(J_mh_e)'*J_mh_e)+(Mt*(J_mt_e)'*J_mt_e)+(Mh*J_m1_e'*J_mh_e)+(Mh*J_m2_e'*J_mh_e);

E = simplify(jacobian(P2e,    qe));

% Partial of any point on biped, hip chosen in this case
dY_dq = jacobian(pMh_e,q);

% Write impact map to a file
% Inputs:
%       q
%       dq
%       params
%
% Outputs: Matrices needed to compile impact map
%       De: Extended inertia matrix
%       E:
%       dY_dq:
%       
write_symbolic_term_to_mfile(q,dq,params,De,E,dY_dq)


%-------------------------------------------------------------------------%
%%%% For controller

% Vector fields
fx = 0;
gx = 0;

% Bezier poly - needed for output function
syms s delq
%s = (q1 - q1_plus)/delq; 
%delq = q1_minus - q1_plus 
%ds/dt = dq1/delq; ds/dq1 = 1/delq;

syms a21 a22 a23 a24 a25 
syms a31 a32 a33 a34 a35

a2 = [a21 a22 a23 a24 a25];
a3 = [a31 a32 a33 a34 a35];
M = 4;

b2 = 0; b3 = 0;
for k = 0:M
    b2 = b2 + a2(1,k+1)*(factorial(M)/(factorial(k)*factorial(M-k)))*s^k*(1-s)^(M-k);
end

for k = 0:M
    b3 = b3 + a3(1,k+1)*(factorial(M)/(factorial(k)*factorial(M-k)))*s^k*(1-s)^(M-k);
end

% Defining outputs

h = [q2 - b2; q3 - b3];

% y_dot = Lfh = dh/dx*fx - independent of gx*u since relative degree is 2
% However, h is a function of (s,q2,q3), not q1 directly, so the following
% is used:
% dh/dq1 = dh/ds*ds/dq1 = dh/ds*1/delq
%
% Temporary variable that multiples the 1st column with 1/delq
temp = sym(eye(6)); temp(1)  = 1/delq;

Lfh = jacobian(h,[s; q2; q3; dq])*temp*fx;

dLfh = jacobian(Lfh,[s; q2; q3; dq])*temp;

% Write matrix used in feedback linearization - d/dx(Lfh) to file
% Inputs:
%       s = (q1 - q1_plus)/delq: gait timing variable
%       delq = q1_minus - q1_plus: difference in cyclic variable during gait 
%       dq1
%       params: 
%       a2: bezier coefficents (1st - 5th) for q2
%       a3: bezier coefficents (1st - 5th) for q3
%
% Outputs:
%       dLfh: partial of Lfh, to be used to compute L2fh and LgLfh
%
write_symbolic_term_to_mfile([s,delq],dq1,[a2,a3],dLfh);



%-------------------------------------------------------------------------%
%%%% For Zero Dynamics


db_ds2 = 0;
for k = 0:M-1
    db_ds2 = db_ds2 + (a2(1,k+2)-a2(1,k+1))*(factorial(M)/(factorial(k)*factorial(M-k-1)))*s^k*(1-s)^(M-k-1);
end

db_ds3 = 0;
for k = 0:M-1
    db_ds3 = db_ds3 + (a3(1,k+2)-a3(1,k+1))*(factorial(M)/(factorial(k)*factorial(M-k-1)))*s^k*(1-s)^(M-k-1);
end

partial_db_ds2 = jacobian(db_ds2,s)*dq1/delq;

partial_db_ds3 = jacobian(db_ds3,s)*dq1/delq;

beta1 = [partial_db_ds2; partial_db_ds3]*dq1/delq;

eta2 = jacobian(K,dq1);

write_symbolic_term_to_mfile(s,[dq1, delq],[a2, a3],beta1)

write_symbolic_term_to_mfile(q,dq,params,eta2)

