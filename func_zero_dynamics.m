% Simulation zero dynamics of 3 link biped in the sagittal plane
%
% Inputs:
%       t: time
%       z: cyclic variables [q1, dq1]
%       a: bezier coefficient for q2 - a(1:5), and q3 - a(6:10)
%       s_params: [z_min, z_max] - max and min angles in gait
%       
% Outputs:       
%       dz = [dq1, ddq1]
%
function dz = func_zero_dynamics(t,z,a,s_params)

z_min = s_params(1); 
z_max = s_params(2); 
delq = z_max - z_min;

q1 = z(1); dq1 = z(2);

% Computes gait timing variable s
% 
% Inputs: 
%       q1: cyclic variable
%       q1_min: minimum value for q1 during gait
%       q1_max: max value for gait
%
% Output: 
%       s: gait timing variable
%
s = func_gait_timing(q1,z_min,z_max);

alpha2 = a(1:5);
alpha3 = a(6:10);

%Need bezier formulation
function b = control_points(s, a2)
    M = 4;
    b = 0;
    for k = 0:M
        b = b + a2(1,k+1) * (factorial(M)/(factorial(k)*factorial(M-k))) * s^k * (1-s)^(M-k);
    end
end

b2 = control_points(s, alpha2);
b3 = control_points(s, alpha3);
h_desired = [q1; b2; b3];
%need to change
q_converted = eye(3)*h_desired;
%dq_converted = jacobian(h_desired, s) * (1/delq);

% find based on s and bezier definition


q2 = q_converted(2);
%q2_desired = q_converted(2);
q3 = q_converted(3);
%q3_desired = q_converted(3);

function db_ds = compute_partial_bezier(M, alpha, s)

    db_ds = 0;
    for k = 0:M-1
        db_ds = db_ds + M * ...
            (alpha(k+2) - alpha(k+1)) * ...
            (factorial(M-1)/(factorial(k)*factorial(M-1-k))) * ...
            s^k * (1-s)^(M-1-k);
    end

end

M = 4;
db_ds2 = compute_partial_bezier(M,alpha2,s) 
db_ds3 = compute_partial_bezier(M,alpha3,s)

beta2 = [db_ds2; db_ds3]/delq;

%Flagged look at 6.13 in grizzle book. 
dq2 = db_ds2 * (dq1/delq) %dq_converted(2);
dq3 = db_ds3 * (dq1/delq) %dq_converted(3);

q = [q1, q2, q3];
dq = [dq1, dq2, dq3];

% Get model paramters
[r,m,Mh,Mt,l,g] = func_model_params;
params = [r,m,Mh,Mt,l,g];

% Compute D,C,G,B matrices
% Inputs:
%       q
%       dq
%       params
%
% Output: [D,C,G,B]
%
[D,C,G,B] = func_compute_D_C_G_B(q,dq,params);


%%%%%%%%%%%%%%%%%% Apply dynamics partitioning to compute the zero dynamics
%%%%%%%%%%%%%%%%%% equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% E.g., 

%Flag
% D_inv = D \ eye(3);
% Omega = -C*dq'-G;
% fx_under = D_inv*(Omega);
% fx = [dq'; fx_under'];
% 
% gx_top = [0, 0, 0]';
% gx_under = D_inv*B;
% gx = [gx_top; gx_under];

%Flag
%H1 = eye(3);

%Flag
D1 = D(1,1);
D2 = D(1,2:3);
H1 = C(1,:)*dq' + G(1);

%------------------------------------------------------------------------%

% Get beta1 term to compute dz
%
% Inputs: s, [dq1, delq], [alpha2, alpha3]
%       s: gait timing variable
%       dq1: velocity of cyclic variable = z2
%       delq: difference between z_max and z_min in s, needed when
%               computing db/ds
%       alpha2: Bezier coefficient (1st to 5th) for q2
%       alpha3: Bezier coefficient (1st to 5th) for q3
%
% Output:
%       beta1: one of 2 matrices obtained when computing ddq_b (acceleration
%               of body coordinates) 
%       	y = h(x) = q_b - b(s), when y = 0,  q_b = b(s)
%           taking two time derivatives we get:
%       	ddq_b = d/ds(db(s)/ds*ds/dt)*ds/dt + db(s)/ds * d^2(s)/dt^2
%                 = beta1 + beta2*ddq_1
%           Note: d^2(s)/dt^2 = 1/delq*d^2(q1)/dt^2
%
%       	beta1 = d/ds(db(s)/ds*ds/dt)*ds/dt
%        	beta2 = db(s)/ds*1/delq
%
beta1 = func_compute_beta1(s, [dq1, delq], [alpha2, alpha3]);

% Compute 1st partial derivative of Bezier polynomial
% Inputs:
%       s: gait timing variable
%       M: order of Bezier polynomial, in this case 4th order
%       alpha: Coefficients for the polynomial, need M+1
%
% Outputs:
%       db/ds
%

% function db_ds = compute_partial_bezier(M, alpha, s)
% 
%     db_ds = 0;
%     for k = 0:M-1
%         db_ds = db_ds + M * ...
%             (alpha(k+2) - alpha(k+1)) * ...
%             (factorial(M-1)/(factorial(k)*factorial(M-1-k))) * ...
%             s^k * (1-s)^(M-1-k);
%     end
% 
% end
% 
% db_ds2 = compute_partial_bezier(M,alpha2,s) 
% db_ds3 = compute_partial_bezier(M,alpha3,s)
% 
% beta2 = [db_ds2; db_ds3]/delq;

%------------------------------------------------------------------------%

dz(1) = z(2);

%Flag
dz(2) = -(D2*beta1 + H1)/(D1 + D2*beta2);

dz = [dz(1), dz(2)]';
end