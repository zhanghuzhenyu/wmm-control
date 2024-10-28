function [sys,x0,str,ts] = spacemodel(t,x,u,flag)
switch flag,
case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;
case 1,
    sys=mdlDerivatives(t,x,u);
case 3,
    sys=mdlOutputs(t,x,u);
case {2,4,9}
    sys=[];
otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end
function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 13;  % X dX [the_R;the_L;phy;th1;th2]   [dthe_R;dthe_L;dphy;dth1;dth2]
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 13 + 10;  % 4+4 + 5个角度 + q(5) + dq(5)
sizes.NumInputs      = 4;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);

qr = 0;
ql = 0;
phai = 0; % -pi/10
q1 = pi/3 ; % pi/2 + pi/6
q2 = -pi/3; %  -pi/2 - pi/6

l1 = 0.514;
l2 = 0.362; 
d = 0.116; % 0.116

global X_init dX_init
X_init =  [d*cos(phai)+l1*cos(phai+q1)+l2*cos(phai+q1+q2);d*sin(phai)+l1*sin(phai+q1)+l2*sin(phai+q1+q2); d*cos(phai);  d*sin(phai)];
dX_init = [0; 0; 0; 0];
global E_bar C_bar M_bar 
E_bar = 0*eye(4);
C_bar = 0*ones(4);
M_bar = 0*ones(4);

x0 = [X_init(1);X_init(2);X_init(3);X_init(4); dX_init(1);dX_init(2);dX_init(3);dX_init(4);
    qr;ql;phai;q1;q2
    ]; 
global velocity_to_next_loop
velocity_to_next_loop = [0;0;0;0;0];
global q dq
q = [0;0;0;0;0];
dq = [0;0;0;0;0];
str = [];
ts  = [0 0];

function sys=mdlDerivatives(t,x,u)
X = x(1:4);
dX = x(5:8);
qr=x(9);ql=x(10); phai=x(11); q1=x(12);q2=x(13);

global velocity_to_next_loop
qrdot = velocity_to_next_loop(1);
qldot = velocity_to_next_loop(2);
phai_dot = velocity_to_next_loop(3);
q1dot = velocity_to_next_loop(4);
q2dot = velocity_to_next_loop(5);
global q dq
q = [mod(qr,2*pi);mod(ql,2*pi);phai;q1;q2];
dq = [qrdot;qldot;phai_dot;q1dot;q2dot];
utol=[u(1);u(2);u(3);u(4)];

mp=17.25;   % 17.25
m1=2.56; m2=1.07; mw=0.159; 
l1 = 0.514; l11 = 0.252;
l2 = 0.362; l22 = 0.243;
b = 0.182;
r = 0.0508; 
d = 0.116;  % 0.116
Iphai = 0.297;
I1 = 0.148;
I2 = 0.0228; 
Iw = 0.0002; % Iw = 0.0002;

% 动力学矩阵
if 1>0
M11 = m1 + m2 + mp + 2*mw;
M12 = 0;
M13 = d*mp*sin(phai) + 2*d*mw*sin(phai) - l22*m2*sin(phai + q1 + q2) - l1*m2*sin(phai + q1) - l11*m1*sin(phai + q1);
M14 = - (m2*(2*l1*sin(phai + q1) + 2*l22*sin(phai + q1 + q2)))/2 - l11*m1*sin(phai + q1);
M15 = -l22*m2*sin(phai + q1 + q2);
M16 = 0;
M17 = 0;
M21 = 0;
M22 = m1 + m2 + mp + 2*mw;
M23 = l22*m2*cos(phai + q1 + q2) - 2*d*mw*cos(phai) - d*mp*cos(phai) + l1*m2*cos(phai + q1) + l11*m1*cos(phai + q1);
M24 = (m2*(2*l1*cos(phai + q1) + 2*l22*cos(phai + q1 + q2)))/2 + l11*m1*cos(phai + q1);
M25 = l22*m2*cos(phai + q1 + q2);
M26 = 0;
M27 = 0;
M31 = d*mp*sin(phai) + 2*d*mw*sin(phai) - l22*m2*sin(phai + q1 + q2) - l1*m2*sin(phai + q1) - l11*m1*sin(phai + q1);
M32 = l22*m2*cos(phai + q1 + q2) - 2*d*mw*cos(phai) - d*mp*cos(phai) + l1*m2*cos(phai + q1) + l11*m1*cos(phai + q1);
M33 = I1 + I2 + Iphai + 2*b^2*mw + d^2*mp + 2*d^2*mw + l1^2*m2 + l11^2*m1 + l22^2*m2 + 2*l1*l22*m2*cos(q2);
M34 = I1 + I2 + l1^2*m2 + l11^2*m1 + l22^2*m2 + 2*l1*l22*m2*cos(q2);
M35 = I2 + l22^2*m2 + l1*l22*m2*cos(q2);
M36 = 0;
M37 = 0;
M41 = - (m2*(2*l1*sin(phai + q1) + 2*l22*sin(phai + q1 + q2)))/2 - l11*m1*sin(phai + q1);
M42 = (m2*(2*l1*cos(phai + q1) + 2*l22*cos(phai + q1 + q2)))/2 + l11*m1*cos(phai + q1);
M43 = I1 + I2 + l1^2*m2 + l11^2*m1 + l22^2*m2 + 2*l1*l22*m2*cos(q2);
M44 = I1 + I2 + l1^2*m2 + l11^2*m1 + l22^2*m2 + 2*l1*l22*m2*cos(q2);
M45 = I2 + l22^2*m2 + l1*l22*m2*cos(q2);
M46 = 0;
M47 = 0;
M51 = -l22*m2*sin(phai + q1 + q2);
M52 = l22*m2*cos(phai + q1 + q2);
M53 = I2 + l22^2*m2 + l1*l22*m2*cos(q2);
M54 = I2 + l22^2*m2 + l1*l22*m2*cos(q2);
M55 = I2 + l22^2*m2;
M56 = 0;
M57 = 0;
M61 = 0;
M62 = 0;
M63 = 0;
M64 = 0;
M65 = 0;
M66 = Iw;
M67 = 0;
M71 = 0;
M72 = 0;
M73 = 0;
M74 = 0;
M75 = 0;
M76 = 0;
M77 = Iw;

M = [M11 M12 M13 M14 M15 M16 M17;
     M21 M22 M23 M24 M25 M26 M27;
     M31 M32 M33 M34 M35 M36 M37;
     M41 M42 M43 M44 M45 M46 M47;
     M51 M52 M53 M54 M55 M56 M57;
     M61 M62 M63 M64 M65 M66 M67;
     M71 M72 M73 M74 M75 M76 M77;
     ];

C11 = 0;
C12 = 0;
C13 = d*mp*phai_dot*cos(phai) - l11*m1*phai_dot*cos(phai + q1) - l1*m2*q1dot*cos(phai + q1) - l11*m1*q1dot*cos(phai + q1) - l1*m2*phai_dot*cos(phai + q1) + 2*d*mw*phai_dot*cos(phai) - l22*m2*phai_dot*cos(phai + q1 + q2) - l22*m2*q1dot*cos(phai + q1 + q2) - l22*m2*q2dot*cos(phai + q1 + q2);
C14 = - l1*m2*phai_dot*cos(phai + q1) - l11*m1*phai_dot*cos(phai + q1) - l1*m2*q1dot*cos(phai + q1) - l11*m1*q1dot*cos(phai + q1) - l22*m2*phai_dot*cos(phai + q1 + q2) - l22*m2*q1dot*cos(phai + q1 + q2) - l22*m2*q2dot*cos(phai + q1 + q2);
C15 = -l22*m2*cos(phai + q1 + q2)*(phai_dot + q1dot + q2dot);
C16 = 0;
C17 = 0;
C21 = 0;
C22 = 0;
C23 = d*mp*phai_dot*sin(phai) - l11*m1*phai_dot*sin(phai + q1) - l1*m2*q1dot*sin(phai + q1) - l11*m1*q1dot*sin(phai + q1) - l1*m2*phai_dot*sin(phai + q1) + 2*d*mw*phai_dot*sin(phai) - l22*m2*phai_dot*sin(phai + q1 + q2) - l22*m2*q1dot*sin(phai + q1 + q2) - l22*m2*q2dot*sin(phai + q1 + q2);
C24 = - l1*m2*phai_dot*sin(phai + q1) - l11*m1*phai_dot*sin(phai + q1) - l1*m2*q1dot*sin(phai + q1) - l11*m1*q1dot*sin(phai + q1) - l22*m2*phai_dot*sin(phai + q1 + q2) - l22*m2*q1dot*sin(phai + q1 + q2) - l22*m2*q2dot*sin(phai + q1 + q2);
C25 = -l22*m2*sin(phai + q1 + q2)*(phai_dot + q1dot + q2dot);
C26 = 0;
C27 = 0;
C31 = 0;
C32 = 0;
C33 = -l1*l22*m2*q2dot*sin(q2);
C34 = -l1*l22*m2*q2dot*sin(q2);
C35 = -l1*l22*m2*sin(q2)*(phai_dot + q1dot + q2dot);
C36 = 0;
C37 = 0;
C41 = 0;
C42 = 0;
C43 = -l1*l22*m2*q2dot*sin(q2);
C44 = -l1*l22*m2*q2dot*sin(q2);
C45 = -l1*l22*m2*sin(q2)*(phai_dot + q1dot + q2dot);
C46 = 0;
C47 = 0;
C51 = 0;
C52 = 0;
C53 = l1*l22*m2*sin(q2)*(phai_dot + q1dot);
C54 = l1*l22*m2*sin(q2)*(phai_dot + q1dot);
C55 = 0;
C56 = 0;
C57 = 0;
C61 = 0;
C62 = 0;
C63 = 0;
C64 = 0;
C65 = 0;
C66 = 0;
C67 = 0;
C71 = 0;
C72 = 0;
C73 = 0;
C74 = 0;
C75 = 0;
C76 = 0;
C77 = 0;

C = [C11 C12 C13 C14 C15 C16 C17;
     C21 C22 C23 C24 C25 C26 C27;
     C31 C32 C33 C34 C35 C36 C37;
     C41 C42 C43 C44 C45 C46 C47;
     C51 C52 C53 C54 C55 C56 C57;
     C61 C62 C63 C64 C65 C66 C67;
     C71 C72 C73 C74 C75 C76 C77;
     ];

 
% E=[cos(phai)/r cos(phai)/r 0 0;
%     sin(phai)/r sin(phai)/r 0 0;
%     b/r -b/r 0 0;
%     0 0 1 0;
%     0 0 0 1;
%     1 0 0 0;
%     0 1 0 0
%     ];

J11 = (r*cos(phai))/2 - (r*(l1*sin(phai + q1) + l2*sin(phai + q1 + q2)))/(2*b) - (d*r*sin(phai))/(2*b);
J12 = (r*cos(phai))/2 + (r*(l1*sin(phai + q1) + l2*sin(phai + q1 + q2)))/(2*b) + (d*r*sin(phai))/(2*b);
J13 = - l1*sin(phai + q1) - l2*sin(phai + q1 + q2);
J14 = -l2*sin(phai + q1 + q2);
J21 = (r*sin(phai))/2 + (r*(l1*cos(phai + q1) + l2*cos(phai + q1 + q2)))/(2*b) + (d*r*cos(phai))/(2*b);
J22 = (r*sin(phai))/2 - (r*(l1*cos(phai + q1) + l2*cos(phai + q1 + q2)))/(2*b) - (d*r*cos(phai))/(2*b);
J23 = l1*cos(phai + q1) + l2*cos(phai + q1 + q2);
J24 = l2*cos(phai + q1 + q2);
J31 = (r*cos(phai))/2 - (d*r*sin(phai))/(2*b);
J32 = (r*cos(phai))/2 + (d*r*sin(phai))/(2*b);
J33 = 0;
J34 = 0;
J41 = (r*sin(phai))/2 + (d*r*cos(phai))/(2*b);
J42 = (r*sin(phai))/2 - (d*r*cos(phai))/(2*b);
J43 = 0;
J44 = 0;

J = [J11 J12 J13 J14;
     J21 J22 J23 J24;
     J31 J32 0 0;
     J41 J42 0 0];
 
dJ11 = (r^2*(qldot - qrdot)*(l1*cos(phai + q1) + d*cos(phai) + b*sin(phai) + l2*cos(phai + q1 + q2)))/(4*b^2) - (q1dot*r*(l1*cos(phai + q1) + l2*cos(phai + q1 + q2)))/(2*b) - (l2*q2dot*r*cos(phai + q1 + q2))/(2*b);
dJ12 = (q1dot*r*(l1*cos(phai + q1) + l2*cos(phai + q1 + q2)))/(2*b) - (r^2*(qldot - qrdot)*(l1*cos(phai + q1) + d*cos(phai) - b*sin(phai) + l2*cos(phai + q1 + q2)))/(4*b^2) + (l2*q2dot*r*cos(phai + q1 + q2))/(2*b);
dJ13 = (r*(qldot - qrdot)*(l1*cos(phai + q1) + l2*cos(phai + q1 + q2)))/(2*b) - l2*q2dot*cos(phai + q1 + q2) - q1dot*(l1*cos(phai + q1) + l2*cos(phai + q1 + q2));
dJ14 = -(l2*cos(phai + q1 + q2)*(2*b*q1dot + 2*b*q2dot - qldot*r + qrdot*r))/(2*b);
dJ21 = (r^2*(qldot - qrdot)*(l1*sin(phai + q1) - b*cos(phai) + d*sin(phai) + l2*sin(phai + q1 + q2)))/(4*b^2) - (q1dot*r*(l1*sin(phai + q1) + l2*sin(phai + q1 + q2)))/(2*b) - (l2*q2dot*r*sin(phai + q1 + q2))/(2*b);
dJ22 = (q1dot*r*(l1*sin(phai + q1) + l2*sin(phai + q1 + q2)))/(2*b) - (r^2*(qldot - qrdot)*(l1*sin(phai + q1) + b*cos(phai) + d*sin(phai) + l2*sin(phai + q1 + q2)))/(4*b^2) + (l2*q2dot*r*sin(phai + q1 + q2))/(2*b);
dJ23 = (r*(qldot - qrdot)*(l1*sin(phai + q1) + l2*sin(phai + q1 + q2)))/(2*b) - l2*q2dot*sin(phai + q1 + q2) - q1dot*(l1*sin(phai + q1) + l2*sin(phai + q1 + q2));
dJ24 = -(l2*sin(phai + q1 + q2)*(2*b*q1dot + 2*b*q2dot - qldot*r + qrdot*r))/(2*b);
dJ31 = (r^2*(d*cos(phai) + b*sin(phai))*(qldot - qrdot))/(4*b^2);
dJ32 = -(r^2*(d*cos(phai) - b*sin(phai))*(qldot - qrdot))/(4*b^2);
dJ33 = 0;
dJ34 = 0;
dJ41 = -(r^2*(b*cos(phai) - d*sin(phai))*(qldot - qrdot))/(4*b^2);
dJ42 = -(r^2*(b*cos(phai) + d*sin(phai))*(qldot - qrdot))/(4*b^2);
dJ43 = 0;
dJ44 = 0;
dJ = [dJ11 dJ12 dJ13 dJ14;
     dJ21 dJ22 dJ23 dJ24;
     dJ31 dJ32 0 0;
     dJ41 dJ42 0 0;
     ];


S_q = [(r/2)*cos(phai)-d*r/(2*b)*sin(phai) (r/2)*cos(phai)+d*r/(2*b)*sin(phai) 0 0;
     (r/2)*sin(phai)+d*r/(2*b)*cos(phai) (r/2)*sin(phai)-d*r/(2*b)*cos(phai) 0 0;
     r/(2*b) -r/(2*b) 0 0;
     0 0 1 0;
     0 0 0 1;
     1 0 0 0;
     0 1 0 0;
     ];
dS_q = zeros(7,4);
dS_q(1,1) = (r^2*(d*cos(phai) + b*sin(phai))*(qldot - qrdot))/(4*b^2);
dS_q(1,2) =  -(r^2*(d*cos(phai) - b*sin(phai))*(qldot - qrdot))/(4*b^2);
dS_q(2,1) = -(r^2*(b*cos(phai) - d*sin(phai))*(qldot - qrdot))/(4*b^2);
dS_q(2,2) =   -(r^2*(b*cos(phai) + d*sin(phai))*(qldot - qrdot))/(4*b^2);
% dS_q =[(r^2*(d*cos(phai) + b*sin(phai))*(qldot - qrdot))/(4*b^2)  -(r^2*(d*cos(phai) - b*sin(phai))*(qldot - qrdot))/(4*b^2) 0 0;
%       (-(r^2*(b*cos(phai) - d*sin(phai))*(qldot - qrdot))/(4*b^2) -(r^2*(b*cos(phai) + d*sin(phai))*(qldot - qrdot))/(4*b^2) 0 0;
%      0 0 0 0;
%      0 0 0 0;
%      0 0 0 0;
%      0 0 0 0;
%      0 0 0 0
%      ];
end

J_inv = inv(J);
% cond_J = cond(J)
% det_J = det(J)

global E_bar C_bar M_bar 
M_bar = S_q' * M * S_q * J_inv;

% C_bar = S_q' * (M * (dS_q *  pinv(J) + S_q * (-pinv(J) * dJ * pinv(J))) + C * S_q * J_inv); 
C_bar = S_q' * (M * (dS_q *  J_inv + S_q * (-J_inv * dJ * J_inv)) + C * S_q * J_inv); 

% E_bar = S_q' * E
E_bar = eye(4);

% if t>10
%     dtol = [5; 5; 2; 2];
% else
%     dtol = [0; 0; 0; 0];
% end

dtol = [2*sin(3*t);2*cos(3*t);5*sin(1.5*t)-5;5*cos(1.5*t)+5];

% dtol = [0.2*sin(3*t);0.2*cos(3*t);0.5*sin(1.5*t);0.5*cos(1.5*t)];


% dtol = [0;0;5*sin(1.5*t);5*cos(1.5*t)];

% dtol = [0; 0; 0; 0];

ddX=pinv(M_bar)*(E_bar*utol - C_bar*dX + dtol); %  - G_bar 

v = J_inv * dX;

global velocity_to_next_loop
% velocity_to_next_loop = [v(1);v(2); (r/(2*b))*(v(1)-v(2));v(3);v(4)]
velocity_to_next_loop = [v(1);v(2); r*(v(1)-v(2))/(2*b);v(3);v(4)];
% velocity_to_next_loop = [v(1);v(2); r*(qrdot-qldot)/(2*b);v(3);v(4)];
d_theta = velocity_to_next_loop;

% dXdt = [dX; ddX; v;];
sys(1:4) = dX;
sys(5:8) = ddX;
sys(9:13) = d_theta;  % [the_R;the_L;phai;q1;q2]


function sys=mdlOutputs(t,x,u)
sys(1:8)=x(1:8);
sys(9:13) = x(9:13);
global q dq
sys(14:18) = q(1:5);
sys(19:23) = dq(1:5);