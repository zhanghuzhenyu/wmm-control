function [sys,x0,str,ts]=s_function(t,x,u,flag)
switch flag,
case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;
case 1,
    sys=mdlDerivatives(t,x,u);
case 3,
    sys=mdlOutputs(t,x,u);
case {2,4,9 }
    sys = [];
otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end
function [sys,x0,str,ts]=mdlInitializeSizes

global lra lra1 lra2 lra3 lra4 
global c_a1 c_a2 c_a3 c_a4 width_a1 width_a2 width_a3 width_a4
global Wa1_limit Wa2_limit Wa3_limit Wa4_limit
global Node_a1 Node_a2 Node_a3 Node_a4
%%%%%%%%%%%% Actor %%%%%%%%%%%%
% a = [-0.5 -0.4 -0.3 -0.2 -0.15 -0.1 -0.05 0.0 0.05 0.1 0.15 0.2 0.3 0.4 0.5].*10;
a = -6:0.5:6;
% aij1 = -5:0.1:5; aij2 = -5:0.1:5; aij3 = -5:0.1:5; aij4 = -5:0.1:5;
aij1 =a; aij2 = a; aij3 = a; aij4 = a;
% c_a1 =[aij1;aij1]; c_a2 =[aij2;aij2]; c_a3 =[aij3;aij3]; c_a4 =[aij4;aij4];
c_a1 =aij1; c_a2 =aij2; c_a3 =aij3; c_a4 =aij4;
width_a1 = 1; width_a2 = 1; width_a3 = 1; width_a4 = 1;
Node_a1 = length(c_a1);
Node_a2 = length(c_a2);
Node_a3 = length(c_a3);
Node_a4 = length(c_a4);
Wa1_limit=5*ones(Node_a1,1); 
Wa2_limit=5*ones(Node_a2,1); 
Wa3_limit=5*ones(Node_a3,1); 
Wa4_limit=5*ones(Node_a4,1); 
lra=1; % 初始值 4
lra1=1; lra2=1; lra3=1.5; lra4=2.5;

global cij Node_c Sc Wc Wc1 lrc width_c phi D R
cij = -0.5:0.01:0.5; % -0.15:0.03:0.15;
Node_c = length(cij);
Sc=zeros(Node_c,1);
Wc=zeros(Node_c,1);
Wc1=5*ones(Node_c,1);
lrc=0.1;  % 初始值 0.2
width_c=0.1;
phi=0.2;  % 初始值 0.2
D=500*eye(4);   % 500  0.1 是Actor-only的参数
R=0.1*eye(4);

global Ki
Ki=1;  % 初始值 0.2

sizes = simsizes;
sizes.NumContStates  = Node_a1 + Node_a2 + Node_a3 +Node_a4 +Node_c; % 4*7+7
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 38;% 其中前四个元素是控制律utol，输入进sys中
sizes.NumInputs      = 13; % 8 + 5agl
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;
sys=simsizes(sizes);
% x0=[];
str=[];
ts=[0 0];
x0  = 0*ones(1,Node_a1 + Node_a2 + Node_a3 +Node_a4 +Node_c);


function sys=mdlDerivatives(t,x,u)

if (t>=-1)

global X_init 
% Xd = [0.2*t+X_init(1) - 0.1; -0.10*sin(0.2*pi*t)+X_init(2) + 0.3; 0.2*t+ X_init(3) + 0.3; X_init(4) + 0.1];
Xd = [0.2*t+0.65; -0.10*sin(0.2*pi*t)+0.75; 0.2*t+ 0.4; 0.1];
dXd = [0.2; -0.040*pi*cos(0.2*pi*t); 0.2; 0];
ddXd = [0; +0.008*pi^2*sin(0.2*pi*t); 0; 0];

X = u(1:4);
dX = u(5:8);

% 求e2 --> 无法插入 只能用e1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e1=X-Xd;
de1=dX-dXd;

k1 = [10;10;10;10]*5; 

up_bound = [0.5;0.5;0.5;0.5]; % 初始上边界
down_bound = [0.1;0.1;0.1;0.1]; % 最低边界
rate = [2;2;2;2];

miu = zeros(4,1);
dmiu = zeros(4,1);

for i = 1:4
    kb(i) = (up_bound(i)-down_bound(i))*exp(-rate(i)*t)+down_bound(i);
    kb_dot(i) = -rate(i)*(up_bound(i)-down_bound(i))*exp(-rate(i)*t);
    kb_ddot(i) = rate(i)^2*(up_bound(i)-down_bound(i))*exp(-rate(i)*t);
    miu(i) = -k1(i)*e1(i)+kb_dot(i)/kb(i)*e1(i)+dXd(i);
    dmiu(i) = ddXd(i) - k1(i)*de1(i) + (kb_ddot(i)*kb(i)-kb_dot(i)^2)/(kb(i)^2)*e1(i) + kb_dot(i)/kb(i)*de1(i);
end

e2=dX-miu;

%% Critic NN
global cij Node_c Sc Wc Wc1 lrc width_c phi D R
global utol hatI
Wc=[x(1:Node_c)];
for i=1:1:Node_c
    Sc(i,1)=exp(-(e1-cij(:,i))'*(e1-cij(:,i))/width_c^2);
end
hatI=Wc'*Sc;

for i=1:1:Node_c
    A(i,1)=-(Sc(i)/phi)+(-2*Sc(i)*(e1-cij(1,i))'/(width_c)^2)*de1;
end
%     dWc=-lrc*(e1'*D*e1+tau0'*R*tau0+Wc'*A)*A;

if norm(Wc)<=norm(Wc1)||norm(Wc)==norm(Wc1)&&Wc'*((e1'*D*e1+utol'*R*utol+Wc'*A)*A)>0
    dWc=-lrc*(e1'*D*e1+utol'*R*utol+Wc'*A)*A;
elseif norm(Wc)==norm(Wc1)&&Wc'*((e1'*D*e1+utol'*R*utol+Wc'*A)*A)<=0
    dWc=-lrc*(e1'*D*e1+utol'*R*utol+Wc'*A)*A+lrc*(Wc'*((e1'*D*e1+utol'*R*utol+Wc'*A)*A)/(norm(Wc))^2)*Wc;
else
    disp("erorr!")
end
global inst_reward
inst_reward1 = e1'*D*e1
inst_reward2 = utol'*R*utol
inst_reward = e1'*D*e1+utol'*R*utol;
%% Actor NN
global Node_a Wa lra 
global Ki
global c_a1 c_a2 c_a3 c_a4 width_a1 width_a2 width_a3 width_a4
global Wa1_limit Wa2_limit Wa3_limit Wa4_limit
global Node_a1 Node_a2 Node_a3 Node_a4
global W1 W2 W3 W4
start_W1 = Node_c + 1; end_W1 = Node_c + Node_a1;
start_W2 = end_W1 + 1; end_W2 = Node_c + Node_a1 + Node_a2;
start_W3 = end_W2 + 1; end_W3 = Node_c + Node_a1 + Node_a2 + Node_a3;
start_W4 = end_W3 + 1; end_W4 = Node_c + Node_a1 + Node_a2 + Node_a3 + Node_a4;
W1 = x(start_W1:end_W1);
W2 = x(start_W2:end_W2);
W3 = x(start_W3:end_W3);
W4 = x(start_W4:end_W4);

% % % Wa = [W1 W2 W3 W4];

% za1=[e1(1);de1(1)];
% za2=[e1(2);de1(2)];
% za3=[e1(3);de1(3)];
% za4=[e1(4);de1(4)];

za1=e2(1);
za2=e2(2);
za3=e2(3);
za4=e2(4);

for i=1:Node_a1
    Sa1(i,1)=exp(-(za1-c_a1(:,i))'*(za1-c_a1(:,i))/(width_a1)^2);
end
for i=1:Node_a2
    Sa2(i,1)=exp(-(za2-c_a2(:,i))'*(za2-c_a2(:,i))/(width_a2)^2);
end
for i=1:Node_a3
    Sa3(i,1)=exp(-(za3-c_a3(:,i))'*(za3-c_a3(:,i))/(width_a3)^2);
end
for i=1:Node_a4
    Sa4(i,1)=exp(-(za4-c_a4(:,i))'*(za4-c_a4(:,i))/(width_a4)^2);
end  
    
global Fnn
Fnn1 = W1'*Sa1;
Fnn2 = W2'*Sa2;
Fnn3 = W3'*Sa3;
Fnn4 = W4'*Sa4;

miu_a = Fnn1 + Fnn2 + Fnn3 + Fnn4;
Fnn = [Fnn1;Fnn2;Fnn3;Fnn4];

% updata_num = miu_a+Ki*hatI
% tanh_update_num = tanh(miu_a+Ki*hatI)

Ki = 1;
temp1 = miu_a
temp2 = Ki*hatI

% prho1=Sa1.*(miu_a+Ki*hatI);
% prho2=Sa2.*(miu_a+Ki*hatI);
% prho3=Sa3.*(miu_a+Ki*hatI);
% prho4=Sa4.*(miu_a+Ki*hatI);

prho1=Sa1.*za1;
prho2=Sa2.*za2;
prho3=Sa3.*za3;
prho4=Sa4.*za4;
% 
% prho1=Sa1.*(miu_a+Ki*hatI+ Fnn1+za1);
% prho2=Sa2.*(miu_a+Ki*hatI+ Fnn2+za2);
% prho3=Sa3.*(miu_a+Ki*hatI+ Fnn3+za3);
% prho4=Sa4.*(miu_a+Ki*hatI+ Fnn4+za4);

% a = 0.2;  % 简单堆叠Actor的误差函数
% prho1=Sa1.*(a*(miu_a+hatI) + (1-a)*za1);
% prho2=Sa2.*(a*(miu_a+hatI) + (1-a)*za2);
% prho3=Sa3.*(a*(miu_a+hatI) + (1-a)*za3);
% prho4=Sa4.*(a*(miu_a+hatI) + (1-a)*za4);

global lra1 lra2 lra3 lra4 

if norm(W1)<=norm(Wa1_limit)||norm(W1)==norm(Wa1_limit)&&W1'*prho1>0
    dW1=-lra1*prho1;
elseif norm(W1)==norm(Wa1_limit)&&W1'*prho1<=0
    dW1=-lra1*prho1+lra1*((W1'*prho1)/(norm(W1))^2)*W1;
end

if norm(W2)<=norm(Wa2_limit)||norm(W2)==norm(Wa2_limit)&&W2'*prho2>0
    dW2=-lra2*prho2;
elseif norm(W2)==norm(Wa2_limit)&&W2'*prho2<=0
    dW2=-lra2*prho2+lra2*((W2'*prho2)/(norm(W2))^2)*W2;
end

if norm(W3)<=norm(Wa3_limit)||norm(W3)==norm(Wa3_limit)&&W3'*prho3>0
    dW3=-lra3*prho3;
elseif norm(W3)==norm(Wa3_limit)&&W3'*prho3<=0
    dW3=-lra3*prho3+lra3*((W3'*prho3)/(norm(W3))^2)*W3;
end

if norm(W4)<=norm(Wa4_limit)||norm(W4)==norm(Wa4_limit)&&W4'*prho4>0
    dW4=-lra4*prho4;
elseif norm(W4)==norm(Wa4_limit)&&W4'*prho4<=0
    dW4=-lra4*prho4+lra4*((W4'*prho4)/(norm(W4))^2)*W4;
end

sys(1:Node_c) = dWc;

sys(start_W1:end_W1) = dW1;
sys(start_W2:end_W2) = dW2;
sys(start_W3:end_W3) = dW3;
sys(start_W4:end_W4) = dW4;


else
    global Node_a Node_c
    global Node_a1 Node_a2 Node_a3 Node_a4
    sys(1:Node_a1 + Node_a2 + Node_a3 +Node_a4 +Node_c) = 0;

end

function sys=mdlOutputs(t,x,u)

global X_init 

% Xd = [0.2*t+X_init(1) - 0.1; -0.10*sin(0.2*pi*t)+X_init(2) + 0.3; 0.2*t+ X_init(3) + 0.3; X_init(4) + 0.1];  
Xd = [0.2*t+0.65; -0.10*sin(0.2*pi*t)+0.75; 0.2*t+ 0.4; 0.1];
dXd = [0.2; -0.040*pi*cos(0.2*pi*t); 0.2; 0];
ddXd = [0; +0.008*pi^2*sin(0.2*pi*t); 0; 0];

X = u(1:4);
dX = u(5:8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k1 = [10;10;10;10]*5; 
k2 = 5;

up_bound = [0.5;0.5;0.5;0.5]; % 初始上边界
down_bound = [0.1;0.1;0.1;0.1]; % 最低边界
rate = [2;2;2;2];

e1=X-Xd;
de1=dX-dXd;
miu = zeros(4,1);
dmiu = zeros(4,1);

for i = 1:4
    kb(i) = (up_bound(i)-down_bound(i))*exp(-rate(i)*t)+down_bound(i);
    kb_dot(i) = -rate(i)*(up_bound(i)-down_bound(i))*exp(-rate(i)*t);
    kb_ddot(i) = rate(i)^2*(up_bound(i)-down_bound(i))*exp(-rate(i)*t);
    miu(i) = -k1(i)*e1(i)+kb_dot(i)/kb(i)*e1(i)+dXd(i);
    dmiu(i) = ddXd(i) - k1(i)*de1(i) + (kb_ddot(i)*kb(i)-kb_dot(i)^2)/(kb(i)^2)*e1(i) + kb_dot(i)/kb(i)*de1(i);
end

e2=dX-miu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global E_bar  M_bar  C_bar

term1 = zeros(4,1);
for i = 1:4
    term1(i)=- e1(i)/(kb(i)*kb(i)-e1(i)*e1(i));
end
disp('')

global utol
global Fnn
if isempty(Fnn)
    Fnn = zeros(4,1);
end

% utol=pinv(E_bar) * (-e1 - k2*M_bar*e2 + C_bar*miu + M_bar*dmiu); % 什么也没有
% utol=pinv(E_bar) * (term1 - k2*M_bar*e2 + C_bar*miu + M_bar*dmiu); % BLF
% utol=pinv(E_bar) * (-e1 - k2*M_bar*e2 + C_bar*miu + M_bar*dmiu + Fnn); % RBF
% utol=pinv(E_bar) * (term1 - k2*M_bar*e2 + C_bar*miu + M_bar*dmiu + Fnn); % RBF+BLF

utol=pinv(E_bar) * (term1 - k2*M_bar*e2 + C_bar*miu + M_bar*dmiu + Fnn); % 测试Critic  + Fnn


% utol=pinv(E_bar) * (- k2*e2 - e1  + C_bar*miu + M_bar*dmiu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global hatI
global inst_reward
global W1 W2 W3 W4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sys(1:4)=utol(1:4);
sys(5:8)=e1(1:4);
sys(9:12)=de1(1:4);

sys(15:18) = Xd(1:4);
sys(19:22) = dXd(1:4);
sys(23) = kb(1);
sys(24) = kb(3);
sys(25:28) = e2(1:4);  % e2 或 z


sys(29:32) = Fnn(1:4);
if isempty(hatI)
    hatI = 0;
end
sys(33) = hatI;
if isempty(inst_reward)
    inst_reward = 0;
end
sys(34) = inst_reward;

sys(35) = norm(W1);
sys(36) = norm(W2);
sys(37) = norm(W3);
sys(38) = norm(W4);




