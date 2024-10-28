close all;
% 
figure(1);
subplot(411);
plot(t,y(:,1),'b','linewidth',1.2);
hold on;
plot(t, y1(:,11), 'r--', 'linewidth', 1.2); % 红色虚线表示期望的X分量
% ylim([0 5]);
ylabel('X_1 (m)');
subplot(412);
plot(t,y(:,2),'b','linewidth',1.2);
hold on;
plot(t, y1(:,12), 'r--', 'linewidth', 1.2); 
% ylim([0.45 0.85]);
ylabel('X_2 (m)');
subplot(413);
plot(t,y(:,3),'b','linewidth',1.2);
hold on;
plot(t, y1(:,13), 'r--', 'linewidth', 1.2); 
% ylim([0 5]);
ylabel('X_3 (m)');
subplot(414);
plot(t,y(:,4),'b','linewidth',1.2);
hold on;
plot(t, y1(:,14), 'r--', 'linewidth', 1.2); 
% ylim([-0.1 0.1]);
xlabel('t (s)');ylabel('X_4 (m)');

figure(2);
subplot(411);
plot(t,y(:,5),'b','linewidth',1.2);
hold on;
plot(t, y1(:,15), 'r--', 'linewidth', 1.2); % 红色虚线表示期望的X分量
ylabel('dX_1 (m/s)');
subplot(412);
plot(t,y(:,6),'b','linewidth',1.2);
hold on;
plot(t, y1(:,16), 'r--', 'linewidth', 1.2); 
ylabel('dX_2 (m/s)');
subplot(413);
plot(t,y(:,7),'b','linewidth',1.2);
hold on;
plot(t, y1(:,17), 'r--', 'linewidth', 1.2); 
ylabel('dX_3 (m/s)');
subplot(414);
plot(t,y(:,8),'b','linewidth',1.2);
hold on;
plot(t, y1(:,18), 'r--', 'linewidth', 1.2);
xlabel('t (s)');ylabel('dX_4 (m/s)');
% 
figure(3);
subplot(411);
plot(t,utol(:,1),'b','linewidth',1.2);
ylabel('u_1 (Nm)');
subplot(412);
plot(t,utol(:,2),'b','linewidth',1.2);
ylabel('u_2 (Nm)');
subplot(413);
plot(t,utol(:,3),'b','linewidth',1.2);
ylabel('u_3 (Nm)');
subplot(414);
plot(t,utol(:,4),'b','linewidth',1.2);
xlabel('t (s)');ylabel('u_4 (Nm)');
% 
figure(4);
subplot(411);
plot(t,y1(:,1),'b','linewidth',1.2);
hold on;
plot(t, y1(:,19), 'r--', 'linewidth', 1.2); % BLF 上下界
plot(t, -y1(:,19), 'r--', 'linewidth', 1.2); 
t_common = linspace(0, 10, 1000);
e_interp = interp1(t, y1(:,1), t_common, 'linear');
RMSE = sqrt(mean(e_interp.^2));
IEA = trapz(t_common, abs(e_interp));
ITEA = trapz(t_common, t_common .*abs(e_interp));
e_interp = y1(:,1);
RMSE_2 = sqrt(mean(e_interp.^2));
IEA_2 = trapz(t, abs(e_interp));
ITEA_2 = trapz(t, t .*abs(e_interp));
text(t(end)-5.5, 0.33, sprintf('RMSE/IEA/ITEA: %.4f %.4f %.4f', RMSE_2, IEA_2, ITEA_2), 'VerticalAlignment', 'top'); % 将 MSE 值添加到图像上
text(t(end)-5.5, 0.43, sprintf('RMSE/IEA/ITEA: %.4f %.4f %.4f', RMSE, IEA, ITEA), 'VerticalAlignment', 'top'); % 将 MSE 值添加到图像上
ylabel('e_1 (m)');
% ITEA_3= 0; dt=0.01;
% for i = 1:length(e_interp)
%     ITEA_3 = ITEA_3 + t(i)*abs(e_interp(i))*dt;
% end

subplot(412);
plot(t,y1(:,2),'b','linewidth',1.2);
hold on;
plot(t, y1(:,19), 'r--', 'linewidth', 1.2);
plot(t, -y1(:,19), 'r--', 'linewidth', 1.2);
e_interp = interp1(t, y1(:,2), t_common, 'linear');
RMSE = sqrt(mean(e_interp.^2));
IEA = trapz(t_common, abs(e_interp));
ITEA = trapz(t_common, t_common .*abs(e_interp));
e_interp = y1(:,2);
RMSE_2 = sqrt(mean(e_interp.^2));
IEA_2 = trapz(t, abs(e_interp));
ITEA_2 = trapz(t, t .*abs(e_interp));
text(t(end)-5.5, 0.33, sprintf('RMSE/IEA/ITEA: %.4f %.4f %.4f', RMSE_2, IEA_2, ITEA_2), 'VerticalAlignment', 'top'); % 将 MSE 值添加到图像上
text(t(end)-5.5, 0.43, sprintf('RMSE/IEA/ITEA: %.4f %.4f %.4f', RMSE, IEA, ITEA), 'VerticalAlignment', 'top'); % 将 MSE 值添加到图像上
ylabel('e_2 (m)');
subplot(413);
plot(t,y1(:,3),'b','linewidth',1.2);
hold on;
plot(t, y1(:,20), 'r--', 'linewidth', 1.2);
plot(t, -y1(:,20), 'r--', 'linewidth', 1.2);
e_interp = interp1(t, y1(:,3), t_common, 'linear');
RMSE = sqrt(mean(e_interp.^2));
IEA = trapz(t_common, abs(e_interp));
ITEA = trapz(t_common, t_common .*abs(e_interp));
e_interp = y1(:,3);
RMSE_2 = sqrt(mean(e_interp.^2));
IEA_2 = trapz(t, abs(e_interp));
ITEA_2 = trapz(t, t .*abs(e_interp));
text(t(end)-5.5, 0.33, sprintf('RMSE/IEA/ITEA: %.4f %.4f %.4f', RMSE_2, IEA_2, ITEA_2), 'VerticalAlignment', 'top'); % 将 MSE 值添加到图像上
text(t(end)-5.5, 0.43, sprintf('RMSE/IEA/ITEA: %.4f %.4f %.4f', RMSE, IEA, ITEA), 'VerticalAlignment', 'top'); % 将 MSE 值添加到图像上
ylabel('e_3 (m)');
subplot(414);
plot(t,y1(:,4),'b','linewidth',1.2);
hold on;
plot(t, y1(:,20), 'r--', 'linewidth', 1.2);
plot(t, -y1(:,20), 'r--', 'linewidth', 1.2);
e_interp = interp1(t, y1(:,4), t_common, 'linear');
RMSE = sqrt(mean(e_interp.^2));
IEA = trapz(t_common, abs(e_interp));
ITEA = trapz(t_common, t_common .*abs(e_interp));
e_interp = y1(:,4);
RMSE_2 = sqrt(mean(e_interp.^2));
IEA_2 = trapz(t, abs(e_interp));
ITEA_2 = trapz(t, t .*abs(e_interp));
text(t(end)-5.5, 0.33, sprintf('RMSE/IEA/ITEA: %.4f %.4f %.4f', RMSE_2, IEA_2, ITEA_2), 'VerticalAlignment', 'top'); % 将 MSE 值添加到图像上
text(t(end)-5.5, 0.43, sprintf('RMSE/IEA/ITEA: %.4f %.4f %.4f', RMSE, IEA, ITEA), 'VerticalAlignment', 'top'); % 将 MSE 值添加到图像上
xlabel('t (s)');ylabel('e_4 (m)');
% 
figure(5);
subplot(411);
plot(t,y1(:,21),'b','linewidth',1.2);
ylabel('e2_1');
subplot(412);
plot(t,y1(:,22),'b','linewidth',1.2);
ylabel('e2_2');
subplot(413);
plot(t,y1(:,23),'b','linewidth',1.2);
ylabel('e2_3');
subplot(414);
plot(t,y1(:,24),'b','linewidth',1.2);
xlabel('t (s)');ylabel('e2_4');
% 
figure(6);
subplot(411);
plot(t,y1(:,25),'b','linewidth',1.2);
ylabel('fn_1');
subplot(412);
plot(t,y1(:,26),'b','linewidth',1.2);
ylabel('fn_2');
subplot(413);
plot(t,y1(:,27),'b','linewidth',1.2);
ylabel('fn_3');
subplot(414);
plot(t,y1(:,28),'b','linewidth',1.2);
ylabel('fn_4');
xlabel('t (s)');

figure(7);
subplot(411);
plot(t,y1(:,31),'b','linewidth',1.2);
ylabel('||W_{a1}||');
subplot(412);
plot(t,y1(:,32),'b','linewidth',1.2);
ylabel('||W_{a2}||');
subplot(413);
plot(t,y1(:,33),'b','linewidth',1.2);
ylabel('||W_{a3}||');
subplot(414);
plot(t,y1(:,34),'b','linewidth',1.2);
ylabel('||W_{a4}||');
xlabel('t (s)');

figure(8);
subplot(511);
plot(t,y1(:,29),'b','linewidth',1.2);
ylabel('hatQ');

subplot(512);
% 计算真实Q
phi_real = y1(:,30);
real_Q = zeros(1,length(phi_real)); psi=0.2; dt=0.01;

for i = 1:length(phi_real)
    for j = 0:length(phi_real) - i
        dis_factor = exp(-dt*j/psi);
        real_Q(i) = real_Q(i) + dis_factor * phi_real(i+j)* dt;
    end
end
plot(t, real_Q, 'r--', 'linewidth', 1.2);
ylabel('Real cost function');

subplot(513);
plot(t, phi_real, 'r--', 'linewidth', 1.2);
ylabel('reward');

Q_estimation = y1(:,29);
subplot(514); % 计算TD误差
TD_error(1)=0;
for i=2:1:length(phi_real)
    TD_error(i) = phi_real(i) - 1/psi*Q_estimation(i) + (Q_estimation(i)-Q_estimation(i-1))/dt;
end
plot(t,TD_error,'b','linewidth',1.2);
ylabel('TD_error');

subplot(515);
plot(t,y1(:,30),'b','linewidth',1.2);
ylabel('||W_c||');
xlabel('t (s)');


