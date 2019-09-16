%run after main.m // LTV1.mod

%mortgage loans b_m
figure('Name','mortgage loans');
subplot(3,2,1);
plot(b_m_obs_epsiA,'color','black','lineWidth',3);
legend('productivity shock');
subplot(3,2,2);
plot(b_m_obs_epsiJ,'color','black','lineWidth',3);
legend('housing preference shock');
subplot(3,2,3);
plot(b_m_obs_epsiHd,'color','black','lineWidth',3);
legend('housing depr. shock');
subplot(3,2,4);
plot(b_m_obs_epsiSm,'color','black','lineWidth',3);
legend('housing risk shock');
subplot(3,2,5);
plot(b_m_obs_epsiWb,'color','black','lineWidth',3);
legend('banker net worth');
subplot(3,2,6);
plot(b_m_obs_epsiSH,'color','black','lineWidth',3);
legend('mortgage bank risk shock');

%output growth dy_data
figure('Name','output');
subplot(3,2,1);
plot(Y_net_obs_epsiA,'color','black','lineWidth',3);
legend('productivity shock');
subplot(3,2,2);
plot(Y_net_obs_epsiJ,'color','black','lineWidth',3);
legend('housing preference shock');
subplot(3,2,3);
plot(Y_net_obs_epsiHd,'color','black','lineWidth',3);
legend('housing depr. shock');
subplot(3,2,4);
plot(Y_net_obs_epsiSm,'color','black','lineWidth',3);
legend('housing risk shock');
subplot(3,2,5);
plot(Y_net_obs_epsiWb,'color','black','lineWidth',3);
legend('banker net worth');
subplot(3,2,6);
plot(Y_net_obs_epsiSH,'color','black','lineWidth',3);
legend('mortgage bank risk shock');


%household credit spread bsp_H_data
figure('Name','HH credit spread');
subplot(3,2,1);
plot(bsp_H_data_epsiA,'color','black','lineWidth',3);
legend('productivity shock');
subplot(3,2,2);
plot(bsp_H_data_epsiJ,'color','black','lineWidth',3);
legend('housing preference shock');
subplot(3,2,3);
plot(bsp_H_data_epsiHd,'color','black','lineWidth',3);
legend('housing depr. shock');
subplot(3,2,4);
plot(bsp_H_data_epsiSm,'color','black','lineWidth',3);
legend('housing risk shock');
subplot(3,2,5);
plot(bsp_H_data_epsiWb,'color','black','lineWidth',3);
legend('banker net worth');
subplot(3,2,6);
plot(bsp_H_data_epsiSH,'color','black','lineWidth',3);
legend('mortgage bank risk shock');