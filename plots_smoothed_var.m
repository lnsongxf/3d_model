%============================================
%SMOOTHED SHOCKS
figure;
plot(oo_.SmoothedVariables.A,'color','black','lineWidth',5);
title('A');

figure;
plot(oo_.SmoothedVariables.EJ,'color','black','lineWidth',5);
title('EJ');

figure;
plot(oo_.SmoothedVariables.EK,'color','black','lineWidth',5);
title('EK');

figure;
plot(oo_.SmoothedVariables.EH,'color','black','lineWidth',5);
title('EH');

figure;
plot(oo_.SmoothedVariables.ESe,'color','black','lineWidth',5);
title('ESe');

figure;
plot(oo_.SmoothedVariables.ESm,'color','black','lineWidth',5);
title('ESm');

figure;
plot(oo_.SmoothedVariables.ESF,'color','black','lineWidth',5);
title('ESF');

figure;
plot(oo_.SmoothedVariables.ESH,'color','black','lineWidth',5);
title('ESH');

figure;
plot(oo_.SmoothedVariables.EWe,'color','black','lineWidth',5);
title('EWe');

figure;
plot(oo_.SmoothedVariables.EWb,'color','black','lineWidth',5);
title('EWb');

figure;
plot(oo_.SmoothedVariables.EdH,'color','black','lineWidth',5);
title('EdH');

figure;
plot(oo_.SmoothedVariables.EdK,'color','black','lineWidth',5);
title('EdK');

%==========================================================
%KEY IRFS
figure;
subplot(4,3,1);
plot(dy_data_epsiA,'lineWidth',3,'color','black');
subplot(4,3,2);
plot(dy_data_epsiH,'lineWidth',3,'color','black');
subplot(4,3,3);
plot(dy_data_epsiHd,'lineWidth',3,'color','black');
subplot(4,3,4);
plot(dy_data_epsiHk,'lineWidth',3,'color','black');
subplot(4,3,5);
plot(dy_data_epsiJ,'lineWidth',3,'color','black');
subplot(4,3,6);
plot(dy_data_epsiK,'lineWidth',3,'color','black');
subplot(4,3,7);
plot(dy_data_epsiSe,'lineWidth',3,'color','black');
subplot(4,3,8);
plot(dy_data_epsiSF,'lineWidth',3,'color','black');
subplot(4,3,9);
plot(dy_data_epsiSH,'lineWidth',3,'color','black');
subplot(4,3,10);
plot(dy_data_epsiSm,'lineWidth',3,'color','black');
subplot(4,3,11);
plot(dy_data_epsiWb,'lineWidth',3,'color','black');
subplot(4,3,12);
plot(dy_data_epsiWe,'lineWidth',3,'color','black');
