%close all
xo = [1.0000    1.0000    1.0000    1.0000    0.9107    0.6861    0.4276    0.2205    0.0940    0.0332    0.0097    0.0023    0.0005    0.0001    0.0000    0.0000...
    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000];
delta = 0:0.01:0.25;


xc = [ 1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    0.9646    0.8381    0.6549    0.4602    0.2909    0.1654    0.0845    0.0389...
    0.0161    0.0060    0.0020    0.0006    0.0002    0.0000    0.0000    0.0000    0.0000    0.0000];



xo_phi2 = [ 1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    0.8552    0.4409    0.1881    0.0664    0.0194    0.0047    0.0009    0.0002    0.0000    0.0000    0.0000    0.0000 ...
 0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000];

xc_phi2 = [1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    0.9205    0.5818    0.3307    0.1691    0.0777    0.0321    0.0120...
 0.0040    0.0012    0.0003    0.0001    0.0000    0.0000    0.0000    0.0000];

c=[59	76	192;
   180	4	38];


figure
hold on
plot(delta(1:21),xo(1:21),'LineWidth',3.5,'Color', [255 192 203]/255);
plot(delta(1:21),xc(1:21),'LineWidth',3.5,'Color', [173 216 230]/255);

plot(delta(1:21),xo_phi2(1:21),'LineWidth',3.5,'Color', [170, 227, 60]/255);
plot(delta(1:21),xc_phi2(1:21),'LineWidth',3.5,'Color', [255,165,0]/255);

xlim([0,0.2])
ylim([0,1]) 
grid on
xlabel('\delta');
%ylabel('\phi_1(x,T^x_\gamma,\delta)');
%legend({'x^o','x^c'})

ylabel('\phi')
legend({'\phi_1(x^{o},T^{o}_\gamma,\delta)',...
    '\phi_1(x^{*},T^{*}_\gamma,\delta)',...
    '\phi_2(x^{o},T^{o}_\gamma,\delta)',...
    '\phi_2(x^{*},T^{*}_\gamma,\delta)'},'Orientation','horizontal','Location', 'northoutside')

set(gca,'FontSize',18)
set(gca,'YTick',0:0.2:1.0)
set(gca,'LineWidth',1.5)