

maxDelta = 0.20;
p_ub{1} = p_ub{1}{1}(upper_bound_comp_opts.delta < maxDelta);
stat_probs_x_o = stat_probs_x_o(upper_bound_comp_opts.delta < maxDelta);

p_ub{2} = p_ub{2}{1}(upper_bound_comp_opts.delta < maxDelta);
stat_probs_x_c = stat_probs_x_c(upper_bound_comp_opts.delta < maxDelta);
delta = upper_bound_comp_opts.delta(upper_bound_comp_opts.delta < maxDelta);


figure
hold on
plot(delta,p_ub{1},'LineWidth',3.5,'Color', [173 216 230]/255);
plot(delta,stat_probs_x_o,'LineWidth',3.5,'Color', [173 216 230]/255,'LineStyle','--');

plot(delta,p_ub{2},'LineWidth',3.5,'Color', [255,165,0]/255);
plot(delta,stat_probs_x_c,'LineWidth',3.5,'Color', [255,165,0]/255,'LineStyle','--');

xlim([0,maxDelta])
ylim([0,1]) 
grid on
xlabel('\delta');


ylabel('\phi')
if upper_bound_comp_opts.invariance
    subscript = '2';
else
    subscript = '1';
end
legend({['$\hat{\phi}_',subscript,'(x^{o},T^{o}_\gamma,\delta)$'],...
    ['$\bar{\phi}_',subscript,'(x^{o},T^{o}_\gamma,\delta)$'],...
    ['$\hat{\phi}_',subscript,'(x^{*},T^{*}_\gamma,\delta)$'],...
    ['$\bar{\phi}_',subscript,'(x^{*},T^{*}_\gamma,\delta)$']},'Interpreter','latex')

set(gca,'FontSize',20)
set(gca,'YTick',0:0.2:1.0)
set(gca,'LineWidth',1.5)

