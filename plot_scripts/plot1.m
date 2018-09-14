result_dir = '/results/mega_script/';

depth = 3;


colors = {[255 192 203]/255, [173 216 230]/255, [170, 227, 60]/255, [255,165,0]/255,[147, 112, 216]/255 };

[~,~,test_data,test_norm,~,~,~,...
    ~,~,~,~,~,~,~]=main_trainGP('mnist',depth);

close all
files = ls(['results/mega_script/p_ub_ReLU_mnist_',int2str(depth),'*']);
files = strsplit(files);
files = files(1:end-1);
%for ii = 1:length(files)
for ii = 4:4
    load(files{ii});
    epsilons = T.epsilons;
    feats = T.pixel2Modify;
    sigma_b_2 = T.theta;
    sigma_w_2 = T.sigma_prior;
    %respampling with bigger delta
    delta = 0:0.01:1.0;
    p_ub = cell(size(T.p_ub));
    test_point_init = test_data(T.testPointIdx,:);
    normTestPoint =  test_norm(T.testPointIdx);
    hs = cell(size(epsilons));
    for hh = 1:length(epsilons)
        hs{hh} = figure();
        title(['epsilons = ', num2str(epsilons(hh)) ]);
    end
    for jj = 1:length(feats)

        p_ub{jj} = zeros(length(delta),length(epsilons));
        
        for hh = 1:length(epsilons)
            testPoint = prioritize_pixel2mod(test_point_init,feats{jj} );
            pix2Mod =  1:length(feats{jj});
            [x_L, x_U] = compute_hyper_rectangle(epsilons(hh),testPoint*normTestPoint, pix2Mod,'mnist');
            %testPoint = get_spherical_coordinates(testPoint,[]);
            [x_L,x_U] = compute_angle_space_hyper_rectangle(x_L,x_U);
            [K,alpha_L,alpha_U] =  compute_K_constant_for_ReLU(x_L,x_U,T.theta,T.sigma_prior,length(x_L)+1,depth);
            ub_mu = T.sup_mu((jj-1)*length(epsilons) + hh);
            ub_xi = T.sup_xi((jj-1)*length(epsilons) + hh);
            sup_d = 2.0*ub_xi;
            fun = @(z)( sqrt( log( 2*K*alpha_U./z + 1 )  ) );
            eta_int = 12*integral(fun,0,0.5*sup_d);
            for ll = 1:length(delta)
                eta = delta(ll) - (ub_mu + eta_int);
                if eta <=0
                    p_ub{jj}(ll,hh) = 1;
                else
                    p_ub{jj}(ll,hh) = exp( - eta^2 / (2*ub_xi) );
                end
            end
            
            %1d plotting of the results
            figure(hs{hh})
            hold on
            plot(delta,p_ub{jj}(:,hh),'LineWidth',3.5,'Color',colors{jj});
            %
            
        end
    end
    for hh = 1:length(epsilons)
        figure(hs{hh});
        hold on
        grid on
        ylim([0,1.0])
        xlim([0,0.3])
        xlabel('\delta');
        ylabel('\phi_1');
        %legend({'Feature 1','Feature 2','Feature 3','Feature 4','Feature 5'})
        set(gca,'FontSize',20)
        set(gca,'LineWidth',2.0)
        set(gca,'YTick',0:0.2:1.0)
        width = 519;
        height = 270;
        set(gcf,'units','points','position',[10,10,width,height])
    end
    %
end
figure
hold on
for ii = 1:5
    plot(0,0,'LineWidth',3.5,'Color',colors{ii})
end
legend({'Feature 1','Feature 2','Feature 3','Feature 4','Feature 5'},'Orientation','horizontal','Location', 'northoutside')
set(gca,'FontSize',20)
set(gca,'LineWidth',2.0)
set(gca, 'visible', 'off')