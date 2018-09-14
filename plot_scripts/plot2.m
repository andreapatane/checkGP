%close all
clear all
clc


result_dir = 'results/variance/';


testPointIdxs = [4];

depth_analysis =  true;
%depths = [1:2:11];
currMin = +inf;
currMax = -inf;
for ii = 1:length(testPointIdxs)
    
    curr_idx = testPointIdxs(ii);
    if depth_analysis
        load([result_dir,'variance_and_depth_ReLU_mnist_',int2str(curr_idx),'.mat']);
        if currMax < max(max(T.lls./T.refs))
            currMax =  max(max(T.lls./T.refs));
        end
        if currMin > min(min(T.lls./T.refs))
            currMin =  min(min(T.lls./T.refs));
        end
    end
end


for ii = 1:length(testPointIdxs)
    
    curr_idx = testPointIdxs(ii);
    if depth_analysis
        load([result_dir,'variance_and_depth_ReLU_mnist_',int2str(curr_idx),'.mat']); 
        depths = T.depth;
    else
        load([result_dir,'variance_ReLU_mnist_',int2str(curr_idx),'.mat']);
    end
    %figure
    %hold on
    if depth_analysis
        T.lls = reshape(T.lls(1,:,:),size(T.lls,2),size(T.lls,3));
        T.uus = reshape(T.uus(1,:,:),size(T.uus,2),size(T.uus,3));
        T.refs = reshape(T.refs(1,:,:),size(T.refs,2),size(T.refs,3));
    end
    %if depth_analysis
    %    surf(depths,T.actualTrainingSamples_vec,T.lls)
    %    ylim([T.actualTrainingSamples_vec(1),T.actualTrainingSamples_vec(end)])
    %    xlim([depths(1),depths(end)])
    %else
    %    surf(T.actualTrainingSamples_vec,T.sigma_w_2_vec,T.lls)
    %    xlim([T.actualTrainingSamples_vec(1),T.actualTrainingSamples_vec(end)])
    %    ylim([T.sigma_w_2_vec(1),T.sigma_w_2_vec(end)])
    %end
    %colormap('coolwarm')
    %colorbar
    %hold off
    
    figure
    hold on
    grid on
    if depth_analysis
        surf(depths,T.actualTrainingSamples_vec,T.lls./T.refs)
        ylim([T.actualTrainingSamples_vec(1),T.actualTrainingSamples_vec(end)])
        xlim([depths(1),depths(end)])
        zlim([currMin,currMax])
        xlabel('L');
        ylabel('|D|');
        
    else
        surf(T.actualTrainingSamples_vec,T.sigma_w_2_vec,T.lls./T.refs)
        xlim([T.actualTrainingSamples_vec(1),T.actualTrainingSamples_vec(end)])
        ylim([T.sigma_w_2_vec(1),T.sigma_w_2_vec(end)])
    end
    zlabel('Normalised Variance')
    caxis([currMin,currMax])
    set(gca,'FontSize',20)
    set(gca,'LineWidth',2.0)
    colormap('coolwarm')
    colorbar
    
    hold off
end