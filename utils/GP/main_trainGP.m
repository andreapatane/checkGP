function [sigma_w_2,sigma_b_2,test_data,test_norm,test_labels,mu_test,kernel,...
    Kstar,training_data,trainedSystem,K_inv,actualTrainingSamples,R,std_Y,sigma_test]=main_trainGP(dataSet,depth,sigma_w_2,actualTrainingSamples)




if strcmp(dataSet,'x_products')
    kernel = 'sqe';
    mle_FLag = true;
elseif strcmp(dataSet,'mnist')
    kernel = 'N-ReLU';
    mle_FLag = false;
end

trainSamplesFile = 2000;
if strcmp(dataSet,'mnist')
    if nargin < 4
        actualTrainingSamples = 1000;
    end
    testSamples = 200;
    plotPosteriori = false;
elseif strcmp(dataSet,'x_products')
    plotPosteriori = false;
    actualTrainingSamples = 128;
    if plotPosteriori
        testSamples = 625;
    else
        testSamples = 2;
    end
end
pool = true;
cache_K_inv = false;





data_folder = get_data_folder();

if strcmp(dataSet,'mnist')
    training_data = csvread([data_folder,'x_train',int2str(trainSamplesFile),'.csv']);
    training_labels = csvread([data_folder,'y_train',int2str(trainSamplesFile),'.csv']);
    training_data = training_data(1:actualTrainingSamples,:);
    training_labels = training_labels(1:actualTrainingSamples,:);
    [training_data,nPixels]  = averagePooling(training_data);
    input_dim = nPixels^2;
    output_dim = 10;
    test_data = csvread([data_folder,'x_test.csv']);
    test_data = test_data(1:testSamples,:);
    test_data = averagePooling(test_data);
    test_labels = csvread([data_folder,'y_test.csv']);
    test_labels = test_labels(1:testSamples);
    std_Y = [];
elseif strcmp(dataSet,'x_products')
    [training_data,training_labels,test_data,test_labels,input_dim,~,std_Y] =  x_product_data(actualTrainingSamples,testSamples,plotPosteriori);
end

if nargin < 3
    sigma_w_2 = 3.19;
end
sigma_b_2 = 0.00;

sig_eps = 1e-10;





if strcmp(kernel,'sqe')
    sigma_b_2 =  sigma_b_2 * ones(input_dim,1);
    test_norm = [];
end

if strcmp(kernel,'N-ReLU')
    training_data = to_norm1(training_data);
    [test_data,test_norm] = to_norm1(test_data);
end

disp('training the GP...')
[trainedSystem,R,diags_K,K_inv,sigma_w_2,sigma_b_2] = trainGP(training_data,training_labels,sigma_w_2,sigma_b_2,depth,kernel,sig_eps,cache_K_inv,dataSet,mle_FLag);

disp('Computing kernels on test set')
[Kstar,Kstarstar] = compute_kernel_on_test(test_data,training_data,kernel,depth,sigma_w_2,sigma_b_2,diags_K);

disp('finally making predictions')
[mu_test, sigma_test,y_test_hat] = make_predictions_on_test(Kstar,trainedSystem,Kstarstar,R);

if strcmp(dataSet,'mnist')
    test_accuracy = mean(y_test_hat ==  test_labels);
    disp(test_accuracy)
    goodIdxs = find(y_test_hat ==  test_labels);
    test_labels = test_labels(goodIdxs);
    test_data = test_data(goodIdxs,:);
    test_norm = test_norm(goodIdxs);
    mu_test = mu_test(goodIdxs,:);
    Kstar = Kstar(goodIdxs,:);
    
    
elseif strcmp(dataSet,'x_products')
    rmse = sqrt(  mean( (mu_test - test_labels).^2   )  );
    disp(rmse)
end

if strcmp(dataSet,'x_products')
    if plotPosteriori
        tt = linspace(-3.5,3.5,sqrt(testSamples));
        mu_4_plot = zeros(length(tt),length(tt));
        var_4_plot = zeros(length(tt),length(tt));
        
        for ii = 1:length(tt)
            for jj = 1:length(tt)
                mu_4_plot(ii,jj) = mu_test(  (ii-1)*length(tt) + jj  );
            end
        end
        for ii = 1:length(tt)
            for jj = 1:length(tt)
                var_4_plot(ii,jj) = sigma_test(  (ii-1)*length(tt) + jj, (ii-1)*length(tt) + jj );
            end
        end
        figure
        hold on
        h = surf(tt,tt,mu_4_plot);
        xlim([-3.5,3.5])
        ylim([-3.5,3.5])
        %set(h,'edgecolor','none')
        colormap(coolwarm(64))
        colorbar
        scatter3(training_data(:,1),training_data(:,2),max(max(mu_4_plot))*ones(size(training_labels)),'k','x','LineWidth',1.5)
        set(gca,'FontSize',16)
        xlabel('x_1')
        ylabel('x_2')
        hold off
        
        figure
        hold on
        surf(tt,tt,var_4_plot);
        xlim([-3.5,3.5])
        ylim([-3.5,3.5])
        colormap(coolwarm(64))
        scatter3(training_data(:,1),training_data(:,2),max(max(var_4_plot))*ones(size(training_labels)),'k','x','LineWidth',1.5)
        colorbar
        xlabel('x_1')
        ylabel('x_2')
        set(gca,'FontSize',16)
        
        hold off
        
        figure
        surf(tt,tt,var_4_plot./mu_4_plot);
        colormap(coolwarm(64))
        colorbar
    end
end


end



