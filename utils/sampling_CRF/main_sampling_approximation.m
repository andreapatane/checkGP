function stat_probs = main_sampling_approximation(nsamples,test_data,sigma_test,gp_training_opts,delta,invariance,mu_test)
%main_sampling_approximation compute the values of \bar{\phi}_1 and
%\bar{\phi}_2, that is sampling based approximation of safety and
%invariance. The function relies on the  Random Field Simulation  Toolbox that has to be in the
%working path for the function to work (can be found at https://uk.mathworks.com/matlabcentral/fileexchange/27613-random-field-simulation)





test_data = test_data(:,1:2);

corr.C = sigma_test;
corr.name = 'gauss';
corr.c0 = 1;
corr.sigma = 1;
[F,~] = randomfield(corr,test_data,'mean',mu_test,'nsamples',nsamples);


tt_x = linspace(gp_training_opts.extremes(1,1),gp_training_opts.extremes(1,2),sqrt( size(test_data,1) ));
tt_y = linspace(gp_training_opts.extremes(2,1),gp_training_opts.extremes(2,2),sqrt( size(test_data,1) ));




stat_probs = zeros(size(delta));
for kk = 1:nsamples
    samp_4_plot = zeros(length(tt_x),length(tt_y));
    
    
    for ii = 1:length(tt_x)
        for jj = 1:length(tt_y)
            samp_4_plot(ii,jj) = real(F(  (ii-1)*length(tt_y) + jj ,kk ));
        end
    end
    
    

    initialValue = samp_4_plot(ceil(length(tt_x)/2),ceil(length(tt_y)/2));
    for jj = 1:length(delta)
        valVec = initialValue - samp_4_plot(:);
        if invariance
            valVec = abs(valVec);
        end
        if any( valVec  >  delta(jj))
            stat_probs(jj) = stat_probs(jj) + 1;
        end
    end
    
    
end





end