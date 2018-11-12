function [mu_test, sigma_test, y_test_hat] = make_predictions_on_test(Kstar,trainedSystem,Kstarstar,R,output_mode)

if ~isempty(trainedSystem)
    mu_test = Kstar*trainedSystem;
else
    mu_test = [];
end
sigma_test = Kstarstar - Kstar * ( R\(R'\(Kstar')));

if strcmp(output_mode,'class_via_regress')
    [~, y_test_hat] = max(mu_test,[],2);
    y_test_hat = y_test_hat - 1;
else
    y_test_hat = [];
end

end