function saveResults(p_ub,sup_mu,sup_xi,pixel2Modify_cell,upper_bound_comp_opts,startIdx,endIdx,nameTempl,method,result_folder)

T.upper_bound_comp_opts = upper_bound_comp_opts;
T.testPointIdxs = startIdx:endIdx;
T.pixel2Modify_cell = pixel2Modify_cell;
T.p_ub = p_ub;
T.sup_mu = sup_mu;
T.sup_xi = sup_xi;
if strcmp(method,'formal')
    nameTempl = ['hat_',nameTempl];
elseif strcmp(method,'approximated')
    nameTempl = ['bar_',nameTempl];
else
    nameTempl = ['---',nameTempl];
end

if ~upper_bound_comp_opts.invariance
    save([result_folder,'phi_1_',nameTempl,'.mat'],'T');
else
    save([result_folder,'phi_2_',nameTempl,'.mat'],'T');
end

end