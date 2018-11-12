function [lbstar,ubstar] = b_n_b_for_lb_computation(mult_coefficients,x_L,x_U,...
    training_data,outputIdx,testPoint,mode,test_prediction,r_star,R,upper_bound_comp_opts)
%main function for branch and bound. mode string indicates what quantity I
%want to optimize for. This will decide what function to use in order to
%compute bounds, and will decide if it is a maximization or minimization
%problem

debug = true;
tic

%deciding if I want maximum or minimum of this functions...
if strcmp(mode,'mu')
    min_or_max = 'min';
    maxIterations = upper_bound_comp_opts.maxIterations_mu;
elseif strcmp(mode,'xi')
    min_or_max = 'max';
    maxIterations = upper_bound_comp_opts.maxIterations_xi;
elseif strcmp(mode,'sigma')
    min_or_max = 'max';
    maxIterations = upper_bound_comp_opts.maxIterations_sigma;
end


%epsilon here defined is the global error tolerance on the optimum value
if strcmp(mode,'mu')
    epsilon = upper_bound_comp_opts.tollerance_mu;
elseif strcmp(mode,'xi')
    epsilon = upper_bound_comp_opts.tollerance_xi;
elseif strcmp(mode,'sigma')
    % in this case I use a relative tolerance. This is because variance
    % order of magnitude varies a lot depending on the number of training
    % points used. So this will ensure a fair comparison in Figure 4
    epsilon  = upper_bound_comp_opts.tollerance_sigma_percentage*upper_bound_comp_opts.ref_value;
end

%If I am using relu kernel, than I need both the angular and cartesian coordinates of
%training and test point.
if strcmp(upper_bound_comp_opts.kernel,'ReLU')
    training_data_temp.spherical = training_data;
    testPoint_temp.spherical = testPoint;
    training_data_temp.cartesian = get_cartesian_coordinates(training_data_temp.spherical,[]);
    testPoint_temp.cartesian = get_cartesian_coordinates(testPoint_temp.spherical,[]);
    
    training_data = training_data_temp;
    testPoint = testPoint_temp;
    clear training_data_temp testPoint_temp
elseif strcmp(upper_bound_comp_opts.kernel,'sqe')
    training_data_temp.cartesian  = training_data;
    testPoint_temp.cartesian  = testPoint;
    
    training_data = training_data_temp;
    testPoint = testPoint_temp;
    clear training_data_temp testPoint_temp
end


%the pixels I want to modify are found by checking which hyper-rectangle
%side is not trivial
pixels2modify = find(x_U > x_L);
%and this is used to introduce compact representation for the regions (as trivial sides are easily reconstructed from x_L)
region = { [x_L(pixels2modify);x_U(pixels2modify)] };

%The initial region of branch and bound is trivially the whole region [x_L,x_U]
[c_x_l,c_x_u] = expand_domain_variables(region{1},pixels2modify,x_L);

%I compute lower and upper bounds on the current region, as well as a point
%solution of the over-approximated problem (this can be iteratively used to start off search again in the same region...)
[lbstar,ubstar,x_over_star] = get_lb_and_ub_in_current_region(mode,upper_bound_comp_opts,mult_coefficients,c_x_l,...
    c_x_u,training_data,outputIdx,testPoint,r_star,R,[]);


%I save current lower annd upper bounds
ubs = {ubstar};
lbs = {lbstar};
x_overs = {x_over_star};

%if lower and upper bounds are epsilon-close, then the search is over!
if lbstar > (ubstar - epsilon)
    region = {};
end



%here I start the while loop exploring the tree space associated to branch
%and bound...
iteration_count = 1;
while ~isempty(region) %until it is still worth to look inside a region...
    
    %I pick the next next region to split, and I split it in two
    %sub-regions
    %picking and branching
    [flagBreak,c_r_i_cell,ubs,lbs,region,x_overs,curr_x_over,lb_b4_split,ub_b4_split] = pick_next_region(lbstar,ubstar,region,lbs,ubs,x_overs,min_or_max,debug); 
    if flagBreak
        break;
    end
    
    %performing lower and upper bounding on the freshly generated two
    %subregions in c_r_i_cell
    for ii = 1:length(c_r_i_cell)
        iteration_count = iteration_count + 1;
        [c_x_l,c_x_u] = expand_domain_variables(c_r_i_cell{ii},pixels2modify,x_L);
        
        %computing region diameter for stats purposes...
        rs = max(c_r_i_cell{ii}(2,:) - c_r_i_cell{ii}(1,:));
        
        %bounding
        [c_lb,c_ub,x_over] = get_lb_and_ub_in_current_region(mode,upper_bound_comp_opts,mult_coefficients,c_x_l,...
            c_x_u,training_data,outputIdx,testPoint,r_star,R,curr_x_over);
        
        
        if debug
            if ~mod(iteration_count-1,50)
                strig_c_lb = sprintf('%.6f',c_lb);
                strig_c_ub = sprintf('%.6f',c_ub);
                strig_lbstar = sprintf('%.6f',lbstar);
                strig_ubstar = sprintf('%.6f',ubstar);
                string_lb_b4_split = sprintf('%.6f',lb_b4_split);
                string_ub_b4_split = sprintf('%.6f',ub_b4_split);
                string_rs = sprintf('%.6f',rs);
                
                disp(['c_lb: ', strig_c_lb, '; c_ub: ', strig_c_ub, ...
                    '; lb_b4_split: ', string_lb_b4_split ,'; ub_b4_split: ', string_ub_b4_split , '; rs: ', string_rs ...
                    '; lbstar: ',strig_lbstar,'; ubstar: ', strig_ubstar,' ; currently in stack: ', int2str(length(region))])
            end
        end
        %fathoming...
        open_region = fathoming(min_or_max,c_lb,ubstar,c_ub,lbstar,debug);
        if open_region
            %closing region if lower and upper bound are within tolerance
            if c_lb > (c_ub - epsilon)
                open_region = false;
            else
                open_region = true;
            end
            if strcmp(min_or_max,'min')
                if c_ub <= ubstar
                    ubstar = c_ub;
                    lbstar = c_lb;
                end
            elseif strcmp(min_or_max,'max')
                if c_lb >= lbstar
                    ubstar = c_ub;
                    lbstar = c_lb;
                end
            end
        end
        %if region is still worth looking at, I stack it.
        if open_region
            region{end+1} = c_r_i_cell{ii};
            ubs{end+1} = c_ub;
            lbs{end+1} = c_lb;
            x_overs{end+1} = x_over;
        end
    end
    if strcmp(min_or_max,'min')
        if  (lb_b4_split > (ubstar - epsilon))
            lbstar = max(lbstar,lb_b4_split);
            break
        else
            if ~mod(iteration_count-1,50)
                mar = sprintf('%.6f',ubstar-lb_b4_split); 
                disp(['Current margin: ', mar])
            end
        end
    end
    if strcmp(min_or_max,'max')
        if  (lbstar > (ub_b4_split - epsilon))
            ubstar = min(ubstar,ub_b4_split);
            break
        else
            if ~mod(iteration_count-1,50)
                mar = sprintf('%.6f',ub_b4_split-lbstar);
                disp(['Current margin: ', mar])
            end
        end
    end
    if iteration_count > maxIterations
        warning('Branch and bound exceeded maximum number of iteration specified. The bound provided is still valid, but most-likely exceeds the tolerance thereshold specified...')
        break;
    end
end

%doing some traslation for mean computation...
if strcmp(mode,'mu')
    temp = lbstar;
    lbstar = test_prediction(outputIdx) - ubstar;
    ubstar = test_prediction(outputIdx) - temp;
end

disp(['iteration_count: ',int2str(iteration_count)])
toc

end


function [c_x_l,c_x_u] = expand_domain_variables(reg,pixels2modify,x_L)
%expands compact representations of region to full representations

c_x_l = x_L;
c_x_u = x_L;
c_x_l(pixels2modify) = reg(1,:);
c_x_u(pixels2modify) = reg(2,:);


end


