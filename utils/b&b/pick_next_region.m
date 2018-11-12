function [flagBreak,c_r_i_cell,ubs,lbs,region,x_overs,curr_x_over,lb_b4_split,ub_b4_split] = pick_next_region(lbstar,ubstar,region,lbs,ubs,x_overs,min_or_max,debug)

%best-bound-first
count = 0;
if strcmp(min_or_max,'max')
    aux = -inf;
    while aux < lbstar
        if isempty(region)
            flagBreak = true;
            break;
        else
            flagBreak = false;
        end
        %[~,idx] = max(cell2mat(lbs));
        [~,idx] = max(cell2mat(ubs));
        lb_b4_split = lbs{idx};
        ub_b4_split = ubs{idx};
        c_r_i_cell = splitregion(region{idx});
        curr_x_over = x_overs{idx};
        aux = ubs{idx};
        %disp(lbstar - aux)
        region(idx) = [];
        ubs(idx) = [];
        lbs(idx) = [];
        x_overs(idx) = [];
        if debug && (count > 1)
            disp('killed in pick')
        end
        count = count + 1;
        
    end
elseif strcmp(min_or_max,'min')
    aux = inf;
    while aux > ubstar
        if isempty(region)
            flagBreak = true;
            break;
        else
            flagBreak = false;
        end
        [~,idx] = min(cell2mat(lbs));
        lb_b4_split = lbs{idx};
        ub_b4_split = ubs{idx};
        %[~,idx] = max(cell2mat(lbs));
        %[~,idx] = min(cell2mat(ubs));
        %disp(['lbs: ', num2str((cell2mat(lbs)))])
        %disp(['ubs: ', num2str((cell2mat(ubs)))])
        c_r_i_cell = splitregion(region{idx});
        curr_x_over = x_overs{idx};
        aux = lbs{idx};
        region(idx) = [];
        ubs(idx) = [];
        lbs(idx) = [];
        %disp('After picking')
        %disp(['lbs: ', num2str((cell2mat(lbs)))])
        %disp(['ubs: ', num2str((cell2mat(ubs)))])
        x_overs(idx) = [];
        if debug && (count > 1)
            disp('killed in pick')
        end
        count = count + 1;
    end
    
    
    
end




%depth-first
%c_r_i_cell = splitregion(region{end});
%region = region(1:(end-1));

%breadth-first
%c_r_i_cell = splitregion(region{1});
%region = region(2:end);



end


function c_r_i_cell = splitregion(reg)

[~,split_idx] =  max(reg(2,:) - reg(1,:));

%split_idx = randi(size(reg,2));
mb = 0.5*(reg(1,split_idx) + reg(2,split_idx));

c_r_i_cell{1} = reg;
c_r_i_cell{2} = reg;
c_r_i_cell{1}(2,split_idx) = mb;
c_r_i_cell{2}(1,split_idx) = mb;


end