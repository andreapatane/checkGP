function open_region = fathoming(min_or_max,c_lb,ubstar,c_ub,lbstar,debug)
if strcmp(min_or_max,'min')
    if c_lb > ubstar
        open_region = false;
%         if debug
%             disp('killed in fathoming')
%         end
    else
        open_region = true;
    end
    
elseif strcmp(min_or_max,'max')
    if c_ub < lbstar
        open_region = false;
%         if debug
%             disp('killed in fathoming')
%         end
    else
        open_region = true;
    end
end
end