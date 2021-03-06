clear; 

c = [1,9,1];
a_i = [1,2,3,1,0;
    3,2,2,0,1;];
b = [9;
    15];

% initial basic variables, the first basic variables are
% the slack variables ie x_3 and x_4
basic = [4,5];

% initial non basic variables are the regular decision variables
non_basic = [1,2,3];

% Calculate basis matrix





% function optTest will have multiple returns. 
% [1 if an optimal sol'n is found (if not 0),
% the coeffs of the non_basic variables,
% the entering variable]
function [optimal,coeffs, entering] = optTest(a,c_basis,b_inv, b, c)
    non_basic_coeffs = (c_basis*b_inv*a) -c ;
    min_coeff = min(non_basic_coeffs);
    min_coeff_index = find(non_basic_coeffs == min_coeff);
    if all(non_basic_coeffs >= 0)
        optimal = 1;
    else
        optimal = 0;
    end
    coeffs = non_basic_coeffs;
    entering = min_coeff_index;
    
end

% function minRatTest will return the leaving basic variable, type
% 1x1 integer. 
function leaving_var = minRatTest(entering_basic_vec,b)
    ratios = [];
    for i = 1:size(b):
        if b(i) > 0 && entering_basic_vec(i) > 0
            ratios(i) = b(i)/entering_basic_vec(i);
        end
    end
    leaving_var = find(entering_basic_vec==min(
    
            

end
