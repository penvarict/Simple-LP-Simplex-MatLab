clear; 
clc;

c_aug = [1,9,1,0,0,0];
c = [1,9,1];
a_aug = [1,2,3,1,0;
    3,2,2,0,1;];
a= [1,2,3;
    3,2,2];
b = [9;
    15];
basic = [4,5];
non_basic = [1,2,3];

op_soln_found = 0;
iterations = 0
while op_soln_found ~= 1
    %generate the basis matrix with the initial set of non_basic vars
    basis = generateBasisMatrix(a_aug, non_basic);
    
    %invert b matrix and calculate x_b
    inv_b = inv(basis);
    x_b = inv_b*b;
    
    %calculate c_basis, the objective func coeffs of the basis vectors
    c_b = extractCBasis(basic,c_aug);
    %check for optimal soln, find entering variable
    [op_soln_found, entering_var_index] = optTest(a,c_b,inv_b,c);
    
    %if optimal soln was not found above, calculate the exiting variable
    entering_basis_vector = extractEnteringBasicVec(a_aug, entering_var_index);
    exiting_var_index = minRatTest(entering_basis_vector,x_b);
    
    %update variable set
    basic_last = basic
    basic = updateBasicSet(basic,exiting_var_index,non_basic(entering_var_index))
    non_basic = updateNonBasicSet(non_basic,entering_var_index,basic_last(exiting_var_index))
    leaving_var_val = non_basic(entering_var_index)
    exiting_var_val = basic(exiting_var_index)
    iterations =iterations +1;

end
z = c_b*inv_b*b;
fprintf("Optimal value is Z = %f\n",z)
fprintf("The value occurs at the variable x_%.0f = %f respectively\n",[basic_last.',x_b].')
fprintf("All other variables are presumed to equal zero\n")


% function handle solve
function handleSolve()


end



% function optTest will have multiple returns. 
% [1 if an optimal sol'n is found (if not 0),
% the coeffs of the non_basic variables,
% the entering variable]
function [optimal, entering_var_ind] = optTest(a,c_basis,inv_b, c)
    non_basic_coeffs = (c_basis*inv_b*a) -c ;
    other_coeffs = c_basis*inv_b;
    min_coeff = min(non_basic_coeffs);
    min_coeff_index = find(non_basic_coeffs == min_coeff);
    if all(non_basic_coeffs >= 0) && all(other_coeffs >=0)
        optimal = 1;
    else
        optimal = 0;
    end
    
    entering_var_ind = min_coeff_index;
    
end

% function minRatTest will return the index of the leaving basic variable, type
% 1x1 integer. 
function leaving_var_ind = minRatTest(entering_basis_vec,b)
    ratios = [];
    for i = 1:size(b,1)
        if b(i) > 0 && entering_basis_vec(i) > 0
            ratios(i,1) = b(i)/entering_basis_vec(i)
        end
    end
    %leaving variable ind
    leaving_var_ind = find(ratios==min(ratios))
    
end


function basis = generateBasisMatrix(a_aug, non_basic)
    temp_b = a_aug;
    % strip off columns from a_aug of non_basic vectors
    non_basic = sort(non_basic);
    temp_b(:,non_basic)= []
    basis = temp_b;
end

function c_basis = extractCBasis(basic_vars,c)
    temp_c = zeros(1,length(basic_vars));
    for i = 1:length(temp_c)
        temp_c(i) = c(basic_vars(i));
    end
    c_basis = temp_c;

end

function entering_basis_vec = extractEnteringBasicVec(a, entering_column)
    entering_basis_vec = a(:,entering_column)
end

function basic_set = updateBasicSet(old_basic_set,replacement_ind,new_val)
    temp = old_basic_set;
    temp(replacement_ind) = new_val;
    basic_set = temp

end

function non_basic_set = updateNonBasicSet(old_non_basic_set, replacement_ind, new_val)
    temp = old_non_basic_set;
    temp(replacement_ind) = new_val;
    non_basic_set = temp
end