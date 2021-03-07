clear; 
clc;

% define initial matrix
c = [1,9,1];
a= [1,2,3;
    3,2,2];
b = [9;
    15];

% define the initial basic and non basic variables. The integers represent
% the sub script of the decision variable in a typical lp model.

basic = [4,5]; % basic starts as the slack variables x_4 and x_5
non_basic = [1,2,3]; % non basic starts as the original decision variables.

handleSolve(a,b,c,basic,non_basic)

% function handleSolve will take in the lp model and the initial basic and
% non basic variables and will print out the optimal solution.
function handleSolve(a, b, c, basic, non_basic)
    % create augmented matrix of a.
    i = eye(size(a, 1));
    a_aug = [a i];
    
    % create augmented row vector of c.
    c_aug = [c zeros(1,size(a, 1))];
    
    %initialize copies of the set of basic/non-basic variables.
    basic_local = basic;
    non_basic_local = non_basic;
    
    op_soln_found = 0; % bool type, default to 0 (false)
    
    %loop until optimal solution is found. 
    while op_soln_found ~= 1
        basis = generateBasisMatrix(a_aug, non_basic_local);
    
        %invert b matrix and calculate x_b
        inv_b = inv(basis);
        x_b = inv_b*b;
        
        %calculate c_basis, the objective func coeffs of the basis vectors
        c_b = extractCBasis(basic_local, c_aug);
        
        %check for optimal soln, find entering variable. 
        [op_soln_found, entering_var_index] = optTest(a, c_b, inv_b, c);
    
        %if optimal soln was not found above, calculate the exiting variable
        entering_basis_vector = extractEnteringBasisVec(a_aug, entering_var_index);
        exiting_var_index = minRatTest(entering_basis_vector, x_b);
        
        %update variable set
        basic_last_local = basic_local;
        basic_local = updateBasicSet(basic_local, exiting_var_index, non_basic_local(entering_var_index));
        non_basic_local = updateNonBasicSet(non_basic_local, entering_var_index, basic_last_local(exiting_var_index));

    end
    % calculate z
    z = c_b*inv_b*b;
    fprintf("Optimal value is Z = %f\n",z)
    % format-print the variables and values of the optimal solution.
    fprintf("The value occurs at the variable x_%.0f = %f respectively\n",[basic_last_local.',x_b].')
    fprintf("All other variables are presumed to equal zero\n")

    
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
            ratios(i,1) = b(i)/entering_basis_vec(i);
        end
    end
    %leaving variable ind
    leaving_var_ind = find(ratios==min(ratios));
    
end

% function generateBasisMatrix will generate the basis matrix which is
% essentially the augmented a matrix without the non_basic vectors.
% Returns type size(a) x size(a-non_basic cols)
function basis = generateBasisMatrix(a_aug, non_basic)
    temp_b = a_aug;
    % strip off columns from a_aug of non_basic vectors
    non_basic = sort(non_basic);
    temp_b(:,non_basic)= [];
    basis = temp_b;
end

% function extractCBasis will extract the coeffs from the objective
% function (c matrix) of the basic variables. Return type row matrix.
function c_basis = extractCBasis(basic_vars,c)
    temp_c = zeros(1,length(basic_vars));
    for i = 1:length(temp_c)
        temp_c(i) = c(basic_vars(i));
    end
    c_basis = temp_c;

end

% function extractEnteringBasisVec will return a the column from the a 
% matrix given an entering column. Type 2x1 vec. 
function entering_basis_vec = extractEnteringBasisVec(a, entering_column)
    entering_basis_vec = a(:,entering_column);
end

% function updateBasicSet will replace the exiting variable with the
% entering variable. Type row vector.
function basic_set = updateBasicSet(old_basic_set,replacement_ind,new_val)
    temp = old_basic_set;
    temp(replacement_ind) = new_val;
    basic_set = temp;

end

% function updateNonBasicSet will take out the entering variable of the 
% non basic set and will put the exiting variable back into the non basic
% set.
function non_basic_set = updateNonBasicSet(old_non_basic_set, replacement_ind, new_val)
    temp = old_non_basic_set;
    temp(replacement_ind) = new_val;
    non_basic_set = temp;
end