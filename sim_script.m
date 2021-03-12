clear; 
clc;
% Charlie Penvari 
% DASE 4000 - Project 1
% Due 13th of March
% email: cpenvari@uccs.edu




% Hardcode problem from project description
c = [1,9,1];
a= [1,2,3;
    3,2,2];
b = [9;
    15];

% define the initial basic and non basic variables. The integers represent
% the sub script of the decision variable in a typical lp model.

basic = [4,5]; % basic starts as the slack variables x_4 and x_5
non_basic = [1,2,3]; % non basic starts as the original decision variables.

% solve the problem.
fprintf("Project problem\n------------------------------------------------")
handleSolve(a,b,c,basic,non_basic)

% solve another, wyndor problem
c_1 = [3,5];
a_1 = [1,0;
        0,2;
       3,2];
b_1 = [4;12;18];

basic_1 = [3,4,5];
non_basic_1 = [1,2];

fprintf("Wyndor Problem\n------------------------------------------------")
handleSolve(a_1,b_1,c_1,basic_1,non_basic_1)
%------------------------------------------------------


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
        basis = generateBasisMatrix(a_aug, basic_local,non_basic_local);
    
        %invert basis matrix and calculate x_b
        inv_basis = inv(basis);
        x_b = inv_basis*b;
        
        %calculate c_basis, the objective func coeffs of the basis vectors
        c_b = extractCBasis(basic_local, c_aug);
        
        %check for optimal soln, find entering variable. 
        [op_soln_found, entering_var_index] = optTest(a, c_b, inv_basis, c);
        %Note, entering_var_index refers to the index of the entering
        %variable in the non-basic set
    
        %if optimal soln was not found above, extract the basis vec from
        %'a'and calculate the exiting variable
        entering_basis_vector = extractEnteringBasisVec(a_aug, entering_var_index);
        exiting_var_index = minRatTest(entering_basis_vector, x_b);
        %Note that the exiting_var_index refers to the index of the basic
        %set.
        
        %update variable set
        basic_loc_string =sprintf('%.0f ',basic_local);
        n_basic_loc_string = sprintf('%.0f ',non_basic_local);
        fprintf("\nBasic variables: x_%s Non-basic variables: x_%s", ...
                basic_loc_string,n_basic_loc_string)
        basic_last_local = basic_local;
        basic_local = updateBasicSet(basic_local, exiting_var_index, non_basic_local(entering_var_index));
        non_basic_local = updateNonBasicSet(non_basic_local, entering_var_index, basic_last_local(exiting_var_index));
        z = c_b*inv_basis*b;
        fprintf("\nEntering basic variable: x_%.0f Leaving basic variable: x_%.0f\n", ...
                [basic_local(exiting_var_index),...
                non_basic_local(entering_var_index)])
        fprintf("Current cpf %.2f\n",[z])
        
    end
    fprintf("\nEnd simplex loop\n")
    fprintf("Optimal value is Z = %.3f\n",z)
    % format-print the variables and values of the optimal solution.
    fprintf("The value occurs at the variable x_%.0f = %.3f\n",[basic_last_local.',x_b].')
    fprintf("All other variables are presumed to equal zero\n")

    
end



% function optTest will have multiple returns. 
% [1 if an optimal sol'n is found (if not 0),
% the coeffs of the non_basic variables,
% the entering variable]
function [optimal, entering_var_ind] = optTest(a,c_basis,inv_b, c)
    non_basic_coeffs = (c_basis*inv_b*a) -c ;
    other_coeffs = c_basis*inv_b;
    
    % the find the min coeff and find its index in the non basic set.
    min_coeff = min(non_basic_coeffs);
    min_coeff_index = find(non_basic_coeffs == min_coeff);
    
    % tie break if applicable
    if size(min_coeff_index,2) >1
        min_coeff_index = min_coeff_index(1);
    end
    % optimality test
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
        else
            ratios(i,1) =0;
        end
    end
    
    %tie break 
    leaving_var_ind_temp = find(ratios==min(ratios(ratios>0)));
    if size(leaving_var_ind_temp,1) >1
        leaving_var_ind = leaving_var_ind_temp(1)
    else
    leaving_var_ind = leaving_var_ind_temp;
    end
end

% function generateBasisMatrix will generate the basis matrix which is
% essentially the augmented a matrix without the non_basic vectors.
% Returns type size(a) x size(a-non_basic cols)
function basis = generateBasisMatrix(a_aug, basic,non_basic)
    
    temp_new_basis = zeros(size(a_aug,1),size(a_aug,2)-(size(non_basic,2)));
    
    for i= 1:size(basic,2)
        temp_new_basis(:,i) = a_aug(:,basic(i));
    end

    basis = temp_new_basis;
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

% function updateBasicSet will replace the exiting variable from the basic
% set with the entering variable from the non-basic set.
function basic_set = updateBasicSet(old_basic_set,replacement_ind,new_val)
    temp = old_basic_set;
    temp(replacement_ind) = new_val;
    basic_set = temp;

end

% function updateNonBasicSet will replace the entering variable that is in
% the non basic set and replace it with the exiting variable.
function non_basic_set = updateNonBasicSet(old_non_basic_set, replacement_ind, new_val)
    temp = old_non_basic_set;
    temp(replacement_ind) = new_val;
    non_basic_set = temp;
end