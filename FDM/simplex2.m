clear
clc
close all

% coefficients of decision variables for less than equal to and equal to constraints 
a_matrix = [0 0.1 ; 0.25 0; 3 2];
% positive rhs
b_vector = [200 ; 800; 12000];
% row vector containing the indices of the equations which are equal to constraint
row_index_for_equality_constraint = [-1, -1] ;
% coefficients of the z = f(x)
objective_function_coefficients_row_vector = [0.2 0.1] ;



function [] = twoPhaseSimplex(a_matrix, b_vector, objective_function_coefficients_row_vector, row_index_for_equality_constraint)
    [coefficients_matrix, rhs_vector, z1, z2] = convertToAugmentedForm(a_matrix, b_vector, objective_function_coefficients_row_vector, row_index_for_equality_constraint);
    
    [] = simplex(coefficients_matrix, rhs_vector);


end


function [basis_variables_index_vector, rhs_vector] = simplex(coefficients_matrix, rhs_vector)

%     disp(coefficients_matrix);
%     disp(rhs_vector);
%     disp(basis_variables_index_vector) ;

    number_of_equations = size(coefficients_matrix, 1);
    number_of_variables = size(coefficients_matrix, 2); % it does not includes Z as a variable
    

    
    while(true)
    
        % finding the entering basic variable
        [a, j] = min(coefficients_matrix(1, :));
    
    
        % check for no better solution
        if(a>=0)
            break;
        end
    
        


        % minimum ratio test
        c = inf;
        for i=2:number_of_equations
            if(coefficients_matrix(i, j)>0)
                if(rhs_vector(i, 1)/coefficients_matrix(i, j) < c) 
                    c = rhs_vector(i, 1)/coefficients_matrix(i, j);
                end
            end
        end
    
        % check for no leaving basic variable (unbounded solution)
        if(c == inf)
            printf('unbounded solution')
            break ;
        end
    
        % index_of_leaving_bv_in_vector is the index of the leaving basic  variable in the vector of basic variables
        index_of_leaving_bv_in_vector = find(rhs_vector(2:number_of_equations, 1)./coefficients_matrix(2:number_of_equations, j) == c, 1) ;
        % index_of_leaving_bv is the subscript of leaving basic variable
        index_of_leaving_bv = basis_variables_index_vector(index_of_leaving_bv_in_vector) ;
        % index_of_entering_bv is the subscript of the entering basic variable
        index_of_entering_bv = j ;
    
    
    
        pivot_element = coefficients_matrix(1+index_of_leaving_bv_in_vector, index_of_entering_bv) ;

        coefficients_matrix(1+index_of_leaving_bv_in_vector, :) = coefficients_matrix(1+index_of_leaving_bv_in_vector, :)/pivot_element;
        rhs_vector(1+index_of_leaving_bv_in_vector, 1) = rhs_vector(1+index_of_leaving_bv_in_vector, 1)/pivot_element ;
    
        % row transformation
        for i = 1:number_of_equations
            if(i == 1+index_of_leaving_bv_in_vector)
                continue ;
            end
    
            rhs_vector(i, 1) = rhs_vector(i, 1) - coefficients_matrix(i, index_of_entering_bv)*rhs_vector(1+index_of_leaving_bv_in_vector, 1) ;
            coefficients_matrix(i, :) = coefficients_matrix(i, :) - coefficients_matrix(i, index_of_entering_bv)*coefficients_matrix(1+index_of_leaving_bv_in_vector, :) ;
    
        end
    
        % updating the list of basic variables
        basis_variables_index_vector(index_of_leaving_bv_in_vector) = index_of_entering_bv;
        
    end
    
%     disp(basis_variables_index_vector)
%     disp(rhs_vector)

end




%  handles the equal to constraints and converts the problem to augmented form
function [cm_1, cm_2, rhs_vector] = convertToAugmentedForm(a_matrix, b_vector, objective_function_coefficients_row_vector, row_index_for_equality_constraint)

    number_of_inequalities = size(a_matrix, 1);
    number_of_dv = size(a_matrix, 2) ; % number of decision variables
    M = 10000000 ;


    % this step completes the coeffcient_matrix except the z row (first row)
    cm_1 = [a_matrix, eye(number_of_inequalities)];
    cm_2 = cm_1 ;

    % a is the first row of coefficient matrix
    a1 = -[objective_function_coefficients_row_vector, zeros(1, number_of_inequalities)] ;
    a2 = a1 ;
    cm_1 = [a1; a_matrix] ;

    % making the coefficient in first row as M for equality constraints
    for i = 1:size(row_index_for_equality_constraint, 2)
        if(row_index_for_equality_constraint(i)>0)
            a(row_index_for_equality_constraint(i) + number_of_dv) = M ;
        end
    end 
    
 

    for i = 1:size(z2, 2)
        if(a(i) == M)
            a2(i) = 1;
        end
    end
    cm_2 = [a2; a_matrix] ;

    
    rhs_vector = [0; b_vector];

    % making the coefficient of basic variable as zero in first row by row transformations
    for i = 1:size(row_index_for_equality_constraint, 2)
        if(row_index_for_equality_constraint(i)>0)
            coeffiecent_matrix(1, :) = coeffiecent_matrix(1, :) - M*coeffiecent_matrix(row_index_for_equality_constraint(i) + 1, :);
            rhs_vector(1) = rhs_vector(1) - M*rhs_vector(row_index_for_equality_constraint(i) + 1) ;
        end
    end 
    
    % index of basic variables in the basis 
    basis_variables_index_vector = ((size(a_matrix, 2)+1):(size(a_matrix, 2) + number_of_inequalities))' ;

end