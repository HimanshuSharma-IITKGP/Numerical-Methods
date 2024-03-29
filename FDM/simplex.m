clear
clc
close all




% coefficients of decision variables for less than equal to and equal to constraints 
a_matrix = [1 -1 1 3 2;-1 1 0 1 -1; 2 -2 1 0 0;1 -1 2 1 2];
% positive rhs
b_vector = [4; 1; 2; 2];
% row vector containing the indices of the equations which are equal to constraint
row_index_for_equality_constraint = [4, -1] ;
% coefficients of the z = f(x)
objective_function_coefficients_row_vector = [-2 2 1 -4 3] ;



[basis, rhs] = optimize( a_matrix, b_vector, objective_function_coefficients_row_vector, row_index_for_equality_constraint ) ;
disp(rhs)
disp(basis)


function [basis_variables_index_vector, rhs_vector] = optimize(a_matrix, b_vector, objective_function_coefficients_row_vector, row_index_for_equality_constraint)

    [coeffiecients_matrix, rhs_vector, basis_variables_index_vector] = convertToAugmentedForm(a_matrix, b_vector, objective_function_coefficients_row_vector, row_index_for_equality_constraint);
%     disp(coeffiecients_matrix);
%     disp(rhs_vector);
%     disp(basis_variables_index_vector) ;

    number_of_equations = size(coeffiecients_matrix, 1);
    number_of_variables = size(coeffiecients_matrix, 2); % it does not includes Z as a variable
    

    
    while(true)
    
        % finding the entering basic variable
        [a, j] = min(coeffiecients_matrix(1, :));
    
    
        % check for no better solution
        if(a>=0)
            break;
        end
    
        


        % minimum ratio test
        c = inf;
        for i=2:number_of_equations
            if(coeffiecients_matrix(i, j)>0)
                if(rhs_vector(i, 1)/coeffiecients_matrix(i, j) < c) 
                    c = rhs_vector(i, 1)/coeffiecients_matrix(i, j);
                end
            end
        end
    
        % check for no leaving basic variable (unbounded solution)
        if(c == inf)
            printf('unbounded solution')
            break ;
        end
    
        % index_of_leaving_bv_in_vector is the index of the leaving basic  variable in the vector of basic variables
        index_of_leaving_bv_in_vector = find(rhs_vector(2:number_of_equations, 1)./coeffiecients_matrix(2:number_of_equations, j) == c, 1) ;
        % index_of_leaving_bv is the subscript of leaving basic variable
        index_of_leaving_bv = basis_variables_index_vector(index_of_leaving_bv_in_vector) ;
        % index_of_entering_bv is the subscript of the entering basic variable
        index_of_entering_bv = j ;
    
    
    
        pivot_element = coeffiecients_matrix(1+index_of_leaving_bv_in_vector, index_of_entering_bv) ;

        coeffiecients_matrix(1+index_of_leaving_bv_in_vector, :) = coeffiecients_matrix(1+index_of_leaving_bv_in_vector, :)/pivot_element;
        rhs_vector(1+index_of_leaving_bv_in_vector, 1) = rhs_vector(1+index_of_leaving_bv_in_vector, 1)/pivot_element ;
    
        % row transformation
        for i = 1:number_of_equations
            if(i == 1+index_of_leaving_bv_in_vector)
                continue ;
            end
    
            rhs_vector(i, 1) = rhs_vector(i, 1) - coeffiecients_matrix(i, index_of_entering_bv)*rhs_vector(1+index_of_leaving_bv_in_vector, 1) ;
            coeffiecients_matrix(i, :) = coeffiecients_matrix(i, :) - coeffiecients_matrix(i, index_of_entering_bv)*coeffiecients_matrix(1+index_of_leaving_bv_in_vector, :) ;
    
        end
    
        % updating the list of basic variables
        basis_variables_index_vector(index_of_leaving_bv_in_vector) = index_of_entering_bv;
        
    end
    
%     disp(basis_variables_index_vector)
%     disp(rhs_vector)

end




%  handles the equal to constraints and converts the problem to augmented form
function [coeffiecent_matrix, rhs_vector, basis_variables_index_vector] = convertToAugmentedForm(a_matrix, b_vector, objective_function_coefficients_row_vector, row_index_for_equality_constraint)

    number_of_inequalities = size(a_matrix, 1);
    number_of_dv = size(a_matrix, 2) ; % number of decision variables
    M = 10000000 ;


    % this step completes the coeffcient_matrix except the z row (first row)
    coeffiecent_matrix = [a_matrix, eye(number_of_inequalities)];

    % a is the first row of coefficient matrix
    a = -[objective_function_coefficients_row_vector, zeros(1, number_of_inequalities)] ;

    % making the coefficient in first row as M for equality constraints
    for i = 1:size(row_index_for_equality_constraint, 2)
        if(row_index_for_equality_constraint(i)>0)
            a(row_index_for_equality_constraint(i) + number_of_dv) = M ;
        end
    end 
    coeffiecent_matrix = [a; coeffiecent_matrix] ;
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