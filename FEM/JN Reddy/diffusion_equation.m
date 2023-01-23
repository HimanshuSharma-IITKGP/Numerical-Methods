

syms x;

num_nodes_per_element = 2;
num_elements = 4;
num_dof_per_node = 1;
num_total_nodes = 5;

L = (1/num_elements)*ones(num_elements) ;
alpha = 1/2 ; % crank-Nicolson scheme
dt = 0.001;
f_ext = [1 0; num_total_nodes 0] ; % forcing term (boundary condition)
u = 0*ones(num_total_nodes*num_dof_per_node, 1) ; % initial condition 
u_new = 0*ones(num_total_nodes*num_dof_per_node, 1) ;






for e = 1:num_elements
    
    start_node_of_this_element = e     ;
    end_node_of_this_element   = e + 1 ;

    l = L(e);
    psi_array = [1-x/l, x/l] ; % hardcoded

    M1 = calculateM1(num_dof_per_node, num_nodes_per_element, psi_array, x, l);
    K  = calculateK(num_dof_per_node, num_nodes_per_element, psi_array, x, l);

    a1 = alpha*dt ;
    a2 = (1-alpha)*dt ;

    K_cap_s  = M1 - a2*K ;
    K_cap_s1 = M1 + a1*K ;


    F_sORs1 = zeros(num_nodes_per_element*num_dof_per_node, 1);


    for i = 1:size(f_ext, 1)
        if f_ext(i, 1) == start_node_of_this_element
            F_sORs1(1) = f_ext(i, 2) ;
        elseif f_ext(i, 1) == end_node_of_this_element
            F_sORs1(2) = f_ext(i, 2 ) ;
        end
    end


    F_s_s1 = dt*(alpha*F_sORs1 + (1-alpha)*F_sORs1);



end


function [M1] = calculateM1(num_dof_per_node, num_nodes_per_element, psi_array, x, l)
    M1 = zeros(num_dof_per_node*num_nodes_per_element, num_dof_per_node*num_nodes_per_element );

    for i = 1:size(M1, 1)
        for j = 1:size(M1, 1)
            M1(i, j) = int(psi_array(i)*psi_array(j), x, 0, l) ;
        end
    end
end




function [K] = calculateK(num_dof_per_node, num_nodes_per_element, psi_array, x, l)
    K = zeros(num_dof_per_node*num_nodes_per_element, num_dof_per_node*num_nodes_per_element );

    for i = 1:size(K, 1)
        for j = 1:size(K, 1)
            didx = diff(psi_array(i), x);
            djdx = diff(psi_array(j), x);
            K(i, j) = int(didx*djdx, x, 0, l) ;
        end
    end
end