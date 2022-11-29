% solution of axial bar problem -d/dx(AE*du/dx) - q = 0 with boundary condition
% u(0) = 0 and Q(L) = P


syms x;



num_elements = 3;
num_nodes_per_element = 2;
num_nodes = num_elements*num_nodes_per_element + 1;





function [K_G] = assemble_K(num_element)


function [Ue] = solveE(Be, Fe, Qe )
    Ue = pinv(Be)*(Fe-Qe);
end



function [phiE_array] = createPhiE(num_nodes_per_element, x, Le)

    if(num_nodes_per_element == 2)
        phiE_array(1) = (1-(x/Le)) ;
        phiE_array(2) = x/Le       ;
    end
    
end


function [Fe] = createFe(phiE_array, num_nodes_per_element, Le, q)
    
    Fe = zeros(num_nodes_per_element, 1);

    for i = 1:num_nodes_per_element
        Fe(i) = linearFunctional(phiE_array, i, Le, q); 
    end
end

function [Ke] = createKe(phiE_array, num_nodes_per_element, x, Le, A, E)
    Ke = zeros(num_nodes_per_element, 1);

    for i = 1:num_nodes_per_element
        for j =1:num_nodes_per_element
            Ke(i, j) = bilinearFunctional(phiE_array, i, j, x, Le, A, E);
        end
    end
end


function [Fe_i] = linearFunctional(phiE_array, i, Le, q)
    Fe_i = int(q*phiE_array(i), 0, Le);
end

function [Ke_ij] = bilinearFunctional(phiE_array, i, j, x, Le, A, E)
   
    phiE_i = phiE_array(i);
    phiE_j = phiE_array(j);
    
    ddx_i = diff(phiE_i, x);
    ddx_j = diff(phiE_j, x);

    Ke_ij = int(A*E*ddx_j*ddx_i, 0, Le) ;
end





