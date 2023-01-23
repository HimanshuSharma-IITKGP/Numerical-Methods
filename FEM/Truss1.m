clear ;
clc; 
close all;

%% input

% Nodal Coordinates
node = [1/2, sqrt(3)/2; 1, 0; 0, 0];

% Members Connectivity
element = [3, 1; 2, 1; 3, 2] ;

% numbers
num_node = size(node, 1);
num_ele = size(element, 1);
num_dof = 2*num_node ;

% Member Properties
E(1:num_ele) = 1;
A(1:num_ele) = 1;

% Boundary conditions
bcs = [2 1 0, 3 0 0];   % Node number , x-direction constraint, y-direction constraint
                        % 0: constraint, 1: no constraint (free to move)


% External Forces
F_ext = [1, 0, -1];     % Node number, x-direction Component, y-direction Component

%% 

% initialization of Displacement vector and stiffness matrix
U = zeros(num_dof, 1);
K_g = zeros(num_dof, num_dof); 

% Applied loads at DOFs
F = zeros(num_dof, 1);

for i = 1:size(F_ext, 1)
    k = F_ext(i);
    F(2*k-1: 2*k) = F_ext(i, 2: 3) ;
end

% indices of known DOFs due to the boundary conditions
Id_U_k = [] ;
for i=1:size(bcs, 1)
    i1 = 2*bcs(i, 1) - 1;
    i2 = 2*bcs(i, 1);

    if(bcs(i, 2) == 0)
        Id_U_k = [Id_U_k, i1];
    end
    if(bcs(i, 3) == 0)
        Id_U_k = [Id_U_k, i2];
    end
end

% indices of unknown DOFs 
Id_U_u = [];

for i = 1:num_dof
    flat = find(Id_U_k == i);
    if(flag>0)
    else 
        Id_U_u = [Id_U_u, i] ;
    end
end

% initalization of the partitioned stiffness matrix
K_ku = zeros(length(Id_U_u), length(Id_U_u));
K_kk = zeros(length(Id_U_k), length(Id_U_u));


%%

% Computation of the element level stiffness matrix and populating the global stifness matrix

for e = 1:num_ele
    startNode = element(e, 1); % i1
    endNode = element(e, 2); % j1

    % length of the member
    L(e) = sqrt(( node(startNode, 1) - node(endNode, 1) )^2 + ( node(startNode, 2) - node(endNode, 2 )^2));

    % lamda_x and lamda_y
    C = (node(endNode, 1) - node(startNode, 1))/L(e); % lamda_x
    S = (node(endNode, 2) - node(startNode, 2))/L(e); % lamda_y 


    % member stiffness matrix 
    K_e = (A(e)*E(e)/L(e))*[C*C C*S -C*C -C*S; C*S S*S -C*S -S*S; -C*C -C*S C*C C*S, -C*S -S*S C*S S*S];

    % element displacemnt vector
    e_dof_vec = [2*startNode-1, 2*startNode, 2*endNode-1, 2*endNode] ;


    % populating the global stifness matrix

    for i = 1:4
        for j =1:4
            K_g(e_dof_vec(1, i), e_dof_vec(1, j)) = K_g(e_dof_vec(1, i), e_dof_vec(1, j)) + k_e(i, j) ;
        end
    end


end


