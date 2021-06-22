function [phi_i,phi_e]=assign(phi_i,phi_e,Data,femregion)

% shape functions
[shape_basis] = basis_lagrange(femregion.fem);

% quadrature nodes and weights for integrals
[nodes_1D, w_1D, nodes_2D, w_2D] = quadrature(Data.nqn);

% evaluation of shape functions on quadrature poiint
[dphiq, Grad, B_edge, G_edge] = evalshape(shape_basis,nodes_2D,nodes_1D,femregion.nln);


for ie = 1:femregion.ne
    
    % Local to global map --> To be used in the assembly phase
    index = (ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    
    % Index of the current edges
    index_element = femregion.nedges*(ie-1).*ones(femregion.nedges,1) + [1:1:femregion.nedges]';
    
    % Coordinates of the verteces of the current triangle
    coords_elem = femregion.coords_element(index_element, :);
    
    [BJ, BJinv, pphys_2D] = get_jacobian_physical_points(coords_elem, nodes_2D);
    
    c=0;
    d=0;
    % =====================================================================
    % Compute integrals over triangles
    % =====================================================================
    for k = 1:length(w_2D) % loop over 2D quadrature nodes
        % scaled weight for the quadrature formula
        dx = w_2D(k)*det(BJ);
 
        % evaluation of the load term       
        x = pphys_2D(k,1);
        y = pphys_2D(k,2);
       
        for i = 1 : femregion.nln
            c=c+phi_i(index(i))*dphiq(1,k,i).*dx;
            d=d+phi_e(index(i))*dphiq(1,k,i).*dx;
        end
    end
    phi_i(index(i))=phi_i(index(i))-c;
    phi_e(index(i))=phi_e(index(i))-d;
end