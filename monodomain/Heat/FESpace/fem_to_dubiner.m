function [u0] = fem_to_dubiner (uh, femregion, Data)

%FUNCTION TO CONVERT THE COMPONENTS OF A VECTOR W.R.T. FEM BASIS TO COMPONONENTS W.R.T. DUBINER BASIS

% shape functions
[shape_basis] = basis_legendre_dubiner(femregion.fem);

% quadrature nodes and weights for integrals
[nodes_1D, w_1D, nodes_2D, w_2D] = quadrature(Data.nqn);

% evaluation of shape functions on quadrature poiint
[dphiq, Grad, B_edge, G_edge] = evalshape_tria_dubiner(shape_basis,nodes_1D, nodes_2D,Data.nqn,femregion.nln);


u0 = zeros(femregion.ne*femregion.nln,1);

for ie = 1:femregion.ne
    index = (ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    for k = 1:length(w_2D)
        for i = 1 : femregion.nln
            u0(index(i)) = u0(index(i)) + uh(index(i))*dphiq(1,k,i).*w_2D(k);
        end
    end    
end