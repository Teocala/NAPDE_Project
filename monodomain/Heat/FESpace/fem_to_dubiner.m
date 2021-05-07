function [u0] = fem_to_dubiner (uh, femregion, Data)

%FUNCTION TO CONVERT THE COMPONENTS OF A VECTOR W.R.T. FEM BASIS TO COMPONONENTS W.R.T. DUBINER BASIS

% shape functions


% quadrature nodes and weights for integrals
[nodes_1D, w_1D, nodes_2D, w_2D] = quadrature(Data.nqn);

% evaluation of shape functions on quadrature poiint
[shape_basis] = basis_legendre_dubiner(femregion.fem);
[phi_dub, Grad, B_edge, G_edge] = evalshape_tria_dubiner(shape_basis,nodes_1D, nodes_2D,Data.nqn,femregion.nln);
[shape_basis] = basis_lagrange(append("P", femregion.fem(2)));
[phi_fem, Grad, B_edge, G_edge] = evalshape(shape_basis,nodes_2D,nodes_1D,femregion.nln);

u0 = zeros(femregion.ne*femregion.nln,1);



for ie = 1:femregion.ne
    index = (ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    for i = 1 : femregion.nln
        for k = 1:length(w_2D) 
            uh_eval_k = 0;
            for j = 1:femregion.nln
                uh_eval_k = uh_eval_k + uh(index(j))*phi_fem(1,k,j);
            end
            u0(index(i)) = u0(index(i)) + uh_eval_k*phi_dub(1,k,i).*w_2D(k);
        end
    end    
end