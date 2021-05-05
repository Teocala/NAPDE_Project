%--------------------------------------------------------------------
% PURPOSE:
%
function [errors]= compute_errors(Dati,femregion,solutions,S,t)
%% [errors]= compute_errors(Dati,femregion,solutions,S)
%==========================================================================
% Compute L2, H1, DG and L-inf errors
%==========================================================================
%    called in main2D.m
%
%    INPUT:
%          Dati        : (struct)  see dati.m
%          femregion   : (struct)  see create_dof.m
%          solutions   : (struct)  see postprocessing.m
%          S           : (array real) matrix associated to the penalty term
%    OUTPUT:
%          errors      : (struct)  


nln=femregion.nln;
ne=femregion.ne;

% initialization
E_L2=0;
E_SEMI_H1=0;
E_DG=0;

% scalar shape functions
if (Dati.fem(1) == 'P')
[shape_basis]= basis_lagrange(Dati.fem);
else
[shape_basis]= basis_legendre_dubiner(Dati.fem);
end

% 1D and 2D quadrature nodes and weights 
[nodes_1D, w_1D, nodes_2D, w_2D]=quadrature(Dati.nqn);

% evaluation of shape functions
if (Dati.fem(1) == 'P')
[dphiq, Grad, B_edge, G_edge] = evalshape(shape_basis,nodes_2D,nodes_1D,femregion.nln);
else
[dphiq, Grad, B_edge, G_edge] = evalshape_tria_dubiner(shape_basis,nodes_1D, nodes_2D,Dati.nqn,femregion.nln);
end

for ie=1:ne % loop over elements

    index=(ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    index_element=femregion.nedges*(ie-1).*ones(femregion.nedges,1) + [1:1:femregion.nedges]';

    coords_elem=femregion.coords_element(index_element, :);

    local_uh=solutions.u_h(index);
    
    [BJ, BJinv, pphys_2D] = get_jacobian_physical_points(coords_elem, nodes_2D); %  Jacobian and physical coordinates at the quadrature points
    
    for k=1:length(w_2D) % loop over quadrature nodes

        dx = det(BJ).*w_2D(k);

        x=pphys_2D(k,1);
        y=pphys_2D(k,2);
        local_exact=eval(Dati.exact_sol)';
        local_grad_exact=[eval(Dati.grad_exact_1),eval(Dati.grad_exact_2)];
        local_grad_aprox=zeros(1,2);
        local_aprox=0;
        for s=1:nln  % reconstruct the discrete solution and its gradient at the quadrature nodes
            local_aprox =local_aprox  + dphiq(1,k,s).*local_uh(s); 
            local_grad_aprox =local_grad_aprox + Grad(k,:,s).*local_uh(s);
        end

        pointwise_diff =(local_grad_exact)-(local_grad_aprox*BJinv);

        E_SEMI_H1= E_SEMI_H1 + (pointwise_diff * transpose(pointwise_diff)).*dx;
        E_L2=E_L2 + ((local_aprox-local_exact).^2).*dx;
    end
end



E_DG = E_SEMI_H1 + (solutions.u_h-solutions.u_ex)'* S * (solutions.u_h-solutions.u_ex); 

E_L2=sqrt(E_L2);
E_SEMI_H1=sqrt(E_SEMI_H1);
E_DG = sqrt(E_DG);


errors=struct('E_L2',E_L2,...
    'E_SEMI_H1',E_SEMI_H1,...
    'E_H1', sqrt((E_L2.^2) + (E_SEMI_H1.^2)),...
    'E_inf',norm(solutions.u_h-solutions.u_ex,inf),...
    'E_DG', E_DG);




