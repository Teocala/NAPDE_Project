function [ C ] = assemble_nonlinear(femregion,Data,u0)
%ASSEMBLE_NONLINEAR Summary of this function goes here
%   Detailed explanation goes here
addpath FESpace
addpath Assembly



% shape functions
[shape_basis] = basis_lagrange(femregion.fem);

% quadrature nodes and weights for integrals
[nodes_1D, w_1D, nodes_2D, w_2D] = quadrature(Data.nqn);
nqn_1D = length(w_1D);

% evaluation of shape functions on quadrature poiint
[dphiq, Grad, B_edge, G_edge] = evalshape(shape_basis,nodes_2D,nodes_1D,femregion.nln);

C=sparse(femregion.ndof,femregion.ndof);%int_omega((u-a)(u-1)u*vdx)


a=Data.a;
kappa= Data.kappa;
ChiM= Data.ChiM;
for ie = 1:femregion.ne
    
    % Local to global map --> To be used in the assembly phase
    index = (ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    
    % Index of the current edges
    index_element = femregion.nedges*(ie-1).*ones(femregion.nedges,1) + [1:1:femregion.nedges]';
    
    
    % Coordinates of the verteces of the current triangle
    coords_elem = femregion.coords_element(index_element, :);
    
    % BJ        = Jacobian of the elemental map
    % BJinv     = Inverse Jacobian of the elemental map
    % pphys_2D = vertex coordinates in the physical domain
    [BJ, BJinv, pphys_2D] = get_jacobian_physical_points(coords_elem, nodes_2D);
    
 
    
    % =====================================================================
    % Compute integrals over triangles
    % =====================================================================
    for k = 1:length(w_2D) % loop over 2D quadrature nodes
        
        % scaled weight for the quadrature formula
        dx = w_2D(k)*det(BJ);
        nonlu=0;
        for s=1:femregion.nln
        nonlu= nonlu+ u0(index(s)).*dphiq(1,k,s);
        end
        
        nonl= (nonlu-1)*(nonlu-a);
       
        
        for i = 1 : femregion.nln    
            
            for j = 1 : femregion.nln
             %Compute non linear matrix
                C(index(i),index(j))= C(index(i),index(j)) +...
                   nonl * dphiq(1,k,i)*(dphiq(1,k,j)).*dx;
            end
        end
    end
end
C=-kappa*ChiM*C;
end

