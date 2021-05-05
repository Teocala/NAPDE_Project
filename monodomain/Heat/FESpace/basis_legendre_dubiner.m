%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Genera le funzioni di base di legendre Qp--Pp sull'elemento di riferimento
% [-1,1]X[-1,1]
%
% INPUT: tipo di funzioni di base da utilizzare 
% OUTPUT: la struttura contenente il numero delle funzioni di 
%       base e la loro espressione analitica
%       e le componenti del gradiente
%
% funzioni di base --> Espressione analitica.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [shape_basis]=basis_legendre_dubiner(fem)


fem=char(fem);

type_elem=fem(1);
degree=sscanf(fem(2:end),'%f');
switch type_elem
    case{'P','D'}, 
        nln=0.5.*(degree+1).*(degree+2);
        n_edge=3;
    case{'Q'}, 
        nln=(degree+1).^2;
        n_edge=4;
end

shape_basis=struct('nln',nln,...
             'n_edge',n_edge,...
             'deg',degree);   



