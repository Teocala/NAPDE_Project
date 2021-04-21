%=======================================================================================================
% This contain all the information for running main
% TEMPLATE OF THE STRUCT DATI
%=======================================================================================================
%%
%  DATI= struct( 'name',              % set the name of the test  
%                'method',            % (string) e.g. 'SIP','NIP'or 'IIP'
%                'Domain',            % set the domain [x1,x2;y1,y2]
%                'T',                 % set the final time
%                't',                 % set the time step
%                'theta',             % set the Theta-method 
%                'initialcond_i',     % set the initial condition intracellular
%                'initialcond_e',     % set the initial condition extracellular
%                'exact_sol_i',       % set the exact intracellular solution
%                'exact_sol_e',       % set the exact extracellular solution
%                'exact_sol_Vm',      % set the exact transmembrane solution
%                'source_e',          % set the extracellular forcing term
%                'source_i',          % set the intracellular forcing term
%                'Neumann_i",         % set the intracellular Boundary condition
%                'Neumann_e",         % set the extracellular Boundary condition
%                'grad_exact_Vm_x',   % set the first componenet of the gradient of the exact solution
%                'grad_exact_Vm_y',   % set the first componenet of the gradient of the exact solution
%                'ChiM',              % set the parameter monodomain equation
%                'Sigma_i',           % set the diffusion scalar intracellular parameter
%                'Sigma_e',           % set the diffusion scalar extracellular parameter
%                'Cm',                % set the membrane capacity in monodomain/bidomain equation
%                'kappa',             % set the factor for the nonlinear reaction in Fitzhug Nagumo model
%                'epsilon',           % set the parameter ODE
%                'gamma',             % set the parameter ODE
%                'a',                 % set the parameter ODE
%                'initialw',          % set the initial condition ODE
%                'exact_w',           % set the exact solution ODE'(epsilon/(epsilon*gamma-5))*sin(2*pi*x).*sin(2*pi*y).*exp((-5).*t)',...
%                'grad_w_x',          % set the first componenet of the gradient of the ODE
%                'grad_w_y',          % set the second componenet of the gradient of the ODE
%                'fem',               % set finite element space (e.g.,'P1', 'P2', 'P3')
%                'penalty_coeff'      % (real) penalty pearameter
%                'nqn',               % (integer) number of 1D Gauss-Ledendre quadrature nodes in 1 
%                       dimension [1,2,3,4,5,6,.....]
%                'snapshot',          % snapshot of the solution
%                'leap',                % set the number of time steps between one snapshot and the successive              
%========================================================================================================

function [DATA] = dati(test)

if test=='Test1' % test con parametri FitzHugh-Nagumo fisiologici 
    DATA = struct( 'name',          test,...
               ... % Test name
               'method',           'SIP',...  
               ... % Set DG discretization
               'domain',           [0,1;0,1],...
               ... % Reaction term
                'T',                0.00001, ...
               ... % Final time 
               'dt',                0.000001, ...
               ... % Time step 
               'theta',             1, ...
               ... % Theta-method ...
               'initialcond_i',      '2*sin(2*pi*x).*sin(2*pi*y)', ...
               ... % Initial condition intracellular 
               'initialcond_e',      'sin(2*pi*x).*sin(2*pi*y)', ...
               ... % Initial condition extracellular
               'exact_sol_i',        '2*sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)',...
               ... % Definition of exact solution intracellular
               'exact_sol_e',        'sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)',...
               ... % Definition of exact solution extracellular
               'exact_sol_Vm',        'sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)',...
               ... % Definition of exact solution 
               'source_i',           '(-5*ChiM*Cm + 16*pi*pi*sigma_i - kappa*ChiM*(sin(2*pi*x).*sin(2*pi*y).*exp(-5*t) - 1)*(sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)- a) - ChiM*(epsilon/(epsilon*gamma-5))) * sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)',...
               ... % Forcing term in time intracelluar            
               'source_e',           '((+5*ChiM*Cm + 8*pi*pi*sigma_e + kappa*ChiM*(sin(2*pi*x).*sin(2*pi*y).*exp(-5*t) - 1)*(sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)- a) + ChiM*(epsilon/(epsilon*gamma-5))) * sin(2*pi*x).*sin(2*pi*y).*exp(-5*t))',...
               ... % Forcing term in time extracellular
               'Neumann_i',           '2*(sigma_i*(-2*pi*sin(2*pi*x).*cos(2*pi*y).*exp(-5*t).*(y==0) + 2*pi*cos(2*pi*x).*sin(2*pi*y).*exp(-5*t).*(x==1) + 2*pi*sin(2*pi*x).*cos(2*pi*y).*exp(-5*t).*(y==1) -  2*pi*cos(2*pi*x).*sin(2*pi*y).*exp(-5*t).*(x==0)))',...
               ... % Boundary condition intracellular
               'Neumann_e',           'sigma_e*(-2*pi*sin(2*pi*x).*cos(2*pi*y).*exp(-5*t).*(y==0) + 2*pi*cos(2*pi*x).*sin(2*pi*y).*exp(-5*t).*(x==1) + 2*pi*sin(2*pi*x).*cos(2*pi*y).*exp(-5*t).*(y==1) -  2*pi*cos(2*pi*x).*sin(2*pi*y).*exp(-5*t).*(x==0))',...
               ... % Boundary condition extracellular (taken with "-")
               'grad_exact_Vm_x',     '2*pi*cos(2*pi*x).*sin(2*pi*y).*exp(-5*t)',... 
               ... % Definition of exact gradient Vm wrt x 
               'grad_exact_Vm_y',     '2*pi*sin(2*pi*x).*cos(2*pi*y).*exp(-5*t)',...    
               ... % Definition of exact gradient Vm wrt y 
               'ChiM',              1,...
               ... % Parameter monodomain equation
               'Sigma_i',           1,...
               ... % Diffusion scalar parameter intracellular
               'Sigma_e',           1,...
               ... % Diffusion scalar parameter extracellular
               'Cm',                1,...
               ... % Membrane capacity in monodomain equation
               'kappa',             1.5*13,...
               ... % Factor for the nonlinear reaction in Fitzhug Nagumo model
               'epsilon',          0.012*100,...
               ... % Parameter ODE
               'gamma',            0.1,...
               ... % Parameter ODE
               'a',                13e-3,...
               ... % Parameter ODE 
               'initialw',         '(epsilon/(epsilon*gamma-5))*sin(2*pi*x).*sin(2*pi*y)',...
               ... % Initial condition ODE
               'exact_w',          '(epsilon/(epsilon*gamma-5))*sin(2*pi*x).*sin(2*pi*y).*exp((-5).*t)',...
               ... % Grad_exact_w_x 
               'grad_w_x',          '2*pi*(epsilon/(epsilon*gamma-5))*cos(2*pi*x).*sin(2*pi*y).*exp((-5).*t)',...
                ... % Grad_exact_w 
               'grad_w_y',          '2*pi*(epsilon/(epsilon*gamma-5))*sin(2*pi*x).*cos(2*pi*y).*exp((-5).*t)',...  
                ... % Exact solution of ODE 
               'fem',               'P2',...   
               ... % Finite element space (other choices 'P2', 'P3')'
               'penalty_coeff',     10,... 
               ... % Penalty coefficient
               'nqn',               4, ...
               ... % Number of 1d GL quadrature nodes
               'snapshot',          'N',...
               ... % Snapshot of the solution
               'leap',               40 ...
               ... % Number of time steps between one snapshot and the successive
               );
           
elseif test=='Test2' % test con parametri FitzHugh-Nagumo unitari 
DATA = struct( 'name',             test,...
               ... % Test name
               'method',           'SIP',...  
               ... % Set DG discretization
               'domain',           [0,1;0,1],...
               ... % Reaction term
                'T',                0.001, ...
               ... % Final time 
               'dt',                0.0001, ...
               ... % Time step 
               'theta',             1, ...
               ... % Theta-method ...
               'initialcond_i',      '2*sin(2*pi*x).*sin(2*pi*y)', ...
               ... % Initial condition intracellular 
               'initialcond_e',      'sin(2*pi*x).*sin(2*pi*y)', ...
               ... % Initial condition extracellular
               'exact_sol_i',        '2*sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)',...
               ... % Definition of exact solution intracellular
               'exact_sol_e',        'sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)',...
               ... % Definition of exact solution extracellular
               'exact_sol_Vm',        'sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)',...
               ... % Definition of exact solution 
               'source_i',           '(-5*ChiM*Cm + 16*pi*pi*sigma_i - kappa*ChiM*(sin(2*pi*x).*sin(2*pi*y).*exp(-5*t) - 1)*(sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)- a) - ChiM*(epsilon/(epsilon*gamma-5))) * sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)',...
               ... % Forcing term in time intracelluar            
               'source_e',           '((+5*ChiM*Cm + 8*pi*pi*sigma_e + kappa*ChiM*(sin(2*pi*x).*sin(2*pi*y).*exp(-5*t) - 1)*(sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)- a) + ChiM*(epsilon/(epsilon*gamma-5))) * sin(2*pi*x).*sin(2*pi*y).*exp(-5*t))',...
               ... % Forcing term in time extracellular
               'Neumann_i',           '2*(sigma_i*(-2*pi*sin(2*pi*x).*cos(2*pi*y).*exp(-5*t).*(y==0) + 2*pi*cos(2*pi*x).*sin(2*pi*y).*exp(-5*t).*(x==1) + 2*pi*sin(2*pi*x).*cos(2*pi*y).*exp(-5*t).*(y==1) -  2*pi*cos(2*pi*x).*sin(2*pi*y).*exp(-5*t).*(x==0)))',...
               ... % Boundary condition intracellular
               'Neumann_e',           'sigma_e*(-2*pi*sin(2*pi*x).*cos(2*pi*y).*exp(-5*t).*(y==0) + 2*pi*cos(2*pi*x).*sin(2*pi*y).*exp(-5*t).*(x==1) + 2*pi*sin(2*pi*x).*cos(2*pi*y).*exp(-5*t).*(y==1) -  2*pi*cos(2*pi*x).*sin(2*pi*y).*exp(-5*t).*(x==0))',...
               ... % Boundary condition extracellular (taken with "-")
               'grad_exact_Vm_x',     '2*pi*cos(2*pi*x).*sin(2*pi*y).*exp(-5*t)',... 
               ... % Definition of exact gradient Vm wrt x 
               'grad_exact_Vm_y',     '2*pi*sin(2*pi*x).*cos(2*pi*y).*exp(-5*t)',...    
               ... % Definition of exact gradient Vm wrt y 
               'ChiM',              1,...
               ... % Parameter monodomain equation
               'Sigma_i',           1,...
               ... % Diffusion scalar parameter intracellular
               'Sigma_e',           1,...
               ... % Diffusion scalar parameter extracellular
               'Cm',                1,...
               ... % Membrane capacity in monodomain equation
               'kappa',             1,...
               ... % Factor for the nonlinear reaction in Fitzhug Nagumo model
               'epsilon',          1,...
               ... % Parameter ODE
               'gamma',            1,...
               ... % Parameter ODE
               'a',                1,...
               ... % Parameter ODE 
               'initialw',         '(epsilon/(epsilon*gamma-5))*sin(2*pi*x).*sin(2*pi*y)',...
               ... % Initial condition ODE
               'exact_w',          '(epsilon/(epsilon*gamma-5))*sin(2*pi*x).*sin(2*pi*y).*exp((-5).*t)',...
               ... % Grad_exact_w_x 
               'grad_w_x',          '2*pi*(epsilon/(epsilon*gamma-5))*cos(2*pi*x).*sin(2*pi*y).*exp((-5).*t)',...
                ... % Grad_exact_w 
               'grad_w_y',          '2*pi*(epsilon/(epsilon*gamma-5))*sin(2*pi*x).*cos(2*pi*y).*exp((-5).*t)',...  
                ... % Exact solution of ODE 
               'fem',               'P3',...   
               ... % Finite element space (other choices 'P2', 'P3')'
               'penalty_coeff',     10,... 
               ... % Penalty coefficient
               'nqn',               4, ...
               ... % Number of 1d GL quadrature nodes
               'snapshot',          'N',...
               ... % Snapshot of the solution
               'leap',               40 ...
               ... % Number of time steps between one snapshot and the successive
               );
end;