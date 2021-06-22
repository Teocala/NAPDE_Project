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

if test=='Test1' % Test for wave front simulation (no exact solution) 
    DATA = struct( 'name',          test,...
               ... % Test name
               'method',              'SI',...
               ... % Set time discretization (SI, OS, GO)
               'DG_method',           'SIP',...  
               ... % Set DG discretization
               'domain',           [0,0.1;0,0.1],...
               ... % Reaction term
                'T',                0.005, ...
               ... % Final time 
               'dt',                0.0001, ...
               ... % Time step 
               'theta',             1, ...
               ... % Theta-method ...
               'initialcond',      '0.*x.*y', ...
               ... % Initial condition for trans-membrane potentialz
               'exact_sol_i',        '0.*x.*y',...
               ... % Definition of exact solution intracellular
               'exact_sol_e',        '0.*x.*y',...
               ... % Definition of exact solution extracellular
               'exact_sol_Vm',       '0.*x.*y',...
               ... % Definition of exact solution 
               'source_i',           '100*10^3*(t>=0.001).*(t<0.002).*(x>0.04).*(x<0.06).*(y>0.04).*(y<0.06)',...
               ... % Forcing term in time intracelluar            
               'source_e',           '-100*10^3*(t>=0.001).*(t<0.002).*(x>0.04).*(x<0.06).*(y>0.04).*(y<0.06)',...    
               ... % Forcing term in time extracellular
               'Neumann_i',           '0.*x.*y',...
               ... % Boundary condition intracellular
               'Neumann_e',           '0.*x.*y',...
               ... % Boundary condition extracellular (taken with "-")
               'grad_exact_Vm_x',     '0.*x.*y',... 
               ... % Definition of exact gradient Vm wrt x 
               'grad_exact_Vm_y',     '0.*x.*y',...    
               ... % Definition of exact gradient Vm wrt y 
               'ChiM',              1e5,...
               ... % Parameter monodomain equation
               'Sigma_i',           [0.17 0; 0 0.019],...
               ... % Diffusion scalar parameter intracellular
               'Sigma_e',           [0.62 0; 0 0.24],...
               ... % Diffusion scalar parameter extracellular
               'Cm',                1e-2,...
               ... % Membrane capacity in monodomain equation
               'kappa',             1.5*13,...
               ... % Factor for the nonlinear reaction in Fitzhug Nagumo model
               'epsilon',          0.012*100,...
               ... % Parameter ODE
               'gamma',            0.1,...
               ... % Parameter ODE
               'a',                13e-3,...
               ... % Parameter ODE 
               'initialw',         '0.*x.*y',...
               ... % Initial condition ODE
               'exact_w',          '0.*x.*y',...
               ... % Exact solution of ODE 
               'grad_w_x',          '0.*x.*y',...
                ... % Grad_exact_w_x
               'grad_w_y',          '0.*x.*y',...  
                ... % Grad_exact_w_y  
               'fem',               'D1',...   
               ... % Finite element space (choices 'P1,'D1','P2','D2', 'P3','D3')
               'penalty_coeff',     10,... 
               ... % Penalty coefficient
               'nqn',               4, ...
               ... % Number of 1d GL quadrature nodes
               'snapshot',          'Y',...
               ... % Snapshot of the solution
               'leap',               10, ...
               ... % Number of time steps between one snapshot and the successive
               'assign',             2 ...
               ... % Type of phi_i imposition (0 = no imposition, 1 = value in a point, 2 = zero average)
               );
           

elseif test=='Test2' % Test for error analysis
DATA = struct( 'name',             test,...
               ... % Test name
               'method',              'SI',...
               ... % Set time discretization (SI, OS, GO)
               'DG_method',           'SIP',...  
               ... % Set DG discretization
               'domain',           [0,1;0,1],...
               ... % Reaction term
                'T',                0.001, ...
               ... % Final time 
               'dt',                0.0001, ...
               ... % Time step 
               'theta',             1, ...
               ... % Theta-method ...
               'initialcond',      'sin(2*pi*x).*sin(2*pi*y)', ...
               ... % Initial condition for trans-membrane potential
               'exact_sol_i',        '2*sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)',...
               ... % Definition of exact solution intracellular
               'exact_sol_e',        'sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)',...
               ... % Definition of exact solution extracellular
               'exact_sol_Vm',        'sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)',...
               ... % Definition of exact solution 
               'source_i',           '(-5*ChiM*Cm + 16*pi*pi*0.12 + kappa*ChiM*(sin(2*pi*x).*sin(2*pi*y).*exp(-5*t) - 1)*(sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)- a) + ChiM*(epsilon/(epsilon*gamma-5))) * sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)',...
               ... % Forcing term in time intracelluar            
               'source_e',           '(-5*ChiM*Cm - 8*pi*pi*0.12 + kappa*ChiM*(sin(2*pi*x).*sin(2*pi*y).*exp(-5*t) - 1)*(sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)- a) + ChiM*(epsilon/(epsilon*gamma-5))) * sin(2*pi*x).*sin(2*pi*y).*exp(-5*t)',...
               ... % Forcing term in time extracellular
               'Neumann_i',           '2*(0.12*(-2*pi*sin(2*pi*x).*cos(2*pi*y).*exp(-5*t).*(y==0) + 2*pi*cos(2*pi*x).*sin(2*pi*y).*exp(-5*t).*(x==1) + 2*pi*sin(2*pi*x).*cos(2*pi*y).*exp(-5*t).*(y==1) -  2*pi*cos(2*pi*x).*sin(2*pi*y).*exp(-5*t).*(x==0)))',...
               ... % Boundary condition intracellular
               'Neumann_e',           '0.12*(-2*pi*sin(2*pi*x).*cos(2*pi*y).*exp(-5*t).*(y==0) + 2*pi*cos(2*pi*x).*sin(2*pi*y).*exp(-5*t).*(x==1) + 2*pi*sin(2*pi*x).*cos(2*pi*y).*exp(-5*t).*(y==1) -  2*pi*cos(2*pi*x).*sin(2*pi*y).*exp(-5*t).*(x==0))',...
               ... % Boundary condition extracellular (taken with "-")
               'grad_exact_Vm_x',     '2*pi*cos(2*pi*x).*sin(2*pi*y).*exp(-5*t)',... 
               ... % Definition of exact gradient Vm wrt x 
               'grad_exact_Vm_y',     '2*pi*sin(2*pi*x).*cos(2*pi*y).*exp(-5*t)',...    
               ... % Definition of exact gradient Vm wrt y 
               'ChiM',              1e5,...
               ... % Parameter monodomain equation
               'Sigma_i',           0.12.*[1 0; 0 1],...
               ... % Diffusion scalar parameter intracellular
               'Sigma_e',           0.12.*[1 0; 0 1],...
               ... % Diffusion scalar parameter extracellular
               'Cm',                1e-2,...
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
               ... % Exact solution of ODE
               'grad_w_x',          '2*pi*(epsilon/(epsilon*gamma-5))*cos(2*pi*x).*sin(2*pi*y).*exp((-5).*t)',...
                ... % Grad_exact_w_x 
               'grad_w_y',          '2*pi*(epsilon/(epsilon*gamma-5))*sin(2*pi*x).*cos(2*pi*y).*exp((-5).*t)',...  
                ... % Grad_exact_w_y
               'fem',               'D1',...   
               ... % Finite element space (choices 'P1,'D1','P2','D2', 'P3','D3')
               'penalty_coeff',     10,... 
               ... % Penalty coefficient
               'nqn',               4, ...
               ... % Number of 1d GL quadrature nodes
               'snapshot',          'N',...
               ... % Snapshot of the solution
               'leap',               40, ...
               ... % Number of time steps between one snapshot and the successive
               'assign',             0 ...
               ... % Type of phi_i imposition (0 = no imposition, 1 = value in a point, 2 = zero average)
               );

end