function [errors,solutions,femregion,Data]= main2D(TestName,nRef)
%==========================================================================
% Solution of the Poisson's problem with DG finite elements
% (non homogeneous Dirichlet boundary conditions)
%==========================================================================
%
%    INPUT:
%          Data        : (struct)  see dati.m
%          nRef        : (int)     refinement level
%
%    OUTPUT:
%          errors      : (struct) contains the computed errors
%          solutions   : (sparse) nodal values of the computed and exact
%                        solution
%          femregion   : (struct) infos about finite elements
%                        discretization
%          Data        : (struct)  see dati.m
%          
% Usage: 
%    [errors,solutions,femregion,Dati] = main2D('Test1',3)



addpath Assembly
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath PostProcessing

close all

%==========================================================================
% LOAD DATA FOR TEST CASE
%==========================================================================

Data = dati(TestName);

%==========================================================================
% MESH GENERATION
%==========================================================================

[region] = generate_mesh(Data,nRef);

%==========================================================================
% FINITE ELEMENT REGION
%==========================================================================

[femregion] = create_dof(Data,region);


%==========================================================================
% CONNECTIVITY FOR NEIGHBOURING ELEMENTS
%==========================================================================

[neighbour] = neighbours(femregion);


%==========================================================================
% BUILD FINITE ELEMENT MATRICES and RIGHT-HAND SIDE
%==========================================================================
if (Data.fem(1) == 'D')
   [Matrices] = matrix2D_dubiner(femregion,neighbour,Data,0);
else
   [Matrices] = matrix2D(femregion,neighbour,Data,0);
end
%==========================================================================
% SOLVE THE LINEAR SYSTEM
%==========================================================================

A=Matrices.A;
f0=Matrices.f;
M=Matrices.M;

%time integration parameters
t=0;
T=Data.T;
dt=Data.dt;
theta=Data.theta;


x=femregion.dof(:,1);
y=femregion.dof(:,2);


u0 = eval(Data.initialcond);

% in the case of Dubiner basis transforms the initial condition wrt Dubiner
% basis
if (Data.fem(1)=='D')
   u0 = fem_to_dubiner (u0, femregion,Data);
end


for t=dt:dt:T
    
    f1=assemble_rhs(femregion,neighbour,Data,t);

    r=dt*(theta*f1+(1-theta)*f0)+(M-dt*(1-theta)*(A))*u0;
    
    u1=(M+dt*theta*(A))\r;
    if (Data.snapshot=='Y' && ( (mod(round(t/dt),Data.leap)==0) || (t/dt)<=5))
        DG_Par_Snapshot(femregion, Data, u0,t);
    end
    f0=f1;
    u0=u1;
end


% pass to fem coordinates before post-processing
if (Data.fem(1)=='D')
   u0 = dubiner_to_fem (u0, femregion,Data);
end

%==========================================================================
% POST-PROCESSING OF THE SOLUTION
%=========================================================================

[solutions]= postprocessing(femregion,Data,u0,T);

%==========================================================================
% ERROR ANALYSIS
%==========================================================================

if (Data.fem(1) == 'D')
    S = matrix_S(append("P", femregion.fem(2)),femregion,neighbour,Data,0);
    [errors]= compute_errors(Data,femregion,solutions,S,T);
else
    [errors]= compute_errors(Data,femregion,solutions,Matrices.S,T);
end







