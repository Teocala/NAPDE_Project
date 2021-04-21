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

[Matrices] = matrix2D(femregion,neighbour,Data,0);

%==========================================================================
% SOLVE THE LINEAR SYSTEM
%==========================================================================

A=Matrices.A;
f0=Matrices.f;
M=Matrices.M;
% M=masslumping(M);

%time integration parameters
t=0;
T=Data.T;
dt=Data.dt;
theta=Data.theta;
epsilon = Data.epsilon;
gamma = Data.gamma;
ChiM= Data.ChiM;
Cm= Data.Cm;

x=femregion.dof(:,1);
y=femregion.dof(:,2);


w = eval(Data.initialw);

%figure(1)
u = eval(Data.initialcond);

for t=dt:dt:T
    
    f1=assemble_rhs(femregion,neighbour,Data,t);
    [C] = assemble_nonlinear(femregion,Data,u);
   
    w=1/(1+epsilon*gamma*dt)*(w+epsilon*dt*u);
   
    r=dt*(theta*f1+(1-theta)*f0)+(ChiM*Cm*M-dt*(1-theta)*(A+C))*u+dt*ChiM*M*w;
    
    u=(ChiM*Cm*M+dt*theta*(A+C))\r;
    if (Data.snapshot=='Y' && ( (mod(round(t/dt),Data.leap)==0) || (t/dt)<=5))
%         ODE_Snapshot(femregion,Data,w1,t)
        DG_Par_Snapshot(femregion, Data, u,t);
    end
    f0=f1;
end


%==========================================================================
% POST-PROCESSING OF THE SOLUTION
%=========================================================================

[solutions]= postprocessing(femregion,Data,u,T);

%==========================================================================
% ERROR ANALYSIS
%==========================================================================
 [errors]= compute_errors(Data,femregion,solutions,Matrices.S,T);

%  solutionsW=struct('u_h',w0,'u_ex',eval(Data.exact_w));
%  [errorsW]= compute_errors(Data,femregion,solutionsW,Matrices.S,T);
end