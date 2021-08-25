TestName='Test3';
cond_vector=zeros(5,1);

for nRef=2:6

Data = dati(TestName);
[region] = generate_mesh(Data,nRef);
[femregion] = create_dof(Data,region);
[neighbour] = neighbours(femregion);

if (Data.fem(1) == 'D')
   [Matrices] = matrix2D_dubiner(femregion,neighbour,Data,0);
else
   [Matrices] = matrix2D(femregion,neighbour,Data,0);
end


A=Matrices.A;
f0=Matrices.f;
M=Matrices.M;

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

u = eval(Data.initialcond);


if (Data.fem(1)=='D')
   u = fem_to_dubiner (u, femregion,Data);
   w = fem_to_dubiner(w,femregion,Data);
end

t=dt;
    
    f1=assemble_rhs(femregion,neighbour,Data,t);
    [C] = assemble_nonlinear(femregion,Data,u);
   
    w=1/(1+epsilon*gamma*dt)*(w+epsilon*dt*u);
   
    cond_vector(nRef-1)= condest( ChiM*Cm*M+dt*theta*(A+C));
end