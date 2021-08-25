TestName='Test2';
cond_vector=zeros(4,3);

for nRef=2:5

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

x=femregion.dof(:,1);
y=femregion.dof(:,2);


u0 = eval(Data.initialcond);


if (Data.fem(1)=='D')
   u0 = fem_to_dubiner (u0, femregion,Data);
end

t=dt;
    
    f1=assemble_rhs(femregion,neighbour,Data,t);

    cond_vector(nRef-1,1)= cond( M,2);
    cond_vector(nRef-1,2)= cond(A,2);
    cond_vector(nRef-1,3)= cond( M+dt*theta*(A),2);
end