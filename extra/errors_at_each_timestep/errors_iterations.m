function [error] = errors_iterations(femregion,Data,u_h,t)

[neighbour] = neighbours(femregion);


if (Data.fem(1)=='D')   
   u_h = dubiner_to_fem (u_h,femregion,Data);
end

x=femregion.dof(:,1);
y=femregion.dof(:,2);
u_ex=eval(Data.exact_sol_Vm);
solutions=struct('u_h',u_h,'u_ex',u_ex);
S = matrix_S(append("P", femregion.fem(2)),femregion,neighbour,Data,0);
error = compute_errors_Vm(Data,femregion,solutions,S,t); 