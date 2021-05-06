function [u0] = dubiner_to_fem (uh, femregion, Data)

%FUNCTION TO CONVERT THE COMPONENTS OF A VECTOR W.R.T. DUBINER BASIS TO COMPONONENTS W.R.T. FEM BASIS

         
deg=sscanf(Data.fem(2:end),'%f');
s=0;

if (deg==1)  % for the time being, implemented only for D1
    csi = [0;1;0];
    eta = [0;0;1];
    a   = [-1; 1; -1];
    b   = [-1; -1; 1];
end

for j=0:(deg)
    for i=0:(deg)
        if (i+j) <= deg
           s=s+1;
           [pi] = eval_jacobi_polynomial(i,0,0,a);
           [pj] = eval_jacobi_polynomial(j,2.*i+1,0,b);
           cij=sqrt((2.*i +1).*2.*(i+j+1)./4.^i);
           phi(1,:,s)=cij.*(2.^i).*((1-eta).^i).*pi.*pj;
        end
    end
end

u0 = zeros(femregion.ne*femregion.nln,1);

for ie = 1:femregion.ne
    
   index = (ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
   
   for i = 1 : femregion.nln
       for j = 1: femregion.nln
         u0(index(i)) = u0(index(i)) +  uh(index(j))*phi(1,i,j);
       end
   end
    
end