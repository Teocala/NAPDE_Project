function convergence_testW(TestName,nRef)
%              usage:   convergence_test('Test1',[2 3 4 5]) 

%     num_tests=length(nRef);
%     for i=1:num_tests
%         fprintf('Stage %d on %d \n', i, num_tests);
%         [errors,solutions,femregion,Data]= main2D(TestName,nRef(i));
%         
%         err_L2(i)=errors.E_L2;
%    
%         h(i)=femregion.h;
% 
%     end

    close all
    figure()
    
    
    loglog(h, err_L2, 'o-', 'linewidth', 2)
    hold on
    loglog(h , err_L2(1)*h/h(1), 'k-', 'linewidth', 2)
    loglog(h , err_L2(1)*(h/h(1)).^2, 'k--', 'linewidth', 2)
    loglog(h , err_L2(1)*(h/h(1)).^3, 'k-.', 'linewidth', 2)
    xlabel('h')
    legend('err L^2', 'h', 'h^2', 'h^3')
    
   
   % txt=['dt= ', num2str(Data.dt), ', T= ', num2str(Data.T), ', Penalty ', Data.method, ', Degree ', Data.fem ];
    %sgtitle(txt)
    
    

    
end

