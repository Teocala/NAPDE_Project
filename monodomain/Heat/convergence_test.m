function convergence_test(TestName,nRef)
%              usage:   convergence_test('Test2',[2 3 4 5]) 

    num_tests=length(nRef);
    for i=1:num_tests
        fprintf('Stage %d on %d \n', i, num_tests);
        [errors,solutions,femregion,Data]= main2D(TestName,nRef(i));
        err_L2(i)=errors.E_L2;
        err_H1(i)=errors.E_H1;
        err_DG(i)=errors.E_DG;
        h(i)=femregion.h;
    end
    close all
    figure()
    
    subplot(1,3,1)
    loglog(h, err_L2, 'o-', 'linewidth', 2)
    hold on
   % loglog(h , err_L2(1)*h/h(1), 'k-', 'linewidth', 2)
    loglog(h , err_L2(1)*(h/h(1)).^2, 'k--', 'linewidth', 2)
    loglog(h , err_L2(1)*(h/h(1)).^3, 'k-.', 'linewidth', 2)
    loglog(h , err_L2(1)*(h/h(1)).^4, 'k-.', 'linewidth', 2)
    xlabel('h')
    legend('err L^2', 'h^2', 'h^3' , 'h^4')
    
    subplot(1,3,2)
    loglog(h, err_H1, 'o-', 'linewidth', 2)
    hold on
    loglog(h , err_H1(1)*h/h(1), 'k-', 'linewidth', 2)
    loglog(h , err_H1(1)*(h/h(1)).^2, 'k--', 'linewidth', 2)
    loglog(h , err_H1(1)*(h/h(1)).^3, 'k-.', 'linewidth', 2)
    xlabel('h')
    legend('err H^1', 'h', 'h^2', 'h^3')
    
    subplot(1,3,3)
    loglog(h, err_DG, 'o-','linewidth',2)
    hold on
    loglog(h , err_DG(1)*h/h(1), 'k-', 'linewidth', 2)
    loglog(h , err_DG(1)*(h/h(1)).^2, 'k--', 'linewidth', 2)
    loglog(h , err_DG(1)*(h/h(1)).^3, 'k-.', 'linewidth', 2)
    xlabel('h')
    legend('err DG', 'h', 'h^2', 'h^3')

    
end