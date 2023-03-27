clear all
clc
close all

% f = @(x) x.^3-25
% a = 2
% b = 3
% tolerance = 0.2


% f = @(x) x.^3; 
% a = 1; 
% b = 3;

% f = @(x) x.^3; 
% a = -3; 
% b = 1;

% f = @(x) sin(x)
% a = -pi/2;
% b = 3*pi/4;

f = @(x) x-3;
a=2;
b=5;

tolerance = 10^-5;
iteration = 0;

x = a:0.1:b;
xaxis = zeros(length(x),1);

figure(1)
plot(x, f(x), 'b-s', x, xaxis, 'k-')
grid on
xlabel('x')
ylabel('y')
%bisection method
if (f(a)*f(b))<0
    p_old = (a);
    error1 = (b-a)/2
    while error1 > tolerance || error2 > tolerance
%         if iteration == 50
%             break
%         end
        iteration = iteration + 1;
        p_new = (a+b)/2;
        if (f(a)*f(p_new)) < 0
            b = p_new;
            p_new
            error1 = abs(f(p_new))
            iteration;
        elseif (f(a)*f(p_new)) > 0
            a = p_new;
            p_new
            error1 = abs(f(p_new))
            iteration;
        end
        format long e
        p_old
        p_new
        error2 = abs(p_new-p_old)
        error(iteration) = error2;
        figure(2)
        num_iter = 1:1:iteration;
        loglog(num_iter,error,'g-s')
        grid on
        xlabel('number of iterations')
        ylabel('error')
        p_old = p_new
        if error2 == 0
            break
        end
    end
    root = p_new
else
    fprintf('\n No root found %\n')
end


fixed_point()
function fixed_point
    %%No interval possible, |g'(x)| > 1
    % g = @(x) x.^5+5*x.^3-x.^2+1+x 
    
    %%No interval guaranteed by convergence theorem
    % g = @(x) (-5*x^3+x^2-1)^1/5
    
    %%No interval guaranteed by convergence theorem
    % g = @(x) -x^7-6*x^5+x^4-5*x^3+x-1
    
    
    %%Converges on [1.213,2] guaranteed by convergence theorem
    % g = @(x) x^4-8*x^3+24*x^2-32*x+16+x 
    
    %%guaranteed by convergence theorem CONVERGES [-1,2]
    % g = @(x) ((x^4+24*x^2-32*x+16)/8)^(1/3) 
    
    %%CONVERGES on [-1,2] guaranteed by convergence theorem
    g = @(x) -16/(x.^3-8*x.^2+24*x-32) 
    
    
    x0 = 1;
    tol = 10^-6;
    max_iter = 50;
    
    fixed_pt_method(g,x0, tol, max_iter);
    function root = fixed_pt_method(g, x0, tol, max_iter)
        iter = 0;
        while ( iter < max_iter);
            x_next = g(x0);
            err = abs(x_next-x0)
            iter = iter + 1
            if (err < tol)
                root = x_next
                break
            else
                x0 = x_next
    
            error(iter) = err;
            figure(3)
            num_iter = 1:1:iter;
            loglog(num_iter,error,'g-s')
            grid on
            xlabel('number of iterations')
            ylabel('error')
            end
        end
    end
end




