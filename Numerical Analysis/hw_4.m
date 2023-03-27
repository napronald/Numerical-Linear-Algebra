clear all
clc
close all

% f = @(x) cos(x+sqrt(2))+x*((x/2)+sqrt(2));
% derf = @(x) x-sin(x+sqrt(2))+sqrt(2);

% f = @(x) exp(6*x)+3*(log(2))^2*exp(2*x)-(log(8))*exp(4*x)-(log(2))^3
% derf = @(x) 6*exp(2*x)*(log(2))^2-4*exp(4*x)*(log(8))+6*exp(6*x)

f = @(x) exp(6*x)+1.441*exp(2*x)-2.079*exp(4*x)-0.333;
derf = @(x) 6*exp(6*x)+2.882*exp(2*x)-8.316*exp(4*x);
derf_2 = @(x) 36*exp(6*x)+5.764*exp(2*x)-33.264*exp(4*x);

start = -0.5;
tol = 10^-5;
max_iter = 50;

[newton_root, iter] = newton(f,derf,start,tol,max_iter)
%--------------------------------------------------------------------

p0 = -2;
p1 = -1;

[secant_root,iter] = secant(f,p0,p1,tol,max_iter)

%--------------------------------------------------------------------

start = 0;

[multi_root,iter] = multi(f,derf,derf_2,start,tol,max_iter)

%--------------------------------------------------------------------

function [root, iter] = newton(f,derf,start,tol,max_iter)
    iter = 1;
    while (iter < max_iter)
        p = start - (f(start)/derf(start))
        
        err(iter) = abs(p-start);

        if (err(iter) < tol)
            root = p;
            break
        else
            iter = iter + 1;
            start = p;
        end
    end
    x = err(1:end-1);
    y = err(2:end);
    logx = log(x);
    logy = log(y);
    
    figure(1)
    plot(logx,logy,'r-o')
    grid on
    
    slope =  (logy(end) - logy(end-1)) / (logx(end) - logx(end-1))
end 


function [root, iter] = secant(f,p0,p1,tol,max_iter)
    q0 = f(p0);
    q1 = f(p1);
    iter = 1;
    while iter < max_iter
       p = p1-q1*(p1-p0)/(q1-q0)
       error(iter) = abs(p-p1);
       if error(iter) < tol
           root = p
           break
       else
           iter = iter + 1
           p0 = p1;
           q0 = q1;
           p1 = p;
           q1 = f(p);
       end
    end
    x = error(1:end-1);
    y = error(2:end);
    logx = log(x);
    logy = log(y);
    
    figure(2)
    plot(logx,logy,'r-o')
    grid on
    slope =  (logy(end) - logy(end-1)) / (logx(end) - logx(end-1))
end

function [root,iter] = multi(f,derf,derf_2,start,tol,max_iter)
    iter = 1;
    while (iter < max_iter)
        p = start - (f(start)*derf(start))/(derf(start)^2-f(start)*derf_2(start));
        % For cubic convergence
%         p = start - (2*f(start)*derf(start))/(2*derf(start)^2-f(start)*derf_2(start));

        err(iter) = abs(p-start);
    
        if (err(iter) < tol)
            root = p
            break
        else
            iter = iter + 1
            start = p;
        end 
    end 
    x = err(1:end-1);
    y = err(2:end);
    logx = log(x);
    logy = log(y);
    
    figure(3)
    plot(logx,logy,'r-o')
    grid on
    
    slope =  (logy(end) - logy(end-1)) / (logx(end) - logx(end-1))
end