clear all
clc
close all

f = @(x) tanh(x);
x0 = -1;
x1 = 1;

h = 0.01767;
x = x0:h:x1
d = zeros(1,length(x));
f_prime = @(x) sech(x).^2;
exact = f_prime(x)

[d] = five_point_method(f,x,h,d)
% [d] = forward_diff_method(f,x,h,d)
% [d] = centered_diff_method(f,x,h,d)
[two_norm] = validation(f,x,h,d,exact,f_prime)

subplot(2,1,1)
plot(x,exact,'bo-',x,d,'rx-');
xlabel('x');
ylabel('f(x)')
legend('Real Function','Predicted function')

subplot(2,1,2)
plot(x,abs(exact-d),'.-');
xlabel('x')
ylabel('error')
title('Absolute Error = |f(x)-f(x^*)|')
legend('Error')

function [two_norm] = validation(f,x,h,d,exact,f_prime)
    two_norm = zeros(1,4);
    for j=1:4
        for i = 1:length(x)  
            d(i) = 1./(12.*h).*(-25.*f(x(i))+48.*f(x(i)+h)-36.*f(x(i)+2.*h)+16.*f(x(i)+3.*h)-3.*f(x(i)+4.*h));
        end
        two_norm(j) = norm(abs(exact-d));
        h = h/2;
        x = -1:h:1;
        d = zeros(1,length(x));
        exact = f_prime(x);
    end
end

function [d] = five_point_method(f,x,h,d)
    for i = 1:length(x)  
        d(i) = 1./(12.*h).*(-25.*f(x(i))+48.*f(x(i)+h)-36.*f(x(i)+2.*h)+16.*f(x(i)+3.*h)-3.*f(x(i)+4.*h));
    end
end

function [d] = forward_diff_method(f,x,h,d)
    for i = 1:length(x)  
        d(i) = (f(x(i)+h)-f(x(i)))/h;
    end
end

function [d] = centered_diff_method(f,x,h,d)
    for i = 1:length(x)  
        d(i) = (f(x(i)+h)-f(x(i)-h))/(2*h);
    end
end