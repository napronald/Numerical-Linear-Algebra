clear all
clc
close all

f = @(x) sin(x);
x0 = 0;
xn = 2*pi;

% f = @(x) 1./(x.^2+1);
% x0 =-2;
% xn=2;

N = 5; % going to need to zoom in to see real function as N increases
num_predicted_points = 60;

data_x = linspace(x0,xn,N);
data_y = f(data_x);
data_z = linspace(x0,xn,num_predicted_points);

[coefficents] = inter_poly(data_x,data_y)

[pn] = predicted_points(data_x,data_z,coefficents)

subplot(2,1,1)
fplot(f,[x0,xn])
hold on
plot(data_x,f(data_x),'bo',data_z,pn,'rx-');
xlabel('z');
ylabel('f(z), P_n(z)')
titletxt = sprintf('f(z) = --,\t P_n(z) = xxx,\t n+1 = %d',N);
title(titletxt)
legend('Real Function','Actual Data Points','Predicted points')

subplot(2,1,2)
plot(data_z,abs(f(data_z)-pn),'.-');
xlabel('z')
ylabel('|f(z)-P_n(z)|')
title('Absolute Error = |f(z)-P_n(z)|')
legend('Error')

function [coefficents] = inter_poly(data_x,data_y)
    n = length(data_x);
    N = zeros(n,n);
    N(:,1)=data_y;
    for i = 2:n
        for j = i:n
            N(j,i) = (N(j,i-1) - N(j-1,i-1)) / (data_x(j) - data_x(j-i+1));
        end
    end
    coefficents = diag(N);
end

function [pn] = predicted_points(data_x,data_z,coefficents)
    pn = coefficents(1);
    syms x 
    x = data_z;
    n = length(data_x);
    for j=2:n
        ai=1;
        for k=1:j-1
            ai=ai.*(x-data_x(k));
        end
        pn =pn+coefficents(j).*ai;
    end
    pn;
end