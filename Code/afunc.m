function [ y ] = afunc( x)

%Ackley's
a = 20;
b = 0.2;
c = 2*pi;
%d = 1;
%xrange = 100
%n = xrange^d;
n = 100;
%y= -a*exp(-b*sqrt(1/n*sum(x.^2))) - exp(1/n*sum(cos(c*x))) + a + exp(1);

%De Jong's
%y = sum(x.^2);

%Rastringin's
y = 10*n + sum(x.^2 - 10*cos(2*pi*x));
end