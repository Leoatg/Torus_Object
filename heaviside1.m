
%redefine the heaviside function

function Y = heaviside1(X)
%HEAVISIDE    Step function.
Y = zeros(size(X));
Y(X > 0) = 1;
Y(X == 0) = .5;
end