function [Cost, Power] = f1(R)
%%  Sphere function (f1)
Cost = sum(R .^ 2, 2);
Power = 1/Cost;
end