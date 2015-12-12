function [Frames, X, Y, tau] = Demo_Heat(GridType, approx, fixedRank)

%For solving the two dimensional heat equation dA/dt - laplace(A) = F(A)
%with dirichlet boundary conditions on an arbitrary domain D we get the
%solution A by first solving -laplace(A) = F(A) with finite differences and
%then propagate the solution in time by calculating dA/dt = F(A) using
%Dynamic low rank approximation (DLR)

close all hidden
if nargin <= 1
    approx = 0.33; %default approximation level, keep only one third of data
end
if nargin <= 2
    approxRankFixed = '';
elseif strcmp(fixedRank, 'fixed')
    approxRankFixed = '!';
end
approxAsRank = strcat('', approxRankFixed);
if approx > 1
    approxAsRank = strcat('asrank', approxRankFixed);
end
%grid inside region [0,1]x[0,1]
T = 15;
tau = 0.05;
slowDownFactor = 1;

n = 128;

function value = F(t, x, y)
    frac = t / T;
    r = 0.05;
    d1 = sqrt((x-frac)^2 + (y-0.5)^2);
    d2 = sqrt((x-0.5)^2 + (y-frac)^2);
    value = 0;
    if (d1 <= r)
        value = value + 200000;
    end
    if (d2 <= r)
        value = value + 200000;
    end
end
[Frames, X, Y] = DLR_Heat_Integrator(GridType, n, tau, T, @F, approx, approxAsRank);

Animated_Surf(Frames, X, Y, tau, slowDownFactor);

end

