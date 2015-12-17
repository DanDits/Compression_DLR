function [Frames, X, Y, tau] = Demo_Heat(GridType, varargin)

%For solving the two dimensional heat equation dA/dt - laplace(A) = F(A)
%with dirichlet boundary conditions on an arbitrary domain D we get the
%solution A by first solving -laplace(A) = F(A) with finite differences and
%then propagate the solution in time by calculating dA/dt = F(A) using
%Dynamic low rank approximation (DLR)

%% Setup and parse parameters to Demo function

%set default values
availableGrids = {'S', 'L', 'C', 'D', 'A', 'H', 'B', 'N', 'E'}; %see numgrid and numgridEllipse

inp = inputParser;
addRequired(inp, 'GridType', @(x) any(validatestring(x, availableGrids)));
addParameter(inp,'approx',0.33,@isnumeric);
addParameter(inp,'rank',0,@isnumeric);
addParameter(inp,'fixed',true,@islogical);

parse(inp,GridType, varargin{:});
approx = inp.Results.approx;
approxAsRank = approx > 2;
if inp.Results.rank > 0
    %explicitly wants to use given rank for approximation
    approx = inp.Results.rank;
    approxAsRank = true;
end
approxRankFixed = inp.Results.fixed;

%% Start setup of functions for calculation
close all hidden

%grid inside region [0,1]x[0,1]
T = 15; %Time goes from 0 to T
tau = 0.05; % timestep
slowDownFactor = 1; %for visualization: higher to slow down movie
n = 128; %discretization of grid, how many nodes in one dimension

%Some inhomogenity functions F for the right hand side of the heat equation
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

function value = F2(t, x, y)
    frac = t / T;
    maxR = 0.25;
    d1 = sqrt((x-0.25)^2 + (y-0.5)^2);
    d2 = sqrt((x-0.5)^2 + (y-0.5)^2);
    value = 0;
    if d1 <= maxR * sin(frac * 2 * pi)
        value = 200000;
    end
    if d2 <= maxR * sin(frac * 2 * pi + pi/2)
        value = value + 200000;
    end
end

approxParam = 'approx';
if approxAsRank
    approxParam = 'rank';
end
[Frames, X, Y] = DLR_Heat_Integrator(GridType, n, tau, T, @F, approxParam, approx, 'fixed', approxRankFixed);

Animated_Surf(Frames, X, Y, tau, slowDownFactor);

end

