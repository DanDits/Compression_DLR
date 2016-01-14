function [Frames, X, Y] = DLR_Heat_Integrator(GridType, GridSize, tau, T, F, varargin)
%Solves dA/dt -laplace(A) = F(A(t)) with A(0)=AStart for discrete timesteps tau
% on Grid defined by GridType (see numgrid or 'E' for ellipse) and GridSize (>0)

%% Input parameter setup and parsing
availableGrids = {'S', 'L', 'C', 'D', 'A', 'H', 'B', 'N', 'E'}; %see numgrid and numgridEllipse

inp = inputParser;
addRequired(inp, 'GridType', @(x) any(validatestring(x, availableGrids)));
addRequired(inp, 'GridSize', @isnumeric);
addRequired(inp, 'tau', @isnumeric);
addRequired(inp, 'T', @isnumeric);
addRequired(inp, 'F', @ (f) isa(f, 'function_handle'));
addParameter(inp,'approx',0.33,@isnumeric);
addParameter(inp,'rank',0,@isnumeric);
addParameter(inp,'fixed',true,@islogical);

parse(inp,GridType, GridSize, tau, T, F, varargin{:});
approx = inp.Results.approx;
approxAsRank = approx > 2;
if inp.Results.rank > 0
    %explicitly wants to use given rank for approximation
    approx = inp.Results.rank;
    approxAsRank = true;
end
approxRankFixed = inp.Results.fixed;

%% Integrator parameters and setup
timeSteps = floor(T / tau); % assume T, tau > 0

% Create square grid with indexed points in the set. Points not in the set 
% are zero.
if GridType == 'E'
    ellipseA = 0.5;
    ellipseB = 0.3;
    Grid = numgridEllipse(GridSize, ellipseA, ellipseB);
else
    Grid = numgrid(GridType, GridSize);
end
%Solve using finite differences on the area [0,1]x[0,1]
%n=amount of inner points in one direction
n=max(size(Grid))-2;
h=1/(n+1);%grid size

X = 0:h:1;
Y = 0:h:1;

%Calculate the solution! Of problem with 0 dirichlet boundary data
M = delsq(Grid); %5-point-star matrix
innerIndicies = Grid>0;
%use 1/(n+1.5)^2 instead of h^2 for better results?!
hCorrected = 1/(n+1.5);

function DeltaA = Get_Delta(Current, stepIndex)
    %for each inner point of the region compute F and fill into f vector
    % to get the right hand side of the poisson equation
    gridSize = size(Grid);
    fVector = zeros(size(M,1),1);
    for i = 1:gridSize(1)
        for j = 1:gridSize(2)
            index = Grid(i, j);
            if index > 0
               %inner point of region
               %evaluate F for this point
               fVector(index) = F((stepIndex - 1) * tau, (j - 1) * h, (i - 1) * h);
            end
        end
    end
    DeltaA = tau * (Solve_Laplace(M, fVector, innerIndicies, hCorrected, Grid) - Current);
    
%    DeltaA = tau * (Get_Laplace(M, Current, innerIndicies, hCorrected, Grid) - Current);
end
approxParam = 'approx';
if approxAsRank
    approxParam = 'rank';
end

StartA = zeros(size(Grid));

frameCount = timeSteps + 1;
Frames = zeros(size(StartA, 1), size(StartA, 2), frameCount);
function Frame = Make_Simple_Frame(U, S, V, frameIndex)
    Frame = U * S * V';
    Frames(:,:,frameIndex) = Frame;
end

DLR(StartA, frameCount, @Get_Delta, @Make_Simple_Frame, approxParam, approx, 'fixed', approxRankFixed);
end

function A = Solve_Laplace(M, fVector, innerIndicies, h, Grid)
u = M\fVector;
A = zeros(size(Grid));
A(innerIndicies) = h^2 * full(u(Grid(innerIndicies))); 

end

function A = Get_Laplace(M, CurrA, innerIndicies, h, Grid)
currU = zeros(length(innerIndicies), 1);
currU(Grid(innerIndicies))=CurrA(innerIndicies);

currU = M * currU; %laplace of current
A = zeros(size(Grid));
A(innerIndicies) = h^2 * full(currU(Grid(innerIndicies))); 
values='VALUES'
norm(A-CurrA)
max(max(A))
max(max(CurrA))
end
