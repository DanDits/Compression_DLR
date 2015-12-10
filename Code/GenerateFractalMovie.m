function [Movie, Frames] = GenerateFractalMovie()
%% Setup parameters
fractalValueBound = 30; %When iterations exeeds this, point not in fractal, >= 2
maxIt = 50; %Max amount of iterations to check if inside fractal

%Which part and what stepsize to see of the fractal region
h=0.01;
coordRange = -1:h:1-h;
cRange = -1.5:0.05:-0.5;

fracDim = length(coordRange);
Frames = zeros(fracDim, fracDim, length(cRange));
index = 1; 
dispstat('','init')
dispstat(sprintf('Frame generation started...'),'keepthis','timestamp');
for c = cRange
    fractalFunc = @(z) z*z + c; % Mandelbrot fractal function, takes complex values z 
    Frames(:,:,index) = generateFractalImage(fractalFunc, fractalValueBound, maxIt, coordRange);
    index = index + 1;
    dispstat(sprintf('Progress %d%%',min(100, floor(index / length(cRange) * 100))),'timestamp');
end
dispstat('Finished frame generation. Now generating movie...','keepprev');
Movie = FramesToMovie(Frames);
disp('Movie generated');
end

function Frac = generateFractalImage(fractalFunc, fractalValueBound, maxIt, coordRange)
fracDim = length(coordRange);
Frac = zeros(fracDim);
%% Do the fractal iteration for each point in the coordinate range
for i = 1:fracDim
    x = coordRange(i);
    for j = 1:fracDim
        y = coordRange(j);
        z = x+1i*y;
        N = fractalIterations(fractalFunc, fractalValueBound, maxIt, z);
        Frac(i, j)=N;
    end
end
% points in the fractal get the highest number
Frac(Frac==-1)= maxIt + 1; 

end

%% helper functions
function [N] = fractalIterations(fracFunc, valueBound, maxIt, z)
    N=0;
    while abs(z) < valueBound && N < maxIt
       N = N+1;
       z = fracFunc(z);
    end
    if abs(z) < valueBound
        N = -1; % Element of the fractal as value is still in bound after maxIt
    end
end