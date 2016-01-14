function [ ] = Demo_Integrator(varargin )
%% Setup and parse parameters to Demo function

%set default values
approx = 0.33; %default approximation level, keep only one third of data
approxRankFixed = true; %default if approximation level is fixed during approximation

inp = inputParser;
addParameter(inp,'approx',approx,@isnumeric);
addParameter(inp,'rank',0,@isnumeric);
addParameter(inp,'fixed',approxRankFixed,@islogical);

parse(inp, varargin{:});
approx = inp.Results.approx;
approxAsRank = approx > 2;
if inp.Results.rank > 0
    %explicitly wants to use given rank for approximation
    approx = inp.Results.rank;
    approxAsRank = true;
end
approxRankFixed = inp.Results.fixed;

%% Parameters of integrator itself

close all hidden
n = 100;
h = 1e-4;
endT = 1;%start time is zero
timeSteps = round(endT / h);
pertubation = 0;% 1e-6;

A2 =  Generate_Low_Rank(n, pertubation);

T1 = Generate_Skew_Symmetric(n); T1I = T1+eye(n);
T2 = Generate_Skew_Symmetric(n);
Q1Start = eye(n); %Generate_Orthogonal(n);
Q2Start = eye(n); %Generate_Orthogonal(n);

function DeltaA = Get_Delta(Current, stepIndex)
    DeltaA = h*(T1I * Current + Current * T2');
end
function Frame = Make_Frame(U, S, V, frameIndex)
   Frame = U * S * V';
   Frames(:,:,frameIndex) = Frame;   
end
StartA = Q1Start * (A2) * Q2Start';
Frames = zeros(size(StartA,1), size(StartA,2), timeSteps);

approxParam = 'approx';
if approxAsRank
    approxParam = 'rank';
end
r=DLR(StartA, timeSteps, @Get_Delta, @Make_Frame, approxParam, approx, 'fixed', approxRankFixed);

function A = Get_A_At_Time(time)
    Q1 = Get_Q_At_Time(T1, Q1Start, time, h);
    Q2T = Get_Q_At_Time(T2, Q2Start, time, h)';
    A =  Q1 * (exp(time) * A2) * Q2T;
end

%% Calculate and visualize errors
figure();
errors = zeros(timeSteps, 1);
for i = 1:timeSteps
    CurrA = Get_A_At_Time(h * (i - 1)); 
    errors(i) = norm(Frames(:,:,i) - CurrA, Inf);  
end
times = (1:timeSteps) * h;
semilogy(times, errors);
disp(strcat('Max error=', num2str(max(errors)), ' min error=', num2str(min(errors))));

end

function Q = Get_Q_At_Time(T, Q0, time, dt)
    if time == 0
        Q = Q0;
        return
    end
    %dQ/dt = T*Q and Q(t) = Q0 at start t = 0
    Q = expm(T*time)*Q0;
end

function A = Generate_Low_Rank(n, pertubation)
    A = zeros(n);
    r = round(min(n / 2, 10));
    %Generate matrix A of rank r
    A(1:r, 1:r) = eye(r) + rand(r, r) * 0.5;
    
    %Add some pertubation to total A
    A = A + rand(n, n) * pertubation;
end
function T = Generate_Skew_Symmetric(n)
    T = rand(n,n);
    
    %T = zeros(n);
    %T = T + diag(ones(n-1,1), -1);
    
    T = T + T'; % now symmetric
    for i = 1:n %rows
        for j = 1:n %columns
            if i == j
                T(i, j) = 0;
            elseif i < j
                T(i, j) = -T(i, j);
            end
        end
    end
end

function Q = Generate_Orthogonal(n) 
    Q = [];
    %Random matrices mostly are of full rank, but still we ensure Q to be
    %nxn
    while min(size(Q)) < n
        Q = orth(rand(n, n));
    end
end