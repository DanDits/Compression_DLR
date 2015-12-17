function [ ] = Demo_Paper(varargin)
%Demo for the example in the paper involving the matrix function
%A(t)=Q1(t)*(A1+exp(t)A2)*Q2(t)' solved by choosing two different DeltaA
%for the DLR routine: traditional with DeltaA = A(t)-A(t-1) or my method
%with DeltaA = A(t) - Y_{t-1}.
%Hereby A1 and A2 and therefore A are low rank (at most rank 10) except for
%a little pertubation that can be controlled and (sightly) increased the
%rank. Q1 and Q2 is a series of orthogonal matrices. Test using
%approximation of rank 10 or 20 and pertubation of 1e-6, 1e-3 or 0 for
%validating the results in the paper on page 184.

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

%% Start the setup and calculation of the function A(t)
close all hidden

n = 100;
h = 1e-3;
endT = 1;%start time is zero
pertubation = 1e-6;
A1 = Generate_Low_Rank(n, pertubation);
A2 =  Generate_Low_Rank(n, pertubation);

T1 = Generate_Skew_Symmetric(n);
T2 = Generate_Skew_Symmetric(n);
Q1Start = eye(n); %Generate_Orthogonal(n);
Q2Start = eye(n); %Generate_Orthogonal(n);

function DeltaA = Get_Delta_With_Current(Current, stepIndex)
    time = stepIndex * h;
    DeltaA = Get_A_At_Time(time) - Current;

end
function DeltaA = Get_Delta_Traditional(Current, stepIndex)
    time = stepIndex * h;
    DeltaA = Get_A_At_Time(time) - Get_A_At_Time((stepIndex - 1) * h);
end

function A = Get_A_At_Time(time)
    Q1 = Get_Q_At_Time(T1, Q1Start, time, h);
    Q2T = Get_Q_At_Time(T2, Q2Start, time, h)';
    A =  Q1 * (A1 + exp(time) * A2) * Q2T;
end

timeSteps = round(endT / h);
StartA = Q1Start * (A1 + A2) * Q2Start';

approxParam = '';
if approxAsRank
    approxParam = 'asrank';
end
if approxRankFixed
    approxParam = strcat(approxParam, '!');
end
FramesTraditional = zeros(size(StartA, 1), size(StartA, 2), timeSteps);
function Frame = Make_Traditional_Frame(U, S, V, frameIndex)
   Frame = U * S * V';
   FramesTraditional(:,:,frameIndex) = Frame;   
end
FramesWithCurrent = zeros(size(StartA, 1), size(StartA, 2), timeSteps);
function Frame = Make_Frame_With_Current(U, S, V, frameIndex)
   Frame = U * S * V';
   FramesWithCurrent(:,:,frameIndex) = Frame;   
end
approxParam = 'approx';
if approxAsRank
    approxParam = 'rank';
end
DLR(StartA, timeSteps, @Get_Delta_Traditional, @Make_Traditional_Frame, approxParam, approx, 'fixed', approxRankFixed);
DLR(StartA, timeSteps, @Get_Delta_With_Current, @Make_Frame_With_Current, approxParam, approx, 'fixed', approxRankFixed);

%% Plot the error of approximated frames to true frames(A at time steps)
figure();
errorsTraditional = zeros(timeSteps, 1);
errorsWithCurrent = zeros(timeSteps, 1);
for i = 1:timeSteps
    CurrA = Get_A_At_Time(h * (i - 1));
    errorsTraditional(i) = norm(FramesTraditional(:,:,i) - CurrA, 'fro');  
    errorsWithCurrent(i) = norm(FramesWithCurrent(:,:,i) - CurrA, 'fro');  
end
disp(strcat('Traditional Max error=', num2str(max(errorsTraditional)), ' min error=', num2str(min(errorsTraditional))));
disp(strcat('WithCurrent Max error=', num2str(max(errorsWithCurrent)), ' min error=', num2str(min(errorsWithCurrent))));

times = (1:timeSteps) * h;
semilogy(times, errorsTraditional, times, errorsWithCurrent);
title(strcat('Errors of DLR in frobenius norm, dt=', num2str(h), ', approx=', num2str(approx), ', pertubation=', num2str(pertubation)));
legend('With DeltaA=A(t)-A(t-1)', 'With DeltaA=A(t)-Y_t-1');
xlabel('time');
ylabel('error');
end


function Q = Get_Q_At_Time(T, Q0, time, dt)
    if time == 0
        Q = Q0;
        return
    end
    %dQ/dt = T*Q and Q(t) = Q0 at start t = 0
    Q = expm(T*time)*Q0;
    
%     Q= zeros(size(Q0));
%     function result = Ode_Func(t, y)
%         result = T * y;
%     end
%     for i = 1:size(Q,2)
%        [~, QResult] = ode15s(@Ode_Func, [0, time], Q0(:, i));
%        Q(:,i) = QResult(end, :)';
%     end
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