function [] = Demo_GameOfLife(approx, startRandom)
%Versuchsreihe für n=100, startRandom=false (=von 1/3 bis 2/3 rechteck
%lebendig starten), evolutionSteps = 100
%(31,false) ->Exakte Berechnung
%(30,false) ->Ziemlich gut, minimale Mutation
%(29,false) ->Größere Mutation
%(20,false) ->Komplett verwuchert, Regeln erkennbar verletzt

%Versuchsreihe wie oben nur mit startRandom=true
%(80,true) -> Teilweise Regeln erkennbar, viel Chaos
%(90,true) -> Einzelne Mutationen, im groben Augemaß richtig
%(95,true) -> Exakte Berechnung
close all hidden
n = 100;
map = [0.2 1. 0.1; 0.9 0.1 0.1];
evolutionSteps = 100;
if approx < 1
    approx = 1;
end
approx = round(approx);

if startRandom
    StartA = rand(n);StartA(StartA>=0.5)=1;StartA(StartA<0.5)=0;%Random, dies very fast
else
    StartA = zeros(n); StartA(floor(n/3):floor(2/3*n), floor(n/3):floor(2/3*n))=1;
end
TrueSolutions = zeros(size(StartA,1), size(StartA,2), evolutionSteps);
LastSolutionBase = StartA;
TrueSolutions(:,:,1)=Make_Display_Frame(StartA);
Frames = zeros(size(StartA,1), size(StartA,2), evolutionSteps);
function DeltaA = Get_Delta(Current, stepIndex)
    DeltaA = zeros(size(Current));
    for i = 1:size(Current,1)
        for j = 1:size(Current,2)
            livingCount = Get_Living_Neighbors(Current, i, j);
            if Current(i,j) >= 0.5 %current is alive
                if livingCount > 3 || livingCount < 2
                    DeltaA(i, j)=-Current(i,j); %DIE!
                else %Stay alive
                    DeltaA(i, j)=1-Current(i,j);
                end
            else %current is dead
               if livingCount == 3 
                  DeltaA(i, j)=1-Current(i,j); %SPAWN! 
               else %Stay dead
                   DeltaA(i, j)=-Current(i,j);
               end
            end
        end
    end
    NextSol = LastSolutionBase + DeltaA;
    TrueSolutions(:,:,stepIndex+1)=Make_Display_Frame(NextSol);
    LastSolutionBase = NextSol;
end

function Frame = Make_Frame(U, S, V, frameIndex)
   Frame = U * S * V';
   Frames(:,:,frameIndex) = Make_Display_Frame(Frame);
end


r=DLR(StartA, evolutionSteps, @Get_Delta, @Make_Frame, 'rank', approx, 'fixed', true);
M = FramesToMovie(Frames, map); 
MW=implay(M);

MT = FramesToMovie(TrueSolutions, map);
MWT=implay(MT);

errors = zeros(length(evolutionSteps), 1);
ranksT = zeros(length(evolutionSteps), 1);
for i = 1:evolutionSteps
   errors(i)= norm(Frames(:,:,i)-TrueSolutions(:,:,i),'fro');
   ranksT(i)=rank(TrueSolutions(:,:,i));
end
figure;
plot(1:evolutionSteps, errors);
title('Estimated error in Frobeniusnorm');

figure;
plot(1:evolutionSteps, r, 1:evolutionSteps, ranksT);
legend('Approx rank', 'True rank');
end % Demo_GameOfLife

function count = Get_Living_Neighbors(A, i, j)
    count = 0;
    lastRow = size(A,1);
    lastColumn = size(A,2);
    for di = [-1 0 1]
        %Exlude (i,j), not its own neighbor
       currI = i + di;
       if currI == 0 %top neighbor of first row is bottom row
           currI = lastRow;
       elseif currI == lastRow + 1
           currI = 1;
       end
       
        for dj = [-1 0 1]
            if di ~= dj || di ~= 0
               currJ = j + dj;
               if currJ == 0 %left neighbor of first column is right column
                   currJ = lastColumn;
               elseif currJ == lastColumn + 1
                   currJ = 1;
               end
               
               if A(currI, currJ) >= 0.5
                  count = count + 1; 
               end
            end
        end
    end
end

function Frame = Make_Display_Frame(Frame)
   Frame(Frame>=0.5)=1;
   Frame(Frame<0.5)=2;
end
