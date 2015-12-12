function [] = Demo(type, approx, fixedRank)
%Interesting parameters: 'gif3.gif', 0.33 and then 0.65 -> großes Bild
%   'gif2.gif' -> unstetige indizierte Farben
%   'gif1.gif' -> maximaler Rang
global M
global Frames
global map
global autoReverse
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
if ~isempty(type)
    if strcmp(type,'fractal')
        [M,Frames]=GenerateFractalMovie();MW1=implay(M); map=[]; autoReverse=3;% Generate sequence of frames of a changing fractal
    elseif strcmp(type, 'simple')
        [M,Frames]=GenerateSimpleMovie();MW1=implay(M); map=[]; autoReverse=3;% Generate sequence of frames of a moving rectangle
    elseif length(type)>4 && strcmp(type(end-3:end),'.gif')
        [M,Frames, map]=GenerateMovieOfGif(type);MW1=implay(M);autoReverse=1; %Extract sequence of frames of a gif file
    end
else
    MW1 = implay(M); %reopen
end

FramesC = DLR_Frames(Frames, approx, approxAsRank, 'errors', 'ranks!'); M2 = FramesToMovie(FramesC, map); MW2=implay(M2);

%% UI fine tuning
%Automatically set figure windows to left and right side of screen and
set(MW1.Parent,'units','normalized','outerposition',[0 0 0.5 1])
set(MW2.Parent,'units','normalized','outerposition',[0.5 0 0.5 1])
% bring them to the front
figure(MW1.Parent)
figure(MW2.Parent)
%and start playing at the same time (using set auto reverse mode)
StartPlayer(MW1, autoReverse);
StartPlayer(MW2, autoReverse);