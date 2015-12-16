function [] = Demo(varargin)
%Some globals so that calculation and loading of data does not have to
%happen all the time Demo is restarted
global M
global Frames
global map
global autoReverse
%% Setup and parse parameters to Demo function

availableTypes = {'fractal', 'simple', 'expand'}; %and everything ending with '.gif'

inp = inputParser;
addOptional(inp,'type','',@(x) (length(x)>4 && strcmp(x(end-3:end),'.gif'))...
    || any(validatestring(x,availableTypes)));
addParameter(inp,'approx',0.33,@isnumeric);
addParameter(inp,'rank',0,@isnumeric);
addParameter(inp,'fixed',true,@islogical);
addParameter(inp,'autoplay',true,@islogical);
addParameter(inp,'showbestapprox',false,@islogical);

parse(inp,varargin{:});
approx = inp.Results.approx;
approxAsRank = approx > 2;
if inp.Results.rank > 0
    %explicitly wants to use given rank for approximation
    approx = inp.Results.rank;
    approxAsRank = true;
end
approxRankFixed = inp.Results.fixed;
type = inp.Results.type;
autoPlay = inp.Results.autoplay;
showBestApprox = inp.Results.showbestapprox;


%% Execute Demo with given parameters
close all hidden

if ~isempty(type)
    if strcmp(type,'fractal')
        [M,Frames]=GenerateFractalMovie();MW1=implay(M); map=[]; autoReverse=3;% Generate sequence of frames of a changing fractal
    elseif strcmp(type, 'simple')
        [M,Frames]=GenerateSimpleMovie();MW1=implay(M); map=[]; autoReverse=3;% Generate sequence of frames of a moving rectangle
    elseif strcmp(type, 'expand')
        [M,Frames]=GenerateExpandingMovie();MW1=implay(M); map=[]; autoReverse=1;
    elseif length(type)>4 && strcmp(type(end-3:end),'.gif')
        [M,Frames, map]=GenerateMovieOfGif(type);MW1=implay(M);autoReverse=1; %Extract sequence of frames of a gif file
    end
else
    MW1 = implay(M); %reopen
end

approxParam = 'approx';
if approxAsRank
    approxParam = 'rank';
end
[FramesC, FramesBestApprox] = DLR_Frames(Frames, 'approx', approx, approxParam, approx, 'fixed', approxRankFixed,...
    'ploterror', true, 'plotranks', true, 'showbestapprox', showBestApprox); 
M2 = FramesToMovie(FramesC, map); MW2=implay(M2);

if showBestApprox
    M3 = FramesToMovie(FramesBestApprox, map); MW3=implay(M3);
    set(MW3.Parent,'units','normalized','outerposition',[0.25 0 0.5 1])
    set(MW3.Parent, 'Name', strcat('Best approximation data, approx=', num2str(approx)));
end

%% UI fine tuning
%Automatically set figure windows to left and right side of screen and
set(MW1.Parent,'units','normalized','outerposition',[0 0 0.5 1])
set(MW2.Parent,'units','normalized','outerposition',[0.5 0 0.5 1])
set(MW1.Parent, 'Name', 'Source data');
set(MW2.Parent, 'Name', strcat('Compressed data, approx=', num2str(approx)));
% bring them to the front
figure(MW1.Parent)
figure(MW2.Parent)

if autoPlay
    %and start playing at the same time (using set auto reverse mode)
    StartPlayer(MW1, autoReverse);
    StartPlayer(MW2, autoReverse);
    if showBestApprox
       StartPlayer(MW3, autoReverse); 
    end
end