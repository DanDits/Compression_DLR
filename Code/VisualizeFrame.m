function [] = VisualizeFrame(Frame, map)
%% Plot the fractal with default color map
figure;
image(Frame)
if nargin > 1
    colormap(map);
end