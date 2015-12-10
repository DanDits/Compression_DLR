function [] = VisualizeFrame(Frame)
%% Plot the fractal with default color map
figure;
image(Frame)
colormap('jet')
colorbar