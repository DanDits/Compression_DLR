function [] = Animated_Contour(Frames, X, Y, tau, slowDownFactor)
%% Animate frames in a surface plot
N = 20; %amount of contour lines

T = tau * (size(Frames,3) - 1);
for j = 1:size(Frames,3)
   clf
   mesh(X, Y, Frames(:,:,j));
   title(strcat('Time passed:  ', num2str((j - 1) * tau, '%.1f'), '/', num2str(T, '%.1f')));
   pause(tau * slowDownFactor);
end

end

