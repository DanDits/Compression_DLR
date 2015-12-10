function [] = Demo_Heat(GridType, approx)

%For solving the two dimensional heat equation dA/dt - laplace(A) = F(A)
%with dirichlet boundary conditions on an arbitrary domain D we get the
%solution A by first solving -laplace(A) = F(A) with finite differences and
%then propagate the solution in time by calculating dA/dt = F(A) using
%Dynamic low rank approximation (DLR)

close all hidden
if nargin <= 1
    approx = 0.33; %default approximation level, keep only one third of data
end
approxAsRank = '';
if approx > 1
    approxAsRank = 'asrank';
end
%grid inside region [0,1]x[0,1]
T = 5;
tau = 0.05;
slowDownFactor = 5;

n = 128;

function value = F(t, x, y)

    d1 = sqrt((x-0.5)^2 + (y-0.25)^2);
    d2 = sqrt((x-0.5)^2 + (y-0.75)^2);
    r = 0.05;
    if ((d1 <= r && t < 0.33*T) || (d2 <= r && t < 0.66 * T))
       value = min(d1, d2) * 1000000;
    else 
        value = 0;
    end
end
[Frames, X, Y] = DLR_Heat_Integrator(GridType, n, tau, T, @F, approx, approxAsRank);

%% Animate frames in a surface plot
s = surf(X, Y, Frames(:,:,1),'EdgeColor','none','LineStyle','none');
colormap('hot')
colorbar
totalMin = min(min(min(Frames)));
totalMax =  max(totalMin + 1, max(max(max(Frames))));
caxis([totalMin totalMax] ) ; 
mode manual;
zlim([totalMin totalMax]);
view(45,45);
for j = 2:size(Frames,3)
   s.ZData = Frames(:,:,j); 
   title(strcat('Time passed:  ', num2str((j - 1) * tau, '%.1f'), '/', num2str(T, '%.1f')));
   pause(tau * slowDownFactor);
end

end

