function G = numgridEllipse(n, a, b) 
%Code adapted from https://github.com/AustinChou/Computational_Mathematics_HT2014/blob/master/project2/ellipse_ev.m

%Curvature of ellipse is k(x)=ab/(a^2sin(x)^2+b^2sin(x)^2)^(3/2) 
%so max and min curvature are a/b^2 and b/a^2
%(for a=0.5 and b=0.25 max(k)=8 and min(k)=1)
  %% make a rectangular grid that will contain the ellipse
  x1 = -a;  x2 = a;
  y1 = -b;  y2 = b;
  % now expand it a bit: add a margin of roughly 10%
  expandFactor = 0.1;
  w = min(x2-x1, y2-y1);
  x1 = x1-expandFactor*w;
  x2 = x2+expandFactor*w;
  y1 = y1-expandFactor*w;
  y2 = y2+expandFactor*w;
  %x1 = floor(x1-expandFactor*w)    % or could round to integers
  %x2 = ceil(x2+expandFactor*w)

  dx = (x2-x1) / n;
  x1d = x1:dx:x2;
  y1d = y1:dx:y2;  % note last y pt isn't y2
  [x,y] = meshgrid(x1d,y1d);


  %% An implicit representation of the shape
  % negative inside and positive outside.
  % ellipse: x^2/a^2 + y^2/b^2 = 1
  A = x.^2 / a^2 + y.^2 / b^2 - 1;
  % or a rectangle:
  %A = max(abs(x)-a, abs(y)-b);

  %% Build a grid
  % this is from numgrid.m: we choose an ordering for the points
  % inside the ellipse (all other points are labeled 0).
  G = A < 0;             % first, label points inside with a 1.
  k = find(G);           % now find the linear index of these.
  G = zeros(size(x));    % new all zero matrix
  G(k) = (1:length(k))'; % label them inside ones from 1 upwards

  %ensure G is square, center G in temporary H that is square
  n = max(size(G));
  H = zeros(n,n);
  columnsG = size(G,2);
  startColumn = floor((n-columnsG) / 2) + 1;
  rowsG = size(G,1);
  startRow = floor((n-rowsG) / 2) + 1;
  H(startRow : startRow + rowsG - 1, startColumn : startColumn + columnsG - 1)=G;
  G=H;
  
  %spy(G)
  %axis equal
