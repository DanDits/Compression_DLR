function [U,S,V] = Get_Rank_Approx(A, r)
% Calculate singular value decomposition and use spectral cutoff to get
% rank r

[U,S,V] = svd(A, 'econ');

%Sort to cut off greatest eigenvalues (absolute value) 
% and make matrix rank <= r, eigenvalues of svd are already sorted
%[S,SortInd]=sort(abs(diag(S)), 'descend');
%SortInd = SortInd(1:r);

S = diag(S);
SortInd = 1:r;

%Spectral cut off and set new smaller U,S,V
U = U(:,SortInd);
S = diag(S(SortInd));
V = V(:,SortInd);

end

