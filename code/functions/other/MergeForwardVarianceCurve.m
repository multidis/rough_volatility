function idxGrp = MergeForwardVarianceCurve(Txi,Tobs)
% Description: Merges a piecewise flat sections of a forward variance curve 
% to obtain a curve with longer flat sections all depending on the expirations 
% that are observed. 
%
% Parameters:
%   Txi:  [Nx1 real] Forward variance curve grid points. I.e. the curve is 
%         (initially) assumed flat from the maturity 0 to Txi(1) and again
%         from Txi(1) to Txi(2) ... etc.
%   Tobs: [Mx1 real] Unique observed expirations.
%
% Output:
%   idxGrp: [Nx1 integer] Indices of the Txi vector specifying which
%           sections are grouped together.
% 

% Determine flat sections:
idxGrp = NaN(size(Txi));
idxGrp(end) = 1;
T = [0;Txi];
N = size(T,1)-1;
for i=N-1:-1:1
    if ~any(Tobs > T(i) & Tobs <= T(i+1))
        idxGrp(i) = idxGrp(i+1);
    else
        idxGrp(i) = idxGrp(i+1) + 1;
    end
end

% Reverse group values:
n_grps = max(idxGrp);
idxGrp = n_grps - idxGrp + 1;

end