function s = s_matrix_full(n)
% ---------------------------------------------------- %
% State space for continuous-time Markov chain
%
% State encoding:
%   0 = no base pair
%   1 = correct base pair
%   2 = incorrect base pair
%
% Position 1:
%   - Primer
%   - Always 1 (non-dynamic)
%
% Positions 2...n:
%   - Can be 0, 1, or 2
% ---------------------------------------------------- %

if n < 1
    error('n must be >= 1');
end

% Number of dynamic positions
nd = n - 1;

% Total number of states
numStates = 3^nd;

% Generate ternary representations for positions 2...n
s_rest = dec2base(0:numStates-1, 3, nd) - '0';

% Prepend primer = 1
s = [ones(numStates,1), s_rest];

end
