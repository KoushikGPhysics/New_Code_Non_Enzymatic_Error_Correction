function [prob_all, prob_ter, steady_state] = steady_state_probability(t, initial_state, final_state, Q)

% --------- PROBABILITY EVOLUTION OF MARKOV CHAIN ----------- %
% INPUT:
% t = time vector
% initial_state = index of starting state
% final_state = indices of absorbing states
% Q = transition rate matrix

% OUTPUT:
% prob_all = probability evolution of all states
% prob_ter = probability evolution of final (absorbing) states
% steady_state = steady-state probabilities (t ->inf)

pi_0 = zeros(1, size(Q,1));
pi_0(initial_state) = 1;

% ------ Matrix Exponential Evolution -------- %
P = zeros(length(t), size(Q,1), size(Q,2));
for i = 1:length(t)
    P(i,:,:) = expm(Q * t(i));
end

% ------ Probability evolution -------- %
pi_evolution = zeros(length(t), size(Q,1));
for i = 1:length(t)
    pi_evolution(i,:) = pi_0 * squeeze(P(i,:,:));
end

prob_all = pi_evolution';  
prob_ter = prob_all(final_state,:);  

% ------ Compute steady state (t -> inf) -------- %
transient = setdiff(1:size(Q,1), final_state);
Q_TT = Q(transient, transient);
Q_TA = Q(transient, final_state);
N = -inv(Q_TT);
B = N * Q_TA;

steady_state = zeros(1, size(Q,1));
if ismember(initial_state, transient)
    steady_state(final_state) = B(initial_state == transient, :);
else
    steady_state(initial_state) = 1;  
end

end
