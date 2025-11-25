
function[prob_all,prob_ter]=P_Evolution(t,initial_state,final_state,Q)

% % --------- PROBABILITY EVOLUTION OF MARKOV CHAIN ----------- %

% INPUT: t= Time , R=Row number of tergate states, Q=Transiiton Matrix

% OUTPUT: prob_all = probability evolution of all states, prob_ter
% prob_ter= probability evolution of terminal states.



pi_0=zeros(1,size(Q,1));
pi_0(1,initial_state)=1;

% ------Matrix Exponential --------%
P = zeros(length(t), size(Q, 1), size(Q, 2));
for i = 1:length(t)
    P(i, :, :) = expm(Q * t(i));
end

% -------- Multplying by initial probability distribution -----------%
pi_evolution = zeros(length(t), size(Q, 1));
for i = 1:length(t)
    pi_evolution(i, :) = pi_0 * squeeze(P(i, :, :));
end

prob_all=pi_evolution';                             % Probability Evolution of all States

prob_ter=prob_all(final_state,:);              % Probability Evolution of Final States


end

