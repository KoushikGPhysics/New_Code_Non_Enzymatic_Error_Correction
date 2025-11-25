

function[s]=unidirectional_s_matrix(n,E_position)

% ---------------------------------%
%%% ===Function for making state matrix of the time continuous markov
%%% chain===%%%
%%%  0 = No base pair    o-----o------o
%%%  1 = correct base pair 
%%%  2 = incorrect base pair

%%% === Input ===%%%
%%% n = Number of Base pair
%%% E_position = Position of Error 

%%% === Output ===%%%
%%% s= state matrix
% ----------------------------------%

%%% initial state =00000
%%% final states=11111/11211
% -------------------------------%

s_position = 1;

s = dec2bin(0:2^n-1) - '0';


if E_position==s_position
    disp('Not allowed');
else
    [row,~]=find(s(:,E_position)==1);
    s_temp=s;
    s_temp(row,E_position)=2;
    s_temp=s_temp(row,:);
    s=[s;s_temp];
end

% % % % % % s=[s;
% % % % % %     5,5,5,5,5;
% % % % % %     6,6,6,6,6];

end