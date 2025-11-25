



function[s]=covalent_s_matrix(n,E_position)

%%% Function for state matrix including the covalent bond formation

% 0 = No base pair
% 1 = correct base pair
%  2 = incorrect base pair

%  state cosisting all 5  is the covalent bonded state corresponding to
%  the right state
% state consisting 5 and one 6 in error position is the covalent bonded
% state corresponding to the wrong state

%%% ==== Input =====%%%

%  n= number of base pairs;
%  E_position= error position;

%%% ==== Output =====%%%
% s = state matrix

%%% -------------------------%%%
%%% initial state=10000
%%% final states = 55555/55655
%%% -----------------------%%%

%%% ====== Reversible States =====%%%

s_position = 1;

s = dec2bin(0:2^n-1) - '0';

f = s(:,s_position) == 1;
s = s(f,:);

if E_position==s_position
    disp('Not allowed');
else
    [row,~]=find(s(:,E_position)==1);
    s_temp=s;
    s_temp(row,E_position)=2;
    s_temp=s_temp(row,:);
    s=[s;s_temp];
end

%%% ====== Absorbing States =====%%%

final_state = 5 * ones(2, n);     
final_state(2,E_position)=6;     

s = [s; final_state];

end






