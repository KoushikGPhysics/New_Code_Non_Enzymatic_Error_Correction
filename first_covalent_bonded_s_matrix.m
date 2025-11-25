

function[matrix]=first_covalent_bonded_s_matrix()

C1 = 1;                        % First column always 1
C2 = [0 1];                   % Possible values for column 2
C3 = [0 1 2];                % Possible values for column 3 (only column where 2 allowed)
C4 = [0 1];                   % Possible values for column 4
C5 = [0 1];                   % Possible values for column 5
C6 = [0 1];  


matrix = [];

% Generate all combinations
for c2 = C2
    for c3 = C3
        for c4 = C4
            for c5 = C5
                for c6=C6
                      matrix = [matrix; C1, c2, c3, c4, c5,c6];
                end
            end
        end
    end
end


matrix_1=[5,5,5,5,5,5; 6,6,6,6,6,6];


matrix=[matrix;matrix_1];


end


