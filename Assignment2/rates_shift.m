function [ratesSet_shift] = rates_shift(ratesSet)
% This functions shifts the whole rates curve by 1bp (0.0001)
%
% INPUT
% ratesSet: struct containing bid and ask rates of the instruments
% OUTPUT
% ratesSet_shift: new rates data shifted by 1bp

ratesSet_shift = ratesSet;

% Get the number of fields in the struct
num_fields = numel(fieldnames(ratesSet));

% Iterate over each field of the struct
for field_index = 1:num_fields
    % Extract the name of the field
    field_name = fieldnames(ratesSet);
    
    % Extract the matrix associated with that field
    matrix = ratesSet.(field_name{field_index});
    
    % Shift each value of the matrix by adding +0.0001 (1bp)
    shifted_matrix = matrix + 0.0001;
    
    % Assign the shifted matrix to the corresponding field in the output struct
    ratesSet_shift.(field_name{field_index}) = shifted_matrix;
end

end     % function rates_shift