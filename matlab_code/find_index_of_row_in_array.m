function [index] = find_index_of_row_in_array(row_2_index_in_array, array_with_row)
%{
Given an [nx3] array, return the first index where row matches array.
%}
index = find(ismember(array_with_row,row_2_index_in_array,'rows'));
index = index(1,1);
end