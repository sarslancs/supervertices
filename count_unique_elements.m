function [ elems, counts ] = count_unique_elements( v )
%COUNT_UNIQUE_ELEMENTS Counts the number of unique elements in an array
%   Counts the number of unique elements in the given array V and sorts the
%   quantity of each element in descending order. Returns them in COUNTS, 
%   with the unique elements in ELEMS
  
    y = sort(v);
    p = find([true;diff(y)~=0;true]);
    elems = y(p(1:end-1));
    counts = diff(p);
    [counts, idx] = sort(counts, 'descend');
    elems = elems(idx);   
end

