function unique = UniquifyDimension(array, dim)
if dim == 1
    array = sortrows(array);
    
    unique = zeros(1,size(array,2));
    currentRow = array(1,:);

    for rowIndex = 2:size(array,1)
       if (array(rowIndex,1) - array(rowIndex-1) > 0)
           unique = [unique; currentRow];
           currentRow = zeros(1,size(array,2));
           currentRow = currentRow + array(rowIndex,:);
       else
           currentRow = CombineRows(currentRow, [0, array(rowIndex,2:end)]);
       end
    end

    % pesky ends
    unique = [unique; array(end,:)];
    unique = unique(2:end,:);
    
elseif dim == 2
    array = array.';
    unique = UniquifyDimension(array, 1);
    unique = unique.';
end

end