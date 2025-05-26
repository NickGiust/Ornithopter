function combined = CombineRows(A, B)
combined = ((A + B) .* (((A ~= 0) + (B ~= 0)) < 2)) + (((A + B) .* (((A ~= 0) + (B ~= 0)) == 2))/2);
end