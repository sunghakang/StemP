function [value] = wobble_pairing(current_point,comparison)
    value = 0;
    if current_point == 'g' && comparison == 'u', value = 1; end
    if current_point == 'u' && comparison == 'g', value = 1; end
    if current_point == 'U' && comparison == 'G', value = 1; end
    if current_point == 'G' && comparison == 'U', value = 1; end
    if current_point == 'g' && comparison == 't', value = 1; end
    if current_point == 't' && comparison == 'g', value = 1; end
    if current_point == 'T' && comparison == 'G', value = 1; end
    if current_point == 'G' && comparison == 'T', value = 1; end
    %if (current_point ~= 'A' && current_point ~= 'U' &&... 
        %current_point ~= 'G' && current_point ~= 'C') ||...
        %(comparison ~= 'A' && comparison ~= 'U' &&... 
       %comparison ~= 'G' && comparison ~= 'C'), value = 1, end