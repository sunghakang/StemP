function [B] = find_basepairingture(V,m,n)

[r,c] = find (V == m); 
B = 0;

if (~isempty(r)) && (V(r(1),3-c(1)) == n)
    B = 1;
end

    
end

