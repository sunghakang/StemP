function value = base_pairing(current_point,comparison)

% base_pair
% input:  current_point and comparison that represents two bases 
% output: 0/1 whether two input bases are canonical base pairs

value = 0;
if (current_point=='A')&&(comparison=='T'), value = 1; end
if (current_point=='T')&&(comparison=='A'), value = 1; end
if (current_point=='A')&&(comparison=='U'), value = 1; end

if (current_point=='U')&&(comparison=='A'), value = 1; end
if (current_point=='C')&&(comparison=='G'), value = 1; end
if (current_point=='G')&&(comparison=='C'), value = 1; end
if (current_point=='a')&&(comparison=='t'), value = 1; end
if (current_point=='t')&&(comparison=='a'), value = 1; end
if (current_point=='a')&&(comparison=='u'), value = 1; end
if (current_point=='u')&&(comparison=='a'), value = 1; end
if (current_point=='c')&&(comparison=='g'), value = 1; end
if (current_point=='g')&&(comparison=='c'), value = 1; end


end