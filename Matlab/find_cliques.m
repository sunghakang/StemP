function cliques = find_cliques(M)

cliques = {};
column= length(M(1,:));
row = length(M(:,1));

    for i = 1 : column
        x = [];
        for j = 1 : row
            if M(j,i) == 1
                x = [x,j];
            end
        end
        cliques{end+1,1} = x;
    end

