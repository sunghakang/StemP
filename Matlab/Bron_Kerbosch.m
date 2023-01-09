function [MC] = Bron_Kerbosch(adjacency_matrix,version)

% Author: Kumbit Hwang 
% last update: 04/29/2017

% the Bron?Kerbosch algorithm is an algorithm for finding maximal cliques in an undirected graph.
% Find maximal cliques using Bron-Kerbosch algorithm

%   Output: MC is the output matrix that contains the maximal cliques in its 
%   columns.

%   Reference: Bron, Coen and Kerbosch, Joep, "Algorithm 457: finding all cliques
%   of an undirected graph", Communications of the ACM, vol. 16, no. 9, 
%   pp: 575?577, September 1973.
%
%   Reference: Cazals, F. and Karande, C., "A note on the problem of reporting 
%   maximal cliques", Theoretical Computer Science (Elsevier), vol. 407,
%   no. 1-3, pp: 564-568, November 2008.

n = size(adjacency_matrix,2);      % number of vertices
MC = [];            % storage for maximal cliques
R = [];             % currently growing clique
P = 1:n;            % pro spective nodes connected to all nodes in R
X = [];             % nodes already processed

if strcmp(version, '1')
    BK1(R,P,X);
else
    BK2(R,P,X);
end
    
    function [] = BK1(R,P,X)
        if isempty(P) && isempty(X)
            % report R as a maximal clique
            newMC = zeros(1,n);
            newMC(R) = 1;
            MC = [MC newMC'];
        else
            for v = P
                P = setxor(P,v);
                Rnew = [R v];
                Nv = find(adjacency_matrix(v,:));
                Pnew = intersect(P,Nv);
                Xnew = intersect(X,Nv);
                BK1(Rnew, Pnew, Xnew);
                X = [X v];
            end
        end
    end
        
    function [] = BK2 ( R, P, X )

        ignore = [];                        
        if (isempty(P) && isempty(X))
            % report R as a maximal clique
            newMC = zeros(1,n);
            newMC(R) = 1; 
            MC = [MC newMC.'];
        else
            % choose pivot
            ppivots = union(P,X);           % potential pivots
            binP = zeros(1,n);
            binP(P) = 1;                    % binP contains ones at indices equal to the values in P          
            % rows of A(ppivots,:) contain ones at the neighbors of ppivots
            pcounts = adjacency_matrix(ppivots,:)*binP.';  % cardinalities of the sets of neighbors of each ppivots intersected with P
            [ignore,ind] = max(pcounts);
            v_p = ppivots(ind);             % select one of the ppivots with the largest count
            
            for v = intersect(find(~adjacency_matrix(v_p,:)),P)   % all prospective nodes who are not neighbors of the pivot
                P = setxor(P,v);
                Rnew = [R v];
                Nv = find(adjacency_matrix(v,:));
                Pnew = intersect(P,Nv);
                Xnew = intersect(X,Nv);
                BK2(Rnew, Pnew, Xnew);
           
                X = [X v];
            end
        end
        
    end % BKv2
%{
    function [] = BK3 ( R, P, X )
        
    end

%}

end
