function [G] = genRandomGraph(N,p,q)

G = zeros(N);

% We construct the graph iteratively using the following operations:
% * New edge in existing node set (p), pick random nodes and edge direction
% (q)
% * New edge not in existing edge set, pick edge direction (r) and node in
% existing graph, and add new node.

k = 1;
while(k < N)
    
    % New edge spans the current node set?
    if(rand > p && k > 1)
        i = randi(k);
        nodes = 1:k;
        nodes(i) = [];
        j = nodes(randi(numel(nodes)));      
        G(i,j) = 1;             
    else
        i = randi(k);
        
        % New edge direction
        if(rand > q)
            G(i,k+1) = 1;
        else
            G(k+1,i) = 1;
        end
        
        k = k+1;
    end
end