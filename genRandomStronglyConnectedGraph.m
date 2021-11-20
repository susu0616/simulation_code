function [ G ] = genRandomStronglyConnectedGraph( N )

G = genRandomGraph(N,rand,rand);

% Add directed edges or reorient edges until strongly connected
while(1)
    [S, ~] = graphconncomp(sparse(G));
    if(S == 1)
        break;
    else
        idx = find(G+eye(N) == 0);
        if(isempty(idx))
            G = genRandomGraph(N,rand,rand);
        else
            G(idx(randi(numel(idx)))) = 1;
        end
    end        
end

end