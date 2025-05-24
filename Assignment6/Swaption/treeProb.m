function tree = treeProb(tree, p_u, p_m, p_d)
% treeProb  Assigns transition probabilities p to each node of a trinomial tree
%
%   tree = treeProb(tree, p_u, p_m, p_d)
%
%   Inputs:
%     tree  –  a cell array {1…N} of node-arrays, each node with field .x
%     p_u   –  “up” probability
%     p_m   –  “middle” probability
%     p_d   –  “down” probability
%
%   Output:
%     tree  –  same structure but with an added field .p on each node
%
%   The root node (time 0) is seeded with .p = 1, and thereafter
%   each parent’s .p is split among its (up,mid,down) children.
%   If the tree has stopped growing (barrier reached), we simply
%   copy probabilities forward.
% GRANDE CHAT

    N = numel(tree);

    % seed the root
    tree{1}(1).p = 1;

    for i = 2:N
        prev = tree{i-1};
        curr = tree{i};
        nPrev = numel(prev);
        nCurr = numel(curr);

        % initialize all probabilities on this layer to zero
        [curr.p] = deal(0);

        if nPrev == nCurr
            % barrier case: no new nodes, just carry forward
            for k = 1:nPrev
                curr(k).p = prev(k).p;
            end
        else
            % normal branching: each parent spreads p to 3 children
            for j = 1:nPrev
                pj = prev(j).p;
                curr(j    ).p = curr(j    ).p + pj * p_u;
                curr(j + 1).p = curr(j + 1).p + pj * p_m;
                curr(j + 2).p = curr(j + 2).p + pj * p_d;
            end
        end

        % (optional) sanity‐check: layer sum must still be 1
        % total = sum([curr.p]);
        % assert(abs(total-1) < 1e-12, 'Layer %d prob sum = %f ≠ 1', i, total);

        tree{i} = curr;
    end
end
