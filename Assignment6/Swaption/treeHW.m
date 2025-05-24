function [tree,Lmax,N_max] = treeHW(a, sigma, ttm, deltaT)

mu = 1-exp(-a*deltaT);

Lmax = floor((1-sqrt(2/3))/mu); 

N_max = 2*Lmax+1; %max number of nodes

sigma_hat = sigma*sqrt((1-exp(-2*a*deltaT))/(2*a));
deltaX = sqrt(3)*sigma_hat;

T = 0:deltaT:ttm;
N = length(T);
tree = cell(N, 1);
    
% At time 0, there's only one node
tree{1}(1).x = 0;

% Build the tree
for i = 2:N
    prevNodes = tree{i-1};

    if length(prevNodes)+2 > N_max
        %fprintf("%d hit the barrier\n",i);
        tree{i} = prevNodes;
    else

        M = 2*(i-1)+1;  % number of nodes at this level
        tree{i} = repmat(struct('x', 0), M, 1);
       
        for j=1:length(prevNodes)
                x_curr = prevNodes(j).x;
    
                tree{i}(j).x = x_curr + deltaX;  % up
                tree{i}(j+1).x = x_curr;             % middle
                tree{i}(j+2).x = x_curr - deltaX;   % down
        end
    end
end

p_u = 0.5*(1/3-mu*Lmax+mu^2*Lmax^2);
p_m = 2/3-mu^2*Lmax^2;
p_d = 0.5*(1/3+mu*Lmax+mu^2*Lmax^2);

p_u_top = 0.5*(7/3-3*Lmax*mu+(Lmax*mu)^2);
p_m_top = -1/3+2*Lmax*mu-(Lmax*mu)^2;
p_d_top = 0.5*(1/3-Lmax*mu+(Lmax*mu)^2);

p_u_bot = 0.5*(1/3+Lmax*mu+(Lmax*mu)^2);
p_m_bot = -1/3-2*Lmax*mu-(Lmax*mu)^2;
p_d_bot = 0.5*(7/3+3*Lmax*mu+(Lmax*mu)^2);

tree = treeProb(tree,p_u,p_m,p_d);
%tree = treeProb2(tree,p_u,p_m,p_d,max_branch,p_u_top,p_m_top,p_d_top,p_u_bot,p_m_bot,p_d_bot);
end