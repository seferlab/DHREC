function [ids, dist1, u, l] = find_culprits_laplacian_multi_conve(infected, k, beta)

G = [0,1,0,1;1,0,1,0;0,1,0,0;1,0,0,0]
infected = [1,2,3]
k=5
beta = 0.1

new_infected = infected;
ids= zeros(k,1);
dist1= zeros(k,1);
curr_dist = 999999999; 
u = [];
l = [];
for i=1:k
		i
		[ids(i), u, l] = find_culprits_laplacian(G, new_infected, 1, beta);
		ids(i)
                exit(1)
		pause
		new_infected = setdiff(infected, ids(1:i));
		[dist1(i) path fsize temp] = mdl_full_code_test2(G, infected, ids(1:i), beta);
	%	if dist1(i) > curr_dist
				% stop no need to go further
				display 'here'
%				ids = ids(1:i-1);
	%			break;
%		end
		curr_dist = dist1(i);
end
