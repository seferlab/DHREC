function [dist1, dist2, dist3, ids1, orig, ids3] = analyze_one_run(output_arr, beta, GOIX)

%graph = 'oix.20061201.v3.1.20.100';
%output_arr = 'OUTPUT/oix-.si-seeds-final-state.0.2.3.out';
%output_arr = 'OUTPUT/oix-.si-seeds-final-state.0.01.3.out';
%output_arr = 'OUTPUT/oix-.si-seeds-final-state.0.01.2.out';
%beta = 0.2;
T = load(output_arr);
%T = [ 2     2     4     0     1     1     1     0     0     0     0];
%cpp_nodeMap = load(['OUTPUT/oix-node-mapping']); %'OUTPUT/' graph '.node_mapping']);
%cpp_nodeMap = load(node_mapping);
%cpp_nodeMap = load(['OUTPUT/' 'chain-' 'node-mapping']);
%n = max(cpp_nodeMap(:, 1));
%nM = zeros(n, 1);
%nM(cpp_nodeMap(:, 1)) = cpp_nodeMap(:, 2)+1;

%graph_file = ['../DATA/topologies/' graph '.topology'];
%graph_file = ['chain.txt'];
%[GOIX] = sparse(load_graph_topology_nodeMap(graph_file, nM));

infected_ids = find(T(T(1)+2:end));
infected_graph = GOIX(infected_ids, infected_ids);
[S, C] = graphconncomp(infected_graph);

dist1 = -99;
dist2 = -99;
dist3 = -99;
ids1 = [];
orig = [];
ids3 = [];
if S == 1
		[ids1, u1, l1] = find_culprits_laplacian_multi(GOIX, infected_ids, T(1), beta);
	%	[ids2, u1, l1] = find_culprits_laplacian_test2(GOIX, infected_ids, T(1), beta);
%		[ids2, u2, l2] = find_culprits_adjacency(GOIX, infected_ids, T(1), beta);
        % random guess now
        tosses = randsample(length(infected_ids), T(1));
        ids3 = infected_ids(tosses);
       % [ids3, u3, l3] = find_culprits_laplacian_local(GOIX, infected_ids, T(1), beta);
		%ids1'
		%ids2'
        %ids3'
		orig = T(2:2-1+T(1));
       
		%pause
        %[vdist1 IND1] = find_hop_distance(infected_graph, find(infected_ids == ids1), find(infected_ids == orig));
        %[vdist2 IND2] = find_hop_distance(infected_graph, find(infected_ids == ids2), find(infected_ids == orig));
        %[vdist3 IND3] = find_hop_distance(infected_graph, find(infected_ids == ids3), find(infected_ids == orig));
        [dist1 path fsize temp] = mdl_full_code_test2(GOIX, infected_ids, ids1, beta); %find_hop_distance(infected_graph, find(infected_ids == ids1), find(infected_ids == orig));
        [dist2 path fsize temp] = mdl_full_code_test2(GOIX, infected_ids, orig, beta); %find_hop_distance(infected_graph, find(infected_ids == ids1), find(infected_ids == orig));
        [dist3 path fsize temp] = mdl_full_code_test2(GOIX, infected_ids, ids3, beta); %find_hop_distance(infected_graph, find(infected_ids == ids1), find(infected_ids == orig));
  %      [vdist2 IND2] = find_hop_distance(infected_graph, find(infected_ids == ids2), find(infected_ids == orig));
  %      [vdist3 IND3] = find_hop_distance(infected_graph, find(infected_ids == ids3), find(infected_ids == orig));
        
        %dist1 = mean(vdist1);
  %      dist2 = mean(vdist2);
   %     dist3 = mean(vdist3);
		%u1full = zeros(length(GOIX), 1);
		%u1full(infected_ids) = u1;
		%ts = sum(u1);
		%sum(u1full(ids1)')*ts - sum(u1full(ids1))^2
		%sum(u1full(orig)')*ts - sum(u1full(orig))^2
		%{
		TD = sum(GOIX, 2);
		inds = setdiff(infected_ids, ids1);
		A = GOIX(inds, inds);
		TD = diag(TD(inds));
		LA = TD - A;
		size(LA)
		[u, l] = eigs(LA, 1, 'SM');
		[sum(abs(u)) abs(l)]
		pause
		orig = [1701 1663]
		inds = setdiff(infected_ids, orig);
		A = GOIX(inds, inds);
		LA = TD - A;
		size(LA)
		[u, l] = eigs(LA, 1, 'SM');
		[sum(abs(u)) abs(l)]
		%}
		


	%	u1full(orig)'

		%u2full = u2;
		%u2full(ids2)'
		%u2full(orig)'
        
      %  u3full = zeros(length(GOIX), 1);
	%	u3full(infected_ids) = u3;
	%	u3full(ids3)'
	%	u3full(orig)'
end
