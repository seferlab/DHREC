function [ids, u, l] = find_culprits_laplacian(G, infected, k, beta, ...
                                               outfile)

size(G)

%G = beta * G;
F = [];
ninfec = length(infected);
%nnotinfec = length(F);

A = G(infected, infected);

DI = diag(sum(A));
factor_i = 1; %/sum(sum(A));
DI = factor_i * DI;

TD =  sum(G);
TD = diag(TD(infected));

DNI = TD - DI;
nnotinfec = 2*length(F); %sum(diag(DNI));
%factor = 1/nnotinfec;
factor = 1;
DNI = factor * DNI;
%{
SING = find(diag(TD) == 1);
for i = 1:length(SING)
	DNI(SNI(i), SNI(i)) = 1;
end
%}
%in = find(infected == 23)
%DNI(in, in) = 0;
%DI(in, in)

LA = DNI + DI - A;


% full laplacian
%L = diag(sum(G))/tA - G;

% take submatrix
%LA = L(infected, infected);


% find smallest eig 
[u, l] = eigs(LA, 1, 'SM');
u = abs(u);
l = abs(l);

% find top-k
% using O(kn) (w/o heap etc)
u_copy = u;
ind = zeros(k,1);
for i=1:k
    [val, ind(i)] = max(u_copy);
%	eq_ind = find(u_copy == val);
%	infected(eq_ind)
    u_copy(ind(i)) = -99;
	%DNI(ind(i), ind(i))
	%DI(ind(i), ind(i))
end
 
ids = infected(ind(1:k));

save(outfile, ids)
