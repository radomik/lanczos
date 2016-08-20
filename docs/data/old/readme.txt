% A_n-<N>_g-<G>.csv
% startvec_n-<N>.csv
%	g = <G>
%	n = <N>
	G = numgrid('L', g); 
	A = delsq(G); 
	Aa = full(A);
	n = length(A);
	csvwrite('A_n-<N>_g-<G>.csv', Aa);

	z = randn(n,1);
	startvec = z/norm(z);
	csvwrite('startvec_n-<N>.csv', startvec);
