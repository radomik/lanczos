function init_data(work_dir, g, a_prefix, startvec_prefix)

% generate data
t = cputime;
G = numgrid('L', g); 
A = delsq(G); 
Aa = full(A);
n = length(A);
fprintf('init_data: n = %d, g = %d\n', n, g)

z = randn(n,1);
startvec = z/norm(z);

t1 = cputime-t;
fprintf('init_data: Input data generation completed. (t = %f)\n', t1)

% output data paths
out_A = strcat(work_dir,'/',a_prefix,'_n-',num2str(n),'_g-',num2str(g),'.csv');
out_startvec = strcat(work_dir,'/',startvec_prefix,'_n-',num2str(n),'.csv');

% save output data
fprintf('init_data: Writing A -> "%s"\n', out_A)
csvwrite(out_A, Aa);
fprintf('init_data: Writing startvec -> "%s"\n', out_startvec)
csvwrite(out_startvec, startvec);
