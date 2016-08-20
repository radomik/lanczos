function lancno_init(work_dir, m, a_name, startvec_name, a_vec_prefix, b_vec_prefix, anorm_prefix)

% input data paths
in_csv_A = strcat(work_dir, '/', a_name, '.csv');
in_csv_startvec = strcat(work_dir, '/', startvec_name, '.csv');

% read input data
fprintf('lancno_init: Reading A <- "%s"\n', in_csv_A)
A = csvread(in_csv_A);
fprintf('lancno_init: Reading startvec <- "%s"\n', in_csv_startvec)
startvec = csvread(in_csv_startvec);

% data size
n = length(A);
fprintf('lancno_init: n = %d, m = %d\n', n, m)

% output data paths
out_a = strcat(work_dir,'/',a_vec_prefix,'_n-',num2str(n),'_m-',num2str(m),'.csv');
out_b = strcat(work_dir,'/',b_vec_prefix,'_n-',num2str(n),'_m-',num2str(m),'.csv');
out_an = strcat(work_dir,'/',anorm_prefix,'_n-',num2str(n),'_m-',num2str(m),'.csv');

% initial Lanczos iteration, m steps
t = cputime;
v = startvec;
a = zeros(m,1); b = zeros(m-1,1);

for k = 1:m
    if k == 1
        r = A*v;
    else
        r = A*v - b(k-1)*v2;
    end
    
    a(k) = v'*r;
    r = r - a(k)*v;
    b(k) = norm(r);
    v2 = v;
    v = 1/b(k)*r;
    
    % estimate |A|_2 by |T|_1
    if k == 1
        anorm = abs( a(1)+b(1) );
    else
        anorm = max( anorm, b(k-1)+abs(a(k))+b(k) );
    end;
end

t1 = cputime-t;
fprintf('lancno_init: Initial Lanczos iteration completed. (t = %f)\n', t1)

% save output data
fprintf('lancno_init: Writing a -> "%s"\n', out_a)
csvwrite(out_a, a);
fprintf('lancno_init: Writing b -> "%s"\n', out_b)
csvwrite(out_b, b);
fprintf('lancno_init: Writing anorm -> "%s"\n', out_an)
csvwrite(out_an, anorm);
