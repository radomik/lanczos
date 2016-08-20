% Lanczos with no reorthogonalization
% Input - already defined variables:
%   A        - symmetric, square matrix
%   startvec - start vector of size (length(A),1) and unit length
%   m        - count of lanczos iterations

n = length(A);
v = startvec;

% initial Lanczos iteration, m steps
t = cputime;
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
fprintf('Initial Lanczos iteration completed. (t = %2.2f)\n', t1)

% tridiagonal matrix
[ritz,S] = mrrr(a,b); % calls the LAPACK routine DSTEGR
lres = abs(b(m)*S(m,:))'; % residual estimation
cul = abs(S(1,:))';

% remove non-converged and spurious Ritz values
tol = eps*anorm;
idx = (lres < tol) & (cul > tol);
ritz = ritz(idx);
lres = lres(idx);
cul = cul(idx);
S = S(:,idx);
mm = sum(idx);

% distinguish the clusters of eigenvalues
tol = 1e-8;
i = 1; ci = {};
for k = 2:mm
    if abs(ritz(k)-ritz(i)) > tol
        ci{end+1} = (i:k-1)';
        i = k; k = k + 1;
    end
end
ci{end+1} = (i:mm)';

% reflect using Householder
e = []; c = []; S2 = zeros(m,length(ci));
for i = 1:length(ci)
    cii = ci{i};
    [y,j] = max(cul(cii));
    ji = cii(1)+j-1;
    
    x = S(1,cii)'; u = x;
    u(j) = u(j) + sign(x(j))*norm(x);
    s = S(:,ji)-2/(u'*u)*(S(:,cii)*u)*u(j);
    S2(:,i) = s;
    
    idx(cii) = false; idx(ji) = true;
    e = [e; ritz(ji)];
    c = [c; lres(ji)];
end

t2 = cputime-t-t1;
fprintf('Eigenvalue computation complete. (t = %2.2f)\n',t2)

% compute eigenvectors
v = startvec;
X = zeros(n,size(S2,2));
for k = 1:m-1
    X = X + v*S2(k,:);
    if k == 1
    	r = A*v;
    else
        r = A*v - b2*v2;
    end
    a = v'*r;
    r = r - a*v;
    b = norm(r);
    v2 = v; a2 = a; b2 = b;
    v = 1/b*r;
end
X = X + v*S2(m,:);
X = X*diag(1./normc(X));

t3 = cputime-t-t1-t2;
fprintf('Eigenvector computation complete. (t = %2.2f)\n', t3)
t = cputime-t;

% output
fprintf(['%d/%d eigenvalues converged.\n' ...
'||AX-XD|| = %e\n||X’’X-I|| = %e\n' ...
'CPU time used: %2.2f\n'], length(e), n, ...
norm(A*X-X*diag(e))/anorm, norm(X'*X-eye(size(X,2))), t)
