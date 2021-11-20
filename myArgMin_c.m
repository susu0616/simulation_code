function x = myArgMin_c(a_c, alpha, x_t,n)
cvx_begin quiet
    variable x(n)
    l = zeros(n,1);
    u = 5*ones(n,1);
    minimize( alpha*x'*a_c + alpha* (norm(x,1)+5*x'*x) +BregDiv(x,x_t,1))
    subject to
        l <= x <= u
cvx_end