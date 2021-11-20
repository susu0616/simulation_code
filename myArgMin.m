function x = myArgMin(Df, tilde_mu, H_i, d, rho, tilde_x, n, zeta, E_ini, E_ref, E_min, E_max)
cvx_begin quiet
    variable x(n)
    variable e(n)

    minimize( Df'*x+ tilde_mu'*(H_i*x-d) +BregDiv(x,tilde_x,1)/rho)
    subject to
        0.1*ones(n,1) <= x <= 1.5*ones(n,1);
        e(1) == E_ini;
        e(n) == E_ref;
        for i=1:n-1
            e(i+1) == e(i)+ zeta*x(i);
            E_min <= e(i+1) <= E_max;
        end
cvx_end