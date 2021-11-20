function y = myArgMax(tilde_mu, g, delta, rho, m)
cvx_begin quiet
    variable y(m)
    maximize( y'*(g+delta) -BregDiv(y,tilde_mu,2)/rho)
    subject to
        y>=0
cvx_end