function [y] = solver_ADMM(A, Hn, zn, rho_n)

b = rho_n(1)*Hn{1}.applyAdjoint(zn{1});
for n = 2:length(Hn)
    b = b+rho_n(n)*Hn{n}.applyAdjoint(zn{n});
end
y = A.applyInverse(b);

end

