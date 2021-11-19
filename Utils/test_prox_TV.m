sz = [100,100,10];

y = randn(sz); alpha = 1e-1;
D = LinOpGrad(sz, 3, 'mirror');
L1 = CostL1(D.sizeout);
L2 = CostL2(sz, y);
CostTV = CostTV1D(sz, 3);
y_prox = CostTV.applyProx(y, alpha);
ADMM = OptiADMM(L2, {alpha*L1}, {D});
Outop = OutputOpti(true); ADMM.OutOp = Outop; ADMM.ItUpOut = 1;
ADMM.maxiter = 1e2;
ADMM.run(y_prox);
y_ADMM = ADMM.xopt;
cost_prox = L2.apply(y_prox) + alpha * CostTV.apply(y_prox);
cost_ADMM = L2.apply(y_ADMM) + alpha * CostTV.apply(y_ADMM);
cost_ADMM - cost_prox
figure; plot(ADMM.OutOp.evolcost);
