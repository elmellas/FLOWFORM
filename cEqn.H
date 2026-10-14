fvScalarMatrix cEqn
(
    fvm::ddt(alphaC, c)
  + fvm::div(phic, c)
  - fvm::laplacian(Dc*alphaFace, c)
 ==
   -fvm::Sp(ScImp, c)
   -ScExp
);

cEqn.relax();
cEqn.solve();

c.max(0.0);
c.correctBoundaryConditions();