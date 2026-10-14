// Solid-phase volume-fraction equation

fvScalarMatrix psiEqn
(
    fvm::ddt(psi)
  + fvm::div(phipsi, psi)
  - fvm::laplacian(Dpsi, psi)
 ==
   -fvm::Sp(SpsiImp, psi)
   -SpsiExp
);

psiEqn.relax();
psiEqn.solve();


// ---------------------------------------------------------
// Diagnostic: check whether bounding modifies psi
// ---------------------------------------------------------

Info<< "psi BEFORE bounds: min = "
    << min(psi).value()
    << ", max = "
    << max(psi).value()
    << endl;

const scalar solidVolumeBeforeBounds =
    gSum(psi.primitiveField()*mesh.V());


// Physical bounds
psi.max(0.0);
psi.min(1.0);

psi.correctBoundaryConditions();


const scalar solidVolumeAfterBounds =
    gSum(psi.primitiveField()*mesh.V());

Info<< "psi AFTER bounds: min = "
    << min(psi).value()
    << ", max = "
    << max(psi).value()
    << endl;

Info<< "Solid-volume change from bounding = "
    << solidVolumeAfterBounds - solidVolumeBeforeBounds
    << endl;