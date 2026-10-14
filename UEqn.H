/* Extra viscosity caused by deposited solid.
 The base molecular viscosity is already included in
 turbulence->divDevReff(U).*/
const volScalarField deltaNu  //nu_app-nu0
(
    "deltaNu",
    nuApp - turbulence->nu()
);
// Momentum predictor
tmp<fvVectorMatrix> tUEqn
(
   fvm::ddt(U)

   +
   turbulence->divDevReff(U)   //(nu0+nut)
   +
   fvm::div(phi, U)
   -
   fvm::laplacian(deltaNu, U)  
   -
   fvc::div(deltaNu*dev2(fvc::grad(U)().T()))  
   + fvm::Sp(betaBrinkman, U)
   ==
   fvOptions(U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}
