const volScalarField ceq
(
    "ceq",
    crif + cT*(T - Trif)
);

// epsH (smoothing width of the Heaviside switch) is read in createFields.H

const volScalarField H
(
    "H",
    0.5*(1.0 + tanh((c - ceq)/epsH))
);


// ---------------------------------------------------------
// Concentration reaction coefficients
//
// deposition:   requires fluid volume alphaFluid
// dissolution:  requires exposed solid, alphaDiss*psi
// ---------------------------------------------------------

const volScalarField ScImp
(
    "ScImp",
      H*alphaFluid/tau1
    + (1.0 - H)*alphaDiss*psi/tau2
);

const volScalarField ScExp
(
    "ScExp",
    -ScImp*ceq
);


// ---------------------------------------------------------
// Solid-fraction reaction coefficients
// ---------------------------------------------------------

const volScalarField SpsiImp
(
    "SpsiImp",
    -(1.0 - H)
    *(c - ceq)/tau2
    *psifact
    *alphaDiss
);

const volScalarField SpsiExp
(
    "SpsiExp",
    -H
    *(c - ceq)/tau1
    *psifact
    *alphaFluid
);