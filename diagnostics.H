// Per-timestep diagnostics of flowForm (solid/transition/fluid regions)

{
// ---------------------------------------------------------------------
// Velocity in partially deposited and fully solid regions
//
// Transition region: 0.5 < psi <= 0.99
// Fully solid region: psi > 0.99
// ---------------------------------------------------------------------

scalar maxUTransition = 0.0;
scalar sumUTransition = 0.0;
scalar volTransition  = 0.0;

scalar maxUSolid = 0.0;
scalar sumUSolid = 0.0;
scalar volSolid  = 0.0;

forAll(psi, celli)
{
    const scalar Vc   = mesh.V()[celli];
    const scalar Umag = mag(U[celli]);

    if (psi[celli] > 0.99)
    {
        maxUSolid = max(maxUSolid, Umag);
        sumUSolid += Umag*Vc;
        volSolid  += Vc;
    }
    else if (psi[celli] > 0.5)
    {
        maxUTransition = max(maxUTransition, Umag);
        sumUTransition += Umag*Vc;
        volTransition  += Vc;
    }
}

reduce(maxUTransition, maxOp<scalar>());
reduce(sumUTransition, sumOp<scalar>());
reduce(volTransition,  sumOp<scalar>());

reduce(maxUSolid, maxOp<scalar>());
reduce(sumUSolid, sumOp<scalar>());
reduce(volSolid,  sumOp<scalar>());


// ---------------------------------------------------------------------
// Concentration ranges in fluid and solid regions
//
// Fluid/transition: psi <= 0.5
// Deposited/solid:  psi > 0.5
// ---------------------------------------------------------------------

scalar maxCFluid = -GREAT;
scalar minCFluid =  GREAT;

scalar maxCSolid = -GREAT;
scalar minCSolid =  GREAT;

forAll(c, celli)
{
    if (psi[celli] > 0.5)
    {
        maxCSolid = max(maxCSolid, c[celli]);
        minCSolid = min(minCSolid, c[celli]);
    }
    else
    {
        maxCFluid = max(maxCFluid, c[celli]);
        minCFluid = min(minCFluid, c[celli]);
    }
}

reduce(maxCFluid, maxOp<scalar>());
reduce(minCFluid, minOp<scalar>());
reduce(maxCSolid, maxOp<scalar>());
reduce(minCSolid, minOp<scalar>());


// ---------------------------------------------------------------------
// Global solid volume
// ---------------------------------------------------------------------

const scalar totalSolidVolume =
    gSum(psi.primitiveField()*mesh.V());


// ---------------------------------------------------------------------
// Mean velocities
//
// Avoid division problems if one of the regions temporarily contains
// no cells.
// ---------------------------------------------------------------------

const scalar meanUTransition =
    sumUTransition/(volTransition + SMALL);

const scalar meanUSolid =
    sumUSolid/(volSolid + SMALL);


// ---------------------------------------------------------------------
// Print timestep diagnostics
// ---------------------------------------------------------------------

Info<< nl
    << "---------------- flowForm diagnostics ----------------" << nl

    << "psi:           min = " << min(psi).value()
    << ", max = " << max(psi).value() << nl

    << "betaBrinkman: min = " << min(betaBrinkman).value()
    << ", max = " << max(betaBrinkman).value() << nl

    << "alphaFluid:    min = " << min(alphaFluid).value()
    << ", max = " << max(alphaFluid).value() << nl

    << "alphaC:        min = " << min(alphaC).value()
    << ", max = " << max(alphaC).value() << nl

    << "alphaDiss:     min = " << min(alphaDiss).value()
    << ", max = " << max(alphaDiss).value() << nl

    << "c fluid:       min = " << minCFluid
    << ", max = " << maxCFluid << nl

    << "c solid:       min = " << minCSolid
    << ", max = " << maxCSolid << nl

    << "Transition |U| (0.5<psi<=0.99): max = "
    << maxUTransition
    << ", mean = " << meanUTransition << nl

    << "Solid |U| (psi>0.99):            max = "
    << maxUSolid
    << ", mean = " << meanUSolid << nl

    << "Transition volume: "
    << volTransition << nl

    << "Solid volume (psi>0.99): "
    << volSolid << nl

    << "Total SolidVolume: "
    << totalSolidVolume << nl

    << "T:             min = " << min(T).value()
    << ", max = " << max(T).value() << nl

    << "------------------------------------------------------"
    << nl << endl;
}
