/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
    for more details.

Application

    flowForm

    Ismail El Mellas

Description

    Transient incompressible solver for transport and deposition of a
    dispersed/solid phase coupled to evolving rheology and Brinkman
    penalisation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "turbulentTransportModel.H"
#include "singlePhaseTransportModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createTimeControls.H"
    #include "createMesh.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    pimpleControl pimple(mesh);

    turbulence->validate();


    Info<< "\nStarting time loop\n" << endl;


    // =====================================================================
    // Physical time loop
    // =====================================================================

    while (runTime.run())
    {
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;


        // =================================================================
        // PIMPLE outer-corrector loop
        // =================================================================

        while (pimple.loop())
        {
            // -------------------------------------------------------------
            // Phase fractions derived from current psi
            // -------------------------------------------------------------

            alphaFluid = scalar(1) - psi;

            // Numerical regularisation for concentration storage
            alphaC = max(alphaFluid, alphaMin);


            // -------------------------------------------------------------
            // Fluid accessibility for dissolution
            //
            // alphaDiss > 0 in solid cells neighbouring fluid.
            // Buried solid remains chemically inaccessible.
            // -------------------------------------------------------------

            alphaDiss = fvc::average(fvc::interpolate(alphaFluid));
            alphaDiss = max(alphaDiss, alphaFluid);

            alphaDiss.max(0.0);
            alphaDiss.min(1.0);


            // -------------------------------------------------------------
            // Apparent viscosity (Maron-Pierce type, pole at psi = psim;
            // psim > 1 keeps it finite for psi <= 1)
            // -------------------------------------------------------------

            nuApp = turbulence->nu()/sqr(1.0 - psi/psim);


            // -------------------------------------------------------------
            // Brinkman resistance
            //
            // beta = betaMax * psi^n
            // -------------------------------------------------------------

            betaBrinkman =
                betaMax*pow(psi, brinkmanExponent.value());

            
            // -------------------------------------------------------------
            // Temperature-dependent material density
            // -------------------------------------------------------------

            rho = rho0 - as*(T - Trif);


            // =============================================================
            // Momentum equation
            // =============================================================

            Info<< "Solving for fluid velocity" << endl;

            #include "UEqn.H"


            // =============================================================
            // Pressure-velocity correction
            // =============================================================

            while (pimple.correct())
            {
                Info<< "Solving for mixture pressure" << endl;

                #include "pEqn.H"

                Info<< "Solving for dispersed phase slip" << endl;

                #include "updatePhiSlip.H"
            }


            // =============================================================
            // Advective fluxes for dispersed and dissolved phases
            // =============================================================

            const surfaceScalarField phipsi
            (
                "phi",
                phi + phislip
            );


            // Harmonic interpolation of the regularised fluid fraction.
            //
            // This prevents large dissolved-species fluxes from entering
            // cells with alphaFluid = 0 and only alphaMin numerical storage.
            const surfaceScalarField alphaFace
            (
                "alphaFace",
                1.0/fvc::interpolate(1.0/alphaC)
            );


            const surfaceScalarField phic
            (
                "phi",
                phi*alphaFace
            );


            // =============================================================
            // Thermal diffusivity
            // =============================================================

            const volScalarField alphaT
            (
                "alphaT",
                kTermal/(rho*cp)
            );


            const volScalarField alphaEff
            (
                "alphaEff",
                alphaT + turbulence->nut()/Prt
            );


            // =============================================================
            // Temperature equation
            // =============================================================

            Info<< "Solving for temperature" << endl;

            #include "TEqn.H"


            // Update rho using the newly solved temperature
            rho = rho0 - as*(T - Trif);


            // -------------------------------------------------------------
            // Concentration-to-solid-volume conversion
            // -------------------------------------------------------------

            const volScalarField psifact
            (
                "psifact",
                MM/rho
            );


            // =============================================================
            // Reaction source terms
            // =============================================================

            Info<< "Computing source terms" << endl;

            #include "SEqn.H"

            Info<< "c - ceq: min = "
            << min(c - ceq).value()
            << ", max = "
            << max(c - ceq).value()
            << endl;


            // =============================================================
            // Solid/dispersed phase equation
            // =============================================================

            Info<< "Solving for dispersed phase" << endl;

            #include "psiEqn.H"


            // =============================================================
            // Dissolved concentration equation
            //
            // IMPORTANT:
            // alphaFluid, alphaC, alphaDiss and alphaFace are intentionally
            // NOT updated between psiEqn and cEqn. Both equations therefore
            // use the same lagged phase state and reaction coefficients.
            // =============================================================

            Info<< "Solving for concentration" << endl;

            #include "cEqn.H"


            // =============================================================
            // Refresh derived fields after updated psi
            //
            // These values are used by the next PIMPLE iteration and are
            // consistent with the fields written at the end of the timestep.
            // =============================================================

            alphaFluid = scalar(1) - psi;

            alphaC = max(alphaFluid, alphaMin);

            alphaDiss = fvc::average(fvc::interpolate(alphaFluid));
            alphaDiss = max(alphaDiss, alphaFluid);

            alphaDiss.max(0.0);
            alphaDiss.min(1.0);

            betaBrinkman =
                betaMax*pow(psi, brinkmanExponent.value());


            // -------------------------------------------------------------
            // Update transport/turbulence models
            // -------------------------------------------------------------

            laminarTransport.correct();
            turbulence->correct();
        }


        // Diagnostics: one report per physical timestep
        #include "diagnostics.H"


        // =================================================================
        // Write fields
        // =================================================================

        runTime.write();


        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = "
            << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //