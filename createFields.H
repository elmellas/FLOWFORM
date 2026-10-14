
    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // concentration diffusivity
    dimensionedScalar Dc
    (
        transportProperties.lookup("Dc")
    );

    // wax diffusivity
    dimensionedScalar Dpsi
    (
        transportProperties.lookup("Dpsi")
    );

    dimensionedScalar rho0
    (
        transportProperties.lookup("rho0")
    );
    // viscosity-law regularisation parameter
    dimensionedScalar psim
    (
        transportProperties.lookup("psim")
    );
    // deposition relaxation time
     dimensionedScalar tau1
    (
        transportProperties.lookup("tau1")
    );
    // redissolution relaxation time
     dimensionedScalar tau2
    (
        transportProperties.lookup("tau2")
    );
    // particle relaxation time
    dimensionedScalar tauP
    (
        transportProperties.lookup("tauP")
    );
    // reference temperature
    dimensionedScalar Trif
    (
        transportProperties.lookup("Trif")
    );

    // reference concentration
    dimensionedScalar crif
    (
        transportProperties.lookup("crif")
    );
    // concentration per unit temperature
    dimensionedScalar cT
    (
        transportProperties.lookup("cT")
    );
    //Molecular mass
    dimensionedScalar MM
    (
        transportProperties.lookup("MM")
    );
    //constant for solve density
    dimensionedScalar as
    (
        transportProperties.lookup("as")
    );

    //specific heat
    dimensionedScalar cp
    (
        transportProperties.lookup("cp")
    );
    //thermal conductivity
    dimensionedScalar kTermal
    (
        transportProperties.lookup("kTermal")
    );
    // turbulent Prandtl number
    dimensionedScalar Prt
    (
        transportProperties.lookup("Prt")
    );

    // Maximum Brinkman penalisation coefficient
    dimensionedScalar betaMax
    (
        transportProperties.lookup("betaMax")
    );
    
    // Brinkman interpolation exponent
    dimensionedScalar brinkmanExponent
    (
        transportProperties.lookup("brinkmanExponent")
    );

    dimensionedScalar alphaMin
    (
        transportProperties.lookup("alphaMin")
    );

    // Smoothing width of the precipitation/dissolution switch
    // H = 0.5*(1 + tanh((c - ceq)/epsH)).  Default = value of the Fig. 2 runs.
    const dimensionedScalar epsH
    (
        dimensionedScalar::getOrDefault
        (
            "epsH",
            transportProperties,
            dimensionSet(0,-3,0,0,0,0,0),
            0.05
        )
    );
    

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );



    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


   
    


    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field c\n" << endl;
    volScalarField c
    (
        IOobject
        (
            "c",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    Info<< "Reading field psi\n" << endl;
    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    volScalarField alphaFluid
    (
    IOobject
    (
        "alphaFluid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    scalar(1) - psi
    );
    volScalarField alphaC
    (
        IOobject
        (
            "alphaC",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        max(alphaFluid, alphaMin)
    );
    volScalarField alphaDiss
    (
        IOobject
        (
            "alphaDiss",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphaFluid
    );
    
    alphaC.oldTime();
    alphaC.oldTime().oldTime();
    volScalarField betaBrinkman
    (
        IOobject
        (
            "betaBrinkman",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        betaMax*pow(psi, brinkmanExponent.value())
    );

// Required by fvm::ddt(alphaFluid,c)
alphaFluid.oldTime();
alphaFluid.oldTime().oldTime();
#include "createPhi.H"

volVectorField Uslip
(
    IOobject
    (
        "Uslip",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "zero",
        dimVelocity,
        vector::zero
    )
);

    //volVectorField Uslip("Uslip",U*0);
     surfaceScalarField phislip
     (
         "phislip",
         phi*0
     );


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);
    mesh.setFluxRequired(p.name());

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );



    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho0-as*(T - Trif)
    );

    volScalarField nuApp
(
    IOobject
    (
        "nuApp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    turbulence->nu()
);
//

    #include "createFvOptions.H" //added

    

/////////////////////////////////////////////////////////////////
