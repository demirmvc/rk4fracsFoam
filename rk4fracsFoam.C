/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Runge Kutta 4th order
-------------------------------------------------------------------------------
LORA
Implementation of fractional step algorithm with Runge-Kutta method in OpenFOAM for Large Eddy Simulations in hydraulic engineering

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow,"
        " using the PISO algorithm."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    #include "pEqn.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        Uold=U; Uc=U;
        phi = (fvc::interpolate(U) & mesh.Sf());
        dU = runTime.deltaT()*(fvc::laplacian(turbulence->nuEff(),U) - fvc::div(phi,U) );
        Uc = Uc + (1.0/6.0)*dU; U = Uold + 0.5*dU;
        #include "pCorrection.H"
        phi = (fvc::interpolate(U) & mesh.Sf());
        dU = runTime.deltaT()*(fvc::laplacian(turbulence->nuEff(),U) - fvc::div(phi,U) );
        Uc = Uc + (1.0/3.0)*dU; U = Uold + 0.5*dU;
        #include "pCorrection.H"
        phi = (fvc::interpolate(U) & mesh.Sf());
        dU = runTime.deltaT()*(fvc::laplacian(turbulence->nuEff(),U) - fvc::div(phi,U) );
        Uc = Uc + (1.0/3.0)*dU; U = Uold + dU;
        #include "pCorrection.H"
        phi = (fvc::interpolate(U) & mesh.Sf());
        dU = runTime.deltaT()*(fvc::laplacian(turbulence->nuEff(),U) - fvc::div(phi,U) );
        Uc = Uc + (1.0/6.0)*dU; U = Uc;
        #include "pCorrection.H"
        phi = (fvc::interpolate(U) & mesh.Sf());
        turbulence->correct();

        laminarTransport.correct();
        turbulence->correct();

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //