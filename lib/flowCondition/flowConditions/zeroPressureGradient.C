
#include "zeroPressureGradient.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
// #include "cyclicPolyPatch.H"
// #include "GeometricField.H"
// #include "volMesh.H"
// // #include "IOdictionary.H"
// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zeroPressureGradient, 0);
    addToRunTimeSelectionTable(flowCondition, zeroPressureGradient, dictionary);
}
// ------------------------ Private Member Funcstion ----------------------- //
void Foam::zeroPressureGradient::initiateDPdx()
{
    
}

// ------------------------------- Constructor ----------------------------- //
Foam::zeroPressureGradient::zeroPressureGradient
(
    eulerMesh& m,
    const dictionary& dict
)
:
    flowCondition(m, dict), // bat buoc phai mang ten cua class ke thua cua no s
    emesh_(m)
{}

// -----------------------------Member functions---------------------------- //

void Foam::zeroPressureGradient::correctDPdx()
{
   // Info<< "day la ham cua duc anh dz"<<endl;
}

Foam::volVectorField Foam::zeroPressureGradient::pressureField()
{
    //- Update dPdx
    volVectorField gradP
    (
        IOobject
        (
            "gradP",
            emesh_.mesh().time().timeName(),
            emesh_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        emesh_.mesh(),
        dimensionedVector
        (
            "gradP", 
            dimensionSet(0,1,-2,0,0,0,0), 
            vector::zero
        )
    );
    // Info<< "day la ham cua duc anh dz"<<endl;
  	return gradP;
    
}
