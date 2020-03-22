
#include "fixedDP.H"
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
    defineTypeNameAndDebug(fixedDP, 0);
    addToRunTimeSelectionTable(flowCondition, fixedDP, dictionary);
}

// ------------------------------- Constructor ----------------------------- //
Foam::fixedDP::fixedDP
(
    dynamicFvMesh& mesh,
    const dictionary& dict
)
:
    flowCondition(mesh, dict),
    mesh_(mesh),
    dPdx_(readScalar(dict.lookup("dPdx"))),
    flowDir_(dict.lookup("flowDirection"))
{}

// -----------------------------Member functions---------------------------- //
void Foam::fixedDP::correctDPdx()
{
	// do nothing
}

Foam::volVectorField Foam::fixedDP::pressureField()
{
    volVectorField gradP
    (
        IOobject
        (
            "gradP",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "gradP", 
            dimensionSet(0,1,-2,0,0,0,0), 
            vector::zero
        )
    );

    Info<<"Adding field dP/dx = "<<dPdx_<<nl<<endl;

    forAll(gradP, celli)
    {
        if (flowDir_ == "horizontal")
            gradP[celli] = vector(dPdx_, 0, 0);
        else if (flowDir_ == "vertical")
            gradP[celli] = vector(0, dPdx_, 0);
        else if (flowDir_ == "curvedChannel")
        {
            //- Assuming center of curve channel is (0 0 0)
            vector position = mesh_.C()[celli];
            vector flowDir = vector(-position.y(), position.x(), 0);
            flowDir /= mag(flowDir);
            gradP[celli] = dPdx_ * flowDir;
        }
        else
        FatalErrorIn
        (
            "IBTechnique::pressureGradField()"
        )	<< "Unknown direction of external pressure gradient field "
            << flowDir_ <<endl
            << "Valid approaching method are : " <<endl
            << "(" <<endl
            << "horizontal" <<endl 
            << "vertical" <<endl
            << "curvedChannel"<<endl
            <<")"<<endl
            << exit(FatalError);
    }

  	return gradP;
}
