
#include "fixedRe.H"
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
    defineTypeNameAndDebug(fixedRe, 0);
    addToRunTimeSelectionTable(flowCondition, fixedRe, dictionary);
}
// ------------------------ Private Member Funcstion ----------------------- //
void Foam::fixedRe::initiateDPdx()
{
    scalar avgU = Reynolds_ * nu_ / Dh_;

    dPdx_ = 12.0*avgU*nu_*1.0 / (Dh_*Dh_); 
    
    Info<<"Initialize dP/dx = "<<dPdx_<<" for Re = "<<Reynolds_<<nl<<endl;
}

// ------------------------------- Constructor ----------------------------- //
Foam::fixedRe::fixedRe
(
    dynamicFvMesh& mesh,
    const dictionary& dict
)
:
    flowCondition(mesh, dict),
    mesh_(mesh),
    Reynolds_(readScalar(dict.lookup("value"))),
    dPdx_(0),
    Dh_(readScalar(dict.lookup("Dh"))),
    tol_(readScalar(dict.lookup("tolerance"))),
    prevReynolds_(SMALL),
    flowDir_(dict.lookup("flowDirection")),
    outletName_(dict.lookup("outletName")),
    addPressureGradientField_(false)
{
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nutmp
    (
        "nu",
        dimensionSet(0,2,-1,0,0,0,0),
        transportProperties.lookup("nu")
    );

    nu_ = nutmp.value();

    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];
        if (isA<cyclicPolyPatch>(pp))
        {
            addPressureGradientField_ = true;
            break;
        }
    }

    if (dict.found("dPdx"))
    {
        dPdx_ = readScalar(dict.lookup("dPdx"));
    }
    else initiateDPdx();
}

// -----------------------------Member functions---------------------------- //

void Foam::fixedRe::correctDPdx()
{
    Info<<"Correcting dP/dx"<<endl;
    //- Caluclate flowrate at the outlet
    label outletPatchID = mesh_.boundaryMesh().findPatchID(outletName_);
    const polyPatch& outletpp = mesh_.boundaryMesh()[outletPatchID];
    const surfaceScalarField& phi = mesh_.lookupObject<surfaceScalarField>("phi");

    //- Mass flux at the outlet
    scalar massFluxCur = gSum(phi.boundaryField()[outletPatchID]);

    //- Calculate outlet Area
    scalar outletArea = gSum(mag(outletpp.faceAreas()));

    //- Calculate average velocity at outlet
    scalar avgU = mag(massFluxCur)/outletArea;

    //- Calculate current Reynolds number
    scalar curReynolds = avgU * Dh_ / nu_;
    Info<<"  Current Reynolds = "<<curReynolds<<","
        <<"  Previous Reynolds = "<<prevReynolds_<<",";

    scalar changeRate = mag(curReynolds - prevReynolds_)*100/prevReynolds_;   
    Info<<"  changeRate = "<<changeRate<<"%"<<endl;

    if (changeRate > 0.05)
    {
        Info<<"  => Remain dP/dx = "<<dPdx_<<nl<<endl;
    }
    else
    {
        //- Calculate the discrepancy of Reynolds in percent
        scalar deltaRe = mag(curReynolds - Reynolds_)*100/Reynolds_;
        
        Info<<"  Desired Reynolds = "<<Reynolds_<<","
            <<"  Current Reynolds = "<<curReynolds<<","
            <<"  DeltaRe = "<<deltaRe<<"%"<<nl<<endl;

        //- Adjust dPdx when the error more than 0.5%
        if (deltaRe > 0.01)
        {
            Info<<"  => Adjust dPdx from "<<dPdx_;
            dPdx_ = Reynolds_ / curReynolds * dPdx_;
            Info<<" to "<<dPdx_<<nl<<endl;
        }
    }
    prevReynolds_ = curReynolds;
}

Foam::volVectorField Foam::fixedRe::pressureField()
{
    //- Update dPdx
    if (mesh_.time().timeIndex() > 2) correctDPdx();

    volVectorField gradP
    (
        IOobject
        (
            "gradP",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "gradP", 
            dimensionSet(0,1,-2,0,0,0,0), 
            vector::zero
        )
    );

    if (addPressureGradientField_)
    {
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
    }

  	return gradP;
}
