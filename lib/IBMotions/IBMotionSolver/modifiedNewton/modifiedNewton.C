/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - modifiedNewton class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "modifiedNewton.H"
#include "wallPolyPatch.H"
#include "IBParticle.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(modifiedNewton, 0);
    addToRunTimeSelectionTable(IBMotionSolver, modifiedNewton, dictionary);
}

//---------------------------------Constructors------------------------------//
Foam::modifiedNewton::modifiedNewton
(
    word typeName,
	dynamicFvMesh& mesh,
    eulerMesh& emesh,
    PtrList<IBObject>& ibo,
    autoPtr<IBModel>& ibmodelPtr
)
:
    IBMotionSolver(typeName, mesh, emesh, ibo, ibmodelPtr),
	mesh_(mesh),
    emesh_(emesh),
    ibo_(ibo),
    ibmodelPtr_(ibmodelPtr)
{}

// -------------------------------Member Functions----------------------------//

void  Foam::modifiedNewton::moveObjects()
{
    for(int i=0; i<ibo_.size(); i++)
    {
        if (ibo_[i].movable())
        {
            Info<< "Moving object: "<<ibo_[i].name()<<endl;
            
            vector repulsiveForce(vector::zero);
            
            for(int j=0; j<ibo_.size(); j++)
            {
                if (i == j)
                    continue;
                else
                {
                    repulsiveForce += objMutualRepulsive(i, j);
                }
            }
            
            const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
            const volVectorField& Uold = U.oldTime();
            const scalar dT = mesh_.time().deltaTValue(); 
            vector volIntegralU_ 
            = 
                (
                    volIntegralU(ibo_[i], U) 
                  - volIntegralU(ibo_[i], Uold)
                ) / dT;

            vector volIntegralRxU_ 
            = 
                (
                    volIntegralRxU(ibo_[i], U) 
                  - volIntegralRxU(ibo_[i], Uold)
                ) / dT;
        
            ibo_[i].updateObjectMotionTobias
            (
                ibmodelPtr_->FLagrange()[i],
                repulsiveForce,
                emesh_.rhoFluid(),
                emesh_.g(),
                volIntegralU_,
                volIntegralRxU_
            );

            ibo_[i].updateVelocity();
            ibo_[i].movePoints();
        }
    }
}

Foam::vector Foam::modifiedNewton::volIntegralU
(
    IBObject& ibobj,
    const volVectorField& U
)
{   
    vector IU(vector::zero);

    forAll(ibobj.solidCellsExt(), cellI)
    {
        scalar vf = volFraction(ibobj, ibobj.solidCellsExt()[cellI]);

        IU += U[ibobj.solidCellsExt()[cellI]] * emesh_.dV() * vf;
    }

    return IU;
}

Foam::vector Foam::modifiedNewton::volIntegralRxU
(
    IBObject& ibobj,
    const volVectorField& U
)
{
    const volVectorField& cc = mesh_.C();

    vector IRxU(vector::zero);

    forAll(ibobj.solidCellsExt(), cellI)
    {
        scalar vf = volFraction(ibobj, ibobj.solidCellsExt()[cellI]);

        vector r = cc[ibobj.solidCellsExt()[cellI]] - ibobj.CG();

        IRxU +=
        ( 
            (r ^ U[ibobj.solidCellsExt()[cellI]]) * emesh_.dV() * vf
        );
    } 

    return IRxU;
}

Foam::scalar Foam::modifiedNewton::volFraction(IBObject& ibobj, label cellID)
{
    const faceList& ff = mesh_.faces();
    const pointField& pp = mesh_.points();
    const cell& cc = mesh_.cells()[cellID];
    pointField cellVertices = cc.points(ff, pp);

    IBParticle& ibp = refCast<IBParticle>(ibobj);
    
    scalar alphaIJK(0.0);
    scalar sumPhi(0.0);
    
    forAll(cellVertices, pointI)
    {
        scalar phi_m = LevelSetFunc(ibp.R(), cellVertices[pointI], ibp.CG());
        alphaIJK += -phi_m * HeavisideFunc(-phi_m);
        sumPhi += mag(phi_m);
    }

    return (alphaIJK/sumPhi);
}

Foam::scalar Foam::modifiedNewton::HeavisideFunc(scalar phi)
{
    if (phi<=0) 
	{
		return 0;
	}
	else
		return 1;
}

Foam::scalar Foam::modifiedNewton::LevelSetFunc
(
	scalar r, 
	point p, 
	point c
)
{
	if(mesh_.nGeometricD() == 2)
	{
	    return 
	    (
	    	sqrt
			(
				pow(p.x()-c.x(),2) 
			  + pow(p.y()-c.y(),2)
			)/r - 1.0
	    );
	}
	else 
	{
		return ( mag(p-c)/r - 1.0);
	}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //