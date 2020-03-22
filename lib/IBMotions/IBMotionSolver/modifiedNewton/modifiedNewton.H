/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - modifiedNewton class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
modifiedNewton class:
    - Solve newtons equations for 6DoF motion.
    - Integrate the motion of fluids inside objects.
	- Can handle objects with density close to fluids density
    - Kempe, T. J. Comput. Phys. 231, 3663–3684 (2012).
\*---------------------------------------------------------------------------*/
#ifndef modifiedNewton_H
#define modifiedNewton_H

// #include "IOmanip.H"
#include "autoPtr.H"
#include "dynamicFvMesh.H"
#include "eulerMesh.H"
#include "IBObject.H"
#include "IBModel.H"
#include "IBMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
							Class modifiedNewton declaration
\*---------------------------------------------------------------------------*/
class IBModel;

class modifiedNewton
:
	public IBMotionSolver
{
	dynamicFvMesh& mesh_;

	eulerMesh& emesh_;

	PtrList<IBObject>& ibo_;
	
	autoPtr<IBModel>& ibmodelPtr_;

public:
	
	TypeName("modifiedNewton");
	
	modifiedNewton
	(
		word typeName,
		dynamicFvMesh& mesh, 
		eulerMesh& emesh,
		PtrList<IBObject>& ibo,
		autoPtr<IBModel>& ibmodelPtr
	);

	virtual ~modifiedNewton()
	{}
	
	//- Move objects
	void moveObjects();
	
	vector volIntegralU
	(
	    IBObject& ibobj,
	    const volVectorField& U
	);
	
	vector volIntegralRxU
	(
	    IBObject& ibobj,
	    const volVectorField& U
	);

	scalar volFraction
	(
		IBObject& ibobj, 
		label cellID
	);
	
	scalar HeavisideFunc
	(
		scalar phi
	);
	
	scalar LevelSetFunc
	(
		scalar r, 
		point p, 
		point c
	);

}; //End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //