/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBObject class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/
#include "error.H"
// #include "IBMotion.H"
#include "IBObject.H"
#include "GeometricField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IBObject, 0);
    defineRunTimeSelectionTable(IBObject, dictionary);
}
//---------------------------------Constructors------------------------------//
Foam::IBObject::IBObject
(
	const word& typeName,
	dynamicFvMesh& mesh,
	eulerMesh& emesh,
	const dictionary& dict
)
{}

//------------------------------------Selectors------------------------------//

Foam::autoPtr<Foam::IBObject> Foam::IBObject::New
(
	const word& typeName,
	dynamicFvMesh& mesh,
	eulerMesh& emesh,
	const dictionary& dict
)
{
	dictionaryConstructorTable::iterator dictCstrIter = 
		dictionaryConstructorTablePtr_->find(typeName);

	if(dictCstrIter == dictionaryConstructorTablePtr_->end())
	{
		FatalErrorIn
		(
			"IBObject::New(const dynamicFvMesh& mesh, "
			"const word& name, const dictionary& dict)"
		)	<< "Unknown object type "<< typeName 
			<< endl
			<< "Valid object type are: "<<endl
			<<dictionaryConstructorTablePtr_->toc()
			<<exit(FatalError);
	}

	return autoPtr<IBObject>(dictCstrIter()(typeName, mesh, emesh, dict));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
