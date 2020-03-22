/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBModel class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "flowCondition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam 
{
    defineTypeNameAndDebug(flowCondition, 0);
    defineRunTimeSelectionTable(flowCondition, dictionary);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::flowCondition::flowCondition
(
    dynamicFvMesh& mesh,
    const dictionary& dict 
)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::flowCondition> Foam::flowCondition::New
(
	dynamicFvMesh& mesh,
	const dictionary& dict
)
{
	word typeName = dict.lookup("type");
	
	Info<< "Manage flow condition: "<<typeName<< endl;
	
	dictionaryConstructorTable::iterator cstrIter = 
		dictionaryConstructorTablePtr_->find(typeName);
		
	if (cstrIter == dictionaryConstructorTablePtr_->end())
	{
		FatalErrorIn
		(
			"flowCondition::New(dynamicFvMesh& mesh, const dictionary& dict)"
		)	<< "Unknown approaching method "<< typeName <<endl <<endl
			<< "Valid approaching method are : " <<endl
			<< dictionaryConstructorTablePtr_->toc()
			<< exit(FatalError);
	}
	
	return autoPtr<flowCondition>(cstrIter()(mesh, dict));
}

