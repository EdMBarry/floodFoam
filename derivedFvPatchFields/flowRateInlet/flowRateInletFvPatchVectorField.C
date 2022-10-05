/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Author
    KM-Turbulenz GmbH, 2009
\*---------------------------------------------------------------------------*/

#include "flowRateInletFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
flowRateInletFvPatchVectorField::flowRateInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    flowRate_(Zero),
    HName_("H")
{}


Foam::flowRateInletFvPatchVectorField::flowRateInletFvPatchVectorField
(
    const flowRateInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    flowRate_(ptf.flowRate_),
    HName_(ptf.HName_)
{
}


Foam::flowRateInletFvPatchVectorField::flowRateInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    flowRate_(readScalar(dict.lookup("flowRate"))),
    HName_("H")
{

    if (dict.found("H"))
    {
        dict.lookup("H") >> HName_;
    }
}


Foam::flowRateInletFvPatchVectorField::flowRateInletFvPatchVectorField
(
    const flowRateInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    flowRate_(ptf.flowRate_),
    HName_(ptf.HName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::flowRateInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    vectorField n = patch().nf();

    const volScalarField& H = db().lookupObject<volScalarField>(HName_);

    scalar surfaceSumH = gSum(H.boundaryField()[patch().index()]*patch().magSf());
    scalar Umean = flowRate_/surfaceSumH;
    scalarField phip = patch().lookupPatchField<surfaceScalarField, scalar>("phi");
    scalar surfaceSumPhip = gSum(phip*H.boundaryField()[patch().index()]);
    Info << "Inflow rate - Q: " <<  surfaceSumPhip << endl;   

    operator==(-n*Umean*H.boundaryField()[patch().index()]);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::flowRateInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);

    os.writeKeyword("flowRate") << flowRate_
        << token::END_STATEMENT << nl;

    if (HName_ != "H")
    {
        os.writeKeyword("H") << HName_ << token::END_STATEMENT << nl;
    }

    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       flowRateInletFvPatchVectorField
   );
}


// ************************************************************************* //
