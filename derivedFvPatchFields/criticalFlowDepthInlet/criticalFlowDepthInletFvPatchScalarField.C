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

\*---------------------------------------------------------------------------*/

#include "criticalFlowDepthInletFvPatchScalarField.H"
#include "freestreamFvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

criticalFlowDepthInletFvPatchScalarField::criticalFlowDepthInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    criticalH_(Zero)
{}


criticalFlowDepthInletFvPatchScalarField::criticalFlowDepthInletFvPatchScalarField
(
    const criticalFlowDepthInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    criticalH_(ptf.criticalH_)    
{}


criticalFlowDepthInletFvPatchScalarField::criticalFlowDepthInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    criticalH_(readScalar(dict.lookup("criticalH")))
{}


criticalFlowDepthInletFvPatchScalarField::criticalFlowDepthInletFvPatchScalarField
(
    const criticalFlowDepthInletFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wbppsf, iF),
    criticalH_(wbppsf.criticalH_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void criticalFlowDepthInletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
 
//  give internal values of H next to boundary patch
    scalarField HInt = patch().lookupPatchField<volScalarField, scalar>("H").patchInternalField();
    scalarField Sp = patch().lookupPatchField<volScalarField, scalar>("S");
    
//    operator==((Hint-criticalH_)*pos(Hint-criticalH_) + criticalH_);
    operator==((HInt - (criticalH_ - Sp))*pos(HInt-(criticalH_-Sp)) + (criticalH_-Sp));

    fixedValueFvPatchScalarField::updateCoeffs();
}

void criticalFlowDepthInletFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    os.writeKeyword("criticalH") << criticalH_
        << token::END_STATEMENT << nl;

    writeEntry(os, "value", *this);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, criticalFlowDepthInletFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
