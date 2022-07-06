/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    suspendedConcentrationFoam

Description
    Solves the steady or transient concentration equation of suspended sediment.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
		#include "createFields.H"
		#include "createTimeControls.H"
		#include "readTimeControls.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    Info<< "\nSolving concentration equation\n" << endl;
           
    while (runTime.loop())
    {	
		
				#include "readTimeControls.H"
		Info<< nl << "Time = " << runTime.timeName() << endl;
		
				
				
					volScalarField   EpsS = alpha*kappa*ustar*(y-(y*y)/h);
					volScalarField  dEpsSdy = alpha*kappa*ustar*(1.0 - (2.0*y)/h);
				    
				    
					volVectorField gradC = 	fvc::grad(C);
					volScalarField dcdx = gradC.component(vector::X);
					volScalarField dcdy = gradC.component(vector::Y);
					volVectorField laplacianC = fvc::grad(dcdy);
					volScalarField d2cdy2 = laplacianC.component(vector::Y);
					
					solve
					(
						fvm::ddt(C)
						+ 
						(alpha*kappa*ustar)*dcdx
						- 
						w0*dcdy		 
						==
					    dEpsSdy*dcdy
					    +
					    EpsS*d2cdy2
					  
					);
											
				runTime.write();
				
			}
				
			Info << "End\n" << endl;

			return 0;

    
}

// ************************************************************************* //
