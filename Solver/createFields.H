    /*------------------------------------------------------------------------------*/
    Info<< "Reading field C\n" << endl;

    volScalarField C
    (
        IOobject
        (
            "C",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );    
	/*------------------------------------------------------------------------------*/   
	Info<< "Reading vertical depth field y\n" << endl;

	volScalarField y
	(
		IOobject
		(
			"y",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh.C().component(vector::Y)
	); 	
	/*------------------------------------------------------------------------------*/       
    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    /*------------------------------------------------------------------------------*/ 
    Info<< "Reading settling Velocity w0\n" << endl;
    dimensionedScalar w0
    (
        transportProperties.lookup("w0")
    );
    
    
    Info<< "Reading alpha\n" << endl;
    dimensionedScalar alpha
    (
        transportProperties.lookup("alpha")
    );
    
    Info<< "Reading von-karman constant kappa\n" << endl;
    dimensionedScalar kappa
    (
        transportProperties.lookup("kappa")
    );
    
    Info<< "Reading shear velocity ustar\n" << endl;
    dimensionedScalar ustar
    (
        transportProperties.lookup("ustar")
    );
    
    Info<< "Reading maximum flow depth h\n" << endl;
    dimensionedScalar h
    (
        transportProperties.lookup("h")
    );
	/*------------------------------------------------------------------------------*/

