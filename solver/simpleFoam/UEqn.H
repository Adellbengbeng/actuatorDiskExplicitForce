    // Momentum predictor

    MRF.correctBoundaryVelocity(U);

	//Added by me
    //Clear  old  values  in  ioVolumeForce
    for(label k = 0; k < VolumeForce.size(); k++) {
        VolumeForce[k] = vector::zero;
    }

    // Calculate  volume  force  from  actuator  disk
    actuatorDisk.CalcActuatorDiskVolForce(mesh, VolumeForce);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U) + VolumeForce
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }


