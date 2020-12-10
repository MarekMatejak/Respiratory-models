within ;
model LungsCompartment
   import Modelica.SIunits.*;
   parameter Boolean useSigmoidCompliance = true "sigmoid compliance e.g. lungs"
   annotation(Evaluate=true, HideResult=true, choices(checkBox=true),Dialog(group="Computational model"));

   parameter Volume ResidualVolume = 0.00123  "Allways presented residual volume if useSigmoidCompliance";
   parameter Volume FunctionalResidualCapacity = 0.00231 "Zero pressure volume for linear compliance model";
   parameter Volume VitalCapacity = 0.00493  "Relative volume capacity if useSigmoidCompliance";
   parameter Volume BaseTidalVolume = 0.000543 "Base value of tidal volume";
   parameter Physiolib.Types.HydraulicCompliance Compliance = 0.000543/1000 "Compliance e.g. TidalVolume/TidalPressureGradient";

   Volume volume;
   Pressure pressure_sigmoid;
   Pressure pressure_linear;



   parameter Physiolib.Types.Pressure MinimalCollapsingPressure=0;
   parameter Physiolib.Types.Pressure a=MinimalCollapsingPressure/log(Modelica.Constants.eps);


   Volume BaseMeanVolume = FunctionalResidualCapacity + BaseTidalVolume/2  "Point of equality with linear presentation such as (FunctionalResidualCapacity + TidalVolume/2)";
   Pressure d = (BaseMeanVolume-ResidualVolume) * (VitalCapacity-(BaseMeanVolume-ResidualVolume)) / (Compliance*VitalCapacity);
   Pressure c = (BaseMeanVolume-FunctionalResidualCapacity)/Compliance + d*log((VitalCapacity/(BaseMeanVolume-ResidualVolume) - 1));

   //  Compliance * (-d*VC) = ((RV.-SelectedVolume).*(VC+RV.-SelectedVolume));  //rovnost derivace
  //(SelectedVolume.-ZeroPressureVolume)./Compliance = (-d.*log((VC./(SelectedVolume.-RV)).-1)+c); //rovnost hodnoty
   /*
   constant Integer n=3;
   parameter Volume SelectedVolume[n] = { FRC + TV/2, RV + VC/2, FRC};
   Pressure c[n](start = {0,0,0}),d[n](start = {1000,1000,1000});
   Pressure Pcl[n]=c-2*d, Pcu[n]=c+2*d;
   Pressure p[n];
   Volume zeroCrossingVolume[n];
  */
equation

  // p = 1*time - 5000;
  volume = 0.001*time + 0.00123+0.00001;

  pressure_linear =
      if ( useSigmoidCompliance) then
        smooth(0,
          if noEvent(volume>ResidualVolume) then
            (max( 0, volume - FunctionalResidualCapacity)/Compliance)
          else
            (a*log(max(Modelica.Constants.eps,volume/ResidualVolume))))
     else   -d*log((VitalCapacity/(volume-ResidualVolume))-1)+c;

  pressure_sigmoid =
      if ( not useSigmoidCompliance) then
        smooth(0,
          if noEvent(volume>ResidualVolume) then
            (max( 0, volume - FunctionalResidualCapacity)/Compliance)
          else
            (a*log(max(Modelica.Constants.eps,volume/ResidualVolume))))
     else   -d*log((VitalCapacity/(volume-ResidualVolume))-1)+c;


/*
  //vstupy zpv, compliance
  //nove vstupy RV = 0, VC = volume_start*2 

   p = -d*log((VC/(volume-RV))-1)+c
   //the same:  
   //V*ones(3) = RV + ( VC ./ (ones(3) + exp(-(p-c)./d));



  //vypocet d a c z Compliance a SelectedVolume
  Compliance * (-d*VC) = ((RV.-SelectedVolume).*(VC+RV.-SelectedVolume));  //rovnost derivace
  (SelectedVolume.-ZeroPressureVolume)./Compliance = (-d.*log((VC./(SelectedVolume.-RV)).-1)+c); //rovnost hodnoty 


  //prechod krivky cez 0:
  zeros(n) = (-d.*log((VC./(zeroCrossingVolume.-RV)).-1)+c);
*/

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    experiment(StopTime=15000, __Dymola_Algorithm="Dassl"));
end LungsCompartment;
