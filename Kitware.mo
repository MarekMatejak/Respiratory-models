within ;
package Kitware
  package Air_IdealGas
    "Gases O2,CO2,H2O,N2 with constant specific heat capacity"

    extends Modelica.Media.Interfaces.PartialMedium(
      final mediumName="Air (IdealGas)",
      final singleState=true,
      final reducedX=false,
      final fixedX=false,
      substanceNames = {"O2", "CO2", "H2O", "N2"},
      ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.pTX,
      reference_T = 310.15,
      reference_p = 101325,
      reference_X= x_default .* substanceData.MolarWeight ./ (x_default*substanceData.MolarWeight),
      Temperature(
        min=273,
        max=350,
        start=310.15));

    package stateOfMatter = Chemical.Interfaces.IdealGas
    "Substances model to translate data into substance properties";
    constant Modelica.SIunits.MoleFraction x_default[nCS]=
     {0.21,
      0.0004,
      0.02,
      0.7696} "Initial mole fractions of all substances (Please check: x_default*ones(nCS) = 1)";

    constant Integer nCS=4 "Number of chemical substances";

    constant stateOfMatter.SubstanceData substanceData[nCS] = {
      Chemical.Substances.Oxygen_gas(),
      Chemical.Substances.CarbonDioxide_gas(),
      Chemical.Substances.Water_gas(),
      Chemical.Substances.Nitrogen_gas()}
       "Definition of the substances";

    redeclare replaceable record extends ThermodynamicState
      "Thermodynamic state variables"
      AbsolutePressure p "Absolute pressure of medium";
      Temperature T "Temperature of medium";
      MassFraction[4] X(start=reference_X)
        "Mass fractions (= (component mass)/total mass  m_i/m)";
    end ThermodynamicState;

    redeclare model extends BaseProperties(final standardOrderComponents=true)
      "Base properties of medium"

    equation
      d = (X*ones(nCS)) / (X*(stateOfMatter.molarVolume(substanceData,T=T,p=p)./substanceData.MolarWeight));
      h = X*(stateOfMatter.molarEnthalpy(substanceData,T=T,p=p) ./ substanceData.MolarWeight);
      u = h - p/d;
      MM = (X*ones(nCS)) / ((X./substanceData.MolarWeight) * ones(nCS));
      R = 8.3144/MM;
      state.X = X;
      state.p = p;
      state.T = T;
    end BaseProperties;

    redeclare function extends setState_pTX
      "Return thermodynamic state as function of p, T and composition X or Xi"
    algorithm
      state.p := p;
      state.T := T;
      state.X := X;
    end setState_pTX;

    redeclare function extends setState_phX
      "Return thermodynamic state as function of p, h and composition X or Xi"
      input MassFraction X[Medium.nX] = reference_X;
    protected
      MolarMass MM=(X*ones(nCS)) / ((X./substanceData.MolarWeight) * ones(nCS));
      MoleFraction x[nCS]=(X./substanceData.MolarWeight) ./ (X*(ones(nCS)./substanceData.MolarWeight));
    algorithm
      state.p := p;
      state.T := stateOfMatter.solution_temperature(substanceData,h*MM,x,p);
      state.X := X;
    end setState_phX;

    redeclare function extends specificEnthalpy "Return specific enthalpy"
    algorithm
      h := state.X*(stateOfMatter.molarEnthalpy(substanceData,T=state.T,p=state.p) ./ substanceData.MolarWeight);
    end specificEnthalpy;

    redeclare function extends specificEntropy "Return specific entropy"
    protected
      Real a[nCS] "activity of substance";
      Modelica.SIunits.MolarEnergy u[nCS] "electro-chemical potential of substances in the solution";
    algorithm
      a := stateOfMatter.activityCoefficient(substanceData);

      u := stateOfMatter.chemicalPotentialPure(
          substanceData,
          state.T,
          state.p) + Modelica.Constants.R*state.T*log(a);

      s := state.X*(stateOfMatter.molarEntropy(u,substanceData,T=state.T,p=state.p) ./ substanceData.MolarWeight);
      annotation (Documentation(info="<html>

</html>"));
    end specificEntropy;

    redeclare function extends specificHeatCapacityCp
      "Return specific heat capacity at constant pressure"
    algorithm
      cp := state.X*(stateOfMatter.molarHeatCapacityCp(substanceData,T=state.T,p=state.p) ./ substanceData.MolarWeight);
      annotation (Documentation(info="<html>

</html>"));
    end specificHeatCapacityCp;

    redeclare function extends specificHeatCapacityCv
      "Return specific heat capacity at constant volume"
    algorithm
      cv := state.X*(stateOfMatter.molarHeatCapacityCv(substanceData,T=state.T,p=state.p) ./ substanceData.MolarWeight);
      annotation (Documentation(info="<html>

</html>"));
    end specificHeatCapacityCv;

    redeclare function extends density
    algorithm
      d := (state.X*ones(nCS)) / (state.X*(stateOfMatter.molarVolume(substanceData,T=state.T,p=state.p)./substanceData.MolarWeight));
    end density;

    redeclare function extends temperature
    algorithm
      T := state.T;
    end temperature;

    redeclare function extends pressure
    algorithm
      p := state.p;
    end pressure;

    function C_outflow
     input Modelica.SIunits.MassFraction x_mass[nCS];
     output Real C_outflow[nC];
    algorithm
      C_outflow := C_default;
      annotation(Inline=true);
    end C_outflow;

    function Xi_outflow
      input Modelica.SIunits.MassFraction x_mass[nCS];
      output Modelica.SIunits.MassFraction Xi[nXi];
    algorithm
      Xi := x_mass;
      annotation(Inline=true);
    end Xi_outflow;

    function x_mass
      input Modelica.SIunits.MassFraction actualStream_Xi[nXi];
      input Real actualStream_C[nC];
      output Modelica.SIunits.MassFraction x_mass[nCS];
    algorithm
      x_mass := actualStream_Xi;
      annotation(Inline=true);
    end x_mass;

    function concentration "Concentration of base substance molecules from Xi and C"
      input ThermodynamicState state;
      input Modelica.SIunits.MassFraction Xi[nXi];
      input Real C[nC];
      output Modelica.SIunits.Concentration concentration[nCS];
    algorithm
      concentration := (Xi./substanceData.MolarWeight)/(Xi*(stateOfMatter.molarVolume(substanceData,T=state.T,p=state.p)./substanceData.MolarWeight));
    end concentration;

    function specificEnthalpyOffsets "Difference between chemical substance enthalpy and medium substance enthalpy at temperature 298.15 K and 100kPa"
      input Modelica.SIunits.ElectricPotential v=0;
      input Real I=0;
      output SpecificEnthalpy h[nCS];
    algorithm
      h := zeros(nCS);
    end specificEnthalpyOffsets;
    annotation (Documentation(info="<html>
<p>
This package is a <strong>template</strong> for <strong>new medium</strong> models. For a new
medium model just make a copy of this package, remove the
\"partial\" keyword from the package and provide
the information that is requested in the comments of the
Modelica source.
</p>
</html>"));
  end Air_IdealGas;

  model RespirationMuscle "Respiration model"
    extends Modelica.Icons.Example;

    import Modelica.SIunits.*;

    replaceable package Air = Chemical.Media.Air_MixtureGasNasa;

    parameter Frequency RespirationRate(displayUnit="1/s")=8/60                                              "Respiration rate";
    parameter Volume ResidualVolume=0.0013                                                      "Lungs residual volume";

    parameter Volume FunctionalResidualCapacity=0.00231                                                      "Functional residual capacity";
    parameter Physiolibrary.Types.HydraulicResistance TotalResistance=147099.75
      "Total lungs pathways conductance";
    parameter Physiolibrary.Types.HydraulicCompliance TotalCompliance=1.0197162129779e-06
                                                   "Total lungs compliance";

    parameter Pressure Pmin=-3432.3275                                      "Relative external lungs pressure minimum caused by respiratory muscles";
    parameter Pressure Pmax=1078.7315                "Relative external lungs pressure maximum";
    parameter Real RespiratoryMusclePressureCycle[:,3] = {
    {0,0,0},
  {0.02/8,0,Pmin/(0.2/8)},
  {0.4/8,0.9*Pmin,0.1*Pmin/(0.3/8)},
  {1/8,Pmin,0},
  {2/8,Pmin,0},
  {2.2/8,0.5*Pmin,-Pmin/(0.2/8)},
  {2.3/8,0.1*Pmin,-0.1*(Pmin)/(0.2/8)},
  {3/8,0,0},
  {4/8,0,0},
  {4.02/8,0.1*Pmax,0.1*Pmax/(0.4/8)},
    {4.4/8,0.9*Pmax,0.1*Pmax/(0.3/8)},
    {4.8/8,Pmax,0},
    {6/8,Pmax,0},
    {6.02/8,Pmax,-Pmax/(0.2/8)},
    {6.3/8,0.1*Pmax,-0.1*Pmax/(0.2/8)},
    {7/8,0,0},
    {1,0,0}}  "Absolute external lungs pressure during respiration cycle scaled to time period (0,1)";
         /* {0,Pmax,0},
        {3/8,Pmin,0},
        {1,Pmax,0}}*/

    parameter Temperature CoreTemperature=310.15                       "body temperature";
    parameter Temperature EnvironmentTemperature=298.15                       "external air temperature";

    parameter Volume LungsAirVolume_initial = FunctionalResidualCapacity;

    parameter Density d = Air.density(Air.setState_pTX(system.p_ambient+Pmax,CoreTemperature));

    Physiolibrary.Blocks.Source.PeriodicCurveSource respiratoryMusclePressureCycle(data=
          RespiratoryMusclePressureCycle, nu=4*9*5)
      "Relative position in respiratory cycle (0,1) to absolute external lungs pressure"
      annotation (Placement(transformation(extent={{18,54},{38,74}})));
                                                                                                                //0.0133,

    inner Modelica.Fluid.System system(T_ambient=CoreTemperature)
                                       "Human body system setting"
      annotation (Placement(transformation(extent={{60,66},{80,86}})));

    Physiolibrary.Types.Constants.FrequencyConst frequency(k=RespirationRate)
      annotation (Placement(transformation(extent={{-46,62},{-38,70}})));
      Pressure transMusclePressue;
  equation
     transMusclePressue = respiratoryMusclePressureCycle.val;
    connect(frequency.y, respiratoryMusclePressureCycle.frequence)
      annotation (Line(points={{-37,66},{-8,66},{-8,64},{18,64}},    color={0,0,127}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}})),                                        Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
      experiment(StopTime=16, __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>References:</p>
<p><br>Mecklenburgh, J. S., and W. W. Mapleson. &quot;Ventilatory assistance and respiratory muscle activity. 1: Interaction in healthy volunteers.&quot; <i>British journal of anaesthesia</i> 80.4 (1998): 422-433.</p>
</html>"));
  end RespirationMuscle;

  model LungsCompliance "Respiration model"
    extends Modelica.Icons.Example;
    Physiolibrary.Fluid.Components.ElasticVessel lungs(
      redeclare package Medium = Air,
      volume_start=ResidualVolume + VitalCapacity - 0.00023,
      ZeroPressureVolume=FunctionalResidualCapacity,
      ResidualVolume(displayUnit="l") = ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity(displayUnit="l") = VitalCapacity,
      BaseTidalVolume=0.0005,
      EnthalpyNotUsed=true,
      nPorts=1) "Lungs"
      annotation (Placement(transformation(extent={{36,-28},{56,-8}})));

    import Modelica.SIunits.*;
    replaceable package Air = Kitware.Air_IdealGas;
  //  Chemical.Media.SimpleAir_C;
   //  Chemical.Media.Air_MixtureGasNasa;
                              //Chemical.Media.Water_Incompressible; //

    parameter Frequency RespirationRate(displayUnit="1/s")=8/60                                              "Respiration rate";
    parameter Volume ResidualVolume=0.0013                                                      "Lungs residual volume";

    parameter Volume FunctionalResidualCapacity=0.00231                                                      "Functional residual capacity";
    parameter Volume VitalCapacity=0.00493                                                       "Functional residual capacity";

    parameter Physiolibrary.Types.HydraulicResistance TotalResistance=1
      "Total lungs pathways conductance";
    parameter Physiolibrary.Types.HydraulicCompliance TotalCompliance(
        displayUnit="l/cmH2O")=1.0197162129779e-06 "Total lungs compliance";

    parameter Pressure Pmin=-3432.3275                                      "Relative external lungs pressure minimum caused by respiratory muscles";
    parameter Pressure Pmax=1078.7315                "Relative external lungs pressure maximum";


    parameter Temperature CoreTemperature=310.15                       "body temperature";
    parameter Temperature EnvironmentTemperature=298.15                       "external air temperature";

    parameter Volume LungsAirVolume_initial = FunctionalResidualCapacity;

  //  parameter Density d = Air.density(Air.setState_pTX(system.p_ambient+Pmax,CoreTemperature));
                                                                                                                //0.0133,

    inner Modelica.Fluid.System system(T_ambient=CoreTemperature)
                                       "Human body system setting"
      annotation (Placement(transformation(extent={{60,66},{80,86}})));

    Modelica.Blocks.Math.Add add annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={54,28})));
    Physiolibrary.Types.Constants.PressureConst ambient_pressure(k=system.p_ambient)
      annotation (Placement(transformation(extent={{84,50},{74,58}})));
      Pressure pressue;
    Modelica.Blocks.Sources.Clock clock(offset=-4000)
      annotation (Placement(transformation(extent={{-38,50},{-18,70}})));
    Physiolibrary.Fluid.Sources.PressureSource environment(redeclare package
        Medium = Air, T=EnvironmentTemperature)
      "External environment"
      annotation (Placement(transformation(extent={{-76,-30},{-56,-10}})));
    Physiolibrary.Fluid.Sensors.FlowMeasure airflowMeasure(redeclare package
        Medium = Air, EnthalpyNotUsed=true)
                              "Lungs pathway airflow"
      annotation (Placement(transformation(extent={{-40,-30},{-20,-10}})));
    Physiolibrary.Fluid.Components.Resistor resistor(redeclare package Medium =
          Air,
      EnthalpyNotUsed=true,
      Resistance(displayUnit="(Pa.s)/m3") = 1)
      annotation (Placement(transformation(extent={{-6,-30},{14,-10}})));
  equation
     pressue = clock.y;
    connect(add.y, lungs.externalPressure)
      annotation (Line(points={{54,17},{54,4},{54,-9},{52,-9}},
                                                 color={0,0,127}));
    connect(ambient_pressure.y, add.u1)
      annotation (Line(points={{72.75,54},{60,54},{60,40}}, color={0,0,127}));
    connect(clock.y, add.u2)
      annotation (Line(points={{-17,60},{48,60},{48,40}}, color={0,0,127}));
    connect(airflowMeasure.q_in,environment. y) annotation (Line(
        points={{-40,-20},{-56,-20}},
        color={127,0,0},
        thickness=0.5));
    connect(airflowMeasure.q_out,resistor. q_in) annotation (Line(
        points={{-20,-20},{-6,-20}},
        color={127,0,0},
        thickness=0.5));
    connect(resistor.q_out, lungs.q_in[1]) annotation (Line(
        points={{14,-20},{14,-18},{45.9,-18}},
        color={127,0,0},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}})),                                        Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
      experiment(StopTime=8000, __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>References:</p>
<p><br>Mecklenburgh, J. S., and W. W. Mapleson. &quot;Ventilatory assistance and respiratory muscle activity. 1: Interaction in healthy volunteers.&quot; <i>British journal of anaesthesia</i> 80.4 (1998): 422-433.</p>
</html>"));
  end LungsCompliance;

  model LungsCompliance2 "Respiration model"
    extends Modelica.Icons.Example;
    Physiolibrary.Fluid.Components.ElasticVessel lungs(
      redeclare package Medium = Air,
      volume_start=ResidualVolume + 1e-5,
      ZeroPressureVolume=FunctionalResidualCapacity,
      ResidualVolume(displayUnit="l") = ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=false,
      useSigmoidCompliance=true,
      VitalCapacity(displayUnit="l") = VitalCapacity,
      BaseTidalVolume=0.0005,
      EnthalpyNotUsed=true,
      nPorts=2) "Lungs"
      annotation (Placement(transformation(extent={{36,-28},{56,-8}})));

    import Modelica.SIunits.*;
    replaceable package Air = Chemical.Media.Water_Incompressible;// Chemical.Media.SimpleAir_C;
     //Chemical.Media.Air_MixtureGasNasa;

    parameter Frequency RespirationRate(displayUnit="1/s")=8/60                                              "Respiration rate";
    parameter Volume ResidualVolume=0.0013                                                      "Lungs residual volume";

    parameter Volume FunctionalResidualCapacity=0.00231                                                      "Functional residual capacity";
    parameter Volume VitalCapacity=0.00493                                                       "Functional residual capacity";

    parameter Physiolibrary.Types.HydraulicResistance TotalResistance=1
      "Total lungs pathways conductance";
    parameter Physiolibrary.Types.HydraulicCompliance TotalCompliance(
        displayUnit="l/cmH2O")=1.0197162129779e-06 "Total lungs compliance";

    parameter Pressure Pmin=-3432.3275                                      "Relative external lungs pressure minimum caused by respiratory muscles";
    parameter Pressure Pmax=1078.7315                "Relative external lungs pressure maximum";

    parameter Temperature CoreTemperature=310.15                       "body temperature";
    parameter Temperature EnvironmentTemperature=298.15                       "external air temperature";

    parameter Volume LungsAirVolume_initial = FunctionalResidualCapacity;

    parameter Density d = Air.density(Air.setState_pTX(system.p_ambient+Pmax,CoreTemperature));
                                                                                                                //0.0133,

    Physiolibrary.Fluid.Sensors.PressureMeasure lungsPressureMeasure(
        redeclare package Medium = Air) "Lungs pressure"
      annotation (Placement(transformation(extent={{70,-20},{90,0}})));

    inner Modelica.Fluid.System system(T_ambient=CoreTemperature)
                                       "Human body system setting"
      annotation (Placement(transformation(extent={{60,66},{80,86}})));

    Modelica.Blocks.Math.Add add annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={54,28})));
    Physiolibrary.Types.Constants.PressureConst ambient_pressure(k=system.p_ambient)
      annotation (Placement(transformation(extent={{84,50},{74,58}})));
      Pressure pressue;
    Modelica.Blocks.Sources.Clock clock(offset=-4000)
      annotation (Placement(transformation(extent={{-38,50},{-18,70}})));
    Physiolibrary.Fluid.Sources.PressureSource environment(redeclare package
        Medium = Air,
      usePressureInput=true,
                      T=EnvironmentTemperature)
      "External environment"
      annotation (Placement(transformation(extent={{-76,-30},{-56,-10}})));
    Physiolibrary.Fluid.Components.Conductor conductor(redeclare package Medium =
          Air,
      EnthalpyNotUsed=true,
      Conductance(displayUnit="l/(cmH2O.s)") = 1.0197162129779)
      annotation (Placement(transformation(extent={{-24,-30},{-4,-10}})));
  equation
     pressue = clock.y;
    connect(lungsPressureMeasure.q_in, lungs.q_in[1]) annotation (Line(
        points={{76,-16},{76,-16.7},{45.9,-16.7}},
        color={127,0,0},
        thickness=0.5));
    connect(ambient_pressure.y, add.u1)
      annotation (Line(points={{72.75,54},{60,54},{60,40}}, color={0,0,127}));
    connect(clock.y, add.u2)
      annotation (Line(points={{-17,60},{48,60},{48,40}}, color={0,0,127}));
    connect(add.y, environment.pressure) annotation (Line(points={{54,17},{54,0},{-82,
            0},{-82,-20},{-76,-20}}, color={0,0,127}));
    connect(environment.y, conductor.q_in) annotation (Line(
        points={{-56,-20},{-24,-20}},
        color={127,0,0},
        thickness=0.5));
    connect(conductor.q_out, lungs.q_in[2]) annotation (Line(
        points={{-4,-20},{10,-20},{10,-19.3},{45.9,-19.3}},
        color={127,0,0},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}})),                                        Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
      experiment(StopTime=8000, __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>References:</p>
<p><br>Mecklenburgh, J. S., and W. W. Mapleson. &quot;Ventilatory assistance and respiratory muscle activity. 1: Interaction in healthy volunteers.&quot; <i>British journal of anaesthesia</i> 80.4 (1998): 422-433.</p>
</html>"));
  end LungsCompliance2;

  model Respiration "Respiration model"
    extends Modelica.Icons.Example;
    Physiolibrary.Fluid.Components.ElasticVessel lungs(
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      ZeroPressureVolume=FunctionalResidualCapacity,
      ResidualVolume(displayUnit="l") = ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity(displayUnit="l") = TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=TidalVolume,
      nPorts=2) "Lungs"
      annotation (Placement(transformation(extent={{36,-28},{56,-8}})));

    import Modelica.SIunits.*;

    replaceable package Air =  Chemical.Media.Water_Incompressible; //Chemical.Media.Air_MixtureGasNasa;

    parameter Frequency RespirationRate(displayUnit="1/s")=8/60                                              "Respiration rate";
    parameter Volume ResidualVolume=0.0013                                                      "Lungs residual volume";
    parameter Volume TotalLungCapacity=0.00623                                                      "Total Lung Capacity";
    parameter Volume TidalVolume=0.0005                                                      "Base Tidal Volume";



    parameter Volume FunctionalResidualCapacity=0.00231                                                      "Functional residual capacity";
    parameter Physiolibrary.Types.HydraulicResistance TotalResistance=147099.75
      "Total lungs pathways conductance";
    parameter Physiolibrary.Types.HydraulicCompliance TotalCompliance=1.0197162129779e-06
                                                   "Total lungs compliance";

    parameter Pressure Pmin=-3432.3275                                      "Relative external lungs pressure minimum caused by respiratory muscles";
    parameter Pressure Pmax=1078.7315                "Relative external lungs pressure maximum";
    parameter Real RespiratoryMusclePressureCycle[:,3] = {
    {0,0,0},
  {0.02/8,0,Pmin/(0.2/8)},
  {0.4/8,0.9*Pmin,0.1*Pmin/(0.3/8)},
  {1/8,Pmin,0},
  {2/8,Pmin,0},
  {2.2/8,0.5*Pmin,-Pmin/(0.2/8)},
  {2.3/8,0.1*Pmin,-0.1*(Pmin)/(0.2/8)},
  {3/8,0,0},
  {4/8,0,0},
  {4.02/8,0.1*Pmax,0.1*Pmax/(0.4/8)},
    {4.4/8,0.9*Pmax,0.1*Pmax/(0.3/8)},
    {4.8/8,Pmax,0},
    {6/8,Pmax,0},
    {6.02/8,Pmax,-Pmax/(0.2/8)},
    {6.3/8,0.1*Pmax,-0.1*Pmax/(0.2/8)},
    {7/8,0,0},
    {1,0,0}}  "External lungs pressure during respiration cycle scaled to time period (0,1)";
  //{4.4/8,0.9*Pmax,0.1*Pmax/(0.3/8)},
  //{4.8/8,Pmax,0},
  //{6/8,Pmax,0},
  //{6.02/8,Pmax,-Pmax/(0.2/8)},
  //{6.3/8,0.1*Pmax,-0.1*Pmax/(0.2/8)},
  //{7/8,0,0},
         /* {0,Pmax,0},
        {3/8,Pmin,0},
        {1,Pmax,0}}*/

    parameter Temperature CoreTemperature=310.15                       "body temperature";
    parameter Temperature EnvironmentTemperature=298.15                       "external air temperature";

    parameter Volume LungsAirVolume_initial = FunctionalResidualCapacity;

    parameter Density d = Air.density(Air.setState_pTX(system.p_ambient+Pmax,CoreTemperature));

    Physiolibrary.Blocks.Source.PeriodicCurveSource respiratoryMusclePressureCycle(data=
          RespiratoryMusclePressureCycle)
      "Relative position in respiratory cycle (0,1) to absolute external lungs pressure"
      annotation (Placement(transformation(extent={{18,54},{38,74}})));
                                                                                                                //0.0133,

    Physiolibrary.Fluid.Sensors.PressureMeasure lungsPressureMeasure(
        redeclare package Medium = Air) "Lungs pressure"
      annotation (Placement(transformation(extent={{70,-20},{90,0}})));

    inner Modelica.Fluid.System system(T_ambient=CoreTemperature)
                                       "Human body system setting"
      annotation (Placement(transformation(extent={{60,66},{80,86}})));

    Physiolibrary.Fluid.Sources.PressureSource environment(redeclare package
        Medium =         Air, T=EnvironmentTemperature)
      "External environment"
      annotation (Placement(transformation(extent={{-76,-30},{-56,-10}})));

    Physiolibrary.Fluid.Sensors.FlowMeasure airflowMeasure(redeclare package
        Medium =         Air) "Lungs pathway airflow"
      annotation (Placement(transformation(extent={{-40,-30},{-20,-10}})));

    Physiolibrary.Fluid.Components.Resistor resistor(redeclare package Medium =
          Air, Resistance=TotalResistance)
      annotation (Placement(transformation(extent={{-6,-30},{14,-10}})));

    Physiolibrary.Types.Constants.FrequencyConst frequency(k=RespirationRate)
      annotation (Placement(transformation(extent={{-46,62},{-38,70}})));
    Modelica.Blocks.Math.Add add annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={54,28})));
    Physiolibrary.Types.Constants.PressureConst ambient_pressure(k=system.p_ambient)
      annotation (Placement(transformation(extent={{84,50},{74,58}})));
      Pressure transMusclePressue;
  equation
     transMusclePressue = respiratoryMusclePressureCycle.val;
    connect(lungsPressureMeasure.q_in, lungs.q_in[1]) annotation (Line(
        points={{76,-16},{76,-16.7},{45.9,-16.7}},
        color={127,0,0},
        thickness=0.5));
    connect(airflowMeasure.q_in, environment.y) annotation (Line(
        points={{-40,-20},{-56,-20}},
        color={127,0,0},
        thickness=0.5));
    connect(airflowMeasure.q_out, resistor.q_in) annotation (Line(
        points={{-20,-20},{-6,-20}},
        color={127,0,0},
        thickness=0.5));
    connect(resistor.q_out, lungs.q_in[2]) annotation (Line(
        points={{14,-20},{28,-20},{28,-19.3},{45.9,-19.3}},
        color={127,0,0},
        thickness=0.5));
    connect(frequency.y, respiratoryMusclePressureCycle.frequence)
      annotation (Line(points={{-37,66},{-8,66},{-8,64},{18,64}},    color={0,0,127}));
    connect(add.y, lungs.externalPressure)
      annotation (Line(points={{54,17},{54,4},{54,-9},{52,-9}},
                                                 color={0,0,127}));
    connect(add.u2, respiratoryMusclePressureCycle.val)
      annotation (Line(points={{48,40},{48,64},{38,64}}, color={0,0,127}));
    connect(ambient_pressure.y, add.u1)
      annotation (Line(points={{72.75,54},{60,54},{60,40}}, color={0,0,127}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}})),                                        Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
      experiment(StopTime=16, __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>References:</p>
<p><br>Mecklenburgh, J. S., and W. W. Mapleson. &quot;Ventilatory assistance and respiratory muscle activity. 1: Interaction in healthy volunteers.&quot; <i>British journal of anaesthesia</i> 80.4 (1998): 422-433.</p>
</html>"));
  end Respiration;

  model MinimalRespiration "Minimal respiration model"
    extends Modelica.Icons.Example;

    import Modelica.SIunits.*;

    replaceable package Air = Chemical.Media.Air_MixtureGasNasa;

    parameter Frequency RespirationRate=0.133                                                                "Respiration rate";
    parameter Volume ResidualVolume=0.0013                                                      "Lungs residual volume";
    parameter Volume TotalLungCapacity=0.00623                                                      "Total Lung Capacity";
    parameter Volume BaseTidalVolume=0.0005                                                      "Base Tidal Volume";
    parameter Volume LungsAirVolume_initial = FunctionalResidualCapacity;

    parameter Volume FunctionalResidualCapacity=0.00231                                                      "Functional residual capacity";
    parameter Physiolibrary.Types.HydraulicResistance TotalResistance(displayUnit=
         "(cmH2O.s)/l")=147099.75
      "Total lungs pathways conductance";
    parameter Physiolibrary.Types.HydraulicCompliance TotalCompliance(displayUnit=
         "l/cmH2O")=1.0197162129779e-06            "Total lungs compliance";

    parameter Pressure Pmin=-533.28954966                                   "Relative external lungs pressure minimum caused by respiratory muscles";
    parameter Pressure Pmax(displayUnit="mmHg")=0    "Relative external lungs pressure maximum";
    parameter Real RespiratoryMusclePressureCycle[:,3] = {
          {0,system.p_ambient + Pmax,0},
          {3/8,system.p_ambient + Pmin,0},
          {1,system.p_ambient + Pmax,0}}
            "Absolute external lungs pressure during respiration cycle scaled to time period (0,1)";

    parameter Temperature CoreTemperature=310.15                       "body temperature";
    parameter Temperature EnvironmentTemperature=298.15                       "external air temperature";



    parameter Density d = Air.density(Air.setState_pTX(system.p_ambient+Pmax,CoreTemperature));

    Physiolibrary.Blocks.Source.PeriodicCurveSource respiratoryMusclePressureCycle(data=
          RespiratoryMusclePressureCycle)
      "Relative position in respiratory cycle (0,1) to absolute external lungs pressure"
      annotation (Placement(transformation(extent={{18,54},{38,74}})));

   // parameter Mass m_initial = LungsAirVolume_initial*Air.density(Air.setState_pTX(system.p_ambient
   //        + Pmax, CoreTemperature));

    Physiolibrary.Fluid.Components.ElasticVessel lungs(
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      ZeroPressureVolume=FunctionalResidualCapacity,
      ResidualVolume=ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=true,
      nPorts=2) "Lungs"
      annotation (Placement(transformation(extent={{36,-28},{56,-8}})));                                        //0.0133,

    Physiolibrary.Fluid.Sensors.PressureMeasure lungsPressureMeasure(
        redeclare package Medium = Air) "Lungs pressure"
      annotation (Placement(transformation(extent={{70,-20},{90,0}})));

    inner Modelica.Fluid.System system(T_ambient=CoreTemperature)
                                       "Human body system setting"
      annotation (Placement(transformation(extent={{60,66},{80,86}})));

    Physiolibrary.Fluid.Sources.PressureSource environment(redeclare package
        Medium =         Air, T=EnvironmentTemperature)
      "External environment"
      annotation (Placement(transformation(extent={{-78,-32},{-58,-12}})));

    Physiolibrary.Fluid.Sensors.FlowMeasure airflowMeasure(redeclare package
        Medium =         Air) "Lungs pathway airflow"
      annotation (Placement(transformation(extent={{-42,-32},{-20,-12}})));

    Physiolibrary.Fluid.Components.Resistor resistor(redeclare package Medium =
          Air, Resistance=TotalResistance)
      annotation (Placement(transformation(extent={{-6,-32},{14,-12}})));

    Physiolibrary.Types.Constants.FrequencyConst frequency(k=RespirationRate)
      annotation (Placement(transformation(extent={{-46,62},{-38,70}})));
    Chemical.Sources.SubstanceOutflow O2(SubstanceFlow(displayUnit="mmol/min") = 0.000257)
      annotation (Placement(transformation(extent={{46,-70},{66,-50}})));
    Chemical.Sources.SubstanceInflow CO2(SubstanceFlow(displayUnit="mmol/min") = 0.00020566666666667)
      annotation (Placement(transformation(extent={{0,-70},{20,-50}})));
  equation

    connect(lungsPressureMeasure.q_in, lungs.q_in[1]) annotation (Line(
        points={{76,-16},{76,-16.7},{45.9,-16.7}},
        color={127,0,0},
        thickness=0.5));
    connect(resistor.q_out, lungs.q_in[2]) annotation (Line(
        points={{14,-22},{28,-22},{28,-19.3},{45.9,-19.3}},
        color={127,0,0},
        thickness=0.5));
    connect(frequency.y, respiratoryMusclePressureCycle.frequence)
      annotation (Line(points={{-37,66},{-8,66},{-8,64},{18,64}},    color={0,0,127}));
    connect(respiratoryMusclePressureCycle.val, lungs.externalPressure)
      annotation (Line(points={{38,64},{52,64},{52,-9}},color={0,0,127}));
    connect(lungs.substances[1], O2.port_a) annotation (Line(points={{36,-18},{36,-60},
            {46,-60}},      color={158,66,200}));
    connect(CO2.port_b, lungs.substances[2]) annotation (Line(points={{20,-60},{30,
            -60},{30,-18},{36,-18}},   color={158,66,200}));
    connect(airflowMeasure.q_in, environment.y) annotation (Line(
        points={{-42,-22},{-58,-22}},
        color={127,0,0},
        thickness=0.5));
    connect(airflowMeasure.q_out, resistor.q_in) annotation (Line(
        points={{-20,-22},{-6,-22}},
        color={127,0,0},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}})),                                        Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
      experiment(StopTime=16, __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>References:</p>
<p><br>Mecklenburgh, J. S., and W. W. Mapleson. &quot;Ventilatory assistance and respiratory muscle activity. 1: Interaction in healthy volunteers.&quot; <i>British journal of anaesthesia</i> 80.4 (1998): 422-433.</p>
</html>"));
  end MinimalRespiration;

  model Respiration2 "Minimal respiration model"
    extends Modelica.Icons.Example;

    import Modelica.SIunits.*;

    replaceable package Air = Chemical.Media.Air_MixtureGasNasa;
    replaceable package PleuralFluid = Chemical.Media.Water_Incompressible;


    parameter Pressure IntrathoraxPressure = system.p_ambient - 700;
    parameter Frequency RespirationRate=0.133                                                                "Respiration rate";
    parameter Volume ResidualVolume=0.0013                                                      "Lungs residual volume";
    parameter Volume TotalLungCapacity=0.00623                                                      "Total Lung Capacity";
    parameter Volume BaseTidalVolume=0.0005                                                      "Base Tidal Volume";
    parameter Volume LungsAirVolume_initial = FunctionalResidualCapacity;
    parameter Volume PleuralVolume_initial = 0.0001 "Initial volume of pleaural fluid";
     parameter Physiolibrary.Types.HydraulicCompliance PleuralCompliance(displayUnit=
         "l/cmH2O")=1.0197162129779e-06            "Total lungs compliance";

    parameter Volume FunctionalResidualCapacity=0.00231                                                      "Functional residual capacity";
    parameter Physiolibrary.Types.HydraulicResistance TotalResistance(displayUnit=
         "(cmH2O.s)/l")=147099.75
      "Total lungs pathways conductance";
    parameter Physiolibrary.Types.HydraulicCompliance TotalCompliance(displayUnit=
         "l/cmH2O")=1.0197162129779e-06            "Total lungs compliance";

    parameter Pressure Pmin=-533.28954966                                   "Relative external lungs pressure minimum caused by respiratory muscles";
    parameter Pressure Pmax(displayUnit="mmHg")=0    "Relative external lungs pressure maximum";
    parameter Real RespiratoryMusclePressureCycle[:,3] = {
          {0,Pmax,0},
          {3/8,Pmin,0},
          {1,Pmax,0}}
            "Absolute external lungs pressure during respiration cycle scaled to time period (0,1)";

    parameter Temperature CoreTemperature=310.15                       "body temperature";
    parameter Temperature EnvironmentTemperature=CoreTemperature                       "external air temperature";

    parameter Density d = Air.density(Air.setState_pTX(system.p_ambient+Pmax,CoreTemperature));

    Physiolibrary.Blocks.Source.PeriodicCurveSource respiratoryMusclePressureCycle(data=
          RespiratoryMusclePressureCycle)
      "Relative position in respiratory cycle (0,1) to absolute external lungs pressure"
      annotation (Placement(transformation(extent={{-62,80},{-42,100}})));

    parameter Mass m_initial = LungsAirVolume_initial*Air.density(Air.setState_pTX(system.p_ambient
           + Pmax, CoreTemperature));

    Physiolibrary.Fluid.Components.ElasticVessel alveoli(
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      ZeroPressureVolume=FunctionalResidualCapacity,
      ResidualVolume=ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=true,
      nPorts=4) "Alveolar space of lungs"
      annotation (Placement(transformation(extent={{48,-26},{68,-6}})));                                        //0.0133,

    Physiolibrary.Fluid.Sensors.PressureMeasure alveolarPressure(redeclare
        package Medium = Air) "Alveolar pressure"
      annotation (Placement(transformation(extent={{80,-20},{100,0}})));

    inner Modelica.Fluid.System system(T_ambient=CoreTemperature)
                                       "Human body system setting"
      annotation (Placement(transformation(extent={{60,66},{80,86}})));

    Physiolibrary.Fluid.Sources.PressureSource environment(redeclare package
        Medium =         Air, T=EnvironmentTemperature)
      "External environment"
      annotation (Placement(transformation(extent={{-98,-28},{-78,-8}})));

    Physiolibrary.Fluid.Sensors.FlowMeasure airflowMeasure(redeclare package
        Medium =         Air) "Lungs pathway airflow"
      annotation (Placement(transformation(extent={{-72,-28},{-50,-8}})));

    Physiolibrary.Fluid.Components.Resistor resistor(redeclare package Medium =
          Air, Resistance=TotalResistance)
      annotation (Placement(transformation(extent={{-42,-28},{-22,-8}})));

    Physiolibrary.Types.Constants.FrequencyConst frequency(k=RespirationRate)
      annotation (Placement(transformation(extent={{-96,86},{-88,94}})));
    Chemical.Sources.SubstanceOutflow O2(SubstanceFlow(displayUnit="mmol/min") = 0.000257)
      annotation (Placement(transformation(extent={{66,-80},{86,-60}})));
    Chemical.Sources.SubstanceInflow CO2(SubstanceFlow(displayUnit="mmol/min") = 0.00020566666666667)
      annotation (Placement(transformation(extent={{-8,-92},{12,-72}})));
    Physiolibrary.Fluid.Sensors.PressureMeasure pleuralPressure(redeclare
        package
        Medium = PleuralFluid, GetAbsolutePressure=true) "Pleural pressure"
      annotation (Placement(transformation(extent={{-2,42},{18,62}})));
    Modelica.Blocks.Math.Add add annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-4,72})));
    Physiolibrary.Types.Constants.PressureConst ambient_pressure(k=
          IntrathoraxPressure)
      annotation (Placement(transformation(extent={{30,88},{20,96}})));



    Physiolibrary.Fluid.Sensors.FlowMeasure AlvelarFlowMeasure(redeclare package
        Medium = Air)
      annotation (Placement(transformation(extent={{-8,-30},{14,-10}})));
    Modelica.Fluid.Sensors.MassFractions alveolarO2(redeclare package Medium =
          Air, substanceName="O2")
      annotation (Placement(transformation(extent={{12,18},{32,38}})));
    Modelica.Fluid.Sensors.MassFractions alveolarCO2(redeclare package Medium =
          Air, substanceName="CO2")
      annotation (Placement(transformation(extent={{14,-10},{34,10}})));
    Physiolibrary.Fluid.Components.ElasticVessel pleural(
      redeclare package Medium = PleuralFluid,
      volume_start=PleuralVolume_initial,
      ZeroPressureVolume=PleuralVolume_initial,
      Compliance=PleuralCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=false,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=false,
      nPorts=1) "Pleaural space with lungs"
      annotation (Placement(transformation(extent={{-20,36},{0,56}})));
  equation



    connect(alveolarPressure.q_in,alveoli. q_in[1]) annotation (Line(
        points={{86,-16},{86,-14.05},{57.9,-14.05}},
        color={127,0,0},
        thickness=0.5));
    connect(frequency.y, respiratoryMusclePressureCycle.frequence)
      annotation (Line(points={{-87,90},{-62,90}},                   color={0,0,127}));
    connect(alveoli.substances[1], O2.port_a) annotation (Line(points={{48,-16},{36,
            -16},{36,-70},{66,-70}},                       color={158,66,200}));
    connect(CO2.port_b,alveoli. substances[2]) annotation (Line(points={{12,-82},{36,
            -82},{36,-16},{48,-16}},      color={158,66,200}));
    connect(airflowMeasure.q_in, environment.y) annotation (Line(
        points={{-72,-18},{-78,-18}},
        color={127,0,0},
        thickness=0.5));
    connect(airflowMeasure.q_out, resistor.q_in) annotation (Line(
        points={{-50,-18},{-42,-18}},
        color={127,0,0},
        thickness=0.5));
    connect(ambient_pressure.y,add. u1)
      annotation (Line(points={{18.75,92},{2,92},{2,84}},   color={0,0,127}));
    connect(add.u2, respiratoryMusclePressureCycle.val)
      annotation (Line(points={{-10,84},{-8,84},{-8,90},{-42,90}},
                                                         color={0,0,127}));
    connect(AlvelarFlowMeasure.q_out, alveoli.q_in[2]) annotation (Line(
        points={{14,-20},{34,-20},{34,-14},{40,-14},{40,-15.35},{57.9,-15.35}},
        color={127,0,0},
        thickness=0.5));
    connect(alveolarO2.port,alveoli. q_in[3]) annotation (Line(points={{22,18},{40,
            18},{40,-16.65},{57.9,-16.65}},     color={0,127,255}));
    connect(alveolarCO2.port,alveoli. q_in[4]) annotation (Line(points={{24,-10},{42,
            -10},{42,-17.95},{57.9,-17.95}},     color={0,127,255}));
    connect(resistor.q_out, AlvelarFlowMeasure.q_in) annotation (Line(
        points={{-22,-18},{-22,-20},{-8,-20}},
        color={127,0,0},
        thickness=0.5));
    connect(add.y, pleural.externalPressure)
      annotation (Line(points={{-4,61},{-4,55}}, color={0,0,127}));
    connect(pleural.q_in[1], pleuralPressure.q_in) annotation (Line(
        points={{-10.1,46},{4,46}},
        color={127,0,0},
        thickness=0.5));
    connect(pleuralPressure.pressure, alveoli.externalPressure)
      annotation (Line(points={{14,48},{64,48},{64,-7}}, color={0,0,127}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}})),                                        Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
      experiment(
        StopTime=5000,
        __Dymola_NumberOfIntervals=5000,
        Tolerance=1e-06,
        __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>References:</p>
<p><br>Mecklenburgh, J. S., and W. W. Mapleson. &quot;Ventilatory assistance and respiratory muscle activity. 1: Interaction in healthy volunteers.&quot; <i>British journal of anaesthesia</i> 80.4 (1998): 422-433.</p>
</html>"));
  end Respiration2;

  model Respiration3 "Respiration model"
    extends Modelica.Icons.Example;

    import Modelica.SIunits.*;

    replaceable package Air = Chemical.Media.Air_MixtureGasNasa;
    replaceable package PleuralFluid = Chemical.Media.Water_Incompressible;

    parameter Pressure IntrathoraxPressure = system.p_ambient - 700;
    parameter Frequency RespirationRate=0.133                                                                "Respiration rate";
    parameter Volume ResidualVolume=0.0013                                                      "Lungs residual volume";
    parameter Volume TotalLungCapacity=0.00623                                                      "Total Lung Capacity";
    parameter Volume BaseTidalVolume=0.0005                                                      "Base Tidal Volume";
    parameter Volume LungsAirVolume_initial = FunctionalResidualCapacity "Initial volume of alveolar space";
    parameter Volume pleuralVolume_initial = 0.0001 "Initial volume of pleural fluid";

    parameter Volume FunctionalResidualCapacity=0.00231                                                      "Functional residual capacity";

    parameter Physiolibrary.Types.HydraulicResistance TotalResistance=147099.75
      "Total lungs pathways resistance";

    parameter Real BronchiResistanceFraction = 0.3;
    parameter Real AlveoliDuctResistanceFraction = 0.2;
    parameter Real TracheaResistanceFraction = 1 - (BronchiResistanceFraction+AlveoliDuctResistanceFraction)/2;

    parameter Physiolibrary.Types.HydraulicResistance TracheaResistance=TotalResistance*TracheaResistanceFraction
      "Left Bronchi Resistance";
    parameter Physiolibrary.Types.HydraulicResistance LeftBronchiResistance=TotalResistance*BronchiResistanceFraction
      "Left Bronchi Resistance";
    parameter Physiolibrary.Types.HydraulicResistance LeftAlveoliResistance=TotalResistance*AlveoliDuctResistanceFraction
      "Left Alveoli Resistance";
      parameter Physiolibrary.Types.HydraulicResistance RightBronchiResistance=TotalResistance*BronchiResistanceFraction
      "Right Bronchi Resistance";
      parameter Physiolibrary.Types.HydraulicResistance RightAlveoliResistance=TotalResistance*AlveoliDuctResistanceFraction
      "Right Alveoli Resistance";

    parameter Physiolibrary.Types.HydraulicCompliance TotalCompliance(displayUnit=
         "l/cmH2O")=1.0197162129779e-06            "Total lungs compliance";

    parameter Pressure Pmin=-533.28954966                                   "Relative external lungs pressure minimum caused by respiratory muscles";
    parameter Pressure Pmax(displayUnit="mmHg")=0    "Relative external lungs pressure maximum";
    parameter Real RespiratoryMusclePressureCycle[:,3] = {
          {0,Pmax,0},
          {3/8,Pmin,0},
          {1,Pmax,0}}
            "Absolute external lungs pressure during respiration cycle scaled to time period (0,1)";

    parameter Temperature CoreTemperature=310.15                       "body temperature";
    parameter Temperature EnvironmentTemperature=298.15                       "external air temperature";

    parameter Mass m_initial = LungsAirVolume_initial*Air.density(Air.setState_pTX(system.p_ambient
           + Pmax, CoreTemperature));

  //  parameter Density d = Air.density(Air.setState_pTX(system.p_ambient+Pmax,CoreTemperature));

    Physiolibrary.Blocks.Source.PeriodicCurveSource respiratoryMusclePressureCycle(data=
          RespiratoryMusclePressureCycle)
      "Relative position in respiratory cycle (0,1) to absolute external lungs pressure"
      annotation (Placement(transformation(extent={{2,48},{22,68}})));

    Physiolibrary.Fluid.Components.ElasticVessel leftAlveoli(
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      ZeroPressureVolume=FunctionalResidualCapacity,
      ResidualVolume=ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=true,
      EnthalpyNotUsed=true,
      nPorts=2) "Left alveolar space"
      annotation (Placement(transformation(extent={{-162,16},{-142,36}})));                                     //0.0133,

    Physiolibrary.Fluid.Sensors.PressureMeasure leftPleauralPressure(redeclare
        package Medium = PleuralFluid, GetAbsolutePressure=true)
      "Left Pleaural pressure" annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=0,
          origin={-70,64})));

    inner Modelica.Fluid.System system(T_ambient=CoreTemperature)
                                       "Human body system setting"
      annotation (Placement(transformation(extent={{60,66},{80,86}})));

    Physiolibrary.Fluid.Sources.PressureSource environment(redeclare package
        Medium =         Air, T=EnvironmentTemperature)
      "External environment"
      annotation (Placement(transformation(extent={{-356,-10},{-336,10}})));

    Physiolibrary.Fluid.Components.Resistor leftBronchi(redeclare package Medium =
          Air,
      EnthalpyNotUsed=true,
      Resistance=LeftBronchiResistance + LeftAlveoliResistance)
      annotation (Placement(transformation(extent={{-252,24},{-232,44}})));

    Physiolibrary.Types.Constants.FrequencyConst frequency(k=RespirationRate)
      annotation (Placement(transformation(extent={{-18,54},{-10,62}})));
    Chemical.Sources.SubstanceOutflow O2_left(SubstanceFlow(displayUnit="mmol/min")=
           0.0001285)
      annotation (Placement(transformation(extent={{-166,-16},{-146,4}})));
    Chemical.Sources.SubstanceInflow CO2_left(SubstanceFlow(displayUnit="mmol/min")=
           0.00010283333333333)
      annotation (Placement(transformation(extent={{-218,-16},{-198,4}})));
    Physiolibrary.Fluid.Sensors.PressureMeasure leftAlveolarPressure(redeclare
        package Medium = Air) "Left Alveolar pressure"
      annotation (Placement(transformation(extent={{-124,22},{-104,42}})));
    Modelica.Blocks.Math.Add musclePressure annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={42,20})));
    Physiolibrary.Types.Constants.PressureConst ambient_pressure(k=
          IntrathoraxPressure)
      annotation (Placement(transformation(extent={{84,44},{74,52}})));
    Physiolibrary.Fluid.Components.Resistor rightBronchi(redeclare package
        Medium =
          Air,
      EnthalpyNotUsed=true,
      Resistance=RightBronchiResistance + RightAlveoliResistance)
      annotation (Placement(transformation(extent={{-252,-54},{-232,-34}})));
    Physiolibrary.Fluid.Components.ElasticVessel rightAlveoli(
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      ZeroPressureVolume=FunctionalResidualCapacity,
      ResidualVolume=ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=true,
      EnthalpyNotUsed=true,
      nPorts=2) "Right alveolar space"
      annotation (Placement(transformation(extent={{-156,-58},{-136,-38}})));
    Chemical.Sources.SubstanceInflow CO2_right(SubstanceFlow(displayUnit="mmol/min")=
           0.00010283333333333)
      annotation (Placement(transformation(extent={{-220,-96},{-200,-76}})));
    Chemical.Sources.SubstanceOutflow O2_right(SubstanceFlow(displayUnit="mmol/min")=
           0.0001285)
      annotation (Placement(transformation(extent={{-158,-96},{-138,-76}})));
    Physiolibrary.Fluid.Components.ElasticVessel leftPlearalSpace(
      redeclare package Medium = PleuralFluid,
      volume_start=pleuralVolume_initial,
      ZeroPressureVolume=pleuralVolume_initial,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=false,
      useSubstances=false,
      nPorts=1) "Left Plearal space"
      annotation (Placement(transformation(extent={{-76,18},{-56,38}})));
    Physiolibrary.Fluid.Sensors.PressureMeasure rightPleauralPressure(redeclare
        package Medium = PleuralFluid, GetAbsolutePressure=true)
      "Right pleaural pressure" annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=0,
          origin={-70,-12})));
    Physiolibrary.Fluid.Components.ElasticVessel rightPleuralSpace(
      redeclare package Medium = PleuralFluid,
      volume_start=pleuralVolume_initial,
      ZeroPressureVolume=pleuralVolume_initial,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=false,
      useSubstances=false,
      nPorts=1) "Right Plearal space"
      annotation (Placement(transformation(extent={{-76,-58},{-56,-38}})));
    Physiolibrary.Fluid.Sensors.PressureMeasure rightAlveolarPressure(redeclare
        package Medium = Air) "Right Alveolar pressure"
      annotation (Placement(transformation(extent={{-122,-54},{-102,-34}})));
  equation

    connect(frequency.y, respiratoryMusclePressureCycle.frequence)
      annotation (Line(points={{-9,58},{2,58}},                      color={0,0,127}));
    connect(leftAlveoli.substances[1], O2_left.port_a) annotation (Line(points={{-162,26},
            {-180,26},{-180,-6},{-166,-6}},     color={158,66,200}));
    connect(CO2_left.port_b, leftAlveoli.substances[2]) annotation (Line(points={{-198,-6},
            {-182,-6},{-182,26},{-162,26}},     color={158,66,200}));
    connect(ambient_pressure.y, musclePressure.u1)
      annotation (Line(points={{72.75,48},{48,48},{48,32}}, color={0,0,127}));
    connect(musclePressure.u2, respiratoryMusclePressureCycle.val)
      annotation (Line(points={{36,32},{36,58},{22,58}}, color={0,0,127}));
    connect(CO2_right.port_b, rightAlveoli.substances[2]) annotation (Line(points={{-200,
            -86},{-178,-86},{-178,-48},{-156,-48}},       color={158,66,200}));
    connect(O2_right.port_a, rightAlveoli.substances[1]) annotation (Line(points={{-158,
            -86},{-176,-86},{-176,-48},{-156,-48}},      color={158,66,200}));
    connect(leftAlveolarPressure.q_in, leftAlveoli.q_in[1]) annotation (Line(
        points={{-118,26},{-152,26},{-152,24},{-152.1,24},{-152.1,27.3}},
        color={127,0,0},
        thickness=0.5));
    connect(leftPlearalSpace.q_in[1], leftPleauralPressure.q_in) annotation (Line(
        points={{-66.1,28},{-66.1,48},{-66,48},{-66,58}},
        color={127,0,0},
        thickness=0.5));
    connect(leftPleauralPressure.pressure, leftAlveoli.externalPressure)
      annotation (Line(points={{-76,60},{-146,60},{-146,35}}, color={0,0,127}));
    connect(musclePressure.y, leftPlearalSpace.externalPressure) annotation (Line(
          points={{42,9},{42,-4},{-36,-4},{-36,46},{-60,46},{-60,37}}, color={0,0,127}));
    connect(rightPleuralSpace.q_in[1], rightPleauralPressure.q_in) annotation (
        Line(
        points={{-66.1,-48},{-66.1,-28},{-66,-28},{-66,-18}},
        color={127,0,0},
        thickness=0.5));
    connect(musclePressure.y,rightPleuralSpace.externalPressure)  annotation (Line(
          points={{42,9},{42,-4},{-36,-4},{-36,-32},{-60,-32},{-60,-39}}, color={0,
            0,127}));
    connect(rightAlveoli.externalPressure, rightPleauralPressure.pressure)
      annotation (Line(points={{-140,-39},{-140,-16},{-76,-16}}, color={0,0,127}));
    connect(rightAlveoli.q_in[1], rightAlveolarPressure.q_in) annotation (Line(
        points={{-146.1,-46.7},{-146.1,-48},{-146,-48},{-146,-50},{-116,-50}},
        color={127,0,0},
        thickness=0.5));
    connect(leftBronchi.q_out, leftAlveoli.q_in[2]) annotation (Line(
        points={{-232,34},{-152,34},{-152,30},{-152.1,30},{-152.1,24.7}},
        color={127,0,0},
        thickness=0.5));
    connect(rightBronchi.q_out, rightAlveoli.q_in[2]) annotation (Line(
        points={{-232,-44},{-148,-44},{-148,-40},{-146,-40},{-146,-49.3},{-146.1,
            -49.3}},
        color={127,0,0},
        thickness=0.5));
    connect(environment.y, leftBronchi.q_in) annotation (Line(
        points={{-336,0},{-258,0},{-258,34},{-252,34}},
        color={127,0,0},
        thickness=0.5));
    connect(environment.y, rightBronchi.q_in) annotation (Line(
        points={{-336,0},{-258,0},{-258,-44},{-252,-44}},
        color={127,0,0},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-360,-100},{100,100}})),
      experiment(StopTime=16, __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>References:</p>
<p><br>Mecklenburgh, J. S., and W. W. Mapleson. &quot;Ventilatory assistance and respiratory muscle activity. 1: Interaction in healthy volunteers.&quot; <i>British journal of anaesthesia</i> 80.4 (1998): 422-433.</p>
</html>"));
  end Respiration3;

  model PulseRespiration "Minimal respiration model"
    extends Modelica.Icons.Example;

    import Modelica.SIunits.*;

    replaceable package Air = Kitware.Air_IdealGas; //Chemical.Media.SimpleAir_C; //Chemical.Media.Air_MixtureGasNasa;
    replaceable package PleuralFluid = Chemical.Media.Water_Incompressible;

    parameter Boolean EnthalpyNotUsed = false;

    parameter Pressure IntrathoraxPressure = system.p_ambient - 700;
    parameter Frequency RespirationRate=0.133                                                                "Respiration rate";
    parameter Volume ResidualVolume=0.0013                                                      "Lungs residual volume";
    parameter Volume TotalLungCapacity=0.00623                                                      "Total Lung Capacity";
    parameter Volume BaseTidalVolume=0.0005                                                      "Base Tidal Volume";
    parameter Volume LungsAirVolume_initial = FunctionalResidualCapacity "Initial volume of alveolar space";
    parameter Volume pleuralVolume_initial = 0.0001 "Initial volume of pleural fluid";

    parameter Volume FunctionalResidualCapacity=0.00231                                                      "Functional residual capacity";


    parameter Physiolibrary.Types.HydraulicResistance TotalResistance=147099.75
      "Total lungs pathways resistance";


    parameter Real BronchiResistanceFraction = 0.3;
    parameter Real AlveoliDuctResistanceFraction = 0.2;
    parameter Real TracheaResistanceFraction = 1 - (BronchiResistanceFraction+AlveoliDuctResistanceFraction)/2;


    parameter Physiolibrary.Types.HydraulicResistance TracheaResistance=TotalResistance*TracheaResistanceFraction
      "Left Bronchi Resistance";
    parameter Physiolibrary.Types.HydraulicResistance LeftBronchiResistance=TotalResistance*BronchiResistanceFraction
      "Left Bronchi Resistance";
    parameter Physiolibrary.Types.HydraulicResistance LeftAlveoliResistance=TotalResistance*AlveoliDuctResistanceFraction
      "Left Alveoli Resistance";
      parameter Physiolibrary.Types.HydraulicResistance RightBronchiResistance=TotalResistance*BronchiResistanceFraction
      "Right Bronchi Resistance";
      parameter Physiolibrary.Types.HydraulicResistance RightAlveoliResistance=TotalResistance*AlveoliDuctResistanceFraction
      "Right Alveoli Resistance";

    parameter Physiolibrary.Types.HydraulicCompliance TotalCompliance(displayUnit=
         "l/cmH2O")=1.0197162129779e-06            "Total lungs compliance";

    parameter Pressure Pmin=-533.28954966                                   "Relative external lungs pressure minimum caused by respiratory muscles";
    parameter Pressure Pmax(displayUnit="mmHg")=0    "Relative external lungs pressure maximum";
    parameter Real RespiratoryMusclePressureCycle[:,3] = {
          {0,Pmax,0},
          {3/8,Pmin,0},
          {1,Pmax,0}}
            "Absolute external lungs pressure during respiration cycle scaled to time period (0,1)";

    parameter Temperature CoreTemperature=310.15                       "body temperature";
    parameter Temperature EnvironmentTemperature=298.15                       "external air temperature";

    parameter Mass m_initial = LungsAirVolume_initial*Air.density(Air.setState_pTX(system.p_ambient
           + Pmax, CoreTemperature,Air.reference_X));

  //  parameter Density d = Air.density(Air.setState_pTX(system.p_ambient+Pmax,CoreTemperature));

    Physiolibrary.Blocks.Source.PeriodicCurveSource respiratoryMusclePressureCycle(data=
          RespiratoryMusclePressureCycle)
      "Relative position in respiratory cycle (0,1) to absolute external lungs pressure"
      annotation (Placement(transformation(extent={{2,48},{22,68}})));



    Physiolibrary.Fluid.Components.ElasticVessel leftAlveoli(
      ConstantTemperature=true,
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      ZeroPressureVolume=FunctionalResidualCapacity,
      ResidualVolume=ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=true,
      EnthalpyNotUsed=EnthalpyNotUsed,
      nPorts=4) "Left alveolar space"
      annotation (Placement(transformation(extent={{-162,16},{-142,36}})));                                     //0.0133,

    Physiolibrary.Fluid.Sensors.PressureMeasure leftPleauralPressure(redeclare
        package Medium = PleuralFluid, GetAbsolutePressure=true)
      "Left Pleaural pressure" annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=0,
          origin={-70,64})));

    inner Modelica.Fluid.System system(T_ambient=CoreTemperature)
                                       "Human body system setting"
      annotation (Placement(transformation(extent={{60,66},{80,86}})));

    Physiolibrary.Fluid.Sources.PressureSource environment(redeclare package
        Medium =         Air, T=EnvironmentTemperature)
      "External environment"
      annotation (Placement(transformation(extent={{-356,-10},{-336,10}})));

    Physiolibrary.Fluid.Components.Resistor leftBronchi(redeclare package Medium =
          Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=LeftBronchiResistance)
      annotation (Placement(transformation(extent={{-252,24},{-232,44}})));

    Physiolibrary.Types.Constants.FrequencyConst frequency(k=RespirationRate)
      annotation (Placement(transformation(extent={{-18,54},{-10,62}})));
    Chemical.Sources.SubstanceOutflow O2_left(SubstanceFlow(displayUnit="mmol/min")=
           0.0001285)
      annotation (Placement(transformation(extent={{-166,-16},{-146,4}})));
    Chemical.Sources.SubstanceInflow CO2_left(SubstanceFlow(displayUnit="mmol/min")=
           0.00010283333333333)
      annotation (Placement(transformation(extent={{-218,-16},{-198,4}})));
    Physiolibrary.Fluid.Sensors.PressureMeasure leftAlveolarPressure(redeclare
        package Medium = Air) "Left Alveolar pressure"
      annotation (Placement(transformation(extent={{-124,22},{-104,42}})));
    Modelica.Blocks.Math.Add musclePressure annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={42,20})));
    Physiolibrary.Types.Constants.PressureConst ambient_pressure(k=
          IntrathoraxPressure)
      annotation (Placement(transformation(extent={{84,44},{74,52}})));
    Physiolibrary.Fluid.Components.Resistor rightBronchi(redeclare package Medium =
          Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=RightBronchiResistance)
      annotation (Placement(transformation(extent={{-252,-54},{-232,-34}})));
    Physiolibrary.Fluid.Components.ElasticVessel rightAlveoli(
      ConstantTemperature=true,
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      ZeroPressureVolume=FunctionalResidualCapacity,
      ResidualVolume=ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=true,
      EnthalpyNotUsed=EnthalpyNotUsed,
      nPorts=2) "Right alveolar space"
      annotation (Placement(transformation(extent={{-156,-58},{-136,-38}})));
    Chemical.Sources.SubstanceInflow CO2_right(SubstanceFlow(displayUnit="mmol/min")=
           0.00010283333333333)
      annotation (Placement(transformation(extent={{-220,-96},{-200,-76}})));
    Chemical.Sources.SubstanceOutflow O2_right(SubstanceFlow(displayUnit="mmol/min")=
           0.0001285)
      annotation (Placement(transformation(extent={{-158,-96},{-138,-76}})));
    Physiolibrary.Fluid.Components.ElasticVessel leftPlearalSpace(
      redeclare package Medium = PleuralFluid,
      volume_start=pleuralVolume_initial,
      ZeroPressureVolume=pleuralVolume_initial,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=false,
      useSubstances=false,
      nPorts=1) "Left Plearal space"
      annotation (Placement(transformation(extent={{-76,18},{-56,38}})));
    Physiolibrary.Fluid.Sensors.PressureMeasure rightPleauralPressure(redeclare
        package Medium = PleuralFluid, GetAbsolutePressure=true)
      "Right pleaural pressure" annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=0,
          origin={-70,-12})));
    Physiolibrary.Fluid.Components.ElasticVessel rightPleuralSpace(
      redeclare package Medium = PleuralFluid,
      volume_start=pleuralVolume_initial,
      ZeroPressureVolume=pleuralVolume_initial,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=false,
      useSubstances=false,
      nPorts=1) "Right Plearal space"
      annotation (Placement(transformation(extent={{-76,-58},{-56,-38}})));
    Physiolibrary.Fluid.Sensors.PressureMeasure rightAlveolarPressure(redeclare
        package Medium = Air) "Right Alveolar pressure"
      annotation (Placement(transformation(extent={{-122,-54},{-102,-34}})));
    Physiolibrary.Fluid.Components.Resistor trachea(
      redeclare package Medium = Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=TracheaResistance)
      annotation (Placement(transformation(extent={{-298,-10},{-278,10}})));
    Physiolibrary.Fluid.Components.Resistor leftAlveolarDuct(
      redeclare package Medium = Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=LeftAlveoliResistance)
      annotation (Placement(transformation(extent={{-218,24},{-198,44}})));
    Physiolibrary.Fluid.Sensors.FlowMeasure flowMeasure(redeclare package Medium =
          Air)
      annotation (Placement(transformation(extent={{-328,-10},{-308,10}})));
    Physiolibrary.Fluid.Components.Resistor rightAlveolarDuct(
      redeclare package Medium = Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=RightAlveoliResistance)
      annotation (Placement(transformation(extent={{-218,-54},{-198,-34}})));
    Physiolibrary.Fluid.Sensors.MassFractions O2_massFraction(redeclare package
        Medium = Air, substanceName="O2")
      annotation (Placement(transformation(extent={{-186,72},{-166,92}})));
    Physiolibrary.Fluid.Sensors.MassFractions CO2_massFraction(redeclare package
        Medium = Air, substanceName="CO2")
      annotation (Placement(transformation(extent={{-150,72},{-130,92}})));
  equation

    connect(frequency.y, respiratoryMusclePressureCycle.frequence)
      annotation (Line(points={{-9,58},{2,58}},                      color={0,0,127}));
    connect(leftAlveoli.substances[1], O2_left.port_a) annotation (Line(points={{-162,26},
            {-180,26},{-180,-6},{-166,-6}},     color={158,66,200}));
    connect(CO2_left.port_b, leftAlveoli.substances[2]) annotation (Line(points={{-198,-6},
            {-182,-6},{-182,26},{-162,26}},     color={158,66,200}));
    connect(ambient_pressure.y, musclePressure.u1)
      annotation (Line(points={{72.75,48},{48,48},{48,32}}, color={0,0,127}));
    connect(musclePressure.u2, respiratoryMusclePressureCycle.val)
      annotation (Line(points={{36,32},{36,58},{22,58}}, color={0,0,127}));
    connect(CO2_right.port_b, rightAlveoli.substances[2]) annotation (Line(points={{-200,
            -86},{-178,-86},{-178,-48},{-156,-48}},       color={158,66,200}));
    connect(O2_right.port_a, rightAlveoli.substances[1]) annotation (Line(points={{-158,
            -86},{-176,-86},{-176,-48},{-156,-48}},      color={158,66,200}));
    connect(leftAlveolarPressure.q_in, leftAlveoli.q_in[1]) annotation (Line(
        points={{-118,26},{-152,26},{-152,24},{-152.1,24},{-152.1,27.95}},
        color={127,0,0},
        thickness=0.5));
    connect(leftPlearalSpace.q_in[1], leftPleauralPressure.q_in) annotation (Line(
        points={{-66.1,28},{-66.1,48},{-66,48},{-66,58}},
        color={127,0,0},
        thickness=0.5));
    connect(leftPleauralPressure.pressure, leftAlveoli.externalPressure)
      annotation (Line(points={{-76,60},{-146,60},{-146,35}}, color={0,0,127}));
    connect(musclePressure.y, leftPlearalSpace.externalPressure) annotation (Line(
          points={{42,9},{42,-4},{-36,-4},{-36,46},{-60,46},{-60,37}}, color={0,0,127}));
    connect(rightPleuralSpace.q_in[1], rightPleauralPressure.q_in) annotation (
        Line(
        points={{-66.1,-48},{-66.1,-28},{-66,-28},{-66,-18}},
        color={127,0,0},
        thickness=0.5));
    connect(musclePressure.y,rightPleuralSpace.externalPressure)  annotation (Line(
          points={{42,9},{42,-4},{-36,-4},{-36,-32},{-60,-32},{-60,-39}}, color={0,
            0,127}));
    connect(rightAlveoli.externalPressure, rightPleauralPressure.pressure)
      annotation (Line(points={{-140,-39},{-140,-16},{-76,-16}}, color={0,0,127}));
    connect(rightAlveoli.q_in[1], rightAlveolarPressure.q_in) annotation (Line(
        points={{-146.1,-46.7},{-146.1,-48},{-146,-48},{-146,-50},{-116,-50}},
        color={127,0,0},
        thickness=0.5));

    connect(leftBronchi.q_out, leftAlveolarDuct.q_in) annotation (Line(
        points={{-232,34},{-218,34}},
        color={127,0,0},
        thickness=0.5));
    connect(leftAlveolarDuct.q_out, leftAlveoli.q_in[2]) annotation (Line(
        points={{-198,34},{-152,34},{-152,30},{-152.1,30},{-152.1,26.65}},
        color={127,0,0},
        thickness=0.5));
    connect(trachea.q_out, leftBronchi.q_in) annotation (Line(
        points={{-278,0},{-268,0},{-268,34},{-252,34}},
        color={127,0,0},
        thickness=0.5));
    connect(trachea.q_out, rightBronchi.q_in) annotation (Line(
        points={{-278,0},{-268,0},{-268,-44},{-252,-44}},
        color={127,0,0},
        thickness=0.5));
    connect(environment.y, flowMeasure.q_in) annotation (Line(
        points={{-336,0},{-328,0}},
        color={127,0,0},
        thickness=0.5));
    connect(flowMeasure.q_out, trachea.q_in) annotation (Line(
        points={{-308,0},{-298,0}},
        color={127,0,0},
        thickness=0.5));
    connect(rightBronchi.q_out, rightAlveolarDuct.q_in) annotation (Line(
        points={{-232,-44},{-218,-44}},
        color={127,0,0},
        thickness=0.5));
    connect(rightAlveolarDuct.q_out, rightAlveoli.q_in[2]) annotation (Line(
        points={{-198,-44},{-148,-44},{-148,-40},{-146,-40},{-146,-49.3},{-146.1,-49.3}},
        color={127,0,0},
        thickness=0.5));

    connect(O2_massFraction.port, leftAlveoli.q_in[3]) annotation (Line(points={{-176,
            72},{-176,34},{-152.1,34},{-152.1,25.35}}, color={0,127,255}));
    connect(CO2_massFraction.port, leftAlveoli.q_in[4]) annotation (Line(points={{-140,
            72},{-140,24.05},{-152.1,24.05}}, color={0,127,255}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-360,-100},{100,100}})),
      experiment(StopTime=200, __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>References:</p>
<p><br>Mecklenburgh, J. S., and W. W. Mapleson. &quot;Ventilatory assistance and respiratory muscle activity. 1: Interaction in healthy volunteers.&quot; <i>British journal of anaesthesia</i> 80.4 (1998): 422-433.</p>
</html>"));
  end PulseRespiration;

  model check_Medium
    import Modelica.SIunits.*;

    package Medium =
    Kitware.Air_IdealGas;
    //Chemical.Media.Air_MixtureGasNasa;
    //Chemical.Media.SimpleAir_C;
    //Chemical.Media.Water_Incompressible;
    //Chemical.Media.SimpleO2Gas_C;
   // Chemical.Media.EthanolInWater_C;
   // Chemical.Media.StandardWater_C;
   // Physiolibrary.Media.SimpleBodyFluid;

    //constants
    constant Integer nCS = Medium.nCS;
    constant Integer nX = Medium.nX;
    constant Integer nXi = Medium.nXi;
    constant Integer nC = Medium.nC;
    constant Medium.stateOfMatter.SubstanceData substanceData[Medium.nCS]=Medium.substanceData;

    //state variables
    Pressure p;
    Temperature T;
    MassFraction X[Medium.nX] "mass fraction", sumX = ones(Medium.nX)*X;

    MoleFraction x[Medium.nCS] = (X./molarMass) ./ (X*(ones(Medium.nCS)./molarMass)) "mole fraction", sumx=ones(Medium.nCS)*x;

    //functions
    Medium.ThermodynamicState state_pTX = Medium.setState_pTX(p,T,X);
    Medium.ThermodynamicState state_phX = Medium.setState_phX(p,specificEnthalpy,X);

    SpecificEnthalpy specificEnthalpy = Medium.specificEnthalpy(state_pTX);
    SpecificEntropy specificEntropy = Medium.specificEntropy(state_pTX);
    SpecificHeatCapacity specificHeatCapacityCp = Medium.specificHeatCapacityCp(state_pTX);
    SpecificHeatCapacity specificHeatCapacityCv = Medium.specificHeatCapacityCv(state_pTX);
    Density density=Medium.density(state_pTX);
    Temperature temperature = Medium.temperature(state_pTX);
    Pressure pressure = Medium.pressure(state_pTX);
    //MassFraction x_mass[Medium.nCS] = Medium.x_mass( Medium.Xi_outflow(X), Medium.C_outflow(X));
   // Concentration concentration[Medium.nCS] = Medium.concentration(state_pTX, Medium.Xi_outflow(X), Medium.C_outflow(X));

    //by susbstance
    MolarEnthalpy molarEnthalpy[Medium.nCS] = Medium.stateOfMatter.molarEnthalpy(substanceData,T=T,p=p);
    MolarMass molarMass[Medium.nCS] = Medium.stateOfMatter.molarMass(substanceData);
    SpecificEnthalpy _specificEnthalpy[Medium.nCS] = molarEnthalpy ./ molarMass;

    MolarEnthalpy _molarEnthalpy = x*molarEnthalpy;
    MolarEnthalpy _solution_h_base = x*Medium.stateOfMatter.molarEnthalpy(substanceData,298.15,p);
    MolarHeatCapacity _solution_Cp = sum(x[i]*substanceData[i].Cp for i in 1:size(x,1));
    MolarEnthalpy _h_diff;
    Temperature _solution_temperature,_T,_T2;

    MolarMass _MM= (ones(nCS)*X) / ((X./substanceData.MolarWeight) * ones(nCS));
    MoleFraction _x[nCS]=(X./substanceData.MolarWeight) ./ (X*(ones(nCS)./substanceData.MolarWeight));

    MolarEnthalpy _molarEnthalpy2;
  equation
    T=310.15;
    p=101325;
    X=Medium.reference_X;

    assert(abs(sumX-1) < 1e-2, "Mass fractions are inconsistent!");
    assert(abs(sumx-1) < 1e-2, "Mole fractions are inconsistent!");

    assert(abs(state_pTX.T-state_phX.T) < 1e-2, "Enthalpy-temperature relationship is inconsistent!");


    _h_diff = _molarEnthalpy - _solution_h_base;
    _solution_temperature = _h_diff/_solution_Cp +  298.15;

    _molarEnthalpy2 = specificEnthalpy*_MM;

    // h := state.X*(stateOfMatter.molarEnthalpy(substanceData,T=state.T,p=state.p) ./ substanceData.MolarWeight);

    //X*(molarEnthalpy ./ substanceData.MolarWeight);

    _T = Medium.stateOfMatter.solution_temperature(substanceData,_molarEnthalpy,x,p);
    _T2 = Medium.stateOfMatter.solution_temperature(substanceData,_molarEnthalpy2,x,p);


    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end check_Medium;
  annotation (uses(
      Modelica(version="3.2.3"),
      Chemical(version="1.4.0-alpha2"),
      Physiolibrary(version="3.0.0-alpha2")));
end Kitware;
