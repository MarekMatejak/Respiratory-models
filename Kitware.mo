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
    constant Modelica.Units.SI.MoleFraction x_default[nCS]={0.21,0.0004,
        0.02,0.7696}
      "Initial mole fractions of all substances (Please check: x_default*ones(nCS) = 1)";

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
      R_s = 8.3144/MM;
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
      Modelica.Units.SI.MolarEnergy u[nCS]
        "electro-chemical potential of substances in the solution";
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
      input Modelica.Units.SI.MassFraction x_mass[nCS];
     output Real C_outflow[nC];
    algorithm
      C_outflow := C_default;
      annotation(Inline=true);
    end C_outflow;

    function Xi_outflow
      input Modelica.Units.SI.MassFraction x_mass[nCS];
      output Modelica.Units.SI.MassFraction Xi[nXi];
    algorithm
      Xi := x_mass;
      annotation(Inline=true);
    end Xi_outflow;

    function x_mass
      input Modelica.Units.SI.MassFraction actualStream_Xi[nXi];
      input Real actualStream_C[nC];
      output Modelica.Units.SI.MassFraction x_mass[nCS];
    algorithm
      x_mass := actualStream_Xi;
      annotation(Inline=true);
    end x_mass;

    function concentration "Concentration of base substance molecules from Xi and C"
      input ThermodynamicState state;
      input Modelica.Units.SI.MassFraction Xi[nXi];
      input Real C[nC];
      output Modelica.Units.SI.Concentration concentration[nCS];
    algorithm
      concentration := (Xi./substanceData.MolarWeight)/(Xi*(stateOfMatter.molarVolume(substanceData,T=state.T,p=state.p)./substanceData.MolarWeight));
    end concentration;

    function specificEnthalpyOffsets "Difference between chemical substance enthalpy and medium substance enthalpy at temperature 298.15 K and 100kPa"
      input Modelica.Units.SI.ElectricPotential v=0;
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

    import Modelica.Units.SI.*;

   replaceable package Air = Physiolibrary.Media.Air;

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

    parameter Boolean EnthalpyNotUsed = false;

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
      EnthalpyNotUsed=EnthalpyNotUsed,
      nPorts=2) "Lungs"
      annotation (Placement(transformation(extent={{36,-28},{56,-8}})));

    import Modelica.Units.SI.*;
    replaceable package Air = Physiolibrary.Media.Air;//Kitware.Air_IdealGas;
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
    Modelica.Blocks.Sources.ContinuousClock clock(offset=-4000)
      annotation (Placement(transformation(extent={{-38,50},{-18,70}})));
    Physiolibrary.Fluid.Sources.PressureSource environment(redeclare package
        Medium = Air, T=EnvironmentTemperature)
      "External environment"
      annotation (Placement(transformation(extent={{-76,-30},{-56,-10}})));
    Physiolibrary.Fluid.Sensors.FlowMeasure airflowMeasure(redeclare package
        Medium = Air, EnthalpyNotUsed=EnthalpyNotUsed)
                              "Lungs pathway airflow"
      annotation (Placement(transformation(extent={{-40,-30},{-20,-10}})));
    Physiolibrary.Fluid.Components.Resistor resistor(redeclare package Medium =
          Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance(displayUnit="(mmHg.min)/l") = 79993432.449)
      annotation (Placement(transformation(extent={{-6,-30},{14,-10}})));
    Physiolibrary.Fluid.Sensors.Temperature temperature(redeclare package
        Medium = Air)
      annotation (Placement(transformation(extent={{16,-4},{36,16}})));
    Physiolibrary.Fluid.Sensors.Temperature temperature1(redeclare package
        Medium = Air)
      annotation (Placement(transformation(extent={{-60,10},{-40,30}})));
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
        points={{14,-20},{14,-16.7},{45.9,-16.7}},
        color={127,0,0},
        thickness=0.5));
    connect(temperature.port, lungs.q_in[2]) annotation (Line(points={{26,-4},
            {26,-19.3},{45.9,-19.3}}, color={0,127,255}));
    connect(environment.y, temperature1.port) annotation (Line(
        points={{-56,-20},{-50,-20},{-50,10}},
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

  model LungsComplianceWithCooling "Respiration model"
    extends Modelica.Icons.Example;

    parameter Boolean EnthalpyNotUsed = false;

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
      EnthalpyNotUsed=EnthalpyNotUsed,
      nPorts=3) "Lungs"
      annotation (Placement(transformation(extent={{36,-28},{56,-8}})));

    import Modelica.Units.SI.*;
    replaceable package Air = Physiolibrary.Media.Air;// Physiolibrary.Media.Water_Incompressible;
                                                                  // Chemical.Media.SimpleAir_C;
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
    Modelica.Blocks.Sources.ContinuousClock clock(offset=-4000)
      annotation (Placement(transformation(extent={{-38,50},{-18,70}})));
    Physiolibrary.Fluid.Sources.PressureSource environment(redeclare package
        Medium = Air,
      usePressureInput=true,
                      T=EnvironmentTemperature)
      "External environment"
      annotation (Placement(transformation(extent={{-76,-30},{-56,-10}})));
    Physiolibrary.Fluid.Components.Conductor conductor(redeclare package
        Medium =
          Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Conductance(displayUnit="l/(mmHg.min)") = 1.2501026264094e-08)
      annotation (Placement(transformation(extent={{-24,-30},{-4,-10}})));
    Physiolibrary.Fluid.Sensors.Temperature temperature(redeclare package
        Medium = Air)
      annotation (Placement(transformation(extent={{-56,10},{-36,30}})));
    Physiolibrary.Fluid.Sensors.Temperature temperature1(redeclare package
        Medium = Air)
      annotation (Placement(transformation(extent={{8,12},{28,32}})));
  equation
     pressue = clock.y;
    connect(lungsPressureMeasure.q_in, lungs.q_in[1]) annotation (Line(
        points={{76,-16},{76,-16.2667},{45.9,-16.2667}},
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
        points={{-4,-20},{10,-20},{10,-18},{45.9,-18}},
        color={127,0,0},
        thickness=0.5));
    connect(temperature.port, environment.y) annotation (Line(points={{-46,10},
            {-46,-20},{-56,-20}}, color={0,127,255}));
    connect(temperature1.port, lungs.q_in[3]) annotation (Line(points={{18,12},
            {18,-19.3},{45.9,-19.3},{45.9,-19.7333}}, color={0,127,255}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}})),                                        Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
      experiment(StopTime=8000, __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>References:</p>
<p><br>Mecklenburgh, J. S., and W. W. Mapleson. &quot;Ventilatory assistance and respiratory muscle activity. 1: Interaction in healthy volunteers.&quot; <i>British journal of anaesthesia</i> 80.4 (1998): 422-433.</p>
</html>"));
  end LungsComplianceWithCooling;

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
      nPorts=3) "Lungs"
      annotation (Placement(transformation(extent={{36,-28},{56,-8}})));

    import Modelica.Units.SI.*;

    replaceable package Air = Physiolibrary.Media.Air;// Physiolibrary.Media.Water_Incompressible;
                                                                    //Chemical.Media.Air_MixtureGasNasa;

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

    Physiolibrary.Fluid.Components.Resistor resistor(redeclare package
        Medium =
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
    Physiolibrary.Fluid.Sensors.Temperature temperature(redeclare package
        Medium = Air)
      annotation (Placement(transformation(extent={{-56,16},{-36,36}})));
    Physiolibrary.Fluid.Sensors.Temperature temperature1(redeclare package
        Medium = Air)
      annotation (Placement(transformation(extent={{8,16},{28,36}})));
  equation
     transMusclePressue = respiratoryMusclePressureCycle.val;
    connect(lungsPressureMeasure.q_in, lungs.q_in[1]) annotation (Line(
        points={{76,-16},{76,-16.2667},{45.9,-16.2667}},
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
        points={{14,-20},{28,-20},{28,-18},{45.9,-18}},
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
    connect(temperature.port, airflowMeasure.q_in) annotation (Line(points={{
            -46,16},{-48,16},{-48,-20},{-40,-20}}, color={0,127,255}));
    connect(temperature1.port, lungs.q_in[3]) annotation (Line(points={{18,16},
            {18,-19.7333},{45.9,-19.7333}},
                                  color={0,127,255}));
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

    import Modelica.Units.SI.*;

    replaceable package Air = Physiolibrary.Media.Air;//Chemical.Media.SimpleAir_C; //Chemical.Media.Air_MixtureGasNasa;

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
      nPorts=3) "Lungs"
      annotation (Placement(transformation(extent={{36,-28},{56,-8}})));                                        //0.0133,

    Physiolibrary.Fluid.Sensors.PressureMeasure lungsPressureMeasure(
        redeclare package Medium = Air) "Lungs pressure"
      annotation (Placement(transformation(extent={{70,-20},{90,0}})));

    inner Modelica.Fluid.System system(T_ambient=CoreTemperature)
                                       "Human body system setting"
      annotation (Placement(transformation(extent={{60,66},{80,86}})));

    Physiolibrary.Fluid.Sources.PressureSource environment(redeclare package
                Medium = Air, T=EnvironmentTemperature)
      "External environment"
      annotation (Placement(transformation(extent={{-78,-32},{-58,-12}})));

    Physiolibrary.Fluid.Sensors.FlowMeasure airflowMeasure(redeclare package
                Medium = Air) "Lungs pathway airflow"
      annotation (Placement(transformation(extent={{-42,-32},{-20,-12}})));

    Physiolibrary.Fluid.Components.Resistor resistor(redeclare package
        Medium =
          Air, Resistance=TotalResistance)
      annotation (Placement(transformation(extent={{-6,-32},{14,-12}})));

    Physiolibrary.Types.Constants.FrequencyConst frequency(k=RespirationRate)
      annotation (Placement(transformation(extent={{-46,62},{-38,70}})));
    Chemical.Sources.SubstanceOutflow O2(SubstanceFlow(displayUnit="mmol/min") = 0.000257)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=180,
          origin={10,-44})));
    Chemical.Sources.SubstanceInflowT CO2(
                                         SubstanceFlow(displayUnit="mmol/min") = 0.00020566666666667,
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.CarbonDioxide_gas())
      annotation (Placement(transformation(extent={{0,-70},{20,-50}})));
    Physiolibrary.Fluid.Sensors.Temperature temperature(redeclare package
        Medium = Air)
      annotation (Placement(transformation(extent={{-58,10},{-38,32}})));
    Physiolibrary.Fluid.Sensors.Temperature temperature1(redeclare package
        Medium = Air)
      annotation (Placement(transformation(extent={{12,12},{32,32}})));
  equation

    connect(lungsPressureMeasure.q_in, lungs.q_in[1]) annotation (Line(
        points={{76,-16},{76,-16.2667},{45.9,-16.2667}},
        color={127,0,0},
        thickness=0.5));
    connect(resistor.q_out, lungs.q_in[2]) annotation (Line(
        points={{14,-22},{28,-22},{28,-18},{45.9,-18}},
        color={127,0,0},
        thickness=0.5));
    connect(frequency.y, respiratoryMusclePressureCycle.frequence)
      annotation (Line(points={{-37,66},{-8,66},{-8,64},{18,64}},    color={0,0,127}));
    connect(respiratoryMusclePressureCycle.val, lungs.externalPressure)
      annotation (Line(points={{38,64},{52,64},{52,-9}},color={0,0,127}));
    connect(CO2.port_b, lungs.substances[2]) annotation (Line(points={{20,-60},
            {32,-60},{32,-18},{36,-18}},
                                       color={158,66,200}));
    connect(airflowMeasure.q_in, environment.y) annotation (Line(
        points={{-42,-22},{-58,-22}},
        color={127,0,0},
        thickness=0.5));
    connect(airflowMeasure.q_out, resistor.q_in) annotation (Line(
        points={{-20,-22},{-6,-22}},
        color={127,0,0},
        thickness=0.5));
    connect(temperature.port, environment.y) annotation (Line(points={{-48,10},
            {-48,-22},{-58,-22}}, color={0,127,255}));
    connect(temperature1.port, lungs.q_in[3]) annotation (Line(points={{22,12},
            {22,-22},{28,-22},{28,-19.7333},{45.9,-19.7333}}, color={0,127,255}));
    connect(O2.port_a, lungs.substances[1]) annotation (Line(points={{20,
            -44},{30,-44},{30,-18},{36,-18}}, color={158,66,200}));
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

   import Modelica.Units.SI.*;

   replaceable package Air = Physiolibrary.Media.Air;//Chemical.Media.SimpleAir_C; //Chemical.Media.Air_MixtureGasNasa;
   replaceable package PleuralFluid =
       Physiolibrary.Media.Water;


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
   parameter Temperature EnvironmentTemperature=298.15                       "external air temperature";

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
     nPorts=6) "Alveolar space of lungs"
     annotation (Placement(transformation(extent={{48,-26},{68,-6}})));                                        //0.0133,

   Physiolibrary.Fluid.Sensors.PressureMeasure alveolarPressure(redeclare package
                Medium =Air) "Alveolar pressure"
     annotation (Placement(transformation(extent={{80,-20},{100,0}})));

   inner Modelica.Fluid.System system(T_ambient=CoreTemperature)
                                      "Human body system setting"
     annotation (Placement(transformation(extent={{60,66},{80,86}})));

   Physiolibrary.Fluid.Sources.PressureSource environment(redeclare package
        Medium =        Air, T=EnvironmentTemperature)
     "External environment"
     annotation (Placement(transformation(extent={{-98,-28},{-78,-8}})));

   Physiolibrary.Fluid.Sensors.FlowMeasure airflowMeasure(redeclare package
        Medium =        Air) "Lungs pathway airflow"
     annotation (Placement(transformation(extent={{-72,-28},{-50,-8}})));

   Physiolibrary.Fluid.Components.Resistor resistor(redeclare package Medium =
         Air, Resistance=TotalResistance)
     annotation (Placement(transformation(extent={{-42,-28},{-22,-8}})));

   Physiolibrary.Types.Constants.FrequencyConst frequency(k=RespirationRate)
     annotation (Placement(transformation(extent={{-96,86},{-88,94}})));
   Physiolibrary.Fluid.Sensors.PressureMeasure pleuralPressure(redeclare package
                Medium =
                PleuralFluid, GetAbsolutePressure=true) "Pleural pressure"
     annotation (Placement(transformation(extent={{-2,42},{18,62}})));
   Modelica.Blocks.Math.Add add annotation (Placement(transformation(
         extent={{-10,-10},{10,10}},
         rotation=270,
         origin={-4,72})));
   Physiolibrary.Types.Constants.PressureConst ambient_pressure(k=
         IntrathoraxPressure)
     annotation (Placement(transformation(extent={{30,88},{20,96}})));



   Physiolibrary.Fluid.Sensors.FlowMeasure AlvelarFlowMeasure(redeclare package
                Medium =
                Air)
     annotation (Placement(transformation(extent={{-8,-30},{14,-10}})));
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
    Physiolibrary.Fluid.Sensors.PartialPressure pCO2_alveolar(
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.CarbonDioxide_gas(),
      redeclare package Medium = Air)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=180,
          origin={66,-74})));
    Physiolibrary.Fluid.Sensors.PartialPressure pO2_alveolar(
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.Oxygen_gas(),
      redeclare package Medium = Air)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=180,
          origin={52,-46})));
    Chemical.Sources.SubstanceOutflow O2(
      SubstanceFlow(displayUnit="mmol/min") = 0.000257)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=180,
          origin={22,-46})));
    Chemical.Sources.SubstanceInflowT CO2(
      SubstanceFlow(displayUnit="mmol/min") = 0.00020566666666667,
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.CarbonDioxide_gas())
      annotation (Placement(transformation(extent={{22,-84},{42,-64}})));
    Physiolibrary.Fluid.Sensors.Temperature Temperature_alveolar(redeclare package
                Medium = Air)
      annotation (Placement(transformation(extent={{14,10},{34,30}})));
    Physiolibrary.Fluid.Sensors.PartialPressure pH2O_alveolar(
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.Water_gas(),
      redeclare package Medium = Air) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={82,-92})));
    Physiolibrary.Fluid.Sensors.Temperature Temperature_mouth(redeclare package
                Medium = Air)
      annotation (Placement(transformation(extent={{-70,8},{-50,28}})));
  equation



   connect(alveolarPressure.q_in,alveoli. q_in[1]) annotation (Line(
       points={{86,-16},{86,-13.8333},{57.9,-13.8333}},
       color={127,0,0},
       thickness=0.5));
   connect(frequency.y, respiratoryMusclePressureCycle.frequence)
     annotation (Line(points={{-87,90},{-62,90}},                   color={0,0,127}));
   connect(alveoli.substances[1], O2.port_a) annotation (Line(points={{48,-16},
            {38,-16},{38,-46},{32,-46}},                  color={158,66,200}));
   connect(CO2.port_b,alveoli. substances[2]) annotation (Line(points={{42,-74},
            {50,-74},{50,-16},{48,-16}}, color={158,66,200}));
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
       points={{14,-20},{34,-20},{34,-14},{40,-14},{40,-14.7},{57.9,-14.7}},
       color={127,0,0},
       thickness=0.5));
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
    connect(pCO2_alveolar.referenceFluidPort, alveoli.q_in[3]) annotation (
        Line(
        points={{66,-64.2},{66,-15.5667},{57.9,-15.5667}},
        color={127,0,0},
        thickness=0.5));
    connect(pCO2_alveolar.port_a, CO2.port_b) annotation (Line(points={{56,-74},
            {42,-74}},                   color={158,66,200}));
    connect(pO2_alveolar.referenceFluidPort, alveoli.q_in[4]) annotation (
        Line(
        points={{52,-36.2},{52,-16.4333},{57.9,-16.4333}},
        color={127,0,0},
        thickness=0.5));
    connect(pO2_alveolar.port_a, O2.port_a) annotation (Line(points={{42,-46},
            {32,-46}},                            color={158,66,200}));
    connect(Temperature_alveolar.port, alveoli.q_in[5]) annotation (Line(
          points={{24,10},{26,10},{26,-20},{34,-20},{34,-14},{40,-14},{40,-15.35},
            {57.9,-15.35},{57.9,-17.3}}, color={0,127,255}));
    connect(pH2O_alveolar.referenceFluidPort, alveoli.q_in[6]) annotation (
        Line(
        points={{82,-82.2},{82,-32},{66,-32},{66,-18.1667},{57.9,-18.1667}},
        color={127,0,0},
        thickness=0.5));

    connect(pH2O_alveolar.port_a, alveoli.substances[3]) annotation (Line(
          points={{72,-92},{48,-92},{48,-18},{56,-18},{56,-16},{48,-16}},
          color={158,66,200}));
    connect(Temperature_mouth.port, airflowMeasure.q_in) annotation (Line(
          points={{-60,8},{-60,0},{-72,0},{-72,-18}}, color={0,127,255}));
   annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
             {100,100}})),                                        Diagram(
         coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
     experiment(
        StopTime=500,
        __Dymola_NumberOfIntervals=5000,
        Tolerance=1e-05,
        __Dymola_Algorithm="Dassl"),
     Documentation(info="<html>
<p>References:</p>
<p><br>Mecklenburgh, J. S., and W. W. Mapleson. &quot;Ventilatory assistance and respiratory muscle activity. 1: Interaction in healthy volunteers.&quot; <i>British journal of anaesthesia</i> 80.4 (1998): 422-433.</p>
</html>"));
  end Respiration2;

  model Respiration3 "Respiration model"
    extends Modelica.Icons.Example;

    import Modelica.Units.SI.*;

    replaceable package Air = Physiolibrary.Media.Air; //Chemical.Media.Air_MixtureGasNasa;
    replaceable package PleuralFluid =
        Physiolibrary.Media.Water;

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
      EnthalpyNotUsed=EnthalpyNotUsed,
      nPorts=3) "Left alveolar space"
      annotation (Placement(transformation(extent={{-162,16},{-142,36}})));                                     //0.0133,

    Physiolibrary.Fluid.Sensors.PressureMeasure leftPleauralPressure(redeclare package
                Medium = PleuralFluid, GetAbsolutePressure=true)
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
      Resistance=LeftBronchiResistance + LeftAlveoliResistance)
      annotation (Placement(transformation(extent={{-252,24},{-232,44}})));

    Physiolibrary.Types.Constants.FrequencyConst frequency(k=RespirationRate)
      annotation (Placement(transformation(extent={{-18,54},{-10,62}})));
    Chemical.Sources.SubstanceOutflow O2_left(SubstanceFlow(displayUnit="mmol/min")=
           0.0001285)
      annotation (Placement(transformation(extent={{-166,-16},{-146,4}})));
    Chemical.Sources.SubstanceInflowT CO2_left(
                                              SubstanceFlow(displayUnit="mmol/min")=
           0.00010283333333333,
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.CarbonDioxide_gas())
      annotation (Placement(transformation(extent={{-218,-16},{-198,4}})));
    Physiolibrary.Fluid.Sensors.PressureMeasure leftAlveolarPressure(redeclare package
                Medium = Air) "Left Alveolar pressure"
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
      EnthalpyNotUsed=EnthalpyNotUsed,
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
      EnthalpyNotUsed=EnthalpyNotUsed,
      nPorts=4) "Right alveolar space"
      annotation (Placement(transformation(extent={{-156,-58},{-136,-38}})));
    Chemical.Sources.SubstanceInflowT CO2_right(
                                               SubstanceFlow(displayUnit="mmol/min")=
           0.00010283333333333,
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.CarbonDioxide_gas())
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
    Physiolibrary.Fluid.Sensors.Temperature alveolarTemperature(redeclare package
                Medium = Air)
      annotation (Placement(transformation(extent={{-174,58},{-154,78}})));
    Physiolibrary.Fluid.Sensors.PartialPressure pCO2(
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.CarbonDioxide_gas(),
      redeclare package Medium = Air)
      annotation (Placement(transformation(extent={{-218,-64},{-198,-84}})));
    Physiolibrary.Fluid.Sensors.PartialPressure pO2(
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.Oxygen_gas(),
      redeclare package Medium = Air)
      annotation (Placement(transformation(extent={{-138,-64},{-158,-84}})));
    Physiolibrary.Fluid.Sensors.Temperature Temperature_mount(redeclare package
                Medium = Air)
      annotation (Placement(transformation(extent={{-316,40},{-296,60}})));
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
        points={{-118,26},{-152,26},{-152,24},{-152.1,24},{-152.1,27.7333}},
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
        points={{-146.1,-46.05},{-146.1,-48},{-146,-48},{-146,-50},{-116,-50}},
        color={127,0,0},
        thickness=0.5));
    connect(leftBronchi.q_out, leftAlveoli.q_in[2]) annotation (Line(
        points={{-232,34},{-152,34},{-152,30},{-152.1,30},{-152.1,26}},
        color={127,0,0},
        thickness=0.5));
    connect(rightBronchi.q_out, rightAlveoli.q_in[2]) annotation (Line(
        points={{-232,-44},{-148,-44},{-148,-40},{-146,-40},{-146,-47.35},{-146.1,
            -47.35}},
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
    connect(alveolarTemperature.port, leftAlveoli.q_in[3]) annotation (Line(
          points={{-164,58},{-164,34},{-152.1,34},{-152.1,24.2667}}, color={0,
            127,255}));
    connect(pCO2.referenceFluidPort, rightAlveoli.q_in[3]) annotation (Line(
        points={{-208,-64.2},{-208,-44},{-148,-44},{-148,-40},{-146,-40},{-146,
            -49.3},{-146.1,-49.3},{-146.1,-48.65}},
        color={127,0,0},
        thickness=0.5));
    connect(CO2_right.port_b, pCO2.port_a) annotation (Line(points={{-200,-86},
            {-192,-86},{-192,-74},{-198,-74}}, color={158,66,200}));
    connect(O2_right.port_a, pO2.port_a) annotation (Line(points={{-158,-86},
            {-164,-86},{-164,-74},{-158,-74}}, color={158,66,200}));
    connect(pO2.referenceFluidPort, rightAlveoli.q_in[4]) annotation (Line(
        points={{-148,-64.2},{-146,-64.2},{-146,-44},{-146.1,-44},{-146.1,-49.95}},
        color={127,0,0},
        thickness=0.5));

    connect(Temperature_mount.port, leftBronchi.q_in) annotation (Line(
          points={{-306,40},{-308,40},{-308,0},{-258,0},{-258,34},{-252,34}},
          color={0,127,255}));
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

    import Modelica.Units.SI.*;

    replaceable package Air = Physiolibrary.Media.Air; //Chemical.Media.SimpleAir_C; //Kitware.Air_IdealGas; //Chemical.Media.SimpleAir_C; //Chemical.Media.Air_MixtureGasNasa;
    replaceable package PleuralFluid =
        Physiolibrary.Media.Water;

    parameter Boolean EnthalpyNotUsed = false;

    parameter Pressure IntrathoraxPressure = system.p_ambient - 700;
    parameter Frequency RespirationRate=0.08                                                               "Respiration rate";
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
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      useThermalPort=true,
      ZeroPressureVolume=FunctionalResidualCapacity,
      ResidualVolume=ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=true,
      EnthalpyNotUsed=EnthalpyNotUsed,
      nPorts=2) "Left alveolar space"
      annotation (Placement(transformation(extent={{-162,16},{-142,36}})));                                     //0.0133,

    Physiolibrary.Fluid.Sensors.PressureMeasure leftPleauralPressure(redeclare package
                Medium = PleuralFluid, GetAbsolutePressure=true)
      "Left Pleaural pressure" annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=0,
          origin={-70,64})));

    inner Modelica.Fluid.System system(T_ambient=CoreTemperature)
                                       "Human body system setting"
      annotation (Placement(transformation(extent={{60,66},{80,86}})));

    Physiolibrary.Fluid.Sources.PressureSource environment(redeclare package
                Medium = Air, T=EnvironmentTemperature)
      "External environment"
      annotation (Placement(transformation(extent={{-360,78},{-340,98}})));

    Physiolibrary.Fluid.Components.Resistor leftBronchi(redeclare package
        Medium =
          Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=LeftBronchiResistance)
      annotation (Placement(transformation(extent={{-252,24},{-232,44}})));

    Physiolibrary.Types.Constants.FrequencyConst frequency(k=RespirationRate)
      annotation (Placement(transformation(extent={{-18,54},{-10,62}})));
    Chemical.Sources.SubstanceOutflow O2_left(SubstanceFlow(displayUnit="mmol/min")=
           0.0001285)
      annotation (Placement(transformation(extent={{-166,-16},{-146,4}})));
    Chemical.Sources.SubstanceInflowT CO2_left(
                                              SubstanceFlow(displayUnit="mmol/min")=
           0.00010283333333333,
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.CarbonDioxide_gas())
      annotation (Placement(transformation(extent={{-218,-16},{-198,4}})));
    Physiolibrary.Fluid.Sensors.PressureMeasure leftAlveolarPressure(redeclare package
                Medium = Air) "Left Alveolar pressure"
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
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=RightBronchiResistance)
      annotation (Placement(transformation(extent={{-252,-54},{-232,-34}})));
    Physiolibrary.Fluid.Components.ElasticVessel rightAlveoli(
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      useThermalPort=true,
      ZeroPressureVolume=FunctionalResidualCapacity,
      ResidualVolume=ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=true,
      EnthalpyNotUsed=EnthalpyNotUsed,
      nPorts=6) "Right alveolar space"
      annotation (Placement(transformation(extent={{-156,-58},{-136,-38}})));
    Chemical.Sources.SubstanceInflowT CO2_right(
                                               SubstanceFlow(displayUnit="mmol/min")=
           0.00010283333333333,
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.CarbonDioxide_gas())
      annotation (Placement(transformation(extent={{-220,-96},{-200,-76}})));
    Chemical.Sources.SubstanceOutflow O2_right(SubstanceFlow(displayUnit="mmol/min")=
           0.0001285)
      annotation (Placement(transformation(extent={{-158,-96},{-138,-76}})));
    Physiolibrary.Fluid.Components.ElasticVessel leftPlearalSpace(
      redeclare package Medium = PleuralFluid,
      volume_start=pleuralVolume_initial,
      useThermalPort=false,
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
      annotation (Placement(transformation(extent={{-134,-38},{-114,-18}})));
    Physiolibrary.Fluid.Components.Resistor trachea(
      redeclare package Medium = Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=0.5*TracheaResistance)
      annotation (Placement(transformation(extent={{-298,-10},{-278,10}})));
    Physiolibrary.Fluid.Components.Resistor leftAlveolarDuct(
      redeclare package Medium = Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=LeftAlveoliResistance)
      annotation (Placement(transformation(extent={{-210,24},{-190,44}})));
    Physiolibrary.Fluid.Sensors.FlowMeasure flowMeasure(redeclare package
        Medium =
          Air)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-318,66})));
    Physiolibrary.Fluid.Components.Resistor rightAlveolarDuct(
      redeclare package Medium = Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=RightAlveoliResistance)
      annotation (Placement(transformation(extent={{-212,-54},{-192,-34}})));
    Physiolibrary.Fluid.Components.ElasticVessel upperRespiratoryTract(
      redeclare package Medium = Air,
      useSubstances=true,
      volume_start=0.0001,
      useThermalPort=true,
      Compliance=TotalCompliance/100,
      ZeroPressureVolume(displayUnit="ml") = 0.0001,
      ExternalPressure=system.p_ambient,
      ResidualVolume(displayUnit="ml") = 0.0001,
      nPorts=4)
      annotation (Placement(transformation(extent={{-328,-10},{-308,10}})));
    Physiolibrary.Fluid.Components.Resistor upperRespiratoryTractResistance(
      redeclare package Medium = Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=0.5*TracheaResistance) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-318,32})));
    Chemical.Sources.PureSubstance water(redeclare package stateOfMatter =
          Chemical.Interfaces.Incompressible, substanceData=
          Chemical.Substances.Water_liquid()) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-314,-68})));
    Chemical.Components.GasSolubility gasSolubility1(KC=1e-7) annotation (
        Placement(transformation(extent={{-362,-48},{-342,-28}})));
    Physiolibrary.Fluid.Sensors.PartialPressure pCO2(
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.CarbonDioxide_gas(),
      redeclare package Medium = Air)
      annotation (Placement(transformation(extent={{-218,-62},{-198,-82}})));
    Physiolibrary.Fluid.Sensors.PartialPressure pO2(
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.Oxygen_gas(),
      redeclare package Medium = Air)
      annotation (Placement(transformation(extent={{-138,-64},{-158,-84}})));
    Physiolibrary.Fluid.Sensors.Temperature Temperature_alveolar(redeclare package
                Medium = Air)
      annotation (Placement(transformation(extent={{-110,-40},{-90,-20}})));
    Physiolibrary.Fluid.Sensors.PartialPressure pH2O_alveolar(
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.Water_gas(),
      redeclare package Medium = Air)
      annotation (Placement(transformation(extent={{-122,-66},{-102,-86}})));
    Physiolibrary.Fluid.Sensors.PartialPressure pH2O_upperRespiratory(
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.Water_gas(),
      redeclare package Medium = Air)
      annotation (Placement(transformation(extent={{-364,34},{-344,14}})));
    Physiolibrary.Fluid.Sensors.Temperature Temperature_upperRespiratory(
        redeclare package Medium = Air)
      annotation (Placement(transformation(extent={{-298,30},{-278,50}})));
    Physiolibrary.Fluid.Sensors.Temperature Temperature_mouth(redeclare package
                Medium = Air)
      annotation (Placement(transformation(extent={{-296,72},{-276,92}})));
    Physiolibrary.Thermal.Components.Conductor conductor(Conductance(
          displayUnit="W/K") = 10)
      annotation (Placement(transformation(extent={{-302,-44},{-322,-24}})));
    Physiolibrary.Thermal.Sources.UnlimitedHeat coreHeat(T=system.T_ambient)
      annotation (Placement(transformation(extent={{-274,-44},{-294,-24}})));
    Physiolibrary.Thermal.Components.Conductor conductor1(Conductance(
          displayUnit="W/K") = 10)
      annotation (Placement(transformation(extent={{-210,6},{-230,26}})));
    Physiolibrary.Thermal.Components.Conductor conductor2(Conductance(
          displayUnit="W/K") = 10)
      annotation (Placement(transformation(extent={{-208,-40},{-228,-20}})));
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
        points={{-146.1,-45.8333},{-148,-45.8333},{-148,-60},{-128,-60},{
            -128,-34}},
        color={127,0,0},
        thickness=0.5));

    connect(leftBronchi.q_out, leftAlveolarDuct.q_in) annotation (Line(
        points={{-232,34},{-210,34}},
        color={127,0,0},
        thickness=0.5));
    connect(leftAlveolarDuct.q_out, leftAlveoli.q_in[2]) annotation (Line(
        points={{-190,34},{-152,34},{-152,30},{-152.1,30},{-152.1,24.7}},
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
        points={{-340,88},{-318,88},{-318,76}},
        color={127,0,0},
        thickness=0.5));
    connect(rightBronchi.q_out, rightAlveolarDuct.q_in) annotation (Line(
        points={{-232,-44},{-212,-44}},
        color={127,0,0},
        thickness=0.5));
    connect(rightAlveolarDuct.q_out, rightAlveoli.q_in[2]) annotation (Line(
        points={{-192,-44},{-146,-44},{-146,-46.7},{-146.1,-46.7}},
        color={127,0,0},
        thickness=0.5));

    connect(upperRespiratoryTract.q_in[1], trachea.q_in) annotation (Line(
        points={{-318.1,1.95},{-318.1,-6},{-318,-6},{-318,-2},{-298,-2},{
            -298,0}},
        color={127,0,0},
        thickness=0.5));
    connect(flowMeasure.q_out, upperRespiratoryTractResistance.q_out)
      annotation (Line(
        points={{-318,56},{-318,42}},
        color={127,0,0},
        thickness=0.5));
    connect(upperRespiratoryTractResistance.q_in, upperRespiratoryTract.q_in[
      2]) annotation (Line(
        points={{-318,22},{-318,0.65},{-318.1,0.65}},
        color={127,0,0},
        thickness=0.5));
    connect(water.port_a, gasSolubility1.liquid_port) annotation (Line(
          points={{-324,-68},{-352,-68},{-352,-48}}, color={158,66,200}));
    connect(gasSolubility1.gas_port, upperRespiratoryTract.substances[3])
      annotation (Line(points={{-352,-28},{-352,2},{-328,2},{-328,0}},
          color={158,66,200}));
    connect(pCO2.port_a, CO2_right.port_b) annotation (Line(points={{-198,-72},
            {-190,-72},{-190,-86},{-200,-86}}, color={158,66,200}));
    connect(pCO2.referenceFluidPort, rightAlveoli.q_in[3]) annotation (Line(
        points={{-208,-62.2},{-184,-62.2},{-184,-44},{-148,-44},{-148,
            -49.3},{-146.1,-49.3},{-146.1,-47.5667}},
        color={127,0,0},
        thickness=0.5));
    connect(pO2.port_a, O2_right.port_a) annotation (Line(points={{-158,-74},
            {-164,-74},{-164,-86},{-158,-86}}, color={158,66,200}));
    connect(pO2.referenceFluidPort, rightAlveoli.q_in[4]) annotation (Line(
        points={{-148,-64.2},{-148,-44},{-146.1,-44},{-146.1,-48.4333}},
        color={127,0,0},
        thickness=0.5));
    connect(pH2O_alveolar.port_a, rightAlveoli.substances[3]) annotation (
        Line(points={{-102,-76},{-94,-76},{-94,-48},{-156,-48}}, color={158,66,
            200}));
    connect(pH2O_alveolar.referenceFluidPort, rightAlveoli.q_in[5])
      annotation (Line(
        points={{-112,-66.2},{-112,-50},{-146,-50},{-146,-44},{-146.1,-44},
            {-146.1,-49.3}},
        color={127,0,0},
        thickness=0.5));
    connect(pH2O_upperRespiratory.port_a, upperRespiratoryTract.substances[3])
      annotation (Line(points={{-344,24},{-334,24},{-334,0},{-328,0}}, color=
           {158,66,200}));
    connect(pH2O_upperRespiratory.referenceFluidPort, upperRespiratoryTract.q_in[
      3]) annotation (Line(
        points={{-354,33.8},{-354,48},{-330,48},{-330,14},{-318.1,14},{
            -318.1,-0.65}},
        color={127,0,0},
        thickness=0.5));
    connect(conductor.q_out, upperRespiratoryTract.heatPort) annotation (
        Line(
        points={{-322,-34},{-324,-34},{-324,-10.2}},
        color={191,0,0},
        thickness=0.5));
    connect(conductor1.q_out, coreHeat.port) annotation (Line(
        points={{-230,16},{-238,16},{-238,-16},{-298,-16},{-298,-34},{-294,
            -34}},
        color={191,0,0},
        thickness=0.5));
    connect(coreHeat.port, conductor.q_in) annotation (Line(
        points={{-294,-34},{-302,-34}},
        color={191,0,0},
        thickness=0.5));
    connect(conductor2.q_out, coreHeat.port) annotation (Line(
        points={{-228,-30},{-238,-30},{-238,-16},{-298,-16},{-298,-34},{
            -294,-34}},
        color={191,0,0},
        thickness=0.5));
    connect(conductor1.q_in, leftAlveoli.heatPort) annotation (Line(
        points={{-210,16},{-158,15.8}},
        color={191,0,0},
        thickness=0.5));
    connect(conductor2.q_in, rightAlveoli.heatPort) annotation (Line(
        points={{-208,-30},{-168,-30},{-168,-58.2},{-152,-58.2}},
        color={191,0,0},
        thickness=0.5));
    connect(flowMeasure.q_in, Temperature_mouth.port) annotation (Line(
        points={{-318,76},{-318,82},{-298,82},{-298,72},{-286,72}},
        color={127,0,0},
        thickness=0.5));
    connect(upperRespiratoryTract.q_in[4], Temperature_upperRespiratory.port)
      annotation (Line(
        points={{-318.1,-1.95},{-318.1,10},{-318,10},{-318,8},{-288,8},{
            -288,30}},
        color={127,0,0},
        thickness=0.5));
    connect(rightAlveoli.q_in[6], Temperature_alveolar.port) annotation (
        Line(
        points={{-146.1,-50.1667},{-100,-50.1667},{-100,-40}},
        color={127,0,0},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-360,-100},{100,100}})),
      experiment(StopTime=200, __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>References:</p>
<p><br>Mecklenburgh, J. S., and W. W. Mapleson. &quot;Ventilatory assistance and respiratory muscle activity. 1: Interaction in healthy volunteers.&quot; <i>British journal of anaesthesia</i> 80.4 (1998): 422-433.</p>
</html>"));
  end PulseRespiration;

  model check_Medium
    import Modelica.Units.SI.*;

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

  model PulseRespiration2 "Minimal respiration model"
    extends Modelica.Icons.Example;

    import Modelica.Units.SI.*;

    replaceable package Air = Physiolibrary.Media.Air; //Kitware.Air_IdealGas; //Chemical.Media.SimpleAir_C; //Chemical.Media.Air_MixtureGasNasa;
    replaceable package PleuralFluid =
        Physiolibrary.Media.Water;

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
                                                                                                                //0.0133,

    inner Modelica.Fluid.System system(T_ambient=310.15)
                                       "Human body system setting"
      annotation (Placement(transformation(extent={{60,66},{80,86}})));

    Physiolibrary.Fluid.Components.ElasticVessel upperRespiratoryTract(
      redeclare package Medium = Air,
      useSubstances=true,
      volume_start=0.0001,
      useThermalPort=true,
      Compliance=TotalCompliance/100,
      ZeroPressureVolume(displayUnit="ml") = 0.0001,
      ExternalPressure=system.p_ambient,
      ResidualVolume(displayUnit="ml") = 0.0001,
      nPorts=3) annotation (Placement(transformation(extent={{-328,-10},{
              -308,10}})));
    Chemical.Sources.PureSubstance water(redeclare package stateOfMatter =
          Chemical.Interfaces.Incompressible, substanceData=
          Chemical.Substances.Water_liquid()) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-314,-68})));
    Chemical.Components.GasSolubility gasSolubility1(KC=1e-7)
                                                     annotation (
        Placement(transformation(extent={{-362,-48},{-342,-28}})));
    Physiolibrary.Fluid.Components.Resistor upperRespiratoryTractResistance(
      redeclare package Medium = Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=0.5*TracheaResistance) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-330,30})));
    Physiolibrary.Fluid.Sources.PressureSource environment(redeclare package
                Medium = Air, T=EnvironmentTemperature)
      "External environment"
      annotation (Placement(transformation(extent={{-372,76},{-352,96}})));
    Physiolibrary.Fluid.Sensors.FlowMeasure flowMeasure(redeclare package
        Medium = Air)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-330,64})));
    Chemical.Sources.ExternalMoleFraction h2o_gas(
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.Water_gas(),
      MoleFraction=0.12182580804342462373550456452011) annotation (
        Placement(transformation(extent={{-202,-22},{-182,-2}})));
    Physiolibrary.Fluid.Sensors.PartialPressure pH2O(
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.Water_gas(),
      redeclare package Medium = Air)
      annotation (Placement(transformation(extent={{-264,66},{-244,86}})));
    Physiolibrary.Fluid.Sensors.Temperature temperature(redeclare package
        Medium = Air)
      annotation (Placement(transformation(extent={{-272,8},{-252,28}})));
    Physiolibrary.Thermal.Components.Conductor conductor(Conductance(
          displayUnit="W/K") = 1) annotation (Placement(transformation(
            extent={{-308,-50},{-288,-30}})));
    Physiolibrary.Thermal.Sources.UnlimitedHeat unlimitedHeat(T=310.15)
      annotation (Placement(transformation(extent={{-202,-56},{-222,-36}})));
  equation

    connect(water.port_a, gasSolubility1.liquid_port) annotation (Line(
          points={{-324,-68},{-352,-68},{-352,-48}}, color={158,66,200}));
    connect(gasSolubility1.gas_port, upperRespiratoryTract.substances[3])
      annotation (Line(points={{-352,-28},{-352,2},{-328,2},{-328,0}},
          color={158,66,200}));
    connect(environment.y,flowMeasure. q_in) annotation (Line(
        points={{-352,86},{-330,86},{-330,74}},
        color={127,0,0},
        thickness=0.5));
    connect(flowMeasure.q_out, upperRespiratoryTractResistance.q_out)
      annotation (Line(
        points={{-330,54},{-330,40}},
        color={127,0,0},
        thickness=0.5));
    connect(upperRespiratoryTractResistance.q_in, upperRespiratoryTract.q_in[1])
          annotation (Line(
        points={{-330,20},{-318.1,20},{-318.1,1.73333}},
        color={127,0,0},
        thickness=0.5));
    connect(temperature.port, upperRespiratoryTract.q_in[2]) annotation (
        Line(points={{-262,8},{-262,0},{-318.1,0},{-318.1,1.11022e-16}},
          color={0,127,255}));
    connect(pH2O.referenceFluidPort, upperRespiratoryTract.q_in[3])
      annotation (Line(
        points={{-254,66.2},{-254,34},{-292,34},{-292,32},{-318.1,32},{
            -318.1,-1.73333}},
        color={127,0,0},
        thickness=0.5));
    connect(pH2O.port_a, upperRespiratoryTract.substances[3]) annotation (
        Line(points={{-244,76},{-236,76},{-236,-16},{-334,-16},{-334,0},{
            -328,0}}, color={158,66,200}));
    connect(upperRespiratoryTract.heatPort, conductor.q_in) annotation (
        Line(points={{-324,-10.2},{-324,-40},{-308,-40}}, color={191,0,0}));
    connect(unlimitedHeat.port, conductor.q_out) annotation (Line(
        points={{-222,-46},{-278,-46},{-278,-40},{-288,-40}},
        color={191,0,0},
        thickness=1));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-360,-100},{100,100}})),
      experiment(StopTime=200, __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>References:</p>
<p><br>Mecklenburgh, J. S., and W. W. Mapleson. &quot;Ventilatory assistance and respiratory muscle activity. 1: Interaction in healthy volunteers.&quot; <i>British journal of anaesthesia</i> 80.4 (1998): 422-433.</p>
</html>"));
  end PulseRespiration2;

  package Media

    model BloodGases_
        input Real C[12]=
        {0.44, 8.197, 20.9, 1e-8,   8.4, 0.042, 0.042, 0.66, 28, 0.153, 5.4, 37.67}
        "Volume, amount of substance or mass of substance per total volume of solution";

        input Modelica.Units.SI.Pressure p = 101325 "Pressure";
        input Modelica.Units.SI.Temperature T = 310.15 "Temperature";
        input Modelica.Units.SI.ElectricPotential electricPotential=0 "Electric potential";
        input Modelica.Units.SI.MoleFraction moleFractionBasedIonicStrength=0 "Mole-fraction based ionic strength";
        output Modelica.Units.SI.Pressure pO2( start = 101325*70/760) "Oxygen partial pressure";
        output Modelica.Units.SI.Pressure pCO2(start = 101325*40/760) "Carbon dioxide partial pressure";
        output Modelica.Units.SI.Pressure pCO(start = 101325*0.001/760) "Carbon monoxide partial pressure";

    protected
        Real Hct=C[1] "haematocrit",
             tO2=C[2] "oxygen content per volume of blood",
             tCO2=C[3] "carbon dioxide content per volume of blood",
             tCO=C[4] "carbon monoxide content per volume of blood",
             tHb=C[5] "haemoglobin content per volume of blood",
             FMetHb=C[6]/C[5] "fraction of metheamoglobin",
             FHbF=C[7]/C[5] "fraction of foetalheamoglobin",
             ctHb_ery=C[5]/C[1] "haemoglobin concentration in red cells",
             tAlb=C[8] "albumin concentration in blood plasma",
             tGlb=C[9] "globulin concentration in blood plasma",
             tPO4=C[10] "inorganic phosphates concentration in blood plasma",
             cDPG=C[11] "DPG concentration in blood plasma",
             SID=C[12] "strong ion difference of blood";


        constant Physiolibrary.Types.Temperature T0 = 273.15+37 "normal temperature";
        constant Physiolibrary.Types.pH pH0 = 7.4 "normal pH";
        constant Physiolibrary.Types.pH pH_ery0 = 7.19 "normal pH in erythrocyte";
        constant Physiolibrary.Types.Pressure pCO20 = (40/760)*101325 "normal CO2 partial pressure";

        Real NSIDP, NSIDE, NSID, BEox, pH,pH_ery, cdCO2;

        Physiolibrary.Types.GasSolubility aCO2N = 0.00023 "solubility of CO2 in blood plasma at 37 degC";
        Physiolibrary.Types.GasSolubility aCO2 = 0.00023 * 10^(-0.0092*(T-310.15)) "solubility of CO2 in blood plasma";
        Physiolibrary.Types.GasSolubility aCO2_ery( displayUnit="mmol/l/mmHg")=0.000195 "solubility 0.23 (mmol/l)/kPa at 25degC";
        Physiolibrary.Types.GasSolubility aO2= exp(log(0.0105) + (-0.0115*(T - T0)) + 0.5*0.00042*(T - T0)^2)/1000
                                                                          "oxygen solubility in blood";
        Physiolibrary.Types.GasSolubility aCO=(0.00099/0.0013)*aO2 "carbon monoxide solubility in blood";

        Real pK = 6.1 + (-0.0026)*(T-310.15) "Henderson-Hasselbalch";
        Real pK_ery = 6.125 - log10(1+10^(pH_ery-7.84-0.06*sO2));

        parameter Real pKa1=2.1, pKa2=6.8, pKa3=12.7 "H2PO4 dissociation";

        parameter Real betaOxyHb = 3.1 "buffer value for oxygenated Hb without CO2";
        parameter Real pIo=7.13  "isoelectric pH for oxygenated Hb without CO2";

        parameter Real pKzD=7.73 "pKa for NH3+ end of deoxygenated haemoglobin chain";
        parameter Real pKzO=7.25 "pKa for NH3+ end of oxygenated haemoglobin chain";
        parameter Real pKcD=7.54 "10^(pH-pKcR) is the dissociation constatnt for HbNH2 + CO2 <-> HbNHCOO- + H+ ";
        parameter Real pKcO=8.35 "10^(pH-pKcO) is the dissociation constatnt for O2HbNH2 + CO2 <-> O2HbNHCOO- + H+ ";
        parameter Real pKhD=7.52 "10^(pH-pKhD) is the dissociation constatnt for HbAH <-> HbA- + H+ ";
        parameter Real pKhO=6.89 "10^(pH-pKhO) is the dissociation constatnt for O2HbAH <-> O2HbA- + H+ ";

       Real cdCO2N;
       Real sCO2N,fzcON;

      Physiolibrary.Types.Concentration beta, cHCO3;

      Physiolibrary.Types.Fraction sO2CO(start=0.9);
      Physiolibrary.Types.Fraction sCO( start=1e-15), sO2, FCOHb;
      Physiolibrary.Types.Concentration ceHb "effective hemoglobin";

      Physiolibrary.Types.Concentration tCO2_P(start=24, displayUnit="mmol/l");
      Physiolibrary.Types.Concentration tCO2_ery( displayUnit="mmol/l");

    equation
        cdCO2N = aCO2N * pCO20 "free disolved CO2 concentration at pCO2=40mmHg and T=37degC";

        NSIDP = -(tAlb*66.463)*( 0.123 * pH0 - 0.631)
                - tGlb*(2.5/28)
                - tPO4*(10^(pKa2-pH0)+2+3*10^(pH0-pKa3))/(10^(pKa1+pKa2-2*pH0)+10^(pKa2-pH0)+1+10^(pH0-pKa3))
                - cdCO2N*10^(pH0-pK) "strong ion difference of blood plasma at pH=7.4, pCO2=40mmHg, T=37degC and sO2=1";

        fzcON = 1/(1+ 10^(pKzO-pH_ery0) + cdCO2N * 10^(pH_ery0-pKcO)) "fraction of heamoglobin units with HN2 form of amino-terminus at pH=7.4 (pH_ery=7.19), pCO2=40mmHg, T=37degC and sO2=1";
        sCO2N = 10^(pH_ery0-pKcO) * fzcON * cdCO2N "CO2 saturation of hemoglobin amino-termini at pH=7.4 (pH_ery=7.19), pCO2=40mmHg, T=37degC and sO2=1";
        NSIDE = - cdCO2N*10^(pH_ery0-pK)
                - ctHb_ery*(betaOxyHb * (pH_ery0-pIo) +  sCO2N*(1+2*10^(pKzO-pH_ery0))/(1+10^(pKzO-pH_ery0)) + 0.82)
                "strong ion difference of red cells at pH=7.4 (pH_ery=7.19), pCO2=40mmHg, T=37degC and sO2=1";

        NSID = Hct * NSIDE + (1-Hct) * NSIDP "strong ion difference of blood at pH=7.4 (pH_ery=7.19), pCO2=40mmHg, T=37degC and sO2=1";

        BEox = (-SID) - NSID "base excess of oxygenated blood";

        beta = 2.3*tHb + 8*tAlb + 0.075*tGlb + 0.309*tPO4 "buffer value of blood";

        pH = pH0 + (1/beta)*(((BEox + 0.3*(1-sO2CO))/(1-tHb/43)) - (cHCO3-24.5)) "Van Slyke";
        pH_ery = 7.19 + 0.77*(pH-7.4) + 0.035*(1-sO2);

        sO2CO =Physiolibrary.Media.BloodBySiggaardAndersen.sO2CO(
            pH,
            pO2,
            pCO2,
            pCO,
            T,
            tHb,
            cDPG,
            FMetHb,
            FHbF);

        sCO*(pO2 + 218*pCO)= 218*sO2CO* (pCO);
        FCOHb= sCO*(1-FMetHb);
        tCO = aCO*pCO + FCOHb*tHb;

        ceHb =  tHb * (1-FCOHb-FMetHb);
        sO2 = (sO2CO*(tHb*(1-FMetHb)) - tHb*FCOHb)/ceHb;
        tO2 = aO2*pO2 + ceHb*sO2;

      cdCO2 = aCO2*pCO2;
      cdCO2 * 10^(pH-pK) = cHCO3;

      tCO2_P = cHCO3 + cdCO2;
      tCO2_ery=aCO2_ery*pCO2*(1+10^(pH_ery-pK_ery));
      tCO2 = tCO2_ery*Hct + tCO2_P*(1-Hct);

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(
              preserveAspectRatio=false)));
    end BloodGases_;

    function totalO2
        input Real pH, pO2, pCO2, pCO, T, tHb, cDPG, FMetHb, FHbF;
        output Real tO2;
    protected
       constant Physiolibrary.Types.Temperature T0 = 273.15+37 "normal temperature";
        constant Physiolibrary.Types.pH pH0 = 7.4 "normal pH";
        constant Physiolibrary.Types.pH pH_ery0 = 7.19 "normal pH in erythrocyte";
        constant Physiolibrary.Types.Pressure pCO20 = (40/760)*101325 "normal CO2 partial pressure";

        parameter Physiolibrary.Types.Concentration cDPG0 = 5
        "normal DPG,used by a";
        parameter Real dadcDPG0 = 0.3 "used by a";
        parameter Real dadcDPGxHbF = -0.1 "or perhabs -0.125";
        parameter Real dadpH = -0.88 "used by a";
        parameter Real dadlnpCO2 = 0.048 "used by a";
        parameter Real dadxMetHb = -0.7 "used by a";
        parameter Real dadxHbF = -0.25 "used by a";

        Real aO2, cdO2, sO2;
        Physiolibrary.Types.Fraction sO2CO;
        Physiolibrary.Types.Pressure pO2CO;
        Physiolibrary.Types.Concentration cO2Hb;
        Physiolibrary.Types.Fraction sCO;
        Physiolibrary.Types.Concentration ceHb;
        Real a;
        Real k;
        Real x;
        Real y;
        Real h;
        Physiolibrary.Types.Fraction FCOHb;
    algorithm
        aO2 :=exp(log(0.0105) + (-0.0115*(T - T0)) + 0.5*0.00042*(T - T0)^2)/1000
                                                                          "oxygen solubility in blood";
        cdO2 :=aO2*pO2 "free disolved oxygen in blood";
        a:=dadpH*(pH-pH0)+dadlnpCO2*log(max(1e-15+1e-22*pCO2,pCO2/pCO20)) +dadxMetHb*FMetHb+(dadcDPG0 + dadcDPGxHbF*FHbF)*(cDPG/cDPG0 - 1);
        k:=0.5342857;
        h:=3.5 + a;

        pO2CO:= pO2 + 218*pCO;
        x:=log(pO2CO/7000) - a - 0.055*(T-T0);
        y:=1.8747+x+h*tanh(k*x);

        sO2CO:= exp(y)/(1+exp(y));

        sCO:=218*sO2CO* (pCO/pO2CO);
        FCOHb:=sCO*(1-FMetHb);

        ceHb := tHb * (1-FCOHb-FMetHb);

        cO2Hb:= sO2CO*(tHb*(1-FMetHb)) - tHb*FCOHb "sO2CO=(cO2Hb + tHb*FCOHb)/(tHb*(1-FMetHb))";

        sO2 := cO2Hb/ceHb;

        tO2 :=aO2*pO2 + ceHb*sO2;

    end totalO2;

    model SiggardF
        Real X[12] =    {0.44, 8.197, 21.265, 1e-6, 8.4, 0.042, 0.042,   0.66, 28, 0.153, 5.4, 37.67};
        Real MM[12] =      {1, 0.032,0.044,0.028,   1,      1,    1, 66.463,  1,     1,   1,   1};

        Modelica.Units.SI.Pressure p = 101325 "pressure";
        Modelica.Units.SI.Temperature T = 310.15 "temperature";
        input Modelica.Units.SI.ElectricPotential electricPotential=0;
        input Modelica.Units.SI.MoleFraction moleFractionBasedIonicStrength=0;
      //  output Modelica.Units.SI.ChemicalPotential u[12];
    //  protected
        Real Hct=X[1] "haematocrit",
             tO2=X[2] "oxygen content per volume of blood",
             tCO2=X[3] "carbon dioxide content per volume of blood",
             tCO=X[4] "carbon monoxide content per volume of blood",
             tHb=X[5] "haemoglobin content per volume of blood",
             FMetHb=X[6]/X[5] "fraction of metheamoglobin",
             FHbF=X[7]/X[5] "fraction of foetalheamoglobin",
             ctHb_ery=X[5]/X[1] "haemoglobin concentration in red cells",
             tAlb=X[8] "albumin concentration in blood plasma",
             tGlb=X[9] "globulin concentration in blood plasma",
             tPO4=X[10] "inorganic phosphates concentration in blood plasma",
             cDPG=X[11] "DPG concentration in blood plasma",
             SID=X[12] "strong ion difference of blood";

            // ceHbe=X[5]-X[6]-X[7] "effective haemoglobin content per volume of blood",

        constant Physiolibrary.Types.Temperature T0 = 273.15+37 "normal temperature";
        constant Physiolibrary.Types.pH pH0 = 7.4 "normal pH";
        constant Physiolibrary.Types.pH pH_ery0 = 7.19 "normal pH in erythrocyte";
        constant Physiolibrary.Types.Pressure pCO20 = (40/760)*101325 "normal CO2 partial pressure";

        Real NSIDP, NSIDE, NSID, BEox, pH,pH_ery;
        Modelica.Units.SI.Pressure pO2, pCO2, pCO;

        Physiolibrary.Types.GasSolubility aCO2N = 0.00023 "solubility of CO2 in blood plasma at 37 degC";
        Physiolibrary.Types.GasSolubility aCO2 = 0.00023 * 10^(-0.0092*(T-310.15)) "solubility of CO2 in blood plasma";
        Physiolibrary.Types.GasSolubility aCO2_ery( displayUnit="mmol/l/mmHg")=0.000195 "solubility 0.23 (mmol/l)/kPa at 25degC";
        Physiolibrary.Types.GasSolubility aO2= exp(log(0.0105) + (-0.0115*(T - T0)) + 0.5*0.00042*(T - T0)^2)/1000
                                                                          "oxygen solubility in blood";
        Physiolibrary.Types.GasSolubility aCO=(0.00099/0.0013)*aO2 "carbon monoxide solubility in blood";

        Real pK = 6.1 + (-0.0026)*(T-310.15) "Henderson-Hasselbalch";
        Real pK_ery = 6.125 - log10(1+10^(pH_ery-7.84-0.06*sO2));

        parameter Real pKa1=2.1, pKa2=6.8, pKa3=12.7 "H2PO4 dissociation";

        parameter Real betaOxyHb = 3.1 "buffer value for oxygenated Hb without CO2";
        parameter Real pIo=7.13  "isoelectric pH for oxygenated Hb without CO2";

        parameter Real pKzD=7.73 "pKa for NH3+ end of deoxygenated haemoglobin chain";
        parameter Real pKzO=7.25 "pKa for NH3+ end of oxygenated haemoglobin chain";
        parameter Real pKcD=7.54 "10^(pH-pKcR) is the dissociation constatnt for HbNH2 + CO2 <-> HbNHCOO- + H+ ";
        parameter Real pKcO=8.35 "10^(pH-pKcO) is the dissociation constatnt for O2HbNH2 + CO2 <-> O2HbNHCOO- + H+ ";
        parameter Real pKhD=7.52 "10^(pH-pKhD) is the dissociation constatnt for HbAH <-> HbA- + H+ ";
        parameter Real pKhO=6.89 "10^(pH-pKhO) is the dissociation constatnt for O2HbAH <-> O2HbA- + H+ ";

       Real cdCO2N;
       Real sCO2N,fzcON;

      Physiolibrary.Types.Concentration beta, cHCO3;

      Physiolibrary.Types.Fraction sO2CO(start=0.75);
      Physiolibrary.Types.Fraction sCO, sO2, FCOHb;
      Physiolibrary.Types.Concentration ceHb "effective hemoglobin";



    equation
        cdCO2N = aCO2N * pCO20 "free disolved CO2 concentration at pCO2=40mmHg and T=37degC";

        NSIDP = -(tAlb*66.463)*( 0.123 * pH0 - 0.631)
                - tGlb*(2.5/28)
                - tPO4*(10^(pKa2-pH0)+2+3*10^(pH0-pKa3))/(10^(pKa1+pKa2-2*pH0)+10^(pKa2-pH0)+1+10^(pH0-pKa3))
                - cdCO2N*10^(pH0-pK) "strong ion difference of blood plasma at pH=7.4, pCO2=40mmHg, T=37degC and sO2=1";

        fzcON = 1/(1+ 10^(pKzO-pH_ery0) + cdCO2N * 10^(pH_ery0-pKcO)) "fraction of heamoglobin units with HN2 form of amino-terminus at pH=7.4 (pH_ery=7.19), pCO2=40mmHg, T=37degC and sO2=1";
        sCO2N = 10^(pH_ery0-pKcO) * fzcON * cdCO2N "CO2 saturation of hemoglobin amino-termini at pH=7.4 (pH_ery=7.19), pCO2=40mmHg, T=37degC and sO2=1";
        NSIDE = - cdCO2N*10^(pH_ery0-pK)
                - ctHb_ery*(betaOxyHb * (pH_ery0-pIo) +  sCO2N*(1+2*10^(pKzO-pH_ery0))/(1+10^(pKzO-pH_ery0)) + 0.82)
                "strong ion difference of red cells at pH=7.4 (pH_ery=7.19), pCO2=40mmHg, T=37degC and sO2=1";

        NSID = Hct * NSIDE + (1-Hct) * NSIDP "strong ion difference of blood at pH=7.4 (pH_ery=7.19), pCO2=40mmHg, T=37degC and sO2=1";

        BEox = (-SID) - NSID "base excess of oxygenated blood";

        beta = 2.3*tHb + 8*tAlb + 0.075*tGlb + 0.309*tPO4 "buffer value of blood";

        pH = pH0 + (1/beta)*(((BEox + 0.3*(1-sO2CO))/(1-tHb/43)) - (cHCO3-24.5)) "Van Slyke";
        pH_ery = 7.19 + 0.77*(pH-7.4) + 0.035*(1-sO2);

        pCO2 = tCO2/(Hct*aCO2_ery*(1+10^(pH_ery-pK_ery)) + (1-Hct)*aCO2*(1+10^(pH-pK)));
        cHCO3 = aCO2*pCO2*10^(pH-pK);

        sO2CO =Physiolibrary.Media.BloodBySiggaardAndersen.sO2CO(
            pH,
            pO2,
            pCO2,
            pCO,
            T,
            tHb,
            cDPG,
            FMetHb,
            FHbF);

        (tCO - aCO*pCO)*(pO2 + 218*pCO)= 218*sO2CO* (pCO)*(1-FMetHb)*tHb;

        ceHb =  tHb * (1-FCOHb-FMetHb);
        sO2 = (sO2CO*(tHb*(1-FMetHb)) - tHb*FCOHb)/ceHb;
        tO2 = aO2*pO2 + ceHb*sO2;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(
              preserveAspectRatio=false)));
    end SiggardF;

    model Blood_test
      import Modelica.Units.SI.*;

      package Medium = Physiolibrary.Media.BloodBySiggaardAndersen;

      parameter Pressure p = 101325;
      parameter Temperature T = 310.15;

      parameter Medium.ThermodynamicState state_default = Medium.setState_pTX(p,T,Medium.X_default);

      parameter SpecificEnthalpy h = Physiolibrary.Media.BloodBySiggaardAndersen.specificEnthalpy(
                                                            state_default);
      parameter Density d = Physiolibrary.Media.BloodBySiggaardAndersen.density(
                                          state_default);

      parameter Concentration C[Physiolibrary.Media.BloodBySiggaardAndersen.nCS]=Physiolibrary.Media.BloodBySiggaardAndersen.C_default;

      parameter MolarMass MM[Medium.nCS]= Medium.molarMasses();
      parameter MassFraction x_mass_start[Medium.nCS]= (C .* MM)./d;
     // parameter Modelica.Units.SI.Mass m_start[Medium.nCS]=mass_start*x_mass_start;

      parameter MassFraction X[Physiolibrary.Media.BloodBySiggaardAndersen.nCS]=(C .* Medium.molarMasses())
           ./ Medium.density(state_default);
     // parameter MassFraction C[Blood.nCS] = {0.44, 8.197, 20.9, 1e-6, 8.4, 0.042, 0.042,   0.66, 28, 0.153, 5.4, 37.67}; // x_mass_start;

      Physiolibrary.Media.BloodBySiggaardAndersen.ChemicalSolution chemicalSolution(
        p=p,
        h=h,
        X=X) annotation (Placement(transformation(extent={{-64,6},{-44,26}})));
      Physiolibrary.Fluid.Components.ElasticVessel elasticVessel(
        redeclare package Medium =Medium,
        useSubstances=true,
        use_mass_start=false,
        mass_start=1,
        EnthalpyNotUsed=true,
        Compliance=1) annotation (Placement(transformation(extent={{8,12},{28,32}})));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(
              preserveAspectRatio=false)));

    end Blood_test;

    package Incompressible "Incompressible medium"

      extends Physiolibrary.Media.Interfaces.PartialMedium(
        ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pTX,
        singleState=true,
        reducedX=false,
        substanceNames={"H2O"},
        reference_T=298.15,
        reference_p=101325,
        reference_X={1},
        SpecificEnthalpy(start=0, nominal=1.0e5),
        Density(start=1e3, nominal=1e3),
        AbsolutePressure(start=1.0e5, nominal=1.0e5),
        Temperature(
          min=273,
          max=350,
          start=298.15));

      constant Chemical.Interfaces.Incompressible.SubstanceData substanceData[nS] = {Chemical.Substances.Water_liquid()};
      constant MolarMass MM[nS]={0.018015};
      constant SpecificHeatCapacity Cp=4180;

      replaceable model extends ChemicalSolution
        Modelica.Blocks.Interfaces.RealOutput T "temperature";
      protected
        ThermodynamicState state = setState_phX(p,h,X);

      equation
        T=state.T;

        actualStreamMolarEnthalpies =
          if EnthalpyNotUsed then zeros(nS)
          else actualStream(substances.h_outflow) "molar enthalpy in stream";

        substances.u = Chemical.Interfaces.Incompressible.chemicalPotentialPure(Water, T, p, v, I)
             "electro-chemical potential of pure water";

        substances.h_outflow =
          if EnthalpyNotUsed then zeros(nS)
          else Chemical.Interfaces.Incompressible.molarEnthalpy(Water,T,p,v,I) "molar enthalphy of the substances";

        molarFlows = substances.q;
      end ChemicalSolution;

      redeclare model extends BaseProperties(final standardOrderComponents=true)
        "Base properties of medium"

      equation
        d = 1000;
        h = stateOfMatter.molarEnthalpy(substanceData[1],T=T,p=p) / stateOfMatter.molarMass(substanceData[1]);
        u = h - p/d;
        MM = stateOfMatter.molarMass(substanceData[1]);
        R_s = 8.3144/MM;
        state.p = p;
        state.T = T;
      end BaseProperties;

      redeclare replaceable record ThermodynamicState
        "A selection of variables that uniquely defines the thermodynamic state"
        extends Modelica.Icons.Record;
        AbsolutePressure p "Absolute pressure of medium";
        Temperature T "Temperature of medium";
        /*   
    Modelica.Units.SI.ElectricPotential electricPotential=0
      "Electric potential of chemical solution";
    Modelica.Units.SI.MoleFraction moleFractionBasedIonicStrength=0
      "Ionic strangth of chemical solution";*/
        annotation (Documentation(info="<html>

</html>"));
      end ThermodynamicState;

      replaceable function molarMasses
       output Modelica.Units.SI.MolarMass molarMasses[nCS];
      algorithm
        molarMasses := {Chemical.Interfaces.Incompressible.molarMass(Water)};
        annotation(Inline=true);
      end molarMasses;

      redeclare function extends setState_pTX
        "Return thermodynamic state as function of p, T and composition X or Xi"
      algorithm
        state.p := p;
        state.T := T;
      end setState_pTX;

      redeclare function extends setState_phX
        "Return thermodynamic state as function of p, h and composition X or Xi"
      algorithm
        state.p := p;
        state.T := Chemical.Interfaces.Incompressible.solution_temperature(Water,h*Chemical.Interfaces.Incompressible.molarMass(Water),{1},p);
      end setState_phX;

      redeclare function extends specificEnthalpy "Return specific enthalpy"
      algorithm
        h := stateOfMatter.molarEnthalpy(substanceData[1],T=state.T,p=state.p) / stateOfMatter.molarMass(substanceData[1]);
      end specificEnthalpy;

      redeclare function extends specificEntropy "Return specific entropy"
      protected
        Real a "activity of substance";
        Modelica.Units.SI.MolarEnergy u
          "electro-chemical potential of substances in the solution";
      algorithm
        a := stateOfMatter.activityCoefficient(substanceData[1]);

        u := stateOfMatter.chemicalPotentialPure(
            substanceData[1],
            state.T,
            state.p) + Modelica.Constants.R*state.T*log(a);

        s := stateOfMatter.molarEntropy(u,substanceData[1],T=state.T,p=state.p) / stateOfMatter.molarMass(substanceData[1]);
        annotation (Documentation(info="<html>

</html>"));
      end specificEntropy;

      redeclare function extends specificHeatCapacityCp
        "Return specific heat capacity at constant pressure"
      algorithm
        cp := stateOfMatter.molarHeatCapacityCp(substanceData[1],T=state.T,p=state.p) / stateOfMatter.molarMass(substanceData[1]);
        annotation (Documentation(info="<html>

</html>"));
      end specificHeatCapacityCp;

      redeclare function extends specificHeatCapacityCv
        "Return specific heat capacity at constant volume"
      algorithm
        cv := stateOfMatter.molarHeatCapacityCv(substanceData[1],T=state.T,p=state.p) / stateOfMatter.molarMass(substanceData[1]);
        annotation (Documentation(info="<html>

</html>"));
      end specificHeatCapacityCv;

      redeclare function extends density
      algorithm
        d := 1000; //stateOfMatter.molarMass(substanceData[1]) / stateOfMatter.molarVolume(substanceData[1]);
      end density;

      redeclare replaceable function density_pTX "Return density from p, T, and X or Xi"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input Temperature T "Temperature";
          input MassFraction X[:] "Mass fractions";
          output Density d "Density";
      algorithm
        d := 1000;//stateOfMatter.molarMass(substanceData[1]) / stateOfMatter.molarVolume(substanceData[1]);
      end density_pTX;

      redeclare replaceable function density_phX "Return density from p, h, and X or Xi"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEnthalpy h "Specific enthalpy";
          input MassFraction X[:]=reference_X "Mass fractions";
          output Density d "Density";
      algorithm
        d := 1000; //stateOfMatter.molarMass(substanceData[1]) / stateOfMatter.molarVolume(substanceData[1]);
      end density_phX;

      redeclare function extends temperature
      algorithm
        T := state.T;
      end temperature;

      redeclare function extends pressure
      algorithm
        p := state.p;
      end pressure;
      annotation (Documentation(info="<html>
<p>
This package is a <strong>template</strong> for <strong>new medium</strong> models. For a new
medium model just make a copy of this package, remove the
\"partial\" keyword from the package and provide
the information that is requested in the comments of the
Modelica source.
</p>
</html>"));
    end Incompressible;

      package IdealGasMSL "Air as mixture of O2,CO2,H2O and N2"

        extends Modelica.Media.IdealGases.Common.MixtureGasNasa(
          mediumName="Air",
          data={
            Modelica.Media.IdealGases.Common.SingleGasesData.O2,
            Modelica.Media.IdealGases.Common.SingleGasesData.CO2,
            Modelica.Media.IdealGases.Common.SingleGasesData.H2O,
            Modelica.Media.IdealGases.Common.SingleGasesData.N2},
          fluidConstants={
            Modelica.Media.IdealGases.Common.FluidData.O2,
            Modelica.Media.IdealGases.Common.FluidData.CO2,
            Modelica.Media.IdealGases.Common.FluidData.H2O,
            Modelica.Media.IdealGases.Common.FluidData.N2},
          substanceNames = {"O2", "CO2", "H2O", "N2"},
          excludeEnthalpyOfFormation = false,
          referenceChoice = Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C,
          h_offset = 0,
          reference_X=x_default .* stateOfMatter.molarMass(substanceData) ./ (x_default*stateOfMatter.molarMass(substanceData)),
          T_default=310.15);

        package stateOfMatter = Chemical.Interfaces.IdealGasMSL
        "Substances model to translate data into substance properties";
    protected
        constant Modelica.Units.SI.MoleFraction x_default[nCS]={0.21,0.0004,0.02,
            0.7696}
          "Initial mole fractions of all substances (Please check: x_default*ones(nCS) = 1)";

    public
        constant Integer nCS=nX "Number of chemical substances";

        constant stateOfMatter.SubstanceData substanceData[nCS] = {
          stateOfMatter.SubstanceData(                data=data[i]) for i in 1:nX}
           "Definition of the substances";

        constant Modelica.Units.SI.MolarMass MM[nCS] = stateOfMatter.molarMass(substanceData);

        replaceable model ChemicalSolution
          Chemical.Interfaces.SubstancePorts_a substances[nCS];
          input Modelica.Units.SI.Pressure p "pressure";
          input Modelica.Units.SI.SpecificEnthalpy h "specific enthalpy";
          input Modelica.Units.SI.MassFraction X[nCS] "mass fractions of substances";
          input Modelica.Units.SI.ElectricPotential v=0 "electric potential";
          input Modelica.Units.SI.MoleFraction I=0 "mole fraction based ionic strength";

          Modelica.Blocks.Interfaces.RealOutput molarFlows[nCS](each unit="mol/s") "molar flows of substances";
          Modelica.Blocks.Interfaces.RealOutput actualStreamMolarEnthalpies[nCS](each unit="J/mol")
            "molar enthalpies in streams";

          parameter Boolean EnthalpyNotUsed=false annotation (
            Evaluate=true,
            HideResult=true,
            choices(checkBox=true),
            Dialog(tab="Advanced", group="Performance"));

          Modelica.Blocks.Interfaces.RealOutput T "temperature";
      protected
          ThermodynamicState state=setState_phX(
                p,
                h,
                X);
        equation
          T = state.T;

          actualStreamMolarEnthalpies = if EnthalpyNotUsed then zeros(nCS) else actualStream(substances.h_outflow)
            "molar enthalpy in stream";

          substances.u = electrochemicalPotentials_pTXvI(
              p,
              T,
              X,
              v,
              I);

          substances.h_outflow = molarEnthalpies_pTvI(
              p,
              T,
              v,
              I);

          molarFlows = substances.q;
        end ChemicalSolution;

        replaceable function electrochemicalPotentials_pTXvI
          input Modelica.Units.SI.Pressure p;
          input Modelica.Units.SI.Temperature T;
          input Modelica.Units.SI.MassFraction X[nCS];
          input Modelica.Units.SI.ElectricPotential electricPotential=0;
          input Modelica.Units.SI.MoleFraction moleFractionBasedIonicStrength=0;
          output Modelica.Units.SI.ChemicalPotential u[nCS];
      protected
          parameter Modelica.Units.SI.MolarMass MM[nCS]=stateOfMatter.molarMass(
               substanceData);
          Real a[nCS];
          Modelica.Units.SI.ChargeNumberOfIon z[nCS];
        algorithm
          a := stateOfMatter.activityCoefficient(substanceData, T, p, electricPotential, moleFractionBasedIonicStrength) .*
            (((X ./ MM)) / ((X ./ MM) * ones(nCS)));
          z := stateOfMatter.chargeNumberOfIon(substanceData, T, p, electricPotential, moleFractionBasedIonicStrength);
          u:= stateOfMatter.chemicalPotentialPure(substanceData, T, p, electricPotential, moleFractionBasedIonicStrength)
             .+ Modelica.Constants.R*T*log(a)
             .+ z*Modelica.Constants.F*electricPotential;
        end electrochemicalPotentials_pTXvI;

        replaceable function molarEnthalpies_pTvI
          input Modelica.Units.SI.Pressure p;
          input Modelica.Units.SI.Temperature T;
          input Modelica.Units.SI.ElectricPotential electricPotential=0;
          input Modelica.Units.SI.MoleFraction moleFractionBasedIonicStrength=0;
          output Modelica.Units.SI.MolarEnthalpy h[nCS];
        algorithm
          h:= stateOfMatter.molarEnthalpy(
              substanceData, T, p, electricPotential, moleFractionBasedIonicStrength);
        end molarEnthalpies_pTvI;

      end IdealGasMSL;

    model FluidAdapter
      "Adapter between chemical substances of one homogenous chemical solution and Modelica.Fluid package components of MSL 3.2, where substances are stored as molarities in expraProperties"
      import Chemical;

      outer Modelica.Fluid.System system "System wide properties";

      replaceable package Medium =
          Physiolibrary.Media.Interfaces.PartialMedium          constrainedby
        Physiolibrary.Media.Interfaces.PartialMedium
        "Medium model" annotation (choicesAllMatching=true);

      // Fluid Port definitions
      parameter Integer nFluidPorts=0 "Number of fluid ports"
        annotation (Evaluate=true, Dialog(
          connectorSizing=true,
          tab="General",
          group="Ports"));

      Modelica.Fluid.Vessels.BaseClasses.VesselFluidPorts_b fluidPorts[nFluidPorts](redeclare
          each package Medium =
                   Medium) "Fluid inlets and outlets" annotation (Placement(transformation(
            extent={{-40,-10},{40,10}},
            origin={100,0},
            rotation=90)));

      Chemical.Interfaces.SubstanceMassPorts_a substances[Medium.nCS]
        "All chemical substances of the solution" annotation (Placement(
            transformation(
            extent={{-10,-40},{10,40}},
            rotation=180,
            origin={-100,0}), iconTransformation(
            extent={{-10,-40},{10,40}},
            rotation=180,
            origin={-100,0})));
      Chemical.Interfaces.SolutionPort solution "Chemical solution"
        annotation (Placement(transformation(extent={{-50,-40},{-30,-20}}),
            iconTransformation(extent={{-50,-40},{-30,-20}})));

      Modelica.Units.SI.MassFraction x_mass[Medium.nCS]
        "Mass fraction of the substance from Chemical.Substance[]";
      Modelica.Units.SI.MassFraction xx_mass[nFluidPorts,Medium.nCS]
        "Mass fraction of the substance per actual stream in fluid port";
      Modelica.Units.SI.MassFlowRate m_flow[nFluidPorts,Medium.nCS] "Mass flow rate from fluid ports";
      Modelica.Units.SI.MassFlowRate m_flow_sum[Medium.nCS] "Mass flow rate of substance";
      Modelica.Units.SI.Concentration actualC_outflow[nFluidPorts,Medium.nC] "Actual C at fluid ports";
      Modelica.Units.SI.Concentration actualC_outflow_sum[nFluidPorts] "Sum of all C at fluid port";
      Modelica.Units.SI.MassFraction actualXi_outflow[nFluidPorts,Medium.nXi] "Actual Xi at fluid ports";
      Modelica.Units.SI.MassFraction actualXi_outflow_sum[nFluidPorts] "Sum of all Xi at fluid port";

      Modelica.Units.SI.MolarMass molarMass[Medium.nCS] "Molar mass of the substance";

      Modelica.Units.SI.Temperature temperature "Temperature of the solution";

      Modelica.Units.SI.Pressure pressure "Pressure of the solution";

      Modelica.Units.SI.ElectricPotential electricPotential(start=0) "Electric potential of the solution";

      Medium.ThermodynamicState state;

      parameter Boolean EnthalpyNotUsed=false annotation (
        Evaluate=true,
        HideResult=true,
        choices(checkBox=true),
        Dialog(tab="Advanced", group="Performance"));

    equation
      //fluid connectors
      for i in 1:nFluidPorts loop

        fluidPorts[i].p = pressure;

        //tok smerom zo substancii
        fluidPorts[i].C_outflow = Medium.C_outflow(x_mass);
        //e.g. (x_mass ./ molarMass);
        fluidPorts[i].Xi_outflow = Medium.Xi_outflow(x_mass);
        //e.g. Medium.X_default[1:Medium.nXi];

        //molarne frakce v jednotlivych fluid portoch smerom zo i do substancii
        actualC_outflow[i, :] = actualStream(fluidPorts[i].C_outflow);
        actualC_outflow_sum[i] = actualStream(fluidPorts[i].C_outflow)*ones(Medium.nC);
        actualXi_outflow[i, :] = actualStream(fluidPorts[i].Xi_outflow);
        actualXi_outflow_sum[i] = actualStream(fluidPorts[i].Xi_outflow)*ones(Medium.nXi);

        if (Medium.nXi > 0 or Medium.nC > 0) then
          xx_mass[i, :] = Medium.x_mass(actualStream_Xi=actualStream(fluidPorts[i].Xi_outflow),
            actualStream_C=actualStream(fluidPorts[i].C_outflow));
        else
          xx_mass[i, :] = ones(Medium.nCS);
        end if;

        //molarne toky v jednotlivych fluid portoch smerom zo i do substancii
        for s in 1:Medium.nCS loop
          m_flow[i, s] = (xx_mass[i, s]*fluidPorts[i].m_flow);
          // / molarMass[s];
        end for;

        //energy balance
        if (EnthalpyNotUsed) then
          fluidPorts[i].h_outflow = Medium.h_default;
        else
          fluidPorts[i].h_outflow = Medium.specificEnthalpy(state);
        end if;

      end for;

      state = Medium.setState_pTX(
        pressure,
        temperature,
        x_mass);
      //, electricPotential, solution.I);

      //substance flow balances

      for s in 1:Medium.nCS loop
        //fluidPorts.m_flow * actualStream(fluidPorts.C_outflow[s])/density + substances[s].q = 0;
        m_flow[:, s]*ones(nFluidPorts) + substances[s].m_flow = 0;
        m_flow_sum[s] = m_flow[:, s]*ones(nFluidPorts);
      end for;

      //substances
      molarMass = Medium.molarMasses();

      x_mass = substances.x_mass;

      //solution aliasses
      temperature = solution.T;
      pressure = solution.p;
      electricPotential = solution.v;

      if (EnthalpyNotUsed) then
        solution.dH = 0;
      else
        solution.dH = -fluidPorts.m_flow*(actualStream(fluidPorts.h_outflow) + (xx_mass*
          Medium.specificEnthalpyOffsets(electricPotential, solution.I)));
      end if;

      //do not affect solution at port?
      solution.i = 0;
      solution.dV = 0;
      solution.Gj = 0;
      solution.nj = 0;
      solution.mj = 0;
      solution.Qj = 0;
      solution.Ij = 0;
      solution.Vj = 0;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
              Line(
              points={{-90,0},{90,0}},
              color={158,66,200},
              thickness=1)}));
    end FluidAdapter;

    package ChemicalExamples
      model FluidAdapter
       extends Modelica.Icons.Example;

       replaceable package Medium = Physiolibrary.Media.Water;
                                                          //Chemical.Examples.Media.StandardWater_C;

        inner Modelica.Fluid.System system
          annotation (Placement(transformation(extent={{-82,66},{-62,86}})));
        Kitware.Media.FluidAdapter fluidConversion1(redeclare package
            Medium = Medium, nFluidPorts=1) annotation (Placement(
              transformation(extent={{-50,-2},{-30,18}})));
        Chemical.Components.Solution leftSolution(ConstantTemperature=true,
            BasePressure=110000)
          annotation (Placement(transformation(extent={{-96,-20},{-26,40}})));
        Chemical.Components.Substance H2O_left(substanceData=
              Chemical.Substances.Water_liquid_without_selfClustering(),
            mass_start=1) annotation (Placement(transformation(extent={{-80,
                  -2},{-60,18}})));
        Chemical.Components.Solution rightSolution(ConstantTemperature=false,
            temperature_start=299.15)
          annotation (Placement(transformation(extent={{24,-20},{98,42}})));
        Chemical.Components.Substance H2O_right(substanceData=
              Chemical.Substances.Water_liquid_without_selfClustering(),
            mass_start=1)
          annotation (Placement(transformation(extent={{84,-2},{64,18}})));
        Kitware.Media.FluidAdapter fluidConversion2(redeclare package
            Medium = Medium, nFluidPorts=1)
          annotation (Placement(transformation(extent={{56,-2},{36,18}})));
        Modelica.Fluid.Pipes.StaticPipe pipe1(
          length=1,
          diameter=0.005,
          redeclare package Medium = Medium)
          annotation (Placement(transformation(extent={{-10,-2},{10,18}})));
      equation
        connect(fluidConversion1.solution, leftSolution.solution) annotation (Line(
              points={{-44,5},{-44,-8},{-40,-8},{-40,-19.4}}, color={127,127,0}));
        connect(H2O_left.solution, leftSolution.solution) annotation (Line(points={
                {-76,-2},{-76,-8},{-40,-8},{-40,-19.4}}, color={127,127,0}));
        connect(fluidConversion2.solution, rightSolution.solution) annotation (Line(
              points={{50,5},{50,-8},{80,-8},{80,-20},{83.2,-20},{83.2,-19.38}},
              color={127,127,0}));
        connect(H2O_right.solution, rightSolution.solution) annotation (Line(points=
               {{80,-2},{80,-19.38},{83.2,-19.38}}, color={127,127,0}));
        connect(fluidConversion1.fluidPorts[1], pipe1.port_a) annotation (Line(points=
               {{-30,8},{-20,8},{-20,8},{-10,8}}, color={0,127,255}));
        connect(pipe1.port_b, fluidConversion2.fluidPorts[1])
          annotation (Line(points={{10,8},{36,8}}, color={0,127,255}));
        connect(H2O_left.port_m, fluidConversion1.substances[1]) annotation (Line(
              points={{-59.8,-2},{-54,-2},{-54,8},{-50,8}}, color={0,0,0}));
        connect(fluidConversion2.substances[1], H2O_right.port_m) annotation (Line(
              points={{56,8},{60,8},{60,-2},{63.8,-2}}, color={0,0,0}));
        annotation (    experiment(StopTime=30, __Dymola_Algorithm="Dassl"));
      end FluidAdapter;

      model FluidAdapter2
       extends Modelica.Icons.Example;

       package Medium = Chemical.Media.EthanolInWater_C;

        inner Modelica.Fluid.System system
          annotation (Placement(transformation(extent={{-82,66},{-62,86}})));
        Kitware.Media.FluidAdapter fluidConversion1(redeclare package
            Medium = Medium, nFluidPorts=1) annotation (Placement(
              transformation(extent={{-50,-2},{-30,18}})));

        Chemical.Components.Solution simpleSolution1(BasePressure=110000)
          annotation (Placement(transformation(extent={{-96,-20},{-26,40}})));
        Chemical.Components.Substance H2O(substanceData=
              Chemical.Substances.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{-90,-2},{-70,18}})));
        Chemical.Components.Solution simpleSolution2
          annotation (Placement(transformation(extent={{24,-20},{98,42}})));
        Chemical.Components.Substance H2O_(substanceData=
              Chemical.Substances.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{88,-6},{68,14}})));
        Kitware.Media.FluidAdapter fluidConversion2(redeclare package
            Medium = Medium, nFluidPorts=1)
          annotation (Placement(transformation(extent={{56,-2},{36,18}})));

        Modelica.Fluid.Pipes.StaticPipe pipe1(
          length=1,
        diameter=0.005,
        redeclare package Medium =Medium)
                      annotation (Placement(transformation(extent={{-10,-2},{10,
                18}})));
        Chemical.Components.Substance C2H5OH(
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          substanceData=Chemical.Substances.Ethanol_liquid(),
          use_mass_start=false,
          amountOfSubstance_start=100)
          annotation (Placement(transformation(extent={{-70,18},{-50,38}})));

        Chemical.Components.Substance C2H5OH_(
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          substanceData=Chemical.Substances.Ethanol_liquid(),
          use_mass_start=false,
          amountOfSubstance_start=55)
          annotation (Placement(transformation(extent={{90,20},{70,40}})));

      equation
      connect(fluidConversion1.solution, simpleSolution1.solution) annotation (
          Line(
          points={{-44,5},{-44,-8},{-40,-8},{-40,-19.4}},
          color={127,127,0}));
        connect(H2O.solution, simpleSolution1.solution) annotation (Line(
            points={{-86,-2},{-86,-8},{-40,-8},{-40,-19.4}},
            color={127,127,0}));
      connect(fluidConversion2.solution, simpleSolution2.solution) annotation (
          Line(
          points={{50,5},{50,-6},{82,-6},{82,-19.38},{83.2,-19.38}},
          color={127,127,0}));
      connect(H2O_.solution, simpleSolution2.solution) annotation (Line(
          points={{84,-6},{82,-6},{82,-19.38},{83.2,-19.38}},
          color={127,127,0}));
      connect(C2H5OH.solution, simpleSolution1.solution) annotation (Line(
          points={{-66,18},{-66,-8},{-40,-8},{-40,-19.4}},
          color={127,127,0}));
      connect(C2H5OH_.solution, simpleSolution2.solution) annotation (Line(
          points={{86,20},{86,-6},{82,-6},{82,-19.38},{83.2,-19.38}},
          color={127,127,0}));
        connect(fluidConversion1.fluidPorts[1], pipe1.port_a)
          annotation (Line(points={{-30,8},{-10,8}}, color={0,127,255}));
        connect(pipe1.port_b, fluidConversion2.fluidPorts[1])
          annotation (Line(points={{10,8},{36,8}}, color={0,127,255}));
        connect(H2O.port_m, fluidConversion1.substances[1]) annotation (Line(points=
               {{-69.8,-2},{-60,-2},{-60,8},{-50,8}}, color={0,0,0}));
        connect(C2H5OH.port_m, fluidConversion1.substances[2])
          annotation (Line(points={{-49.8,18},{-50,18},{-50,8}}, color={0,0,0}));
        connect(H2O_.port_m, fluidConversion2.substances[1]) annotation (Line(
              points={{67.8,-6},{62,-6},{62,8},{56,8}}, color={0,0,0}));
        connect(C2H5OH_.port_m, fluidConversion2.substances[2]) annotation (Line(
              points={{69.8,20},{62,20},{62,8},{56,8}}, color={0,0,0}));
        annotation (    experiment(
            StopTime=18.4),
          Documentation(info="<html>
<p>Demonstration of compatibility with FluidPort from Modelica Standard Library.</p>
</html>"));
      end FluidAdapter2;

      model FluidAdapter2_0
       extends Modelica.Icons.Example;

       package Medium = Chemical.Media.EthanolInWater_C;

        inner Modelica.Fluid.System system
          annotation (Placement(transformation(extent={{-82,66},{-62,86}})));
        Kitware.Media.FluidAdapter fluidConversion1(redeclare package
            Medium = Medium, nFluidPorts=1) annotation (Placement(
              transformation(extent={{-50,-2},{-30,18}})));

        Chemical.Components.Solution simpleSolution1(BasePressure=110000)
          annotation (Placement(transformation(extent={{-96,-20},{-26,40}})));
        Chemical.Components.Substance H2O(substanceData=
              Chemical.Substances.Water_liquid(), mass_start=1)
          annotation (Placement(transformation(extent={{-90,-2},{-70,18}})));

        Chemical.Components.Substance C2H5OH(
          redeclare package stateOfMatter = Chemical.Interfaces.Incompressible,
          substanceData=Chemical.Substances.Ethanol_liquid(),
          use_mass_start=false,
          amountOfSubstance_start=100)
          annotation (Placement(transformation(extent={{-70,18},{-50,38}})));

        Modelica.Fluid.Sources.MassFlowSource_T boundary(m_flow=-1,           redeclare
            package Medium = Medium,
          nPorts=1)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={84,58})));
        Modelica.Fluid.Sensors.TraceSubstancesTwoPort waterFlow(substanceName="H2O",
            redeclare package Medium = Medium)
          annotation (Placement(transformation(extent={{-14,48},{6,68}})));
        Modelica.Fluid.Sensors.TraceSubstancesTwoPort etchanolFlow(substanceName="C2H5OH",
            redeclare package Medium = Medium)
          annotation (Placement(transformation(extent={{18,48},{38,68}})));
        Modelica.Fluid.Sensors.MassFlowRate massFlowRate(redeclare package
            Medium =
              Medium)
          annotation (Placement(transformation(extent={{48,48},{68,68}})));
      equation
      connect(fluidConversion1.solution, simpleSolution1.solution) annotation (
          Line(
          points={{-44,5},{-44,-8},{-40,-8},{-40,-19.4}},
          color={127,127,0}));
        connect(H2O.solution, simpleSolution1.solution) annotation (Line(
            points={{-86,-2},{-86,-8},{-40,-8},{-40,-19.4}},
            color={127,127,0}));
      connect(C2H5OH.solution, simpleSolution1.solution) annotation (Line(
          points={{-66,18},{-66,-8},{-40,-8},{-40,-19.4}},
          color={127,127,0}));
        connect(fluidConversion1.fluidPorts[1], waterFlow.port_a) annotation (Line(
              points={{-30,8},{-22,8},{-22,58},{-14,58}}, color={0,127,255}));
        connect(waterFlow.port_b, etchanolFlow.port_a)
          annotation (Line(points={{6,58},{18,58}}, color={0,127,255}));
        connect(etchanolFlow.port_b, massFlowRate.port_a)
          annotation (Line(points={{38,58},{48,58}}, color={0,127,255}));
        connect(massFlowRate.port_b, boundary.ports[1])
          annotation (Line(points={{68,58},{74,58}}, color={0,127,255}));
        connect(H2O.port_m, fluidConversion1.substances[1]) annotation (Line(points=
               {{-69.8,-2},{-60,-2},{-60,8},{-50,8}}, color={0,0,0}));
        connect(C2H5OH.port_m, fluidConversion1.substances[2])
          annotation (Line(points={{-49.8,18},{-50,18},{-50,8}}, color={0,0,0}));
        annotation (    experiment(StopTime=5.5, __Dymola_Algorithm="Dassl"),
          Documentation(info="<html>
<p>Demonstration of compatibility with FluidPort from Modelica Standard Library.</p>
</html>"));
      end FluidAdapter2_0;

      model FluidAdapter_Gas
       extends Modelica.Icons.Example;

       replaceable package Medium = Kitware.Media.IdealGasMSL;

        inner Modelica.Fluid.System system
          annotation (Placement(transformation(extent={{-82,66},{-62,86}})));
        Kitware.Media.FluidAdapter fluidConversion1(redeclare package
            Medium = Medium, nFluidPorts=1) annotation (Placement(
              transformation(extent={{-50,-2},{-30,18}})));
        Chemical.Components.Solution leftSolution(redeclare package
            stateOfMatter = Medium.stateOfMatter, BasePressure=200000)
          annotation (Placement(transformation(extent={{-96,-20},{-26,40}})));
        Chemical.Components.Substance leftSubstance[Medium.nCS](
          redeclare package stateOfMatter = Medium.stateOfMatter,
          substanceData=Medium.substanceData,
          each mass_start=1) annotation (Placement(transformation(extent={
                  {-80,-2},{-60,18}})));
        Chemical.Components.Solution rightSolution(
            redeclare package stateOfMatter = Medium.stateOfMatter,
          ConstantTemperature=false,
          temperature_start=299.15)
          annotation (Placement(transformation(extent={{24,-20},{98,42}})));
        Chemical.Components.Substance rightSubstance[Medium.nCS](
          redeclare package stateOfMatter = Medium.stateOfMatter,
          substanceData=Medium.substanceData,
          each mass_start=1)
          annotation (Placement(transformation(extent={{84,-2},{64,18}})));
        Kitware.Media.FluidAdapter fluidConversion2(redeclare package
            Medium = Medium, nFluidPorts=1)
          annotation (Placement(transformation(extent={{56,-2},{36,18}})));
        Modelica.Fluid.Pipes.StaticPipe pipe1(
          length=1,
          diameter=0.005,
          redeclare package Medium = Medium)
          annotation (Placement(transformation(extent={{-10,-2},{10,18}})));
      equation
        connect(fluidConversion1.solution, leftSolution.solution) annotation (Line(
              points={{-44,5},{-44,-8},{-40,-8},{-40,-19.4}}, color={127,127,0}));
        for i in 1:Medium.nCS loop
          connect(leftSubstance[i].solution, leftSolution.solution) annotation (Line(points={
                {-76,-2},{-76,-8},{-40,-8},{-40,-19.4}}, color={127,127,0}));
          connect(rightSubstance[i].solution, rightSolution.solution) annotation (Line(points=
               {{80,-2},{80,-19.38},{83.2,-19.38}}, color={127,127,0}));

        end for;
        connect(leftSubstance.port_m, fluidConversion1.substances) annotation (Line(
              points={{-59.8,-2},{-54,-2},{-54,8},{-50,8}}, color={0,0,0}));

        connect(fluidConversion2.solution, rightSolution.solution) annotation (Line(
              points={{50,5},{50,-8},{80,-8},{80,-20},{83.2,-20},{83.2,-19.38}},
              color={127,127,0}));
        connect(fluidConversion1.fluidPorts[1], pipe1.port_a) annotation (Line(points=
               {{-30,8},{-20,8},{-20,8},{-10,8}}, color={0,127,255}));
        connect(pipe1.port_b, fluidConversion2.fluidPorts[1])
          annotation (Line(points={{10,8},{36,8}}, color={0,127,255}));
        connect(fluidConversion2.substances, rightSubstance.port_m) annotation (Line(
              points={{56,8},{60,8},{60,-2},{63.8,-2}}, color={0,0,0}));
        annotation (    experiment(StopTime=31));
      end FluidAdapter_Gas;
    end ChemicalExamples;

    package Examples
      "Integrative examples of Physiolibrary cross-domain usage"
      extends Modelica.Icons.ExamplesPackage;

    end Examples;
  end Media;

  model MinimalCirculation "Minimal circulation models driven by cardiac output"
    extends Modelica.Icons.Example;

    replaceable package Blood =
        Physiolibrary.Media.BloodBySiggaardAndersen;

    Physiolibrary.Fluid.Components.MassPump heart(
    redeclare package Medium = Blood,
    useSolutionFlowInput=
       true)
   annotation (Placement(transformation(extent={{-20,60},{0,80}})));
  Physiolibrary.Fluid.Components.ElasticVessel arteries(
   redeclare package Medium = Blood,
   volume_start(displayUnit="l") = 1e-3,
   nPorts=3,
   Compliance(displayUnit="ml/mmHg") = 1.1625954425608e-08,
   ZeroPressureVolume(displayUnit="ml") = 0.00085)
   annotation (Placement(transformation(extent={{54,40},{78,62}})));

  Physiolibrary.Fluid.Components.ElasticVessel veins(
   redeclare package Medium = Blood,
   volume_start(displayUnit="l") = 3.2e-3,
   nPorts=2,
   ZeroPressureVolume(displayUnit="ml") = 0.00295,
   Compliance(displayUnit="ml/mmHg") = 6.1880080007267e-07)
   annotation (Placement(transformation(extent={{-58,40},{-38,60}})));

    Modelica.Blocks.Sources.Pulse pulse(
   width=25,
   period=60/75,
   amplitude=3.3e-1)
   annotation (Placement(transformation(extent={{-94,74},{-74,94}})));
    Physiolibrary.Fluid.Sensors.PressureMeasure pressureMeasure(
    redeclare package Medium = Blood)
      annotation (Placement(transformation(extent={{82,68},{102,88}})));
    Physiolibrary.Fluid.Sensors.FlowMeasure flowMeasure(redeclare package Medium = Blood) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={22,32})));
    Physiolibrary.Fluid.Components.Conductor resistance(redeclare package Medium = Blood,Conductance=6.2755151845753e-09)
      annotation (Placement(transformation(extent={{-18,22},{2,42}})));
  equation
    connect(
     pulse.y, heart.solutionFlow) annotation (Line(points={{-73,84},{-10,
         84},{-10,77}}, color={0,0,127}));
    connect(
     veins.q_in[1], heart.q_in) annotation (Line(
     points={{-48.1,51.3},{-46,51.3},{-46,70},{-20,70}},
     color={127,0,0},
     thickness=0.5));
    connect(
     pressureMeasure.q_in, arteries.q_in[1]) annotation (Line(
     points={{88,72},{88,52.9067},{65.88,52.9067}},
     color={127,0,0},
     thickness=0.5));
    connect(
     resistance.q_in, veins.q_in[2]) annotation (Line(
     points={{-18,32},{-32,32},{-32,48.7},{-48.1,48.7}},
     color={127,0,0},
     thickness=0.5));
    connect(
     heart.q_out, arteries.q_in[2]) annotation (Line(
     points={{0,70},{65.88,70},{65.88,51}},
     color={127,0,0},
     thickness=0.5));
    connect(
     resistance.q_out, flowMeasure.q_out) annotation (Line(
     points={{2,32},{12,32}},
     color={127,0,0},
     thickness=0.5));
    connect(
     flowMeasure.q_in, arteries.q_in[3]) annotation (Line(
     points={{32,32},{38,32},{38,49.0933},{65.88,49.0933}},
     color={127,0,0},
     thickness=0.5));
    annotation (
   Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
           {100,100}}), graphics={Text(
         extent={{-54,98},{66,88}},
         lineColor={175,175,175},
         textString="Minimal circulation driven by cardiac output")}),
   Documentation(revisions=
                         "<html>
        <p><i>2014-2018</i></p>
        <p>Marek Matejak, marek@matfyz.cz </p>
        </html>"),
      experiment(StopTime=10));
  end MinimalCirculation;

  model PulseRespiration3 "Minimal respiration model"
    extends Modelica.Icons.Example;

    import Modelica.Units.SI.*;

    replaceable package Air = Physiolibrary.Media.Air; // Chemical.Media.SimpleAir_C; //Kitware.Air_IdealGas; //Chemical.Media.SimpleAir_C; //Chemical.Media.Air_MixtureGasNasa;
    replaceable package Blood =
        Physiolibrary.Media.BloodBySiggaardAndersen;
    replaceable package PleuralFluid =
        Physiolibrary.Media.Water;

    parameter Boolean EnthalpyNotUsed = true; // false;

    parameter Pressure IntrathoraxPressure = system.p_ambient - 700;
    parameter Frequency RespirationRate=0.08                                                               "Respiration rate";
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

    Physiolibrary.Fluid.Sensors.PressureMeasure leftPleauralPressure(redeclare package
                Medium = PleuralFluid, GetAbsolutePressure=true)
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
      annotation (Placement(transformation(extent={{-360,78},{-340,98}})));

    Physiolibrary.Fluid.Components.Resistor leftBronchi(redeclare package
        Medium =
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
    Physiolibrary.Fluid.Sensors.PressureMeasure leftAlveolarPressure(redeclare package
                Medium = Air) "Left Alveolar pressure"
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
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=RightBronchiResistance)
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
      useThermalPort=false,
      ZeroPressureVolume=pleuralVolume_initial,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=false,
      useSubstances=false,
      EnthalpyNotUsed=EnthalpyNotUsed,
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
      EnthalpyNotUsed=EnthalpyNotUsed,
      nPorts=1) "Right Plearal space"
      annotation (Placement(transformation(extent={{-76,-58},{-56,-38}})));
    Physiolibrary.Fluid.Sensors.PressureMeasure rightAlveolarPressure(redeclare
        package Medium = Air) "Right Alveolar pressure"
      annotation (Placement(transformation(extent={{-122,-54},{-102,-34}})));
    Physiolibrary.Fluid.Components.Resistor trachea(
      redeclare package Medium = Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=0.5*TracheaResistance)
      annotation (Placement(transformation(extent={{-298,-10},{-278,10}})));
    Physiolibrary.Fluid.Components.Resistor leftAlveolarDuct(
      redeclare package Medium = Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=LeftAlveoliResistance)
      annotation (Placement(transformation(extent={{-218,24},{-198,44}})));
    Physiolibrary.Fluid.Sensors.FlowMeasure flowMeasure(redeclare package
        Medium =
          Air)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-318,66})));
    Physiolibrary.Fluid.Components.Resistor rightAlveolarDuct(
      redeclare package Medium = Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=RightAlveoliResistance)
      annotation (Placement(transformation(extent={{-218,-54},{-198,-34}})));
    Physiolibrary.Fluid.Sensors.MassFractions O2_massFraction(redeclare package
                Medium =
                 Air, substanceName="O2")
      annotation (Placement(transformation(extent={{-186,72},{-166,92}})));
    Physiolibrary.Fluid.Sensors.MassFractions CO2_massFraction(redeclare package
                Medium =
                 Air, substanceName="CO2")
      annotation (Placement(transformation(extent={{-150,72},{-130,92}})));
    Physiolibrary.Fluid.Components.ElasticVessel upperRespiratoryTract(
      EnthalpyNotUsed=EnthalpyNotUsed,
      redeclare package Medium = Air,
      useSubstances=true,
      volume_start=0.0001,
      Compliance=TotalCompliance/100,
      ZeroPressureVolume(displayUnit="ml") = 0.0001,
      ExternalPressure=system.p_ambient,
      ResidualVolume(displayUnit="ml") = 0.0001,
      nPorts=3)
      annotation (Placement(transformation(extent={{-328,-10},{-308,10}})));
    Physiolibrary.Fluid.Components.Resistor upperRespiratoryTractResistance(
      redeclare package Medium = Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=0.5*TracheaResistance) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-318,32})));
    Chemical.Sources.PureSubstance water(redeclare package stateOfMatter =
          Chemical.Interfaces.Incompressible, substanceData=
          Chemical.Substances.Water_liquid()) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-314,-68})));
    Chemical.Components.GasSolubility gasSolubility1(KC=1e-5) annotation (
        Placement(transformation(extent={{-362,-48},{-342,-28}})));
    Physiolibrary.Fluid.Sensors.MassFractions H2O_massFraction(redeclare package
                Medium = Air, substanceName="H2O")
      annotation (Placement(transformation(extent={{-296,46},{-276,66}})));
    Physiolibrary.Fluid.Components.ElasticVessel blood(
      redeclare package Medium = Blood,
      useSubstances=true,
      use_mass_start=true,
      mass_start=1,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Compliance=1) annotation (Placement(transformation(extent={{-158,-176},{-138,-156}})));
    Chemical.Components.GasSolubility gasSolubility
      annotation (Placement(transformation(extent={{-212,-144},{-192,-124}})));
    Chemical.Components.GasSolubility gasSolubility2
      annotation (Placement(transformation(extent={{-194,-144},{-174,-124}})));
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
        points={{-340,88},{-318,88},{-318,76}},
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
    connect(upperRespiratoryTract.q_in[1], trachea.q_in) annotation (Line(
        points={{-318.1,1.73333},{-318.1,-6},{-318,-6},{-318,-2},{-298,-2},{
            -298,0}},
        color={127,0,0},
        thickness=0.5));
    connect(flowMeasure.q_out, upperRespiratoryTractResistance.q_out)
      annotation (Line(
        points={{-318,56},{-318,42}},
        color={127,0,0},
        thickness=0.5));
    connect(upperRespiratoryTractResistance.q_in, upperRespiratoryTract.q_in[
      2]) annotation (Line(
        points={{-318,22},{-318,2.22045e-16},{-318.1,2.22045e-16}},
        color={127,0,0},
        thickness=0.5));
    connect(water.port_a, gasSolubility1.liquid_port) annotation (Line(
          points={{-324,-68},{-352,-68},{-352,-48}}, color={158,66,200}));
    connect(gasSolubility1.gas_port, upperRespiratoryTract.substances[3])
      annotation (Line(points={{-352,-28},{-352,2},{-328,2},{-328,0}},
          color={158,66,200}));
    connect(upperRespiratoryTract.q_in[3], H2O_massFraction.port)
      annotation (Line(
        points={{-318.1,-1.73333},{-318.1,16},{-318,16},{-318,14},{-286,14},
            {-286,46}},
        color={127,0,0},
        thickness=0.5));
    connect(gasSolubility.liquid_port, blood.substances[2])
      annotation (Line(points={{-202,-144},{-202,-166},{-158,-166}}, color={158,66,200}));
    connect(gasSolubility2.liquid_port, blood.substances[3])
      annotation (Line(points={{-184,-144},{-184,-172},{-158,-172},{-158,-166}}, color={158,66,200}));
    connect(gasSolubility.gas_port, O2_right.port_a) annotation (Line(points={{-202,-124},{-202,-102},{-168,
            -102},{-168,-86},{-158,-86}}, color={158,66,200}));
    connect(CO2_right.port_b, gasSolubility2.gas_port)
      annotation (Line(points={{-200,-86},{-186,-86},{-186,-124},{-184,-124}}, color={158,66,200}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-360,-260},{100,100}})),
                                                                   Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-360,-260},{100,100}})),
      experiment(StopTime=200, __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>References:</p>
<p><br>Mecklenburgh, J. S., and W. W. Mapleson. &quot;Ventilatory assistance and respiratory muscle activity. 1: Interaction in healthy volunteers.&quot; <i>British journal of anaesthesia</i> 80.4 (1998): 422-433.</p>
</html>"));
  end PulseRespiration3;

  model PulseRespiration_T "Minimal respiration model"
    extends Modelica.Icons.Example;

    import Modelica.Units.SI.*;

    replaceable package Air = Physiolibrary.Media.Air; //Chemical.Media.SimpleAir_C; //Kitware.Air_IdealGas; //Chemical.Media.SimpleAir_C; //Chemical.Media.Air_MixtureGasNasa;
    replaceable package PleuralFluid =
        Physiolibrary.Media.Water;

    parameter Boolean EnthalpyNotUsed = false;

    parameter Pressure IntrathoraxPressure = system.p_ambient - 700;
    parameter Frequency RespirationRate=0.08                                                               "Respiration rate";
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
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      useThermalPort=true,
      ZeroPressureVolume=FunctionalResidualCapacity,
      ResidualVolume=ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=false,
      EnthalpyNotUsed=EnthalpyNotUsed,
      nPorts=2) "Left alveolar space"
      annotation (Placement(transformation(extent={{-162,16},{-142,36}})));                                     //0.0133,

    Physiolibrary.Fluid.Sensors.PressureMeasure leftPleauralPressure(redeclare package
                Medium = PleuralFluid, GetAbsolutePressure=true)
      "Left Pleaural pressure" annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=0,
          origin={-70,64})));

    inner Modelica.Fluid.System system(T_ambient=CoreTemperature)
                                       "Human body system setting"
      annotation (Placement(transformation(extent={{60,66},{80,86}})));

    Physiolibrary.Fluid.Sources.PressureSource environment(redeclare package
                Medium = Air, T=EnvironmentTemperature)
      "External environment"
      annotation (Placement(transformation(extent={{-360,78},{-340,98}})));

    Physiolibrary.Fluid.Components.Resistor leftBronchi(redeclare package
        Medium =
          Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=LeftBronchiResistance)
      annotation (Placement(transformation(extent={{-252,24},{-232,44}})));

    Physiolibrary.Types.Constants.FrequencyConst frequency(k=RespirationRate)
      annotation (Placement(transformation(extent={{-18,54},{-10,62}})));
    Physiolibrary.Fluid.Sensors.PressureMeasure leftAlveolarPressure(redeclare package
                Medium = Air) "Left Alveolar pressure"
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
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=RightBronchiResistance)
      annotation (Placement(transformation(extent={{-252,-54},{-232,-34}})));
    Physiolibrary.Fluid.Components.ElasticVessel rightAlveoli(
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      useThermalPort=true,
      ZeroPressureVolume=FunctionalResidualCapacity,
      ResidualVolume=ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=false,
      EnthalpyNotUsed=EnthalpyNotUsed,
      nPorts=3) "Right alveolar space"
      annotation (Placement(transformation(extent={{-156,-58},{-136,-38}})));
    Physiolibrary.Fluid.Components.ElasticVessel leftPlearalSpace(
      redeclare package Medium = PleuralFluid,
      volume_start=pleuralVolume_initial,
      useThermalPort=false,
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
      annotation (Placement(transformation(extent={{-130,-48},{-110,-28}})));
    Physiolibrary.Fluid.Components.Resistor trachea(
      redeclare package Medium = Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=0.5*TracheaResistance)
      annotation (Placement(transformation(extent={{-298,-10},{-278,10}})));
    Physiolibrary.Fluid.Components.Resistor leftAlveolarDuct(
      redeclare package Medium = Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=LeftAlveoliResistance)
      annotation (Placement(transformation(extent={{-218,24},{-198,44}})));
    Physiolibrary.Fluid.Sensors.FlowMeasure flowMeasure(redeclare package
        Medium =
          Air)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-318,66})));
    Physiolibrary.Fluid.Components.Resistor rightAlveolarDuct(
      redeclare package Medium = Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=RightAlveoliResistance)
      annotation (Placement(transformation(extent={{-218,-54},{-198,-34}})));
    Physiolibrary.Fluid.Components.ElasticVessel upperRespiratoryTract(
      redeclare package Medium = Air,
      useSubstances=true,
      volume_start=0.0001,
      useThermalPort=true,
      Compliance=TotalCompliance/100,
      ZeroPressureVolume(displayUnit="ml") = 0.0001,
      ExternalPressure=system.p_ambient,
      ResidualVolume(displayUnit="ml") = 0.0001,
      nPorts=4)
      annotation (Placement(transformation(extent={{-328,-10},{-308,10}})));
    Physiolibrary.Fluid.Components.Resistor upperRespiratoryTractResistance(
      redeclare package Medium = Air,
      EnthalpyNotUsed=EnthalpyNotUsed,
      Resistance=0.5*TracheaResistance) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-318,32})));
    Chemical.Sources.PureSubstance water(redeclare package stateOfMatter =
          Chemical.Interfaces.Incompressible, substanceData=
          Chemical.Substances.Water_liquid()) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-314,-68})));
    Chemical.Components.GasSolubility gasSolubility1(KC=1e-7) annotation (
        Placement(transformation(extent={{-362,-48},{-342,-28}})));
    Physiolibrary.Fluid.Sensors.Temperature Temperature_alveolar(redeclare package
                Medium = Air)
      annotation (Placement(transformation(extent={{-176,-22},{-156,-2}})));
    Physiolibrary.Fluid.Sensors.PartialPressure pH2O_upperRespiratory(
      redeclare package stateOfMatter = Chemical.Interfaces.IdealGas,
      substanceData=Chemical.Substances.Water_gas(),
      redeclare package Medium = Air)
      annotation (Placement(transformation(extent={{-364,34},{-344,14}})));
    Physiolibrary.Fluid.Sensors.Temperature Temperature_upperRespiratory(
        redeclare package Medium = Air)
      annotation (Placement(transformation(extent={{-298,30},{-278,50}})));
    Physiolibrary.Fluid.Sensors.Temperature Temperature_mouth(redeclare package
                Medium = Air)
      annotation (Placement(transformation(extent={{-296,72},{-276,92}})));
    Physiolibrary.Thermal.Components.Conductor conductor(Conductance(
          displayUnit="W/K") = 10)
      annotation (Placement(transformation(extent={{-322,-44},{-302,-24}})));
    Physiolibrary.Thermal.Sources.UnlimitedHeat coreHeat(T=system.T_ambient)
      annotation (Placement(transformation(extent={{-274,-44},{-294,-24}})));
    Physiolibrary.Thermal.Components.Conductor conductor1(Conductance(
          displayUnit="W/K") = 10)
      annotation (Placement(transformation(extent={{-212,-18},{-232,2}})));
    Physiolibrary.Thermal.Components.Conductor conductor2(Conductance(
          displayUnit="W/K") = 10)
      annotation (Placement(transformation(extent={{-212,-28},{-232,-8}})));
  equation

    connect(frequency.y, respiratoryMusclePressureCycle.frequence)
      annotation (Line(points={{-9,58},{2,58}},                      color={0,0,127}));
    connect(ambient_pressure.y, musclePressure.u1)
      annotation (Line(points={{72.75,48},{48,48},{48,32}}, color={0,0,127}));
    connect(musclePressure.u2, respiratoryMusclePressureCycle.val)
      annotation (Line(points={{36,32},{36,58},{22,58}}, color={0,0,127}));
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
        points={{-146.1,-46.2667},{-148,-46.2667},{-148,-44},{-124,-44}},
        color={127,0,0},
        thickness=0.5));

    connect(leftBronchi.q_out, leftAlveolarDuct.q_in) annotation (Line(
        points={{-232,34},{-218,34}},
        color={127,0,0},
        thickness=0.5));
    connect(leftAlveolarDuct.q_out, leftAlveoli.q_in[2]) annotation (Line(
        points={{-198,34},{-152,34},{-152,30},{-152.1,30},{-152.1,24.7}},
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
        points={{-340,88},{-318,88},{-318,76}},
        color={127,0,0},
        thickness=0.5));
    connect(rightBronchi.q_out, rightAlveolarDuct.q_in) annotation (Line(
        points={{-232,-44},{-218,-44}},
        color={127,0,0},
        thickness=0.5));
    connect(rightAlveolarDuct.q_out, rightAlveoli.q_in[2]) annotation (Line(
        points={{-198,-44},{-146,-44},{-146,-48},{-146.1,-48}},
        color={127,0,0},
        thickness=0.5));

    connect(upperRespiratoryTract.q_in[1], trachea.q_in) annotation (Line(
        points={{-318.1,1.95},{-318.1,-6},{-318,-6},{-318,-2},{-298,-2},{-298,
            0}},
        color={127,0,0},
        thickness=0.5));
    connect(flowMeasure.q_out, upperRespiratoryTractResistance.q_out)
      annotation (Line(
        points={{-318,56},{-318,42}},
        color={127,0,0},
        thickness=0.5));
    connect(upperRespiratoryTractResistance.q_in, upperRespiratoryTract.q_in[
      2]) annotation (Line(
        points={{-318,22},{-318,0.65},{-318.1,0.65}},
        color={127,0,0},
        thickness=0.5));
    connect(water.port_a, gasSolubility1.liquid_port) annotation (Line(
          points={{-324,-68},{-352,-68},{-352,-48}}, color={158,66,200}));
    connect(gasSolubility1.gas_port, upperRespiratoryTract.substances[3])
      annotation (Line(points={{-352,-28},{-352,2},{-328,2},{-328,0}},
          color={158,66,200}));
    connect(pH2O_upperRespiratory.port_a, upperRespiratoryTract.substances[3])
      annotation (Line(points={{-344,24},{-334,24},{-334,0},{-328,0}}, color=
           {158,66,200}));
    connect(pH2O_upperRespiratory.referenceFluidPort, upperRespiratoryTract.q_in[
      3]) annotation (Line(
        points={{-354,33.8},{-354,48},{-330,48},{-330,14},{-318.1,14},{-318.1,
            -0.65}},
        color={127,0,0},
        thickness=0.5));
    connect(Temperature_upperRespiratory.port, upperRespiratoryTract.q_in[4])
      annotation (Line(points={{-288,30},{-304,30},{-304,6},{-318,6},{-318,1.11022e-16},
            {-318.1,1.11022e-16},{-318.1,-1.95}}, color={0,127,255}));
    connect(Temperature_mouth.port, flowMeasure.q_in) annotation (Line(
          points={{-286,72},{-302,72},{-302,82},{-318,82},{-318,76}}, color={
            0,127,255}));
    connect(upperRespiratoryTract.heatPort, conductor.q_in) annotation (Line(
          points={{-324,-10.2},{-324,-34},{-322,-34}}, color={191,0,0}));
    connect(coreHeat.port, conductor.q_out) annotation (Line(
        points={{-294,-34},{-302,-34}},
        color={191,0,0},
        thickness=1));
    connect(Temperature_alveolar.port, rightAlveoli.q_in[3]) annotation (
        Line(points={{-166,-22},{-166,-36},{-134,-36},{-134,-44},{-146.1,
            -44},{-146.1,-49.7333}}, color={0,127,255}));
    connect(conductor1.q_in, leftAlveoli.heatPort) annotation (Line(
        points={{-212,-8},{-182,-8},{-182,6},{-158,6},{-158,15.8}},
        color={191,0,0},
        thickness=1));
    connect(conductor2.q_in, rightAlveoli.heatPort) annotation (Line(
        points={{-212,-18},{-182,-18},{-182,-64},{-152,-64},{-152,-58.2}},
        color={191,0,0},
        thickness=1));

    connect(conductor2.q_out, coreHeat.port) annotation (Line(
        points={{-232,-18},{-298,-18},{-298,-34},{-294,-34}},
        color={191,0,0},
        thickness=1));
    connect(conductor1.q_out, coreHeat.port) annotation (Line(
        points={{-232,-8},{-236,-8},{-236,-18},{-298,-18},{-298,-34},{-294,
            -34}},
        color={191,0,0},
        thickness=1));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-360,-100},{100,100}})),
      experiment(StopTime=200, __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>References:</p>
<p><br>Mecklenburgh, J. S., and W. W. Mapleson. &quot;Ventilatory assistance and respiratory muscle activity. 1: Interaction in healthy volunteers.&quot; <i>British journal of anaesthesia</i> 80.4 (1998): 422-433.</p>
</html>"));
  end PulseRespiration_T;
  annotation (uses(
      Modelica(version="4.0.0"),
      Chemical(version="1.4.0-alpha7"),
      Physiolibrary(version="3.0.0-alpha7")),
    version="1",
    conversion(from(version="", script=
            "modelica://Kitware/ConvertFromKitware_.mos")));
end Kitware;
