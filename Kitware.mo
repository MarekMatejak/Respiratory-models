within ;
package Kitware
  model RespirationMuscle "Respiration model"
    extends Modelica.Icons.Example;

    import Modelica.SIunits.*;

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
    Physiolibrary.Fluid.Components.ElasticVessel lungs(
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      ZeroPressureVolume=FunctionalResidualCapacity,
      CollapsingPressureVolume(displayUnit="l") = ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity(displayUnit="l") = VitalCapacity,
      BaseTidalVolume=5e-07,
      nHydraulicPorts=2) "Lungs"
      annotation (Placement(transformation(extent={{36,-28},{56,-8}})));

    import Modelica.SIunits.*;

    replaceable package Air = Physiolibrary.Media.Air;

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
        Medium = Air, T=EnvironmentTemperature)
      "External environment"
      annotation (Placement(transformation(extent={{-76,-30},{-56,-10}})));
    Physiolibrary.Fluid.Sensors.FlowMeasure airflowMeasure(redeclare package
        Medium = Air)         "Lungs pathway airflow"
      annotation (Placement(transformation(extent={{-40,-30},{-20,-10}})));
    Physiolibrary.Fluid.Components.Resistor resistor(redeclare package Medium =
          Air, Resistance=TotalResistance)
      annotation (Placement(transformation(extent={{-6,-30},{14,-10}})));
  equation
     pressue = clock.y;
    connect(lungsPressureMeasure.q_in, lungs.q_in[1]) annotation (Line(
        points={{76,-16},{76,-16.7},{45.7,-16.7}},
        color={127,0,0},
        thickness=0.5));
    connect(add.y, lungs.externalPressure)
      annotation (Line(points={{54,17},{54,-8}}, color={0,0,127}));
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
    connect(resistor.q_out, lungs.q_in[2]) annotation (Line(
        points={{14,-20},{14,-19.3},{45.7,-19.3}},
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

  model Respiration "Respiration model"
    extends Modelica.Icons.Example;
    Physiolibrary.Fluid.Components.ElasticVessel lungs(
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      ZeroPressureVolume=FunctionalResidualCapacity,
      CollapsingPressureVolume(displayUnit="l") = ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity(displayUnit="l") = TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=TidalVolume,
      nHydraulicPorts=2) "Lungs"
      annotation (Placement(transformation(extent={{36,-28},{56,-8}})));

    import Modelica.SIunits.*;

    replaceable package Air = Physiolibrary.Media.Air;

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
        points={{76,-16},{76,-16.7},{45.7,-16.7}},
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
        points={{14,-20},{28,-20},{28,-19.3},{45.7,-19.3}},
        color={127,0,0},
        thickness=0.5));
    connect(frequency.y, respiratoryMusclePressureCycle.frequence)
      annotation (Line(points={{-37,66},{-8,66},{-8,64},{18,64}},    color={0,0,127}));
    connect(add.y, lungs.externalPressure)
      annotation (Line(points={{54,17},{54,-8}}, color={0,0,127}));
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

    replaceable package Air = Physiolibrary.Media.Air;

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
      CollapsingPressureVolume=ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=true,
      nHydraulicPorts=2) "Lungs"
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
        points={{76,-16},{76,-16.7},{45.7,-16.7}},
        color={127,0,0},
        thickness=0.5));
    connect(resistor.q_out, lungs.q_in[2]) annotation (Line(
        points={{14,-22},{28,-22},{28,-19.3},{45.7,-19.3}},
        color={127,0,0},
        thickness=0.5));
    connect(frequency.y, respiratoryMusclePressureCycle.frequence)
      annotation (Line(points={{-37,66},{-8,66},{-8,64},{18,64}},    color={0,0,127}));
    connect(respiratoryMusclePressureCycle.val, lungs.externalPressure)
      annotation (Line(points={{38,64},{54,64},{54,-8}},color={0,0,127}));
    connect(lungs.substances[1], O2.port_a) annotation (Line(points={{36.6,-18},{36.6,
            -60},{46,-60}}, color={158,66,200}));
    connect(CO2.port_b, lungs.substances[2]) annotation (Line(points={{20,-60},{30,
            -60},{30,-18},{36.6,-18}}, color={158,66,200}));
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

    replaceable package Air = Physiolibrary.Media.Air;

    parameter Pressure IntrathoraxPressure = system.p_ambient - 700;
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
      annotation (Placement(transformation(extent={{18,54},{38,74}})));

    parameter Mass m_initial = LungsAirVolume_initial*Air.density(Air.setState_pTX(system.p_ambient
           + Pmax, CoreTemperature));

    Physiolibrary.Fluid.Components.ElasticVessel pleural(
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      ZeroPressureVolume=FunctionalResidualCapacity,
      CollapsingPressureVolume=ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=true,
      nHydraulicPorts=2) "Pleaural space with lungs"
      annotation (Placement(transformation(extent={{36,-28},{56,-8}})));                                        //0.0133,

    Physiolibrary.Fluid.Sensors.PressureMeasure pleauralPressure(redeclare
        package Medium = Air) "Pleaural pressure"
      annotation (Placement(transformation(extent={{70,-20},{90,0}})));

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
      annotation (Placement(transformation(extent={{-46,62},{-38,70}})));
    Chemical.Sources.SubstanceOutflow O2(SubstanceFlow(displayUnit="mmol/min") = 0.000257)
      annotation (Placement(transformation(extent={{66,-80},{86,-60}})));
    Chemical.Sources.SubstanceInflow CO2(SubstanceFlow(displayUnit="mmol/min") = 0.00020566666666667)
      annotation (Placement(transformation(extent={{-8,-92},{12,-72}})));
    Physiolibrary.Fluid.Components.ElasticMembrane alveoli(redeclare package
        Medium = Air, Compliance=TotalCompliance)
      annotation (Placement(transformation(extent={{-2,-28},{18,-8}})));
    Physiolibrary.Fluid.Sensors.PressureMeasure alveolarPressure(redeclare
        package Medium = Air) "Alveolar pressure"
      annotation (Placement(transformation(extent={{-10,10},{10,30}})));
    Modelica.Blocks.Math.Add add annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={54,22})));
    Physiolibrary.Types.Constants.PressureConst ambient_pressure(k=
          IntrathoraxPressure)
      annotation (Placement(transformation(extent={{84,44},{74,52}})));
  equation

    connect(pleauralPressure.q_in, pleural.q_in[1]) annotation (Line(
        points={{76,-16},{76,-16.7},{45.7,-16.7}},
        color={127,0,0},
        thickness=0.5));
    connect(frequency.y, respiratoryMusclePressureCycle.frequence)
      annotation (Line(points={{-37,66},{-8,66},{-8,64},{18,64}},    color={0,0,127}));
    connect(pleural.substances[1], O2.port_a) annotation (Line(points={{36.6,-18},
            {36,-18},{36,-70},{66,-70}},                   color={158,66,200}));
    connect(CO2.port_b, pleural.substances[2]) annotation (Line(points={{12,-82},{
            36,-82},{36,-18},{36.6,-18}}, color={158,66,200}));
    connect(airflowMeasure.q_in, environment.y) annotation (Line(
        points={{-72,-18},{-78,-18}},
        color={127,0,0},
        thickness=0.5));
    connect(airflowMeasure.q_out, resistor.q_in) annotation (Line(
        points={{-50,-18},{-42,-18}},
        color={127,0,0},
        thickness=0.5));
    connect(resistor.q_out, alveoli.q_in) annotation (Line(
        points={{-22,-18},{-2,-18}},
        color={127,0,0},
        thickness=0.5));
    connect(alveoli.q_out, pleural.q_in[2]) annotation (Line(
        points={{18,-18},{28,-18},{28,-19.3},{45.7,-19.3}},
        color={127,0,0},
        thickness=0.5));
    connect(alveolarPressure.q_in, alveoli.q_in) annotation (Line(
        points={{-4,14},{-10,14},{-10,-18},{-2,-18}},
        color={127,0,0},
        thickness=0.5));
    connect(ambient_pressure.y,add. u1)
      annotation (Line(points={{72.75,48},{60,48},{60,34}}, color={0,0,127}));
    connect(pleural.externalPressure, add.y)
      annotation (Line(points={{54,-8},{54,11},{54,11}}, color={0,0,127}));
    connect(add.u2, respiratoryMusclePressureCycle.val)
      annotation (Line(points={{48,34},{48,64},{38,64}}, color={0,0,127}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}})),                                        Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
      experiment(StopTime=16, __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>References:</p>
<p><br>Mecklenburgh, J. S., and W. W. Mapleson. &quot;Ventilatory assistance and respiratory muscle activity. 1: Interaction in healthy volunteers.&quot; <i>British journal of anaesthesia</i> 80.4 (1998): 422-433.</p>
</html>"));
  end Respiration2;

  model PulseRespiration "Minimal respiration model"
    extends Modelica.Icons.Example;

    import Modelica.SIunits.*;

    replaceable package Air = Physiolibrary.Media.Air;

    parameter Pressure IntrathoraxPressure = system.p_ambient - 700;
    parameter Frequency RespirationRate=0.133                                                                "Respiration rate";
    parameter Volume ResidualVolume=0.0013                                                      "Lungs residual volume";
    parameter Volume TotalLungCapacity=0.00623                                                      "Total Lung Capacity";
    parameter Volume BaseTidalVolume=0.0005                                                      "Base Tidal Volume";
    parameter Volume LungsAirVolume_initial = FunctionalResidualCapacity;

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



    Physiolibrary.Fluid.Components.ElasticVessel leftPleuralSpaceWithLungs(
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      ZeroPressureVolume=FunctionalResidualCapacity,
      CollapsingPressureVolume=ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=true,
      nHydraulicPorts=2) "Pleaural space with lungs"
      annotation (Placement(transformation(extent={{-88,26},{-68,46}})));                                       //0.0133,

    Physiolibrary.Fluid.Sensors.PressureMeasure pleauralPressure(redeclare
        package Medium = Air) "Pleaural pressure"
      annotation (Placement(transformation(extent={{-74,54},{-54,74}})));

    inner Modelica.Fluid.System system(T_ambient=CoreTemperature)
                                       "Human body system setting"
      annotation (Placement(transformation(extent={{60,66},{80,86}})));

    Physiolibrary.Fluid.Sources.PressureSource environment(redeclare package
        Medium =         Air, T=EnvironmentTemperature)
      "External environment"
      annotation (Placement(transformation(extent={{-356,-10},{-336,10}})));

    Physiolibrary.Fluid.Sensors.FlowMeasure airflowMeasure(redeclare package
        Medium =         Air) "Lungs pathway airflow"
      annotation (Placement(transformation(extent={{-328,-10},{-306,10}})));

    Physiolibrary.Fluid.Components.Resistor leftBronchi(redeclare package Medium =
          Air, Resistance=LeftBronchiResistance)
      annotation (Placement(transformation(extent={{-252,24},{-232,44}})));

    Physiolibrary.Types.Constants.FrequencyConst frequency(k=RespirationRate)
      annotation (Placement(transformation(extent={{-18,54},{-10,62}})));
    Chemical.Sources.SubstanceOutflow O2_left(SubstanceFlow(displayUnit="mmol/min")=
           0.0001285)
      annotation (Placement(transformation(extent={{-104,-8},{-84,12}})));
    Chemical.Sources.SubstanceInflow CO2_left(SubstanceFlow(displayUnit="mmol/min")=
           0.00010283333333333)
      annotation (Placement(transformation(extent={{-156,-8},{-136,12}})));
    Physiolibrary.Fluid.Components.ElasticMembrane leftAlveoli(redeclare
        package
        Medium = Air, Compliance=TotalCompliance)
      annotation (Placement(transformation(extent={{-156,24},{-136,44}})));
    Physiolibrary.Fluid.Sensors.PressureMeasure alveolarPressure(redeclare
        package Medium = Air) "Alveolar pressure"
      annotation (Placement(transformation(extent={{-178,52},{-158,72}})));
    Modelica.Blocks.Math.Add musclePressure annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={42,10})));
    Physiolibrary.Types.Constants.PressureConst ambient_pressure(k=
          IntrathoraxPressure)
      annotation (Placement(transformation(extent={{84,44},{74,52}})));
    Physiolibrary.Fluid.Components.Resistor leftAlveoliDuct(redeclare package
        Medium = Air, Resistance=LeftAlveoliResistance)
      annotation (Placement(transformation(extent={{-212,24},{-192,44}})));
    Physiolibrary.Fluid.Components.Resistor rightBronchi(redeclare package Medium =
          Air, Resistance=RightBronchiResistance)
      annotation (Placement(transformation(extent={{-252,-54},{-232,-34}})));
    Physiolibrary.Fluid.Components.Resistor rightAlveoliDuct(redeclare package
        Medium = Air, Resistance=RightAlveoliResistance)
      annotation (Placement(transformation(extent={{-212,-54},{-192,-34}})));
    Physiolibrary.Fluid.Components.Resistor trachea(redeclare package Medium =
          Air, Resistance=TracheaResistance)
      annotation (Placement(transformation(extent={{-294,-10},{-274,10}})));
    Physiolibrary.Fluid.Components.ElasticMembrane rightAlveoli(redeclare
        package
        Medium = Air, Compliance=TotalCompliance)
      annotation (Placement(transformation(extent={{-158,-54},{-138,-34}})));
    Physiolibrary.Fluid.Components.ElasticVessel rightPlearalSpaceWithLungs(
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      ZeroPressureVolume=FunctionalResidualCapacity,
      CollapsingPressureVolume=ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=true,
      nHydraulicPorts=1) "Pleaural space with lungs"
      annotation (Placement(transformation(extent={{-88,-52},{-68,-32}})));
    Chemical.Sources.SubstanceInflow CO2_right(SubstanceFlow(displayUnit="mmol/min")=
           0.00010283333333333)
      annotation (Placement(transformation(extent={{-162,-94},{-142,-74}})));
    Chemical.Sources.SubstanceOutflow O2_right(SubstanceFlow(displayUnit="mmol/min")=
           0.0001285)
      annotation (Placement(transformation(extent={{-100,-94},{-80,-74}})));
  equation

    connect(pleauralPressure.q_in, leftPleuralSpaceWithLungs.q_in[1]) annotation (
       Line(
        points={{-68,58},{-68,37.3},{-78.3,37.3}},
        color={127,0,0},
        thickness=0.5));
    connect(frequency.y, respiratoryMusclePressureCycle.frequence)
      annotation (Line(points={{-9,58},{2,58}},                      color={0,0,127}));
    connect(leftPleuralSpaceWithLungs.substances[1], O2_left.port_a) annotation (
        Line(points={{-87.4,36},{-116,36},{-116,2},{-104,2}}, color={158,66,200}));
    connect(CO2_left.port_b, leftPleuralSpaceWithLungs.substances[2]) annotation (
       Line(points={{-136,2},{-118,2},{-118,36},{-87.4,36}}, color={158,66,200}));
    connect(airflowMeasure.q_in, environment.y) annotation (Line(
        points={{-328,0},{-336,0}},
        color={127,0,0},
        thickness=0.5));
    connect(leftAlveoli.q_out, leftPleuralSpaceWithLungs.q_in[2]) annotation (
        Line(
        points={{-136,34},{-90,34},{-90,34.7},{-78.3,34.7}},
        color={127,0,0},
        thickness=0.5));
    connect(alveolarPressure.q_in, leftAlveoli.q_in) annotation (Line(
        points={{-172,56},{-172,34},{-156,34}},
        color={127,0,0},
        thickness=0.5));
    connect(ambient_pressure.y, musclePressure.u1)
      annotation (Line(points={{72.75,48},{48,48},{48,22}}, color={0,0,127}));
    connect(leftPleuralSpaceWithLungs.externalPressure, musclePressure.y)
      annotation (Line(points={{-70,46},{-70,-1},{42,-1}}, color={0,0,127}));
    connect(musclePressure.u2, respiratoryMusclePressureCycle.val)
      annotation (Line(points={{36,22},{36,58},{22,58}}, color={0,0,127}));
    connect(leftBronchi.q_out, leftAlveoliDuct.q_in) annotation (Line(
        points={{-232,34},{-212,34}},
        color={127,0,0},
        thickness=0.5));
    connect(leftAlveoliDuct.q_out, leftAlveoli.q_in) annotation (Line(
        points={{-192,34},{-156,34}},
        color={127,0,0},
        thickness=0.5));
    connect(rightBronchi.q_out, rightAlveoliDuct.q_in) annotation (Line(
        points={{-232,-44},{-212,-44}},
        color={127,0,0},
        thickness=0.5));
    connect(trachea.q_out, leftBronchi.q_in) annotation (Line(
        points={{-274,2.22045e-16},{-268,2.22045e-16},{-268,34},{-252,34}},
        color={127,0,0},
        thickness=0.5));
    connect(trachea.q_out, rightBronchi.q_in) annotation (Line(
        points={{-274,2.22045e-16},{-268,2.22045e-16},{-268,-44},{-252,-44}},
        color={127,0,0},
        thickness=0.5));
    connect(airflowMeasure.q_out, trachea.q_in) annotation (Line(
        points={{-306,0},{-294,0}},
        color={127,0,0},
        thickness=0.5));
    connect(rightAlveoliDuct.q_out, rightAlveoli.q_in) annotation (Line(
        points={{-192,-44},{-158,-44}},
        color={127,0,0},
        thickness=0.5));
    connect(rightAlveoli.q_out, rightPlearalSpaceWithLungs.q_in[1]) annotation (
        Line(
        points={{-138,-44},{-88,-44},{-88,-42},{-78.3,-42}},
        color={127,0,0},
        thickness=0.5));
    connect(musclePressure.y, rightPlearalSpaceWithLungs.externalPressure)
      annotation (Line(points={{42,-1},{42,-12},{-70,-12},{-70,-32}}, color={0,0,127}));
    connect(CO2_right.port_b, rightPlearalSpaceWithLungs.substances[2])
      annotation (Line(points={{-142,-84},{-120,-84},{-120,-42},{-87.4,-42}},
          color={158,66,200}));
    connect(O2_right.port_a, rightPlearalSpaceWithLungs.substances[1])
      annotation (Line(points={{-100,-84},{-118,-84},{-118,-42},{-87.4,-42}},
          color={158,66,200}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-360,-100},{100,100}})),
      experiment(StopTime=16, __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>References:</p>
<p><br>Mecklenburgh, J. S., and W. W. Mapleson. &quot;Ventilatory assistance and respiratory muscle activity. 1: Interaction in healthy volunteers.&quot; <i>British journal of anaesthesia</i> 80.4 (1998): 422-433.</p>
</html>"));
  end PulseRespiration;

  model PulseRespiration2 "Minimal respiration model"
    extends Modelica.Icons.Example;

    import Modelica.SIunits.*;

    replaceable package Air = Chemical.Examples.Media.SimpleAir_C;

    parameter Pressure IntrathoraxPressure = system.p_ambient - 700;
    parameter Frequency RespirationRate=0.133                                                                "Respiration rate";
    parameter Volume ResidualVolume=0.0013                                                      "Lungs residual volume";
    parameter Volume TotalLungCapacity=0.00623                                                      "Total Lung Capacity";
    parameter Volume BaseTidalVolume=0.0005                                                      "Base Tidal Volume";
    parameter Volume LungsAirVolume_initial = FunctionalResidualCapacity;

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

    Physiolibrary.Fluid.Components.ElasticVessel leftPleuralSpaceWithLungs(
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      ZeroPressureVolume=FunctionalResidualCapacity,
      CollapsingPressureVolume=ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=true,
      nHydraulicPorts=2) "Pleaural space with lungs"
      annotation (Placement(transformation(extent={{-88,26},{-68,46}})));                                       //0.0133,

    Physiolibrary.Fluid.Sensors.PressureMeasure pleauralPressure(redeclare
        package Medium = Air) "Pleaural pressure"
      annotation (Placement(transformation(extent={{-74,54},{-54,74}})));

    inner Modelica.Fluid.System system(T_ambient=CoreTemperature)
                                       "Human body system setting"
      annotation (Placement(transformation(extent={{60,66},{80,86}})));

    Physiolibrary.Fluid.Sources.PressureSource environment(redeclare package
        Medium =         Air, T=EnvironmentTemperature)
      "External environment"
      annotation (Placement(transformation(extent={{-356,-10},{-336,10}})));

    Physiolibrary.Fluid.Sensors.FlowMeasure airflowMeasure(redeclare package
        Medium =         Air) "Lungs pathway airflow"
      annotation (Placement(transformation(extent={{-328,-10},{-306,10}})));

    Physiolibrary.Fluid.Components.Resistor leftBronchi(redeclare package Medium =
          Air, Resistance=LeftBronchiResistance)
      annotation (Placement(transformation(extent={{-252,24},{-232,44}})));

    Physiolibrary.Types.Constants.FrequencyConst frequency(k=RespirationRate)
      annotation (Placement(transformation(extent={{-18,54},{-10,62}})));
    Chemical.Sources.SubstanceOutflow O2_left(SubstanceFlow(displayUnit="mmol/min")=
           0.0001285)
      annotation (Placement(transformation(extent={{-104,-8},{-84,12}})));
    Chemical.Sources.SubstanceInflow CO2_left(SubstanceFlow(displayUnit="mmol/min")=
           0.00010283333333333)
      annotation (Placement(transformation(extent={{-156,-8},{-136,12}})));
    Physiolibrary.Fluid.Components.ElasticMembrane leftAlveoli(redeclare
        package
        Medium = Air, Compliance=TotalCompliance)
      annotation (Placement(transformation(extent={{-156,24},{-136,44}})));
    Physiolibrary.Fluid.Sensors.PressureMeasure alveolarPressure(redeclare
        package Medium = Air) "Alveolar pressure"
      annotation (Placement(transformation(extent={{-178,52},{-158,72}})));
    Modelica.Blocks.Math.Add musclePressure annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={42,10})));
    Physiolibrary.Types.Constants.PressureConst ambient_pressure(k=
          IntrathoraxPressure)
      annotation (Placement(transformation(extent={{84,44},{74,52}})));
    Physiolibrary.Fluid.Components.Resistor leftAlveoliDuct(redeclare package
        Medium = Air, Resistance=LeftAlveoliResistance)
      annotation (Placement(transformation(extent={{-212,24},{-192,44}})));
    Physiolibrary.Fluid.Components.Resistor rightBronchi(redeclare package Medium =
          Air, Resistance=RightBronchiResistance)
      annotation (Placement(transformation(extent={{-252,-54},{-232,-34}})));
    Physiolibrary.Fluid.Components.Resistor rightAlveoliDuct(redeclare package
        Medium = Air, Resistance=RightAlveoliResistance)
      annotation (Placement(transformation(extent={{-212,-54},{-192,-34}})));
    Physiolibrary.Fluid.Components.Resistor trachea(redeclare package Medium =
          Air, Resistance=TracheaResistance)
      annotation (Placement(transformation(extent={{-294,-10},{-274,10}})));
    Physiolibrary.Fluid.Components.ElasticMembrane rightAlveoli(redeclare
        package
        Medium = Air, Compliance=TotalCompliance)
      annotation (Placement(transformation(extent={{-158,-54},{-138,-34}})));
    Physiolibrary.Fluid.Components.ElasticVessel rightPlearalSpaceWithLungs(
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      ZeroPressureVolume=FunctionalResidualCapacity,
      CollapsingPressureVolume=ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=true,
      nHydraulicPorts=1) "Pleaural space with lungs"
      annotation (Placement(transformation(extent={{-88,-52},{-68,-32}})));
    Chemical.Sources.SubstanceInflow CO2_right(SubstanceFlow(displayUnit="mmol/min")=
           0.00010283333333333)
      annotation (Placement(transformation(extent={{-162,-94},{-142,-74}})));
    Chemical.Sources.SubstanceOutflow O2_right(SubstanceFlow(displayUnit="mmol/min")=
           0.0001285)
      annotation (Placement(transformation(extent={{-100,-94},{-80,-74}})));
    Physiolibrary.Fluid.Components.ElasticVessel leftAnatomicalDeadSpace(
      redeclare package Medium = Air,
      volume_start=2e-05,
      ZeroPressureVolume=1e-05,
      Compliance=TotalCompliance/100,
      useExternalPressureInput=true,
      nHydraulicPorts=1)
      annotation (Placement(transformation(extent={{-232,66},{-212,86}})));
    Physiolibrary.Fluid.Components.ElasticVessel rightAnatomicalDeadSpace(
      redeclare package Medium = Air,
      volume_start=2e-05,
      ZeroPressureVolume=1e-05,
      Compliance=TotalCompliance/100,
      useExternalPressureInput=true,
      nHydraulicPorts=1)
      annotation (Placement(transformation(extent={{-234,-80},{-214,-60}})));
  equation

    connect(pleauralPressure.q_in, leftPleuralSpaceWithLungs.q_in[1]) annotation (
       Line(
        points={{-68,58},{-68,37.3},{-78.3,37.3}},
        color={127,0,0},
        thickness=0.5));
    connect(frequency.y, respiratoryMusclePressureCycle.frequence)
      annotation (Line(points={{-9,58},{2,58}},                      color={0,0,127}));
    connect(leftPleuralSpaceWithLungs.substances[1], O2_left.port_a) annotation (
        Line(points={{-87.4,36},{-116,36},{-116,2},{-104,2}}, color={158,66,200}));
    connect(CO2_left.port_b, leftPleuralSpaceWithLungs.substances[2]) annotation (
       Line(points={{-136,2},{-118,2},{-118,36},{-87.4,36}}, color={158,66,200}));
    connect(airflowMeasure.q_in, environment.y) annotation (Line(
        points={{-328,0},{-336,0}},
        color={127,0,0},
        thickness=0.5));
    connect(leftAlveoli.q_out, leftPleuralSpaceWithLungs.q_in[2]) annotation (
        Line(
        points={{-136,34},{-90,34},{-90,34.7},{-78.3,34.7}},
        color={127,0,0},
        thickness=0.5));
    connect(alveolarPressure.q_in, leftAlveoli.q_in) annotation (Line(
        points={{-172,56},{-172,34},{-156,34}},
        color={127,0,0},
        thickness=0.5));
    connect(ambient_pressure.y, musclePressure.u1)
      annotation (Line(points={{72.75,48},{48,48},{48,22}}, color={0,0,127}));
    connect(leftPleuralSpaceWithLungs.externalPressure, musclePressure.y)
      annotation (Line(points={{-70,46},{-70,-1},{42,-1}}, color={0,0,127}));
    connect(musclePressure.u2, respiratoryMusclePressureCycle.val)
      annotation (Line(points={{36,22},{36,58},{22,58}}, color={0,0,127}));
    connect(leftBronchi.q_out, leftAlveoliDuct.q_in) annotation (Line(
        points={{-232,34},{-212,34}},
        color={127,0,0},
        thickness=0.5));
    connect(leftAlveoliDuct.q_out, leftAlveoli.q_in) annotation (Line(
        points={{-192,34},{-156,34}},
        color={127,0,0},
        thickness=0.5));
    connect(rightBronchi.q_out, rightAlveoliDuct.q_in) annotation (Line(
        points={{-232,-44},{-212,-44}},
        color={127,0,0},
        thickness=0.5));
    connect(trachea.q_out, leftBronchi.q_in) annotation (Line(
        points={{-274,2.22045e-16},{-268,2.22045e-16},{-268,34},{-252,34}},
        color={127,0,0},
        thickness=0.5));
    connect(trachea.q_out, rightBronchi.q_in) annotation (Line(
        points={{-274,2.22045e-16},{-268,2.22045e-16},{-268,-44},{-252,-44}},
        color={127,0,0},
        thickness=0.5));
    connect(airflowMeasure.q_out, trachea.q_in) annotation (Line(
        points={{-306,0},{-294,0}},
        color={127,0,0},
        thickness=0.5));
    connect(rightAlveoliDuct.q_out, rightAlveoli.q_in) annotation (Line(
        points={{-192,-44},{-158,-44}},
        color={127,0,0},
        thickness=0.5));
    connect(rightAlveoli.q_out, rightPlearalSpaceWithLungs.q_in[1]) annotation (
        Line(
        points={{-138,-44},{-88,-44},{-88,-42},{-78.3,-42}},
        color={127,0,0},
        thickness=0.5));
    connect(musclePressure.y, rightPlearalSpaceWithLungs.externalPressure)
      annotation (Line(points={{42,-1},{42,-12},{-70,-12},{-70,-32}}, color={0,0,127}));
    connect(CO2_right.port_b, rightPlearalSpaceWithLungs.substances[2])
      annotation (Line(points={{-142,-84},{-120,-84},{-120,-42},{-87.4,-42}},
          color={158,66,200}));
    connect(O2_right.port_a, rightPlearalSpaceWithLungs.substances[1])
      annotation (Line(points={{-100,-84},{-118,-84},{-118,-42},{-87.4,-42}},
          color={158,66,200}));
    connect(leftAnatomicalDeadSpace.q_in[1], leftAlveoliDuct.q_in) annotation (
        Line(
        points={{-222.3,76},{-224,76},{-224,34},{-212,34}},
        color={127,0,0},
        thickness=0.5));
    connect(rightAnatomicalDeadSpace.q_in[1], rightAlveoliDuct.q_in) annotation (
        Line(
        points={{-224.3,-70},{-224.3,-44},{-212,-44}},
        color={127,0,0},
        thickness=0.5));
    connect(ambient_pressure.y, leftAnatomicalDeadSpace.externalPressure)
      annotation (Line(points={{72.75,48},{-70,48},{-70,86},{-214,86}}, color={0,0,
            127}));
    connect(ambient_pressure.y, rightAnatomicalDeadSpace.externalPressure)
      annotation (Line(points={{72.75,48},{-72,48},{-72,-60},{-216,-60}}, color={0,
            0,127}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-360,-100},{100,100}})),
      experiment(StopTime=16, __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>References:</p>
<p><br>Mecklenburgh, J. S., and W. W. Mapleson. &quot;Ventilatory assistance and respiratory muscle activity. 1: Interaction in healthy volunteers.&quot; <i>British journal of anaesthesia</i> 80.4 (1998): 422-433.</p>
</html>"));
  end PulseRespiration2;

  model PulseRespiration3 "Minimal respiration model"
    extends Modelica.Icons.Example;

    import Modelica.SIunits.*;

    replaceable package Air = Chemical.Examples.Media.SimpleAir_C;

    parameter Pressure IntrathoraxPressure = system.p_ambient - 700;
    parameter Frequency RespirationRate=0.133                                                                "Respiration rate";
    parameter Volume ResidualVolume=0.0013                                                      "Lungs residual volume";
    parameter Volume TotalLungCapacity=0.00623                                                      "Total Lung Capacity";
    parameter Volume BaseTidalVolume=0.0005                                                      "Base Tidal Volume";
    parameter Volume LungsAirVolume_initial = FunctionalResidualCapacity;

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

    Physiolibrary.Fluid.Components.ElasticVessel leftPleuralSpaceWithLungs(
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      ZeroPressureVolume=FunctionalResidualCapacity,
      CollapsingPressureVolume=ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=true,
      nHydraulicPorts=2) "Pleaural space with lungs"
      annotation (Placement(transformation(extent={{-88,26},{-68,46}})));                                       //0.0133,

    Physiolibrary.Fluid.Sensors.PressureMeasure pleauralPressure(redeclare
        package Medium = Air) "Pleaural pressure"
      annotation (Placement(transformation(extent={{-50,26},{-30,46}})));

    inner Modelica.Fluid.System system(T_ambient=CoreTemperature)
                                       "Human body system setting"
      annotation (Placement(transformation(extent={{-334,-84},{-314,-64}})));

    Physiolibrary.Fluid.Sources.PressureSource environment(redeclare package
        Medium =         Air, T=EnvironmentTemperature)
      "External environment"
      annotation (Placement(transformation(extent={{-356,-10},{-336,10}})));

    Physiolibrary.Fluid.Sensors.FlowMeasure airflowMeasure(redeclare package
        Medium =         Air) "Lungs pathway airflow"
      annotation (Placement(transformation(extent={{-328,-10},{-306,10}})));

    Physiolibrary.Fluid.Components.Resistor leftBronchi(redeclare package Medium =
          Air, Resistance=LeftBronchiResistance)
      annotation (Placement(transformation(extent={{-252,26},{-232,46}})));

    Physiolibrary.Types.Constants.FrequencyConst frequency(k=RespirationRate)
      annotation (Placement(transformation(extent={{-18,54},{-10,62}})));
    Chemical.Sources.SubstanceOutflow O2_left(SubstanceFlow(displayUnit="mmol/min")=
           0.0001285)
      annotation (Placement(transformation(extent={{-98,-6},{-78,14}})));
    Chemical.Sources.SubstanceInflow CO2_left(SubstanceFlow(displayUnit="mmol/min")=
           0.00010283333333333)
      annotation (Placement(transformation(extent={{-144,-6},{-124,14}})));
    Physiolibrary.Fluid.Components.ElasticMembrane leftAlveoli(redeclare
        package
        Medium = Air, Compliance=TotalCompliance)
      annotation (Placement(transformation(extent={{-140,26},{-120,46}})));
    Physiolibrary.Fluid.Sensors.PressureMeasure alveolarPressure(redeclare
        package Medium = Air) "Alveolar pressure"
      annotation (Placement(transformation(extent={{-164,36},{-144,56}})));
    Modelica.Blocks.Math.Add musclePressure annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={42,10})));
    Physiolibrary.Types.Constants.PressureConst ambient_pressure(k=
          IntrathoraxPressure)
      annotation (Placement(transformation(extent={{98,82},{88,90}})));
    Physiolibrary.Fluid.Components.Resistor leftAlveoliDuct(redeclare package
        Medium = Air, Resistance=LeftAlveoliResistance)
      annotation (Placement(transformation(extent={{-192,26},{-172,46}})));
    Physiolibrary.Fluid.Components.Resistor rightBronchi(redeclare package Medium =
          Air, Resistance=RightBronchiResistance)
      annotation (Placement(transformation(extent={{-252,-54},{-232,-34}})));
    Physiolibrary.Fluid.Components.Resistor rightAlveoliDuct(redeclare package
        Medium = Air, Resistance=RightAlveoliResistance)
      annotation (Placement(transformation(extent={{-188,-54},{-168,-34}})));
    Physiolibrary.Fluid.Components.Resistor trachea(redeclare package Medium =
          Air, Resistance=TracheaResistance)
      annotation (Placement(transformation(extent={{-294,-10},{-274,10}})));
    Physiolibrary.Fluid.Components.ElasticMembrane rightAlveoli(redeclare
        package
        Medium = Air, Compliance=TotalCompliance)
      annotation (Placement(transformation(extent={{-142,-54},{-122,-34}})));
    Physiolibrary.Fluid.Components.ElasticVessel rightPlearalSpaceWithLungs(
      redeclare package Medium = Air,
      volume_start=LungsAirVolume_initial,
      ZeroPressureVolume=FunctionalResidualCapacity,
      CollapsingPressureVolume=ResidualVolume,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      useSigmoidCompliance=true,
      VitalCapacity=TotalLungCapacity - ResidualVolume,
      BaseTidalVolume=BaseTidalVolume,
      useSubstances=true,
      nHydraulicPorts=1) "Pleaural space with lungs"
      annotation (Placement(transformation(extent={{-88,-52},{-68,-32}})));
    Chemical.Sources.SubstanceInflow CO2_right(SubstanceFlow(displayUnit="mmol/min")=
           0.00010283333333333)
      annotation (Placement(transformation(extent={{-146,-28},{-126,-8}})));
    Chemical.Sources.SubstanceOutflow O2_right(SubstanceFlow(displayUnit="mmol/min")=
           0.0001285)
      annotation (Placement(transformation(extent={{-98,-28},{-78,-8}})));
    Physiolibrary.Fluid.Components.ElasticVessel leftAnatomicalDeadSpace(
      redeclare package Medium = Air,
      volume_start=2e-05,
      ZeroPressureVolume=1e-05,
      Compliance=TotalCompliance/100,
      useExternalPressureInput=true,
      nHydraulicPorts=1)
      annotation (Placement(transformation(extent={{-232,6},{-212,26}})));
    Physiolibrary.Fluid.Components.ElasticVessel rightAnatomicalDeadSpace(
      redeclare package Medium = Air,
      volume_start=2e-05,
      ZeroPressureVolume=1e-05,
      Compliance=TotalCompliance/100,
      useExternalPressureInput=true,
      nHydraulicPorts=1)
      annotation (Placement(transformation(extent={{-232,-82},{-212,-62}})));
    Physiolibrary.Fluid.Components.ElasticVessel leftAlveolarDeadSpace(
      redeclare package Medium = Air,
      volume_start=2e-05,
      ZeroPressureVolume=0,
      Compliance=TotalCompliance,
      useExternalPressureInput=true,
      nHydraulicPorts=1)
      annotation (Placement(transformation(extent={{-140,66},{-120,86}})));
    Physiolibrary.Fluid.Components.Resistor leftAlveoliDeadSpaceDuct(redeclare
        package Medium = Air, Resistance=LeftAlveoliResistance)
      annotation (Placement(transformation(extent={{-192,66},{-172,86}})));
    Physiolibrary.Fluid.Components.Resistor rightAlveoliDeadSpaceDuct(redeclare
        package Medium = Air, Resistance=RightAlveoliResistance)
      annotation (Placement(transformation(extent={{-188,-82},{-168,-62}})));
    Physiolibrary.Fluid.Components.ElasticVessel rightAlveolarDeadSpace(
      redeclare package Medium = Air,
      volume_start=2e-05,
      ZeroPressureVolume=0,
      Compliance=TotalCompliance/100,
      useExternalPressureInput=true,
      nHydraulicPorts=1)
      annotation (Placement(transformation(extent={{-142,-82},{-122,-62}})));
  equation

    connect(pleauralPressure.q_in, leftPleuralSpaceWithLungs.q_in[1]) annotation (
       Line(
        points={{-44,30},{-44,37.3},{-78.3,37.3}},
        color={127,0,0},
        thickness=0.5));
    connect(frequency.y, respiratoryMusclePressureCycle.frequence)
      annotation (Line(points={{-9,58},{2,58}},                      color={0,0,127}));
    connect(leftPleuralSpaceWithLungs.substances[1], O2_left.port_a) annotation (
        Line(points={{-87.4,36},{-110,36},{-110,4},{-98,4}}, color={158,66,200}));
    connect(CO2_left.port_b, leftPleuralSpaceWithLungs.substances[2]) annotation (
       Line(points={{-124,4},{-112,4},{-112,36},{-87.4,36}}, color={158,66,200}));
    connect(airflowMeasure.q_in, environment.y) annotation (Line(
        points={{-328,0},{-336,0}},
        color={127,0,0},
        thickness=0.5));
    connect(leftAlveoli.q_out, leftPleuralSpaceWithLungs.q_in[2]) annotation (
        Line(
        points={{-120,36},{-90,36},{-90,34.7},{-78.3,34.7}},
        color={127,0,0},
        thickness=0.5));
    connect(alveolarPressure.q_in, leftAlveoli.q_in) annotation (Line(
        points={{-158,40},{-158,36},{-140,36}},
        color={127,0,0},
        thickness=0.5));
    connect(ambient_pressure.y, musclePressure.u1)
      annotation (Line(points={{86.75,86},{48,86},{48,22}}, color={0,0,127}));
    connect(leftPleuralSpaceWithLungs.externalPressure, musclePressure.y)
      annotation (Line(points={{-70,46},{-70,-1},{42,-1}}, color={0,0,127}));
    connect(musclePressure.u2, respiratoryMusclePressureCycle.val)
      annotation (Line(points={{36,22},{36,58},{22,58}}, color={0,0,127}));
    connect(leftBronchi.q_out, leftAlveoliDuct.q_in) annotation (Line(
        points={{-232,36},{-192,36}},
        color={127,0,0},
        thickness=0.5));
    connect(leftAlveoliDuct.q_out, leftAlveoli.q_in) annotation (Line(
        points={{-172,36},{-140,36}},
        color={127,0,0},
        thickness=0.5));
    connect(rightBronchi.q_out, rightAlveoliDuct.q_in) annotation (Line(
        points={{-232,-44},{-188,-44}},
        color={127,0,0},
        thickness=0.5));
    connect(trachea.q_out, leftBronchi.q_in) annotation (Line(
        points={{-274,2.22045e-16},{-268,2.22045e-16},{-268,36},{-252,36}},
        color={127,0,0},
        thickness=0.5));
    connect(trachea.q_out, rightBronchi.q_in) annotation (Line(
        points={{-274,2.22045e-16},{-268,2.22045e-16},{-268,-44},{-252,-44}},
        color={127,0,0},
        thickness=0.5));
    connect(airflowMeasure.q_out, trachea.q_in) annotation (Line(
        points={{-306,0},{-294,0}},
        color={127,0,0},
        thickness=0.5));
    connect(rightAlveoliDuct.q_out, rightAlveoli.q_in) annotation (Line(
        points={{-168,-44},{-142,-44}},
        color={127,0,0},
        thickness=0.5));
    connect(rightAlveoli.q_out, rightPlearalSpaceWithLungs.q_in[1]) annotation (
        Line(
        points={{-122,-44},{-88,-44},{-88,-42},{-78.3,-42}},
        color={127,0,0},
        thickness=0.5));
    connect(musclePressure.y, rightPlearalSpaceWithLungs.externalPressure)
      annotation (Line(points={{42,-1},{42,-12},{-70,-12},{-70,-32}}, color={0,0,127}));
    connect(CO2_right.port_b, rightPlearalSpaceWithLungs.substances[2])
      annotation (Line(points={{-126,-18},{-110,-18},{-110,-42},{-87.4,-42}},
          color={158,66,200}));
    connect(O2_right.port_a, rightPlearalSpaceWithLungs.substances[1])
      annotation (Line(points={{-98,-18},{-108,-18},{-108,-42},{-87.4,-42}},
          color={158,66,200}));
    connect(leftAnatomicalDeadSpace.q_in[1], leftAlveoliDuct.q_in) annotation (
        Line(
        points={{-222.3,16},{-222,16},{-222,36},{-192,36}},
        color={127,0,0},
        thickness=0.5));
    connect(rightAnatomicalDeadSpace.q_in[1], rightAlveoliDuct.q_in) annotation (
        Line(
        points={{-222.3,-72},{-222.3,-44},{-188,-44}},
        color={127,0,0},
        thickness=0.5));
    connect(ambient_pressure.y, leftAnatomicalDeadSpace.externalPressure)
      annotation (Line(points={{86.75,86},{-214,86},{-214,26}}, color={0,0,127}));
    connect(ambient_pressure.y, rightAnatomicalDeadSpace.externalPressure)
      annotation (Line(points={{86.75,86},{-20,86},{-20,-62},{-214,-62}}, color={0,
            0,127}));
    connect(leftAlveoliDeadSpaceDuct.q_in, leftAlveoliDuct.q_in) annotation (Line(
        points={{-192,76},{-192,78},{-208,78},{-208,36},{-192,36}},
        color={127,0,0},
        thickness=0.5));
    connect(leftAlveoliDeadSpaceDuct.q_out, leftAlveolarDeadSpace.q_in[1])
      annotation (Line(
        points={{-172,76},{-130.3,76}},
        color={127,0,0},
        thickness=0.5));
    connect(ambient_pressure.y, leftAlveolarDeadSpace.externalPressure)
      annotation (Line(points={{86.75,86},{-122,86}}, color={0,0,127}));
    connect(rightAlveoliDeadSpaceDuct.q_out, rightAlveolarDeadSpace.q_in[1])
      annotation (Line(
        points={{-168,-72},{-132.3,-72}},
        color={127,0,0},
        thickness=0.5));
    connect(rightAlveolarDeadSpace.externalPressure, rightAnatomicalDeadSpace.externalPressure)
      annotation (Line(points={{-124,-62},{-214,-62}}, color={0,0,127}));
    connect(rightAlveoliDeadSpaceDuct.q_in, rightAlveoliDuct.q_in) annotation (
        Line(
        points={{-188,-72},{-206,-72},{-206,-44},{-188,-44}},
        color={127,0,0},
        thickness=0.5));
    connect(ambient_pressure.y, rightAlveolarDeadSpace.externalPressure)
      annotation (Line(points={{86.75,86},{-20,86},{-20,-62},{-124,-62}}, color={0,
            0,127}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-360,-100},{100,100}})),
      experiment(StopTime=16, __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
<p>References:</p>
<p><br>Mecklenburgh, J. S., and W. W. Mapleson. &quot;Ventilatory assistance and respiratory muscle activity. 1: Interaction in healthy volunteers.&quot; <i>British journal of anaesthesia</i> 80.4 (1998): 422-433.</p>
</html>"));
  end PulseRespiration3;
  annotation (uses(
      Modelica(version="3.2.3"),
      Chemical(version="1.3.1"),
      Physiolibrary(version="3.0.0")));
end Kitware;
