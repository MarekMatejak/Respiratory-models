within ;
model AlveolarVentilationRegulationParameters
  Real Gp(start=30.24),Gc(start=1.44),Ip(start=35.5),Ic(start=35.5);
  inner Modelica.Fluid.System system
    annotation (Placement(transformation(extent={{-2,16},{18,36}})));
equation

  3.88  = Gp*exp(-0.05*96.54)*(40-Ip) + Gc*(40-Ic);
  11.31  = Gp*exp(-0.05*96.01)*(42.12-Ip) + Gc*(42.12-Ic);
  31.61 = Gp*exp(-0.05*95.45)*(43.18-Ip) + Gc*(43.18-Ic);
  58.3  = Gp*exp(-0.05*95.18)*(45.21-Ip) + Gc*(45.21-Ic);

  annotation (uses(Modelica(version="3.2.3")));
end AlveolarVentilationRegulationParameters;
