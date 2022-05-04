L1 = 0.3; % meters
L2 = 0.1; % meters
g = 9.81; % meters/s^2
woodDensity = 500; % kg/m^3

A1 = 0.01^2;
A2 = 0.01^2;

motorMass1 = 0.2;
motorMass2 = 0.2;

waterBottleMass = 0.025 + 0.15; % kg (bottle + water)

Tground = (L1/2)*A1*woodDensity*g + (L2/2+L1)*A2*woodDensity*g + L1*motorMass1*g + (L1+L2)*motorMass2*g + (L1+L2+0.05)*waterBottleMass*g;

motorTorque1 = (1/100)*9.8*(30);