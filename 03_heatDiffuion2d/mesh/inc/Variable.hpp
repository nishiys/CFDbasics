#pragma once

class Variable
{
public:
    Variable();
    ~Variable();

    // Temperature
    double temperature;
    // Heat Source Per Unit Volume (W/m3)
    double volumetric_source;
    // Thermal Conducivity (W/mK)
    double thermal_cond;
    // Thickness of the plate (m)
    double thickness;
};