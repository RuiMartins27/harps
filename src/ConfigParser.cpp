#include "ConfigParser.h"

#include <fstream>
#include <sstream>
#include <regex>
#include <iostream>
#include <algorithm>

HARPSConfig ConfigParser::parseFile(const std::string& filename) {
    HARPSConfig config;
    std::map<std::string, std::string> values;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open config file: " + filename);
    }
    
    std::string line, section;
    while (std::getline(file, line)) {
        // Remove comments and trim whitespace
        size_t commentPos = line.find('#');
        if (commentPos != std::string::npos) {
            line = line.substr(0, commentPos);
        }
        line = trim(line);
        
        if (line.empty()) continue;
        
        if (line.find("Start(") == 0) {
            section = line.substr(6, line.find(')') - 6);
        } else if (line.find("End(") == 0) {
            section = "";
        } else if (line.find('=') != std::string::npos) {
            auto parts = split(line, '=');
            if (parts.size() == 2) {
                std::string key = trim(parts[0]);
                std::string value = trim(parts[1]);
                values[section + "." + key] = value;
            }
        }
    }
    
    // Parameters are only parsed if they are on the input file, this function won't try every possible one
    // Parse grid parameters
    
    if (values.find("Grid.Nx") != values.end()){
        config.n_x = string_to_int(values["Grid.Nx"],"Grid.Nx",'p');
        if (config.n_x == 2) throw std::runtime_error("Grid.Nx cannot be 2, use 1 for assuming reduced dimension or >2 for considering this dimension.");
    } else throw std::runtime_error("Grid.Nx needed");

    if (values.find("Grid.Ny") != values.end()){
        config.n_y = string_to_int(values["Grid.Ny"],"Grid.Ny",'p');
        if (config.n_y == 2) throw std::runtime_error("Grid.Ny cannot be 2, use 1 for assuming reduced dimension or >2 for considering this dimension.");
    }else throw std::runtime_error("Grid.Ny needed");

    if (values.find("Grid.Nz") != values.end()){
        config.n_z = string_to_int(values["Grid.Nz"],"Grid.Nz",'p');
        if (config.n_z == 2) throw std::runtime_error("Grid.Nz cannot be 2, use 1 for assuming reduced dimension or >2 for considering this dimension.");
    }else throw std::runtime_error("Grid.Nz needed");

    if (values.find("Grid.LengthX") != values.end()) config.lengthX = string_to_double(values["Grid.LengthX"],"Grid.LengthX",'p');
    else throw std::runtime_error("Grid.LengthX needed");

    if (values.find("Grid.LengthY") != values.end()) config.lengthY = string_to_double(values["Grid.LengthY"],"Grid.LengthY",'p');
    else throw std::runtime_error("Grid.LengthY needed");

    if (values.find("Grid.LengthZ") != values.end()) config.lengthZ = string_to_double(values["Grid.LengthZ"],"Grid.LengthX",'p');
    else throw std::runtime_error("Grid.LengthZ needed");

    // Non uniform grid refinement
    if (values.find("Grid.RefinementFactorX") != values.end()) {
        config.refinementFactor_x = string_to_double(values["Grid.RefinementFactorX"],"Grid.RefinementFactorX");

        if(values.find("Grid.x_BL") != values.end()) config.x_BL = string_to_double(values["Grid.x_BL"],"Grid.x_BL",'p');
        else throw std::runtime_error("Grid.x_BL needed for Grid refinement");
        if(values.find("Grid.x_BR") != values.end()) config.x_BR = string_to_double(values["Grid.x_BR"],"Grid.x_BR",'p');
        else throw std::runtime_error("Grid.x_BR needed for Grid refinement");
        if(values.find("Grid.x_RL") != values.end()) config.x_RL = string_to_double(values["Grid.x_RL"],"Grid.x_RL",'p');
        else throw std::runtime_error("Grid.x_RL needed for Grid refinement");
        if(values.find("Grid.x_RR") != values.end()) config.x_RR = string_to_double(values["Grid.x_RR"],"Grid.x_RR",'p');
        else throw std::runtime_error("Grid.x_RR needed for Grid refinement");
    }
    if (values.find("Grid.RefinementFactorY") != values.end()){
        config.refinementFactor_y = string_to_double(values["Grid.RefinementFactorY"],"Grid.RefinementFactorY");

        if(values.find("Grid.y_BL") != values.end()) config.y_BL = string_to_double(values["Grid.y_BL"],"Grid.y_BL",'p');
        else throw std::runtime_error("Grid.y_BL needed for Grid refinement");
        if(values.find("Grid.y_BR") != values.end()) config.y_BR = string_to_double(values["Grid.y_BR"],"Grid.y_BR",'p');
        else throw std::runtime_error("Grid.y_BR needed for Grid refinement");
        if(values.find("Grid.y_RL") != values.end()) config.y_RL = string_to_double(values["Grid.y_RL"],"Grid.y_RL",'p');
        else throw std::runtime_error("Grid.y_RL needed for Grid refinement");
        if(values.find("Grid.y_RR") != values.end()) config.y_RR = string_to_double(values["Grid.y_RR"],"Grid.y_RR",'p');
        else throw std::runtime_error("Grid.y_RR needed for Grid refinement");
    }
    if (values.find("Grid.RefinementFactorZ") != values.end()){
        config.refinementFactor_z = string_to_double(values["Grid.RefinementFactorZ"],"Grid.RefinementFactorZ");

        if(values.find("Grid.z_BL") != values.end()) config.z_BL = string_to_double(values["Grid.z_BL"],"Grid.z_BL",'p');
        else throw std::runtime_error("Grid.z_BL needed for Grid refinement");
        if(values.find("Grid.z_BR") != values.end()) config.z_BR = string_to_double(values["Grid.z_BR"],"Grid.z_BR",'p');
        else throw std::runtime_error("Grid.z_BR needed for Grid refinement");
        if(values.find("Grid.z_RL") != values.end()) config.z_RL = string_to_double(values["Grid.z_RL"],"Grid.z_RL",'p');
        else throw std::runtime_error("Grid.z_RL needed for Grid refinement");
        if(values.find("Grid.z_RR") != values.end()) config.z_RR = string_to_double(values["Grid.z_RR"],"Grid.z_RR",'p');
        else throw std::runtime_error("Grid.z_RR needed for Grid refinement");
    }

    // Directy providing the grid points
    if (values.find("Grid.x_Grid") != values.end()){
        if (config.refinementFactor_x != -1) throw std::runtime_error("Grid.x_Grid and Grid.RefinementFactorX cannot be used at the same time.");
        config.x_grid = parseDoubleValueArray(values["Grid.x_Grid"], config);
        if (static_cast<int>(config.x_grid.size()) != config.n_x) {
            throw std::runtime_error("Number of x_Grid points (" + std::to_string(config.x_grid.size()) +  ") does not match n_x (" + std::to_string(config.n_x) + ").");
        }
    }
    if (values.find("Grid.y_Grid") != values.end()){
        if (config.refinementFactor_y != -1) throw std::runtime_error("Grid.y_Grid and Grid.RefinementFactorY cannot be used at the same time.");
        config.y_grid = parseDoubleValueArray(values["Grid.y_Grid"], config);
        if (static_cast<int>(config.y_grid.size()) != config.n_y) {
            throw std::runtime_error("Number of y_Grid points (" + std::to_string(config.y_grid.size()) +  ") does not match n_y (" + std::to_string(config.n_y) + ").");
        }
    }
    if (values.find("Grid.z_Grid") != values.end()){
        if (config.refinementFactor_z != -1) throw std::runtime_error("Grid.z_Grid and Grid.RefinementFactorZ cannot be used at the same time.");
        config.z_grid = parseDoubleValueArray(values["Grid.z_Grid"], config);
        if (static_cast<int>(config.z_grid.size()) != config.n_z) {
            throw std::runtime_error("Number of z_Grid points (" + std::to_string(config.z_grid.size()) +  ") does not match n_z (" + std::to_string(config.n_z) + ").");
        }
    }


    // Parse physical parameters
    if (values.find("PhysicalParameters.Frequency") != values.end()) config.frequency = string_to_double(values["PhysicalParameters.Frequency"],"PhysicalParameters.Frequency",'p');
    else throw std::runtime_error("PhysicalParameters.Frequency needed");

    if (values.find("PhysicalParameters.WaveguideNumber") != values.end()) config.waveguide_number = string_to_double(values["PhysicalParameters.WaveguideNumber"],"PhysicalParameters.WaveguideNumber",'p');

    if (values.find("PhysicalParameters.yCenter") != values.end()) config.yCenter = string_to_double(values["PhysicalParameters.yCenter"],"PhysicalParameters.yCenter",'p');
    
    if (values.find("PhysicalParameters.FlagCylindricalPlasma") != values.end()) config.flag_cylindrical_plasma = parseBoolean(values["PhysicalParameters.FlagCylindricalPlasma"], "PhysicalParameters.FlagCyllindricalPlasma");

    if (values.find("PhysicalParameters.RealMobility") != values.end()) config.realMobility = string_to_double(values["PhysicalParameters.RealMobility"],"PhysicalParameters.RealMobility",'a');
    else throw std::runtime_error("PhysicalParameters.RealMobility needed");

    if (values.find("PhysicalParameters.ImagMobility") != values.end()) config.imagMobility = string_to_double(values["PhysicalParameters.ImagMobility"],"PhysicalParameters.ImagMobility",'a');
    else throw std::runtime_error("PhysicalParameters.ImagMobility needed");

    if (values.find("PhysicalParameters.RealPermittivity") != values.end()) config.realPermittivity = string_to_double(values["PhysicalParameters.RealPermittivity"],"PhysicalParameters.RealPermittivity",'p');
    else throw std::runtime_error("PhysicalParameters.RealPermittivity needed");

    if (values.find("PhysicalParameters.ElectronDensity") != values.end()) config.electronDensity = string_to_double(values["PhysicalParameters.ElectronDensity"],"PhysicalParameters.ElectronDensity",'p'); 
    else throw std::runtime_error("PhysicalParameters.ElectronDensity needed");

    if (values.find("PhysicalParameters.RealPermittivityLocations") != values.end()) config.realPermittivityLocations = parseLocationArray(values["PhysicalParameters.RealPermittivityLocations"], config);
    if (values.find("PhysicalParameters.RealPermittivityValues") != values.end()) config.realPermittivityValues = parseDoubleValueArray(values["PhysicalParameters.RealPermittivityValues"], config);
    if (config.realPermittivityLocations.size() != config.realPermittivityValues.size()) {
        throw std::runtime_error("Number of real Permittivity locations (" + std::to_string(config.realPermittivityLocations.size()) + 
                    ") does not match number of real Permittivity values (" + std::to_string(config.realPermittivityValues.size()) + ").");
    }

    if (values.find("PhysicalParameters.ElectronDensityLocations") != values.end()) config.electronDensityLocations = parseLocationArray(values["PhysicalParameters.ElectronDensityLocations"], config);
    if (values.find("PhysicalParameters.ElectronDensityValues") != values.end()) config.electronDensityValues = parseDoubleValueArray(values["PhysicalParameters.ElectronDensityValues"], config);
    if (config.electronDensityLocations.size() != config.electronDensityValues.size()) {
        throw std::runtime_error("Number of electron Density locations (" + std::to_string(config.electronDensityLocations.size()) + 
                    ") does not match number of electron Density values (" + std::to_string(config.electronDensityValues.size()) + ").");
    }

    if (values.find("PhysicalParameters.RealMobilityLocations") != values.end()) config.realMobilityLocations = parseLocationArray(values["PhysicalParameters.RealMobilityLocations"], config);
    if (values.find("PhysicalParameters.RealMobilityValues") != values.end()) config.realMobilityValues = parseDoubleValueArray(values["PhysicalParameters.RealMobilityValues"], config);
    if (config.realMobilityLocations.size() != config.realMobilityValues.size()) {
        throw std::runtime_error("Number of real mobility locations (" + std::to_string(config.realMobilityLocations.size()) + 
                    ") does not match number of real mobility values (" + std::to_string(config.realMobilityValues.size()) + ").");
    }

    if (values.find("PhysicalParameters.ImagMobilityLocations") != values.end()) config.imagMobilityLocations = parseLocationArray(values["PhysicalParameters.ImagMobilityLocations"], config);
    if (values.find("PhysicalParameters.ImagMobilityValues") != values.end()) config.imagMobilityValues = parseDoubleValueArray(values["PhysicalParameters.ImagMobilityValues"], config);
    if (config.imagMobilityLocations.size() != config.imagMobilityValues.size()) {
        throw std::runtime_error("Number of imag mobility locations (" + std::to_string(config.imagMobilityLocations.size()) + 
                    ") does not match number of imag mobility values (" + std::to_string(config.imagMobilityValues.size()) + ").");
    }

    
    // Parse solver parameters
    if (values.find("Solver.Method") != values.end()) config.solverMethod = values["Solver.Method"]; // in the future catch non valid solver methods here
    else throw std::runtime_error("Solver.Method needed");

    if (values.find("Solver.PC") != values.end()) config.PC = values["Solver.PC"]; // preconditioner, can be empty for direct solvers

    if (values.find("Solver.Scalar") != values.end()) config.scalar = string_to_int(values["Solver.Scalar"],"Solver.Scalar");

    if (values.find("Solver.Tolerance") != values.end()) config.tolerance = string_to_double(values["Solver.Tolerance"],"Solver.Tolerance",'p');
    else if(config.solverMethod != "direct_LU") throw std::runtime_error("Solver.Tolerance needed for iterative solvers");

    if (values.find("Solver.MaxIterations") != values.end()) config.maxIterations = string_to_int(values["Solver.MaxIterations"],"Solver.MaxIterations",'p');
    else if(config.solverMethod != "direct_LU") throw std::runtime_error("Solver.MaxIterations needed for iterative solvers");

    if (values.find("Solver.InitialGuessFile") != values.end()) config.initial_guess_file = values["Solver.InitialGuessFile"];
    
    // Parse boundary conditions
    if (values.find("BoundaryConditions.X_lower_Boundary") != values.end()) config.xLowerBC = values["BoundaryConditions.X_lower_Boundary"];
    if (values.find("BoundaryConditions.X_upper_Boundary") != values.end()) config.xUpperBC = values["BoundaryConditions.X_upper_Boundary"];
    if (values.find("BoundaryConditions.Y_lower_Boundary") != values.end()) config.yLowerBC = values["BoundaryConditions.Y_lower_Boundary"];
    if (values.find("BoundaryConditions.Y_upper_Boundary") != values.end()) config.yUpperBC = values["BoundaryConditions.Y_upper_Boundary"];
    if (values.find("BoundaryConditions.Z_lower_Boundary") != values.end()) config.zLowerBC = values["BoundaryConditions.Z_lower_Boundary"];
    if (values.find("BoundaryConditions.Z_upper_Boundary") != values.end()) config.zUpperBC = values["BoundaryConditions.Z_upper_Boundary"];
    
    const bool hasExcitation_BC = 
    config.xLowerBC == "Excitation" || config.yLowerBC == "Excitation" || config.zLowerBC == "Excitation" ||
    config.xUpperBC == "Excitation" || config.yUpperBC == "Excitation" || config.zUpperBC == "Excitation";  
    if (values.find("BoundaryConditions.Excitation_Value") != values.end()) config.excitationValue = parseComplexValueRow(values["BoundaryConditions.Excitation_Value"]);
    else if (hasExcitation_BC) throw std::runtime_error("Excitation boundary condition requires Fields values (BoundaryConditions.Excitation_Value) to be specified");

    const bool hasRobin_BC = 
    config.xLowerBC == "Robin" || config.yLowerBC == "Robin" || config.zLowerBC == "Robin" ||
    config.xUpperBC == "Robin" || config.yUpperBC == "Robin" || config.zUpperBC == "Robin";  
    if (values.find("BoundaryConditions.Robin_Values") != values.end()){
        config.injectionValues = parseComplexValueArray(values["BoundaryConditions.Robin_Values"], config);
        
        std::vector<std::vector<Complex>> transposed(3, std::vector<Complex>(config.injectionValues.size()));
        for(unsigned int p = 0; p < config.injectionValues.size(); ++p)  for(int c = 0; c < 3; ++c) transposed[c][p] = config.injectionValues[p][c];
        config.injectionValues = transposed;
    } else if(hasRobin_BC) throw std::runtime_error("Robin boundary condition requires Fields values (BoundaryConditions.Robin_Values) to be specified");

    if (values.find("BoundaryConditions.SourceLocations") != values.end()) config.sourceLocations = parseLocationArray(values["BoundaryConditions.SourceLocations"], config);
    if (values.find("BoundaryConditions.SourceValues") != values.end()) config.sourceValues = parseComplexValueArray(values["BoundaryConditions.SourceValues"], config);
    if (config.sourceLocations.size() != config.sourceValues.size()) {
        throw std::runtime_error("Number of source locations (" + std::to_string(config.sourceLocations.size()) + 
                    ") doesn't match number of source values (" + std::to_string(config.sourceValues.size()) + ").");
    }
    
    if (values.find("BoundaryConditions.ExcitationLocations") != values.end()) config.excitationLocations = parseLocationArray(values["BoundaryConditions.ExcitationLocations"], config);
    if (values.find("BoundaryConditions.ExcitationValues") != values.end()) config.excitationValues = parseComplexValueArray(values["BoundaryConditions.ExcitationValues"], config);
    if (config.excitationLocations.size() != config.excitationValues.size()) {
        throw std::runtime_error("Number of excitation locations (" + std::to_string(config.excitationLocations.size()) + 
                    ") doesn't match number of excitation values (" + std::to_string(config.excitationValues.size()) + ").");
    }

    if (values.find("BoundaryConditions.MetalLocations") != values.end()) config.metalLocations = parseLocationArray(values["BoundaryConditions.MetalLocations"], config);
    if (values.find("BoundaryConditions.yReflector") != values.end()) config.y_reflector = string_to_double(values["BoundaryConditions.yReflector"], "BoundaryConditions.yReflector", 'p');

    // Parse PML parameters
    bool PMLExists = std::any_of(values.begin(), values.end(), [](const auto& pair) { return pair.first.find("PML.") == 0;});  
    if (PMLExists){
        if (values.find("PML.Enable_X_lower") != values.end()) config.enableXLowerPML = parseBoolean(values["PML.Enable_X_lower"], "PML.Enable_X_lower");
        if (values.find("PML.Enable_X_upper") != values.end()) config.enableXUpperPML = parseBoolean(values["PML.Enable_X_upper"], "PML.Enable_X_upper");
        if (values.find("PML.Enable_Y_lower") != values.end()) config.enableYLowerPML = parseBoolean(values["PML.Enable_Y_lower"], "PML.Enable_Y_lower");
        if (values.find("PML.Enable_Z_upper") != values.end()) config.enableYUpperPML = parseBoolean(values["PML.Enable_Y_upper"], "PML.Enable_Y_upper");
        if (values.find("PML.Enable_Z_lower") != values.end()) config.enableZLowerPML = parseBoolean(values["PML.Enable_Z_lower"], "PML.Enable_Z_lower");
        if (values.find("PML.Enable_Z_upper") != values.end()) config.enableZUpperPML = parseBoolean(values["PML.Enable_Z_upper"], "PML.Enable_Z_upper");

        if(config.enableXLowerPML) config.xLowerLayers = string_to_int(values["PML.X_lower_Layers"],"PML.X_lower_Layers",'p');
        if(config.enableXUpperPML) config.xUpperLayers = string_to_int(values["PML.X_upper_Layers"],"PML.X_upper_Layers",'p');
        if(config.enableYLowerPML) config.yLowerLayers = string_to_int(values["PML.Y_lower_Layers"],"PML.Y_lower_Layers",'p');
        if(config.enableYUpperPML) config.yUpperLayers = string_to_int(values["PML.Y_upper_Layers"],"PML.Y_upper_Layers",'p');
        if(config.enableZLowerPML) config.zLowerLayers = string_to_int(values["PML.Z_lower_Layers"],"PML.Z_lower_Layers",'p');
        if(config.enableZUpperPML) config.zUpperLayers = string_to_int(values["PML.Z_upper_Layers"],"PML.X_upper_Layers",'p');

        config.orderPML = string_to_int(values["PML.Order"],"PML.Order",'p');
        config.sigma_0 = string_to_double(values["PML.Sigma_0"],"PML.Sigma_0",'p');
    } 

    // Parse output parameters
    if (values.find("Output.PrintMatrix") != values.end()) config.printMatrix = parseBoolean(values["Output.PrintMatrix"], "Output.PrintMatrix");
    if (values.find("Output.PrintProgress") != values.end()) config.printProgress = parseBoolean(values["Output.PrintProgress"], "Output.PrintProgress");
    if (values.find("Output.PrintSolution") != values.end()) config.printSolution = parseBoolean(values["Output.PrintSolution"], "Output.PrintSolution");
    if (values.find("Output.StoreResults") != values.end()) config.storeResults = parseBoolean(values["Output.StoreResults"], "Output.StoreResults");
    if (values.find("Output.PrintConfig") != values.end()) config.printConfig = parseBoolean(values["Output.PrintConfig"], "Output.PrintConfig");
    if (values.find("Output.PrintGrid") != values.end()) config.printGrid = parseBoolean(values["Output.PrintGrid"], "Output.PrintGrid");

    if (values.find("Output.OutputDirectory") != values.end()) config.outputDirectory = values["Output.OutputDirectory"];
    else if(config.storeResults) throw std::runtime_error("Output.OutputDirectory needed for storing results");

    if (values.find("Output.Fields") != values.end()) config.outputFields = split(values["Output.Fields"], ',');
    else if(config.storeResults) throw std::runtime_error("Output.Fields needed for storing results (Options: AbsorbedPower, FieldAmplitudes, FieldRealPart, FieldComponents, Inputs)");

    return config;
}
    

std::string ConfigParser::trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, (last - first + 1));
}

std::vector<std::string> ConfigParser::split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);

    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(trim(token));
    }
    return tokens;
}
// Optimized parseLocationArray with file streaming
std::vector<std::tuple<int, int, int>> ConfigParser::parseLocationArray(const std::string& str, const HARPSConfig& config) {
    std::vector<std::tuple<int, int, int>> locations;
    
    // Check if it's a file reference
    if (str.find(".dat") != std::string::npos) {
        return parseLocationArrayFromFile(str, config);
    }
    
    // Parse inline string (original logic for small strings)
    if (str.empty()) return locations;
    
    std::string content = str.substr(1, str.length() - 2);
    
    std::vector<std::string> pointStrings;
    size_t start = 0;
    while (start < content.length()) {
        size_t openParen = content.find('(', start);
        size_t closeParen = content.find(')', openParen);
        if (openParen == std::string::npos || closeParen == std::string::npos) break;
        
        std::string pointStr = content.substr(openParen + 1, closeParen - openParen - 1);
        pointStrings.push_back(pointStr);
        start = closeParen + 1;
    }
    
    // Parse each point
    for (const auto& pointStr : pointStrings) {
        locations.push_back(parseLocationTuple(pointStr, config));
    }

    return locations;
}

std::vector<std::tuple<int, int, int>> ConfigParser::parseLocationArrayFromFile(const std::string& filename, const HARPSConfig& config) {
    std::ifstream file(harps_dir + filename);
    if (!file.is_open()) {
        throw std::invalid_argument("File(" + harps_dir + filename + "): Cannot be opened");
    }

    std::vector<std::tuple<int, int, int>> locations;
    locations.reserve(10000);
    
    std::string line;
    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty()) continue;
        if (!std::any_of(line.begin(), line.end(), ::isdigit)) continue;

        // Each line should be a tuple: (x,y,z)
        size_t openParen = line.find('(');
        size_t closeParen = line.find(')');
        
        if (openParen == std::string::npos || closeParen == std::string::npos) continue;
        
        std::string pointStr = line.substr(openParen + 1, closeParen - openParen - 1);
        locations.push_back(parseLocationTuple(pointStr, config));
    }
    
    return locations;
}

std::tuple<int, int, int> ConfigParser::parseLocationTuple(const std::string& pointStr, const HARPSConfig& config) {
    std::vector<std::string> coords = split(pointStr, ',');
    if (coords.size() != 3) {
        throw std::runtime_error("Invalid point format: " + pointStr);
    }

    int coord_x = parseExpression(trim(coords[0]), config);
    int coord_y = parseExpression(trim(coords[1]), config);
    int coord_z = parseExpression(trim(coords[2]), config);

    if (coord_x < 0 || coord_x > (config.n_x - 1)) {
        throw std::runtime_error("Coordinate Value i needs to be smaller than Nx(" + 
                               std::to_string(config.n_x) + "): " + std::to_string(coord_x));
    }
    if (coord_y < 0 || coord_y > (config.n_y - 1)) {
        throw std::runtime_error("Coordinate Value j needs to be smaller than Ny(" + 
                               std::to_string(config.n_y) + "): " + std::to_string(coord_y));
    }
    if (coord_z < 0 || coord_z > (config.n_z - 1)) {
        throw std::runtime_error("Coordinate Value k needs to be smaller than Nz(" + 
                               std::to_string(config.n_z) + "): " + std::to_string(coord_z));
    }

    return std::make_tuple(coord_x, coord_y, coord_z);
}

std::vector<double> ConfigParser::parseDoubleValueArray(const std::string& str, const HARPSConfig& config) {
    std::vector<double> doubleValues;

    // Check if it's a file reference
    if (str.find(".dat") != std::string::npos) {
        return parseDoubleValueArrayFromFile(str);
    }
    
    // Parse inline string (original logic for small strings)
    if (str.empty()) return doubleValues;

    std::string content = str.substr(1, str.length() - 2);
    
    std::vector<std::string> valueStrings;
    size_t start = 0;
    while (start < content.length()) {
        size_t openParen = content.find('(', start);
        size_t closeParen = content.find(')', openParen);
        if (openParen == std::string::npos || closeParen == std::string::npos) break;
        
        std::string pointStr = content.substr(openParen + 1, closeParen - openParen - 1);
        valueStrings.push_back(pointStr);
        start = closeParen + 1;
    }

    for (size_t i = 0; i < valueStrings.size(); i++) {
        doubleValues.push_back(string_to_double(valueStrings[i], "double from array", 'a'));
    }

    return doubleValues;
}

std::vector<double> ConfigParser::parseDoubleValueArrayFromFile(const std::string& filename) {
    std::ifstream file(harps_dir + filename);
    if (!file.is_open()) {
        throw std::invalid_argument("File(" + harps_dir + filename + "): Cannot be opened");
    }

    std::vector<double> doubleValues;
    doubleValues.reserve(10000);
    
    std::string line;
    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty()) continue;
        if (!std::any_of(line.begin(), line.end(), ::isdigit)) continue;
        
        // Each line should be a value: (value)
        size_t openParen = line.find('(');
        size_t closeParen = line.find(')');
        
        if (openParen == std::string::npos || closeParen == std::string::npos) continue;
        
        std::string valueStr = line.substr(openParen + 1, closeParen - openParen - 1);
        doubleValues.push_back(string_to_double(trim(valueStr), "double from array", 'a'));
    }
    
    return doubleValues;
}

std::vector<Complex> ConfigParser::parseComplexValueRow(const std::string& str) {
    std::vector<Complex> row;
    
    size_t pos = 0;
    size_t len = str.length();
    
    while (pos < len && str[pos] != '[') pos++;
    if (pos >= len) return row;
    pos++;
    
    while (pos < len) {
        // Skip whitespace and separators
        while (pos < len && (str[pos] == ' ' || str[pos] == ';')) pos++;
        if (pos >= len || str[pos] == ']') break;

        // 1. Handle "F"
        if (str[pos] == 'F') {
            row.emplace_back(Constants::FREE);  // or whatever FREE means
            pos++;
            continue;
        }

        // 2. Handle Complex(...)
        if (pos + 8 < len && str.compare(pos, 8, "Complex(") == 0) {
            pos += 8; // move past "Complex("

            char* end_ptr;

            double real = std::strtod(str.c_str() + pos, &end_ptr);
            pos = end_ptr - str.c_str();

            // skip spaces and comma
            while (pos < len && (str[pos] == ' ' || str[pos] == ',')) pos++;

            double imag = std::strtod(str.c_str() + pos, &end_ptr);
            pos = end_ptr - str.c_str();

            row.emplace_back(real, imag);

            // skip to closing ')'
            while (pos < len && str[pos] != ')') pos++;
            if (pos < len) pos++;

            continue;
        }

        // Unknown token → skip or break
        break;
    }

    return row;
}

std::vector<std::vector<Complex>> ConfigParser::parseComplexValueArray(const std::string& str, const HARPSConfig& config) {
    
    std::vector<std::vector<Complex>> complexValues;

    if (str.find(".dat") != std::string::npos) {
        return parseComplexValueArrayFromFile(str);
    }
    
    if (str.empty()) return complexValues;
    
    std::string content = str.substr(1, str.length() - 2);
    
    size_t start = 0;
    while (start < content.length()) {
        size_t openParen = content.find('[', start);
        size_t closeParen = content.find(']', openParen);
        if (openParen == std::string::npos || closeParen == std::string::npos) break;
        
        std::string pointStr = content.substr(openParen, closeParen - openParen + 1);
        complexValues.push_back(parseComplexValueRow(pointStr));
        start = closeParen + 1;
    }
    
    return complexValues;
}

std::vector<std::vector<Complex>> ConfigParser::parseComplexValueArrayFromFile(const std::string& filename) {
    std::ifstream file(harps_dir + filename);
    if (!file.is_open()) {
        throw std::invalid_argument("File(" + harps_dir + filename + "): Cannot be opened");
    }

    std::vector<std::vector<Complex>> complexValues;
    complexValues.reserve(10000);
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (!std::any_of(line.begin(), line.end(), ::isdigit)) continue;

        complexValues.push_back(parseComplexValueRow(line));
    }
    
    return complexValues;
}

int ConfigParser::parseExpression(const std::string& expr, const HARPSConfig& config) {
    std::string expression = trim(expr);
    
    if (expression.find('/') == std::string::npos) {
        return string_to_int(expression, "Point Coordinates", 'p');
    }
    
    std::vector<std::string> parts = split(expression, '/');
    if (parts.size() != 2) {
        throw std::runtime_error("Invalid expression format: " + expression);
    }
    
    std::string var = trim(parts[0]);
    std::string rightPart = trim(parts[1]);
    
    int baseValue;
    if (var == "Nz") baseValue = config.n_z;
    else if (var == "Nx") baseValue = config.n_x;
    else if (var == "Ny") baseValue = config.n_y;
    else throw std::runtime_error("Unknown variable in expression: " + var);
    
    size_t plusPos = rightPart.find('+');
    size_t minusPos = rightPart.find('-');
    
    if (plusPos == std::string::npos && minusPos == std::string::npos) {
        return baseValue / string_to_int(rightPart, "denominator of N/integer expression");
    } else if (plusPos != std::string::npos) {
        int divisor = string_to_int(trim(rightPart.substr(0, plusPos)), 
                                    "denominator of N/integer+integer expression");
        int addend = string_to_int(trim(rightPart.substr(plusPos + 1)), 
                                   "addition of N/integer+integer expression");
        return (baseValue / divisor) + addend;
    } else {
        int divisor = string_to_int(trim(rightPart.substr(0, minusPos)), 
                                    "denominator of N/integer-integer expression");
        int addend = string_to_int(trim(rightPart.substr(minusPos + 1)), 
                                   "denominator of N/integer-integer expression");
        return (baseValue / divisor) - addend;
    }
}

int ConfigParser::string_to_int(std::string str, std::string key, char cond){
    try {
        std::size_t pos;
        int value = std::stoi(str, &pos);
        if (pos != str.length())  throw std::invalid_argument("Invalid characters in integer (" + key + "): " + str);
        if(cond == 'p' && value<0) throw std::invalid_argument("This integer (" + key + ") should be positive: " + std::to_string(value));
        else if (cond == 'n' && value>0) throw std::invalid_argument("This integer (" + key + ") should be negative: " + std::to_string(value));
        return value;
    } catch (const std::invalid_argument& e) {
        if (std::string(e.what()).find(key) != std::string::npos) throw;
        throw std::invalid_argument("Invalid integer format (" + key + "): " + str);
    } catch (const std::out_of_range& e) {
        throw std::invalid_argument("Integer out of range (" + key + "): " + str);
    }
}

double ConfigParser::string_to_double(std::string str, std::string key, char cond){
    try {
        std::size_t pos;
        double value = std::stod(str, &pos);
        if (pos != str.length())  throw std::invalid_argument("Invalid characters in float (" + key + "): " + str);
        if(cond == 'p' && value<0) throw std::invalid_argument("This float (" + key + ") should be positive: " + std::to_string(value));
        else if (cond == 'n' && value>0) throw std::invalid_argument("This float (" + key + ") should be negative: " + std::to_string(value));
        return value;
    } catch (const std::invalid_argument& e) {
        if (std::string(e.what()).find(key) != std::string::npos) throw;
        throw std::invalid_argument("Invalid float format (" + key + "): " + str);
    } catch (const std::out_of_range& e) {
        throw std::invalid_argument("Float out of range (" + key + "): " + str);
    }
}

bool ConfigParser::parseBoolean(const std::string& value, std::string key) {
    if (value == "TRUE" || value == "True" || value == "true" || value == "1" || value == "yes"  || value == "Yes") {
        return true;
    } else if (value == "FALSE" || value == "False" || value == "false" || value == "0" || value == "-1" || value == "no"  || value == "No") {
        return false;
    } else {
        throw std::invalid_argument("Invalid boolean value (" + key + ", please try true or false): " + value);
    }
}


    
std::ostream& operator<<(std::ostream& os, const std::tuple<int, int, int>& t) {
    os << "(" << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << ")";
    return os;
}
    
// Operator<< for HARPSConfig
std::ostream& operator<<(std::ostream& os, const HARPSConfig& config) {
    os << "HARPSConfig {\n";
    os << "  Grid: n_x=" << config.n_x << ", n_y=" << config.n_y << ", n_z=" << config.n_z << "\n";
    os << "  Dimensions: (" << config.lengthX << ", " << config.lengthY << ", " << config.lengthZ << ")\n";
    os << "  Refinement: (" << config.refinementFactor_x << ", " << config.refinementFactor_y << ", " << config.refinementFactor_z << ")\n";
    if(config.refinementFactor_x > 0) os << "  Refinement Positions X: (" << config.x_BL << ", " << config.x_RL << ", " << config.x_RR <<  ", " << config.x_BR << ")\n";
    if(config.refinementFactor_y > 0) os << "  Refinement Positions Y: (" << config.y_BL << ", " << config.y_RL << ", " << config.y_RR <<  ", " << config.y_BR << ")\n";
    if(config.refinementFactor_z > 0) os << "  Refinement Positions Z: (" << config.z_BL << ", " << config.z_RL << ", " << config.z_RR <<  ", " << config.z_BR << ")\n";
    if(config.x_grid.size() > 0) os << "  X grid: " << config.x_grid << "\n";
    if(config.y_grid.size() > 0) os << "  Y grid: " << config.y_grid << "\n";
    if(config.z_grid.size() > 0) os << "  Z grid: " << config.z_grid << "\n";
    os << "  Frequency: " << config.frequency << " Hz\n";
    os << "  Waveguide Number: " << config.waveguide_number << "\n";
    os << "  yCenter: " << config.yCenter << "\n";
    os << "  Cylindrical Plasma Flag: " << config.flag_cylindrical_plasma << "\n";
    os << "  Real Mobility: " << config.realMobility << "\n";
    os << "  Real Mobility Locations: " << config.realMobilityLocations << "\n";
    os << "  Real Mobility Values: " << config.realMobilityValues << "\n";
    os << "  Imag Mobility: " << config.imagMobility << "\n";
    os << "  Imag Mobility Locations: " << config.imagMobilityLocations << "\n";
    os << "  Imag Mobility Values: " << config.imagMobilityValues << "\n";
    os << "  Real Permittivity: " << config.realPermittivity << "\n";
    os << "  Real Permittivity Locations: " << config.realPermittivityLocations << "\n";
    os << "  Real Permittivity Values: " << config.realPermittivityValues << "\n";
    os << "  Electron Density: " << config.electronDensity << "\n";
    os << "  Electron Density Locations: " << config.electronDensityLocations << "\n";
    os << "  Electron Density Values: " << config.electronDensityValues << "\n";
    os << "  Solver: " << config.solverMethod << ", PC=" << config.PC << ", scalar=" << config.scalar << ", tolerance=" << config.tolerance << ", maxIterations=" << config.maxIterations << "Initial Guess File=" << config.initial_guess_file << "\n";
    os << "  Boundary Conditions: X(" << config.xLowerBC << ", " << config.xUpperBC << "), Y(" << config.yLowerBC << ", " << config.yUpperBC << "), Z(" << config.zLowerBC << ", " << config.zUpperBC << ")\n";
    os << "  Excitation BC Value: " << config.excitationValue << "\n";
    os << "  Source Locations: " << config.sourceLocations << "\n";
    os << "  Source Values: " << config.sourceValues << "\n";
    os << "  Excitation Locations: " << config.excitationLocations << "\n";
    os << "  Excitation Values: " << config.excitationValues << "\n";
    os << "  Metal (in volume) Locations: " << config.metalLocations << "\n";
    os << "  y Reflector: " << config.y_reflector << "\n";
    os << "  PML Enabled: X(" << config.enableXLowerPML << ", " << config.enableXUpperPML << "), Y(" << config.enableYLowerPML << ", " << config.enableYUpperPML << "), Z(" << config.enableZLowerPML << ", " << config.enableZUpperPML << ")\n";
    os << "  PML Layers: X(" << config.xLowerLayers << ", " << config.xUpperLayers << "), Y(" << config.yLowerLayers << ", " << config.yUpperLayers << "), Z(" << config.zLowerLayers << ", " << config.zUpperLayers << ")\n";
    os << "  PML Order: " << config.orderPML << ", PML Sigma_0: " << config.sigma_0 << "\n";
    os << "  Output: printMatrix=" << config.printMatrix << ", printProgress=" << config.printProgress << ", printSolution=" << config.printSolution << ", storeResults=" << config.storeResults << ", printConfig=" << config.printConfig << ", printGrid=" << config.printGrid << "\n";
    os << "  Output Directory: " << config.outputDirectory << "\n";
    os << "  Output Fields: " << config.outputFields << "\n";
    os << "}";
    return os;
}









