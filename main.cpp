#include <iostream>
#include "header.h"

int main() {
    //Configuracion variables
    const double alphaAbsorption = 0.2;
    const double deltaDiffusion = 0.2;
    const int numberOfRays = 12;
    const double initialRayEnergy = 100/numberOfRays;
    const int roomEdgeSize=10;
    // Configuración de la fuente de rayos
    Point sourcePosition(-1.5, -1.5, -1.5);

    Source source(sourcePosition, numberOfRays);
    source.createRays();


    // Configuración de la habitación (cubo)
    Room room(roomEdgeSize);

    // Realizar el trazado de rayos
    room.TraceRays(source, alphaAbsorption, deltaDiffusion, initialRayEnergy);

    // Mostrar los resultados
    room.showResults();
    room.saveResultsToFile();

    return 0;
}
