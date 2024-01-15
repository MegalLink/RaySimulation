#include <iostream>
#include "header.h"

int main() {
    //Configuracion variables
    const double alphaAbsorption = 0.2;
    const double deltaDiffusion = 0.2;
    const double initialRayEnergy = 100;

    // Configuración de la fuente de rayos
    Point sourcePosition(-1.5, -1.5, -1.5);
    int numberOfRays = 12;
    Source source(sourcePosition, numberOfRays);
    source.createRays();  // Suponiendo que esta función inicializa los rayos correctamente

    // Configuración de la habitación (cubo)
    Room room;
    // Definir los planos que forman el cubo.
    // Cada plano se define con cuatro puntos (puede requerir triangulación si es necesario)
    room.planes[0] = Plane(Point(-2, -2, -2), Point(-2, -2, 2), Point(-2, 2, 2), Point(-2, 2, -2));
    room.planes[1] = Plane(Point(2, -2, -2), Point(2, -2, 2), Point(2, 2, 2), Point(2, 2, -2));
    room.planes[2] = Plane(Point(-2, -2, -2), Point(-2, -2, 2), Point(2, -2, 2), Point(2, -2, -2));
    room.planes[3] = Plane(Point(-2, 2, -2), Point(-2, 2, 2), Point(2, 2, 2), Point(2, 2, -2));
    room.planes[4] = Plane(Point(-2, -2, -2), Point(-2, 2, -2), Point(2, 2, -2), Point(2, -2, -2));
    room.planes[5] = Plane(Point(-2, -2, 2), Point(-2, 2, 2), Point(2, 2, 2), Point(2, -2, 2));

    // Realizar el trazado de rayos
    room.TraceRays(source, alphaAbsorption, deltaDiffusion, initialRayEnergy);

    // Mostrar los resultados
    room.showResults();
    room.saveResultsToFile();

    return 0;
}
