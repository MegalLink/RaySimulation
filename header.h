#ifndef PROYECTOSIM_HEADER_H
#define PROYECTOSIM_HEADER_H

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

//La clase Point representa un punto en un espacio tridimensional:
class Point {
public:
    double x, y, z;

    Point() : x(0), y(0), z(0) {}

    Point(double x, double y, double z) : x(x), y(y), z(z) {}

    // Sobrecarga de operadores si es necesario
    Point operator-(const Point &other) const {
        return Point(x - other.x, y - other.y, z - other.z);
    }
};


//La clase Vector representa un vector en el espacio tridimensional y
//ofrece operaciones comunes como el cálculo del módulo, producto punto, etc.
class Vector {
public:
    double x, y, z;

    Vector() : x(0), y(0), z(0) {}

    Vector(double x, double y, double z) : x(x), y(y), z(z) {}

    double modulo() const {
        return sqrt(x * x + y * y + z * z);
    }

    Vector operator+(const Vector &v) const {
        return Vector(x + v.x, y + v.y, z + v.z);
    }

    Vector operator-(const Vector &v) const {
        return Vector(x - v.x, y - v.y, z - v.z);
    }

    double dot(const Vector &v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    Vector cross(const Vector &v) const {
        return Vector(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }

    Vector &normalize() {
        double m = modulo();
        if (m > 0) {
            x /= m;
            y /= m;
            z /= m;
        }
        return *this;
    }

    Vector operator*(double scalar) const {
        return Vector(x * scalar, y * scalar, z * scalar);
    }

    Vector &operator+=(const Vector &v) {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    // Constructor que toma un Point
    Vector(const Point &p) : x(p.x), y(p.y), z(p.z) {}
};

//La clase Triangle representa un triángulo en el espacio 3D.
//Aquí reutilizo la definición que proporcionaste, adaptándola a un estilo más convencional de C++.
class Triangle {
public:
    Point p0, p1, p2;
    Vector normal;
    int ID;

    Triangle() : p0(), p1(), p2(), normal(), ID(0) {}

    Triangle(const Point &p0, const Point &p1, const Point &p2, int ID)
            : p0(p0), p1(p1), p2(p2), ID(ID) {
        calculateNormal();
    }

    void calculateNormal() {
        Vector v0 = p1 - p0;
        Vector v1 = p2 - p0;
        normal = v0.cross(v1).normalize();
    }
};

//La clase Plane representa un plano definido por puntos y una normal.
class Plane {
public:
    Point points[4]; // Asumiendo que el plano se define con cuatro puntos
    Vector normal;
    vector<Triangle> triangles;

    Plane() {}

    Plane(const Point &p0, const Point &p1, const Point &p2, const Point &p3) {
        points[0] = p0;
        points[1] = p1;
        points[2] = p2;
        points[3] = p3;
        calculateNormal();
        // Aquí también puedes inicializar los triángulos si es necesario
    }

    void calculateNormal() {
        Vector v0 = points[1] - points[0];
        Vector v1 = points[2] - points[0];
        normal = v0.cross(v1).normalize();
    }
};

class ResidualEnergy {
public:
    double energy;
    int planeIndex;
    int triangleIndex;
    Vector reflectedDirection;

    ResidualEnergy(double energy, int planeIndex, int triangleIndex, const Vector &reflectedDirection)
            : energy(energy), planeIndex(planeIndex), triangleIndex(triangleIndex),
              reflectedDirection(reflectedDirection) {}
};

// La clase Source representa una fuente de rayos.
class Source {
public:
    Point p;
    Vector *Rays;
    int NRays;

    Source() : p(), Rays(nullptr), NRays(0) {}

    Source(const Point &p, int nRays) : p(p), NRays(nRays) {
        Rays = new Vector[nRays];
        // Aquí debes implementar la lógica para inicializar los rayos
    }

    ~Source() {
        delete[] Rays;
    }

    // Método para crear rayos en la direccion de un icosaedro
    void createRays() {
        const double phi = M_PI * (3.0 - sqrt(5.0)); // Golden Ratio

        for (int i = 0; i < NRays; ++i) {
            double y = 1 - (i / static_cast<double>(NRays - 1)) * 2;  // Range from -1 to 1
            double radius = sqrt(1 - y * y);
            double theta = phi * i;

            double x = cos(theta) * radius;
            double z = sin(theta) * radius;

            Rays[i] = Vector(x, y, z);
            Rays[i].normalize(); // Make sure the rays are normalized
        }
    }
};


// La clase Room representa la habitación, es decir, el espacio en el que se realizará el trazado de rayos.
class Room {
public:
    vector<vector<ResidualEnergy>> residualEnergies; // Almacenamiento de energías residuales
    Plane planes[6]; // Asumiendo que la habitación es un cubo

public:
    //Plane planes[6]; // Asumiendo que la habitación es un cubo

    Room() {}

    Room(const Plane &p0, const Plane &p1, const Plane &p2,
         const Plane &p3, const Plane &p4, const Plane &p5) {
        planes[0] = p0;
        planes[1] = p1;
        planes[2] = p2;
        planes[3] = p3;
        planes[4] = p4;
        planes[5] = p5;
    }


    void TraceRays(Source &source, double alpha, double delta, double initialEnergy) {
        const int maxReflections = 50;
        residualEnergies.resize(source.NRays); // Asegurarse de que tiene el tamaño adecuado


        // Ejemplo de cómo inicializar los triángulos para cada plano:
        for (int i = 0; i < 6; ++i) {
            // Suponiendo que cada plano se divide en dos triángulos para simplificar
            planes[i].triangles.push_back(
                    Triangle(planes[i].points[0], planes[i].points[1], planes[i].points[2], i * 2));
            planes[i].triangles.push_back(
                    Triangle(planes[i].points[2], planes[i].points[3], planes[i].points[0], i * 2 + 1));
        }

        for (int r = 0; r < source.NRays; ++r) {
            Vector rayDirection = source.Rays[r]; // Asumiendo que 'Rays' es un array de 'Vector'
            Point rayOrigin = source.p;
            double energy = initialEnergy;
            int reflections = 0;

            while (reflections < maxReflections && energy > 0) {
                double closestDistance = numeric_limits<double>::max();
                Point intersectionPoint;
                int hitPlaneIndex = -1;
                int hitTriangleIndex = -1;

                // Buscar la intersección más cercana
                for (int p = 0; p < 6; ++p) {
                    Plane &plane = planes[p];
                    for (int t = 0; t < plane.triangles.size(); ++t) {
                        double distance;
                        Point tempIntersectionPoint;
                        if (RayIntersectsTriangle(rayOrigin, rayDirection, plane.triangles[t], tempIntersectionPoint,
                                                  distance)) {
                            if (distance < closestDistance) {
                                closestDistance = distance;
                                intersectionPoint = tempIntersectionPoint;
                                hitPlaneIndex = p;
                                hitTriangleIndex = t; // Store the triangle index
                            }
                        }
                    }
                }

                if (hitPlaneIndex != -1) {
                    // Actualizar origen del rayo y calcular la dirección reflejada
                    rayOrigin = intersectionPoint;
                    rayDirection = ReflectRay(rayDirection, planes[hitPlaneIndex].normal);

                    // Calcular energía residual
                    energy *= (1 - alpha) * (1 - delta); // Asumiendo que alpha y delta son conocidos
                    residualEnergies[r].push_back(
                            ResidualEnergy(energy, hitPlaneIndex, hitTriangleIndex, rayDirection));

                    ++reflections;
                } else {
                    // Rayo perdido o sin energía
                    break;
                }
            }
        }
    }

    // Opcional: Manejar las energías residuales como se requiera
    const vector<vector<ResidualEnergy>> &getResidualEnergies() const {
        return residualEnergies;
    }

    void showResults() const {
        for (size_t rayIndex = 0; rayIndex < getResidualEnergies().size(); ++rayIndex) {
            cout << "===========================Energias residuales para el rayo: " << rayIndex + 1
                 << "========================:\n";
            for (size_t bounceIndex = 0; bounceIndex < getResidualEnergies()[rayIndex].size(); ++bounceIndex) {
                ResidualEnergy residualEnergy = getResidualEnergies()[rayIndex][bounceIndex];
                cout << " ---->Reflection number: " << bounceIndex + 1 << " Residual Energy:" << residualEnergy.energy
                     << "<-----\n";
                cout << " Plane Index:" << residualEnergy.planeIndex << " Triangle Index:"
                     << residualEnergy.triangleIndex << "\n";
                cout << " Ray Reflex direction:" << printf("x:%f,y:%f,z:%f", residualEnergy.reflectedDirection.x,
                                                           residualEnergy.reflectedDirection.y,
                                                           residualEnergy.reflectedDirection.z) << "\n";
            }
        }
    }

    void saveResultsToFile() const {
        // Residual Energies Files
        for (size_t rayIndex = 0; rayIndex < getResidualEnergies().size(); ++rayIndex) {
            // Create a stringstream to format the file name
            std::ostringstream fileNameStream;
            fileNameStream << R"(D:\Repositories\ProyectoSim\output\ResidualEnergies_Ray_)" << rayIndex << ".txt";
            std::string fileName = fileNameStream.str();
            std::ofstream outputFile(fileName, std::ios::out);

            // Check if the file was opened successfully
            if (!outputFile.is_open()) {
                std::cerr << "Error opening or creating file: " << fileName << std::endl;
                return; // Exit the function due to an error
            }

            for (size_t bounceIndex = 0; bounceIndex < getResidualEnergies()[rayIndex].size(); ++bounceIndex) {
                if (bounceIndex == 0) {
                    outputFile
                            << "index,residualEnergy,planeIndex,triangleIndex,reflectionDirectionX,reflectionDirectionY,reflectionDirectionZ";
                }
                ResidualEnergy residualEnergy = getResidualEnergies()[rayIndex][bounceIndex];
                Vector direction = residualEnergy.reflectedDirection;
                outputFile << bounceIndex << "," << residualEnergy.energy << ","
                           << residualEnergy.planeIndex << "," << residualEnergy.triangleIndex << ","
                           << residualEnergy.reflectedDirection.x << "," << residualEnergy.reflectedDirection.y << ","
                           << residualEnergy.reflectedDirection.z << "\t";
                outputFile << '\n';
            }
            outputFile.close();
            std::cout << "File saved successfully: " << std::endl;
        }
    }

public:
    bool RayIntersectsTriangle(const Point &rayOrigin, const Vector &rayDirection,
                               const Triangle &triangle, Point &intersectionPoint, double &distance) {
        const double EPSILON = 1e-6;

        Vector edge1 = triangle.p1 - triangle.p0;
        Vector edge2 = triangle.p2 - triangle.p0;

        Vector h = rayDirection.cross(edge2);
        double a = edge1.dot(h);

        if (a > -EPSILON && a < EPSILON)
            return false; // Ray is parallel to the triangle plane

        double f = 1.0 / a;
        Vector s = rayOrigin - triangle.p0;
        double u = f * s.dot(h);

        if (u < 0.0 || u > 1.0)
            return false;

        Vector q = s.cross(edge1);
        double v = f * rayDirection.dot(q);

        if (v < 0.0 || u + v > 1.0)
            return false;

        distance = f * edge2.dot(q);
        if (distance > EPSILON) {
            intersectionPoint.x = rayOrigin.x + rayDirection.x * distance;
            intersectionPoint.y = rayOrigin.y + rayDirection.y * distance;
            intersectionPoint.z = rayOrigin.z + rayDirection.z * distance;
            return true;
        }

        return false; // Intersection behind the ray's origin
    }

    Vector ReflectRay(const Vector &ray, const Vector &normal) {
        // Implementar la reflexión del rayo
        // La fórmula es R = D - 2 * (D . N) * N, donde R es el rayo reflejado, D el rayo incidente, y N la normal
        return ray - normal * 2 * ray.dot(normal);
    }
};

#endif //PROYECTOSIM_HEADER_H
