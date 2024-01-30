#ifndef PROYECTOSIM_HEADER_H
#define PROYECTOSIM_HEADER_H

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

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


//La clase Vector representa un vector en el espacio tridimensional
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

    // producto punto de vectores
    double dot(const Vector &v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    // producto cruz entre vectores
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

class SubPlane {
private:
    void generateTriangles() {
        // Generar automáticamente 2 triángulos en función de los puntos
        triangles.clear(); // Limpiar cualquier triángulo existente

        // Supongamos que los puntos de cada triángulo se generan de alguna manera
        // Aquí, simplemente se crean puntos de ejemplo
        Point p0 = points[0];
        Point p1 = points[1];
        Point p2 = points[2];
        triangles.push_back(Triangle(p0, p1, p2, 0));

        Point p3 = points[2];
        Point p4 = points[3];
        triangles.push_back(Triangle(p2, p3, p4, 1));
    }


public:
    Point points[4];
    Vector normal;
    int subPlaneID;
    vector<Triangle> triangles;

    SubPlane(int subPlaneID, const Point &p0, const Point &p1, const Point &p2, const Point &p3) {
        points[0] = p0;
        points[1] = p1;
        points[2] = p2;
        points[3] = p3;
        subPlaneID = subPlaneID;
        generateTriangles();
    }
};

//La clase Plane representa un plano definido por puntos y una normal.
class Plane {
public:
    Point points[4];
    Vector normal;
    int planeID;
    std::vector<SubPlane> subPlanes;

    Plane(int planeID, const Point &p0, const Point &p1, const Point &p2, const Point &p3) {
        points[0] = p0;
        points[1] = p1;
        points[2] = p2;
        points[3] = p3;
        planeID = planeID;
        calculateNormal();
        createSubPlanes(1);
    }

    void calculateNormal() {
        Vector v0 = points[1] - points[0];
        Vector v1 = points[2] - points[0];
        normal = v0.cross(v1).normalize();
    }

    void createSubPlanes(int horizontalVerticalSlices) {
        // Calculate the distance from p0 to p1
        double distance = Vector(points[1] - points[0]).modulo();

        // Calculate the step size for slices
        double stepSize = distance / horizontalVerticalSlices;

        // Number of subplanes = slices * slices
        int subplaneCount = horizontalVerticalSlices * horizontalVerticalSlices;

        // Generate subplanes
        for (int i = 0; i < horizontalVerticalSlices; ++i) {
            for (int j = 0; j < horizontalVerticalSlices; ++j) {
                double offsetX = static_cast<double>(i) * stepSize / distance;
                double offsetY = static_cast<double>(j) * stepSize / distance;

                Point subplanePoints[4];
                for (int k = 0; k < 4; ++k) {
                    Vector subplaneVector = Vector(points[k]) + normal * offsetX + normal * offsetY;
                    subplanePoints[k] = Point(subplaneVector.x, subplaneVector.y, subplaneVector.z);
                }

                subPlanes.push_back(SubPlane(subplaneCount, subplanePoints[0], subplanePoints[1], subplanePoints[2],
                                             subplanePoints[3]));
                subplaneCount++;
            }
        }
    }
};

class ResidualEnergy {
public:
    double energy;
    int planeIndex;
    int subPlaneIndex;
    int triangleIndex;
    Vector reflectedDirection;

    ResidualEnergy(double energy, int planeIndex, int subPlaneIndex, int triangleIndex,
                   const Vector &reflectedDirection)
            : energy(energy), planeIndex(planeIndex), subPlaneIndex(subPlaneIndex), triangleIndex(triangleIndex),
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
    vector<vector<ResidualEnergy>> residualEnergies;
    std::vector<Plane> planes;

private:
    std::vector<Plane> createCube(double edgeSize) {
        std::vector<Plane> result;

        // Assuming that the cube is centered at the origin
        Point vertices[8] = {
                Point(-edgeSize / 2, -edgeSize / 2, -edgeSize / 2),
                Point(edgeSize / 2, -edgeSize / 2, -edgeSize / 2),
                Point(edgeSize / 2, -edgeSize / 2, edgeSize / 2),
                Point(-edgeSize / 2, -edgeSize / 2, edgeSize / 2),
                Point(-edgeSize / 2, edgeSize / 2, -edgeSize / 2),
                Point(edgeSize / 2, edgeSize / 2, -edgeSize / 2),
                Point(edgeSize / 2, edgeSize / 2, edgeSize / 2),
                Point(-edgeSize / 2, edgeSize / 2, edgeSize / 2)
        };

        // Create the planes of the cube
        result.push_back(Plane(0, vertices[0], vertices[1], vertices[2], vertices[3]));
        result.push_back(Plane(1, vertices[4], vertices[5], vertices[6], vertices[7]));
        result.push_back(Plane(2, vertices[0], vertices[3], vertices[7], vertices[4]));
        result.push_back(Plane(3, vertices[1], vertices[2], vertices[6], vertices[5]));
        result.push_back(Plane(4, vertices[0], vertices[1], vertices[5], vertices[4]));
        result.push_back(Plane(5, vertices[2], vertices[3], vertices[7], vertices[6]));

        return result;
    }

public:
    Room(double edgeSize)
            : residualEnergies(), planes(createCube(edgeSize)) {
    }

    void TraceRays(Source &source, double alpha, double delta, double initialEnergy) {
        const int maxReflections = 50; //TODO THIS MIGH BE ON DEFINE OR PASS IT TO METHOD
        residualEnergies.resize(source.NRays); // Asegurarse de que tiene el tamaño adecuado

        for (int r = 0; r < source.NRays; ++r) {
            Vector rayDirection = source.Rays[r];
            Point rayOrigin = source.p;
            double energy = initialEnergy;
            int reflections = 0;

            while (reflections < maxReflections && energy > 0) {
                double closestDistance = numeric_limits<double>::max();
                Point intersectionPoint;
                int hitPlaneIndex = -1;
                int hitSubPlaneIndex = -1;
                int hitTriangleIndex = -1;

                // Buscar la intersección más cercana
                for (int planeID = 0; planeID < planes.size(); ++planeID) {
                    Plane &plane = planes[planeID];
                    //cout << "Plane ID:" << planeID<< endl;
                    for (int subPlaneID = 0; subPlaneID < plane.subPlanes.size(); ++subPlaneID) {
                        SubPlane &subPlane = plane.subPlanes[subPlaneID];
                        // cout << "Sub Plane ID:" << subPlaneID<< endl;
                        for (int triangleID = 0; triangleID < subPlane.triangles.size(); ++triangleID) {
                            //   cout << "Triangle ID:" << triangleID << endl;
                            double distance;
                            Point tempIntersectionPoint;
                            if (RayIntersectsTriangle(rayOrigin, rayDirection, subPlane.triangles[triangleID],
                                                      tempIntersectionPoint,
                                                      distance)) {
                                if (distance < closestDistance) {
                                    closestDistance = distance;
                                    intersectionPoint = tempIntersectionPoint;
                                    hitPlaneIndex = planeID;
                                    hitSubPlaneIndex = subPlaneID;
                                    hitTriangleIndex = triangleID; // Store the triangle index
                                }
                            }
                        }
                    }
                }

                if (hitTriangleIndex != -1) {
                    // Actualizar origen del rayo y calcular la dirección reflejada
                    rayOrigin = intersectionPoint;
                    rayDirection = ReflectRay(rayDirection, planes[hitSubPlaneIndex].normal);

                    // Calcular energía residual
                    energy *= (1 - alpha) * (1 - delta); // Asumiendo que alpha y delta son conocidos
                    residualEnergies[r].push_back(
                            ResidualEnergy(energy, hitPlaneIndex, hitSubPlaneIndex, hitTriangleIndex, rayDirection));

                    ++reflections;
                } else {
                    // Rayo perdido o sin energía
                    break;
                }
            }
        }
    }

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
                cout << " Plane Index:" << residualEnergy.planeIndex << " SubPlane Index:"
                     << residualEnergy.subPlaneIndex << " Triangle Index:" << residualEnergy.triangleIndex << "\n";
                cout << " Ray Reflex direction:" << printf("x:%f,y:%f,z:%f", residualEnergy.reflectedDirection.x,
                                                           residualEnergy.reflectedDirection.y,
                                                           residualEnergy.reflectedDirection.z) << "\n";
            }
        }
    }

    void saveResultsToFile() const {
        // Get the current working directory on your computer uncoment this line
        //std::filesystem::path currentPath = std::filesystem::current_path();
        const string OUTPUT_PATH = "D:\\Repositories\\ProyectoSim\\output\\";

        // Specify the output folder name
        std::string outputFolder = "output";
        // Residual Energies Files
        for (size_t rayIndex = 0; rayIndex < getResidualEnergies().size(); ++rayIndex) {
            std::ostringstream fileNameStream;
            // uncomment this line on your pc
            //fileNameStream << outputFolder << "/ResidualEnergies_Ray_" << rayIndex << ".txt";
            fileNameStream << OUTPUT_PATH << R"(ResidualEnergies_Ray_)" << rayIndex << ".txt";
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
                            << "index,residualEnergy,planeIndex,subPlaneIndex,triangleIndex,reflectionDirectionX,reflectionDirectionY,reflectionDirectionZ"<<endl;
                }
                ResidualEnergy residualEnergy = getResidualEnergies()[rayIndex][bounceIndex];
                Vector direction = residualEnergy.reflectedDirection;
                outputFile << bounceIndex << "," << residualEnergy.energy << ","
                           << residualEnergy.planeIndex << "," << residualEnergy.subPlaneIndex << ","
                           << residualEnergy.triangleIndex << ","
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

        // Calculate vectors representing edges of the triangle
        Vector triangleEdge1 = triangle.p1 - triangle.p0;
        Vector triangleEdge2 = triangle.p2 - triangle.p0;

        // Calculate the cross product of the ray direction and triangle edge2
        Vector h = rayDirection.cross(triangleEdge2);
        double determinant = triangleEdge1.dot(h);

        // Check if the ray is parallel to the triangle plane
        if (determinant > -EPSILON && determinant < EPSILON)
            return false;

        // Calculate the reciprocal of the determinant
        double determinantReciprocal = 1.0 / determinant;

        // Calculate vectors used for barycentric coordinates
        Vector rayOriginToVertex = rayOrigin - triangle.p0;
        double u = determinantReciprocal * rayOriginToVertex.dot(h);

        // Check if the intersection point is outside the triangle
        if (u < 0.0 || u > 1.0)
            return false;

        Vector q = rayOriginToVertex.cross(triangleEdge1);
        double v = determinantReciprocal * rayDirection.dot(q);

        // Check if the intersection point is outside the triangle
        if (v < 0.0 || u + v > 1.0)
            return false;

        // Calculate the distance to the intersection point
        distance = determinantReciprocal * triangleEdge2.dot(q);

        // Check if the intersection point is in front of the ray
        if (distance > EPSILON) {
            // Calculate the intersection point coordinates
            intersectionPoint.x = rayOrigin.x + rayDirection.x * distance;
            intersectionPoint.y = rayOrigin.y + rayDirection.y * distance;
            intersectionPoint.z = rayOrigin.z + rayDirection.z * distance;
            return true;
        }

        // Intersection behind the ray's origin
        return false;
    }

    Vector ReflectRay(const Vector &ray, const Vector &normal) {
        // Implementar la reflexión del rayo
        // La fórmula es R = D - 2 * (D . N) * N, donde R es el rayo reflejado, D el rayo incidente, y N la normal
        return ray - normal * 2 * ray.dot(normal);
    }
};

#endif //PROYECTOSIM_HEADER_H
