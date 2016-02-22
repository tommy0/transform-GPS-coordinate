#ifndef COORDINATES_H
#define COORDINATES_H
#include <vector>
#include <cmath>
#define PI 3.14159265359

enum CorrType{RectWGS, GeodesicWGS, RectPE, GeodesicPE};

typedef std::vector<double> Coor;

class Coordinates
{
    public:
        Coordinates(double first, double second, double third, CorrType Type);
        Coor getRectCoor(CorrType type);
        Coor getGeodesicCoor(CorrType type);
        double getDistRect();
    private:
        double X_WGS;
        double Y_WGS;
        double Z_WGS;
        double B_WGS;
        double L_WGS;
        double H_WGS;
        double X_PE;
        double Y_PE;
        double Z_PE;
        double B_PE;
        double L_PE;
        double H_PE;
        void WGS_To_PE();
        void PE_To_WGS();
        void RectToGeodesic(CorrType type);
        void GeodesicToRect(CorrType type);
};

#endif // COORDINATES_H
