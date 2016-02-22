#include <iostream>
#include "Coordinates.h"

using namespace std;

int main()
{
    Coordinates coordinat(127,-170,30000,GeodesicPE);
    Coor result = coordinat.getGeodesicCoor(GeodesicWGS);
    double dr=coordinat.getDistRect();
    std::cout<<dr<<std::endl;
    for(unsigned int i=0; i<result.size();++i)
    {
        std::cout<<result[i]<<" ";
    }
    return 0;
}
