#include "Coordinates.h"

Coordinates::Coordinates(double first, double second, double third, CorrType Type)
{
switch (Type)
{
    case RectWGS:
    {
        X_WGS=first;
        Y_WGS=second;
        Z_WGS=third;
        RectToGeodesic(RectWGS);
        WGS_To_PE();
        RectToGeodesic(RectPE);
        break;
    }
    case RectPE:
    {
        X_PE=first;
        Y_PE=second;
        Z_PE=third;
        RectToGeodesic(RectPE);
        PE_To_WGS();
        RectToGeodesic(RectWGS);
        break;
    }
    case GeodesicWGS:
    {
        B_WGS=(first/180)*PI;
        L_WGS=(second/180)*PI;
        H_WGS=third;
        GeodesicToRect(GeodesicWGS);
        WGS_To_PE();
        RectToGeodesic(RectPE);
        break;
    }
    case GeodesicPE:
    {
        B_PE=(first/180)*PI;
        L_PE=(second/180)*PI;
        H_PE=third;
        GeodesicToRect(GeodesicPE);
        PE_To_WGS();
        RectToGeodesic(RectWGS);
        break;
    }
}
}

void Coordinates::GeodesicToRect(CorrType type)
{
    switch (type)
    {
    case GeodesicPE:
        {
            double a=6378136;
            double alpha=1/298.25784;
            double eSquare=2*alpha-alpha*alpha;
            double N=a/sqrt(1-eSquare*(sin(B_PE)*sin(B_PE)));
            X_PE=(N+H_PE)*cos(B_PE)*cos(L_PE);
            Y_PE=(N+H_PE)*cos(B_PE)*sin(L_PE);
            Z_PE=((1-eSquare)*N + H_PE)*sin(B_PE);
            break;
        }
    case GeodesicWGS:
        {
            double a=6378137;
            double alpha=1/298.257223563;
            double eSquare=2*alpha-alpha*alpha;
            double N=a/sqrt(1-eSquare*(sin(B_WGS)*sin(B_WGS)));
            X_WGS=(N+H_WGS)*cos(B_WGS)*cos(L_WGS);
            Y_WGS=(N+H_WGS)*cos(B_WGS)*sin(L_WGS);
            Z_WGS=((1-eSquare)*N + H_WGS)*sin(B_WGS);
            break;
        }
    }
}

void Coordinates::RectToGeodesic(CorrType type)
{
    switch(type)
    {
    case RectWGS:
        {
            double a=6378137;
            double alpha=1/298.257223563;
            double eSquare=2*alpha-alpha*alpha;
            double D=sqrt(X_WGS*X_WGS + Y_WGS*Y_WGS);
            if (D==0)
            {
                B_WGS=(PI/2)*Z_WGS/fabs(Z_WGS);
                L_WGS=0;
                H_WGS=Z_WGS*sin(B_WGS)-a*sqrt(1-eSquare*sin(B_WGS)*sin(B_WGS));
            }
            if(D>0)
            {
                double La=asin(Y_WGS/D);
                if(Y_WGS<0 && X_WGS>0) L_WGS=2*PI-La;
                if(Y_WGS<0 && X_WGS<0) L_WGS=2*PI+La;
                if(Y_WGS>0 && X_WGS<0) L_WGS=PI-La;
                if(Y_WGS>0 && X_WGS>0) L_WGS=La;
            }
            if(Z_WGS==0)
            {
                B_WGS=0;
                H_WGS=D-a;
            }
            else
            {
                double r=sqrt(X_WGS*X_WGS + Y_WGS*Y_WGS + Z_WGS*Z_WGS);
                double c=asin(Z_WGS/r);
                double p=eSquare*a/(2*r);
                double s1=0, s2=0, d=1,b=0;
                while(d>=(2.78e-8/180)*PI)
                {
                    s1=s2;
                    b=c+s1;
                    s2=asin(p*sin(2*b)/(sqrt(1-eSquare*sin(b)*sin(b))));
                    d=fabs(s2-s1);
                }
                B_WGS=b;
                H_WGS=D*cos(B_WGS)+Z_WGS*sin(B_WGS)-a*sqrt(1-eSquare*sin(B_WGS)*sin(B_WGS));
            }
            break;
        }
    case RectPE:
        {
            double a=6378136;
            double alpha=1/298.25784;
            double eSquare=2*alpha-alpha*alpha;
            double D=sqrt(X_PE*X_PE + Y_PE*Y_PE);
            if (D=0)
            {
                B_PE=(PI/2)*Z_PE/fabs(Z_PE);
                L_PE=0;
                H_PE=Z_PE*sin(B_PE)-a*sqrt(1-eSquare*sin(B_PE)*sin(B_PE));
            }
            if(D>0)
            {
                double La=asin(Y_PE/D);
                if(Y_PE<0 && X_PE>0) L_PE=2*PI-La;
                if(Y_PE<0 && X_PE<0) L_PE=2*PI+La;
                if(Y_PE>0 && X_PE<0) L_PE=PI-La;
                if(Y_PE>0 && X_PE>0) L_PE=La;
            }
            if(Z_PE==0)
            {
                B_PE=0;
                H_PE=D-a;
            }
            else
            {
                double r=sqrt(X_PE*X_PE + Y_PE*Y_PE + Z_PE*Z_PE);
                double c=asin(Z_PE/r);
                double p=eSquare*a/(2*r);
                double s1=0.0, s2=0.0, d=1.0,b=0;
                while(d>=(2.78e-8/180)*PI)
                {
                    s1=s2;
                    b=c+s1;
                    s2=asin(p*sin(2*b)/(sqrt(1-eSquare*sin(b)*sin(b))));
                    d=fabs(s2-s1);
                }
                B_PE=b;
                H_PE=D*cos(B_PE)+Z_PE*sin(B_PE)-a*sqrt(1-eSquare*sin(B_PE)*sin(B_PE));
            }
            break;
        }
    }
}

void Coordinates::WGS_To_PE()
{
    double mat[3][3] ={{1, 0.82e-6,0},{-0.82e-6,1,0},{0, 0, 1}};
    double col[3]={1.1,0.3,0.9};
    double m=(1+0.12e-6);
    double vecWGS[3] = {X_WGS,Y_WGS,Z_WGS};
    double vecPE[3];
    for(int i=0; i<3; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            mat[i][j]*=m;
        }
    }
    for(int i=0; i<3; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            vecPE[i]=mat[i][j]*vecWGS[j];
        }
        vecPE[i]+=col[i];
    }
    X_PE=vecPE[0];
    Y_PE=vecPE[1];
    Z_PE=vecPE[2];
}

void Coordinates::PE_To_WGS()
{
    double mat[3][3] ={{1, -0.82e-6,0},{0.82e-6,1,0},{0, 0, 1}};
    double col[3]={-1.1,-0.3,-0.9};
    double m=(1-0.12e-6);
    double vecPE[3] = {X_PE,Y_PE,Z_PE};
    double vecWGS[3]={0,0,0};
    for(int i=0; i<3; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            mat[i][j]*=m;
        }
    }
    for(int i=0; i<3; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            vecWGS[i]+=mat[i][j]*vecPE[j];
        }
        vecWGS[i]+=col[i];
    }
    X_WGS=vecWGS[0];
    Y_WGS=vecWGS[1];
    Z_WGS=vecWGS[2];
}

Coor Coordinates::getGeodesicCoor(CorrType type)
{
    switch (type)
    {
    case GeodesicPE:
        {
            Coor temp={B_PE*180/PI,L_PE*180/PI,H_PE};
            return temp;
        }
    case GeodesicWGS:
        {
            Coor temp={B_WGS*180/PI,L_WGS*180/PI,H_WGS};
            return temp;
        }
    }
}


Coor Coordinates::getRectCoor(CorrType type)
{
    switch (type)
    {
    case RectPE:
        {
            Coor temp={X_PE,Y_PE,Z_PE};
            return temp;
        }
    case RectWGS:
        {
            Coor temp={X_WGS,Y_WGS,Z_WGS};
            return temp;
        }
    }
}

double Coordinates::getDistRect()
{
    return sqrt((X_WGS-X_PE)*(X_WGS-X_PE) + (Y_WGS-Y_PE)*(Y_WGS-Y_PE) + (Z_WGS-Z_PE)*(Z_WGS-Z_PE));
}
