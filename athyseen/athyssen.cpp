#include "stdlib.h"
#include <string>
#include <vector>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <assert.h>
using namespace std;

#define SUBMISSION
//#define FULL_SUBMISSION

#ifndef SUBMISSION
//#define ALLOW_SAVEBMP
#endif

#define PANELWIDTH  (0.22f)
#define PANELHEIGHT (0.32f)

//=========================================================================

class RobonautEye
{
public:
std::vector<string> recognizeObjects( std::vector<int>& leftEyeImage, std::vector<int>& rightEyeImage);
};

//=========================================================================

typedef unsigned char         byte;
typedef unsigned short        u16;
typedef short                 s16;
typedef unsigned int          u32;
typedef unsigned long long    u64;
typedef  long long            i64;

#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))

//=========================================================================
// RANDOM
//=========================================================================
static       i64 rand_root = 13;
static const i64 rand_mask =(i64)0xFFFFFFFFFFFFLL;
static const i64 rand_add = (i64)0xBLL;
static const i64 rand_mul = (i64)0x5DEECE66DLL;
inline int myrand()
{
    i64 n = (rand_root*rand_mul + rand_add) & rand_mask;  
    rand_root = n; 
    return int(n>>(48-31));
}

u32 rand0=13;
u32 rand1=13;
#define fastrand() (rand0 = 36969 * (rand0 & 65535) + (rand0 >> 16),rand1 = 18000 * (rand1 & 65535) + (rand1 >> 16),((rand0 << 16) + rand1)&0x7FFFFFFF)


//=========================================================================
// CLOCK
//=========================================================================
static clock_t clock_start;
static double SecsUsed = 0;
static bool Paused = false;
void StartClock() { SecsUsed = 0; Paused=false; clock_start = clock();}
void PauseClock() { if(!Paused) {SecsUsed += (clock()-clock_start) / (double)CLOCKS_PER_SEC;} Paused=true;}
void UnpauseClock() { clock_start = clock(); Paused=false;}
double SecondsRunning() { double T =SecsUsed; if( !Paused ) {T+=(clock()-clock_start) / (double)CLOCKS_PER_SEC;} return T;}

//=========================================================================
// LOGGING
//=========================================================================
#ifdef FULL_SUBMISSION
inline void log(char * format, ...) {(void)format;}
#else
void log(char * format, ...)
{
#ifndef FULL_SUBMISSION
      char buffer[512];
      va_list args;
      va_start (args, format);
      vsprintf (buffer,format, args);
  
    #ifdef SUBMISSION
      fprintf(stdout,"%s",buffer);
      fflush(stdout);
    #else
      static FILE* fp = NULL;
      if( fp==NULL )
      {
          fp = fopen("logtext.txt","wt");
          fclose(fp);
          fp = fopen("logtext.txt","at");
      }
      fprintf(fp,buffer);
      fflush(fp);

      fprintf(stderr,buffer);
      fflush(stderr);
    #endif

      va_end (args);
#endif
}
#endif

void logclock(const char* Str)
{
    static double PrevT = 0;
    double T = SecondsRunning();
    log("<%8.3f,%8.3f> - %s",(float)T,T-PrevT,Str);
    PrevT = T;
}

struct logtime
{
public:    
    double StartT;
    char Title[32];
    logtime( const char* T ) 
    { 
        strncpy(Title,T,31);
        StartT = SecondsRunning(); 
        log(">>>>>>>>>>>>>>>>>>>>>> LOGTIME ENTRY <%8.3f,%8.3f> - %s\n",(float)StartT,(float)0,Title); 
    }
    ~logtime() 
    { 
        double T2=SecondsRunning(); double Secs = T2-StartT; 
        log(">>>>>>>>>>>>>>>>>>>>>> LOGTIME EXIT <%8.3f,%8.3f> - %s\n",(float)T2,(float)Secs,Title); 
    }
};

//=========================================================================
//=========================================================================
//=========================================================================
//=========================================================================

enum ENUM_ENVIRONMENTS
{
    ENVIRON_ISS=0,
    ENVIRON_LAB,
    ENVIRON_LAB2,
    ENVIRON_LAB3,
    ENVIRON_SIM,
};

char* EnvironmentNames[] = 
{
    "ISS",
    "Lab",
    "Lab2",
    "Lab3",
    "Sim"
};

struct EnvironmentInfo
{
    int Directory;
    int Side;
    char* Name;

    float Focus;
    float SX;
    float SY;
    float CX;
    float CY;
    float Kappa;

    float PwrCapRedRatio;

    float CamPosTransfer[3];
};

EnvironmentInfo EnvironmentInfos[5*2] = 
{
    // ISS
    {ENVIRON_ISS,0,"ISS_LEFT",        0.008205f,3.445e-006f,3.450e-006f,1196.4345f,965.3666f,-1186.1950f,  1.40f,-0.140060846356575,0.000755406445758013,0.00437628070103844},
    {ENVIRON_ISS,1,"ISS_RIGHT",        0.008277f,3.452e-006f,3.450e-006f,1191.7501f,984.7040f,-1050.1765f,  1.40f,0.140060846356575,-0.000755406445758013,-0.00437628070103844},
    // Lab
    {ENVIRON_LAB,0,"Lab_LEFT",        0.008205f,3.445e-006f,3.450e-006f,1196.4345f,965.3666f,-1186.1950f,  2.50f,-0.140060846356575,0.000755406445758013,0.00437628070103844},
    {ENVIRON_LAB,1,"Lab_RIGHT",        0.008277f,3.452e-006f,3.450e-006f,1191.7501f,984.7040f,-1050.1765f,  2.50f,0.140060846356575,-0.000755406445758013,-0.00437628070103844},
    // Lab2
    {ENVIRON_LAB2,0,"Lab2_LEFT",    0.008277f,4.90e-006f,4.900e-006f,800.00f,600.000f,0.0f,  2.00f,-0.145246983427213,-0.00143344082412253,-0.000722509379433324},
    {ENVIRON_LAB2,1,"Lab2_RIGHT",    0.008277f,4.90e-006f,4.900e-006f,800.00f,600.000f,0.0f,  2.00f,0.145246983427213,0.00143344082412253,0.000722509379433324},
    // Lab3
    {ENVIRON_LAB3,0,"Lab3_LEFT",    0.008277f,4.90e-006f,4.900e-006f,800.00f,600.000f,0.0f,  1.00f,-0.145789224949895,-0.00171330596683407,0.00268476925603585},
    {ENVIRON_LAB3,1,"Lab3_RIGHT",    0.008277f,4.90e-006f,4.900e-006f,800.00f,600.000f,0.0f,  1.00f,0.145789224949895,0.00171330596683407,-0.00268476925603585},
    // Sim
    {ENVIRON_SIM,0,"Sim_LEFT",        0.00838024245713582,1.04942396718597E-05,1.04942396718597E-05,800.00f,615.126,0.0f,  1.50f,-0.0989986070366767,2.26194914525757E-05,0.000134341771614058},
    {ENVIRON_SIM,1,"Sim_RIGHT",        0.00838024245713582,1.04942396718597E-05,1.04942396718597E-05,800.00f,615.126,0.0f,  1.50f,0.0989986070366765,-2.26194914525762E-05,-0.000134341771614058},
};

//=========================================================================
// 2D GEOMETRY
//=========================================================================
struct Vec2
{
    float X;
    float Y;

    Vec2() {}

    inline Vec2( int x, int y) { X=(float)x; Y=(float)y; }
    inline Vec2( float x, float y) { X=x; Y=y; }
    inline bool        Same    ( const Vec2& P ) const { return (P.X==X) && (P.Y==Y); }
    inline float    Dist2    ( const Vec2& P ) const { return (P.X-X)*(P.X-X) + (P.Y-Y)*(P.Y-Y); }
    inline float    Dist    ( const Vec2& P ) const { return sqrtf( (P.X-X)*(P.X-X) + (P.Y-Y)*(P.Y-Y) ); }
    inline float    Len2    ( void ) const { return (X*X) + (Y*Y); }
    inline float    Len      ( void ) const { return sqrtf( (X*X) + (Y*Y) ); }
    inline Vec2        DeltaTo    ( const Vec2& P ) const { return Vec2(P.X-X,P.Y-Y); }
    inline Vec2        Minus    ( const Vec2& P ) const { return Vec2(X-P.X,Y-P.Y); }
    inline Vec2        Mult    ( float S ) const { return Vec2(X*S,Y*S); }
    inline float    Dot        ( const Vec2& P ) const { return P.X*X + P.Y*Y; }
    inline float    Cross    ( const Vec2& P ) const { return X*P.Y - P.X*Y; }
    inline Vec2        Perp    ( void ) const { return Vec2(-Y,X); }
};

struct Vec3
{
    float X;
    float Y;
    float Z;

    Vec3() {}

    inline Vec3( float x, float y, float z) { X=x; Y=y; Z=z;}
    inline Vec3( Vec2 V) { X=V.X; Y=V.Y; Z=0;}
    inline bool        Same    ( const Vec3& P ) const { return (P.X==X) && (P.Y==Y) && (P.Z==Z); }
    inline float    Dist2    ( const Vec3& P ) const { return (P.X-X)*(P.X-X) + (P.Y-Y)*(P.Y-Y) + (P.Z-Z)*(P.Z-Z); }
    inline float    Dist    ( const Vec3& P ) const { return sqrtf( (P.X-X)*(P.X-X) + (P.Y-Y)*(P.Y-Y) + (P.Z-Z)*(P.Z-Z) ); }
    inline float    Len2    ( void ) const { return (X*X) + (Y*Y) + (Z*Z); }
    inline float    Len      ( void ) const { return sqrtf( (X*X) + (Y*Y) + (Z*Z) ); }
    inline Vec3        DeltaTo    ( const Vec3& P ) const { return Vec3(P.X-X,P.Y-Y,P.Z-Z); }
    inline Vec3        Minus    ( const Vec3& P ) const { return Vec3(X-P.X,Y-P.Y,Z-P.Z); }
    inline Vec3        Mult    ( float S ) const { return Vec3(X*S,Y*S,Z*S); }
    inline float    Dot        ( const Vec3& P ) const { return P.X*X + P.Y*Y + P.Z*Z; }
    inline Vec3        Cross    ( const Vec3& V ) const { return Vec3( (Y*V.Z) - (Z*V.Y), (Z*V.X) - (X*V.Z), (X*V.Y) - (Y*V.X) ); }
    inline Vec3        Normalize    () const { float L = Len(); return Vec3(X/L,Y/L,Z/L); }
};

Vec2 Find2DLineIntersection( Vec2 P1, Vec2 P2, Vec2 P3, Vec2 P4 )
{
    float bx = P2.X - P1.X;
    float by = P2.Y - P1.Y;
    float dx = P4.X - P3.X;
    float dy = P4.Y - P3.Y; 
    float bd = bx*dy - by*dx;
    if(bd == 0) return Vec2(0,0);
    float cx = P3.X-P1.X; 
    float cy = P3.Y-P1.Y;
    float t = (cx*dy - cy*dx) / bd; 
    float X = P1.X+t*bx;
    float Y = P1.Y+t*by;
    return Vec2(X,Y);
}



struct Quad
{
    //  0---3
    //  |   |
    //  1---2
    Vec2 Pts[4];

    void Expand( float T, float B, float L, float R )
    {
        Pts[0].X -= L;
        Pts[0].Y -= T;
        Pts[1].X -= L;
        Pts[1].Y += B;
        Pts[2].X += R;
        Pts[2].Y += B;
        Pts[3].X += R;
        Pts[3].Y -= T;
    }

    Vec2 CalcBilinearPos( float XT, float YT )
    {
        float W00 = (1-XT)*(1-YT);
        float W01 = (1-XT)*(YT);
        float W11 = (XT*YT);
        float W10 = (XT)*(1-YT);

        float X = Pts[0].X*W00 + Pts[1].X*W01 + Pts[2].X*W11 + Pts[3].X*W10;
        float Y = Pts[0].Y*W00 + Pts[1].Y*W01 + Pts[2].Y*W11 + Pts[3].Y*W10;

        return Vec2(X,Y);
    }
};


void CalcBary( Vec2 P, Vec2 A, Vec2 B, Vec2 C, float& wa, float& wb, float& wc )
{
    Vec2 V0( B.X-A.X, B.Y-A.Y);
    Vec2 V1( C.X-A.X, C.Y-A.Y);
    Vec2 V2( P.X-A.X, P.Y-A.Y);
    float d00 = V0.Dot(V0);
    float d01 = V0.Dot(V1);
    float d11 = V1.Dot(V1);
    float d20 = V2.Dot(V0);
    float d21 = V2.Dot(V1);
    float invDenom = 1.0f / (d00 * d11 - d01 * d01);
    wb = (d11 * d20 - d01 * d21) * invDenom;
    wc = (d00 * d21 - d01 * d20) * invDenom;
    wa = 1.0f - wb - wc;
}

void CalcBaryXY( Vec3 P, Vec3 A, Vec3 B, Vec3 C, float& wa, float& wb, float& wc )
{
    Vec2 V0( B.X-A.X, B.Y-A.Y);
    Vec2 V1( C.X-A.X, C.Y-A.Y);
    Vec2 V2( P.X-A.X, P.Y-A.Y);
    float d00 = V0.Dot(V0);
    float d01 = V0.Dot(V1);
    float d11 = V1.Dot(V1);
    float d20 = V2.Dot(V0);
    float d21 = V2.Dot(V1);
    float invDenom = 1.0f / (d00 * d11 - d01 * d01);
    wb = (d11 * d20 - d01 * d21) * invDenom;
    wc = (d00 * d21 - d01 * d20) * invDenom;
    wa = 1.0f - wb - wc;
}

Vec2 CalcBaryPosition( Vec2 A, Vec2 B, Vec2 C, float wa, float wb, float wc )
{
    Vec2 P;
    P.X = A.X*wa + B.X*wb + C.X*wc;
    P.Y = A.Y*wa + B.Y*wb + C.Y*wc;
    return P;
}

Vec3 CalcBaryPosition( Vec3 A, Vec3 B, Vec3 C, float wa, float wb, float wc )
{
    Vec3 P;
    P.X = A.X*wa + B.X*wb + C.X*wc;
    P.Y = A.Y*wa + B.Y*wb + C.Y*wc;
    P.Z = A.Z*wa + B.Z*wb + C.Z*wc;
    return P;
}

Vec2 CalcBaryPosition( int TriIndex, Vec2* pTriPts, float* Weights )
{
    return CalcBaryPosition( pTriPts[TriIndex*3+0],
                             pTriPts[TriIndex*3+1],
                             pTriPts[TriIndex*3+2],
                             Weights[0], Weights[1], Weights[2] );
}


void ImageToRay( EnvironmentInfo& EI, float ImgX, float ImgY, float& WXWZ, float& WYWZ )
{
    float Up = (ImgX-EI.CX)*EI.SX;
    float Vp = (ImgY-EI.CY)*EI.SY;
    float U = Up / (1+EI.Kappa*(Up*Up+Vp*Vp));
    float V = Vp / (1+EI.Kappa*(Up*Up+Vp*Vp));
    WXWZ = U / EI.Focus;
    WYWZ = V / EI.Focus;

    log("IMAGE_TO_RAY %f %f %f %f %f %f %f %f\n",ImgX,ImgY,Up,Vp,U,V,WXWZ,WYWZ);
}

void WorldToImage( EnvironmentInfo& EI, float WX, float WY, float WZ, float& ImgX, float& ImgY )
{
    float U = EI.Focus * WX / WZ;
    float V = EI.Focus * WY / WZ;
    float Up = (2*U) / (1+sqrtf(1-4*EI.Kappa*(U*U+V*V)));
    float Vp = (2*V) / (1+sqrtf(1-4*EI.Kappa*(U*U+V*V)));
    ImgX = Up/EI.SX + EI.CX;
    ImgY = Vp/EI.SY + EI.CY;
}

float ScoreQ3DFIRST(  EnvironmentInfo& EI, Quad Q, float* Q3D )
{
    float WorldLen[6] = {0.32,0.22,0.32,0.22,0.388,0.388};
    int WorldLenIA[6] = {0,1,2,3,0,1};
    int WorldLenIB[6] = {1,2,3,0,2,3};
    float ProjTotalErr=0;
    float EdgeLenTotalErr=0;

    for( int i=0; i<4; i++ )
    {
        float ImgX=0;
        float ImgY=0;
        WorldToImage(EI,Q3D[i*3+0],Q3D[i*3+1],Q3D[i*3+2],ImgX,ImgY);
        float DX = Q.Pts[i].X - ImgX;
        float DY = Q.Pts[i].Y - ImgY;
        ProjTotalErr += (DX*DX + DY*DY);
    }
    log("ProjTotalErr: %f\n",ProjTotalErr);
    for( int i=0; i<6; i++ )
    {
        int IA = WorldLenIA[i];
        int IB = WorldLenIB[i];
        float WL = WorldLen[i];
        float WDX = Q3D[IA*3+0] - Q3D[IB*3+0];
        float WDY = Q3D[IA*3+1] - Q3D[IB*3+1];
        float WDZ = Q3D[IA*3+2] - Q3D[IB*3+2];
        float WLen = sqrtf(WDX*WDX + WDY*WDY + WDZ*WDZ);
        EdgeLenTotalErr += (WLen-WL)*(WLen-WL);
    }
    log("EdgeLenTotalErr: %f\n",EdgeLenTotalErr);
    return ProjTotalErr*100 + EdgeLenTotalErr*1000;
}

    bool RobustCameraRaySphereIntersect( Vec3 Ray, Vec3 SphereCenter, double Radius, Vec3& P0, Vec3& P1 )
    {
        // return P0,P1 if two intersections
        // return P0=P1 if only one intersection
        // return P0=P1 of closest approach

        //Compute A, B and C coefficients
        // Ray origin is (0,0,0)
        Vec3 RayOrigin(-SphereCenter.X,-SphereCenter.Y,-SphereCenter.Z);
        double a = Ray.Dot(Ray);//dot(ray.d, ray.d);
        double b = 2 * Ray.Dot(RayOrigin);
        double c = RayOrigin.Dot(RayOrigin) - Radius * Radius;//dot(ray.o, ray.o) - (r * r);

        //Find discriminant
        double disc = b * b - 4 * a * c;

        //log("a:%f b:%f c:%f disc:%f\n",a,b,c,disc);
    
        // if discriminant is negative there are no real roots, so ray misses sphere
        if (disc < 0)
        {
            // solve closest pt on Ray to SphereCenter
            Vec3 NormRay = Ray.Normalize();
            double DistAlongRay = SphereCenter.Dot(NormRay);
            P0.X = NormRay.X * DistAlongRay;
            P0.Y = NormRay.Y * DistAlongRay;
            P0.Z = NormRay.Z * DistAlongRay;
            P1 = P0;
            return true;
        }

        // compute q as described above
        double distSqrt = sqrt(disc);
        double q;
        if (b < 0)  q = (-b - distSqrt)/2.0;
        else        q = (-b + distSqrt)/2.0;

        // compute t0 and t1
        double t0 = q / a;
        double t1 = c / q;

        // make sure t0 is smaller than t1
        if (t0 > t1) { double temp = t0; t0 = t1; t1 = temp; }

        // if t1 is less than zero, the object is in the ray's negative direction
        // and consequently the ray misses the sphere
        if (t1 < 0) return false;

        //log("t0:%f  t1:%f\n",t0,t1);
        // if t0 is less than zero, the intersection point is at t1
        if (t0 < 0) { t0 = t1; }

        P0.X = Ray.X * t0;
        P0.Y = Ray.Y * t0;
        P0.Z = Ray.Z * t0;
        P1.X = Ray.X * t1;
        P1.Y = Ray.Y * t1;
        P1.Z = Ray.Z * t1;

        //log("raysph %f %f %f %f, %f %f %f %f\n", t0, P0.X, P0.Y, P0.Z, t1, P1.X, P1.Y, P1.Z);
        return true;
    }


    float Score4PointPose( Vec3 P0, Vec3 P1, Vec3 P2, Vec3 P3, double* MDist, double WA, double WB, double WC )
    {
        Vec3 PlanarP3 = CalcBaryPosition(P0, P1, P2, WA, WB, WC);
        double DP = PlanarP3.Dist(P3);
        double D01 = MDist[0] - P0.Dist(P1);
        double D02 = MDist[1] - P0.Dist(P2);
        double D03 = MDist[2] - P0.Dist(P3);
        double D12 = MDist[3] - P1.Dist(P2);
        double D13 = MDist[4] - P1.Dist(P3);
        double D23 = MDist[5] - P2.Dist(P3);

        //log("%f,%f,%f,%f,%f,%f,%f\n",DP,D01,D02,D03,D12,D13,D23);

        double PairErr = DP*DP + D01*D01 + D02*D02 + D03*D03 + D12*D12 + D13*D13 + D23*D23;
        return sqrtf(PairErr/7);
    }

    float Solve4PointCoplanarPose( EnvironmentInfo& EI, vector<Vec3> RayDir, vector<Vec3> ModelPos, vector<Vec3>& CamPos )
    {
        Vec3* pRayDir = &RayDir[0];
        Vec3* pModel = &ModelPos[0];
        Vec3  NRayDir[4];
        double Zs[4];
        double ModelEdgeDist[4];
        double ZeroToPtDist[4];

        for( int i=0; i<4; i++ )
        {
            NRayDir[i] = RayDir[i].Normalize();
            ModelEdgeDist[i] = ModelPos[i].Dist(ModelPos[(i+1)%4]);
            ZeroToPtDist[i] = ModelPos[0].Dist(ModelPos[i]);

            log(" RayDir %d %f %f %f\n",i,RayDir[i].X,RayDir[i].Y,RayDir[i].Z);
            log("NRayDir %d %f %f %f\n",i,NRayDir[i].X,NRayDir[i].Y,NRayDir[i].Z);
        }

        double MDist[6];
        MDist[0] = (ModelPos[0].Dist(ModelPos[1]));
        MDist[1] = (ModelPos[0].Dist(ModelPos[2]));
        MDist[2] = (ModelPos[0].Dist(ModelPos[3]));
        MDist[3] = (ModelPos[1].Dist(ModelPos[2]));
        MDist[4] = (ModelPos[1].Dist(ModelPos[3]));
        MDist[5] = (ModelPos[2].Dist(ModelPos[3]));

        double ModelDist02 = ModelPos[0].Dist(ModelPos[2]);
        double ModelDist13 = ModelPos[1].Dist(ModelPos[3]);
        double BestErr = 1000000.0;

        float P3WA = 0;
        float P3WB = 0;
        float P3WC = 0;
        CalcBaryXY(ModelPos[3], ModelPos[0], ModelPos[1], ModelPos[2], P3WA, P3WB, P3WC);


        // Loop through range of Camera Zs for point 0
        double StartZ = 0.100;
        double EndZ = 1.500;
        double StepZ = 0.001;
        for( double Z0 = StartZ; Z0<=EndZ; Z0+=StepZ )
        {
            Vec3 P0;
            P0.X = pRayDir[0].X*Z0;
            P0.Y = pRayDir[0].Y*Z0;
            P0.Z = pRayDir[0].Z*Z0;

            // P1 should lie on a sphere of radius ZeroToPtDist[1] from P0
            Vec3 P1s[2];
            if( !RobustCameraRaySphereIntersect( NRayDir[1], P0, ZeroToPtDist[1], P1s[0], P1s[1] ) )
                continue;

            // P2 should lie on a sphere of radius ZeroToPtDist[2] from P0
            Vec3 P2s[2];
            if( !RobustCameraRaySphereIntersect( NRayDir[2], P0, ZeroToPtDist[2], P2s[0], P2s[1] ) )
                continue;

            // P3 should lie on a sphere of radius ZeroToPtDist[3] from P0
            Vec3 P3s[2];
            if( !RobustCameraRaySphereIntersect( NRayDir[3], P0, ZeroToPtDist[3], P3s[0], P3s[1] ) )
                continue;

            // There should be a P1 to P2 dist = ModelEdgeDist[1]
            Vec3 Pts[4];
            Pts[0] = P0;
            double BestP1P2P3Error = 10000000;
            {
                for( int i=0; i<2; i++ )
                {
                    Vec3 P1 = P1s[i];
                    for( int j=0; j<2; j++ )
                    {
                        Vec3 P2 = P2s[j];
                        for (int k = 0; k < 2; k++)
                        {
                            Vec3 P3 = P3s[k];

                            double PairErr = Score4PointPose(P0,P1,P2,P3,MDist,P3WA, P3WB, P3WC);

                            if (PairErr < BestP1P2P3Error)
                            {
                                BestP1P2P3Error = PairErr;
                                Pts[1] = P1;
                                Pts[2] = P2;
                                Pts[3] = P3;
                            }
                        }
                    }
                }
            }

            // Add up total error for this solution
            double Err = Score4PointPose(Pts[0],Pts[1],Pts[2],Pts[3],MDist,P3WA, P3WB, P3WC);

            //log("POSE %8.3f %1.6f\n",Z0,Err);

            if( Err < BestErr )
            {
                BestErr = Err;
                CamPos[0] = Pts[0];
                CamPos[1] = Pts[1];
                CamPos[2] = Pts[2];
                CamPos[3] = Pts[3];
            }
        }

        return BestErr;
    }


float ScoreQ3D(  EnvironmentInfo& EI, Quad Q, float* Q3D )
{
    float WorldLen[6] = {0.32f,0.22f,0.32f,0.22f,0.388f,0.388f};
    int WorldLenIA[6] = {0,1,2,3,0,1};
    int WorldLenIB[6] = {1,2,3,0,2,3};
    float ProjTotalErr=0;
    float EdgeLenTotalErr=0;

    for( int i=0; i<4; i++ )
    {
        float ImgX=0;
        float ImgY=0;
        WorldToImage(EI,Q3D[i*3+0],Q3D[i*3+1],Q3D[i*3+2],ImgX,ImgY);
        float DX = Q.Pts[i].X - ImgX;
        float DY = Q.Pts[i].Y - ImgY;
        ProjTotalErr += (DX*DX + DY*DY);
    }

    for( int i=0; i<6; i++ )
    {
        int IA = WorldLenIA[i];
        int IB = WorldLenIB[i];
        float WL = WorldLen[i];
        float WDX = Q3D[IA*3+0] - Q3D[IB*3+0];
        float WDY = Q3D[IA*3+1] - Q3D[IB*3+1];
        float WDZ = Q3D[IA*3+2] - Q3D[IB*3+2];
        float WLen = sqrtf(WDX*WDX + WDY*WDY + WDZ*WDZ);
        EdgeLenTotalErr += (WLen-WL)*(WLen-WL);
    }
    return ProjTotalErr*100 + EdgeLenTotalErr*1000;
}

void Solve3DQuad( EnvironmentInfo& EI, Quad Q, float* Q3D )
{

    float Dirs[8];
    float WorldLen[4] = {0.32,0.22,0.32,0.22};
    float ImgLen[4] = {0};
    //float AvgZ[4] = {0};

    log("Solve3DQuad -------------------\n");

    log("EI.Dir:  %d\n",EI.Directory);
    log("EI.Side: %d\n",EI.Side);
    log("EI.Focus:%f EI.Kappa:%f\n",EI.Focus, EI.Kappa);
    log("EI.CX:%f EI.CY:%f\n",EI.CX,EI.CY);
    log("EI.SX:%10.8f EI.SI:%10.8f\n",EI.SX,EI.SY);
    
    //// Compute camera rays
    //for( int i=0; i<4; i++ )
    //{
    //    ImageToRay( EI,Q.Pts[i].X, Q.Pts[i].Y, Dirs[i*2+0], Dirs[i*2+1] );
    //    ImgLen[i] = Q.Pts[i].Dist(Q.Pts[(i+1)%4]);
    //    AvgZ[i] = (EI.Focus/EI.SX) * WorldLen[i] / ImgLen[i];

    //    log("[%d] Dir:(%f,%f) AvgZ:%f\n",i,Dirs[i*2+0],Dirs[i*2+1],AvgZ[i]);
    //}

    vector<Vec3> RayDir;
    vector<Vec3> ModelPos;
    vector<Vec3> CamPos;
    Vec2 PanelRefCorners[4] = {Vec2(0,0), Vec2(0.0f,PANELHEIGHT), Vec2(PANELWIDTH,PANELHEIGHT), Vec2(PANELWIDTH,0.0f)};
    for( int i=0; i<4; i++ )
    {
        Vec3 RDir(0,0,1);
        ImageToRay(EI,Q.Pts[i].X,Q.Pts[i].Y,RDir.X,RDir.Y);
        RayDir.push_back(RDir);
        ModelPos.push_back( PanelRefCorners[i]);
        CamPos.push_back(Vec3(0,0,0));
    }

    Solve4PointCoplanarPose(EI,RayDir,ModelPos,CamPos);

    for( int i=0; i<4; i++ )
    {
        Q3D[i*3+0] = CamPos[i].X;
        Q3D[i*3+1] = CamPos[i].Y;
        Q3D[i*3+2] = CamPos[i].Z;
    }

/*
    // Solve Zs that minimize projection error
    for( int i=0; i<4; i++ )
    {
        float Z = (AvgZ[i] + AvgZ[(i+3)%4])*0.5f;
        Q3D[i*3+0] = Dirs[i*2+0] * Z;
        Q3D[i*3+1] = Dirs[i*2+1] * Z;
        Q3D[i*3+2] = Z;
    }

    float BestDirs[8];
    memcpy(BestDirs,Dirs,sizeof(BestDirs));

    float BestZs[4] = {Q3D[0*3+2],Q3D[1*3+2],Q3D[2*3+2],Q3D[3*3+2]};

    for( int i=0; i<4; i++ ) log("Q   [%d]  (%f,%f)\n",i, Q.Pts[i].X, Q.Pts[i].Y);
    for( int i=0; i<4; i++ ) log("Q3D [%d]  (%f,%f,%f)\n",i, Q3D[i*3+0],Q3D[i*3+1],Q3D[i*3+2]);

    float BestErr = ScoreQ3DFIRST(EI, Q, Q3D );


    for( int nTrials=0; nTrials<5000; nTrials++ )
    {
        if( (nTrials%500)==0 )
            log("[%d] BestErr %f\n",nTrials,BestErr);

        float NewDirs[8];
        float NewZs[4];
        memcpy(NewDirs,BestDirs,sizeof(NewDirs));
        memcpy(NewZs,BestZs,sizeof(NewZs));

        //if( (nTrials<5000) || (myrand()%1000 > 500) )
        {
            // Bump Z
            int I = myrand()%4;
            //float Amt = -0.010f + (float)((myrand()%1000)/1000.0f)*(0.020f);
            float Amt = -0.020f + (float)((myrand()%1000)/1000.0f)*(0.040f);
            NewZs[I] += Amt;

        }
        //else
        //{
        //    // Bump ImgXY
        //    int I = myrand()%4;
        //    float AmtX = -0.10f + (float)((myrand()%1000)/1000.0f)*(0.20f);
        //    float AmtY = -0.10f + (float)((myrand()%1000)/1000.0f)*(0.20f);
        //    Q.Pts[I].X += AmtX;
        //    Q.Pts[I].Y += AmtY;
        //    ImageToRay( Q.Pts[I].X, Q.Pts[I].Y, NewDirs[I*2+0], NewDirs[I*2+1] );
        //}

        for( int i=0; i<4; i++ )
        {
            Q3D[i*3+0] = NewDirs[i*2+0] * NewZs[i];
            Q3D[i*3+1] = NewDirs[i*2+1] * NewZs[i];
            Q3D[i*3+2] = NewZs[i];
        }

        float Err = ScoreQ3D(EI, Q, Q3D );
        if( Err < BestErr )
        {
            BestErr = Err;
            memcpy(BestDirs,NewDirs,sizeof(NewDirs));
            memcpy(BestZs,NewZs,sizeof(NewZs));
            //log("[%6d] BestErr: %f\n",nTrials,BestErr);
        }
    }

    for( int i=0; i<4; i++ )
    {
        Q3D[i*3+0] = BestDirs[i*2+0] * BestZs[i];
        Q3D[i*3+1] = BestDirs[i*2+1] * BestZs[i];
        Q3D[i*3+2] = BestZs[i];
    }
    */

    for( int i=0; i<4; i++ ) log("Final Q3D [%d]  (%f,%f,%f)\n",i, Q3D[i*3+0],Q3D[i*3+1],Q3D[i*3+2]);
}

Vec2 CalcPerspectiveAwareBilinearPos( EnvironmentInfo& EI, float* Q3D, float XT, float YT, float ZOffset )
{
    float W00 = (1-XT)*(1-YT);
    float W01 = (1-XT)*(YT);
    float W11 = (XT*YT);
    float W10 = (XT)*(1-YT);

    float WX = Q3D[0*3+0]*W00 + Q3D[1*3+0]*W01 + Q3D[2*3+0]*W11 + Q3D[3*3+0]*W10;
    float WY = Q3D[0*3+1]*W00 + Q3D[1*3+1]*W01 + Q3D[2*3+1]*W11 + Q3D[3*3+1]*W10;
    float WZ = Q3D[0*3+2]*W00 + Q3D[1*3+2]*W01 + Q3D[2*3+2]*W11 + Q3D[3*3+2]*W10;


    // Compute normal
    Vec3 V01 = Vec3( Q3D[1*3+0]-Q3D[0*3+0],  Q3D[1*3+1]-Q3D[0*3+1],  Q3D[1*3+2]-Q3D[0*3+2] );
    Vec3 V03 = Vec3( Q3D[3*3+0]-Q3D[0*3+0],  Q3D[3*3+1]-Q3D[0*3+1],  Q3D[3*3+2]-Q3D[0*3+2] );
    Vec3 Dir = V01.Cross(V03);
    Vec3 Normal = Dir.Normalize();

    WX += Normal.X * ZOffset;
    WY += Normal.Y * ZOffset;
    WZ += Normal.Z * ZOffset;

    float ImgX;
    float ImgY;
    WorldToImage(EI,WX,WY,WZ,ImgX,ImgY);

    return Vec2(ImgX,ImgY);
}

//=========================================================================
//=========================================================================
//=========================================================================
//=========================================================================
enum ENUM_ITEMS
{
    ITEM_PANEL_POWER_SWITCH=0,        
    ITEM_PANEL_POWER_COVER,        
    ITEM_PANEL_POWER_LED,        
    ITEM_A01_ROCKER_SWITCH,        
    ITEM_A01_ROCKER_LED_TOP,        
    ITEM_A01_ROCKER_LED_BOTTOM,  
    ITEM_A02_LED_NUM_PAD_A1,        
    ITEM_A02_LED_NUM_PAD_A2,        
    ITEM_A02_LED_NUM_PAD_A3,        
    ITEM_A02_LED_NUM_PAD_B1,        
    ITEM_A02_LED_NUM_PAD_B2,        
    ITEM_A02_LED_NUM_PAD_B3,        
    ITEM_A02_LED_NUM_PAD_C1,        
    ITEM_A02_LED_NUM_PAD_C2,        
    ITEM_A02_LED_NUM_PAD_C3,        
    ITEM_A03_TOGGLE,                
    ITEM_A03_LED,                
    ITEM_A04_TOGGLE,                
    ITEM_A04_LED_TOP,            
    ITEM_A04_LED_BOTTOM,            
    ITEM_A05_TOGGLE,                
    ITEM_A05_LED,             

    // Additional reference points
    ITEM_PANEL_POWER_COVER_UP,
    ITEM_PANEL_POWER_COVER_SIM,
    ITEM_PANEL_POWER_COVER_UP_SIM,
};

char* StateNames[] = 
{
    "UP",
    "DOWN",
    "ON",
    "OFF",
    "CENTER",
    "OUT",
    "HIDDEN",
};

enum ENUM_STATES
{
    STATE_UP=0,
    STATE_DOWN,
    STATE_ON,
    STATE_OFF,
    STATE_CENTER,
    STATE_OUT,
    STATE_HIDDEN,
};

struct ItemInfo
{
    //----- start of default structures
    int Index;
    char* Name;
    Vec2 RefPos;
    float ZOffset;
    int StateMask[7];
    int DefaultState;
    //----- end of default structures

    Vec3 PanelPos;

    float StateConfidence[2][7];
    Vec2 Pos[2];

    int FinalState;
};

ItemInfo Items[] = 
{
    {ITEM_PANEL_POWER_SWITCH,    "PANEL_POWER_SWITCH",        Vec2(352,112),0.004,        1,1,0,0,0,0,1,  STATE_DOWN},
    {ITEM_PANEL_POWER_COVER,     "PANEL_POWER_COVER",        Vec2(350,135),0.027,        1,1,0,0,0,0,1,  STATE_DOWN},
    {ITEM_PANEL_POWER_LED,       "PANEL_POWER_LED",            Vec2(396,112),0.002,        0,0,1,1,0,0,1,  STATE_OFF},
    {ITEM_A01_ROCKER_SWITCH,     "A01_ROCKER_SWITCH",        Vec2(396,352),0.012,    1,1,0,0,1,0,1,  STATE_CENTER},
    {ITEM_A01_ROCKER_LED_TOP,    "A01_ROCKER_LED_TOP",        Vec2(397,291),0.002,        0,0,1,1,0,0,1,  STATE_OFF},
    {ITEM_A01_ROCKER_LED_BOTTOM,"A01_ROCKER_LED_BOTTOM",    Vec2(397,411),0.002,        0,0,1,1,0,0,1,  STATE_OFF},
    {ITEM_A02_LED_NUM_PAD_A1,    "A02_LED_NUM_PAD_A1",        Vec2( 99,288),0.004,    0,0,1,1,0,0,1,  STATE_OFF},
    {ITEM_A02_LED_NUM_PAD_A2,    "A02_LED_NUM_PAD_A2",        Vec2(162,288),0.004,    0,0,1,1,0,0,1,  STATE_OFF},
    {ITEM_A02_LED_NUM_PAD_A3,    "A02_LED_NUM_PAD_A3",        Vec2(226,288),0.004,    0,0,1,1,0,0,1,  STATE_OFF},
    {ITEM_A02_LED_NUM_PAD_B1,    "A02_LED_NUM_PAD_B1",        Vec2( 98,351),0.004,    0,0,1,1,0,0,1,  STATE_OFF},
    {ITEM_A02_LED_NUM_PAD_B2,    "A02_LED_NUM_PAD_B2",        Vec2(163,351),0.004,    0,0,1,1,0,0,1,  STATE_OFF},
    {ITEM_A02_LED_NUM_PAD_B3,    "A02_LED_NUM_PAD_B3",        Vec2(226,351),0.004,    0,0,1,1,0,0,1,  STATE_OFF},
    {ITEM_A02_LED_NUM_PAD_C1,    "A02_LED_NUM_PAD_C1",        Vec2( 98,415),0.004,    0,0,1,1,0,0,1,  STATE_OFF},
    {ITEM_A02_LED_NUM_PAD_C2,    "A02_LED_NUM_PAD_C2",        Vec2(162,415),0.004,    0,0,1,1,0,0,1,  STATE_OFF},
    {ITEM_A02_LED_NUM_PAD_C3,    "A02_LED_NUM_PAD_C3",        Vec2(226,414),0.004,    0,0,1,1,0,0,1,  STATE_OFF},
    {ITEM_A03_TOGGLE,            "A03_TOGGLE",                Vec2(109,585),0.004,    1,1,0,0,0,1,1,  STATE_UP},
    {ITEM_A03_LED,               "A03_LED",                    Vec2(110,524),0.002,        0,0,1,1,0,0,1,  STATE_OFF},
    {ITEM_A04_TOGGLE,            "A04_TOGGLE",                Vec2(262,584),0.004,    1,1,0,0,1,1,1,  STATE_UP},
    {ITEM_A04_LED_TOP,           "A04_LED_TOP",                Vec2(262,524),0.002,        0,0,1,1,0,0,1,  STATE_OFF},
    {ITEM_A04_LED_BOTTOM,        "A04_LED_BOTTOM",            Vec2(262,644),0.002,        0,0,1,1,0,0,1,  STATE_OFF},
    {ITEM_A05_TOGGLE,            "A05_TOGGLE",                Vec2(412,584),0.004,    1,1,0,0,0,1,1,  STATE_UP},
    {ITEM_A05_LED,               "A05_LED",                    Vec2(412,524),0.002,        0,0,1,1,0,0,1,  STATE_OFF},


    {ITEM_PANEL_POWER_COVER_UP,     "PANEL_POWER_COVER_UP",     Vec2(348,23),0.040,        1,1,0,0,0,0,1,  STATE_DOWN},
    {ITEM_PANEL_POWER_COVER_SIM,    "PANEL_POWER_COVER_SIM",    Vec2(351,170),0.025,        1,1,0,0,0,0,1,  STATE_DOWN},
    {ITEM_PANEL_POWER_COVER_UP_SIM, "PANEL_POWER_COVER_UP_SIM",    Vec2(345,54),0.040,        1,1,0,0,0,0,1,  STATE_DOWN},
};
            
void InitItems( void )
{
    log("n item references %d\n",sizeof(Items)/sizeof(Items[0]));
    // Put any other default information into the Items
    for( int i=0; i<sizeof(Items)/sizeof(Items[0]); i++ )
    {
        Items[i].RefPos.X /= 518.0f;
        Items[i].RefPos.Y /= 750.0f;

        Items[i].PanelPos.X = Items[i].RefPos.X * 0.22f;
        Items[i].PanelPos.Y = Items[i].RefPos.Y * 0.32f;
        Items[i].PanelPos.Z = Items[i].ZOffset;

        memset(Items[i].StateConfidence,0,sizeof(Items[i].StateConfidence));
    }
}

//=========================================================================


//=========================================================================

string ImageName;

void SaveBMP( const char* prefix,const char* filename, float* Pix, int W, int H)
{
#ifndef ALLOW_SAVEBMP
    (void)filename;
    (void)Pix;
    (void)W;
    (void)H;
    return;
#else
 
    int DstStride = W*3;
    while( (DstStride%4) !=0 ) DstStride++;
 
    int HeaderSize = 54;
    int FileSize = (DstStride*H) + HeaderSize;
    byte* Buffer = (byte*)malloc(FileSize);
    memset(Buffer,0,FileSize);
    byte* Pixels = Buffer + HeaderSize;
    if( !Buffer )
        return;
 
    float* pSrcBaseRow = Pix;
 
    // Copy memory into pixel buffer
    for( int y=0; y<H; y++ )
    {
        float* pSrc = Pix + y*W;
        byte* pDst = Pixels + ((H-1-y)*DstStride);
 
        for( int x=0; x<W; x++ )
        {
            float V = pSrc[0];
            byte B = (byte)(V*255);
            pDst[0] = B;
            pDst[1] = B;
            pDst[2] = B;
            pSrc++;
            pDst+=3;
        }
    }

    // Now write windows bitmap header
    int I=0;
    *((byte*) (Buffer+I)) = 'B';                I+=1;
    *((byte*) (Buffer+I)) = 'M';                I+=1;
    *((int*)  (Buffer+I)) = (int)FileSize;        I+=4;
    *((int*)  (Buffer+I)) = (int)0;                I+=4;
    *((int*)  (Buffer+I)) = (int)HeaderSize;    I+=4;
    *((int*)  (Buffer+I)) = (int)40;            I+=4;
    *((int*)  (Buffer+I)) = (int)W;                I+=4;
    *((int*)  (Buffer+I)) = (int)H;                I+=4;
    *((short*)(Buffer+I)) = (int)1;                I+=2;
    *((short*)(Buffer+I)) = (int)24;            I+=2;
    *((int*)  (Buffer+I)) = (int)0;                I+=4;
    *((int*)  (Buffer+I)) = (int)FileSize;        I+=4;
    *((int*)  (Buffer+I)) = (int)2835;            I+=4;
    *((int*)  (Buffer+I)) = (int)2835;            I+=4;
    *((int*)  (Buffer+I)) = (int)0;                I+=4;
    *((int*)  (Buffer+I)) = (int)0;                I+=4;
 
    // Write to disk
    char fname[256];
    sprintf(fname,"%s_%s",prefix,filename);
    FILE* fp = fopen(fname,"wb");
    if( fp )
    {
        fwrite(Buffer,FileSize,1,fp);
        fclose(fp);
    }
 
    free(Buffer);

#endif
}

#ifdef SUBMISSION
inline void DrawLine( float* Pix, int W, int H, float V, float X0, float Y0, float X1, float Y1 ) {}
inline void DrawCircle( float* Pix, int W, int H, float V, float CX, float CY, float R ) {}
inline void DrawRect( float* Pix, int W, int H, float V, float CX, float CY, float RW, float RH ) {}

#else
void DrawLine( float* Pix, int W, int H, float V, float X0, float Y0, float X1, float Y1 )
{
    float DX = X1-X0;
    float DY = Y1-Y0;
    int nSteps = (int)(2 + (2*sqrtf(DX*DX + DY*DY)));
    for( int i=0; i<nSteps; i++ )
    {
        float T = (float)i / (float)(nSteps-1);
        int X = (int)(X0 + T*(DX) + 0.5f);
        int Y = (int)(Y0 + T*(DY) + 0.5f);
        if( (X>=0) && (X<W) && (Y>=0) && (Y<H) )
            Pix[X+Y*W] = V;
    }
}

void DrawCircle( float* Pix, int W, int H, float V, float CX, float CY, float R )
{
    int nSteps = (int)((R*2*3.141)*2);
    for( int i=0; i<nSteps; i++ )
    {
        float T = (float)i / (float)(nSteps-1);
        float DX = R*cosf( T*3.141592f*2);
        float DY = R*sinf( T*3.141592f*2);
        int X = (int)(CX + DX + 0.5f);
        int Y = (int)(CY + DY + 0.5f);
        if( (X>=0) && (X<W) && (Y>=0) && (Y<H) )
            Pix[X+Y*W] = V;
    }
}

void DrawRect( float* Pix, int W, int H, float V, float CX, float CY, float RW, float RH )
{
    DrawLine(Pix,W,H,V,CX,CY,CX,CY+RH);
    DrawLine(Pix,W,H,V,CX,CY+RH,CX+RW,CY+RH);
    DrawLine(Pix,W,H,V,CX+RW,CY+RH,CX+RW,CY);
    DrawLine(Pix,W,H,V,CX+RW,CY,CX,CY);
}

#endif

void ComputeEdges( float* PixDst, float* PixSrc, int W, int H )
{
    memset(PixDst,0,W*H*sizeof(float));

    float MinG = +10000000;
    float MaxG = -10000000;
    for( int Y=1; Y<H-1; Y++ )
        for( int X=1; X<W-1; X++ )
        {
            int O = X+Y*W;

            float GX =      PixSrc[O-1-W] -   PixSrc[O+1-W] +
                        2*PixSrc[O-1  ] - 2*PixSrc[O+1  ] +
                          PixSrc[O-1+W] -   PixSrc[O+1+W];

            float GY =      PixSrc[O-1-W] + 2*PixSrc[O-W] + PixSrc[O+1-W]
                          -PixSrc[O-1+W] - 2*PixSrc[O+W] - PixSrc[O+1+W];

            float G = sqrtf(GX*GX + GY*GY);
            if( G<MinG ) MinG = G;
            if( G>MaxG ) MaxG = G;

            PixDst[O] = G;
        }

    for( int i=0; i<W*H; i++ )
    {
        float T = (PixDst[i]-MinG)/(MaxG-MinG);
        if( T<0 ) T=0;
        if( T>1 ) T=1;
        PixDst[i] = T;
    }
}

void ComputeEdges2( float* PixDst, float* PixSrc, int W, int H )
{
    memset(PixDst,0,W*H*sizeof(float));

    float MinEV = +10000000;
    float MaxEV = -10000000;

    for( int Y=12; Y<H-12; Y++ )
        for( int X=12; X<W-12; X++ )
        {
            int O = X+Y*W;

            float LV = (PixSrc[O-4] + PixSrc[O-8] + PixSrc[O-12])/3;
            float RV = (PixSrc[O+4] + PixSrc[O+8] + PixSrc[O+12])/3;
            float Diff = RV-LV;
            

            float EV = 0;
            if( Diff > 0 )
                EV = fabs(Diff);
            if( EV<MinEV ) MinEV = EV;
            if( EV>MaxEV ) MaxEV = EV;

            PixDst[O] = EV;
        }

    for( int i=0; i<W*H; i++ )
    {
        float T = (PixDst[i]-MinEV)/(MaxEV-MinEV);
        if( T<0 ) T=0;
        if( T>1 ) T=1;
        PixDst[i] = T;
    }
}

void ComputeEdges3( float* PixDst, float* PixSrc, int W, int H )
{
    memset(PixDst,0,W*H*sizeof(float));

    for( int Y=16; Y<H-16; Y++ )
        for( int X=16; X<W-16; X++ )
        {
            int O = X+Y*W;

            int LWC = 0;
            int LBC = 0;
            int RWC = 0;
            int RBC = 0;

            if( PixSrc[O-16]==0 ) LBC++; else LWC++;
            if( PixSrc[O-12]==0 ) LBC++; else LWC++;
            if( PixSrc[O- 8]==0 ) LBC++; else LWC++;
            if( PixSrc[O- 4]==0 ) LBC++; else LWC++;
            if( PixSrc[O- 2]==0 ) LBC++; else LWC++;
            if( PixSrc[O+16]==0 ) RBC++; else RWC++;
            if( PixSrc[O+12]==0 ) RBC++; else RWC++;
            if( PixSrc[O+ 8]==0 ) RBC++; else RWC++;
            if( PixSrc[O+ 4]==0 ) RBC++; else RWC++;
            if( PixSrc[O+ 2]==0 ) RBC++; else RWC++;

            //if( PixSrc[O-16-W*4]==0 ) LBC++; else LWC++;
            //if( PixSrc[O-12-W*4]==0 ) LBC++; else LWC++;
            //if( PixSrc[O- 8-W*4]==0 ) LBC++; else LWC++;
            //if( PixSrc[O- 4-W*4]==0 ) LBC++; else LWC++;
            //if( PixSrc[O- 2-W*4]==0 ) LBC++; else LWC++;
            //if( PixSrc[O+16-W*4]==0 ) RBC++; else RWC++;
            //if( PixSrc[O+12-W*4]==0 ) RBC++; else RWC++;
            //if( PixSrc[O+ 8-W*4]==0 ) RBC++; else RWC++;
            //if( PixSrc[O+ 4-W*4]==0 ) RBC++; else RWC++;
            //if( PixSrc[O+ 2-W*4]==0 ) RBC++; else RWC++;

            //if( PixSrc[O-16+W*4]==0 ) LBC++; else LWC++;
            //if( PixSrc[O-12+W*4]==0 ) LBC++; else LWC++;
            //if( PixSrc[O- 8+W*4]==0 ) LBC++; else LWC++;
            //if( PixSrc[O- 4+W*4]==0 ) LBC++; else LWC++;
            //if( PixSrc[O- 2+W*4]==0 ) LBC++; else LWC++;
            //if( PixSrc[O+16+W*4]==0 ) RBC++; else RWC++;
            //if( PixSrc[O+12+W*4]==0 ) RBC++; else RWC++;
            //if( PixSrc[O+ 8+W*4]==0 ) RBC++; else RWC++;
            //if( PixSrc[O+ 4+W*4]==0 ) RBC++; else RWC++;
            //if( PixSrc[O+ 2+W*4]==0 ) RBC++; else RWC++;

            int TWC = 0;
            int TBC = 0;
            int BWC = 0;
            int BBC = 0;

            //if( PixSrc[O-16*W]==0 ) TBC++; else TWC++;
            //if( PixSrc[O-12*W]==0 ) TBC++; else TWC++;
            //if( PixSrc[O- 8*W]==0 ) TBC++; else TWC++;
            //if( PixSrc[O- 4*W]==0 ) TBC++; else TWC++;
            //if( PixSrc[O- 2*W]==0 ) TBC++; else TWC++;
            //if( PixSrc[O+16*W]==0 ) BBC++; else BWC++;
            //if( PixSrc[O+12*W]==0 ) BBC++; else BWC++;
            //if( PixSrc[O+ 8*W]==0 ) BBC++; else BWC++;
            //if( PixSrc[O+ 4*W]==0 ) BBC++; else BWC++;
            //if( PixSrc[O+ 2*W]==0 ) BBC++; else BWC++;


            if((LBC==5) && (RWC==5)) PixDst[O] = 0.25f;
            if((LWC==5) && (RBC==5)) PixDst[O] = 0.50f;
            if((TBC==5) && (BWC==5)) PixDst[O] = 0.75f;
            if((TWC==5) && (BBC==5)) PixDst[O] = 1.00f;
        }
}


float* ComputeHistogram( int nBins, float* Pix, int W, int H )
{
    float* Bins = (float*)malloc(sizeof(float)*nBins);
    memset(Bins,0,sizeof(float)*nBins);

    for( int i=0; i<W*H; i++ )
    {
        float V = Pix[i];
        int B = (V*(nBins-1)+0.5f);
        if( B<0 ) B=0;
        if( B>=nBins ) B=nBins-1;
        Bins[B]++;
    }

    return Bins;
}

float FindOtsuThreshold( float* Hist, int nBins )
{
    //for( int i=0; i<nBins; i++ ) log("%3d %f\n",i,Hist[i]);

    int BestK = 0;
    double BestScore = -10000000;
    for (int k = 1; k < nBins; k++)
    {
        double p1 = 0;
        double p2 = 0;
        double wsum1 = 0;
        double wsum2 = 0;
        for (int i = 0; i < nBins; i++)
        {
            if (i <= k)
            {
                p1 += Hist[i];
                wsum1 += i * Hist[i];
            }
            else
            {
                p2 += Hist[i];
                wsum2 += i * Hist[i];
            }
        }
        double p1p2 = p1 * p2;
        if (p1p2 == 0) p1p2 = 1;

        //double diff = wsum1 * p2 - wsum2 * p1;
        double diff = wsum1 * p2*2 - wsum2 * p1;    // Expecting 2X white pixels over black
        double score = (diff * diff) / (p1p2);

        if (score > BestScore)
        {
            BestScore = score;
            BestK = k;
        }
    }
    return (float)BestK/(nBins-1);
}


float* AllocPix( int W, int H )
{
    float* Dst = (float*)malloc(sizeof(float)*W*H);
    memset(Dst,0,sizeof(float)*W*H);
    return Dst;
}
void CopyPix( float* Dst, float* Src, int W, int H )
{
    memcpy(Dst,Src,W*H*4);
}
float* ClonePix( float* Pix, int W, int H )
{
    float* Dst = (float*)malloc(sizeof(float)*W*H);
    memcpy(Dst,Pix,sizeof(float)*W*H);
    return Dst;
}

void ClearOutsideRange( float Value, float* Pix, int W, int H, float Min, float Max )
{
    for( int i=0; i<W*H; i++ ) if( (Pix[i]<Min) || (Pix[i]>Max) ) Pix[i] = Value;
}
void ClearInsideRange( float Value, float* Pix, int W, int H, float Min, float Max )
{
    for( int i=0; i<W*H; i++ ) if( (Pix[i]>=Min) && (Pix[i]<=Max) ) Pix[i] = Value;
}

void ScalePixels( float S, float* Pix, int W, int H)
{
    for( int i=0; i<W*H; i++ ) Pix[i] *= S;
}

void MedianFilter( float* SrcPix, float* DstPix, int W, int H )
{
    float V[9];

    for( int CY=2; CY<H-2; CY++ )
    for( int CX=2; CX<W-2; CX++ )
    {
        int O = CX+CY*W;
        V[0] = SrcPix[O-1-W];
        V[1] = SrcPix[O-0-W];
        V[2] = SrcPix[O+1-W];
        V[3] = SrcPix[O-1-0];
        V[4] = SrcPix[O-0-0];
        V[5] = SrcPix[O+1-0];
        V[6] = SrcPix[O-1+W];
        V[7] = SrcPix[O-0+W];
        V[8] = SrcPix[O+1+W];

        for( int i=0; i<9; i++ )
        {
            int B=i;
            for( int j=i+1; j<9; j++ ) if( V[j] < V[B] ) B=j;
            float TV = V[i];
            V[i] = V[B];
            V[B] = TV;
        }

        DstPix[O] = V[4];
    }

}

void AvgFilter( float* SrcPix, float* DstPix, int W, int H )
{
    for( int CY=1; CY<H-1; CY++ )
    for( int CX=1; CX<W-1; CX++ )
    {
        int O = CX+CY*W;
        float Sum = SrcPix[O-1-W] + SrcPix[O-0-W] + SrcPix[O+1-W] +
                    SrcPix[O-1-0] + SrcPix[O-0-0] + SrcPix[O+1-0] +
                    SrcPix[O-1+W] + SrcPix[O-0+W] + SrcPix[O+1+W];
        DstPix[O] = Sum * (1.0f/9.0f);
    }
}

void ShrinkPix( float* Pix, int W, int H, float Threshold )
{
    float* Dst = ClonePix(Pix,W,H);

    for( int Y=1; Y<H-1; Y++ )
        for( int X=1; X<W-1; X++ )
        {
            int O = X+Y*W;
            if( Pix[O] >= Threshold )
            {
                int C=0;
                if( Pix[O-1] > Threshold ) C++;
                if( Pix[O+1] > Threshold ) C++;
                if( Pix[O-W] > Threshold ) C++;
                if( Pix[O+W] > Threshold ) C++;
                
                if( Pix[O-1-W] > Threshold ) C++;
                if( Pix[O+1-W] > Threshold ) C++;
                if( Pix[O-1+W] > Threshold ) C++;
                if( Pix[O+1+W] > Threshold ) C++;

                if( C<6 ) Dst[O] = 0;
                else Dst[O] = 1;
            }
        }
    CopyPix(Pix,Dst,W,H);
    free(Dst);
}

void SinglePixelNoiseRemoval( float* Pix, int W, int H )
{
    float* Dst = ClonePix(Pix,W,H);

    for( int Y=1; Y<H-1; Y++ )
        for( int X=1; X<W-1; X++ )
        {
            int O = X+Y*W;
            if( Pix[O] != 0 )
            {
                if( ((Pix[O-1]==0) && (Pix[O+1]==0)) ||
                    ((Pix[O-W]==0) && (Pix[O+W]==0)) )
                    Dst[O] = 0;
            }
        }
    CopyPix(Pix,Dst,W,H);
    free(Dst);
}

void GrowPix( float* Pix, int W, int H, float Threshold )
{
    float* Dst = ClonePix(Pix,W,H);

    for( int Y=1; Y<H-1; Y++ )
        for( int X=1; X<W-1; X++ )
        {
            int O = X+Y*W;
            if( Pix[O] < Threshold )
            {
                int C=0;
                if( Pix[O-1] > Threshold ) C++;
                if( Pix[O+1] > Threshold ) C++;
                if( Pix[O-W] > Threshold ) C++;
                if( Pix[O+W] > Threshold ) C++;
                
                if( Pix[O-1-W] > Threshold ) C++;
                if( Pix[O+1-W] > Threshold ) C++;
                if( Pix[O-1+W] > Threshold ) C++;
                if( Pix[O+1+W] > Threshold ) C++;

                if( C>3 ) Dst[O] = 1;
                else Dst[O] = 0;
            }
        }
    CopyPix(Pix,Dst,W,H);
    free(Dst);
}


bool DoRecursiveFloodFill(int StartX, int StartY, int ID, float* Pix, int* IDs, int W, int H, int& MinX, int& MaxX, int& MinY, int& MaxY)
{
    int YOffset = StartY*W;
    int StartI = StartX + YOffset;
    if (IDs[StartI] != -1) return false;
    if (Pix[StartI] == 0) return false;

    MinY = min(MinY,StartY);
    MaxY = max(MaxY,StartY);

    IDs[StartI] = ID;

    // Find left and right extents 
    int LeftX = StartX;
    while ((LeftX > 0) && (IDs[LeftX - 1 + YOffset]==-1) && (Pix[LeftX - 1 + YOffset] > 0))
    {
        LeftX--;
        IDs[LeftX + YOffset] = ID;
    }

    int RightX = StartX;
    while ((RightX < W-1) && (IDs[RightX + 1 + YOffset]==-1) && (Pix[RightX + 1 + YOffset] > 0))
    {
        RightX++;
        IDs[RightX + YOffset] = ID;
    }

    // Update extents
    MinX = min(MinX,LeftX);
    MaxX = max(MaxX,RightX);

    // Loop from left to right and start recursive fills on line below
    int AboveY = StartY - 1;
    int BelowY = StartY + 1;
    for (int X = LeftX; X <= RightX; X++)
    {
        // Search above
        if (IDs[X + AboveY*W] == -1)
        {
            DoRecursiveFloodFill(X,AboveY,ID,Pix,IDs,W,H,MinX,MaxX,MinY,MaxY);
        }

        // Search Below
        if (IDs[X + BelowY*W] == -1)
        {
            DoRecursiveFloodFill(X,BelowY,ID,Pix,IDs,W,H,MinX,MaxX,MinY,MaxY);
        }
    }

    return true;
}

void DoFloodFill2(int StartX, int StartY, int ID, float* Pix, int* IDs, int W, int H, int& MinX, int& MaxX, int& MinY, int& MaxY)
{
    int StartI = StartX + StartY*W;
    if (IDs[StartI] != -1) return;
    if (Pix[StartI] == 0) return;

    int* SeedX = (int*)malloc(sizeof(int)*W*H);
    int* SeedY = (int*)malloc(sizeof(int)*W*H);
    int nSeeds = 0;

    IDs[StartI] = ID;
    SeedX[nSeeds] = StartX;
    SeedY[nSeeds] = StartY;
    nSeeds++;

    while( nSeeds>0 )
    {
        int CX = SeedX[nSeeds-1];
        int CY = SeedY[nSeeds-1];
        nSeeds--;

        if( (CX>0) && (CX<W-1) && (CY>0) && (CY<H-1) )
        {
            MinY = min(MinY,CY);
            MaxY = max(MaxY,CY);
            MinX = min(MinX,CX);
            MaxX = max(MaxX,CX);

            int I = CX+CY*W;
            if( (IDs[I-1]==-1) && (Pix[I-1]>0) ) {SeedX[nSeeds]=(CX-1);    SeedY[nSeeds]=(CY);        nSeeds++;    IDs[I-1] = ID;}
            if( (IDs[I+1]==-1) && (Pix[I+1]>0) ) {SeedX[nSeeds]=(CX+1);    SeedY[nSeeds]=(CY);        nSeeds++;    IDs[I+1] = ID;}
            if( (IDs[I-W]==-1) && (Pix[I-W]>0) ) {SeedX[nSeeds]=(CX);        SeedY[nSeeds]=(CY-1);    nSeeds++;    IDs[I-W] = ID;}
            if( (IDs[I+W]==-1) && (Pix[I+W]>0) ) {SeedX[nSeeds]=(CX);        SeedY[nSeeds]=(CY+1);    nSeeds++;    IDs[I+W] = ID;}
        }
    }

    free(SeedX);
    free(SeedY);
}

int FloodFill( float* Pix, int W, int H, float* DstPix, int* BlobIDs )
{
    logtime LT("FloodFill");

    int nPixels = W*H;
    //int* BlobIDs = (int*)malloc(sizeof(int)*W*H);
    memset(BlobIDs,0xFFFFFFFF,sizeof(int)*W*H);

    int BlobSequence = 0;
    vector<int> KeepIDs;
    
    float SegYT[] = {0.50f,0.60f,0.70f};
    for( int SegY=0; SegY<3; SegY++)
    {
        int StartY = (int)(SegYT[SegY]*H);

        log("floodfill %d\n",StartY);

        for (int StartX = 0; StartX < W; StartX+=16)
        {
            int StartI = StartX + StartY*W;

            if ((Pix[StartI]>0) && (BlobIDs[StartI]==-1))
            {
                int MinX = StartX;
                int MaxX = StartX;
                int MinY = StartY;
                int MaxY = StartY;

                log("====> %d  %d %d \n",BlobSequence,StartX,StartY);

                // Call recursive filling routine
                //if( DoRecursiveFloodFill(StartX,StartY,BlobSequence,Pix,BlobIDs,W,H,MinX,MaxX,MinY,MaxY) )
                DoFloodFill2(StartX,StartY,BlobSequence,Pix,BlobIDs,W,H,MinX,MaxX,MinY,MaxY);

                int BlobW = MaxX-MinX+1;
                int BlobH = MaxY-MinY+1;
                log("BlobW %d   BlobH %d\n",BlobW,BlobH);
                if( (BlobW>=(W/8)) && (BlobH>=(H/4)) )
                {
                    log("keeping %d\n",BlobSequence);
                    KeepIDs.push_back(BlobSequence);
                }

                
                BlobSequence++;
            }
        }
    }

    log("finished floodfill\n");

    if( KeepIDs.size() > 0 )
    {
        int nKeeps = KeepIDs.size();
        int* pKeeps = &KeepIDs[0];
        for( int i=0; i<W*H; i++ )
        {
            bool Keep = false;
            int KeepIndex = -1;
            int ID = BlobIDs[i];
            if( ID != -1 )
            {
                for( int j=0; j<nKeeps; j++ )
                    if( ID==pKeeps[j] ) {Keep = true; KeepIndex=j; break; }
            }
            BlobIDs[i] = KeepIndex;
            if( Keep ) DstPix[i] = 1;
            else DstPix[i] = 0;
        }
    }

    return KeepIDs.size();
}

//=========================================================================

int IllumFloodFill( float* Pix, int W, int H, float* DstPix, int* BlobIDs, vector<Vec2>& TL, vector<Vec2>& BR )
{
    int nPixels = W*H;
    //int* BlobIDs = (int*)malloc(sizeof(int)*W*H);
    memset(BlobIDs,0xFFFFFFFF,sizeof(int)*W*H);

    int BlobSequence = 0;
    vector<int> KeepIDs;
    
    for( int StartY=50; StartY<H-180; StartY+=4)
    {
        for (int StartX = 0; StartX < W-50; StartX+=4)
        {
            int StartI = StartX + StartY*W;

            if ((Pix[StartI]>0) && (BlobIDs[StartI]==-1))
            {
                int MinX = StartX;
                int MaxX = StartX;
                int MinY = StartY;
                int MaxY = StartY;

                //log("====> %d  %d %d \n",BlobSequence,StartX,StartY);

                DoFloodFill2(StartX,StartY,BlobSequence,Pix,BlobIDs,W,H,MinX,MaxX,MinY,MaxY);

                int BlobW = MaxX-MinX+1;
                int BlobH = MaxY-MinY+1;
                //log("BlobW %d   BlobH %d\n",BlobW,BlobH);
                if( (BlobW>=(6)) && (BlobH>=(6) && (BlobW<=128)) && (BlobH<=128))
                {
                    //log("keeping %d\n",BlobSequence);
                    KeepIDs.push_back(BlobSequence);
                    TL.push_back(Vec2((float)MinX,(float)MinY));
                    BR.push_back(Vec2((float)MaxX,(float)MaxY));
                }
                
                BlobSequence++;
            }
        }
    }

    log("finished floodfill\n");

    if( KeepIDs.size() > 0 )
    {
        int nKeeps = KeepIDs.size();
        int* pKeeps = &KeepIDs[0];
        for( int i=0; i<W*H; i++ )
        {
            bool Keep = false;
            int KeepIndex = -1;
            int ID = BlobIDs[i];
            if( ID != -1 )
            {
                for( int j=0; j<nKeeps; j++ )
                    if( ID==pKeeps[j] ) {Keep = true; KeepIndex=j; break; }
            }
            BlobIDs[i] = KeepIndex;
            if( Keep ) DstPix[i] = 1;
            else DstPix[i] = 0;
        }
    }

    return KeepIDs.size();
}


int DarkFloodFill( float* Pix, int W, int H, float* DstPix, int* BlobIDs, vector<Vec2>& TL, vector<Vec2>& BR )
{
    int nPixels = W*H;
    //int* BlobIDs = (int*)malloc(sizeof(int)*W*H);
    memset(BlobIDs,0xFFFFFFFF,sizeof(int)*W*H);

    int BlobSequence = 0;
    vector<int> KeepIDs;
    
    for( int StartY=50; StartY<H-100; StartY+=4)
    {
        for (int StartX = 0; StartX < W; StartX+=4)
        {
            int StartI = StartX + StartY*W;

            if ((Pix[StartI]>0) && (BlobIDs[StartI]==-1))
            {
                int MinX = StartX;
                int MaxX = StartX;
                int MinY = StartY;
                int MaxY = StartY;

                //log("====> %d  %d %d \n",BlobSequence,StartX,StartY);

                DoFloodFill2(StartX,StartY,BlobSequence,Pix,BlobIDs,W,H,MinX,MaxX,MinY,MaxY);

                int BlobW = MaxX-MinX+1;
                int BlobH = MaxY-MinY+1;
                //log("BlobW %d   BlobH %d\n",BlobW,BlobH);
                if( (BlobW>=(10)) && (BlobH>=(16) && (BlobW<=40)) && (BlobH<=40) && (abs(BlobW-BlobH)<10))
                {
                    //log("keeping %d\n",BlobSequence);
                    KeepIDs.push_back(BlobSequence);
                    TL.push_back(Vec2((float)MinX,(float)MinY));
                    BR.push_back(Vec2((float)MaxX,(float)MaxY));
                }
                
                BlobSequence++;
            }
        }
    }

    log("finished floodfill\n");

    if( KeepIDs.size() > 0 )
    {
        int nKeeps = KeepIDs.size();
        int* pKeeps = &KeepIDs[0];
        for( int i=0; i<W*H; i++ )
        {
            bool Keep = false;
            int KeepIndex = -1;
            int ID = BlobIDs[i];
            if( ID != -1 )
            {
                for( int j=0; j<nKeeps; j++ )
                    if( ID==pKeeps[j] ) {Keep = true; KeepIndex=j; break; }
            }
            BlobIDs[i] = KeepIndex;
            if( Keep ) DstPix[i] = 1;
            else DstPix[i] = 0;
        }
    }

    return KeepIDs.size();
}

//=========================================================================
//=========================================================================
//=========================================================================

struct Eye
{
    string Name;
    int W;
    int H;
    int nPix;
    float* PixR;
    float* PixG;
    float* PixB;
    float* PixI;
    float* PixIEdges;
    float* PixIThresh;
    float* PixDebug;
    float* PixIllum;
    float* PixDark;
    float* PixFFDark;
    int* BlobIDs;
    int nBlobs;
    int Side;
    float* PixIHistogram;
    int Environment;
    EnvironmentInfo EnvInfo;

    Quad PanelQuad;
    float Panel3D[4*3];
    bool PanelFound;

    void Init( string name, vector<int>& Image, int side )
    {
        logtime LT("Init");

        Side = side;
        Name = name;
        H = Image[0];
        W = Image[1];
        nPix = W*H;

        log("Width:  %d\n",W);
        log("Height: %d\n",H);

        PixR = (float*)malloc(sizeof(float)*nPix);
        PixG = (float*)malloc(sizeof(float)*nPix);
        PixB = (float*)malloc(sizeof(float)*nPix);
        PixI = (float*)malloc(sizeof(float)*nPix);
        PixDebug = (float*)malloc(sizeof(float)*nPix);
        PixIEdges = (float*)malloc(sizeof(float)*nPix);
        BlobIDs = (int*)malloc(sizeof(int)*W*H);
        PixIllum = (float*)malloc(sizeof(float)*nPix);
        int* pImage = &Image[0];

        for( int i=0; i<nPix; i++ )
        {
            u32 C = pImage[i+2];
            float R = ((C>>16)&0xFF)*0.00392156862745f;
            float G = ((C>>8)&0xFF)*0.00392156862745f;
            float B = ((C>>0)&0xFF)*0.00392156862745f;
            PixR[i] = R;
            PixG[i] = G;
            PixB[i] = B;
            PixI[i] = (R+G+B) * 0.33333f;
        }

   //     {
            //logtime LT("MedianFilter");
   //         float* P = ClonePix(PixI,W,H);
   //         MedianFilter(P,PixI,W,H);
   //         free(P);
   //     }

        {
            logtime LT("AverageFilter");
            float* P = ClonePix(PixI,W,H);
            AvgFilter(P,PixI,W,H);
            free(P);
        }

        PixIHistogram = ComputeHistogram(256,PixI,W,H);

        PixDebug = ClonePix(PixI,W,H);
        ScalePixels( 0.5f, PixDebug, W, H );

        DetermineEnvironment();
        log("=======================\n");
        log("ENVIRONMENT: %s\n",EnvInfo.Name);
        log("=======================\n");

        CalcIlluminationMask();

        //SaveBMP(EnvInfo.Name,"planeR.bmp",PixR,W,H);
        //SaveBMP(EnvInfo.Name,"planeG.bmp",PixG,W,H);
        //SaveBMP(EnvInfo.Name,"planeB.bmp",PixB,W,H);
        SaveBMP(EnvInfo.Name,"greyscale.bmp",PixI,W,H);
        FindPanel();

        if( PanelFound==false )
            log("!!! PANEL NOT FOUND\n");

    }

    void CalcIlluminationMask()
    {
        logtime LT("CalcIlluminationMask");

        if( EnvInfo.Directory == ENVIRON_ISS )
        {
            for( int i=0; i<nPix; i++ )
            {
                float R = PixR[i];
                float G = PixG[i];
                float B = PixB[i];

                //if( (B>(240.0f/255.0f)) && (max(R,G)<(180.0f/255.0f)) ) PixIllum[i] = 1.0f;
                //else
                if( (G>(245.0f/255.0f)) && (min(R,B)<(180.0f/255.0f)) ) PixIllum[i] = 1.0f;
                else PixIllum[i] = 0;
            }
        }
        else
        if( EnvInfo.Directory == ENVIRON_SIM )
        {
            for( int i=0; i<nPix; i++ )
            {
                float R = PixR[i];
                float G = PixG[i];
                float B = PixB[i];

                //if( (B>(200.0f/255.0f)) && (max(R,G)<(170.0f/255.0f)) ) PixIllum[i] = 1.0f;
                //else
                //if( (G>(230.0f/255.0f)) && (max(R,B)<(175.0f/255.0f)) ) PixIllum[i] = 1.0f;
                //else PixIllum[i] = 0;

                //if( (B>(200.0f/255.0f)) && (max(R,G)<(200.0f/255.0f)) ) PixIllum[i] = 1.0f;
                //else
                //if( (G>(230.0f/255.0f)) && (max(R,B)<(200.0f/255.0f)) ) PixIllum[i] = 1.0f;
                //else PixIllum[i] = 0;

                if( (B>(200.0f/255.0f)) && (max(R,G)<B) ) PixIllum[i] = 1.0f;
                else
                if( (G>(230.0f/255.0f)) && (max(R,B)<G) ) PixIllum[i] = 1.0f;
                else PixIllum[i] = 0;

            }
        }
        else
        if( EnvInfo.Directory == ENVIRON_LAB )
        {
            for( int i=0; i<nPix; i++ )
            {
                float R = PixR[i];
                float G = PixG[i];
                float B = PixB[i];

                if( ((G-R) > 0.20) && ((G-B) > 0.20) ) PixIllum[i] = 1.0f;
                else
                    if( ((min(B,G)-R) > 0.20) ) PixIllum[i] = 1.0f;
                    else PixIllum[i]=0;
            }

        }
        else
        if( EnvInfo.Directory == ENVIRON_LAB3 )
        {
            for( int Y=8; Y<H-8; Y++ )
            for( int X=8; X<W-8; X++ )
            {
                int O = X+Y*W;

                float R = PixR[O];
                float G = PixG[O];
                float B = PixB[O];

                PixIllum[O] = 0;

                if( (G > (250/255.0f)) && (min(R,B) < (200/255.0f)) ) PixIllum[O] = 1;
                else
                    if( (G > (190/255.0f)) && (min(R,B) < (80/255.0f)) ) PixIllum[O] = 1;
                else
                if( ((G-R) > 0.10) && ((G-B) > 0.20) ) PixIllum[O] = 1;
                else
                    if( ((min(B,G)-R) > 0.20) ) PixIllum[O] = 1;
                    else PixIllum[O]=0;

                if( PixIllum[O]==0 )
                {
                    if( (R>0.78f) && (G>0.78f) && (B>0.78f) && ((PixIllum[O-W]>0) && (PixIllum[O-1]>0)) )
                        PixIllum[O] = 1;
                }

                ////if( ((X==834)&&(Y==427)) ||
                ////    ((X==835)&&(Y==427)) ||
                ////    ((X==838)&&(Y==424)) ||
                ////    ((X==837)&&(Y==426)) ||
                ////    ((X==837)&&(Y==425)) ||
                ////    ((X==836)&&(Y==426)) ||
                ////    ((X==834)&&(Y==428)) )
                //if( ((X==833)&&(Y==432))
                //    )
                //{
                //    log("!!!! %d,%d  %1.4f,%1.4f,%1.4f    %1.1f,%1.1f,%1.1f  %1.1f %1.1f %1.1f \n",X,Y,R,G,B,R*255,G*255,B*255,PixIllum[O],PixIllum[O-W],PixIllum[O-1]);
                //}

            }
        }
        else
        {
            for( int i=0; i<nPix; i++ )
            {
                float R = PixR[i];
                float G = PixG[i];
                float B = PixB[i];

                if( ((G-R) > 0.20) && ((G-B) > 0.20) ) PixIllum[i] = 1.0f;
                else
                    if( ((min(B,G)-R) > 0.20) ) PixIllum[i] = 1.0f;
                    else PixIllum[i]=0;
            }
        }

        //ShrinkPix(PixIllum,W,H,0.5f);

        SaveBMP(EnvInfo.Name,"illum.bmp",PixIllum,W,H);
    }

    double AvgBright;
    void DetermineEnvironment()
    {
        logtime LT("DetermineEnvironment");
    
        float WhiteSum = 0;
        for( int i=245; i<=255; i++ ) WhiteSum += PixIHistogram[i];
        float WhiteRatio = WhiteSum / nPix;

        float BlackSum = 0;
        for( int i=0; i<=5; i++ ) BlackSum += PixIHistogram[i];
        float BlackRatio = BlackSum / nPix;

        float TopAvg = 0;
        for( int i=0; i<W; i++ ) TopAvg += PixI[8*W+i];
        TopAvg /= W;

        AvgBright = 0;
        for( int i=0; i<W*H; i++ )
            AvgBright += PixI[i];
        AvgBright /= W*H;

        float BBarRatio=0;
        {
            int MaxBBC = 0;
            int Threshold = (W/8)*3/4;
            for( int Y=0; Y<H; Y+=8 )
            {
                int C=0;
                for( int X=0; X<W; X+=8 )
                    if( PixI[X+Y*W] < 0.10 ) C++;
                if( C > MaxBBC ) MaxBBC = C;
            }
            BBarRatio = MaxBBC / (float)Threshold;
        }

        log("W:           %d\n",W);
        log("WhiteRatio:  %f\n",WhiteRatio );
        log("BlackRatio:  %f\n",BlackRatio );
        log("TopAvg:      %f\n",TopAvg );
        log("AvgBright:   %f\n",(float)AvgBright);
        log("BlackBar :   %f\n",BBarRatio);

        Environment = 0;
        if( W>1600 )
        {
            if( AvgBright < 0.20 ) Environment = ENVIRON_ISS;
            else
            {
                if( WhiteRatio > 0.05 ) Environment=ENVIRON_ISS;
                else Environment = ENVIRON_LAB;
            }
        }
        else
        {   
            
            if( BlackRatio<0.01 )
            {
                Environment = ENVIRON_LAB3;
            }
            else
            {
                if( TopAvg < 0.25 ) Environment = ENVIRON_SIM;
                else
                if( TopAvg > 0.30 ) Environment = ENVIRON_LAB2;
                else
                {
                    if(WhiteRatio > 0.005 ) Environment = ENVIRON_SIM;
                    else Environment = ENVIRON_LAB2;
                }
            }

            //// Changed due to WYOF
            //if( (TopAvg < 0.225) || ((TopAvg<0.30) && (BlackRatio>0.1f) && (BlackRatio<0.35)) ) Environment = ENVIRON_SIM;
            ////if( (TopAvg < 0.25) || ((TopAvg<0.30) && (BlackRatio>0.1f) && (BlackRatio<0.35)) ) Environment = ENVIRON_SIM;
   //         else
   //         {
   //             if( BBarRatio > 1.0f ) Environment = ENVIRON_LAB2;
   //             else Environment = ENVIRON_LAB3;
   //         }
        }

        //{
        //    FILE* fp = fopen("C:\\Projects\\TC_Robonaut\\Environments.csv","at");
        //    fprintf(fp,"%d,%d,%f,%f,%f,%f,%f,\n",(int)Environment,W,WhiteRatio,BlackRatio,TopAvg,(float)AvgBright,BBarRatio);
        //    fclose(fp);
        //}

        EnvInfo = EnvironmentInfos[ Environment*2 + Side ];

        if( Environment == ENVIRON_SIM )
        {
            // Patch power cover positions
            Items[ITEM_PANEL_POWER_COVER].RefPos =  Items[ITEM_PANEL_POWER_COVER_SIM].RefPos;
            Items[ITEM_PANEL_POWER_COVER].PanelPos =  Items[ITEM_PANEL_POWER_COVER_SIM].PanelPos;
            Items[ITEM_PANEL_POWER_COVER].ZOffset =  Items[ITEM_PANEL_POWER_COVER_SIM].ZOffset;
            Items[ITEM_PANEL_POWER_COVER_UP].RefPos =  Items[ITEM_PANEL_POWER_COVER_UP_SIM].RefPos;
            Items[ITEM_PANEL_POWER_COVER_UP].PanelPos =  Items[ITEM_PANEL_POWER_COVER_UP_SIM].PanelPos;
            Items[ITEM_PANEL_POWER_COVER_UP].ZOffset =  Items[ITEM_PANEL_POWER_COVER_UP_SIM].ZOffset;
        }
    }

    bool IsQuadShapedLikePanel( Quad Q )
    {
        int PanelShapErr = 0;
        float YBorder = -500;
        float XBorder = -500;
        for( int i=0; i<4; i++ )
        {
            if( (Q.Pts[i].X<XBorder) || 
                (Q.Pts[i].X>W-XBorder) || 
                (Q.Pts[i].Y<YBorder) || 
                (Q.Pts[i].Y>H-YBorder))
                PanelShapErr = 1;
        }

        float L0 = Q.Pts[0].Dist(Q.Pts[1]);
        float L1 = Q.Pts[3].Dist(Q.Pts[2]);
        if( (L0<50) || (L1<50) ) PanelShapErr = 4;
        float Ratio = max(L1,L0) / min(L1,L0);
        log("panelshape LR ratio: %f\n",Ratio);
        if( Ratio > 1.8f ) PanelShapErr = 2;

        L0 = Q.Pts[0].Dist(Q.Pts[3]);
        L1 = Q.Pts[1].Dist(Q.Pts[2]);
        if( (L0<50) || (L1<50) ) PanelShapErr = 5;
        Ratio = max(L1,L0) / min(L1,L0);
        log("panelshape TB ratio: %f\n",Ratio);
        if( Ratio > 2.5f ) PanelShapErr = 3;

        log("PanelShapeErr: %d\n",PanelShapErr);
        return (PanelShapErr==0);
    }

    int CalcColorMatchScore(EnvironmentInfo& EI, float* Q3D )
    {
        //log("CalcColorMatchScore\n");
        int PwrCapScore = 0;
        {
            Vec2 Center = CalcPerspectiveAwareBilinearPos(EnvInfo,Q3D,Items[ITEM_PANEL_POWER_SWITCH].RefPos.X,Items[ITEM_PANEL_POWER_SWITCH].RefPos.Y,Items[ITEM_PANEL_POWER_SWITCH].ZOffset);
            int Radius = 64;
            int Step = 4;
            int CX = (int)Center.X;
            int CY = (int)Center.Y;

            DrawLine(PixDebug,W,H,1,CX-Radius,CY-Radius,CX+Radius,CY+Radius);
            DrawLine(PixDebug,W,H,1,CX-Radius,CY+Radius,CX+Radius,CY-Radius);

            int YC=0;
            int NC=0;
            for( int Y=CY-Radius; Y<=CY+Radius; Y+=Step )
                for( int X=CX-Radius; X<=CX+Radius; X+=Step )
                {
                    if( (Y>=0) && (Y<H) && (X>=0) && (X<W) )
                    {
                        int O = X+Y*W;
                        float R = PixR[O];
                        float G = PixG[O];
                        float B = PixB[O];
                        float GB = (G+B)*0.5f;
                        float Ratio = 0;
                        if( GB>0 ) Ratio = R/GB;
                    
                        if( (R>G) && (R>B) && (Ratio > EI.PwrCapRedRatio) ) YC++;
                        else NC++;
                    }
                }

            if( (YC+NC) > 0 )
                PwrCapScore = (YC*100/(YC+NC));
            //log("PwrCapScore %f  (%f,%f)  %d %d %d\n",EI.PwrCapRedRatio,Center.X,Center.Y,YC,NC,PwrCapScore);
        }
    
        int IllumBtnScore = 0;
        {
            int nChecks = 0;
            for( int i=0; i<22; i++ )
            {
                // Check LEDs
                if( Items[i].StateMask[STATE_ON] == 1 )
                {
                    nChecks++;

                    float SearchRadius = 0;
                    if((i>=ITEM_A02_LED_NUM_PAD_A1) && (i<=ITEM_A02_LED_NUM_PAD_C3) ) SearchRadius = 20.0f / 518;
                    else SearchRadius = 8.0f / 518;

                    Vec2 Center = CalcPerspectiveAwareBilinearPos(EnvInfo,Q3D,Items[i].RefPos.X,Items[i].RefPos.Y,Items[i].ZOffset);
                    Vec2 RRef = CalcPerspectiveAwareBilinearPos(EnvInfo,Q3D,Items[i].RefPos.X,Items[i].RefPos.Y-SearchRadius,Items[i].ZOffset);
                    int ImgRadius = (int)RRef.Dist(Center);
                    if( ImgRadius > 100 ) ImgRadius = 100;
                    if( ImgRadius < 4 ) ImgRadius = 4;
                    int ImgRadius2 = ImgRadius*ImgRadius;

                    //DrawCircle(PixDebug,W,H,1.0f,Center.X,Center.Y,ImgRadius);

                    int NC=0;
                    int YC=0;

                    for( int DY=-ImgRadius; DY<=ImgRadius; DY++ )
                        for( int DX=-ImgRadius; DX<=ImgRadius; DX++ )
                        {
                            float D2 = DX*DX + DY*DY;
                            if( D2 <= ImgRadius2 )
                            {
                                int X = (int)(Center.X + DX);
                                int Y = (int)(Center.Y + DY);
                                if( (X>=0) && (X<W) && (Y>=0) && (Y<H) )
                                {
                                    float V = PixIllum[X+Y*W];
                                    if( V > 0.5f ) YC++; else NC++;
                                }
                            }
                        }

                    if( (YC+NC) > 0 )
                        IllumBtnScore += (YC*100/(YC+NC));
                }
            }
            //IllumBtnScore /= nChecks;
        }

        int FinalScore = PwrCapScore*10 + IllumBtnScore;
        log("============= CalcColorMatchScore %d %d %d\n",PwrCapScore,IllumBtnScore,FinalScore);
        return FinalScore;
    }

    bool DoesQuadLookLikePanel( Quad Q )
    {

        return true;
    }

    float SolveThresholdFromPanelLights()
    {
        //float* PixIllumShrunk = ClonePix(PixIllum,W,H);
        //ShrinkPix( PixIllumShrunk, W, H, 0.5f );
        
        //int nBins = 255;
        //float* Bins = (float*)malloc(sizeof(float)*nBins);
        //memset(Bins,0,sizeof(float)*nBins);

        float Avg = 0;
        float Max = 0;
        float Min = 10000;
        float AvgC = 0;
        for( int Y=4; Y<=H-4; Y+=2 )
            for( int X=4; X<=W-4; X+=2 )
            {
                int O = X+Y*W;
                if( (PixIllum[O]>0) && (PixIllum[O-1]>0) && (PixIllum[O+1]>0) )
                {
                    int R = 32;
                    for( int HY=Y-R; HY<=Y+R; HY+=2 )
                        for( int HX=X-R; HX<=X+R; HX+=2 )
                        {
                            if( (HX>=0) && (HX<W) && (HY>=0) && (HY<H) && (PixIllum[HX+HY*W]==0))
                            {
                                float V = PixI[HX+HY*W];

                                Avg += V;
                                AvgC++;
                                Max = max(Max,V);
                                Min = min(Min,V);

                                //int B = (V*(nBins-1)+0.5f);
                                //if( B<0 ) B=0;
                                //if( B>=nBins ) B=nBins-1;
                                //Bins[B]++;
                            }
                        }
                }
            }

        if( AvgC>0 ) Avg /= AvgC;

        log("SolveThresholdFromPanelLights %f %f %f %f\n",Avg,Max,Min,AvgC);

        //free(PixIllumShrunk);
        //free(Bins);

        return Avg;
    }

    ////float SolveThresholdFromPanelLights()
    ////{
    ////    //float* PixIllumShrunk = ClonePix(PixIllum,W,H);
    ////    //ShrinkPix( PixIllumShrunk, W, H, 0.5f );
    ////    
    ////    //int nBins = 255;
    ////    //float* Bins = (float*)malloc(sizeof(float)*nBins);
    ////    //memset(Bins,0,sizeof(float)*nBins);

    ////    float Avg = 0;
    ////    float Max = 0;
    ////    float Min = 10000;
    ////    float AvgC = 0;
    ////    for( int Y=175; Y<=1000; Y+=2 )
    ////        for( int X=100; X<=1500; X+=2 )
    ////        {
    ////            int O = X+Y*W;
    ////            if( (PixIllum[O]>0) && 
                ////    (PixIllum[O-2]>0) && (PixIllum[O+2]>0) &&
                ////    (PixIllum[O-W*2]>0) && (PixIllum[O+W*2]>0) )
    ////            {
    ////                int R = 32;
    ////                for( int HY=Y-R; HY<=Y+R; HY+=2 )
    ////                    for( int HX=X-R; HX<=X+R; HX+=2 )
    ////                    {
    ////                        if( (HX>=0) && (HX<W) && (HY>=0) && (HY<H) && (PixIllum[HX+HY*W]==0))
    ////                        {
    ////                            float V = PixI[HX+HY*W];

    ////                            Avg += V;
    ////                            AvgC++;
    ////                            Max = max(Max,V);
    ////                            Min = min(Min,V);

    ////                            //int B = (V*(nBins-1)+0.5f);
    ////                            //if( B<0 ) B=0;
    ////                            //if( B>=nBins ) B=nBins-1;
    ////                            //Bins[B]++;
    ////                        }
    ////                    }
    ////            }
    ////        }

    ////    if( AvgC>0 ) Avg /= AvgC;

    ////    log("SolveThresholdFromPanelLights %f %f %f %f\n",Avg,Max,Min,AvgC);

    ////    //free(PixIllumShrunk);
    ////    //free(Bins);

    ////    return Avg;
    ////}

    int IllumButtonWidth;
    int IllumButtonHeight;

    void FindLab3AnchorsAndInliers( vector<Vec2>& ButtonIllums, vector<Vec2>& LEDIllums, vector<Vec2>& LEDDarks )
    {
        SinglePixelNoiseRemoval(PixIllum,W,H);
        SinglePixelNoiseRemoval(PixIllum,W,H);

        vector<Vec2> BlobMinPos;
        vector<Vec2> BlobMaxPos;
        float* PixFF = ClonePix(PixIllum,W,H);
        nBlobs = IllumFloodFill( PixIllum, W, H, PixFF,BlobIDs,BlobMinPos,BlobMaxPos );

        if( nBlobs > 10 ) nBlobs = 10;

        

        //vector<Vec2> BlobImgCenters;
        for( int i=0; i<nBlobs; i++ )
        {
            int BlobW = BlobMaxPos[i].X-BlobMinPos[i].X;
            int BlobH = BlobMaxPos[i].Y-BlobMinPos[i].Y;
            float CX = (BlobMinPos[i].X + BlobMaxPos[i].X)*0.5f;
            float CY = (BlobMinPos[i].Y + BlobMaxPos[i].Y)*0.5f;
            //BlobImgCenters.push_back( Vec2(CX,CY) );

            DrawLine(PixFF,W,H,1,CX-32,CY,CX+32,CY);
            DrawLine(PixFF,W,H,1,CX,CY-32,CX,CY+32);

            if( BlobW < 16 ) DrawCircle(PixFF,W,H,1,CX,CY,32);

            if( BlobW < 16 ) LEDIllums.push_back(Vec2(CX,CY));
            else
            if( PixIllum[ (int)CX + ((int)CY)*W ] == 0 )
            {
                ButtonIllums.push_back(Vec2(CX,CY));
                IllumButtonWidth = BlobW; 
                IllumButtonHeight = BlobH;
            }
        }

        SaveBMP(EnvInfo.Name,"ff.bmp",PixFF,W,H);

        //float* PixEdges = ClonePix(PixI,W,H);
        //ComputeEdges(PixEdges,PixI,W,H);
        //SaveBMP(EnvInfo.Name,"edges1.bmp",PixEdges,W,H);

        // If we don't have enough anchors...game over
        if( ButtonIllums.size()==0 )
            return;

        // Check for dark inliers
        {
            //AvgFilter(PixI,PixDark,W,H);

            float AvgIllum = 0;
            float AvgC = 0;
            {
                for( int i=0; i<ButtonIllums.size(); i++ )
                {
                    int CX = (int)ButtonIllums[i].X;
                    int CY = (int)ButtonIllums[i].Y;
                    for( int Y=CY-35; Y<=CY+35; Y++ )
                    for( int X=CX-35; X<=CX+35; X++ )
                    {
                        if( (Y>=0) && (Y<H) && (X>=0) && (X<W) )
                        {
                            AvgIllum += PixI[X+Y*W];
                            AvgC++;
                        }
                    }
                }

                AvgIllum /= AvgC;
            }

            log("AvgIllum %f\n",AvgIllum);

            // Find walls guided by anchors
            int WallL = 100000;
            int WallR = 0;
            for( int i=0; i<ButtonIllums.size(); i++ )
            {
                WallL = min((int)ButtonIllums[i].X-150,WallL);
                WallR = max((int)ButtonIllums[i].X+350,WallR);
            }
            for( int i=0; i<LEDIllums.size(); i++ )
            {
                WallL = min((int)LEDIllums[i].X-50,WallL);
                WallR = max((int)LEDIllums[i].X+50,WallR);
            }

            if( WallL<0 ) WallL = 0;
            if( WallR>W-32 ) WallR = W-32;

            log("DarkWalls %d %d\n",WallL,WallR);

            PixDark = AllocPix(W,H);

            //float Threshold[] = {50/255.0f,40/255.0f,20/255.0f};
            bool GoodThreshold = false;
            //for( int iThresh=0; iThresh<3; iThresh++ )
            {
                //float Thresh = Threshold[iThresh];
                float Thresh = AvgIllum * 0.5f;

                int nPixels = 0;
                for( int Y=32; Y<H-32; Y++ )
                for( int X=WallL; X<=WallR; X++ )
                {
                    int O = X+Y*W;

                    float R = PixR[O];
                    float G = PixG[O];
                    float B = PixB[O];

                    PixDark[O] = 0;

                    //if( (R<Thresh) && (G<Thresh) && (B<Thresh) ) PixDark[O] = 1;
                    if( ((R+G+B)*0.3333f)<Thresh)
                    {
                        int A = O+24;
                        int B = O-24;
                        int C = O-W*24;
                        int D = O+W*24;

                        if( (((PixR[A] + PixG[A] + PixB[A])*0.3333f) > Thresh) &&
                            (((PixR[B] + PixG[B] + PixB[B])*0.3333f) > Thresh) &&
                            (((PixR[C] + PixG[C] + PixB[C])*0.3333f) > Thresh) &&
                            (((PixR[D] + PixG[D] + PixB[D])*0.3333f) > Thresh)
                            )
                            PixDark[O] = 1;
                    }

                    nPixels++;
                }

                SinglePixelNoiseRemoval(PixDark,W,H);
                SinglePixelNoiseRemoval(PixDark,W,H);
                GrowPix(PixDark,W,H,0.5f);
                GrowPix(PixDark,W,H,0.5f);
                SaveBMP(EnvInfo.Name,"dark.bmp",PixDark,W,H);

                int YC=0;
                int NC=0;
                for( int Y=8; Y<H-8; Y++ )
                for( int X=WallL; X<=WallR; X++ )
                {
                    if( PixDark[X+Y*W] != 0 ) YC++;
                    else NC++;
                }

                float Ratio = (float)YC/NC;
                log("PixC %d %d %f\n",YC,NC,Ratio);

                //if( (Ratio>=0.002f) && (Ratio<0.10f) )
                //if( ((iThresh==0) && (YC < ((YC+NC)/5))) ||
                //    ((iThresh==1) && (YC < ((YC+NC)/5))) ||
                //    ((iThresh==2) && (YC < ((YC+NC)/5))))
                {
                    GoodThreshold=true;
                    //break;
                }
            }

            //if( GoodThreshold )
            {
                vector<Vec2> BlobMinPos;
                vector<Vec2> BlobMaxPos;
                PixFFDark = ClonePix(PixDark,W,H);
                nBlobs = DarkFloodFill( PixDark, W, H, PixFFDark,BlobIDs,BlobMinPos,BlobMaxPos );

                if( nBlobs > 25 ) nBlobs = 25;

                vector<Vec2> BlobImgCenters;
                for( int i=0; i<nBlobs; i++ )
                {
                    int BlobW = BlobMaxPos[i].X-BlobMinPos[i].X;
                    int BlobH = BlobMaxPos[i].Y-BlobMinPos[i].Y;
                    float CX = (BlobMinPos[i].X + BlobMaxPos[i].X)*0.5f;
                    float CY = (BlobMinPos[i].Y + BlobMaxPos[i].Y)*0.5f;
                    BlobImgCenters.push_back( Vec2(CX,CY) );

                    //DrawLine(PixFFDark,W,H,1,CX-32,CY,CX+32,CY);
                    //DrawLine(PixFFDark,W,H,1,CX,CY-32,CX,CY+32);

                    LEDDarks.push_back(Vec2(CX,CY));
                }
                SaveBMP(EnvInfo.Name,"ffdark.bmp",PixFFDark,W,H);
            }

        }
    }



    void FindPanelLab3()
    {
        PanelFound = false;
        vector<Vec2> ButtonIllums;
        vector<Vec2> LEDIllums;
        vector<Vec2> LEDDarks;

        
        FindLab3AnchorsAndInliers( ButtonIllums, LEDIllums, LEDDarks );
        log("FindPanelLab3 BI:%d LI:%d LD:%d\n",ButtonIllums.size(),LEDIllums.size(),LEDDarks.size());


        // Find the three toggle switch LEDs
        bool ToggleLEDsFound = false;
        Vec2 TLEDPos[3];
        float ButtonXBounds[2] = {0,0};
        float ButtonYBounds[2] = {0,0};
        {
            vector<Vec2> Candidates;
            vector<bool> ValidForB;

            // find lowest Y of ButtonIllums
            float BIllumMaxY = 0;
            float BIllumMaxX = 0;
            for( int i=0; i<ButtonIllums.size(); i++ ) BIllumMaxY = max(BIllumMaxY,ButtonIllums[i].Y);
            for( int i=0; i<ButtonIllums.size(); i++ ) BIllumMaxX = max(BIllumMaxX,ButtonIllums[i].X);
            for( int i=0; i<LEDDarks.size(); i++ ) if( LEDDarks[i].Y > BIllumMaxY ) Candidates.push_back(LEDDarks[i]);
            for( int i=0; i<LEDIllums.size(); i++ ) if( LEDIllums[i].Y > BIllumMaxY ) Candidates.push_back(LEDIllums[i]);

            int nCandidates = Candidates.size();
            if( nCandidates<3 ) return;


            
            {
                float BestErr = 10000000;

                int B=0;
                {
                    float BestDist=10000000;
                    for( int i=0; i<nCandidates; i++ )
                    {
                        if( (Candidates[i].X > BIllumMaxX) && (Candidates[i].Y > (BIllumMaxY+80)))
                        {
                            ValidForB.push_back(true);
                            float Dst = (Candidates[i].Y - BIllumMaxY)*(Candidates[i].Y - BIllumMaxY) + 
                                        (Candidates[i].X - BIllumMaxX)*(Candidates[i].X - BIllumMaxX);
                            if( Dst<BestDist )
                            {
                                BestDist = Dst;
                                B = i;
                            }
                        }
                        else
                            ValidForB.push_back(false);
                    }
                }

                float ACMaxXDist = IllumButtonWidth*7.2 + 150;
                log("ACMaxXDist %f\n",ACMaxXDist);


                // Loop through all combinations
                for( int A=0; A<nCandidates; A++ )
                    //for( int B=0; B<nCandidates; B++ )
                        for( int C=0; C<nCandidates; C++ )
                            if( ValidForB[B] && (A!=B) && (B!=C) && (C!=A) )
                            {
                                if( (Candidates[A].X < Candidates[B].X) && 
                                    (Candidates[B].X < Candidates[C].X) &&
                                    ((Candidates[C].X-Candidates[A].X)<ACMaxXDist))
                                {
                                    // Check if the true B is between A and B
                                    {
                                        int X = (int)(Candidates[A].X + Candidates[B].X)/2;
                                        int Y = (int)(Candidates[A].Y + Candidates[B].Y)/2;
                                        if( PixFFDark[X+Y*W] == 1 )
                                            continue;
                                    }

                                    float DstAB = Candidates[A].Dist(Candidates[B]);
                                    float DstBC = Candidates[B].Dist(Candidates[C]);
                                    float DstAC = Candidates[A].Dist(Candidates[C]);
                                    float Err = fabs(DstBC-DstAB) + fabs(DstAC-(DstAB+DstBC)) + fabs(Candidates[C].Y-Candidates[A].Y) + fabs(Candidates[B].Y-Candidates[A].Y);
                                    if( Err < BestErr )
                                    {
                                        TLEDPos[0] = Candidates[A];
                                        TLEDPos[1] = Candidates[B];
                                        TLEDPos[2] = Candidates[C];
                                        BestErr = Err;
                                        ToggleLEDsFound = true;
                                    }
                                }
                            }

                if( !ToggleLEDsFound )
                {
                    log("!!! ToggleLEDs Not Found\n");
                    return;
                }

                DrawLine(PixDebug,W,H,1,TLEDPos[0].X,TLEDPos[0].Y-2,TLEDPos[2].X,TLEDPos[2].Y-2);
                DrawLine(PixDebug,W,H,1,TLEDPos[0].X,TLEDPos[0].Y+2,TLEDPos[2].X,TLEDPos[2].Y+2);

    //            float Dst01 = TLEDPos[0].Dist(TLEDPos[1]);
                //float Dst02X = TLEDPos[2].X - TLEDPos[0].X;
    //            ButtonXBounds[0] = TLEDPos[0].X + Dst01*0.2f;
    //            ButtonXBounds[1] = TLEDPos[0].X + Dst01*0.55f;
    //            //ButtonYBounds[0] = TLEDPos[0].Y - Dst01*1.39f;
    //            //ButtonYBounds[1] = TLEDPos[0].Y - Dst01*0.92f;
    //            //ButtonYBounds[0] = TLEDPos[0].Y - 140 - 64;//Dst01*1.45f;
    //            //ButtonYBounds[1] = TLEDPos[0].Y - 140;
    //            ButtonYBounds[0] = TLEDPos[0].Y - Dst02X*0.97f;//1.00f;
    //            ButtonYBounds[1] = TLEDPos[0].Y - Dst02X*0.67f;//0.70f;


                float Dst01 = TLEDPos[0].Dist(TLEDPos[1]);
                float Dst02X = TLEDPos[2].X - TLEDPos[0].X;
                ButtonXBounds[0] = TLEDPos[0].X + Dst01*0.20f;
                ButtonXBounds[1] = TLEDPos[0].X + Dst01*0.55f;
                ButtonYBounds[0] = TLEDPos[0].Y - IllumButtonHeight*5.6f;
                ButtonYBounds[1] = TLEDPos[0].Y - IllumButtonHeight*3.9f;

                //IllumButtonWidth = BlobW; 
                //IllumButtonHeight = BlobH;

                DrawLine(PixDebug,W,H,1,ButtonXBounds[0],0,ButtonXBounds[0],H);
                DrawLine(PixDebug,W,H,1,ButtonXBounds[1],0,ButtonXBounds[1],H);
                DrawLine(PixDebug,W,H,1,0,ButtonYBounds[0],W,ButtonYBounds[0]);
                DrawLine(PixDebug,W,H,1,0,ButtonYBounds[1],W,ButtonYBounds[1]);
            }
        }


        // Fill out positions of the various reference points
        vector<Vec2> ItemImagePos;
        vector<Vec3> ItemPanelPos;
        vector<int> ItemIndex;
        {
            // Add power led
            {
                Vec2 Pos;
                int Count=0;
                for( int i=0; i<LEDIllums.size(); i++ )
                {
                    float X = LEDIllums[i].X;
                    float Y = LEDIllums[i].Y;
                    if( (X>(ButtonXBounds[1]+75)) && (Y<(ButtonYBounds[0]-150)) ) 
                    {
                        Pos = LEDIllums[i];
                        Count++;
                    }
                }

                if( Count==1 )
                {
                    ItemImagePos.push_back(Pos);
                    ItemIndex.push_back(ITEM_PANEL_POWER_LED);
                    log("pushed power LED\n");
                }
            }

            //if( ItemIndex.size() ==0 )
            {
                // Add illuminated buttons
                for( int i=0; i<ButtonIllums.size(); i++ )
                {
                    float BX = ButtonIllums[i].X;
                    float BY = ButtonIllums[i].Y;
                    int ItemI = -1;
                    if( (BX<ButtonXBounds[0]) && (BY<ButtonYBounds[0]) ) ItemI = ITEM_A02_LED_NUM_PAD_A1; else
                    if( (BX<ButtonXBounds[1]) && (BY<ButtonYBounds[0]) ) ItemI = ITEM_A02_LED_NUM_PAD_A2; else
                    if(                             (BY<ButtonYBounds[0]) ) ItemI = ITEM_A02_LED_NUM_PAD_A3; else
                    if( (BX<ButtonXBounds[0]) && (BY<ButtonYBounds[1]) ) ItemI = ITEM_A02_LED_NUM_PAD_B1; else
                    if( (BX<ButtonXBounds[1]) && (BY<ButtonYBounds[1]) ) ItemI = ITEM_A02_LED_NUM_PAD_B2; else
                    if(                             (BY<ButtonYBounds[1]) ) ItemI = ITEM_A02_LED_NUM_PAD_B3; else
                    if( (BX<ButtonXBounds[0])                          ) ItemI = ITEM_A02_LED_NUM_PAD_C1; else
                    if( (BX<ButtonXBounds[1])                          ) ItemI = ITEM_A02_LED_NUM_PAD_C2; else
                                                                         ItemI = ITEM_A02_LED_NUM_PAD_C3;

                    ItemImagePos.push_back(ButtonIllums[i]);
                    ItemIndex.push_back(ItemI);

                    //break;
                }
            }

            // Add toggle LEDS
            ItemIndex.push_back(ITEM_A03_LED);
            ItemIndex.push_back(ITEM_A04_LED_TOP);
            ItemIndex.push_back(ITEM_A05_LED);
            ItemImagePos.push_back(TLEDPos[0]);
            ItemImagePos.push_back(TLEDPos[1]);
            ItemImagePos.push_back(TLEDPos[2]);


            // Fill out panel positions
            for( int i=0; i<ItemIndex.size(); i++ )
            {
                ItemPanelPos.push_back( Items[ItemIndex[i]].PanelPos );
            }
        }

        float Q3D[4*3];
        if( !SolvePanelPoseFromItemsNew( EnvInfo, ItemImagePos, ItemPanelPos, Q3D ) )
        {
            log("SolvePanelPoseFromItemsNew...FAILED\n");
            if( !SolvePanelPoseFromItemsOLD( EnvInfo, ItemImagePos, ItemPanelPos, Q3D ) )
            {
                log("SolvePanelPoseFromItemsOLD...FAILED\n");
                return;
            }
        }
        log("SolvedPanelPoseFromItems...\n");

        memcpy(Panel3D,Q3D,sizeof(float)*12);
        PanelFound = true;

        for( int j=0; j<22; j++ )
        {
            Vec2 Center = CalcPerspectiveAwareBilinearPos(EnvInfo,Q3D,Items[j].RefPos.X,Items[j].RefPos.Y,Items[j].ZOffset);
            DrawCircle(PixDebug,W,H,1,Center.X,Center.Y,8);
        }

        Vec2 Corners[4] = {Vec2(0,0), Vec2(0,1), Vec2(1,1), Vec2(1,0)};
        for( int j=0; j<4; j++ )
        {
            Vec2 Center = CalcPerspectiveAwareBilinearPos(EnvInfo,Q3D,Corners[j].X,Corners[j].Y,0);
            DrawCircle(PixDebug,W,H,1,Center.X,Center.Y,8);
            DrawCircle(PixDebug,W,H,1,Center.X,Center.Y,9);
            DrawCircle(PixDebug,W,H,1,Center.X,Center.Y,10);
            DrawCircle(PixDebug,W,H,1,Center.X,Center.Y,11);
            DrawLine(PixDebug,W,H,1,Center.X-8,Center.Y,Center.X+8,Center.Y);
            DrawLine(PixDebug,W,H,1,Center.X,Center.Y-8,Center.X,Center.Y+8);
        }

        log("FindPanelLab3 completed.\n");
    }

    bool SolvePanelPoseFromItemsOLD( EnvironmentInfo& EI, vector<Vec2> ItemImagePos, vector<Vec3> ItemPanelPos, float* Q3D )
    {
        vector<Vec3> RayDir;
        vector<Vec3> CamPos;
        vector<Vec3> PanelPos;
        vector<Vec2> ImagePos;

        // Choose the 3 items that create the triangle with the largest
        // area in image space
        int BestI[3];
        {
            float BestArea = 0;
            for( int A=0; A<ItemImagePos.size(); A++ )
                for( int B=0; B<ItemImagePos.size(); B++ )
                    for( int C=0; C<ItemImagePos.size(); C++ )
                        if( (A!=B) && (B!=C) && (C!=A) )
                        {
                            float a = ItemImagePos[A].Dist( ItemImagePos[B] );
                            float b = ItemImagePos[B].Dist( ItemImagePos[C] );
                            float c = ItemImagePos[C].Dist( ItemImagePos[A] );
                            float p = (a+b+c)*0.5f;
                            float Area = sqrtf(p*(p-a)*(p-b)*(p-c));
                            if( Area > BestArea )
                            {
                                BestArea = Area;
                                BestI[0] = A;
                                BestI[1] = B;
                                BestI[2] = C;
                            }
                        }
        }

        for( int i=0; i<3; i++ )
            DrawLine(PixDebug,W,H,1,ItemImagePos[BestI[i]].X,ItemImagePos[BestI[i]].Y,ItemImagePos[BestI[(i+1)%3]].X,ItemImagePos[BestI[(i+1)%3]].Y);


        int nItems = 3;
        for( int i=0; i<3; i++ )
        {
            Vec3 Dir;
            ImageToRay( EI,ItemImagePos[BestI[i]].X, ItemImagePos[BestI[i]].Y, Dir.X, Dir.Y );
            RayDir.push_back(Dir);
            CamPos.push_back(Vec3(0,0,0));
            PanelPos.push_back(ItemPanelPos[BestI[i]]);
            ImagePos.push_back(ItemImagePos[BestI[i]]);
        }

        // Compute Avg internal distance in panel
        double AvgInternalPanelDist = 0;
        double AvgInternalImageDist = 0;
        {
            int C = 0;
            for( int i=0; i<nItems; i++ )
                for (int j = i+1; j < nItems; j++)
                {
                    double Dist = PanelPos[i].Dist(PanelPos[j]);
                    AvgInternalPanelDist += Dist;

                    double ImgDist = ImagePos[i].Dist(ImagePos[j]);
                    AvgInternalImageDist += ImgDist;

                    C++;
                }

            if (C == 0) return false;
            AvgInternalPanelDist /= C;
            AvgInternalImageDist /= C;
        }

        // Initial Z estimate 
        double InitialZ = (EI.Focus / EI.SX) * (AvgInternalPanelDist / AvgInternalImageDist);
        log("InitialZ %f\n",InitialZ);
        for (int i = 0; i < nItems; i++)
        {
            CamPos[i].X = RayDir[i].X * InitialZ;
            CamPos[i].Y = RayDir[i].Y * InitialZ;
            CamPos[i].Z = InitialZ;
        }

        if( nItems>3 ) nItems=3; 

        double BestErr = 100000000.0;
        vector<Vec3> NewCamPos(nItems);
        //vector<Vec3> NewImgPos(nItems);
        for (int nTrials = 0; nTrials < 10000; nTrials++)
        {
            for (int i = 0; i < nItems; i++)
                NewCamPos[i].Z = CamPos[i].Z;

            // Randomize NewZs
            int I = myrand() % nItems;
            double Amt = -0.020f + ((myrand()%10000)/10000.0) * (0.040f);
            NewCamPos[I].Z += Amt;

            ////!!! Force all the same Zs
            //for( int i=0; i<nItems; i++ ) NewCamPos[i].Z = NewCamPos[I].Z;

            // Calc new positions
            for (int i = 0; i < nItems; i++)
            {
                NewCamPos[i].X = RayDir[i].X * NewCamPos[i].Z;
                NewCamPos[i].Y = RayDir[i].Y * NewCamPos[i].Z;
                //WorldToImage(EI,NewCamPos[i].X,NewCamPos[i].Y,NewCamPos[i].Z,NewImgPos[i].X,NewImgPos[i].Y);
            }

            // Sum internal position errors
            double CamPosErr = 0;
            double nChecks = 0;
            {
                for (int i = 0; i < nItems; i++)
                for (int j = i+1; j < nItems; j++)
                {
                    double PanelDist = PanelPos[i].Dist(PanelPos[j]);
                    double CamDist = NewCamPos[i].Dist(NewCamPos[j]);
                    CamPosErr += (CamDist - PanelDist) * (CamDist - PanelDist);
                    nChecks++;
                }
            }
            //CamPosErr /= nChecks;
            CamPosErr = sqrtf(CamPosErr / nChecks);

            if (CamPosErr < BestErr)
            {
                BestErr = CamPosErr;
                for (int i = 0; i < nItems; i++)
                    CamPos[i] = NewCamPos[i];
            }
            //if( (nTrials%100)==0 )
                //log("[%05d] SolvePanelPoseFromItems %1.4f\n",nTrials,BestErr);
        }

        for( int i=0; i<nItems; i++ )
            log("CamPos %d   %f,%f,%f\n",i,CamPos[i].X,CamPos[i].Y,CamPos[i].Z);

        // CamPos contains best estimates of item positions
        memset(Q3D,0,sizeof(float)*12);
        Vec2 PanelRefCorners[4] = {Vec2(0,0), Vec2(0.0f,PANELHEIGHT), Vec2(PANELWIDTH,PANELHEIGHT), Vec2(PANELWIDTH,0.0f)};

        int Triangles = 0;
        for( int A=0; A<nItems; A++ )
        for( int B=A+1; B<nItems; B++ )
        for( int C=B+1; C<nItems; C++ )
        {
            for( int i=0; i<4; i++ )
            {
                float WA,WB,WC;
                CalcBaryXY( PanelRefCorners[i], PanelPos[A], PanelPos[B], PanelPos[C], WA, WB, WC );

                Q3D[i*3+0] += CamPos[A].X*WA + CamPos[B].X*WB + CamPos[C].X*WC;
                Q3D[i*3+1] += CamPos[A].Y*WA + CamPos[B].Y*WB + CamPos[C].Y*WC;
                Q3D[i*3+2] += CamPos[A].Z*WA + CamPos[B].Z*WB + CamPos[C].Z*WC;
            }

            Triangles++;
            break;
        }

        for( int i=0; i<4; i++ )
        {
            Q3D[i*3+0] /= (float)Triangles;
            Q3D[i*3+1] /= (float)Triangles;
            Q3D[i*3+2] /= (float)Triangles;
        }

        return true;
    }


    bool SolvePanelPoseFromItemsNew( EnvironmentInfo& EI, vector<Vec2> ItemImagePos, vector<Vec3> ItemPanelPos, float* Q3D )
    {
        // Choose the 4 items such that
        // A,B,C create a large area triangle
        // D has the largest negative WB barycentric weight
        int BestI[4] = {-1,-1,-1,-1};
        //{
        //    float BestArea = 0;
        //    for( int A=0; A<ItemImagePos.size(); A++ )
        //        for( int B=0; B<ItemImagePos.size(); B++ )
        //            for( int C=0; C<ItemImagePos.size(); C++ )
        //                for( int D=0; D<ItemImagePos.size(); D++ )
        //                if( (A!=B) && (A!=C) && (A!=D) &&
        //                    (B!=C) && (B!=D) && (C!=D) )
        //                {
        //                    float WA,WB,WC;
        //                    CalcBary(ItemImagePos[D],ItemImagePos[A],ItemImagePos[B],ItemImagePos[C],WA,WB,WC);
        //                    if( (WB < 0) && (WA>=0) && (WA<=1) && (WC>=0) && (WC<=1) )
        //                    {
        //                        float a = ItemImagePos[A].Dist( ItemImagePos[B] );
        //                        float b = ItemImagePos[B].Dist( ItemImagePos[C] );
        //                        float c = ItemImagePos[C].Dist( ItemImagePos[A] );
        //                        float p = (a+b+c)*0.5f;
        //                        float Area = sqrtf(p*(p-a)*(p-b)*(p-c));
        //                        if( Area > BestArea )
        //                        {
        //                            BestArea = Area;
        //                            BestI[0] = A;
        //                            BestI[1] = B;
        //                            BestI[2] = C;
        //                            BestI[3] = D;
        //                        }
        //                    }
        //                }
        //}

        log("total itemimagePos options: %d\n",ItemImagePos.size());
        {
            float BestArea = 0;
            float BestD = 0;
            for( int A=0; A<ItemImagePos.size(); A++ )
                for( int B=0; B<ItemImagePos.size(); B++ )
                    for( int C=0; C<ItemImagePos.size(); C++ )
                        for( int D=0; D<ItemImagePos.size(); D++ )
                        if( (A!=B) && (A!=C) && (A!=D) && (B!=C) && (B!=D) && (C!=D) )
                        {
                            {
                                float a = ItemImagePos[A].Dist( ItemImagePos[B] );
                                float b = ItemImagePos[B].Dist( ItemImagePos[C] );
                                float c = ItemImagePos[C].Dist( ItemImagePos[A] );
                                float p = (a+b+c)*0.5f;
                                float Area = sqrtf(p*(p-a)*(p-b)*(p-c));
                                if( (Area+1) >= BestArea )
                                {
                                    float WA,WB,WC;
                                    CalcBary(ItemImagePos[D],ItemImagePos[A],ItemImagePos[B],ItemImagePos[C],WA,WB,WC);
                                    //if( (WB < 0) && (WA>=0) && (WA<=1) && (WC>=0) && (WC<=1) )

                                    if( (Area>BestArea) || ((fabs(Area-BestArea)<0.001) && (WB<BestD)) )
                                    {
                                        BestArea = Area;
                                        BestI[0] = A;
                                        BestI[1] = B;
                                        BestI[2] = C;
                                        BestI[3] = D;
                                        BestD = WB;

                                        log("Area: %f  WA:%f  WB:%f  WC:%f\n",Area,WA,WB,WC);    
                                    }
                                }
                            }
                        }
        }
        for( int i=0; i<4; i++ )
            if( BestI[i]==-1 )
            {
                memset(Q3D,0,sizeof(float)*12);
                return false;
            }

        for( int i=0; i<4; i++ )
        {
            DrawLine(PixDebug,W,H,1,ItemImagePos[BestI[i]].X,ItemImagePos[BestI[i]].Y,ItemImagePos[BestI[(i+1)%4]].X,ItemImagePos[BestI[(i+1)%4]].Y);
        }


        // PREPARE FOR 4POINT POSE ESTIMATE
        vector<Vec3> RayDir;
        vector<Vec3> PanelPos;
        vector<Vec3> CamPos;
        for( int i=0; i<4; i++ )
        {
            Vec3 RDir(0,0,1);
            ImageToRay(EI,ItemImagePos[BestI[i]].X,ItemImagePos[BestI[i]].Y,RDir.X,RDir.Y);
            RayDir.push_back(RDir);
            PanelPos.push_back( ItemPanelPos[BestI[i]]);
            CamPos.push_back(Vec3(0,0,0));
        }

        Solve4PointCoplanarPose(EI,RayDir,PanelPos,CamPos);

        for( int i=0; i<4; i++ )
            log("CamPos %d   %f,%f,%f\n",i,CamPos[i].X,CamPos[i].Y,CamPos[i].Z);

        // CamPos contains best estimates of item positions
        memset(Q3D,0,sizeof(float)*12);
        Vec2 PanelRefCorners[4] = {Vec2(0,0), Vec2(0.0f,PANELHEIGHT), Vec2(PANELWIDTH,PANELHEIGHT), Vec2(PANELWIDTH,0.0f)};

        for( int i=0; i<4; i++ )
        {
            float WA,WB,WC;
            CalcBaryXY( PanelRefCorners[i], PanelPos[0], PanelPos[1], PanelPos[2], WA, WB, WC );

            Q3D[i*3+0] = CamPos[0].X*WA + CamPos[1].X*WB + CamPos[2].X*WC;
            Q3D[i*3+1] = CamPos[0].Y*WA + CamPos[1].Y*WB + CamPos[2].Y*WC;
            Q3D[i*3+2] = CamPos[0].Z*WA + CamPos[1].Z*WB + CamPos[2].Z*WC;
        }
        return true;
    }




    void FindPanel()
    {
        logtime LT("FindPanel");

        PanelFound = false;
        if( EnvInfo.Directory == ENVIRON_LAB3 )
        {
            FindPanelLab3();
            return;
        }
        else
        if( EnvInfo.Directory == ENVIRON_SIM )
        {
            PixIThresh = ClonePix(PixI,W,H);
            ClearOutsideRange(0,PixIThresh,W,H,32.0f/255.0f,250.0f/255.0f);
            ClearInsideRange(1,PixIThresh,W,H,32.0f/255.0f,250.0f/255.0f);
            for( int i=0; i<5; i++ ) ShrinkPix(PixIThresh,W,H,0.5f);
            for( int i=0; i<5; i++ ) GrowPix(PixIThresh,W,H,0.5f);
        }
        else
        if( EnvInfo.Directory == ENVIRON_LAB2 )
        {
            PixIThresh = ClonePix(PixI,W,H);

            float Threshold = SolveThresholdFromPanelLights();
            AvgFilter(PixI,PixIThresh,W,H);
            float Low = Threshold * 0.30f;
            //float High = 150.0f / 255.0f;
            float High = Threshold + 0.10f;
            ClearOutsideRange(0,PixIThresh,W,H,Low,High);
            ClearInsideRange(1,PixIThresh,W,H,Low,High);

            //AvgFilter(PixI,PixIThresh,W,H);
            //float Low = 15.0f / 255.0f;
            //float High = 150.0f / 255.0f;
            //ClearOutsideRange(0,PixIThresh,W,H,Low,High);
            //ClearInsideRange(1,PixIThresh,W,H,Low,High);
        }
        else
        {
            // If image is too dark...shift up and recalc histogram
            ////double AvgBright = 0;
            ////for( int i=0; i<W*H; i++ ) AvgBright += PixI[i];
            ////AvgBright /= W*H;
            ////if( AvgBright < 0.20 )
            ////{
            ////    float Offset = 0.7f - AvgBright;
            ////    for( int i=0; i<W*H; i++ ) 
            ////    {
            ////        float V = PixI[i] + Offset;
            ////        if( V>1 ) V=1;
            ////        if( V<0 ) V=0;
            ////        PixI[i] = V;
            ////    }

            ////    PixIHistogram = ComputeHistogram(256,PixI,W,H);
            ////}

            float Threshold = FindOtsuThreshold(PixIHistogram,256);
            log("Threshold %f\n",Threshold);

            float Th = Threshold;
            PixIThresh = ClonePix(PixI,W,H);
        
            // Try to mask out Robonaut glove
            if( Environment == ENVIRON_ISS )
            {
                if( AvgBright < 0.20 )
                {
                    //for( int i=0; i<W*H; i++ ) if( PixIThresh[i] > (0.156) ) PixIThresh[i] = 0;
                    for( int i=0; i<W*H; i++ ) 
                        if( (PixIThresh[i] < 0.047) || (PixIThresh[i] > 0.156) ) PixIThresh[i] = 0;
                        else PixIThresh[i] = 1;

                    for( int i=0; i<10; i++ ) ShrinkPix(PixIThresh,W,H,0.5f);
                }
                else
                {
                    for( int i=0; i<W*H; i++ ) if( PixIThresh[i] > (0.96) ) PixIThresh[i] = 0;
                    ClearInsideRange(0,PixIThresh,W,H,0,Th);
                    ClearOutsideRange(1,PixIThresh,W,H,0,Th);
                    for( int i=0; i<10; i++ ) ShrinkPix(PixIThresh,W,H,0.5f);
                }
            }
            else
            {
                ClearInsideRange(0,PixIThresh,W,H,0,Th);
                ClearOutsideRange(1,PixIThresh,W,H,0,Th);
                for( int i=0; i<10; i++ ) ShrinkPix(PixIThresh,W,H,0.5f);
            }

        }

       
        SaveBMP(EnvInfo.Name,"thresh.bmp",PixIThresh,W,H);

        float* PixFF = ClonePix(PixI,W,H);
        nBlobs = FloodFill( PixIThresh, W, H, PixFF,BlobIDs );
        SaveBMP(EnvInfo.Name,"ff.bmp",PixFF,W,H);

        int BestPowerCoverScore = 0;
        for( int i=0; i<nBlobs; i++ )
        {
            Quad Q = FindRANSACQuad(i);

            // Does quad look like a panel?
            if( IsQuadShapedLikePanel(Q) && DoesQuadLookLikePanel(Q))
            {
                if( EnvInfo.Directory != ENVIRON_SIM )
                    Q.Expand(-8,5,5,5);

                DrawLine(PixDebug,W,H,1.0f,Q.Pts[0].X,Q.Pts[0].Y,Q.Pts[1].X,Q.Pts[1].Y);
                DrawLine(PixDebug,W,H,1.0f,Q.Pts[1].X,Q.Pts[1].Y,Q.Pts[2].X,Q.Pts[2].Y);
                DrawLine(PixDebug,W,H,1.0f,Q.Pts[2].X,Q.Pts[2].Y,Q.Pts[3].X,Q.Pts[3].Y);
                DrawLine(PixDebug,W,H,1.0f,Q.Pts[3].X,Q.Pts[3].Y,Q.Pts[0].X,Q.Pts[0].Y);
                //DrawLine(PixDebug,W,H,1.0f,Q.Pts[0].X,Q.Pts[0].Y,Q.Pts[2].X,Q.Pts[2].Y);
                //DrawLine(PixDebug,W,H,1.0f,Q.Pts[1].X,Q.Pts[1].Y,Q.Pts[3].X,Q.Pts[3].Y);
                //DrawLine(PixDebug,W,H,1.0f,(Q.Pts[0].X+Q.Pts[1].X)/2,(Q.Pts[0].Y+Q.Pts[1].Y)/2,(Q.Pts[2].X+Q.Pts[3].X)/2,(Q.Pts[2].Y+Q.Pts[3].Y)/2);

                float Q3D[4*3];
                Solve3DQuad(EnvInfo,Q,Q3D);


                for( int j=0; j<22; j++ )
                {
                    Vec2 Center = CalcPerspectiveAwareBilinearPos(EnvInfo,Q3D,Items[j].RefPos.X,Items[j].RefPos.Y,Items[j].ZOffset);
                    DrawCircle(PixDebug,W,H,1,Center.X,Center.Y,8);
                }

                int PCScore = CalcColorMatchScore(EnvInfo,Q3D);

                if( (PCScore>0) && (PCScore > BestPowerCoverScore) )
                {
                    BestPowerCoverScore = PCScore;

                    PanelQuad = Q;
                    memcpy(Panel3D,Q3D,sizeof(float)*12);
                    PanelFound = true;

                    Vec2 Corners[4] = {Vec2(0,0), Vec2(0,1), Vec2(1,1), Vec2(1,0)};
                    for( int j=0; j<4; j++ )
                    {
                        Vec2 Center = CalcPerspectiveAwareBilinearPos(EnvInfo,Q3D,Corners[j].X,Corners[j].Y,0);
                        DrawCircle(PixDebug,W,H,1,Center.X,Center.Y,8);
                        DrawCircle(PixDebug,W,H,1,Center.X,Center.Y,9);
                        DrawCircle(PixDebug,W,H,1,Center.X,Center.Y,10);
                        DrawCircle(PixDebug,W,H,1,Center.X,Center.Y,11);
                        DrawLine(PixDebug,W,H,1,Center.X-8,Center.Y,Center.X+8,Center.Y);
                        DrawLine(PixDebug,W,H,1,Center.X,Center.Y-8,Center.X,Center.Y+8);
                    }
                }

            }
        }

        if( PanelFound )
        {
            DrawLine(PixDebug,W,H,1.0f,PanelQuad.Pts[0].X,PanelQuad.Pts[0].Y,PanelQuad.Pts[2].X,PanelQuad.Pts[2].Y);
            DrawLine(PixDebug,W,H,1.0f,PanelQuad.Pts[1].X,PanelQuad.Pts[1].Y,PanelQuad.Pts[3].X,PanelQuad.Pts[3].Y);
        }
    }

    Vec2 RunBlobEdgeFinder( int BlobID, int StartX, int StartY, int DirX, int DirY)
    {
        int BigStep = 16;
        int PrevX = StartX;
        int PrevY = StartY;
        int NewX = StartX;
        int NewY = StartY;
        int StepX = DirX * BigStep;
        int StepY = DirY * BigStep;
        bool SmallStepping = false;

        while( true )
        {
            PrevX = NewX;
            PrevY = NewY;
            NewX += StepX;
            NewY += StepY;
            if( (NewX>=0) && (NewX<W) && (NewY>=0) && (NewY<H) )
            {
                if( BlobIDs[ NewX + NewY*W ] == BlobID )
                {
                    if( SmallStepping )
                        return Vec2(NewX,NewY);
                    else
                    {
                        SmallStepping = true;
                        StepX = DirX;
                        StepY = DirY;
                        NewX = PrevX;
                        NewY = PrevY;
                    }
                }
            }
            else
                return Vec2(0,0);
        }
    }

    bool FindRANSACLine( vector<Vec2>& Pts, Vec2 ExpectedDir, Vec2& FinalPt1, Vec2& FinalPt2 )
    {
        int nPts = (int)Pts.size();

        log( "FindRANSACLine nPts: %d\n",nPts);
        if( nPts<2 )
        {
            FinalPt1.X = FinalPt1.Y = 0;
            FinalPt2.X = FinalPt2.Y = 0;
            return false;
        }

        Vec2* pPts = &Pts[0];

        float InlierThreshold = 4;
        float BestInliers = 0;
        float BestErr = 1000000.0f;
        Vec2 BestPt1;
        Vec2 BestPt2;

        for( int Trial=0; Trial<500; Trial++ )
        {
            int A = myrand() % nPts;
            int B = myrand() % nPts;

            Vec2 Pt = pPts[A];
            Vec2 Dir = Pt.DeltaTo(pPts[B]);
            float L = Dir.Len2();
            if( L < 64*64 ) continue;
            L = sqrtf(L);
            Dir = Dir.Mult( 1.0f / L );

            float ExpectedDirDot = Dir.Dot(ExpectedDir);
            if( fabs(ExpectedDirDot) < 0.7f )
                continue;

            Vec2 Perp = Dir.Perp();

            int nInliers = 0;
            float Err = 0;
            for( int i=0; i<nPts; i++ )
            {
                float DX = pPts[i].X - Pt.X;
                float DY = pPts[i].Y - Pt.Y;
                float PerpDist = DX*Perp.X + DY*Perp.Y;
                if( fabs(PerpDist) <= InlierThreshold )
                {
                    nInliers++;
                    Err += PerpDist * PerpDist;
                }
            }

            if( (nInliers>BestInliers) || ((nInliers==BestInliers) && (Err<BestErr)) )
            {
                BestErr = Err;
                BestInliers = nInliers;
                BestPt1 = pPts[A];
                BestPt2 = pPts[B];
            }
        }

        FinalPt1 = BestPt1;
        FinalPt2 = BestPt2;
        return true;
    }


    Quad FindRANSACQuad( int BlobID )
    {
        logtime LT("FindRANSACQuad");

        // Find edges
        vector<Vec2> TPts;
        vector<Vec2> BPts;
        vector<Vec2> LPts;
        vector<Vec2> RPts;
        Vec2 Zero(0,0);

        for( int X=0; X<W; X+=8 )
        {
            Vec2 TP = RunBlobEdgeFinder( BlobID, X, 0, 0, 1 );
            Vec2 BP = RunBlobEdgeFinder( BlobID, X, H-1, 0, -1 );
            if( !TP.Same(Zero) && (TP.Y>8) ) TPts.push_back(TP);
            if( !BP.Same(Zero) && (BP.Y<(H-8)) ) BPts.push_back(BP);
        }

        for( int Y=0; Y<H; Y+=8 )
        {
            Vec2 LP = RunBlobEdgeFinder( BlobID, 0, Y, 1, 0 );
            Vec2 RP = RunBlobEdgeFinder( BlobID, W-1, Y, -1, 0 );
            if( !LP.Same(Zero) && (LP.X>8) ) LPts.push_back(LP);
            if( !RP.Same(Zero) && (RP.X<(W-8))) RPts.push_back(RP);
        }

        Vec2 LinePts1[4]; // T,B,L,R
        Vec2 LinePts2[4]; // T,B,L,R

        bool RANSACSuccessful = true;
        RANSACSuccessful &= FindRANSACLine( TPts, Vec2(1,0), LinePts1[0], LinePts2[0] );
        RANSACSuccessful &= FindRANSACLine( BPts, Vec2(1,0), LinePts1[1], LinePts2[1] );
        RANSACSuccessful &= FindRANSACLine( LPts, Vec2(0,1), LinePts1[2], LinePts2[2] );
        RANSACSuccessful &= FindRANSACLine( RPts, Vec2(0,1), LinePts1[3], LinePts2[3] );

        if( !RANSACSuccessful )
        {
            log("FindRANSACLine FAILED!\n");
            Quad Q;
            memset(&Q,0,sizeof(Q));
            return Q;
        }

        for( int i=0; i<4; i++ )
        {
            DrawCircle(PixDebug,W,H,1.0f,LinePts1[i].X,LinePts1[i].Y,6);
            DrawCircle(PixDebug,W,H,1.0f,LinePts2[i].X,LinePts2[i].Y,6);

            log("LPS [%d] (%f,%f) (%f,%f)\n",i,LinePts1[i].X,LinePts1[i].Y,LinePts2[i].X,LinePts2[i].Y);
        }

        Quad Q;

        Q.Pts[0] = Find2DLineIntersection( LinePts1[0], LinePts2[0], LinePts1[2], LinePts2[2] ); // TL
        Q.Pts[1] = Find2DLineIntersection( LinePts1[1], LinePts2[1], LinePts1[2], LinePts2[2] ); // BL
        Q.Pts[2] = Find2DLineIntersection( LinePts1[1], LinePts2[1], LinePts1[3], LinePts2[3] ); // BR
        Q.Pts[3] = Find2DLineIntersection( LinePts1[0], LinePts2[0], LinePts1[3], LinePts2[3] ); // TR

        for( int i=0; i<4; i++ )
        {
            log("QPS [%d] (%f,%f)\n",i,Q.Pts[i].X,Q.Pts[i].Y);
        }

        return Q;
    }

    void TransferQ3D( Eye& E )
    {
        for( int i=0; i<4; i++ )
        {
            E.Panel3D[i*3+0] = Panel3D[i*3+0] + EnvInfo.CamPosTransfer[0];
            E.Panel3D[i*3+1] = Panel3D[i*3+1] + EnvInfo.CamPosTransfer[1];
            E.Panel3D[i*3+2] = Panel3D[i*3+2] + EnvInfo.CamPosTransfer[2];
        }
        E.PanelFound = true;
    }

    void FinalizePositions( void )
    {
        log("FinalizePositions\n");
        if( !PanelFound )
        {
            for( int i=0; i<22; i++ )
            {
                Items[i].Pos[Side].X = W/2;
                Items[i].Pos[Side].Y = H/2;
            }
            return;
        }

        for( int i=0; i<22; i++ )
        {
            Vec2 CourseCenter = CalcPerspectiveAwareBilinearPos(EnvInfo,Panel3D,Items[i].RefPos.X,Items[i].RefPos.Y,Items[i].ZOffset);
            Items[i].Pos[Side].X = CourseCenter.X;
            Items[i].Pos[Side].Y = CourseCenter.Y;
        }


        // Check cover position
        int PwrCapScore[2] = {0};
        {

            bool CoverIsUp = false;

            if( Environment == ENVIRON_LAB3 )
            {
                Vec2 CenterDown = CalcPerspectiveAwareBilinearPos(EnvInfo,Panel3D,Items[ITEM_PANEL_POWER_COVER].RefPos.X,Items[ITEM_PANEL_POWER_COVER].RefPos.Y,Items[ITEM_PANEL_POWER_COVER].ZOffset);
                Vec2 CenterUp = CalcPerspectiveAwareBilinearPos(EnvInfo,Panel3D,Items[ITEM_PANEL_POWER_COVER_UP].RefPos.X,Items[ITEM_PANEL_POWER_COVER_UP].RefPos.Y+0.004f,Items[ITEM_PANEL_POWER_COVER_UP].ZOffset*0.75f);
                int Radius = 8;
                int Step = 1;
                int CX[2] = {(int)CenterDown.X,(int)CenterUp.X};
                int CY[2] = {(int)CenterDown.Y,(int)CenterUp.Y};
                DrawLine(PixDebug,W,H,1,CX[0]-Radius,CY[0]-Radius,CX[0]+Radius,CY[0]+Radius);
                DrawLine(PixDebug,W,H,1,CX[0]-Radius,CY[0]+Radius,CX[0]+Radius,CY[0]-Radius);
                DrawLine(PixDebug,W,H,1,CX[1]-Radius,CY[1]-Radius,CX[1]+Radius,CY[1]+Radius);
                DrawLine(PixDebug,W,H,1,CX[1]-Radius,CY[1]+Radius,CX[1]+Radius,CY[1]-Radius);

                for( int i=0; i<2; i++ )
                {
                    double Avg=0;
                    double AvgC=0;
                    for( int Y=CY[i]-Radius; Y<=CY[i]+Radius; Y+=Step )
                        for( int X=CX[i]-Radius; X<=CX[i]+Radius; X+=Step )
                            if( (Y>=0) && (Y<H) && (X>=0) && (X<W) )
                            {
                                Avg += PixI[X+Y*W];
                                AvgC++;
                            }
                    if( AvgC>0 ) Avg /= AvgC;
                    PwrCapScore[i] = (int)(Avg*100);
                }

                log("PwrCapScore Scores: DOWN:%d UP:%d\n",PwrCapScore[0],PwrCapScore[1]);
                CoverIsUp = PwrCapScore[1] < PwrCapScore[0];
            }
            else
            {
                Vec2 CenterDown = CalcPerspectiveAwareBilinearPos(EnvInfo,Panel3D,Items[ITEM_PANEL_POWER_COVER].RefPos.X,Items[ITEM_PANEL_POWER_COVER].RefPos.Y,Items[ITEM_PANEL_POWER_COVER].ZOffset);
                Vec2 CenterUp = CalcPerspectiveAwareBilinearPos(EnvInfo,Panel3D,Items[ITEM_PANEL_POWER_COVER_UP].RefPos.X,Items[ITEM_PANEL_POWER_COVER_UP].RefPos.Y,Items[ITEM_PANEL_POWER_COVER_UP].ZOffset);
                int Radius = 8;
                int Step = 1;
                int CX[2] = {(int)CenterDown.X,(int)CenterUp.X};
                int CY[2] = {(int)CenterDown.Y,(int)CenterUp.Y};
                DrawLine(PixDebug,W,H,1,CX[0]-Radius,CY[0]-Radius,CX[0]+Radius,CY[0]+Radius);
                DrawLine(PixDebug,W,H,1,CX[0]-Radius,CY[0]+Radius,CX[0]+Radius,CY[0]-Radius);
                DrawLine(PixDebug,W,H,1,CX[1]-Radius,CY[1]-Radius,CX[1]+Radius,CY[1]+Radius);
                DrawLine(PixDebug,W,H,1,CX[1]-Radius,CY[1]+Radius,CX[1]+Radius,CY[1]-Radius);

                for( int i=0; i<2; i++ )
                {
                    int YC=0;
                    int NC=0;
                    for( int Y=CY[i]-Radius; Y<=CY[i]+Radius; Y+=Step )
                        for( int X=CX[i]-Radius; X<=CX[i]+Radius; X+=Step )
                        {
                            if( (Y>=0) && (Y<H) && (X>=0) && (X<W) )
                            {
                                int O = X+Y*W;

                                float R = PixR[O];
                                float G = PixG[O];
                                float B = PixB[O];
                                float GB = (G+B)*0.5f;
                                float Ratio = 0;
                                if( GB>0 ) Ratio = R/GB;
                    
                                if( (R>G) && (R>B) && (Ratio > EnvInfo.PwrCapRedRatio) ) YC++;
                                else NC++;
                            }
                        }

                    if( YC+NC > 0 ) 
                        PwrCapScore[i] = (YC*100/(YC+NC));
                }

                log("PwrCapScore Scores: DOWN:%d UP:%d\n",PwrCapScore[0],PwrCapScore[1]);
                CoverIsUp = (PwrCapScore[0]<25) && (PwrCapScore[1] > PwrCapScore[0]);
            }


            if( CoverIsUp )
            {
                log("Cover Up\n");
                Vec2 Pos = CalcPerspectiveAwareBilinearPos(EnvInfo,Panel3D,Items[ITEM_PANEL_POWER_COVER_UP].RefPos.X,Items[ITEM_PANEL_POWER_COVER_UP].RefPos.Y,Items[ITEM_PANEL_POWER_COVER_UP].ZOffset);
                Items[ITEM_PANEL_POWER_COVER].Pos[Side] = Pos;
                Items[ITEM_PANEL_POWER_COVER].StateConfidence[Side][STATE_UP] = 1;
            }
            else
            {
                log("Cover Down\n");
                Vec2 Pos = CalcPerspectiveAwareBilinearPos(EnvInfo,Panel3D,Items[ITEM_PANEL_POWER_COVER].RefPos.X,Items[ITEM_PANEL_POWER_COVER].RefPos.Y,Items[ITEM_PANEL_POWER_COVER].ZOffset);
                Items[ITEM_PANEL_POWER_COVER].Pos[Side] = Pos;
                // Use default state
                //Items[ITEM_PANEL_POWER_COVER].StateConfidence[Side][STATE_UP] = 1;
            }
        }
        log("FinalizePositions completed\n");
    }

    void SolveStates( void )
    {
        log("SolveStates..started\n");
        if( !PanelFound ) return;
        for( int i=0; i<22; i++ )
        {
            // Check LEDs
            if( Items[i].StateMask[STATE_ON] == 1 )
            {
                float XSearchRadius = 12.0f / 518;
                float YSearchRadius = 9.0f / 518;
                if((i>=ITEM_A02_LED_NUM_PAD_A1) && (i<=ITEM_A02_LED_NUM_PAD_C3) )
                {        
                    XSearchRadius = 25.0f / 518;
                    YSearchRadius = 16.0f / 518;
                }

                Vec2 Center = Items[i].Pos[Side];
                Vec2 XRRef = CalcPerspectiveAwareBilinearPos(EnvInfo,Panel3D,Items[i].RefPos.X-XSearchRadius,Items[i].RefPos.Y,Items[i].ZOffset);
                Vec2 YRRef = CalcPerspectiveAwareBilinearPos(EnvInfo,Panel3D,Items[i].RefPos.X,Items[i].RefPos.Y-YSearchRadius,Items[i].ZOffset);
                int XRadius = (int)XRRef.Dist(Center);
                int YRadius = (int)YRRef.Dist(Center);
                if( XRadius > 100 ) XRadius = 100;
                if( XRadius < 8 ) XRadius = 8;
                if( YRadius > 100 ) YRadius = 100;
                if( YRadius < 8 ) YRadius = 8;

                DrawRect(PixDebug,W,H,1,Center.X-XRadius,Center.Y-YRadius,XRadius*2,YRadius*2);

                //DrawCircle(PixDebug,W,H,1.0f,Center.X,Center.Y,ImgRadius);

                int NC=0;
                int YC=0;

                for( int X=(int)(Center.X-XRadius); X<=(int)(Center.X+XRadius); X++ )
                    for( int Y=(int)(Center.Y-YRadius); Y<=(int)(Center.Y+YRadius); Y++ )
                        if( (X>=0) && (X<W) && (Y>=0) && (Y<H) )
                        {
                            float V = PixIllum[X+Y*W];
                            if( V > 0.5f ) YC++; else NC++;
                        }


                int Threshold = 10;
                if( Environment == ENVIRON_LAB ) Threshold = 5;
                if( Environment == ENVIRON_ISS ) Threshold = 12;

                if( YC >= Threshold ) Items[i].StateConfidence[Side][STATE_ON] = 1.0f;
                else         Items[i].StateConfidence[Side][STATE_OFF] = 1.0f;

                log("%2d - ImgRadius %d,%d  YC:%d  NC:%d  %f\n",i,XRadius,YRadius,YC,NC,Items[i].StateConfidence[Side][STATE_ON]);
            }
        }

        if( Environment == ENVIRON_LAB3 )
        {
            Items[ITEM_PANEL_POWER_SWITCH].StateConfidence[Side][STATE_UP] = 1;
            Items[ITEM_PANEL_POWER_COVER].StateConfidence[Side][STATE_UP] = 1;
            Items[ITEM_PANEL_POWER_LED].StateConfidence[Side][STATE_ON] = 1;
        }

        log("SolveStates..finished\n");
    }

    //void SolveStates( void )
    //{
    //    log("SolveStates\n");
    //    if( !PanelFound ) return;
    //    for( int i=0; i<22; i++ )
    //    {
    //        // Check LEDs
    //        if( Items[i].StateMask[STATE_ON] == 1 )
    //        {
    //            float SearchRadius = 0;
    //            if((i>=ITEM_A02_LED_NUM_PAD_A1) && (i<=ITEM_A02_LED_NUM_PAD_C3) ) SearchRadius = 20.0f / 518;
    //            else SearchRadius = 12.0f / 518;

    //            Vec2 Center = Items[i].Pos[Side];
    //            Vec2 RRef = CalcPerspectiveAwareBilinearPos(EnvInfo,Panel3D,Items[i].RefPos.X,Items[i].RefPos.Y-SearchRadius,Items[i].ZOffset);
    //            int ImgRadius = (int)RRef.Dist(Center);
    //            if( ImgRadius > 100 ) ImgRadius = 100;
    //            if( ImgRadius < 8 ) ImgRadius = 8;
    //            int ImgRadius2 = ImgRadius*ImgRadius;

    //            
    //            DrawCircle(PixDebug,W,H,1.0f,Center.X,Center.Y,ImgRadius);

    //            int NC=0;
    //            int YC=0;

    //            for( int DY=-ImgRadius; DY<=ImgRadius; DY++ )
    //                for( int DX=-ImgRadius; DX<=ImgRadius; DX++ )
    //                {
    //                    float D2 = DX*DX + DY*DY;
    //                    if( D2 <= ImgRadius2 )
    //                    {
    //                        int X = (int)(Center.X + DX);
    //                        int Y = (int)(Center.Y + DY);
    //                        if( (X>=0) && (X<W) && (Y>=0) && (Y<H) )
    //                        {
    //                            float V = PixIllum[X+Y*W];
    //                            if( V > 0.5f ) YC++; else NC++;
    //                        }
    //                    }
    //                }

    //            log("%2d - ImgRadius %d  YC:%d  NC:%d\n",i,ImgRadius,YC,NC);

    //            int Threshold = 10;
    //            if( Environment == ENVIRON_LAB ) Threshold = 5;
    //            if( Environment == ENVIRON_ISS ) Threshold = 12;

    //            if( YC >= Threshold ) Items[i].StateConfidence[Side][STATE_ON] = 1.0f;
    //            else         Items[i].StateConfidence[Side][STATE_OFF] = 1.0f;
    //        }
    //    }

    //    if( Environment == ENVIRON_LAB3 )
    //    {
    //        Items[ITEM_PANEL_POWER_SWITCH].StateConfidence[Side][STATE_UP] = 1;
    //        Items[ITEM_PANEL_POWER_COVER].StateConfidence[Side][STATE_UP] = 1;
    //        Items[ITEM_PANEL_POWER_LED].StateConfidence[Side][STATE_ON] = 1;
    //    }
    //}


    void Kill()
    {
        //if(PixR) free(PixR);
        //if(PixG) free(PixG);
        //if(PixB) free(PixB);
        //if(PixI) free(PixI);
        //if(PixIEdges) free(PixIEdges);
        //if(PixIThresh) free(PixIThresh);
        //if(BlobIDs) free(BlobIDs);
        //if(PixDebug) free(PixDebug);
        //if(PixIllum) free(PixIllum);
    }
};


//=========================================================================
//=========================================================================
//=========================================================================

vector<string> RobonautEye::recognizeObjects( vector<int>& leftEyeImage, vector<int>& rightEyeImage)
{
    StartClock();

    log("\n\nRobonautEye::recognizeObjects ENTER\n");

    InitItems();

    // Fill out some lame answers
    int Width = leftEyeImage[0];
    int Height = leftEyeImage[1];
    for( int i=0; i<sizeof(Items)/sizeof(ItemInfo); i++ )
    {
        Items[i].FinalState = STATE_HIDDEN;//Items[i].DefaultState;
        Items[i].Pos[0].X= Width/2;
        Items[i].Pos[0].Y= Height/2;
        Items[i].Pos[1].X = Width/2;
        Items[i].Pos[1].Y = Height/2;
    }
    
    Eye LeftEye;
    Eye RightEye;


  ////  {
  ////      logtime LT("LeftEYE");

  ////      log("\n\n");
  ////      LeftEye.Init( "Left", leftEyeImage, 0 );
  ////      LeftEye.FinalizePositions();
  ////      LeftEye.SolveStates();
  ////      SaveBMP(LeftEye.EnvInfo.Name,"debug.bmp",LeftEye.PixDebug,LeftEye.W,LeftEye.H);
  ////  }

  ////  {    
        ////logtime LT("RightEYE");

  ////      log("\n\n");
  ////      RightEye.Init( "Right", rightEyeImage, 1 );
  ////      RightEye.FinalizePositions();
  ////      RightEye.SolveStates();
  ////      SaveBMP(RightEye.EnvInfo.Name,"debug.bmp",RightEye.PixDebug,RightEye.W,RightEye.H);
  ////  }


    {
        log("\n\n");
        LeftEye.Init( "Left", leftEyeImage, 0 );
        log("\n\n");
        RightEye.Init( "Right", rightEyeImage, 1 );

        if( !LeftEye.PanelFound && RightEye.PanelFound )
        {
            RightEye.TransferQ3D(LeftEye);
            LeftEye.FinalizePositions();
            RightEye.FinalizePositions();
            LeftEye.SolveStates();
            RightEye.SolveStates();
            SaveBMP(LeftEye.EnvInfo.Name,"debug.bmp",LeftEye.PixDebug,LeftEye.W,LeftEye.H);
            SaveBMP(RightEye.EnvInfo.Name,"debug.bmp",RightEye.PixDebug,RightEye.W,RightEye.H);
        }
        else
        if( LeftEye.PanelFound && !RightEye.PanelFound )
        {
            LeftEye.TransferQ3D(RightEye);
            LeftEye.FinalizePositions();
            RightEye.FinalizePositions();
            LeftEye.SolveStates();
            RightEye.SolveStates();
            SaveBMP(LeftEye.EnvInfo.Name,"debug.bmp",LeftEye.PixDebug,LeftEye.W,LeftEye.H);
            SaveBMP(RightEye.EnvInfo.Name,"debug.bmp",RightEye.PixDebug,RightEye.W,RightEye.H);
        }
        else
        if( LeftEye.PanelFound && RightEye.PanelFound )
        {
            LeftEye.FinalizePositions();
            RightEye.FinalizePositions();
            LeftEye.SolveStates();
            RightEye.SolveStates();
            SaveBMP(LeftEye.EnvInfo.Name,"debug.bmp",LeftEye.PixDebug,LeftEye.W,LeftEye.H);
            SaveBMP(RightEye.EnvInfo.Name,"debug.bmp",RightEye.PixDebug,RightEye.W,RightEye.H);
        }

    }

    // Resolve final LED states
    for( int i=0; i<22; i++ )
    {
        if( Items[i].StateMask[STATE_ON] == 1 )
        {
            log("[%2d] - ItemState %f %f\n",i,Items[i].StateConfidence[0][STATE_ON],Items[i].StateConfidence[1][STATE_ON]);

            if( (Items[i].StateConfidence[0][STATE_ON]>0) || (Items[i].StateConfidence[1][STATE_ON]>0) )
            {
                Items[i].FinalState = STATE_ON;
                Items[ITEM_PANEL_POWER_LED].StateConfidence[0][STATE_ON] = 1;
                Items[ITEM_PANEL_POWER_LED].StateConfidence[1][STATE_ON] = 1;
                Items[ITEM_PANEL_POWER_LED].FinalState = STATE_ON;
            }
            else
                Items[i].FinalState = STATE_OFF;
        }
    }

    if( Items[ITEM_PANEL_POWER_LED].FinalState == STATE_ON ) Items[ITEM_PANEL_POWER_SWITCH].FinalState = STATE_UP;
    else Items[ITEM_PANEL_POWER_SWITCH].FinalState = STATE_DOWN;

    if( Items[ITEM_A03_LED].FinalState == STATE_ON ) Items[ITEM_A03_TOGGLE].FinalState = STATE_UP;
    else Items[ITEM_A03_TOGGLE].FinalState = STATE_DOWN;

    if( Items[ITEM_A04_LED_TOP].FinalState == STATE_ON ) Items[ITEM_A04_TOGGLE].FinalState = STATE_UP;
    else
    if( Items[ITEM_A04_LED_BOTTOM].FinalState == STATE_ON ) Items[ITEM_A04_TOGGLE].FinalState = STATE_DOWN;
    else Items[ITEM_A04_TOGGLE].FinalState = STATE_CENTER;

    if( Items[ITEM_A05_LED].FinalState == STATE_ON ) Items[ITEM_A05_TOGGLE].FinalState = STATE_UP;
    else Items[ITEM_A05_TOGGLE].FinalState = STATE_DOWN;

    
    if( Items[ITEM_A01_ROCKER_LED_TOP].FinalState == STATE_ON ) Items[ITEM_A01_ROCKER_SWITCH].FinalState = STATE_UP;
    else
    if( Items[ITEM_A01_ROCKER_LED_BOTTOM].FinalState == STATE_ON ) Items[ITEM_A01_ROCKER_SWITCH].FinalState = STATE_DOWN;
    else Items[ITEM_A01_ROCKER_SWITCH].FinalState = STATE_CENTER;
    

    if( (Items[ITEM_PANEL_POWER_COVER].StateConfidence[0][STATE_UP]==1) ||
        (Items[ITEM_PANEL_POWER_COVER].StateConfidence[1][STATE_UP]==1) )
        Items[ITEM_PANEL_POWER_COVER].FinalState = STATE_UP;
    else
        Items[ITEM_PANEL_POWER_COVER].FinalState = STATE_DOWN;

    double Runtime = SecondsRunning();
    log("Seconds: %8.3f\n",Runtime);

    // Build response
    vector<string> Response;
    ////if( Runtime > 4.9 )//LeftEye.Environment != ENVIRON_LAB3 )
    ////{
    ////    for( int i=0; i<22; i++ )
    ////    {
    ////        char Buff[256];
    ////        sprintf(Buff,"%s,%d,%d,%d,%d","IGNORE",0,0,0,0);
    ////        Response.push_back( string(Buff) );

    ////        log("[%2d] %s\n",i,Buff);
    ////    }
    ////}
    ////else
    {
        for( int i=0; i<22; i++ )
        {
            char Buff[256];
            sprintf(Buff,"%s,%d,%d,%d,%d",StateNames[Items[i].FinalState],(int)Items[i].Pos[0].X,(int)Items[i].Pos[0].Y,(int)Items[i].Pos[1].X,(int)Items[i].Pos[1].Y);
            Response.push_back( string(Buff) );

            log("[%2d] %s\n",i,Buff);
        }
    }

    LeftEye.Kill();
    RightEye.Kill();

    log("RobonautEye::recognizeObjects EXIT\n");

    fflush(stdout);
    fflush(stderr);
    return Response;
}


//=========================================================================
//=========================================================================
//=========================================================================
//=========================================================================
//#ifndef FULL_SUBMISSION

#if !defined(FULL_SUBMISSION) && !defined(SUBMISSION)

#include "stdafx.h"
#include <iostream>
#include <stdio.h>

inline int readint(void)
{
    int I;
    fscanf(stdin,"%d",&I);
    return I;
}

int _tmain(int argc, char* argv[])
{
    log("_tmain\n");

    RobonautEye REye;

    int len = readint();

    vector<int> LeftEye(len);
    vector<int> RightEye(len);
    
    for( int i=0; i<len; i++ ) LeftEye[i] = readint();
    for( int i=0; i<len; i++ ) RightEye[i] = readint();

    vector<string> resp = REye.recognizeObjects( LeftEye, RightEye );

    for( int i=0; i<22; i++ )
    {
        fprintf(stdout,resp[i].c_str());
        fprintf(stdout, "\n");
    }
    fflush(stdout);

}

#endif
