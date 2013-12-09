#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include <assert.h>

#include <vector>
#include <list>
#include <map>
#include <set>
#include <queue>
#include <deque>
#include <stack>
#include <bitset>
#include <algorithm>
#include <functional>
#include <numeric>
#include <utility>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <complex>
#include <fstream>

using namespace std;

#define FOR(i, b, e)    for(int i = (b); i <= (e); i++)
#define FORL(i, b, e)    for(int i = (b); i < (e); i++)
#define FORD(i, e, b)    for(int i = (e); i >= (b); i--)
#define FOR0(i, e)        FORL(i, 0, e)

#define min(a, b)        (((a) < (b)) ? (a) : (b))
#define max(a, b)        (((a) > (b)) ? (a) : (b))
#define MINA(a, b)        do { if ((a) > (b)) (a) = (b); } while(0)
#define MAXA(a, b)        do { if ((a) < (b)) (a) = (b); } while(0)
#define MINA2(a, b, i, j)        do { if ((a) > (b)) { (a) = (b); (i) = (j); } } while(0)
#define MAXA2(a, b, i, j)        do { if ((a) < (b)) { (a) = (b); (i) = (j); } } while(0)

#define SWAP(a, b)        do { int _t = a; a = b; b = _t; } while(0)
#define SWAPT(a, b, t)    do { t _t = a; a = b; b = _t; } while(0)
#define SQR(a)            ((a) * (a))
#define MSET(a, b)        memset(a, b, sizeof(a))

#define INT                int
#define INT_CAP            0x3F3F3F3F

typedef long long int    LI;

typedef pair<int, int>    II;
typedef vector<int>       VI;
typedef vector<double>    VD;
typedef vector<string>    VS;
#define ALL(c)            c.begin(), c.end()
#define SZ(c)            (static_cast<int>(c.size()))
#define FORALL(it, c)     for(it = c.begin(); it != c.end(); ++it)
#define FORALLR(it, c)    for(it = c.rbegin(); it != c.rend(); ++it)
#define PB                push_back
#define MP                make_pair
#define P1                first
#define P2                second

FILE *debugf = stdout;
char *testcase = "";
char *testmode = "";
int tester, testdump, testparam, testcommand;
int timerstep, timerstepcap, localtimeout;
LI globalstart, localstart;

#define DEBUG(f, ...)    do { fprintf(debugf, "DEBUG:" f "\n", ##__VA_ARGS__); fflush(debugf); } while(0) 

#ifdef _DEBUG
#define DEBUG1(f, ...)    do { fprintf(debugf, "DEBUG:" f "\n", ##__VA_ARGS__); fflush(debugf); } while(0) 
#define ASSERT    assert
#else
#define DEBUG1(f, ...)      do {} while(0)
#define ASSERT
#endif

#ifndef LOCAL

#include <sys/time.h>
 
LI gettimems() 
{
    timeval tv; 
    gettimeofday(&tv, 0); 
    return tv.tv_sec * 1000 + tv.tv_usec / 1000;
}
#else

#define sprintf sprintf_s

#include <time.h>

LI gettimems() 
{
    return (LI) clock() * 1000 / CLOCKS_PER_SEC;
}

#endif

int randn0(int n)
{
/*    DEBUG("RANDDDD!");
    if (n <= 1)
        return 0; */
//    return 1;  
    return rand() % n;                 /* fast, but incorrect */
}

/*---*/

#define COORD_DOUBLE
#define GETDBL(x)        x = COORD_CAP

#define PI        3.141592653589793L
#define RAD(x)    (PI * (x) / 180)
#define EPS        1e-4 // must be set
#define ROUND(x) ((x)>=0?(int)((x)+0.5):(int)((x)-0.5))
#define DBL_CAP    (DBL_MAX / 2)
#define INT_CAP    0x3F3F3F3F

#ifdef COORD_DOUBLE
#define COORD        double
#define COORDX        double
#define CABS        fabs
#define CABS0(x)    (fabs(x) < EPS)
#define CLE(x, y)    ((x) < (y) + EPS)
#define CL(x, y)    ((x) < (y) - EPS)
#define CHALF(x)    (0.5 * (x))
#define COORD_CAP    (DBL_MAX / 2)
#define GETCOORD    GETDBL
#else
#define COORD        int
#define COORDX        long long int
#define CABS(x)        ((x>0)?(x):-(x))
#define CABS0(x)    ((x) == 0)
#define CLE(x, y)    ((x) <= (y))
#define CL(x, y)    ((x) < (y))
#define CHALF(x)    ((x) >> 1)
#define COORD_CAP    INT_CAP
#define GETCOORD    GETI
#endif

#define SQRC(x)        SQR((COORDX) x)

COORD lenc2(COORD x, COORD y)
{
    return SQRC(x) + SQRC(y);
}

double lenc(COORD x, COORD y)
{
    return sqrt((double) lenc2(x, y));
}

double distc(COORD x1, COORD y1, COORD x2, COORD y2)
{
    return lenc(x1 - x2, y1 - y2);
}

COORDX distc2(COORD x1, COORD y1, COORD x2, COORD y2)
{
    return lenc2(x1 - x2, y1 - y2);
}

double anglec(double a, double b, double c) /* angle opposite of side a */
{
    return acos((SQR(b) + SQR(c) - SQR(a)) / (2.0 * b * c));
}

#ifdef __cplusplus 
struct pv {
    COORD x;
    COORD y;
    pv(COORD ax = COORD_CAP, COORD ay = COORD_CAP): x(ax), y(ay) {}
};
struct pvp {
    pv p[2];
};
#else
typedef struct pv {
    COORD x;
    COORD y;
} pv;
typedef struct pvp { // point pair (just for convenience)
    pv p[2];
} pvp;
#endif


bool isvalid(pv p)
{
    return p.x != COORD_CAP;
}

pv invalidpv()
{
    pv p;
#ifndef __cplusplus
    p.x = COORD_CAP;
    p.y = COORD_CAP;
#endif
    return p;
}

double len(pv p)
{
    return lenc(p.x, p.y);
}

double len2(pv p)
{
    return lenc2(p.x, p.y);
}


double dist(pv a, pv b)
{
    return distc(a.x, a.y, b.x, b.y);
}

COORDX dist2(pv a, pv b)
{
    return distc2(a.x, a.y, b.x, b.y);
}

pv add(pv a, pv b)
{
    pv p;
    p.x = a.x + b.x;
    p.y = a.y + b.y;
    return p;
}

pv sub(pv a, pv b)
{
    pv p;
    p.x = a.x - b.x;
    p.y = a.y - b.y;
    return p;
}

COORDX dot(pv a, pv b)
{
    return (COORDX) a.x * b.x + (COORDX) a.y * b.y;
}

COORDX cross(pv a, pv b)
{
    return (COORDX) a.x * b.y - (COORDX) a.y * b.x;
}

pv rotl(pv a)
{
    pv p;
    p.x = -a.y;
    p.y = a.x;
    return p;
}

pv rotr(pv a)
{
    pv p;
    p.x = a.y;
    p.y = -a.x;
    return p;
}

pv inv(pv a)
{
    pv p;
    p.x = -a.x;
    p.y = -a.y;
    return p;
}

bool eq(pv a, pv b)
{
    return CABS0(a.x-b.x) && CABS0(a.y-b.y); // fast, but not perfectly EPS, as it checks coords separately, not distance
}

int cmp(const void *arg1, const void *arg2)
{
    pv *p1, *p2;
    p1 = (pv *) arg1;
    p2 = (pv *) arg2;
    if (CABS0(p1->x - p2->x)) {
        if (CABS0(p1->y - p2->y))
            return 0;
        else if (p1->y > p2->y)
            return 1;
        return -1;
    }
    else if (p1->x > p2->x)
        return 1;
    return -1;
}

int side(pv p, pv a, pv b) /* 0: on, 1:right, -1:left */
{
    COORDX s;
    s = cross(sub(a, p), sub(b, p));
    if (CABS0(s)) 
        return 0;
    if (s < 0)
        return -1;
    return 1;
}

COORDX trap2(pv a, pv b) /* (doubled and can be negative) for area calc */
{
    return (COORDX) (a.x - b.x) * (a.y + b.y);
}

pv getp()
{
    pv p;
    GETCOORD(p.x);
    GETCOORD(p.y);
    return p;
}

#ifdef __cplusplus

pv operator+(pv a, pv b)
{
    return add(a, b);
}

pv operator-(pv a, pv b)
{
    return add(a, b);
}

pv operator==(pv a, pv b)
{
    return eq(a, b);
}
#endif

#ifdef COORD_DOUBLE
pv mul(double m, pv a)
{
    pv p;
    p.x = (COORD) (m * a.x);
    p.y = (COORD) (m * a.y);
    return p;
}

#ifdef __cplusplus
pv operator*(pv p, double m)
{
    return mul(m, p);
}

pv operator*(double m, pv p)
{
    return mul(m, p);
}
#endif

pv norm(pv p)
{
    double d;
    d = lenc(p.x, p.y);
    if (CABS0(d))
        return p;        /* error */
    return mul(1 / d, p);
}

pv rot(pv p, pv b, double al)
{
    pv pd, pf;
    pd = sub(p, b);
    pf.x = (COORD) (cos(al) * pd.x - sin(al) * pd.y);
    pf.y = (COORD) (sin(al) * pd.x + cos(al) * pd.y);
    return add(b, pf);
}

double angle(pv a, pv b, pv c)    /* alpha: angle at a */
{
    pv p1, p3;
    p1 = norm(sub(b, a));
    p3 = sub(c, a);
    return atan2((double) cross(p1, p3), (double) dot(p3, p1));
}

pv mirror(pv p, pv a, pv b)
{
    if (eq(p, a))
        return p;
    return rot(p, a, -2.0 * angle(p, a, b));
}

pv segp(pv a, pv b, double m)
{
    return add(a, mul(m, sub(b, a)));
}

pv segpd2(pv a, pv b, double d2)
{
    return segp(a, b, d2 / dist2(a, b));
}

pv segpd(pv a, pv b, double d)
{
    return segp(a, b, d / dist(a, b));
}

pv closest(pv p, pv a, pv b)
{
    return segpd2(a, b, (double) dot(sub(p, a), sub(b, a)));
}
#else /* not COORD_DOUBLE */

double angle(pv a, pv b, pv c)    /* alpha: angle at a */
{
    double dab, dac, dbc;
    dab = dist(a, b);
    dac = dist(a, c);
    dbc = dist(b, c);
    return anglec(dbc, dac, dab);
}

#endif

void normiv(COORD *a, COORD *b)
{
    if (*a > *b)
        SWAPT(*a, *b, COORD);
}

bool oniv(COORD a, COORD b1, COORD b2)
{
    return CLE(b1, a) && CLE(a, b2);
}

bool iniv(COORD a, COORD b1, COORD b2)
{
    return CL(b1, a) && CL(a, b2);
}

bool ivoniv(COORD a1, COORD a2, COORD b1, COORD b2)
{
    return CLE(b1, a1) && CLE(a2, b2);
}

bool iviniv(COORD a1, COORD a2, COORD b1, COORD b2)
{
    return CL(b1, a1) && CL(a2, b2);
}

bool ivtouchiv(COORD a1, COORD a2, COORD b1, COORD b2)
{
    if (a1 < b1)
        return CABS0(a2 - b1);
    else
        return CABS0(b2 - a1);
}

bool ivoverlapiv(COORD a1, COORD a2, COORD b1, COORD b2)
{
    return CL(max(a1, b1), min(a2, b2));
}

bool ivtouchoverlapiv(COORD a1, COORD a2, COORD b1, COORD b2)
{
    return CLE(max(a1, b1), min(a2, b2));
}

bool online(pv p, pv a, pv b)
{
    return side(p, a, b) == 0;
}

bool onseg(pv p, pv a, pv b)
{
    if (side(p, a, b) != 0)
        return false;
    if (CABS0(a.x - b.x)) {
        SWAPT(a.x, a.y, COORD);
        SWAPT(b.x, b.y, COORD);
        SWAPT(p.x, p.y, COORD);
    }
    normiv(&a.x, &b.x);
    return oniv(p.x, a.x, b.x);
}

bool inseg(pv p, pv a, pv b)
{
    if (side(p, a, b) != 0)
        return false;
    if (CABS0(a.x - b.x)) {
        SWAPT(a.x, a.y, COORD);
        SWAPT(b.x, b.y, COORD);
        SWAPT(p.x, p.y, COORD);
    }
    normiv(&a.x, &b.x);
    return iniv(p.x, a.x, b.x);
}

pv closesttoseg(pv p, pv a, pv b) 
{
    COORDX d2;
    d2 = dot(sub(p, a), sub(b, a));
    if (CABS0(d2))
        return a;
    else if (CABS0(dot(sub(p, b), sub(a, b))))
        return b;
    pv pr = segpd2(a, b, d2);
    if (onseg(pr, a, b))
        return pr;
    if (dist2(p, a) < dist2(p, b))
        return a;
    return b;
}

pv intersect(pv a1, pv a2, pv b1, pv b2)
{
    pv p; 
    COORDX d;
    p = sub(b2, b1);
    d = cross(sub(a2, a1), p);
    if (CABS0(d))
        return invalidpv();
    return segp(a1, a2, (double) cross(sub(b1, a1), p) / d);
}

bool isparallel(pv a1, pv a2, pv b1, pv b2)
{
    pv av = norm(sub(a2, a1));
    pv bv = norm(sub(b2, b1));
    double d2 = dist2(av, bv);
    return CABS0(d2);
}

bool hassegintersect(pv a1, pv a2, pv b1, pv b2)
{
    pv p; 
    COORDX d, n1, n2;
    p = sub(b2, b1);
    d = cross(sub(a2, a1), p);
    n1 = cross(sub(b1, a1), p);
    n2 = cross(sub(b1, a1), sub(a2, a1));
    if (CABS0(d)) 
        return CABS0(n1) && CABS0(n2);
    if (d < 0)
        return CLE(d, n1) && CLE(n1, 0) && CLE(d, n2) && CLE(n2, 0);
    else
        return CLE(0, n1) && CLE(n1, d) && CLE(0, n2) && CLE(n2, d);
} 

pv segintersect(pv a1, pv a2, pv b1, pv b2)
{
    if (!hassegintersect(a1, a2, b1, b2))
        return invalidpv();
    return intersect(a1, a2, b1, b2);
}

int abovepolyseg(pv p, pv a, pv b)
{
    COORDX h, l;
    COORD xa, xb;
    xa = a.x;
    xb = b.x;
    normiv(&xa, &xb);
    if (CL(xa, p.x) && CLE(p.x, xb)) {
        l = (COORDX) (b.x - a.x) * (b.y - p.y);
        h = (COORDX) (b.y - a.y) * (b.x - p.x);
        if ((a.x - b.x) > 0)
            return l < h;
        else 
            return h < l;
    }
    return 0;
}

bool inpoly0(pv p, pv *pp, int n)
{
    int k;
    k = 0;
    FOR0(i, n) {
        if (abovepolyseg(p, pp[i], pp[i + 1]))
            k++;
    }
    return k & 1;
}

pv mid(pv a, pv b)
{
    return add(a, mul(0.5, sub(b, a)));
}

//----------------------------------

enum {
    PPS,
    PPC,
    PPLED,
    A01RS,
    A01RLEDT,
    A01RLEDB,
    A02LEDA1,
    A02LEDA2,
    A02LEDA3,
    A02LEDB1,
    A02LEDB2,
    A02LEDB3,
    A02LEDC1,
    A02LEDC2,
    A02LEDC3,
    A03T,
    A03LED,
    A04T,
    A04LEDT,
    A04LEDB,
    A05T,
    A05LED,
    MAX_ELEMENT
};

enum {
    TOPLEFT,
    TOPRIGHT,
    BOTTOMRIGHT,
    BOTTOMLEFT,
    MAX_CORNER
};

enum {
    UNKNOWN,
    ON,
    OFF,
    UP,
    DOWN,
    CENTER,
    OUT,
    MAX_STATE
};

string elementstr[MAX_ELEMENT] = {"PPS", "PPC", "PPLED", "A01RS", "A01RLEDT", 
 "A01RLEDB", "A02LEDA1", "A02LEDA2", "A02LEDA3", "A02LEDB1", "A02LEDB2", "A02LEDB3", "A02LEDC1", "A02LEDC2", "A02LEDC3",
 "A03T", "A03LED", "A04T", "A04LEDT", "A04LEDB", "A05T", "A05LED"};

string statestr[MAX_STATE] = {"UNKNOWN", "ON", "OFF", "UP", "DOWN", "CENTER", "OUT"};

enum {
    SIM,
    LAB,
    LAB2,
    LAB3,
    ISS,
    MAX_MODE
};

#define MAX_H        2500
#define MAX_W        2500
#define MAX_XY        2500
#define MAX_DOTS    1000

char modestr[MAX_MODE][5] = {"Sim","Lab","Lab2","Lab3","ISS"};

struct pv4 {
    pv p[MAX_CORNER]; 
};

enum {
    BLACK = 0x0,
    RED = 0xFF0000,
    GREEN = 0x00FF00,
    BLUE = 0x0000FF,
    MIDGREY = 0x3F3F3F,
    GREY = 0x7F7F7F,
    WHITE = 0xFFFFFF
};

enum {
    RB_CMIN,
    RB_RMIN,
    RB_HIT128,
    LBL_CREL,
    LBL_RMAX,

    CLINE_GRP,
    CLINE_CNT,
    CLINE_EPS,
    CLINE_DIFF,
    CLINE_XMIN,
    CLINE_XMAX,
    CLINE_YMIN,
    CLINE_YMAX,
    CLINE_BMIN,
    CLINE_BHIT128,

    A01RS_MAXXY,
    A01RS_CIRCLER,
    A01RS_CNT,

    BIGLED_MAXXY,
    BIGLED_CIRCLER,
    BIGLED_CNT,
    BIGLED_MINCIRCLER,
    BIGLED_GREYCNT,

    LED_MAXXY,
    LED_CIRCLER,
    LED_CNT,
    LED_BLACKMAXXY,
    LED_BLACKCIRCLER,
    LED_BLACKCNT,

    MAX_PARAM
};

int param[2][MAX_PARAM] = {{20, 40, 110, 400, 400, 
 5, 19, 2, 40, 50, 600, 50, 800, 20, 110,
 20, 20, 20, 
 15, 24, 5, 10, 10, 
 20, 6, 7, 15, 8, 5}};
bool rparam[MAX_PARAM] = {0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0};
bool cparam[MAX_PARAM] = {1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0};
bool rateparam[MAX_PARAM] = {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
bool panelparam[MAX_PARAM] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; 

pv oripanel[5] = {pv(490,118),pv(992,141),pv(999,880),pv(456,870)};
pv oripos[MAX_ELEMENT] = 
{
    pv(831,238),pv(842,121),pv(873,242),pv(874,470),pv(873,414),pv(873,533),pv(577,398),pv(639,401),pv(704,403),pv(575,461),pv(638,462),
    pv(703,466),pv(572,523),pv(637,526),pv(701,527),pv(578,692),pv(580,635),pv(734,699),pv(734,639),pv(732,762),pv(887,702),pv(888,643)
};

pv elements[MAX_ELEMENT];
pv pos[2][MAX_ELEMENT];
double oridist[5];
pv panelsidenorm[2][5];
pvp paneldir[2][5];
double panelsidedist[2][5], panelsidepersp[2][5];
pv perspect[2][2];
bool parallel[2][2];
pv perspectlineend[2][2];
double oneperoridist[5];
pv center[2];
pv redblock[2][2];
pv panel[2][5];
pvp pleft[2], pright[2], ptop[2], pbottom[2];
bool ppledfix[2];
bool upperpos[2];

bool bigok[2], redblockdown[2], bigrecheck[2];
int h, w, hw, mode;
int npanelside[2], naddpanelcorner[2];
int img[2][MAX_H * MAX_W];
int saveimg[2][MAX_H * MAX_W];
int tmpimg[2][MAX_H * MAX_W];
int sum[2][MAX_H * MAX_W];
int state[MAX_ELEMENT];
int ravg[2], bavg[2], gavg[2], davg[2], avg[2], wbvavg[2], greyavg[2], blackavg[2], baseblack[2];
double oriw, orih;
pv dots[MAX_DOTS], nodots[MAX_DOTS];
int ndot[2], nodot[2];
int leftrightadjust;
bool dotpanel[2];
int nuseddots[2];
bool settled[2][MAX_ELEMENT];
int blackshift[2];
bool lightdowncheck[2];
int minred[2], maxred[2];
bool worseside;
bool bigledfound[2];
int blackpanelside[2];
int minr[2], maxr[2], minc[2], maxc[2];

void debugp4(char *s, pv p1, pv p2, pv p3, pv p4)
{
    DEBUG("%s %d:%d %d:%d %d:%d %d:%d", s, (int) p1.x, (int) p1.y, (int) p2.x, (int) p2.y, (int) p3.x, (int) p3.y, (int) p4.x, (int) p4.y);
}

void setp(int sd, int r, int c, int v)
{
    img[sd][r * w + c] = v;
}

bool isin(int r, int c)
{
    return r >= 0 && r < h && c >= 0 && c < w;
}

int getp(int sd, int r, int c)
{
    return img[sd][r * w + c];
}

int red(int rgb)
{
    return rgb >> 16;
}

int green(int rgb)
{
    return (rgb >> 8) & 0xFF;
}

int blue(int rgb)
{
    return rgb & 0xFF;
}

bool islow(int v, int avg)
{
    return v <= 5 * avg / 8;
}

bool isverylow(int v, int avg)
{
    return v < 20 || v <= 3 * avg / 8;
}

bool ishigh(int v, int avg)
{
    return v >= 3 * avg / 2;
}

bool isveryhigh(int v, int avg)
{
    return v >= 230;
}

bool isavg(int v, int avg)
{
    return v > avg / 2 && v < 5 * avg / 2;
}

bool isdarkred(int sd, int rgb)
{
    if (npanelside[sd] > 0 && mode != SIM && red(rgb) > ravg[sd] / 4)
        return false;
    return red(rgb) > (mode == LAB2 ? 40 : 20) && (mode == SIM || islow(green(rgb), gavg[sd]) && islow(blue(rgb), bavg[sd])) && red(rgb) > 3 * (green(rgb) + blue(rgb)) / 4;
}

bool isblack(int sd, int rgb)
{
    return red(rgb) + green(rgb) + blue(rgb) < max(3 * (12 - blackshift[sd]), blackavg[sd]);
}

bool iswhite(int sd, int rgb)
{
    return green(rgb) + blue(rgb) > 2 * 230; // we don't care about red
}

bool isgreen(int sd, int rgb)
{
    return green(rgb) + blue(rgb) < 2 * 230 && ishigh(green(rgb), gavg[sd]) && green(rgb) > blue(rgb) && green(rgb) > 5 * red(rgb) / 4;
}

bool isblue(int sd, int rgb)
{

    return green(rgb) + blue(rgb) < 2 * 230  && (isveryhigh(blue(rgb), bavg[sd]) || mode == SIM && ishigh(blue(rgb), bavg[sd])) && blue(rgb) > green(rgb) && blue(rgb) > 5 * red(rgb) / 4;
}

bool ismidgrey(int sd, int rgb)
{
    return red(rgb) + green(rgb) + blue(rgb) < 2 * greyavg[sd];
}

double perspectdist(int sd, double mul, double dist, double dist1, double dist2, int iter)
{
    if (iter == 0)
        return mul * dist;
/*    if (dist1 > dist2 * 1.1)
        dist1 *= 1.1;
    else if (dist2 > dist1 * 1.1)
        dist2 *= 1.1; */
    double midd = dist1 / (dist1 + dist2);
    double d = midd * dist;
    double mind = min(dist1, dist2);
    double maxd = max(dist1, dist2);
    double newd = mind + mind * (maxd - mind) / (mind + maxd);
    if (mul > 0.5)
        return d + perspectdist(sd, 2 * (mul - 0.5), dist - d, newd, dist2, iter - 1);
    return perspectdist(sd, 2 * mul, d, dist1, newd, iter - 1);
}

void setparallel(int sd, pv *pp, double *d, int i)
{
//    pp[i] = add(panel[sd][i], mul(d[i] * oneperoridist[i] * panelsidedist[sd][i], panelsidenorm[sd][i])); 
    double based = d[i] * oneperoridist[i] * panelsidedist[sd][i];
    double specd = perspectdist(sd, d[i] * oneperoridist[i], panelsidedist[sd][i],
        panelsidedist[sd][(i + 3) % 4], panelsidedist[sd][i + 1], 5);
    pp[i] = add(panel[sd][i], mul(specd, panelsidenorm[sd][i])); 
}

pv transform(int sd, pv orip)
{
    double d[4];
    FOR0(i, 4) {
        pv cp = closest(orip, oripanel[i], oripanel[i + 1]);
        d[(i + 1) % 4] = dist(orip, cp);
    }
    pv pp[4];
    setparallel(sd, pp, d, 0);
    setparallel(sd, pp, d, 2);
    setparallel(sd, pp, d, 1);
    setparallel(sd, pp, d, 3);
    pv p = intersect(pp[0], pp[2], pp[1], pp[3]);
    if (!parallel[sd][0]) {
        if (perspect[sd][0].x < 0)
            p.x -= 12;
        else
            p.x += 12;
    }
    return p;
}

#define MAX_PERSPECT_DIST    SQR(15000.0)

void transformall()
{
    FOR0(sd, 2) {
        panel[sd][4] = panel[sd][0];
        FOR0(i, 4) {
            panelsidedist[sd][i] = dist(panel[sd][i], panel[sd][i + 1]);
            panelsidenorm[sd][i] = norm(sub(panel[sd][i + 1], panel[sd][i]));
        }
        panelsidedist[sd][4] = panelsidedist[sd][0];
        parallel[sd][0] = true;
        parallel[sd][1] = true;
        perspect[sd][0] = intersect(panel[sd][0], panel[sd][1], panel[sd][2], panel[sd][3]);
        parallel[sd][0] = !isvalid(perspect[sd][0]) || dist2(pv(0, 0), perspect[sd][0]) > MAX_PERSPECT_DIST;
/*        if (!parallel[sd][0]) {
            double d = dist(perspect[sd][0], panel[sd][0]);
            perspectlineend[sd][0] = segpd(perspect[sd][0], panel[sd][3], d);
        }
        perspect[sd][1] = intersect(panel[sd][1], panel[sd][2], panel[sd][3], panel[sd][4]);
        parallel[sd][1] = !isvalid(perspect[sd][1]) || dist2(pv(0, 0), perspect[sd][1]) > MAX_PERSPECT_DIST;
        if (!parallel[sd][0]) {
            double d = dist(perspect[sd][1], panel[sd][0]);
            perspectlineend[sd][1] = segpd(perspect[sd][1], panel[sd][1], d);
        } */
        FOR0(i, MAX_ELEMENT) {
            pos[sd][i] = transform(sd, oripos[i]);
        }
        if (isvalid(redblock[sd][0]) && (ndot[sd] > 0 || naddpanelcorner[sd] > 0)) {
            if (redblock[sd][0].x < pos[sd][PPLED].x - 300 || redblock[sd][1].x > pos[sd][PPLED].x + 200 ||
                redblock[sd][0].y < pos[sd][PPLED].y - 300 || redblock[sd][1].y > pos[sd][PPLED].y + 200)
                redblock[sd][0] = pv();
        }
    }
}


bool checkredline(int sd, int r, int c, int &minc, int& maxc)
{
    int cnt = 1;
    int cd;
    MINA2(minc, c, minred[sd], r);
    for(cd = 1; cd < param[sd][RB_CMIN] * 4; ++cd) {
        if (c + cd >= w)
            break;
        if (getp(sd, r, c + cd) == RED) {
            cnt++;
            MINA2(minc, c + cd, minred[sd], r);
            MAXA2(maxc, c + cd, maxred[sd], r);
        }
        if (cd > 2 + param[sd][RB_CMIN] / 4 && cnt < param[sd][RB_HIT128] * cd / 128)
            break;
    }
    if (cd < param[sd][RB_CMIN])
        return false;
    return true;
}

bool checkredblock(int sd, int r, int c, int& minc, int& maxc, int& maxr)
{
    int cnt = 1;
    int rd;
    for(rd = 1; rd < param[sd][RB_RMIN] * 8; ++rd) {
        if (r + rd >= h)
            break;
        if (isin(r + rd, c - 1) && getp(sd, r + rd, c - 1) == RED) {
            --c;
            if (isin(r + rd, c - 1) && getp(sd, r + rd, c - 1) == RED)
                --c;
        }
        else if (isin(r + rd, c) && getp(sd, r + rd, c) != RED) {
            ++c;
            if (isin(r + rd, c) && getp(sd, r + rd, c) != RED)
                ++c;
        }
        if (checkredline(sd, r + rd, c, minc, maxc)) {
            maxr = r + rd;        
            ++cnt;
        } else if (cnt < 3 * rd / 4)
            break;
    }
    if (rd < param[sd][RB_RMIN]) 
        return false;
    return true;
}

void findredblock(int sd)
{
    FOR(r, param[sd][RB_RMIN], h - param[sd][RB_RMIN]) {
        FOR(c, param[sd][RB_CMIN], w - param[sd][RB_CMIN]) {
            if (getp(sd, r, c) != RED)
                continue;
            int minc = INT_CAP;
            int maxc = 0;
            int maxr = 0;
            if (!checkredline(sd, r, c, minc, maxc))
                continue;
            if (checkredblock(sd, r, c, minc, maxc, maxr)) {
                redblock[sd][0] = pv(minc, r);
                redblock[sd][1] = pv(maxc, maxr);
                return;
            }
            c += 5; 
        }
    }
}

void calcavg(int sd)
{
    int rtotal, gtotal, btotal, dtotal;
    rtotal = gtotal = btotal = dtotal = 0;
    FOR0(r, h) {
        int rgbprev = 0;
        FOR0(c, w) {
            int rgb = getp(sd, r, c);
            rtotal += red(rgb);
            gtotal += green(rgb);
            btotal += blue(rgb);
            if (c) {
                dtotal += (abs(red(rgbprev) - red(rgb)) + abs(green(rgbprev) - green(rgb)) + abs(blue(rgbprev) - blue(rgb)));
            }
            if (r < h - 1) {
                int rgbrow = getp(sd, r + 1, c);
                dtotal += (abs(red(rgbrow) - red(rgb)) + abs(green(rgbrow) - green(rgb)) + abs(blue(rgbrow) - blue(rgb)));
            }
            if (isblack(sd, rgb))
                baseblack[sd]++;
            rgbprev = rgb;
        }
    }
    ravg[sd] = rtotal / hw;
    gavg[sd] = gtotal / hw;
    bavg[sd] = btotal / hw;
    avg[sd] = int(((LI)rtotal + gtotal + btotal) / hw / 3);
    davg[sd] = dtotal / hw;
    wbvavg[sd] = avg[sd];
    greyavg[sd] = avg[sd];
    blackavg[sd] = avg[sd];
}

void setmode()
{
    int da = (davg[0] + davg[1]) / 2;
    int aa = (avg[0] + avg[1]) / 2;
    if (h == 1200) {
        if (da <= 12)
            mode = SIM;
        else if (baseblack[0] < hw / 8)
            mode = LAB3;
        else
            mode = LAB2;
    } else {
        if (da <= 17 && aa >= 30)
            mode = LAB;
        else
            mode = ISS;
    }
}

void adjustparams()
{
    // some params can be totally different depending on situations
    FOR0(sd, 2) {
        FOR0(i, MAX_PARAM) {
            if (rparam[i])
                param[sd][i] = h * param[sd][i] / 1200;
            else if (cparam[i])
                param[sd][i] = w * param[sd][i] / 1600;
            if (davg[0] >= 20 && rateparam[i])
                param[sd][i] = 10 * param[sd][i] / 16;
        }
    }
}

int smoothv(int v1, int v2, int v3)
{
    return (2 * v1 + v2 + v3) / 4;
}

void smooth(int sd)
{
    FOR(r, 1, h - 1) {
        FOR(c, 1, w - 1) {
            int rgb1 = getp(sd, r, c);
            int rgb2 = getp(sd, r - 1, c);
            int rgb3 = getp(sd, r, c - 1);
            setp(sd, r, c, (smoothv(red(rgb1), red(rgb2), red(rgb3)) << 16) +
                (smoothv(green(rgb1), green(rgb2), green(rgb3)) << 8) +
                smoothv(blue(rgb1), blue(rgb2), blue(rgb3)));
        }
    }
}

void init()
{
    hw = h * w;
    leftrightadjust = 340 * w / 1600;
    oriw = ((oripanel[1].x - oripanel[0].x) + (oripanel[2].x - oripanel[3].x)) / 2;
    orih = ((oripanel[2].y - oripanel[1].y) + (oripanel[3].y - oripanel[0].y)) / 2;
    FOR0(sd, 2)
        memcpy(saveimg[sd], img[sd], sizeof(int) * h * w);
    FOR0(i, MAX_PARAM)
        param[1][i] = param[0][i];
    oripanel[4] = oripanel[0];
    oneperoridist[0] = 1.0 / oriw;
    oneperoridist[1] = 1.0 / orih;
    oneperoridist[2] = 1.0 / oriw;
    oneperoridist[3] = 1.0 / orih;
    FOR0(sd, 2) {
        minc[sd] = 0;
        minr[sd] = 0;
        maxc[sd] = w - 1;
        maxr[sd] = h - 1;
    }
}

void filldefault()
{
    FOR0(i, MAX_ELEMENT) {
        state[i] = UNKNOWN;
        FOR0(sd, 2)
            pos[sd][i] = pv(200 + 200 * sd, 200 + 50 * i);
    }
}

int wbv(int sd, int v, int cavg)
{
//    return (v + 3 * v * avg[sd] / cavg) / 4;
    return min(0xFF, v * wbvavg[sd] / cavg);
}

void lightup(int sd)
{
    static bool lit[2] = {false, false};
    if (avg[sd] >= 25)
        return;
    FOR0(r, h) {
        FOR0(c, w) {
            int rgb = getp(sd, r, c);
            setp(sd, r, c, (min(3 * red(rgb) / 2, 0xFF) << 16) +
                (min(3 * green(rgb) / 2, 0xFF) << 8) + min(3 * blue(rgb) / 2, 0xFF));
        }
    }
    if (!lit[sd]) {
        avg[sd] *= 2;
        ravg[sd] *= 2;
        bavg[sd] *= 2;
        gavg[sd] *= 2;
        davg[sd] *= 2; 
        lit[sd] = true;
    }
}

void recolor(int sd)
{
    int midgrey = 0;
    int black = 0;
    int first = (mode == SIM || mode == LAB) ? 1 : 0;
    memcpy(tmpimg[sd], img[sd], sizeof(int) * hw);
    FOR(doit, first, 1) {
        if (doit && !first)
            memcpy(img[sd], tmpimg[sd], sizeof(int) * hw);
        FOR(r, minr[sd], maxr[sd]) {
            FOR(c, minc[sd], maxc[sd]) {
                int rgb = getp(sd, r, c);
                if (mode != SIM) {
                    rgb = (wbv(sd, red(rgb), ravg[sd]) << 16) +
                        (wbv(sd, green(rgb), gavg[sd]) << 8) +
                        wbv(sd, blue(rgb), bavg[sd]); 
                } 
                if (isdarkred(sd, rgb))
                    rgb = RED;
                else if (isgreen(sd, rgb))
                    rgb = GREEN;
                else if (isblue(sd, rgb))
                    rgb = BLUE;
                else if (isblack(sd, rgb)) {
                    rgb = BLACK;
                    black++;
                }
                else if (iswhite(sd, rgb))
                    rgb = WHITE;
                else if (ismidgrey(sd, rgb)) {
                    rgb = MIDGREY; 
                    midgrey++;
                } else
                    rgb = GREY;
                setp(sd, r, c, rgb);
            }
        }
        bool change = false;
        if (mode == LAB2) {
            if (sd == 1)
                greyavg[sd] = 3 * avg[sd] / 5;  
            else
                greyavg[sd] = 3 * avg[sd] / 4;  
            blackavg[sd] = 3 * avg[sd] / 4;
            blackshift[sd] = 4;
            change = true;
        } else {
            if (black > hw / 3 || avg[sd] < 30) {
                blackavg[sd] = 3 * avg[sd] / 4;
                blackshift[sd] = 4;
                change = true;
            } 
            if (midgrey > hw / 4) {
                greyavg[sd] = avg[sd] = 3 * avg[sd] / 4;  //!!! argh avg shouldnt be changed, but its too late to retest
                change = true;
            }
        }
        if (!change)
            break;
    }
    if (mode == LAB) { // last minute hack away that face
        FOR0(r, 150) {
            FOR0(c, 150) {
                int rgb = getp(sd, r, c);
                if (rgb == RED)
                    setp(sd, r, c, GREY);
                c += 1;
            }
            r += 1;
        }
    } 
}

pv guess_center_redblock(int sd)
{
    double rh = redblock[sd][1].y - redblock[sd][0].y;
    double rw = redblock[sd][1].x - redblock[sd][0].x;
    pv p;
    if (rh > 4 * rw / 5) { // old/buggy...but got optimization for it
        redblockdown[sd] = true;
        p = pv(redblock[sd][0].x - min(rw, 50), redblock[sd][1].y + rh);
    } else
        p = pv(redblock[sd][0].x, redblock[sd][1].y + min(200, 2 * rh));
    if (p.x < w / 4)
        p.x = max(200, redblock[sd][1].x);
    else if (p.x > 3 * w / 4)
        p.x = min(w - 300, redblock[sd][0].x - rw);
    if (p.y < h / 6 || p.y > 5 * h / 6)
        p.y = min(h - 200, redblock[sd][1].y);
    return p;
}

bool isblacksegment(int sd, int *buff, int rr, int d, int xydmax)
{
    int cnt = 1;
    FOR(qq, 1, param[sd][CLINE_BMIN]) {
        int rgb = buff[xydmax + d * (rr + qq)];
        if (rgb == BLACK)
            cnt++;
        else {
            if (qq >= 4 && cnt < param[sd][CLINE_BHIT128] * qq / 128)
                return false;
        }
    }
    return true;
}

void getpaneledge(int sd, int *buff, int d, int& edge, int xydmin, int xydmax)
{
    edge = -1;
    if (isblacksegment(sd, buff, xydmin, d, xydmax))
        return;
    FOR(rr, xydmin + 1, xydmax) {
        int rgb = buff[xydmax + d * rr];
        if (rgb & 0xFF000000) // stopper
            return;
        if (rgb == BLACK) {
            if (isblacksegment(sd, buff, rr, d, xydmax)) {
                edge = xydmax + d * rr;
                return;
            }
        }
    }
}

void getpaneledge(int sd, bool horizontal, int relpos, int *edge, int xydmin, int xydmax)
{
    int buff[2 * MAX_XY];
    memset(buff, 0xFF, sizeof(buff));
    edge[0] = edge[1] = -1;
    if (horizontal) {
        int r = (int) center[sd].y + relpos;
        if (r < 0 || r >= h)
            return;
        int rc = (int) center[sd].x - xydmax;
        int c = max(0, rc);
        FOR(dx, c - rc, 2 * xydmax) {
            if (c >= w)
                break;
            buff[dx] = getp(sd, r, c);
            ++c;
        }
    } else {
        int c = (int) center[sd].x + relpos;
        if (c < 0 || c >= w)
            return;
        int rr = (int) center[sd].y - xydmax;
        int r = max(0, rr);
        FOR(dy, r - rr, 2 * xydmax) {
            if (r >= h)
                break;
            buff[dy] = getp(sd, r, c);
            ++r;
        }
    }
    getpaneledge(sd, buff, -1, edge[0], xydmin, xydmax);
//    getpaneledge(sd, buff, 1, edge[1], xydmin, xydmax);
    int saveval = param[sd][CLINE_BHIT128];
    if (!isvalid(redblock[sd][0]) && horizontal)
        param[sd][CLINE_BHIT128] = 128; //need all blacks to avoid darkredblock
    getpaneledge(sd, buff, 1, edge[1], !horizontal ? xydmin : max(xydmin, 
        (int) (redblock[sd][0].x + redblock[sd][1].x) / 2 - (int) center[sd].x), xydmax);
    param[sd][CLINE_BHIT128] = saveval;
    FOR0(s, 2) {
        if (edge[s] >= 0)
            edge[s] += (int) (horizontal ? center[sd].x : center[sd].y) - xydmax;
    }
}

pvp getpanelline(int sd, bool horizontal, int relpos)
{
    pvp res;
    int cnt[2];
    cnt[0] = cnt[1] = 0;
    int sedge[2];
    sedge[0] = sedge[1] = 0;
    FOR0(i, param[sd][CLINE_GRP]) {
        int edge[2];
        getpaneledge(sd, horizontal, relpos + i - param[sd][CLINE_GRP] / 2, edge,
            horizontal ? param[sd][CLINE_XMIN] : param[sd][CLINE_YMIN], horizontal ? param[sd][CLINE_XMAX] : param[sd][CLINE_YMAX]);
        FOR0(s, 2) {
            if (edge[s] >= 0 && (cnt[s] == 0 || abs(sedge[s] / cnt[s] - edge[s]) <= 10)) {
                sedge[s] += edge[s];
                cnt[s]++;
            }
        }
    }
    FOR0(s, 2) {
        if (!isvalid(redblock[sd][0]) && horizontal && relpos  < -200 && s == 1) // ignore upper part of right side (redblock)
            continue;
        if (cnt[s] >= 3) {
            if (horizontal) {
                res.p[s].x = sedge[s] / cnt[s];
                res.p[s].y = center[sd].y + relpos;
            } else {
                res.p[s].y = sedge[s] / cnt[s];
                res.p[s].x = center[sd].x + relpos;
            }
        }
    }
    return res;
}

bool almostonline(pv p, pv p1, pv p2, int cap)
{
    pv tp = closest(p, p1, p2);
    if (!isvalid(p))
        return false;
    double d = dist(p, tp);
    return d < cap;
}

bool badintersect(pv a1, pv a2, pv b1, pv b2)
{
    pv p = segintersect(a1, a2, b1, b2);
    return isvalid(p);
/*    if (!isvalid(p))
        return true;
    return onseg(p, b1, b2); */
}

bool getpanelside(int sd, bool horizontal, pvp *clines, int s, pvp& res, int markpos, int cap)
{
    if (!horizontal) {
        FOR0(i, param[sd][CLINE_CNT]) {
            if (isvalid(pleft[sd].p[0]) && almostonline(clines[i].p[s], pleft[sd].p[0], pleft[sd].p[1], 5) ||
              isvalid(pright[sd].p[0]) && almostonline(clines[i].p[s], pright[sd].p[0], pright[sd].p[1], 5) ||
              isvalid(pleft[sd].p[0]) && side(clines[i].p[s], pleft[sd].p[0], pleft[sd].p[1]) == 1) // single on the left is bad
                clines[i].p[s] = pv();
        }
    } 
    if (testparam == 2 && !horizontal && s == 1 || testparam == 3 && horizontal && s == 1 ||
        testparam == 4 && horizontal && s == 0) {
        int dmark = MAX_ELEMENT - 1;
        FOR0(i, param[sd][CLINE_CNT]) {
            if (isvalid(clines[i].p[s])) {
                pos[sd][dmark] = clines[i].p[s];
                --dmark;
            }
        }
    } 
    FOR0(i, param[sd][CLINE_CNT] - 3) {
        if (!isvalid(clines[i].p[s]))
            continue;
        FORD(j, param[sd][CLINE_CNT] - 1, i + 1) {
            if (!isvalid(clines[j].p[s]))
                continue;
            int cnt = 0;
            FOR(k, i + 1, j - 1) {
                if (!isvalid(clines[k].p[s]))
                    continue;
                if (almostonline(clines[k].p[s], clines[i].p[s], clines[j].p[s], param[sd][CLINE_EPS])) {
                    cnt++;
                    if (cnt >= cap)
                        break;
                }
            }
            if (cnt >= cap) {
                if (!horizontal) {
                    if (isvalid(pright[sd].p[0])) {
                        if (cnt <= 2 && badintersect(clines[i].p[s], clines[j].p[s], pright[sd].p[0], pright[sd].p[1]))
                            continue;
                        if (side(clines[i].p[s], pright[sd].p[0], pright[sd].p[1]) == -1 && 
                            side(clines[j].p[s], pright[sd].p[0], pright[sd].p[1]) == -1) {
                            continue;
                        }
                    }
                }
                res.p[0] = clines[i].p[s];
                res.p[1] = clines[j].p[s];
                pos[sd][markpos] = res.p[0];
                pos[sd][markpos+1] = res.p[1];
                npanelside[sd]++;
                return true;
            }
        }
    }
    //!! tovabb lehet engedni, mar eleg megbizhatoak vagyunk, hogy akar 2 pontra is illesszuk
    // ?ISS CIZA, SIM PJEG
    // !!! de nar nem biztos, hogy kell a jobb-bal igazitas miatt
    if (!horizontal && s == 1 && (mode == LAB || mode == ISS || mode == SIM)) {
        if (cap <= 1)
            return false;
    } else
        if (cap <= 2)
            return false;
    return getpanelside(sd, horizontal, clines, s, res, markpos, cap - 1);
}

bool collectpanelside(int sd, bool horizontal, pvp& left, pvp& right)
{
    pvp clines[20];
    FOR0(i, param[sd][CLINE_CNT])
        clines[i] = getpanelline(sd, horizontal, (i - param[sd][CLINE_CNT] / 2) * param[sd][CLINE_DIFF]);
    bool ok = getpanelside(sd, horizontal, clines, 0, left, 3 + horizontal * 4, 5);
    ok = getpanelside(sd, horizontal, clines, 1, right, 3 + horizontal * 4 + 2, 5) && ok;
    return ok;
}

void adjustpaneltop(int sd)
{
    if (npanelside[sd] < 4 || dotpanel[0] || dotpanel[1])
        return;
    double topdist = dist(panel[sd][0], panel[sd][1]);
    double bottomdist = dist(panel[sd][2], panel[sd][3]);
    if (topdist <= bottomdist * 1.05)
        return;
    panel[sd][0] = segp(panel[sd][0], panel[sd][3], 0.02);
    panel[sd][1] = segp(panel[sd][1], panel[sd][2], 0.02);
}

int  specledv(int rgb)
{
    switch(rgb) {
    case GREEN: return 12;
    case BLUE: return 6;
    case BLACK: return 3;
    case WHITE: return 2;
    case RED: case MIDGREY: return 1;
    }
    return 0;
}

pv specledpos(int& bestcmpcnt, int sd, pv inp, int maxxy, int circler, int mincnt)
{
    bestcmpcnt = mincnt;
    int bestdiff = INT_CAP;
    int circler2 = SQR(circler);
    pv bestp;
    FOR(r, (int) inp.y - maxxy, (int) inp.y + maxxy) {
        FOR(c, (int) inp.x - maxxy, (int) inp.x + maxxy) {
            int cnt = 0;
            int outcnt = 0;
            int left = 0;
            int right = 0;
            int top = 0;
            int bottom = 0;
            int primary = 0;
            int greycnt = 0;
            FOR(rr, r - circler, r + circler) {
                FOR(cc, c - circler, c + circler) {
                    if (!isin(rr, cc))
                        continue;
                    int d2 = SQR(r - rr) + SQR(c - cc);
                    int rgb = getp(sd, rr, cc);
                    if (rgb != GREY) {
                        if (d2 <= circler2) {
                            cnt += specledv(rgb);
                            if (cc - c < 0)
                                left++;
                            else if (cc - c > 0)
                                right++;
                            if (rr - r < 0)
                                top++;
                            else if (rr - r > 0)
                                bottom++;
                        } else {
                            outcnt += specledv(rgb);
                            if (rgb == GREEN || rgb == BLUE)
                                outcnt += 4 * specledv(rgb); // hard, to avoid useing bigleds
                        }
                    } else if (d2 <= circler2)
                        --greycnt;
                }
                if (rr > r && cnt < bestcmpcnt / 3)
                    break;
            }
            outcnt = max(0, outcnt - greycnt);
            int cmpcnt = cnt;
            cmpcnt -= 3 * outcnt;
            cmpcnt -= (abs(right - left) + abs(top - bottom)) / 2;
            if (cmpcnt > bestcmpcnt) {
                bestcmpcnt = cmpcnt;
                bestp.y = r;
                bestp.x = c;
            }
            c++;
        }
        r++;
    }
    return bestp;
}



pv bestledpos(int& bestcnt, int sd, pv inp, int maxxy, int circler, int color, int mincnt, int mincircle = 0, bool nocircle = false, bool singlecolor = false)
{
    maxxy = min(150, maxxy);
    bestcnt = mincnt;
    int bestcmpcnt = mincnt;
    int bestdiff = INT_CAP;
    int circler2 = SQR(circler);
    int mincircler2 = SQR(mincircle);
    int color2 = color;
    int color3 = color;
    if (!singlecolor) {
        if (color == BLACK)
            color2 = MIDGREY;
        else if (color == MIDGREY) {
            color2 = GREEN;
            color3 = BLUE;
        } else if (color == RED)
            color2 = BLACK;
        else {
            color2 = WHITE;
//            if (mincircle > 0 && color == BLUE)
//                color3 = GREEN;
//            else 
            if (mode != SIM && !mincircle)
                color3 = BLACK;
        }
    }
    pv bestp;
    int my = maxxy;
    if (lightdowncheck[sd])
        my = 3 * my / 4;
    FOR(r, (int) inp.y - my, (int) inp.y + my) {
        FOR(c, (int) inp.x - maxxy, (int) inp.x + maxxy) {
            int cnt = 0;
            int outcnt = 0;
            int left = 0;
            int right = 0;
            int top = 0;
            int bottom = 0;
            int primary = 0;
            FOR(rr, r - circler, r + circler) {
                FOR(cc, c - circler, c + circler) {
                    if (!isin(rr, cc))
                        continue;
                    int d2 = SQR(r - rr) + SQR(c - cc);
                    int rgb = getp(sd, rr, cc);
                    if (rgb == color || rgb == color2 || rgb == color3) {
                        if (d2 <= circler2 && d2 >= mincircler2 || nocircle) {
                            cnt++;
                            if (cc - c < 0)
                                left++;
                            else if (cc - c > 0)
                                right++;
                            if (rr - r < 0)
                                top++;
                            else if (rr - r > 0)
                                bottom++;
                            if (rgb == color)
                                primary++;
                        } else
                            outcnt++;
                    }
                    else if (mincircle > 0 && rgb == RED) // human hand
                        outcnt++;
                }
                if (rr > r && cnt < bestcnt / 3)
                    break;
            }
            int cmpcnt = cnt;
            cmpcnt -= outcnt / 2;
            cmpcnt -= (abs(right - left) + abs(top - bottom)) / 2;
            cmpcnt += primary / 4;
            if (primary >= cnt / 8) {
                if (cmpcnt > bestcmpcnt) {
                    bestcnt = cnt;
                    bestcmpcnt = cmpcnt;
                    bestp.y = r;
                    bestp.x = c;
                }
            }
            c += (h == 1200 ? 1 : 2);
        }
        r += (h == 1200 ? 1 : 2);
    }
    return bestp;
}

void adjustforppled()
{
//    if (mode != SIM)
//        return;
    FOR0(sd, 2) {
        upperpos[sd] = dist(panel[sd][0], panel[sd][1]) > 1.1 * dist(panel[sd][2], panel[sd][3]);
        int cnt;
        pv p = bestledpos(cnt, sd, pos[sd][PPLED], 3 * param[sd][LED_MAXXY], param[sd][LED_CIRCLER], GREEN, SQR(param[sd][LED_CNT]));
        if (!isvalid(p)) {
            if (isvalid(redblock[sd][0]) && pos[sd][PPLED].x > redblock[sd][1].x && pos[sd][PPLED].x < w) {
                int adjustr = 0;
                if (upperpos[sd] && mode != SIM)
                    adjustr = 3 * h / 1200;
                pv basepos = pos[sd][PPLED];
                basepos.x += param[sd][LED_BLACKMAXXY] + adjustr;
                p = bestledpos(cnt, sd, basepos, param[sd][LED_BLACKMAXXY] + adjustr, param[sd][LED_BLACKCIRCLER] + adjustr, BLACK, SQR(param[sd][LED_BLACKCNT]));
            }
        } else
            state[PPLED] = true;
        if (isvalid(p)) {
            pv dv = sub(p, pos[sd][PPLED]);
            dv.x = min(dv.x, 3 * param[sd][LED_MAXXY] / 2);
            dv.y = min(dv.y, 3 * param[sd][LED_MAXXY] / 2);
//            if (dv.x > 0)
//                dv.x /= 2;
//            dv.y /= 2;
            DEBUG("adjustppled: %d %d:%d cnt %d", sd, (int) p.x, (int) p.y, cnt);
            FOR0(i, MAX_ELEMENT) {
                if (i == PPC && isvalid(redblock[sd][0]))
                    continue;
                if (i == PPLED)
                    pos[sd][i] = p; 
                else
                    pos[sd][i] = add(pos[sd][i], dv); 
            } 
            ppledfix[sd] = true;
        }
    }
}
/*
void collectbigleds(int sd)
{
    pv p;
    FOR(r, 10, h - 10) {
        int lastrgb = -1;
        FOR(c, 10, w - 10) {
            if (isin(r, c)) {
                int rgb = getp(sd, r, c);
                if ((rgb == BLUE || rgb == GREEN) && lastrgb == rgb && getp(sd, r + 1, c) == rgb) {
                    int cnt;
                    pv checkpos;
                    checkpos.y = r + 20;
                    checkpos.x = c;
                    p = bestledpos(cnt, sd, checkpos, 5, 20, rgb, 80, 15, false, true);
                    if (isvalid(p)) {
                        bool ok = true;
                        FOR0(i, nbigleds[sd]) {
                            if (dist2(bigleds[i], tp) < SQR(50)) {
                                ok = false;
                                break;
                            }
                        }
                        if (ok) {
                            DEBUG("bigled: %4.0f:%4.0f", p.x, p.y);
                            bigledcolor[nbigled[sd]++] = rgb;
                            bigleds[nbigled[sd]++] = p;
                        }
                    }
                }
                lastrgb = rgb;
            }
        }
    }
    return p;
}
*/

pv findabigled(int sd)
{
    pv p;
    LI start = gettimems();
    FOR(r, 10, h - 10) {
        int lastrgb = -1;
        FOR(c, 10, w - 10) {
            if (isin(r, c)) {
                int rgb = getp(sd, r, c);
                if ((rgb == BLUE || rgb == GREEN) && lastrgb == rgb && getp(sd, r + 1, c) == rgb) {
                    int cnt;
                    pv checkpos;
                    checkpos.y = r + param[sd][BIGLED_CIRCLER];
                    checkpos.x = c;
                    p = bestledpos(cnt, sd, checkpos, 5, param[sd][BIGLED_CIRCLER], rgb, SQR(15), param[sd][BIGLED_MINCIRCLER]);
                    if (isvalid(p)) {
                        if (cnt < 1000) { // weird counts are possible
                            int sqwhite = 0;
                            int sqcnt = 0;
                            int dd = param[sd][BIGLED_CIRCLER] + 4;
                            int circler2 = SQR(dd);
                            FOR(dr, -dd, dd) {
                                FOR(dc, -dd, dd) {
                                    if (isin(r + dr, c + dc)) {
                                        int d2 = SQR(dr) + SQR(dc);
                                        if (d2 > circler2) {
                                            int rgb2 = getp(sd, (int)p.y + dr, (int)p.x + dc);
                                            if (rgb2 == WHITE)
                                                sqwhite++;
                                            else if (rgb2 == rgb)
                                                sqcnt++;
                                        }
                                    }
                                }
                            }
                            if (sqcnt + sqwhite < 50) {
                                bigledfound[sd] = true;
                                DEBUG("findabig success %d time: %d", sd, (int) (gettimems() - start));
                                return p;
                            }
                        }
                        r += 2 * param[sd][BIGLED_CIRCLER];
                    }
                }
                lastrgb = rgb;
            }
        }
    }
    DEBUG("findabig failure %d time: %d", sd, (int) (gettimems() - start));
    return p;
}


void calcsums(int sd)
{
    FOR0(r, h) {
        FOR0(c, w) {
            int rgb = getp(sd, r, c);
            int v = specledv(rgb);
            if (r > 0)
                v += sum[sd][(r - 1) * w + c];
            if (c > 0)
                v += sum[sd][r * w + c - 1];
            if (r > 0 && c > 0)
                v -= sum[sd][(r - 1) * w + c - 1];
            sum[sd][r * w + c] = v;
        }
    }
}

void collectdots3(int sd)
{
    const int bd = 30;
    const int dd = 20;
    const int mincnt = 400;
    const int maxcnt = 3000;
    int dmark = MAX_ELEMENT - 1;
    FOR0(i, MAX_ELEMENT)
        pos[sd][i] = pv(200 + 200 * sd, 200 + 50 * i);
    if (bigledfound[sd])
        pos[sd][4] = ::center[sd];
    for(int r = 0; r < h - bd; r += 6) {
        int lastrgb = -1;
        for(int c = 0; c < w - bd; c += 6) {
            int rgb = getp(sd, r, c);
            if (rgb == GREY && getp(sd, r + 1, c) == GREY && getp(sd, r, c + 1) == GREY) {
                int v = sum[sd][r * w + c] + sum[sd][(r + bd) * w + c + bd] -
                    sum[sd][(r + bd) * w + c] - sum[sd][r * w + c + bd];
                if (v >= mincnt && v <= maxcnt) {
                    pv p;
                    p.y = r + (dd + bd) / 2 / 2;
                    p.x = c + (dd + bd) / 2 / 2;
                    bool ok = true;
                    FOR0(i, ndot[sd]) {
                        double d = dist(dots[i], p);
                        if (d < bd) {
                            ok = false;
                            break;
                        }
                    }
                    FOR0(i, nodot[sd]) {
                        double d = dist(nodots[i], p);
                        if (d < bd) {
                            ok = false;
                            break;
                        }
                    }
                    if (ok) {
                        pv tp;
                        if (sd == 1 && c >= 537 && r >= 700) //!!!!
                            p.y--;
                        int cnt = 0;
                        int cc = (dd + bd) / 2 / 2;
                        tp = specledpos(cnt, sd, p, 4, cc, mincnt);
                        if (!isvalid(tp) && v > 800) {
                            cc = bd / 2;
                            p.y = r + cc;
                            p.x = c + cc;
                            tp = specledpos(cnt, sd, p, 2, cc, 3 * mincnt / 2);
                        }
                        if (!isvalid(tp) && v < 1200) {
                            cc = dd / 2;
                            p.y = r + cc;
                            p.x = c + cc;
                            tp = specledpos(cnt, sd, p, 6, cc, 3 * mincnt / 4);
                        }
                        if (isvalid(tp)) {
                            if (v > 1400) {
                                pv sp = bestledpos(cnt, sd, tp, 0, cc, GREEN, SQR(param[sd][LED_CNT]));
                                if (!isvalid(sp))
                                    sp = bestledpos(cnt, sd, tp, 0, cc, BLUE, SQR(param[sd][LED_CNT]));
                                ok = isvalid(sp);
                            } 
                            if (ok) {
                                DEBUG("ndot: %d %d %d:%d v %d bestcmpcnt %d", sd, dmark, (int) tp.x, (int) tp.y, v, cnt);
                                dots[ndot[sd]++] = tp;
                                if (testparam == 2) {
                                    if (dmark >= 5) {
                                        pos[sd][dmark] = tp;
                                        --dmark;
                                    }
                                }
                            }
/*                            else
                                nodots[nodot[sd]++] = tp;  */
                        }
                    }
                }
            }
        }
    }
}

void collectdots2(int sd)
{
    const int bd = 30;
    const int dd = 20;
    const int mincnt = 400;
    const int maxcnt = 4000;
    int dmark = MAX_ELEMENT - 1;
    FOR0(i, MAX_ELEMENT)
        pos[sd][i] = pv(200 + 200 * sd, 200 + 50 * i);
    if (bigledfound[sd])
        pos[sd][4] = ::center[sd];
    for(int r = 0; r < h - bd; r += 6) {
        int lastrgb = -1;
        for(int c = 0; c < w - bd; c += 6) {
            int rgb = getp(sd, r, c);
            if (rgb == GREY && lastrgb == GREY && getp(sd, r + 1, c) == GREY) {
                int v = sum[sd][r * w + c] + sum[sd][(r + bd) * w + c + bd] -
                    sum[sd][(r + bd) * w + c] - sum[sd][r * w + c + bd];
                if (v >= mincnt && v <= maxcnt) {
                    pv p;
                    p.y = r + (dd + bd) / 2 / 2;
                    p.x = c + (dd + bd) / 2 / 2;
                    bool ok = true;
                    FOR0(i, ndot[sd]) {
                        double d = dist(dots[i], p);
                        if (d < bd) {
                            ok = false;
                            break;
                        }
                    }
                    FOR0(i, nodot[sd]) {
                        double d = dist(nodots[i], p);
                        if (d < bd) {
                            ok = false;
                            break;
                        }
                    }
                    if (ok) {
                        pv tp;
                        if (sd == 1 && c >= 537 && r >= 700) //!!!!
                            p.y--;
                        int cnt = 0;
                        int cc = (dd + bd) / 2 / 2;
                        tp = specledpos(cnt, sd, p, 4, cc, mincnt);
                        if (!isvalid(tp) && v > 600) {
                            cc = bd / 2;
                            p.y = r + cc;
                            p.x = c + cc;
                            tp = specledpos(cnt, sd, p, 2, cc, 3 * mincnt / 2);
                        }
                        if (!isvalid(tp) && v < 1500) {
                            cc = dd / 2;
                            p.y = r + cc;
                            p.x = c + cc;
                            tp = specledpos(cnt, sd, p, 6, cc, 3 * mincnt / 4);
                        }
                        if (isvalid(tp)) {
                            if (v > 1000) {
                                pv sp = bestledpos(cnt, sd, tp, 0, cc, GREEN, SQR(param[sd][LED_CNT]));
                                if (!isvalid(sp))
                                    sp = bestledpos(cnt, sd, tp, 0, cc, BLUE, SQR(param[sd][LED_CNT]));
                                ok = isvalid(sp);
                            } 
                            if (ok) {
                                DEBUG("ndot: %d %d %d:%d v %d bestcmpcnt %d", sd, dmark, (int) tp.x, (int) tp.y, v, cnt);
                                dots[ndot[sd]++] = tp;
                                if (testparam == 2) {
                                    if (dmark >= 5) {
                                        pos[sd][dmark] = tp;
                                        --dmark;
                                    }
                                }
                            }
/*                            else
                                nodots[nodot[sd]++] = tp;  */
                        }
                    }
                }
            }
            lastrgb = rgb;
        }
    }
}

void collectdots(int sd)
{
    if (mode == LAB2)
        collectdots2(sd);
    else
        collectdots3(sd);
}


void collectminidots(int sd)
{
    const int bd = 20;
    const int dd = 12;
    const int mincnt = 200;
    const int maxcnt = 800;
    int dmark = MAX_ELEMENT - 1;
    FOR0(i, MAX_ELEMENT)
        pos[sd][i] = pv(200 + 200 * sd, 200 + 50 * i);
    for(int r = 0; r < h - bd; r += 6) {
        int lastrgb = -1;
        for(int c = 0; c < w - bd; c += 6) {
            int rgb = getp(sd, r, c);
            if (rgb == GREY && lastrgb == GREY && getp(sd, r + 1, c) == GREY) {
                int v = sum[sd][r * w + c] + sum[sd][(r + bd) * w + c + bd] -
                    sum[sd][(r + bd) * w + c] - sum[sd][r * w + c + bd];
                if (v >= mincnt && v <= maxcnt) {
                    pv p;
                    p.y = r + (dd + bd) / 2 / 2;
                    p.x = c + (dd + bd) / 2 / 2;
                    bool ok = true;
                    FOR0(i, ndot[sd]) {
                        double d = dist(dots[i], p);
                        if (d < bd) {
                            ok = false;
                            break;
                        }
                    }
                    FOR0(i, nodot[sd]) {
                        double d = dist(nodots[i], p);
                        if (d < bd) {
                            ok = false;
                            break;
                        }
                    }
                    if (ok) {
                        pv tp;
                        if (sd == 1 && c >= 537 && r >= 700) //!!!!
                            p.y--;
                        int cnt = 0;
                        int cc = (dd + bd) / 2 / 2;
                        tp = specledpos(cnt, sd, p, 4, cc, mincnt);
                        if (!isvalid(tp) && v > 300) {
                            cc = bd / 2;
                            p.y = r + cc;
                            p.x = c + cc;
                            tp = specledpos(cnt, sd, p, 2, cc, 3 * mincnt / 2);
                        }
                        if (!isvalid(tp) && v < 500) {
                            cc = dd / 2;
                            p.y = r + cc;
                            p.x = c + cc;
                            tp = specledpos(cnt, sd, p, 6, cc, 3 * mincnt / 4);
                        }
                        if (isvalid(tp)) {
                            if (v > 1400) {
                                pv sp = bestledpos(cnt, sd, tp, 0, cc, GREEN, SQR(param[sd][LED_CNT]));
                                if (!isvalid(sp))
                                    sp = bestledpos(cnt, sd, tp, 0, cc, BLUE, SQR(param[sd][LED_CNT]));
                                ok = isvalid(sp);
                            } 
                            if (ok) {
                                DEBUG("ndot: %d %d %d:%d v %d bestcmpcnt %d", sd, dmark, (int) tp.x, (int) tp.y, v, cnt);
                                dots[ndot[sd]++] = tp;
                                if (testparam == 2) {
                                    if (dmark >= 5) {
                                        pos[sd][dmark] = tp;
                                        --dmark;
                                    }
                                }
                            }
/*                            else
                                nodots[nodot[sd]++] = tp;  */
                        }
                    }
                }
            }
            lastrgb = rgb;
        }
    }
}


void transformdots(int sd, pv left, pv center, pv right, pv bottom, pv top, 
 int leftelem, int centerelem, int rightelem, int bottomelem, int topelem)
{
    double ldis, rdis, bdis;

    if (!isvalid(right))
        right = add(center, sub(center, left));
    else if (!isvalid(left))
        left = add(center, sub(center, right));
    rdis = dist(right, center);
    ldis = dist(left, center);
    double lravg = (ldis + rdis) / 2;
    bool specshift = false;
    if (!isvalid(bottom)) {
        double dx = max(ldis, rdis) + (max(ldis, rdis) - min(ldis, rdis));
        double dd = dist(oripos[centerelem], oripos[bottomelem]) * dx / dist(oripos[leftelem], oripos[centerelem]);
        pv ip = add(center, rotl(sub(right, center)));
        pv bp1 = segpd(center, ip, dd);
        pv bp2 = center;
        bp2.y += dd;
        if (mode == SIM)
            bottom = bp1;
        else if (mode == LAB3)
            bottom = mid(bp1, bp2);
        else {// assume iss/lab/lab2 style 
            bottom = center;
            bottom.y += 0.6 * dd;
            if (left.y > right.y)
                bottom.x -= 25;
            else
                bottom.x += 25;
        }
    }
    bdis = dist(bottom, center);

    double brate, lrate, rrate, trate;
    brate = bdis / dist(oripos[centerelem], oripos[bottomelem]);
    lrate = ldis / dist(oripos[leftelem], oripos[centerelem]);
    rrate = rdis / dist(oripos[centerelem], oripos[rightelem]);
    trate = brate;

    pv obp = closest(oripos[centerelem], oripanel[3], oripanel[2]);
    double odb = dist(oripos[centerelem], obp);
    pv olp = closest(oripos[centerelem], oripanel[0], oripanel[3]);
    double odl = dist(oripos[centerelem], olp);
    pv orp = closest(oripos[centerelem], oripanel[1], oripanel[2]);
    double odr = dist(oripos[centerelem], orp);
    pv otp = closest(oripos[centerelem], oripanel[0], oripanel[1]);
    double odt = dist(oripos[centerelem], otp);

    pv actbp = segpd(center, bottom, brate * odb);
    pv actlp = segpd(center, left, lrate * odl);
    pv actrp = segpd(center, right, rrate * odr);
/*    panel[sd][3] = add(actbp, sub(actlp, center));
    panel[sd][2] = add(actbp, sub(actrp, center)); */
    double upmul = mode == LAB3 ? (bottom.y > 3 * h / 4 ? 1.05 : 1.02) : 1.2;

    panel[sd][3] = segpd(actlp, add(actlp, sub(actbp, center)), brate * odb * ldis / lravg) ;
    panel[sd][2] = segpd(actrp, add(actrp, sub(actbp, center)), brate * odb * rdis / lravg) ;
    if (mode != LAB3 && mode != SIM) {
        panel[sd][2].x -= 20;
        panel[sd][3].x += 20;
    }
    panel[sd][1] = segpd(panel[sd][2], actrp, trate * (odt + odb) * upmul * rdis / lravg);
    panel[sd][0] = segpd(panel[sd][3], actlp, trate * (odt + odb) * upmul * ldis / lravg); 
    npanelside[sd] = 4;
    dotpanel[sd] = true;
    pos[sd][0] = left;
    pos[sd][1] = bottom;
    pos[sd][2] = center;
    pos[sd][3] = right;
    if (bigledfound[sd])
        pos[sd][4] = ::center[sd];

    if (bigledfound[sd] && !inpoly0(::center[sd], panel[sd], 4) ||
        dist(panel[sd][2], panel[sd][3]) > w / 2 ||
        dist(panel[sd][0], panel[sd][1]) + dist(panel[sd][3], panel[sd][2]) >
        dist(panel[sd][0], panel[sd][3]) + dist(panel[sd][1], panel[sd][2])) {
        DEBUG("gotcha1");
        npanelside[sd] = 0;
        nuseddots[sd] = 0;
    }
    if (bigledfound[sd]) {
        pv olp = closest(::center[sd], panel[sd][0], panel[sd][3]);
        pv orp = closest(::center[sd], panel[sd][1], panel[sd][2]);
        if (dist(::center[sd], olp) > dist(::center[sd], orp)) {
            DEBUG("gotcha2");
            npanelside[sd] = 0;
            nuseddots[sd] = 0;
        }
    }
}

void dotstopanel(int sd)
{
    int cap = 8;
    FORD(i, ndot[sd] - 3, 0) {
        if (dots[i + 2].y - dots[i].y > 100 * h / 1200 || !almostonline(dots[i], dots[i + 1], dots[i + 2], cap))
            continue;
        if (dots[i].x > dots[i + 1].x)
            swap(dots[i], dots[i + 1]);
        if (dots[i].x > dots[i + 2].x)
            swap(dots[i], dots[i + 2]);
        if (dots[i + 1].x > dots[i + 2].x)
            swap(dots[i + 1], dots[i + 2]);
        double d1 = dist(dots[i], dots[i + 1]);
        double d2 = dist(dots[i + 1], dots[i + 2]);
        if (min(d1, d2) / max(d1, d2) < 0.7)
            continue;
        if (d1 < w / 20 || d1 > w / 6 || abs(dots[i].x - dots[i + 1].x) < 2 * abs(dots[i].y - dots[i + 1].y))
            continue;
        if (bigledfound[sd] && 
            (dots[i].y < ::center[sd].y || dots[i + 1].x < ::center[sd].x - 100 || dots[i].x > ::center[sd].x + 100))
            break;
        bool ok = true;
        FORD(k, i - 1, 0) {
            if (dist(dots[k], dots[i + 1]) < min(d1, d2)) {
                ok = false;
                break;
            }
        }
        if (!ok)
            continue;
        pv bottompv;
        nuseddots[sd] = 3;
        FORD(j, ndot[sd] - 1, i + 3) {
            double db = dist(dots[i + 1], dots[j]);
            if (db > 5 * d1 / 4 || d1 > 2 * db)
                continue;
            if (almostonline(dots[j], dots[i + 1], add(dots[i + 1], rotr(sub(dots[i + 1], dots[i]))), 2 * cap) ||
                abs(dots[i + 1].x - dots[j].x) < 30) {
                bottompv = dots[j];
                nuseddots[sd] = 4;
                break;
            }
        }
        transformdots(sd, dots[i], dots[i + 1], dots[i + 2], bottompv, pv(), 
            A03LED, A04LEDT, A05LED, A04LEDB, -1);
        return;
    }
    FORD(i, ndot[sd] - 2, 0) {
        pv left = dots[i];
        pv center = dots[i + 1];
        if (left.x > center.x)
            swap(left, center);
        double dd = dist(left, center);
        if (dd < w / 20 || dd > w / 6 || abs(dots[i].x - dots[i + 1].x) < 2 * abs(dots[i].y - dots[i + 1].y))
            continue;
        if (bigledfound[sd] && 
            (dots[i].y < ::center[sd].y || center.x < ::center[sd].x - 100 || left.x > ::center[sd].x + 100))
            break;
        bool ok = true;
        FORD(k, i - 1, 0) {
            if (abs(left.x - dots[k].x) < 30 && dist(left, dots[k]) < 3 * dd / 4) {
                ok = false;
                break;
            }
        }
        if (!ok)
            continue;
        FOR(j, i + 2, ndot[sd] - 1) {
            double db = dist(center, dots[j]);
            if (abs(center.x - dots[j].x) > 30 || db > 5 * dd / 4 || dd > 2 * db)
    //            !almostonline(center, add(center, rotr(sub(left, center))), dots[i + 2], 20))
                continue;
            bool ok = true;
            FOR(k, j + 1, ndot[sd] - 1) {
                if (abs(center.x - dots[k].x) < 30) {
                    ok = false;
                    break;
                }
            }
            if (!ok)
                break;
            nuseddots[sd] = 3;
            transformdots(sd, left, center, pv(), dots[j], pv(), 
                A03LED, A04LEDT, A05LED, A04LEDB, -1);
            return;
        }
    }
    FORD(i, ndot[sd] - 2, 0) {
        pv center = dots[i];
        pv right = dots[i + 1];
        if (center.x > right.x)
            swap(center, right);
        double dd = dist(right, center);
        if (dd < w / 20 || dd > w / 6 || abs(dots[i].x - dots[i + 1].x) < 2 * abs(dots[i].y - dots[i + 1].y))
            continue;
        if (bigledfound[sd] && 
            (dots[i].y < ::center[sd].y || center.x < ::center[sd].x - 100 || center.x > ::center[sd].x + 300))
            break;
        FOR(j, i + 2, ndot[sd] - 1) {
            double db = dist(center, dots[j]);
            if (abs(center.x - dots[j].x) > 30 || db > 5 * dd / 4 || dd > 2 * db)
    //            !almostonline(center, add(center, rotr(sub(center, right))), dots[i + 2], 20))
                continue;
            nuseddots[sd] = 3;
            transformdots(sd, pv(), center, right, dots[j], pv(), 
                A03LED, A04LEDT, A05LED, A04LEDB, -1);
            return;
        }
    }
    FORD(i, ndot[sd] - 2, 0) {
        pv center = mid(dots[i], dots[i + 1]);
        pv right = dots[i + 1];
        if (center.x > right.x)
            swap(center, right);
        double dd = dist(right, center);
        if (dd < w / 20 || dd > w / 6 || abs(dots[i].x - dots[i + 1].x) < 2 * abs(dots[i].y - dots[i + 1].y))
            continue;
        if (bigledfound[sd] && 
            (dots[i].y < ::center[sd].y || center.x < ::center[sd].x  - 100 || center.x > ::center[sd].x + 300))
            break;
        bool ok = true;
        FORD(k, i - 1, 0) {
            if (abs(center.x - dots[k].x) < 30 && dist(center, dots[k]) < 3 * dd / 4) {
                ok = false;
                break;
            }
        }
        if (!ok)
            continue;
        FOR(j, i + 2, ndot[sd] - 1) {
            double db = dist(center, dots[j]);
            if (abs(center.x - dots[j].x) > 30 || db > 5 * dd / 4 || dd > 2 * db)
    //            !almostonline(center, add(center, rotr(sub(center, right))), dots[i + 2], 20))
                continue;
            nuseddots[sd] = 3;
            transformdots(sd, pv(), center, right, dots[j], pv(), 
                A03LED, A04LEDT, A05LED, A04LEDB, -1);
            return;
        }
    }
}


void transformtwo(int sd, pv center, pv bottom)
{
    int centerelem = PPLED;
    int bottomelem = A01RLEDT;

    double rate = dist(center, bottom) / dist(oripos[centerelem], oripos[bottomelem]);
    pv obp = closest(oripos[centerelem], oripanel[3], oripanel[2]);
    double odb = dist(oripos[centerelem], obp);
    pv olp = closest(oripos[centerelem], oripanel[0], oripanel[3]);
    double odl = dist(oripos[centerelem], olp);
    pv orp = closest(oripos[centerelem], oripanel[1], oripanel[2]);
    double odr = dist(oripos[centerelem], orp);
    pv otp = closest(oripos[centerelem], oripanel[0], oripanel[1]);
    double odt = dist(oripos[centerelem], otp);

    pv actbp = segpd(center, bottom, rate * odb);
    pv sb = sub(bottom, center);
    pv actlp = segpd(center, add(center, rotl(sb)), (center.x < w / 2 ? 0.9 : 0.8) * rate * odl);
    pv actrp = segpd(center, add(center, rotr(sb)), 0.8 * rate * odr);
    panel[sd][3] = segpd(actlp, add(actlp, sub(actbp, center)), rate * odb) ;
    panel[sd][2] = segpd(actrp, add(actrp, sub(actbp, center)), rate * odb) ;
    panel[sd][1] = segpd(panel[sd][2], actrp, rate * (odt + odb));
    panel[sd][0] = segpd(panel[sd][3], actlp, rate * (odt + odb)); 
    npanelside[sd] = 4;
    dotpanel[sd] = true;
    pos[sd][0] = bottom;
    pos[sd][2] = center;
    if (bigledfound[sd])
        pos[sd][4] = ::center[sd];

    if (bigledfound[sd] && !inpoly0(::center[sd], panel[sd], 4) ||
        dist(panel[sd][2], panel[sd][3]) > w / 2 ||
        dist(panel[sd][0], panel[sd][1]) + dist(panel[sd][3], panel[sd][2]) >
        dist(panel[sd][0], panel[sd][3]) + dist(panel[sd][1], panel[sd][2])) {
            DEBUG("gotcha1a");
            npanelside[sd] = 0;
            nuseddots[sd] = 0;
    }
    if (bigledfound[sd]) {
        pv olp = closest(::center[sd], panel[sd][0], panel[sd][3]);
        pv orp = closest(::center[sd], panel[sd][1], panel[sd][2]);
        if (dist(::center[sd], olp) > dist(::center[sd], orp)) {
            DEBUG("gotcha2a");
            npanelside[sd] = 0;
            nuseddots[sd] = 0;
        }
    }

}

void dotstopanellast(int sd)
{
    if (ndot[sd] < 1)
        return;
    int topi = 0;
    bool topiok = false;
    FOR(i, 0, ndot[sd] - 1) {
        int cnt;
        pv p = bestledpos(cnt, sd, dots[i], 3, param[sd][LED_CIRCLER], GREEN, SQR(param[sd][LED_CNT]));
        if (isvalid(p)) {
            topi = i;
            topiok = true;
            break;
        }
    }
    if (bigledfound[sd] && dots[topi].x < ::center[sd].x)
        return;
    FOR(j, topi + 1, ndot[sd] - 1) {
        if (abs(dots[topi].x - dots[j].x) <= 50 && dist(dots[topi], dots[j]) < h / 5) {
            nuseddots[sd] = 2;
            transformtwo(sd, dots[topi], dots[j]);
            return;
        }
    }
    if (topiok) {
        FOR0(i, 5) {
            pv testp = dots[topi];
            testp.x += 10;
            testp.y += 60 + i * 50;
            int cnt;
            pv p = bestledpos(cnt, sd, testp, param[sd][LED_MAXXY], param[sd][LED_CIRCLER], BLACK, SQR(param[sd][LED_CNT]));
            if (isvalid(p)) {
                nuseddots[sd] = 1;
                transformtwo(sd, dots[topi], p);
                return;
            }
        } 
    }
}

void assesspanel(int sd)
{
    calcsums(sd);
    collectdots(sd);
    dotstopanel(sd);
}

void assessmini(int sd)
{
    ndot[sd] = 0;
    collectminidots(sd);
    dotstopanel(sd);
}

void assesslast(int sd)
{
    ndot[sd] = 0;
    collectdots(sd);
    dotstopanellast(sd);
}

void assesslastmini(int sd)
{
    ndot[sd] = 0;
    collectminidots(sd);
    dotstopanellast(sd);
}

void buildonbigled(int sd)
{
    if (!bigledfound[sd])
        return;
    dotpanel[sd] = true;
    nuseddots[sd] = 1;
    naddpanelcorner[sd] = 4;
    if (mode == LAB3) {
        pv bottom = center[sd];
        bottom.x += 250;
        bottom.y -= 40;
        pv top = bottom;
        top.y -= 190;
        transformtwo(sd, top, bottom);
    } else {
        panel[sd][3].x = ::center[sd].x - 250;
        panel[sd][3].y = ::center[sd].y + 350;
        panel[sd][2].x = ::center[sd].x + 290;
        panel[sd][2].y = ::center[sd].y + 260;
        panel[sd][0].x = ::center[sd].x - 190;
        panel[sd][0].y = ::center[sd].y - 280;
        panel[sd][1].x = ::center[sd].x + 460;
        panel[sd][1].y = ::center[sd].y - 340;
        npanelside[sd] = 4;
        if (bigledfound[sd])
            pos[sd][4] = ::center[sd];
    }
}

void calcpos()
{    
    LI start = gettimems();
    FOR0(sd, 2)
        if (mode != SIM)
            smooth(sd);
    FOR0(sd, 2)
        calcavg(sd);
    setmode();
    adjustparams();
    bool nothingelse = false;
    FOR0(sd, 2) {
//        lightup(sd);
        recolor(sd);
#if 1
        if (mode == LAB3 || mode == LAB2)
            center[sd] = findabigled(sd);
        if (!isvalid(center[sd])) {
            findredblock(sd);
            if (isvalid(redblock[sd][0])) {
                FOR0(i, 2)
                    pos[sd][i] = redblock[sd][i];
                center[sd] = guess_center_redblock(sd);
            }
            else
                center[sd] = findabigled(sd);
        }
#endif
    }
    FOR0(sd, 2) {
        if (!isvalid(center[sd]))  {
            if (isvalid(center[1 - sd])) {
                center[sd] = center[1-sd];
                if (sd == 0)
                    center[sd].x += leftrightadjust;
                else
                    center[sd].x -= leftrightadjust;
            } 
            else {
                center[sd].y = h / 2;
                if (mode == LAB2 || mode == LAB3)
                    center[sd].x = 5 * w / 8;
                else if (mode == LAB)
                    center[sd].x = w / 3;
                else
                    center[sd].x = w / 4;
                if (sd == 0)
                    center[sd].x += leftrightadjust;
                else
                    center[sd].x -= leftrightadjust;
            }
        }
        if (isvalid(center[sd]))
            pos[sd][2] = center[sd]; 
    }
#if 1 //!!!!!!!!!!!!!!
    FOR0(sd, 2) {
        if (isvalid(center[sd])) {
            bool ok = collectpanelside(sd, true, pleft[sd], pright[sd]);
            ok = collectpanelside(sd, false, ptop[sd], pbottom[sd]) && ok;
            if (mode == LAB2 && isvalid(pright[sd].p[0]) && isvalid(pleft[sd].p[0]) && sd == 1 && npanelside[0] == 4) {
                pv pr = closest(center[sd], pright[sd].p[0], pright[sd].p[1]);
                double dr = dist(center[sd], pr);
                pv pl = closest(center[sd], pleft[sd].p[0], pleft[sd].p[1]);
                double dl = dist(center[sd], pl);
                if (dr < dl * 0.8)
                    pright[sd].p[0] = pv(); // invalidate
            }
            if (isvalid(ptop[sd].p[0]) && isvalid(pleft[sd].p[0]))
                panel[sd][0] = intersect(ptop[sd].p[0], ptop[sd].p[1], pleft[sd].p[0], pleft[sd].p[1]);
            if (isvalid(ptop[sd].p[0]) && isvalid(pright[sd].p[0]))
                panel[sd][1] = intersect(ptop[sd].p[0], ptop[sd].p[1], pright[sd].p[0], pright[sd].p[1]);
            if (isvalid(pbottom[sd].p[0]) && isvalid(pright[sd].p[0]))
                panel[sd][2] = intersect(pbottom[sd].p[0], pbottom[sd].p[1], pright[sd].p[0], pright[sd].p[1]);
            if (isvalid(pbottom[sd].p[0]) && isvalid(pleft[sd].p[0]))
                panel[sd][3] = intersect(pbottom[sd].p[0], pbottom[sd].p[1], pleft[sd].p[0], pleft[sd].p[1]);
            paneldir[sd][0] = ptop[sd];
            paneldir[sd][1] = pright[sd];
            paneldir[sd][2] = pbottom[sd];
            paneldir[sd][3] = pleft[sd];
            blackpanelside[sd] = npanelside[sd];
        }
    }
    if (npanelside[0] < 4 && npanelside[1] < 4) {
        FOR0(sd, 2) {
            if (npanelside[sd] == 3) {
                FOR(dir, -1, 0) {
                    int bpdir = dir == 0 ? 1 : -1;
                    FOR0(i, 4) {
                        pv p0 = paneldir[sd][(i + 4 + dir) %4].p[0];
                        pv p1 = paneldir[sd][(i + 4 + dir) %4].p[1];
                        if (!isvalid(panel[sd][i]) && 
                            isvalid(panel[sd][(i + 4 + bpdir) % 4]) && isvalid(panel[sd][(i + 4 + 2 * bpdir) % 4]) &&
                            isvalid(p0) && isvalid(p1)) {
                            double xd =    dist(panel[sd][(i + 4 + bpdir) %4], panel[sd][(i + 4 + 2 * bpdir) %4]);
                            double xmul = i <= 1 ? 1.25 : 0.95;
                            double d = (((i + 1 + dir) & 1) ? (xmul * oriw / orih) : orih / oriw) * xd;
                            pv xp;
                            if (mode == SIM || mode == LAB) { // with good visibility we assume the panel partly off screen (otherwise we would have detected the edge)
                                if (i == 0 && dir == -1 || i == 1 && dir == 0)
                                    xp = intersect(pv(0.0, 0.0), pv(w - 1, 0.0), p0, p1);
                                else  if (i == 1 && dir == -1 || i == 2 && dir == 0)
                                    xp = intersect(pv(w - 1, 0.0), pv(w - 1, h - 1), p0, p1);
                                else if (i == 2 && dir == -1 || i == 3 && dir == 0)
                                    xp = intersect(pv(0.0, h - 1), pv(w - 1, h - 1), p0, p1);
                                else
                                    xp = intersect(pv(0.0, h - 1), pv(0.0, 0.0), p0, p1);
                                if (isvalid(xp)) {
                                    double specd = dist(panel[sd][(i + 4 + bpdir) %4], xp);
                                    if (specd > d)
                                        d = specd;
                                }
                            }
                            if ((i + 4 + dir) % 4 >= 2)
                                d = -d;
                            if (dir == 0)
                                d = -d;
                            pv pn = norm(sub(p1, p0));
                            panel[sd][i] = add(panel[sd][(i + 4 + bpdir) %4], mul(d, pn));
                        }
                    }
                }
                naddpanelcorner[sd] += 2;
                npanelside[sd]++;
                break;
            }
        }
    }
    FOR0(sd, 2) {
        if (npanelside[sd] < 4) {
            int cnt = 0;
            FOR(dir, -1, 0) {
                int bpdir = dir == 0 ? 1 : -1;
                FOR0(i, 4) {
                    if (!isvalid(panel[sd][i])) {
                        if (isvalid(panel[1 - sd][i]) && isvalid(panel[1 - sd][(i + 4 + bpdir) %4]) &&
                            isvalid(panel[sd][(i + 4 + bpdir) % 4]) && isvalid(paneldir[sd][(i + 4 + dir) %4].p[0])) {
                            double d = dist(panel[1 - sd][i], panel[1 - sd][(i + 4 + bpdir) %4]);
                            if ((i + 4 + dir) % 4 >= 2)
                                d = -d;
                            if (dir == 0)
                                d = -d;
                            pv pn = norm(sub(paneldir[sd][(i + 4 + dir) %4].p[1], paneldir[sd][(i + 4 + dir) %4].p[0]));
                            panel[sd][i] = add(panel[sd][(i + 4 + bpdir) %4], mul(d, pn));
                            naddpanelcorner[sd]++;
                            cnt++;
                        }
                    }
                }
            }
            npanelside[sd] += cnt / 2;
            if (!isvalid(panel[sd][0]))
                panel[sd][0] = add(panel[sd][1], mul(1.5, sub(panel[sd][3], panel[sd][2])));
        }
     }
#endif
#if 1
    FOR0(sd, 2)
        if (npanelside[sd] < 4)
            assesspanel(sd); 
    if (npanelside[0] < 4 && npanelside[1] < 4)
        FOR0(sd, 2)
            assessmini(sd); 
    if (npanelside[0] < 4 && npanelside[1] < 4)
        FOR0(sd, 2)
            assesslast(sd);  
    if (npanelside[0] < 4 && npanelside[1] < 4)
        FOR0(sd, 2)
            assesslastmini(sd);  
#endif
    if (npanelside[0] < 4 && npanelside[1] < 4)
        FOR0(sd, 2) 
            buildonbigled(sd);

    FOR0(sd, 2) {
        int cnt = 0;
        FOR0(i, 4) {
            if (!isvalid(panel[sd][i]) && isvalid(panel[1 - sd][i])) {
                double adjust;
                if (mode == LAB3 || mode == SIM) {
                    adjust = 340 * oneperoridist[1];
                    if (i == 0 || i == 3)
                        adjust *= dist(panel[1 - sd][0], panel[1 - sd][3]);
                    else
                        adjust *= dist(panel[1 - sd][1], panel[1 - sd][2]);
                    adjust = max(adjust, 200);
                } else {
                    adjust = 340 * h / 1200;
                    if (i <= 1)
                        adjust *= 1.2;
                }
                panel[sd][i] = panel[1 - sd][i];
                if (sd == 1)
                    panel[sd][i].x -= adjust;
                else
                    panel[sd][i].x += adjust;
                naddpanelcorner[sd]++;
                cnt++;
            }
        }
        if (cnt == 4)
            npanelside[sd] = 4;
        else if (cnt == 3)
            npanelside[sd] = 2;
    }
    FOR0(sd, 2)
        adjustpaneltop(sd);
    FOR0(sd, 2) 
        FOR0(i, 4)
            pos[sd][11 + i] = panel[sd][i];
    DEBUG("calcpos total: %d", (int) (gettimems() - start));
}

void adjustbigleds(int sd, pv oldpos, pv newpos, int cnt, int element)
{
//    if (newpos.x < 20 || newpos.x > w - 20 || newpos.y < 20 || newpos.y > h - 20)
//        return;
    if (!bigok[sd]) {
        DEBUG("adjustbig: %d %d:%d %d elem: %d cnt: %d", sd, (int) newpos.x, (int) newpos.y, bigrecheck[sd], element, cnt);
        pv dv = sub(newpos, oldpos);
        if (!ppledfix[sd] && (naddpanelcorner[sd] || dotpanel[sd])) { 
            param[sd][LED_MAXXY] /= 2;
            param[sd][LED_BLACKMAXXY] /= 2;
            param[sd][BIGLED_MAXXY] /= 2;
//            dv = sub(segp(oldpos, newpos, 0.5), oldpos);
        }
        int start = ppledfix[sd] ? (mode == SIM ? A01RLEDT : A02LEDA1) : PPS;
        FOR(i, start, MAX_ELEMENT - 1)
            pos[sd][i] = add(pos[sd][i], dv);
        bigok[sd] = true;
    }
}

void checkled(int sd, int element, bool skipblack = false, int color = GREEN, int toggle = -1)
{
    if (mode == SIM)
        skipblack = true; // no black, and darkgrey is too risky
    pv oldpos = pos[sd][element];
    int cnt;
    int expectcnt = param[sd][LED_CNT];
    if (element == A04LEDB)
        expectcnt = expectcnt / 2;
    pv p = bestledpos(cnt, sd, pos[sd][element], param[sd][LED_MAXXY], param[sd][LED_CIRCLER], color, SQR(expectcnt));
    if (!isvalid(p) && color == BLUE && bigok[sd])
        p = bestledpos(cnt, sd, pos[sd][element], param[sd][LED_MAXXY] / 2, param[sd][LED_CIRCLER], GREEN, SQR(expectcnt));
    if (isvalid(p)) {
        if (!bigok[sd] && (element == A03LED || element == A04LEDT) && p.y < pos[sd][element].y - 10 && 
            (mode == LAB3 || mode == LAB2)) {
            pv bp = bestledpos(cnt, sd, p, 10, 
                param[sd][BIGLED_CIRCLER], GREEN, SQR(param[sd][BIGLED_CNT]), param[sd][BIGLED_MINCIRCLER]);
            if (isvalid(bp))
                return; // not good!
        }
        state[element] = ON;
    }
    else if (!skipblack) {
        p = bestledpos(cnt, sd, pos[sd][element], param[sd][LED_BLACKMAXXY], param[sd][LED_BLACKCIRCLER], BLACK, SQR(param[sd][LED_BLACKCNT]));
        if (isvalid(p))  {
            if (!bigok[sd] && (element == A03LED || element == A04LEDT) && 
                p.y < pos[sd][element].y - 10 && (mode == LAB2 || mode == LAB3)) {
                pv bp = bestledpos(cnt, sd, p, 10, param[sd][BIGLED_CIRCLER], BLACK, SQR(param[sd][BIGLED_CNT] + 3));
                if (!isvalid(bp)) 
                    bp = bestledpos(cnt, sd, p, 10, 
                        param[sd][BIGLED_CIRCLER], GREEN, SQR(param[sd][BIGLED_CNT]), param[sd][BIGLED_MINCIRCLER]);
                if (!isvalid(bp)) 
                    bp = bestledpos(cnt, sd, p, 10, 
                        param[sd][BIGLED_CIRCLER], BLUE, SQR(param[sd][BIGLED_CNT]), param[sd][BIGLED_MINCIRCLER]);
                if (isvalid(bp))
                    return; // not good!
            }
            if (state[element] != ON)
                state[element] = OFF;
        }
    }
    if (isvalid(p)) {
        pos[sd][element] = p;
        if (!bigok[sd]) {
            if (element == A03LED || element == A04LEDT) 
                adjustbigleds(sd, oldpos, p, cnt, element);
        } else if (toggle >= 0) {
            pv dv = sub(p, oldpos);
            if (mode == LAB2 && toggle == A03T) { // dont know why
                dv.x /= 2;
                dv.y /= 2;
            }
            pos[sd][toggle] = add(pos[sd][toggle], dv);
            if (element == A04LEDT)
                pos[sd][A04LEDB] = add(pos[sd][A04LEDB], dv);
            if (element == A05LED && mode == SIM) {
                if (state[A01RLEDT] != ON)
                    pos[sd][A01RLEDT] = add(pos[sd][A01RLEDT], dv);
                if (state[A01RLEDB] != ON)
                    pos[sd][A01RLEDB] = add(pos[sd][A01RLEDB], dv);
            }
        }
    }
}

void checkbigled(int sd, int element, int color)
{
    int cnt;
    bool bluegreenspecial = false;
    pv p;
    if (settled[sd][element])
        return;
    int mxy = param[sd][BIGLED_MAXXY];
//    if (mode == LAB3)
//        mxy /= 2;
    if (!bigok[sd])
        mxy += mxy / 2;
    if (!bigrecheck[sd]) {
        if (worseside && naddpanelcorner[sd] >= 3 && state[element] != ON)
            return; // avoid being wrong
        if (!ppledfix[sd] && !bigok[sd] && element <= A02LEDA3)
            mxy += mxy / 2;
        pv xp = pos[sd][element];
        if (!bigok[sd] && (mode == LAB2 || mode == LAB3 || naddpanelcorner[sd] || dotpanel[sd])) {
            int dd = 10;
            if (naddpanelcorner[sd] || dotpanel[sd])
                dd = 25;
            if (element <= A02LEDA3)
                xp.y -= dd;
            else if (element >= A02LEDC1)
                xp.y += dd;
            if (element == A02LEDA1 || element == A02LEDB1 || element == A02LEDC1)
                xp.x -= dd;
            else if (element == A02LEDA3 || element == A02LEDB3 || element == A02LEDC3)
                xp.x += dd;
        }
        p = bestledpos(cnt, sd, xp, bigok[sd] && !dotpanel[sd] && !naddpanelcorner[sd] ? mxy / 2 : mxy, 
            param[sd][BIGLED_CIRCLER], color, SQR(param[sd][BIGLED_CNT]), param[sd][BIGLED_MINCIRCLER]); // danger of overlap
    }
    if (!isvalid(p)) {
        if (bigrecheck[sd] && !settled[sd][element]) {
            if (worseside && naddpanelcorner[sd] >= 3 && state[element] == ON)
                return; // avoid being wrong
            if (color == BLUE && bigok[sd]) {
                p = bestledpos(cnt, sd, pos[sd][element], mxy / 2, 
                    param[sd][BIGLED_CIRCLER], GREEN, SQR(param[sd][BIGLED_CNT]), param[sd][BIGLED_MINCIRCLER]);
                if (isvalid(p))
                    bluegreenspecial = true;
            }
            if (!isvalid(p)) {
                if (mode == LAB2 || mode == LAB3) 
                    p = bestledpos(cnt, sd, pos[sd][element], mxy, 
                        param[sd][BIGLED_CIRCLER], MIDGREY, SQR(param[sd][BIGLED_CNT]), param[sd][BIGLED_MINCIRCLER]);
                if (!isvalid(p) && mode != SIM)
                    p = bestledpos(cnt, sd, pos[sd][element], mxy, 
                        param[sd][BIGLED_CIRCLER], BLACK, SQR(param[sd][BIGLED_CNT]), param[sd][BIGLED_MINCIRCLER]);
                if (!isvalid(p)) {
                    int mul = mode == SIM ? 2 : 3;
                    p = bestledpos(cnt, sd, pos[sd][element], mxy, mul * param[sd][BIGLED_CIRCLER] / 2, MIDGREY, SQR(param[sd][BIGLED_GREYCNT]));
                }
            }
        } else if (mode == SIM)
            p = bestledpos(cnt, sd, pos[sd][element], mxy, param[sd][BIGLED_CIRCLER], WHITE, 20);
    }
    if (isvalid(p)) {
        if (!bigrecheck[sd] || bluegreenspecial) {
            if (worseside && state[element] != ON)
                return; // we cant believe it
            if (mode == SIM || cnt >= 60) {
                state[element] = ON;
                settled[sd][element] = true; //!!!! fura, nem pozicionaljuk????
            }
        }
        if (cnt >= 60) {
            adjustbigleds(sd, pos[sd][element], p, cnt, element);
            if (mode != SIM || state[element] == ON)
                pos[sd][element] = p;
            settled[sd][element] = true;
        }
    }
}

void recalcpanelavg()
{
    FOR0(sd, 2) {
        if (npanelside[sd] == 4) {
            memcpy(img[sd], saveimg[sd], sizeof(int) * h * w);
            lightup(sd);
            if (mode != SIM && mode != LAB)
                smooth(sd);
            minc[sd] = (int) max(min(panel[sd][0].x, panel[sd][3].x) - 200, 0);
            minr[sd] = (int) max(min(panel[sd][0].y, panel[sd][1].y) - 200, 0);
            maxc[sd] = (int) min(max(panel[sd][1].x, panel[sd][2].x) + 200, w - 1);
            maxr[sd] = (int) min(max(panel[sd][2].y, panel[sd][3].y) + 200, h - 1);
            int rtot, gtot, btot, cnt;
            rtot = gtot = btot = cnt = 0;
            int dy = (maxr[sd] - minr[sd]) / 40;
            int dx = (maxc[sd] - minc[sd]) / 40;
            for(int r = minr[sd]; r <= maxr[sd]; r += dy) {
                for(int c = minc[sd]; c <= maxc[sd]; c += dx) { 
                    pv p;
                    p.y = r;
                    p.x = c; 
                    if (isin((int) p.y, (int) p.x) && inpoly0(p, panel[sd], 4)) {
                        int rgb = getp(sd, (int) p.y, (int) p.x);
                        if (!iswhite(sd, rgb) && !isblack(sd, rgb)) { // skip iss hand
                            rtot += red(rgb);
                            gtot += green(rgb);
                            btot += blue(rgb);
                            cnt++;
                        }
                    }
                }
            }
            avg[sd] = (rtot + gtot + btot) / (cnt + 1) / 3;
            greyavg[sd] = avg[sd];
            blackavg[sd] = avg[sd];
            recolor(sd);
        }
    }
}

void checkleds()
{
    LI start = gettimems();
    int sd;
    FOR0(sdd, 2) {
        if (sdd == 0)
            sd = naddpanelcorner[0] == 0 ? 0 : 1;
        else {
            sd = 1 - sd;
            worseside = naddpanelcorner[sd] > 0 && naddpanelcorner[1 - sd] == 0;
/*            if (naddpanelcorner[sd] > 2) {
                FOR0(i, MAX_ELEMENT) {
                    pos[sd][i] = pos[1 - sd][i];
                    if (sd == 1)
                        pos[sd][i].x -= leftrightadjust;
                    else
                        pos[sd][i].x += leftrightadjust;
                }
            } */
        }
        if (!ppledfix[sd])
            checkled(sd, PPLED);
        if (naddpanelcorner[sd] || dotpanel[sd]) { 
            if (!ppledfix[sd]) {
                param[sd][LED_MAXXY] *= 3;
                param[sd][LED_BLACKMAXXY] *= 3;
                param[sd][BIGLED_MAXXY] *= 2;
            } else {
                param[sd][LED_MAXXY] *= 2;
                param[sd][LED_BLACKMAXXY] *= 2;
            }
            if (mode == LAB3 || mode == LAB2) { 
                if (!ppledfix[sd])
                    lightdowncheck[sd] = true;
            }
        }
        if (mode == LAB3 || mode == LAB2) { // cant rely on bigs with blue-green problems
            checkled(sd, A03LED, false, GREEN, A03T);
            checkled(sd, A04LEDT, false, GREEN, A04T);
            lightdowncheck[sd] = false;
        }
        while(1) {
            if (mode == LAB2 || mode == LAB3 || naddpanelcorner[sd] || dotpanel[sd]) { 
                checkbigled(sd, A02LEDC1, BLUE);
                checkbigled(sd, A02LEDC3, BLUE);
                checkbigled(sd, A02LEDA3, BLUE);
                checkbigled(sd, A02LEDA1, BLUE);
                if (!bigrecheck[sd] && bigok[sd]) {
                    pos[sd][A02LEDA2] = mid(pos[sd][A02LEDA1], pos[sd][A02LEDA3]);
                    pos[sd][A02LEDC2] = mid(pos[sd][A02LEDC1], pos[sd][A02LEDC3]);
                    pos[sd][A02LEDB1] = mid(pos[sd][A02LEDA1], pos[sd][A02LEDC1]);
                    pos[sd][A02LEDB3] = mid(pos[sd][A02LEDA3], pos[sd][A02LEDC3]);
                }
                checkbigled(sd, A02LEDA2, GREEN);
                checkbigled(sd, A02LEDB1, GREEN);
                checkbigled(sd, A02LEDB3, GREEN);
                checkbigled(sd, A02LEDC2, GREEN);
                if (!bigrecheck[sd] && bigok[sd])
                    pos[sd][A02LEDB2] = mid(pos[sd][A02LEDB1], pos[sd][A02LEDB3]);
                checkbigled(sd, A02LEDB2, BLUE);
            } else {
                checkbigled(sd, A02LEDA1, BLUE);
                checkbigled(sd, A02LEDA2, GREEN);
                checkbigled(sd, A02LEDA3, BLUE);
                checkbigled(sd, A02LEDB1, GREEN);
                checkbigled(sd, A02LEDB2, BLUE);
                checkbigled(sd, A02LEDB3, GREEN);
                checkbigled(sd, A02LEDC1, BLUE);
                checkbigled(sd, A02LEDC2, GREEN);
                checkbigled(sd, A02LEDC3, BLUE);
            }
            if (bigrecheck[sd])
                break;
            if (bigok[sd] && !naddpanelcorner[sd] && !dotpanel[sd])
                break;
            bigrecheck[sd] = true;
            checkled(sd, A03LED, false, GREEN, A03T);
            checkled(sd, A04LEDT, false, GREEN, A04T);
        }
        checkled(sd, A01RLEDT);
        if (!upperpos[sd] || mode == SIM)
            checkled(sd, A01RLEDB, mode != LAB3, BLUE);
        if (!bigrecheck[sd]) {
            checkled(sd, A03LED, false, GREEN, A03T);
            checkled(sd, A04LEDT, false, GREEN, A04T);
        }
        checkled(sd, A04LEDB, true, BLUE);
        checkled(sd, A05LED, false, GREEN, A05T);
    }
    DEBUG("checkleds total: %d", (int) (gettimems() - start));
}

void applyrules()
{
    if (state[PPLED] == OFF) {
        if (isvalid(redblock[0][0]) || isvalid(redblock[1][0])) { // otherwise it could be bad
            FOR0(i, MAX_ELEMENT) {
                switch (i) {
                case PPS: case PPC: case A01RS: case A03T: case A04T: case A05T:
                    break;
                default:
                    state[i] = OFF;
                    break;
                }
            }
        } else {
            int cnt = 0;
            FOR0(i, MAX_ELEMENT) {
                if (state[i] == ON) 
                    cnt++;
            }
            if (cnt >= 3)
                state[PPLED] = ON;
        }
    }
    if (state[PPLED] == UNKNOWN) {
        int cnt = 0;
        FOR0(i, MAX_ELEMENT) {
            if (state[i] == ON) 
                cnt++;
        }
        if (cnt >= 2 || cnt == 1 && mode == LAB3)
            state[PPLED] = ON;
    }
    FOR0(i, MAX_ELEMENT) {
        if (state[i] == UNKNOWN) {
            switch (i) {
            case PPS: case PPC: case A01RS: case A03T: case A04T: case A05T:
                state[i] = DOWN;
                break;
            default:
                state[i] = OFF;
                break;
            }
        }
    } 
    if (state[PPLED] == ON)
        state[PPS] = UP;
    if (state[A01RLEDT] == ON)
        state[A01RS] = UP;
    else if (state[A01RLEDB] == ON)
        state[A01RS] = DOWN;
    else
        state[A01RS] = CENTER;
    if (state[A03LED] == ON)
        state[A03T] = UP;
    if (state[A04LEDT] == ON)
        state[A04T] = UP;
    else if (state[A04LEDB] == ON)
        state[A04T] = DOWN;
    else
        state[A04T] = CENTER;
    if (state[A05LED] == ON)
        state[A05T] = UP; 

}

void overridewithori()
{
    FOR0(sd, 2) {
        FOR0(i, MAX_ELEMENT) {
            pos[sd][i] = oripos[i];
        }
    }
}

bool isredblockleft(int sd)
{
    if (mode == SIM)
        return redblock[sd][0].x < w / 2;
    if (redblock[sd][0].x < w / 4)
        return true;
    if (redblock[sd][0].x > 3 * w / 4)
        return false;
    return panel[sd][3].y < panel[sd][2].y;
}

void fixspecials()
{
    FOR0(sd, 2) {
        if (isvalid(redblock[sd][0])) {
            double rh = redblock[sd][1].y - redblock[sd][0].y;
            double rw = redblock[sd][1].x - redblock[sd][0].x;
            double rate = dist(panel[sd][0], panel[sd][1]) * oneperoridist[0];
            if (rh > 2 * rw || redblock[sd][1].y >= pos[sd][PPS].y + rate * 50 ||
                rh > 3 * rw / 2 && redblock[sd][1].y >= pos[sd][PPS].y + rate * 30)
                state[PPC] = DOWN;
            else if (rw > rh)
                state[PPC] = UP;
        }
    }
    if (state[PPC] == UNKNOWN) {
        if (mode == LAB3 || mode == SIM)
            state[PPC] = UP;
        else {
            state[PPC] = UP;
            FOR0(sd, 2) {
                if (!isvalid(panel[sd][0]) || !isvalid(panel[sd][1]))
                    continue;
                double rate = dist(panel[sd][0], panel[sd][1]) * oneperoridist[0];
                int cnt;
                pv xp = pos[sd][PPLED];
                if (isredblockleft(sd)) 
                    xp.x -= rate * 30;
                xp.y += rate * 70;
                pv p = bestledpos(cnt, sd, xp, param[sd][A01RS_MAXXY], param[sd][A01RS_CIRCLER], RED, SQR(param[sd][A01RS_CNT]), 0, true);
                if (!isvalid(p))
                    p = bestledpos(cnt, sd, xp, (int) (rate * 10), param[sd][A01RS_CIRCLER], BLACK, SQR(param[sd][A01RS_CNT]), 0, true);
                if (isvalid(p)) {
                    state[PPC] = DOWN;
                    break;
                }
            }
        }
    }
    FOR0(sd, 2) {
        if (isvalid(redblock[sd][0])) {
            if (state[PPC] == DOWN) {
                pos[sd][PPC].y = redblock[sd][1].y;
                pos[sd][PPC].x = (redblock[sd][0].x + redblock[sd][1].x) / 2;
            } else {
                if (upperpos[sd]) {
                    double rate = dist(panel[sd][0], panel[sd][1]) * oneperoridist[0];
                    if (isredblockleft(sd)) {
                        pos[sd][PPC].y = minred[sd] + 10;
                        pos[sd][PPC].x = redblock[sd][0].x + rate * 40 / 2;
                    } else {
                        pos[sd][PPC].y = maxred[sd] + 10;
                        pos[sd][PPC].x = redblock[sd][1].x - rate * 40 / 2;
                    }
                } else {
                    pos[sd][PPC].y = (3 * redblock[sd][0].y + redblock[sd][1].y) / 4;
                    pos[sd][PPC].x = (redblock[sd][0].x + redblock[sd][1].x) / 2;
                }
            }
        } else {
            if (state[PPC] == DOWN) {
                double lrate = dist(panel[sd][0], panel[sd][1]) * oneperoridist[0];
                double rrate = dist(panel[sd][1], panel[sd][2]) * oneperoridist[1];
                pos[sd][PPC].y = pos[sd][PPLED].y + rrate * 70 + sd * 15;
                if (!isredblockleft(sd))
                    pos[sd][PPC].x = pos[sd][PPLED].x - 20 + (1 - sd) * 35;
                else
                    pos[sd][PPC].x = pos[sd][PPLED].x - lrate * 40;
            } else if (mode != LAB3) {
                    double rate = dist(panel[sd][0], panel[sd][1]) * oneperoridist[0];
                    pos[sd][PPC].y = pos[sd][PPLED].y - 20;
                    if (!isredblockleft(sd))
                        pos[sd][PPC].x = pos[sd][PPLED].x + rate * 20;
                    else
                        pos[sd][PPC].x = pos[sd][PPLED].x - rate * 60;
            }
        }
        pos[sd][A01RS].y += 10;
        int cnt;
        pv p = bestledpos(cnt, sd, pos[sd][A01RS], param[sd][A01RS_MAXXY], param[sd][A01RS_CIRCLER], BLACK, SQR(param[sd][A01RS_CNT]), 0, true);
        if (isvalid(p)) {
            pos[sd][A01RS] = p;
            if (upperpos[sd])
                pos[sd][A01RS].y += 10;
        }
    }
}

void adjustpanelparams()
{
    //!!! poziciofuggo, lefele csokken a meret (felso nezopontbol)
    FOR0(sd, 2) {
        if (npanelside[sd] == 4) {
            int minx = (int) min(panel[sd][0].x, panel[sd][3].x);
            int miny = (int) min(panel[sd][0].y, panel[sd][1].y);
            int maxx = (int) max(panel[sd][1].x, panel[sd][2].x);
            int maxy = (int) max(panel[sd][2].y, panel[sd][3].y);
            double avgratio = ((maxx - minx) * oneperoridist[0] * 1600 / w + (maxy - miny) * oneperoridist[1] * 1200 / h) / 2;
            FOR0(i, MAX_PARAM)
                if (panelparam[i])
                    param[sd][i] = (int) (avgratio * param[sd][i]);
        }
    }
}

#ifdef LOCAL

void savetga(string sfn, int sd)
{
    unsigned char buf[4 * MAX_W];
    ofstream ofs(sfn.c_str(), ios::out | ios::binary);
    ofs.put(0);
    ofs.put(0);
    ofs.put(2);
    FOR0(i, 9)
        ofs.put(0);         
    ofs.put(w & 0xFF);
    ofs.put(w / 256);
    ofs.put(h & 0xFF);
    ofs.put(h / 256);
    ofs.put(32);
    ofs.put(0);
    FORD(y, h -1, 0) {
        unsigned char *p = buf;
        FOR0(x, w) {
            int rgb = getp(sd, y, x);
            *p++ = blue(rgb);
            *p++ = green(rgb);
            *p++ = red(rgb);
            *p++ = 0;
        }
        ofs.write((char *) buf, 4 * w);
    }
}

void saveimgs(string post)
{
    FOR0(sd, 2) {
        string sfn = testcase;
        if (sd == 0)
            sfn += "_left";
        else
            sfn += "_right";
        sfn += "_";
        sfn += post + ".tga";
        savetga(sfn, sd);
    }
}
#endif

VS recognize()
{
    if (!globalstart)
        globalstart = gettimems();
    VS res;
    init();
    filldefault();
    LI saveimgtime = 0;
    try {
        calcpos();
#ifdef LOCAL
        LI preimgsavetime = gettimems();
        if (testcommand & 8)
            saveimgs("calc");
        saveimgtime = gettimems() - preimgsavetime;
#endif
        if (!testparam) {
            transformall();
            adjustpanelparams();
            recalcpanelavg();
            adjustforppled();
            checkleds();
            fixspecials();
        }
        applyrules(); 
    }
    catch (...) {
        DEBUG("exception");
        filldefault();
    }
    LI globalend = gettimems();
    FOR0(sd, 2)
        DEBUG("result:,%s,%s,%d,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%5.2f", testcase, testmode, sd, modestr[mode], davg[sd], avg[sd], bigledfound[sd], isvalid(redblock[sd][0]), dotpanel[sd], (int) (globalend - globalstart - saveimgtime), npanelside[sd], blackpanelside[sd], naddpanelcorner[sd], nuseddots[sd], (double) hw / baseblack[0]);
    int threshold = 10;
    if (mode == ISS || mode == LAB)
        threshold = 15;

    FOR0(sd, 2) {
        FOR0(i, MAX_ELEMENT) {
            if (state[i] == UNKNOWN) {
                switch (i) {
                case PPS: case PPC: case A01RS: case A03T: case A04T: case A05T:
                    state[i] = DOWN;
                    break;
                default:
                    state[i] = OFF;
                    break;
                }
            }
            int x =    (int) pos[sd][i].x;
            int y =    (int) pos[sd][i].y;
            if (x < 0)
                pos[sd][i].x = threshold / 3;
            if (y < 0)
                pos[sd][i].y = threshold / 3;
            if (x >= w)
                pos[sd][i].x = w - threshold / 3;
            if (y >= h)
                pos[sd][i].y = h - threshold / 3;
        }
    }

    FOR0(i, MAX_ELEMENT) {
        char buf[100];
        sprintf(buf, "%s,%d,%d,%d,%d", statestr[state[i]].c_str(), (int) pos[0][i].x, (int) pos[0][i].y, (int) pos[1][i].x, (int) pos[1][i].y);
        res.PB(buf);
    }
    return res;
}

class RobonautEye {
public:
    VS recognizeObjects(VI& leftEyeImage, VI& rightEyeImage) 
    {
        globalstart = gettimems();
        h = leftEyeImage[0];
        w = leftEyeImage[1];
        FOR0(r, h)
            FOR0(c, w)
                setp(0, r, c, leftEyeImage[2 + r * w + c]);
        FOR0(r, h)
            FOR0(c, w)
                setp(1, r, c, rightEyeImage[2 + r * w + c]);
        DEBUG("array handle time %d", (int) (gettimems() - globalstart));
        return recognize();
    }
};

#ifdef LOCAL

VS split(string s, char c)
{
    VS res;
    char *p = (char *) s.c_str();
    do {
        char *bp = p;
        while(*p && *p != c)
            p++;
        res.PB(string(bp, p));
    } while (*p++);
    return res;
}

void score()
{
    pv realpos[2][MAX_ELEMENT];
    string realstate[MAX_ELEMENT];
    double score[MAX_ELEMENT][3];
    ifstream ifs("model.csv");
    string s;
    string selector(testcase);
    bool ok = false;
    ifs >> s;
    do {
        ifs >> s;
        if (s.empty())
            break;
        VS val = split(s, ',');
        string sel = val[0] + val[1];
        sel = sel.substr(0, SZ(sel) - 4);
        if (sel == selector) {
            FOR0(i, MAX_ELEMENT) {
                realstate[i] = val[2 + i * 5];
                realpos[0][i].x = atoi(val[2 + i * 5 + 1].c_str());
                realpos[0][i].y = atoi(val[2 + i * 5 + 2].c_str());
                realpos[1][i].x = atoi(val[2 + i * 5 + 3].c_str());
                realpos[1][i].y = atoi(val[2 + i * 5 + 4].c_str());
            }
            ok = true;
            break;
        }
    } while(!ifs.eof());
    if (!ok) {
        DEBUG("testcase %s not found in model.csv", selector.c_str());
        return;
    }
    int threshold;
    if (selector[0] == 'I' || (selector[0] == 'L' && selector[3] != '2' && selector[3] != '3'))
        threshold = 15;
    else
        threshold = 10;
    FOR0(i, MAX_ELEMENT) {
        FOR0(j, 3)
            score[i][j] = 0;
        if (realstate[i] != "HIDDEN") {
            FOR0(sd, 2) {
                if (realpos[sd][i].x == -1 && realpos[sd][i].y == -1) {
                    score[i][sd] = 1;
                } else {
                    int X1 = (int) realpos[sd][i].x, X2 = (int) pos[sd][i].x;
                    int Y1 = (int) realpos[sd][i].y, Y2 = (int) pos[sd][i].y;
                    double distanceLeft = sqrt((double) (X1 - X2) * (X1 - X2) + (Y1 - Y2) * (Y1 - Y2));
                    if (distanceLeft <= 2 * threshold) {
                        if (distanceLeft <= threshold) {
                            score[i][sd] = 1;
                        } else {
                            score[i][sd] = 1 - pow((distanceLeft - threshold) / threshold, 1.5);
                        }
                    }
                }
            }
            if (realstate[i] == statestr[state[i]])
                score[i][2] = 1;
        } else {
            FOR0(j, 3)
                score[i][j] = 1;
        }
    }

    double tot = 0;
    DEBUG("Scores for each object (left image location, right image location, state):");
    FOR0(i, MAX_ELEMENT) {
        FOR0(j, 3) {
            tot += score[i][j];
        }
        DEBUG("%s %g %g %g", elementstr[i].c_str(), score[i][0], score[i][1], score[i][2]);
    }

    tot /= MAX_ELEMENT * 3;
    DEBUG("Score : %f", tot);
}

int main(int argc, char *argv[])
{
    srand(1);
    if (argc >= 2)
        testcase = argv[1];
    if (argc >= 3)
        testparam = atoi(argv[2]);
    else
        testparam = 0;
    if (argc >= 4)
        testcommand = atoi(argv[3]);
    else
        testcommand = 0;
    if (argc >= 5)
        testmode = argv[4];
    debugf = stderr;

    VS res;
    string fn = "gen_";
    fn += testcase;
    fn += ".dat";
    if (testcommand & 1) {
        ifstream fgen(fn.c_str(), ios::binary);
        fgen.read((char *)&h, sizeof(h));
        fgen.read((char *)&w, sizeof(w));
        FOR0(sd, 2)
            fgen.read((char *)img[sd], h * w * sizeof(int));
        res = recognize();
    } else {
        string bs = "blank.dat";
        if (testcommand & 2)
            bs = fn;
        ofstream fgen(fn.c_str(), ios::binary);
        int size;
        VI vi[2];
        cin >> size;
        FOR0(sd, 2) {
            cin >> h >> w;
            if (testcommand & 2 && sd == 0) {
                fgen.write((char *)&h, sizeof(h));
                fgen.write((char *)&w, sizeof(w));
            }
            vi[sd].PB(h);
            vi[sd].PB(w);
            FOR0(r, h) {
                FOR0(c, w) {
                    int v;
                    cin >> v;
                    vi[sd].PB(v);
                    if (testcommand & 2)
                        setp(sd, r, c, v);
                }
            }
            if (testcommand & 2)
                fgen.write((char *)img[sd], h * w * sizeof(int));
            DEBUG("img %d done", sd);
        }
//        if (testcommand == 1) {
            RobonautEye re;
            res = re.recognizeObjects(vi[0], vi[1]);
//        }
    }
    if (testcommand & 4)
        saveimgs("final");
    score();
    DEBUG("out"); // keep this for now
    FOR0(i, MAX_ELEMENT) {
        if (SZ(res) <= i)
            cout << endl;
        else {
            DEBUG("%s %s", elementstr[i].c_str(), res[i].c_str());
            cout << res[i] << endl;
        }
    }
    cout.flush();
}
    
#endif
