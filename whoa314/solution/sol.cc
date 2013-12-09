#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

typedef long long int64;
typedef unsigned long long uint64;
typedef int int32;
typedef unsigned int uint32;

using namespace std;

class RobonautEye {
 public:
  vector<string> recognizeObjects(
      const vector<int>& leftEyeImage, const vector<int>& rightEyeImage);
};

static int getHue(int R, int G, int B) {
  if (R >= G && G >= B) {
    if (R - B == 0) return 0;
    return 60 * (G - B) / (R - B);
  } else if (G > R && R >= B) {
    return 60 * (2 - (R - B) / (G - B));
  } else if (G >= B && B > R) {
    return 60 * (2 + (B - R) / (G - R));
  } else if (B > G && G > R) {
    return 60 * (4 - (G - R) / (B - R));
  } else if (B > R && R >= G) {
    return 60 * (4 + (R - G) / (B - G));
  } else if (R >= B && B > G) {
    return 60 * (6 - (B - G) / (R - G));
  } else {
    return 0;
  }
}

struct SumImage {
  vector<int> sumR;
  vector<int> sumG;
  vector<int> sumB;
  vector<int64> sumSqR;
  vector<int64> sumSqG;
  vector<int64> sumSqB;
  vector<bool> isForeground;
  vector<bool> isBackground;
  vector<bool> isRed;
  vector<bool> isGreen;
  vector<bool> isBlue;
  SumImage(const vector<int>& arr) {
    int width = arr[1];
    int height = arr[0];
    sumR.resize(width * height);
    sumG.resize(width * height);
    sumB.resize(width * height);
    sumSqR.resize(width * height);
    sumSqG.resize(width * height);
    sumSqB.resize(width * height);
    isForeground.resize(width * height);
    isBackground.resize(width * height);
    isRed.resize(width * height);
    isGreen.resize(width * height);
    isBlue.resize(width * height);
    int* ptrR = &sumR[0];
    int* ptrG = &sumG[0];
    int* ptrB = &sumB[0];
    int64* ptrSqR = &sumSqR[0];
    int64* ptrSqG = &sumSqG[0];
    int64* ptrSqB = &sumSqB[0];
    const int* ptrCur = &arr[2];
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        int rgb = *ptrCur++;
        int r = (rgb >> 16) & 0x00ff;
        int g = (rgb >> 8) & 0x00ff;
        int b = (rgb) & 0x00ff;
        *ptrR = r;
        *ptrG = g;
        *ptrB = b;
        *ptrSqR = r * r;
        *ptrSqG = g * g;
        *ptrSqB = b * b;
        if (y - 1 >= 0) {
          *ptrR += *(ptrR - width);
          *ptrG += *(ptrG - width);
          *ptrB += *(ptrB - width);
          *ptrSqR += *(ptrSqR - width);
          *ptrSqG += *(ptrSqG - width);
          *ptrSqB += *(ptrSqB - width);
        }
        if (y - 1 >= 0 && x - 1 >= 0) {
          *ptrR -= *(ptrR - width - 1);
          *ptrG -= *(ptrG - width - 1);
          *ptrB -= *(ptrB - width - 1);
          *ptrSqR -= *(ptrSqR - width - 1);
          *ptrSqG -= *(ptrSqG - width - 1);
          *ptrSqB -= *(ptrSqB - width - 1);
        }
        if (x - 1 >= 0) {
          *ptrR += *(ptrR - 1);
          *ptrG += *(ptrG - 1);
          *ptrB += *(ptrB - 1);
          *ptrSqR += *(ptrSqR - 1);
          *ptrSqG += *(ptrSqG - 1);
          *ptrSqB += *(ptrSqB - 1);
        }
        ++ptrR;
        ++ptrG;
        ++ptrB;
        ++ptrSqR;
        ++ptrSqG;
        ++ptrSqB;
      }
    }
  }
};

struct SccProperty {
  int xmin, xmax, ymin, ymax;
  int xminy, xmaxy;
  int yminx, ymaxx;
  int x, y;  // average.
  int size;
  int r, g, b;
  int propNum;
  bool isGreenLed;
  bool isBlueLed;
  bool isRockerSwitch;
  bool containsCenter;
  bool isFirstFromLeft;
};

struct ModelObject {
  bool leftFound, rightFound;
  string state;
  int leftX, leftY;
  int rightX, rightY;
  SccProperty* leftProp;
  SccProperty* rightProp;
  ModelObject() {
    leftFound = false;
    rightFound = false;
    state = "HIDDEN";
    leftX = 0;
    leftY = 0;
    rightX = 0;
    rightY = 0;
    leftProp = NULL;
    rightProp = NULL;
  }
};

struct ModelSummary {
  ModelObject objects[22];
  int width, height;
  ModelSummary(int w, int h) : width(w), height(h) {
  }
};

static double compRatios9[22][2] = {
    {0.02,0.15},{0.07,0.19},{0.00,0.00},{1.00,0.00},
    {0.71,-0.01},{1.22,0.02},{0.78,1.02},{0.77,0.81},
    {0.75,0.59},{1.00,1.00},{0.99,0.78},{0.98,0.57},
    {1.22,0.98},{1.21,0.77},{1.20,0.56},{1.76,0.89},
    {1.57,0.90},{1.73,0.42},{1.54,0.42},{1.91,0.41},
    {1.70,-0.05},{1.51,-0.06}};
static double compRatios10[22][2] = {
    {0.02,0.20},{0.06,0.24},{0.00,0.00},{1.00,0.00},
    {0.72,-0.01},{1.22,0.03},{0.79,1.30},{0.78,1.03},
    {0.75,0.75},{1.01,1.27},{1.00,1.00},{0.99,0.73},
    {1.23,1.24},{1.22,0.98},{1.20,0.71},{1.78,1.13},
    {1.58,1.15},{1.74,0.53},{1.55,0.54},{1.92,0.52},
    {1.71,-0.07},{1.51,-0.07}};

// Newer one, small only.
static double compRatiosOrthoSmallLeft[22][2] = {
    {-0.01,-0.17},{-0.21,-0.10},{0.00,0.00},{1.00,0.00},
    {0.75,-0.02},{1.22,-0.04},{0.71,-1.21},{0.71,-0.95},
    {0.72,-0.69},{0.96,-1.20},{0.97,-0.94},{0.97,-0.69},
    {1.21,-1.19},{1.21,-0.94},{1.21,-0.68},{1.84,-1.12},
    {1.62,-1.13},{1.84,-0.55},{1.62,-0.55},{2.07,-0.55},
    {1.84,0.00},{1.63,0.01}};
static double compRatiosOrthoSmallRight[22][2] = {
    {0.03,-0.19},{-0.09,-0.10},{0.00,0.00},{1.00,0.00},
    {0.73,-0.01},{1.22,-0.01},{0.88,-1.30},{0.85,-1.02},
    {0.81,-0.74},{1.12,-1.29},{1.09,-1.00},{1.05,-0.73},
    {1.36,-1.27},{1.31,-0.99},{1.28,-0.72},{1.92,-1.16},
    {1.72,-1.19},{1.83,-0.57},{1.64,-0.58},{2.02,-0.57},
    {1.75,0.01},{1.55,0.02}};
static double compRatiosOrthoLargeLeft[22][2] = {
    {-0.02,-0.19},{0.15,-0.31},{0.00,0.00},{1.00,0.00},
    {0.71,0.02},{1.22,-0.03},{0.55,-1.27},{0.59,-1.01},
    {0.61,-0.74},{0.77,-1.24},{0.80,-0.97},{0.84,-0.71},
    {0.98,-1.20},{1.01,-0.94},{1.04,-0.69},{1.49,-1.07},
    {1.31,-1.09},{1.56,-0.49},{1.39,-0.50},{1.73,-0.47},
    {1.66,0.10},{1.48,0.10}};
static double compRatiosOrthoLargeRight[22][2] = {
    {-0.04,-0.18},{0.09,-0.39},{0.00,0.00},{1.00,0.00},
    {0.72,0.04},{1.22,-0.03},{0.39,-1.13},{0.46,-0.89},
    {0.52,-0.65},{0.61,-1.09},{0.68,-0.86},{0.75,-0.62},
    {0.83,-1.05},{0.89,-0.82},{0.96,-0.59},{1.35,-0.92},
    {1.17,-0.94},{1.52,-0.39},{1.33,-0.41},{1.62,-0.28},
    {1.70,0.14},{1.51,0.14}};
static double compRatiosOrthoLeft[22][2] = {
    {-0.01,-0.18},{-0.02,-0.21},{0.00,0.00},{1.00,0.00},
    {0.73,0.00},{1.22,-0.04},{0.63,-1.24},{0.65,-0.98},
    {0.66,-0.72},{0.86,-1.22},{0.88,-0.96},{0.90,-0.70},
    {1.09,-1.19},{1.11,-0.94},{1.13,-0.69},{1.66,-1.09},
    {1.46,-1.11},{1.70,-0.52},{1.50,-0.53},{1.88,-0.51},
    {1.75,0.05},{1.55,0.06}};
static double compRatiosOrthoRight[22][2] = {
    {0.01,-0.18},{-0.02,-0.21},{0.00,0.00},{1.00,0.00},
    {0.73,0.01},{1.22,-0.01},{0.69,-1.24},{0.70,-0.97},
    {0.70,-0.71},{0.92,-1.21},{0.93,-0.95},{0.94,-0.69},
    {1.16,-1.19},{1.15,-0.93},{1.16,-0.67},{1.70,-1.07},
    {1.51,-1.09},{1.71,-0.50},{1.52,-0.51},{1.89,-0.48},
    {1.73,0.06},{1.53,0.06}};

static bool canReach(int width, int height,
    int x0, int y0, int x1, int y1, vector<vector<int> >& scc,
    map<int, SccProperty>& sccPropLarge) {
    if (x0 < 0) x0 = 0;
    if (x0 >= width) x0 = width - 1;
    if (y0 < 0) y0 = 0;
    if (y0 >= height) y0 = height - 1;
    if (x1 < 0) x1 = 0;
    if (x1 >= width) x1 = width - 1;
    if (y1 < 0) y1 = 0;
    if (y1 >= height) y1 = height - 1;
    int dx = x1 - x0;
    int dy = y1 - y0;
  int dist = (int) sqrt(dx * dx + dy * dy);
  if (dist == 0) return true;
  for (int i = 0; i < dist; ++i) {
    int x = dx * i / dist + x0;
    int y = dy * i / dist + y0;
    if (sccPropLarge.find(scc[x][y]) != sccPropLarge.end()) {
      return false;
    }
  }
  return true;
}

// Reachability as defined by the straight line connecting the two not intersecting a "large" component.
static bool canReach(int width, int height,
    SccProperty& prop0, SccProperty& prop1, vector<vector<int> >& scc,
    map<int, SccProperty>& sccPropLarge) {
  if (prop0.x < 0 || prop0.x >= width ||
      prop1.x < 0 || prop1.x >= width ||
      prop0.y < 0 || prop0.y >= height ||
      prop1.y < 0 || prop1.y >= height) {
    cerr << "Invalid prop out of bounds!! " << prop0.x << "," << prop0.y << " " <<
        prop1.x << "," << prop1.y << endl;
    return false;
  }
  return canReach(width, height, prop0.x, prop0.y, prop1.x, prop1.y, scc, sccPropLarge);
}

static bool mightBeRocker(
    int width, int height, vector<vector<int> >& scc, int ymin, SccProperty& prop) {
  if (prop.xmax - prop.xmin < width / 15 &&
      prop.ymax - prop.ymin < height / 10 &&
      max(prop.r, max(prop.g, prop.b)) < 100 &&
      //prop.ymax - prop.ymin > prop.xmax - prop.xmin &&
      prop.y > ymin) {
    //maxSize = prop.size;
    //maxProp = prop;
    int numInside = 0;
    int count = 0;
    for (int dy = max(0, prop.y - 20); dy < min(height, prop.y + 20); ++dy) {
      for (int dx = max(0, prop.x - 10); dx < min(width, prop.x + 10); ++dx) {
        ++count;
        if (scc[dx][dy] == prop.propNum) {
          ++numInside;
        }
      }
    }
    double fillFraction = numInside / (double)count;
    if (fillFraction > .8) {
      cerr <<"Candidate rocker switch with fill fraction: " << fillFraction << endl;
      return true;
    }
  }
  return false;
}

static double UPPER_NIBLACK = 1.5;

static inline int getAvgInt(const vector<int>& sumImg, int x, int y, int radius, int width, int height) {
  // Inclusive.
  int xmin = x - radius;
  int ymin = y - radius;
  if (xmin < 0) xmin = 0;
  if (ymin < 0) ymin = 0;

  // Exclusive.
  int xmax = x + radius;
  int ymax = y + radius;
  if (xmax > width) xmax = width;
  if (ymax > height) ymax = height;

  int64 tally = sumImg[(xmax - 1) + width * (ymax - 1)];
  if (ymin - 1 >= 0) {
    // Top.
    tally -= sumImg[(xmax - 1) + width * (ymin - 1)];
  }
  if (xmin - 1 >= 0) {
    // Left.
    tally -= sumImg[(xmin - 1) + width * (ymax - 1)];
  }
  if (ymin - 1 >= 0 && xmin - 1 >= 0) {
    // Top left.
    tally += sumImg[(xmin - 1) + width * (ymin - 1)];
  }
  tally /= (xmax - xmin) * (ymax - ymin);
  return (int) tally;
}

static inline int getAvg(const vector<int64>& sumImg, int x, int y, int radius, int width, int height) {
  // Inclusive.
  int xmin = x - radius;
  int ymin = y - radius;
  if (xmin < 0) xmin = 0;
  if (ymin < 0) ymin = 0;

  // Exclusive.
  int xmax = x + radius;
  int ymax = y + radius;
  if (xmax > width) xmax = width;
  if (ymax > height) ymax = height;

  int64 tally = sumImg[(xmax - 1) + width * (ymax - 1)];
  if (ymin - 1 >= 0) {
    // Top.
    tally -= sumImg[(xmax - 1) + width * (ymin - 1)];
  }
  if (xmin - 1 >= 0) {
    // Left.
    tally -= sumImg[(xmin - 1) + width * (ymax - 1)];
  }
  if (ymin - 1 >= 0 && xmin - 1 >= 0) {
    // Top left.
    tally += sumImg[(xmin - 1) + width * (ymin - 1)];
  }
  tally /= (xmax - xmin) * (ymax - ymin);
  return (int) tally;
}

static void doNiblackMulti(const vector<int>& arr, SumImage& sumImage, int kernel) {
  int width = arr[1];
  int height = arr[0];
  bool isSmall = (width < 1700 && height < 1700);
  double UPPER_LIM = UPPER_NIBLACK * UPPER_NIBLACK;
  double LOWER_LIM = -.75 * .75;
  const int* curPtr = &arr[2];
  int fullKernelSize = 4 * kernel * kernel;
  for (int y = 0, cur = 0; y < height; ++y) {
    int avgR = getAvgInt(sumImage.sumR, 0, y, kernel, width, height);
    int avgG = getAvgInt(sumImage.sumG, 0, y, kernel, width, height);
    int avgB = getAvgInt(sumImage.sumB, 0, y, kernel, width, height);
    int avgSqR = getAvg(sumImage.sumSqR, 0, y, kernel, width, height);
    int avgSqG = getAvg(sumImage.sumSqG, 0, y, kernel, width, height);
    int avgSqB = getAvg(sumImage.sumSqB, 0, y, kernel, width, height);
    int sigmaSquaredR = avgSqR - avgR * avgR;
    int sigmaSquaredG = avgSqG - avgG * avgG;
    int sigmaSquaredB = avgSqB - avgB * avgB;
    for (int x = 0; x < width; ++x, ++cur) {
      if (isSmall || ((x + y) & 0x01) == 0) {
        if (x > kernel && y > kernel && x < width - kernel && y < height - kernel) {
          int topY = width * (y + kernel - 1);
          int botY = width * (y - kernel - 1);
          avgR = (int)((sumImage.sumR[x + kernel - 1 + topY] -
                 sumImage.sumR[x + kernel - 1 + botY] -
                 sumImage.sumR[x - kernel - 1 + topY] +
                 sumImage.sumR[x - kernel - 1 + botY]) /
                  fullKernelSize);
          avgG = (int)((sumImage.sumG[x + kernel - 1 + topY] -
                 sumImage.sumG[x + kernel - 1 + botY] -
                 sumImage.sumG[x - kernel - 1 + topY] +
                 sumImage.sumG[x - kernel - 1 + botY]) /
                  fullKernelSize);
          avgB = (int)((sumImage.sumB[x + kernel - 1 + topY] -
                 sumImage.sumB[x + kernel - 1 + botY] -
                 sumImage.sumB[x - kernel - 1 + topY] +
                 sumImage.sumB[x - kernel - 1 + botY]) /
                  fullKernelSize);
          avgSqR = (int)((sumImage.sumSqR[x + kernel - 1 + topY] -
                 sumImage.sumSqR[x + kernel - 1 + botY] -
                 sumImage.sumSqR[x - kernel - 1 + topY] +
                 sumImage.sumSqR[x - kernel - 1 + botY]) /
                  fullKernelSize);
          avgSqG = (int)((sumImage.sumSqG[x + kernel - 1 + topY] -
                 sumImage.sumSqG[x + kernel - 1 + botY] -
                 sumImage.sumSqG[x - kernel - 1 + topY] +
                 sumImage.sumSqG[x - kernel - 1 + botY]) /
                  fullKernelSize);
          avgSqB = (int)((sumImage.sumSqB[x + kernel - 1 + topY] -
                 sumImage.sumSqB[x + kernel - 1 + botY] -
                 sumImage.sumSqB[x - kernel - 1 + topY] +
                 sumImage.sumSqB[x - kernel - 1 + botY]) /
                  fullKernelSize);
        } else {
          avgR = getAvgInt(sumImage.sumR, x, y, kernel, width, height);
          avgG = getAvgInt(sumImage.sumG, x, y, kernel, width, height);
          avgB = getAvgInt(sumImage.sumB, x, y, kernel, width, height);
          avgSqR = getAvg(sumImage.sumSqR, x, y, kernel, width, height);
          avgSqG = getAvg(sumImage.sumSqG, x, y, kernel, width, height);
          avgSqB = getAvg(sumImage.sumSqB, x, y, kernel, width, height);
        }
        sigmaSquaredR = avgSqR - avgR * avgR;
        sigmaSquaredG = avgSqG - avgG * avgG;
        sigmaSquaredB = avgSqB - avgB * avgB;
      }
      // 0 is b, 1 is g, 2 is r.
      int minChan = 999;
      double minDev = 999;
      int maxChan = -1;
      double maxDev = -1;
      int rgb = *curPtr++;
      int R = (rgb >> 16) & 0x00ff;
      int G = (rgb >> 8) & 0x00ff;
      int B = (rgb) & 0x00ff;
      int dR = R - avgR;
      dR *= abs(R - avgR);
      int dG = G - avgG;
      dG *= abs(G - avgG);
      int dB = B - avgB;
      dB *= abs(B - avgB);
      double devR = dR / (double) sigmaSquaredR;
      double devG = dG / (double) sigmaSquaredG;
      double devB = dB / (double) sigmaSquaredB;
      if (devR <= devG && devR <= devB) {
        minChan = 2;
        minDev = devR;
      } else if (devG <= devR && devG <= devB) {
        minChan = 1;
        minDev = devG;
      } else {
        minChan = 0;
        minDev = devB;
      }
      if (devR >= devG && devR >= devB) {
        maxChan = 2;
        maxDev = devR;
      } else if (devG >= devR && devG >= devB) {
        maxChan = 1;
        maxDev = devG;
      } else {
        maxChan = 0;
        maxDev = devB;
      }
      if (abs(maxDev) > abs(minDev) && maxDev > UPPER_LIM) {
        sumImage.isForeground[cur] = true;
        if (maxChan == 2) {
          sumImage.isRed[cur] = true;
        } else if (maxChan == 1) {
          sumImage.isGreen[cur] = true;
        } else {
          sumImage.isBlue[cur] = true;
        }

      } else if (abs(maxDev) < abs(minDev) && minDev < LOWER_LIM) {
        sumImage.isBackground[cur] = true;
        if (maxChan == 2) {
          sumImage.isRed[cur] = true;
        } else if (maxChan == 1) {
          sumImage.isGreen[cur] = true;
        } else {
          sumImage.isBlue[cur] = true;
        }
      }
      /*int avgCol = (avgR + avgG + avgB) / 3;

      int col = (R + G + B) / 3;
      int dR = R - avgR;
      int dG = G - avgG;
      int dB = B - avgB;
      double devR = sigmaSquaredR > 0 ? dR * dR / (double) sigmaSquaredR : 0;
      double devG = sigmaSquaredG > 0 ? dG * dG / (double) sigmaSquaredG : 0;
      double devB = sigmaSquaredB > 0 ? dB * dB / (double) sigmaSquaredB : 0;*/
      /*double avgDev = devR + devG + devB;
      avgDev /= 3;
      avgDev = sqrt(avgDev);
      if (col > avgCol) {
        if (avgDev > 1) {
          sumImage.isForeground[cur] = true;
        }
      } else if (col < avgCol) {
        if (avgDev > .75) {
          sumImage.isBackground[cur] = true;
        }
      }*/
    }
  }
}

static int64 getAverageContrast(const vector<int>& arr) {
  int width = arr[1];
  int height = arr[0];
  int count = 0;
  int64 diffSum = 0;
  int64 maxDiff = 0;
  for (int y = 1; y < height - 1; ++y) {
    for (int x = 1; x < width - 1; ++x) {
      int cur = 2 + y * width + x;
      int rs = (arr[cur] >> 16) & 0x00ff;
      int gs = (arr[cur] >> 8) & 0x00ff;
      int bs = (arr[cur]) & 0x00ff;
      for (int dy = y - 1; dy <= y + 1; ++dy) {
        for (int dx = x - 1; dx <= x + 1; ++dx) {
          if (dy == y && dx == x) continue;
          int dcur = 2 + dy * width + dx;
          int drs = (arr[dcur] >> 16) & 0x00ff;
          int dgs = (arr[dcur] >> 8) & 0x00ff;
          int dbs = (arr[dcur]) & 0x00ff;

          // Piecewise diff.
          int64 dr = drs - rs;
          int64 dg = dgs - gs;
          int64 db = dbs - bs;
          dr *= dr;
          dg *= dg;
          db *= db;
          int64 diff = dr + dg + db;

          diffSum += diff;
          if (diff > maxDiff) {
            maxDiff = diff;
          }
          ++count;
        }
      }
    }
  }
  diffSum /= count;
  cerr << "Average contrast " << diffSum << endl;
  cerr << "Max contrast " << maxDiff << endl;
  return diffSum;
}
static int64 getAverageLum(const vector<int>& arr) {
  int width = arr[1];
  int height = arr[0];
  int64 totalLum = 0;
  for (int y = 0, cur = 2; y < height; ++y) {
    for (int x = 0; x < width; ++x, ++cur) {
      int r = (arr[cur] >> 16) & 0x00ff;
      int g = (arr[cur] >> 8) & 0x00ff;
      int b = (arr[cur]) & 0x00ff;
      totalLum += r * r + g * g + b * b;
    }
  }
  totalLum /= width * height;
  return totalLum;
}

// Merge the two into prop0.
static void mergeScc(SccProperty& prop0, SccProperty& prop1) {
  cerr << "Merging " << prop1.x << "," << prop1.y << " into " << prop0.x << "," << prop0.y << endl;
  prop0.x = (int)((((int64)prop0.x) * prop0.size + ((int64)prop1.x) * prop1.size) /
      (prop0.size + prop1.size));
  prop0.y = (int)((((int64)prop0.y) * prop0.size + ((int64)prop1.y) * prop1.size) /
      (prop0.size + prop1.size));
  prop0.r = (int)((((int64)prop0.r) * prop0.size + ((int64)prop1.r) * prop1.size) /
      (prop0.size + prop1.size));
  prop0.g = (int)((((int64)prop0.g) * prop0.size + ((int64)prop1.g) * prop1.size) /
      (prop0.size + prop1.size));
  prop0.b = (int)((((int64)prop0.b) * prop0.size + ((int64)prop1.b) * prop1.size) /
      (prop0.size + prop1.size));
  if (prop1.xmin < prop0.xmin) {
    prop0.xmin = prop1.xmin;
    prop0.xminy = prop1.xminy;
  }
  if (prop1.xmax > prop0.xmax) {
    prop0.xmax = prop1.xmax;
    prop0.xmaxy = prop1.xmaxy;
  }
  if (prop1.ymin < prop0.ymin) {
    prop0.ymin = prop1.ymin;
    prop0.yminx = prop1.yminx;
  }
  if (prop1.ymax > prop0.ymax) {
    prop0.ymax = prop1.ymax;
    prop0.ymaxx = prop1.ymaxx;
  }
  prop0.isGreenLed = prop0.isGreenLed || prop1.isGreenLed;
  prop0.isBlueLed = prop0.isBlueLed || prop1.isBlueLed;
  prop0.isRockerSwitch = prop0.isRockerSwitch || prop1.isRockerSwitch;
  prop0.containsCenter = prop0.containsCenter || prop1.containsCenter;
  prop0.isFirstFromLeft = prop0.isFirstFromLeft || prop1.isFirstFromLeft;
  prop0.size += prop1.size;
}

// Find the strongly connected components.
static void findScc(
    const vector<int>& arr, vector<vector<int> >& scc,
    map<int, SccProperty>& sccProp,
    map<int, SccProperty>& sccPropLarge,
    SumImage& sumImage) {
  int width = arr[1];
  int height = arr[0];
  int maxComponentSize = width * height / 50;
  int minComponentSize = width * height / 45000;

  //int64 diffSum = getAverageContrast(arr);
  //int64 averageLum = getAverageLum(arr);
  //cerr << "Average lum " << averageLum << endl;

  int curComponent = 1;
  int numInterestingComponents = 0;
  vector<int> greenLedKeys;
  vector<SccProperty> greenLedValues;
  vector<int> blueLedKeys;
  vector<SccProperty> blueLedValues;
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      if (scc[x][y] != 0) {
        continue;
      }
      if (!sumImage.isForeground[x + width * y] && !sumImage.isBackground[x + width * y]) continue;
      bool localIsForeground = sumImage.isForeground[x + width * y];
      // Start a new fill.
      deque<pair<int, int> > q;
      scc[x][y] = curComponent;
      q.push_back(make_pair(x, y));
      int idx = 2 + y * width + x;
      int r = (arr[idx] >> 16) & 0x00ff;
      int g = (arr[idx] >> 8) & 0x00ff;
      int b = (arr[idx]) & 0x00ff;
      int localHue = getHue(r, g, b);
      int64 localLum = r * r + g * g + b * b;
      int componentSize = 1;
      int compXmin = width + 1;
      int compXmax = 0;
      int compYmin = height + 1;
      int compYmax = 0;
      int compXminY = 0;
      int compXmaxY = 0;
      int compYminX = 0;
      int compYmaxX = 0;
      int64 sumSquaredR = r * r;
      int64 sumSquaredG = g * g;
      int64 sumSquaredB = b * b;
      int64 compR = r;
      int64 compG = g;
      int64 compB = b;
      int64 compX = x;
      int64 compY = y;
      while (!q.empty()) {
        pair<int, int> cur = q.front();
        q.pop_front();
        if (cur.first > compXmax) {
          compXmax = cur.first;
          compXmaxY = cur.second;
        }
        if (cur.first < compXmin) {
          compXmin = cur.first;
          compXminY = cur.second;
        }
        if (cur.second > compYmax) {
          compYmax = cur.second;
          compYmaxX = cur.first;
        }
        if (cur.second < compYmin) {
          compYmin = cur.second;
          compYminX = cur.first;
        }
        /*int tmpidx = 2 + cur.second * width + cur.first;
        int tmpr = (arr[tmpidx] >> 16) & 0x00ff;
        int tmpg = (arr[tmpidx] >> 8) & 0x00ff;
        int tmpb = (arr[tmpidx]) & 0x00ff;
        int neighHue = getHue(tmpr, tmpg, tmpb);*/
        int dyStart = max(cur.second - 1, 0);
        int dyEnd = min(height, cur.second + 2);
        int dxStart = max(cur.first - 1, 0);
        int dxEnd = min(width, cur.first + 2);
        for (int dy = dyStart; dy < dyEnd; ++dy) {
          for (int dx = dxStart; dx < dxEnd; ++dx) {
            if (scc[dx][dy] != 0) {
              continue;
            }
            if (!sumImage.isForeground[dx + width * dy] && !sumImage.isBackground[dx + width * dy]) continue;
            if (localIsForeground != sumImage.isForeground[dx + width * dy]) continue;
            // Unvisited neighbor.
            int didx = 2 + dy * width + dx;
            int64 dr = (arr[didx] >> 16) & 0x00ff;
            int64 dg = (arr[didx] >> 8) & 0x00ff;
            int64 db = (arr[didx]) & 0x00ff;
            //int curHue = getHue((int)dr, (int)dg, (int)db);
            int64 dr0 = dr;
            int64 dg0 = dg;
            int64 db0 = db;
            /*dr -= r;
            dg -= g;
            db -= b;
            dr *= dr;
            dg *= dg;
            db *= db;
            int64 diff = dr + dg + db;*/
            //int hueDiff = abs(curHue - neighHue);
            //if (hueDiff > 180) {
            //  hueDiff = 360 - hueDiff;
            //}
            //if (hueDiff < 150) {
            //if (abs(curHue - localHue) < 60) {
            //if (diff < max(localLum, 10000LL)) {
              sumSquaredR += dr0 * dr0;
              sumSquaredG += dg0 * dg0;
              sumSquaredB += db0 * db0;
              compR += dr0;
              compG += dg0;
              compB += db0;
              compX += dx;
              compY += dy;
              scc[dx][dy] = curComponent;
              ++componentSize;
              q.push_back(make_pair(dx, dy));
            //}
          }
        }
      }  // while
      if (componentSize < maxComponentSize) {
        if (componentSize > 20) {
          sumSquaredR /= componentSize;
          sumSquaredG /= componentSize;
          sumSquaredB /= componentSize;
          bool isGreenLed =
              (sumSquaredG - min(sumSquaredR, sumSquaredB) > 14500) &&
              (sumSquaredG - max(sumSquaredR, sumSquaredB) > 2500) && localIsForeground;
          bool isBlueLed = componentSize > 30 && componentSize < height * 3 / 2 &&
              (sumSquaredR < 11000 || sumSquaredB - sumSquaredR > 43000) &&
              sumSquaredB > 40000 &&
              sumSquaredB - max(sumSquaredR, sumSquaredG) > 5000 && localIsForeground;
          if (isGreenLed || isBlueLed || componentSize > minComponentSize) {
            SccProperty prop;
            prop.xmin = compXmin;
            prop.xmax = compXmax;
            prop.ymin = compYmin;
            prop.ymax = compYmax;
            prop.xminy = compXminY;
            prop.xmaxy = compXmaxY;
            prop.yminx = compYminX;
            prop.ymaxx = compYmaxX;
            prop.size = componentSize;
            prop.r = (int) (compR / componentSize);
            prop.g = (int) (compG / componentSize);
            prop.b = (int) (compB / componentSize);
            prop.x = (int) (compX / componentSize);
            prop.y = (int) (compY / componentSize);
            prop.containsCenter = true;
            prop.isFirstFromLeft = false;
            prop.propNum = curComponent;
            if (isGreenLed) {
              prop.isGreenLed = true;
              prop.isBlueLed = false;
              greenLedKeys.push_back(curComponent);
              greenLedValues.push_back(prop);
              cerr << prop.x << "," << prop.y << " (g): " << componentSize << " " <<
                  sumSquaredR << "/" << sumSquaredG << "/" << sumSquaredB << endl;
            } else if (isBlueLed) {
              prop.isGreenLed = false;
              prop.isBlueLed = true;
              blueLedKeys.push_back(curComponent);
              blueLedValues.push_back(prop);
              cerr << prop.x << "," << prop.y << " (b): " << componentSize << " " <<
                  sumSquaredR << "/" << sumSquaredG << "/" << sumSquaredB << endl;
            } else {
              prop.isGreenLed = false;
              prop.isBlueLed = false;
              sccProp[curComponent] = prop;
            }
            ++numInterestingComponents;
          }
        }
      } else {
        // Large component
        SccProperty prop;
        prop.xmin = compXmin;
        prop.xmax = compXmax;
        prop.ymin = compYmin;
        prop.ymax = compYmax;
        prop.xminy = compXminY;
        prop.xmaxy = compXmaxY;
        prop.yminx = compYminX;
        prop.ymaxx = compYmaxX;
        prop.size = componentSize;
        prop.r = (int) (compR / componentSize);
        prop.g = (int) (compG / componentSize);
        prop.b = (int) (compB / componentSize);
        prop.x = (int) (compX / componentSize);
        prop.y = (int) (compY / componentSize);
        prop.containsCenter = true;
        prop.isFirstFromLeft = false;
        prop.propNum = curComponent;
        sccPropLarge[curComponent] = prop;
      }
      ++curComponent;
    }
  }
  int64 maxDist = (width + height) / 2 * 25LL / 1000L;
  maxDist *= maxDist;  // Compare in the squared space to avoid sqrt.
  for (int i = 0; i < greenLedKeys.size(); ++i) {
    if (sccProp.find(greenLedKeys[i]) != sccProp.end()) continue;
    SccProperty& prop0 = greenLedValues[i];
    vector<int> mergedComponents;
    mergedComponents.push_back(greenLedKeys[i]);
    for (int j = i + 1; j < greenLedKeys.size(); ++j) {
      if (sccProp.find(greenLedKeys[j]) != sccProp.end()) continue;
      SccProperty prop1 = greenLedValues[j];
      int dx = prop0.x - prop1.x;
      int dy = prop0.y - prop1.y;
      if (dx * dx + dy * dy < maxDist) {
        mergeScc(prop0, prop1);
        mergedComponents.push_back(greenLedKeys[j]);
      }
    }
    for (int k = 0; k < mergedComponents.size(); ++k) {
      sccProp[mergedComponents[k]] = prop0;
    }
  }
  for (int i = 0; i < blueLedKeys.size(); ++i) {
    if (sccProp.find(blueLedKeys[i]) != sccProp.end()) continue;
    SccProperty& prop0 = blueLedValues[i];
    vector<int> mergedComponents;
    mergedComponents.push_back(blueLedKeys[i]);
    for (int j = i + 1; j < blueLedKeys.size(); ++j) {
      if (sccProp.find(blueLedKeys[j]) != sccProp.end()) continue;
      SccProperty prop1 = blueLedValues[j];
      int dx = prop0.x - prop1.x;
      int dy = prop0.y - prop1.y;
      if (dx * dx + dy * dy < maxDist) {
        mergeScc(prop0, prop1);
        mergedComponents.push_back(blueLedKeys[j]);
      }
    }
    for (int k = 0; k < mergedComponents.size(); ++k) {
      sccProp[mergedComponents[k]] = prop0;
    }
  }

  // Now merge *everything*.
  /*int64 maxDistForMerge = (width + height) / 2 * 25L / 1000L;
  maxDistForMerge *= maxDistForMerge;
  vector<SccProperty*> allGridScc;
  set<int> allGridSccNum;
  for (map<int, SccProperty>::iterator it = sccProp.begin(); it != sccProp.end(); ++it) {
    SccProperty& prop = it->second;
    if (prop.isGreenLed || prop.isBlueLed) continue;
    if (allGridSccNum.find(prop.propNum) == allGridSccNum.end()) {
      allGridSccNum.insert(prop.propNum);
      allGridScc.push_back(&prop);
    }
  }
  vector<SccProperty*> mergedGrid;
  bool hadMerge = true;
  while (hadMerge) {
    set<int> alreadyMergedGrid;
    hadMerge = false;
    for (int i = 0; i < allGridScc.size(); ++i) {
      SccProperty& prop0 = *allGridScc[i];
      if (alreadyMergedGrid.find(prop0.propNum) != alreadyMergedGrid.end()) continue;
      vector<int> mergedComponents;
      alreadyMergedGrid.insert(prop0.propNum);
      mergedComponents.push_back(prop0.propNum);
      for (int j = i + 1; j < allGridScc.size(); ++j) {
        SccProperty& prop1 = *allGridScc[j];
        if (alreadyMergedGrid.find(prop1.propNum) != alreadyMergedGrid.end()) continue;
        int dx = prop0.x - prop1.x;
        int dy = prop0.y - prop1.y;
        if (dx * dx + dy * dy < maxDistForMerge) {
          mergeScc(prop0, prop1);
          alreadyMergedGrid.insert(prop1.propNum);
          mergedComponents.push_back(prop1.propNum);
          hadMerge = true;
        }
      }
      for (int k = 0; k < mergedComponents.size(); ++k) {
        sccProp[mergedComponents[k]] = prop0;
      }
      mergedGrid.push_back(&prop0);
    }
    if (hadMerge) {
      allGridScc = mergedGrid;
      mergedGrid.clear();
    }
  }*/

  for (map<int, SccProperty>::iterator it = sccProp.begin(); it != sccProp.end(); ++it) {
    SccProperty& prop = it->second;
    if (mightBeRocker(width, height, scc, 0, prop)) {
      // Might be a rocker, we'll say it is if it can reach all LEDs.
      bool reachAll = true;
      for (map<int, SccProperty>::iterator it2 = sccProp.begin(); it2 != sccProp.end(); ++it2) {
        SccProperty& prop2 = it2->second;
        if (prop2.isGreenLed || prop2.isBlueLed) {
          if (!canReach(width, height, prop, prop2, scc, sccPropLarge)) {
            reachAll = false;
            break;
          }
        }
      }
      if (reachAll) {
        cerr << "Candidate rocker post-reachability: " << prop.x << "," << prop.y << endl;
        prop.isRockerSwitch = true;
      }
    }
  }
  cerr << "Num components: " << (curComponent - 1) << endl;
  cerr << "Num interesting components: " << numInterestingComponents << endl;
  cerr << "Num large components: " << sccPropLarge.size() << endl;
}

static void findRed(const vector<int>& arr, ModelSummary& model, bool isLeft,
    vector<vector<int> >& scc, map<int, SccProperty>& sccProp, map<int, SccProperty>& sccPropLarge,
    int xmin, int ymin, int xmax, int ymax, double requiredFraction) {
  int width = arr[1];
  int height = arr[0];
  xmin = max(0, xmin);
  ymin = max(0, ymin);
  xmax = min(xmax, width);
  ymax = min(ymax, height);
  int numRed = 0;
  vector<int> redSum(width * height);
  int* redSumPtr = &redSum[0];
  vector<bool> isRed(width * height);
  for (int y = 0, cur = 2; y < height; ++y) {
    for (int x = 0; x < width; ++x, ++cur, ++redSumPtr) {
      if (x < xmin || x > xmax || y < ymin || y > ymax) continue;
      if (sccPropLarge.find(scc[x][y]) != sccPropLarge.end()) {
        if (y - 1 >= 0) {
          *redSumPtr += *(redSumPtr - width);
        }
        if (x - 1 >= 0) {
          *redSumPtr += *(redSumPtr - 1);
        }
        if (y - 1 >= 0 && x - 1 >= 0) {
          *redSumPtr -= *(redSumPtr - width - 1);
        }
        continue;
      }
      int64 r = (arr[cur] >> 16) & 0x00ff;
      int64 g = (arr[cur] >> 8) & 0x00ff;
      int64 b = (arr[cur]) & 0x00ff;
      r *= r;
      r *= r;
      g *= g;
      g *= g;
      b *= b;
      b *= b;
      if ((r > 0 && max(g, b) == 0) || r / (double)max(g, b) >= 8) {
        ++numRed;
        isRed[y * width + x] = true;
        *redSumPtr = 1;
      }
      if (y - 1 >= 0) {
        *redSumPtr += *(redSumPtr - width);
      }
      if (x - 1 >= 0) {
        *redSumPtr += *(redSumPtr - 1);
      }
      if (y - 1 >= 0 && x - 1 >= 0) {
        *redSumPtr -= *(redSumPtr - width - 1);
      }
    }
  }
  cerr << "Num red points: " << numRed << endl;
  int kerX = width / 80;
  int kerY = height / 40;
  double samp = (kerX * 2 + 1) * (kerY * 2 + 1);
  double maxFrac = 0;
  int bestX = 0;
  int bestY = 0;
  for (int x = kerX + xmin + 1; x < xmax - kerX - 1; ++x) {
    for (int y = kerY + ymin + 1; y < ymax - kerY - 1; ++y) {
      if (x < xmin || x > xmax || y < ymin || y > ymax) continue;
      if (!isRed[y * width + x]) continue;
      int count = redSum[x + kerX + width * (y + kerY)] -
                  redSum[x + kerX + width * (y - kerY - 1)] -
                  redSum[x - kerX - 1 + width * (y + kerY)] +
                  redSum[x - kerX - 1 + width * (y - kerY - 1)];
      double redFrac = count / samp;
      if (redFrac > maxFrac) {
        maxFrac = redFrac;
        bestX = x;
        bestY = y;
      }
    }
  }
  if (maxFrac > requiredFraction) {
    cerr << "Found red cover in range: " << xmin << "," << ymin << "," << xmax << "," << ymax << endl;
    map<int, int> propCount;
    for (int x = bestX - kerX; x <= bestX + kerX; ++x) {
      for (int y = bestY - kerY; y <= bestY + kerY; ++y) {
        if (sccProp.find(scc[x][y]) != sccProp.end()) {
          propCount[scc[x][y]]++;
        }
      }
    }
    SccProperty* bestProp = NULL;
    int bestPropValue = -1;
    for (map<int, int>::iterator it = propCount.begin(); it != propCount.end(); ++it) {
      if (it->second > bestPropValue && sccProp.find(it->first) != sccProp.end()) {
        bestPropValue = it->second;
        bestProp = &sccProp[it->first];
      }
    }
    if (bestProp != NULL) {
      cerr << "Found SCC corresponding to cover with size " << bestProp->size << endl;
      bestX = bestProp->x;
      bestY = bestProp->y;
    } else {
      cerr << "No SCC found corresponding to cover, just using center of kernel." << endl;
    }
    if (isLeft) {
      model.objects[1].leftX = bestX;
      model.objects[1].leftY = bestY;
      model.objects[1].leftFound = true;
      model.objects[1].leftProp = bestProp;
    } else {
      model.objects[1].rightX = bestX;
      model.objects[1].rightY = bestY;
      model.objects[1].rightFound = true;
      model.objects[1].rightProp = bestProp;
    }
  } else {
    cerr << "Failed to find red cover in range: " << xmin << "," << ymin << "," << xmax << ","
         << ymax << endl;
  }
}

static void findPanelLed(const vector<int>& arr,
    vector<vector<int> >& scc, map<int, SccProperty>& sccProp, map<int, SccProperty>& sccPropLarge,
    SumImage& sumImage,
    ModelSummary& model, bool isLeft) {
  int width = arr[1];
  int height = arr[0];
  int64 xpos = 0;
  int64 ypos = 0;
  int count = 0;
  int coverX, coverY;
  int disallowedProp = -1;
  if (isLeft) {
    coverX = model.objects[1].leftX;
    coverY = model.objects[1].leftY;
    if (model.objects[1].leftProp != NULL) {
      disallowedProp = model.objects[1].leftProp->propNum;
    }
  } else {
    coverX = model.objects[1].rightX;
    coverY = model.objects[1].rightY;
    if (model.objects[1].rightProp != NULL) {
      disallowedProp = model.objects[1].rightProp->propNum;
    }
  }
  for (int y = max(coverY - height / 20, 0); y < coverY + height / 15 && y < height; ++y) {
    for (int x = coverX; x < coverX + width / 15 && x < width; ++x) {
      int cur = y * width + x + 2;
      int r = (arr[cur] >> 16) & 0x00ff;
      int g = (arr[cur] >> 8) & 0x00ff;
      int b = (arr[cur]) & 0x00ff;
     // if (g > 150 && max(r, b) < min(180, g)) {
      if (g > 150 && g > max(r, b) && g - max(r, b) > 30) {
          //&& scc[x][y] != disallowedProp && sccProp.find(scc[x][y]) != sccProp.end()) {
      //if (sumImage.isForeground[cur - 2] && g > max(r, b)) {//sumImage.isGreen[cur - 2]) {
        xpos += x;
        ypos += y;
        ++count;
      }
    }
  }
  if (count > (width < 1700 ? 10 : 20)) {
    cerr << "Found panel LED with " << count << " green pixels" << endl;
    xpos /= count;
    ypos /= count;
    if (isLeft) {
      model.objects[2].leftX = (int) xpos;
      model.objects[2].leftY = (int) ypos;
      model.objects[2].leftFound = true;
    } else {
      model.objects[2].rightX = (int) xpos;
      model.objects[2].rightY = (int) ypos;
      model.objects[2].rightFound = true;
    }
  } else {
    cerr << "Failed to find panel LED" << endl;
  }
}

static void findRockerSwitch(const vector<int>& arr,
    vector<vector<int> >& scc, map<int, SccProperty>& sccProp, map<int, SccProperty>& sccPropLarge,
    SumImage& sumImage,
    ModelSummary& model, bool isLeft) {
  int width = arr[1];
  int height = arr[0];
  int coverX, coverY;
  if (isLeft) {
    coverX = model.objects[1].leftX;
    coverY = model.objects[1].leftY;
  } else {
    coverX = model.objects[1].rightX;
    coverY = model.objects[1].rightY;
  }
  int ymin = min(coverY + height / 15, height);
  int ymax = min(coverY + height / 3, height);
  int xmin = max(coverX - width / 20, 0);
  int xmax = min(coverX + width / 5, width);
  int maxSize = 0;
  SccProperty* maxProp = NULL;
  for (int y = ymin; y < ymax; ++y) {
    if (sumImage.isBackground[xmin + width * y] &&
        sccPropLarge.find(scc[xmin][y]) != sccPropLarge.end()) {
      // Stop search if we've gone down far enough to hit a large component.
      break;
    }
    for (int x = xmin; x < xmax; ++x) {
      if (sumImage.isBackground[x + width * y] &&
          sccPropLarge.find(scc[x][y]) != sccPropLarge.end()) {
        // Go to next row if we've gone to the right far enough to hit a large component.
        break;
      }
      map<int, SccProperty>::iterator it = sccProp.find(scc[x][y]);
      if (it != sccProp.end()) {
        SccProperty& prop = it->second;
        if (prop.size > maxSize &&
            mightBeRocker(width, height, scc, coverY + height / 10, prop)) {
          maxSize = prop.size;
          maxProp = &prop;
        }
      }
    }
  }
  if (maxProp != NULL) {
    cerr << "Found Rocker!" << endl;
    if (isLeft) {
      model.objects[3].leftX = maxProp->x;
      model.objects[3].leftY = maxProp->y;
      model.objects[3].leftFound = true;
    } else {
      model.objects[3].rightX = maxProp->x;
      model.objects[3].rightY = maxProp->y;
      model.objects[3].rightFound = true;
    }
  } else {
    cerr << "Failed to find rocker." << endl;
  }
}

struct DiffVect {
  int x, y;
  double length;
  DiffVect(int x1, int y1, double length1) : x(x1), y(y1), length(length1) {
  }
};

// matchup to *index* of possibleComps.
// matchup.length == greenLeds.size().
// used.length == possibleComps.length.
static bool dfsMatchStrict = false;
static SccProperty* panelLedForDfs = NULL;
static SccProperty* rockerForDfs = NULL;
static void dfsMatchRect(vector<int>& matchup, int idx, vector<int>& possibleComps, vector<bool>& used,
    vector<SccProperty*>& leds,
    int ledX, int ledY,
    double a1, double b1, double a2, double b2,
    vector<int>& bestMatchup, double* bestDist,
    double compRatiosOrtho[22][2]) {

  int rockerX = (int) (ledX + a1);
  int rockerY = (int) (ledY + b1);
  if (idx == matchup.size()) {
    // All matched up. Compute sum of dist squareds.
    vector<DiffVect> diffVectors;
    double distSum = 0;
    double minScale = 999;
    double maxScale = -999;
    double minTheta = 999;
    double maxTheta = -999;
    for (int i = 0; i < matchup.size(); ++i) {
      SccProperty& prop = *leds[i];

      int comp = possibleComps[matchup[i]];
      int expX = (int) (compRatiosOrtho[comp][0] * a1 + compRatiosOrtho[comp][1] * a2);
      int expY = (int) (compRatiosOrtho[comp][0] * b1 + compRatiosOrtho[comp][1] * b2);
      expX += ledX;
      expY += ledY;
      int dx2 = expX - prop.x;
      int dy2 = expY - prop.y;
      double dist = dx2 * dx2 + dy2 * dy2;
      distSum += dist;
      double length = sqrt(dist);
      if (length > 10) {
        diffVectors.push_back(DiffVect(dx2, dy2, length));
      }
      expX -= rockerX;
      expY -= rockerY;
      double expLength = sqrt(expX * expX + expY * expY);
      int actX = prop.x - rockerX;
      int actY = prop.y - rockerY;
      double actLength = sqrt(actX * actX + actY * actY);
      double scaleConv = actLength / expLength;
      double crossProduct = expX * actY - expY * actX;
      crossProduct /= actLength;
      crossProduct /= expLength;
      double thetaConv = asin(crossProduct);
      if (scaleConv < minScale) minScale = scaleConv;
      if (scaleConv > maxScale) maxScale = scaleConv;
      if (thetaConv < minTheta) minTheta = thetaConv;
      if (thetaConv > maxTheta) maxTheta = thetaConv;
    }
    distSum /= matchup.size();
    if (diffVectors.size() > 1) {
      distSum *= (maxScale - minScale) + 1;
      distSum *= (maxTheta - minTheta) + 1;
      double pairwiseAngleSum = 0;
      int numSum = 0;
      for (int i = 0; i < diffVectors.size(); ++i) {
        DiffVect& v1 = diffVectors[i];
        for (int j = i + 1; j < diffVectors.size(); ++j) {
          DiffVect& v2 = diffVectors[j];
          double dotProduct = v1.x * v2.x + v1.y * v2.y;
          dotProduct /= v1.length;
          dotProduct /= v2.length;
          pairwiseAngleSum += abs(1.0 - dotProduct) + .5;
          ++numSum;
        }
      }
      if (numSum > 0) {
        pairwiseAngleSum /= numSum;
      } else {
        pairwiseAngleSum = 1.0;
      }
      distSum *= pairwiseAngleSum;
    }
    if (*bestDist < 0 || distSum < *bestDist) {
      //cerr << "Got better matchup of dist: " << distSum << endl;
      *bestDist = distSum;
      for (int i = 0; i < matchup.size(); ++i) {
       // cerr << leds[i]->x << "," << leds[i]->y << " -> " <<
       //     possibleComps[matchup[i]] << endl;
        bestMatchup[i] = matchup[i];
      }
    }
    return;
  }
  for (int i = 0; i < used.size(); ++i) {
    if (used[i]) continue;
    bool slotIsGreen = true;
    switch (possibleComps[i]) {
      case 6:
      case 8:
      case 10:
      case 12:
      case 14:
      case 19:
        slotIsGreen = false;
        break;
      default:
        break;
    }
    /*if ((possibleComps[i] < 6 || possibleComps[i] > 14) &&
        !leds[idx]->containsCenter) continue;*/
    if (leds[idx]->isGreenLed && !slotIsGreen) continue;
    if (leds[idx]->isBlueLed && slotIsGreen) continue;
    int topoY[22] = {0, 0, 0, 2, 1, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 5, 4, 5, 4, 6, 5, 4};
    int topoX[22] = {3, 3, 3, 3, 3, 3, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 0, 2, 2, 2, 3, 3};
    //if (leds[idx]->isGreenLed || leds[idx]->isBlueLed) {
    {
      vector<SccProperty*> curMatchup(22);
      for (int j = 0; j < idx; ++j) {
        curMatchup[possibleComps[matchup[j]]] = leds[j];
      }
      int ymin = 0;
      int ymax = 9999;
      int xmin = 0;
      int xmax = 9999;
      if (panelLedForDfs != NULL) {
        ymin = panelLedForDfs->y;
      }
      if (rockerForDfs != NULL) {
        xmax = rockerForDfs->xmax + 200; 
      }
      int topoClassY = topoY[possibleComps[i]];
      int topoClassX = topoX[possibleComps[i]];
      for (int j = 0; j < 22; ++j) {
        if (curMatchup[j] != NULL) {
          if (topoY[j] < topoClassY) {
            if (dfsMatchStrict) {
              ymin = max(ymin, curMatchup[j]->ymax);
            } else {
              ymin = max(ymin, curMatchup[j]->y);
            }
          } else if (topoY[j] > topoClassY) {
            if (dfsMatchStrict) {
              ymax = min(ymax, curMatchup[j]->ymin);
            } else {
              ymax = min(ymax, curMatchup[j]->y);
            }
          }
          if (topoX[j] < topoClassX) {
            if (dfsMatchStrict) {
              xmin = max(xmin, curMatchup[j]->xmax);
            } else {
              xmin = max(xmin, curMatchup[j]->x);
            }
          } else if (topoX[j] > topoClassX) {
            if (dfsMatchStrict) {
              xmax = min(xmax, curMatchup[j]->xmin);
            } else {
              xmax = min(xmax, curMatchup[j]->x);
            }
          }
        }
      }
      if (leds[idx]->y < ymin) {
        continue;
      }
      if (leds[idx]->y >= ymax) {
        continue;
      }
      if (leds[idx]->x < xmin) {
        continue;
      }
      if (leds[idx]->x >= xmax) {
        continue;
      }
    }
    used[i] = true;
    matchup[idx] = i;
    dfsMatchRect(matchup, idx + 1, possibleComps, used, leds, ledX, ledY, a1, b1, a2, b2,
        bestMatchup, bestDist, compRatiosOrtho);
    used[i] = false;
  }
}

static bool matchSccToBothLedsDfsRect(const vector<int>& arr, vector<vector<int> >& scc,
    map<int, SccProperty>& sccProp, map<int, SccProperty>& sccPropLarge,
    ModelSummary& model, bool isLeft) {
  int width = arr[1];
  int height = arr[0];
  int ledX, ledY;
  int rockerX, rockerY;
  if (isLeft) {
    ledX = model.objects[2].leftX;
    ledY = model.objects[2].leftY;
    rockerX = model.objects[3].leftX;
    rockerY = model.objects[3].leftY;
  } else {
    ledX = model.objects[2].rightX;
    ledY = model.objects[2].rightY;
    rockerX = model.objects[3].rightX;
    rockerY = model.objects[3].rightY;
  }
  int dx = rockerX - ledX;
  int dy = rockerY - ledY;
  double a1, b1;
  double a2, b2;
  {
    a1 = dx;
    b1 = dy;
    a2 = b1;
    b2 = -a1;
    int length = dx * dx + dy * dy;
    if (length < 100) return false;
    double sinVal = dx / sqrt(length);
    //double theta = asin(sinVal) * 2;
    double theta = asin(sinVal);
    double a3 = cos(theta) * a2 - sin(theta) * b2;
    double b3 = sin(theta) * a2 + cos(theta) * b2;
    a2 = a3;
    b2 = b3;
  }
  vector<SccProperty*> blueLeds;
  {
    set<pair<int, int> > blueLedsSet;
    int blueLedsCount = 0;
    for (map<int, SccProperty>::iterator it = sccProp.begin(); it != sccProp.end(); ++it) {
      SccProperty& prop = it->second;
      if (!prop.isBlueLed) continue;
      if (blueLedsSet.find(make_pair(prop.x, prop.y)) != blueLedsSet.end()) continue;
      if (!canReach(width, height, prop.x, prop.y, ledX, ledY, scc, sccPropLarge)) continue;
      int compX = prop.x;
      int compY = prop.y;
      int dx2 = compX - ledX;
      int dy2 = compY - ledY;
      int length2 = dx2 * dx2 + dy2 * dy2;
      if (length2 < 100) {
        cerr << "Skipping a component of length2: " << length2 << endl;
        continue;
      }
      ++blueLedsCount;
      // Assuming panelLED is correct, throw out obviously wrong ones.
      if (prop.x > ledX) continue;

      blueLedsSet.insert(make_pair(prop.x, prop.y));
      blueLeds.push_back(&prop);
    }
    if (blueLedsCount > 12) {
      cerr << "Way too many blueLeds: " << blueLeds.size() << " bailing out." << endl;
      blueLeds.clear();
    } else if (blueLeds.size() > 6) {
      cerr << "Error, got more blueLeds than possible: " << blueLeds.size()
           << " trimming..." << endl;
      blueLeds.resize(6);
    }
  }
  vector<SccProperty*> greenLeds;
  {
    set<pair<int, int> > greenLedsSet;
    int greenLedsCount = 0;
    for (map<int, SccProperty>::iterator it = sccProp.begin(); it != sccProp.end(); ++it) {
      SccProperty& prop = it->second;
      if (!prop.isGreenLed) continue;
      if (greenLedsSet.find(make_pair(prop.x, prop.y)) != greenLedsSet.end()) continue;
      if (!canReach(width, height, prop.x, prop.y, ledX, ledY, scc, sccPropLarge)) continue;
      int compX = prop.x;
      int compY = prop.y;
      int dx2 = compX - ledX;
      int dy2 = compY - ledY;
      int length2 = dx2 * dx2 + dy2 * dy2;
      if (length2 < 100) {
        cerr << "Skipping a component of length2: " << length2 << endl;
        continue;
      }
      ++greenLedsCount;

      greenLedsSet.insert(make_pair(prop.x, prop.y));
      greenLeds.push_back(&prop);
    }
    if (greenLeds.size() > 8) {
      cerr << "Error, got more greenLeds than possible: " << greenLeds.size()
           << " trimming..." << endl;
      greenLeds.resize(8);
    }
  }
  vector<SccProperty*> allLeds;
  allLeds.insert(allLeds.end(), blueLeds.begin(), blueLeds.end());
  allLeds.insert(allLeds.end(), greenLeds.begin(), greenLeds.end());
  int possibleCompsArr[14] = { 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 19, 21 };
  vector<int> possibleComps(14);
  for (int i = 0; i < 14; ++i) {
    possibleComps[i] = possibleCompsArr[i];
  }
  vector<int> matchup(allLeds.size());
  vector<bool> used(14);
  for (int i = 0; i < 14; ++i) {
    used[i] = false;
  }
  vector<int> bestMatchup(allLeds.size());
  double bestDist = -1;
  dfsMatchRect(matchup, 0, possibleComps, used, allLeds, ledX, ledY, a1, b1, a2, b2,
      bestMatchup, &bestDist, isLeft ? compRatiosOrthoLeft : compRatiosOrthoRight);
  if (bestDist < 0) return false;
  for (int i = 0; i < allLeds.size(); ++i) {
    SccProperty& prop = *allLeds[i];
    int matchedComp = possibleComps[bestMatchup[i]];
    cerr << "Dfs Matched prop " << prop.x << "," << prop.y << " to comp " << matchedComp << endl;
    if (isLeft) {
      model.objects[matchedComp].leftX = prop.x;
      model.objects[matchedComp].leftY = prop.y;
      model.objects[matchedComp].leftFound = true;
      model.objects[matchedComp].state = "ON";
    } else {
      model.objects[matchedComp].rightX = prop.x;
      model.objects[matchedComp].rightY = prop.y;
      model.objects[matchedComp].rightFound = true;
      model.objects[matchedComp].state = "ON";
    }
  }
  return true;
}

static void matchSccToGreenLedsDfsRect(const vector<int>& arr, vector<vector<int> >& scc,
    map<int, SccProperty>& sccProp, map<int, SccProperty>& sccPropLarge,
    ModelSummary& model, bool isLeft) {
  int width = arr[1];
  int height = arr[0];
  int ledX, ledY;
  int rockerX, rockerY;
  if (isLeft) {
    ledX = model.objects[2].leftX;
    ledY = model.objects[2].leftY;
    rockerX = model.objects[3].leftX;
    rockerY = model.objects[3].leftY;
  } else {
    ledX = model.objects[2].rightX;
    ledY = model.objects[2].rightY;
    rockerX = model.objects[3].rightX;
    rockerY = model.objects[3].rightY;
  }

  vector<SccProperty*> greenLeds;
  set<pair<int, int> > greenLedsSet;
  for (map<int, SccProperty>::iterator it = sccProp.begin(); it != sccProp.end(); ++it) {
    SccProperty& prop = it->second;
    if (!prop.isGreenLed) continue;
    if (greenLedsSet.find(make_pair(prop.x, prop.y)) != greenLedsSet.end()) continue;
    if (!canReach(width, height, prop.x, prop.y, ledX, ledY, scc, sccPropLarge)) continue;
    int compX = prop.x;
    int compY = prop.y;
    int dx2 = compX - ledX;
    int dy2 = compY - ledY;
    int length2 = dx2 * dx2 + dy2 * dy2;
    if (length2 < 100) {
      cerr << "Skipping a component of length2: " << length2 << endl;
      continue;
    }
    int propCircleArea = (prop.xmax - prop.xmin) * (prop.ymax - prop.ymin) * 3 / 4;
    if (propCircleArea == 0) propCircleArea = 1;
    double fillRatio = prop.size / (double) (propCircleArea);
    bool centerEmpty = false;
    if (sccProp.find(scc[prop.x][prop.y]) == sccProp.end()) {
      centerEmpty = true;
      prop.containsCenter = false;
    } else {
      SccProperty& centerProp = sccProp[scc[prop.x][prop.y]];
      if (centerProp.x != prop.x || centerProp.y != prop.y) {
        centerEmpty = true;
      }
    }

//    if (fillRatio < 0.9 && centerEmpty) {
//      greenDonuts.push_back(&prop);
//    } else {
//      greenDots.push_back(&prop);
//    }
    greenLeds.push_back(&prop);
    greenLedsSet.insert(make_pair(prop.x, prop.y));
  }
  {
    int possibleCompsArr[8] = { 4, 7, 9, 11, 13, 16, 18, 21 };
    vector<int> possibleCompsDonut(8);
    for (int i = 0; i < 8; ++i) {
      possibleCompsDonut[i] = possibleCompsArr[i];
    }
    if (greenLeds.size() > possibleCompsDonut.size()) {
      cerr << "Error, got more greenLeds than possible: " << greenLeds.size()
           << " trimming..." << endl;
      greenLeds.resize(possibleCompsDonut.size());
    }
    int dx = rockerX - ledX;
    int dy = rockerY - ledY;
    double a1, b1;
    double a2, b2;
    {
      a1 = dx;
      b1 = dy;
      a2 = b1;
      b2 = -a1;
      int length = dx * dx + dy * dy;
      if (length < 100) return;
      double sinVal = dx / sqrt(length);
      //double theta = asin(sinVal) * 2;
      double theta = asin(sinVal);
      double a3 = cos(theta) * a2 - sin(theta) * b2;
      double b3 = sin(theta) * a2 + cos(theta) * b2;
      a2 = a3;
      b2 = b3;
    }
    vector<int> matchup(greenLeds.size());
    vector<bool> used(possibleCompsDonut.size());
    for (int i = 0; i < possibleCompsDonut.size(); ++i) {
      used[i] = false;
    }
    vector<int> bestMatchup(greenLeds.size());
    double bestDist = -1;
    dfsMatchRect(matchup, 0, possibleCompsDonut, used, greenLeds, ledX, ledY, a1, b1, a2, b2,
        bestMatchup, &bestDist, isLeft ? compRatiosOrthoLeft : compRatiosOrthoRight);
    if (bestDist < 0) return;
    for (int i = 0; i < greenLeds.size(); ++i) {
      SccProperty& prop = *greenLeds[i];
      int matchedComp = possibleCompsDonut[bestMatchup[i]];
      cerr << "Dfs Matched prop " << prop.x << "," << prop.y << " to comp " << matchedComp << endl;
      if (isLeft) {
        model.objects[matchedComp].leftX = prop.x;
        model.objects[matchedComp].leftY = prop.y;
        model.objects[matchedComp].leftFound = true;
        model.objects[matchedComp].state = "ON";
      } else {
        model.objects[matchedComp].rightX = prop.x;
        model.objects[matchedComp].rightY = prop.y;
        model.objects[matchedComp].rightFound = true;
        model.objects[matchedComp].state = "ON";
      }
    }
  }
}

static void matchSccToBlueLedsDfsRect(const vector<int>& arr, vector<vector<int> >& scc,
    map<int, SccProperty>& sccProp, map<int, SccProperty>& sccPropLarge,
    ModelSummary& model, bool isLeft) {
  int possibleCompsArr[6] = { 6, 8, 10, 12, 14, 19 };
  vector<int> possibleComps(6);
  for (int i = 0; i < 6; ++i) {
    possibleComps[i] = possibleCompsArr[i];
  }
  int width = arr[1];
  int height = arr[0];
  int ledX, ledY;
  int rockerX, rockerY;
  if (isLeft) {
    ledX = model.objects[2].leftX;
    ledY = model.objects[2].leftY;
    rockerX = model.objects[3].leftX;
    rockerY = model.objects[3].leftY;
  } else {
    ledX = model.objects[2].rightX;
    ledY = model.objects[2].rightY;
    rockerX = model.objects[3].rightX;
    rockerY = model.objects[3].rightY;
  }
  int dx = rockerX - ledX;
  int dy = rockerY - ledY;
  double a1, b1;
  double a2, b2;
  {
    a1 = dx;
    b1 = dy;
    a2 = b1;
    b2 = -a1;
    int length = dx * dx + dy * dy;
    if (length < 100) return;
    double sinVal = dx / sqrt(length);
    //double theta = asin(sinVal) * 2;
    double theta = asin(sinVal);
    double a3 = cos(theta) * a2 - sin(theta) * b2;
    double b3 = sin(theta) * a2 + cos(theta) * b2;
    a2 = a3;
    b2 = b3;
  }

  vector<SccProperty*> blueLeds;
  set<pair<int, int> > blueLedsSet;
  int blueLedsCount = 0;
  for (map<int, SccProperty>::iterator it = sccProp.begin(); it != sccProp.end(); ++it) {
    SccProperty& prop = it->second;
    if (!prop.isBlueLed) continue;
    if (blueLedsSet.find(make_pair(prop.x, prop.y)) != blueLedsSet.end()) continue;
    if (!canReach(width, height, prop.x, prop.y, ledX, ledY, scc, sccPropLarge)) continue;
    int compX = prop.x;
    int compY = prop.y;
    int dx2 = compX - ledX;
    int dy2 = compY - ledY;
    int length2 = dx2 * dx2 + dy2 * dy2;
    if (length2 < 100) {
      cerr << "Skipping a component of length2: " << length2 << endl;
      continue;
    }
    ++blueLedsCount;
    // Assuming panelLED is correct, throw out obviously wrong ones.
    if (prop.x > ledX) continue;

    blueLedsSet.insert(make_pair(prop.x, prop.y));
    blueLeds.push_back(&prop);
  }
  if (blueLedsCount > 12) {
    cerr << "Way too many blueLeds: " << blueLeds.size() << " bailing out." << endl;
    return;
  } else if (blueLeds.size() > 6) {
    cerr << "Error, got more blueLeds than possible: " << blueLeds.size()
         << " trimming..." << endl;
    blueLeds.resize(6);
  }
  vector<int> matchup(blueLeds.size());
  vector<bool> used(6);
  for (int i = 0; i < 6; ++i) {
    used[i] = false;
  }
  vector<int> bestMatchup(blueLeds.size());
  double bestDist = -1;
  dfsMatchRect(matchup, 0, possibleComps, used, blueLeds, ledX, ledY, a1, b1, a2, b2,
      bestMatchup, &bestDist, isLeft ? compRatiosOrthoLeft : compRatiosOrthoRight);
  if (bestDist < 0) return;
  for (int i = 0; i < blueLeds.size(); ++i) {
    SccProperty& prop = *blueLeds[i];
    int matchedComp = possibleComps[bestMatchup[i]];
    cerr << "Dfs Matched prop " << prop.x << "," << prop.y << " to comp " << matchedComp << endl;
    if (isLeft) {
      model.objects[matchedComp].leftX = prop.x;
      model.objects[matchedComp].leftY = prop.y;
      model.objects[matchedComp].leftFound = true;
      model.objects[matchedComp].state = "ON";
    } else {
      model.objects[matchedComp].rightX = prop.x;
      model.objects[matchedComp].rightY = prop.y;
      model.objects[matchedComp].rightFound = true;
      model.objects[matchedComp].state = "ON";
    }
  }
}

static bool gridHasObj(ModelSummary& model, bool isLeft, int obj) {
  if (isLeft) {
    return model.objects[obj].leftFound;
  } else {
    return model.objects[obj].rightFound;
  }
}

static bool compareSccBySize(SccProperty* prop0, SccProperty* prop1) {
  int size0 = (prop0->xmax - prop0->xmin) * (prop0->ymax - prop0->ymin);
  int size1 = (prop1->xmax - prop1->xmin) * (prop1->ymax - prop1->ymin);
  return size0 < size1;
  //return prop0->size < prop1->size;
}
static bool compareSccByY(SccProperty* prop0, SccProperty* prop1) {
  return prop0->y < prop1->y;
}

static void doBasisDfsMatch(const vector<int>& arr, vector<vector<int> >& scc,
    map<int, SccProperty>& sccProp, map<int, SccProperty>& sccPropLarge,
    ModelSummary& model, bool isLeft,
    double a1, double b1, double a2, double b2, double compRatios[22][2]) {
  int width = arr[1];
  int height = arr[0];
  int ledX, ledY;
  if (isLeft) {
    ledX = model.objects[2].leftX;
    ledY = model.objects[2].leftY;
  } else {
    ledX = model.objects[2].rightX;
    ledY = model.objects[2].rightY;
  }
  vector<SccProperty*> allLeds;
  for (map<int, SccProperty>::iterator it = sccProp.begin(); it != sccProp.end(); ++it) {
    SccProperty& prop = it->second;
    if (prop.isGreenLed || prop.isBlueLed) {
      allLeds.push_back(&prop);
    }
  }

  // For each not-found grid component.
  double maxDistFromPredicted = height / 27;
  maxDistFromPredicted *= maxDistFromPredicted;
  int64 maxDistForMerge = (width + height) / 2 * 25L / 1000L;
  maxDistForMerge *= maxDistForMerge;
  set<int> disallowedScc;  // Used by other already-found ones.
  for (int i = 6; i <= 14; ++i) {
    int takenX = 0, takenY = 0;
    if (isLeft && model.objects[i].leftFound) {
      takenX = model.objects[i].leftX;
      takenY = model.objects[i].leftY;
    } else if (!isLeft && model.objects[i].rightFound) {
      takenX = model.objects[i].rightX;
      takenY = model.objects[i].rightY;
    } else {
      continue;
    }
    for (map<int, SccProperty>::iterator it = sccProp.begin(); it != sccProp.end(); ++it) {
      SccProperty& prop = it->second;
      if (prop.isGreenLed || prop.isBlueLed) continue;
      int dx = prop.x - takenX;
      int dy = prop.y - takenY;
      double dist = dx * dx + dy * dy;
      if (dist < maxDistFromPredicted) {
        disallowedScc.insert(prop.propNum);
      }
    }
  }
  vector<SccProperty*> allGridScc;
  vector<int> componentsMissing;
  for (int i = 6; i <= 14; ++i) {
    if (isLeft && model.objects[i].leftFound) continue;
    if (!isLeft && model.objects[i].rightFound) continue;
    componentsMissing.push_back(i);

    // Find its predicted position.
    int newX = (int) (compRatios[i][0] * a1 + compRatios[i][1] * a2);
    int newY = (int) (compRatios[i][0] * b1 + compRatios[i][1] * b2);
    newX += ledX;
    newY += ledY;

    // Then find all SCC in its vicinity.
    for (map<int, SccProperty>::iterator it = sccProp.begin(); it != sccProp.end(); ++it) {
      SccProperty& prop = it->second;
      if (prop.isGreenLed || prop.isBlueLed) continue;
      if (disallowedScc.find(prop.propNum) != disallowedScc.end()) continue;
      //if (!canReach(width, height, prop.x, prop.y, ledX, ledY, scc, sccPropLarge)) continue;
      int dx = prop.x - newX;
      int dy = prop.y - newY;
      double dist = dx * dx + dy * dy;
      if (dist < maxDistFromPredicted) {
        // As long as it's not just part of another LED.
        bool isEdgeOfLed = false;
        for (int i = 0; i < allLeds.size(); ++i) {
          SccProperty& led = *allLeds[i];
          dx = prop.x - led.x;
          dy = prop.y - led.y;
          if (dx * dx + dy * dy < maxDistForMerge) {
            cerr << "Skipping SCC " << prop.x << "," << prop.y << " due to proximity " <<
                "to led: " << led.x << "," << led.y << endl;
            isEdgeOfLed = true;
            break;
          }
        }
        if (!isEdgeOfLed) {
          allGridScc.push_back(&prop);
        }
      }
    }
  }
  if (componentsMissing.size() == 0) return;
  if (allGridScc.size() == 0) return;

  // Merge gridScc which are close to each other.
  vector<SccProperty*> mergedGrid;
  bool hadMerge = true;
  while (hadMerge) {
    set<int> alreadyMergedGrid;
    hadMerge = false;
    for (int i = 0; i < allGridScc.size(); ++i) {
      SccProperty& prop0 = *allGridScc[i];
      if (alreadyMergedGrid.find(prop0.propNum) != alreadyMergedGrid.end()) continue;
      vector<int> mergedComponents;
      alreadyMergedGrid.insert(prop0.propNum);
      mergedComponents.push_back(prop0.propNum);
      for (int j = i + 1; j < allGridScc.size(); ++j) {
        SccProperty& prop1 = *allGridScc[j];
        if (alreadyMergedGrid.find(prop1.propNum) != alreadyMergedGrid.end()) continue;
        int dx = prop0.x - prop1.x;
        int dy = prop0.y - prop1.y;
        if (dx * dx + dy * dy < maxDistForMerge) {
          mergeScc(prop0, prop1);
          alreadyMergedGrid.insert(prop1.propNum);
          mergedComponents.push_back(prop1.propNum);
          hadMerge = true;
        }
      }
      for (int k = 0; k < mergedComponents.size(); ++k) {
        sccProp[mergedComponents[k]] = prop0;
      }
      mergedGrid.push_back(&prop0);
    }
    if (hadMerge) {
      allGridScc = mergedGrid;
      mergedGrid.clear();
    }
  }

  if (mergedGrid.size() > componentsMissing.size()) {
    cerr << "Found " << mergedGrid.size() << " SCC for " << componentsMissing.size() <<
        " possible components; bailing out." << endl;
    //sort(mergedGrid.begin(), mergedGrid.end(), compareSccBySize);
    mergedGrid.resize(componentsMissing.size());
    return;
  }

  // Even though we may have used a different basis to get "nearby" components, for now
  // dfsMatchRect requires the rotated-ortho version.
  int rockerX, rockerY;
  if (isLeft) {
    rockerX = model.objects[3].leftX;
    rockerY = model.objects[3].leftY;

  } else {
    rockerX = model.objects[3].rightX;
    rockerY = model.objects[3].rightY;
  }
  {
    int dx = rockerX - ledX;
    int dy = rockerY - ledY;
    a1 = dx;
    b1 = dy;
    a2 = b1;
    b2 = -a1;

    double length = dx * dx + dy * dy;
    double sinVal = dx / sqrt(length);
    //double theta = asin(sinVal) * 2;
    double theta = asin(sinVal);
    double a3 = cos(theta) * a2 - sin(theta) * b2;
    double b3 = sin(theta) * a2 + cos(theta) * b2;
    a2 = a3;
    b2 = b3;
  }

  vector<int> possibleComps = componentsMissing;
  vector<int> matchup(mergedGrid.size());
  vector<bool> used(possibleComps.size());
  vector<int> bestMatchup(mergedGrid.size());
  double bestDist = -1;
  dfsMatchRect(matchup, 0, possibleComps, used, mergedGrid, ledX, ledY, a1, b1, a2, b2,
      bestMatchup, &bestDist, isLeft ? compRatiosOrthoLeft : compRatiosOrthoRight);
  if (bestDist < 0) return;
  for (int i = 0; i < mergedGrid.size(); ++i) {
    SccProperty& prop = *mergedGrid[i];
    int matchedComp = possibleComps[bestMatchup[i]];
    cerr << "Dfs Matched prop " << prop.x << "," << prop.y << " to comp " << matchedComp << endl;
    if (isLeft) {
      model.objects[matchedComp].leftX = prop.x;
      model.objects[matchedComp].leftY = prop.y;
      model.objects[matchedComp].leftFound = true;
    } else {
      model.objects[matchedComp].rightX = prop.x;
      model.objects[matchedComp].rightY = prop.y;
      model.objects[matchedComp].rightFound = true;
    }
  }
}

static void basisDfsMatch(const vector<int>& arr, vector<vector<int> >& scc,
    map<int, SccProperty>& sccProp, map<int, SccProperty>& sccPropLarge,
    ModelSummary& model, bool isLeft, int gridComp, double compRatios[22][2]) {
  int width = arr[1];
  int height = arr[0];
  int ledX, ledY;
  int rockerX, rockerY;
  int gridX, gridY;
  if (isLeft) {
    ledX = model.objects[2].leftX;
    ledY = model.objects[2].leftY;
    rockerX = model.objects[3].leftX;
    rockerY = model.objects[3].leftY;
    gridX = model.objects[gridComp].leftX;
    gridY = model.objects[gridComp].leftY;

  } else {
    ledX = model.objects[2].rightX;
    ledY = model.objects[2].rightY;
    rockerX = model.objects[3].rightX;
    rockerY = model.objects[3].rightY;
    gridX = model.objects[gridComp].rightX;
    gridY = model.objects[gridComp].rightY;
  }
  double a1, b1;
  double a2, b2;
  {
    int dx = rockerX - ledX;
    int dy = rockerY - ledY;
    a1 = dx;
    b1 = dy;
  }
  {
    int dx = gridX - rockerX;
    int dy = gridY - rockerY;
    a2 = dx;
    b2 = dy;
  }
  doBasisDfsMatch(arr, scc, sccProp, sccPropLarge, model, isLeft, a1, b1, a2, b2, compRatios);
}

static void basisDfsMatchWeak(const vector<int>& arr, vector<vector<int> >& scc,
    map<int, SccProperty>& sccProp, map<int, SccProperty>& sccPropLarge,
    ModelSummary& model, bool isLeft) {
  int width = arr[1];
  int height = arr[0];
  int ledX, ledY;
  int rockerX, rockerY;
  if (isLeft) {
    ledX = model.objects[2].leftX;
    ledY = model.objects[2].leftY;
    rockerX = model.objects[3].leftX;
    rockerY = model.objects[3].leftY;

  } else {
    ledX = model.objects[2].rightX;
    ledY = model.objects[2].rightY;
    rockerX = model.objects[3].rightX;
    rockerY = model.objects[3].rightY;
  }
  cerr << "Doing basisDfsMatchWeak with " << ledX << "," << ledY << " "
       << rockerY << "," << rockerY << endl;
  double a1, b1;
  double a2, b2;
  {
    int dx = rockerX - ledX;
    int dy = rockerY - ledY;
    a1 = dx;
    b1 = dy;
    a2 = b1;
    b2 = -a1;

    double length = dx * dx + dy * dy;
    double sinVal = dx / sqrt(length);
    //double theta = asin(sinVal) * 2;
    double theta = asin(sinVal);
    double a3 = cos(theta) * a2 - sin(theta) * b2;
    double b3 = sin(theta) * a2 + cos(theta) * b2;
    a2 = a3;
    b2 = b3;
  }

  // TODO: maybe include 20 and 21?
  double avgScale = 0;
  double avgTheta = 0;
  int numSamples = 0;
  // For each component, find the scale and rotation correction relative to the vector
  // made by that component to the rocker as compared to the predicted position, and
  // average the values in.
  vector<double> scaleSamples;
  vector<double> thetaSamples;
  for (int comp = 6; comp <= 18; ++comp) {
    int actX = -1, actY = -1;
    if (isLeft) {
      if (model.objects[comp].leftFound) {
        actX = model.objects[comp].leftX;
        actY = model.objects[comp].leftY;
      }
    } else {
      if (model.objects[comp].rightFound) {
        actX = model.objects[comp].rightX;
        actY = model.objects[comp].rightY;
      }
    }
    if (actX == -1 || actY == -1) continue;
    actX -= rockerX;
    actY -= rockerY;
    int predX, predY;
    if (isLeft) {
      predX = (int) (compRatiosOrthoLeft[comp][0] * a1 + compRatiosOrthoLeft[comp][1] * a2) + ledX;
      predY = (int) (compRatiosOrthoLeft[comp][0] * b1 + compRatiosOrthoLeft[comp][1] * b2) + ledY;
    } else {
      predX = (int) (compRatiosOrthoRight[comp][0] * a1 + compRatiosOrthoRight[comp][1] * a2) + ledX;
      predY = (int) (compRatiosOrthoRight[comp][0] * b1 + compRatiosOrthoRight[comp][1] * b2) + ledY;
    }
    predX -= rockerX;
    predY -= rockerY;
    double predLength = sqrt(predX * predX + predY * predY);
    double actLength = sqrt(actX * actX + actY * actY);
    double scaleConv = actLength / predLength;
    double crossProduct = predX * actY - predY * actX;
    crossProduct /= actLength;
    crossProduct /= predLength;
    double thetaConv = asin(crossProduct);
    avgScale += scaleConv - 1;
    avgTheta += thetaConv;
    ++numSamples;
    scaleSamples.push_back(scaleConv - 1);
    thetaSamples.push_back(thetaConv);
  }
  /*if (numSamples > 0) {
    avgScale = (scaleSamples[(numSamples - 1) / 2] + scaleSamples[numSamples / 2]) / 2;
    avgTheta = (thetaSamples[(numSamples - 1) / 2] + thetaSamples[numSamples / 2]) / 2;
    avgScale += 1;
    // Correct the normal vector by the average correction factors.
    double a3 = cos(avgTheta) * a2 - sin(avgTheta) * b2;
    double b3 = sin(avgTheta) * a2 + cos(avgTheta) * b2;
    a3 *= avgScale;
    b3 *= avgScale;
    a2 = a3;
    b2 = b3;
  }*/
  if (numSamples > 0) {
    avgScale /= numSamples;
    avgTheta /= numSamples;
    avgScale += 1;
    // Correct the normal vector by the average correction factors.
    double a3 = cos(avgTheta) * a2 - sin(avgTheta) * b2;
    double b3 = sin(avgTheta) * a2 + cos(avgTheta) * b2;
    a3 *= avgScale;
    b3 *= avgScale;
    a2 = a3;
    b2 = b3;
  }
  doBasisDfsMatch(arr, scc, sccProp, sccPropLarge, model, isLeft, a1, b1, a2, b2,
      isLeft ? compRatiosOrthoLeft : compRatiosOrthoRight);
}

static bool blockBasisExtrapolate[22] = { 0 };
static bool allowBasisExtrapolateRaw = true;
static void doBasisExtrapolate(const vector<int>& arr, map<int, SccProperty>& sccProp,
    ModelSummary& model, bool isLeft,
    double a1, double b1, double a2, double b2, double compRatios[22][2]) {
  int width = arr[1];
  int height = arr[0];
  int ledX, ledY;
  if (isLeft) {
    ledX = model.objects[2].leftX;
    ledY = model.objects[2].leftY;
  } else {
    ledX = model.objects[2].rightX;
    ledY = model.objects[2].rightY;
  }
  for (int i = 0; i < 22; ++i) {
    if (blockBasisExtrapolate[i]) continue;
    if (isLeft && model.objects[i].leftFound) continue;
    if (!isLeft && model.objects[i].rightFound) continue;
    int newX = (int) (compRatios[i][0] * a1 + compRatios[i][1] * a2);
    int newY = (int) (compRatios[i][0] * b1 + compRatios[i][1] * b2);
    newX += ledX;
    newY += ledY;
    double bestDist = height / 27; //(width + height);
    bestDist *= bestDist;
    SccProperty* bestProp = NULL;
    for (map<int, SccProperty>::iterator it = sccProp.begin(); it != sccProp.end(); ++it) {
      SccProperty& prop = it->second;
      int dx = it->second.x - newX;
      int dy = it->second.y - newY;
      double dist = dx * dx + dy * dy;
      if (dist < bestDist) {
        bestProp = &it->second;
        bestDist = dist;
      }
    }
    if (bestProp == NULL) {
      if (allowBasisExtrapolateRaw) {
        if (isLeft) {
          model.objects[i].leftX = newX;
          model.objects[i].leftY = newY;
          model.objects[i].leftFound = true;
          cerr << "Basis extrapolation matched no SCC, using raw for comp " << i << endl;
        } else {
          model.objects[i].rightX = newX;
          model.objects[i].rightY = newY;
          model.objects[i].rightFound = true;
          cerr << "Basis extrapolation matched no SCC, using raw for comp " << i << endl;
        }
      }
      continue;
    }
    if (i == 15 || i == 17 || i == 20) {
      // Toggles use topmost point instead of center of mass.
      newX = bestProp->yminx;
      newY = bestProp->ymin;
    } else {
      newX = bestProp->x;
      newY = bestProp->y;
    }
    if (isLeft) {
      model.objects[i].leftX = newX;
      model.objects[i].leftY = newY;
      model.objects[i].leftFound = true;
      cerr << "Basis extrapolation matched comp " << i << " to scc: " <<
          newX << "," << newY << endl;
    } else {
      model.objects[i].rightX = newX;
      model.objects[i].rightY = newY;
      model.objects[i].rightFound = true;
      cerr << "Basis extrapolation matched comp " << i << " to scc: " <<
          newX << "," << newY << endl;
    }
  }
}

static void basisExtrapolate(const vector<int>& arr, map<int, SccProperty>& sccProp,
    ModelSummary& model, bool isLeft, int gridComp, double compRatios[22][2]) {
  int width = arr[1];
  int height = arr[0];
  int ledX, ledY;
  int rockerX, rockerY;
  int gridX, gridY;
  if (isLeft) {
    ledX = model.objects[2].leftX;
    ledY = model.objects[2].leftY;
    rockerX = model.objects[3].leftX;
    rockerY = model.objects[3].leftY;
    gridX = model.objects[gridComp].leftX;
    gridY = model.objects[gridComp].leftY;

  } else {
    ledX = model.objects[2].rightX;
    ledY = model.objects[2].rightY;
    rockerX = model.objects[3].rightX;
    rockerY = model.objects[3].rightY;
    gridX = model.objects[gridComp].rightX;
    gridY = model.objects[gridComp].rightY;
  }
  double a1, b1;
  double a2, b2;
  {
    int dx = rockerX - ledX;
    int dy = rockerY - ledY;
    a1 = dx;
    b1 = dy;
  }
  {
    int dx = gridX - rockerX;
    int dy = gridY - rockerY;
    a2 = dx;
    b2 = dy;
  }
  doBasisExtrapolate(arr, sccProp, model, isLeft, a1, b1, a2, b2, compRatios);
}

static void basisExtrapolateWeak(const vector<int>& arr, map<int, SccProperty>& sccProp,
    ModelSummary& model, bool isLeft) {
  int width = arr[1];
  int height = arr[0];
  int ledX, ledY;
  int rockerX, rockerY;
  if (isLeft) {
    ledX = model.objects[2].leftX;
    ledY = model.objects[2].leftY;
    rockerX = model.objects[3].leftX;
    rockerY = model.objects[3].leftY;

  } else {
    ledX = model.objects[2].rightX;
    ledY = model.objects[2].rightY;
    rockerX = model.objects[3].rightX;
    rockerY = model.objects[3].rightY;
  }
  cerr << "Doing basisExtrapolateWeak with " << ledX << "," << ledY << " "
       << rockerY << "," << rockerY << endl;
  double a1, b1;
  double a2, b2;
  {
    int dx = rockerX - ledX;
    int dy = rockerY - ledY;
    a1 = dx;
    b1 = dy;
    a2 = b1;
    b2 = -a1;

    double length = dx * dx + dy * dy;
    double sinVal = dx / sqrt(length);
    //double theta = asin(sinVal) * 2;
    double theta = asin(sinVal);
    double a3 = cos(theta) * a2 - sin(theta) * b2;
    double b3 = sin(theta) * a2 + cos(theta) * b2;
    a2 = a3;
    b2 = b3;
  }

  // TODO: maybe include 20 and 21?
  double avgScale = 0;
  double avgTheta = 0;
  int numSamples = 0;
  // For each component, find the scale and rotation correction relative to the vector
  // made by that component to the rocker as compared to the predicted position, and
  // average the values in.
  vector<double> scaleSamples;
  vector<double> thetaSamples;
  for (int comp = 6; comp <= 18; ++comp) {
    int actX = -1, actY = -1;
    if (isLeft) {
      if (model.objects[comp].leftFound) {
        actX = model.objects[comp].leftX;
        actY = model.objects[comp].leftY;
      }
    } else {
      if (model.objects[comp].rightFound) {
        actX = model.objects[comp].rightX;
        actY = model.objects[comp].rightY;
      }
    }
    if (actX == -1 || actY == -1) continue;
    actX -= rockerX;
    actY -= rockerY;
    int predX, predY;
    if (isLeft) {
      predX = (int) (compRatiosOrthoLeft[comp][0] * a1 + compRatiosOrthoLeft[comp][1] * a2) + ledX;
      predY = (int) (compRatiosOrthoLeft[comp][0] * b1 + compRatiosOrthoLeft[comp][1] * b2) + ledY;
    } else {
      predX = (int) (compRatiosOrthoRight[comp][0] * a1 + compRatiosOrthoRight[comp][1] * a2) + ledX;
      predY = (int) (compRatiosOrthoRight[comp][0] * b1 + compRatiosOrthoRight[comp][1] * b2) + ledY;
    }
    predX -= rockerX;
    predY -= rockerY;
    double predLength = sqrt(predX * predX + predY * predY);
    double actLength = sqrt(actX * actX + actY * actY);
    double scaleConv = actLength / predLength;
    double crossProduct = predX * actY - predY * actX;
    crossProduct /= actLength;
    crossProduct /= predLength;
    double thetaConv = asin(crossProduct);
    avgScale += scaleConv - 1;
    avgTheta += thetaConv;
    ++numSamples;
    scaleSamples.push_back(scaleConv - 1);
    thetaSamples.push_back(thetaConv);
  }
  /*if (numSamples > 0) {
    avgScale = (scaleSamples[(numSamples - 1) / 2] + scaleSamples[numSamples / 2]) / 2;
    avgTheta = (thetaSamples[(numSamples - 1) / 2] + thetaSamples[numSamples / 2]) / 2;
    avgScale += 1;
    // Correct the normal vector by the average correction factors.
    double a3 = cos(avgTheta) * a2 - sin(avgTheta) * b2;
    double b3 = sin(avgTheta) * a2 + cos(avgTheta) * b2;
    a3 *= avgScale;
    b3 *= avgScale;
    a2 = a3;
    b2 = b3;
  }*/
  if (numSamples > 0) {
    avgScale /= numSamples;
    avgTheta /= numSamples;
    avgScale += 1;
    // Correct the normal vector by the average correction factors.
    double a3 = cos(avgTheta) * a2 - sin(avgTheta) * b2;
    double b3 = sin(avgTheta) * a2 + cos(avgTheta) * b2;
    a3 *= avgScale;
    b3 *= avgScale;
    a2 = a3;
    b2 = b3;
  }
  doBasisExtrapolate(arr, sccProp, model, isLeft, a1, b1, a2, b2,
      isLeft ? compRatiosOrthoLeft : compRatiosOrthoRight);
}

static void extractExtrapolation(ModelSummary& model, bool isLeft, int target,
    vector<vector<pair<int, int> > >& proposedGrid) {
  /*if (isLeft) {
    proposedGrid[target].push_back(make_pair(
        model.objects[target].leftX, model.objects[target].leftY));
    model.objects[target].leftFound = false;
  } else {
    proposedGrid[target].push_back(make_pair(
        model.objects[target].rightX, model.objects[target].rightY));
    model.objects[target].rightFound = false;
  }*/
}

static void centerObj(ModelSummary& model, bool isLeft, int target, int objA, int objB) {
  if (isLeft) {
    model.objects[target].leftX = (model.objects[objA].leftX + model.objects[objB].leftX) / 2;
    model.objects[target].leftY = (model.objects[objA].leftY + model.objects[objB].leftY) / 2;
    if (model.objects[target].state != "ON") {
      model.objects[target].state = "OFF";
    }
    model.objects[target].leftFound = true;
  } else {
    model.objects[target].rightX = (model.objects[objA].rightX + model.objects[objB].rightX) / 2;
    model.objects[target].rightY = (model.objects[objA].rightY + model.objects[objB].rightY) / 2;
    if (model.objects[target].state != "ON") {
      model.objects[target].state = "OFF";
    }
    model.objects[target].rightFound = true;
  }
}

static void extrapolateObj(ModelSummary& model, bool isLeft, int target, int objA, int objB) {
  if (isLeft) {
    model.objects[target].leftX = (model.objects[objA].leftX +
        2 * (model.objects[objB].leftX - model.objects[objA].leftX));
    model.objects[target].leftY = (model.objects[objA].leftY +
        2 * (model.objects[objB].leftY - model.objects[objA].leftY));
    if (model.objects[target].state != "ON") {
      model.objects[target].state = "OFF";
    }
    model.objects[target].leftFound = true;
  } else {
    model.objects[target].rightX = (model.objects[objA].rightX +
        2 * (model.objects[objB].rightX - model.objects[objA].rightX));
    model.objects[target].rightY = (model.objects[objA].rightY +
        2 * (model.objects[objB].rightY - model.objects[objA].rightY));
    if (model.objects[target].state != "ON") {
      model.objects[target].state = "OFF";
    }
    model.objects[target].rightFound = true;
  }
}

static void fillGrid(const vector<int>& arr,
    map<int, SccProperty>& sccProp,
    ModelSummary& model, bool isLeft) {
  // 6  7  8
  // 9  10 11
  // 12 13 14
  // Edge centers.
  bool foundExtrapolation = true;
  while (foundExtrapolation) {
    foundExtrapolation = false;
    vector<vector<pair<int, int> > > proposedGrid(22);
    for (int target = 9; target <= 11; ++target) {
      if (!gridHasObj(model, isLeft, target) &&
          gridHasObj(model, isLeft, target - 3) && gridHasObj(model, isLeft, target + 3)) {
        centerObj(model, isLeft, target, target - 3, target + 3);
        foundExtrapolation = true;
        extractExtrapolation(model, isLeft, target, proposedGrid);
      }
    }
    for (int target = 7; target <= 13; target += 3) {
      if (!gridHasObj(model, isLeft, target) &&
          gridHasObj(model, isLeft, target - 1) && gridHasObj(model, isLeft, target + 1)) {
        centerObj(model, isLeft, target, target - 1, target + 1);
        foundExtrapolation = true;
        extractExtrapolation(model, isLeft, target, proposedGrid);
      }
    }
    // Dead center.
    if (!gridHasObj(model, isLeft, 10) &&
        gridHasObj(model, isLeft, 6) && gridHasObj(model, isLeft, 14)) {
      centerObj(model, isLeft, 10, 6, 14);
      foundExtrapolation = true;
      extractExtrapolation(model, isLeft, 10, proposedGrid);
    }
    if (!gridHasObj(model, isLeft, 10) &&
        gridHasObj(model, isLeft, 12) && gridHasObj(model, isLeft, 8)) {
      centerObj(model, isLeft, 10, 12, 8);
      foundExtrapolation = true;
      extractExtrapolation(model, isLeft, 10, proposedGrid);
    }
    // Straight extrapolates.
    for (int target = 6; target <= 8; ++target) {
      if (!gridHasObj(model, isLeft, target) &&
          gridHasObj(model, isLeft, target + 3) && gridHasObj(model, isLeft, target + 6)) {
        extrapolateObj(model, isLeft, target, target + 6, target + 3);
        foundExtrapolation = true;
        extractExtrapolation(model, isLeft, target, proposedGrid);
      }
    }
    for (int target = 12; target <= 14; ++target) {
      if (!gridHasObj(model, isLeft, target) &&
          gridHasObj(model, isLeft, target - 3) && gridHasObj(model, isLeft, target - 6)) {
        extrapolateObj(model, isLeft, target, target - 6, target - 3);
        foundExtrapolation = true;
        extractExtrapolation(model, isLeft, target, proposedGrid);
      }
    }
    for (int target = 6; target <= 12; target += 3) {
      if (!gridHasObj(model, isLeft, target) &&
          gridHasObj(model, isLeft, target + 1) && gridHasObj(model, isLeft, target + 2)) {
        extrapolateObj(model, isLeft, target, target + 2, target + 1);
        foundExtrapolation = true;
        extractExtrapolation(model, isLeft, target, proposedGrid);
      }
    }
    for (int target = 8; target <= 14; target += 3) {
      if (!gridHasObj(model, isLeft, target) &&
          gridHasObj(model, isLeft, target - 1) && gridHasObj(model, isLeft, target - 2)) {
        extrapolateObj(model, isLeft, target, target - 2, target - 1);
        foundExtrapolation = true;
        extractExtrapolation(model, isLeft, target, proposedGrid);
      }
    }
    // Diagonal extrapolation.
    if (!gridHasObj(model, isLeft, 6) &&
        gridHasObj(model, isLeft, 10) && gridHasObj(model, isLeft, 14)) {
      extrapolateObj(model, isLeft, 6, 14, 10);
      foundExtrapolation = true;
      extractExtrapolation(model, isLeft, 6, proposedGrid);
    }
    if (!gridHasObj(model, isLeft, 8) &&
        gridHasObj(model, isLeft, 10) && gridHasObj(model, isLeft, 12)) {
      extrapolateObj(model, isLeft, 8, 12, 10);
      foundExtrapolation = true;
      extractExtrapolation(model, isLeft, 8, proposedGrid);
    }
    if (!gridHasObj(model, isLeft, 12) &&
        gridHasObj(model, isLeft, 10) && gridHasObj(model, isLeft, 8)) {
      extrapolateObj(model, isLeft, 12, 8, 10);
      foundExtrapolation = true;
      extractExtrapolation(model, isLeft, 12, proposedGrid);
    }
    if (!gridHasObj(model, isLeft, 14) &&
        gridHasObj(model, isLeft, 10) && gridHasObj(model, isLeft, 6)) {
      extrapolateObj(model, isLeft, 14, 6, 10);
      foundExtrapolation = true;
      extractExtrapolation(model, isLeft, 14, proposedGrid);
    }

    // Corner extrapolation.
    if (!gridHasObj(model, isLeft, 6) &&
        gridHasObj(model, isLeft, 7) &&
        gridHasObj(model, isLeft, 9) &&
        gridHasObj(model, isLeft, 10)) {
      if (isLeft) {
        model.objects[6].leftX = model.objects[7].leftX +
            (model.objects[9].leftX - model.objects[10].leftX);
        model.objects[6].leftY = model.objects[9].leftY +
            (model.objects[7].leftY - model.objects[10].leftY);
        if (model.objects[6].state != "ON") {
          model.objects[6].state = "OFF";
        }
        model.objects[6].leftFound = true;
      } else {
        model.objects[6].rightX = model.objects[7].rightX +
            (model.objects[9].rightX - model.objects[10].rightX);
        model.objects[6].rightY = model.objects[9].rightY +
            (model.objects[7].rightY - model.objects[10].rightY);
        if (model.objects[6].state != "ON") {
          model.objects[6].state = "OFF";
        }
        model.objects[6].rightFound = true;
      }
      foundExtrapolation = true;
      extractExtrapolation(model, isLeft, 6, proposedGrid);
    }
    /*if (foundExtrapolation) {
      for(int i = 0; i < 22; ++i) {
        if (proposedGrid[i].size() == 0) continue;
        int avgX = 0;
        int avgY = 0;
        for (int j = 0; j < proposedGrid[i].size(); ++j) {
          avgX += proposedGrid[i][j].first;
          avgY += proposedGrid[i][j].second;
        }
        avgX /= proposedGrid[i].size();
        avgY /= proposedGrid[i].size();
        if (isLeft) {
          model.objects[i].leftX = avgX;
          model.objects[i].leftY = avgY;
          model.objects[i].leftFound = true;
        } else {
          model.objects[i].rightX = avgX;
          model.objects[i].rightY = avgY;
          model.objects[i].rightFound = true;
        }
      }
    }*/

    // Straight-line L-extrapolation; only after all others have been exhausted.
    /*if (!foundExtrapolation) {
      if (!gridHasObj(model, isLeft, 13) &&
          gridHasObj(model, isLeft, 9) &&
          gridHasObj(model, isLeft, 10) &&
          gridHasObj(model, isLeft, 11)) {
        if (isLeft) {
          int dx = model.objects[11].leftX - model.objects[10].leftX;
          int dy = model.objects[11].leftY - model.objects[10].leftY;
          int newX = -dy + model.objects[10].leftX;
          int newY = dx + model.objects[10].leftY;
          SccProperty* bestProp = NULL;
          double bestDist = 99999;
          for (map<int, SccProperty>::iterator it = sccProp.begin(); it != sccProp.end(); ++it) {
            SccProperty& prop = it->second;
            if (newX >= prop.xmin && newX <= prop.xmax &&
                newY >= prop.ymin && newY <= prop.ymax) {
              double diffX = newX - prop.x;
              double diffY = newY - prop.y;
              if (diffX * diffX + diffY * diffY < bestDist) {
                bestDist = diffX * diffX + diffY * diffY;
                bestProp = &prop;
              }
            }
          }
          if (bestProp != NULL) {
            newX = bestProp->x;
            newY = bestProp->y;
          }
          model.objects[13].leftX = newX;
          model.objects[13].leftY = newY;
          if (model.objects[13].state != "ON") {
            model.objects[13].state = "OFF";
          }
          model.objects[13].leftFound = true;
        } else {
          int dx = model.objects[11].rightX - model.objects[10].rightX;
          int dy = model.objects[11].rightY - model.objects[10].rightY;
          int newX = -dy + model.objects[10].rightX;
          int newY = dx + model.objects[10].rightY;
          SccProperty* bestProp = NULL;
          double bestDist = 99999;
          for (map<int, SccProperty>::iterator it = sccProp.begin(); it != sccProp.end(); ++it) {
            SccProperty& prop = it->second;
            if (newX >= prop.xmin && newX <= prop.xmax &&
                newY >= prop.ymin && newY <= prop.ymax) {
              double diffX = newX - prop.x;
              double diffY = newY - prop.y;
              if (diffX * diffX + diffY * diffY < bestDist) {
                bestDist = diffX * diffX + diffY * diffY;
                bestProp = &prop;
              }
            }
          }
          if (bestProp != NULL) {
            newX = bestProp->x;
            newY = bestProp->y;
          }
          model.objects[13].rightX = newX;
          model.objects[13].rightY = newY;
          if (model.objects[13].state != "ON") {
            model.objects[13].state = "OFF";
          }
          model.objects[13].rightFound = true;
        }
        foundExtrapolation = true;
      }
    }*/
  }

  // Bottom three LEDs.
  if (!gridHasObj(model, isLeft, 16) &&
      gridHasObj(model, isLeft, 18) && gridHasObj(model, isLeft, 21)) {
    extrapolateObj(model, isLeft, 16, 21, 18);
  }
  if (!gridHasObj(model, isLeft, 21) &&
      gridHasObj(model, isLeft, 18) && gridHasObj(model, isLeft, 16)) {
    extrapolateObj(model, isLeft, 21, 16, 18);
  }
  if (!gridHasObj(model, isLeft, 18) &&
      gridHasObj(model, isLeft, 16) && gridHasObj(model, isLeft, 21)) {
    centerObj(model, isLeft, 18, 16, 21);
  }

//  for (int i = 6; i <= 14; ++i) {
//    int objX, objY;
//    if (isLeft) {
//      objX = model.objects[i].leftX;
//      objY = model.objects[i].leftY;
//    } else {
//      objX = model.objects[i].rightX;
//      objY = model.objects[i].rightY;
//    }
//    SccProperty* bestProp = NULL;
//    int bestDist = 400;  // Max distance we can try to adjust the position.
//    int64 avgX = 0;
//    int64 avgY = 0;
//    int numSamples = 0;
//    for (map<int, SccProperty>::iterator it = sccProp.begin(); it != sccProp.end(); ++it) {
//      SccProperty& prop = it->second;
//      int dx = objX - prop.x;
//      int dy = objY - prop.y;
//      int dist = dx * dx + dy * dy;
//      if (objX >= prop.xmin && objX <= prop.xmax &&
//          objY >= prop.ymin && objY <= prop.ymax &&
//          dist < bestDist) {
//        bestDist = dist;
//        bestProp = &prop;
//        avgX += prop.x;
//        avgY += prop.y;
//        ++numSamples;
//      }
//    }
//    //if (bestProp != NULL && (bestProp->x != objX || bestProp->y != objY)) {
//    //  cerr << "Extrapolated point moved adjusted from " << objX << "," << objY << " to "
//    //      << bestProp->x << "," << bestProp->y << " for comp " << i << endl;
//    if (numSamples > 0) {
//      avgX /= numSamples;
//      avgY /= numSamples;
//      if (isLeft) {
//        //model.objects[i].leftX = bestProp->x;
//        //model.objects[i].leftY = bestProp->y;
//        model.objects[i].leftX = avgX;
//        model.objects[i].leftY = avgY;
//      } else {
//        //model.objects[i].rightX = bestProp->x;
//        //model.objects[i].rightY = bestProp->y;
//        model.objects[i].rightX = avgX;
//        model.objects[i].rightY = avgY;
//      }
//    }
//  }
}

static void interpolateRockerLeds(
    const vector<int>& arr, vector<vector<int> >& scc, map<int, SccProperty>& sccProp,
    ModelSummary& model, bool isLeft) {
  int width = arr[1];
  int height = arr[0];
  int ledX, ledY;
  int rockerX, rockerY;
  if (isLeft) {
    ledX = model.objects[2].leftX;
    ledY = model.objects[2].leftY;
    rockerX = model.objects[3].leftX;
    rockerY = model.objects[3].leftY;
  } else {
    ledX = model.objects[2].rightX;
    ledY = model.objects[2].rightY;
    rockerX = model.objects[3].rightX;
    rockerY = model.objects[3].rightY;
  }
  int vecX = rockerX - ledX;
  int vecY = rockerY - ledY;
  int ledTopX = vecX * 73 / 100 + ledX;
  if (ledTopX < 0) ledTopX = 0;
  if (ledTopX >= width) ledTopX = width - 1;
  int ledTopY = vecY * 73 / 100 + ledY;
  if (ledTopY < 0) ledTopY = 0;
  if (ledTopY >= height) ledTopY = height - 1;
  int ledBottomX = vecX * 122 / 100 + ledX;
  if (ledBottomX < 0) ledBottomX = 0;
  if (ledBottomX >= width) ledBottomX = width - 1;
  int ledBottomY = vecY * 122 / 100 + ledY;
  if (ledBottomY < 0) ledBottomY = 0;
  if (ledBottomY >= height) ledBottomY = height - 1;

  int64 maxDist = (width + height) / 2 * 25L / 1000L;
  maxDist *= maxDist;
  SccProperty* topProp = NULL;
  int64 bestDist = maxDist;
  for (map<int, SccProperty>::iterator it = sccProp.begin(); it != sccProp.end(); ++it) {
    SccProperty& prop = it->second;
    int dx = prop.x - ledTopX;
    int dy = prop.y - ledTopY;
    if (dx * dx + dy * dy < bestDist) {
      bestDist = dx * dx + dy * dy;
      topProp = &prop;
    }
  }

  // Just naive position-guessing, so don't mark as found just yet.
  if (isLeft) {
    // Only update if not found/matched-up already.
    if (!model.objects[4].leftFound) {
      if (topProp != NULL) {
        model.objects[4].leftX = topProp->x;
        model.objects[4].leftY = topProp->y;
        model.objects[4].leftProp = topProp;
      } else {
        model.objects[4].leftX = ledTopX;
        model.objects[4].leftY = ledTopY;
      }
      model.objects[4].leftFound = true;
    }
    if (!model.objects[5].leftFound) {
      model.objects[5].leftX = ledBottomX;
      model.objects[5].leftY = ledBottomY;
      model.objects[5].leftFound = true;
    }
  } else {
    if (!model.objects[4].rightFound) {
      if (topProp != NULL) {
        model.objects[4].rightX = topProp->x;
        model.objects[4].rightY = topProp->y;
        model.objects[4].rightProp = topProp;
      } else {
        model.objects[4].rightX = ledTopX;
        model.objects[4].rightY = ledTopY;
      }
      model.objects[4].rightFound = true;
    }
    if (!model.objects[5].rightFound) {
      model.objects[5].rightX = ledBottomX;
      model.objects[5].rightY = ledBottomY;
      model.objects[5].rightFound = true;
    }
  }
}

static void estimateGrid(double grid[9][2], ModelSummary& model, bool isLeft) {
  int width = model.width;
  int height = model.height;
  int ledX, ledY, rockerX, rockerY;
  if (isLeft) {
    if (!model.objects[2].leftFound || !model.objects[3].leftFound) return;
    ledX = model.objects[2].leftX;
    ledY = model.objects[2].leftY;
    rockerX = model.objects[3].leftX;
    rockerY = model.objects[3].leftY;
  } else {
    if (!model.objects[2].rightFound || !model.objects[3].rightFound) return;
    ledX = model.objects[2].rightX;
    ledY = model.objects[2].rightY;
    rockerX = model.objects[3].rightX;
    rockerY = model.objects[3].rightY;
  }
  double a1 = rockerX - ledX;
  double b1 = rockerY - ledY;
  double a2 = b1;
  double b2 = -a1;
  for (int i = 0; i < 9; ++i) {
    int a3 = (int) (grid[i][0] * a1 + grid[i][1] * a2);
    int b3 = (int) (grid[i][0] * b1 + grid[i][1] * b2);
    int x = a3 + ledX;
    int y = b3 + ledY;
    if (x < 0) x = 0;
    if (x >= width) x = width - 1;
    if (y < 0) y = 0;
    if (y >= height) y = height - 1;
    if (isLeft) {
      model.objects[i + 6].leftX = x;
      model.objects[i + 6].leftY = y;
    } else {
      model.objects[i + 6].rightX = x;
      model.objects[i + 6].rightY = y;
    }
  }
}

static void reconcileMissing(ModelSummary& model) {
  // Extract all the ones present in both.
  int dx = 0;
  int dy = 0;
  int count = 0;
  // For now, the red switch is the most reliable.
  //for (int i = 0; i < 22; ++i) {
  // Panel LED and rocker are most reliable.
  for (int i = 2; i <= 3; ++i) {
    if (model.objects[i].leftFound && model.objects[i].rightFound) {
      dx += model.objects[i].leftX - model.objects[i].rightX;
      dy += model.objects[i].leftY - model.objects[i].rightY;
      ++count;
    }
  }
  if (count == 0) {
    cerr << "No objects present in both; skipping reconciliation." << endl;
    return;
  }
  dx /= count;
  dy /= count;
  for (int i = 0; i < 22; ++i) {
    if (!model.objects[i].leftFound && model.objects[i].rightFound) {
      cerr << "Found right but not left for object " << i << "; reconciling." << endl;
      model.objects[i].leftX = model.objects[i].rightX + dx;
      model.objects[i].leftY = model.objects[i].rightY + dy;
      model.objects[i].leftFound = true;
    } else if (model.objects[i].leftFound && !model.objects[i].rightFound) {
      cerr << "Found left but not right for object " << i << "; reconciling." << endl;
      model.objects[i].rightX = model.objects[i].leftX - dx;
      model.objects[i].rightY = model.objects[i].leftY - dy;
      model.objects[i].rightFound = true;
    }
  }
}

#ifdef DEBUG
int64 CurrentTimeMillis();
#endif

vector<string> RobonautEye::recognizeObjects(
    const vector<int>& leftEyeImageOrig, const vector<int>& rightEyeImageOrig) {
  string ans[] = { "UP", "0", "0", "0", "0",        // PANEL_POWER_SWITCH (0) (0)
                   "UP", "0", "0", "0", "0",        // PANEL_POWER_COVER (5) (1)
                   "ON", "0", "0", "0", "0",        // PANEL_POWER_LED (10) (2)
                   "CENTER", "0", "0", "0", "0",    // A01_ROCKER_SWITCH (15) (3)
                   "OFF", "0", "0", "0", "0",       // A01_ROCKER_LED_TOP (20) (4)
                   "OFF", "0", "0", "0", "0",       // A01_ROCKER_LED_BOTTOM (25) (5)
                   "OFF", "0", "0", "0", "0",       // A02_LED_NUM_PAD_A1 (30) (6)
                   "OFF", "0", "0", "0", "0",       // A02_LED_NUM_PAD_A2 (35) (7)
                   "OFF", "0", "0", "0", "0",       // A02_LED_NUM_PAD_A3 (40) (8)
                   "OFF", "0", "0", "0", "0",       // A02_LED_NUM_PAD_B1 (45) (9)
                   "OFF", "0", "0", "0", "0",       // A02_LED_NUM_PAD_B2 (50) (10)
                   "OFF", "0", "0", "0", "0",       // A02_LED_NUM_PAD_B3 (55) (11)
                   "OFF", "0", "0", "0", "0",       // A02_LED_NUM_PAD_C1 (60) (12)
                   "OFF", "0", "0", "0", "0",       // A02_LED_NUM_PAD_C2 (65) (13)
                   "OFF", "0", "0", "0", "0",       // A02_LED_NUM_PAD_C3 (70) (14)
                   "DOWN", "0", "0", "0", "0",      // A03_TOGGLE (75) (15)
                   "OFF", "0", "0", "0", "0",       // A03_LED (80) (16)
                   "CENTER", "0", "0", "0", "0",    // A04_TOGGLE (85) (17)
                   "OFF", "0", "0", "0", "0",       // A04_LED_TOP (90) (18)
                   "OFF", "0", "0", "0", "0",       // A04_LED_BOTTOM (95) (19)
                   "DOWN", "0", "0", "0", "0",      // A05_TOGGLE (100) (20)
                   "OFF", "0", "0", "0", "0" };     // A05_LED (105) (21)

  /*int origHeight = leftEyeImageOrig[0];
  int origWidth = leftEyeImageOrig[1];

  int newHeight = origHeight;
  int newWidth = origWidth;
  vector<int> leftEyeImage(2 + newHeight * newWidth);
  vector<int> rightEyeImage(2 + newHeight * newWidth);
  leftEyeImage[0] = newHeight;
  leftEyeImage[1] = newWidth;
  rightEyeImage[0] = newHeight;
  rightEyeImage[1] = newWidth;
  for (int y = 0; y < newHeight; ++y) {
    for (int x = 0; x < newWidth; ++x) {
      if (x == 0 || x == newWidth - 1 || y == 0 || y == newHeight - 1) {
        leftEyeImage[2 + x + y * newWidth] = leftEyeImageOrig[2 + x + y * origWidth];
        rightEyeImage[2 + x + y * newWidth] = rightEyeImageOrig[2 + x + y * origWidth];
      } else {
        {
          int rs = 0, gs = 0, bs = 0;
          for (int dx = x - 1; dx <= x + 1; ++dx) {
            for (int dy = y - 1; dy <= y + 1; ++dy) {
              int rgb = leftEyeImageOrig[2 + dx + dy * origWidth];
              int r = (rgb >> 16) & 0x00ff;
              int g = (rgb >> 8) & 0x00ff;
              int b = (rgb) & 0x00ff;
              rs += r;
              gs += g;
              bs += b;
            }
          }
          rs /= 9;
          gs /= 9;
          bs /= 9;
          leftEyeImage[2 + x + y * newWidth] = (rs << 16) | (gs << 8) | bs;
        }
        {
          int rs = 0, gs = 0, bs = 0;
          for (int dx = x - 1; dx <= x + 1; ++dx) {
            for (int dy = y - 1; dy <= y + 1; ++dy) {
              int rgb = rightEyeImageOrig[2 + dx + dy * origWidth];
              int r = (rgb >> 16) & 0x00ff;
              int g = (rgb >> 8) & 0x00ff;
              int b = (rgb) & 0x00ff;
              rs += r;
              gs += g;
              bs += b;
            }
          }
          rs /= 9;
          gs /= 9;
          bs /= 9;
          rightEyeImage[2 + x + y * newWidth] = (rs << 16) | (gs << 8) | bs;
        }
      }
    }
  }*/
  const vector<int>& leftEyeImage = leftEyeImageOrig;
  const vector<int>& rightEyeImage = rightEyeImageOrig;

  int height = leftEyeImage[0];
  int width = leftEyeImage[1];
  bool isSmall = false;

  if (height > 2000 && width > 2000) {
    for (int i = 0; i < 22; ++i) {
      compRatiosOrthoLeft[i][0] = compRatiosOrthoLargeLeft[i][0];
      compRatiosOrthoLeft[i][1] = compRatiosOrthoLargeLeft[i][1];
      compRatiosOrthoRight[i][0] = compRatiosOrthoLargeRight[i][0];
      compRatiosOrthoRight[i][1] = compRatiosOrthoLargeRight[i][1];
    }
    UPPER_NIBLACK = 1.5;
  } else if (height < 1700 && width < 1700) {
    isSmall = true;
    for (int i = 0; i < 22; ++i) {
      compRatiosOrthoLeft[i][0] = compRatiosOrthoSmallLeft[i][0];
      compRatiosOrthoLeft[i][1] = compRatiosOrthoSmallLeft[i][1];
      compRatiosOrthoRight[i][0] = compRatiosOrthoSmallRight[i][0];
      compRatiosOrthoRight[i][1] = compRatiosOrthoSmallRight[i][1];
    }
    UPPER_NIBLACK = 1.5;
  }

#ifdef DEBUG
  int64 startTime = CurrentTimeMillis();
#endif

  SumImage sumLeft(leftEyeImage);
  SumImage sumRight(rightEyeImage);

#ifdef DEBUG
  cerr << "Took " << CurrentTimeMillis() - startTime << " ms to sumImage" << endl;
  startTime = CurrentTimeMillis();
#endif
  int kernel = width / 5;
  doNiblackMulti(leftEyeImage, sumLeft, kernel);
  doNiblackMulti(rightEyeImage, sumRight, kernel);
#ifdef DEBUG
  cerr << "Took " << CurrentTimeMillis() - startTime << " ms to niblack" << endl;
  startTime = CurrentTimeMillis();
  //if (startTime > 42) exit(1);
#endif

  vector<vector<int> > sccLeft(width);
  vector<vector<int> > sccRight(width);
  for (int i = 0; i < width; ++i) {
    sccLeft[i].resize(height);
    sccRight[i].resize(height);
  }
  map<int, SccProperty> sccPropLeft;
  map<int, SccProperty> sccPropLargeLeft;
  map<int, SccProperty> sccPropRight;
  map<int, SccProperty> sccPropLargeRight;
  findScc(leftEyeImage, sccLeft, sccPropLeft, sccPropLargeLeft, sumLeft);
  findScc(rightEyeImage, sccRight, sccPropRight, sccPropLargeRight, sumRight);
#ifdef DEBUG
  cerr << "Took " << CurrentTimeMillis() - startTime << " ms to SCC" << endl;
  startTime = CurrentTimeMillis();
#endif
  bool foundLedLeft = false;
  SccProperty* candidatePanelLedLeft = NULL;
  {
    vector<SccProperty*> allGreenLeds;
    set<int> allGreenLedsNum;
    for (map<int, SccProperty>::iterator it = sccPropLeft.begin(); it != sccPropLeft.end(); ++it) {
      SccProperty& prop = it->second;
      if (prop.isGreenLed && allGreenLedsNum.find(prop.propNum) == allGreenLedsNum.end()) {
        allGreenLedsNum.insert(prop.propNum);
        allGreenLeds.push_back(&prop);
      }
    }
    if (allGreenLeds.size() > 1) {
      sort(allGreenLeds.begin(), allGreenLeds.end(), compareSccByY);
      if (allGreenLeds[1]->y - allGreenLeds[0]->y > height / 15) {
        foundLedLeft = true;
        candidatePanelLedLeft = allGreenLeds[0];
      }
    } else if (allGreenLeds.size() == 1) {
      foundLedLeft = true;
      candidatePanelLedLeft = allGreenLeds[0];
    }
  }
  bool foundLedRight = false;
  SccProperty* candidatePanelLedRight = NULL;
  {
    vector<SccProperty*> allGreenLeds;
    set<int> allGreenLedsNum;
    for (map<int, SccProperty>::iterator it = sccPropRight.begin(); it != sccPropRight.end(); ++it) {
      SccProperty& prop = it->second;
      if (prop.isGreenLed && allGreenLedsNum.find(prop.propNum) == allGreenLedsNum.end()) {
        allGreenLedsNum.insert(prop.propNum);
        allGreenLeds.push_back(&prop);
      }
    }
    if (allGreenLeds.size() > 1) {
      sort(allGreenLeds.begin(), allGreenLeds.end(), compareSccByY);
      if (allGreenLeds[1]->y - allGreenLeds[0]->y > height / 15) {
        foundLedRight = true;
        candidatePanelLedRight = allGreenLeds[0];
      }
    } else if (allGreenLeds.size() == 1) {
      foundLedRight = true;
      candidatePanelLedRight = allGreenLeds[0];
    }
  }
  ModelSummary model(width, height);
  for (int i = 0; i < 22; ++i) {
    model.objects[i].state = ans[5 * i];
  }

  {
    if (foundLedLeft && candidatePanelLedLeft != NULL) {
      findRed(leftEyeImage, model, true,
          sccLeft, sccPropLeft, sccPropLargeLeft,
          candidatePanelLedLeft->x - width / 10, candidatePanelLedLeft->y - height / 8,
          candidatePanelLedLeft->x + width / 80, candidatePanelLedLeft->y + height / 8,
          .05);
      if (model.objects[1].leftFound) {
        // If we found it, commit both the cover as well as the led.
        model.objects[2].leftFound = true;
        model.objects[2].leftX = candidatePanelLedLeft->x;
        model.objects[2].leftY = candidatePanelLedLeft->y;
        model.objects[2].leftProp = candidatePanelLedLeft;
      }
    }
    if (foundLedRight && candidatePanelLedRight != NULL) {
      findRed(rightEyeImage, model, false,
          sccRight, sccPropRight, sccPropLargeRight,
          candidatePanelLedRight->x - width / 10, candidatePanelLedRight->y - height / 8,
          candidatePanelLedRight->x + width / 80, candidatePanelLedRight->y + height / 8,
          .05);
      if (model.objects[1].rightFound) {
        // If we found it, commit both the cover as well as the led.
        model.objects[2].rightFound = true;
        model.objects[2].rightX = candidatePanelLedRight->x;
        model.objects[2].rightY = candidatePanelLedRight->y;
        model.objects[2].rightProp = candidatePanelLedRight;
      }
    }

    // Looking for the PANEL_POWER_COVER.
    if (!model.objects[1].leftFound) {
      findRed(leftEyeImage, model, true,
          sccLeft, sccPropLeft, sccPropLargeLeft,
          width / 6, 0, width, height * 3 / 4, 0);
    }
    if (!model.objects[1].rightFound) {
      findRed(rightEyeImage, model, false,
          sccRight, sccPropRight, sccPropLargeRight,
          width / 6, 0, width, height * 3 / 4, 0);
    }
#ifdef DEBUG
  cerr << "Took " << CurrentTimeMillis() - startTime << " ms to look for cover" << endl;
  startTime = CurrentTimeMillis();
#endif

    // Guess the switch is straight down.
    model.objects[0].leftX = model.objects[1].leftX;
    model.objects[0].leftY = model.objects[1].leftY + height /20;
    model.objects[0].rightX = model.objects[1].rightX;
    model.objects[0].rightY = model.objects[1].rightY + height /20;

    // PANEL_POWER_LED.
    if (!model.objects[2].leftFound) {
      findPanelLed(leftEyeImage, sccLeft, sccPropLeft, sccPropLargeLeft, sumLeft, model, true);
    }
    if (!model.objects[2].rightFound) {
      findPanelLed(rightEyeImage, sccRight, sccPropRight, sccPropLargeRight, sumRight, model, false);
    }

    // Panel state.
    bool panelIsProbablyOff = false;
    if (!foundLedLeft && !foundLedRight &&
        !model.objects[2].leftFound && !model.objects[2].rightFound) {
      panelIsProbablyOff = true;
      ans[10] = "OFF";
      ans[0] = "DOWN";
    }

    // A01_ROCKER_SWITCH
    findRockerSwitch(leftEyeImage, sccLeft, sccPropLeft, sccPropLargeLeft, sumLeft, model, true);
    findRockerSwitch(rightEyeImage, sccRight, sccPropRight, sccPropLargeRight, sumRight, model, false);

#ifdef DEBUG
  cerr << "Took " << CurrentTimeMillis() - startTime << " ms to look for LED/switch" << endl;
  startTime = CurrentTimeMillis();
#endif

    // Rocker might have a pretty good guess. This solves Lab2/MPUI.
    /*if (!model.objects[3].leftFound) {
      vector<SccProperty*> allPossibleRockerSwitches;
      for (map<int, SccProperty>::iterator it = sccPropLeft.begin(); it != sccPropLeft.end(); ++it) {
        SccProperty& prop = it->second;
        if (prop.isRockerSwitch) {
          allPossibleRockerSwitches.push_back(&prop);
        }
      }
      if (allPossibleRockerSwitches.size() == 1) {
        SccProperty& prop = *allPossibleRockerSwitches[0];
        model.objects[3].leftX = prop.x;
        model.objects[3].leftY = prop.y;
        model.objects[3].leftFound = true;

        if (!model.objects[2].leftFound) {
          model.objects[2].leftX = model.objects[3].leftX;
          model.objects[2].leftY = model.objects[3].leftY - height / 5;
          model.objects[2].leftFound = true;
        }
      } 
    }
    if (!model.objects[3].rightFound) {
      vector<SccProperty*> allPossibleRockerSwitches;
      for (map<int, SccProperty>::iterator it = sccPropRight.begin(); it != sccPropRight.end(); ++it) {
        SccProperty& prop = it->second;
        if (prop.isRockerSwitch) {
          allPossibleRockerSwitches.push_back(&prop);
        }
      }
      if (allPossibleRockerSwitches.size() == 1) {
        SccProperty& prop = *allPossibleRockerSwitches[0];
        model.objects[3].rightX = prop.x;
        model.objects[3].rightY = prop.y;
        model.objects[3].rightFound = true;

        if (!model.objects[2].rightFound) {
          model.objects[2].rightX = model.objects[3].rightX;
          model.objects[2].rightY = model.objects[3].rightY - height / 5;
          model.objects[2].rightFound = true;
        }
      } 
    }*/

    if (!model.objects[2].leftFound && !model.objects[3].leftFound &&
        (model.objects[2].rightFound || model.objects[3].rightFound)) {
      int dx = model.objects[1].leftX - (model.objects[1].rightX + width / 4);
      int dy = model.objects[1].leftY - model.objects[1].rightY;
      if (dx * dx + dy * dy > width * width / 36) {
        // Special case; trust the right eye for red cover in this situation.
        model.objects[1].leftX = model.objects[1].rightX + width / 4;
        if (model.objects[1].leftX >= width) {
          model.objects[1].leftX = width - 1;
        }
        model.objects[1].leftY = model.objects[1].rightY;
        model.objects[1].leftProp = NULL;

        // Redo all the things which depended on the red switch.
        findPanelLed(leftEyeImage, sccLeft, sccPropLeft, sccPropLargeLeft, sumLeft, model, true);
        findRockerSwitch(leftEyeImage, sccLeft, sccPropLeft, sccPropLargeLeft, sumLeft, model, true);
      }
    }
    if (!model.objects[2].rightFound && !model.objects[3].rightFound &&
        (model.objects[2].leftFound || model.objects[3].leftFound)) {
      int dx = model.objects[1].rightX - (model.objects[1].leftX - width / 4);
      int dy = model.objects[1].rightY - model.objects[1].leftY;
      if (dx * dx + dy * dy > width * width / 36) {
        // Special case; trust the left eye for red cover in this situation.
        model.objects[1].rightX = model.objects[1].leftX - width / 4;
        if (model.objects[1].rightX < 0) {
          model.objects[1].rightX = 0;
        }
        model.objects[1].rightY = model.objects[1].leftY;
        model.objects[1].rightProp = NULL;

        // Redo all the things which depended on the red switch.
        findPanelLed(rightEyeImage, sccRight, sccPropRight, sccPropLargeRight, sumRight, model, false);
        findRockerSwitch(rightEyeImage, sccRight, sccPropRight, sccPropLargeRight, sumRight, model, false);
      }
    }

    // Wild guess for 2.
    if (model.objects[1].leftFound && !model.objects[2].leftFound) {
      int predX;
      if (isSmall) {
        predX = model.objects[1].leftX + width / 40;
      } else {
        predX = model.objects[1].leftX + width / 33;
      }
      int predY = model.objects[1].leftY;
      SccProperty* bestProp = NULL;
      if (panelIsProbablyOff) {
        int ymin = 0;
        int ymax = height;
        if (model.objects[1].leftProp != NULL) {
          ymin = model.objects[1].leftProp->ymin;
          ymax = model.objects[1].leftProp->ymax;
        }
        int64 bestDist = (width + height) / 2 / 20;
        bestDist *= bestDist;
        for (map<int, SccProperty>::iterator it = sccPropLeft.begin(); it != sccPropLeft.end(); ++it) {
          SccProperty& prop = it->second;
          if (prop.x > model.objects[1].leftX &&
              prop.y >= ymin && prop.y < ymax &&
              prop.xmax - prop.xmin < width / 40 &&
              prop.ymax - prop.ymin < height / 30 &&
              max(prop.r, max(prop.g, prop.b)) < 100 &&
              canReach(width, height, prop.x, prop.y, predX, predY, sccLeft, sccPropLargeLeft)) {
            int dx = prop.x - predX;
            int dy = prop.y - predY;
            if (dx * dx + dy * dy < bestDist) {
              bestDist = dx * dx + dy * dy;
              bestProp = &prop;
            }
          }
        }
      }
      if (bestProp != NULL) {
        predX = bestProp->x;
        predY = bestProp->y;
      }
      model.objects[2].leftX = predX;
      model.objects[2].leftY = predY;
      model.objects[2].leftProp = bestProp;
      model.objects[2].leftFound = true;
    }
    if (model.objects[1].rightFound && !model.objects[2].rightFound) {
      int predX;
      if (isSmall) {
        predX = model.objects[1].rightX + width / 40;
      } else {
        predX = model.objects[1].rightX + width / 33;
      }
      int predY = model.objects[1].rightY;
      SccProperty* bestProp = NULL;
      if (panelIsProbablyOff) {
        int ymin = 0;
        int ymax = height;
        if (model.objects[1].rightProp != NULL) {
          ymin = model.objects[1].rightProp->ymin;
          ymax = model.objects[1].rightProp->ymax;
        }
        int64 bestDist = (width + height) / 2 / 20;
        bestDist *= bestDist;
        for (map<int, SccProperty>::iterator it = sccPropRight.begin(); it != sccPropRight.end(); ++it) {
          SccProperty& prop = it->second;
          if (prop.x > model.objects[1].rightX &&
              prop.y >= ymin && prop.y < ymax &&
              prop.xmax - prop.xmin < width / 40 &&
              prop.ymax - prop.ymin < height / 30 &&
              max(prop.r, max(prop.g, prop.b)) < 100 &&
              canReach(width, height, prop.x, prop.y, predX, predY, sccRight, sccPropLargeRight)) {
            int dx = prop.x - predX;
            int dy = prop.y - predY;
            if (dx * dx + dy * dy < bestDist) {
              bestDist = dx * dx + dy * dy;
              bestProp = &prop;
            }
          }
        }
      }
      if (bestProp != NULL) {
        predX = bestProp->x;
        predY = bestProp->y;
      }
      model.objects[2].rightX = predX;
      model.objects[2].rightY = predY;
      model.objects[2].rightProp = bestProp;
      model.objects[2].rightFound = true;
    }

    // Wild guess for 3.
    if (model.objects[2].leftFound && !model.objects[3].leftFound) {
      model.objects[3].leftX = model.objects[2].leftX;
      model.objects[3].leftY = model.objects[2].leftY + height / 5;
      model.objects[3].leftFound = true;
    }
    if (model.objects[2].rightFound && !model.objects[3].rightFound) {
      model.objects[3].rightX = model.objects[2].rightX;
      model.objects[3].rightY = model.objects[2].rightY + height / 5;
      model.objects[3].rightFound = true;
    }

    // Possibly all leds.
    /*if (model.objects[2].leftFound && model.objects[3].leftFound) {
      dfsMatchStrict = true;
      if (!matchSccToBothLedsDfsRect(
              leftEyeImage, sccLeft, sccPropLeft, sccPropLargeLeft, model, true)) {
        dfsMatchStrict = false;
        if (model.objects[2].leftFound && model.objects[3].leftFound) {
          matchSccToGreenLedsDfsRect(
              leftEyeImage, sccLeft, sccPropLeft, sccPropLargeLeft, model, true);
        }
        if (model.objects[2].leftFound && model.objects[3].leftFound) {
          matchSccToBlueLedsDfsRect(
              leftEyeImage, sccLeft, sccPropLeft, sccPropLargeLeft, model, true);
        }
      }
    }
    if (model.objects[2].rightFound && model.objects[3].rightFound) {
      dfsMatchStrict = true;
      if (!matchSccToBothLedsDfsRect(
              rightEyeImage, sccRight, sccPropRight, sccPropLargeRight, model, false)) {
        dfsMatchStrict = false;
        if (model.objects[2].rightFound && model.objects[3].rightFound) {
          matchSccToGreenLedsDfsRect(
              rightEyeImage, sccRight, sccPropRight, sccPropLargeRight, model, false);
        }
        if (model.objects[2].rightFound && model.objects[3].rightFound) {
          matchSccToBlueLedsDfsRect(
              rightEyeImage, sccRight, sccPropRight, sccPropLargeRight, model, false);
        }
      }
    }
    dfsMatchStrict = false;*/
    if (model.objects[2].leftFound && model.objects[3].leftFound) {
      panelLedForDfs = model.objects[2].leftProp;
      rockerForDfs = model.objects[3].leftProp;
      matchSccToGreenLedsDfsRect(
          leftEyeImage, sccLeft, sccPropLeft, sccPropLargeLeft, model, true);
    }
    if (model.objects[2].rightFound && model.objects[3].rightFound) {
      panelLedForDfs = model.objects[2].rightProp;
      rockerForDfs = model.objects[3].rightProp;
      matchSccToGreenLedsDfsRect(
          rightEyeImage, sccRight, sccPropRight, sccPropLargeRight, model, false);
    }
    if (model.objects[2].leftFound && model.objects[3].leftFound) {
      panelLedForDfs = model.objects[2].leftProp;
      rockerForDfs = model.objects[3].leftProp;
      matchSccToBlueLedsDfsRect(
          leftEyeImage, sccLeft, sccPropLeft, sccPropLargeLeft, model, true);
    }
    if (model.objects[2].rightFound && model.objects[3].rightFound) {
      panelLedForDfs = model.objects[2].rightProp;
      rockerForDfs = model.objects[3].rightProp;
      matchSccToBlueLedsDfsRect(
          rightEyeImage, sccRight, sccPropRight, sccPropLargeRight, model, false);
    }
    panelLedForDfs = NULL;
    rockerForDfs = NULL;
#ifdef DEBUG
  cerr << "Took " << CurrentTimeMillis() - startTime << " ms to look for dfs match LEDs" << endl;
  startTime = CurrentTimeMillis();
#endif

    fillGrid(leftEyeImage, sccPropLeft, model, true);
    fillGrid(rightEyeImage, sccPropRight, model, false);

    // A01_ROCKER_LED_TOP and A01_ROCKER_LED_BOTTOM
    if (model.objects[2].leftFound && model.objects[3].leftFound) {
      interpolateRockerLeds(leftEyeImage, sccLeft, sccPropLeft, model, true);
    }
    if (model.objects[2].rightFound && model.objects[3].rightFound) {
      interpolateRockerLeds(rightEyeImage, sccRight, sccPropRight, model, false);
    }

    // For now, correct for red switch position just for scoring purposes.
    bool leftCoverDown = false;
    if (model.objects[1].leftProp != NULL) {
      SccProperty* prop = model.objects[1].leftProp;
      bool isTall = false;
      if (prop->xmax - prop->xmin < (prop->ymax - prop->ymin) * 4 / 5) {
        isTall = true;
      }
      bool isClose = false;
      for (int x = max(prop->x - 10, 0); x < min(prop->x + 10, width); ++x) {
        for (int y = max(prop->ymax - 10, 0); y < min(prop->ymax + 10, height); ++y) {
          map<int, SccProperty>::iterator it = sccPropLeft.find(sccLeft[x][y]);
          if (it != sccPropLeft.end() &&
              it->second.x == prop->x &&
              it->second.y == prop->y &&
              it->second.xmin == prop->xmin &&
              it->second.xmax == prop->xmax &&
              it->second.ymin == prop->ymin &&
              it->second.ymax == prop->ymax) {
            isClose = true;
            break;
          }
        }
        if (isClose) break;
      }
      int markX, markY;
      if (isTall && isClose) {
        cerr << "Cover appears to be down for calculating its position" << endl;
        markX = (prop->ymaxx + prop->yminx) / 2;
        markY = prop->ymax;
      } else {
        cerr << "Cover appears to be up for calculating its position" << endl;
        if (prop->xminy < prop->xmaxy) {
          markX = prop->xmin + (prop->yminx - prop->xmin) / 4;
          markY = prop->xminy;
        } else {
          markX = prop->xmax - (prop->xmax - prop->yminx) / 4;
          markY = prop->xmaxy;
        }
      }
      model.objects[1].leftX = markX;
      model.objects[1].leftY = markY;
      leftCoverDown = isTall && isClose;
    }
    bool rightCoverDown = false;
    if (model.objects[1].rightProp != NULL) {
      SccProperty* prop = model.objects[1].rightProp;
      bool isTall = false;
      if (prop->xmax - prop->xmin < (prop->ymax - prop->ymin) * 4 / 5) {
        isTall = true;
      }
      bool isClose = false;
      for (int x = max(prop->x - 10, 0); x < min(prop->x + 10, width); ++x) {
        for (int y = max(prop->ymax - 10, 0); y < min(prop->ymax + 10, height); ++y) {
          map<int, SccProperty>::iterator it = sccPropRight.find(sccRight[x][y]);
          if (it != sccPropRight.end() &&
              it->second.x == prop->x &&
              it->second.y == prop->y &&
              it->second.xmin == prop->xmin &&
              it->second.xmax == prop->xmax &&
              it->second.ymin == prop->ymin &&
              it->second.ymax == prop->ymax) {
            isClose = true;
            break;
          }
        }
        if (isClose) break;
      }
      int markX, markY;
      if (isTall && isClose) {
        cerr << "Cover appears to be down for calculating its position" << endl;
        markX = (prop->ymaxx + prop->yminx) / 2;
        markY = prop->ymax;
      } else {
        cerr << "Cover appears to be up for calculating its position" << endl;
        if (prop->xminy < prop->xmaxy) {
          markX = prop->xmin + (prop->yminx - prop->xmin) / 4;
          markY = prop->xminy;
        } else {
          markX = prop->xmax - (prop->xmax - prop->yminx) / 4;
          markY = prop->xmaxy;
        }
      }
      model.objects[1].rightX = markX;
      model.objects[1].rightY = markY;
      rightCoverDown = isTall && isClose;
    }
    if (leftCoverDown && rightCoverDown) {
      ans[5] = "DOWN";
    }

    // Basis extrapolate is pretty good for 21.
    for (int i = 0; i < 22; ++i) {
      blockBasisExtrapolate[i] = true;
    }
    blockBasisExtrapolate[21] = false;
    if (model.objects[2].leftFound && model.objects[3].leftFound) {
      if (model.objects[10].leftFound) {
        basisExtrapolate(leftEyeImage, sccPropLeft, model, true, 10, compRatios10);
      } else {
        basisExtrapolateWeak(leftEyeImage, sccPropLeft, model, true);
      }
    }
    if (model.objects[2].rightFound && model.objects[3].rightFound) {
      if (model.objects[10].rightFound) {
        basisExtrapolate(rightEyeImage, sccPropRight, model, false, 10, compRatios10);
      } else {
        basisExtrapolateWeak(rightEyeImage, sccPropRight, model, false);
      }
    }
    for (int i = 0; i < 22; ++i) {
      blockBasisExtrapolate[i] = false;
    }

    // Refill to maybe get 16/18.
    fillGrid(leftEyeImage, sccPropLeft, model, true);
    fillGrid(rightEyeImage, sccPropRight, model, false);

    // Special logic to find straight-line extrapolation of 16/18.
    /*if (model.objects[21].leftFound && !model.objects[18].leftFound && !model.objects[16].leftFound) {
      int objX = model.objects[21].leftX;
      int objY = model.objects[21].leftY;
      vector<SccProperty*> candidates;
      for (map<int, SccProperty>::iterator it = sccPropLeft.begin(); it != sccPropLeft.end(); ++it) {
        SccProperty& prop = it->second;
        if (prop.x < objX - width / 25 &&
            prop.size > width / 4 &&
            prop.ymax - prop.ymin > height / 200 &&
            abs(prop.y - objY) < height / 21 &&
            !prop.isGreenLed && !prop.isBlueLed &&
            canReach(width, height, prop.x, prop.y, objX, objY, sccLeft, sccPropLargeLeft)) {
          candidates.push_back(&prop);
        }
      }
      double bestScore = 999;
      SccProperty* cand16 = NULL;
      SccProperty* cand18 = NULL;
      for (int i = 0; i < candidates.size(); ++i) {
        SccProperty* prop16 = candidates[i];
        for (int j = 0; j < candidates.size(); ++j) {
          if (i == j) continue;
          SccProperty* prop18 = candidates[j];
          if (prop16->x >= prop18->x) continue;
          int dx16 = prop18->x - prop16->x;
          int dy16 = prop18->y - prop16->y;
          double diff16 = sqrt(dx16 * dx16 + dy16 * dy16);
          int dx18 = objX - prop18->x;
          int dy18 = objY - prop18->y;
          double diff18 = sqrt(dx18 * dx18 + dy18 * dy18);
          double pairDiff = abs(diff16 - diff18);
          double theta = (dx16 * dx18 + dy16 * dy18) / diff16 / diff18;
          theta = abs(1.0 - theta);
          double score = (theta + 1) * (pairDiff + 1);
          if (score < bestScore) {
            bestScore = score;
            cand16 = prop16;
            cand18 = prop18;
          }
        }
      }
      if (cand16 != NULL && cand18 != NULL) {
        cerr <<"Found 16 and 18 via linear matching: " << cand16->x << "," << cand16->y << " "
            << cand18->x << "," << cand18->y << " with score " << bestScore << endl;
        model.objects[16].leftX = cand16->x;
        model.objects[16].leftY = cand16->y;
        model.objects[16].leftFound = true;
        model.objects[16].leftProp = cand16;
        model.objects[18].leftX = cand18->x;
        model.objects[18].leftY = cand18->y;
        model.objects[18].leftFound = true;
        model.objects[18].leftProp = cand18;
      }
    }
    if (model.objects[21].rightFound && !model.objects[18].rightFound && !model.objects[16].rightFound) {
      int objX = model.objects[21].rightX;
      int objY = model.objects[21].rightY;
      vector<SccProperty*> candidates;
      for (map<int, SccProperty>::iterator it = sccPropRight.begin(); it != sccPropRight.end(); ++it) {
        SccProperty& prop = it->second;
        if (prop.x < objX - width / 25 &&
            prop.size > width / 4 &&
            prop.ymax - prop.ymin > height / 200 &&
            abs(prop.y - objY) < height / 21 &&
            !prop.isGreenLed && !prop.isBlueLed &&
            canReach(width, height, prop.x, prop.y, objX, objY, sccRight, sccPropLargeRight)) {
          candidates.push_back(&prop);
        }
      }
      double bestScore = 999;
      SccProperty* cand16 = NULL;
      SccProperty* cand18 = NULL;
      for (int i = 0; i < candidates.size(); ++i) {
        SccProperty* prop16 = candidates[i];
        for (int j = 0; j < candidates.size(); ++j) {
          if (i == j) continue;
          SccProperty* prop18 = candidates[j];
          if (prop16->x >= prop18->x) continue;
          int dx16 = prop18->x - prop16->x;
          int dy16 = prop18->y - prop16->y;
          double diff16 = sqrt(dx16 * dx16 + dy16 * dy16);
          int dx18 = objX - prop18->x;
          int dy18 = objY - prop18->y;
          double diff18 = sqrt(dx18 * dx18 + dy18 * dy18);
          double pairDiff = abs(diff16 - diff18);
          double theta = (dx16 * dx18 + dy16 * dy18) / diff16 / diff18;
          theta = abs(1.0 - theta);
          double score = (theta + 1) * (pairDiff + 1);
          if (score < bestScore) {
            bestScore = score;
            cand16 = prop16;
            cand18 = prop18;
          }
        }
      }
      if (cand16 != NULL && cand18 != NULL) {
        cerr <<"Found 16 and 18 via linear matching: " << cand16->x << "," << cand16->y << " "
            << cand18->x << "," << cand18->y << " with score " << bestScore << endl;
        model.objects[16].rightX = cand16->x;
        model.objects[16].rightY = cand16->y;
        model.objects[16].rightFound = true;
        model.objects[16].rightProp = cand16;
        model.objects[18].rightX = cand18->x;
        model.objects[18].rightY = cand18->y;
        model.objects[18].rightFound = true;
        model.objects[18].rightProp = cand18;
      }
    }*/

    if (model.objects[2].leftFound && model.objects[3].leftFound) {
      if (model.objects[10].leftFound) {
        basisDfsMatch(leftEyeImage, sccLeft, sccPropLeft, sccPropLargeLeft, model, true, 10, compRatios10);
      } /*else if (model.objects[9].leftFound) {
        basisDfsMatch(leftEyeImage, sccPropLeft, model, true, 9, compRatios9);
      } */else {
        basisDfsMatchWeak(leftEyeImage, sccLeft, sccPropLeft, sccPropLargeLeft, model, true);
      }
    }
    if (model.objects[2].rightFound && model.objects[3].rightFound) {
      if (model.objects[10].rightFound) {
        basisDfsMatch(rightEyeImage, sccRight, sccPropRight, sccPropLargeRight, model, false, 10, compRatios10);
      } /*else if (model.objects[9].rightFound) {
        basisDfsMatch(rightEyeImage, sccPropRight, model, false, 9, compRatios9);
      } */else {
        basisDfsMatchWeak(rightEyeImage, sccRight, sccPropRight, sccPropLargeRight, model, false);
      }
    }

    if (model.objects[2].leftFound && model.objects[3].leftFound) {
      if (model.objects[10].leftFound) {
        basisExtrapolate(leftEyeImage, sccPropLeft, model, true, 10, compRatios10);
      } /*else if (model.objects[9].leftFound) {
        basisExtrapolate(leftEyeImage, sccPropLeft, model, true, 9, compRatios9);
      } */else {
        basisExtrapolateWeak(leftEyeImage, sccPropLeft, model, true);
      }
    }
    if (model.objects[2].rightFound && model.objects[3].rightFound) {
      if (model.objects[10].rightFound) {
        basisExtrapolate(rightEyeImage, sccPropRight, model, false, 10, compRatios10);
      } /*else if (model.objects[9].rightFound) {
        basisExtrapolate(rightEyeImage, sccPropRight, model, false, 9, compRatios9);
      } */else {
        basisExtrapolateWeak(rightEyeImage, sccPropRight, model, false);
      }
    }
#ifdef DEBUG
  cerr << "Took " << CurrentTimeMillis() - startTime << " ms to look for dfs match all" << endl;
  startTime = CurrentTimeMillis();
#endif
    reconcileMissing(model);
  }

  // Flip inverted orderings.
  /*for (int i = 6; i <= 8; ++i) {
    int j = i + 6;
    if (model.objects[i].leftY > model.objects[j].leftY) {
      int tmpX = model.objects[i].leftX;
      int tmpY = model.objects[i].leftY;
      model.objects[i].leftX = model.objects[j].leftX;
      model.objects[i].leftY = model.objects[j].leftY;
      model.objects[j].leftX = tmpX;
      model.objects[j].leftY = tmpY;
    }
    if (model.objects[i].rightY > model.objects[j].rightY) {
      int tmpX = model.objects[i].rightX;
      int tmpY = model.objects[i].rightY;
      model.objects[i].rightX = model.objects[j].rightX;
      model.objects[i].rightY = model.objects[j].rightY;
      model.objects[j].rightX = tmpX;
      model.objects[j].rightY = tmpY;
    }
  }*/

  /*for (int i = 0; i < 22; ++i) {
    model.objects[i].leftX *= 2;
    model.objects[i].leftY *= 2;
    model.objects[i].rightX *= 2;
    model.objects[i].rightY *= 2;
  }*/
  for (int i = 0; i < 22; ++i) {
    model.objects[i].leftX = max(model.objects[i].leftX, 0);
    model.objects[i].leftX = min(model.objects[i].leftX, width - 1);
    model.objects[i].leftY = max(model.objects[i].leftY, 0);
    model.objects[i].leftY = min(model.objects[i].leftY, height - 1);
    model.objects[i].rightX = max(model.objects[i].rightX, 0);
    model.objects[i].rightX = min(model.objects[i].rightX, width - 1);
    model.objects[i].rightY = max(model.objects[i].rightY, 0);
    model.objects[i].rightY = min(model.objects[i].rightY, height - 1);
  }
  {
    char buf[50];
    for (int obj = 0; obj <= 3; ++obj) {
      int base = obj * 5;
      sprintf(buf, "%d", model.objects[obj].leftX);
      ans[base + 1] = buf;
      sprintf(buf, "%d", model.objects[obj].leftY);
      ans[base + 2] = buf;
      sprintf(buf, "%d", model.objects[obj].rightX);
      ans[base + 3] = buf;
      sprintf(buf, "%d", model.objects[obj].rightY);
      ans[base + 4] = buf;
    }
    {
      // Green leds.
      int possibleComps[9] = { 4, 5, 7, 9, 11, 13, 16, 18, 21 };
      for (int i = 0; i < 9; ++i) {
        int obj = possibleComps[i];
        int base = obj * 5;
        if (model.objects[obj].leftFound) {
          if (i > 1) ans[base] = model.objects[obj].state;
          sprintf(buf, "%d", model.objects[obj].leftX);
          ans[base + 1] = buf;
          sprintf(buf, "%d", model.objects[obj].leftY);
          ans[base + 2] = buf;
        }
        if (model.objects[obj].rightFound) {
          if (i > 1) ans[base] = model.objects[obj].state;
          sprintf(buf, "%d", model.objects[obj].rightX);
          ans[base + 3] = buf;
          sprintf(buf, "%d", model.objects[obj].rightY);
          ans[base + 4] = buf;
        }
      }
    }
    {
      // blue leds.
      int possibleComps[6] = { 6, 8, 10, 12, 14, 19 };
      for (int i = 0; i < 6; ++i) {
        int obj = possibleComps[i];
        int base = obj * 5;
        if (model.objects[obj].leftFound) {
          ans[base] = model.objects[obj].state;
          sprintf(buf, "%d", model.objects[obj].leftX);
          ans[base + 1] = buf;
          sprintf(buf, "%d", model.objects[obj].leftY);
          ans[base + 2] = buf;
        }
        if (model.objects[obj].rightFound) {
          ans[base] = model.objects[obj].state;
          sprintf(buf, "%d", model.objects[obj].rightX);
          ans[base + 3] = buf;
          sprintf(buf, "%d", model.objects[obj].rightY);
          ans[base + 4] = buf;
        }
      }
    }
    {
      // All the rest.
      int possibleComps[3] = { 15, 17, 20 };
      for (int i = 0; i < 3; ++i) {
        int obj = possibleComps[i];
        int base = obj * 5;
        if (model.objects[obj].leftFound) {
          ans[base] = model.objects[obj].state;
          sprintf(buf, "%d", model.objects[obj].leftX);
          ans[base + 1] = buf;
          sprintf(buf, "%d", model.objects[obj].leftY);
          ans[base + 2] = buf;
        }
        if (model.objects[obj].rightFound) {
          ans[base] = model.objects[obj].state;
          sprintf(buf, "%d", model.objects[obj].rightX);
          ans[base + 3] = buf;
          sprintf(buf, "%d", model.objects[obj].rightY);
          ans[base + 4] = buf;
        }
      }
    }

    // Let the three bottom LEDs propagate state for the toggle switches.
    if (model.objects[16].state == "ON") {
      ans[75] = "UP";
    }
    if (model.objects[18].state == "ON") {
      ans[85] = "UP";
    } else if (model.objects[19].state == "ON") {
      ans[85] = "DOWN";
    }
    if (model.objects[21].state == "ON") {
      ans[100] = "UP";
    }
  }

  vector<string> ret(22);
  for (int i = 0, cur = 0; i < 22; ++i) {
    string buf;
    for (int j = 0; j < 5; ++j, ++cur) {
      buf += ans[cur];
      if (j != 4) {
        buf += ",";
      }
    }
#ifdef DEBUG
    cerr << buf << endl;
#endif
    ret[i] = buf;
  }
  return ret;
}

#include <sys/time.h>
int64 CurrentTimeMillis() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (int64) (tv.tv_sec * 1000 + tv.tv_usec / 1000);
}

int ntoh(int orig) {
  return ((orig >> 24) & 0x000000ff) |
         ((orig >> 8) & 0x0000ff00) |
         ((orig << 8) & 0x00ff0000) |
         ((orig << 24) & 0xff000000);
}

int main() {
  int64 start_time = CurrentTimeMillis();
  int length;
  cin.read((char*)(&length), 4);
  length = ntoh(length);
  vector<int> left(length);
  vector<int> right(length);
  cin.read((char*)(&left[0]), 4 * length);
  cin.read((char*)(&right[0]), 4 * length);
  int64 done_read = CurrentTimeMillis();

  RobonautEye sol;
  vector<string> ret = sol.recognizeObjects(left, right);
  int64 done_compute = CurrentTimeMillis();

  for (int i = 0; i < 22; ++i) {
    cout << ret[i] << endl;
  }
  cout.flush();

  cerr << "Read time: " << done_read - start_time << endl;
  cerr << "Compute time: " << done_compute - done_read << endl;
  return 0;
}
