#ifndef UTILITY_H
#define UTILITY_H

#include <vector>
#include <QImage>
#include <QFile>
#include <QStringList>
#include <map>
#include <QPainter>
#include <QColor>
#include <QDebug>
#include <QDir>
#include <QTime>
#include <QDate>
#include <QElapsedTimer>
#include "robonauteye.h"
#include "imageutility.h"
#include "image.h"
#include "geometry.h"

using namespace std;

struct SecondaryObjectPosition {
  int x;
  int y;
  int imageType; // 0: screws, 1: numbers and letters
};

struct RealObjectInformations {
  QString state;
  int leftX;
  int leftY;
  int rightX;
  int rightY;
};

struct DetailedScore {
  float hiddenScore;
  float proximityScore;
  float stateScore;
  float overallScore;
  int nbFiles;
};

class Utility
{
public:
  Utility();

  static vector<int> loadImage(QString path);
  static map <QString, vector <RealObjectInformations> > loadCSV(QString path);
  static void drawInformations(QPainter & p, vector <RealObjectInformations> infos, bool isLeft = true, float scale = 1);
  static void drawSecondaryObjects(QPainter & p, vector <SecondaryObjectPosition> secondaryPositions, float scale = 1);
  static void drawSquares(QPainter & p, vector <ObjectInformations> squares, float scale = 1);
  static map <QString, vector <RealObjectInformations> > loadSecondaryObjects(QString csvPath, QString secondaryPath);
  static void saveSecondaryObjectsFile(QString secondaryObjectsPath, vector <SecondaryObjectPosition> secondaryPositions);
  static vector <SecondaryObjectPosition> loadSecondaryObjectsFile(QString secondaryObjectsPath);
  static void scaleInformations(vector <RealObjectInformations> & infos, float scale = 1);
  static vector <ObjectInformations> loadObjectsInformations(QString path);
  static void saveObjectsInformations(QString path, vector<ObjectInformations> & infos);
  static void exportObjectsInformationsToImages(QString pathObjectsInformations, QString pathTo);
  static void drawPositions(QPainter & p, vector<Point> & positions, float scale);
  static void drawPositions(QPainter & p, Point * positions, size_t nb, float scale);
  static void drawPositionsWithIndex(QPainter & p, vector<Point> & positions, float scale);
  static void drawPolygon(QPainter & p, Point * positions, size_t nb, float scale);
  static void drawPrimaryObject(QPainter & p, PrimaryObject & primaryObject, float scale);
  static void exportPrimaryObjectsToImages(QString pathTo);
  static vector<Point> assignSecondaryInformations(vector <RealObjectInformations> & imageSecondaryInfos, vector<Point> & projections, bool isLeft);
  static vector<Point> getImageSecondaryPoint(vector <RealObjectInformations> & imageSecondaryInfos, bool isLeft);
  static vector< vector<float> > * getOrthographicHeuristicStatisticsItem(vector<Point> & assignments);
  static void generateOrthographicHeuristicStatistics();
  static vector <Point> getPrimaryObjectsPositions(vector <RealObjectInformations> & infos, bool isLeft);
  static void upgradeSecondaryPositions();
  static void selectSecondaryObjects(vector <RealObjectInformations> & primaryObjects);
  static void getScoreStatistics();
  static DetailedScore evaluateScoreFile(QString file);
  static DetailedScore evaluateScoreItem(RealObjectInformations * primaryInfos, FinalObject * finalObjects, float threshold = 15);
  static void exportDecals(QString pathTo);
  static void getMostFrequentState(QString modelCsvPath);
};

#endif // UTILITY_H
