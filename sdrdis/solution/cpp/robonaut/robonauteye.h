#ifndef ROBONAUTEYE_H
#define ROBONAUTEYE_H


#include <vector>
#include <string>
#include "image.h"
#include "imagesmoother.h"
#include "statistics.h"
#include "imagedisjointset.h"
#include "felzenhuttensegmentation.h"
#include "standardselector.h"
#include "standardmerger.h"
#include "secondaryobjectpca.h"
#include "secondaryobjectneuralnetwork.h"
#include "primaryobjectpca.h"
#include "primaryobjectneuralnetwork.h"
#include "geometry.h"
#include "quadrupleindexes.h"
#include <sstream>

using namespace std;

#define PREPROCESS_BLUR 1.1
#define PREPROCESS_CONTRAST 100

#define SEGMENT_THRESHOLD 200
#define SEGMENT_MIN_SIZE 30

#define SELECT_1_MAX_SIZE 600
#define SELECT_1_MIN_WIDTH_HEIGHT_RATIO 0.25
#define SELECT_1_MAX_WIDTH_HEIGHT_RATIO 4
#define SELECT_1_MIN_SIZE_BOX_RATIO 0.15

#define MERGE_MAX_COLOR_DIFF 50
#define MERGE_MIN_SHARED_THRESHOLD 0.2

#define SELECT_2_MAX_SIZE 250
#define SELECT_2_MIN_WIDTH_HEIGHT_RATIO 0.5
#define SELECT_2_MAX_WIDTH_HEIGHT_RATIO 2
#define SELECT_2_MIN_SIZE_BOX_RATIO 0.5

#define AREA_MARGIN 0.3
#define AREA_MIN_SIZE 30
#define AREA_LOCATION_BUCKET_SIZE 10

#define SECONDARY_OBJECT_CLASSIFIER_IMAGE_SIZE 20
#define PRIMARY_OBJECT_CLASSIFIER_IMAGE_SIZE 30

#define SECONDARY_OBJECT_CLASSIFIER_SCORE_THRESHOLD 0.35
#define SECONDARY_OBJECT_HOMOGRAPHY_MIN_SIZE_QUADRUPLE 250.0
#define SECONDARY_OBJECT_HOMOGRAPHY_SUCCESS_THRESHOLD 0.6
#define SECONDARY_OBJECT_HOMOGRAPHY_SUCCESS_LOWERBOUND 550.0
#define SECONDARY_OBJECT_HOMOGRAPHY_SUCCESS_UPPERBOUND 850.0
#define SECONDARY_OBJECT_HOMOGRAPHY_PURIFY_MIN_FACTOR 1.2
#define SECONDARY_OBJECT_HOMOGRAPHY_PURIFY_MAX_POINTS 2

#define SECONDARY_OBJECT_SEARCH_AREA 30

struct ObjectInformations {
  size_t index;
  string state;
  int x;
  int y;
  int width;
  int height;
  int centerX;
  int centerY;
  int sizeWeight;
  float secondaryObjectClassifierScore;
  float imageTypeScores[4];

  ObjectInformations() {
    imageTypeScores[0] = 0;
    imageTypeScores[1] = 0;
    imageTypeScores[2] = 0;
    imageTypeScores[3] = 0;
  }

  ObjectInformations(float centerX, float centerY) {
    this->centerX = centerX;
    this->centerY = centerY;
    ObjectInformations();
  }
};


struct ObjectInformationsMapping {
  int * objectInformationsMap;
  vector<ObjectInformations> positives;
};

struct ProjectionVector {
  float angle;
  float distance;
};

struct SecondaryObject {
  Point position;
  int imageType;

  SecondaryObject(float x, float y, int imageType) {
    this->position.x = x;
    this->position.y = y;
    this->imageType = imageType;
  }
};

struct PrimaryObject {
  Point position;
  int size;
  float stateScores[3];
  bool isCorrect;

  PrimaryObject() {
    stateScores[0] = 0;
    stateScores[1] = 0;
    stateScores[2] = 0;
    isCorrect = false;
  }

  PrimaryObject(float x, float y, int size) {
    stateScores[0] = 0;
    stateScores[1] = 0;
    stateScores[2] = 0;
    isCorrect = false;
    this->position.x = x;
    this->position.y = y;
    this->size = size;
  }
};

struct FinalObject {
  int state;
  Point leftPosition;
  Point rightPosition;
};

static const SecondaryObject TASKBOARD_SECONDARY_OBJECTS_POSITIONS[] = {
  SecondaryObject(54, 112, 1),  // 0
  SecondaryObject(129, 52, 1),  // 1
  SecondaryObject(182, 70, 1),  // 2
  SecondaryObject(314, 70, 1),  // 3
  SecondaryObject(429, 52, 1),  // 4
  SecondaryObject(504, 112, 1), // 5
  SecondaryObject(415, 133, 2), // 6
  SecondaryObject(71, 239, 1),  // 7
  SecondaryObject(200, 239, 1), // 8
  SecondaryObject(314, 239, 1), // 9
  SecondaryObject(415, 313, 2), // 10
  SecondaryObject(415, 432, 2), // 11
  SecondaryObject(129, 545, 2), // 12
  SecondaryObject(279, 545, 2), // 13
  SecondaryObject(429, 545, 2), // 14
  SecondaryObject(279, 664, 2), // 15
  SecondaryObject(54, 677, 1),  // 16
  SecondaryObject(129, 736, 1), // 17
  SecondaryObject(429, 736, 1), // 18
  SecondaryObject(504, 677, 1), // 19
  SecondaryObject(118, 261, 3), // 20
  SecondaryObject(180, 261, 3), // 21
  SecondaryObject(244, 261, 3), // 22
  SecondaryObject(68, 310, 3),  // 23
  SecondaryObject(68, 372, 3),  // 24
  SecondaryObject(68, 435, 3),  // 25
};

#define TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS 26

static const PrimaryObject TASKBOARD_PRIMARY_OBJECTS_POSITIONS[] = {
  PrimaryObject(369, 133, 180), // 0
  PrimaryObject(415, 133, 40), // 1
  PrimaryObject(415, 373, 140), // 2
  PrimaryObject(415, 313, 40), // 3
  PrimaryObject(415, 432, 40), // 4
  PrimaryObject(116, 309, 50), // 5
  PrimaryObject(180, 309, 50), // 6
  PrimaryObject(244, 309, 50), // 7
  PrimaryObject(116, 372, 50), // 8
  PrimaryObject(180, 372, 50), // 9
  PrimaryObject(244, 372, 50), // 10
  PrimaryObject(116, 436, 50), // 11
  PrimaryObject(180, 436, 50), // 12
  PrimaryObject(244, 436, 50), // 13
  PrimaryObject(128, 602, 140), // 14
  PrimaryObject(128, 545, 40), // 15
  PrimaryObject(279, 602, 140), // 16
  PrimaryObject(279, 545, 40), // 17
  PrimaryObject(279, 664, 40), // 18
  PrimaryObject(429, 602, 140), // 19
  PrimaryObject(429, 545, 40), // 20
};

static const size_t TASKBOARD_PRIMARY_OBJECTS_DEFAULT_STATES[] = {
  1, // 0
  1, // 1
  2, // 2
  0, // 3
  0, // 4
  0, // 5
  0, // 6
  0, // 7
  0, // 8
  0, // 9
  0, // 10
  0, // 11
  0, // 12
  0, // 13
  0, // 14
  0, // 15
  2, // 16
  0, // 17
  0, // 18
  0, // 19
  0, // 20
};

static const size_t TASKBOARD_PRIMARY_OBJECTS_NB_POSITIONS = 21;

static const size_t TASKBOARD_NB_FINAL_OBJECTS = 22;

static const bool IS_FINAL_OBJECT_ON_OFF[] = {
  false, // 0
  false, // 1
  true,  // 2
  false, // 3
  true,  // 4
  true,  // 5
  true,  // 6
  true,  // 7
  true,  // 8
  true,  // 9
  true,  // 10
  true,  // 11
  true,  // 12
  true,  // 13
  true,  // 14
  false, // 15
  true,  // 16
  false, // 17
  true,  // 18
  true,  // 19
  false, // 20
  true   // 21
};

struct SecondaryObjectAssociations {
  int associations[TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS];
  float score;
};

class RobonautEye
{
public:
  RobonautEye();
  vector<string> recognizeObjects(vector<int> leftEyeImage, vector<int> rightEyeImage);
  FinalObject * getFinalObjects(Image<float> * leftImage, Image<float> * rightImage);
  PrimaryObject * processImage(Image<float> * image);
  PrimaryObject * processImage(Image<float> * image, Array2D<double> * &H);
  static Image<float> * convertVectorToImage(vector<int> & imageVector);
  ObjectInformationsMapping * getSecondaryObjectsInformationsFromImage(Image<float> * realImage);
  static Image<float> * preprocessImage(Image<float> * from, float blur, float contrast);
  static vector < Set<float> * > * segmentImage(Image<float> * image, float splitThreshold, float splitMinSize);
  void computeSecondaryObjectClassifierScore(Image<float> * image, ObjectInformations & objectInformations);
  static PrimaryObject getPerspectiveProjection(Array2D<double> * H, const PrimaryObject & primaryObject);
  static Image<float> * getPrimaryObjectImage(Image<float> * fromImage, PrimaryObject & primaryObject);
  void computePrimaryObjectState(Image<float> * image, PrimaryObject & primaryObject);
  static ProjectionVector getOrthographicProjectionVector(OrthographicProjectionEquation & eq);
  static float * getOrthographicProjectionInformations(Point & fromA, Point & fromB, Point & fromC, Point & toA, Point & toB, Point & toC);
  static float * getOrthographicProjectionInformations(Point & fromA, Point & fromB, Point & fromC, Point & fromD, Point & toA, Point & toB, Point & toC, Point & toD);
  static bool isSetCorrect(Point & fromA, Point & fromB, Point & fromC, Point & toA, Point & toB, Point & toC);
  static bool isSetCorrect(Point & fromA, Point & fromB, Point & fromC, Point & fromD, Point & toA, Point & toB, Point & toC, Point & toD);
  Array2D<double> * getHomography(Image<float> * realImage, ObjectInformationsMapping * infos, bool & error);
  static multiset<QuadrupleIndexes> * getQuadruples(vector <ObjectInformations> & points);
  SecondaryObjectAssociations findClosestQuadruples(Image<float> * realImage, const QuadrupleIndexes & quadrupleFrom, ObjectInformationsMapping * infos);
  SecondaryObjectAssociations evaluateHomography(Image<float> * realImage, Array2D<double> * H, ObjectInformationsMapping * infos);
  int searchForSecondaryObjectOnArea(Image<float> * realImage, Point & point, int imageType, vector <bool> * blacklist, ObjectInformationsMapping * infos);
  static void getPrimaryItemInformations(size_t id, Array2D<double> * H, float &length1, float &length2, float &angle1, float &angle2, Point &toPosition);
  static Point getSwitchPosition(Array2D<double> * H);
  static Point getSwitchDecal(float length1, float length2, float angle1, float angle2);
  static Point getPowerSwitchPositionDown(Array2D<double> * H);
  static Point getPowerSwitchDecalDown(float length1, float length2, float angle1, float angle2);
  static Point getPowerSwitchPositionUp(Array2D<double> * H);
  static Point getPowerSwitchDecalUp(float length1, float length2, float angle1, float angle2);
  Point * getPrimaryPoints();
  Point * getSecondaryPoints();
  Array2D<double> * getPurifiedHomography(Image<float> * realImage, ObjectInformationsMapping * infos, SecondaryObjectAssociations association, float scoreThreshold, size_t maxDepth, float &score, bool &satisfied);
  Array2D<double> * getPurifiedHomography(Image<float> * realImage, ObjectInformationsMapping * infos, vector<Point> &pointsFrom, vector<Point> &pointsTo, float scoreThreshold, size_t maxDepth, float &score, bool &satisfied);
  void addObjectToMapping(ObjectInformations& obj, ObjectInformationsMapping * mapping, Image<float> * image);
  void refreshMapping(int bucketX, int bucketY, ObjectInformationsMapping * mapping, Image<float> * image);
  //void searchForUndetectedSecondaryObjects(Image<float> * realImage, Point & point, vector<bool> * blacklist, ObjectInformationsMapping * mapping);
  //void locateUndetectedSecondaryObject(Image<float> * realImage, int bucketX, int bucketY, int meanSize, bool alreadyDetected, ObjectInformationsMapping * mapping, vector<Point> & points);
  void correctObject(ObjectInformations & obj, Image<float> * image);

protected:

  StandardSelectorSettings firstSelectSettings;
  StandardMergerSettings mergeSettings;
  StandardSelectorSettings secondSelectSettings;

  SecondaryObjectPCA secondaryObjectPCA;
  PrimaryObjectPCA primaryObjectPCA;
  SecondaryObjectNeuralNetwork secondaryObjectNeuralNetwork;
  PrimaryObjectNeuralNetwork primaryObjectNeuralNetwork;
};

#endif // ROBONAUTEYE_H
