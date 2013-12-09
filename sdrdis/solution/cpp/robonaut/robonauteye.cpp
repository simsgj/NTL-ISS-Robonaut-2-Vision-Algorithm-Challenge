#include "robonauteye.h"


#define bestIndex(scores) (scores[0] > scores[1] ? (scores[0] > scores[2] ? 0 : 2) : (scores[1] > scores[2] ? 1 : 2))
#define lowerBoundScore(infos) (max(SECONDARY_OBJECT_HOMOGRAPHY_SUCCESS_LOWERBOUND, min(SECONDARY_OBJECT_HOMOGRAPHY_SUCCESS_UPPERBOUND, SECONDARY_OBJECT_HOMOGRAPHY_SUCCESS_THRESHOLD * infos->positives.size() * 100.0)))


static Point * TASKBOARD_PRIMARY_OBJECTS_POSITIONS_POINTS;
static Point * TASKBOARD_SECONDARY_OBJECTS_POSITIONS_POINTS;

RobonautEye::RobonautEye()
{
  firstSelectSettings.minSize = 0;
  firstSelectSettings.maxSize = SELECT_1_MAX_SIZE;
  firstSelectSettings.minWidthHeightRatio = SELECT_1_MIN_WIDTH_HEIGHT_RATIO;
  firstSelectSettings.maxWidthHeightRatio = SELECT_1_MAX_WIDTH_HEIGHT_RATIO;
  firstSelectSettings.minSizeBoxRatio = SELECT_1_MIN_SIZE_BOX_RATIO;

  mergeSettings.maxColorDiff = MERGE_MAX_COLOR_DIFF;
  mergeSettings.minSharedThreshold = MERGE_MIN_SHARED_THRESHOLD;

  secondSelectSettings.minSize = 0;
  secondSelectSettings.maxSize = SELECT_2_MAX_SIZE;
  secondSelectSettings.minWidthHeightRatio = SELECT_2_MIN_WIDTH_HEIGHT_RATIO;
  secondSelectSettings.maxWidthHeightRatio = SELECT_2_MAX_WIDTH_HEIGHT_RATIO;
  secondSelectSettings.minSizeBoxRatio = SELECT_2_MIN_SIZE_BOX_RATIO;


  TASKBOARD_PRIMARY_OBJECTS_POSITIONS_POINTS = new Point[TASKBOARD_PRIMARY_OBJECTS_NB_POSITIONS];
  for (size_t i = 0; i < TASKBOARD_PRIMARY_OBJECTS_NB_POSITIONS; i++) {
    TASKBOARD_PRIMARY_OBJECTS_POSITIONS_POINTS[i] = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[i].position;
  }

  TASKBOARD_SECONDARY_OBJECTS_POSITIONS_POINTS = new Point[TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS];
  for (size_t i = 0; i < TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS; i++) {
    TASKBOARD_SECONDARY_OBJECTS_POSITIONS_POINTS[i] = TASKBOARD_SECONDARY_OBJECTS_POSITIONS[i].position;
  }
}

Point * RobonautEye::getPrimaryPoints() {
  return TASKBOARD_PRIMARY_OBJECTS_POSITIONS_POINTS;
}

Point * RobonautEye::getSecondaryPoints() {
  return TASKBOARD_SECONDARY_OBJECTS_POSITIONS_POINTS;
}

vector<string> RobonautEye::recognizeObjects(vector<int> leftEyeImage, vector<int> rightEyeImage) {
  vector<string> results;

  Image <float> * leftImage = convertVectorToImage(leftEyeImage);
  Image <float> * rightImage = convertVectorToImage(rightEyeImage);

  FinalObject * finalObjects = getFinalObjects(leftImage, rightImage);

  ostringstream convert;
  string resultItem;
  for (size_t i = 0; i < TASKBOARD_NB_FINAL_OBJECTS; i++) {
    resultItem = "";
    if (finalObjects[i].state == 0) {
      resultItem = IS_FINAL_OBJECT_ON_OFF[i] ? "OFF," : "DOWN,";
    } else if (finalObjects[i].state == 1) {
      resultItem = IS_FINAL_OBJECT_ON_OFF[i] ? "ON," : "UP,";
    } else {
      resultItem = "CENTER,";
    }
    convert.str("");
    convert.clear();
    convert << round(finalObjects[i].leftPosition.x);
    resultItem += convert.str() + ",";
    convert.str("");
    convert.clear();
    convert << round(finalObjects[i].leftPosition.y);
    resultItem += convert.str() + ",";
    convert.str("");
    convert.clear();
    convert << round(finalObjects[i].rightPosition.x);
    resultItem += convert.str() + ",";
    convert.str("");
    convert.clear();
    convert << round(finalObjects[i].rightPosition.y);
    resultItem += convert.str();


    results.push_back(resultItem);
  }

  delete leftImage;
  delete rightImage;


  return results;
}

FinalObject * RobonautEye::getFinalObjects(Image<float> * leftImage, Image<float> * rightImage) {
  // initialization here to prevent warning, but doesn't that create a memory leak ?
  Array2D<double> * leftHomography = new Array2D<double>();
  Array2D<double> * rightHomography = new Array2D<double>();
  PrimaryObject * leftPrimaryObjects = processImage(leftImage, leftHomography);
  PrimaryObject * rightPrimaryObjects = processImage(rightImage, rightHomography);
  FinalObject * finalObjects = new FinalObject[TASKBOARD_NB_FINAL_OBJECTS];
  float localScores[3];
  float stateScores[TASKBOARD_PRIMARY_OBJECTS_NB_POSITIONS][3];
  for (size_t i = 0; i < TASKBOARD_PRIMARY_OBJECTS_NB_POSITIONS; i++) {
    for (size_t j = 0; j < 3; j++) {
      stateScores[i][j] = leftPrimaryObjects[i].stateScores[j] + rightPrimaryObjects[i].stateScores[j];
    }
    if (!leftPrimaryObjects[i].isCorrect && !rightPrimaryObjects[i].isCorrect) {
      stateScores[i][TASKBOARD_PRIMARY_OBJECTS_DEFAULT_STATES[i]] = 1;
    }
  }

  finalObjects[0].leftPosition.x = leftPrimaryObjects[0].position.x;
  finalObjects[0].leftPosition.y = leftPrimaryObjects[0].position.y;
  finalObjects[0].rightPosition.x = rightPrimaryObjects[0].position.x;
  finalObjects[0].rightPosition.y = rightPrimaryObjects[0].position.y;
  finalObjects[1].leftPosition.x = -1;
  finalObjects[1].leftPosition.y = -1;
  finalObjects[1].rightPosition.x = -1;
  finalObjects[1].rightPosition.y = -1;
  for (size_t i = 2; i < TASKBOARD_NB_FINAL_OBJECTS; i++) {
    finalObjects[i].leftPosition.x = leftPrimaryObjects[i - 1].position.x;
    finalObjects[i].leftPosition.y = leftPrimaryObjects[i - 1].position.y;
    finalObjects[i].rightPosition.x = rightPrimaryObjects[i - 1].position.x;
    finalObjects[i].rightPosition.y = rightPrimaryObjects[i - 1].position.y;
  }

  // PANEL_POWER_SWITCH + PANEL_POWER_LEDs
  finalObjects[0].state = stateScores[1][0] > stateScores[1][1] ? 0 : 1;

  // NUM_PAD
  for (size_t i = 6; i <= 14; i++) {
    finalObjects[i].state = (stateScores[i - 1][0] - stateScores[i -  1][1]) > 0 ? 0 : 1;
    if (finalObjects[i].state == 1) {
      finalObjects[0].state = 1;
    }
  }

  // PANEL_POWER_SWITCH + PANEL_POWER_LEDs
  if (stateScores[3][1] > stateScores[3][0]) finalObjects[0].state = 1;
  if (stateScores[15][1] > stateScores[15][0]) finalObjects[0].state = 1;
  if (stateScores[17][1] > stateScores[17][0]) finalObjects[0].state = 1;
  if (stateScores[20][1] > stateScores[20][0]) finalObjects[0].state = 1;


  // PANEL_POWER_COVER
  finalObjects[1].state = stateScores[0][0] > stateScores[0][1] ? 0 : 1;

  // POST PROCESSING POSITION
  if (finalObjects[1].state == 0) {
    if (leftHomography->dim1() > 0) {
      finalObjects[1].leftPosition = getPowerSwitchPositionDown(leftHomography);
    }
    if (rightHomography->dim1() > 0) {
      finalObjects[1].rightPosition = getPowerSwitchPositionDown(rightHomography);
    }
  } else {
    if (leftHomography->dim1() > 0) {
      finalObjects[1].leftPosition = getPowerSwitchPositionUp(leftHomography);
    }
    if (rightHomography->dim1() > 0) {
      finalObjects[1].rightPosition = getPowerSwitchPositionUp(rightHomography);
    }
  }

  // A01_ROCKER_SWITCH + A01_ROCKER_LED_TOP + A01_ROCKER_BOTTOM
  if (finalObjects[0].state == 1) {
    localScores[0] = 0.5 * stateScores[2][0] + 1 * stateScores[3][0] + 0.2 * stateScores[4][1];
    localScores[1] = 0.5 * stateScores[2][1] + 1 * stateScores[3][1] + 0.2 * stateScores[4][0];
    localScores[2] = 0.5 * stateScores[2][2] + 1 * stateScores[3][0] + 0.2 * stateScores[4][0];
    finalObjects[3].state = bestIndex(localScores);
    finalObjects[4].state = finalObjects[3].state == 1 ? 1 : 0;
    finalObjects[5].state = finalObjects[3].state == 0 ? 1 : 0;
  } else {
    finalObjects[3].state = bestIndex(stateScores[2]);
    finalObjects[4].state = 0;
    finalObjects[5].state = 0;
  }

  // A03_TOGGLE + A03_LEFT
  if (finalObjects[0].state == 1) {
    localScores[0] = 0.5 * stateScores[14][0] + 1 * stateScores[15][0];
    localScores[1] = 0.5 * stateScores[14][1] + 1 * stateScores[15][1];
    finalObjects[15].state = localScores[0] > localScores[1] ? 0 : 1;
    finalObjects[16].state = finalObjects[15].state;
  } else {
    finalObjects[15].state = stateScores[14][0] > stateScores[14][1] ? 0 : 1;
    finalObjects[16].state = 0;
  }

  // A04_TOGGLE + A04_LED
  if (finalObjects[0].state == 1) {
    localScores[0] = 0.5 * stateScores[16][0] + 1 * stateScores[17][0] + 0.2 * stateScores[18][1];
    localScores[1] = 0.5 * stateScores[16][1] + 1 * stateScores[17][1] + 0.2 * stateScores[18][0];
    localScores[2] = 0.5 * stateScores[16][2] + 1 * stateScores[17][0] + 0.2 * stateScores[18][0];
    finalObjects[17].state = bestIndex(localScores);
    finalObjects[18].state = finalObjects[17].state == 1 ? 1 : 0;
    finalObjects[19].state = finalObjects[17].state == 0 ? 1 : 0;
  } else {
    finalObjects[17].state = bestIndex(stateScores[16]);
    finalObjects[18].state = 0;
    finalObjects[19].state = 0;
  }


  // A05_TOGGLE + A05_LED
  if (finalObjects[0].state == 1) {
    localScores[0] = 0.5 * stateScores[19][0] + 1 * stateScores[20][0];
    localScores[1] = 0.5 * stateScores[19][1] + 1 * stateScores[20][1];
    finalObjects[20].state = localScores[0] > localScores[1] ? 0 : 1;
    finalObjects[21].state = finalObjects[20].state;
  } else {
    finalObjects[20].state = stateScores[19][0] > stateScores[19][1] ? 0 : 1;
    finalObjects[21].state = 0;
  }

  finalObjects[2].state = finalObjects[0].state;

  for (size_t i = 0; i < TASKBOARD_NB_FINAL_OBJECTS; i++) {
    if (i != 1 && i != 3) {
      if (finalObjects[i].leftPosition.y != -1) {
        finalObjects[i].leftPosition.y += 5;
      }
      if (finalObjects[i].rightPosition.y != -1) {
        finalObjects[i].rightPosition.y += 5;
      }
    }

    if (finalObjects[i].leftPosition.x < 0) {
      finalObjects[i].leftPosition.x = -1;
    }
    if (finalObjects[i].rightPosition.x < 0) {
      finalObjects[i].rightPosition.x = -1;
    }
    if (finalObjects[i].leftPosition.x > leftImage->getWidth()) {
      finalObjects[i].leftPosition.x = leftImage->getWidth() - 1;
    }
    if (finalObjects[i].rightPosition.x > rightImage->getWidth()) {
      finalObjects[i].rightPosition.x = rightImage->getWidth() - 1;
    }
    if (finalObjects[i].leftPosition.y < 0) {
      finalObjects[i].leftPosition.y = -1;
    }
    if (finalObjects[i].rightPosition.y < 0) {
      finalObjects[i].rightPosition.y = -1;
    }
    if (finalObjects[i].leftPosition.y > leftImage->getHeight()) {
      finalObjects[i].leftPosition.y = leftImage->getHeight() - 1;
    }
    if (finalObjects[i].rightPosition.y > rightImage->getHeight()) {
      finalObjects[i].rightPosition.y = rightImage->getHeight() - 1;
    }
  }

  delete leftPrimaryObjects;
  delete rightPrimaryObjects;
  delete leftHomography;
  delete rightHomography;
  return finalObjects;
}

PrimaryObject * RobonautEye::processImage(Image<float> * image) {
  Array2D<double> * H = new Array2D<double>();
  PrimaryObject * res = processImage(image, H);
  delete H;
  return res;
}

PrimaryObject * RobonautEye::processImage(Image<float> * image, Array2D<double> * &H) {
  PrimaryObject * primaryObjects = new PrimaryObject[TASKBOARD_PRIMARY_OBJECTS_NB_POSITIONS];
  ObjectInformationsMapping * secondaryObjectsInformations = getSecondaryObjectsInformationsFromImage(image);
  bool error;
  if (secondaryObjectsInformations->positives.size() >= 4) {
    H = getHomography(image, secondaryObjectsInformations, error);
    if (!error) {
      for (size_t i = 0; i < TASKBOARD_PRIMARY_OBJECTS_NB_POSITIONS; i++) {
        primaryObjects[i] = RobonautEye::getPerspectiveProjection(H, TASKBOARD_PRIMARY_OBJECTS_POSITIONS[i]);
        if (primaryObjects[i].position.x < 0 || primaryObjects[i].position.y < 0 || primaryObjects[i].position.x > image->getWidth() || primaryObjects[i].position.y > image->getHeight()) {
          primaryObjects[i].position = Point(-1 , -1);
          primaryObjects[i].isCorrect = false;
        } else {
          computePrimaryObjectState(image, primaryObjects[i]);
          primaryObjects[i].isCorrect = true;
        }
      }
      // POST PROCESSING SWITCH POSITION (since it is not on a plan...)
      primaryObjects[2].position = getSwitchPosition(H);
    }
  } else {
    error = true;
  }

  if (error) {
    for (size_t i = 0; i < TASKBOARD_PRIMARY_OBJECTS_NB_POSITIONS; i++) {
      primaryObjects[i].position = Point(-1 , -1);
      primaryObjects[i].isCorrect = false;
    }
  }

  delete secondaryObjectsInformations->objectInformationsMap;
  delete secondaryObjectsInformations;
  return primaryObjects;
}



Image<float> * RobonautEye::convertVectorToImage(vector<int> & imageVector) {
  int height = imageVector.at(0);
  int width = imageVector.at(1);
  size_t length = width * height;
  float * imageData = new float[length];
  size_t i = 0;
  for (vector<int>::iterator it = imageVector.begin() + 2 ; it != imageVector.end(); ++it) {
    imageData[i] = (*it) * 255.0 / 16777215.0;
    i++;
  }

  return new Image<float>(height, width, imageData);
}

ObjectInformationsMapping * RobonautEye::getSecondaryObjectsInformationsFromImage(Image<float> * realImage) {
  ObjectInformationsMapping * mapping = new ObjectInformationsMapping;
  Image<float> * reducedImage = Image<float>::resize(realImage, 800);
  float scale = reducedImage->getWidth() * 1.0 / realImage->getWidth();
  Image<float> * preprocessedImage = preprocessImage(reducedImage, PREPROCESS_BLUR, PREPROCESS_CONTRAST);
  vector < Set<float> * > * sets = segmentImage(preprocessedImage, SEGMENT_THRESHOLD, SEGMENT_MIN_SIZE);

  StandardSelector<float> firstSelector(firstSelectSettings, preprocessedImage->getWidth(), preprocessedImage->getHeight());
  firstSelector.apply(sets);

  StandardMerger<float> merger(mergeSettings, preprocessedImage->getWidth(), preprocessedImage->getHeight());
  merger.apply(sets);

  StandardSelector<float> secondSelector(secondSelectSettings, preprocessedImage->getWidth(), preprocessedImage->getHeight());
  secondSelector.apply(sets);

  int bucketWidth = realImage->getWidth() / AREA_LOCATION_BUCKET_SIZE;
  int bucketHeight = realImage->getHeight() / AREA_LOCATION_BUCKET_SIZE;

  map<int, ObjectInformations> preObjectInformationsMap;

  for (size_t i = 0; i < sets->size(); i++) {
    Set <float> * currentSet = sets->at(i);
    int currentWidth = currentSet->maxX - currentSet->minX;
    int currentHeight = currentSet->maxY - currentSet->minY;
    int currentSize = currentWidth > currentHeight ? currentWidth : currentHeight;

    int centerX = (currentSet->minX + currentWidth / 2.0) / scale;
    int centerY = (currentSet->minY + currentHeight / 2.0) / scale;
    int size = (1 + AREA_MARGIN) * currentSize / scale;
    if (size < AREA_MIN_SIZE) {
      size = AREA_MIN_SIZE;
    }
    ObjectInformations infos;
    infos.centerX = centerX;
    infos.centerY = centerY;
    infos.x = centerX - size / 2;
    infos.y = centerY - size / 2;
    infos.width = size;
    infos.height = size;
    infos.sizeWeight = size;

    correctObject(infos, realImage);

    int fromXBucket = infos.centerX / AREA_LOCATION_BUCKET_SIZE;
    int fromYBucket = infos.centerY / AREA_LOCATION_BUCKET_SIZE;

    for (int x = max(0, fromXBucket - 1); x < min(bucketWidth, fromXBucket + 2); x++) {
      for (int y = max(0, fromYBucket - 1); y < min(bucketHeight, fromYBucket + 2); y++) {
        map<int, ObjectInformations>::iterator it = preObjectInformationsMap.find(x + y * bucketWidth);
        if (it != preObjectInformationsMap.end()) {
          ObjectInformations other = it->second;

          int sumWeight = infos.sizeWeight + other.sizeWeight;
          int overallMinX = infos.x < other.x ? infos.x : other.x;
          int overallMinY = infos.y < other.y ? infos.y : other.y;
          int overallMaxX = infos.x + infos.width > other.x + other.width ? infos.x + infos.width : other.x + other.width;
          int overallMaxY = infos.y + infos.height > other.y + other.height ? infos.y + infos.height : other.y + other.height;
          int overallWidth = overallMaxX - overallMinX;
          int overallHeight = overallMaxY - overallMinY;
          int overallSize = overallWidth > overallHeight ? overallWidth : overallHeight;

          infos.centerX = (infos.centerX * infos.sizeWeight + other.centerX * other.sizeWeight) / sumWeight;
          infos.centerY = (infos.centerY * infos.sizeWeight + other.centerY * other.sizeWeight) / sumWeight;
          infos.x = infos.centerX - overallSize / 2;
          infos.y = infos.centerY - overallSize / 2;
          infos.width = overallSize;
          infos.height = overallSize;
          infos.sizeWeight = sumWeight;

          correctObject(infos, realImage);

          preObjectInformationsMap.erase(it);
        }
      }
    }

    int fromIdBucket = infos.centerX / AREA_LOCATION_BUCKET_SIZE + (infos.centerY / AREA_LOCATION_BUCKET_SIZE) * bucketWidth;
    preObjectInformationsMap.insert(pair<int, ObjectInformations> (fromIdBucket, infos));
  }

  for (size_t i = 0; i < sets->size(); i++) {
    delete sets->at(i);
  }
  sets->clear();

  mapping->objectInformationsMap = new int[bucketWidth * bucketHeight];
  for (int i = 0; i < bucketWidth * bucketHeight; i++) {
    mapping->objectInformationsMap[i] = -2;
  }

  for (map<int, ObjectInformations>::iterator it=preObjectInformationsMap.begin(); it!=preObjectInformationsMap.end(); ++it) {
    addObjectToMapping(it->second, mapping, realImage);
  }


  delete sets;
  delete preprocessedImage;
  delete reducedImage;

  return mapping;
}

void RobonautEye::correctObject(ObjectInformations & infos, Image<float> * realImage) {
  if (infos.x < 0) {
    infos.x = 0;
  }
  if (infos.x + infos.width > realImage->getWidth()) {
    infos.x = realImage->getWidth() - infos.width;
  }
  if (infos.y < 0) {
    infos.y = 0;
  }
  if (infos.y + infos.height > realImage->getHeight()) {
    infos.y = realImage->getHeight() - infos.height;
  }
}

void RobonautEye::addObjectToMapping(ObjectInformations& obj, ObjectInformationsMapping * mapping, Image<float> * image) {

  int bucketX = obj.centerX / AREA_LOCATION_BUCKET_SIZE;
  int bucketY = obj.centerY / AREA_LOCATION_BUCKET_SIZE;
  int bucketWidth = image->getWidth() / AREA_LOCATION_BUCKET_SIZE;
  computeSecondaryObjectClassifierScore(image, obj);
  if (obj.secondaryObjectClassifierScore > SECONDARY_OBJECT_CLASSIFIER_SCORE_THRESHOLD) {
    obj.index = mapping->positives.size();
    mapping->positives.push_back(obj);
    refreshMapping(bucketX, bucketY, mapping, image);
  } else {
    obj.index = -1;
    mapping->objectInformationsMap[bucketX + bucketY * bucketWidth] = -1;
  }

  //mapping->objectInformationsMap[fromIdBucket] = obj.index;
}

void RobonautEye::refreshMapping(int bucketX, int bucketY, ObjectInformationsMapping * mapping, Image<float> * image) {
  int bucketWidth = image->getWidth() / AREA_LOCATION_BUCKET_SIZE;
  int bucketHeight = image->getHeight() / AREA_LOCATION_BUCKET_SIZE;
  int bucketArea = ceil(SECONDARY_OBJECT_SEARCH_AREA / AREA_LOCATION_BUCKET_SIZE) + 1;

  int fromBucketX = max(0, bucketX - bucketArea);
  int toBucketX = min(bucketWidth, bucketX + bucketArea);

  int fromBucketY = max(0, bucketY - bucketArea);
  int toBucketY = min(bucketHeight, bucketY + bucketArea);


  size_t positivesSize = mapping->positives.size();
  for (int x = fromBucketX; x < toBucketX; x++) {
    for (int y = fromBucketY; y < toBucketY; y++) {
      int chosen = -1;
      int bestDistance = -1;
      int realX = x * AREA_LOCATION_BUCKET_SIZE + AREA_LOCATION_BUCKET_SIZE * 0.5;
      int realY = y * AREA_LOCATION_BUCKET_SIZE + AREA_LOCATION_BUCKET_SIZE * 0.5;
      for (size_t i = 0; i < positivesSize; i++) {
        ObjectInformations obj = mapping->positives.at(i);
        float distance = sqrt((realX - obj.centerX) * (realX - obj.centerX) + (realY - obj.centerY) * (realY - obj.centerY));
        if (distance < SECONDARY_OBJECT_SEARCH_AREA && (chosen == -1 || distance < bestDistance)) {
          bestDistance = distance;
          chosen = i;
        }
      }

      if (chosen >= 0) {
        mapping->objectInformationsMap[x + y * bucketWidth] = chosen;
      }
    }
  }
}

Image<float> * RobonautEye::preprocessImage(Image<float> * fromImage, float blur, float contrast) {
  Image<float> * im = ImageSmoother<float>::apply(fromImage, blur);
  Statistics<float>::normalize(im->img, im->getWidth() * im->getHeight(), 128, contrast);
  Statistics<float>::applyLowerBound(im->img, im->getWidth() * im->getHeight(), 0);
  Statistics<float>::applyUpperBound(im->img, im->getWidth() * im->getHeight(), 255);
  return im;
}

vector < Set<float> * > * RobonautEye::segmentImage(Image<float> * image, float segmentThreshold, float segmentMinSize) {
  FelzenHuttenSegmentation<float> segmenter(segmentThreshold, segmentMinSize);
  ImageDisjointSet<float> * disjointSet = segmenter.segmentImage(image);
  vector< Set<float> * > * sets = disjointSet->getSets(image);
  delete disjointSet;

  return sets;
}

void RobonautEye::computeSecondaryObjectClassifierScore(Image<float> * image, ObjectInformations & objectInformations) {
  Image<float> * extract = image->reshape(
        image,
        objectInformations.x,
        objectInformations.y,
        objectInformations.width,
        objectInformations.height,
        SECONDARY_OBJECT_CLASSIFIER_IMAGE_SIZE,
        SECONDARY_OBJECT_CLASSIFIER_IMAGE_SIZE
  );
  float * reduced = secondaryObjectPCA.reduce(extract->img);
  float * classifierScores = secondaryObjectNeuralNetwork.evaluate(reduced);
  objectInformations.imageTypeScores[0] = classifierScores[0];
  objectInformations.imageTypeScores[1] = classifierScores[1];
  objectInformations.imageTypeScores[2] = classifierScores[2];
  objectInformations.imageTypeScores[3] = classifierScores[3];
  objectInformations.secondaryObjectClassifierScore = max(max(classifierScores[1], classifierScores[2]), classifierScores[3]); //
  delete classifierScores;
  delete reduced;
  delete extract;
}

void RobonautEye::computePrimaryObjectState(Image<float> * image, PrimaryObject & primaryObject) {
  Image<float> * extract = getPrimaryObjectImage(image, primaryObject);
  float * reduced = primaryObjectPCA.reduce(extract->img);
  float * evaluation = primaryObjectNeuralNetwork.evaluate(reduced);
  primaryObject.stateScores[0] = evaluation[0];
  primaryObject.stateScores[1] = evaluation[1];
  primaryObject.stateScores[2] = evaluation[2];
  delete evaluation;
  delete reduced;
  delete extract;
}

PrimaryObject RobonautEye::getPerspectiveProjection(Array2D<double> * H, const PrimaryObject & primaryObject) {
  Point points[5] = {
    primaryObject.position, // center
    Point(
      primaryObject.position.x - primaryObject.size / 2,
      primaryObject.position.y - primaryObject.size / 2
    ),
    Point(
      primaryObject.position.x + primaryObject.size / 2,
      primaryObject.position.y - primaryObject.size / 2
    ),
    Point(
      primaryObject.position.x + primaryObject.size / 2,
      primaryObject.position.y + primaryObject.size / 2
    ),
    Point(
      primaryObject.position.x - primaryObject.size / 2,
      primaryObject.position.y + primaryObject.size / 2
    )
  };

  Point * projectedPoints = Geometry::getPerspectiveProjectionPoints(H, points, 5);
  double x = projectedPoints[0].x;
  double y = projectedPoints[0].y;

  double squaredSize = 0;
  for (size_t i = 1; i < 5; i++) {
    double squaredSizeItem = (projectedPoints[i].x - x) * (projectedPoints[i].x - x) + (projectedPoints[i].y - y) * (projectedPoints[i].y - y);
    if (squaredSizeItem > squaredSize) {
      squaredSize = squaredSizeItem;
    }
  }

  delete projectedPoints;
  return PrimaryObject(x, y, max(30.0, sqrt(squaredSize) * 2.0 / sqrt(2)));
}

Image<float> * RobonautEye::getPrimaryObjectImage(Image<float> * image, PrimaryObject & primaryObject) {
  int x = primaryObject.position.x - primaryObject.size / 2;
  int y = primaryObject.position.y - primaryObject.size / 2;
  int size = primaryObject.size;
  if (x < 0) {
    x = 0;
  }
  if (x + size > image->getWidth()) {
    x = image->getWidth() - size;
  }
  if (y < 0) {
    y = 0;
  }
  if (y + size > image->getHeight()) {
    y = image->getHeight() - size;
  }
  Image<float> * extract = image->reshape(
        image,
        x,
        y,
        size,
        size,
        PRIMARY_OBJECT_CLASSIFIER_IMAGE_SIZE,
        PRIMARY_OBJECT_CLASSIFIER_IMAGE_SIZE
  );
  return extract;
}


ProjectionVector RobonautEye::getOrthographicProjectionVector(OrthographicProjectionEquation & eq) {
  ProjectionVector v;
  Point * points = new Point[2];
  points[0] = Point(0, 0);
  points[1] = Point(1, 0);
  Point * projections = Geometry::getOrthographicProjectionPoints(eq, points, 2);
  v.angle = atan2(projections[1].y - projections[0].y, projections[1].x - projections[0].x);
  v.distance = sqrt((projections[1].y - projections[0].y) * (projections[1].y - projections[0].y) + (projections[1].x - projections[0].x) * (projections[1].x - projections[0].x));
  delete projections;
  delete points;

  return v;
}

bool RobonautEye::isSetCorrect(Point & fromA, Point & fromB, Point & fromC, Point & toA, Point & toB, Point & toC) {
  float * infos = getOrthographicProjectionInformations(fromA, fromB, fromC, toA, toB, toC);
  bool isCorrect = true;

  if (infos[8] < 0) {
    isCorrect = false;
  }

  if (infos[0] < 0.5 || infos[0] > 3.25) {
    isCorrect = false;
  }

  if (infos[1] < -1 || infos[1] > 1) {
    isCorrect = false;
  }

  if (infos[2] < 0.25 || infos[2] > 2.5) {
    isCorrect = false;
  }

  if (infos[3] < -0.5 || infos[3] > 0.75) {
    isCorrect = false;
  }

  /*
  if (infos[4] < -50 || infos[4] > 1150) {
    isCorrect = false;
  }
  */

  if (infos[5] < 0.25 || infos[5] > 2) {
    isCorrect = false;
  }

  if (infos[6] < -1.5 || infos[6] > 2) {
    isCorrect = false;
  }

  /*
  if (infos[7] < -100 || infos[7] > 1300) {
    isCorrect = false;
  }
  */

  delete infos;
  return isCorrect;
}

bool RobonautEye::isSetCorrect(Point & fromA, Point & fromB, Point & fromC, Point & fromD, Point & toA, Point & toB, Point & toC, Point & toD) {
  float * infos = getOrthographicProjectionInformations(fromA, fromB, fromC, fromD, toA, toB, toC, toD);

  bool isCorrect = true;

  if (infos[8] < 0) {
    isCorrect = false;
  }

  if (infos[0] > 1.4) {
    isCorrect = false;
  }

  if (infos[1] > 0.8) {
    isCorrect = false;
  }

  if (infos[2] > 1.4) {
    isCorrect = false;
  }

  if (infos[3] > 0.8) {
    isCorrect = false;
  }

  if (infos[4] > 350) {
    isCorrect = false;
  }

  if (infos[5] > 0.8) {
    isCorrect = false;
  }

  if (infos[6] > 2.2) {
    isCorrect = false;
  }

  if (infos[7] > 750) {
    isCorrect = false;
  }

  delete infos;
  return isCorrect;
}

float * RobonautEye::getOrthographicProjectionInformations(Point & fromA, Point & fromB, Point & fromC, Point & toA, Point & toB, Point & toC) {
  float * infos = new float[9];
  OrthographicProjectionEquation eq = Geometry::getOrthographicProjectionEquation(fromA, fromB, fromC, toA, toB, toC);
  if (!eq.xEq.correct || !eq.yEq.correct) {
    infos[8] = -1;
    return infos;
  }
  ProjectionVector v = RobonautEye::getOrthographicProjectionVector(eq);
  infos[0] = v.distance;
  infos[1] = v.angle;
  infos[2] = eq.xEq.a;
  infos[3] = eq.xEq.b;
  infos[4] = eq.xEq.c;
  infos[5] = eq.yEq.a;
  infos[6] = eq.yEq.b;
  infos[7] = eq.yEq.c;
  infos[8] = 1;
  return infos;
}

float * RobonautEye::getOrthographicProjectionInformations(Point & fromA, Point & fromB, Point & fromC, Point & fromD, Point & toA, Point & toB, Point & toC, Point & toD) {
  float * infos = new float[9];
  float ** triplesInfos = new float *[4];
  triplesInfos[0] = getOrthographicProjectionInformations(fromA, fromB, fromC, toA, toB, toC);
  triplesInfos[1] = getOrthographicProjectionInformations(fromA, fromB, fromD, toA, toB, toD);
  triplesInfos[2] = getOrthographicProjectionInformations(fromA, fromC, fromD, toA, toC, toD);
  triplesInfos[3] = getOrthographicProjectionInformations(fromB, fromC, fromD, toB, toC, toD);

  if (triplesInfos[0][8] < 0 || triplesInfos[1][8] < 0 || triplesInfos[2][8] < 0 || triplesInfos[3][8] < 0) {
    infos[8] = -1;
    return infos;
  }
  for (size_t i = 0; i < 8; i++) {
    float * list = new float[4];
    list[0] = 0;
    list[1] = 0;
    list[2] = 0;
    list[3] = 0;
    for (size_t j = 0; j < 4; j++) {
      list[j] += triplesInfos[j][i];
    }
    infos[i] = sqrt(Statistics<float>::variance(list, 4));
    delete list;
  }
  infos[8] = 1;

  delete triplesInfos[0];
  delete triplesInfos[1];
  delete triplesInfos[2];
  delete triplesInfos[3];
  delete triplesInfos;
  return infos;
}


Array2D<double> * RobonautEye::getHomography(Image<float> * realImage, ObjectInformationsMapping * infos, bool & error) {
  // HEAP
  error = false;
  multiset<QuadrupleIndexes> * quadruples = getQuadruples(infos->positives);

  SecondaryObjectAssociations bestAssociation;
  bool started = false;
  while (quadruples->size() > 0) {
    multiset<QuadrupleIndexes>::iterator it = quadruples->begin();
    QuadrupleIndexes quadruple = (*it);
    quadruples->erase(quadruples->begin());
    SecondaryObjectAssociations association = findClosestQuadruples(realImage, quadruple, infos);
    if (!started || association.score > bestAssociation.score) {
      started = true;
      bestAssociation = association;
    }
    if (bestAssociation.score > lowerBoundScore(infos)) {
      break;
    }
  }

  if (bestAssociation.score < 400) {
    error = true;
    return new Array2D<double>();
  }

  // todo here
  bool satisfied = false;
  float score = 0;
  Array2D<double> * H = getPurifiedHomography(realImage, infos, bestAssociation, bestAssociation.score * SECONDARY_OBJECT_HOMOGRAPHY_PURIFY_MIN_FACTOR, SECONDARY_OBJECT_HOMOGRAPHY_PURIFY_MAX_POINTS, score, satisfied);

  return H;
}

Array2D<double> * RobonautEye::getPurifiedHomography(Image<float> * realImage, ObjectInformationsMapping * infos, SecondaryObjectAssociations association, float scoreThreshold, size_t maxDepth, float &score, bool &satisfied) {
  vector<Point> pointsFrom;
  vector<Point> pointsTo;
  for (size_t i = 0; i < TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS; i++) {
    if (association.associations[i] >= 0) {
      ObjectInformations objectInfos = infos->positives[association.associations[i]];
      pointsFrom.push_back(TASKBOARD_SECONDARY_OBJECTS_POSITIONS[i].position);
      pointsTo.push_back(Point(objectInfos.centerX, objectInfos.centerY));
    }
  }
  return getPurifiedHomography(realImage, infos, pointsFrom, pointsTo, scoreThreshold, maxDepth, score, satisfied);
}

Array2D<double> * RobonautEye::getPurifiedHomography(Image<float> * realImage, ObjectInformationsMapping * infos, vector<Point> &pointsFrom, vector<Point> &pointsTo, float scoreThreshold, size_t maxDepth, float &score, bool &satisfied) {
  Array2D<double> *H = Geometry::getPerspectiveProjectionHomography(&pointsFrom[0], &pointsTo[0], pointsFrom.size());
  SecondaryObjectAssociations resultAssociation = evaluateHomography(realImage, H, infos);
  score = resultAssociation.score;
  if (score < scoreThreshold) {
    satisfied = false;
    if (pointsFrom.size() < 5 || maxDepth == 0) {
      return H;
    }
    for (size_t firstLevel = 0; firstLevel < 2; firstLevel++) {
      for (size_t i = 0; i < pointsFrom.size(); i++) {
        vector<Point> pointsFrom2 = pointsFrom;
        vector<Point> pointsTo2 = pointsTo;
        pointsFrom2.erase(pointsFrom2.begin() + i);
        pointsTo2.erase(pointsTo2.begin() + i);
        bool resSatisfied = false;
        float resScore = 0;
        Array2D<double> *res = getPurifiedHomography(realImage, infos, pointsFrom2, pointsTo2, scoreThreshold, firstLevel == 0 ? 0 : maxDepth - 1, resScore, resSatisfied);
        if (resScore > score) {
          score = resScore;
          delete H;
          H = res;
          if (resSatisfied) {
            return H;
          }
        } else {
          delete res;
        }
      }
    }
  } else {
    satisfied = true;
  }
  return H;
}

multiset<QuadrupleIndexes> * RobonautEye::getQuadruples(vector <ObjectInformations> & points) {
  // HEAP ?
  multiset<QuadrupleIndexes> * quadruples = new multiset<QuadrupleIndexes>();

  ObjectInformations * pointsList = &points[0];
  size_t nbPoints = points.size();

  for (size_t i = 0; i < nbPoints - 3; i++) {
    for (size_t j = i + 1; j < nbPoints - 2; j++) {
      for (size_t k = j + 1; k < nbPoints - 1; k++) {
        for (size_t l = k + 1; l < nbPoints; l++) {
          float centerX = (pointsList[i].centerX + pointsList[j].centerX + pointsList[k].centerX + pointsList[l].centerX) / 4;
          float centerY = (pointsList[i].centerY + pointsList[j].centerY + pointsList[k].centerY + pointsList[l].centerY) / 4;
          QuadrupleIndexes quadruple;
          quadruple.i = i;
          quadruple.j = j;
          quadruple.k = k;
          quadruple.l = l;
          quadruple.score = sqrt((centerX - pointsList[i].centerX) * (centerX - pointsList[i].centerX) + (centerY - pointsList[i].y) * (centerY - pointsList[i].y) +
              (centerX - pointsList[j].centerX) * (centerX - pointsList[j].centerX) + (centerY - pointsList[j].centerY) * (centerY - pointsList[j].centerY) +
              (centerX - pointsList[k].centerX) * (centerX - pointsList[k].centerX) + (centerY - pointsList[k].centerY) * (centerY - pointsList[k].centerY) +
              (centerX - pointsList[l].centerX) * (centerX - pointsList[l].centerX) + (centerY - pointsList[l].centerY) * (centerY - pointsList[l].centerY));
          if (quadruple.score > SECONDARY_OBJECT_HOMOGRAPHY_MIN_SIZE_QUADRUPLE) {
            quadruples->insert(quadruple);
          }
        }
      }
    }
  }

  return quadruples;
}

SecondaryObjectAssociations RobonautEye::findClosestQuadruples(Image<float> * realImage, const QuadrupleIndexes & quadrupleFrom, ObjectInformationsMapping * infos) {
  SecondaryObjectAssociations bestAssociation;
  bool started = false;

  Point to[4];
  Point from[4];
  to[0] = Point(infos->positives[quadrupleFrom.i].centerX, infos->positives[quadrupleFrom.i].centerY);
  to[1] = Point(infos->positives[quadrupleFrom.j].centerX, infos->positives[quadrupleFrom.j].centerY);
  to[2] = Point(infos->positives[quadrupleFrom.k].centerX, infos->positives[quadrupleFrom.k].centerY);
  to[3] = Point(infos->positives[quadrupleFrom.l].centerX, infos->positives[quadrupleFrom.l].centerY);

  size_t total = 0;
  for (size_t i = 0; i < TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS; i++) {
    if (infos->positives[quadrupleFrom.i].imageTypeScores[TASKBOARD_SECONDARY_OBJECTS_POSITIONS[i].imageType] < SECONDARY_OBJECT_CLASSIFIER_SCORE_THRESHOLD) {
      continue;
    }
    from[0] = TASKBOARD_SECONDARY_OBJECTS_POSITIONS[i].position;
    for (size_t j = 0; j < TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS; j++) {
      if (j == i) continue;
      if (infos->positives[quadrupleFrom.j].imageTypeScores[TASKBOARD_SECONDARY_OBJECTS_POSITIONS[j].imageType] < SECONDARY_OBJECT_CLASSIFIER_SCORE_THRESHOLD) {
        continue;
      }
      from[1] = TASKBOARD_SECONDARY_OBJECTS_POSITIONS[j].position;
      for (size_t k = 0; k < TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS; k++) {
        if (k == j || k == i) continue;
        if (infos->positives[quadrupleFrom.k].imageTypeScores[TASKBOARD_SECONDARY_OBJECTS_POSITIONS[k].imageType] < SECONDARY_OBJECT_CLASSIFIER_SCORE_THRESHOLD) {
          continue;
        }
        from[2] = TASKBOARD_SECONDARY_OBJECTS_POSITIONS[k].position;
        if (isSetCorrect(from[0], from[1], from[2], to[0], to[1], to[2])) {
          for (size_t l = 0; l < TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS; l++) {
            if (l == k || l == j || l == i) continue;
            if (infos->positives[quadrupleFrom.l].imageTypeScores[TASKBOARD_SECONDARY_OBJECTS_POSITIONS[l].imageType] < SECONDARY_OBJECT_CLASSIFIER_SCORE_THRESHOLD) {
              continue;
            }
            from[3] = TASKBOARD_SECONDARY_OBJECTS_POSITIONS[l].position;
            if (
                isSetCorrect(from[0], from[1], from[3], to[0], to[1], to[3]),
                isSetCorrect(from[0], from[2], from[3], to[0], to[2], to[3]),
                isSetCorrect(from[1], from[2], from[3], to[1], to[2], to[3]),
                isSetCorrect(from[0], from[1], from[2], from[3], to[0], to[1], to[2], to[3])
                ) {
              total++;

              SecondaryObjectAssociations association;
              Array2D<double> * H = Geometry::getPerspectiveProjectionHomography(&from[0], &to[0], 4);
              association = evaluateHomography(realImage, H, infos);
              if (association.score > SECONDARY_OBJECT_HOMOGRAPHY_SUCCESS_LOWERBOUND) {
                float score = 0;
                bool satisfied = false;
                Array2D<double> * purifiedH = getPurifiedHomography(realImage, infos, association, 10000, SECONDARY_OBJECT_HOMOGRAPHY_PURIFY_MAX_POINTS, score, satisfied);
                association = evaluateHomography(realImage, purifiedH, infos);
                delete purifiedH;
              }
              delete H;

              if (!started || association.score > bestAssociation.score) {
                bestAssociation = association;
                started = true;
              }

              if (bestAssociation.score > lowerBoundScore(infos)) {
                return bestAssociation;
              }
            }
          }
        }
      }
    }
  }

  if (!started) {
    bestAssociation.score = 0;
  }

  return bestAssociation;
}

SecondaryObjectAssociations RobonautEye::evaluateHomography(Image<float> * realImage, Array2D<double> * H, ObjectInformationsMapping * infos) {
  Point * projected = Geometry::getPerspectiveProjectionPoints(H, &TASKBOARD_SECONDARY_OBJECTS_POSITIONS_POINTS[0], TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS);
  SecondaryObjectAssociations associations;
  associations.score = 0;
  vector <bool> * blacklist = new vector<bool>(infos->positives.size(), false);
  for (size_t i = 0; i < TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS; i++) {
    int selectedId = searchForSecondaryObjectOnArea(realImage, projected[i], TASKBOARD_SECONDARY_OBJECTS_POSITIONS[i].imageType, blacklist, infos);
    associations.associations[i] = selectedId;

    if (selectedId >= 0) {
      (*blacklist)[selectedId] = true;
      ObjectInformations selectedObject = infos->positives.at(selectedId);
      float distance = sqrt((projected[i].x - selectedObject.centerX) * (projected[i].x - selectedObject.centerX) + (projected[i].y - selectedObject.centerY) * (projected[i].y - selectedObject.centerY));
      float score = (SECONDARY_OBJECT_SEARCH_AREA - distance) * 1.0 / SECONDARY_OBJECT_SEARCH_AREA;
      associations.score += score * score * 100.0;
    }
  }

  /*
  for (size_t i = 0; i < TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS; i++) {
    for (size_t j = i + 1; j < TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS; j++) {
      float distance = sqrt((projected[i].x - projected[j].x) * (projected[i].x - projected[j].x) + (projected[i].y - projected[j].y) * (projected[i].y - projected[j].y));
      if (distance < 20) {
        associations.score /= 2;
      }
    }
  }
  */
  delete blacklist;

  delete projected;
  return associations;
}

int RobonautEye::searchForSecondaryObjectOnArea(Image<float> * realImage, Point & point, int imageType, vector<bool> * blacklist, ObjectInformationsMapping * infos) {
  int bucketX = point.x / AREA_LOCATION_BUCKET_SIZE;
  int bucketY = point.y / AREA_LOCATION_BUCKET_SIZE;
  int bucketImageWidth = realImage->getWidth() / AREA_LOCATION_BUCKET_SIZE;
  //qDebug() << "BUCKET WIDTH SECONDARY OBJECT" << bucketImageWidth;
  int bucketImageHeight = realImage->getHeight() / AREA_LOCATION_BUCKET_SIZE;
  int searchArea = ceil(SECONDARY_OBJECT_SEARCH_AREA * 1.0 / AREA_LOCATION_BUCKET_SIZE) + 1;
/*
  int chosen = -1;
  float bestDistance = -1;
  for (size_t i = 0; i < infos->positives.size(); i++) {
    if ((*blacklist)[i]) {
      continue;
    }
    ObjectInformations positive = infos->positives.at(i);
    if (positive.imageTypeScores[imageType] < SECONDARY_OBJECT_CLASSIFIER_SCORE_THRESHOLD) {
      continue;
    }
    float distance = sqrt((point.x - positive.centerX) * (point.x - positive.centerX) + (point.y - positive.centerY) * (point.y - positive.centerY));
    if (distance < SECONDARY_OBJECT_SEARCH_AREA && (chosen == -1 || distance < bestDistance)) {
      bestDistance = distance;
      chosen = i;
    }
  }
  //return chosen;
*/
  int chosenSearchArea = -1;
  float bestDistance = -1;

  for (int i = 0; i < searchArea; i++) {
    for (int x = bucketX - i; x < bucketX + i + 1; x++) {
      for (int y = bucketY - i; y < bucketY + i + 1; y++) {
        if (y == bucketY - i || y == bucketY + i || x == bucketX - i || x == bucketX + i) {
          if (x >= 0 && y >= 0 && x < bucketImageWidth && y < bucketImageHeight) {
            int fromIdBucket = x + y * bucketImageWidth;

            int id = infos->objectInformationsMap[fromIdBucket];
            if (id < 0) {
              continue;
            }
            ObjectInformations positive = infos->positives[id];
            float distance = sqrt((point.x - positive.centerX) * (point.x - positive.centerX) + (point.y - positive.centerY) * (point.y - positive.centerY));
            if (positive.imageTypeScores[imageType] > SECONDARY_OBJECT_CLASSIFIER_SCORE_THRESHOLD && (!(*blacklist)[positive.index]) && distance < SECONDARY_OBJECT_SEARCH_AREA && (chosenSearchArea == -1 || distance < bestDistance)) {
              bestDistance = distance;
              chosenSearchArea = positive.index;
            }
          }
        }
      }
    }
  }

/* For this to work Secondary Object neural networks and PCA need to be trained to process any extract on the image...
  if (chosenSearchArea == -1) {
    size_t sizeBefore = infos->positives.size();
    searchForUndetectedSecondaryObjects(realImage, point, blacklist, infos);
    if (infos->positives.size() != sizeBefore) {
      qDebug() << "relaunched";
      searchForSecondaryObjectOnArea(realImage, point, imageType, blacklist, infos);
    }
  }
*/
  /*
  if (chosen != chosenSearchArea) {
    float distance = sqrt((point.x - infos->positives.at(chosen).centerX) * (point.x - infos->positives.at(chosen).centerX) + (point.y - infos->positives.at(chosen).centerY) * (point.y - infos->positives.at(chosen).centerY));
    qDebug() << "DIFF: " << point.x << point.y << bucketX << bucketY << chosen << infos->positives.at(chosen).centerX << infos->positives.at(chosen).centerY << chosenSearchArea << distance;
  }
*/

  return chosenSearchArea;
}

/*
void RobonautEye::searchForUndetectedSecondaryObjects(Image<float> * realImage, Point & point, vector<bool> * blacklist, ObjectInformationsMapping * mapping) {
  int bucketX = point.x / AREA_LOCATION_BUCKET_SIZE;
  int bucketY = point.y / AREA_LOCATION_BUCKET_SIZE;
  int bucketImageWidth = realImage->getWidth() / AREA_LOCATION_BUCKET_SIZE;
  int bucketImageHeight = realImage->getHeight() / AREA_LOCATION_BUCKET_SIZE;
  int searchArea = ceil(SECONDARY_OBJECT_SEARCH_AREA * 1.0 / AREA_LOCATION_BUCKET_SIZE);

  float meanSize = 0;
  for (vector<ObjectInformations>::iterator it = mapping->positives.begin() ; it != mapping->positives.end(); ++it) {
    meanSize += ((*it).height + (*it).width) / (2.0 * mapping->positives.size());
  }



  for (int i = 0; i < searchArea; i++) {
    for (int x = bucketX - i; x < bucketX + i + 1; x++) {
      for (int y = bucketY - i; y < bucketY + i + 1; y++) {
        if (y == bucketY - i || y == bucketY + i || x == bucketX - i || x == bucketX + i) {
          if (x >= 0 && y >= 0 && x < bucketImageWidth && y < bucketImageHeight) {
            int fromIdBucket = x + y * bucketImageWidth;
            if (mapping->objectInformationsMap[fromIdBucket] == -2) {
              ObjectInformations info;
              info.centerX = bucketX * AREA_LOCATION_BUCKET_SIZE + AREA_LOCATION_BUCKET_SIZE * 0.5;
              info.centerY = bucketY * AREA_LOCATION_BUCKET_SIZE + AREA_LOCATION_BUCKET_SIZE * 0.5;
              info.width = meanSize;
              info.height = meanSize;
              info.sizeWeight = 1;
              correctObject(info, realImage);
              computeSecondaryObjectClassifierScore(realImage, info);

              if (info.secondaryObjectClassifierScore > SECONDARY_OBJECT_CLASSIFIER_SCORE_THRESHOLD) {
                // addObjectToMapping(info, mapping, realImage);
                vector<Point> points;
                locateUndetectedSecondaryObject(realImage, x, y, meanSize, true, mapping, points);
                Point meanPoint;
                meanPoint.x = 0;
                meanPoint.y = 0;
                for (vector<Point>::iterator it = points.begin() ; it != points.end(); ++it) {
                  meanPoint.x += (*it).x / points.size();
                  meanPoint.y += (*it).y / points.size();
                  mapping->objectInformationsMap[((int) ((*it).x / AREA_LOCATION_BUCKET_SIZE)) + ((int)((*it).y / AREA_LOCATION_BUCKET_SIZE)) * bucketImageWidth] = -1;
                }
                ObjectInformations info2;
                info2.centerX = meanPoint.x;
                info2.centerY = meanPoint.y;
                info2.width = meanSize;
                info2.height = meanSize;
                info2.sizeWeight = 1;
                correctObject(info2, realImage);

                qDebug() << "before";
                addObjectToMapping(info2, mapping, realImage);
                qDebug() << "POINTS" << points.size() << meanPoint.x << meanPoint.y << mapping->positives.size();
              } else {
                mapping->objectInformationsMap[x + y * bucketImageWidth] = -1;
              }
            }
          }
        }
      }
    }
  }
}

void RobonautEye::locateUndetectedSecondaryObject(Image<float> * realImage, int bucketX, int bucketY, int meanSize, bool alreadyDetected, ObjectInformationsMapping * mapping, vector<Point> & points) {
  int bucketImageWidth = realImage->getWidth() / AREA_LOCATION_BUCKET_SIZE;
  int bucketImageHeight = realImage->getHeight() / AREA_LOCATION_BUCKET_SIZE;
  if (mapping->objectInformationsMap[bucketX + bucketY * bucketImageWidth] == -2) {
    if (!alreadyDetected) {
      ObjectInformations info;
      info.centerX = bucketX * AREA_LOCATION_BUCKET_SIZE + AREA_LOCATION_BUCKET_SIZE * 0.5;
      info.centerY = bucketY * AREA_LOCATION_BUCKET_SIZE + AREA_LOCATION_BUCKET_SIZE * 0.5;
      info.width = meanSize;
      info.height = meanSize;
      info.sizeWeight = 1;
      correctObject(info, realImage);
      computeSecondaryObjectClassifierScore(realImage, info);
      if (info.secondaryObjectClassifierScore < SECONDARY_OBJECT_CLASSIFIER_SCORE_THRESHOLD) {
        mapping->objectInformationsMap[bucketX + bucketY * bucketImageWidth] = -1;
        return;
      }
    }
    points.push_back(Point(bucketX * AREA_LOCATION_BUCKET_SIZE + AREA_LOCATION_BUCKET_SIZE * 0.5, bucketY * AREA_LOCATION_BUCKET_SIZE + AREA_LOCATION_BUCKET_SIZE * 0.5));
    mapping->objectInformationsMap[bucketX + bucketY * bucketImageWidth] = -3;

    if (bucketX > 0) {
      locateUndetectedSecondaryObject(realImage, bucketX - 1, bucketY, meanSize, false, mapping, points);
    }

    if (bucketY > 0) {
      locateUndetectedSecondaryObject(realImage, bucketX, bucketY - 1, meanSize, false, mapping, points);
    }

    if (bucketX < bucketImageWidth - 1) {
      locateUndetectedSecondaryObject(realImage, bucketX + 1, bucketY, meanSize, false, mapping, points);
    }

    if (bucketY < bucketImageHeight - 1) {
      locateUndetectedSecondaryObject(realImage, bucketX, bucketY + 1, meanSize, false, mapping, points);
    }
  }
}

*/

void RobonautEye::getPrimaryItemInformations(size_t id, Array2D<double> * H, float &length1, float &length2, float &angle1, float &angle2, Point &toPosition) {
  Point * pointsFrom = new Point[3];
  pointsFrom[0].x = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[id].position.x + 1;
  pointsFrom[0].y = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[id].position.y;
  pointsFrom[1].x = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[id].position.x;
  pointsFrom[1].y = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[id].position.y + 1;
  pointsFrom[2].x = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[id].position.x;
  pointsFrom[2].y = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[id].position.y;
  Point * pointsTo = Geometry::getPerspectiveProjectionPoints(H, pointsFrom, 3);
  length1 = sqrt((pointsTo[0].x - pointsTo[2].x) * (pointsTo[0].x - pointsTo[2].x) + (pointsTo[0].y - pointsTo[2].y) * (pointsTo[0].y - pointsTo[2].y));
  length2 = sqrt((pointsTo[1].x - pointsTo[2].x) * (pointsTo[1].x - pointsTo[2].x) + (pointsTo[1].y - pointsTo[2].y) * (pointsTo[1].y - pointsTo[2].y));
  angle1 = atan2(pointsTo[0].y - pointsTo[2].y, pointsTo[0].x - pointsTo[2].x);
  angle2 = atan2(pointsTo[1].y - pointsTo[2].y, pointsTo[1].x - pointsTo[2].x);
  toPosition = pointsTo[2];
  delete pointsTo;
  delete pointsFrom;
}

Point RobonautEye::getSwitchPosition(Array2D<double> * H) {
  float length1;
  float length2;
  float angle1;
  float angle2;
  Point toPosition;
  getPrimaryItemInformations(2, H, length1, length2, angle1, angle2, toPosition);

  Point newPosition = getSwitchDecal(length1, length2, angle1, angle2);
  newPosition.x = toPosition.x + newPosition.x;
  newPosition.y = toPosition.y + newPosition.y;
  return newPosition;
}

Point RobonautEye::getSwitchDecal(float length1, float length2, float angle1, float angle2) {
  float v1[3];
  float v2[5];
  v1[0] = 1;
  v1[1] = length1;
  v1[2] = length2;
  v2[0] = 1;
  v2[1] = cos(angle1);
  v2[2] = sin(angle1);
  v2[3] = cos(angle2);
  v2[4] = sin(angle2);
  Point p;
  p.x = -1162.530163 * v1[0] * v2[0] + 1338.110485 * v1[0] * v2[1] + 20.627228 * v1[0] * v2[2] + 24.101029 * v1[0] * v2[3] + -164.312769 * v1[0] * v2[4] + -803.866034 * v1[1] * v2[0] + -107.277511 * v1[1] * v2[1] + -72.779681 * v1[1] * v2[2] + -65.497817 * v1[1] * v2[3] + 888.796113 * v1[1] * v2[4] + 2369.954027 * v1[2] * v2[0] + -1574.047412 * v1[2] * v2[1] + 38.832119 * v1[2] * v2[2] + 15.243031 * v1[2] * v2[3] + -780.077884 * v1[2] * v2[4];
  p.y = 307.688958 * v1[0] * v2[0] + -357.018157 * v1[0] * v2[1] + -14.131346 * v1[0] * v2[2] + -10.046833 * v1[0] * v2[3] + 38.359596 * v1[0] * v2[4] + -344.779204 * v1[1] * v2[0] + 68.560833 * v1[1] * v2[1] + -42.566229 * v1[1] * v2[2] + -8.407327 * v1[1] * v2[3] + 316.312812 * v1[1] * v2[4] + 1051.428428 * v1[2] * v2[0] + -519.855443 * v1[2] * v2[1] + 51.830521 * v1[2] * v2[2] + 12.563200 * v1[2] * v2[3] + -558.599816 * v1[2] * v2[4];
  return p;
}

Point RobonautEye::getPowerSwitchPositionDown(Array2D<double> * H) {
  float length1;
  float length2;
  float angle1;
  float angle2;
  Point toPosition;
  getPrimaryItemInformations(0, H, length1, length2, angle1, angle2, toPosition);

  Point newPosition = getPowerSwitchDecalDown(length1, length2, angle1, angle2);
  newPosition.x = toPosition.x + newPosition.x;
  newPosition.y = toPosition.y + newPosition.y;

  return newPosition;
}

Point RobonautEye::getPowerSwitchDecalDown(float length1, float length2, float angle1, float angle2) {
  float v1[2];
  float v2[4];
  v1[0] = 1;
  v1[1] = length1;
  v2[0] = 1;
  v2[1] = sin(angle1);
  v2[2] = cos(angle2);
  v2[3] = sin(angle2);
  Point p;
  p.x = 695.297033 * v1[0] * v2[0] + -0.172132 * v1[0] * v2[1] + 119.311219 * v1[0] * v2[2] + -688.254764 * v1[0] * v2[3] + -551.801175 * v1[1] * v2[0] + -103.996969 * v1[1] * v2[1] + -122.301928 * v1[1] * v2[2] + 550.213253 * v1[1] * v2[3];
  p.y = -237.759693 * v1[0] * v2[0] + 13.174071 * v1[0] * v2[1] + 34.315789 * v1[0] * v2[2] + 250.527486 * v1[0] * v2[3] + 392.744412 * v1[1] * v2[0] + -12.576138 * v1[1] * v2[1] + -17.625509 * v1[1] * v2[2] + -336.630563 * v1[1] * v2[3];
  return p;
}


Point RobonautEye::getPowerSwitchPositionUp(Array2D<double> * H) {
  float length1;
  float length2;
  float angle1;
  float angle2;
  Point toPosition;
  getPrimaryItemInformations(0, H, length1, length2, angle1, angle2, toPosition);

  Point newPosition = getPowerSwitchDecalUp(length1, length2, angle1, angle2);
  newPosition.x = toPosition.x + newPosition.x;
  newPosition.y = toPosition.y + newPosition.y;

  return newPosition;
}

Point RobonautEye::getPowerSwitchDecalUp(float length1, float length2, float angle1, float angle2) {
  float v1[3];
  float v2[5];
  v1[0] = 1;
  v1[1] = length1;
  v1[2] = length2;
  v2[0] = 1;
  v2[1] = sin(angle1);
  v2[2] = cos(angle1);
  v2[3] = cos(angle2);
  v2[4] = sin(angle2);
  Point p;
  p.x = 2130.138092 * v1[0] * v2[0] + -62.123833 * v1[0] * v2[1] + 47.104116 * v1[0] * v2[2] + -150.578554 * v1[0] * v2[3] + -2195.233041 * v1[0] * v2[4] + 8099.645522 * v1[1] * v2[0] + -1345.430946 * v1[1] * v2[1] + -5679.577326 * v1[1] * v2[2] + 131.205217 * v1[1] * v2[3] + -2609.461572 * v1[1] * v2[4] + -11234.849320 * v1[2] * v2[0] + 1424.535557 * v1[2] * v2[1] + 6714.490017 * v1[2] * v2[2] + -242.492537 * v1[2] * v2[3] + 4747.174675 * v1[2] * v2[4];
  p.y = 5733.545206 * v1[0] * v2[0] + 139.229432 * v1[0] * v2[1] + -2606.898359 * v1[0] * v2[2] + -291.505887 * v1[0] * v2[3] + -3253.153451 * v1[0] * v2[4] + 6327.662467 * v1[1] * v2[0] + 120.022144 * v1[1] * v2[1] + -3074.397064 * v1[1] * v2[2] + -378.087673 * v1[1] * v2[3] + -3134.510620 * v1[1] * v2[4] + -9555.393618 * v1[2] * v2[0] + -231.000119 * v1[2] * v2[1] + 4378.151795 * v1[2] * v2[2] + 660.937854 * v1[2] * v2[3] + 5086.882153 * v1[2] * v2[4];
  return p;
}
