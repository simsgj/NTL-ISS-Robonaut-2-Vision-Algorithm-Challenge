#include "utility.h"

Utility::Utility()
{
}


vector<int> Utility::loadImage(QString path) {
  QImage img(path);

  int width = img.width();
  int height = img.height();

  vector<int> result(width * height + 2);
  result[0] = height;
  result[1] = width;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      QRgb rgb = img.pixel(x, y);
      result[x + y * width + 2] = qRed(rgb) * 65536 + 256 * qGreen(rgb) + qBlue(rgb);
    }
  }
  return result;
}

map <QString, vector <RealObjectInformations> > Utility::loadCSV(QString path) {
  map <QString, vector <RealObjectInformations> > informations;

  QFile file(path);
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
      return informations;
  }

  bool firstLine = true;
  while (!file.atEnd()) {
      QString line = QString(file.readLine());
      if (firstLine == true) {
        firstLine = false;
        continue;
      }
      QStringList lineList = line.split(',');
      QString filename = lineList.at(0) + "/" + lineList.at(1);
      vector <RealObjectInformations> realInformations;
      for (int i = 0; i < (lineList.size() - 2) / 5; ++i) {
        RealObjectInformations info;
        info.state = lineList.at(i * 5 + 2);
        info.leftX = lineList.at(i * 5 + 3).toInt();
        info.leftY = lineList.at(i * 5 + 4).toInt();
        info.rightX = lineList.at(i * 5 + 5).toInt();
        info.rightY = lineList.at(i * 5 + 6).toInt();
        realInformations.push_back(info);
      }
      informations.insert(pair<QString, vector <RealObjectInformations> >(filename, realInformations));
  }

  return informations;
}

void Utility::selectSecondaryObjects(vector <RealObjectInformations> & primaryObjects) {
  primaryObjects.erase(primaryObjects.begin()); // PANER_POWER_SWITCH
  primaryObjects.erase(primaryObjects.begin()); // PANEL_POWER_COVER
  // PANEL_POWER_LED
  primaryObjects.erase(primaryObjects.begin() + 1); // A01_ROCKER_SWITCH
  // A01_ROCKER_LED_TOP
  // A01_ROCKER_LED_BOTTOM
  primaryObjects.erase(primaryObjects.begin() + 3); // A02_LED_NUM_PAD_A1
  primaryObjects.erase(primaryObjects.begin() + 3); // A02_LED_NUM_PAD_A2
  primaryObjects.erase(primaryObjects.begin() + 3); // A02_LED_NUM_PAD_A3
  primaryObjects.erase(primaryObjects.begin() + 3); // A02_LED_NUM_PAD_B1
  primaryObjects.erase(primaryObjects.begin() + 3); // A02_LED_NUM_PAD_B2
  primaryObjects.erase(primaryObjects.begin() + 3); // A02_LED_NUM_PAD_B3
  primaryObjects.erase(primaryObjects.begin() + 3); // A02_LED_NUM_PAD_C1
  primaryObjects.erase(primaryObjects.begin() + 3); // A02_LED_NUM_PAD_C2
  primaryObjects.erase(primaryObjects.begin() + 3); // A02_LED_NUM_PAD_C3
  primaryObjects.erase(primaryObjects.begin() + 3); // A03_TOGGLE
  // A03_LED
  primaryObjects.erase(primaryObjects.begin() + 4); // A04_TOGGLE
  // A04_LED_TOP
  // A04_LED_BOTTOM
  primaryObjects.erase(primaryObjects.begin() + 6); // A05_TOGGLE
  // A05_LED
}

map <QString, vector <RealObjectInformations> > Utility::loadSecondaryObjects(QString csvPath, QString secondaryPath) {
  map <QString, vector <RealObjectInformations> > informations = loadCSV(csvPath);

  for (map<QString, vector <RealObjectInformations> >::iterator it=informations.begin(); it!=informations.end(); ++it) {
    selectSecondaryObjects(it->second);

    QStringList pathList = it->first.split('/');
    QString initialFilename = pathList[1];
    pathList[1] = "LeftImage_"  + initialFilename;
    QString secondaryLeftPath = secondaryPath + pathList.join('_') + ".txt";
    pathList[1] = "RightImage_"  + initialFilename;
    QString secondaryRightPath = secondaryPath + pathList.join('_') + ".txt";

    vector <SecondaryObjectPosition> leftSecondaryPositions = Utility::loadSecondaryObjectsFile(secondaryLeftPath + "");
    vector <SecondaryObjectPosition> rightSecondaryPositions = Utility::loadSecondaryObjectsFile(secondaryRightPath + "");

    for (size_t i = 0; i < leftSecondaryPositions.size(); i++) {
      RealObjectInformations info;
      info.leftX = leftSecondaryPositions.at(i).x;
      info.leftY = leftSecondaryPositions.at(i).y;
      info.rightX = -1;
      info.rightY = -1;
      info.state = leftSecondaryPositions.at(i).imageType == 0 ? "SECONDARY" : "TERTIARY";
      it->second.push_back(info);
    }

    for (size_t i = 0; i < rightSecondaryPositions.size(); i++) {
      RealObjectInformations info;
      info.rightX = rightSecondaryPositions.at(i).x;
      info.rightY = rightSecondaryPositions.at(i).y;
      info.leftX = -1;
      info.leftY = -1;
      info.state = rightSecondaryPositions.at(i).imageType == 0 ? "SECONDARY" : "TERTIARY";
      it->second.push_back(info);
    }
  }

  return informations;
}


void Utility::saveSecondaryObjectsFile(QString secondaryObjectsPath, vector <SecondaryObjectPosition> secondaryPositions) {
  QFile file(secondaryObjectsPath);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
      return;
  }

  QTextStream out(&file);
  out << secondaryPositions.size() << endl;

  for (vector<SecondaryObjectPosition>::iterator it = secondaryPositions.begin() ; it != secondaryPositions.end(); ++it) {
    out << (*it).x << endl;
    out << (*it).y << endl;
    out << (*it).imageType << endl;
  }
}

vector <SecondaryObjectPosition> Utility::loadSecondaryObjectsFile(QString secondaryObjectsPath) {
  QFile file(secondaryObjectsPath);
  vector <SecondaryObjectPosition> secondaryPositions;
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
      return secondaryPositions;
  }



  QByteArray lineNb = file.readLine();
  lineNb.truncate(lineNb.size() - 1);
  size_t nb = lineNb.toInt();
  for (size_t i = 0; i < nb; i++) {
      QByteArray lineX = file.readLine();
      QByteArray lineY = file.readLine();
      QByteArray lineType = file.readLine();
      lineX.truncate(lineX.size() - 1);
      lineY.truncate(lineY.size() - 1);
      lineType.truncate(lineType.size() - 1);

      SecondaryObjectPosition position;
      position.x = lineX.toInt();
      position.y = lineY.toInt();
      position.imageType = lineType.toInt();

      secondaryPositions.push_back(position);
  }

  return secondaryPositions;
}


void Utility::drawInformations(QPainter & p, vector <RealObjectInformations> infos, bool isLeft, float scale) {
  for (size_t i = 0; i < infos.size(); i++) {
    int x;
    int y;
    if (isLeft) {
      x = infos.at(i).leftX;
      y = infos.at(i).leftY;
    } else {
      x = infos.at(i).rightX;
      y = infos.at(i).rightY;
    }
    if (x < 0 && y < 0) {
      continue;
    }
    x *= scale;
    y *= scale;
    p.drawLine(x - 3, y, x + 3, y);
    p.drawLine(x, y - 3, x, y + 3);
  }
}

void Utility::drawSecondaryObjects(QPainter & p, vector <SecondaryObjectPosition> secondaryPositions, float scale) {
  for (size_t i = 0; i < secondaryPositions.size(); i++) {
    if (secondaryPositions.at(i).imageType == 0) {
      p.setPen(QColor(255, 0, 0));
    } else {
      p.setPen(QColor(255, 255, 0));
    }
    int x = secondaryPositions.at(i).x;
    int y = secondaryPositions.at(i).y;
    x *= scale;
    y *= scale;
    p.drawLine(x - 3, y, x + 3, y);
    p.drawLine(x, y - 3, x, y + 3);
  }
}

void Utility::scaleInformations(vector <RealObjectInformations> & infos, float scale) {
  for (vector<RealObjectInformations>::iterator it = infos.begin() ; it != infos.end(); ++it) {
    (*it).leftX *= scale;
    (*it).leftY *= scale;
    (*it).rightX *= scale;
    (*it).rightY *= scale;
  }
}

void Utility::drawSquares(QPainter & p, vector <ObjectInformations> squares, float scale) {
  for (size_t i = 0; i < squares.size(); i++) {
    p.setPen(squares.at(i).state == "0" ? QColor(255, 0, 0) : QColor(0, 255, 0));
    p.drawRect(squares.at(i).x * scale, squares.at(i).y * scale, squares.at(i).width * scale, squares.at(i).height * scale);
  }
}

vector <ObjectInformations> Utility::loadObjectsInformations(QString path) {
  QFile file(path);
  vector <ObjectInformations> infos;
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
      return infos;
  }



  QByteArray lineNb = file.readLine();
  lineNb.truncate(lineNb.size() - 1);
  size_t nb = lineNb.toInt();
  for (size_t i = 0; i < nb; i++) {
    QByteArray lineState = file.readLine();
    QByteArray lineX = file.readLine();
    QByteArray lineY = file.readLine();
    QByteArray lineWidth = file.readLine();
    QByteArray lineHeight = file.readLine();
    lineState.truncate(lineState.size() - 1);
    lineX.truncate(lineX.size() - 1);
    lineY.truncate(lineY.size() - 1);
    lineWidth.truncate(lineWidth.size() - 1);
    lineHeight.truncate(lineHeight.size() - 1);

    ObjectInformations info;
    info.state = QString(lineState).toStdString();
    info.x = lineX.toInt();
    info.y = lineY.toInt();
    info.width = lineWidth.toInt();
    info.height = lineHeight.toInt();

    infos.push_back(info);
  }

  return infos;
}

void Utility::saveObjectsInformations(QString path, vector<ObjectInformations> & infos) {
  QFile file(path);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
      return;
  }

  QTextStream out(&file);
  out << infos.size() << endl;

  for (vector<ObjectInformations>::iterator it = infos.begin() ; it != infos.end(); ++it) {
    out << QString((*it).state.c_str()) << endl;
    out << (*it).x << endl;
    out << (*it).y << endl;
    out << (*it).width << endl;
    out << (*it).height << endl;
  }
}

void Utility::exportObjectsInformationsToImages(QString pathObjectsInformations, QString pathTo) {
  QDir objectsInformationsDirectory(pathObjectsInformations);
  QStringList list = objectsInformationsDirectory.entryList();
  map <QString, vector <RealObjectInformations> > primaryInfos = Utility::loadCSV("data/model.csv");
  int nbLeds = 0;
  int nbScrews = 0;
  int nbNumpad = 0;
  int nbNegatives = 0;
  for (int i = 0; i < list.size(); i++) {
    qDebug() << i << "/" << list.size();
    QString filename = list.at(i);
    if (filename.at(0) == '.') {
      continue;
    }
    QString pathObjectsInformationsItem = pathObjectsInformations + "/" + filename;
    QString pathSecondaryPositionsItem = "data/secondary_positions/" + filename;
    QStringList filenameSep = filename.split("_");
    QString subPath = filenameSep.at(0);
    filenameSep.removeAt(0);
    QStringList filenameSep2 = filenameSep.join("_").split(".");
    filenameSep2.removeAt(filenameSep2.size() - 1);
    QString pathImage = "data/" + subPath + "/" + filenameSep2.join(".");
    bool isLeft = filenameSep.at(0) == "LeftImage";
    QString pathPrimary = subPath + "/" + filenameSep.at(1).split(".").at(0) + "." + filenameSep.at(1).split(".").at(1);
    vector <RealObjectInformations> primaryObjects = primaryInfos.find(pathPrimary)->second;

    vector <int> imageVector = Utility::loadImage(pathImage);
    Image <float> * img = RobonautEye::convertVectorToImage(imageVector);

    vector <ObjectInformations> infos = Utility::loadObjectsInformations(pathObjectsInformationsItem);
    vector <SecondaryObjectPosition> secondaryObjectPositions = Utility::loadSecondaryObjectsFile(pathSecondaryPositionsItem);

    for (size_t j = 0; j < infos.size(); j++) {
      ObjectInformations infosItem = infos.at(j);
      float centerX = infosItem.x + infosItem.width / 2;
      float centerY = infosItem.y + infosItem.height / 2;
      bool positive = infosItem.state == "1";
      Image <float> * extract = Image<float>::reshape(img, infosItem.x, infosItem.y, infosItem.height, infosItem.width, 20, 20);
      //
      QString extractPath;
      if (positive) {
        QString subpath = "/leds/";
        bool isScrew = true;
        for (size_t k = 0; k < primaryObjects.size(); k++) {
          RealObjectInformations objectPosition = primaryObjects[k];
          float distance;
          if (isLeft) {
            distance = sqrt((objectPosition.leftX - centerX) * (objectPosition.leftX - centerX) + (objectPosition.leftY - centerY) * (objectPosition.leftY - centerY));
          } else {
            distance = sqrt((objectPosition.rightX - centerX) * (objectPosition.rightX - centerX) + (objectPosition.rightY - centerY) * (objectPosition.rightY - centerY));
          }

          if (distance < 15) {
            extractPath = pathTo + "/leds/" + QString().setNum(nbLeds) + ".png";
            nbLeds++;
            isScrew = false;
          }
        }

        for (size_t k = 0; k < secondaryObjectPositions.size(); k++) {
          SecondaryObjectPosition objectPosition = secondaryObjectPositions[k];
          float distance = sqrt((objectPosition.x - centerX) * (objectPosition.x - centerX) + (objectPosition.y - centerY) * (objectPosition.y - centerY));
          if (distance < 15) {
            if (objectPosition.imageType == 0) {
              extractPath = pathTo + "/screws/" + QString().setNum(nbScrews) + ".png";
              nbScrews++;
            } else {
              extractPath = pathTo + "/numpad/" + QString().setNum(nbNumpad) + ".png";
              nbNumpad++;
            }
            isScrew = false;
            break;
          }
        }

        if (isScrew) {
          extractPath = pathTo + "/screws/" + QString().setNum(nbScrews) + ".png";
          nbScrews++;
        }
      } else {
        extractPath = pathTo + "/negatives/" + QString().setNum(nbNegatives) + ".png";
        nbNegatives++;
      }
      ImageUtility::imageFloatToQImage(extract).save(extractPath);
      delete extract;
    }

    delete img;
  }
}


void Utility::drawPositions(QPainter & p, vector<Point> & positions, float scale) {
  for (vector<Point>::iterator it = positions.begin() ; it != positions.end(); ++it) {
    int x = (*it).x;
    int y = (*it).y;
    x *= scale;
    y *= scale;
    p.drawLine(x - 3, y, x + 3, y);
    p.drawLine(x, y - 3, x, y + 3);
  }
}

void Utility::drawPositions(QPainter & p, Point * positions, size_t nb, float scale) {
  for (size_t i = 0; i < nb; i++) {
    int x = positions[i].x;
    int y = positions[i].y;
    x *= scale;
    y *= scale;
    p.drawLine(x - 3, y, x + 3, y);
    p.drawLine(x, y - 3, x, y + 3);
  }
}

void Utility::drawPositionsWithIndex(QPainter & p, vector<Point> & positions, float scale) {
  drawPositions(p, positions, scale);
  int i = 0;
  for (vector<Point>::iterator it = positions.begin() ; it != positions.end(); ++it) {
    int x = (*it).x - 10;
    int y = (*it).y + 40;
    x *= scale;
    y *= scale;
    p.drawText(x, y, QString().setNum(i));
    i++;
  }
}

void Utility::drawPolygon(QPainter & p, Point * positions, size_t nb, float scale) {
  for (size_t j = 0; j < nb; j++) {
    p.drawLine(positions[j].x * scale, positions[j].y * scale, (positions[j + 1 == nb ? 0 : j + 1].x) * scale, (positions[j + 1 == nb ? 0 : j + 1].y) * scale);
  }
}

void Utility::drawPrimaryObject(QPainter & p, PrimaryObject & primaryObject, float scale) {
  Point positions[4] = {
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
    ),
  };
  drawPolygon(p, &positions[0], 4, scale);
}

void Utility::exportPrimaryObjectsToImages(QString pathTo) {
  map <QString, vector <RealObjectInformations> > infos = Utility::loadCSV("data/model.csv");
  vector <QString> filesList;
  for (map<QString, vector <RealObjectInformations> >::iterator it=infos.begin(); it!=infos.end(); ++it) {
    filesList.push_back(it->first);
  }
  int imageNum[3] = {0, 0, 0};
  for (size_t i = 0; i < filesList.size() * 2; i++) {
    if (i == 169 || i == 190 || i == 191 || i == 203 || i == 246) { // blacklisted images...
      continue;
    }
    qDebug() << (i + 1) << "/" << filesList.size() * 2;
    QStringList pathList = filesList.at(i / 2).split('/');
    pathList[1] = ((i % 2 == 0) ? "LeftImage_" : "RightImage_") + pathList[1];
    QString imagePath = "data/" + pathList.join('/');
    QPixmap imFrom(imagePath);
    QPainter p(&imFrom);
    p.setPen(QColor(255, 255, 255));
    vector<int> imageVector = loadImage(imagePath);
    Image<float> * image = RobonautEye::convertVectorToImage(imageVector);
    vector<RealObjectInformations> infosItem = infos.find(filesList.at(i / 2))->second;
    QString bigSwitchState = infosItem.at(1).state;
    infosItem.erase(infosItem.begin() + 1);

    qDebug() << "A";
    vector<PrimaryObject> primaryObjectsFrom;
    vector<Point> pointsFrom;
    vector<Point> pointsTo;
    vector<RealObjectInformations> filteredRealObjectInformations;
    PrimaryObject * realPrimaryObjects = RobonautEye().processImage(image);
    vector<PrimaryObject> realFilteredPrimaryObjects;
    for (size_t j = 0; j < infosItem.size(); j++) {
      RealObjectInformations infosElem = infosItem.at(j);
      Point p;
      if (i % 2 == 0) {
        p.x = infosElem.leftX;
        p.y = infosElem.leftY;
      } else {
        p.x = infosElem.rightX;
        p.y = infosElem.rightY;
      }
      if (infosElem.state != "HIDDEN" && p.x >= 0 && p.y >= 0) {
        pointsFrom.push_back(TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].position);
        pointsTo.push_back(p);
        primaryObjectsFrom.push_back(TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j]);
        filteredRealObjectInformations.push_back(infosElem);
        realFilteredPrimaryObjects.push_back(realPrimaryObjects[j]);
      }
    }
    delete realPrimaryObjects;

    qDebug() << "B";
    Array2D<double> * H = Geometry::getPerspectiveProjectionHomography(&pointsFrom[0], &pointsTo[0], pointsFrom.size());

    for (size_t j = 0; j < primaryObjectsFrom.size(); j++) {
      PrimaryObject primaryObject = RobonautEye::getPerspectiveProjection(H, primaryObjectsFrom[j]);
      //qDebug() << primaryObjectsFrom[j].position.x << primaryObjectsFrom[j].position.y << primaryObjectsFrom[j].size;

      QString state = filteredRealObjectInformations.at(j).state;
      int folderId;
      if (state == "DOWN" || state == "OFF") {
        folderId = 0;
      }
      if (state == "UP" || state == "ON") {
        folderId = 1;
      }
      if (state == "CENTER") {
        folderId = 2;
      }

      Image<float> * extract = RobonautEye::getPrimaryObjectImage(image, primaryObject);
      QImage imExtract = ImageUtility::imageFloatToQImage(extract);
      imExtract.save(pathTo + "/" + QString().setNum(folderId) + "/" + QString().setNum(imageNum[folderId]) +"_s.png");
      imageNum[folderId]++;
      PrimaryObject realPrimaryObject = realFilteredPrimaryObjects.at(j);
      if (realPrimaryObject.isCorrect && sqrt((realPrimaryObject.position.x - primaryObject.position.x) * (realPrimaryObject.position.x - primaryObject.position.x) + (realPrimaryObject.position.y - primaryObject.position.y) * (realPrimaryObject.position.y - primaryObject.position.y)) < 20) {
        Image<float> * realExtract = RobonautEye::getPrimaryObjectImage(image, realPrimaryObject);
        QImage realImExtract = ImageUtility::imageFloatToQImage(realExtract);
        realImExtract.save(pathTo + "/" + QString().setNum(folderId) + "/" + QString().setNum(imageNum[folderId]) +"_r.png");
        imageNum[folderId]++;
        delete realExtract;
      }
      delete extract;
      drawPrimaryObject(p, primaryObject, 1);
    }
    qDebug() << "C";

    if (bigSwitchState == "DOWN") {
      PrimaryObject primaryObject = RobonautEye::getPerspectiveProjection(H, TASKBOARD_PRIMARY_OBJECTS_POSITIONS[0]);
      Image<float> * extract = RobonautEye::getPrimaryObjectImage(image, primaryObject);
      QImage imExtract = ImageUtility::imageFloatToQImage(extract);
      imExtract.save(pathTo + "/0/" + QString().setNum(imageNum[0]) +".png");
      imageNum[0]++;

      drawPrimaryObject(p, primaryObject, 1);

      delete extract;
    }
    qDebug() << "D";

    imFrom.save(pathTo + "/current.png");
    delete H;
    delete image;

  }

}

vector<Point> Utility::assignSecondaryInformations(vector <RealObjectInformations> & imageSecondaryInfos, vector<Point> & projections, bool isLeft) {
  vector<Point> assignments(projections.size(), Point(-1, -1));
  // /!\ When 2 are assigned to the same position, delete the 2... We just need to collect some real data, not all of it :).
  vector<Point> imageSecondaryPoints = getImageSecondaryPoint(imageSecondaryInfos, isLeft);

  vector<int> assignmentsIndexes(projections.size(), -1);
  vector<int> nbAssignments(imageSecondaryPoints.size(), 0);

  for (size_t i = 0; i < projections.size(); i++) {
    Point projection = projections.at(i);
    int chosen = -1;
    float bestDistance = -1;
    for (size_t j = 0; j < imageSecondaryPoints.size(); j++) {
      Point secondaryPoint = imageSecondaryPoints.at(j);
      float distance = sqrt((projection.x - secondaryPoint.x) * (projection.x - secondaryPoint.x) + (projection.y - secondaryPoint.y) * (projection.y - secondaryPoint.y));
      if (distance < 100 && (chosen == -1 || distance < bestDistance)) {
        chosen = j;
        bestDistance = distance;
      }
    }

    if (chosen > -1) {
      assignmentsIndexes[i] = chosen;
      nbAssignments[chosen]++;
    }
  }

  for (size_t i = 0; i < projections.size(); i++) {
    int chosen = assignmentsIndexes[i];
    if (chosen != -1 && nbAssignments[chosen] < 2) {
      assignments[i] = imageSecondaryPoints[chosen];
    }
  }

  return assignments;
}

vector<Point> Utility::getImageSecondaryPoint(vector <RealObjectInformations> & imageSecondaryInfos, bool isLeft) {
  vector<Point> imageSecondaryPoints;
  for (vector<RealObjectInformations>::iterator it = imageSecondaryInfos.begin() ; it != imageSecondaryInfos.end(); ++it) {
    Point p;
    if (isLeft) {
      p.x = (*it).leftX;
      p.y = (*it).leftY;
    } else {
      p.x = (*it).rightX;
      p.y = (*it).rightY;
    }

    imageSecondaryPoints.push_back(p);
  }

  return imageSecondaryPoints;
}

vector< vector<float> > * Utility::getOrthographicHeuristicStatisticsItem(vector<Point> & assignments) {
  vector< vector<float> > * stats = new vector< vector<float> >();
  stats->resize(16);
  size_t assignmentsSize = assignments.size();
  for (size_t i = 0; i < assignmentsSize; i++) {
    Point fromA = TASKBOARD_SECONDARY_OBJECTS_POSITIONS[i].position;
    Point toA = assignments.at(i);
    if (fromA.x < 0 || fromA.y < 0 || toA.x < 0 || toA.y < 0) {
      continue;
    }
    for (size_t j = i + 1; j < assignmentsSize; j++) {
      Point fromB = TASKBOARD_SECONDARY_OBJECTS_POSITIONS[j].position;
      Point toB = assignments.at(j);
      if (fromB.x < 0 || fromB.y < 0 || toB.x < 0 || toB.y < 0) {
        continue;
      }
      for (size_t k = j + 1; k < assignmentsSize; k++) {
        Point fromC = TASKBOARD_SECONDARY_OBJECTS_POSITIONS[k].position;
        Point toC = assignments.at(k);
        if (fromC.x < 0 || fromC.y < 0 || toC.x < 0 || toC.y < 0) {
          continue;
        }
        float * infos = RobonautEye::getOrthographicProjectionInformations(fromA, fromB, fromC, toA, toB, toC);
        if (infos[8] < 0) {
          qDebug() << "incorrect equation!";
          continue;
        }

        for (size_t t = 0; t < 8; t++) {
          (*stats)[t].push_back(infos[t]);
        }

        delete infos;

        for (size_t l = k + 1; l < assignmentsSize; l++) {
          Point fromD = TASKBOARD_SECONDARY_OBJECTS_POSITIONS[l].position;
          Point toD = assignments.at(l);
          if (fromD.x < 0 || fromD.y < 0 || toD.x < 0 || toD.y < 0) {
            continue;
          }

          float * infos4 = RobonautEye::getOrthographicProjectionInformations(fromA, fromB, fromC, fromD, toA, toB, toC, toD);
          if (infos4[8] < 0) {
            qDebug() << "incorrect equation 4!";
            continue;
          }
          for (size_t t = 0; t < 8; t++) {
            (*stats)[t + 8].push_back(infos4[t]);
          }

          delete infos4;
        }
      }
    }
  }

  return stats;
}

void Utility::generateOrthographicHeuristicStatistics() {
  QFile file("data/orthographic_stats.csv");
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
      return;
  }

  QTextStream out(&file);

  map <QString, vector <RealObjectInformations> > primaryInfos = Utility::loadCSV("data/model.csv");
  map <QString, vector <RealObjectInformations> > secondaryInfos = Utility::loadSecondaryObjects("data/model.csv", "data/secondary_positions/");
  vector< vector<float> > overallStats(16);
  vector <QString> filesList;
  for (map<QString, vector <RealObjectInformations> >::iterator it=primaryInfos.begin(); it!=primaryInfos.end(); ++it) {
    filesList.push_back(it->first);
  }
  for (size_t i = 0; i < filesList.size() * 2; i++) {
    qDebug() << (i + 1) << "/" << filesList.size() * 2;
    if (i == 169 || i == 190 || i == 191 || i == 203 || i == 246) { // blacklisted images...
      continue;
    }
    vector <RealObjectInformations> imagePrimaryInfos = primaryInfos.find(filesList.at(i / 2))->second;
    vector <RealObjectInformations> imageSecondaryInfos = secondaryInfos.find(filesList.at(i / 2))->second;

    vector <Point> primaryPositions = getPrimaryObjectsPositions(imagePrimaryInfos, i % 2 == 0);
    primaryPositions.erase(primaryPositions.begin() + 1);

    vector <Point> primaryFilteredPositionsFrom;
    vector <Point> primaryFilteredPositionsTo;
    for (size_t j = 0; j < TASKBOARD_PRIMARY_OBJECTS_NB_POSITIONS; j++) {
      Point pointTo = primaryPositions.at(j);
      if (pointTo.x < 0 || pointTo.y < 0) {
        continue;
      }

      primaryFilteredPositionsFrom.push_back(TASKBOARD_PRIMARY_OBJECTS_POSITIONS[j].position);
      primaryFilteredPositionsTo.push_back(primaryPositions.at(j));
    }

    if (primaryFilteredPositionsFrom.size() > 3) {

      Array2D<double> * H = Geometry::getPerspectiveProjectionHomography(
        &primaryFilteredPositionsFrom[0],
        &primaryFilteredPositionsTo[0],
        primaryFilteredPositionsFrom.size()
      );

      Point * projectionsRaw = Geometry::getPerspectiveProjectionPoints(H, RobonautEye().getSecondaryPoints(), 20);

      vector < Point > projections;
      for (size_t i = 0; i < TASKBOARD_SECONDARY_OBJECTS_NB_POSITIONS; i++) {
        projections.push_back(projectionsRaw[i]);
      }

      vector <Point> assignments = Utility::assignSecondaryInformations(imageSecondaryInfos, projections, i % 2 == 0);
      vector< vector<float> > * stats = getOrthographicHeuristicStatisticsItem(assignments);

      for (size_t j = 0; j < 16; j++) {
        overallStats[j].insert(overallStats[j].end(), (*stats)[j].begin(), (*stats)[j].end());
      }

      delete stats;
    }


    //assignSecondaryInformations()
  }

  float statsPrecisions[16] = {
    0.25,
    0.25,
    0.25,
    0.25,
    50,
    0.25,
    0.25,
    50,
    0.1,
    0.1,
    0.1,
    0.1,
    50,
    0.1,
    0.1,
    50
  };

  QStringList labels;
  labels.push_back("Distance");
  labels.push_back("Angle");
  labels.push_back("Ax");
  labels.push_back("Bx");
  labels.push_back("Cx");
  labels.push_back("Ay");
  labels.push_back("By");
  labels.push_back("Cy");
  labels.push_back("Variance distance");
  labels.push_back("Variance angle");
  labels.push_back("Variance Ax");
  labels.push_back("Variance Bx");
  labels.push_back("Variance Cx");
  labels.push_back("Variance Ay");
  labels.push_back("Variance By");
  labels.push_back("Variance Cy");

  vector < map <int, int> > statsSummary(16);

  for (size_t i = 0; i < 16; i++) {
    for (size_t j = 0; j < overallStats[i].size(); j++) {
      int val = overallStats[i][j] / statsPrecisions[i];
      if (statsSummary[i].find(val) == statsSummary[i].end()) {
        statsSummary[i][val] = 1;
      } else {
        statsSummary[i][val]++;
      }
    }
  }

  for (size_t i = 0; i < 16; i++) {
    out << endl << endl << endl << endl << labels[i] << endl << endl;
    for (map<int,int>::iterator it=statsSummary[i].begin(); it!=statsSummary[i].end(); ++it) {
      out << QString().setNum(it->first * statsPrecisions[i]) + " - " + QString().setNum((it->first + 1) * statsPrecisions[i]) << "," << it->second << endl;
    }
  }

  file.close();
}

vector <Point> Utility::getPrimaryObjectsPositions(vector <RealObjectInformations> & infos, bool isLeft) {
  vector <Point> positions;

  for (vector<RealObjectInformations>::iterator itItem = infos.begin() ; itItem != infos.end(); ++itItem) {
    RealObjectInformations infos = (*itItem);
    Point position;
    if (isLeft) {
      position.x = infos.leftX;
      position.y = infos.leftY;
    } else {
      position.x = infos.rightX;
      position.y = infos.rightY;
    }
    positions.push_back(position);
  }

  return positions;
}


void Utility::upgradeSecondaryPositions() {
  map <QString, vector <RealObjectInformations> > primaryInfos = Utility::loadCSV("data/model.csv");
  vector <QString> filesList;
  for (map<QString, vector <RealObjectInformations> >::iterator it=primaryInfos.begin(); it!=primaryInfos.end(); ++it) {
    filesList.push_back(it->first);
  }

  for (size_t i = 0; i < filesList.size() * 2; i++) {
    map<QString, vector <RealObjectInformations> >::iterator it = primaryInfos.find(filesList.at(i / 2)); //filesList.at(i / 2)
    vector <RealObjectInformations> primaryObjects = it->second;
    selectSecondaryObjects(primaryObjects);
    QStringList pathList = it->first.split('/');
    pathList[1] = ((i % 2 == 0) ? "LeftImage_" : "RightImage_") + pathList[1];
    QString secondaryObjectsPath = "data/secondary_positions/" + pathList.join('_') + ".txt";
    vector <SecondaryObjectPosition> secondaryObjects = loadSecondaryObjectsFile(secondaryObjectsPath);
    for (size_t j = 0; j < secondaryObjects.size(); j++) {
      SecondaryObjectPosition secondaryObject = secondaryObjects[j];
      secondaryObject.imageType = 0;
      secondaryObjects[j] = secondaryObject;
    }
    saveSecondaryObjectsFile(secondaryObjectsPath, secondaryObjects);

  }
}


void Utility::getScoreStatistics() {
  QTime time = QTime::currentTime();
  QDate date = QDate::currentDate();

  QFile file("data/reports/" + QString().setNum(time.hour()) + "-" + QString().setNum(time.minute()) + "_" + QString().setNum(date.day()) + "-" + QString().setNum(date.month()) + ".txt");
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
      return;
  }

  QTextStream out(&file);


  map <QString, vector <RealObjectInformations> > primaryInfos = Utility::loadCSV("data/model.csv");
  vector <QString> filesList;
  for (map<QString, vector <RealObjectInformations> >::iterator it=primaryInfos.begin(); it!=primaryInfos.end(); ++it) {
    filesList.push_back(it->first);
  }


  DetailedScore base;
  base.hiddenScore = 0;
  base.overallScore = 0;
  base.proximityScore = 0;
  base.stateScore = 0;
  base.nbFiles = 0;

  map <QString, DetailedScore> stats;
  stats["ISS"] = base;
  stats["Lab"] = base;
  stats["Lab2"] = base;
  stats["Lab3"] = base;
  stats["Sim"] = base;

  map <QString, float> environmentScores;
  environmentScores["ISS"] = 1.1;
  environmentScores["Lab"] = 1;
  environmentScores["Lab2"] = 0.8;
  environmentScores["Lab3"] = 0.7;
  environmentScores["Sim"] = 0.5;



  int totalNbFiles = 0;
  float totalCoef = 0;
  for (size_t i = 0; i < filesList.size(); i++) {
    qDebug() << (i + 1) << "/" << filesList.size();
    QElapsedTimer timer;
    QString folder = filesList.at(i).split("/").at(0);
    timer.start();
    DetailedScore score = evaluateScoreFile(filesList.at(i));
    int totalTime = timer.elapsed();
    if (score.overallScore < 58) {
      out << "/!\\ ";
    }
    out << i << ", " << filesList.at(i) << ": " << score.overallScore << " (" << totalTime << "ms)" << endl;
    stats[folder].hiddenScore += score.hiddenScore;
    stats[folder].overallScore += score.overallScore;
    stats[folder].proximityScore += score.proximityScore;
    stats[folder].stateScore += score.stateScore;
    stats[folder].nbFiles++;
    totalNbFiles++;
    totalCoef += environmentScores[folder];
  }

  float maxProximityScore = 2.0 / 3;
  float maxStateScore = 1.0 / 3;


  float totalOverallScore = 0;
  float totalHiddenScore = 0;
  float totalStateScore = 0;
  float totalProximityScore = 0;
  float totalProximityScoreEffect = 0;
  float totalStateScoreEffect = 0;

  for (map<QString,float>::iterator it=environmentScores.begin(); it!=environmentScores.end(); ++it) {
    QString folder = it->first;
    DetailedScore score = stats[folder];
    out << endl << endl << "Folder: " << folder << endl << endl;
    out << "Nb files: " << score.nbFiles << endl;
    out << "Environment score: " << it->second << endl;
    if (score.nbFiles > 0) {
      float hiddenScore = score.hiddenScore / (score.nbFiles * 66.0);
      float proximityScore = score.proximityScore / (score.nbFiles * 66.0);
      float stateScore = score.stateScore / (score.nbFiles * 66.0);
      float overallScore = score.overallScore / (score.nbFiles * 66.0);
      float possibleProximityScoreImprovement = maxProximityScore - (proximityScore + (2.0 * hiddenScore / 3));
      float possibleStateScoreImprovement = maxStateScore - (stateScore + (1.0 * hiddenScore / 3));
      float possibleProximityScoreEffect = possibleProximityScoreImprovement * score.nbFiles * environmentScores[folder] / totalCoef;
      float possibleStateScoreEffect = possibleStateScoreImprovement * score.nbFiles * environmentScores[folder] / totalCoef;
      totalOverallScore += overallScore * score.nbFiles * environmentScores[folder] / totalCoef;
      totalHiddenScore += hiddenScore * score.nbFiles * environmentScores[folder] / totalCoef;
      totalStateScore += stateScore * score.nbFiles * environmentScores[folder] / totalCoef;
      totalProximityScore += proximityScore * score.nbFiles * environmentScores[folder] / totalCoef;
      totalProximityScoreEffect += possibleProximityScoreEffect;
      totalStateScoreEffect += possibleStateScoreEffect;
      out << "Hidden score: " << hiddenScore << endl;
      out << "Proximity score: " << proximityScore << "(" << possibleProximityScoreImprovement << ", " << possibleProximityScoreEffect << ")" << endl;
      out << "State score: " << stateScore << "(" << possibleStateScoreImprovement << ", " << possibleStateScoreEffect << ")" << endl;
      out << "Overall score: " << overallScore << "(" << (possibleStateScoreEffect + possibleProximityScoreEffect) << ")" << endl;
    }
  }

  out << endl << endl << "TOTAL" << endl << endl;
  out << "Nb files: " << totalNbFiles << endl;
  out << "Hidden score: " << totalHiddenScore << endl;
  out << "Proximity score: " << totalProximityScore << "(" << totalProximityScoreEffect << ")" << endl;
  out << "State score: " << totalStateScore << "(" << totalStateScoreEffect << ")" << endl;
  out << "Overall score: " << totalOverallScore << "(" << (totalProximityScoreEffect + totalStateScoreEffect) << ")" << endl;


  file.close();
}

DetailedScore Utility::evaluateScoreFile(QString file) {
  map <QString, vector <RealObjectInformations> > primaryInfosList = Utility::loadCSV("data/model.csv");
  map <QString, vector <RealObjectInformations> >::iterator it = primaryInfosList.find(file);
  QString folder = file.split("/").at(0);
  float threshold = 10;
  if (folder == "ISS" || folder == "Lab") {
    threshold = 15;
  }

  QStringList pathList = file.split('/');
  QString pathSecondItem = pathList[1];
  pathList[1] = "LeftImage_"+ pathSecondItem;
  QString leftImagePath = "data/" + pathList.join('/');
  pathList[1] = "RightImage_"+ pathSecondItem;
  QString rightImagePath = "data/" + pathList.join('/');

  vector<int> leftImageVector = Utility::loadImage(leftImagePath);
  Image<float> * leftImage = RobonautEye::convertVectorToImage(leftImageVector);
  vector<int> rightImageVector = Utility::loadImage(rightImagePath);
  Image<float> * rightImage = RobonautEye::convertVectorToImage(rightImageVector);

  RobonautEye eye;
  FinalObject * finalObjects = eye.getFinalObjects(leftImage, rightImage);
  DetailedScore score = evaluateScoreItem(&it->second[0], finalObjects, threshold);
  delete finalObjects;
  delete leftImage;
  delete rightImage;
  return score;
}

DetailedScore Utility::evaluateScoreItem(RealObjectInformations * primaryInfos, FinalObject * finalObjects, float threshold) {
  DetailedScore detailedScore;
  detailedScore.hiddenScore = 0;
  detailedScore.proximityScore = 0;
  detailedScore.stateScore = 0;
  for (size_t i = 0; i < TASKBOARD_NB_FINAL_OBJECTS; i++) {
    RealObjectInformations infosItem = primaryInfos[i];
    if (infosItem.state == "HIDDEN") {
      detailedScore.hiddenScore += 3;
    } else {
      if (infosItem.leftX != -1 || infosItem.leftY != -1) {
        float distance = sqrt((finalObjects[i].leftPosition.x - infosItem.leftX) * (finalObjects[i].leftPosition.x - infosItem.leftX) + (finalObjects[i].leftPosition.y - infosItem.leftY) * (finalObjects[i].leftPosition.y - infosItem.leftY));
        if (distance <= 2 * threshold) {
          if (distance <= threshold) {
            detailedScore.proximityScore += 1;
          } else {
            detailedScore.proximityScore += 1 - pow(((distance - threshold) * 1.0 / threshold), 1.5);
          }
        }
      } else {
        detailedScore.proximityScore += 1;
      }

      if (infosItem.rightX != -1 || infosItem.rightY != -1) {
        float distance = sqrt((finalObjects[i].rightPosition.x - infosItem.rightX) * (finalObjects[i].rightPosition.x - infosItem.rightX) + (finalObjects[i].rightPosition.y - infosItem.rightY) * (finalObjects[i].rightPosition.y - infosItem.rightY));
        if (distance <= 2 * threshold) {
          if (distance <= threshold) {
            detailedScore.proximityScore += 1;
          } else {
            detailedScore.proximityScore += 1 - pow(((distance - threshold) * 1.0 / threshold), 1.5);
          }
        }
      } else {
        detailedScore.proximityScore += 1;
      }

      if (((infosItem.state == "DOWN" || infosItem.state == "OFF") && finalObjects[i].state == 0) ||
          ((infosItem.state == "UP" || infosItem.state == "ON") && finalObjects[i].state == 1) ||
          ((infosItem.state == "CENTER") && finalObjects[i].state == 2)
          ) {
        detailedScore.stateScore += 1;
      }
    }
  }

  detailedScore.overallScore = detailedScore.hiddenScore + detailedScore.proximityScore + detailedScore.stateScore;

  return detailedScore;
}


void Utility::exportDecals(QString pathTo) {
  map <QString, vector <RealObjectInformations> > infos = Utility::loadCSV("data/model.csv");
  vector <QString> filesList;
  for (map<QString, vector <RealObjectInformations> >::iterator it=infos.begin(); it!=infos.end(); ++it) {
    filesList.push_back(it->first);
  }
  int num = 0;
  for (size_t i = 0; i < filesList.size() * 2; i++) {
    if (i == 169 || i == 190 || i == 191 || i == 203 || i == 246) { // blacklisted images...
      continue;
    }


    qDebug() << (i + 1) << "/" << filesList.size() * 2;
    QStringList pathList = filesList.at(i / 2).split('/');
    pathList[1] = ((i % 2 == 0) ? "LeftImage_" : "RightImage_") + pathList[1];
    QString imagePath = "data/" + pathList.join('/');
    QPixmap imFrom(imagePath);
    QPainter p(&imFrom);
    p.setPen(QColor(255, 255, 255));
    vector<int> imageVector = loadImage(imagePath);
    Image<float> * image = RobonautEye::convertVectorToImage(imageVector);
    vector<RealObjectInformations> infosItem = infos.find(filesList.at(i / 2))->second;
    //infosItem.erase(infosItem.begin() + 1); <--- here

    RobonautEye eye;
    PrimaryObject * realPrimaryObjects = new PrimaryObject[TASKBOARD_PRIMARY_OBJECTS_NB_POSITIONS];
    ObjectInformationsMapping * secondaryObjectsInformations = eye.getSecondaryObjectsInformationsFromImage(image);
    bool error;
    if (secondaryObjectsInformations->positives.size() >= 4) {
      Array2D<double> * H = eye.getHomography(image, secondaryObjectsInformations, error);

      // ROCKER SWITCH
      /* add infosItem.erase(infosItem.begin() + 1); at the begining if you want to uncomment this...
      RealObjectInformations infoSwitch = infosItem.at(2);
      Point positionSwitch = (i % 2 == 0) ? Point(infoSwitch.leftX, infoSwitch.leftY) : Point(infoSwitch.rightX, infoSwitch.rightY);
      if (infoSwitch.state != "HIDDEN" && positionSwitch.x >= 0 && positionSwitch.y >= 0) {
        PrimaryObject primaryObject = RobonautEye::getPerspectiveProjection(H, TASKBOARD_PRIMARY_OBJECTS_POSITIONS[2]);
        if (!(primaryObject.position.x < 0 || primaryObject.position.y < 0 || primaryObject.position.x > image->getWidth() || primaryObject.position.y > image->getHeight())) {
          if (sqrt((primaryObject.position.x - positionSwitch.x) * (primaryObject.position.x - positionSwitch.x) + (primaryObject.position.y - positionSwitch.y) * (primaryObject.position.y - positionSwitch.y)) < 150) {
            qDebug() << "here" << positionSwitch.x << positionSwitch.y << primaryObject.position.x << primaryObject.position.y;
            Image<float> * extract = eye.getPrimaryObjectImage(image, primaryObject);
            QImage extractIm = ImageUtility::imageFloatToQImage(extract);
            extractIm.save(pathTo + "/switch/" + QString().setNum(num) + ".png");
            QFile file(pathTo + "/switch/" + QString().setNum(num) + ".txt");
            if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
                return;
            }
            QTextStream out(&file);

            out << positionSwitch.x << endl;
            out << positionSwitch.y << endl;

            Point * pointsFrom = new Point[3];
            pointsFrom[0].x = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[2].position.x;
            pointsFrom[0].y = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[2].position.y + 1;
            pointsFrom[1].x = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[2].position.x + 1;
            pointsFrom[1].y = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[2].position.y;
            pointsFrom[2].x = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[2].position.x;
            pointsFrom[2].y = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[2].position.y + 1;

            Point * pointsTo = Geometry::getPerspectiveProjectionPoints(H, pointsFrom, 3);
            out << pointsTo[0].x << endl;
            out << pointsTo[0].y << endl;
            out << pointsTo[1].x << endl;
            out << pointsTo[1].y << endl;
            out << pointsTo[2].x << endl;
            out << pointsTo[2].y << endl;

            file.close();
            num++;
            delete extract;
            delete pointsFrom;
            delete pointsTo;
          }
        }
        */

      // POWER SWITCH
      RealObjectInformations infoSwitch = infosItem.at(1);
      Point positionSwitch = (i % 2 == 0) ? Point(infoSwitch.leftX, infoSwitch.leftY) : Point(infoSwitch.rightX, infoSwitch.rightY);
      if (infoSwitch.state != "HIDDEN" && positionSwitch.x >= 0 && positionSwitch.y >= 0) {
        PrimaryObject primaryObject = RobonautEye::getPerspectiveProjection(H, TASKBOARD_PRIMARY_OBJECTS_POSITIONS[0]);
        if (!(primaryObject.position.x < 0 || primaryObject.position.y < 0 || primaryObject.position.x > image->getWidth() || primaryObject.position.y > image->getHeight())) {
          if (sqrt((primaryObject.position.x - positionSwitch.x) * (primaryObject.position.x - positionSwitch.x) + (primaryObject.position.y - positionSwitch.y) * (primaryObject.position.y - positionSwitch.y)) < 2000) {
            qDebug() << "here" << positionSwitch.x << positionSwitch.y << primaryObject.position.x << primaryObject.position.y;
            Image<float> * extract = eye.getPrimaryObjectImage(image, primaryObject);
            QImage extractIm = ImageUtility::imageFloatToQImage(extract);

            extractIm.save(pathTo + "/power_switch/" + infoSwitch.state + "/" + QString().setNum(num) + ".png");
            QFile file(pathTo + "/power_switch/" + infoSwitch.state + "/" + QString().setNum(num) + ".txt");
            if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
                return;
            }
            QTextStream out(&file);

            out << positionSwitch.x << endl;
            out << positionSwitch.y << endl;

            Point * pointsFrom = new Point[3];
            pointsFrom[0].x = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[0].position.x;
            pointsFrom[0].y = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[0].position.y;
            pointsFrom[1].x = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[0].position.x + 1;
            pointsFrom[1].y = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[0].position.y;
            pointsFrom[2].x = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[0].position.x;
            pointsFrom[2].y = TASKBOARD_PRIMARY_OBJECTS_POSITIONS[0].position.y + 1;

            Point * pointsTo = Geometry::getPerspectiveProjectionPoints(H, pointsFrom, 3);
            out << pointsTo[0].x << endl;
            out << pointsTo[0].y << endl;
            out << pointsTo[1].x << endl;
            out << pointsTo[1].y << endl;
            out << pointsTo[2].x << endl;
            out << pointsTo[2].y << endl;

            file.close();
            num++;
            delete extract;
            delete pointsFrom;
            delete pointsTo;
          }
        }

        //realPrimaryObjects
      }
      delete H;
    }
    delete secondaryObjectsInformations;
    delete realPrimaryObjects;
    delete image;

  }

}

void Utility::getMostFrequentState(QString modelCsvPath) {
  map<QString, vector <RealObjectInformations> > models = loadCSV(modelCsvPath);
  size_t nbPerItem[TASKBOARD_NB_FINAL_OBJECTS][3];
  for (size_t i = 0; i < TASKBOARD_NB_FINAL_OBJECTS; i++) {
    nbPerItem[i][0] = 0;
    nbPerItem[i][1] = 0;
    nbPerItem[i][2] = 0;
  }

  for (map<QString, vector <RealObjectInformations> >::iterator it=models.begin(); it!=models.end(); ++it) {
    vector <RealObjectInformations> itemInfos = it->second;
    for (size_t i = 0; i < itemInfos.size(); i++) {
      RealObjectInformations realObject = itemInfos.at(i);
      int stateId = 0;
      if (realObject.state == "DOWN" || realObject.state == "OFF") {
        stateId = 0;
      }
      if (realObject.state == "UP" || realObject.state == "ON") {
        stateId = 1;
      }
      if (realObject.state == "CENTER") {
        stateId = 2;
      }

      nbPerItem[i][stateId]++;
    }
  }


  for (size_t i = 0; i < TASKBOARD_NB_FINAL_OBJECTS; i++) {
    size_t bestState = 0;
    //qDebug() << "----" << i;
    for (size_t j = 0; j < 3; j++) {
      if (nbPerItem[i][j] > nbPerItem[i][bestState]) {
        bestState = j;
      }
      //qDebug() << nbPerItem[i][j];
    }
    qDebug() << i << ":" << bestState;
  }
}
