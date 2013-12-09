import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import javax.imageio.*; 
import javax.swing.*;
import java.text.*;

/**
 * This tool loads model.csv and runs through every row of the training data, treating some
 * manually selected triplets of points (e.g. Panel LED, Rocker Switch, and middle of button
 * Grid) as basis vectors for all other components on the board. It then prints the relative
 * coordinates of every component as defined by the basis vectors. It averages all such
 * values across all images.
 */
class TestEntry {
  String dir;
  String filename;
  String[] model;
  public TestEntry(String dir, String filename) {
    this.dir = dir;
    this.filename = filename;
  }
}

public class Trainer {
  static int getLum(int rgb) {
    int r = ((rgb >> 16) & 0x00ff);
    int g = ((rgb >> 8) & 0x00ff);
    int b = ((rgb) & 0x00ff);
    return r + g + b;
  }

  static void loadEntries(Map<String, TestEntry> entries) {
    String[] folders = { "ISS", "Lab", "Lab2", "Lab3", "Sim" };
    for (String folder : folders) {
      File dir = new File(folder);
      Set<String> imgSuffixes = new TreeSet<String>();
      for (String filename : dir.list()) {
        String[] parts = filename.split("_");
        if (parts.length == 2) {
          imgSuffixes.add(parts[1]);
        }
      }

      for (String img : imgSuffixes) {
        entries.put(folder + "/" + img, new TestEntry(folder, img));
      }
    }
  }
  static void loadCsv(Map<String, TestEntry> entries) throws Exception {
    Scanner scan = new Scanner(new File("model.csv"));
    scan.nextLine();
    while (scan.hasNextLine()) {
      String[] line = scan.nextLine().split(",");
      TestEntry entry = entries.get(line[0] + "/" + line[1]);
      if (entry == null) {
        System.out.println("Missing from entries!: " + line[0] + "/" + line[1]);
        continue;
        //System.exit(1);
      }
      entry.model = line;
    }
  }

  // Use three known points to represent vector-space coordinates of the button grid.
  static void calculateBasisCoordinatesForGrid(Map<String, TestEntry> entries) throws Exception {
    RobonautEyeTester.debug = false;
all_entries:    for (TestEntry entry : entries.values()) {
      int[] arrLeft = RobonautEyeTester.imageToArray(
          entry.dir + "/LeftImage_" + entry.filename);
      int[] arrRight = RobonautEyeTester.imageToArray(
          entry.dir + "/RightImage_" + entry.filename);
      int width = arrLeft[1];
      int height = arrLeft[0];

      // Vector led -> rocker.
      double a1, b1;
      DecimalFormat fmt = new DecimalFormat("0.00");
      StringBuilder msg = new StringBuilder();
      {
        if (Integer.parseInt(entry.model[16 + 2]) < 0 || Integer.parseInt(entry.model[11 + 2]) < 0) {
          continue;
        }
        if (Integer.parseInt(entry.model[17 + 2]) < 0 || Integer.parseInt(entry.model[12 + 2]) < 0) {
          continue;
        }
        int dx = Integer.parseInt(entry.model[16 + 2]) - Integer.parseInt(entry.model[11 + 2]);
        int dy = Integer.parseInt(entry.model[17 + 2]) - Integer.parseInt(entry.model[12 + 2]);
        a1 = dx;
        b1 = dy;
      }
      double a2, b2;
      a2 = b1;
      b2 = -a1;
      double det = 1.0 / (a1 * b2 - b1 * a2);
      if (det == 0) System.out.println("IZEROOO!!!");
      int ledX = Integer.parseInt(entry.model[11 + 2]);
      int ledY = Integer.parseInt(entry.model[12 + 2]);
      for (int i = 30; i <= 70; i += 5) {
        if (Integer.parseInt(entry.model[i + 1 + 2]) < 0 || Integer.parseInt(entry.model[i + 2 + 2]) < 0) {
          continue all_entries;
        }
        double a3 = Integer.parseInt(entry.model[i + 1 + 2]) - ledX;
        double b3 = Integer.parseInt(entry.model[i + 2 + 2]) - ledY;
        double x = det * (a3 * b2 - b3 * a2);
        double y = det * (-a3 * b1 + b3 * a1);
        msg.append(fmt.format(x) + "," + fmt.format(y) + " ");
      }
      System.out.println(entry.dir + "/" + entry.filename + ": " + msg.toString());
    }
  }

  // Rectangular coordinates.
  static void calculateBasisVectors(Map<String, TestEntry> entries) throws Exception {
    RobonautEyeTester.debug = false;
    ArrayList<ArrayList<Double>> vectorsA = new ArrayList<ArrayList<Double>>();
    ArrayList<ArrayList<Double>> vectorsB = new ArrayList<ArrayList<Double>>();
    for (int i = 0; i < 22; ++i) {
      vectorsA.add(new ArrayList<Double>());
      vectorsB.add(new ArrayList<Double>());
    }
    for (TestEntry entry : entries.values()) {
      int[] arrLeft = RobonautEyeTester.imageToArray(
          entry.dir + "/LeftImage_" + entry.filename);
      int[] arrRight = RobonautEyeTester.imageToArray(
          entry.dir + "/RightImage_" + entry.filename);
      int width = arrLeft[1];
      int height = arrLeft[0];

      for (int lr = 0; lr <= 1; ++lr) {
        int ledX = Integer.parseInt(entry.model[11 + 2 + 2 * lr]);
        int ledY = Integer.parseInt(entry.model[12 + 2 + 2 * lr]);
        int rockerX = Integer.parseInt(entry.model[16 + 2 + 2 * lr]);
        int rockerY = Integer.parseInt(entry.model[17 + 2 + 2 * lr]);
        //int gridX = Integer.parseInt(entry.model[51 + 2 + 2 * lr]);
        //int gridY = Integer.parseInt(entry.model[52 + 2 + 2 * lr]);
        int gridX = Integer.parseInt(entry.model[46 + 2 + 2 * lr]);
        int gridY = Integer.parseInt(entry.model[47 + 2 + 2 * lr]);
        if (ledX <= 0 || ledY <= 0 ||
            rockerX <= 0 || rockerY <= 0 ||
            gridX <= 0 || gridY <= 0) continue;

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
          //a2 = b1;
          //b2 = -a1;
        }
        double det = 1.0 / (a1 * b2 - b1 * a2);

        DecimalFormat fmt = new DecimalFormat("0.00");
        StringBuilder msg = new StringBuilder();
        for (int comp = 0; comp < 22; ++comp) {
          int compX = Integer.parseInt(entry.model[comp * 5 + 1 + 2 + 2 * lr]);
          int compY = Integer.parseInt(entry.model[comp * 5 + 2 + 2 + 2 * lr]);
          if (compX <= 0 || compY <= 0) {
          if (comp >= 5 && comp <= 14 || comp == 16 || comp == 18 || comp == 19 || comp == 21)
            msg.append("N/A ");
            continue;
          }
          int a3 = compX - ledX;
          int b3 = compY - ledY;
          double x = det * (a3 * b2 - b3 * a2);
          double y = det * (-a3 * b1 + b3 * a1);
          if (Double.isNaN(x) || Double.isNaN(y)) {
            System.out.println("Got NaN for " + entry.dir + "/" + entry.filename);
            System.exit(1);
          }
          if (comp >= 5 && comp <= 14 || comp == 16 || comp == 18 || comp == 19 || comp == 21)
          msg.append(fmt.format(x) + "," + fmt.format(y) + " ");
          vectorsA.get(comp).add(x);
          vectorsB.get(comp).add(y);
        }
        System.out.println(entry.dir + "/" + entry.filename + ": " + msg.toString());
      }
    }  // For each entry.
    DecimalFormat fmt = new DecimalFormat("0.00");
    System.out.print("double compRatios[22][2] = {");
    for (int comp = 0; comp < 22; ++comp) {
      if (comp % 4 == 0) {
        System.out.println();
        System.out.print("    ");
      }
      double sum = 0;
      for (Double d : vectorsA.get(comp)) {
        sum += d;
      }
      sum /= vectorsA.get(comp).size();
      double ang = 0;
      for (Double d : vectorsB.get(comp)) {
        ang += d;
      }
      ang /= vectorsB.get(comp).size();
      System.out.print("{" + fmt.format(sum) + "," + fmt.format(ang) + "}");
      if (comp != 21) {
        System.out.print(",");
      }
    }
    System.out.println("};");
  }

  static void calculateConversions(Map<String, TestEntry> entries) throws Exception {
    double[][] compRatios = {
        {0.01,-0.18},{-0.14,-0.18},{0.00,0.00},{1.00,0.00},
        {0.72,-0.01},{1.22,-0.03},{0.74,-1.23},{0.74,-0.97},
        {0.73,-0.71},{0.98,-1.21},{0.97,-0.95},{0.96,-0.69},
        {1.20,-1.19},{1.19,-0.93},{1.19,-0.68},{1.76,-1.09},
        {1.56,-1.11},{1.74,-0.53},{1.54,-0.53},{1.92,-0.51},
        {1.73,0.03},{1.54,0.04}};
    RobonautEyeTester.debug = false;
    ArrayList<ArrayList<Double>> vectorsA = new ArrayList<ArrayList<Double>>();
    ArrayList<ArrayList<Double>> vectorsB = new ArrayList<ArrayList<Double>>();
    for (int i = 0; i < 22; ++i) {
      vectorsA.add(new ArrayList<Double>());
      vectorsB.add(new ArrayList<Double>());
    }
    for (TestEntry entry : entries.values()) {
      int[] arrLeft = RobonautEyeTester.imageToArray(
          entry.dir + "/LeftImage_" + entry.filename);
      int[] arrRight = RobonautEyeTester.imageToArray(
          entry.dir + "/RightImage_" + entry.filename);
      int width = arrLeft[1];
      int height = arrLeft[0];

      for (int lr = 0; lr <= 1; ++lr) {
        int ledX = Integer.parseInt(entry.model[11 + 2 + 2 * lr]);
        int ledY = Integer.parseInt(entry.model[12 + 2 + 2 * lr]);
        int rockerX = Integer.parseInt(entry.model[16 + 2 + 2 * lr]);
        int rockerY = Integer.parseInt(entry.model[17 + 2 + 2 * lr]);
        if (ledX <= 0 || ledY <= 0 ||
            rockerX <= 0 || rockerY <= 0) continue;

        double a1, b1;
        double a2, b2;
        {
          int dx = rockerX - ledX;
          int dy = rockerY - ledY;
          a1 = dx;
          b1 = dy;
          a2 = b1;
          b2 = -a1;
        }

        DecimalFormat fmt = new DecimalFormat("0.00");
        StringBuilder msg = new StringBuilder();
        for (int comp = 0; comp < 22; ++comp) {
          int compX = Integer.parseInt(entry.model[comp * 5 + 1 + 2 + 2 * lr]);
          int compY = Integer.parseInt(entry.model[comp * 5 + 2 + 2 + 2 * lr]);
          if (compX <= 0 || compY <= 0) {
          if (comp >= 6 && comp <= 14 || comp == 16 || comp == 18 || comp == 19 || comp == 21)
            msg.append("N/A ");
            continue;
          }
          int predX = (int) (compRatios[comp][0] * a1 + compRatios[comp][1] * a2) + ledX;
          int predY = (int) (compRatios[comp][0] * b1 + compRatios[comp][1] * b2) + ledY;
          predX -= rockerX;
          predY -= rockerY;
          double predLength = Math.sqrt(predX * predX + predY * predY);
          compX -= rockerX;
          compY -= rockerY;
          double compLength = Math.sqrt(compX * compX + compY * compY);
          double scaleConv = compLength / predLength;

          double crossProduct = predX * compY - predY * compX;
          crossProduct /= compLength;
          crossProduct /= predLength;
          double thetaConv = Math.asin(crossProduct);
          /*if (Double.isNaN(x) || Double.isNaN(y)) {
            System.out.println("Got NaN for " + entry.dir + "/" + entry.filename);
            System.exit(1);
          }*/
          if (comp >= 6 && comp <= 14 || comp == 16 || comp == 18 || comp == 19 || comp == 21)
          msg.append(fmt.format(scaleConv) + "," + fmt.format(thetaConv) + " ");
          vectorsA.get(comp).add(scaleConv);
          vectorsB.get(comp).add(thetaConv);
        }
        System.out.println(entry.dir + "/" + entry.filename + ": " + msg.toString());
      }
    }  // For each entry.
    DecimalFormat fmt = new DecimalFormat("0.00");
    System.out.print("double compRatios[22][2] = {");
    for (int comp = 0; comp < 22; ++comp) {
      if (comp % 4 == 0) {
        System.out.println();
        System.out.print("    ");
      }
      double sum = 0;
      for (Double d : vectorsA.get(comp)) {
        sum += d;
      }
      sum /= vectorsA.get(comp).size();
      double ang = 0;
      for (Double d : vectorsB.get(comp)) {
        ang += d;
      }
      ang /= vectorsB.get(comp).size();
      System.out.print("{" + fmt.format(sum) + "," + fmt.format(ang) + "}");
      if (comp != 21) {
        System.out.print(",");
      }
    }
    System.out.println("};");
  }

  // Calculate distance ratio against LED/Rocker and angle for each other component.
  static void calculateDistanceVectors(Map<String, TestEntry> entries) throws Exception {
    RobonautEyeTester.debug = false;
    ArrayList<ArrayList<Double>> lengthRatios = new ArrayList<ArrayList<Double>>();
    ArrayList<ArrayList<Double>> angles = new ArrayList<ArrayList<Double>>();
    for (int i = 0; i < 18; ++i) {
      lengthRatios.add(new ArrayList<Double>());
      angles.add(new ArrayList<Double>());
    }
    for (TestEntry entry : entries.values()) {
      int[] arrLeft = RobonautEyeTester.imageToArray(
          entry.dir + "/LeftImage_" + entry.filename);
      int[] arrRight = RobonautEyeTester.imageToArray(
          entry.dir + "/RightImage_" + entry.filename);
      int width = arrLeft[1];
      int height = arrLeft[0];

      for (int lr = 0; lr <= 1; ++lr) {
        int ledX = Integer.parseInt(entry.model[11 + 2 + 2 * lr]);
        int ledY = Integer.parseInt(entry.model[12 + 2 + 2 * lr]);
        int rockerX = Integer.parseInt(entry.model[16 + 2 + 2 * lr]);
        int rockerY = Integer.parseInt(entry.model[17 + 2 + 2 * lr]);
        if (ledX < 0 || ledY < 0 || rockerX < 0 || rockerY < 0) continue;

        int dx = rockerX - ledX;
        int dy = rockerY - ledY;
        double length = Math.sqrt(dx * dx + dy * dy);

        DecimalFormat fmt = new DecimalFormat("0.00");
        StringBuilder msg = new StringBuilder();
        for (int comp = 4; comp <= 21; ++comp) {
          int compX = Integer.parseInt(entry.model[comp * 5 + 1 + 2 + 2 * lr]);
          int compY = Integer.parseInt(entry.model[comp * 5 + 2 + 2 + 2 * lr]);
          if (compX < 0 || compY < 0) {
            msg.append("N/A ");
            continue;
          }
          int dx2 = compX - ledX;
          int dy2 = compY - ledY;
          double length2 = Math.sqrt(dx2 * dx2 + dy2 * dy2);
          double lengthRatio = length2 / length;
          double dot = dx * dx2 + dy * dy2;
          dot /= length;
          dot /= length2;
          if (Double.isNaN(lengthRatio) || Double.isNaN(dot)) {
            msg.append("N/A ");
            continue;
          }
          msg.append(fmt.format(lengthRatio) + "," + fmt.format(dot) + " ");
          lengthRatios.get(comp - 4).add(lengthRatio);
          angles.get(comp - 4).add(dot);
        }
        System.out.println(entry.dir + "/" + entry.filename + ": " + msg.toString());
      }
    }  // For each entry.
    DecimalFormat fmt = new DecimalFormat("0.00");
    System.out.print("double compRatios[18][2] = {");
    for (int comp = 0; comp < 18; ++comp) {
      if (comp % 4 == 0) {
        System.out.println();
        System.out.print("    ");
      }
      double sum = 0;
      for (Double d : lengthRatios.get(comp)) {
        sum += d;
      }
      sum /= lengthRatios.get(comp).size();
      double ang = 0;
      for (Double d : angles.get(comp)) {
        ang += d;
      }
      ang /= angles.get(comp).size();
      System.out.print("{" + fmt.format(sum) + "," + fmt.format(ang) + "}");
      if (comp != 17) {
        System.out.print(",");
      }
    }
    System.out.println("};");
  }

  public static void main(String[] args) throws Exception {
    long startTime = System.currentTimeMillis();
    Map<String, TestEntry> entries = new TreeMap<String, TestEntry>();

    loadEntries(entries);
    loadCsv(entries);
    //calculateBasisCoordinatesForGrid(entries);
    calculateBasisVectors(entries);
    //calculateDistanceVectors(entries);
    //calculateConversions(entries);

    long endTime = System.currentTimeMillis();
    System.out.println("Everything took " + (endTime - startTime) + " ms");
  }
}
