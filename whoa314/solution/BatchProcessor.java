import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import javax.imageio.*; 
import javax.swing.*;
import java.util.concurrent.*;

/**
 * This tool runs through all input files as listed inside Trainer.java, and loads model.csv
 * to test extracted basis-vector coordinates by drawing associated dots for all the taskboard
 * components on a new copy of each input image before saving the marked images to a new directory.
 */
public class BatchProcessor {
  static void processEntry(TestEntry entry) throws Exception {
    int[] arrLeft = RobonautEyeTester.imageToArray(
        entry.dir + "/LeftImage_" + entry.filename);
    int[] arrRight = RobonautEyeTester.imageToArray(
        entry.dir + "/RightImage_" + entry.filename);
    String baseName = entry.filename.substring(0, entry.filename.indexOf('.'));
    String dirOut = entry.dir + "Proc";
    {
      File dirOutFile = new File(dirOut);
      if (!dirOutFile.exists()) {
        dirOutFile.mkdir();
      }
    }

    String leftOut = dirOut + "/LeftImage_" + baseName + ".png";
    String rightOut = dirOut + "/RightImage_" + baseName + ".png";
    int width = arrLeft[1];
    int height = arrLeft[0];
    final BufferedImage imgLeft = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
    final BufferedImage imgRight = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
    for (int y = 0, cur = 2; y < height; ++y) {
      for (int x = 0; x < width; ++x, ++cur) {
        imgLeft.setRGB(x, y, arrLeft[cur] | 0xff000000);
        imgRight.setRGB(x, y, arrRight[cur] | 0xff000000);
      }
    }

    // Led -> Rocker and its normal.
    double[][] compRatios = {
        {0.01,-0.18},{-0.14,-0.18},{0.00,0.00},{1.00,0.00},
        {0.72,-0.01},{1.22,-0.03},{0.74,-1.23},{0.74,-0.97},
        {0.73,-0.71},{0.98,-1.21},{0.97,-0.95},{0.96,-0.69},
        {1.20,-1.19},{1.19,-0.93},{1.19,-0.68},{1.76,-1.09},
        {1.56,-1.11},{1.74,-0.53},{1.54,-0.53},{1.92,-0.51},
        {1.73,0.03},{1.54,0.04}};
    // Led -> Rocker and Rocker -> Grid.
//    double[][] compRatios = {
//        {0.02,0.20},{0.06,0.24},{0.00,0.00},{1.00,0.00},
//        {0.72,-0.01},{1.22,0.03},{0.79,1.30},{0.78,1.03},
//        {0.75,0.75},{1.01,1.27},{1.00,1.00},{0.99,0.73},
//        {1.23,1.24},{1.22,0.98},{1.20,0.71},{1.78,1.13},
//        {1.58,1.15},{1.74,0.53},{1.55,0.54},{1.92,0.52},
//        {1.71,-0.07},{1.51,-0.07}};
    for (int lr = 0; lr <= 1; ++lr) {
      int ledX = Integer.parseInt(entry.model[11 + 2 + 2 * lr]);
      int ledY = Integer.parseInt(entry.model[12 + 2 + 2 * lr]);
      int rockerX = Integer.parseInt(entry.model[16 + 2 + 2 * lr]);
      int rockerY = Integer.parseInt(entry.model[17 + 2 + 2 * lr]);
      int gridX = Integer.parseInt(entry.model[51 + 2 + 2 * lr]);
      int gridY = Integer.parseInt(entry.model[52 + 2 + 2 * lr]);
      if (ledX <= 0 || ledY <= 0 ||
          rockerX <= 0 || rockerY <= 0 ||
          gridX <= 0 || gridY <= 0) continue;
      double a1, b1;
      double a2, b2;
      /*{
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
      }*/
      int dx = rockerX - ledX;
      int dy = rockerY - ledY;
      a1 = dx;
      b1 = dy;
      a2 = b1;
      b2 = -a1;
      int lengthSq = dx * dx + dy * dy;
      if (lengthSq < 100) return;

      // Predicted grid pos.
      int predX = (int) (compRatios[10][0] * a1 + compRatios[10][1] * a2) + ledX;
      int predY = (int) (compRatios[10][0] * b1 + compRatios[10][1] * b2) + ledY;

      // Express predicted as vector from rocker.
      predX -= rockerX;
      predY -= rockerY;
      double predLength = Math.sqrt(predX * predX + predY * predY);

      int actX = gridX - rockerX;
      int actY = gridY - rockerY;
      double actLength = Math.sqrt(actX * actX + actY * actY);

      // We must scale the rocker -> comp vector by this much to meet actual from predicted.
      double scaleConv = actLength / predLength;

      double crossProduct = predX * actY - predY * actX;
      crossProduct /= actLength;
      crossProduct /= predLength;
      double thetaConv = Math.asin(crossProduct);

      /*double sinVal = dx / length;
      double theta = Math.asin(sinVal) * 2;
      double a3 = Math.cos(theta) * a2 - Math.sin(theta) * b2;
      double b3 = Math.sin(theta) * a2 + Math.cos(theta) * b2;
      a2 = a3;
      b2 = b3;*/
      /*{
        double a3 = Math.cos(thetaConv) * a2 - Math.sin(thetaConv) * b2;
        double b3 = Math.sin(thetaConv) * a2 + Math.cos(thetaConv) * b2;
        a3 *= scaleConv;
        b3 *= scaleConv;
        a2 = a3;
        b2 = b3;
      }*/
      for (int i = 0; i < 22; ++i) {
        int comp = i;
        int newX = (int) (compRatios[i][0] * a1 + compRatios[i][1] * a2);
        int newY = (int) (compRatios[i][0] * b1 + compRatios[i][1] * b2);
        newX += ledX;
        newY += ledY;
        /*{
          newX -= rockerX;
          newY -= rockerY;
          newX = (int) (newX * scaleConv);
          newY = (int) (newY * scaleConv);
          double rotX = newX * Math.cos(thetaConv) - newY * Math.sin(thetaConv);
          double rotY = newX * Math.sin(thetaConv) + newY * Math.cos(thetaConv);
          newX = (int) (rotX + rockerX);
          newY = (int) (rotY + rockerY);
        }*/
        IdeaTester.markPos(lr == 0 ? imgLeft : imgRight, newX, newY, "PP" + Integer.toString(comp));
      }
    }  // for left/right
    ImageIO.write(imgLeft, "png", new File(leftOut));
    ImageIO.write(imgRight, "png", new File(rightOut));
  }
  public static void main(String[] args) throws Exception {
    Map<String, TestEntry> entries = new TreeMap<String, TestEntry>();
    Trainer.loadEntries(entries);
    Trainer.loadCsv(entries);

    ScheduledThreadPoolExecutor pool = new ScheduledThreadPoolExecutor(8);
    final CountDownLatch latch = new CountDownLatch(entries.size());
    for (TestEntry entryBase : entries.values()) {
      final TestEntry entry = entryBase;
      pool.execute(new Runnable() {
        @Override
        public void run() {
          try {
            processEntry(entry);
          } catch (Exception e) {
            e.printStackTrace();
          } finally {
            latch.countDown();
          }
        }
      });
    }  // for each entry
    while (!latch.await(3, TimeUnit.SECONDS)) {
      System.out.println("Remaining: " + latch.getCount());
    }
    pool.shutdown();
  }
}
