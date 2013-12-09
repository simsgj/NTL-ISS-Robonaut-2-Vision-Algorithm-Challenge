import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import javax.imageio.*; 
import javax.swing.*;

/**
 * This file contains the logic for performing niblack categorization of pixels, including the
 * initial creation of the per-channel "Sum Images" to efficiently compute average and
 * standard deviation of each pixel. It displays the per-channel niblack results, as well as
 * the combined version ultimately used in the solution.
 */
class SumImage {
  long[][] sumR;
  long[][] sumG;
  long[][] sumB;
  long[][] sumH;
  long[][] sumSqR;
  long[][] sumSqG;
  long[][] sumSqB;
  long[][] sumSqH;
  long[][] sumL;
  long[][] sumSqL;
  boolean[][] isForeground;
  boolean[][] isBackground;
  boolean[][] isRed;
  boolean[][] isGreen;
  boolean[][] isBlue;
  public SumImage(int[] arr) {
    int width = arr[1];
    int height = arr[0];
    sumR = new long[width][height];
    sumG = new long[width][height];
    sumB = new long[width][height];
    sumH = new long[width][height];
    sumL = new long[width][height];
    sumSqR = new long[width][height];
    sumSqG = new long[width][height];
    sumSqB = new long[width][height];
    sumSqH = new long[width][height];
    sumSqL = new long[width][height];
    isForeground = new boolean[width][height];
    isBackground = new boolean[width][height];
    isRed = new boolean[width][height];
    isGreen = new boolean[width][height];
    isBlue = new boolean[width][height];
    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        int cur = 2 + x + y * width;
        int rgb = arr[cur];
        int r = (rgb >> 16) & 0x00ff;
        int g = (rgb >> 8) & 0x00ff;
        int b = (rgb) & 0x00ff;
        int hue = ImageProcessor.getHue(r, g, b);
        int lum = (int)((0.2126 * r) + (0.7152 * g) + (0.0722 * b));
        if (lum > 255) lum = 255;
        sumR[x][y] = r;
        sumG[x][y] = g;
        sumB[x][y] = b;
        sumH[x][y] = hue;
        sumL[x][y] = lum;
        sumSqR[x][y] = r * r;
        sumSqG[x][y] = g * g;
        sumSqB[x][y] = b * b;
        sumSqH[x][y] = hue * hue;
        sumSqL[x][y] = lum * lum;
        if (y - 1 >= 0) {
          sumR[x][y] += sumR[x][y - 1];
          sumG[x][y] += sumG[x][y - 1];
          sumB[x][y] += sumB[x][y - 1];
          sumH[x][y] += sumH[x][y - 1];
          sumL[x][y] += sumL[x][y - 1];
          sumSqR[x][y] += sumSqR[x][y - 1];
          sumSqG[x][y] += sumSqG[x][y - 1];
          sumSqB[x][y] += sumSqB[x][y - 1];
          sumSqH[x][y] += sumSqH[x][y - 1];
          sumSqL[x][y] += sumSqL[x][y - 1];
        }
        if (y - 1 >= 0 && x - 1 >= 0) {
          sumR[x][y] -= sumR[x - 1][y - 1];
          sumG[x][y] -= sumG[x - 1][y - 1];
          sumB[x][y] -= sumB[x - 1][y - 1];
          sumH[x][y] -= sumH[x - 1][y - 1];
          sumL[x][y] -= sumL[x - 1][y - 1];
          sumSqR[x][y] -= sumSqR[x - 1][y - 1];
          sumSqG[x][y] -= sumSqG[x - 1][y - 1];
          sumSqB[x][y] -= sumSqB[x - 1][y - 1];
          sumSqH[x][y] -= sumSqH[x - 1][y - 1];
          sumSqL[x][y] -= sumSqL[x - 1][y - 1];
        }
        if (x - 1 >= 0) {
          sumR[x][y] += sumR[x - 1][y];
          sumG[x][y] += sumG[x - 1][y];
          sumB[x][y] += sumB[x - 1][y];
          sumH[x][y] += sumH[x - 1][y];
          sumL[x][y] += sumL[x - 1][y];
          sumSqR[x][y] += sumSqR[x - 1][y];
          sumSqG[x][y] += sumSqG[x - 1][y];
          sumSqB[x][y] += sumSqB[x - 1][y];
          sumSqH[x][y] += sumSqH[x - 1][y];
          sumSqL[x][y] += sumSqL[x - 1][y];
        }
      }
    }
  }
}

public class ImageProcessor {
  static double UPPER_NIBLACK = 1.5;
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
      System.out.println("Unhandled hue: " + R + "/" + G + "/" + B);
      return 0;
    }
  }
  static int getAvg(long[][] sumImg, int x, int y, int radius) {
    int width = sumImg.length;
    int height = sumImg[0].length;
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

    long tally = sumImg[xmax - 1][ymax - 1];
    if (ymin - 1 >= 0) {
      // Top.
      tally -= sumImg[xmax - 1][ymin - 1];
    }
    if (xmin - 1 >= 0) {
      // Left.
      tally -= sumImg[xmin - 1][ymax - 1];
    }
    if (ymin - 1 >= 0 && xmin - 1 >= 0) {
      // Top left.
      tally += sumImg[xmin - 1][ymin - 1];
    }
    tally /= (xmax - xmin) * (ymax - ymin);
    return (int) tally;
  }
  static void doNiblackMulti(int[] arr, BufferedImage img, SumImage sumImage, int kernel) {
    int width = arr[1];
    int height = arr[0];
    long[][][] sums = new long[3][][];
    sums[2] = sumImage.sumR;
    sums[1] = sumImage.sumG;
    sums[0] = sumImage.sumB;
    long[][][] sumSqs = new long[3][][];
    sumSqs[2] = sumImage.sumSqR;
    sumSqs[1] = sumImage.sumSqG;
    sumSqs[0] = sumImage.sumSqB;
    boolean isSmall = (width < 1700 && height < 1700);
    for (int y = 0; y < height; ++y) {
      int avgR = getAvg(sums[2], 0, y, kernel); //, width, height);
      int avgG = getAvg(sums[1], 0, y, kernel); //, width, height);
      int avgB = getAvg(sums[0], 0, y, kernel); //, width, height);
      int avgSqR = getAvg(sumSqs[2], 0, y, kernel); //, width, height);
      int avgSqG = getAvg(sumSqs[1], 0, y, kernel); //, width, height);
      int avgSqB = getAvg(sumSqs[0], 0, y, kernel); //, width, height);
      int sigmaSquaredR = avgSqR - avgR * avgR;
      int sigmaSquaredG = avgSqG - avgG * avgG;
      int sigmaSquaredB = avgSqB - avgB * avgB;
      double stdDevR = Math.sqrt(sigmaSquaredR);
      double stdDevG = Math.sqrt(sigmaSquaredG);
      double stdDevB = Math.sqrt(sigmaSquaredB);
      for (int x = 0; x < width; ++x) {
        if (isSmall || ((x + y) & 0x01) == 0) {
          avgR = getAvg(sums[2], x, y, kernel); //, width, height);
          avgG = getAvg(sums[1], x, y, kernel); //, width, height);
          avgB = getAvg(sums[0], x, y, kernel); //, width, height);
          avgSqR = getAvg(sumSqs[2], x, y, kernel); //, width, height);
          avgSqG = getAvg(sumSqs[1], x, y, kernel); //, width, height);
          avgSqB = getAvg(sumSqs[0], x, y, kernel); //, width, height);
          sigmaSquaredR = avgSqR - avgR * avgR;
          sigmaSquaredG = avgSqG - avgG * avgG;
          sigmaSquaredB = avgSqB - avgB * avgB;
          stdDevR = Math.sqrt(sigmaSquaredR);
          stdDevG = Math.sqrt(sigmaSquaredG);
          stdDevB = Math.sqrt(sigmaSquaredB);
        }
        // 0 is b, 1 is g, 2 is r.
        int minChan = 999;
        int minCol = 999;
        int minAvg = 999;
        double minStdDev = 999;
        double minDev = 999;
        int maxChan = -1;
        int maxCol = -1;
        int maxAvg = -1;
        double maxStdDev = -1;
        double maxDev = -1;
        int cur = x + y * width;
        int R = (arr[2 + cur] >> 16) & 0x00ff;
        int G = (arr[2 + cur] >> 8) & 0x00ff;
        int B = (arr[2 + cur]) & 0x00ff;
        double devR = (R - avgR) / (double) stdDevR;
        double devG = (G - avgG) / (double) stdDevG;
        double devB = (B - avgB) / (double) stdDevB;
        if (devR <= devG && devR <= devB) {
          minChan = 2;
          minCol = R;
          minAvg = avgR;
          minStdDev = stdDevR;
          minDev = devR;
        } else if (devG <= devR && devG <= devB) {
          minChan = 1;
          minCol = G;
          minAvg = avgG;
          minStdDev = stdDevG;
          minDev = devG;
        } else {
          minChan = 0;
          minCol = B;
          minAvg = avgB;
          minStdDev = stdDevB;
          minDev = devB;
        }
        if (devR >= devG && devR >= devB) {
          maxChan = 2;
          maxCol = R;
          maxAvg = avgR;
          maxStdDev = stdDevR;
          maxDev = devR;
        } else if (devG >= devR && devG >= devB) {
          maxChan = 1;
          maxCol = G;
          maxAvg = avgG;
          maxStdDev = stdDevG;
          maxDev = devG;
        } else {
          maxChan = 0;
          maxCol = B;
          maxAvg = avgB;
          maxStdDev = stdDevB;
          maxDev = devB;
        }
        if (Math.abs(maxDev) > Math.abs(minDev) && maxCol > maxAvg + UPPER_NIBLACK * maxStdDev) {
          int shift = maxChan * 8;
          img.setRGB(x, y, 0xff000000 | (0x000000ff << shift));
          sumImage.isForeground[x][y] = true;
          if (maxChan == 2) {
            sumImage.isRed[x][y] = true;
          } else if (maxChan == 1) {
            sumImage.isGreen[x][y] = true;
          } else {
            sumImage.isBlue[x][y] = true;
          }
        } else if (Math.abs(maxDev) < Math.abs(minDev) && minCol < minAvg - .75 * minStdDev) {
          int shift = maxChan * 8;
          img.setRGB(x, y, 0xff000000 | (0x000000af << shift));
          //img.setRGB(x, y, 0xff000000);
          sumImage.isBackground[x][y] = true;
          if (maxChan == 2) {
            sumImage.isRed[x][y] = true;
          } else if (maxChan == 1) {
            sumImage.isGreen[x][y] = true;
          } else {
            sumImage.isBlue[x][y] = true;
          }
        } else {
          img.setRGB(x, y, 0xffffffff);
        }

        /*int avgCol = (avgR + avgG + avgB) / 3;

        int col = (R + G + B) / 3;
        int dR = R - avgR;
        int dG = G - avgG;
        int dB = B - avgB;
        double devR = sigmaSquaredR > 0 ? dR * dR / (double) sigmaSquaredR : 0;
        double devG = sigmaSquaredG > 0 ? dG * dG / (double) sigmaSquaredG : 0;
        double devB = sigmaSquaredB > 0 ? dB * dB / (double) sigmaSquaredB : 0;
        double avgDev = devR + devG + devB;
        int toPaint = 0xffffffff;
        if (R > Math.max(G, B)) {
          toPaint = 0xffff0000;
        } else if (G > Math.max(R, B)) {
          toPaint = 0xff00ff00;
        } else {
          toPaint = 0xff0000ff;
        }
        avgDev /= 3;
        avgDev = Math.sqrt(avgDev);
        if (col > avgCol && avgDev > 1) {
          img.setRGB(x, y, toPaint);
          sumImage.isForeground[x][y] = true;
        } else if (col < avgCol && avgDev > .5) {
          img.setRGB(x, y, 0xff000000);
          sumImage.isBackground[x][y] = true;
        } else {
          img.setRGB(x, y, 0xffffffff);
        }*/
      }
    }
  }
  static void doNiblackHue(
      int[] arr, BufferedImage img, long[][] sum, long[][] sumSq, int kernel, int shift) {
    int width = sum.length;
    int height = sum[0].length;
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        long avg = getAvg(sum, x, y, kernel);
        long avgSq = getAvg(sumSq, x, y, kernel);
        long sigmaSquared = avgSq - avg * avg;
        if (sigmaSquared < 0) {
          System.out.println("Got negative sigmaSquared! " + sigmaSquared);
          System.exit(1);
        }
        long stdDev = (long) Math.sqrt(sigmaSquared);
        int r = (arr[2 + x + y * width] >> shift) & 0x00ff;
        if (r < avg - .75 * stdDev) {
          //img.setRGB(x, y, 0xff000000 | (0x00ff << shift));
          img.setRGB(x, y, 0xff000000);
        } else if (r > avg + .75 * stdDev) {
          img.setRGB(x, y, 0xff0000ff);
          //img.setRGB(x, y, 0xff000000 | (0x00ff << shift));
         // img.setRGB(x, y, 0xffffffff);
        } else {
         // img.setRGB(x, y, 0xff777777);
          img.setRGB(x, y, 0xffffffff);
        }
      }
    }
  }
  static void doNiblack(
      int[] arr, BufferedImage img, long[][] sum, long[][] sumSq, int kernel, int shift) {
    int width = sum.length;
    int height = sum[0].length;
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        long avg = getAvg(sum, x, y, kernel);
        long avgSq = getAvg(sumSq, x, y, kernel);
        long sigmaSquared = avgSq - avg * avg;
        if (sigmaSquared < 0) {
          System.out.println("Got negative sigmaSquared! " + sigmaSquared);
          System.exit(1);
        }
        long stdDev = (long) Math.sqrt(sigmaSquared);
        int r = (arr[2 + x + y * width] >> shift) & 0x00ff;
        if (r > avg + 3 * stdDev) {
          img.setRGB(x, y, 0xff000000 | (0x00ff << shift));
        } else if (r < avg - .75 * stdDev) {
          //img.setRGB(x, y, 0xff000000 | (0x00ff << shift));
          img.setRGB(x, y, 0xff000000);
        } else if (r > avg + 1.5 * stdDev) {
          img.setRGB(x, y, 0xff000000 | (0x00ff << shift));
         // img.setRGB(x, y, 0xffffffff);
        } else {
         // img.setRGB(x, y, 0xff777777);
          img.setRGB(x, y, 0xffffffff);
        }
      }
    }
  }
  public static void main(String[] args) throws Exception {
    RobonautEyeTester.debug = false;
    int[] arrLeft = RobonautEyeTester.imageToArray(args[0]);
    int[] arrRight = RobonautEyeTester.imageToArray(args[0].replaceFirst("Left", "Right"));
    int width = arrLeft[1];
    int height = arrLeft[0];
    final BufferedImage imgLeftR = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
    final BufferedImage imgLeftG = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
    final BufferedImage imgLeftB = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
    final BufferedImage imgLeftH = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
    final BufferedImage imgLeftL = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
    final BufferedImage imgRightR = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
    final BufferedImage imgRightG = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
    final BufferedImage imgRightB = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
    final BufferedImage imgRightH = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
    final BufferedImage imgRightL = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
    for (int y = 0, cur = 2; y < height; ++y) {
      for (int x = 0; x < width; ++x, ++cur) {
        imgLeftR.setRGB(x, y, (0xffff0000 & (arrLeft[cur] | 0xff000000)));
        imgLeftG.setRGB(x, y, (0xff00ff00 & (arrLeft[cur] | 0xff000000)));
        imgLeftB.setRGB(x, y, (0xff0000ff & (arrLeft[cur] | 0xff000000)));
        imgLeftL.setRGB(x, y, (0xffffffff & (arrLeft[cur] | 0xff000000)));
        imgRightR.setRGB(x, y, (0xffff0000 & (arrRight[cur] | 0xff000000)));
        imgRightG.setRGB(x, y, (0xff00ff00 & (arrRight[cur] | 0xff000000)));
        imgRightB.setRGB(x, y, (0xff0000ff & (arrRight[cur] | 0xff000000)));
        imgRightL.setRGB(x, y, (0xffffffff & (arrRight[cur] | 0xff000000)));
      }
    }


    // Ideas here.
    long startTime = System.currentTimeMillis();
    SumImage sumLeft = new SumImage(arrLeft);
    SumImage sumRight = new SumImage(arrRight);
    /*for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        imgLeftR.setRGB(x, y, 0xff000000 | getAvg(sumLeft.sumR, x, y, width / 10));
        imgLeftG.setRGB(x, y, 0xff000000 | getAvg(sumLeft.sumG, x, y, width / 10));
        imgLeftB.setRGB(x, y, 0xff000000 | getAvg(sumLeft.sumB, x, y, width / 10));
        imgRightR.setRGB(x, y, 0xff000000 | getAvg(sumRight.sumR, x, y, width / 10));
        imgRightG.setRGB(x, y, 0xff000000 | getAvg(sumRight.sumG, x, y, width / 10));
        imgRightB.setRGB(x, y, 0xff000000 | getAvg(sumRight.sumB, x, y, width / 10));
      }
    }*/
    int kernel = width / 5;
    {
      doNiblack(arrLeft, imgLeftR, sumLeft.sumR, sumLeft.sumSqR, kernel, 16);
      doNiblack(arrLeft, imgLeftG, sumLeft.sumG, sumLeft.sumSqG, kernel, 8);
      doNiblack(arrLeft, imgLeftB, sumLeft.sumB, sumLeft.sumSqB, kernel, 0);
      //doNiblackHue(arrLeft, imgLeftB, sumLeft.sumH, sumLeft.sumSqH, kernel, 0);
      //doNiblack(arrLeft, imgLeftB, sumLeft.sumL, sumLeft.sumSqL, kernel, 0);
      doNiblack(arrRight, imgRightR, sumRight.sumR, sumRight.sumSqR, kernel, 16);
      doNiblack(arrRight, imgRightG, sumRight.sumG, sumRight.sumSqG, kernel, 8);
      doNiblack(arrRight, imgRightB, sumRight.sumB, sumRight.sumSqB, kernel, 0);
      //doNiblackHue(arrRight, imgRightB, sumRight.sumH, sumRight.sumSqH, kernel, 0);
      //doNiblack(arrRight, imgRightB, sumRight.sumL, sumRight.sumSqL, kernel, 0);
    }
    {
      doNiblackMulti(arrLeft, imgLeftL, sumLeft, kernel);
      doNiblackMulti(arrRight, imgRightL, sumRight, kernel);
    }

    long finishTime = System.currentTimeMillis();
    System.out.println("Idea took: " + (finishTime - startTime) + " ms");

    JFrame frame = new JFrame("Idea Tester");
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    JPanel panel = new JPanel() {
      @Override
      public void paintComponent(Graphics g) {
        if (g == null) return;
        g.drawImage(imgLeftR.getScaledInstance(
            getWidth() / 4, getHeight() / 2, Image.SCALE_SMOOTH), 0, 0, null);
        g.drawImage(imgLeftG.getScaledInstance(
            getWidth() / 4, getHeight() / 2, Image.SCALE_SMOOTH), getWidth() / 4, 0, null);
        g.drawImage(imgLeftB.getScaledInstance(
            getWidth() / 4, getHeight() / 2, Image.SCALE_SMOOTH), getWidth() * 2 / 4, 0, null);
        g.drawImage(imgLeftL.getScaledInstance(
            getWidth() / 4, getHeight() / 2, Image.SCALE_SMOOTH), getWidth() * 3 / 4, 0, null);

        g.drawImage(imgRightR.getScaledInstance(
            getWidth() / 4, getHeight() / 2, Image.SCALE_SMOOTH), 0, getHeight() / 2, null);
        g.drawImage(imgRightG.getScaledInstance(
            getWidth() / 4, getHeight() / 2, Image.SCALE_SMOOTH), getWidth() / 4, getHeight() / 2, null);
        g.drawImage(imgRightB.getScaledInstance(
            getWidth() / 4, getHeight() / 2, Image.SCALE_SMOOTH), getWidth() * 2 / 4, getHeight() / 2, null);
        g.drawImage(imgRightL.getScaledInstance(
            getWidth() / 4, getHeight() / 2, Image.SCALE_SMOOTH), getWidth() * 3 / 4, getHeight() / 2, null);
      }
      @Override
      public void repaint() {
        //super.repaint();
        paintComponent(getGraphics());
      }
    };
    panel.setPreferredSize(new Dimension(imgLeftR.getWidth() * 4 / 3, imgLeftR.getHeight() * 2 / 3));

    frame.getContentPane().add(panel);
    frame.pack();
    frame.setVisible(true);
  }
}
