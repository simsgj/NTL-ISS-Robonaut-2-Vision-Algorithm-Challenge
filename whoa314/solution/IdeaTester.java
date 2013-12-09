import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import javax.imageio.*; 
import javax.swing.*;

class SccProperty {
  int xmin, xmax, ymin, ymax;
  int xminy, xmaxy;
  int yminx, ymaxx;

  // Average.
  int x, y;
  int size;
  int r, g, b;
  int propNum;

  boolean isGreenLed;
  boolean isBlueLed;
  boolean isRockerSwitch;
  boolean containsCenter;

  boolean isFirstFromLeft;
}

class ModelObject {
  boolean leftFound, rightFound;
  String state;
  int leftX, leftY;
  int rightX, rightY;
  SccProperty leftProp, rightProp;
  public ModelObject() {
    leftFound = false;
    rightFound = false;
    state = "HIDDEN";
    leftX = 0;
    leftY = 0;
    rightX = 0;
    rightY = 0;
    leftProp = null;
    rightProp = null;
  }
}

class ModelSummary {
  ModelObject[] objects = new ModelObject[22];
  public ModelSummary() {
    for (int i = 0; i < objects.length; ++i) {
      objects[i] = new ModelObject();
    }
  }
}

public class IdeaTester {
  static double[][] compRatios10 = {
      {0.02,0.20},{0.06,0.24},{0.00,0.00},{1.00,0.00},
      {0.72,-0.01},{1.22,0.03},{0.79,1.30},{0.78,1.03},
      {0.75,0.75},{1.01,1.27},{1.00,1.00},{0.99,0.73},
      {1.23,1.24},{1.22,0.98},{1.20,0.71},{1.78,1.13},
      {1.58,1.15},{1.74,0.53},{1.55,0.54},{1.92,0.52},
      {1.71,-0.07},{1.51,-0.07}};
  static double[][] compRatiosOrthoSmallLeft = {
      {-0.01,-0.17},{-0.21,-0.10},{0.00,0.00},{1.00,0.00},
      {0.75,-0.02},{1.22,-0.04},{0.71,-1.21},{0.71,-0.95},
      {0.72,-0.69},{0.96,-1.20},{0.97,-0.94},{0.97,-0.69},
      {1.21,-1.19},{1.21,-0.94},{1.21,-0.68},{1.84,-1.12},
      {1.62,-1.13},{1.84,-0.55},{1.62,-0.55},{2.07,-0.55},
      {1.84,0.00},{1.63,0.01}};
  static double[][] compRatiosOrthoSmallRight = {
      {0.03,-0.19},{-0.09,-0.10},{0.00,0.00},{1.00,0.00},
      {0.73,-0.01},{1.22,-0.01},{0.88,-1.30},{0.85,-1.02},
      {0.81,-0.74},{1.12,-1.29},{1.09,-1.00},{1.05,-0.73},
      {1.36,-1.27},{1.31,-0.99},{1.28,-0.72},{1.92,-1.16},
      {1.72,-1.19},{1.83,-0.57},{1.64,-0.58},{2.02,-0.57},
      {1.75,0.01},{1.55,0.02}};
  static double[][] compRatiosOrthoLargeLeft = {
      {-0.02,-0.19},{0.15,-0.31},{0.00,0.00},{1.00,0.00},
      {0.71,0.02},{1.22,-0.03},{0.55,-1.27},{0.59,-1.01},
      {0.61,-0.74},{0.77,-1.24},{0.80,-0.97},{0.84,-0.71},
      {0.98,-1.20},{1.01,-0.94},{1.04,-0.69},{1.49,-1.07},
      {1.31,-1.09},{1.56,-0.49},{1.39,-0.50},{1.73,-0.47},
      {1.66,0.10},{1.48,0.10}};
  static double[][] compRatiosOrthoLargeRight = {
      {-0.04,-0.18},{0.09,-0.39},{0.00,0.00},{1.00,0.00},
      {0.72,0.04},{1.22,-0.03},{0.39,-1.13},{0.46,-0.89},
      {0.52,-0.65},{0.61,-1.09},{0.68,-0.86},{0.75,-0.62},
      {0.83,-1.05},{0.89,-0.82},{0.96,-0.59},{1.35,-0.92},
      {1.17,-0.94},{1.52,-0.39},{1.33,-0.41},{1.62,-0.28},
      {1.70,0.14},{1.51,0.14}};
  static double[][] compRatiosOrthoLeft = {
      {-0.01,-0.18},{-0.02,-0.21},{0.00,0.00},{1.00,0.00},
      {0.73,0.00},{1.22,-0.04},{0.63,-1.24},{0.65,-0.98},
      {0.66,-0.72},{0.86,-1.22},{0.88,-0.96},{0.90,-0.70},
      {1.09,-1.19},{1.11,-0.94},{1.13,-0.69},{1.66,-1.09},
      {1.46,-1.11},{1.70,-0.52},{1.50,-0.53},{1.88,-0.51},
      {1.75,0.05},{1.55,0.06}};
  static double[][] compRatiosOrthoRight = {
      {0.01,-0.18},{-0.02,-0.21},{0.00,0.00},{1.00,0.00},
      {0.73,0.01},{1.22,-0.01},{0.69,-1.24},{0.70,-0.97},
      {0.70,-0.71},{0.92,-1.21},{0.93,-0.95},{0.94,-0.69},
      {1.16,-1.19},{1.15,-0.93},{1.16,-0.67},{1.70,-1.07},
      {1.51,-1.09},{1.71,-0.50},{1.52,-0.51},{1.89,-0.48},
      {1.73,0.06},{1.53,0.06}};
  // Older one.
  static double[][] compRatiosOrthoOld = {
      {0.01,-0.18},{-0.14,-0.18},{0.00,0.00},{1.00,0.00},
      {0.72,-0.01},{1.22,-0.03},{0.74,-1.23},{0.74,-0.97},
      {0.73,-0.71},{0.98,-1.21},{0.97,-0.95},{0.96,-0.69},
      {1.20,-1.19},{1.19,-0.93},{1.19,-0.68},{1.76,-1.09},
      {1.56,-1.11},{1.74,-0.53},{1.54,-0.53},{1.92,-0.51},
      {1.73,0.03},{1.54,0.04}};
  static void drawParallels(BufferedImage imgLeft, BufferedImage imgRight) {
    int width = imgLeft.getWidth();
    int height = imgLeft.getHeight();
    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        if ((x + y) % 2 == 0) {
          int tmp = imgLeft.getRGB(x, y);
          imgLeft.setRGB(x, y, imgRight.getRGB(x, y));
          imgRight.setRGB(x, y, tmp);
        }
      }
    }
  }

  // Reachability as defined by the straight line connecting the two not intersecting a "large" component.
  static boolean canReach(BufferedImage img, int width, int height,
      SccProperty prop0, SccProperty prop1, int[][] scc,
      Map<Integer, SccProperty> sccPropLarge) {
    if (prop0.x < 0 || prop0.x >= width ||
        prop1.x < 0 || prop1.x >= width ||
        prop0.y < 0 || prop0.y >= height ||
        prop1.y < 0 || prop1.y >= height) {
      System.out.println("Invalid prop out of bounds!! " + prop0.x + "," + prop0.y + " " +
          prop1.x + "," + prop1.y);
      System.exit(1);
      return false;
    }
    return canReach(img, width, height, prop0.x, prop0.y, prop1.x, prop1.y, scc, sccPropLarge);
  }
  static boolean canReach(BufferedImage img, int width, int height,
      int x0, int y0, int x1, int y1, int[][] scc,
      Map<Integer, SccProperty> sccPropLarge) {
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
    int dist = (int) Math.sqrt(dx * dx + dy * dy);
    if (dist == 0) return true;
    for (int i = 0; i < dist; ++i) {
      int x = dx * i / dist + x0;
      int y = dy * i / dist + y0;
      //img.setRGB(x, y, Color.green.getRGB());
      if (sccPropLarge.containsKey(scc[x][y])) {
        return false;
      }
    }
    return true;
  }

  static boolean mightBeRocker(BufferedImage img, int width, int height, int[][] scc,
      int ymin, SccProperty prop) {
    if (prop.xmax - prop.xmin < width / 15 &&
        prop.ymax - prop.ymin < height / 10 &&
        Math.max(prop.r, Math.max(prop.g, prop.b)) < 100 &&
        //prop.ymax - prop.ymin > prop.xmax - prop.xmin &&
        //prop.y > coverY + height / 10) {
        prop.y > ymin) {
      //maxSize = prop.size;
      //maxProp = prop;
      int numInside = 0;
      int count = 0;
      for (int dy = Math.max(0, prop.y - 20); dy < Math.min(height, prop.y + 20); ++dy) {
        for (int dx = Math.max(0, prop.x - 10); dx < Math.min(width, prop.x + 10); ++dx) {
          ++count;
          if (scc[dx][dy] == prop.propNum) {
            ++numInside;
          }
          //img.setRGB(dx, dy, Color.yellow.getRGB());
        }
      }
      double fillFraction = numInside / (double)count;
      if (fillFraction > .8) {
        System.out.println("Candidate rocker switch with fill fraction: " + fillFraction +
            " at " + prop.x + "," + prop.y);
        return true;
      }
    }
    return false;
  }

  // Filtering turns out to be useful just for ISS images.
  // Downsides? If none, we should run it anyway.
  static int[] filter(int[] arr) {
    int width = arr[1];
    int height = arr[0];
    int[] ret = new int[arr.length];
    ret[0] = arr[0];
    ret[1] = arr[1];
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        int rs = 255;
        int gs = 255;
        int bs = 255;
        for (int dy = y - 1; dy <= y + 1; ++dy) {
          if (dy < 0) continue;
          if (dy >= height) continue;
          for (int dx = x - 1; dx <= x + 1; ++dx) {
            if (dx < 0) continue;
            if (dx >= width) continue;
            if (dy == y && dx == x) continue;
            int cur = 2 + dy * width + dx;
            rs = Math.min(rs, (arr[cur] >> 16) & 0x00ff);
            gs = Math.min(gs, (arr[cur] >> 8) & 0x00ff);
            bs = Math.min(bs, (arr[cur]) & 0x00ff);
          }
        }
        int cur = 2 + y * width + x;
        ret[cur] = ((rs & 0x00ff) << 16) | ((gs & 0x00ff) << 8) | (bs & 0x00ff);
      }
    }
    return ret;
  }

  static class Pixel {
    int x, y, diff;
    public Pixel(int x, int y, int diff) {
      this.x = x; this.y = y; this.diff = diff;
    }
  }
  static long getAverageContrast(int[] arr) {
    int width = arr[1];
    int height = arr[0];
    int count = 0;
    long diffSum = 0;
    long maxDiff = 0;
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
            long dr = drs - rs;
            long dg = dgs - gs;
            long db = dbs - bs;
            dr *= dr;
            dg *= dg;
            db *= db;
            long diff = dr + dg + db;

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
    System.out.println("Average contrast " + diffSum);
    System.out.println("Max contrast " + maxDiff);
    return diffSum;
  }

  static long getAverageLum(int[] arr) {
    int width = arr[1];
    int height = arr[0];
    long totalLum = 0;
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
  static int[] mixLum(int[] arr) {
    int[] ret = new int[arr.length];
    ret[0] = arr[0];
    ret[1] = arr[1];
    int width = arr[1];
    int height = arr[0];
    long totalLum = getAverageLum(arr);
    for (int y = 0, cur = 2; y < height; ++y) {
      for (int x = 0; x < width; ++x, ++cur) {
        int r = (arr[cur] >> 16) & 0x00ff;
        int g = (arr[cur] >> 8) & 0x00ff;
        int b = (arr[cur]) & 0x00ff;
        int curLum = r + g + b;
        double scale = totalLum / (double) (curLum);
        r = (int)(r * scale);
        g = (int)(g * scale);
        b = (int)(b * scale);
        if (r < 0) r = 0;
        if (r > 255) r = 255;
        if (g < 0) g = 0;
        if (g > 255) g = 255;
        if (b < 0) b = 0;
        if (b > 255) b = 255;
        ret[cur] = (r << 16) | (g << 8) | b;
      }
    }
    return ret;
  }

  static void findHighContrast(int[] arr, BufferedImage img) {
    int width = arr[1];
    int height = arr[0];
    ArrayList<Pixel> pixels = new ArrayList<Pixel>();
    for (int y = 1; y < height - 1; ++y) {
      for (int x = 1; x < width - 1; ++x) {
        int cur = 2 + y * width + x;
        int rs = (arr[cur] >> 16) & 0x00ff;
        int gs = (arr[cur] >> 8) & 0x00ff;
        int bs = (arr[cur]) & 0x00ff;
        int diff = 0;
        for (int dy = y - 1; dy <= y + 1; ++dy) {
          for (int dx = x - 1; dx <= x + 1; ++dx) {
            if (dy == y && dx == x) continue;
            int dcur = 2 + dy * width + dx;
            int drs = (arr[dcur] >> 16) & 0x00ff;
            int dgs = (arr[dcur] >> 8) & 0x00ff;
            int dbs = (arr[dcur]) & 0x00ff;
            int dr = drs - rs;
            int dg = dgs - gs;
            int db = dbs - bs;
            int diff2 = dr * dr + dg * dg + db * db;
            if (diff2 > diff) {
              diff = diff2;
            }
          }
        }
        pixels.add(new Pixel(x, y, diff));
      }
    }
    // Descending order.
    Collections.sort(pixels, new Comparator<Pixel>() {
      @Override
      public int compare(Pixel a, Pixel b) {
        return Integer.valueOf(b.diff).compareTo(Integer.valueOf(a.diff));
      }
      @Override
      public boolean equals(Object o) {
        return false;
      }
    });
    for (int i = 0; i < pixels.size() / 10; ++i) {
      img.setRGB(pixels.get(i).x, pixels.get(i).y, Color.white.getRGB());
    }
  }

  // Doesn't work very well yet.
  static void findWhite(int[] arr, BufferedImage img, int[] xpos, int[] ypos) {
    int width = arr[1];
    int height = arr[0];
    int numWhite = 0;
    boolean[][] isWhite = new boolean[width][height];
    for (int y = 0, cur = 2; y < height; ++y) {
      for (int x = 0; x < width; ++x, ++cur) {
        int r = (arr[cur] >> 16) & 0x00ff;
        int g = (arr[cur] >> 8) & 0x00ff;
        int b = (arr[cur]) & 0x00ff;
        int max = Math.max(r, Math.max(g, b));
        int min = Math.min(r, Math.min(g, b));
        if (min > 100) {// && (max - min) / (double) max < 0.3) {
        //if ((max - min) / (double) max < 0.2) {
          img.setRGB(x, y, Color.red.getRGB());
          ++numWhite;
          isWhite[x][y] = true;
        } else {
          //img.setRGB(x, y, Color.white.getRGB());
        }
      }
    }
    System.out.println("Num red points: " + numWhite);
    int kerX = width / 80;
    int kerY = height / 40;
    double samp = (kerX * 2 + 1) * (kerY * 2 + 1);
    double maxFrac = 0;
    int bestX = 0;
    int bestY = 0;
    for (int x = kerX; x < width - kerX; ++x) {
      for (int y = kerY; y < height - kerY; ++y) {
        if (!isWhite[x][y]) continue;
        int count = 0;
        for (int dx = x - kerX; dx <= x + kerX; ++dx) {
          for (int dy = y - kerY; dy <= y + kerY; ++dy) {
            if (isWhite[dx][dy]) ++count;
          }
        }
        double redFrac = count / samp;
        if (redFrac > maxFrac) {
          maxFrac = redFrac;
          bestX = x;
          bestY = y;
        }
      }
    }
    markPos(img, bestX, bestY, "0");
    xpos[0] = bestX;
    ypos[0] = bestY;
  }

  static boolean allowMarkPos = true;

  static void markPos(BufferedImage img, int xpos, int ypos, String msg) {
    if (!allowMarkPos) return;
    Graphics g = img.getGraphics();
    g.setColor(Color.green);
    g.fillOval(xpos - 10, ypos - 10, 20, 20);
    g.setColor(Color.black);
    g.setFont(new Font("Times", Font.BOLD, 36));
    g.drawString(msg, xpos, ypos);
  }

  static void markPosPink(BufferedImage img, int xpos, int ypos, String msg) {
    if (!allowMarkPos) return;
    Graphics g = img.getGraphics();
    g.setColor(Color.pink);
    g.fillOval(xpos - 8, ypos - 8, 16, 16);
    g.setColor(Color.black);
    g.setFont(new Font("Times", Font.BOLD, 36));
    g.drawString(msg, xpos, ypos);
  }

  static void markPosYellow(BufferedImage img, int xpos, int ypos, String msg) {
    if (!allowMarkPos) return;
    Graphics g = img.getGraphics();
    g.setColor(Color.yellow);
    g.fillOval(xpos - 8, ypos - 8, 16, 16);
    g.setColor(Color.black);
    g.setFont(new Font("Times", Font.BOLD, 36));
    g.drawString(msg, xpos, ypos);
  }

  static void markPosRed(BufferedImage img, int xpos, int ypos, String msg) {
    if (!allowMarkPos) return;
    Graphics g = img.getGraphics();
    g.setColor(Color.red);
    g.fillOval(xpos - 8, ypos - 8, 16, 16);
    g.setColor(Color.black);
    g.setFont(new Font("Times", Font.BOLD, 36));
    g.drawString(msg, xpos, ypos);
  }

  static void markPosBlue(BufferedImage img, int xpos, int ypos, String msg) {
    if (!allowMarkPos) return;
    Graphics g = img.getGraphics();
    g.setColor(Color.blue);
    g.fillOval(xpos - 10, ypos - 10, 20, 20);
    g.setColor(Color.black);
    g.setFont(new Font("Times", Font.BOLD, 36));
    g.drawString(msg, xpos, ypos);
  }

  static void findRed(int[] arr, BufferedImage img, ModelSummary model, boolean isLeft,
      int[][] scc, Map<Integer, SccProperty> sccProp, Map<Integer, SccProperty> sccPropLarge,
      int xmin, int ymin, int xmax, int ymax, double requiredFraction) {
    int width = arr[1];
    int height = arr[0];
    xmin = Math.max(0, xmin);
    ymin = Math.max(0, ymin);
    xmax = Math.min(xmax, width);
    ymax = Math.min(ymax, height);
    int numRed = 0;
    boolean[][] isRed = new boolean[width][height];
    int[][] redSum = new int[width][height];
    // TODO(dhuo): Try to remove the need to limit where to look once there are better
    // ways to backtrack bad guesses at the red cover.
    for (int y = 0, cur = 2; y < height; ++y) {
      for (int x = 0; x < width; ++x, ++cur) {
        if (x < xmin || x >= xmax || y < ymin || y >= ymax) continue;
        if (sccPropLarge.containsKey(scc[x][y])) {
          if (y - 1 >= 0) {
            redSum[x][y] += redSum[x][y - 1];
          }
          if (x - 1 >= 0) {
            redSum[x][y] += redSum[x - 1][y];
          }
          if (y - 1 >= 0 && x - 1 >= 0) {
            redSum[x][y] -= redSum[x - 1][y - 1];
          }
          continue;
        }
        long r = (arr[cur] >> 16) & 0x00ff;
        long g = (arr[cur] >> 8) & 0x00ff;
        long b = (arr[cur]) & 0x00ff;
        r *= r;
        r *= r;
        g *= g;
        g *= g;
        b *= b;
        b *= b;
        if ((r > 0 && Math.max(g, b) == 0) || r / (double)Math.max(g, b) >= 8) {
          img.setRGB(x, y, Color.red.getRGB());
          ++numRed;
          isRed[x][y] = true;
          redSum[x][y] = 1;
        } else {
          //img.setRGB(x, y, Color.white.getRGB());
        }
        if (y - 1 >= 0) {
          redSum[x][y] += redSum[x][y - 1];
        }
        if (x - 1 >= 0) {
          redSum[x][y] += redSum[x - 1][y];
        }
        if (y - 1 >= 0 && x - 1 >= 0) {
          redSum[x][y] -= redSum[x - 1][y - 1];
        }
      }
    }
    System.out.println("Num red points: " + numRed);
    int kerX = width / 80;
    int kerY = height / 40;
    double samp = (kerX * 2 + 1) * (kerY * 2 + 1);
    double maxFrac = 0;
    int bestX = 0;
    int bestY = 0;
    for (int x = kerX + xmin + 1; x < xmax - kerX - 1; ++x) {
      for (int y = kerY + ymin + 1; y < ymax - kerY - 1; ++y) {
        if (!isRed[x][y]) continue;
        int count = redSum[x + kerX][(y + kerY)] -
                    redSum[x + kerX][(y - kerY - 1)] -
                    redSum[x - kerX - 1][(y + kerY)] +
                    redSum[x - kerX - 1][(y - kerY - 1)];
        /*int count = 0;
        for (int dx = x - kerX; dx <= x + kerX; ++dx) {
          for (int dy = y - kerY; dy <= y + kerY; ++dy) {
            if (isRed[dx][dy]) ++count;
          }
        }*/
        double redFrac = count / samp;
        if (redFrac > maxFrac) {
          maxFrac = redFrac;
          bestX = x;
          bestY = y;
        }
      }
    }

    if (maxFrac > requiredFraction) {
      System.out.println(
          "Found red cover in range: " + xmin + "," + ymin + "," + xmax + "," + ymax +
          " with fraction " + maxFrac + " out of " + requiredFraction);
      Map<Integer, Integer> propCount = new HashMap<Integer, Integer>();
      for (int x = bestX - kerX; x <= bestX + kerX; ++x) {
        for (int y = bestY - kerY; y <= bestY + kerY; ++y) {
          SccProperty checkProp = sccProp.get(scc[x][y]);
          if (checkProp != null) {
            if (propCount.get(scc[x][y]) == null) {
              propCount.put(scc[x][y], 0);
            }
            propCount.put(scc[x][y], propCount.get(scc[x][y]) + 1);
          }
        }
      }
      SccProperty bestProp = null;
      int bestPropValue = -1;
      for (Integer prop : propCount.keySet()) {
        if (propCount.get(prop) > bestPropValue) {
          bestPropValue = propCount.get(prop);
          bestProp = sccProp.get(prop);
        }
      }
      /*Graphics g = img.getGraphics();
      g.setColor(Color.green);
      g.fillRect(bestX - kerX, bestY - kerY, 2 * kerX, 2 * kerY);*/
      if (bestProp != null) {
        System.out.println("Found SCC corresponding to cover with size " + bestProp.size);
        bestX = bestProp.x;
        bestY = bestProp.y;
        markPosBlue(img, bestProp.yminx, bestProp.ymin, "");
        markPosBlue(img, bestProp.ymaxx, bestProp.ymax, "");
        markPosBlue(img, bestProp.xmin, bestProp.xminy, "");
        markPosBlue(img, bestProp.xmax, bestProp.xmaxy, "");
        System.out.println("Cover width/height: " +
            (bestProp.xmax - bestProp.xmin) + "/" + (bestProp.ymax - bestProp.ymin));
      } else {
        System.out.println("No SCC found corresponding to cover, just using center of kernel.");
        markPos(img, bestX, bestY, "1");
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
      System.out.println(
          "Failed to find red cover in range: " + xmin + "," + ymin + "," + xmax + "," + ymax +
          " with fraction " + maxFrac + " out of " + requiredFraction);
    }
  }

  static void findPanelLed(int[] arr, BufferedImage img,
      int[][] scc, Map<Integer, SccProperty> sccProp, Map<Integer, SccProperty> sccPropLarge,
      SumImage sumImage,
      ModelSummary model, boolean isLeft) {
    int width = arr[1];
    int height = arr[0];
    long xpos = 0;
    long ypos = 0;
    int count = 0;
    int coverX, coverY;
    if (isLeft) {
      coverX = model.objects[1].leftX;
      coverY = model.objects[1].leftY;
    } else {
      coverX = model.objects[1].rightX;
      coverY = model.objects[1].rightY;
    }
    for (int y = Math.max(coverY - height / 20, 0); y < coverY + height / 15 && y < height; ++y) {
      for (int x = coverX; x < coverX + width / 15 && x < width; ++x) {
        int cur = y * width + x + 2;
        int r = (arr[cur] >> 16) & 0x00ff;
        int g = (arr[cur] >> 8) & 0x00ff;
        int b = (arr[cur]) & 0x00ff;
        //if ((g > 0 && Math.max(r, b) == 0) || g > 150 && Math.max(r, b) < 180) {
        if (g > 150 && g > Math.max(r, b) && g - Math.max(r, b) > 30) {
        //if (sumImage.isForeground[x][y] && g > Math.max(r, b)) {
          xpos += x;
          ypos += y;
          img.setRGB(x, y, Color.blue.getRGB());
          ++count;
        }
      }
    }
    if (count > (width < 1700 ? 10 : 20)) {
      System.out.println("Found LED with " + count + " green pixels");
      xpos /= count;
      ypos /= count;
      markPos(img, (int) xpos, (int) ypos, "2");
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
      System.out.println("Failed to find LED.");
    }
  }

  static void findRockerSwitch(int[] arr, BufferedImage img, int[][] scc,
      Map<Integer, SccProperty> sccProp, Map<Integer, SccProperty> sccPropLarge,
      SumImage sumImage,
      ModelSummary model, boolean isLeft) {
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
    int ymin = Math.min(coverY + height / 15, height);
    int ymax = Math.min(coverY + height / 3, height);
    int xmin = Math.max(coverX - width / 20, 0);
    int xmax = Math.min(coverX + width / 5, width);
    int maxSize = 0;
    SccProperty maxProp = null;
    for (int y = ymin; y < ymax; ++y) {
      if (sumImage.isBackground[xmin][y] && sccPropLarge.containsKey(scc[xmin][y])) {
        break;
      }
      for (int x = xmin; x < xmax; ++x) {
        if (sumImage.isBackground[x][y] && sccPropLarge.containsKey(scc[x][y])) {
          break;
        }
        if (sccProp.containsKey(scc[x][y])) {
          SccProperty prop = sccProp.get(scc[x][y]);
          if (prop.size > maxSize &&
              mightBeRocker(img, width, height, scc, coverY + height / 10, prop)) {
            maxSize = prop.size;
            maxProp = prop;
          }
        }
        if (x % 4 == 0 && y % 4 == 0) {
          img.setRGB(x, y, Color.green.getRGB());
        }
      }
    }
    if (maxProp != null) {
      System.out.println("Found Rocker!");
      markPos(img, maxProp.x, maxProp.y, "3");
      markPosRed(img, (maxProp.xmin + maxProp.xmax) / 2, (maxProp.ymin + maxProp.ymax) / 2, "3");

      if (isLeft) {
        model.objects[3].leftX = maxProp.x;
        model.objects[3].leftY = maxProp.y;
        model.objects[3].leftFound = true;
        model.objects[3].leftProp = maxProp;
      } else {
        model.objects[3].rightX = maxProp.x;
        model.objects[3].rightY = maxProp.y;
        model.objects[3].rightFound = true;
        model.objects[3].rightProp = maxProp;
      }
    } else {
      System.out.println("Failed to find rocker.");
    }
  }

  static void estimateGrid(double[][] grid, int coverX, int coverY, int ledX, int ledY,
      int rockerX, int rockerY, BufferedImage img) {
    int width = img.getWidth();
    int height = img.getHeight();
    double a1 = coverX - ledX;
    double b1 = coverY - ledY;
    double a2 = rockerX - ledX;
    double b2 = rockerY - ledY;
    for (int i = 0; i < 9; ++i) {
      int a3 = (int) (grid[i][0] * a1 + grid[i][1] * a2);
      int b3 = (int) (grid[i][0] * b1 + grid[i][1] * b2);
      int x = a3 + ledX;
      int y = b3 + ledY;
      if (x < 0) x = 0;
      if (x >= width) x = width - 1;
      if (y < 0) y = 0;
      if (y >= height) y = height - 1;
      markPos(img, x, y, "null");
    }
  }
  static void estimateGrid(double[][] grid, BufferedImage img, ModelSummary model, boolean isLeft) {
    int width = img.getWidth();
    int height = img.getHeight();
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
      markPos(img, x, y, "null");
    }
  }

  // matchup to *index* of possibleComps.
  // matchup.length == leds.size().
  // used.length == possibleComps.length.
  static class DiffVect {
    int x, y;
    double length;
    public DiffVect(int x, int y, double length) {
      this.x = x;
      this.y = y;
      this.length = length;
    }
  }
  static long dfsCount = 0;
  static boolean dfsMatchStrict = false;
  static SccProperty panelLedForDfs = null;
  static SccProperty rockerForDfs = null;
  static void dfsMatchRect(int[] matchup, int idx, int[] possibleComps, boolean[] used,
      ArrayList<SccProperty> leds,
      int ledX, int ledY,
      double a1, double b1, double a2, double b2,
      int[] bestMatchup, double[] bestDist, double[][] compRatiosOrtho) {
    ++dfsCount;

    // This works because [3] is (1.0, 0).
    int rockerX = (int) (ledX + a1);
    int rockerY = (int) (ledY + b1);
    if (idx == matchup.length) {
      // All matched up. Compute sum of dist squareds.
      // Stash away the difference vectors.
      ArrayList<DiffVect> diffVectors = new ArrayList<DiffVect>();
      double distSum = 0;
      double minScale = 999;
      double maxScale = -999;
      double minTheta = 999;
      double maxTheta = -999;
     // ArrayList<String> tmpReport = new ArrayList<String>();
      for (int i = 0; i < matchup.length; ++i) {
        SccProperty prop = leds.get(i);

        int comp = possibleComps[matchup[i]];
        int expX = (int) (compRatiosOrtho[comp][0] * a1 + compRatiosOrtho[comp][1] * a2);
        int expY = (int) (compRatiosOrtho[comp][0] * b1 + compRatiosOrtho[comp][1] * b2);
        expX += ledX;
        expY += ledY;
        int dx2 = expX - prop.x;
        int dy2 = expY - prop.y;
        double dist = dx2 * dx2 + dy2 * dy2;
        distSum += dist;
        double length = Math.sqrt(dist);
        if (length > 10) {
          diffVectors.add(new DiffVect(dx2, dy2, length));
        }

        expX -= rockerX;
        expY -= rockerY;
        double expLength = Math.sqrt(expX * expX + expY * expY);
        int actX = prop.x - rockerX;
        int actY = prop.y - rockerY;
        double actLength = Math.sqrt(actX * actX + actY * actY);
        double scaleConv = actLength / expLength;
        double crossProduct = expX * actY - expY * actX;
        crossProduct /= actLength;
        crossProduct /= expLength;
        double thetaConv = Math.asin(crossProduct);
        //tmpReport.add("Comp: " + comp + ": scaleConv: " + scaleConv + " thetaConv: " + thetaConv);
        if (scaleConv < minScale) minScale = scaleConv;
        if (scaleConv > maxScale) maxScale = scaleConv;
        if (thetaConv < minTheta) minTheta = thetaConv;
        if (thetaConv > maxTheta) maxTheta = thetaConv;
      }
      distSum /= matchup.length;
      if (diffVectors.size() > 1) {
        distSum *= (maxScale - minScale) + 1;
        distSum *= (maxTheta - minTheta) + 1;
        double pairwiseAngleSum = 0;
        int numSum = 0;
        for (int i = 0; i < diffVectors.size(); ++i) {
          DiffVect v1 = diffVectors.get(i);
          for (int j = i + 1; j < diffVectors.size(); ++j) {
            DiffVect v2 = diffVectors.get(j);
            double dotProduct = v1.x * v2.x + v1.y * v2.y;
            dotProduct /= v1.length;
            dotProduct /= v2.length;
            pairwiseAngleSum += Math.abs(1.0 - dotProduct) + .5;
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
      if (bestDist[0] < 0 || distSum < bestDist[0]) {
        bestDist[0] = distSum;
        //System.out.println("Got better matchup of dist: " + distSum);
        for (int i = 0; i < matchup.length; ++i) {
          //System.out.println(leds.get(i).x + "," + leds.get(i).y + " -> " +
          //    possibleComps[matchup[i]]);
          bestMatchup[i] = matchup[i];
        }
        /*System.out.println("maxScale: " + maxScale);
        System.out.println("minScale: " + minScale);
        System.out.println("maxTheta: " + maxTheta);
        System.out.println("minTheta: " + minTheta);
        for (String tmp : tmpReport) {
          System.out.println(tmp);
        }
        System.out.println();*/
      }
      return;
    }
    for (int i = 0; i < used.length; ++i) {
      if (used[i]) continue;
      boolean slotIsGreen = true;
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
          !leds.get(idx).containsCenter) continue;*/
      if (leds.get(idx).isGreenLed && !slotIsGreen) continue;
      if (leds.get(idx).isBlueLed && slotIsGreen) continue;
      int[] topoY = {0, 0, 0, 2, 1, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 5, 4, 5, 4, 6, 5, 4};
      int[] topoX = {3, 3, 3, 3, 3, 3, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 0, 2, 2, 2, 3, 3};

      // Look at the current matchup state and eliminate impossible branches.
      //if (leds.get(idx).isGreenLed || leds.get(idx).isBlueLed) {
        {
        SccProperty[] curMatchup = new SccProperty[22];
        for (int j = 0; j < idx; ++j) {
          curMatchup[possibleComps[matchup[j]]] = leds.get(j);
        }
        int ymin = 0;
        int ymax = Integer.MAX_VALUE;
        int xmin = 0;
        int xmax = Integer.MAX_VALUE;
        if (panelLedForDfs != null) {
          ymin = panelLedForDfs.y;
        }
        if (rockerForDfs != null) {
          xmax = rockerForDfs.xmax + 200; 
        }
        int topoClassY = topoY[possibleComps[i]];
        int topoClassX = topoX[possibleComps[i]];
        for (int j = 0; j < 22; ++j) {
          if (curMatchup[j] != null) {
            if (topoY[j] < topoClassY) {
              if (dfsMatchStrict) {
                ymin = Math.max(ymin, curMatchup[j].ymax);
              } else {
                ymin = Math.max(ymin, curMatchup[j].y);
              }
            } else if (topoY[j] > topoClassY) {
              if (dfsMatchStrict) {
                ymax = Math.min(ymax, curMatchup[j].ymin);
              } else {
                ymax = Math.min(ymax, curMatchup[j].y);
              }
            }
            if (topoX[j] < topoClassX) {
              if (dfsMatchStrict) {
                xmin = Math.max(xmin, curMatchup[j].xmax);
              } else {
                xmin = Math.max(xmin, curMatchup[j].x);
              }
            } else if (topoX[j] > topoClassX) {
              if (dfsMatchStrict) {
                xmax = Math.min(xmax, curMatchup[j].xmin);
              } else {
                xmax = Math.min(xmax, curMatchup[j].x);
              }
            }
          }
        }
        if (leds.get(idx).y < ymin) {
          continue;
        }
        if (leds.get(idx).y >= ymax) {
          continue;
        }
        if (leds.get(idx).x < xmin) {
          continue;
        }
        if (leds.get(idx).x >= xmax) {
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

  static boolean matchSccToBothLedsDfsRect(int[] arr, BufferedImage img, int[][] scc,
      Map<Integer, SccProperty> sccProp, Map<Integer, SccProperty> sccPropLarge,
      ModelSummary model, boolean isLeft) {
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
      double sinVal = dx / Math.sqrt(length);
      //double theta = Math.asin(sinVal) * 2;
      double theta = Math.asin(sinVal);
      double a3 = Math.cos(theta) * a2 - Math.sin(theta) * b2;
      double b3 = Math.sin(theta) * a2 + Math.cos(theta) * b2;
      a2 = a3;
      b2 = b3;
    }
    ArrayList<SccProperty> blueLeds = new ArrayList<SccProperty>();
    {
      int blueLedsCount = 0;
      for (SccProperty prop : sccProp.values()) {
        if (!prop.isBlueLed) continue;
        if (blueLeds.contains(prop)) continue;
        if (!canReach(img, width, height, prop.x, prop.y, ledX, ledY, scc, sccPropLarge)) continue;
        // Don't even re-check the blue LED's comp.
        int compX = prop.x;
        int compY = prop.y;
        int dx2 = compX - ledX;
        int dy2 = compY - ledY;
        double length2 = Math.sqrt(dx2 * dx2 + dy2 * dy2);
        if (length2 < 10) {
          System.out.println("Skipping a component of length2: " + length2);
          continue;
        }
        ++blueLedsCount;
        // Assuming panelLED is correct, throw out obviously wrong ones.
        if (prop.x > ledX) continue;

        blueLeds.add(prop);
      }
      if (blueLedsCount > 12) {
        System.out.println("Way too many blueLeds: " + blueLeds.size() + " bailing out.");
        blueLeds.clear();
      } else if (blueLeds.size() > 6) {
        System.out.println("Error, got more blueLeds than possible: " + blueLeds.size() + 
            " trimming...");
        for (int i = blueLeds.size() - 1; i >= 6; --i) {
          blueLeds.remove(i);
        }
      }
    }
    ArrayList<SccProperty> greenLeds = new ArrayList<SccProperty>();
    {
      int greenLedsCount = 0;
      for (SccProperty prop : sccProp.values()) {
        if (!prop.isGreenLed) continue;
        if (greenLeds.contains(prop)) continue;
        if (!canReach(img, width, height, prop.x, prop.y, ledX, ledY, scc, sccPropLarge)) continue;
        // Don't even re-check the green LED's comp.
        int compX = prop.x;
        int compY = prop.y;
        int dx2 = compX - ledX;
        int dy2 = compY - ledY;
        double length2 = Math.sqrt(dx2 * dx2 + dy2 * dy2);
        if (length2 < 10) {
          System.out.println("Skipping a component of length2: " + length2);
          continue;
        }
        ++greenLedsCount;

        greenLeds.add(prop);
      }
      if (greenLeds.size() > 8) {
        System.out.println("Error, got more greenLeds than possible: " + greenLeds.size() + 
            " trimming...");
        for (int i = greenLeds.size() - 1; i >= 8; --i) {
          greenLeds.remove(i);
        }
      }
    }
    int[] possibleComps = { 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 19, 21 };
    ArrayList<SccProperty> allLeds = new ArrayList<SccProperty>();
    allLeds.addAll(blueLeds);
    allLeds.addAll(greenLeds);
    int[] matchup = new int[allLeds.size()];
    boolean[] used = new boolean[possibleComps.length];
    int[] bestMatchup = new int[matchup.length];
    double[] bestDist = new double[1];
    bestDist[0] = -1;
    dfsCount = 0;
    dfsMatchRect(matchup, 0, possibleComps, used, allLeds, ledX, ledY, a1, b1, a2, b2,
        bestMatchup, bestDist, isLeft ? compRatiosOrthoLeft : compRatiosOrthoRight);
    System.out.println("DFS took " + dfsCount + " iterations");
    if (bestDist[0] < 0) {
      System.out.println("MatchSccToBoth found no valid matchings.");
      return false;
    }
    for (int i = 0; i < bestMatchup.length; ++i) {
      SccProperty prop = allLeds.get(i);
      int matchedComp = possibleComps[bestMatchup[i]];
      System.out.println("Matched prop " + prop.x + "," + prop.y + " to comp " + matchedComp +
          "(" + prop.r + "/" + prop.g + "/" + prop.b + ")");
      markPosRed(img, prop.x, prop.y, "M" + Integer.toString(matchedComp));
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

  static boolean matchSccToGreenLedsDfsRect(int[] arr, BufferedImage img, int[][] scc,
      Map<Integer, SccProperty> sccProp, Map<Integer, SccProperty> sccPropLarge,
      ModelSummary model, boolean isLeft) {
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

    ArrayList<SccProperty> greenDots = new ArrayList<SccProperty>();
    ArrayList<SccProperty> greenDonuts = new ArrayList<SccProperty>();
    for (SccProperty prop : sccProp.values()) {
      if (!prop.isGreenLed) continue;
      if (greenDots.contains(prop)) continue;
      if (greenDonuts.contains(prop)) continue;
      if (!canReach(img, width, height, prop.x, prop.y, ledX, ledY, scc, sccPropLarge)) continue;
      // Don't even re-check the green LED's comp.
      int compX = prop.x;
      int compY = prop.y;
      int dx2 = compX - ledX;
      int dy2 = compY - ledY;
      double length2 = dx2 * dx2 + dy2 * dy2;
      if (length2 < 100) {
        System.out.println("Skipping a component of length2: " + length2);
        continue;
      }
      int propCircleArea = (prop.xmax - prop.xmin) * (prop.ymax - prop.ymin) * 3 / 4;
      if (propCircleArea == 0) propCircleArea = 1;
      double fillRatio = prop.size / (double) (propCircleArea);
      System.out.println("Fill ratio: " + fillRatio);
      if (fillRatio < 0.9 &&
          (sccProp.get(scc[prop.x][prop.y]) == null || sccProp.get(scc[prop.x][prop.y]) != prop)) {
        prop.containsCenter = false;
      }
//      else {
//        greenDots.add(prop);
//      }
      greenDonuts.add(prop);
    }
    {
      int[] possibleCompsDonut = { 4, 7, 9, 11, 13, 16, 18, 21 };
      //int[] possibleCompsDonut = { 7, 9, 11, 13 };
      if (greenDonuts.size() > possibleCompsDonut.length) {
        System.out.println("Error, got more greenDonuts than possible: " + greenDonuts.size() + 
            " trimming...");
        for (int i = greenDonuts.size() - 1; i >= possibleCompsDonut.length; --i) {
          greenDonuts.remove(i);
        }
      }

      int[] matchup = new int[greenDonuts.size()];
      boolean[] used = new boolean[possibleCompsDonut.length];
      int[] bestMatchup = new int[matchup.length];
      double[] bestDist = new double[1];
      bestDist[0] = -1;
      //for (double scaleMult = .9; scaleMult < 1.5; scaleMult += .1) {
        //for (double thetaMult = 1; thetaMult < 2; thetaMult += .1) {
          {{
          double thetaMult = 1;
          double scaleMult = 1;
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
            double sinVal = dx / Math.sqrt(length);
            double theta = Math.asin(sinVal) * thetaMult;
            double a3 = Math.cos(theta) * a2 - Math.sin(theta) * b2;
            double b3 = Math.sin(theta) * a2 + Math.cos(theta) * b2;
            a2 = a3;
            b2 = b3;
            a2 *= scaleMult;
            b2 *= scaleMult;
          }
          dfsCount = 0;
          dfsMatchRect(matchup, 0, possibleCompsDonut, used, greenDonuts, ledX, ledY, a1, b1, a2, b2,
              bestMatchup, bestDist, isLeft ? compRatiosOrthoLeft : compRatiosOrthoRight);
          System.out.println("DFS took " + dfsCount + " iterations");
        }
      }
      if (bestDist[0] < 0) return false;
      for (int i = 0; i < bestMatchup.length; ++i) {
        SccProperty prop = greenDonuts.get(i);
        int matchedComp = possibleCompsDonut[bestMatchup[i]];
        System.out.println("Matched prop " + prop.x + "," + prop.y + " to comp " + matchedComp);
        markPosRed(img, prop.x, prop.y, "M" + Integer.toString(matchedComp));
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
    return true;
//    {
//      int[] possibleCompsDot = { 4, 16, 18, 21 };
//      if (greenDots.size() > possibleCompsDot.length) {
//        System.out.println("Error, got more greenDots than possible: " + greenDots.size() + 
//            " trimming...");
//        for (int i = greenDots.size() - 1; i >= possibleCompsDot.length; --i) {
//          greenDots.remove(i);
//        }
//      }
//      int[] matchup = new int[greenDots.size()];
//      boolean[] used = new boolean[possibleCompsDot.length];
//      int[] bestMatchup = new int[matchup.length];
//      double[] bestDist = new double[1];
//      bestDist[0] = -1;
//      dfsMatchRect(matchup, 0, possibleCompsDot, used, greenDots, ledX, ledY, a1, b1, a2, b2,
//          bestMatchup, bestDist, isLeft ? compRatiosOrthoLeft : compRatiosOrthoRight);
//      for (int i = 0; i < bestMatchup.length; ++i) {
//        SccProperty prop = greenDots.get(i);
//        int matchedComp = possibleCompsDot[bestMatchup[i]];
//        System.out.println("Matched prop " + prop.x + "," + prop.y + " to comp " + matchedComp);
//        markPosRed(img, prop.x, prop.y, "M" + Integer.toString(matchedComp));
//        if (isLeft) {
//          model.objects[matchedComp].leftX = prop.x;
//          model.objects[matchedComp].leftY = prop.y;
//          model.objects[matchedComp].leftFound = true;
//          model.objects[matchedComp].state = "ON";
//        } else {
//          model.objects[matchedComp].rightX = prop.x;
//          model.objects[matchedComp].rightY = prop.y;
//          model.objects[matchedComp].rightFound = true;
//          model.objects[matchedComp].state = "ON";
//        }
//      }
//    }
  }

  static void matchSccToBlueLedsDfsRect(int[] arr, BufferedImage img,
      int[][] scc, Map<Integer, SccProperty> sccProp, Map<Integer, SccProperty> sccPropLarge,
      ModelSummary model, boolean isLeft) {
    int[] possibleComps = { 6, 8, 10, 12, 14, 19 };
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
      double sinVal = dx / Math.sqrt(length);
      //double theta = Math.asin(sinVal) * 2;
      double theta = Math.asin(sinVal);
      double a3 = Math.cos(theta) * a2 - Math.sin(theta) * b2;
      double b3 = Math.sin(theta) * a2 + Math.cos(theta) * b2;
      a2 = a3;
      b2 = b3;
      // Draw all the predicted positions.
      //allowMarkPos = true;
      /*for (int i = 0; i < 22; ++i) {
      // int[] pointsOfInterest = {0, 15, 17, 20 };
      //for (int i : pointsOfInterest) {
        int comp = i;
//        {
//          theta = Math.asin(sinVal) * 2;
//          a3 = Math.cos(theta) * a2 - Math.sin(theta) * b2;
//          b3 = Math.sin(theta) * a2 + Math.cos(theta) * b2;
//          int newX = (int) (compRatiosOrthoOld[i][0] * a1 + compRatiosOrthoOld[i][1] * a3);
//          int newY = (int) (compRatiosOrthoOld[i][0] * b1 + compRatiosOrthoOld[i][1] * b3);
//          newX += ledX;
//          newY += ledY;
//          markPosRed(img, newX, newY, "" + Integer.toString(comp));
//        }
        if (isLeft) {
          int newX = (int) (compRatiosOrthoLeft[i][0] * a1 + compRatiosOrthoLeft[i][1] * a2);
          int newY = (int) (compRatiosOrthoLeft[i][0] * b1 + compRatiosOrthoLeft[i][1] * b2);
          newX += ledX;
          newY += ledY;
          markPos(img, newX, newY, "" + Integer.toString(comp));
        } else {
          int newX = (int) (compRatiosOrthoRight[i][0] * a1 + compRatiosOrthoRight[i][1] * a2);
          int newY = (int) (compRatiosOrthoRight[i][0] * b1 + compRatiosOrthoRight[i][1] * b2);
          newX += ledX;
          newY += ledY;
          markPos(img, newX, newY, "" + Integer.toString(comp));
        }
      }*/
      //allowMarkPos = false;
    }

    ArrayList<SccProperty> blueLeds = new ArrayList<SccProperty>();
    int blueLedsCount = 0;
    for (SccProperty prop : sccProp.values()) {
      if (!prop.isBlueLed) continue;
      if (blueLeds.contains(prop)) continue;
      if (!canReach(img, width, height, prop.x, prop.y, ledX, ledY, scc, sccPropLarge)) continue;
      // Don't even re-check the blue LED's comp.
      int compX = prop.x;
      int compY = prop.y;
      int dx2 = compX - ledX;
      int dy2 = compY - ledY;
      double length2 = Math.sqrt(dx2 * dx2 + dy2 * dy2);
      if (length2 < 10) {
        System.out.println("Skipping a component of length2: " + length2);
        continue;
      }
      ++blueLedsCount;
      // Assuming panelLED is correct, throw out obviously wrong ones.
      if (prop.x > ledX) continue;

      blueLeds.add(prop);
      markPosRed(img, prop.x, prop.y, "");
    }
    if (blueLedsCount > 12) {
      System.out.println("Way too many blueLeds: " + blueLeds.size() + " bailing out.");
      return;
    } else if (blueLeds.size() > 6) {
      System.out.println("Error, got more blueLeds than possible: " + blueLeds.size() + 
          " trimming...");
      for (int i = blueLeds.size() - 1; i >= 6; --i) {
        blueLeds.remove(i);
      }
    }
    int[] matchup = new int[blueLeds.size()];
    boolean[] used = new boolean[6];
    int[] bestMatchup = new int[matchup.length];
    double[] bestDist = new double[1];
    bestDist[0] = -1;
    dfsMatchRect(matchup, 0, possibleComps, used, blueLeds, ledX, ledY, a1, b1, a2, b2,
        bestMatchup, bestDist, isLeft ? compRatiosOrthoLeft : compRatiosOrthoRight);
    if (bestDist[0] < 0) return;
    for (int i = 0; i < bestMatchup.length; ++i) {
      SccProperty prop = blueLeds.get(i);
      int matchedComp = possibleComps[bestMatchup[i]];
      System.out.println("Matched prop " + prop.x + "," + prop.y + " to comp " + matchedComp +
          "(" + prop.r + "/" + prop.g + "/" + prop.b + ")");
      markPosRed(img, prop.x, prop.y, "M" + Integer.toString(matchedComp));
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

  static boolean isOverlap(SccProperty prop0, SccProperty prop1) {
    if (prop0.xmin >= prop1.xmin && prop0.xmin <= prop1.xmax &&
        prop0.ymin >= prop1.ymin && prop0.ymin <= prop1.ymax) {
      return true;
    }
    if (prop0.xmin >= prop1.xmin && prop0.xmin <= prop1.xmax &&
        prop0.ymax >= prop1.ymin && prop0.ymax <= prop1.ymax) {
      return true;
    }
    if (prop0.xmax >= prop1.xmin && prop0.xmax <= prop1.xmax &&
        prop0.ymin >= prop1.ymin && prop0.ymin <= prop1.ymax) {
      return true;
    }
    if (prop0.xmax >= prop1.xmin && prop0.xmax <= prop1.xmax &&
        prop0.ymax >= prop1.ymin && prop0.ymax <= prop1.ymax) {
      return true;
    }
    if (prop1.xmin >= prop0.xmin && prop1.xmin <= prop0.xmax &&
        prop1.ymin >= prop0.ymin && prop1.ymin <= prop0.ymax) {
      return true;
    }
    if (prop1.xmin >= prop0.xmin && prop1.xmin <= prop0.xmax &&
        prop1.ymax >= prop0.ymin && prop1.ymax <= prop0.ymax) {
      return true;
    }
    if (prop1.xmax >= prop0.xmin && prop1.xmax <= prop0.xmax &&
        prop1.ymin >= prop0.ymin && prop1.ymin <= prop0.ymax) {
      return true;
    }
    if (prop1.xmax >= prop0.xmin && prop1.xmax <= prop0.xmax &&
        prop1.ymax >= prop0.ymin && prop1.ymax <= prop0.ymax) {
      return true;
    }
    return false;
  }

  static void doBasisDfsMatch(int[] arr, BufferedImage img, int[][] scc,
      Map<Integer, SccProperty> sccProp, Map<Integer, SccProperty> sccPropLarge,
      ModelSummary model, boolean isLeft,
      double a1, double b1, double a2, double b2, double[][] compRatios) {
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

    ArrayList<SccProperty> allLeds = new ArrayList<SccProperty>();
    for (SccProperty prop : sccProp.values()) {
      if (prop.isGreenLed || prop.isBlueLed) {
        allLeds.add(prop);
      }
    }

    // For each not-found grid component.
    double maxDistFromPredicted = height / 27;
    maxDistFromPredicted *= maxDistFromPredicted;
    long maxDistForMerge = (width + height) / 2 * 25L / 1000L;
    maxDistForMerge *= maxDistForMerge;
    HashSet<Integer> disallowedScc = new HashSet<Integer>();  // Used by other already-found ones.
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
      for (SccProperty prop : sccProp.values()) {
        if (prop.isGreenLed || prop.isBlueLed) continue;
        int dx = prop.x - takenX;
        int dy = prop.y - takenY;
        double dist = dx * dx + dy * dy;
        if (dist < maxDistFromPredicted) {
          disallowedScc.add(prop.propNum);
        }
      }
    }
    ArrayList<SccProperty> allGridScc = new ArrayList<SccProperty>();
    ArrayList<Integer> componentsMissing = new ArrayList<Integer>();
    for (int i = 6; i <= 14; ++i) {
      if (isLeft && model.objects[i].leftFound) continue;
      if (!isLeft && model.objects[i].rightFound) continue;
      componentsMissing.add(i);

      // Find its predicted position.
      int newX = (int) (compRatios[i][0] * a1 + compRatios[i][1] * a2);
      int newY = (int) (compRatios[i][0] * b1 + compRatios[i][1] * b2);
      newX += ledX;
      newY += ledY;

      // Then find all SCC in its vicinity.
      for (SccProperty prop : sccProp.values()) {
        if (prop.isGreenLed || prop.isBlueLed) continue;
        if (disallowedScc.contains(prop.propNum)) continue;
        //if (!canReach(img, width, height, prop.x, prop.y, ledX, ledY, scc, sccPropLarge)) continue;
        int dx = prop.x - newX;
        int dy = prop.y - newY;
        double dist = dx * dx + dy * dy;
        if (dist < maxDistFromPredicted) {
          // As long as it's not just part of another LED.
          boolean isEdgeOfLed = false;
          for (SccProperty led : allLeds) {
            dx = prop.x - led.x;
            dy = prop.y - led.y;
            if (dx * dx + dy * dy < maxDistForMerge) {
              System.out.println("Skipping SCC " + prop.x + "," + prop.y + " due to proximity " +
                  "to led: " + led.x + "," + led.y);
              isEdgeOfLed = true;
              break;
            }
          }
          if (!isEdgeOfLed) {
            allGridScc.add(prop);
          }
        }
      }
    }
    if (componentsMissing.size() == 0) return;
    if (allGridScc.size() == 0) return;

    // Merge gridScc which are close to each other.
    ArrayList<SccProperty> mergedGrid = new ArrayList<SccProperty>();
    boolean hadMerge = true;
    while (hadMerge) {
      HashSet<Integer> alreadyMergedGrid = new HashSet<Integer>();
      hadMerge = false;
      for (int i = 0; i < allGridScc.size(); ++i) {
        SccProperty prop0 = allGridScc.get(i);
        if (alreadyMergedGrid.contains(prop0.propNum)) continue;
        ArrayList<Integer> mergedComponents = new ArrayList<Integer>();
        alreadyMergedGrid.add(prop0.propNum);
        mergedComponents.add(prop0.propNum);
        for (int j = i + 1; j < allGridScc.size(); ++j) {
          SccProperty prop1 = allGridScc.get(j);
          if (alreadyMergedGrid.contains(prop1.propNum)) continue;
          int dx = prop0.x - prop1.x;
          int dy = prop0.y - prop1.y;
          if (dx * dx + dy * dy < maxDistForMerge) {
            Scc.mergeScc(prop0, prop1);
            alreadyMergedGrid.add(prop1.propNum);
            mergedComponents.add(prop1.propNum);
            hadMerge = true;
          }
        }
        for (int k = 0; k < mergedComponents.size(); ++k) {
          sccProp.put(mergedComponents.get(k), prop0);
        }
        mergedGrid.add(prop0);
      }
      if (hadMerge) {
        allGridScc = mergedGrid;
        mergedGrid = new ArrayList<SccProperty>();
      }
    }

    for (SccProperty prop : mergedGrid) {
      markPosRed(img, prop.x, prop.y, "");
      Graphics g = img.getGraphics();
      g.setColor(Color.red);
      for (int j = 0; j < 5; ++j) {
        g.drawRect(prop.xmin - j, prop.ymin - j,
            prop.xmax - prop.xmin + 2 * j, prop.ymax - prop.ymin + 2 * j);
      }
    }

    if (mergedGrid.size() > componentsMissing.size()) {
      System.out.println("Found " + mergedGrid.size() + " SCC for " + componentsMissing.size() +
          " possible components; bailing out.");
      for (int i = mergedGrid.size() - 1; i >= componentsMissing.size(); --i) {
        mergedGrid.remove(i);
      }
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
      double sinVal = dx / Math.sqrt(length);
      //double theta = Math.asin(sinVal) * 2;
      double theta = Math.asin(sinVal);
      double a3 = Math.cos(theta) * a2 - Math.sin(theta) * b2;
      double b3 = Math.sin(theta) * a2 + Math.cos(theta) * b2;
      a2 = a3;
      b2 = b3;
    }

    int[] possibleComps = new int[componentsMissing.size()];
    for (int i = 0; i < possibleComps.length; ++i) {
      possibleComps[i] = componentsMissing.get(i);
    }
    int[] matchup = new int[mergedGrid.size()];
    boolean[] used = new boolean[possibleComps.length];
    int[] bestMatchup = new int[matchup.length];
    double[] bestDist = new double[1];
    bestDist[0] = -1;
    dfsMatchRect(matchup, 0, possibleComps, used, mergedGrid, ledX, ledY, a1, b1, a2, b2,
        bestMatchup, bestDist, isLeft ? compRatiosOrthoLeft : compRatiosOrthoRight);
    if (bestDist[0] < 0) {
      return;
    }
    for (int i = 0; i < bestMatchup.length; ++i) {
      SccProperty prop = mergedGrid.get(i);
      int matchedComp = possibleComps[bestMatchup[i]];
      System.out.println("Matched prop " + prop.x + "," + prop.y + " to comp " + matchedComp +
          "(" + prop.r + "/" + prop.g + "/" + prop.b + ")");
      markPosRed(img, prop.x, prop.y, "D" + Integer.toString(matchedComp));
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


  static void basisDfsMatch(int[] arr, BufferedImage img, int[][] scc,
      Map<Integer, SccProperty> sccProp, Map<Integer, SccProperty> sccPropLarge,
      ModelSummary model, boolean isLeft) {
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
      gridX = model.objects[10].leftX;
      gridY = model.objects[10].leftY;

    } else {
      ledX = model.objects[2].rightX;
      ledY = model.objects[2].rightY;
      rockerX = model.objects[3].rightX;
      rockerY = model.objects[3].rightY;
      gridX = model.objects[10].rightX;
      gridY = model.objects[10].rightY;
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
      /*for (int i = 0; i < 22; ++i) {
      // int[] pointsOfInterest = {0, 15, 17, 20 };
      //for (int i : pointsOfInterest) {
        int comp = i;
        {
          int newX = (int) (compRatios10[i][0] * a1 + compRatios10[i][1] * a2);
          int newY = (int) (compRatios10[i][0] * b1 + compRatios10[i][1] * b2);
          newX += ledX;
          newY += ledY;
          markPos(img, newX, newY, "" + Integer.toString(comp));
        }
      }*/
    doBasisDfsMatch(arr, img, scc, sccProp, sccPropLarge, model, isLeft, a1, b1, a2, b2, compRatios10);
  }

  static void basisDfsMatchWeak(int[] arr, BufferedImage img, int[][] scc,
      Map<Integer, SccProperty> sccProp, Map<Integer, SccProperty> sccPropLarge,
      ModelSummary model, boolean isLeft) {
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
    System.out.println("Doing basisExtrapolateWeak with " + ledX + "," + ledY + " " +
        rockerY + "," + rockerY);
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
      double sinVal = dx / Math.sqrt(length);
      //double theta = Math.asin(sinVal) * 2;
      double theta = Math.asin(sinVal);
      double a3 = Math.cos(theta) * a2 - Math.sin(theta) * b2;
      double b3 = Math.sin(theta) * a2 + Math.cos(theta) * b2;
      a2 = a3;
      b2 = b3;
    }

    // TODO: maybe include 20 and 21?
    double avgScale = 0;
    double avgTheta = 0;
    int numSamples = 0;
    ArrayList<Double> scaleSamples = new ArrayList<Double>();
    ArrayList<Double> thetaSamples = new ArrayList<Double>();
    // For each component, find the scale and rotation correction relative to the vector
    // made by that component to the rocker as compared to the predicted position, and
    // average the values in.
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
      double predLength = Math.sqrt(predX * predX + predY * predY);
      double actLength = Math.sqrt(actX * actX + actY * actY);
      double scaleConv = actLength / predLength;
      double crossProduct = predX * actY - predY * actX;
      crossProduct /= actLength;
      crossProduct /= predLength;
      double thetaConv = Math.asin(crossProduct);
      avgScale += scaleConv - 1;
      avgTheta += thetaConv;
      ++numSamples;
      scaleSamples.add(scaleConv - 1);
      thetaSamples.add(thetaConv);
    }
    /*if (numSamples > 0) {
      avgScale = (scaleSamples.get((numSamples - 1) / 2) + scaleSamples.get(numSamples / 2)) / 2;
      avgTheta = (thetaSamples.get((numSamples - 1) / 2) + thetaSamples.get(numSamples / 2)) / 2;
      avgScale += 1;
      // Correct the normal vector by the average correction factors.
      double a3 = Math.cos(avgTheta) * a2 - Math.sin(avgTheta) * b2;
      double b3 = Math.sin(avgTheta) * a2 + Math.cos(avgTheta) * b2;
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
      double a3 = Math.cos(avgTheta) * a2 - Math.sin(avgTheta) * b2;
      double b3 = Math.sin(avgTheta) * a2 + Math.cos(avgTheta) * b2;
      a3 *= avgScale;
      b3 *= avgScale;
      a2 = a3;
      b2 = b3;
    }

    doBasisDfsMatch(arr, img, scc, sccProp, sccPropLarge, model, isLeft, a1, b1, a2, b2,
        isLeft ? compRatiosOrthoLeft : compRatiosOrthoRight);
  }

  static boolean[] blockBasisExtrapolate = new boolean[22];
  static void doBasisExtrapolate(int[] arr, BufferedImage img,
      Map<Integer, SccProperty> sccProp,
      ModelSummary model, boolean isLeft,
      double a1, double b1, double a2, double b2, double[][] compRatios) {
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
      SccProperty bestProp = null;
      for (SccProperty prop : sccProp.values()) {
        int dx = prop.x - newX;
        int dy = prop.y - newY;
        double dist = dx * dx + dy * dy;
        if (dist < bestDist) {
          bestProp = prop;
          bestDist = dist;
        }
      }
      if (bestProp == null) {
        if (isLeft) {
          model.objects[i].leftX = newX;
          model.objects[i].leftY = newY;
          model.objects[i].leftFound = true;
          System.out.println("Basis extrapolation matched no SCC, using raw for comp " + i );
          /*markPos(img, newX, newY, "BR" + Integer.toString(i));
          markPosYellow(img, newX, newY, "");*/
        } else {
          model.objects[i].rightX = newX;
          model.objects[i].rightY = newY;
          model.objects[i].rightFound = true;
          System.out.println("Basis extrapolation matched no SCC, using raw for comp " + i );
          /*markPos(img, newX, newY, "BR" + Integer.toString(i));
          markPosYellow(img, newX, newY, "");*/
        }
        continue;
      }
      int origX = newX;
      int origY = newY;
      if (i == 15 || i == 17 || i == 20) {
        // Toggles use topmost point instead of center of mass.
        newX = bestProp.yminx;
        newY = bestProp.ymin;
      } else {
        newX = bestProp.x;
        newY = bestProp.y;
      }
      if (isLeft) {
        model.objects[i].leftX = newX;
        model.objects[i].leftY = newY;
        model.objects[i].leftFound = true;
        System.out.println("Basis extrapolation matched comp " + i + " to scc: " +
            newX + "," + newY);
        markPos(img, newX, newY, "B" + Integer.toString(i));
        //markPosYellow(img, origX, origY, "");
      } else {
        model.objects[i].rightX = newX;
        model.objects[i].rightY = newY;
        model.objects[i].rightFound = true;
        System.out.println("Basis extrapolation matched comp " + i + " to scc: " +
            newX + "," + newY);
        markPos(img, newX, newY, "B" + Integer.toString(i));
        //markPosYellow(img, origX, origY, "");
      }
    }
  }

  static void basisExtrapolate(int[] arr, BufferedImage img,
      Map<Integer, SccProperty> sccProp,
      ModelSummary model, boolean isLeft) {
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
      gridX = model.objects[10].leftX;
      gridY = model.objects[10].leftY;

    } else {
      ledX = model.objects[2].rightX;
      ledY = model.objects[2].rightY;
      rockerX = model.objects[3].rightX;
      rockerY = model.objects[3].rightY;
      gridX = model.objects[10].rightX;
      gridY = model.objects[10].rightY;
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
    doBasisExtrapolate(arr, img, sccProp, model, isLeft, a1, b1, a2, b2, compRatios10);
  }

  static void basisExtrapolateWeak(int[] arr, BufferedImage img,
      Map<Integer, SccProperty> sccProp,
      ModelSummary model, boolean isLeft) {
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
    System.out.println("Doing basisExtrapolateWeak with " + ledX + "," + ledY + " " +
        rockerY + "," + rockerY);
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
      double sinVal = dx / Math.sqrt(length);
      //double theta = Math.asin(sinVal) * 2;
      double theta = Math.asin(sinVal);
      double a3 = Math.cos(theta) * a2 - Math.sin(theta) * b2;
      double b3 = Math.sin(theta) * a2 + Math.cos(theta) * b2;
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
      double predLength = Math.sqrt(predX * predX + predY * predY);
      double actLength = Math.sqrt(actX * actX + actY * actY);
      double scaleConv = actLength / predLength;
      double crossProduct = predX * actY - predY * actX;
      crossProduct /= actLength;
      crossProduct /= predLength;
      double thetaConv = Math.asin(crossProduct);
      avgScale += scaleConv - 1;
      avgTheta += thetaConv;
      ++numSamples;
    }
    if (numSamples > 0) {
      avgScale /= numSamples;
      avgTheta /= numSamples;
      avgScale += 1;
      // Correct the normal vector by the average correction factors.
      double a3 = Math.cos(avgTheta) * a2 - Math.sin(avgTheta) * b2;
      double b3 = Math.sin(avgTheta) * a2 + Math.cos(avgTheta) * b2;
      a3 *= avgScale;
      b3 *= avgScale;
      a2 = a3;
      b2 = b3;
    }

    doBasisExtrapolate(arr, img, sccProp, model, isLeft, a1, b1, a2, b2,
        isLeft ? compRatiosOrthoLeft : compRatiosOrthoRight);
  }

  static boolean gridHasObj(ModelSummary model, boolean isLeft, int obj) {
    if (isLeft) {
      return model.objects[obj].leftFound;
    } else {
      return model.objects[obj].rightFound;
    }
  }

  static void centerObj(BufferedImage img, ModelSummary model, boolean isLeft,
      int target, int objA, int objB) {
    if (isLeft) {
      model.objects[target].leftX = (model.objects[objA].leftX + model.objects[objB].leftX) / 2;
      model.objects[target].leftY = (model.objects[objA].leftY + model.objects[objB].leftY) / 2;
      if (!model.objects[target].state.equals("ON")) {
        model.objects[target].state = "OFF";
      }
      model.objects[target].leftFound = true;
      markPos(img, model.objects[target].leftX, model.objects[target].leftY,
          "F" + Integer.toString(target));
    } else {
      model.objects[target].rightX = (model.objects[objA].rightX + model.objects[objB].rightX) / 2;
      model.objects[target].rightY = (model.objects[objA].rightY + model.objects[objB].rightY) / 2;
      if (!model.objects[target].state.equals("ON")) {
        model.objects[target].state = "OFF";
      }
      model.objects[target].rightFound = true;
      markPos(img, model.objects[target].rightX, model.objects[target].rightY,
          "F" + Integer.toString(target));
    }
  }

  static void extrapolateObj(BufferedImage img, ModelSummary model, boolean isLeft,
      int target, int objA, int objB) {
    if (isLeft) {
      model.objects[target].leftX = (model.objects[objA].leftX +
          2 * (model.objects[objB].leftX - model.objects[objA].leftX));
      model.objects[target].leftY = (model.objects[objA].leftY +
          2 * (model.objects[objB].leftY - model.objects[objA].leftY));
      if (!model.objects[target].state.equals("ON")) {
        model.objects[target].state = "OFF";
      }
      model.objects[target].leftFound = true;
      markPos(img, model.objects[target].leftX, model.objects[target].leftY,
          "E" + Integer.toString(target));
    } else {
      model.objects[target].rightX = (model.objects[objA].rightX +
          2 * (model.objects[objB].rightX - model.objects[objA].rightX));
      model.objects[target].rightY = (model.objects[objA].rightY +
          2 * (model.objects[objB].rightY - model.objects[objA].rightY));
      if (!model.objects[target].state.equals("ON")) {
        model.objects[target].state = "OFF";
      }
      model.objects[target].rightFound = true;
      markPos(img, model.objects[target].rightX, model.objects[target].rightY,
          "E" + Integer.toString(target));
    }
  }

  static void fillGrid(int[] arr, BufferedImage img,
      Map<Integer, SccProperty> sccProp,
      ModelSummary model, boolean isLeft) {
    // 6  7  8
    // 9  10 11
    // 12 13 14
    // Edge centers.
    boolean foundExtrapolation = true;
    while (foundExtrapolation) {
      foundExtrapolation = false;
      for (int target = 9; target <= 11; ++target) {
        if (!gridHasObj(model, isLeft, target) &&
            gridHasObj(model, isLeft, target - 3) && gridHasObj(model, isLeft, target + 3)) {
          centerObj(img, model, isLeft, target, target - 3, target + 3);
          foundExtrapolation = true;
        }
      }
      for (int target = 7; target <= 13; target += 3) {
        if (!gridHasObj(model, isLeft, target) &&
            gridHasObj(model, isLeft, target - 1) && gridHasObj(model, isLeft, target + 1)) {
          centerObj(img, model, isLeft, target, target - 1, target + 1);
          foundExtrapolation = true;
        }
      }

      // Dead center.
      if (!gridHasObj(model, isLeft, 10) &&
          gridHasObj(model, isLeft, 6) && gridHasObj(model, isLeft, 14)) {
        centerObj(img, model, isLeft, 10, 6, 14);
        foundExtrapolation = true;
      }
      if (!gridHasObj(model, isLeft, 10) &&
          gridHasObj(model, isLeft, 12) && gridHasObj(model, isLeft, 8)) {
        centerObj(img, model, isLeft, 10, 12, 8);
        foundExtrapolation = true;
      }

      // Straight extrapolates.
      for (int target = 6; target <= 8; ++target) {
        if (!gridHasObj(model, isLeft, target) &&
            gridHasObj(model, isLeft, target + 3) && gridHasObj(model, isLeft, target + 6)) {
          extrapolateObj(img, model, isLeft, target, target + 6, target + 3);
          foundExtrapolation = true;
        }
      }
      for (int target = 12; target <= 14; ++target) {
        if (!gridHasObj(model, isLeft, target) &&
            gridHasObj(model, isLeft, target - 3) && gridHasObj(model, isLeft, target - 6)) {
          extrapolateObj(img, model, isLeft, target, target - 6, target - 3);
          foundExtrapolation = true;
        }
      }
      for (int target = 6; target <= 12; target += 3) {
        if (!gridHasObj(model, isLeft, target) &&
            gridHasObj(model, isLeft, target + 1) && gridHasObj(model, isLeft, target + 2)) {
          extrapolateObj(img, model, isLeft, target, target + 2, target + 1);
          foundExtrapolation = true;
        }
      }
      for (int target = 8; target <= 14; target += 3) {
        if (!gridHasObj(model, isLeft, target) &&
            gridHasObj(model, isLeft, target - 1) && gridHasObj(model, isLeft, target - 2)) {
          extrapolateObj(img, model, isLeft, target, target - 2, target - 1);
          foundExtrapolation = true;
        }
      }

      // Diagonal extrapolation.
      if (!gridHasObj(model, isLeft, 6) &&
          gridHasObj(model, isLeft, 10) && gridHasObj(model, isLeft, 14)) {
        extrapolateObj(img, model, isLeft, 6, 14, 10);
        foundExtrapolation = true;
      }
      if (!gridHasObj(model, isLeft, 8) &&
          gridHasObj(model, isLeft, 10) && gridHasObj(model, isLeft, 12)) {
        extrapolateObj(img, model, isLeft, 8, 12, 10);
        foundExtrapolation = true;
      }
      if (!gridHasObj(model, isLeft, 12) &&
          gridHasObj(model, isLeft, 10) && gridHasObj(model, isLeft, 8)) {
        extrapolateObj(img, model, isLeft, 12, 8, 10);
        foundExtrapolation = true;
      }
      if (!gridHasObj(model, isLeft, 14) &&
          gridHasObj(model, isLeft, 10) && gridHasObj(model, isLeft, 6)) {
        extrapolateObj(img, model, isLeft, 14, 6, 10);
        foundExtrapolation = true;
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
          markPos(img, model.objects[6].leftX, model.objects[6].leftY,
              "C" + Integer.toString(6));
        } else {
          model.objects[6].rightX = model.objects[7].rightX +
              (model.objects[9].rightX - model.objects[10].rightX);
          model.objects[6].rightY = model.objects[9].rightY +
              (model.objects[7].rightY - model.objects[10].rightY);
          if (model.objects[6].state != "ON") {
            model.objects[6].state = "OFF";
          }
          model.objects[6].rightFound = true;
          markPos(img, model.objects[6].rightX, model.objects[6].rightY,
              "C" + Integer.toString(6));
        }
        foundExtrapolation = true;
      }

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
            SccProperty bestProp = null;
            double bestDist = 99999;
            for (SccProperty prop : sccProp.values()) {
              if (newX >= prop.xmin && newX <= prop.xmax &&
                  newY >= prop.ymin && newY <= prop.ymax) {
                double diffX = newX - prop.x;
                double diffY = newY - prop.y;
                if (diffX * diffX + diffY * diffY < bestDist) {
                  bestDist = diffX * diffX + diffY * diffY;
                  bestProp = prop;
                }
              }
            }
            if (bestProp != null) {
              newX = bestProp.x;
              newY = bestProp.y;
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
            SccProperty bestProp = null;
            double bestDist = 99999;
            for (SccProperty prop : sccProp.values()) {
              if (newX >= prop.xmin && newX <= prop.xmax &&
                  newY >= prop.ymin && newY <= prop.ymax) {
                double diffX = newX - prop.x;
                double diffY = newY - prop.y;
                if (diffX * diffX + diffY * diffY < bestDist) {
                  bestDist = diffX * diffX + diffY * diffY;
                  bestProp = prop;
                }
              }
            }
            if (bestProp != null) {
              newX = bestProp.x;
              newY = bestProp.y;
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


      // L extrapolation
//      for (int target = 10; target <= 11; ++target) {
//        if (!gridHasObj(model, isLeft, target) &&
//            gridHasObj(model, isLeft, target - 3) && gridHasObj(model, isLeft, target - 4)) {
//
//        }
//      }
    }

    // Bottom three LEDs.
    if (!gridHasObj(model, isLeft, 16) &&
        gridHasObj(model, isLeft, 18) && gridHasObj(model, isLeft, 21)) {
      extrapolateObj(img, model, isLeft, 16, 21, 18);
    }
    if (!gridHasObj(model, isLeft, 21) &&
        gridHasObj(model, isLeft, 18) && gridHasObj(model, isLeft, 16)) {
      extrapolateObj(img, model, isLeft, 21, 16, 18);
    }
    if (!gridHasObj(model, isLeft, 18) &&
        gridHasObj(model, isLeft, 16) && gridHasObj(model, isLeft, 21)) {
      centerObj(img, model, isLeft, 18, 16, 21);
    }

//    for (int i = 6; i <= 14; ++i) {
//      int objX, objY;
//      if (isLeft) {
//        objX = model.objects[i].leftX;
//        objY = model.objects[i].leftY;
//      } else {
//        objX = model.objects[i].rightX;
//        objY = model.objects[i].rightY;
//      }
//      SccProperty bestProp = null;
//      int bestDist = 400;  // Max distance we can try to adjust the position.
//      for (SccProperty prop : sccProp.values()) {
//        int dx = objX - prop.x;
//        int dy = objY - prop.y;
//        int dist = dx * dx + dy * dy;
//        if (objX >= prop.xmin && objX <= prop.xmax &&
//            objY >= prop.ymin && objY <= prop.ymax &&
//            dist < bestDist) {
//          bestDist = dist;
//          bestProp = prop;
//        }
//      }
//      if (bestProp != null && (bestProp.x != objX || bestProp.y != objY)) {
//        markPosPink(img, bestProp.x, bestProp.y, "");
//        System.out.println("Extrapolated point moved adjusted from " + objX + "," + objY + " to "
//            + bestProp.x + "," + bestProp.y + " for comp " + i);
//        if (isLeft) {
//          model.objects[i].leftX = bestProp.x;
//          model.objects[i].leftY = bestProp.y;
//        } else {
//          model.objects[i].rightX = bestProp.x;
//          model.objects[i].rightY = bestProp.y;
//        }
//      }
//    }
  }

  static void interpolateRockerLeds(
      int[] arr, BufferedImage img, int[][] scc, Map<Integer, SccProperty> sccProp,
      ModelSummary model, boolean isLeft) {
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

    long maxDist = (width + height) / 2 * 25L / 1000L;
    maxDist *= maxDist;
    SccProperty topProp = null;
    long bestDist = maxDist;
    for (SccProperty prop : sccProp.values()) {
      int dx = prop.x - ledTopX;
      int dy = prop.y - ledTopY;
      if (dx * dx + dy * dy < bestDist) {
        bestDist = dx * dx + dy * dy;
        topProp = prop;
      }
    }

    // Just naive position-guessing, so maybe don't mark as found just yet?
    if (isLeft) {
      // Only update if not found/matched-up already.
      if (!model.objects[4].leftFound) {
        if (topProp != null) {
          model.objects[4].leftX = topProp.x;
          model.objects[4].leftY = topProp.y;
          model.objects[4].leftProp = topProp;
          markPosYellow(img, ledTopX, ledTopY, "");
          markPos(img, topProp.x, topProp.y, "MI4");
        } else {
          model.objects[4].leftX = ledTopX;
          model.objects[4].leftY = ledTopY;
          markPos(img, ledTopX, ledTopY, "I4");
        }
        model.objects[4].leftFound = true;
      }
      if (!model.objects[5].leftFound) {
        model.objects[5].leftX = ledBottomX;
        model.objects[5].leftY = ledBottomY;
        model.objects[5].leftFound = true;
        markPos(img, ledBottomX, ledBottomY, "I5");
      }
    } else {
      if (!model.objects[4].rightFound) {
        if (topProp != null) {
          model.objects[4].rightX = topProp.x;
          model.objects[4].rightY = topProp.y;
          model.objects[4].rightProp = topProp;
          markPosYellow(img, ledTopX, ledTopY, "");
          markPos(img, topProp.x, topProp.y, "MI4");
        } else {
          model.objects[4].rightX = ledTopX;
          model.objects[4].rightY = ledTopY;
          markPos(img, ledTopX, ledTopY, "I4");
        }
        model.objects[4].rightFound = true;
      }
      if (!model.objects[5].rightFound) {
        model.objects[5].rightX = ledBottomX;
        model.objects[5].rightY = ledBottomY;
        model.objects[5].rightFound = true;
        markPos(img, ledBottomX, ledBottomY, "I5");
      }
    }
  }

  static void reconcileMissing(BufferedImage imgLeft, BufferedImage imgRight, ModelSummary model) {
    // Extract all the ones present in both.
    int dx = 0;
    int dy = 0;
    int count = 0;
    // For now, the red switch is the most reliable.
    for (int i = 2; i < model.objects.length && i <= 3; ++i) {
      if (model.objects[i].leftFound && model.objects[i].rightFound) {
        dx += model.objects[i].leftX - model.objects[i].rightX;
        dy += model.objects[i].leftY - model.objects[i].rightY;
        ++count;
      }
    }
    if (count == 0) {
      System.out.println("No objects present in both; skipping reconciliation.");
      return;
    }
    dx /= count;
    dy /= count;
    for (int i = 0; i < model.objects.length; ++i) {
      if (!model.objects[i].leftFound && model.objects[i].rightFound) {
        System.out.println("Found right but not left for object " + i + "; reconciling.");
        model.objects[i].leftX = model.objects[i].rightX + dx;
        model.objects[i].leftY = model.objects[i].rightY + dy;
        model.objects[i].leftFound = true;
        markPos(imgLeft, model.objects[i].leftX, model.objects[i].leftY, "R" + Integer.toString(i));
      } else if (model.objects[i].leftFound && !model.objects[i].rightFound) {
        System.out.println("Found left but not right for object " + i + "; reconciling.");
        model.objects[i].rightX = model.objects[i].leftX - dx;
        model.objects[i].rightY = model.objects[i].leftY - dy;
        model.objects[i].rightFound = true;
        markPos(imgRight, model.objects[i].rightX, model.objects[i].rightY, "R" + Integer.toString(i));
      }
    }
  }

  public static void main(String[] args) throws Exception {
    RobonautEyeTester.debug = false;
    int[] arrLeftOrig = RobonautEyeTester.imageToArray(args[0]);
    int[] arrRightOrig = RobonautEyeTester.imageToArray(args[0].replaceFirst("Left", "Right"));

    /*int origHeight = arrLeftOrig[0];
    int origWidth = arrLeftOrig[1];

    int newHeight = origHeight / 2;
    int newWidth = origWidth / 2;
    int[] arrLeft = new int[2 + newHeight * newWidth];
    int[] arrRight = new int[2 + newHeight * newWidth];
    arrLeft[0] = newHeight;
    arrLeft[1] = newWidth;
    arrRight[0] = newHeight;
    arrRight[1] = newWidth;
    for (int y = 0; y < newHeight; ++y) {
      for (int x = 0; x < newWidth; ++x) {
        arrLeft[2 + x + y * newWidth] = arrLeftOrig[2 + 2 * x + 2 * y * origWidth];
        arrRight[2 + x + y * newWidth] = arrRightOrig[2 + 2 * x + 2 * y * origWidth];
      }
    }*/

    int[] arrLeft = arrLeftOrig;
    int[] arrRight = arrRightOrig;

    //arrLeft = filter(arrLeft);
    //arrRight = filter(arrRight);
    //arrLeft = mixLum(arrLeft);
    //arrRight = mixLum(arrRight);

    int width = arrLeft[1];
    int height = arrLeft[0];

    boolean isSmall = false;
    if (height > 2000 && width > 2000) {
      compRatiosOrthoLeft = compRatiosOrthoLargeLeft;
      compRatiosOrthoRight = compRatiosOrthoLargeRight;
      ImageProcessor.UPPER_NIBLACK = 1.5;
    } else if (height < 1700 && width < 1700) {
      isSmall = true;
      compRatiosOrthoLeft = compRatiosOrthoSmallLeft;
      compRatiosOrthoRight = compRatiosOrthoSmallRight;
      ImageProcessor.UPPER_NIBLACK = 1.5;
    }
    final BufferedImage imgLeft = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
    final BufferedImage imgRight = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
    for (int y = 0, cur = 2; y < height; ++y) {
      for (int x = 0; x < width; ++x, ++cur) {
        imgLeft.setRGB(x, y, arrLeft[cur] | 0xff000000);
        imgRight.setRGB(x, y, arrRight[cur] | 0xff000000);
      }
    }
    //drawParallels(imgLeft, imgRight);

    if (args.length != 0) {

    // Ideas here.
    long startTime = System.currentTimeMillis();
    SumImage sumLeft = new SumImage(arrLeft);
    SumImage sumRight = new SumImage(arrRight);
    {
      System.out.println("Done creating sumImages in " + (System.currentTimeMillis() - startTime) + " ms");
      int kernel = width / 5;
      ImageProcessor.doNiblackMulti(arrLeft, imgLeft, sumLeft, kernel);
      ImageProcessor.doNiblackMulti(arrRight, imgRight, sumRight, kernel);
      System.out.println("Done creating niBlack in " + (System.currentTimeMillis() - startTime) + " ms");
    }

    // Appears to be good at finding circles, white label.
    int[][] sccLeft = new int[width][height];
    Map<Integer, SccProperty> sccPropLeft = new TreeMap<Integer, SccProperty>();
    Map<Integer, SccProperty> sccPropLargeLeft = new TreeMap<Integer, SccProperty>();
    Scc.findScc(arrLeft, imgLeft, sccLeft, sccPropLeft, sccPropLargeLeft, sumLeft);
    int[][] sccRight = new int[width][height];
    Map<Integer, SccProperty> sccPropRight = new TreeMap<Integer, SccProperty>();
    Map<Integer, SccProperty> sccPropLargeRight = new TreeMap<Integer, SccProperty>();
    Scc.findScc(arrRight, imgRight, sccRight, sccPropRight, sccPropLargeRight, sumRight);

    boolean foundLedLeft = false;
    SccProperty candidatePanelLedLeft = null;
    {
      ArrayList<SccProperty> allGreenLeds = new ArrayList<SccProperty>();
      HashSet<Integer> allGreenLedsNum = new HashSet<Integer>();
      for (SccProperty prop : sccPropLeft.values()) {
        if (prop.isGreenLed && !allGreenLedsNum.contains(prop.propNum)) {
          allGreenLedsNum.add(prop.propNum);
          allGreenLeds.add(prop);
        }
      }
      if (allGreenLeds.size() > 1) {
        Collections.sort(allGreenLeds, new Comparator<SccProperty>() {
          @Override
          public int compare(SccProperty p0, SccProperty p1) {
            return Integer.valueOf(p0.y).compareTo(Integer.valueOf(p1.y));
          }
        });
        // Make sure they're at least a minimum distance away from each other.
        if (allGreenLeds.get(1).y - allGreenLeds.get(0).y > height / 15) {
          foundLedLeft = true;
          candidatePanelLedLeft = allGreenLeds.get(0);
        } else {
          System.out.println("Rejecting LED: " + allGreenLeds.get(0).x + "," + allGreenLeds.get(0).y
              + " because of proximity to " + allGreenLeds.get(1).x + "," + allGreenLeds.get(1).y
              + " propNum0: " + allGreenLeds.get(0).propNum + " propNum1: " + allGreenLeds.get(1).propNum);
        }
      } else if (allGreenLeds.size() == 1) {
        foundLedLeft = true;
        candidatePanelLedLeft = allGreenLeds.get(0);
      }
    }
    boolean foundLedRight = false;
    SccProperty candidatePanelLedRight = null;
    {
      ArrayList<SccProperty> allGreenLeds = new ArrayList<SccProperty>();
      HashSet<Integer> allGreenLedsNum = new HashSet<Integer>();
      for (SccProperty prop : sccPropRight.values()) {
        if (prop.isGreenLed && !allGreenLedsNum.contains(prop.propNum)) {
          allGreenLedsNum.add(prop.propNum);
          allGreenLeds.add(prop);
        }
      }
      if (allGreenLeds.size() > 1) {
        Collections.sort(allGreenLeds, new Comparator<SccProperty>() {
          @Override
          public int compare(SccProperty p0, SccProperty p1) {
            return Integer.valueOf(p0.y).compareTo(Integer.valueOf(p1.y));
          }
        });
        // Make sure they're at least a minimum distance away from each other.
        if (allGreenLeds.get(1).y - allGreenLeds.get(0).y > height / 15) {
          foundLedRight = true;
          candidatePanelLedRight = allGreenLeds.get(0);
        } else {
          System.out.println("Rejecting LED: " + allGreenLeds.get(0).x + "," + allGreenLeds.get(0).y
              + " because of proximity to " + allGreenLeds.get(1).x + "," + allGreenLeds.get(1).y
              + " propNum0: " + allGreenLeds.get(0).propNum + " propNum1: " + allGreenLeds.get(1).propNum);
        }
      } else if (allGreenLeds.size() == 1) {
        foundLedRight = true;
        candidatePanelLedRight = allGreenLeds.get(0);
      }
    }

    ModelSummary model = new ModelSummary();
    if (foundLedLeft && candidatePanelLedLeft != null) {
      findRed(arrLeft, imgLeft, model, true,
          sccLeft, sccPropLeft, sccPropLargeLeft,
          candidatePanelLedLeft.x - width / 10, candidatePanelLedLeft.y - height / 8,
          candidatePanelLedLeft.x + width / 80, candidatePanelLedLeft.y + height / 8,
          .05);
      if (model.objects[1].leftFound) {
        // If we found it, commit both the cover as well as the led.
        model.objects[2].leftFound = true;
        model.objects[2].leftX = candidatePanelLedLeft.x;
        model.objects[2].leftY = candidatePanelLedLeft.y;
        markPos(imgLeft, model.objects[2].leftX, model.objects[2].leftY, "2");
      }
    }
    if (foundLedRight && candidatePanelLedRight != null) {
      findRed(arrRight, imgRight, model, false,
          sccRight, sccPropRight, sccPropLargeRight,
          candidatePanelLedRight.x - width / 10, candidatePanelLedRight.y - height / 8,
          candidatePanelLedRight.x + width / 80, candidatePanelLedRight.y + height / 8,
          .05);
      if (model.objects[1].rightFound) {
        // If we found it, commit both the cover as well as the led.
        model.objects[2].rightFound = true;
        model.objects[2].rightX = candidatePanelLedRight.x;
        model.objects[2].rightY = candidatePanelLedRight.y;
        markPos(imgRight, model.objects[2].rightX, model.objects[2].rightY, "2");
      }
    }

    if (!model.objects[1].leftFound) {
      findRed(arrLeft, imgLeft, model, true,
          sccLeft, sccPropLeft, sccPropLargeLeft,
          width / 6, 0, width, height * 3 / 4, 0);
    }
    if (!model.objects[1].rightFound) {
      findRed(arrRight, imgRight, model, false,
          sccRight, sccPropRight, sccPropLargeRight,
          width / 6, 0, width, height * 3 / 4, 0);
    }

    if (!model.objects[2].leftFound) {
      findPanelLed(arrLeft, imgLeft, sccLeft, sccPropLeft, sccPropLargeLeft, sumLeft, model, true);
    }
    if (!model.objects[2].rightFound) {
      findPanelLed(arrRight, imgRight, sccRight, sccPropRight, sccPropLargeRight, sumRight, model, false);
    }
    
    boolean panelIsProbablyOff = false;
    if (!foundLedLeft && !foundLedRight &&
        !model.objects[2].leftFound && !model.objects[2].rightFound) {
      panelIsProbablyOff = true;
      model.objects[0].state = "DOWN";
      model.objects[2].state = "OFF";
    }

    findRockerSwitch(arrLeft, imgLeft, sccLeft, sccPropLeft, sccPropLargeLeft, sumLeft, model, true);
    findRockerSwitch(arrRight, imgRight, sccRight, sccPropRight, sccPropLargeRight, sumRight, model, false);


    if (model.objects[2].leftFound && model.objects[3].leftFound &&
        model.objects[2].rightFound && model.objects[3].rightFound) {

      {
        double dx = model.objects[3].leftX - model.objects[2].leftX;
        double dy = model.objects[3].leftY - model.objects[2].leftY;
        double leftDist = Math.sqrt(dx * dx + dy * dy);
        System.out.println("leftDist:  " + leftDist);
      }
      {
        double dx = model.objects[3].rightX - model.objects[2].rightX;
        double dy = model.objects[3].rightY - model.objects[2].rightY;
        double rightDist = Math.sqrt(dx * dx + dy * dy);
        System.out.println("rightDist:  " + rightDist);
      }
    }

    // Rocker might have a pretty good guess. This solves Lab2/MPUI.
    /*if (!model.objects[3].leftFound) {
      ArrayList<SccProperty> allPossibleRockerSwitches = new ArrayList<SccProperty>();
      for (SccProperty prop : sccPropLeft.values()) {
        if (prop.isRockerSwitch) {
          allPossibleRockerSwitches.add(prop);
        }
      }
      if (allPossibleRockerSwitches.size() == 1) {
        SccProperty prop = allPossibleRockerSwitches.get(0);
        model.objects[3].leftX = prop.x;
        model.objects[3].leftY = prop.y;
        model.objects[3].leftFound = true;
        markPos(imgLeft, model.objects[3].leftX, model.objects[3].leftY, "EL3");

        if (!model.objects[2].leftFound) {
          model.objects[2].leftX = model.objects[3].leftX;
          model.objects[2].leftY = model.objects[3].leftY - height / 5;
          model.objects[2].leftFound = true;
          markPosRed(imgLeft, model.objects[2].leftX, model.objects[2].leftY, "WEL2");
        }
      } 
    }
    if (!model.objects[3].rightFound) {
      ArrayList<SccProperty> allPossibleRockerSwitches = new ArrayList<SccProperty>();
      for (SccProperty prop : sccPropRight.values()) {
        if (prop.isRockerSwitch) {
          allPossibleRockerSwitches.add(prop);
        }
      }
      if (allPossibleRockerSwitches.size() == 1) {
        SccProperty prop = allPossibleRockerSwitches.get(0);
        model.objects[3].rightX = prop.x;
        model.objects[3].rightY = prop.y;
        model.objects[3].rightFound = true;
        markPosRed(imgRight, model.objects[3].rightX, model.objects[3].rightY, "EL3");
        if (!model.objects[2].rightFound) {
          model.objects[2].rightX = model.objects[3].rightX;
          model.objects[2].rightY = model.objects[3].rightY - height / 5;
          model.objects[2].rightFound = true;
          markPosRed(imgRight, model.objects[2].rightX, model.objects[2].rightY, "WEL2");
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
        model.objects[1].leftProp = null;
        markPosRed(imgLeft, model.objects[1].leftX, model.objects[1].leftY, "ADJ1");

        // Redo all the things which depended on the red switch.
        findPanelLed(arrLeft, imgLeft, sccLeft, sccPropLeft, sccPropLargeLeft, sumLeft, model, true);
        findRockerSwitch(arrLeft, imgLeft, sccLeft, sccPropLeft, sccPropLargeLeft, sumLeft, model, true);
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
        model.objects[1].rightProp = null;
        markPosRed(imgRight, model.objects[1].rightX, model.objects[1].rightY, "ADJ1");

        // Redo all the things which depended on the red switch.
        findPanelLed(arrRight, imgRight, sccRight, sccPropRight, sccPropLargeRight, sumRight, model, false);
        findRockerSwitch(arrRight, imgRight, sccRight, sccPropRight, sccPropLargeRight, sumRight, model, false);
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
      if (predX >= width) predX = width - 1;
      int predY = model.objects[1].leftY;
      SccProperty bestProp = null;
      if (panelIsProbablyOff) {
        int ymin = 0;
        int ymax = height;
        if (model.objects[1].leftProp != null) {
          ymin = model.objects[1].leftProp.ymin;
          ymax = model.objects[1].leftProp.ymax;
        }
        long bestDist = (width + height) / 2 / 20;
        bestDist *= bestDist;
        for (SccProperty prop : sccPropLeft.values()) {
          if (prop.x > model.objects[1].leftX &&
              prop.y >= ymin && prop.y < ymax &&
              prop.xmax - prop.xmin < width / 40 &&
              prop.ymax - prop.ymin < height / 30 &&
              Math.max(prop.r, Math.max(prop.g, prop.b)) < 100 &&
              canReach(imgLeft, width, height, prop.x, prop.y, predX, predY, sccLeft, sccPropLargeLeft)) {
            int dx = prop.x - predX;
            int dy = prop.y - predY;
            if (dx * dx + dy * dy < bestDist) {
              bestDist = dx * dx + dy * dy;
              bestProp = prop;
            }
          }
        }
      }
      if (bestProp != null) {
        predX = bestProp.x;
        predY = bestProp.y;
        markPosRed(imgLeft, predX, predY, "WM2");
      } else {
        markPosRed(imgLeft, predX, predY, "W2");
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
      if (predX >= width) predX = width - 1;
      int predY = model.objects[1].rightY;
      SccProperty bestProp = null;
      if (panelIsProbablyOff) {
        int ymin = 0;
        int ymax = height;
        if (model.objects[1].rightProp != null) {
          ymin = model.objects[1].rightProp.ymin;
          ymax = model.objects[1].rightProp.ymax;
        }
        long bestDist = (width + height) / 2 / 20;
        bestDist *= bestDist;
        for (SccProperty prop : sccPropRight.values()) {
          if (prop.x > model.objects[1].rightX &&
              prop.y >= ymin && prop.y < ymax &&
              prop.xmax - prop.xmin < width / 40 &&
              prop.ymax - prop.ymin < height / 30 &&
              Math.max(prop.r, Math.max(prop.g, prop.b)) < 100 &&
              canReach(imgRight, width, height, prop.x, prop.y, predX, predY, sccRight, sccPropLargeRight)) {
            int dx = prop.x - predX;
            int dy = prop.y - predY;
            if (dx * dx + dy * dy < bestDist) {
              bestDist = dx * dx + dy * dy;
              bestProp = prop;
            }
          }
        }
      }
      if (bestProp != null) {
        predX = bestProp.x;
        predY = bestProp.y;
        markPosRed(imgRight, predX, predY, "WM2");
      } else {
        markPosRed(imgRight, predX, predY, "W2");
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
      markPosRed(imgLeft, model.objects[3].leftX, model.objects[3].leftY, "W3");
    }
    if (model.objects[2].rightFound && !model.objects[3].rightFound) {
      model.objects[3].rightX = model.objects[2].rightX;
      model.objects[3].rightY = model.objects[2].rightY + height / 5;
      model.objects[3].rightFound = true;
      markPosRed(imgRight, model.objects[3].rightX, model.objects[3].rightY, "W3");
    }

    /*if (model.objects[2].leftFound && model.objects[3].leftFound) {
      dfsMatchStrict = true;
      if (!matchSccToBothLedsDfsRect(
              arrLeft, imgLeft, sccLeft, sccPropLeft, sccPropLargeLeft, model, true)) {
        dfsMatchStrict = false;
        if (model.objects[2].leftFound && model.objects[3].leftFound) {
          matchSccToGreenLedsDfsRect(
              arrLeft, imgLeft, sccLeft, sccPropLeft, sccPropLargeLeft, model, true);
        }
        if (model.objects[2].leftFound && model.objects[3].leftFound) {
          matchSccToBlueLedsDfsRect(
              arrLeft, imgLeft, sccLeft, sccPropLeft, sccPropLargeLeft, model, true);
        }
      }
    }
    if (model.objects[2].rightFound && model.objects[3].rightFound) {
      dfsMatchStrict = true;
      if (!matchSccToBothLedsDfsRect(
              arrRight, imgRight, sccRight, sccPropRight, sccPropLargeRight, model, false)) {
        dfsMatchStrict = false;
        if (model.objects[2].rightFound && model.objects[3].rightFound) {
          matchSccToGreenLedsDfsRect(
              arrRight, imgRight, sccRight, sccPropRight, sccPropLargeRight, model, false);
        }
        if (model.objects[2].rightFound && model.objects[3].rightFound) {
          matchSccToBlueLedsDfsRect(
              arrRight, imgRight, sccRight, sccPropRight, sccPropLargeRight, model, false);
        }
      }
    }
    dfsMatchStrict = false;*/

    if (model.objects[2].leftFound && model.objects[3].leftFound) {
      panelLedForDfs = model.objects[2].leftProp;
      rockerForDfs = model.objects[3].leftProp;
      matchSccToGreenLedsDfsRect(
          arrLeft, imgLeft, sccLeft, sccPropLeft, sccPropLargeLeft, model, true);
    }
    if (model.objects[2].rightFound && model.objects[3].rightFound) {
      panelLedForDfs = model.objects[2].rightProp;
      rockerForDfs = model.objects[3].rightProp;
      matchSccToGreenLedsDfsRect(
          arrRight, imgRight, sccRight, sccPropRight, sccPropLargeRight, model, false);
    }
    if (model.objects[2].leftFound && model.objects[3].leftFound) {
      panelLedForDfs = model.objects[2].leftProp;
      rockerForDfs = model.objects[3].leftProp;
      matchSccToBlueLedsDfsRect(
          arrLeft, imgLeft, sccLeft, sccPropLeft, sccPropLargeLeft, model, true);
    }
    if (model.objects[2].rightFound && model.objects[3].rightFound) {
      panelLedForDfs = model.objects[2].rightProp;
      rockerForDfs = model.objects[3].rightProp;
      matchSccToBlueLedsDfsRect(
          arrRight, imgRight, sccRight, sccPropRight, sccPropLargeRight, model, false);
    }
    panelLedForDfs = null;
    rockerForDfs = null;

    fillGrid(arrLeft, imgLeft, sccPropLeft, model, true);
    fillGrid(arrRight, imgRight, sccPropRight, model, false);

    if (model.objects[2].leftFound && model.objects[3].leftFound) {
      interpolateRockerLeds(arrLeft, imgLeft, sccLeft, sccPropLeft, model, true);
    }
    if (model.objects[2].rightFound && model.objects[3].rightFound) {
      interpolateRockerLeds(arrRight, imgRight, sccRight, sccPropRight, model, false);
    }
    // For now, correct for red switch position just for scoring purposes.
    boolean leftCoverDown = false;
    if (model.objects[1].leftProp != null) {
      SccProperty prop = model.objects[1].leftProp;
      boolean isTall = false;
      if (prop.xmax - prop.xmin < (prop.ymax - prop.ymin) * 4 / 5) {
        isTall = true;
      }
      // See if xmiddle/ymax is close to the component; if not, the cover must be up.
      boolean isClose = false;
      for (int x = Math.max(prop.x - 10, 0); x < Math.min(prop.x + 10, width); ++x) {
        for (int y = Math.max(prop.ymax - 10, 0); y < Math.min(prop.ymax + 10, height); ++y) {
          if (sccPropLeft.get(sccLeft[x][y]) != null &&
              sccPropLeft.get(sccLeft[x][y]) == prop) {
            isClose = true;
            break;
          }
        }
        if (isClose) break;
      }
      int markX, markY;
      if (isTall && isClose) {
        System.out.println("Cover appears to be down for calculating its position");
        markX = (prop.ymaxx + prop.yminx) / 2;
        markY = prop.ymax;
        markPos(imgLeft, markX, markY, "1D");
      } else {
        System.out.println("Cover appears to be up for calculating its position");
        if (prop.xminy < prop.xmaxy) {
          markX = prop.xmin + (prop.yminx - prop.xmin) / 4;
          markY = prop.xminy;
          markPos(imgLeft, markX, markY, "1L");
        } else {
          markX = prop.xmax - (prop.xmax - prop.yminx) / 4;
          markY = prop.xmaxy;
          markPos(imgLeft, markX, markY, "1R");
        }
      }
      model.objects[1].leftX = markX;
      model.objects[1].leftY = markY;
      leftCoverDown = isTall && isClose;
    }
    boolean rightCoverDown = false;
    if (model.objects[1].rightProp != null) {
      SccProperty prop = model.objects[1].rightProp;
      boolean isTall = false;
      if (prop.xmax - prop.xmin < (prop.ymax - prop.ymin) * 4 / 5) {
        isTall = true;
      }
      // See if xmiddle/ymax is close to the component; if not, the cover must be up.
      boolean isClose = false;
      for (int x = Math.max(prop.x - 10, 0); x < Math.min(prop.x + 10, width); ++x) {
        for (int y = Math.max(prop.ymax - 10, 0); y < Math.min(prop.ymax + 10, height); ++y) {
          if (sccPropRight.get(sccRight[x][y]) != null &&
              sccPropRight.get(sccRight[x][y]) == prop) {
            isClose = true;
            break;
          }
        }
        if (isClose) break;
      }
      int markX, markY;
      if (isTall && isClose) {
        System.out.println("Cover appears to be down for calculating its position");
        markX = (prop.ymaxx + prop.yminx) / 2;
        markY = prop.ymax;
        markPos(imgRight, markX, markY, "1D");
      } else {
        System.out.println("Cover appears to be up for calculating its position");
        if (prop.xminy < prop.xmaxy) {
          markX = prop.xmin + (prop.yminx - prop.xmin) / 4;
          markY = prop.xminy;
          markPos(imgRight, markX, markY, "1L");
        } else {
          markX = prop.xmax - (prop.xmax - prop.yminx) / 4;
          markY = prop.xmaxy;
          markPos(imgRight, markX, markY, "1R");
        }
      }
      model.objects[1].rightX = markX;
      model.objects[1].rightY = markY;
      rightCoverDown = isTall && isClose;
    }
    if (leftCoverDown && rightCoverDown) {
      System.out.println("Cover state is DOWN");
    } else {
      System.out.println("Cover state is UP");
    }

    // Try to estimate center of grid by taking average of scc points to the left of the rocker.
    /*if (model.objects[1].leftFound && model.objects[1].leftProp != null &&
        model.objects[3].leftFound && model.objects[3].leftProp != null &&
        model.objects[5].leftFound) {
      int ymin = model.objects[1].leftProp.ymax;
      int ymax = model.objects[5].leftY;
      int xmax = model.objects[3].leftProp.xmin;
      if (xmax >= width) xmax = width - 1;
      if (xmax < 0) xmax = 0;
      if (ymin >= height) ymin = height - 1;
      if (ymin < 0) ymin = 0;
      if (ymax >= height) ymax = height - 1;
      if (ymax < 0) ymax = 0;
      long xsum = 0;
      int xsamples = 0;
      for (int y = ymin; y < ymax; ++y) {
        for (int x = xmax; x >= 0; --x) {
          if (sccPropLargeLeft.containsKey(sccLeft[x][y])) {
            break;
          }
          if (sccPropLeft.containsKey(sccLeft[x][y])) {
            xsum += x;
            ++xsamples;
          }
          if (x % 4 == 0 && y % 4 == 0) {
            imgLeft.setRGB(x, y, Color.green.getRGB());
          }
        }
      }
      model.objects[10].leftX = (int) (xsum / xsamples);
      model.objects[10].leftY = model.objects[3].leftProp.y;
      model.objects[10].leftFound = true;
    }*/

    // Basis extrapolate is pretty good for 21.
    for (int i = 0; i < 22; ++i) {
      blockBasisExtrapolate[i] = true;
    }
    blockBasisExtrapolate[21] = false;
    if (model.objects[2].leftFound && model.objects[3].leftFound) {
      if (model.objects[10].leftFound) {
        basisExtrapolate(arrLeft, imgLeft, sccPropLeft, model, true);
      }
      else {
        basisExtrapolateWeak(arrLeft, imgLeft, sccPropLeft, model, true);
      }
    }
    if (model.objects[2].rightFound && model.objects[3].rightFound) {
      if (model.objects[10].rightFound) {
        basisExtrapolate(arrRight, imgRight, sccPropRight, model, false);
      } else {
        basisExtrapolateWeak(arrRight, imgRight, sccPropRight, model, false);
      }
    }
    for (int i = 0; i < 22; ++i) {
      blockBasisExtrapolate[i] = false;
    }
    fillGrid(arrLeft, imgLeft, sccPropLeft, model, true);
    fillGrid(arrRight, imgRight, sccPropRight, model, false);

    /*if (model.objects[21].leftFound && !model.objects[18].leftFound && !model.objects[16].leftFound) {
      int objX = model.objects[21].leftX;
      int objY = model.objects[21].leftY;
      ArrayList<SccProperty> candidates = new ArrayList<SccProperty>();
      for (SccProperty prop : sccPropLeft.values()) {
        if (prop.x < objX - width / 25 &&
            prop.size > width / 4 &&
            prop.ymax - prop.ymin > height / 200 &&
            Math.abs(prop.y - objY) < height / 21 &&
            !prop.isGreenLed && !prop.isBlueLed &&
            canReach(imgLeft, width, height, prop.x, prop.y, objX, objY, sccLeft, sccPropLargeLeft)) {
          //markPosBlue(imgLeft, prop.x, prop.y, "PROP");
          candidates.add(prop);
        }
      }
      double bestScore = 999;
      SccProperty cand16 = null;
      SccProperty cand18 = null;
      for (int i = 0; i < candidates.size(); ++i) {
        SccProperty prop16 = candidates.get(i);
        for (int j = 0; j < candidates.size(); ++j) {
          if (i == j) continue;
          SccProperty prop18 = candidates.get(j);
          if (prop16.x >= prop18.x) continue;
          int dx16 = prop18.x - prop16.x;
          int dy16 = prop18.y - prop16.y;
          double diff16 = Math.sqrt(dx16 * dx16 + dy16 * dy16);
          int dx18 = objX - prop18.x;
          int dy18 = objY - prop18.y;
          double diff18 = Math.sqrt(dx18 * dx18 + dy18 * dy18);
          double pairDiff = Math.abs(diff16 - diff18);
          double theta = (dx16 * dx18 + dy16 * dy18) / diff16 / diff18;
          theta = Math.abs(1.0 - theta);
          double score = (theta + 1) * (pairDiff + 1);
          if (score < bestScore) {
            bestScore = score;
            cand16 = prop16;
            cand18 = prop18;
          }
        }
      }
      if (cand16 != null && cand18 != null) {
        System.out.println("Found 16 and 18 via linear matching: " + cand16.x + "," + cand16.y + " "
            + cand18.x + "," + cand18.y + " with score " + bestScore);
        model.objects[16].leftX = cand16.x;
        model.objects[16].leftY = cand16.y;
        model.objects[16].leftFound = true;
        model.objects[16].leftProp = cand16;
        model.objects[18].leftX = cand18.x;
        model.objects[18].leftY = cand18.y;
        model.objects[18].leftFound = true;
        model.objects[18].leftProp = cand18;
        markPos(imgLeft, cand16.x, cand16.y, "L16");
        markPos(imgLeft, cand18.x, cand18.y, "L18");
      }
    }
    if (model.objects[21].rightFound && !model.objects[18].rightFound && !model.objects[16].rightFound) {
      int objX = model.objects[21].rightX;
      int objY = model.objects[21].rightY;
      ArrayList<SccProperty> candidates = new ArrayList<SccProperty>();
      for (SccProperty prop : sccPropRight.values()) {
        if (prop.x < objX - width / 25 &&
            prop.size > width / 4 &&
            prop.ymax - prop.ymin > height / 200 &&
            Math.abs(prop.y - objY) < height / 21 &&
            !prop.isGreenLed && !prop.isBlueLed &&
            canReach(imgRight, width, height, prop.x, prop.y, objX, objY, sccRight, sccPropLargeRight)) {
          //markPosBlue(imgRight, prop.x, prop.y, "PROP");
          candidates.add(prop);
        }
      }
      double bestScore = 999;
      SccProperty cand16 = null;
      SccProperty cand18 = null;
      for (int i = 0; i < candidates.size(); ++i) {
        SccProperty prop16 = candidates.get(i);
        for (int j = 0; j < candidates.size(); ++j) {
          if (i == j) continue;
          SccProperty prop18 = candidates.get(j);
          if (prop16.x >= prop18.x) continue;
          int dx16 = prop18.x - prop16.x;
          int dy16 = prop18.y - prop16.y;
          double diff16 = Math.sqrt(dx16 * dx16 + dy16 * dy16);
          int dx18 = objX - prop18.x;
          int dy18 = objY - prop18.y;
          double diff18 = Math.sqrt(dx18 * dx18 + dy18 * dy18);
          double pairDiff = Math.abs(diff16 - diff18);
          double theta = (dx16 * dx18 + dy16 * dy18) / diff16 / diff18;
          theta = Math.abs(1.0 - theta);
          double score = (theta + 1) * (pairDiff + 1);
          if (score < bestScore) {
            bestScore = score;
            cand16 = prop16;
            cand18 = prop18;
          }
        }
      }
      if (cand16 != null && cand18 != null) {
        System.out.println("Found 16 and 18 via linear matching: " + cand16.x + "," + cand16.y + " "
            + cand18.x + "," + cand18.y + " with score " + bestScore);
        model.objects[16].rightX = cand16.x;
        model.objects[16].rightY = cand16.y;
        model.objects[16].rightFound = true;
        model.objects[16].rightProp = cand16;
        model.objects[18].rightX = cand18.x;
        model.objects[18].rightY = cand18.y;
        model.objects[18].rightFound = true;
        model.objects[18].rightProp = cand18;
        markPos(imgRight, cand16.x, cand16.y, "L16");
        markPos(imgRight, cand18.x, cand18.y, "L18");
      }
    }*/

    if (model.objects[2].leftFound && model.objects[3].leftFound) {
      if (model.objects[10].leftFound) {
        basisDfsMatch(arrLeft, imgLeft, sccLeft, sccPropLeft, sccPropLargeLeft, model, true);
      }
      else {
        basisDfsMatchWeak(arrLeft, imgLeft, sccLeft, sccPropLeft, sccPropLargeLeft, model, true);
      }
    }
    if (model.objects[2].rightFound && model.objects[3].rightFound) {
      if (model.objects[10].rightFound) {
        basisDfsMatch(arrRight, imgRight, sccRight, sccPropRight, sccPropLargeRight, model, false);
      }
      else {
        basisDfsMatchWeak(arrRight, imgRight, sccRight, sccPropRight, sccPropLargeRight, model, false);
      }
    }

    if (model.objects[2].leftFound && model.objects[3].leftFound) {
      if (model.objects[10].leftFound) {
        basisExtrapolate(arrLeft, imgLeft, sccPropLeft, model, true);
      }
      else {
        basisExtrapolateWeak(arrLeft, imgLeft, sccPropLeft, model, true);
      }
    }
    if (model.objects[2].rightFound && model.objects[3].rightFound) {
      if (model.objects[10].rightFound) {
        basisExtrapolate(arrRight, imgRight, sccPropRight, model, false);
      } else {
        basisExtrapolateWeak(arrRight, imgRight, sccPropRight, model, false);
      }
    }
    fillGrid(arrLeft, imgLeft, sccPropLeft, model, true);
    fillGrid(arrRight, imgRight, sccPropRight, model, false);
    reconcileMissing(imgLeft, imgRight, model);

    long finishTime = System.currentTimeMillis();
    System.out.println("Idea took: " + (finishTime - startTime) + " ms");
    }

    JFrame frame = new JFrame("Idea Tester");
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    JPanel panel = new JPanel() {
      @Override
      public void paintComponent(Graphics g) {
        if (g == null) return;
        g.drawImage(imgLeft.getScaledInstance(
            getWidth() / 2, getHeight(), Image.SCALE_SMOOTH), 0, 0, null);
        g.drawImage(imgRight.getScaledInstance(
            getWidth() / 2, getHeight(), Image.SCALE_SMOOTH), getWidth() / 2, 0, null);
      }
      @Override
      public void repaint() {
        //super.repaint();
        paintComponent(getGraphics());
      }
    };
    panel.setPreferredSize(new Dimension(imgLeft.getWidth() * 2 / 3, imgLeft.getHeight() / 3));

    frame.getContentPane().add(panel);
    frame.pack();
    frame.setVisible(true);

   /* int[][] sccLeft = new int[width][height];
    int[][] sccRight = new int[width][height];
    while (true) {
      Scc.minComponentSize += 1000;
      Scc.findScc(arrLeft, imgLeft);
      Scc.findScc(arrRight, imgRight);
      panel.repaint();
      try {
        Thread.sleep(5000);
      } catch (InterruptedException ie) {
      }
    }*/
  }
}
