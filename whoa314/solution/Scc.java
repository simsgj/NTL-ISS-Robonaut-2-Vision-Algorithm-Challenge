import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import javax.imageio.*; 
import javax.swing.*;

public class Scc {
  static class Pos {
    int x, y;
    public Pos(int x, int y) {
      this.x = x; this.y = y;
    }
  }

  // Prop0 will contain the merge of the two.
  public static void mergeScc(SccProperty prop0, SccProperty prop1) {
    System.out.println(
        "Merging " + prop1.x + "," + prop1.y + " into " + prop0.x + "," + prop0.y);
    prop0.x = (int)((((long)prop0.x) * prop0.size + ((long)prop1.x) * prop1.size) /
        (prop0.size + prop1.size));
    prop0.y = (int)((((long)prop0.y) * prop0.size + ((long)prop1.y) * prop1.size) /
        (prop0.size + prop1.size));
    prop0.r = (int)((((long)prop0.r) * prop0.size + ((long)prop1.r) * prop1.size) /
        (prop0.size + prop1.size));
    prop0.g = (int)((((long)prop0.g) * prop0.size + ((long)prop1.g) * prop1.size) /
        (prop0.size + prop1.size));
    prop0.b = (int)((((long)prop0.b) * prop0.size + ((long)prop1.b) * prop1.size) /
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

  /**
   * Find the strongly connected components.
   */
  public static void findScc(int[] arr, BufferedImage img, int[][] scc,
      Map<Integer, SccProperty> sccProp, Map<Integer, SccProperty> sccPropLarge,
      SumImage sumImage) {
    int width = arr[1];
    int height = arr[0];
    int maxComponentSize = width * height / 50;
    int minComponentSize = width * height / 45000;
    System.out.println("Min component szie: " + minComponentSize);

    int[][] sccRgb = new int[width][height];
    int curComponent = 1;
    int numInterestingComponents = 0;
    Random rand = new Random();
    ArrayList<Integer> greenLedKeys = new ArrayList<Integer>();
    ArrayList<SccProperty> greenLedValues = new ArrayList<SccProperty>();
    ArrayList<Integer> blueLedKeys = new ArrayList<Integer>();
    ArrayList<SccProperty> blueLedValues = new ArrayList<SccProperty>();
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        if (scc[x][y] != 0) {
          continue;
        }
        if (!sumImage.isForeground[x][y] && !sumImage.isBackground[x][y]) continue;
        boolean localIsForeground = sumImage.isForeground[x][y];
        // Start a new fill.
        LinkedList<Pos> q = new LinkedList<Pos>();
        scc[x][y] = curComponent;
        int sharedRgb = (int) (rand.nextInt() | 0xff000000);
        q.add(new Pos(x, y));
        int idx = 2 + y * width + x;
        int r = (arr[idx] >> 16) & 0x00ff;
        int g = (arr[idx] >> 8) & 0x00ff;
        int b = (arr[idx]) & 0x00ff;
        int localHue = ImageProcessor.getHue(r, g, b);
        long localLum = r * r + g * g + b * b;
        int componentSize = 1;
        int compXmin = Integer.MAX_VALUE;
        int compXmax = 0;
        int compYmin = Integer.MAX_VALUE;
        int compYmax = 0;
        int compXminY = 0;
        int compXmaxY = 0;
        int compYminX = 0;
        int compYmaxX = 0;
        long sumSquaredR = r * r;
        long sumSquaredG = g * g;
        long sumSquaredB = b * b;
        long compR = r;
        long compG = g;
        long compB = b;
        long compX = x;
        long compY = y;
        while (!q.isEmpty()) {
          Pos cur = q.poll();
          if (cur.x > compXmax) {
            compXmax = cur.x;
            compXmaxY = cur.y;
          }
          if (cur.x < compXmin) {
            compXmin = cur.x;
            compXminY = cur.y;
          }
          if (cur.y > compYmax) {
            compYmax = cur.y;
            compYmaxX = cur.x;
          }
          if (cur.y < compYmin) {
            compYmin = cur.y;
            compYminX = cur.x;
          }
          /*int tmpidx = 2 + cur.y * width + cur.x;
          int tmpr = (arr[tmpidx] >> 16) & 0x00ff;
          int tmpg = (arr[tmpidx] >> 8) & 0x00ff;
          int tmpb = (arr[tmpidx]) & 0x00ff;
          int neighHue = ImageProcessor.getHue(tmpr, tmpg, tmpb);*/
          for (int dy = Math.max(cur.y - 1, 0); dy < Math.min(height, cur.y + 2); ++dy) {
            for (int dx = Math.max(cur.x - 1, 0); dx < Math.min(width, cur.x + 2); ++dx) {
              if (scc[dx][dy] != 0) {
                continue;
              }
              if (!sumImage.isForeground[dx][dy] && !sumImage.isBackground[dx][dy]) continue;
              if (localIsForeground != sumImage.isForeground[dx][dy]) continue;

              // Unvisited neighbor.
              int didx = 2 + dy * width + dx;
              long dr = (arr[didx] >> 16) & 0x00ff;
              long dg = (arr[didx] >> 8) & 0x00ff;
              long db = (arr[didx]) & 0x00ff;
              int curHue = ImageProcessor.getHue((int)dr, (int)dg, (int)db);
              long dr0 = dr;
              long dg0 = dg;
              long db0 = db;

              dr -= r;
              dg -= g;
              db -= b;
              dr *= dr;
              dg *= dg;
              db *= db;
              long diff = dr + dg + db;
              //int hueDiff = Math.abs(curHue - neighHue);
              //if (hueDiff > 180) {
              //  hueDiff = 360 - hueDiff;
              //}
              //if (hueDiff < 120) {
              //if (Math.abs(curHue - localHue) < 60) {
              //if (diff < Math.max(localLum, 10000)) {
              //{
                sumSquaredR += dr0 * dr0;
                sumSquaredG += dg0 * dg0;
                sumSquaredB += db0 * db0;
                compR += dr0;
                compG += dg0;
                compB += db0;
                compX += dx;
                compY += dy;
                sccRgb[dx][dy] = sharedRgb;
                //img.setRGB(dx, dy, sharedRgb);
                scc[dx][dy] = curComponent;
                ++componentSize;
                q.add(new Pos(dx, dy));
              //}
            }
          }
        }  // while
        if (componentSize < maxComponentSize) {
          sumSquaredR /= componentSize;
          sumSquaredG /= componentSize;
          sumSquaredB /= componentSize;
          boolean isGreenLed = componentSize > 20 &&
              (sumSquaredG - Math.min(sumSquaredR, sumSquaredB) > 14500) &&
              sumSquaredG - Math.max(sumSquaredR, sumSquaredB) > 2500 && localIsForeground;
          boolean isBlueLed = componentSize > 30 && componentSize < height * 3 / 2 &&
              (sumSquaredR < 11000 || sumSquaredB - sumSquaredR > 43000) &&
              sumSquaredB > 40000 &&
              sumSquaredB - Math.max(sumSquaredR, sumSquaredG) > 5000 && localIsForeground;
          /* &&
              Math.abs((compXmax - compXmin) - (compYmax - compYmin)) <
              Math.max(compXmax - compXmin, compYmax - compYmin);*/
          if (isGreenLed || isBlueLed || componentSize > minComponentSize) {
            SccProperty prop = new SccProperty();
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
              greenLedKeys.add(curComponent);
              greenLedValues.add(prop);
              System.out.println(prop.x + "," + prop.y + " (g): " + componentSize + " " +
                  sumSquaredR + "/" + sumSquaredG + "/" + sumSquaredB);
            } else if (isBlueLed) {
              prop.isGreenLed = false;
              prop.isBlueLed = true;
              blueLedKeys.add(curComponent);
              blueLedValues.add(prop);
              System.out.println(prop.x + "," + prop.y + " (b): " + componentSize + " " +
                  sumSquaredR + "/" + sumSquaredG + "/" + sumSquaredB);
            } else {
              prop.isGreenLed = false;
              prop.isBlueLed = false;
              sccProp.put(curComponent, prop);
            }

            ++numInterestingComponents;
          }
        } else {
          // Large component
          SccProperty prop = new SccProperty();
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
          prop.isGreenLed = false;
          prop.isBlueLed = false;
          prop.containsCenter = true;
          prop.isFirstFromLeft = false;
          sccPropLarge.put(curComponent, prop);
        }
        ++curComponent;
      }
    }

    // Merge components which are less than 2.5% of the average linear dimension away from each
    // other.
    // TODO: Bail if there are a crazy number of found leds.
    long maxDist = (width + height) / 2 * 25L / 1000L;
    maxDist *= maxDist;  // Compare in the squared space to avoid sqrt.
    for (int i = 0; i < greenLedKeys.size(); ++i) {
      if (sccProp.containsKey(greenLedKeys.get(i))) continue;
      SccProperty prop0 = greenLedValues.get(i);
      ArrayList<Integer> mergedComponents = new ArrayList<Integer>();
      mergedComponents.add(greenLedKeys.get(i));
      for (int j = i + 1; j < greenLedKeys.size(); ++j) {
        if (sccProp.containsKey(greenLedKeys.get(j))) continue;
        SccProperty prop1 = greenLedValues.get(j);
        int dx = prop0.x - prop1.x;
        int dy = prop0.y - prop1.y;
        if (dx * dx + dy * dy < maxDist) {
          mergeScc(prop0, prop1);
          mergedComponents.add(greenLedKeys.get(j));
        }
      }
      for (int k = 0; k < mergedComponents.size(); ++k) {
        sccProp.put(mergedComponents.get(k), prop0);
      }
    }
    for (int i = 0; i < blueLedKeys.size(); ++i) {
      if (sccProp.containsKey(blueLedKeys.get(i))) continue;
      SccProperty prop0 = blueLedValues.get(i);
      ArrayList<Integer> mergedComponents = new ArrayList<Integer>();
      mergedComponents.add(blueLedKeys.get(i));
      for (int j = i + 1; j < blueLedKeys.size(); ++j) {
        if (sccProp.containsKey(blueLedKeys.get(j))) continue;
        SccProperty prop1 = blueLedValues.get(j);
        int dx = prop0.x - prop1.x;
        int dy = prop0.y - prop1.y;
        if (dx * dx + dy * dy < maxDist) {
          mergeScc(prop0, prop1);
          mergedComponents.add(blueLedKeys.get(j));
        }
      }
      for (int k = 0; k < mergedComponents.size(); ++k) {
        sccProp.put(mergedComponents.get(k), prop0);
      }
    }

    /*long maxDistForMerge = (width + height) / 2 * 25L / 1000L;
    maxDistForMerge *= maxDistForMerge;
    ArrayList<SccProperty> allGridScc = new ArrayList<SccProperty>(sccProp.values());
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
    }*/

    boolean drawInteresting = true;
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        if (drawInteresting) {
          // Draw the "interesting" components.
          if (sccProp.containsKey(scc[x][y])) {
            //img.setRGB(x, y, sccRgb[x][y]);
            SccProperty prop = sccProp.get(scc[x][y]);
            int rgb = 0xff000000 | (prop.r << 16) | (prop.g << 8) | prop.b;
            img.setRGB(x, y, rgb);
          } else {
            img.setRGB(x, y, Color.white.getRGB());
          }
        } else {
          // Draw the "large" components.
          if (sccPropLarge.containsKey(scc[x][y]) && sumImage.isBackground[x][y]) {
            SccProperty prop = sccPropLarge.get(scc[x][y]);
            int rgb = 0xff000000 | (prop.r << 16) | (prop.g << 8) | prop.b;
            img.setRGB(x, y, rgb);
          } else {
            img.setRGB(x, y, Color.white.getRGB());
          }
        }
      }
    }
    /*IdeaTester.allowMarkPos = true;
    for (SccProperty prop : sccProp.values()) {
      IdeaTester.markPos(img, prop.x, prop.y, "");
    }
    IdeaTester.allowMarkPos = false;*/

    /*ArrayList<SccProperty> blackestOnes = new ArrayList<SccProperty>(sccProp.values());
    Collections.sort(blackestOnes, new Comparator<SccProperty>() {
      @Override
      public int compare(SccProperty prop0, SccProperty prop1) {
        int lum0 = prop0.r * prop0.r + prop0.g * prop0.g + prop0.b * prop0.b;
        int lum1 = prop1.r * prop1.r + prop1.g * prop1.g + prop1.b * prop1.b;
        return Integer.valueOf(lum0).compareTo(Integer.valueOf(lum1));
      }
    });
    for (int i = 0; i < blackestOnes.size() / 10; ++i) {
      IdeaTester.markPos(img, blackestOnes.get(i).x, blackestOnes.get(i).y, "Candidate Rocker");
    }*/

    for (SccProperty prop : sccProp.values()) {
      if (IdeaTester.mightBeRocker(img, width, height, scc, 0, prop)) {
        // Might be a rocker, we'll say it is if it can reach all LEDs.
        boolean reachAll = true;
        for (SccProperty prop2 : sccProp.values()) {
          if (prop2.isGreenLed || prop2.isBlueLed) {
            if (!IdeaTester.canReach(img, width, height, prop, prop2, scc, sccPropLarge)) {
              reachAll = false;
              break;
            }
          }
        }
        if (reachAll) {
          //IdeaTester.markPos(img, prop.x, prop.y, "Candidate Rocker");
          //IdeaTester.markPosRed(img, prop.x, prop.y, "Candidate Rocker");
          System.out.println("Candidate rocker post-reachability: " + prop.x + "," + prop.y);
          prop.isRockerSwitch = true;
        }
      }
      /*if (prop.isGreenLed || prop.isBlueLed) {
        // Look to the left.
        boolean isFirstFromLeft = true;
        for (int x = prop.xmin - 1; x >= 0; --x) {
          if (sccPropLarge.containsKey(scc[x][prop.xminy])) {
            break;
          }
          if (sccProp.containsKey(scc[x][prop.xminy]) &&
              sccProp.get(scc[x][prop.xminy]).x < prop.xmin) {
            isFirstFromLeft = false;
            break;
          }
        }
        if (isFirstFromLeft) {
          prop.isFirstFromLeft = true;
          IdeaTester.markPos(img, prop.x, prop.y, "IS FIRST FROM LEFT");
        }
      }*/
      if (prop.isGreenLed) {
        IdeaTester.markPos(img, prop.x, prop.y, "___");
      } else if (prop.isBlueLed) {
        IdeaTester.markPosBlue(img, prop.x, prop.y, "_ _");
      }
    }
    System.out.println("Num compoennts: " + (curComponent - 1));
    System.out.println("Num interesting components: " + numInterestingComponents);
    System.out.println("Num large components: " + sccPropLarge.size());

  }
}
