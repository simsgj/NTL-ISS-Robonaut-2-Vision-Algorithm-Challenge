import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.io.*;
import java.util.Arrays;

public class RobonautEyeTester {
  public static byte[] toBytes(int[] arr) {
    byte[] ret = new byte[arr.length * 4];
    for (int i = 0, j = 0; i < arr.length; ++i) {
      // Little-endian.
      ret[j++] = (byte) (arr[i] & 0x00ff);  // b
      ret[j++] = (byte) ((arr[i] >> 8) & 0x00ff);  // g
      ret[j++] = (byte) ((arr[i] >> 16) & 0x00ff);  // r
      ret[j++] = (byte) ((arr[i] >> 24) & 0x00ff);
    }
    return ret;
  }
    public static int[] imageToArray(String imageFile) throws Exception {
        printMessage("Reading image from " + imageFile + ".");
        if (imageFile.toLowerCase().endsWith("tif")) {
            FileInputStream fis = new FileInputStream(imageFile);
            byte[] contents = new byte[fis.available()];
            fis.read(contents);
            fis.close();

            int H, W;
            if (contents.length < 7000000) {
                H = 1200;
                W = 1600;
            } else {
                H = 2050;
                W = 2448;
            }

            int[] res = new int[H * W + 2];
            int pos = 8;
            res[0] = H;
            res[1] = W;
            for (int i = 2; i < res.length; i++) {
                int R = contents[pos++];
                int G = contents[pos++];
                int B = contents[pos++];
                if (R < 0) R += 256;
                if (G < 0) G += 256;
                if (B < 0) B += 256;
                res[i] = (R << 16) + (G << 8) + B;
            }

            return res;
        } else {
            BufferedImage img = ImageIO.read(new File(imageFile));

            int H = img.getHeight();
            int W = img.getWidth();
            int[] res = new int[H * W + 2];

            int pos = 0;
            res[pos++] = H;
            res[pos++] = W;

            byte[] pixels = ((DataBufferByte) img.getRaster().getDataBuffer()).getData();
            for (int i = 0; i < pixels.length; i += 3) {
                int R = pixels[i+2];
                int G = pixels[i+1];
                int B = pixels[i];
                if (R < 0) R += 256;
                if (G < 0) G += 256;
                if (B < 0) B += 256;
                res[pos++] = (R << 16) + (G << 8) + B;
            }

            return res;
        }
    }

    public void saveRawFile(int[] data, String fileName) throws Exception {
        printMessage("Saving raw image data to " + fileName + ".");
        PrintWriter pw = new PrintWriter(fileName);
        pw.println(data.length);
        for (int elm : data) {
            pw.println(elm);
        }
        pw.flush();
        pw.close();
    }

    static boolean debug = true;
    static String execCommand = null;
    static String convFile = null;
    static String folder = null, image = null;
    static String outputFile = null;
    static String visFile = null;
    static String modelFile = null;

    public static void printMessage(String s) {
        if (debug) {
            System.out.println(s);
        }
    }

    static final String _ = File.separator;

    static final int OBJECT_CNT = 22;
    static final String[] OBJECT_NAME = new String[] {
            "PANEL_POWER_SWITCH", "PANEL_POWER_COVER", "PANEL_POWER_LED", "A01_ROCKER_SWITCH", "A01_ROCKER_LED_TOP",
            "A01_ROCKER_LED_BOTTOM", "A02_LED_NUM_PAD_A1", "A02_LED_NUM_PAD_A2", "A02_LED_NUM_PAD_A3",
            "A02_LED_NUM_PAD_B1", "A02_LED_NUM_PAD_B2", "A02_LED_NUM_PAD_B3", "A02_LED_NUM_PAD_C1",
            "A02_LED_NUM_PAD_C2", "A02_LED_NUM_PAD_C3", "A03_TOGGLE", "A03_LED", "A04_TOGGLE", "A04_LED_TOP", "A04_LED_BOTTOM",
            "A05_TOGGLE", "A05_LED"
    };

    class TaskboardObject {
        String state;
        int leftX, leftY, rightX, rightY;

        public TaskboardObject(String state, int leftX, int leftY, int rightX, int rightY) {
            this.state = state;
            this.leftX = leftX;
            this.leftY = leftY;
            this.rightX = rightX;
            this.rightY = rightY;
        }

        public int getX(int LRFlag) {
            return LRFlag == 0 ? leftX : rightX;
        }

        public int getY(int LRFlag) {
            return LRFlag == 0 ? leftY : rightY;
        }
    }

    public void visualize(int[] imageData, TaskboardObject[] userAns, TaskboardObject[] modelAns, double[][] objectScore, String resFile, int LRFlag) throws Exception {
        int H = imageData[0], W = imageData[1];

        BufferedImage bi = new BufferedImage(W, H, BufferedImage.TYPE_INT_RGB);
        bi.setRGB(0, 0, W, H, imageData, 2, W);

        Graphics2D g = (Graphics2D)bi.getGraphics();

        g.setFont(new Font("Arial", Font.BOLD, 15));

        if (modelFile != null) {
            g.setColor(Color.BLUE);
            for (int i = 0; i < modelAns.length; i++) {
                TaskboardObject obj = modelAns[i];
                if (obj.state.equals("HIDDEN") || obj.getX(LRFlag) == -1 && obj.getY(LRFlag) == -1) {
                    continue;
                }
                g.fillOval(obj.getX(LRFlag) - 5, obj.getY(LRFlag) - 5, 11, 11);
                g.drawString("" + i, obj.getX(LRFlag), obj.getY(LRFlag) + 25);
            }
        }

        if (execCommand != null) {
            for (int i=0; i < userAns.length; i++) {
                TaskboardObject obj = userAns[i];
                if (obj.state.equals("HIDDEN") || obj.getX(LRFlag) == -1 && obj.getY(LRFlag) == -1) {
                    continue;
                }
                if (modelFile == null || objectScore[i][2] == 0.0) {
                    g.setColor(Color.RED);
                } else {
                    g.setColor(Color.GREEN);
                }
                g.fillOval(obj.getX(LRFlag) - 5, obj.getY(LRFlag) - 5, 11, 11);
                g.drawString("" + i, obj.getX(LRFlag), obj.getY(LRFlag) - 15);
            }
        }

        ImageIO.write(bi, "PNG", new File(resFile));
    }

    public void doExec() throws Exception {
      long loadStart = System.currentTimeMillis();
        if (!overrideThreshold) {
            if (folder.equals("ISS") || folder.equals("Lab")) {
                threshold = 15;
            } else {
                threshold = 10;
            }
        }
        if (modelFile != null && execCommand == null && visFile != null) {
            printMessage("Mode: visualize model answer.");
        } else if (modelFile == null && execCommand != null && visFile == null) {
            printMessage("Mode: execute your solution and validate return value.");
        } else if (modelFile != null && execCommand != null && visFile == null) {
            printMessage("Mode: execute your solution and calculate correctness score.");
        } else if (modelFile != null && execCommand != null && visFile != null) {
            printMessage("Mode: execute your solution, calculate correctness score, visualize your and model answers.");
        } else if (modelFile == null && execCommand != null && visFile != null) {
            printMessage("Mode: execute your solution, validate return value, visualize your answer.");
        } else {
            System.out.println("WARNING: nothing to do for this combination of arguments.");
            System.exit(0);
        }

        int[] leftImage = imageToArray(folder + _ + "LeftImage_" + image);
        int[] rightImage = imageToArray(folder + _ + "RightImage_" + image);

        String[] userAns;
      long loadEnd = System.currentTimeMillis();

        TaskboardObject[] userObj = new TaskboardObject[OBJECT_CNT];
        if (execCommand != null) {
            printMessage("Executing your solution: " + execCommand + ".");
            printMessage("Took " + (loadEnd - loadStart) + " ms to load images");

            Process solution = Runtime.getRuntime().exec(execCommand);

            BufferedReader reader = new BufferedReader(new InputStreamReader(solution.getInputStream()));
            DataOutputStream writer = new DataOutputStream(
                new BufferedOutputStream(solution.getOutputStream()));
            new ErrorStreamRedirector(solution.getErrorStream()).start();

            long printStart = System.currentTimeMillis();
            writer.writeInt(leftImage.length);
            writer.write(toBytes(leftImage));
            writer.write(toBytes(rightImage));
            writer.flush();
            long printEnd = System.currentTimeMillis();
            printMessage("Took " + (printEnd - printStart) + " ms to print out");

            userAns = new String[OBJECT_CNT];

            for (int i = 0; i < userAns.length; i++) {
                userAns[i] = reader.readLine();
            }

            solution.destroy();

            for (int i = 0; i < userAns.length; i++) {
                if (userAns[i].length() > 50 || !userAns[i].matches("[A-Z]{1,10},-?[0-9]{1,4},-?[0-9]{1,4},-?[0-9]{1,4},-?[0-9]{1,4}")) {
                    System.out.println("Your return value is not valid. Element " + i + " (0-based) is not properly formatted.");
                    return;
                }
            }

            if (modelFile == null && visFile == null) {
                System.out.println("Your return value is properly formatted. Here it is:");
                for (String elm : userAns) {
                    System.out.println(elm);
                }
                System.exit(0);
            }

            for (int i=0; i < OBJECT_CNT; i++) {
                String[] items = userAns[i].split(",");
                userObj[i] = new TaskboardObject(items[0], Integer.parseInt(items[1]), Integer.parseInt(items[2]), Integer.parseInt(items[3]), Integer.parseInt(items[4]));
            }
        }

        double[][] objectScore = new double[OBJECT_CNT][3];
        TaskboardObject[] modelObj = new TaskboardObject[OBJECT_CNT];
        if (modelFile != null) {
            printMessage("Extracting model answer from " + modelFile + ".");

            BufferedReader br = new BufferedReader(new FileReader(modelFile));
            br.readLine();
            boolean find = false;
            while (true) {
                String s = br.readLine();
                if (s == null) {
                    break;
                }
                String[] items = s.split(",");
                if (items[0].equals(folder) && items[1].equals(image)) {
                    find = true;
                    modelObj = new TaskboardObject[OBJECT_CNT];
                    for (int i = 0; i < OBJECT_CNT; i++) {
                        modelObj[i] = new TaskboardObject(items[5 * i + 2], Integer.parseInt(items[5 * i + 3]), Integer.parseInt(items[5 * i + 4]),
                                Integer.parseInt(items[5 * i + 5]), Integer.parseInt(items[5 * i + 6]));
                    }
                    break;
                }
            }
            br.close();

            if (!find) {
                System.out.println("ERROR: the model file does not contain an answer for folder = " + folder + " and image = " + image + ".");
                System.exit(0);
            }

            if (execCommand != null) {

                for (int i = 0; i < OBJECT_CNT; i++) {
                    if (!modelObj[i].state.equals("HIDDEN")) {
                        if (modelObj[i].leftX == -1 && modelObj[i].leftY == -1) {
                            objectScore[i][0] = 1;
                        } else {
                            int X1 = modelObj[i].leftX, X2 = userObj[i].leftX;
                            int Y1 = modelObj[i].leftY, Y2 = userObj[i].leftY;
                            double distanceLeft = Math.sqrt((X1 - X2) * (X1 - X2) + (Y1 - Y2) * (Y1 - Y2));
                            if (distanceLeft <= 2 * threshold) {
                                if (distanceLeft <= threshold) {
                                    objectScore[i][0] = 1;
                                } else {
                                    objectScore[i][0] = 1 - Math.pow((distanceLeft - threshold) / threshold, 1.5);
                                }
                            }
                        }

                        if(modelObj[i].rightX == -1 && modelObj[i].rightY == -1) {
                            objectScore[i][1] = 1;
                        } else {
                            int X1 = modelObj[i].rightX, X2 = userObj[i].rightX;
                            int Y1 = modelObj[i].rightY, Y2 = userObj[i].rightY;
                            double distanceRight = Math.sqrt((X1 - X2) * (X1 - X2) + (Y1 - Y2) * (Y1 - Y2));
                            if (distanceRight <= 2 * threshold) {
                                if (distanceRight <= threshold) {
                                    objectScore[i][1] = 1;
                                } else {
                                    objectScore[i][1] = 1 - Math.pow((distanceRight - threshold) / threshold, 1.5);
                                }
                            }
                        }

                        if (modelObj[i].state.equals(userObj[i].state)) {
                            objectScore[i][2] = 1;
                        }
                    } else {
                        Arrays.fill(objectScore[i], 1.0);
                    }
                }

                double score = 0;
                printMessage("Scores for each object (left image location, right image location, state):");
                for (int i=0; i < objectScore.length; i++) {
                    for (int j=0; j < objectScore[i].length; j++) {
                        score += objectScore[i][j];
                    }
                    printMessage(OBJECT_NAME[i] + ": " + objectScore[i][0] + " " + objectScore[i][1] + " " + objectScore[i][2]);
                }

                score /= OBJECT_CNT * 3;
                System.out.println("Score = " + score);
            }
        }

        if (visFile != null) {
            printMessage("Generating visualizations. File name: " + visFile + ".");
            visualize(leftImage, userObj, modelObj, objectScore, visFile + "_left.png", 0);
            visualize(rightImage, userObj, modelObj, objectScore, visFile + "_right.png", 1);
        }
    }

    public void doConvert() throws Exception {
        printMessage("Mode: convert an image to a sequence of integers.");
        if (outputFile == null) {
            System.out.println("ERROR: you need to specify where to save the result. Use -out <file name> for that.");
            System.exit(0);
        }
        int[] data = imageToArray(convFile);
        saveRawFile(data, outputFile);
    }

    static boolean overrideThreshold = false;
    static int threshold = 15;

    public static void main(String[] args) throws Exception {
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-conv")) {
                convFile = args[++i];
            } else if (args[i].equals("-out")) {
                outputFile = args[++i];
            } else if (args[i].equals("-folder")) {
                folder = args[++i];
            } else if (args[i].equals("-img")) {
                image = args[++i];
            } else if (args[i].equals("-exec")) {
                execCommand = args[++i];
            } else if (args[i].equals("-vis")) {
                visFile = args[++i];
            } else if (args[i].equals("-model")) {
                modelFile = args[++i];
            } else if (args[i].equals("-silent")) {
                debug = false;
            } else if (args[i].equals("-threshold")) {
                overrideThreshold = true;
                threshold = Integer.parseInt(args[++i]);
            } else {
                System.out.println("WARNING: unknown argument " + args[i] + ".");
            }
        }

        try {
            if (folder != null && image != null) {
                new RobonautEyeTester().doExec();
            } else if (convFile != null) {
                new RobonautEyeTester().doConvert();
            } else {
                System.out.println("WARNING: nothing to do for this combination of arguments.");
            }
        } catch (Exception e) {
            System.out.println("FAILURE: " + e.getMessage());
            e.printStackTrace();
        }
    }

    class ErrorStreamRedirector extends Thread {
        public BufferedReader reader;

        public ErrorStreamRedirector(InputStream is) {
            reader = new BufferedReader(new InputStreamReader(is));
        }

        public void run() {
            while (true) {
                String s;
                try {
                    s = reader.readLine();
                } catch (Exception e) {
                    // e.printStackTrace();
                    return;
                }
                if (s == null) {
                    break;
                }
                System.out.println(s);
            }
        }
    }
}
