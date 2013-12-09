import java.util.*;
import java.io.*;

public class TestAll {
  static class TestEntry {
    String dir;
    String filename;
    public TestEntry(String dir, String filename) {
      this.dir = dir;
      this.filename = filename;
    }
  }

  public static void main(String[] args) throws Exception {
    long startTime = System.currentTimeMillis();
    String[] folders = { "ISS", "Lab", "Lab2", "Lab3", "Sim" };

    final LinkedList<TestEntry> entries = new LinkedList<TestEntry>();
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
        entries.add(new TestEntry(folder, img));
      }
    }

    Thread[] workers = new Thread[8];
    final Vector<String> finalOutput = new Vector<String>();
    for (int i = 0; i < workers.length; ++i) {
      workers[i] = new Thread() {
        @Override
        public void run() {
          try {
            while (true) {
              TestEntry entry;
              synchronized (entries) {
                if (entries.size() == 0) {
                  break;
                }
                entry = entries.poll();
              }

              Process proc = Runtime.getRuntime().exec(
                  "java RobonautEyeTester -folder " + entry.dir + " -img " + entry.filename +
                  " -exec \"./a.exe\" -model model.csv");
              Scanner scan = new Scanner(proc.getInputStream());
              while (scan.hasNextLine()) {
                String line = scan.nextLine();
                if (line.startsWith("Score =")) {
                  finalOutput.add(entry.dir + "/" + entry.filename + ": " + line);
                }
              }
              proc.waitFor();
            }
          } catch (Exception e) {
            e.printStackTrace();
          }
        }
      };
      workers[i].start();
    }
    for (int i = 0; i < workers.length; ++i) {
      workers[i].join();
    }
    double scoreSum = 0;
    Collections.sort(finalOutput);
    for (String s : finalOutput) {
      double score = Double.parseDouble(s.split("=")[1].trim());
      scoreSum += score;
      System.out.println(s);
    }
    long endTime = System.currentTimeMillis();
    System.out.println("Everything took " + (endTime - startTime) + " ms");
    System.out.println("Avg score: " + (scoreSum / finalOutput.size() * 1000000));
  }
}
