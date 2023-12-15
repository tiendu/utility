import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.nio.charset.StandardCharsets;
import java.util.zip.GZIPInputStream;

public class FastaCounter {
    public static void main(String[] args) {
        List<String> files = new ArrayList<>();
        int desiredCores = -1;

        // Parse command-line arguments
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-i") || args[i].equals("--input")) {
                // Extract input files
                for (int j = i + 1; j < args.length; j++) {
                    if (args[j].startsWith("-")) {
                        break;
                    }
                    files.add(args[j]);
                }
            } else if (args[i].equals("-n") || args[i].equals("--cores")) {
                // Extract number of cores
                if (i + 1 < args.length) {
                    try {
                        desiredCores = Integer.parseInt(args[i + 1]);
                    } catch (NumberFormatException e) {
                        System.out.println("Error: Invalid number of cores.");
                        return;
                    }
                }
            }
        }

        // Check if required arguments are provided
        if (files.isEmpty() || desiredCores == -1) {
            System.out.println("Usage: java LineCounter -i <file1> <file2> ... -n <cores>");
            return;
        }

        // Get the number of available cores
        int availableCores = Runtime.getRuntime().availableProcessors();

        // Validate the desired number of cores
        if (desiredCores <= 0 || desiredCores > availableCores) {
            System.out.println("Error: Invalid number of cores. Available cores: " + availableCores);
            return;
        }

        // Create a thread pool with the desired number of cores
        ExecutorService executor = Executors.newFixedThreadPool(desiredCores);

        // Create a list to hold the Future objects
        List<Future<StringBuilder>> futures = new ArrayList<>();

        // Submit tasks to the thread pool
        for (String file : files) {
            Future<StringBuilder> future = executor.submit(() -> countLines(file));
            futures.add(future);
        }

        // Print the line count for each file
        for (int i = 0; i < files.size(); i++) {
            String file = files.get(i);
            StringBuilder result = null;

            try {
                result = futures.get(i).get();
            } catch (InterruptedException | ExecutionException e) {
                e.printStackTrace();
            }
            System.out.println("Output for " + file + ":\n" + result);
        }

        // Shutdown the thread pool
        executor.shutdown();
    }

    private static StringBuilder countLines(String file) {
        int lineCount = 0;
        StringBuilder output = new StringBuilder();
        try {
            BufferedReader reader;
            if (file.endsWith(".gz")) {
                // Handle gzipped file
                InputStream fileStream = new FileInputStream(file);
                InputStream gzipStream = new GZIPInputStream(fileStream);
                reader = new BufferedReader(new InputStreamReader(gzipStream, StandardCharsets.UTF_8));
            } else {
                // Handle regular text file
                reader = new BufferedReader(new FileReader(file, StandardCharsets.UTF_8));
            }

            String line;
            StringBuilder sequence = new StringBuilder();
            String sequenceId = "";

            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (sequence.length() > 0 ) {
                        output.append(file).append(" ").append(sequenceId).append(" ").append(lineCount).append("\n");
                        sequence.setLength(0);
                        lineCount = 0;
                    }
                    sequenceId = line.substring(1);
                } else {
                    sequence.append(line);
                    lineCount += line.length();
                }
            }

            // Add the length of the last sequence if present
            if (sequence.length() > 0) {
                output.append(file).append(" ").append(sequenceId).append(" ").append(lineCount).append("\n");
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return output;
    }
}
