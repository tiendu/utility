import java.io.BufferedReader
import java.io.File
import java.io.FileReader
import java.io.FileInputStream
import java.io.InputStreamReader
import java.io.IOException
import java.util.zip.GZIPInputStream
import java.util.List
import java.util.ArrayList

/**
 * Extracts the lines that come after a line starting with '@' from a given file.
 *
 * @param inputFile The input file to read.
 * @return The lines that come after a line starting with '@'.
 */
List<String> extractSequencesFromFastQ(File inputFile) {
    List<String> sequences = new ArrayList<>()
    try {
        BufferedReader reader
        // Check if the file is gzipped
        if (inputFile.getName().endsWith(".gz")) {
            // Create a buffered reader for gzipped file
            reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFile))))
        } else {
            // Create a buffered reader for regular file
            reader = new BufferedReader(new FileReader(inputFile))
        }

        String line = reader.readLine()
        int count = 1
        while (line != null) {
            if (count.mod(4) == 2) {
                def sequence = line // Read the next line after the line starting with '@' and divisible by 4
                if (sequence != null) {
                    sequences.add(sequence.trim()) // Add the trimmed line
                }
            }
            line = reader.readLine()
            count++
        }

        reader.close()
    } catch (IOException e) {
        e.printStackTrace()
    }

    return sequences // Return the lines that come after a line starting with '@'
}

// Example usage

File inputFile = new File(args[0])
List<String> nextLines = extractSequencesFromFastQ(inputFile)
if (!nextLines.isEmpty()) {
    for (String line : nextLines) {
        System.out.println(line)
    }
}

