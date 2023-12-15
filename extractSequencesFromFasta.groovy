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
 * Extracts the sequences from a FASTA file.
 *
 * @param inputFile The input FASTA file to read.
 * @return The lines that come after a line starting with '@'.
 */
List<String> extractSequencesFromFasta(String filePath) {
    List<String> sequences = new ArrayList<>()
    try {
        BufferedReader reader
        File inputFile = new File(filePath)
        // Check if the file is gzipped
        if (inputFile.getName().endsWith(".gz")) {
            // Create a buffered reader for gzipped file
            reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFile))))
        } else {
            // Create a buffered reader for regular file
            reader = new BufferedReader(new FileReader(inputFile))
        }

        String line = reader.readLine()
        StringBuilder sequenceBuilder = new StringBuilder()
        while (line != null) {
            if (line.startsWith(">")) {
                // Check if there's a previous sequence and add it to the list
                if (sequenceBuilder.length() > 0) {
                    sequences.add(sequenceBuilder.toString())
                    sequenceBuilder = new StringBuilder()
                }
            } else {
                // Append the line to the sequence builder
                sequenceBuilder.append(line.trim())
            }
            line = reader.readLine()
        }

        // Add the last sequence to the list
        if (sequenceBuilder.length() > 0) {
            sequences.add(sequenceBuilder.toString())
        }

        reader.close()
    } catch (IOException e) {
        e.printStackTrace()
    }

    return sequences
}

// Example usage

String filePath = args[0]
List<String> sequences = extractSequencesFromFasta(filePath)
if (!sequences.isEmpty()) {
    for (String sequence : sequences) {
        println sequence
    }
}
