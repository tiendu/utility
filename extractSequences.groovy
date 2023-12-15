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
 * Extracts the sequences from a FASTA or FASTQ file.
 *
 * @param inputPath The input file path to read.
 * @return The lines that come after a line starting with '@' or '>'.
 */
List<String> extractSequencesFromFile(String inputPath) {
    File inputFile = new File(inputPath)
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
        StringBuilder sequenceBuilder = new StringBuilder()
        boolean isFasta = line.startsWith(">")
        boolean isFastq = line.startsWith("@")
        char lineStartChar = isFasta ? '>' : '@'
        int count = 1
        while (line != null) {
            if (line.charAt(0) == lineStartChar) {
                // Line starting character indicates the start of a new sequence
                // Add the previous sequence to the list
                if (sequenceBuilder.length() > 0) {
                    sequences.add(sequenceBuilder.toString())
                    sequenceBuilder = new StringBuilder()
                }
            } else if (isFastq && count.mod(4) == 2) {
                // Append the line to the current sequence if it's a FASTQ file and the line is at position 2 (1-based index)
                sequenceBuilder.append(line.trim())
            } else if (isFasta) {
                // Append the line to the current sequence if it's a FASTA file
                sequenceBuilder.append(line.trim())
            }
            line = reader.readLine()
            count++
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
/*
File inputFile = new File(args[0])
List<String> sequences = extractSequencesFromFile(inputFile)
if (!sequences.isEmpty()) {
    for (String sequence : sequences) {
        println sequence
    }
}*/
