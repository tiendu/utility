// Import necessary libraries
import java.util.concurrent.Callable
import java.util.concurrent.ExecutorService
import java.util.concurrent.CompletionService
import java.util.concurrent.ExecutorCompletionService
import java.util.concurrent.Executors
import java.util.concurrent.Future
import java.util.stream.Stream
import java.util.zip.GZIPInputStream
import java.util.List
import java.util.ArrayList
import java.io.BufferedReader
import java.io.File
import java.io.FileReader
import java.io.FileInputStream
import java.io.InputStreamReader
import java.io.IOException


/**
 * Generate variants for a given DNA sequence by replacing ambiguous nucleotides with their possible values.
 * @param sequence The input DNA sequence.
 * @return A list of all possible variants generated from the input sequence.
 */
def generateVariants(sequence) {
    // Define nucleotide variants
    def nucleotideVariants = [
        'A': ['A'],
        'T': ['T'],
        'G': ['G'],
        'C': ['C'],
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'K': ['G', 'T'],
        'M': ['A', 'C'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'T', 'G', 'C']
    ]
    def totalVariants = BigInteger.ONE
    def nucleotides = sequence.split('')
    def length = nucleotides.size()
    def variantLists = []

    (0..<length).each { i ->
        // Get the list of nucleotide variants for each position in the sequence
        def variants = nucleotideVariants[nucleotides[i]] ?: []
        totalVariants *= variants.size()
        variantLists.add(variants)
    }

    def results = []

    (0..<totalVariants.intValue()).each { i ->
        def variant = ''
        def index = BigInteger.valueOf(i)
        (0..<length).each { j ->
            // Get the list of possible nucleotide variants for the current position
            def variants = variantLists[j] ?: []

            // Get the index within the list of variants based on the current iteration and the list size
            def variantIndex = index.mod(BigInteger.valueOf(variants.size())).intValue()

            // Append the selected nucleotide variant from the list to the result
            variant += variants[variantIndex]

            // Update the index by dividing it by the list size to move to the next position
            index = index.divide(BigInteger.valueOf(variants.size()))
        }

        // Add the generated variant to the list of results
        results.add(variant)
    }

    return results
}

/**
 * Get the reverse complement of a given DNA sequence.
 * @param sequence The input DNA sequence.
 * @return The reverse complement of the input sequence.
 */
def reverseComplement(sequence) {
    def complement = [
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
    ]

    def reverseComplementedSeq = sequence.reverse().collect { complement[it] ?: it }.join('')
    return reverseComplementedSeq
}

/**
 * CheckPrimersTask class represents a task for checking if a sequence contains the given forward and reverse primers.
 */
class CheckPrimersTask implements Callable<String> {
    String sequence
    List<String> forwardPrimerVariants
    List<String> reversePrimerVariants

    public CheckPrimersTask(String sequence, List<String> forwardPrimerVariants, List<String> reversePrimerVariants) {
        this.sequence = sequence
        this.forwardPrimerVariants = forwardPrimerVariants
        this.reversePrimerVariants = reversePrimerVariants
    }

    public String call() throws Exception {
        // Iterate over all forward and reverse primer combinations
        for (String forwardPrimer : forwardPrimerVariants) {
            for (String reversePrimer : reversePrimerVariants) {
                // Check if the sequence contains the forward and reverse primers
                if (sequence.contains(forwardPrimer) && sequence.contains(reversePrimer)) {
                    // If match found, return the forward and reverse primer combination
                    return forwardPrimer + " " + reversePrimer
                }
            }
        }
        // If no match found, return null
        return null
    }
}

/**
 * Check if the primers match the given sequence using multiple threads.
 * @param sequence The DNA sequence to check against.
 * @param forwardPrimer The forward primer sequence.
 * @param reversePrimer The reverse primer sequence.
 * @param numThreads The number of threads to use for parallel processing.
 * @return The number of matching primer variants found.
 */
def checkPrimersMultithreaded(sequence, forwardPrimer, reversePrimer, numThreads) {
    // Generate variants for forward and reverse primers
    def delineatedForwardPrimers = generateVariants(forwardPrimer)
    def delineatedReversePrimers = []
    generateVariants(reversePrimer).each { rp ->
        delineatedReversePrimers.add(reverseComplement(rp))
    }

    List<String> matchingPrimerVariants = [] // Store the matching primer variants

    // Create a thread pool with a fixed number of threads
    ExecutorService executor = Executors.newFixedThreadPool(numThreads)
    CompletionService<String> completionService = new ExecutorCompletionService<>(executor)

    // Submit tasks for each combination of forward and reverse primer variants
    completionService.submit(new CheckPrimersTask(sequence, delineatedForwardPrimers, delineatedReversePrimers))

    try {
        // Process the completed tasks
        Future<String> future = completionService.take()
        String result = future.get()
        if (result) {
            matchingPrimerVariants.add(result)
        }
    } catch (Exception e) {
        e.printStackTrace()
    } finally {
        executor.shutdown()
    }

    // Print the matching primer variants
    matchingPrimerVariants.each { matchedVariant ->
        String[] primers = matchedVariant.split(' ')
        def rcRp = reverseComplement(primers[1])
        println("Sequence: ${sequence}")
        println("Forward Primer: ${primers[0]}")
        println("Reverse Primer: ${rcRp}")
        println()
    }

    return matchingPrimerVariants.size()
}

/**
 * Extract DNA sequences from a file.
 * @param inputPath The path to the input file.
 * @return A list of DNA sequences extracted from the file.
 */
def extractSequencesFromFile(inputPath) {
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

/**
 * Process a file with primers and count the matching sequences.
 * @param inputPath The path to the input file.
 * @param forwardPrimer The forward primer sequence.
 * @param reversePrimer The reverse primer sequence.
 * @param numThreads The number of threads to use for parallel processing.
 */
def processFileWithPrimers(String inputPath, String forwardPrimer, String reversePrimer, int numThreads) {
    List<String> sequences = extractSequencesFromFile(inputPath)

    int matchingSequenceCount = 0

    sequences.each { sequence ->
        int count = checkPrimersMultithreaded(sequence, forwardPrimer, reversePrimer, numThreads)
        if (count > 0) {
            matchingSequenceCount += count
        }
    }

    println("Total matching sequences: ${matchingSequenceCount}")
}

// Example usage

String inputPath = args[0]
String forwardPrimer = args[1]
String reversePrimer = args[2]
int numThreads = Integer.parseInt(args[3])

processFileWithPrimers(inputPath, forwardPrimer, reversePrimer, numThreads)
