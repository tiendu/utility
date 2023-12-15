import java.util.concurrent.Callable
import java.util.concurrent.ExecutorService
import java.util.concurrent.ExecutorCompletionService
import java.util.concurrent.Executors
import java.util.concurrent.Future
import java.math.BigInteger

/**
 * Generate variants for a given DNA sequence by replacing ambiguous nucleotides with their possible values.
 * @param sequence The input DNA sequence.
 * @return A list of all possible variants generated from the input sequence.
 */
def generateVariants(sequence) {
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

/*
class PrimerMatcher {
    static List<CheckPrimersTask> generateTasks(String seq, List<String> forwardPrimers, List<String> reversePrimers) {
        List<CheckPrimersTask> tasks = []
        for (String fp : forwardPrimers) {
            for (String rp : reversePrimers) {
                tasks.add(new CheckPrimersTask(seq, fp, reverseComplement(rp)))
            }
        }
        return tasks
    }
    static String reverseComplement(String seq) {
        def complement = [
            'A': 'T',
            'T': 'A',
            'G': 'C',
            'C': 'G'
        ]
        def rcSeq = seq.reverse().collect { complement[it] ?: it }.join('')
        return rcSeq
    }
}
*/

class CheckPrimersTask implements Callable<Boolean> {
    String sequence
    String forwardPrimer
    String reversePrimer

    public CheckPrimersTask(String sequence, String forwardPrimer, String reversePrimer) {
        this.sequence = sequence
        this.forwardPrimer = forwardPrimer
        this.reversePrimer = reversePrimer
    }

    public Boolean call() throws Exception {
        if (sequence =~ forwardPrimer && sequence =~ reversePrimer) {
            return true
        }
        return false
    }
}

/**
 * Check if the given primers match the sequence using multithreading.
 * @param sequence The input DNA sequence.
 * @param forwardPrimer The forward primer.
 * @param reversePrimer The reverse primer.
 * @param numThreads The number of threads to use for parallel execution.
 * @return true if the primers match the sequence, false otherwise.
 */
def checkPrimersMultithreaded(sequence, forwardPrimer, reversePrimer, numThreads) {
    def delineatedForwardPrimers = generateVariants(forwardPrimer)
    def delineatedReversePrimers = generateVariants(reversePrimer)

    // Create tasks for each combination of forward and reverse primers
//     def tasks = PrimerMatcher.generateTasks(seq, delineatedForwardPrimer, delineatedReversePrimer)
    def tasks = []

    // Add tasks with every variants of primers
    for (fp in delineatedForwardPrimers) {
        for (rp in delineatedReversePrimers) {
            tasks.add(new CheckPrimersTask(sequence, fp, reverseComplement(rp)))
        }
    }

    // Create a thread pool with the specified number of threads
    ExecutorService executor = Executors.newFixedThreadPool(numThreads)

    try {
        // Create an ExecutorCompletionService using the executor
        def completionService = new ExecutorCompletionService<Boolean>(executor)
        
        // Submit tasks to the executor
        tasks.each { task ->
            completionService.submit(task)
        }

        // Process completed tasks as they become available
        int completedTasks = 0
        while (completedTasks < tasks.size()) {
            def completedTask = completionService.take().get()

            if (completedTask) {
                // Matching primers found
                println("Forward Primer: " + tasks[completedTasks].forwardPrimer)
                println("Reverse Primer: " + tasks[completedTasks].reversePrimer)
                return true
            }

            completedTasks++
        }
    } finally {
        executor.shutdown()
    }

    return false // No matching primers found
}

// Example usage
/*
def sequence = args[0]
def forwardPrimer = args[1]
def reversePrimer = args[2]
def numThreads = 8

def primersMatch = checkPrimersMultithreaded(sequence, forwardPrimer, reversePrimer, numThreads)

if (primersMatch) {
    println("Primers match the sequence.")
} else {
    println("Primers do not match the sequence.")
}
*/
