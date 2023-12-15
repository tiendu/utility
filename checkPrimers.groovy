import java.math.BigInteger

/**
 * Generate variants for a given sequence by replacing ambiguous nucleotides with their possible values.
 * @param s The input sequence.
 * @return A list of all possible variants generated from the input sequence.
 */
def generateVariants(s) {
    def hash = [
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
    def solutions = BigInteger.ONE
    def a = s.split('')
    def m = a.size()
    def list = []

    (0..<m).each { i ->
        def b = hash[a[i]] ?: []
        solutions *= b.size()
        list.add(b)
    }

    def results = []

    (0..<solutions.intValue()).each { i ->
        def result = ''
        def idx = BigInteger.valueOf(i)
        (0..<m).each { j ->
            // Get the sublist of possible nucleotides for the current position
            def sublist = list[j] ?: []

            // Get the index within the sublist based on the current iteration and sublist size
            def sublistIndex = idx.mod(BigInteger.valueOf(sublist.size())).intValue()

            // Append the selected nucleotide from the sublist to the result
            result += sublist[sublistIndex]

            // Update the index by dividing it by the sublist size to move to the next position
            idx = idx.divide(BigInteger.valueOf(sublist.size()))
        }

        // Add the generated variant to the list of results
        results.add(result)
    }

    return results
}

/**
 * Get the reverse complement of a given DNA sequence.
 * @param seq The input DNA sequence.
 * @return The reverse complement of the input sequence.
 */
def reverseComplement(seq) {
    def complement = [
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
        ]

    def rcSeq = seq.reverse().collect { complement[it] ?: it }.join('')
    return rcSeq
}

/**
 * Check if the given primers match the sequence.
 * @param seq The input sequence.
 * @param forwardPrimer The forward primer.
 * @param reversePrimer The reverse primer.
 * @return true if the primers match the sequence, false otherwise.
 */
def checkPrimers(seq, forwardPrimer, reversePrimer) {
    def delineatedForwardPrimer = generateVariants(forwardPrimer)
    def delineatedReversePrimer = generateVariants(reversePrimer)
    
    // Iterate over all possible variants of the forward primer
    for (fp in delineatedForwardPrimer) {
        def forwardRegex = fp.replaceAll("[RYWSKMBDHVN]", ".")
         // Check if the forward primer variant matches the sequence
        if (seq =~ forwardRegex) {
            // Iterate over all possible variants of the reverse primer
            for (rp in delineatedReversePrimer) {
                def rcRp = reverseComplement(rp)
                def reverseRegex = rcRp.replaceAll("[RYWSKMBDHVN]", ".")
                // Check if the reverse complement of the reverse primer variant matches the sequence
                if (seq =~ reverseRegex) {
                    println "${seq} ${forwardRegex} ${reverseRegex}"
                    return true
                }
            }
        }
    }
//     for (fp in delineatedForwardPrimer) {
//         if (seq.startsWith(fp)) {
//             for (rp in delineatedReversePrimer) {
//                 def rcRp = reverseComplement(rp)
//                 if (seq.endsWith(rcRp)) {
//                     return true
//                 }
//             }
//         }
//     }

    return false // No matching primers found
}

// Example usage
/*
def seq = args[0]
def primer1 = args[1]
def primer2 = args[2]

def primersMatch = checkPrimers(seq, primer1, primer2)

if (primersMatch) {
    println("Primers match the sequence.")
} else {
    println("Primers do not match the sequence.")
}
*/
