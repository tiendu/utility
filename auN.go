package main

import (
    "bufio"
    "compress/gzip"
    "errors"
    "flag"
    "fmt"
    "io"
    "math"
    "os"
    "runtime"
    "sort"
    "strings"
    "sync"
)

// FileType distinguishes between FASTA and FASTQ files.
type FileType int

const (
    FASTA FileType = iota
    FASTQ
)

// Seq holds a sequence record.
type Seq struct {
    ID       string
    Sequence string
    Quality  *string // nil for FASTA
}

// arrayFlags is a custom flag type for capturing multiple input files.
type arrayFlags []string

func (i *arrayFlags) String() string {
    return strings.Join(*i, ",")
}

// Set supports comma-separated values as well as repeated flags.
func (i *arrayFlags) Set(value string) error {
    parts := strings.Split(value, ",")
    for _, part := range parts {
        trimmed := strings.TrimSpace(part)
        if trimmed != "" {
            *i = append(*i, trimmed)
        }
    }
    return nil
}

// FileTypeFromFilename returns the file type based on the filename extension.
func FileTypeFromFilename(filename string) (FileType, error) {
    fastqExts := []string{".fastq.gz", ".fq.gz", ".fastq", ".fq"}
    fastaExts := []string{
        ".fasta", ".fa", ".faa", ".fna",
        ".fasta.gz", ".fa.gz", ".faa.gz", ".fna.gz",
    }

    for _, ext := range fastqExts {
        if strings.HasSuffix(filename, ext) {
            return FASTQ, nil
        }
    }
    for _, ext := range fastaExts {
        if strings.HasSuffix(filename, ext) {
            return FASTA, nil
        }
    }
    return 0, errors.New(fmt.Sprintf("Unrecognized file extension for '%s'", filename))
}

// readLine is a helper to read a single line from the scanner.
func readLine(scanner *bufio.Scanner) (string, error) {
    if scanner.Scan() {
        return scanner.Text(), nil
    }
    if err := scanner.Err(); err != nil {
        return "", err
    }
    return "", io.EOF
}

// readSequences reads sequences from the file at filePath using the given fileType.
// It supports both uncompressed and gzip-compressed files.
func readSequences(filePath string, fileType FileType) ([]Seq, error) {
    file, err := os.Open(filePath)
    if err != nil {
        return nil, err
    }
    defer file.Close()

    var reader io.Reader
    if strings.HasSuffix(filePath, ".gz") {
        gzReader, err := gzip.NewReader(file)
        if err != nil {
            return nil, err
        }
        defer gzReader.Close()
        reader = gzReader
    } else {
        reader = file
    }

    bufReader := bufio.NewReader(reader)
    scanner := bufio.NewScanner(bufReader)
    var seqs []Seq

    switch fileType {
    case FASTQ:
        // FASTQ records are 4 lines: header, sequence, plus, quality.
        for {
            header, err := readLine(scanner)
            if err == io.EOF {
                break
            } else if err != nil {
                return nil, err
            }
            if !strings.HasPrefix(header, "@") {
                return nil, errors.New("Invalid FASTQ format: header does not start with '@'")
            }

            seqLine, err := readLine(scanner)
            if err != nil {
                return nil, err
            }
            // Skip the '+' line.
            _, err = readLine(scanner)
            if err != nil {
                return nil, err
            }
            qualLine, err := readLine(scanner)
            if err != nil {
                return nil, err
            }
            id := strings.TrimPrefix(header, "@")
            quality := qualLine
            seqs = append(seqs, Seq{
                ID:       id,
                Sequence: seqLine,
                Quality:  &quality,
            })
        }
    case FASTA:
        var currentSeq *Seq = nil
        for scanner.Scan() {
            line := strings.TrimSpace(scanner.Text())
            if line == "" {
                continue
            }
            if strings.HasPrefix(line, ">") {
                if currentSeq != nil {
                    seqs = append(seqs, *currentSeq)
                }
                currentSeq = &Seq{
                    ID:       strings.TrimPrefix(line, ">"),
                    Sequence: "",
                    Quality:  nil,
                }
            } else {
                if currentSeq == nil {
                    // Skip lines before any header
                    continue
                }
                currentSeq.Sequence += strings.TrimSpace(line)
            }
        }
        if currentSeq != nil {
            seqs = append(seqs, *currentSeq)
        }
        if err := scanner.Err(); err != nil {
            return nil, err
        }
    }

    return seqs, nil
}

// Nx computes the Nx metric for a slice of sequence lengths.
func Nx(lengths []int, x int) int {
    totalSum := 0
    for _, l := range lengths {
        totalSum += l
    }
    sort.Ints(lengths)
    cumulativeSum := 0
    for _, l := range lengths {
        cumulativeSum += l
        if cumulativeSum >= totalSum*(100-x)/100 {
            return l
        }
    }
    return 0
}

// MinInt returns the minimum integer in a slice.
func MinInt(array []int) int {
    if len(array) == 0 {
        return 0
    }
    min := array[0]
    for _, v := range array {
        if v < min {
            min = v
        }
    }
    return min
}

// MaxInt returns the maximum integer in a slice.
func MaxInt(array []int) int {
    if len(array) == 0 {
        return 0
    }
    max := array[0]
    for _, v := range array {
        if v > max {
            max = v
        }
    }
    return max
}

// processFile determines the file type, reads sequences, computes statistics,
// and prints the results.
func processFile(filename string) {
    fileType, err := FileTypeFromFilename(filename)
    if err != nil {
        fmt.Fprintf(os.Stderr, "Error determining file type for %s: %v\n", filename, err)
        return
    }

    seqs, err := readSequences(filename, fileType)
    if err != nil {
        fmt.Fprintf(os.Stderr, "Error reading sequences from %s: %v\n", filename, err)
        return
    }

    var lengths []int
    for _, seq := range seqs {
        lengths = append(lengths, len(seq.Sequence))
    }
    if len(lengths) == 0 {
        fmt.Printf("%s: No sequences found\n", filename)
        return
    }
    N50 := Nx(lengths, 50)
    N90 := Nx(lengths, 90)
    var sumLength, sumSquaredLength float64
    for _, l := range lengths {
        sumLength += float64(l)
        sumSquaredLength += math.Pow(float64(l), 2)
    }
    auN := sumSquaredLength / sumLength

    fmt.Printf("%s\n\t- No. of seqs: %d\n\t- Total size (bp): %.0f\n\t- Min: %d\n\t- Max: %d\n\t- auN: %.3f\n\t- N50: %d\n\t- N90: %d\n",
        filename, len(lengths), sumLength, MinInt(lengths), MaxInt(lengths), auN, N50, N90)
}

func main() {
    // Define command-line flags.
    var inputFiles arrayFlags
    flag.Var(&inputFiles, "input", "Input FASTA/FASTQ files (provide multiple times or as a comma-separated list)")
    threadsFlag := flag.Int("threads", 4, "number of threads to use (default is 4)")
    flag.Parse()

    // Combine files provided via --input and positional arguments.
    if len(flag.Args()) > 0 {
        inputFiles = append(inputFiles, flag.Args()...)
    }

    if len(inputFiles) == 0 {
        fmt.Fprintln(os.Stderr, "No input files provided. Use --input to specify files or provide them as positional arguments.")
        os.Exit(1)
    }

    // Determine the number of threads: use the minimum between the user-provided thread count and available CPU cores.
    availableThreads := runtime.NumCPU()
    threads := *threadsFlag
    if threads > availableThreads {
        threads = availableThreads
    }

    sem := make(chan struct{}, threads)
    var wg sync.WaitGroup

    for _, filename := range inputFiles {
        wg.Add(1)
        sem <- struct{}{}
        go func(fn string) {
            defer wg.Done()
            processFile(fn)
            <-sem
        }(filename)
    }
    wg.Wait()
}
