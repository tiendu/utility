package main

import (
    "fmt"
    "flag"
    "os"
    "bufio"
    "strings"
    "regexp"
    "sort"
    "math"
)

func main() {
    flag.Parse()
    filenames := flag.Args()
    var filename string
    for i := 0; i < len(filenames); i++ {
        filename = filenames[i]
        file, err := os.Open(filename)
        if err != nil {
            panic(err)
        }
        defer file.Close()
        header := regexp.MustCompile("^>")
        reader := bufio.NewReader(file)
        var lengths []int
        for {
            line, err := reader.ReadString('\n')
            line = strings.TrimSuffix(line, "\n")
            if err != nil {
                break
            }
            if header.MatchString(line) {
                continue
            } else {
                lengths = append(lengths, len(line))
            }
        }
        N50 := Nx(lengths, 50)
        N90 := Nx(lengths, 90)
        var sum_length, sum_squared_length float64
        for i := 0; i < len(lengths); i++ {
            sum_length += float64(lengths[i])
            sum_squared_length += math.Pow(float64(lengths[i]), 2)
        }
        auN := sum_squared_length / sum_length
        fmt.Printf("%s\n\t- No. of seqs: %d\n\t- Total size (bp): %.0f\n\t- Min: %d\n\t- Max: %d\n\t- auN: %.3f\n\t- N50: %d\n\t- N90: %d\n", filename, len(lengths), sum_length, MinInt(lengths), MaxInt(lengths), auN, N50, N90)
    }
}

func Nx(lengths []int, x int) int {
    var cumulative_sum, total_sum int
    for i := 0; i < len(lengths); i++ {
        total_sum += lengths[i]
    }
    sort.Ints(lengths)
    for i := 0; i < len(lengths); i++ {
        cumulative_sum += lengths[i]
        if cumulative_sum >= total_sum * (100 - x) / 100 {
            return lengths[i]
        }
    }
    return 0
} 

func MinInt(array []int) int {
    var min int
    min = array[0]
    for i := 0; i < len(array); i++ {
        if min > array[i] {
            min = array[i]
        }
    }
    return min
}

func MaxInt(array []int) int {
    var max int
    max = array[0]
    for i := 0; i < len(array); i++ {
        if max < array[i] {
            max = array[i]
        }
    }
    return max
}
