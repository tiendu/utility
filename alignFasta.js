const fs = require('fs');

// Represents a set of scores for sequence alignment
class ScoreSet {
    constructor() {
        this.scoringMatrix = null;
        this.gap = null;
        this.beginGap = null;
        this.endGap = null;
        this.useBeginGapTop = true;
        this.useBeginGapLeft = true;
        this.useEndGapBottom = true;
        this.useEndGapRight = true;
    }

    // Get the score of aligning two residues
    getScore(residue1, residue2) {
        return this.scoringMatrix.getScore(residue1, residue2);
    }

    // Set parameters for the score set
    setScoreSetParams(scoringMatrix, gapPenalty, beginGapPenalty, endGapPenalty) {
        this.scoringMatrix = scoringMatrix;
        this.gap = gapPenalty;
        this.beginGap = beginGapPenalty;
        this.endGap = endGapPenalty;
    }
}

// Represents a scoring matrix for sequence alignment
class ScoringMatrix {
    constructor() {
        this.mismatch = null;
        this.match = null;
    }

    // Get the score for aligning two residues
    getScore(residue1, residue2) {
        residue1 = residue1.toLowerCase();
        residue2 = residue2.toLowerCase();
        if (residue1 === residue2) {
            return this.match;
        } else {
            return this.mismatch;
        }
    }
}

// Represents an identity scoring matrix
class Identity extends ScoringMatrix {
    setMismatch(mismatchScore) {
        this.mismatch = mismatchScore;
    }

    setMatch(matchScore) {
        this.match = matchScore;
    }
}

// Represents a node in dynamic programming matrix
class Node {
    constructor() {
        this.value;
        this.tracebackI;
        this.tracebackJ;
    }
}

// Represents a quadratically aligning pair of sequences
class AlignPairQuad {
    constructor() {
        this.M;
        this.N;
        this.scoreSet;
        this.nodes;
        this.alignedM;
        this.alignedN;
        this.score;
    }

    initializeMatrix(sequenceOne, sequenceTwo, scoreSet) {
        this.scoreSet = scoreSet;

        this.M = sequenceOne;
        this.N = sequenceTwo;
        this.score = 0;

        // Create a two-dimensional array of nodes
        this.nodes = new Array(this.M.length + 1);

        // Row i
        for (let i = 0; i < this.nodes.length; i++) {
            this.nodes[i] = new Array(this.N.length + 1);
            // Column j
            for (let j = 0; j < this.nodes[i].length; j++) {
                this.nodes[i][j] = new Node();
            }
        }

        this.nodes[0][0].value = 0;

        // i rows
        for (let i = 1; i < this.nodes.length; i++) {
            if (this.scoreSet.useBeginGapLeft) {
                this.nodes[i][0].value = this.nodes[i - 1][0].value - this.scoreSet.beginGap;
            } else {
                this.nodes[i][0].value = this.nodes[i - 1][0].value - this.scoreSet.gap;
            }
            this.nodes[i][0].tracebackI = i - 1;
            this.nodes[i][0].tracebackJ = 0;
        }

        // j columns
        for (let j = 1; j < this.nodes[0].length; j++) {
            if (this.scoreSet.useBeginGapTop) {
                this.nodes[0][j].value = this.nodes[0][j - 1].value - this.scoreSet.beginGap;
            } else {
                this.nodes[0][j].value = this.nodes[0][j - 1].value - this.scoreSet.gap;
            }
            this.nodes[0][j].tracebackI = 0;
            this.nodes[0][j].tracebackJ = j - 1;
        }
    }

    dumpMatrix() {
        console.log("Dynamic programming matrix i=" + this.nodes.length + " and j=" + this.nodes[0].length);
        for (let i = 0; i < this.nodes.length; i++) {
            let row = '';
            for (let j = 0; j < this.nodes[i].length; j++) {
                const traceI = this.nodes[i][j].tracebackI !== undefined ? this.nodes[i][j].tracebackI : 'u';
                const traceJ = this.nodes[i][j].tracebackJ !== undefined ? this.nodes[i][j].tracebackJ : 'u';
                row += `(${i},${j})[${traceI},${traceJ}]=${this.nodes[i][j].value} `;
            }
            console.log(row);
        }
    }

    fillMatrix() {
        //i rows
        for (var i = 1; i < this.nodes.length; i++) {
            //j columns
            for (var j = 1; j < this.nodes[0].length; j++) {
            var a;
            var b;
            var c;
        
            //handle end gaps here
        
            if (i == this.nodes.length - 1 && j == this.nodes[0].length - 1) {
                if (this.scoreSet.useEndGapRight) {
                a = this.nodes[i - 1][j].value - this.scoreSet.endGap;
                } else {
                a = this.nodes[i - 1][j].value - this.scoreSet.gap;
                }
        
                if (this.scoreSet.useEndGapBottom) {
                b = this.nodes[i][j - 1].value - this.scoreSet.endGap;
                } else {
                b = this.nodes[i][j - 1].value - this.scoreSet.gap;
                }
            } else if (i == this.nodes.length - 1) {
                a = this.nodes[i - 1][j].value - this.scoreSet.gap;
                if (this.scoreSet.useEndGapBottom) {
                b = this.nodes[i][j - 1].value - this.scoreSet.endGap;
                } else {
                b = this.nodes[i][j - 1].value - this.scoreSet.gap;
                }
            } else if (j == this.nodes[0].length - 1) {
                if (this.scoreSet.useEndGapRight) {
                a = this.nodes[i - 1][j].value - this.scoreSet.endGap;
                } else {
                a = this.nodes[i - 1][j].value - this.scoreSet.gap;
                }
                b = this.nodes[i][j - 1].value - this.scoreSet.gap;
            } else {
                a = this.nodes[i - 1][j].value - this.scoreSet.gap;
                b = this.nodes[i][j - 1].value - this.scoreSet.gap;
            }
        
            c =
                this.nodes[i - 1][j - 1].value +
                this.scoreSet.getScore(this.M[i - 1], this.N[j - 1]);
        
            if (a >= b && a >= c) {
                this.nodes[i][j].value = a;
                this.nodes[i][j].tracebackI = i - 1;
                this.nodes[i][j].tracebackJ = j;
            } else if (b >= c && b >= a) {
                this.nodes[i][j].value = b;
                this.nodes[i][j].tracebackI = i;
                this.nodes[i][j].tracebackJ = j - 1;
            } else {
                this.nodes[i][j].value = c;
                this.nodes[i][j].tracebackI = i - 1;
                this.nodes[i][j].tracebackJ = j - 1;
            }
            }
        }
        this.score = this.nodes[this.nodes.length - 1][
            this.nodes[0].length - 1
        ].value;
    }

    align() {
        this.alignedM = new Array();
        this.alignedN = new Array();
      
        var currentI = this.nodes.length - 1;
        var currentJ = this.nodes[0].length - 1;
      
        var currentNode = this.nodes[this.nodes.length - 1][this.nodes[0].length - 1];
      
        while (
          currentNode.tracebackI != undefined &&
          currentNode.tracebackJ != undefined
        ) {
          if (
            currentNode.tracebackI == currentI - 1 &&
            currentNode.tracebackJ == currentJ - 1
          ) {
            this.alignedM.push(this.M.pop());
            this.alignedN.push(this.N.pop());
          } else if (currentNode.tracebackJ == currentJ - 1) {
            this.alignedM.push("-");
            this.alignedN.push(this.N.pop());
          } else {
            this.alignedM.push(this.M.pop());
            this.alignedN.push("-");
          }
      
          currentI = currentNode.tracebackI;
          currentJ = currentNode.tracebackJ;
      
          currentNode = this.nodes[currentNode.tracebackI][currentNode.tracebackJ];
        }
      
        this.alignedM = this.alignedM.reverse();
        this.alignedN = this.alignedN.reverse();
    }

    getAlignedM() {
        return this.alignedM.join("");
    }

    getAlignedN() {
        return this.alignedN.join("");
    }
}

class AlignPairLinear {
    constructor() {
        this.M;
        this.N;
        this.alignedM;
        this.alignedN;
        this.scoreSet;
        this.Sn;
        this.Sp;
        this.score;
    }

    align() {
        if (this.M.length == 0) {
            for (var j = 1; j <= this.N.length; j++) {
                this.alignedM.push("-");
                this.alignedN.push(this.N[j - 1]);
                this.score = this.score + this.scoreSet.gap;
            }
        } else if (this.N.length == 0) {
            for (var j = 1; j <= this.M.length; j++) {
                this.alignedN.push("-");
                this.alignedM.push(this.M[j - 1]);
                this.score = this.score + this.scoreSet.gap;
            }
        } else if (this.M.length == 0 && this.N.length == 0) {
            // Do nothing
        } else {
            this.path(0, 0, this.M.length, this.N.length);
        }
    }

    path(i1, j1, i2, j2) {
        //alert ("i1, j1, : i2, j2 " + i1 +", " + j1 + ", " + i2 + ", " + j2);
      
        if (i1 + 1 == i2 || j1 == j2) {
          //align using quadratic space alignment
          var subM = new Array();
          var subN = new Array();
      
          for (var i = i1 + 1; i <= i2; i++) {
            subM.push(this.M[i - 1]);
          }
      
          for (var j = j1 + 1; j <= j2; j++) {
            subN.push(this.N[j - 1]);
          }
      
          var alignment = new AlignPairQuad();
      
          subScoreSet = new ScoreSet();
          if (j1 == j2) {
            if (j1 == 0) {
              subScoreSet.setScoreSetParam(
                this.scoreSet.scoringMatrix,
                this.scoreSet.beginGap,
                this.scoreSet.beginGap,
                this.scoreSet.beginGap
              );
            } else if (j1 == this.N.length) {
              subScoreSet.setScoreSetParam(
                this.scoreSet.scoringMatrix,
                this.scoreSet.endGap,
                this.scoreSet.endGap,
                this.scoreSet.endGap
              );
            } else {
              subScoreSet.setScoreSetParam(
                this.scoreSet.scoringMatrix,
                this.scoreSet.gap,
                this.scoreSet.gap,
                this.scoreSet.gap
              );
            }
          } else {
            subScoreSet.setScoreSetParam(
              this.scoreSet.scoringMatrix,
              this.scoreSet.gap,
              this.scoreSet.beginGap,
              this.scoreSet.endGap
            );
            subScoreSet.useBeginGapTop = false;
            subScoreSet.useBeginGapLeft = false;
            subScoreSet.useEndGapBottom = false;
            subScoreSet.useEndGapRight = false;
      
            if (i1 == 0) {
              subScoreSet.useBeginGapTop = true;
            }
      
            if (j1 == 0) {
              subScoreSet.useBeginGapLeft = true;
            }
      
            if (j2 == this.N.length) {
              subScoreSet.useEndGapRight = true;
            }
      
            if (i2 == this.M.length) {
              subScoreSet.useEndGapBottom = true;
            }
          }
      
          alignment.initializeMatrix(subM, subN, subScoreSet);
          alignment.fillMatrix();
          alignment.align();
          //alignment.dumpMatrix();
          this.alignedM.push(alignment.getAlignedM());
          this.alignedN.push(alignment.getAlignedN());
      
          this.score = this.score + alignment.score;
        } else {
          var middle = Math.floor((i1 + i2) / 2);
      
          //linear-space computation of alignment score to middle row
          //forward pass
      
          //gaps along top
      
          this.Sn[j1] = 0;
      
          if (i1 == 0) {
            for (var j = j1 + 1; j <= j2; j++) {
              this.Sn[j] = this.Sn[j - 1] - this.scoreSet.beginGap;
            }
          } else {
            for (var j = j1 + 1; j <= j2; j++) {
              this.Sn[j] = this.Sn[j - 1] - this.scoreSet.gap;
            }
          }
      
          //now continue down rows to middle row
          var diag;
          var left;
          //for (var i = i1 + 1; i <= i2; i++) {
          for (var i = i1 + 1; i <= middle; i++) {
            diag = this.Sn[j1];
            left;
            if (j1 == 0) {
              left = this.Sn[j1] - this.scoreSet.beginGap;
            } else {
              left = this.Sn[j1] - this.scoreSet.gap;
            }
      
            this.Sn[j1] = left;
      
            // we need three values to set the score: diag, left, and above to fill in the row
            for (var j = j1 + 1; j <= j2; j++) {
              // above will be in the this.Sn array, which is holding a mixture of the previous row and the new row
              // var above = this.Sn[j];
      
              //pick max of three and store in next left
              if (j == this.N.length && i == this.M.length) {
                left = Math.max(
                  this.Sn[j] - this.scoreSet.endGap,
                  Math.max(
                    left - this.scoreSet.endGap,
                    diag + this.scoreSet.getScore(this.M[i - 1], this.N[j - 1])
                  )
                );
              } else if (i == this.M.length) {
                left = Math.max(
                  this.Sn[j] - this.scoreSet.gap,
                  Math.max(
                    left - this.scoreSet.endGap,
                    diag + this.scoreSet.getScore(this.M[i - 1], this.N[j - 1])
                  )
                );
              } else if (j == this.N.length) {
                left = Math.max(
                  this.Sn[j] - this.scoreSet.endGap,
                  Math.max(
                    left - this.scoreSet.gap,
                    diag + this.scoreSet.getScore(this.M[i - 1], this.N[j - 1])
                  )
                );
              } else {
                left = Math.max(
                  this.Sn[j] - this.scoreSet.gap,
                  Math.max(
                    left - this.scoreSet.gap,
                    diag + this.scoreSet.getScore(this.M[i - 1], this.N[j - 1])
                  )
                );
              }
              diag = this.Sn[j];
      
              //prepares this.Sn for use in next iteration of i loop
              this.Sn[j] = left;
            }
          }
      
          //linear-space computation of alignment score to middle row
          //reverse pass
      
          //gaps along bottom
      
          this.Sp[j2] = 0;
      
          if (i2 == this.M.length) {
            for (var j = j2 - 1; j >= j1; j--) {
              this.Sp[j] = this.Sp[j + 1] - this.scoreSet.endGap;
            }
          } else {
            for (var j = j2 - 1; j >= j1; j--) {
              this.Sp[j] = this.Sp[j + 1] - this.scoreSet.gap;
            }
          }
      
          //now continue up rows to middle row
          var right;
          //for (var i = i2 - 1; i >= i1; i--) {
          for (var i = i2 - 1; i >= middle; i--) {
            diag = this.Sp[j2];
            if (j2 == this.N.length) {
              right = this.Sp[j2] - this.scoreSet.endGap;
            } else {
              right = this.Sp[j2] - this.scoreSet.gap;
            }
      
            this.Sp[j2] = right;
      
            //we need three values to set the score: diag, right, and below to fill in the row
            for (var j = j2 - 1; j >= j1; j--) {
              //below will be in the this.Sp array, which is holding a mixture of the previous row and the new row
              //var below = this.Sp[j];
      
              //pick max of three and store in next right
              if (j == 0 && i == 0) {
                right = Math.max(
                  this.Sp[j] - this.scoreSet.beginGap,
                  Math.max(
                    right - this.scoreSet.beginGap,
                    diag +
                      this.scoreSet.getScore(this.M[i + 1 - 1], this.N[j + 1 - 1])
                  )
                );
              } else if (j == 0) {
                right = Math.max(
                  this.Sp[j] - this.scoreSet.beginGap,
                  Math.max(
                    right - this.scoreSet.gap,
                    diag +
                      this.scoreSet.getScore(this.M[i + 1 - 1], this.N[j + 1 - 1])
                  )
                );
              } else if (i == 0) {
                right = Math.max(
                  this.Sp[j] - this.scoreSet.gap,
                  Math.max(
                    right - this.scoreSet.beginGap,
                    diag +
                      this.scoreSet.getScore(this.M[i + 1 - 1], this.N[j + 1 - 1])
                  )
                );
              } else {
                right = Math.max(
                  this.Sp[j] - this.scoreSet.gap,
                  Math.max(
                    right - this.scoreSet.gap,
                    diag +
                      this.scoreSet.getScore(this.M[i + 1 - 1], this.N[j + 1 - 1])
                  )
                );
              }
              diag = this.Sp[j];
              this.Sp[j] = right;
            }
          }
      
          // now find the value of j that maximizes this.Sn[j] + this.Sp[j]
          // this point will be in the maximum scoring path in the final alignment.
          // once we have this point we can divide the problem into two new problems,
      
          var maxValue = this.Sn[j1] + this.Sp[j1];
          var maxJ = j1;
      
          for (var j = j1 + 1; j <= j2; j++) {
            if (this.Sn[j] + this.Sp[j] >= maxValue) {
              maxValue = this.Sn[j] + this.Sp[j];
              maxJ = j;
            }
          }
      
          this.path(i1, j1, middle, maxJ);
          this.path(middle, maxJ, i2, j2);
        }
    }

    getAlignedM() {
        return this.alignedM.join("");
    }

    getAlignedN() {
        return this.alignedN.join("");
    }

    setAlignParam(M, N, scoreSet) {
        this.M = M;
        this.N = N;
        this.alignedM = [];
        this.alignedN = [];
        this.scoreSet = scoreSet;
        this.Sn = new Array(this.N.length);
        this.Sp = new Array(this.N.length);
        this.score = 0;
    }
}

// Function to read FASTA file and extract sequence
function getSequenceFromFasta(filePath) {
    const fastaContent = fs.readFileSync(filePath, 'utf8');
    const sequences = fastaContent.split('>');
    sequences.shift(); // Remove empty first element
    const sequenceMap = new Map();

    sequences.forEach(sequence => {
        const [title, ...lines] = sequence.split('\n').filter(line => line.trim() !== '');
        const sequenceString = lines.join('');
        sequenceMap.set(title, sequenceString);
    });

    return sequenceMap;
}

function pairwiseDna(titleOne, newDnaOne, titleTwo, newDnaTwo, matchScore, mismatchScore, gapPenalty, beginGapPenalty, endGapPenalty) {
    const useLinearSpace = true;
    const useQuadraticSpace = false;

    const matrix = new Identity();
    matrix.setMatch(matchScore);
    matrix.setMismatch(mismatchScore);

    const scoreSet = new ScoreSet();
    scoreSet.setScoreSetParam(matrix, gapPenalty, beginGapPenalty, endGapPenalty);

    let alignment;

    if (useLinearSpace) {
        alignment = new AlignPairLinear();
        alignment.setAlignParam(newDnaOne, newDnaTwo, scoreSet);
        alignment.align();

        console.log(">" + titleOne);
        console.log(alignment.getAlignedM());
        console.log();
        console.log(">" + titleTwo);
        console.log(alignment.getAlignedN());
        console.log();
        console.log("Alignment score:", alignment.score);
        console.log();
    }

    if (useQuadraticSpace) {
        alignment = new AlignPairQuad();
        alignment.initializeMatrix(newDnaOne, newDnaTwo, scoreSet);
        alignment.fillMatrix();
        //alignment.dumpMatrix();
        alignment.align();
    
        outputWindow.document.write(">" + titleOne + "\n");
        outputWindow.document.write(addReturns(alignment.getAlignedM()));
        outputWindow.document.write("\n");
        outputWindow.document.write("\n");
        outputWindow.document.write(">" + titleTwo + "\n");
        outputWindow.document.write(addReturns(alignment.getAlignedN()));
        outputWindow.document.write("\n\n");
        outputWindow.document.write("Alignment score: " + alignment.score + "\n\n");
    }
}

// Sample usage
if (process.argv.length < 4) {
    console.error('Insufficient arguments. Usage: node script.js fastaFile1 fastaFile2');
    process.exit(1);
}

const [, , fastaFile1, fastaFile2] = process.argv;
const fastaMap1 = getSequenceFromFasta(fastaFile1);
const fastaMap2 = getSequenceFromFasta(fastaFile2);

const titleOne = [...fastaMap1.keys()][0]; // Assuming there's only one sequence per file
const titleTwo = [...fastaMap2.keys()][0];
const newDnaOne = fastaMap1.get(titleOne);
const newDnaTwo = fastaMap2.get(titleTwo);

pairwiseDna(titleOne, newDnaOne, titleTwo, newDnaTwo, 2, -1, 2, 0, 0); // Default scoring parameters
