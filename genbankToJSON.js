const fs = require('fs');

const MONTHS = [
    'JAN', 'FEB', 'MAR',
    'APR', 'MAY', 'JUN',
    'JUL', 'AUG', 'SEP',
    'OCT', 'NOV', 'DEC',
];

// GenBank specification: https://www.ncbi.nlm.nih.gov/genbank/samplerecord/
// In order: locus name, sequence length, molecule type (e.g. DNA), genbank division, modification date
// GenBank divisions: https://www.ncbi.nlm.nih.gov/genbank/htgs/divisions/
// Locus definition has changed with time, use accession number for a unique identifier

const GENBANK_ANNOTATION_KEY = {
    LOCUS_TAG: 'LOCUS',
    DEFINITION_TAG: 'DEFINITION',
    ACCESSION_TAG: 'ACCESSION',
    VERSION_TAG: 'VERSION',
    KEYWORDS_TAG: 'KEYWORDS',
    SOURCE_TAG: 'SOURCE',
    ORGANISM_TAG: 'ORGANISM',
    REFERENCE_TAG: 'REFERENCE',
    AUTHORS_TAG: 'AUTHORS',
    CONSORTIUM_TAG: 'CONSRTM',
    TITLE_TAG: 'TITLE',
    JOURNAL_TAG: 'JOURNAL',
    PUBMED_TAG: 'PUBMED',
    REMARK_TAG: 'REMARK',
    FEATURES_TAG: 'FEATURES',
    BASE_COUNT_TAG: 'BASE COUNT',
    ORIGIN_TAG: 'ORIGIN',
    CONTIG_TAG: 'CONTIG',
    END_SEQUENCE_TAG: '//',
};

function genbankToJson(sequence) {
    if (typeof sequence !== 'string') {
        throw new TypeError('Input must be a string!');
    }

    const resultsArray = [];
    let result;
    let currentFeatureNote;

    let lines = sequence.split(/\r?\n/);
    let fieldName;
    let subFieldType;
    let featureLocationIndentation;
    let lastLineWasFeaturesTag;
    let lastLineWasLocation;

    let hasFoundLocus = false;

    for (let line of lines) {
        if (line === null) break;
        let lineFieldName = getLineFieldName(line);
        let val = getLineVal(line);
        let isSubKey = isSubKeyword(line);
        let isKey = isKeyword(line);

        if (lineFieldName === GENBANK_ANNOTATION_KEY.END_SEQUENCE_TAG || isKey) {
            fieldName = lineFieldName;
            subFieldType = null;
        } else if (isSubKey) {
            subFieldType = lineFieldName;
        }
        // Ignore these lines
        if (line.trim() === '' || lineFieldName === ';') {
            continue;
        }

        if (!hasFoundLocus && fieldName !== GENBANK_ANNOTATION_KEY.LOCUS_TAG) {
            throw new Error('No LOCUS tag was found!');
        }

        switch (fieldName) {
            case GENBANK_ANNOTATION_KEY.LOCUS_TAG:
                hasFoundLocus = true;
                parseLocus(line);
                break;
            case GENBANK_ANNOTATION_KEY.FEATURES_TAG:
                parseFeatures(line, lineFieldName, val);
                break;
            case GENBANK_ANNOTATION_KEY.ORIGIN_TAG:
                parseOrigin(line, lineFieldName);
                break;
            case GENBANK_ANNOTATION_KEY.CONTIG_TAG:
                parseContig(line, lineFieldName);
                break;
            case GENBANK_ANNOTATION_KEY.DEFINITION_TAG:
            case GENBANK_ANNOTATION_KEY.ACCESSION_TAG:
            case GENBANK_ANNOTATION_KEY.VERSION_TAG:
            case GENBANK_ANNOTATION_KEY.KEYWORDS_TAG:
                parseMultiLineField(fieldName, line, fieldName.toLowerCase());
                break;
            case GENBANK_ANNOTATION_KEY.SOURCE_TAG:
                if (subFieldType === GENBANK_ANNOTATION_KEY.ORGANISM_TAG) {
                    parseMultiLineField(subFieldType, line, 'organism');
                } else {
                    parseMultiLineField(lineFieldName, line, 'source');
                }
                break;
            case GENBANK_ANNOTATION_KEY.REFERENCE_TAG:
                if (lineFieldName === GENBANK_ANNOTATION_KEY.REFERENCE_TAG) {
                    const ref = result.references || [];
                    result.references = ref;
                    ref.push({});
                }
                parseReference(line, subFieldType);
                break;
            case GENBANK_ANNOTATION_KEY.END_SEQUENCE_TAG:
                endSeq();
                break;
            default:
                // Unhandled tag
                break;
        }
    }

    // Catch the case where successfully started a sequence and parsed it, but endSeq isn't called correctly
    if (resultsArray[resultsArray.length - 1] !== result) {
        endSeq();
    }
    return resultsArray;

    function endSeq() {
        // Do some post processing clean-up
        postProcessCurSeq();
        resultsArray.push(result);
    }

    function getCurrentFeature() {
        return result.features[result.features.length - 1];
    }

    function postProcessCurSeq() {
        if (result && result.features) {
            for (let i = 0; i < result.features.length; i++) {
                result.features[i] = postProcessGenbankFeature(result.features[i]);
            }
        }
    }

    function parseOrigin(line, key) {
        if (key !== GENBANK_ANNOTATION_KEY.ORIGIN_TAG) {
            let newLine = line.replace(/[\s]*[0-9]*/g, '');
            result.sequence += newLine;
        }
    }

    function parseContig(line, key) {
        if (key !== GENBANK_ANNOTATION_KEY.CONTIG_TAG) {
            let newLine = line.replace(/[\s]*[0-9]*/g, '');
            result.sequence += newLine;
        }
    }

    function parseLocus(line) {
        // Initialize the result object with default values
        result = {
            features: [],
            name: 'Untitled sequence',
            sequence: '',
            references: [],
        };

        // Remove the LOCUS tag and trim the line
        line = removeFieldName(GENBANK_ANNOTATION_KEY.LOCUS_TAG, line).trim();
        // Extract information using regular expressions
        const m = line.match(/^(\S+)\s+(\d+)\s+bp\s+(\S+)\s+(\S+)\s+(\S+)\s*(\S+)?$/);
        if (!m) {
            throw new Error('Invalid LOCUS line format');
        }

        // Extract individual fields
        const locusName = m[1];
        const size = +m[2];
        const moleculeType = m[3];
        const circular = m[4] === 'circular';
        const seq = result;
        let dateStr = '';

        // Determine the date string
        if (!m[6]) {
            dateStr = m[5];
        } else {
            seq.genbankDivision = m[5];
            dateStr = m[6];
        }

        // Parse the date string
        const dateMatch = dateStr.match(/^(\d{2})-(.{3})-(\d{4})$/);
        if (!dateMatch) {
            throw new Error('Invalid date format');
        }
        const date = new Date(+dateMatch[3], MONTHS.indexOf(dateMatch[2].toUpperCase()), +dateMatch[1], 12, 0, 0, 0);

        // Populate the result object with parsed values
        seq.circular = circular;
        seq.moleculeType = moleculeType;
        seq.date = date.toISOString();
        seq.name = locusName;
        seq.size = size;

        return result;
    }

    function removeFieldName(fName, line) {
        line = line.replace(/^\s*/, '');
        if (line.indexOf(fName) === 0) {
            line = line.replace(fName, '');
        }
        return line.trim();
    }

    function parseReference(line, subType) {
        const refs = result.references;
        let lastRef = refs[refs.length - 1];
        if (!subType) {
            parseMultiLineField(
                GENBANK_ANNOTATION_KEY.REFERENCE_TAG,
                line,
                'description',
                lastRef,
            );
        } else {
            parseMultiLineField(subType, line, subType.toLowerCase(), lastRef);
        }
    }

    function parseFeatures(line, key, val) {
        let strand;
        // For the main features Location/Qualifiers line 
        if (key === GENBANK_ANNOTATION_KEY.FEATURES_TAG) {
            lastLineWasFeaturesTag = true;
            return;
        }

        if (lastLineWasFeaturesTag) {
            // Get the indentation of feature locations
            featureLocationIndentation = getLengthOfWhiteSpaceBeforeStartOfLetters(line, );
            // Reset lastLineWasFeaturesTag
            lastLineWasFeaturesTag = false;
        }
        // For Location/Qualifiers lines
        if (isFeatureLineRunon(line, featureLocationIndentation)) {
            // The line is a continuation of the above line
            if (lastLineWasLocation) {
                // The last line was a location, so the run-on line is expected to be a feature location too
                parseFeatureLocation(line.trim());
                lastLineWasLocation = true;
            } else {
                // The last line was a note
                if (currentFeatureNote) {
                    // Append to the currentFeatureNote
                    currentFeatureNote[
                        currentFeatureNote.length - 1
                    ] += line.match(/translation/) ? 
                    line.trim().replace(/"/g, '') : 
                    ` ${line.trim().replace(/"/g, '')}`;
                }
                lastLineWasLocation = false;
            }
        } else {
            // New Element/Qualifier lines. Not runon lines.
            if (isNote(line)) {
                // Is a new Feature Element (e.g. source, CDS) in the form of  "[\s] KEY  SEQLOCATION"
                // Is a Feature Qualifier in the /KEY="<VALUE>" format; could be multiple per Element
                // Check that feature did not get skipped for missing location
                if (getCurrentFeature()) {
                    parseFeatureNote(line);
                    lastLineWasLocation = false;
                }
            } else {
                // The line is a location, so make a new feature from it
                strand = val.match(/complement/g) ? -1 : 1; // Determine the strand

                newFeature();
                let feat = getCurrentFeature();
                feat.type = key;
                feat.strand = strand;

                parseFeatureLocation(val);
                lastLineWasLocation = true;
            }
        }
    }

    function newFeature() {
        result.features.push({
        notes: {},
        });
    }

    function isNote(line) {
        // Trim leading and trailing whitespace
        line = line.trim();
    
        // Check if the line starts with a slash (/)
        if (line.startsWith('/')) {
            return true;
        }
    
        // Check if the line matches the pattern "/key=value"
        if (line.match(/^\/[\w-]+=/)) {
            return true;
        }
    
        // If none of the above conditions are met, it's not a note
        return false;
    }

    function parseFeatureLocation(locStr) {
        // Trim any leading or trailing whitespace
        locStr = locStr.trim();
    
        // Use a regular expression to extract numbers from the location string
        const matches = locStr.match(/\d+/g);
    
        // If there are no matches, return early
        if (!matches || matches.length === 0) {
            return;
        }
    
        // Parse the start and end positions from the matches
        const start = +matches[0];
        const end = matches.length > 1 ? +matches[1] : start;
    
        // Get the current feature
        const feat = getCurrentFeature();
    
        // Update the feature's start and end positions
        feat.start = start;
        feat.end = end;
    }

    function parseFeatureNote(line) {
        // Trim leading and trailing whitespace and remove leading slashes and quotes
        line = line.trim().replace(/^\/|"$/g, '');
    
        // Split the line into key-value pairs
        const [key, rawValue] = line.split(/="|=/);
    
        // Process the value
        let value = rawValue ? rawValue.replace(/\\/g, ' ').replace(/".*/g, '') : '';
        value = isNaN(value) ? value : +value; // Convert to number if possible
    
        // Get the current feature's notes
        const currentNotes = getCurrentFeature().notes;
    
        // Update the notes object with the key-value pair
        if (currentNotes[key]) {
            currentNotes[key].push(value);
        } else {
            currentNotes[key] = [value];
        }
    
        // Set the current feature note
        currentFeatureNote = currentNotes[key];
    }

    function getLineFieldName(line) {
        let arr;
        arr = line.trim().split(/[\s]+/);

        return arr[0];
    }

    function parseMultiLineField(fName, line, resultKey, r) {
        r = r || result;
        let fieldValue = removeFieldName(fName, line);
        r[resultKey] = r[resultKey] ? `${r[resultKey]} ` : '';
        r[resultKey] += fieldValue;

        return r
    }

    function getLineVal(line) {
        line = line.trim(); // Trim leading and trailing whitespace
    
        // If the line doesn't contain '=', return the trimmed line
        if (line.indexOf('=') < 0) {
            return line;
        } else {
            // Split the line by '='
            const arr = line.split('=');
            // Join all parts after the first '=' as the value
            const value = arr.slice(1).join('=').trim();
            return value;
        }
    }

    function isKeyword(line) {
        let isKey = false;
        if (line.substr(0, 10).match(/^[\S]+/)) {
            isKey = true;
        }
        return isKey;
    }

    function isSubKeyword(line) {
        let isSubKey = false;
        if (line.substr(0, 10).match(/^[\s]+[\S]+/)) {
            isSubKey = true;
        }
        return isSubKey;
    }

    function postProcessGenbankFeature(feat) {
        if (feat.notes.label) {
            feat.name = feat.notes.label[0];
        } else if (feat.notes.gene) {
            feat.name = feat.notes.gene[0];
        } else if (feat.notes.ApEinfo_label) {
            feat.name = feat.notes.ApEinfo_label[0];
        } else if (feat.notes.name) {
            feat.name = feat.notes.name[0];
        } else if (feat.notes.organism) {
            feat.name = feat.notes.organism[0];
        } else if (feat.notes.locus_tag) {
            feat.name = feat.notes.locus_tag[0];
        } else if (feat.notes.note) {
            feat.name = feat.notes.note[0];
        } else {
            feat.name = 'Untitled Feature';
        }
        feat.name = typeof feat.name === 'string' ? feat.name : String(feat.name);
        return feat;
    }
    function isFeatureLineRunon(line, featureLocationIndentation) {
        const indentationOfLine = getLengthOfWhiteSpaceBeforeStartOfLetters(line);
            // Check if the indentation of the line matches the expected feature location indentation
        if (featureLocationIndentation === indentationOfLine) {
            // If the indentations match, it's not a run-on line
            /*
            The feature location indentation calculated right after the feature tag
            Cannot be the same as the indentation of the lineß
            FEATURES             Location/Qualifiers
                 rep_origin      complement(1074..3302)
            01234  <-- this is the indentation we're talking about
            */
            return false;
        }
        
        // Trim the line and check if it starts with a slash (/)
        const trimmedLine = line.trim();
        if (trimmedLine.startsWith('/')) {
            // If the trimmed line starts with a slash, it's not a run-on line
            return false;
        }

        // If none of the above conditions are met, it's a run-on line
        return true;
        // Run-on line example:
        // FEATURES             Location/Qualifiers
        //     rep_origin      complement(1074..3302)
        //                 /label=pSC101**
        //                 /note="REP_ORIGIN REP_ORIGIN pSC101* aka pMPP6, gives plasm
        //                 id number 3 -4 copies per cell, BglII site in pSC101* ori h <--run-on line!
        //                 as been dele ted by quick change agatcT changed to agatcA g <--run-on line!
        //                 iving pSC101* * pSC101* aka pMPP6, gives plasmid number 3-4 <--run-on line!
        //                 copies p er cell, BglII site in pSC101* ori has been delet  <--run-on line!
        //                 ed by quic k change agatcT changed to agatcA giving pSC101* <--run-on line!
        //                 * [pBbS0a-RFP]"                                             <--run-on line!
        //                 /gene="SC101** Ori"
        //                 /note="pSC101* aka pMPP6, gives plasmid number 3-4 copies p
        //                 er cell, BglII site in pSC101* ori has been deleted by qui
        //                 c k change agatcT changed to agatcA giving pSC101**"
        //                 /vntifkey="33"
    }

    function getLengthOfWhiteSpaceBeforeStartOfLetters(string) {
        // Match the whitespace characters before the start of letters
        let match = string.match(/^\s*/);
        
        // Return the length of matched whitespace (or 0 if there's no match)
        return match ? match[0].length : 0;
    }
}

if (process.argv.length < 3) {
    console.error('Insufficient arguments. Usage: node script.js file.gbk');
    process.exit(1);
}

const [, , genbankFile] = process.argv;

fs.readFile(genbankFile, 'utf8', (err, data) => {
    if(err) {
        console.error('Error reading file:', err);
        return;
    }
    
    const result = genbankToJson(data);

    console.log(console.log(JSON.stringify(result)));
});
