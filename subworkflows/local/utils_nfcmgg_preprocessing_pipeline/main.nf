/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Mock subworkflow to please linting
workflow UTILS_NFCMGG_PREPROCESSING_PIPELINE {
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Generate readgroup information from FASTQ header
//
def getReadgroupFromFastq(fastq, SM, LB, CN) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    def line
    fastq.withInputStream { fq ->
        def isGzip = fastq.name.toString().endsWith('.gz')
        def stream = isGzip ? new java.util.zip.GZIPInputStream(fq) as InputStream : fq as InputStream
        def decoder = new InputStreamReader(stream, 'ASCII')
        def buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    def rg = [:]
    rg.LB = LB ?: ''
    rg.CN = CN ?: ''
    rg.PL = 'ILLUMINA'
    rg.SM = SM ?: fastq.name.toString() - ~/_R[0-9]_001.*$/
    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        // def sequencer_serial = fields[0]
        // def run_number       = fields[1]
        def fcid = fields[2]
        def lane = fields[3]
        def index = fields[-1] ==~ /^[GATCN+-]+$/ ? fields[-1] : ''
        rg.ID = [index ?: fcid, lane].join('.')
        rg.PU = [fcid, lane].join('.')
    }
    else if (fields.size() == 5) {
        def fcid = fields[0]
        def lane = fields[1]
        rg.ID = [fcid, lane].join('.')
        rg.PU = [fcid, lane].join('.')
    }
    return rg
}
