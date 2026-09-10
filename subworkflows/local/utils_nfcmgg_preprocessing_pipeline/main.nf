/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow UTILS_NFCMGG_PREPROCESSING_PIPELINE {

    main:
    dummy_emit = true

    emit:
    dummy_emit
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
    rg.PL = 'ILLUMINA'
    rg.SM = SM ?: fastq.name.toString() - ~/_R[0-9]_001.*$/
    if (LB) {
        rg.LB = LB
    }
    if (CN) {
        rg.CN = CN
    }
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

//
// Generate readgroup from bclconvert outputs
//
def getReadgroupsFromBclconvert(ch_fastq_list_csv, ch_fastq) {
    return ch_fastq_list_csv
        .join(ch_fastq, by: [0])
        .map { meta, csv_file, fastq_list ->
            def meta_fastq = []
            csv_file
                .splitCsv(header: true)
                .each { row ->
                    // Create the readgroup tuple
                    // RGID,RGSM,RGLB,Lane,Read1File,Read2File
                    def rg = [:]
                    // row.RGID is index1.index2.lane
                    rg.ID = row.RGID
                    // RGPU is a custom column in the samplesheet containing the flowcell ID
                    rg.PU = row.RGPU ? row.RGPU : meta.id + "." + row.Lane
                    rg.SM = row.RGSM
                    if (row.RGLB) {
                        rg.LB = row.RGLB
                    }
                    rg.PL = "ILLUMINA"

                    // dereference the fastq files in the csv
                    def fastq1 = fastq_list.find { fq -> file(fq).name == file(row.Read1File).name }
                    if (!fastq1) {
                        error("BCL Convert fastq_list.csv Read1File '${row.Read1File}' for sample '${row.RGSM}' was not in the demultiplexed FASTQs")
                    }
                    def fastq2 = row.Read2File ? fastq_list.find { fq -> file(fq).name == file(row.Read2File).name } : null

                    // set fastq metadata
                    def new_meta = meta + [id: fastq1.getSimpleName().toString() - ~/_R[0-9]_001.*$/, readgroup: rg, single_end: !fastq2]

                    meta_fastq << [new_meta, fastq2 ? [fastq1, fastq2] : [fastq1]]
                }
            return meta_fastq
        }
        .flatMap()
}

//
// Pick the sampleinfo row for a demultiplexed FASTQ.
// Sampleinfo is parsed once per flowcell lane, so identical rows repeat and are deduplicated.
// One row per samplename is attached as-is. Multiple rows require a unique
// match of sampleinfo.library to readgroup.LB.
//
def matchSampleinfo(meta, infos) {
    def rows = (infos instanceof Collection ? infos as List : [infos]).unique(false)
    if (rows.size() == 1) {
        return rows[0]
    }
    def samplename = meta.readgroup?.SM ?: meta.samplename
    def lb = meta.readgroup?.LB
    def matches = rows.findAll { row -> row.library && lb && row.library == lb }
    if (matches.size() != 1) {
        error("Multiplexed sample '${samplename}' needs a unique sampleinfo row for library '${lb}'. Set Illumina LibraryName so BCL Convert RGLB matches sampleinfo.library.")
    }
    return matches[0]
}

//
// Associate demultiplexed FASTQs [meta, fastq] with sampleinfo rows.
//
def associateSampleinfo(ch_fastq, ch_sampleinfo) {
    def ch_info = ch_sampleinfo
        .map { sampleinfo -> [sampleinfo.samplename, sampleinfo] }
        .groupTuple()
    return ch_fastq
        .map { meta, fastq -> [meta.readgroup.SM, meta, fastq] }
        .combine(ch_info, by: 0)
        .map { _samplename, meta, fastq, infos ->
            return [meta, fastq, matchSampleinfo(meta, infos)]
        }
}

//
// Output directory for a sample. Multiplexed libraries stay under library/samplename.
//
def samplePublishDir(meta) {
    def samplename = meta.samplename ?: meta.id
    return meta.library ? "${meta.library}/${samplename}" : "${samplename}"
}
