//
// This file holds several functions specific to the preprocessing workflow in the nf-cmgg/preprocessing pipeline
//

class GenomeUtils {
    //
    // Get attribute from genome config file e.g. fasta
    //
    def getGenomeAttribute(genome, attribute) {
        if (genome.containsKey(attribute)) {
            return nextflow.Nextflow.file(genome[attribute], checkIfExists: true)
        } else {
            nextflow.Nextflow.error("Genome config does not contain attribute ${attribute}")
        }
        return null
    }
}
