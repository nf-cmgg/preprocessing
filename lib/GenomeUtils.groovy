//
// This file holds several functions specific to the preprocessing workflow in the nf-cmgg/preprocessing pipeline
//

import nextflow.Nextflow

class GenomeUtils {
    //
    // Get attribute from genome config file e.g. fasta
    //
    public static Object getGenomeAttribute(genome, attribute) {
        if (genome.containsKey(attribute)) {
            return Nextflow.file(genome[attribute], checkIfExists: true)
        } else {
            Nextflow.error("Genome config does not contain attribute ${attribute}")
        }
        return null
    }
}
