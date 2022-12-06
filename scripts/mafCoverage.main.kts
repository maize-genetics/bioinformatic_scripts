@file:DependsOn("org.biokotlin:biokotlin:0.05.01")
@file:DependsOn("org.jetbrains.kotlinx:dataframe:0.8.0-rc-7")

import biokotlin.genome.getCoverageIdentityPercentForMAF
import java.util.*
import java.io.File
import org.jetbrains.kotlinx.dataframe.api.*
import org.jetbrains.kotlinx.dataframe.io.writeCSV

/**
 * Given a directory containing MAF files, determine the % coverage
 * of each chromosome against the reference, storing the results in a CSV
 */

val MAF_LOCATION = args[0]
var covidDf = emptyDataFrame()

File(MAF_LOCATION).walk().filter{it.extension == "maf"}.forEach{
    val mafFilename = it.toString()
    val cultivarName = mafFilename.substring(MAF_LOCATION.length + 1, mafFilename.indexOf(".maf"))

    val mafDf = getCoverageIdentityPercentForMAF(mafFilename)
    var chromosomes = dataFrameOf("cultivar")(cultivarName)

    for (dfRow in mafDf!!) {
        val contigName = dfRow["contig"] as String
        val percentCov = dfRow["percentCov"] as Double
        if (contigName !in chromosomes.columnNames()) {
            chromosomes = chromosomes.add(contigName) { percentCov }
        }
    }

    covidDf = covidDf.concat(chromosomes)
}

covidDf!!.writeCSV("maf-coverage.csv")
