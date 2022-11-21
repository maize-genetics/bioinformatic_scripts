import java.io.*
import kotlin.system.exitProcess

/**
 *  To run this script at the command line, type:
 *   kotlinc -script scinetAnchorwaveScript.kts -- -refFasta <refFasta> -refGFF <refGFF> .... (all other parameters)
 *
 *  Input:
 *  1.  name with path for reference fasta
 *  2.  name with path for reference gff file
 *  3.  reference prefix which will be used to name output files, e.g. B73
 *  4.  A tab delimited file with column headers AssemblyFullName and AssemblyShortName. This file contains the list
 *      of assemblies to be aligned
 *  5.  Local outputDir - folder where the created scripts will be written.
 *  6.  Scinet outputDir - folder on scinet where output is written.  This is used in the created scripts.
 *
 * This script will take input and create from it 2 scripts that can be used to run anchorwave on scinet:
 * These scripts mimic the  steps run in the PHG plugin AssemblyMAFFromAnchorWavePlugin
 *  1.   A script that runs the initial steps to create a CDS fasta from a ref and gff, then
 *      align the resulting CDS fasta to the reference fasta as per documented anchorwave steps.
 *      This script creates files that will be used by all the assemblies that are to to be aligned in the second script.
 *  2.  A script that calls minimap2 and anchorwave proali to align each assembly against the reference
 *      SAM and CDS fasta created in step 1.
 *
 *  These scripts may be used on scinet with a conda anchorwave environment setup.
 *  The first script is run as a regular bash script on a single compute node.  Its output is used in the second script.
 *  The second script should make use of SLURM arrays, to align each assembly on its own compute node.  This script is
 *     run after the first script has successfully completed.
 *
 *  To use these on scinet, you'll need a slurm job that loads miniconda and an anchorwave 
 *  environment.  The latter already exists in the buckler_lab_project space.  Include these
 *  lines in your slurm file:
 *
 *   module load miniconda
 *   source activate /project/buckler_lab_panand/anchorwave
 *
 *  There should be 1 slurm file that runs the refAlignScript.sh NOT in an array.
 *  There should be a second slurm file that runs the assemblyAlignScript using slurm's array option.
 *
 * NOTE: This script does not have "main" in its name, but runs fine.  Perhaps
 * that is needed when there are dependencies to be loaded
 */

// Check for parameters 
if (!(args.contains("-refFasta")) || !(args.contains("-refGFF")) || !(args.contains("-refPrefix")) ||
	!(args.contains("-assemblyListFile")) || !(args.contains("-localOutputDir")) || !(args.contains("-scinetOutputDir")) ) {
    println("This script requires the following arguments: ")
    println("  -refFasta: full path to the reference fasta for aligning")
    println("  -refGFF:  full path to the reference gff file used by anchorwave")
    println("  -refPrefix: short name for the reference, which will be used as a prefix to create the reference sam file, e.g. B73")
    println("  -assemblyListFile: a tab-delimited file with column headers AssemblyFullName and AssemblyShortName.\n")
    println("          The first column contains full path names of the assembly fasta.")
    println("          the second column containing a shortname for the assembly, e.g. CML145.  This becomes the name asscoiated with the .sam file output")
    println("  -localOutputDir: Full path to the location where the anchorwave scripts should be written. ")
    println("  -scinetOutputDir: Full path to the folder on scinet where output files, e.g. the .sam files, will be written")
   println("Please re-run the script with the correct parameters. ")
   exitProcess(1)
}

val refFasta = args[1 + args.indexOf("-refFasta")]
val refGff = args[1 + args.indexOf("-refGFF")]
val refPrefix = args[1 + args.indexOf("-refPrefix")]
val assemblyListFile = args[1 + args.indexOf("-assemblyListFile")]
val localOutputDir = args[1 + args.indexOf("-localOutputDir")]
val scinetOutputDir = args[1 + args.indexOf("-scinetOutputDir")]
var numThreads = 1
if (args.contains("-numThreads")) numThreads = args[1 + args.indexOf("-numThreads")].toInt()

println("Begin processing first script ...")

val refAlignScript = "${localOutputDir}/refAlignScript.sh"

val cdsFasta = "${scinetOutputDir}/cdsFasta.fa"
val refSamOutput = "${scinetOutputDir}/${refPrefix}.sam"
File(refAlignScript).bufferedWriter().use { writer ->

    // First the ref cds fasta is created
    writer.write("anchorwave  gff2seq -r ${refFasta} -i ${refGff} -o ${cdsFasta}\n")

    // second, the ref is aligned via minimap2 against the cds fasta created above
    writer.write("minimap2 -x splice -t ${numThreads} -k 12 -a -p 0.4 -N20 ${refFasta} ${cdsFasta} -o ${refSamOutput}\n")

}
println("Finished creating refAlignScript - output written to ${refAlignScript}")

val assemblyAlignScript = "${localOutputDir}/assemblyAlignScript.sh"
val assemblyLines =  File(assemblyListFile).bufferedReader().readLines()

File(assemblyAlignScript).bufferedWriter().use { writer ->

    // For each assembly in the key file, align assembly and cdsFasta via minimap
    // Then align the assembly .sam and ref.sam via anchorwave
    for (line in assemblyLines) {
        if (line.startsWith("AssemblyFullName")) continue
        val assemblyData = line.split("\t")
        if (assemblyData.size != 2) {
            println("ERROR - assembly line must have 2 tab-delimited columns, one for full assembly name, one for assembly prefix. Error for line ${line}")
            continue
        }
        val assemblySam = "${scinetOutputDir}/${assemblyData[1]}.sam"
        val mafOutputFile = "${scinetOutputDir}/${assemblyData[1]}.maf"
        val anchorsproFile = "${scinetOutputDir}/${assemblyData[1]}_${refPrefix}.anchorspro"

        // This assumes the assemblyFullName is relative to scinet, so should be /project/buck... or /90daydata/...
        // Each line of the output file will contain 2 commands: one for minimap and 1 for anchorwave proali.  Both
        // commands must be run for each assembly.  Commands are separted by "&&"
        writer.write("minimap2 -x splice -t ${numThreads} -k 12 -a -p 0.4 -N20 ${assemblyData[0]} ${cdsFasta} -o ${assemblySam} && " +
                    "anchorwave proali -i ${refGff} -r ${refFasta} -as ${cdsFasta} -a ${assemblySam} " +
                    "-ar ${refSamOutput} -s ${assemblyData[0]} -n ${anchorsproFile} -R 1 -Q 1 -t ${numThreads} -o ${mafOutputFile}\n")
    }
    println("Finished creating assemblyAlignScript - output written to ${assemblyAlignScript}")
}

println("\nDONE - scripts written to folder ${localOutputDir}")


