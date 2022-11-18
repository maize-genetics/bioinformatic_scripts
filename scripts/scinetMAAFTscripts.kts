import java.io.File
import kotlin.system.exitProcess

/**
 * This script created for Aimee Schulz to run batches of MAFFT command via slurm on scinet.  There will
 * be multiple MAFFT commands per line, based on the number provided in the numBatch parameter.  The output is a
 * single file that may be run via a slurm job using arrays.
 *
 * THis script using ";" to separate the commands on a single line.  The ";" will execute subsequent commands when
 * previous ones file.  If this is replaced by "&&", then a failure will stop execution - no subsequent commands will
 * be run.  Change this script accordingly if you care.
 *
 * To run this script on the command line:
 *    kotlinc -script scinetMAAFTscripts.kts -- -fastaListFile <input file with fasta> -scinetInputDir <Folder holding fasta on scinet> .... (all other parameters)
 *
 * It takes as input the following:
 *   1. a file that contains the names of fasta files one per line. Just filename, no path
 *   2. the path on scinet where the inputfiles will live
 *   3. the path on scinet where the output files should be written
 *   4. FUll path and name of the scinet script to be written
 *   5. number of files to create in a batch.
 */

// Check for parameters
if (!(args.contains("-fastaListFile")) || !(args.contains("-scinetInputDir")) || !(args.contains("-scinetOutputDir")) ||
    !(args.contains("-scinetScriptFile")) || !(args.contains("-batchNum"))) {
    println("\nThis script requires the following arguments: ")
    println("  -fastaListFile: full path to the local file that contains the list of fastas to be run through mafft, one file per line, just name, no path")
    println("  -scinetInputDir:  the path on scinet where the input files will live")
    println("  -scinetOutputDir: The path on scinet where the output files will be written.")
    println("  -batchNum: the number of mafft commands to be run per node.  This is the number of files that will be written to each script file that become entries in the slurm input file.")
    println("  -scinetScriptFile: Filename with full path where the scinet script should be written. ")
    println("\nPlease re-run the script with the correct parameters. ")
    exitProcess(1)
}

val fastaFileFile = args[1 + args.indexOf("-fastaListFile")]
val sInDir = args[1 + args.indexOf("-scinetInputDir")]
val sOutDir = args[1 + args.indexOf("-scinetOutputDir")]
val batchNum = args[1 + args.indexOf("-batchNum")].toInt()
val scriptFile = args[1 + args.indexOf("-scinetScriptFile")]

println("Begin processing ...")

// read the list of fasta files to process
val listOfFastas =  File(fastaFileFile).bufferedReader().readLines()

println("Begin... batchNum=${batchNum}, fastaFile size=${listOfFastas.size}\n")
File(scriptFile).bufferedWriter().use { writer ->

    val beginLine = "mafft --ep 0 --genofpair --adjustdirection ${sInDir}/"

    var fastaLine = 0
    for (idx in 0 until (listOfFastas.size-batchNum) step batchNum) {
        val lineSB = StringBuilder()

        // The check in the for loop above guarantees this loop won't attempt to access entries
        // in the fasta file beyond those that exist.
        for (jdx in 0 until batchNum) {
            val fastaIn = listOfFastas.get(idx+jdx)
            val lastDot = fastaIn.lastIndexOf(".")
            var fastaOut = fastaIn.substring(0,lastDot)
            lineSB.append("${beginLine}${listOfFastas.get(idx+jdx)} > ${sOutDir}/${fastaOut}_aligned.fa;")
        }
        // If this must stop when one in the batch fails, then replace the ";" with " && "
        // ";" will keep executing the next command on the line. "&&" will stop executing on error
        //lineSB.setLength(lineSB.length-3) // only if using &&
        lineSB.append("\n")
        println("writing: ${lineSB.toString()}")
        writer.write(lineSB.toString())
        fastaLine = idx
    } // end for idx in 0-fastaFiles.size

    // process remaining lines
    var currentIdx = fastaLine + batchNum
    println("\nEnd of for loop - fastaLine=${fastaLine}, fastaLine+batchNum = ${currentIdx}")
    if (currentIdx < listOfFastas.size) {
        val lineSB = StringBuilder()
        while (currentIdx < listOfFastas.size) {
            val fastaIn = listOfFastas.get(currentIdx)
            val lastDot = fastaIn.lastIndexOf(".")
            var fastaOut = fastaIn.substring(0,lastDot)
            lineSB.append("${beginLine}${listOfFastas.get(currentIdx)} > ${sOutDir}/${fastaOut}_aligned.fa;")
            currentIdx++
        }
        println("\nWriting last fastas to file ${scriptFile}")
        writer.write(lineSB.toString())

    }

} // end writer

println("Finished!! output file is in ${scriptFile}")
