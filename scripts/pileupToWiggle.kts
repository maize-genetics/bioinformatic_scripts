import java.io.File

/**
 * This script takes output from samtools mpileup and converts it to wiggle format
 * so that coverage can be viewed in programs such as igvtools.
 * The script takes one or two parameters.
 * The first (required) parameter is the input mpileup file.
 * The second, optional parameter is the name of the output file.
 * If an output file is not provided, the script will write to standard output.
 *
 * @author ahb232
 *
 */

val helpMessage = "script converts samtools mpileup to wiggle format\n" +
        "Usage: kotlinc -script pileupToWiggle.kts -- [INPUT_FILE] [OUTPUT_FILE]\n" +
        "INPUT_FILE: input file in samtools mpileup format\n" +
        "OUTPUT_FILE: (optional) output file. If not provided, will write to standard output\n"

var fileName = ""
var outFile = ""

if (args.size == 0) {
    System.err.println("ERROR: Input file required.")
    print(helpMessage)
    System.exit(0)
} else {
    if (args[0] == "--help" || args[0] == "-h") {
        print(helpMessage)
        System.exit(0)
    }
    fileName = args[0]

    if (args.size == 2) {
        outFile = args[1]
    } else if (args.size > 2) {
        System.err.println("ERROR: Too many arguments")
        print(helpMessage)
        System.exit(0)
    }

}

var currentChrom = ""
var currentPos = 0

val reader = File(fileName).bufferedReader()
val writer = if (outFile != "") { File(outFile).bufferedWriter() } else { null }

var line = reader.readLine()

while(line != null) {

    val fields = line.split("\t")

    if(fields[0] != currentChrom || fields[1].toInt() != currentPos + 1) {
        if (writer != null) {
            writer.write("fixedStep chrom=${fields[0]} start=${fields[1]} step=1\n")
        } else {
            print("fixedStep chrom=${fields[0]} start=${fields[1]} step=1\n")
        }
        currentChrom = fields[0]
    }
    currentPos=fields[1].toInt()
    if (writer != null) {
        writer.write("${fields[3]}\n")
    } else {
        print("${fields[3]}\n")
    }
    line = reader.readLine()
}

reader.close()
if (writer != null){ writer.close() }

