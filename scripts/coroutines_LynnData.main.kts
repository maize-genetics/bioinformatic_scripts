#!/usr/bin/env kscript

/**
 * This is a test script whose purpose is to demonstrate running Kotlin co-routines
 * in a kotlin script.  All the data processed is dummy data created programmatically
 *
 *  Run this script via: (replace outputDir parameter with your own)
 *  kotlinc -script coroutines_LynnData.main.kts -- -outputDir /workdir/lcj34/kotlin_scripts/coroutinesTestOutput -numThreads 4
 *  NOTE - the "numThreads" parameter is optional, "outputDir" is required
 *
 * kotlinx-coroutines 1.3.2 needs both DependsOn and jcenter
 * kotlinx-coroutines 1.6.4 must only have the coroutines-core, NOT coroutines-core-common
 * kotlinc version 1.7.20 uses the 1.6.4 coroutines version
 *
 * tested on bl01 with kotlinc 1.7.20 installed
 *
 * @author lcj34
 *
 */

//@file:Repository("https://jcenter.bintray.com")
//@file:DependsOn("org.jetbrains.kotlinx:kotlinx-coroutines-core:1.3.2")
//@file:DependsOn("org.jetbrains.kotlinx:kotlinx-coroutines-core-common:1.3.2")
@file:DependsOn("org.jetbrains.kotlinx:kotlinx-coroutines-core:1.6.4")

import kotlin.system.exitProcess
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import java.util.*
import java.io.*

// Setup the output folder
if (!args.contains("-outputDir") ) {
    println("must specify input parameters with -outputDir parameter")
    exitProcess(1)
}
val outputDir = args[1 + args.indexOf("-outputDir")]

// Setup the number of threads for processing
var numThreads = 4;
if (args.contains("-numThreads")) {
   numThreads = args[1 + args.indexOf("-numThreads")].toInt()
}


/**
 * Data class for co-routine Input Channel
 */
data class InputChannelData(val transcript: String, val dataList: List<String>, val fileName: String)


fun processData() {
	println("Creating dataMap of fake data to write")
	val dataMap = createDataMap()

	println("Begin runBlocking ...")
	runBlocking {

	  val inputChannel = Channel<InputChannelData>()

	  launch {
		// We take the lines from the map created above.
		// We create an InputChannelData entry, put it on the
		// input channel.
		// The method that processes this channel will write
		// each value to a file.

		val transcriptEntries = dataMap.entries
		transcriptEntries.forEach { entry ->
		  val ts = entry.key
		  val data = entry.value
		  val outFile = "${outputDir}/${ts}_msa.fa"
		  inputChannel.send(InputChannelData(ts,data,outFile))

		}

		println("Done adding to the input channel")
	        inputChannel.close()

	  }
	
          val workerThreads = (1 .. numThreads).map { threadNum ->
		// there are no other parameters - the result is the file that is written
		// In a fancier example, user might add processed data to a "result" channel,
		// which might then be loaded to a database, or written out, or otherwise processed.
		launch {writeFiles(inputChannel) }
	  }

	  // We have already closed the input channel - nothing else to do here

	} // end runBlocking
	


}

// This creates a dummy object - test data
fun createDataMap(): Map<String, List<String>> {
    val dataMap = mutableMapOf<String,MutableList<String>>()
    val transcriptList = arrayOf("ts1","ts2","ts3","ts4","ts5","ts6","ts7","ts8","ts9","ts10")
    for (transcript in transcriptList) {
        val tsData = arrayOf("${transcript}data1","${transcript}data2","${transcript}data3","${transcript}data4","${transcript}data5","${transcript}data6","${transcript}data7","${transcript}data8")
        val tsDataList:MutableList<String> = tsData.asList() as MutableList<String>
        dataMap.put(transcript,tsDataList)
    }
    return dataMap
}

suspend fun writeFiles(inputChannel: Channel<InputChannelData>) = withContext(Dispatchers.IO) {

	// Take data off queue, write to output
	for (transcriptData in inputChannel ) {
		val ts = transcriptData.transcript
		val entries = transcriptData.dataList
		val outputFile = transcriptData.fileName
	
		File(outputFile).bufferedWriter().use { writer ->
			for (line in entries ) {
				writer.write(line)
			}

		}

	}

}

println("Calling process data")
processData()
