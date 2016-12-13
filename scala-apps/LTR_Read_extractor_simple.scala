import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.SAMFileWriter
import htsjdk.samtools.SAMFileWriterFactory
import htsjdk.samtools.DefaultSAMRecordFactory
import htsjdk.samtools.ValidationStringency
import htsjdk.samtools.SamReader
import htsjdk.samtools.SAMRecord
import java.io.File
import java.io._
import scala.math._
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet

object targetExtract{

	def main (args: Array[String]){
		if (args.contains("-h") || args.contains("--help") || args.size != 2){
			println("JAVA_OPTS=-Xmx4g scala -cp htsjdk-1.139.jar:. findMobileElements <SAMPLEBAMs.txt> <Targets File ch:start-end>")
			System.exit(0)
		}

		val Chroms = Array(("chr1",158337067),("chr2",137060424),("chr3",121430405),("chr4",120829699),("chr5",121191424),("chr6",119458736),("chr7",112638659),("chr8",113384836),("chr9",105708250),("chr10",104305016),("chr11",107310763),("chr12",91163125),("chr13",84240350),("chr14",84648390),("chr15",85296676),("chr16",81724687),("chr17",75158596),("chr18",66004023),("chr19",64057457),("chr20",72042655),("chr21",71599096),("chr22",61435874),("chr23",52530062),("chr24",62714930),("chr25",42904170),("chr26",51681464),("chr27",45407902),("chr28",46312546),("chr29",51505224),("chrX",148823899))
		
		val inputBams = new BufferedReader(new FileReader(new File(args(0))))
		val inputTargets = new BufferedReader(new FileReader(new File(args(1))))
		
		
		//Chrom,Start,End, Breakpoint1, Breakpoint2
		var targetWindows : List[Tuple3[String,Int,Int]] = Nil
		
		while (inputTargets.ready){
			val tmp = inputTargets.readLine.split("\t")
			val chrom = tmp(0).split(":").apply(0).trim
			val start = tmp(0).split(":").apply(1).split("-").apply(0).toInt
			val end = tmp(0).split(":").apply(1).split("-").apply(1).toInt
						
			targetWindows = (chrom, start, end) :: targetWindows
		}
		
		//println(targetWindows)
		
		inputTargets.close

		var allBams : List[String] = Nil
		
		while (inputBams.ready){
			val tmp = inputBams.readLine
			allBams = tmp :: allBams
		}

//		val inputBreak = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(args(1)))

		val outFactory = new htsjdk.samtools.SAMFileWriterFactory
		//Write header too new bam

		for (win <- targetWindows.reverse){
		val inputLocal = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(allBams(0)))
		val output = outFactory.makeSAMOrBAMWriter(inputLocal.getFileHeader,false,new File(s"${win._1 + "_" + win._2}.discordant.bam"))
		inputLocal.close
		/* Begin analysis of Chromosomes for LTRS */

		for (bamFile <- allBams){
		
			val fID = if (bamFile.toString.contains("_")) bamFile.toString.split("/").last.split("_").apply(0) else bamFile.toString.split("/").last.split('.').apply(0)
			val input = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(bamFile))
			val mates = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(bamFile))
					
					var alignments = input.queryOverlapping(win._1,win._2,win._3)
					var typeEvent = new HashMap[String,Tuple2[Int,Int]]
					val chrs = win 
					var candidateWindowStart = win._2
					var candidateWindowEnd = win._3

					/* Process all alignments on a Chromosome */
						while (alignments.hasNext){
							val tmp = alignments.next		

							if (tmp.getReadPairedFlag == true){
								   try { 
									output.addAlignment(tmp)
									output.addAlignment(mates.queryMate(tmp))

								   } catch {
								     case ioe: IOException => System.err.println("Error")
								     case e: Exception => System.err.println("Error")
								   }
							}
							
						} //End of While Alignments
						alignments.close
					/* Have completed a chromosome scan now need to clean up*/
					input.close
					mates.close
					
					} // end for chr from targets
		output.close
	}// End loop through BAMS
	} // def main

} //object