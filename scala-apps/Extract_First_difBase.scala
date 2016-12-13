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

object targetMobileElementsExtract{

	def main (args: Array[String]){
		if (args.contains("-h") || args.contains("--help") || args.size < 3){
			println("JAVA_OPTS=-Xmx4g scala -cp htsjdk-1.139.jar:. findMobileElements <SAMPLEBAMs List> <REPEATS_DB_UCSC> <Targets File ch:start-end>")
			System.exit(0)
		}

		//val Chroms = Array(("chr1",158337067),("chr2",137060424),("chr3",121430405),("chr4",120829699),("chr5",121191424),("chr6",119458736),("chr7",112638659),("chr8",113384836),("chr9",105708250),("chr10",104305016),("chr11",107310763),("chr12",91163125),("chr13",84240350),("chr14",84648390),("chr15",85296676),("chr16",81724687),("chr17",75158596),("chr18",66004023),("chr19",64057457),("chr20",72042655),("chr21",71599096),("chr22",61435874),("chr23",52530062),("chr24",62714930),("chr25",42904170),("chr26",51681464),("chr27",45407902),("chr28",46312546),("chr29",51505224),("chrX",148823899))
		
		val inputBams = new BufferedReader(new FileReader(new File(args(0))))
		val inputTargets = new BufferedReader(new FileReader(new File(args(2))))
		
		
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
		val repeats = new BufferedReader(new FileReader(new File(args(1))))


		//val outFactory = new htsjdk.samtools.SAMFileWriterFactory
		//Write header too new bam

		/* Build repeat structures and load data into them*/		

		val repeatsChrom = new HashMap[String,Array[Tuple6[Int,Int,String,String,String,String]]]
		var repList: List[Tuple6[Int,Int,String,String,String,String]] = Nil
		var prevChrom = ""

		while (repeats.ready){
			val tmp = repeats.readLine.split("\t")
			if (prevChrom == tmp(5)){
				repList = (tmp(6).toInt,tmp(7).toInt,tmp(10),tmp(11),tmp(12),tmp(9)) :: repList
			} else {
				repeatsChrom += prevChrom -> repList.reverse.toArray
				repList = (tmp(6).toInt,tmp(7).toInt,tmp(10),tmp(11),tmp(12),tmp(9)) :: Nil
				prevChrom = tmp(5)
			}
		}//end while repeats

		repeatsChrom += prevChrom -> repList.reverse.toArray
		repList = Nil
		prevChrom = ""
		repeats.close

		def isLTR (pos: Int, curSearch: Array[Tuple6[Int,Int,String,String,String,String]]) : Tuple3[Boolean,String,String] = {
			val half = scala.math.floor(curSearch.length/2).toInt
			if (curSearch.length == 0){
					(false,"","")
			} else {
				if (curSearch(half)._1 <= pos && curSearch(half)._2 >= pos) {
					(true,curSearch(half)._6,s"${curSearch(half)._3}-${curSearch(half)._5}")
					} else {
						if (curSearch.length <= 1) {
							(false,"","")
						} else {
						if (curSearch(half)._1 > pos) isLTR(pos,curSearch.slice(0,half)) else isLTR(pos,curSearch.slice(half + 1,curSearch.length))
					}
				}
			}
		} //end is LTR

		var data = new HashMap[String, HashMap[String,Int]]


		/* Begin analysis of Chromosomes for LTRS */

		for (bamFile <- allBams){
		
			val fID = if (bamFile.toString.contains("_")) bamFile.toString.split("/").last.split("_").apply(0) else bamFile.toString.split("/").last.split(".").apply(0)
			
			//val outEvent = new BufferedWriter(new FileWriter(new File(s"${fID}.event.tab")))
			val input = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(bamFile))
			val mates = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(bamFile))
			//val output = outFactory.makeSAMOrBAMWriter(input.getFileHeader,false,new File(s"${fID}.discordant.bam"))
			//val outLeft = new BufferedWriter(new FileWriter(new File("ERV_5p.tab")))
			//val outRight = new BufferedWriter(new FileWriter(new File("ERV_3p.tab")))



			for (win <- targetWindows.reverse){
			//print(" " + win._1)		
			
				var alignments = input.queryOverlapping(win._1,win._2 - 5,win._3 + 5)
				var typeEvent = new HashMap[String,Tuple2[Int,Int]]
				val chrs = win 
				var updated = false
				//var LTRpRn, LTRpRp, LTRnRp, LTRnRn = 0
				var RpMp, RpMn, RnMp, RnMn = 0
				var fwdP, revP, fwdS, revS, firstEnd, lastStart = 0
				var splitEnd = new HashSet[Int]
				var candidateWindowStart = win._2
				var candidateWindowEnd = win._3
				var windowBoo = false
				var lastfwd = false

				/*  Chr_Pos -> 5'/3' -> Sequence -> Count */

				/* Process all alignments on a Chromosome */
				while (alignments.hasNext){
					val tmp = alignments.next		

					/* PRoper window identification, Check if on good Chromosomes */
					if (repeatsChrom.contains(tmp.getReferenceName) && repeatsChrom.contains(tmp.getMateReferenceName)){

						//if (tmp.getProperPairFlag == false && tmp.getMappingQuality =< 60){

						//	val curMateLTRcheck = isLTR(tmp.getMateAlignmentStart,repeatsChrom(tmp.getMateReferenceName))
						//	val curMainLTRcheck = isLTR(tmp.getAlignmentStart,repeatsChrom(tmp.getReferenceName))
						//	if (curMainLTRcheck._1 == false && curMateLTRcheck._1 == true){

								val cigar = tmp.getCigarString
								val read = tmp.getReadString
								var fwdRef, fwdErv, revRef, revErv, fwdData, revData = ""

								if (cigar.contains("S") && cigar.contains("M") && ! cigar.contains("I") && ! cigar.contains("D") && cigar.map(s => if (s == 'M') 1 else 0).sum == 1){

									if(cigar.indexOf("S") < cigar.indexOf("M")){
										//case Softclip 5'
										val bases = cigar.substring(0,cigar.indexOf("S")).toInt
										if (bases >= 5){
											//outLeft.write( win._1 + ":" + win._2 + "-" + win._3 + "\tERV:" + read.substring(bases -5, bases) + "__" + read.substring(bases,bases + 20) + "\n")
											fwdRef = if (read.size >= bases + 20) read.substring(bases,bases + 20) else read.substring(bases,read.size)
											fwdErv = if (bases -5 >= 0) read.substring(bases -5, bases) else read.substring(0, bases)
											fwdData = "fwd\t" + fwdErv + "\t" + fwdRef
											//println(win._1 + ":" + win._2 + " " + fwdRef + " " + fwdErv + " " + fwdData)
										} else {
											//for (i <- 0 to (scala.math.abs(bases - 10) -1)) outLeft.write("N")
											//outLeft.write(win._1 + ":" + win._2 + "-" + win._3 + "\tERV:" + read.substring(0, bases) + "__" + read.substring(bases,bases + 20) +  "\n")
										}

									} else {
											//case softclip 3'
											val bases = cigar.substring(0,cigar.indexOf("M")).toInt
											val split = cigar.substring(cigar.indexOf("M")+1,cigar.indexOf("S")).toInt
											if (split >= 5){
												//outRight.write(win._1 + ":" + win._2 + "-" + win._3 + "\t" + read.substring(bases - 20, bases) + "__ERV:" + read.substring(bases,bases + 5) + "\n")
												revRef = if (bases - 20 >= 0) read.substring(bases -20, bases) else read.substring(0,bases)
												revErv = if (read.size >= bases +5 ) read.substring(bases,bases + 5) else read.substring(bases,read.size)
												revData = "rev\t" + revRef + "\t" + revErv
												//println(win._1 + ":" + win._2 + " " + revErv + " " + revRef + " " + revData)

											} else {
												//outRight.write(win._1 + ":" + win._2 + "-" + win._3 + "\t" + read.substring(bases - 20, bases) + "__ERV:" + read.substring(bases,bases + split))
												//for (i <- 0 to (scala.math.abs(split - 20) -1)) outRight.write("N")
												//outRight.write("\n")
											}

										}

								}
								/*  Chr_Pos -> 5'/3' -> Sequence -> Count */
								val key = win._1 + ":" + win._2 + "-" + win._3
								//val fwdData = "fwd\t" + fwdErv + "\t" + fwdRef
								//val revData = "rev\t" + revRef + "\t" + revErv

								if (data.contains(key)){
									if (data(key).contains(fwdData)) data(key)(fwdData) = data(key)(fwdData) + 1 else data(key) += fwdData -> 1
									if (data(key).contains(revData)) data(key)(revData) = data(key)(revData) + 1 else data(key) += revData -> 1
								} else {
									data += key -> HashMap(fwdData -> 1)
									data(key) += revData -> 1
								}

							} // is Repeat structure of mates fine?

					//} // end good read mapping and improperly paired

						//if (windowBoo == true) output.addAlignment(tmp)
					//}// End of Good Chromosomes check
				} //End of While Alignments
				alignments.close

			} // end for chr from targets


			mates.close
			//output.close
			//outEvent.close
			input.close
			//outRight.close
			//outLeft.close
		}// End loop through BAMS

				for ( i <- data.keys){
					val revs = data(i).filter(s => s._1.contains("rev"))
					val fwds = data(i).filter(s => s._1.contains("fwd"))

					val sortedRevs = revs.toSeq.sortBy(- _._2).toArray
					val sortedFwds = fwds.toSeq.sortBy(- _._2).toArray

					//println(sortedRevs(0))
					//println(sortedFwds(0))

					if (sortedRevs.size >= 1) if (sortedRevs(0)._1.split("\t").apply(1) == "" && sortedRevs.size >= 2) print(i + "\t" + sortedRevs(1)._1.split("\t")(1) + "\t" + sortedRevs(1)._1.split("\t")(2)) else print(i + "\t" + sortedRevs(0)._1.split("\t")(1) + "\t" + sortedRevs(0)._1.split("\t")(2))
					if (sortedFwds.size >= 1) if (sortedFwds(0)._1.split("\t").apply(1) == "" && sortedFwds.size >= 2) print("\t" + sortedFwds(1)._1.split("\t")(1) + "\t" + sortedFwds(1)._1.split("\t")(2)) else print("\t" + sortedFwds(0)._1.split("\t")(1) + "\t" + sortedFwds(0)._1.split("\t")(2))
					if (sortedRevs.size >= 1 || sortedFwds.size >= 1) print("\n")
				}
	} // def main

} //object