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

object findMobileElements{

	def main (args: Array[String]){
		if (args.contains("-h") || args.contains("--help") || args.size < 3){
			println("JAVA_OPTS=-Xmx4g scala -cp htsjdk-1.139.jar:. findMobileElements <SAMPLENAME> <SAMPLEBAM> <REPEATS_DB_UCSC> <UNIQUE true/false>")
			System.exit(0)
		}

		val Chroms = Array(("chr1",158337067),("chr2",137060424),("chr3",121430405),("chr4",120829699),("chr5",121191424),("chr6",119458736),("chr7",112638659),("chr8",113384836),("chr9",105708250),("chr10",104305016),("chr11",107310763),("chr12",91163125),("chr13",84240350),("chr14",84648390),("chr15",85296676),("chr16",81724687),("chr17",75158596),("chr18",66004023),("chr19",64057457),("chr20",72042655),("chr21",71599096),("chr22",61435874),("chr23",52530062),("chr24",62714930),("chr25",42904170),("chr26",51681464),("chr27",45407902),("chr28",46312546),("chr29",51505224),("chrX",148823899))
		

		val input = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(args(1)))
		//var alignments = input.queryAlignmentStart(Chroms(0)._1,1)
		val inputBreak = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(args(1)))
		//var alignmentsBreak = inputBreak.queryAlignmentStart(Chroms(0)._1,1)
		val repeats = new BufferedReader(new FileReader(new File(args(2))))


		val outFactory = new htsjdk.samtools.SAMFileWriterFactory
		//Write header too new bam
		val output = outFactory.makeSAMOrBAMWriter(input.getFileHeader,true,new File(s"${args(0)}.bam"))
		val outFWD = new BufferedWriter(new FileWriter(new File(s"${args(0)}.bedgraph")))
		val outEvent = new BufferedWriter(new FileWriter(new File(s"${args(0)}.event.tab")))
		input.close

		outFWD.write("track type=bedGraph\n")

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


		/* Begin analysis of Chromosomes for LTRS */


		for (chrs <- Chroms){
		print(" " + chrs._1)
		val input = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(args(1)))
		var alignments = input.query(chrs._1,1,200000000,true)
		var typeEvent = new HashMap[String,Tuple2[Int,Int]]

		var updated = false
		//var LTRpRn, LTRpRp, LTRnRp, LTRnRn = 0
		var RpMp, RpMn, RnMp, RnMn = 0
		var fwdP, revP, fwdS, revS = 0
		var splitEnd = new HashSet[Int]
		var candidateWindowStart, candidateWindowEnd = 0
		var windowBoo = false

		def directionLTR (LTR: String, uniqueNeg: Boolean, mateNeg: Boolean) : Unit = {
			if (uniqueNeg){ 
				// Negative Strand Read
					if (LTR == "+"){
						//LTR is on + strand	
						if (mateNeg) RnMn += 1 else RnMp += 1 
						// Is the mate negative or positive
					} else {
						// LTR is on - strand
						if (mateNeg) RnMp += 1 else RnMn += 1
					}
				} else {
					// Positive Strand Read
					if (LTR == "+"){
						//LTR is on + strand	
						if (mateNeg) RpMn += 1 else RpMp += 1 
						// Is the mate negative or positive
					} else {
						// LTR is on - strand
						if (mateNeg) RpMp += 1 else RpMn += 1
					}
				}
		}

		def updateCounts(breaks: HashSet[Int]) : Unit = {
			if (breaks.size == 2){
				val inputBreak = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(args(1)))
				var breakpoints = breaks.toArray.sorted
				/* Extract reads covering break point if possible */
				//println("Breakpoint refinement")
				val break = inputBreak.queryOverlapping(chrs._1,breakpoints(0),breakpoints(1))
				while (break.hasNext){
					val tmpBreak = break.next
					if ((breakpoints.contains(tmpBreak.getAlignmentStart) || breakpoints.contains(tmpBreak.getAlignmentEnd)) && tmpBreak.getProperPairFlag == true && tmpBreak.getMappingQuality == 60 && tmpBreak.getCigarString.contains("S")) {
						if (tmpBreak.getReadNegativeStrandFlag == false) fwdS += 1 else revS += 1
					}
				}
			}
		}


		/* Process all alignments on a Chromosome */
		while (alignments.hasNext){
			val tmp = alignments.next		

			/* If the next read is outside of the current region then write out region and advance */
			if (candidateWindowEnd != 0 && tmp.getAlignmentStart >= candidateWindowEnd){
				val splitDif = if (splitEnd.size == 2) scala.math.abs(splitEnd.toArray.apply(0) - splitEnd.toArray.apply(1)) else 0
				if (fwdP > 0 && fwdS > 0 && (revP + revS) > 0 && splitDif <= 20) {
					//println(s"Criteria for print ${chrs._1}:${candidateWindowStart}-${candidateWindowEnd}\t${fwdP}\t${revP}\t${fwdS}\t${revS}\t${splitEnd}")
					val targets = splitEnd.toArray.sorted
					if (splitEnd.size == 2 ) updateCounts(splitEnd) else if (splitEnd.size == 1) updateCounts(HashSet(targets(0) - 10, targets(0) + 10))
					//println(s"PASSED ${chrs._1}:${candidateWindowStart}-${candidateWindowEnd}\t${fwdP}\t${revP}\t${fwdS}\t${revS}\t${splitEnd}")
					val ornt = s"${RpMp}:${RpMn}:${RnMp}:${RnMn}"
					outEvent.write(s"${chrs._1}:${candidateWindowStart}-${candidateWindowEnd}\t${fwdP}\t${revP}\t${fwdS}\t${revS}\t${splitEnd}\t${splitDif}\t${ornt}\t")
					typeEvent.foreach(s => outEvent.write(s"${s._1},${s._2._1},${s._2._2}:"))
					outEvent.write("\n")
					outFWD.write(s"${chrs._1}\t${candidateWindowStart}\t${candidateWindowEnd}\t${fwdP + fwdS + revP + revS}\n")
				} // end of write

				candidateWindowEnd = 0
				candidateWindowStart = 0
				windowBoo = false
				fwdP = 0
				fwdS = 0
				revS = 0
				revP = 0
				splitEnd = new HashSet[Int]
				typeEvent = new HashMap[String,Tuple2[Int,Int]]
				RpMn = 0
				RpMp = 0
				RnMp = 0
				RnMn = 0
			} // end of window work


			/* PRoper window identification, Check if on good Chromosomes */
			if (repeatsChrom.contains(tmp.getReferenceName) && repeatsChrom.contains(tmp.getMateReferenceName)){

				if (tmp.getProperPairFlag == false && tmp.getMappingQuality == 60){

					val curMateLTRcheck = isLTR(tmp.getMateAlignmentStart,repeatsChrom(tmp.getMateReferenceName))
					val curMainLTRcheck = isLTR(tmp.getAlignmentStart,repeatsChrom(tmp.getReferenceName))

					if (curMainLTRcheck._1 == false && curMateLTRcheck._1 == true){

						if (windowBoo == false){
							if (tmp.getReadNegativeStrandFlag == false){
								//output.addAlignment(tmp)
								//println(" create new window " + tmp.getAlignmentStart)
								windowBoo = true
								fwdP += 1
								directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
								typeEvent += curMateLTRcheck._3 -> (1,0)
								candidateWindowStart = ((tmp.getAlignmentStart / 1000)*1000)
								// 10300 = /1k*1k = 10,000 - 10300 = -300 vs 200 convert to absolute
								if (abs(candidateWindowStart - tmp.getAlignmentStart) >= abs((candidateWindowStart + 500) - tmp.getAlignmentStart)) candidateWindowStart += 500
								candidateWindowEnd = candidateWindowStart + 1500
								//println(candidateWindowStart + " " + candidateWindowEnd)
							} //end if fwd read

						} else {
							/* Window is True, improperly paired and good mapping quality */
							//output.addAlignment(tmp)
							if (tmp.getStringAttribute("SA") == null ) { //|| tmp.getSupplementaryAlignmentFlag == false){ Not supllemantary alignment
								if (tmp.getReadNegativeStrandFlag == false){ // Not sup and on +ve strand
									if (typeEvent.contains(curMateLTRcheck._3)) typeEvent(curMateLTRcheck._3) = (typeEvent(curMateLTRcheck._3)._1 + 1, typeEvent(curMateLTRcheck._3)._2) else typeEvent += curMateLTRcheck._3 -> (1,0)
									directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
									fwdP += 1
								} else {
									if (typeEvent.contains(curMateLTRcheck._3)) typeEvent(curMateLTRcheck._3) = (typeEvent(curMateLTRcheck._3)._1, typeEvent(curMateLTRcheck._3)._2 + 1) else typeEvent += curMateLTRcheck._3 -> (0,1)
									directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
									revP += 1
								}//end fwd read npp
							} else { // IS Sup Alignment/splitRead
								//println(" is getSupplementaryAlignmentFlag")
								if (tmp.getReadNegativeStrandFlag == false){
										if (typeEvent.contains(curMateLTRcheck._3)) typeEvent(curMateLTRcheck._3) = (typeEvent(curMateLTRcheck._3)._1 + 1, typeEvent(curMateLTRcheck._3)._2) else typeEvent += curMateLTRcheck._3 -> (1,0)
										fwdS += 1
										directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
										splitEnd += tmp.getAlignmentEnd
									} else {
										if (typeEvent.contains(curMateLTRcheck._3))typeEvent(curMateLTRcheck._3) =  (typeEvent(curMateLTRcheck._3)._1, typeEvent(curMateLTRcheck._3)._2 + 1) else typeEvent += curMateLTRcheck._3 -> (0,1)
										revS += 1
										directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
										splitEnd += tmp.getAlignmentStart
									} // end sup is neg read
							} // end Sup Else
						} // end if window true
					} // is Repeat structure of mates fine?

				} // end good read mapping and improperly paired

				//if (windowBoo == true) output.addAlignment(tmp)
			}// End of Good Chromosomes check
		} //End of While Alignments
		/* Have completed a chromosome scan now need to clean up*/
		if (windowBoo){
			val splitDif = if (splitEnd.size == 2) scala.math.abs(splitEnd.toArray.apply(0) - splitEnd.toArray.apply(1)) else 0
			if (fwdP > 0 && fwdS > 0 && (revP + revS) > 0 && splitDif <= 20) {
					val targets = splitEnd.toArray.sorted
					if (splitEnd.size == 2 ) updateCounts(splitEnd) else if (splitEnd.size == 1) updateCounts(HashSet(targets(0) - 10, targets(0) + 10))
					val ornt = s"${RpMp}:${RpMn}:${RnMp}:${RnMn}"
					outEvent.write(s"${chrs._1}:${candidateWindowStart}-${candidateWindowEnd}\t${fwdP}\t${revP}\t${fwdS}\t${revS}\t${splitEnd}\t${splitDif}\t${ornt}\t")
					typeEvent.foreach(s => outEvent.write(s"${s._1},${s._2._1},${s._2._2}:"))
					outEvent.write("\n")
					outFWD.write(s"${chrs._1}\t${candidateWindowStart}\t${candidateWindowEnd}\t${fwdP + fwdS + revP + revS}\n")
			}
		}
		input.close

		} // end for chr from Chroms
		outEvent.close
		outFWD.close
		output.close
	} // def main

} //object