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

object targetMobileElements{

	def main (args: Array[String]){
		if (args.contains("-h") || args.contains("--help") || args.size < 3){
			println("JAVA_OPTS=-Xmx4g scala -cp htsjdk-1.139.jar:. targetMobileElements <SAMPLEBAMs List> <REPEATS_DB_UCSC> <Targets File ch:start-end>")
			System.exit(0)
		}

		val Chroms = Array(("chr1",158337067),("chr2",137060424),("chr3",121430405),("chr4",120829699),("chr5",121191424),("chr6",119458736),("chr7",112638659),("chr8",113384836),("chr9",105708250),("chr10",104305016),("chr11",107310763),("chr12",91163125),("chr13",84240350),("chr14",84648390),("chr15",85296676),("chr16",81724687),("chr17",75158596),("chr18",66004023),("chr19",64057457),("chr20",72042655),("chr21",71599096),("chr22",61435874),("chr23",52530062),("chr24",62714930),("chr25",42904170),("chr26",51681464),("chr27",45407902),("chr28",46312546),("chr29",51505224),("chrX",148823899))
		
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


		val outFactory = new htsjdk.samtools.SAMFileWriterFactory
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


		/* Begin analysis of Chromosomes for LTRS */

		for (bamFile <- allBams){
		
		val fID = if (bamFile.toString.contains("_")) bamFile.toString.split("/").last.split("_").apply(0) else bamFile.toString.split("/").last.split(".").apply(0)
		
		val outEvent = new BufferedWriter(new FileWriter(new File(s"${fID}.event.tab")))
		val input = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(bamFile))

		for (win <- targetWindows.reverse){
		//print(" " + win._1)		
		
			var alignments = input.queryOverlapping(win._1,win._2,win._3)
			var typeEvent = new HashMap[String,Tuple2[Int,Int]]
			val chrs = win 
			var updated = false
			var RpMp, RpMn, RnMp, RnMn = 0
			var fwdP, revP, fwdS, revS, firstEnd, lastStart = 0
			var splitEnd = new HashMap[Int,Int]
			var hetBridge, readDepth = 0
			var windowDepth = 0
			var allIPP = 0
			var candidateWindowStart = win._2
			var candidateWindowEnd = win._3
			var windowBoo = false
			var lastfwd = false
			var c5pSplit = new HashMap[Int,Int]
			var c3pSplit = new HashMap[Int,Int]

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
					val inputBreak = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(bamFile))
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

			def updateBreaks(cur: HashMap[Int,Int], in: Int) : HashMap[Int,Int] ={
				var temp = cur
				if (cur.contains(in)) temp(in) = cur(in) + 1 else temp += in -> 1
				temp
			}

		def delDup(c5p : HashMap[Int,Int], c3p: HashMap[Int,Int]) : String ={
			if (c5p.size >= 1 && c3p.size >= 1){
					val best5p = c5p.toSeq.sortBy(-_._2).toArray.apply(0)._1
					val best3p = c3p.toSeq.sortBy(-_._2).toArray.apply(0)._1
					if (best5p < best3p) {
						"DUP"
					} else {
						if (best5p > best3p) "DEL" else "UNW"
					}
				} else{
					"UNW"
				}
		}


			def findOthers(breaks: HashMap[Int,Int]) : Unit = {
				val inputBreak = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(bamFile))
				var breakpoints = breaks.keys.toArray.sorted
					/* Extract reads covering break point if possible */
					//println("Breakpoint refinement")
				if (breakpoints.size == 2){
					//val keyvals = breakpoints.keys.toArray
					val break = inputBreak.queryOverlapping(chrs._1,breakpoints(0),breakpoints(1))
					while (break.hasNext){
						val tmpBreak = break.next
						if (tmpBreak.getAlignmentStart < breakpoints(0) && tmpBreak.getAlignmentEnd > breakpoints(1)){
							hetBridge += 1
						}
						readDepth += 1
						if (tmpBreak.getMappingQuality >= 20 && tmpBreak.getMappingQuality <= 60 && ( ( tmpBreak.getCigarString.contains("S") && tmpBreak.getStringAttribute("SA") == null) || tmpBreak.getCigarString.contains("H") ) && repeatsChrom.contains(tmpBreak.getReferenceName) && repeatsChrom.contains(tmpBreak.getMateReferenceName)) {
							if (tmpBreak.getProperPairFlag == true){
								if (tmpBreak.getReadNegativeStrandFlag == false) {
									splitEnd = updateBreaks(splitEnd,tmpBreak.getAlignmentStart)
									//if (splitEnd.contains(tmpBreak.getAlignmentStart)) splitEnd(tmpBreak.getAlignmentStart) = splitEnd(tmpBreak.getAlignmentStart) + 1 else splitEnd += tmpBreak.getAlignmentStart -> 1
									fwdS += 1
									if (c5pSplit.contains(tmpBreak.getAlignmentEnd)) c5pSplit(tmpBreak.getAlignmentEnd) = c5pSplit(tmpBreak.getAlignmentEnd) + 1 else c5pSplit += tmpBreak.getAlignmentEnd -> 1
								} else {
									splitEnd = updateBreaks(splitEnd,tmpBreak.getAlignmentEnd)
									//splitEnd += tmpBreak.getAlignmentEnd
									revS += 1
									if (c3pSplit.contains(tmpBreak.getAlignmentStart)) c3pSplit(tmpBreak.getAlignmentStart) = c3pSplit(tmpBreak.getAlignmentStart) + 1 else c3pSplit += tmpBreak.getAlignmentStart -> 1
								}
							} else {
								val curMateLTRcheck = isLTR(tmpBreak.getMateAlignmentStart,repeatsChrom(tmpBreak.getMateReferenceName))
								val curMainLTRcheck = isLTR(tmpBreak.getAlignmentStart,repeatsChrom(tmpBreak.getReferenceName))
								if (curMainLTRcheck._1 == false && curMateLTRcheck._1 == true){							
									if (tmpBreak.getReadNegativeStrandFlag == false) {
										//splitEnd += tmpBreak.getAlignmentEnd
										splitEnd = updateBreaks(splitEnd,tmpBreak.getAlignmentEnd)
										fwdS += 1
										if (c5pSplit.contains(tmpBreak.getAlignmentEnd)) c5pSplit(tmpBreak.getAlignmentEnd) = c5pSplit(tmpBreak.getAlignmentEnd) + 1 else c5pSplit += tmpBreak.getAlignmentEnd -> 1
									} else {
										if (c3pSplit.contains(tmpBreak.getAlignmentStart)) c3pSplit(tmpBreak.getAlignmentStart) = c3pSplit(tmpBreak.getAlignmentStart) + 1 else c3pSplit += tmpBreak.getAlignmentStart -> 1
										splitEnd = updateBreaks(splitEnd,tmpBreak.getAlignmentStart)
										revS += 1
									}
								}
							}
						}
					}
				} else {
					if (breaks.size > 2){
						val sortedBreaks = breaks.toSeq.sortBy(- _._2).toArray
						findOthers(HashMap(sortedBreaks(0), sortedBreaks(1)))
					}
					
				}
				inputBreak.close
			}


			/* Process all alignments on a Chromosome */
			while (alignments.hasNext){
				val tmp = alignments.next		
				windowDepth += 1
				/* PRoper window identification, Check if on good Chromosomes */
				if (repeatsChrom.contains(tmp.getReferenceName) && repeatsChrom.contains(tmp.getMateReferenceName)){

					if (tmp.getProperPairFlag == false && tmp.getMappingQuality >= 20 && tmp.getMappingQuality <= 60){
						allIPP += 1

						val curMateLTRcheck = isLTR(tmp.getMateAlignmentStart,repeatsChrom(tmp.getMateReferenceName))
						val curMainLTRcheck = isLTR(tmp.getAlignmentStart,repeatsChrom(tmp.getReferenceName))
						//if (curMainLTRcheck._1 == false && curMateLTRcheck._1 == true){
						if (curMateLTRcheck._1 == true){
							if (windowBoo == false){
								if (tmp.getReadNegativeStrandFlag == false){
									firstEnd = tmp.getAlignmentEnd
									windowBoo = true
									fwdP += 1
									lastfwd = true
									directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
									typeEvent += curMateLTRcheck._3 -> (1,0)
								} //end if fwd read

							} else {
								/* Window is True, improperly paired and good mapping quality */
								//lastStart = tmp.getAlignmentStart
								if (tmp.getStringAttribute("SA") == null ) { //|| tmp.getSupplementaryAlignmentFlag == false){ Not supllemantary alignment
									if (tmp.getReadNegativeStrandFlag == false){ // Not sup and on +ve strand
										if (lastfwd) firstEnd = tmp.getAlignmentEnd
										if (typeEvent.contains(curMateLTRcheck._3)) typeEvent(curMateLTRcheck._3) = (typeEvent(curMateLTRcheck._3)._1 + 1, typeEvent(curMateLTRcheck._3)._2) else typeEvent += curMateLTRcheck._3 -> (1,0)
										directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
										fwdP += 1
									} else {
										if (lastfwd) {
											lastStart = tmp.getAlignmentStart
											lastfwd = false
										}
										if (typeEvent.contains(curMateLTRcheck._3)) typeEvent(curMateLTRcheck._3) = (typeEvent(curMateLTRcheck._3)._1, typeEvent(curMateLTRcheck._3)._2 + 1) else typeEvent += curMateLTRcheck._3 -> (0,1)
										directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
										revP += 1
									}//end fwd read npp
								} else { // IS Sup Alignment/splitRead
									if (tmp.getReadNegativeStrandFlag == false){
											if (typeEvent.contains(curMateLTRcheck._3)) typeEvent(curMateLTRcheck._3) = (typeEvent(curMateLTRcheck._3)._1 + 1, typeEvent(curMateLTRcheck._3)._2) else typeEvent += curMateLTRcheck._3 -> (1,0)
											fwdS += 1
											fwdP += 1
											if (c5pSplit.contains(tmp.getAlignmentEnd)) c5pSplit(tmp.getAlignmentEnd) = c5pSplit(tmp.getAlignmentEnd) + 1 else c5pSplit += tmp.getAlignmentEnd -> 1
											directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
											splitEnd = updateBreaks(splitEnd,tmp.getAlignmentEnd)
											//splitEnd += tmp.getAlignmentEnd
										} else {
											if (typeEvent.contains(curMateLTRcheck._3))typeEvent(curMateLTRcheck._3) =  (typeEvent(curMateLTRcheck._3)._1, typeEvent(curMateLTRcheck._3)._2 + 1) else typeEvent += curMateLTRcheck._3 -> (0,1)
											revS += 1
											revP += 1
											if (c3pSplit.contains(tmp.getAlignmentStart)) c3pSplit(tmp.getAlignmentStart) = c3pSplit(tmp.getAlignmentStart) + 1 else c3pSplit += tmp.getAlignmentStart -> 1
											directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
											splitEnd = updateBreaks(splitEnd,tmp.getAlignmentStart)
											//splitEnd += tmp.getAlignmentStart
										} // end sup is neg read
								} // end Sup Else
							} // end if window true
						} // is Repeat structure of mates fine?

					} // end good read mapping and improperly paired

					//if (windowBoo == true) output.addAlignment(tmp)
				}// End of Good Chromosomes check
			} //End of While Alignments
			alignments.close
			/* Have completed a chromosome scan now need to clean up*/
			if (windowBoo){
				val targets = splitEnd.toSeq.sortBy(- _._2).toArray
				if (splitEnd.size == 2 && fwdP > 0 && revP > 0) findOthers(splitEnd) else if (splitEnd.size == 1 && fwdP > 0 && revP > 0) findOthers(HashMap((targets(0)._1 - 10,0), (targets(0)._1 + 10,0))) else if (splitEnd.size == 0 && fwdP > 4 && revP > 4) findOthers(HashMap(firstEnd -> 0,lastStart -> 0))
				var splitDif = if (splitEnd.size == 2) scala.math.abs(splitEnd.toArray.apply(0)._1 - splitEnd.toArray.apply(1)._1) else 0
				if (splitDif > 20){
					findOthers(HashMap(win._2 -> 0,win._3 -> 0))
					splitDif = if (splitEnd.size == 2) scala.math.abs(splitEnd.toArray.apply(0)._1 - splitEnd.toArray.apply(1)._1) else 0
				}
				println(s"fwd:${fwdP} rev:${revP} splitDif:${splitDif}")
				if (fwdP >= 1 && revP >= 1 && (fwdP + revP) >= 3 && splitDif <= 20) {
						val best5p = if (c5pSplit.size >= 1) c5pSplit.toSeq.sortBy(-_._2).toArray.apply(0)._1 else -1
						val best3p = if (c3pSplit.size >= 1) c3pSplit.toSeq.sortBy(-_._2).toArray.apply(0)._1 else -1
						val ornt = s"${RpMp}:${RpMn}:${RnMp}:${RnMn}"
						outEvent.write(s"${chrs._1}:${candidateWindowStart}-${candidateWindowEnd}\t${fwdP}\t${revP}\t${fwdS}\t${revS}\t${splitEnd.keys.toSet}\t${splitDif}\t${ornt}\t")
						typeEvent.foreach(s => outEvent.write(s"${s._1},${s._2._1},${s._2._2}:"))
						outEvent.write(s"\t${readDepth}:${hetBridge}\t${windowDepth}\t${allIPP - (fwdP+revP)}\t${(fwdP+revP)/(allIPP.toFloat - (fwdP+revP).toFloat + 0.000001)}\t${delDup(c5pSplit,c3pSplit)}\t${best5p}:${best3p}\n")
				}
			}

		} // end for chr from targets
		outEvent.close
		input.close
	}// End loop through BAMS
	} // def main

} //object