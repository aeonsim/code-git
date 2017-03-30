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
import scala.language.postfixOps

object findMobileElements{

	def main (args: Array[String]){
		if (args.contains("-h") || args.contains("--help") || args.size < 3){
			println("JAVA_OPTS=-Xmx4g scala -cp htsjdk-1.139.jar:. findMobileElements <SAMPLENAME> <SAMPLEBAM> <REPEATS_DB_UCSC>")
			System.exit(0)
		}

		val Chroms = Array(("chr1",158337067),("chr2",137060424),("chr3",121430405),("chr4",120829699),("chr5",121191424),("chr6",119458736),("chr7",112638659),("chr8",113384836),("chr9",105708250),("chr10",104305016),("chr11",107310763),("chr12",91163125),("chr13",84240350),("chr14",84648390),("chr15",85296676),("chr16",81724687),("chr17",75158596),("chr18",66004023),("chr19",64057457),("chr20",72042655),("chr21",71599096),("chr22",61435874),("chr23",52530062),("chr24",62714930),("chr25",42904170),("chr26",51681464),("chr27",45407902),("chr28",46312546),("chr29",51505224),("chrX",148823899))
		

		val input = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(args(1)))
		//var alignments = input.queryAlignmentStart(Chroms(0)._1,1)
		val inputBreak = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(args(1)))
		//var alignmentsBreak = inputBreak.queryAlignmentStart(Chroms(0)._1,1)
		val repeats = new BufferedReader(new FileReader(new File(args(2))))
		val gap_in = new BufferedReader(new FileReader(new File("/home/u/u218757/kos_home/refs/gaps.bed")))

		val outFactory = new htsjdk.samtools.SAMFileWriterFactory
		//Write header too new bam
		//val output = outFactory.makeSAMOrBAMWriter(input.getFileHeader,true,new File(s"${args(0)}.bam"))
		//val outFWD = new BufferedWriter(new FileWriter(new File(s"${args(0)}.bedgraph")))
		val outEvent = new BufferedWriter(new FileWriter(new File(s"${args(0)}.event.tab")))
		input.close

		//outFWD.write("track type=bedGraph\n")

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


		/* Build Gaps structure */

		val gapsChrom = new HashMap[String,Array[Tuple2[Int,Int]]]
		var gapList: List[Tuple2[Int,Int]] = Nil
		//var prevChrom = ""

		while(gap_in.ready){
			val tmp = gap_in.readLine.split("\t")
			if (prevChrom == "") prevChrom = tmp(0)

			if (tmp(0) == prevChrom){
				gapList = (tmp(1).toInt - 1000,tmp(2).toInt + 1000) :: gapList
			} else {
				gapsChrom += prevChrom -> gapList.sortBy(_._1).toArray
				//println( prevChrom + " > " + tmp(0) + " " + gapList.size + " " + gapsChrom(prevChrom).size)
				prevChrom = tmp(0)
				gapList = Nil
				gapList = (tmp(1).toInt - 1000,tmp(2).toInt + 1000) :: gapList
			}
		}
		gapsChrom += prevChrom -> gapList.sortBy(_._1).toArray
		gap_in.close




	val expReads = 4035
	val expIPP = 88
	val expSR = 24

	def bootStrap(reads: Int, evid: Int, geno: Int, avg: Int): Double = {
		val ran = scala.util.Random
		val iterations = 10000
	    val dataPool = Array.fill[String](geno)("R") ++ Array.fill[String](avg)("A")
	    val tAR = evid.toDouble / reads.toDouble
	    var bigOrEq, smallOrEq = 0
	    var tmp = 0
	    if (dataPool.size >= 1) {
	      while (tmp < iterations) {
	        var ref = 0
	        var alt = 0
	        var ctr = 0
	        while (ctr < reads) {
	          if (dataPool(ran.nextInt(dataPool.size)) == "R") ref += 1 else alt += 1
	          ctr += 1
	        }
	        if (alt >= evid) bigOrEq += 1
	        tmp += 1
	      }

	      if (bigOrEq == 0) 1 / iterations.toDouble else bigOrEq / iterations.toDouble
	    } else {
	      1
	    }


	}




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


/* Is Gap? */
		def isGAP (pos: Int, curSearch: Array[Tuple2[Int,Int]]) : Tuple3[Boolean,String,String] = {
			if (pos > 1){
				val half = scala.math.floor(curSearch.length/2).toInt
				if (curSearch.length == 0){
						//println(pos + " no gap")
						(false,"","")
				} else {
					if (curSearch(half)._1 <= pos && curSearch(half)._2 >= pos) {
						//println(pos + " gap " + curSearch(half)._1  + " " + curSearch(half)._2)
						(true,"gap","gap")
						} else {
							if (curSearch.length <= 1) {
								//println(pos + " nogap " + curSearch(half)._1  + " " + curSearch(half)._2)
								(false,"","")
							} else {
							if (curSearch(half)._1 > pos) isGAP(pos,curSearch.slice(0,half)) else isGAP(pos,curSearch.slice(half + 1,curSearch.length))
						}
					}
				}
				} else {
					(false,"","")
				}
		}


		/* Begin analysis of Chromosomes for LTRS */


		for (chrs <- Chroms){
			print(" " + chrs._1)
			val input = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(args(1)))
			val mates = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(args(1)))
			var alignments = input.query(chrs._1,1,200000000,true)
			var typeEvent = new HashMap[String,Tuple2[Int,Int]]

			var updated = false
			//var LTRpRn, LTRpRp, LTRnRp, LTRnRn = 0
			var RpMp, RpMn, RnMp, RnMn = 0
			var fwdP, revP, fwdS, revS, firstEnd, lastStart = 0
			var hetBridge, readDepth = 0
			var splitEnd = new HashSet[Int]
			var windowDepth = 0
			var allIPP = 0
			var candidateWindowStart, candidateWindowEnd = 0
			var windowBoo = false
			var lastfwd = false
			var c5pSplit = new HashMap[Int,Int]
			var c3pSplit = new HashMap[Int,Int]
			var elLength : List[Int] = Nil
			var polyA = false 
			var gapBo = false
			var mateReads : List[String] = Nil
			var gapStart = -10000
			var gapChr = ""
			var posBreaks : List[Int] = Nil
			var numReads, afwdS, arevS = 0

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
				inputBreak.close
			}

			def delDup(c5p : HashMap[Int,Int], c3p: HashMap[Int,Int]) : String ={
				if (c5p.size >= 1 && c3p.size >= 1){
						val best5p = c5p.toSeq.sortBy(-_._2).toArray.apply(0)._1
						val best3p = c3p.toSeq.sortBy(-_._2).toArray.apply(0)._1
						if (best5p < best3p) {
							"DUP"
						} else {
							if (best5p > best3p) "DEL" else if (best5p == best3p) "EQL" else "UNK"
						}
					} else{
						"UNK"
					}
			}


			def polyAT(cig: String, read: String): Double = {
				if (cig.contains("S") && !cig.map(s => List('I','D','X','=').contains(s)).contains(true)) {
					if (cig.indexOf("S") < cig.indexOf("M")) {
						val idx = cig.split("S").apply(0).toInt
						if (idx >= 10) read.substring(0,idx).count(s => s == 'A' || s == 'T')/read.substring(0,idx).length.toDouble else 0.0
					 } else { 
						val idx =  cig.split("M").apply(0).toInt
						if (idx >= 10) read.substring(idx,read.length).count(s => s == 'A' || s == 'T')/read.substring(idx,read.length).length.toDouble else 0.0
					}
				} else {
					0.0
				}

			}

			import sys.process._
			def getLINEMP(reads: List[String]): List[Int] = {
				if (reads.size >= 1){

					val out = new BufferedWriter(new FileWriter(new File(s"/tmp/${args(0)}_tmp.fa")))
					reads.foreach(s => out.write(s">R${scala.util.Random.nextInt(100000)}\n${s}\n"))
					//out.write(s">read\n${read}")
					out.close
					val result = (s"blat L1_BT_seq.fasta /tmp/${args(0)}_tmp.fa stdout" !!).split("\n")
					if (result.size < 6) {
						List(-1)
					} else {
						result.slice(5,result.size).map(s => s.split("\t").apply(15).toInt).toList
					}
					//if (result.size == 6 ) System.err.println(result(5)) else System.err.println(result.size)
					//if (result.size < 6 ) -1 else result(5).split("\t").apply(15).toInt
					} else {
						List(-1)
					}
			}


			def findOthers(breaks: HashSet[Int]) : Unit = {
				val inputBreak = htsjdk.samtools.SamReaderFactory.make.validationStringency(ValidationStringency.SILENT).open(new File(args(1)))
				var breakpoints = breaks.toArray.sorted
					/* Extract reads covering break point if possible */
					//println("Breakpoint refinement")
				if (breakpoints.size == 2){
					val break = inputBreak.queryOverlapping(chrs._1,breakpoints(0),breakpoints(1))
					while (break.hasNext){
						val tmpBreak = break.next
						if(isGAP(tmpBreak.getMateAlignmentStart,gapsChrom(tmpBreak.getReferenceName))._1) gapBo = true
						if (tmpBreak.getAlignmentStart < breakpoints(0) && tmpBreak.getAlignmentEnd > breakpoints(1)){
							hetBridge += 1
						}
						readDepth += 1
						if (tmpBreak.getMappingQuality >= 20 && tmpBreak.getMappingQuality <= 60 && ( ( tmpBreak.getCigarString.contains("S") && tmpBreak.getStringAttribute("SA") == null) || tmpBreak.getCigarString.contains("H") ) && repeatsChrom.contains(tmpBreak.getReferenceName) && repeatsChrom.contains(tmpBreak.getMateReferenceName)) {
							if (tmpBreak.getProperPairFlag == true){
								if (tmpBreak.getReadNegativeStrandFlag == false) {
									splitEnd += tmpBreak.getAlignmentStart
									fwdS += 1
									if (polyAT(tmpBreak.getCigarString,tmpBreak.getReadString) >= 0.8) polyA = true
									if (c5pSplit.contains(tmpBreak.getAlignmentEnd)) c5pSplit(tmpBreak.getAlignmentEnd) = c5pSplit(tmpBreak.getAlignmentEnd) + 1 else c5pSplit += tmpBreak.getAlignmentEnd -> 1
								} else {
									if (c3pSplit.contains(tmpBreak.getAlignmentStart)) c3pSplit(tmpBreak.getAlignmentStart) = c3pSplit(tmpBreak.getAlignmentStart) + 1 else c3pSplit += tmpBreak.getAlignmentStart -> 1
									splitEnd += tmpBreak.getAlignmentEnd
									revS += 1
									if (polyAT(tmpBreak.getCigarString,tmpBreak.getReadString) >= 0.8) polyA = true
								}
							} else {
								val curMateLTRcheck = isLTR(tmpBreak.getMateAlignmentStart,repeatsChrom(tmpBreak.getMateReferenceName))
								val curMainLTRcheck = isLTR(tmpBreak.getAlignmentStart,repeatsChrom(tmpBreak.getReferenceName))
								if (curMainLTRcheck._1 == false && curMateLTRcheck._1 == true && gapBo == false){							
									if (tmpBreak.getReadNegativeStrandFlag == false) {
										splitEnd += tmpBreak.getAlignmentEnd
										fwdS += 1
										if (polyAT(tmpBreak.getCigarString,tmpBreak.getReadString) >= 0.8) polyA = true
										if (c5pSplit.contains(tmpBreak.getAlignmentEnd)) c5pSplit(tmpBreak.getAlignmentEnd) = c5pSplit(tmpBreak.getAlignmentEnd) + 1 else c5pSplit += tmpBreak.getAlignmentEnd -> 1
									} else {
										if (c3pSplit.contains(tmpBreak.getAlignmentStart)) c3pSplit(tmpBreak.getAlignmentStart) = c3pSplit(tmpBreak.getAlignmentStart) + 1 else c3pSplit += tmpBreak.getAlignmentStart -> 1
										splitEnd += tmpBreak.getAlignmentStart
										revS += 1
										if (polyAT(tmpBreak.getCigarString,tmpBreak.getReadString) >= 0.8) polyA = true
									}
								}
							}
						}
					}
				}
				inputBreak.close
			}


			/* Process all alignments on a Chromosome */
			while (alignments.hasNext){
				val tmp = alignments.next
				if (windowBoo) numReads += 1		

				val best5p : Int = if (c5pSplit.size >= 1) c5pSplit.toSeq.sortBy(-_._2).toArray.apply(0)._1 else -1
				val best3p : Int = if (c3pSplit.size >= 1) c3pSplit.toSeq.sortBy(-_._2).toArray.apply(0)._1 else -1
				/* If the next read is outside of the current region then write out region and advance */
				if (candidateWindowEnd != 0 && tmp.getAlignmentStart >= candidateWindowEnd && ((
					if (best5p == -1 ) true else (isGAP(best5p,gapsChrom(tmp.getReferenceName))._1 == false)) && (
					if (best3p == -1 ) true else ( isGAP(best3p,gapsChrom(tmp.getReferenceName))._1 == false)))) {

					val Pipp = bootStrap(numReads,fwdP + revP, expReads, expIPP)

					posBreaks = posBreaks.sorted.slice((0.2*posBreaks.size).toInt,posBreaks.size - (0.2*posBreaks.size).toInt)

					val avgSPDif = if(posBreaks.size > 2) (posBreaks.sorted drop 1, posBreaks.sorted).zipped.map(_-_).sum/posBreaks.size else 0
					println( s"${chrs._1}:${candidateWindowStart}-${candidateWindowEnd} ${numReads} ${fwdP} ${revP} ${Pipp} ${fwdS} ${revS} ${typeEvent.size} avgDif: ${avgSPDif}")


				if (Pipp <= 0.05 ){ //&& bootStrap(numReads,fwdP + revP, expReads, expIPP) <= 0.005){
					fwdS += afwdS
					revS += arevS
				}

					val elementIs = if (elLength.filter(s => s >= 0).size > 0) elLength.filter(s => s >= 0).sorted.head.toString else "NA"
					val targets = splitEnd.toArray.sorted
					if (splitEnd.size == 2 && fwdP > 0 && revP > 0) findOthers(splitEnd) else if (splitEnd.size == 1 && fwdP > 2 && revP > 2) findOthers(HashSet(targets(0) - 10, targets(0) + 10)) else if (splitEnd.size == 0 && fwdP > 4 && revP > 4) findOthers(HashSet(firstEnd,lastStart))
					val splitDif = if (splitEnd.size == 2) scala.math.abs(splitEnd.toArray.apply(0) - splitEnd.toArray.apply(1)) else 0
					
					println(s"${Pipp <= 0.05} ${(fwdS + revS) > 0} ${typeEvent.size == 1}")

					if (Pipp <= 0.05 && fwdP > 0 && revP > 0 && (fwdS + revS) > 0  && typeEvent.size == 1 && avgSPDif < 50) {
					//if (bootStrap(numReads,fwdP + revP + fwdS + revS, expReads, expIPP + expSR) <= 0.05 && (fwdS + revS) > 0 && mateSR >= 4 && typeEvent.size == 1) {
						//println(s"Criteria for print ${chrs._1}:${candidateWindowStart}-${candidateWindowEnd}\t${fwdP}\t${revP}\t${fwdS}\t${revS}\t${splitEnd}")
						//if (splitEnd.size == 2 ) updateCounts(splitEnd) else if (splitEnd.size == 1) updateCounts(HashSet(targets(0) - 10, targets(0) + 10))
						//println(s"PASSED ${chrs._1}:${candidateWindowStart}-${candidateWindowEnd}\t${fwdP}\t${revP}\t${fwdS}\t${revS}\t${splitEnd}")
						val ornt = s"${RpMp}:${RpMn}:${RnMp}:${RnMn}"
						outEvent.write(s"${chrs._1}:${candidateWindowStart}-${candidateWindowEnd}\t${fwdP}\t${revP}\t${fwdS}\t${revS}\t${splitEnd}\t${splitDif}\t${ornt}\t")
						typeEvent.foreach(s => outEvent.write(s"${s._1},${s._2._1},${s._2._2}:"))
						outEvent.write(s"\t${readDepth}:${hetBridge}\t${windowDepth}\t${allIPP - (fwdP+revP)}\t${(fwdP+revP)/(allIPP.toFloat - (fwdP+revP).toFloat + 0.000001)}\t${delDup(c5pSplit,c3pSplit)}\t${best5p}:${best3p}\t${polyA}\t${elementIs}\n")
						//outFWD.write(s"${chrs._1}\t${candidateWindowStart}\t${candidateWindowEnd}\t${fwdP + fwdS + revP + revS}\n")
					} // end of write

					candidateWindowEnd = 0
					candidateWindowStart = 0
					windowBoo = false
					fwdP = 0
					fwdS = 0
					revS = 0
					revP = 0
					arevS = 0
					afwdS = 0
					allIPP = 0
					splitEnd = new HashSet[Int]
					typeEvent = new HashMap[String,Tuple2[Int,Int]]
					RpMn = 0
					RpMp = 0
					RnMp = 0
					RnMn = 0
					c5pSplit = new HashMap[Int,Int]
					c3pSplit = new HashMap[Int,Int]
					firstEnd = 0
					lastStart = 0
					lastfwd = false				
					readDepth = 0
					hetBridge = 0
					windowDepth = 0
					numReads = 0
					elLength = Nil
					polyA = false
					gapBo = false
					mateReads = Nil
					posBreaks = Nil
				} // end of window work

				/* PRoper window identification, Check if on good Chromosomes */
				if (repeatsChrom.contains(tmp.getReferenceName) && repeatsChrom.contains(tmp.getMateReferenceName)){
					if (windowBoo) windowDepth += 1
					if ( gapStart > 0  && (tmp.getAlignmentStart > (gapStart + 2000) || gapChr !=  tmp.getReferenceName )){
						gapStart = -1000
						gapChr = ""
						gapBo = false
					}

					if (tmp.getProperPairFlag == false && tmp.getMappingQuality >= 20 && tmp.getMappingQuality <= 60 && gapStart < 0){
						
						val curMateLTRcheck : Tuple3[Boolean,String,String] = if (isLTR(tmp.getMateAlignmentStart,repeatsChrom(tmp.getMateReferenceName))._1) isLTR(tmp.getMateAlignmentStart,repeatsChrom(tmp.getMateReferenceName)) else isLTR(tmp.getMateAlignmentStart + 100,repeatsChrom(tmp.getMateReferenceName))

						val curMainLTRcheck : Tuple3[Boolean,String,String] = if (isLTR(tmp.getAlignmentStart,repeatsChrom(tmp.getReferenceName))._1 == false ) isLTR(tmp.getAlignmentEnd,repeatsChrom(tmp.getMateReferenceName)) else isLTR(tmp.getAlignmentEnd,repeatsChrom(tmp.getReferenceName))

						if (curMainLTRcheck._1 == false && curMateLTRcheck._1 == true){
							allIPP += 1
							if(isGAP(tmp.getAlignmentStart,gapsChrom(tmp.getReferenceName))._1 == true) {
								gapBo = true
								gapStart = tmp.getAlignmentStart
								gapChr = tmp.getReferenceName
							}

							if (windowBoo == false && gapBo == false){
								if (tmp.getReadNegativeStrandFlag == false){
									firstEnd = tmp.getAlignmentEnd
									//output.addAlignment(tmp)
									println(" create new window " + tmp.getAlignmentStart)
									windowBoo = true
									windowDepth += 1
									fwdP += 1
									lastfwd = true
									directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
									typeEvent += curMateLTRcheck._3 -> (1,0)
									candidateWindowStart = ((tmp.getAlignmentStart / 1000)*1000)
									// 10300 = /1k*1k = 10,000 - 10300 = -300 vs 200 convert to absolute
									if (abs(candidateWindowStart - tmp.getAlignmentStart) > 500 ) candidateWindowStart += 500
									candidateWindowEnd = candidateWindowStart + 1500
									
									//println(candidateWindowStart + " " + candidateWindowEnd)
								} //end if fwd read

							} else {
								/* Window is True, improperly paired and good mapping quality */
								if (gapBo == false) {
									var cMate = tmp 
										try {
											cMate = mates.queryMate(tmp)
											} catch {
												case e: Exception => 
											}
									val useMate = if (cMate == null) false else ! cMate.getCigarString.map(s => List('I','D','X','=').contains(s)).contains(true)
									//if (! tmp.getCigarString.contains("S")) { //|| tmp.getSupplementaryAlignmentFlag == false){ Not supllemantary alignment
									//if (tmp.getStringAttribute("SA") == null || !tmp.getCigarString.contains("S") || !tmp.getCigarString.contains("H")) { 
									//println(tmp.getCigarString + " " + tmp.getCigarString.contains("S") + " " + tmp.getCigarString.contains("H"))

									if (tmp.getCigarString.contains("S") || tmp.getCigarString.contains("H")){
										if (tmp.getCigarString.contains("S")){
											posBreaks = if(tmp.getCigarString.indexOf("S") > tmp.getCigarString.indexOf("M")) tmp.getAlignmentEnd :: posBreaks else tmp.getAlignmentStart :: posBreaks
										} else {
											posBreaks = if(tmp.getCigarString.indexOf("H") > tmp.getCigarString.indexOf("M")) tmp.getAlignmentEnd :: posBreaks else tmp.getAlignmentStart :: posBreaks
										}

									/*
									* Add logic that once we start identifing 3' reads we should not count 5' reads, 
									* if we see 5' reads we should close the window and reopen a new one.
									*/
									println("Split Read")
										//println(" is getSupplementaryAlignmentFlag")
										if (tmp.getReadNegativeStrandFlag == false){
												if (typeEvent.contains(curMateLTRcheck._3)) typeEvent(curMateLTRcheck._3) = (typeEvent(curMateLTRcheck._3)._1 + 1, typeEvent(curMateLTRcheck._3)._2) else typeEvent += curMateLTRcheck._3 -> (1,0)
												fwdS += 1
												fwdP += 1
												//if (useMate ) elLength = getLINEMP(cMate.getReadString) :: elLength
												if (useMate) {
													mateReads = cMate.getReadString :: mateReads
													//if (splitMate(cMate.getCigarString)) mateSR += 1
												}

												if (c5pSplit.contains(tmp.getAlignmentEnd)) c5pSplit(tmp.getAlignmentEnd) = c5pSplit(tmp.getAlignmentEnd) + 1 else c5pSplit += tmp.getAlignmentEnd -> 1
												directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
												splitEnd += tmp.getAlignmentEnd
											} else {
												if (typeEvent.contains(curMateLTRcheck._3))typeEvent(curMateLTRcheck._3) =  (typeEvent(curMateLTRcheck._3)._1, typeEvent(curMateLTRcheck._3)._2 + 1) else typeEvent += curMateLTRcheck._3 -> (0,1)
												revS += 1
												revP += 1
												//if (useMate ) elLength = getLINEMP(cMate.getReadString) :: elLength
												if (useMate) {
													mateReads = cMate.getReadString :: mateReads
													//if (splitMate(cMate.getCigarString)) mateSR += 1
												}
												if (c3pSplit.contains(tmp.getAlignmentStart)) c3pSplit(tmp.getAlignmentStart) = c3pSplit(tmp.getAlignmentStart) + 1 else c3pSplit += tmp.getAlignmentStart -> 1
												directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
												splitEnd += tmp.getAlignmentStart
											} // end sup is neg read

										} else {

											if (tmp.getReadNegativeStrandFlag == false){ // Not sup and on +ve strand
												if (lastfwd) firstEnd = tmp.getAlignmentEnd
												if (typeEvent.contains(curMateLTRcheck._3)) typeEvent(curMateLTRcheck._3) = (typeEvent(curMateLTRcheck._3)._1 + 1, typeEvent(curMateLTRcheck._3)._2) else typeEvent += curMateLTRcheck._3 -> (1,0)
												directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
												//list of positions within the element the mate maps to
												//if (useMate ) elLength = getLINEMP(cMate.getReadString) :: elLength
													if (useMate) {
														mateReads = cMate.getReadString :: mateReads
														//if (splitMate(cMate.getCigarString)) mateSR += 1
													}
												fwdP += 1
											} else {
												if (lastfwd) {
													lastStart = tmp.getAlignmentStart
													lastfwd = false
												}
												if (typeEvent.contains(curMateLTRcheck._3)) typeEvent(curMateLTRcheck._3) = (typeEvent(curMateLTRcheck._3)._1, typeEvent(curMateLTRcheck._3)._2 + 1) else typeEvent += curMateLTRcheck._3 -> (0,1)
													directionLTR(curMateLTRcheck._2,tmp.getReadNegativeStrandFlag,tmp.getMateNegativeStrandFlag)
													revP += 1
													//if (useMate ) elLength = getLINEMP(cMate.getReadString) :: elLength
													if (useMate) {
														mateReads = cMate.getReadString :: mateReads
														//if (splitMate(cMate.getCigarString)) mateSR += 1
													}
											}//end fwd read npp

										}

								}//is gap?

							} // end if window true
						} // is Repeat structure of mates fine?

					} else {
						if (windowBoo && tmp.getMappingQuality == 60 && (tmp.getCigarString.contains("S") || tmp.getCigarString.contains("H"))){
							if (tmp.getReadNegativeStrandFlag) afwdS += 1 else arevS += 1
							if (tmp.getCigarString.contains("S")){
									posBreaks = if(tmp.getCigarString.indexOf("S") > tmp.getCigarString.indexOf("M")) tmp.getAlignmentEnd :: posBreaks else tmp.getAlignmentStart :: posBreaks
								} else {
									posBreaks = if(tmp.getCigarString.indexOf("H") > tmp.getCigarString.indexOf("M")) tmp.getAlignmentEnd :: posBreaks else tmp.getAlignmentStart :: posBreaks
								}
							
						}
					} // end good read mapping and improperly paired

					//if (windowBoo == true) output.addAlignment(tmp)
				}// End of Good Chromosomes check
				//}//top level isGap
			} //End of While Alignments
			/* Have completed a chromosome scan now need to clean up*/
			val best5p = if (c5pSplit.size >= 1) c5pSplit.toSeq.sortBy(-_._2).toArray.apply(0)._1 else -1
			val best3p = if (c3pSplit.size >= 1) c3pSplit.toSeq.sortBy(-_._2).toArray.apply(0)._1 else -1
			//if (windowBoo) {
			if (windowBoo && ((if (best5p != -1 ) isGAP(best5p,gapsChrom(chrs._1))._1 == false else true) && (if (best3p != -1 ) isGAP(best3p,gapsChrom(chrs._1))._1 == false else true))) {
				//println(" no gap: " + isGAP(best5p,gapsChrom(chrs._1)) + " " + isGAP(best3p,gapsChrom(chrs._1)))
				//val elementIs = if (elLength.filter(s => s > 1).sum/elLength.filter(s => s > 1).size.toDouble > 4000) "Short" else if (elLength.filter(s => s > 1).size >= 1) "Long" else "NA"
 					
					val Pipp = bootStrap(numReads,fwdP + revP, expReads, expIPP)

					posBreaks = posBreaks.sorted.slice((0.2*posBreaks.size).toInt,posBreaks.size - (0.2*posBreaks.size).toInt)

					//val avgSPDif = if(posBreaks.size > 2) (posBreaks.sorted drop 1, posBreaks.sorted).zipped.map(_-_).apply(posBreaks.size/2) else 0
					val avgSPDif = if(posBreaks.size > 2) (posBreaks.sorted drop 1, posBreaks.sorted).zipped.map(_-_).sum/posBreaks.size else 0

					println( s"${chrs._1}:${candidateWindowStart}-${candidateWindowEnd} ${numReads} ${fwdP} ${revP} ${Pipp} ${fwdS} ${revS} ${typeEvent.size} avgDif: ${avgSPDif}")

				if (Pipp <= 0.05 ){ //&& bootStrap(numReads,fwdP + revP, expReads, expIPP) <= 0.005){
					fwdS += afwdS
					revS += arevS
				}

				//elLength = getLINEMP(mateReads)
				val elementIs = if (elLength.filter(s => s >= 0).size > 0) elLength.filter(s => s >= 0).sorted.head.toString else "NA"
				val targets = splitEnd.toArray.sorted
				//if (splitEnd.size == 2 && fwdP > 0 && revP > 0) findOthers(splitEnd) else if (splitEnd.size == 1 && fwdP > 0 && revP > 0) findOthers(HashSet(targets(0) - 10, targets(0) + 10)) else if (splitEnd.size == 0 && fwdP > 4 && revP > 4) findOthers(HashSet(firstEnd,lastStart))
				val splitDif = if (splitEnd.size == 2) scala.math.abs(splitEnd.toArray.apply(0) - splitEnd.toArray.apply(1))+1 else 0
				//if (fwdP >= 3 && revP >= 3 && (fwdS + revS) > 0 && mateSR >= 4) {
				//if (bootStrap(numReads,fwdP + revP + fwdS + revS + mateSR, expReads, expIPP + expSR) <= 0.05 && (fwdS + revS) > 0 && mateSR >= 4 && typeEvent.size == 1) {
				if (Pipp <= 0.05 && fwdP > 0 && revP > 0 && (fwdS + revS) > 0 && avgSPDif < 50) {
						val ornt = s"${RpMp}:${RpMn}:${RnMp}:${RnMn}"
						outEvent.write(s"${chrs._1}:${candidateWindowStart}-${candidateWindowEnd}\t${fwdP}\t${revP}\t${fwdS}\t${revS}\t${splitEnd}\t${splitDif}\t${ornt}\t")
						typeEvent.foreach(s => outEvent.write(s"${s._1},${s._2._1},${s._2._2}:"))
						outEvent.write(s"\t${readDepth}:${hetBridge}\t${windowDepth}\t${allIPP - (fwdP+revP)}\t${(fwdP+revP)/(allIPP.toFloat - (fwdP+revP).toFloat + 0.000001)}\t${delDup(c5pSplit,c3pSplit)}\t${best5p}:${best3p}\t${polyA}\t${elementIs}\n")
					}
			}
			input.close
			mates.close
		} // end for chr from Chroms
		outEvent.close
		//outFWD.close
		//output.close
	} // def main

} //object