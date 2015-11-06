import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.SAMFileWriter
import htsjdk.samtools.SAMFileWriterFactory
import htsjdk.samtools.DefaultSAMRecordFactory
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

	val input = htsjdk.samtools.SamReaderFactory.make.open(new File(args(1)))
	//var alignments = input.queryAlignmentStart(Chroms(0)._1,1)

	val inputBreak = htsjdk.samtools.SamReaderFactory.make.open(new File(args(1)))
	//var alignmentsBreak = inputBreak.queryAlignmentStart(Chroms(0)._1,1)

	val outFactory = new htsjdk.samtools.SAMFileWriterFactory
	val output = outFactory.makeSAMOrBAMWriter(input.getFileHeader,true,new File(s"${args(0)}.notPP.sec.bam"))

	val repeats = new BufferedReader(new FileReader(new File(args(2))))

	val outFWD = new BufferedWriter(new FileWriter(new File(s"${args(0)}.bedgraph")))
	val outEvent = new BufferedWriter(new FileWriter(new File(s"${args(0)}.event.tab")))

	outFWD.write("track type=bedGraph\n")

	input.close

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
	}

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
	

//chr1    387557  387760
//isLTR(387760,repeatsChrom("chr1"))


	for (chrs <- Chroms){
		println(chrs._1)
		var LTRpRn, LTRpRp, LTRnRp, LTRnRn = 0
		val input = htsjdk.samtools.SamReaderFactory.make.open(new File(args(1)))
		var alignments = input.query(chrs._1,1,200000000,true)
		var typeEvent = new HashSet[String]

		var fwdP, revP, fwdS, revS = 0
		var splitEnd = new HashSet[Int]
		var candidateWindowStart, candidateWindowEnd = 0
		var windowBoo = false

		def directionLTR (LTR: String, mateNeg: Boolean) : Unit = {
			if (LTR == "+"){
					if (mateNeg) LTRpRn += 1 else LTRpRp += 1 
				} else {
					if (mateNeg) LTRnRn += 1 else LTRnRp += 1
				}
		}

		def updateCounts(tmp: SAMRecord) : Unit = {
			if (tmp.getAlignmentStart >= candidateWindowEnd){
				if (splitEnd.size == 2){
					val inputBreak = htsjdk.samtools.SamReaderFactory.make.open(new File(args(1)))
					var breakpoints = splitEnd.toArray.sorted
					/* Extract reads covering break point if possible */
					val break = inputBreak.queryOverlapping(chrs._1,breakpoints(0),breakpoints(1))
					while (break.hasNext){
						val tmpBreak = break.next
						if ((breakpoints.contains(tmpBreak.getAlignmentStart) || breakpoints.contains(tmpBreak.getAlignmentEnd)) && tmpBreak.getProperPairFlag == true && tmpBreak.getMappingQuality == 60 && tmpBreak.getCigarString.contains("S")) {
							output.addAlignment(tmpBreak)
							if (tmpBreak.getReadNegativeStrandFlag == false) fwdS += 1 else revS += 1
						}
					}
				}
			}
		}

		while (alignments.hasNext){
			val tmp = alignments.next
			
			/* If the next read is outside of the current region then write out region and advance */
			//if (tmp.getAlignmentStart >= end && ( fwdCount >= 1 || revCount >= 1)){
			
			if (tmp.getAlignmentStart >= candidateWindowEnd){
				val splitDif = if (splitEnd.size == 2) scala.math.abs(splitEnd.toArray.apply(0) - splitEnd.toArray.apply(1)) else 0
				if (fwdP > 0 && fwdS > 0 && (revP + revS) > 0) {
					updateCounts(tmp)
					val ornt = s"${LTRpRp}:${LTRpRn}:${LTRnRp}:${LTRnRn}"
					outEvent.write(s"${chrs._1}:${candidateWindowStart}-${candidateWindowEnd}\t${fwdP}\t${revP}\t${fwdS}\t${revS}\t${splitEnd}\t${splitDif}\t${ornt}\n")
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
				typeEvent = new HashSet[String]
				LTRpRn = 0
				LTRpRp = 0
				LTRnRp = 0
				LTRnRn = 0
			} // end of window work

			/* PRoper window identification */
			if (repeatsChrom.contains(tmp.getReferenceName) && repeatsChrom.contains(tmp.getMateReferenceName)){

				if (windowBoo == false && tmp.getProperPairFlag == false && tmp.getMappingQuality == 60 && tmp.getReadNegativeStrandFlag == false){
					val curLTRcheck = isLTR(tmp.getMateAlignmentStart,repeatsChrom(tmp.getMateReferenceName))
					val curMainLTRcheck = isLTR(tmp.getAlignmentStart,repeatsChrom(tmp.getReferenceName))
						if (curMainLTRcheck._1 == false && curLTRcheck._1 == true){
							output.addAlignment(tmp)
							windowBoo = true
							fwdP += 1
							candidateWindowStart = ((tmp.getAlignmentStart / 1000)*1000)
							// 10300 = /1k*1k = 10,000 - 10300 = -300 vs 200 convert to absolute
							if (abs(candidateWindowStart - tmp.getAlignmentStart) >= abs((candidateWindowStart + 500) - tmp.getAlignmentStart)) candidateWindowStart += 500
							candidateWindowEnd = candidateWindowStart + 1500
							//println(candidateWindowStart + " " + candidateWindowEnd)	
						}
				}// end of window creation

				if (windowBoo == true && tmp.getProperPairFlag == false && tmp.getMappingQuality == 60){
				val curLTRcheck = isLTR(tmp.getMateAlignmentStart,repeatsChrom(tmp.getMateReferenceName))
				val curMainLTRcheck = isLTR(tmp.getAlignmentStart,repeatsChrom(tmp.getReferenceName))
					if (curMainLTRcheck._1 == false && curLTRcheck._1 == true){ // Is not in LTR mate is
						if (tmp.getStringAttribute("SA") == null ) { //|| tmp.getSupplementaryAlignmentFlag == false){ Not supllemantary alignment
							if (tmp.getReadNegativeStrandFlag == false){ // Not sup and on +ve strand
								directionLTR(curLTRcheck._2,tmp.getMateNegativeStrandFlag)
								fwdP += 1
							} else {
								directionLTR(curLTRcheck._2,tmp.getMateNegativeStrandFlag)
								revP += 1
							}//end fwd read npp
						} else { // IS Sup Alignment/splitRead
							if (tmp.getReadNegativeStrandFlag == false){
									fwdS += 1
									splitEnd += tmp.getAlignmentEnd
								} else {
									revS += 1
									splitEnd += tmp.getAlignmentStart
								} // end sup is neg read
						} // end Sup Else
					}//is not LTR mate is LTR
				}//appropriate read

				if (windowBoo == true) output.addAlignment(tmp)

			} // end of good chromosomes
		} // end of while for chromosome
		typeEvent = new HashSet[String]
		LTRpRn = 0
		LTRpRp = 0
		LTRnRp = 0
		LTRnRn = 0
		input.close

		} // end of loop for a chrom


		/* Once outside the the Chr prepare things for the next Chr, and write the ends, close the current input*/
}// for loop chroms


} // main

//} //object