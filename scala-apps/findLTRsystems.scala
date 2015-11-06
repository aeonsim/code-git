import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.SAMFileWriter
import htsjdk.samtools.SAMFileWriterFactory
import htsjdk.samtools.DefaultSAMRecordFactory
import htsjdk.samtools.SamReader
import htsjdk.samtools.SAMRecord
import java.io.File
import java.io._
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
var alignments = input.queryAlignmentStart(Chroms(0)._1,1)
val outFactory = new htsjdk.samtools.SAMFileWriterFactory
val output = outFactory.makeSAMOrBAMWriter(input.getFileHeader,true,new File(s"${args(0)}.notPP.sec.bam"))

val repeats = new BufferedReader(new FileReader(new File(args(2))))

val outFWD = new BufferedWriter(new FileWriter(new File(s"${args(0)}.notPP.bedgraph")))
//val outREV = new BufferedWriter(new FileWriter(new File(s"${args(0)}.notPP.rev.bedgraph")))
//val outRSPL = new BufferedWriter(new FileWriter(new File(s"${args(0)}.split.rev.bedgraph")))
val outFSPL = new BufferedWriter(new FileWriter(new File(s"${args(0)}.split.bedgraph")))

outFWD.write("track type=bedGraph\n")
//outREV.write("track type=bedGraph\n")
//outRSPL.write("track type=bedGraph\n")
outFSPL.write("track type=bedGraph\n")

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
	if (curSearch.length == 0) {(false,"","")} else {
		if (curSearch(half)._1 <= pos && curSearch(half)._2 >= pos) {
			(true,curSearch(half)._6,s"${curSearch(half)._3}-${curSearch(half)._5}")
		} else {
			if (curSearch.length <= 1) (false,"","") else {
				if (curSearch(half)._1 > pos) isLTR(pos,curSearch.slice(0,half)) else isLTR(pos,curSearch.slice(half + 1,curSearch.length))
			}
		}
	}
}

//chr1    387557  387760
//isLTR(387760,repeatsChrom("chr1"))


for (chrs <- Chroms){
	println(chrs._1)
	var end = 1000
	var start = 1
	var fwdCount, revCount, splRCount, splFCount, isLTRCount, allCount = 0
	var LTRpRn, LTRpRp, LTRnRp, LTRnRn = 0
	val input = htsjdk.samtools.SamReaderFactory.make.open(new File(args(1)))
	var alignments = input.query(chrs._1,1,200000000,true)
	var typeEvent = new HashSet[String]

	var fwdP, revP, fwdS, revS = 0
	var splitEnd : List[Int] = Nil
	var candidateWindowStart, candidateWindowEnd = 0
	var windowBoo = false

	while (alignments.hasNext){
		val tmp = alignments.next
		
		/* If the next read is outside of the current region then write out region and advance */
		//if (tmp.getAlignmentStart >= end && ( fwdCount >= 1 || revCount >= 1)){
		
		if (tmp.getAlignmentStart >= end){
			val curEvent = if (typeEvent.size >= 1) typeEvent.reduceLeft(_ + "," + _) else ""
			outFWD.write(s"${chrs._1}\t${start}\t${end}\t${fwdCount+revCount}\t${LTRpRp}:${LTRpRn}:${LTRnRp}:${LTRnRn}\t${fwdCount}\t${revCount}\t${curEvent}\n")
			outFSPL.write(s"${chrs._1}\t${start}\t${end}\t${splFCount+splRCount}\t${LTRpRp}:${LTRpRn}:${LTRnRp}:${LTRnRn}\t${splFCount}\t${splRCount}\t${curEvent}\n")
			//outREV.write(s"${chrs._1}\t${start}\t${end}\t${revCount}\t${LTRpRp}:${LTRpRn}:${LTRnRp}:${LTRnRn}\t${curEvent}\n")		
			//outRSPL.write(s"${chrs._1}\t${start}\t${end}\t${splRCount}\t${LTRpRp}:${LTRpRn}:${LTRnRp}:${LTRnRn}\t${curEvent}\n")

			
			start = start + 500
			end += 500
			typeEvent = new HashSet[String]
			fwdCount = 0
			revCount = 0
			splRCount = 0
			splFCount = 0
			allCount =  0
			LTRpRn = 0
			LTRpRp = 0
			LTRnRp = 0
			LTRnRn = 0
		}

		/* PRoper window identification */
		if (repeatsChrom.contains(tmp.getReferenceName) && repeatsChrom.contains(tmp.getMateReferenceName)){

		if (windowBoo == false && tmp.getProperPairFlag == false && tmp.getMappingQuality == 60 && tmp.getReadNegativeStrandFlag == false){
			val curLTRcheck = isLTR(tmp.getMateAlignmentStart,repeatsChrom(tmp.getMateReferenceName))
			val curMainLTRcheck = isLTR(tmp.getAlignmentStart,repeatsChrom(tmp.getReferenceName))
				if (curMainLTRcheck._1 == false && curLTRcheck._1 == true){
					windowBoo = true
					fwdP += 1
					candidateWindowStart = ((tmp.getAlignmentStart / 1000)*1000)
					if ((candidateWindowStart - tmp.getAlignmentStart) > (candidateWindowStart + 500 - tmp.getAlignmentStart)) candidateWindowStart += 500
					candidateWindowEnd = candidateWindowStart + 1000	
				}
		}

		if (windowBoo == true && tmp.getProperPairFlag == false && tmp.getMappingQuality == 60){
			val curLTRcheck = isLTR(tmp.getMateAlignmentStart,repeatsChrom(tmp.getMateReferenceName))
			val curMainLTRcheck = isLTR(tmp.getAlignmentStart,repeatsChrom(tmp.getReferenceName))
			if (curMainLTRcheck._1 == false && curLTRcheck._1 == true){
				if (tmp.getStringAttribute("SA") == null){
					if (tmp.getReadNegativeStrandFlag == false) fwdP += 1 else revP += 1 
				} else {
					if (tmp.getReadNegativeStrandFlag == false){
						splitEnd = tmp.getAlignmentEnd :: splitEnd
						fwdS += 1
					} else{
						splitEnd = tmp.getAlignmentStart :: splitEnd
						revS += 1
					}
				}

			}
		}

	}
		if (tmp.getAlignmentStart >= candidateWindowEnd){
			if (fwdP > 0 && fwdS > 0 && revP > 0 && revS > 0) println(s"${chrs._1}:${candidateWindowStart}-${candidateWindowEnd}\t${fwdP}\t${revP}\t${fwdS}\t${revS}\t${splitEnd}")
			candidateWindowEnd = 0
			candidateWindowStart = 0
			windowBoo = false
			fwdP = 0
			fwdS = 0
			revS = 0
			revP = 0
			splitEnd = Nil
		}


		/* IF your an improper pair or a secondary/split read and your mates chromosome is in the repeats database*/
		
		if (tmp.getProperPairFlag == false && tmp.getStringAttribute("SA") != null){
			if (tmp.getReadNegativeStrandFlag) splRCount +=1 else splFCount += 1
		}
		
		if ((tmp.getProperPairFlag == false || tmp.getStringAttribute("SA") != null || tmp.getSupplementaryAlignmentFlag == true) && repeatsChrom.contains(tmp.getMateReferenceName)){
			
			/*	Check to see if your mate is present in a repeat if it isn't do nothing/loop if it is write it out to the SAM*/
			val curLTRcheck = isLTR(tmp.getMateAlignmentStart,repeatsChrom(tmp.getMateReferenceName))
			val curMainLTRcheck = isLTR(tmp.getAlignmentStart,repeatsChrom(tmp.getReferenceName))
			
			if (curLTRcheck._1 && tmp.getMateUnmappedFlag == false && tmp.getMappingQuality >= 60 && (tmp.getInferredInsertSize == 0 || scala.math.abs(tmp.getInferredInsertSize) >= 1000000) && (if (args(3) == "true") {if (curMainLTRcheck._1 && curLTRcheck._1) false else true} else true)){
				allCount += 1
				output.addAlignment(tmp)
				if (curLTRcheck._1) typeEvent += curLTRcheck._3 else typeEvent += curMainLTRcheck._3
				
				/* If -ve stand record in rev if your improper paired or Split else record in the fwd */
				
				if (tmp.getReadNegativeStrandFlag == true ){
					if (tmp.getProperPairFlag == false){ 
						revCount += 1
						if (curLTRcheck._2 == "+"){
							if (tmp.getMateNegativeStrandFlag == true) LTRnRn += 1 else LTRnRp += 1
						}
					}
					//if (tmp.getSupplementaryAlignmentFlag == true || tmp.getStringAttribute("SA") != null) splRCount += 1
					
				} else {
					if (tmp.getProperPairFlag == false) {
						fwdCount += 1
						if (curLTRcheck._2 == "+"){
							if (tmp.getMateNegativeStrandFlag == true) LTRpRn += 1 else LTRpRp += 1
						}

					}
					//if (tmp.getSupplementaryAlignmentFlag == true || tmp.getStringAttribute("SA") != null) splFCount += 1
				}
			}

		}
		

	}
	
	/* Once outside the the Chr prepare things for the next Chr, and write the ends, close the current input*/
	val curEvent = if (typeEvent.size >= 1) typeEvent.reduceLeft(_ + "," + _) else ""
	outFWD.write(s"${chrs._1}\t${start}\t${end}\t${fwdCount+revCount}\t${LTRpRp}:${LTRpRn}:${LTRnRp}:${LTRnRn}\t${fwdCount}\t${revCount}\t${curEvent}\n")
	outFSPL.write(s"${chrs._1}\t${start}\t${end}\t${splFCount+splRCount}\t${LTRpRp}:${LTRpRn}:${LTRnRp}:${LTRnRn}\t${splFCount}\t${splRCount}\t${curEvent}\n")
	//outFWD.write(s"${chrs._1}\t${start}\t${end}\t${fwdCount}\t${LTRpRp}:${LTRpRn}:${LTRnRp}:${LTRnRn}\t${curEvent}\n")
	//outREV.write(s"${chrs._1}\t${start}\t${end}\t${revCount}\t${LTRpRp}:${LTRpRn}:${LTRnRp}:${LTRnRn}\t${curEvent}\n")
	//outFSPL.write(s"${chrs._1}\t${start}\t${end}\t${splFCount}\t${LTRpRp}:${LTRpRn}:${LTRnRp}:${LTRnRn}\t${curEvent}\n")
	//outRSPL.write(s"${chrs._1}\t${start}\t${end}\t${splRCount}\t${LTRpRp}:${LTRpRn}:${LTRnRp}:${LTRnRn}\t${curEvent}\n")

	fwdCount = 0
	revCount = 0
	splRCount = 0
	splFCount = 0
	typeEvent = new HashSet[String]
	LTRpRn = 0
	LTRpRp = 0
	LTRnRp = 0
	LTRnRn = 0
	input.close
}

output.close
outFWD.close
//outREV.close
outFSPL.close
//outRSPL.close
}

}