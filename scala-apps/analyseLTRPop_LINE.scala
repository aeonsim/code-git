import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet
import org.apache.commons.io.FileUtils._
import java.io._

/* to filter with Awk | awk '{n=split($0,col,"\t"); y=split(col[8],p,":"); if ((p[1]+p[2]+p[3]+p[4]) >= 5 && col[5]/col[7] >= 2 && col[7] >= 2 && col[5] >= 10 && col[4] < 0.2) {print $0};}'*/

val fwd = org.apache.commons.io.FileUtils.listFiles(new File("."),Array("event.tab"),true).iterator
//val sFwd = org.apache.commons.io.FileUtils.listFiles(new File("."),Array("split.fwd.bedgraph"),true).iterator

val depths = new BufferedReader(new FileReader(new File("/home/u/u218757/kos_home/refs/currentDepths.txt")))

/* Chr -> Pos -> Tuple5[Obs,Carriers,ImpropPaired,NumOFBreakpoints(Array[Int]),# SplitReads,BreakDif HashSet[Int], Orination Array, Event HashMap, Breakpoints HashSet, genes hashSet]*/

val data = new HashMap[String,HashMap[Int,Tuple12[Int,List[String],Int,List[Int],Int,HashSet[Int],Array[Int],HashMap[String,Tuple2[Int,Int]],HashMap[Int,Int],HashSet[String],HashMap[String,Array[Int]],HashMap[Int,Int]]]]
//val geneData = new HashMap[String,HashMap[Int,]]

//HashMap[String -> HaspMap[String,Array[Int]]] Pos -> ID -> Reads
var genotypes = new HashMap[String,HashMap[String,Array[Int]]] 

val curDepth = new HashMap[String,Double]
while (depths.ready){
val cline = depths.readLine.split("\t")
curDepth += cline(0) -> cline(1).toDouble
}
depths.close

/* HashMap[Proband -> (Sire, Dam, List[Children])]*/

val families = new HashMap[String,Tuple3[String,String,List[String]]]

val popIn = new BufferedReader(new FileReader(new File("/home/u/u218757/kos_home/refs/Damona-full.ped")))
val annoIn = new BufferedReader(new FileReader(new File("/home/u/u218757/kos_home/refs/Bos_taurus.UMD3.1.82.gtf")))
val annoIn2 = new BufferedReader(new FileReader(new File("/home/u/u218757/kos_home/refs/RefSeq_Genes.gtf")))

//Chr -> Tuple5[start, end, class, name/id, strand]

var anno = new HashMap[String,Array[Tuple5[Int,Int,String,String,String]]]
var annoList : List[Tuple5[Int,Int,String,String,String]] = Nil
		var prevChrom = ""
		
		while (annoIn.ready){
			val tmp = annoIn.readLine.split("\t")
			if (tmp(0).apply(0) != '#'){
				val chr = "chr" + tmp(0).trim
				if (tmp(2) == "transcript") {
					val info = tmp(8).split(" ")
					val name = if (info.contains("gene_name")) info(info.indexOf("gene_name") + 1 ) else info(info.indexOf("gene_id") + 1 )
					if (prevChrom == chr){
						annoList = (tmp(3).toInt,tmp(4).toInt, tmp(2), name, tmp(6)) :: annoList
					} else {
						anno += prevChrom -> annoList.reverse.toArray
						annoList = (tmp(3).toInt,tmp(4).toInt, tmp(2), name, tmp(6)) :: Nil
						prevChrom = chr
					}
				}
			}
		}//end while repeats
		
		anno += prevChrom -> annoList.reverse.toArray
		annoList = Nil
		prevChrom = ""
		annoIn.close
		
		
		
		var anno2 = new HashMap[String,Array[Tuple5[Int,Int,String,String,String]]]
		while (annoIn2.ready){
			val tmp = annoIn2.readLine.split("\t")
			if (tmp(0).apply(0) != '#'){
				val chr = tmp(0)
				if (tmp(2) == "gene") {
					val info = tmp(8).replaceAll("=",";").split(";")
					val name = if (info.contains("symbol_ncbi")) info(info.indexOf("symbol_ncbi") + 1) else if (info.contains("Dbxref")) info(info.indexOf("Dbxref") + 1) else info(info.indexOf("ID") + 1)
					if (prevChrom == chr){
						annoList = (tmp(3).toInt,tmp(4).toInt, tmp(2), name, tmp(6)) :: annoList
					} else {
						anno2 += prevChrom -> annoList.reverse.toArray
						annoList = (tmp(3).toInt,tmp(4).toInt, tmp(2), name, tmp(6)) :: Nil
						prevChrom = chr
					}
				}
			}
		}//end while repeats
		
		anno2 += prevChrom -> annoList.reverse.toArray
		annoList = Nil
		prevChrom = ""
		annoIn.close
		

//boolean, strand, Name?
		def isEle (pos: Int, curSearch: Array[Tuple5[Int,Int,String,String,String]]) : Tuple3[Boolean,String,String] = {
			val half = scala.math.floor(curSearch.length/2).toInt
			if (curSearch.length == 0){
					(false,"","")
			} else {
				if (curSearch(half)._1 <= pos && curSearch(half)._2 >= pos) {
					(true,curSearch(half)._5,s"${curSearch(half)._3}-${curSearch(half)._4}")
					} else {
						if (curSearch.length <= 1) {
							(false,"","")
						} else {
						if (curSearch(half)._1 > pos) isEle(pos,curSearch.slice(0,half)) else isEle(pos,curSearch.slice(half + 1,curSearch.length))
					}
				}
			}
		} //end isEle


val pop = new HashMap[String,Array[String]]
var workingPop = new HashSet[String]

while (popIn.ready){
val current = popIn.readLine.split("\t")
pop += current(1) -> current
}

var populationSize = 0

/* Build pedigrees for Analysis */

	def findChildren (pop: List[Array[String]], parent: String) : List[String]={
		var result: List[String] = Nil
		for (item <- pop){
			if ( (item(2) == parent)||(item(3) == parent) ){
				result = item(1) :: result
			}//eif

		}//efor
		result
	}//edef

/* Recursive function to find all parents & ancestors of nominated individual
* returns list of those present in both the PED & VCF files.
*/


	def itPed (pop: HashMap[String, Array[String]] , rec: String, maxIT: Int): List[String] = {
		if (pop.contains(rec) && maxIT >= 0){
			return rec :: itPed(pop,pop(rec).apply(2),maxIT -1) ::: itPed(pop,pop(rec).apply(3),maxIT -1)
		} else {
			return Nil
		}
	}

var hcTrios = 0

for (indv <- pop){
//if (pop.contains(indv._2(2)) && pop.contains(indv._2(3))){
 val kids = findChildren(pop.values.toList,indv._1)
 //if (kids.size >= 1) families += indv._1 -> (indv._2(2),indv._2(3),kids)
 if (curDepth.contains(indv._2(2)) && curDepth(indv._2(2)) >= 18.0 && curDepth.contains(indv._2(3)) && curDepth(indv._2(3)) >= 18.0) hcTrios += 1 
 families += indv._1 -> (indv._2(2),indv._2(3),kids)
//}
}

println("High Coverage Trios " + hcTrios)

def convert2Pos(input: String): Array[Array[Int]] = {
	//if (input != "Set()") input.substring(4,input.size -1).split(",").map(s => s.trim.toInt) else Array()
	if (input != "Map()") input.slice(4,input.size -1).split(",").map(y => y.split("->").map(_.trim.toInt)) else Array()
}


def addBreaks (breaks: HashSet[Int], newBreaks: String ): HashSet[Int] ={
	val novoBreaks = convert2Pos(newBreaks) // newBreaks.substring(4,newBreaks.size -1).split(",").map(s => s.trim.toInt)
	if (novoBreaks.size <= 2){
		for ( k <- novoBreaks){
			breaks += k(0)
		}
	}
 
	breaks
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

def addBreaksInfo (breaks: HashMap[Int, Int], newBreaks: String ): HashMap[Int,Int] ={

	val novoBreaks = convert2Pos(newBreaks)// newBreaks.substring(4,newBreaks.size -1).split(",").map(s => s.trim.toInt)
	for (split <- novoBreaks){
		if (breaks.contains(split(0))){
			breaks(split(0)) = breaks(split(0)) + split(1)
		} else {
			breaks += split(0) -> split(1)		
		}
	}

	breaks
}


def updateCounts(counts: HashMap[Int, Int], newItem: String ): HashMap[Int,Int] ={

	val novoBreaks = newItem.split(":").map(_.toInt)// newBreaks.substring(4,newBreaks.size -1).split(",").map(s => s.trim.toInt)
	val items = Array(novoBreaks(0),-1 * novoBreaks(1))
	for (split <- items){
		if (counts.contains(split)){
			counts(split) = counts(split) + 1
		} else {
			counts += split -> 1 
		}
	}

	counts
}


//def isEle (pos: Int, curSearch: Array[Tuple5[Int,Int,String,String,String]]) : Tuple3[Boolean,String,String]
//addgenes(curGenesHash,breakpointsString,chrom)
def addGenes (genes: HashSet[String], newPos: String, chr: String): HashSet[String] ={
	val novoGenes = convert2Pos(newPos) // newPos.substring(4,newPos.size -1).split(",").map(s => s.trim.toInt)
	for ( k <- novoGenes){
		val tmp = isEle(k(0),anno(chr))
		val tmp2 = isEle(k(0),anno2(chr))
		if (tmp._1) genes += tmp._3.replaceAll("\"","").replaceAll("-",":") + tmp._2
		if (tmp2._1) genes += tmp2._3.replaceAll("-",":") + ";" + tmp2._2
	}
	genes
}

def pedEventType (posCarriers: List[String]): String = {
	var mend, denovo, posDenovo, unknown, mendLC = 0
	for (individual <- posCarriers){
		if (families.contains(individual)){
			val cfam = families(individual)
			var kids = 0
			for (kid <- cfam._3) if (posCarriers.contains(kid)) kids += 1
			if (curDepth.contains(individual) && curDepth(individual) >= 18.0){// && kids >= 5){
				println(s"Core HC 5K:${individual} ${cfam._1} Car ${posCarriers.contains(cfam._1)} ${cfam._2} Car ${posCarriers.contains(cfam._2)}")
				if(posCarriers.contains(cfam._1) || posCarriers.contains(cfam._2)){
					mend += 1
				} else{
					println(s"popS${workingPop.contains(cfam._1)} popD${workingPop.contains(cfam._2)} ${if (curDepth.contains(cfam._1)) curDepth(cfam._1) else - 0.0 } ${if (curDepth.contains(cfam._2)) curDepth(cfam._2) else -0.0}")
					if (workingPop.contains(cfam._1) && workingPop.contains(cfam._2) && curDepth.contains(cfam._1) && curDepth(cfam._1) >= 18.0 && curDepth.contains(cfam._2) && curDepth(cfam._2) >= 18.0) {
						denovo += 1 
						println(s"Pro:${individual} Sire:${cfam._1} Dam:${cfam._2}")
						} else {
							unknown += 1
						}
				}
			} else {
				if(posCarriers.contains(cfam._1) || posCarriers.contains(cfam._2)){
					mendLC += 1
				} else{
					if (workingPop.contains(cfam._1) && workingPop.contains(cfam._2)) posDenovo += 1 else unknown += 1
				}

			}
		}
	}
	s"${mend}\t${denovo}\t${unknown}\t${posDenovo}\t${mendLC}"

}

/* Loop through all files */

while (fwd.hasNext){

	val F = fwd.next

	val fID = F.toString.split("/").last.split("\\.").apply(0)
	workingPop += fID
	//println(fID)
	val cEvents = new BufferedReader(new FileReader(F))
	populationSize += 1

/* Loop through all lines */

	while (cEvents.ready){
		val cEventLine = cEvents.readLine.split("\t")
		//cEventLine.foreach(s => print(s + " "))
		////HashMap[String -> HaspMap[String,Tuple2[Int,Int]]] Pos -> ID -> Reads

			/* Update Read mate patterns types */
	
	def updatePat (d: Array[Int], chrom: String, pos: Int) : Array[Int] = {
		var curData = d
		for (i <- 0 to 3){
		curData(i) = curData(i) + data(chrom)(pos)._7.apply(i).toInt
		}
		curData
	}
	

	def updateEvents (info: String, curEvents: HashMap[String,Tuple2[Int,Int]]) : HashMap[String,Tuple2[Int,Int]] = {
		/* evs are Event,FWD.int,REV.int */
		var workingEvs = curEvents
		if (info.size >= 1) for (evs <- info.substring(0,info.size -1).split(":")){
			val cevs = evs.split(",")
			if (workingEvs.contains(cevs(0))){
				workingEvs(cevs(0)) = (cevs(1).toInt + workingEvs(cevs(0))._1,cevs(2).toInt + workingEvs(cevs(0))._2) 
			} else {
				if (cevs.size >= 3) workingEvs += cevs(0) -> (cevs(1).toInt,cevs(2).toInt)
			}
		}
		workingEvs
	}
		
		/* If complete line process */
		if (cEventLine.size >= 10){
			/* Input format is: Window, FWD iPP, Rev IPP, FWD Split, REV Split, BreakPoints, Breakpoint Dif, Orientation EventTypes (Event,fwd,rev): */
			/* Store is Chr -> Pos -> Tuple5[Obs,Carriers,ImpropPaired,NumOFBreakpoints(List[Int]),# SplitReads,BreakDif HashSet[Int],MatePairPatArray[Int],MAP[Event(fwd,rev)]]*/
			val NumOFBreakpoints = cEventLine(5).split(",").size
			val chrom = "chr" + cEventLine(0).split(":").apply(0).trim.substring(3,cEventLine(0).split(":").apply(0).trim.size)			
			val start = cEventLine(0).split(":").apply(1).split("-").apply(0).toInt
			var modStart = if (data.contains(chrom) && data(chrom).contains(start.toInt - 500) ) (start - 500) else if (data.contains(chrom) && data(chrom).contains(start.toInt + 500)) (start + 500) else if (data.contains(chrom) && data(chrom).contains(start.toInt - 1000) ) (start - 1000) else if (data.contains(chrom) && data(chrom).contains(start.toInt + 500) ) (start + 500) else start
			val RepMatePat = cEventLine(7).split(":").map(_.toInt)
			
			var gts = Array(0,0)
			try {
			gts = cEventLine(9).split(":").map(_.toInt)
			} catch {
					case e: Exception => System.err.println(e + " " + fID)
				} 

		val depth = if (curDepth.contains(fID)) curDepth(fID) else cEventLine(9).split(":").map(_.toInt).sum
		val iPP = cEventLine(1).toInt + cEventLine(2).toInt
		val iPPGenome = depth / iPP

		if (iPPGenome < 3.0){

			if (data.contains(chrom)){
				if (data(chrom).contains(modStart)){
					/* Update the numbers */
					val wk = data(chrom)(modStart)
					val updatedEvents = updateEvents(cEventLine(8),data(chrom)(modStart)._8)
					data(chrom)(modStart) = (wk._1 + 1, fID :: wk._2, wk._3 + cEventLine(1).toInt + cEventLine(2).toInt, NumOFBreakpoints :: wk._4, wk._5 + cEventLine(3).toInt + cEventLine(4).toInt, wk._6 + cEventLine(6).toInt, updatePat(RepMatePat,chrom,modStart),updatedEvents,addBreaksInfo(wk._9,cEventLine(5)),addGenes(wk._10,cEventLine(5),chrom),wk._11 ++ HashMap(fID -> gts),updateCounts(wk._12,cEventLine(14)))
				} else{
					val updatedEvents = updateEvents(cEventLine(8),HashMap())
					data(chrom) += modStart -> (1, List(fID), cEventLine(1).toInt + cEventLine(2).toInt, List(NumOFBreakpoints), cEventLine(3).toInt + cEventLine(4).toInt, HashSet(cEventLine(6).toInt),RepMatePat,updatedEvents,addBreaksInfo(HashMap(),cEventLine(5)),addGenes(HashSet(),cEventLine(5),chrom),HashMap(fID -> gts),updateCounts(HashMap(0->0),cEventLine(14)))
				}
			} else {
				val updatedEvents = updateEvents(cEventLine(8),HashMap())
				data += chrom -> HashMap(modStart -> (1, List(fID), cEventLine(1).toInt + cEventLine(2).toInt, List(NumOFBreakpoints), cEventLine(3).toInt + cEventLine(4).toInt, HashSet(cEventLine(6).toInt),RepMatePat,updatedEvents,addBreaksInfo(HashMap(),cEventLine(5)),addGenes(HashSet(),cEventLine(5),chrom),HashMap(fID -> gts),updateCounts(HashMap(0->0),cEventLine(14))))
			}

		} // if Coverage

		}
	}

} // end of main while


/* Data should be full loaded so need to Analyse */

val analysis = new BufferedWriter(new FileWriter(new File("Analysis.tab")))
val analysisVCF = new BufferedWriter(new FileWriter(new File("genotypes.vcf")))

val chromOrder = data.keys.filter(s => ! List("chrX","chrM").contains(s)).toList.sortWith(_.slice(3,7).toInt < _.slice(3,7).toInt) ::: List("chrX")

analysis.write(s"CHROM:START-END\tPOP%\tType\tnPP_READS\tsplit_READS\tNUM-CARRIERS\tMendelian\tDenovo\tUnknown\tposDenovo\tLCMend\tR+M+\tR+M-\tR-M+\tR-M-\tMostCommonBreakpointDif\tBreakpointDIFs\tBreakpoints\tGene\tTop_Event\tEvents\tCARRIERS\n")


val refinedData = new HashMap[String,HashMap[Int,Tuple12[Int,List[String],Int,List[Int],Int,HashSet[Int],Array[Int],HashMap[String,Tuple2[Int,Int]],HashMap[Int,Int],HashSet[String],HashMap[String,Array[Int]],HashMap[Int,Int]]]]
println(data.size + " Data Size")
for (chr <- data){
	for (in <- chr._2){
		val tmp = in._2
	if(tmp._9.size > 0){	
		val sortedBreaks = if (tmp._9.size >= 2) tmp._9.toSeq.sortBy(- _._2).toArray else Array(tmp._9.head, (-1,0))
		val startPos = if (sortedBreaks(0)._1 > sortedBreaks(1)._1) sortedBreaks(1)._1 else sortedBreaks(0)._1
		if (refinedData.contains(chr._1)){
			if (refinedData(chr._1).contains(startPos)){
				//Combine windows
				//data(chrom)(modStart) = (wk._1 + 1, fID :: wk._2, wk._3 + cEventLine(1).toInt + cEventLine(2).toInt, NumOFBreakpoints :: wk._4, wk._5 + cEventLine(3).toInt + cEventLine(4).toInt, wk._6 + cEventLine(6).toInt, updatePat(RepMatePat,chrom,modStart),updatedEvents,addBreaksInfo(wk._9,cEventLine(5)),addGenes(wk._10,cEventLine(5),chrom))
				val pd = refinedData(chr._1)(startPos)
				refinedData(chr._1)(startPos) = (pd._1 + tmp._1, pd._2 ++ tmp._2, pd._3 + tmp._3, pd._4 ++ tmp._4, pd._5 + tmp._5, pd._6 ++ tmp._6, pd._7 ++ tmp._7,pd._8 ++ tmp._8, pd._9 ++ tmp._9, pd._10 ++ tmp._10, pd._11 ++ tmp._11, pd._12 ++ tmp._12)
			} else {
				//add window
				refinedData(chr._1) += startPos -> tmp
			}
		} else {
			refinedData += chr._1 -> HashMap(startPos -> tmp)
		}
     }	
	}
}

//println(refinedData.size + " Refined size")

analysisVCF.write("##fileformat=VCFv4.2\n##LocLTR\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
workingPop.toSeq.sorted.foreach(s => analysisVCF.write("\t" + s))
analysisVCF.write("\n")

for (chr <- chromOrder){
 if(refinedData.contains(chr)){
	val tmpDataOrder =  refinedData(chr).keys.toArray.sorted

	for (pos <- tmpDataOrder){
		println(chr + " " + pos)
/* (Count, List[ID's])*/
		
		val tmp = refinedData(chr)(pos)
		
	//	if(tmp._3 >= 5){
		/* Store is Chr -> Pos -> Tuple5[Obs,List[Carriers],ImpropPaired,NumOFBreakpoints(Array[Int]),# SplitReads,BreakDif HashSet[Int]]*/
		if (tmp._4.contains(2)){  // require at least one individual to have both breakpoints
		
		analysisVCF.write(s"${chr}\t${pos}\t.\tN\tNNNNNNNN\t.\t.\tEND=${pos + 8}\tGT:RO:AO")
		for (ind <- workingPop.toSeq.sorted){
			if (! tmp._11.contains(ind)){
				analysisVCF.write("\t0/0:.:.")
			} else {
				val ref = tmp._11(ind).apply(1)
				val alt = tmp._11(ind).apply(0)
				if (ref == 0 && alt >= 1) analysisVCF.write(s"\t1/1:${ref}:${alt}") else if (ref == 0 && alt == 0) analysisVCF.write("\t0/0:.:.") else analysisVCF.write(s"\t0/1:${ref}:${alt}")
			}
		}
		
		analysisVCF.write("\n")

		/*def typeOfEvent(in: HashMap[Int,Int]): String ={
			val fwd = in.toSeq.filter(_._1 >= 1)
			val rev = in.toSeq.filter(_._1 <= -1) //convert revs back to postives

			if (fwd.size >= 1 && rev.size >= 1){
				val best5p = fwd.sortBy(-_._2).toArray.apply(0)._1
				val best3p = rev.sortBy(-_._2).toArray.apply(0)._1 * -1
				if (best5p < best3p) {
					"DUP"
				} else {
					if (best5p > best3p) "DEL" else if (best5p == best3p) "EQL" else "UNK"
				}
			} else{
				"UNK"
			}

		} */
		/* Inverted due to logic bug somewhere else*/

		def typeOfEvent(in: HashMap[Int,Int]): String ={
			val fwd = in.toSeq.filter(_._1 >= 1)
			val rev = in.toSeq.filter(_._1 <= -1) //convert revs back to postives

			if (fwd.size >= 1 && rev.size >= 1){
				val best5p = fwd.sortBy(-_._2).toArray.apply(0)._1
				val best3p = rev.sortBy(-_._2).toArray.apply(0)._1 * -1
				if (best5p < best3p) {
					"DEL"
				} else {
					if (best5p > best3p) "DUP" else if (best5p == best3p) "EQL" else "UNK"
				}
			} else{
				"UNK"
			}

		} 

//		println(tmp._12)
			//val eType = if (tmp._12("DUP") && tmp._12("DEL")) "BOTH" else if (tmp._12("DUP")) "DUP" else if (tmp._12("DEL")) "DEL" else "UNK"

			val sortedBreaks = if (tmp._9.size >= 2) tmp._9.toSeq.sortBy(- _._2).toArray else Array(tmp._9.head, (-1,0))
			analysis.write(s"${chr}:${pos-500}-${pos+500}\t${tmp._1/populationSize.toFloat}\t${typeOfEvent(tmp._12)}\t${tmp._3}\t${tmp._5}\t${tmp._1}")
			analysis.write(s"\t" + pedEventType(tmp._2) + s"\t${tmp._7(0)}\t${tmp._7(1)}\t${tmp._7(2)}\t${tmp._7(3)}\t${scala.math.abs(sortedBreaks(0)._1 - sortedBreaks(1)._1) +1}\t${tmp._6.filterNot(s => s == 0)}\t")
			analysis.write(if (sortedBreaks(0)._1 <  sortedBreaks(1)._1) s"${sortedBreaks(0)._1}-${sortedBreaks(1)._1}" else s"${sortedBreaks(1)._1}-${sortedBreaks(0)._1}")
			analysis.write(s"\t${tmp._10}\t")
			val sortedEvents = tmp._8.toSeq.sortBy(s => - (s._2._1 + s._2._2))
			analysis.write(s"${sortedEvents(0)._1}:${sortedEvents(0)._2._1},${sortedEvents(0)._2._2}\t")
			sortedEvents.foreach(s => analysis.write(s"${s._1},${s._2._1},${s._2._2}:"))
			tmp._2.foreach(da => analysis.write("\t" + da))
			analysis.write("\n")
		}

	//	}
	}
	}

}
analysis.close
analysisVCF.close

//end
