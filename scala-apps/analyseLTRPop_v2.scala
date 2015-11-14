import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet
import org.apache.commons.io.FileUtils._
import java.io._

/* to filter with Awk | awk '{n=split($0,col,"\t"); y=split(col[8],p,":"); if ((p[1]+p[2]+p[3]+p[4]) >= 5 && col[5]/col[7] >= 2 && col[7] >= 2 && col[5] >= 10 && col[4] < 0.2) {print $0};}'*/

val fwd = org.apache.commons.io.FileUtils.listFiles(new File("."),Array("event.tab"),true).iterator
//val sFwd = org.apache.commons.io.FileUtils.listFiles(new File("."),Array("split.fwd.bedgraph"),true).iterator

val depths = new BufferedReader(new FileReader(new File("/scratch/aeonsim/LTR-retros/currentDepths.txt")))

/* Chr -> Pos -> Tuple5[Obs,Carriers,ImpropPaired,NumOFBreakpoints(Array[Int]),# SplitReads,BreakDif HashSet[Int], Orination Array, Event HashMap, Breakpoints HashSet]*/

val data = new HashMap[String,HashMap[Int,Tuple9[Int,List[String],Int,List[Int],Int,HashSet[Int],Array[Int],HashMap[String,Tuple2[Int,Int]],HashSet[Int]]]]

val curDepth = new HashMap[String,Double]
while (depths.ready){
val cline = depths.readLine.split("\t")
curDepth += cline(0) -> cline(1).toDouble
}
depths.close

/* HashMap[Proband -> (Sire, Dam, List[Children])]*/


val families = new HashMap[String,Tuple3[String,String,List[String]]]

val popIn = new BufferedReader(new FileReader(new File("/scratch/aeonsim/vcfs/ref-peds/Damona-Full-8Jan15.ped")))
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

for (indv <- pop){
if (pop.contains(indv._2(2)) && pop.contains(indv._2(3))){
 val kids = findChildren(pop.values.toList,indv._1)
 if (kids.size >= 1) families += indv._1 -> (indv._2(2),indv._2(3),kids)
}
}

def addBreaks (breaks: HashSet[Int], newBreaks: String ): HashSet[Int] ={
	val novoBreaks = newBreaks.substring(4,newBreaks.size -1).split(",").map(s => s.trim.toInt)
	for ( k <- novoBreaks){
		breaks += k
	}
	breaks
}

def pedEventType (posCarriers: List[String]): String = {
	var fullped, partped, denovo, junk = 0
	for (individual <- posCarriers){
	/* if core family member */
		if (families.contains(individual) && curDepth.contains(individual) && curDepth(individual) >= 15){
			val cfam = families(individual)
			if (posCarriers.contains(cfam._1) || posCarriers.contains(cfam._2)){
				var kids = 0
				for (kid <- cfam._3) if (posCarriers.contains(kid)) kids += 1
				if ( kids >= 1 ){
					fullped += 1
				} else {
					partped += 1
				}
			} else {
			/* Not in Parents */
				if (workingPop.contains(cfam._1) && workingPop.contains(cfam._2)){
					var kids = 0
					for (kid <- cfam._3) if (posCarriers.contains(kid) && curDepth.contains(kid) && curDepth(kid) >= 15) kids += 1
					if ( kids >= 1 ){
						denovo += 1
					} else {
						junk += 1
					}
				}
			}

		}


	}
	s"${fullped}\t${partped}\t${denovo}\t${junk}"

}


/* Loop through all files */

while (fwd.hasNext){

	val F = fwd.next
//try {
	//val S = sFwd.next
	val fID = F.toString.split("/").last.split("_").apply(0)
	workingPop += fID

	val cEvents = new BufferedReader(new FileReader(F))
	populationSize += 1

/* Loop through all lines */

	while (cEvents.ready){
		val cEventLine = cEvents.readLine.split("\t")

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
		for (evs <- info.split(":")){
			val cevs = evs.split(",")
			if (workingEvs.contains(cevs(0))){
				workingEvs(cevs(0)) = (cevs(1).toInt + workingEvs(cevs(0))._1,cevs(2).toInt + workingEvs(cevs(0))._2) 
			} else {
				workingEvs += cevs(0) -> (cevs(1).toInt,cevs(2).toInt)
			}
		}
		workingEvs
	}
		
		/* If complete line process */
		if (cEventLine.size >= 8){
			/* Input format is: Window, FWD iPP, Rev IPP, FWD Split, REV Split, BreakPoints, Breakpoint Dif, Orientation EventTypes (Event,fwd,rev): */
			/* Store is Chr -> Pos -> Tuple5[Obs,Carriers,ImpropPaired,NumOFBreakpoints(List[Int]),# SplitReads,BreakDif HashSet[Int],MatePairPatArray[Int],MAP[Event(fwd,rev)]]*/
			val NumOFBreakpoints = cEventLine(5).split(",").size
			val chrom = cEventLine(0).split(":").apply(0).trim

			val start = cEventLine(0).split(":").apply(1).split("-").apply(0).toInt
			val modStart = if (data.contains(chrom) && data(chrom).contains(start.toInt - 500) ) (start - 500) else if (data.contains(chrom) && data(chrom).contains(start.toInt + 500)) (start + 500) else start
			val RepMatePat = cEventLine(7).split(":").map(_.toInt)
			

			if (data.contains(chrom)){
				if (data(chrom).contains(modStart)){
					/* Update the numbers */
					val wk = data(chrom)(modStart)
					val updatedEvents = updateEvents(cEventLine(8),data(chrom)(modStart)._8)
					data(chrom)(modStart) = (wk._1 + 1, fID :: wk._2, wk._3 + cEventLine(1).toInt + cEventLine(2).toInt, NumOFBreakpoints :: wk._4, wk._5 + cEventLine(3).toInt + cEventLine(4).toInt, wk._6 + cEventLine(6).toInt, updatePat(RepMatePat,chrom,modStart),updatedEvents,addBreaks(wk._9,cEventLine(5)))
				} else{
					val updatedEvents = updateEvents(cEventLine(8),HashMap())
					data(chrom) += modStart -> (1, List(fID), cEventLine(1).toInt + cEventLine(2).toInt, List(NumOFBreakpoints), cEventLine(3).toInt + cEventLine(4).toInt, HashSet(cEventLine(6).toInt),RepMatePat,updatedEvents,addBreaks(HashSet(),cEventLine(5)))
				}
			} else {
				val updatedEvents = updateEvents(cEventLine(8),HashMap())
				data += chrom -> HashMap(modStart -> (1, List(fID), cEventLine(1).toInt + cEventLine(2).toInt, List(NumOFBreakpoints), cEventLine(3).toInt + cEventLine(4).toInt, HashSet(cEventLine(6).toInt),RepMatePat,updatedEvents,addBreaks(HashSet(),cEventLine(5))))
			}

		}
	}

/*} catch {
case e: Exception => System.err.println(e + " " + F.toString)

}*/

} // end of main while


/* Data should be full loaded so need to Analyse */

val analysis = new BufferedWriter(new FileWriter(new File("Analysis.tab")))
val chromOrder = data.keys.filter(s => ! List("chrX","chrM").contains(s)).toList.sortWith(_.slice(3,7).toInt < _.slice(3,7).toInt) ::: List("chrX")

analysis.write(s"CHROM:START-END\tPOP%\tnPP_READS\tsplit_READS\tNUM-CARRIERS\t3gen\tTrio\tDenovo\tJunk\tR+M+\tR+M-\tR-M+\tR-M-\tBreakpointDIF\tBreakpoints\tEvents\tCARRIERS\n")

for (chr <- chromOrder){
	val tmpDataOrder = data(chr).keys.toArray.sorted

	for (pos <- tmpDataOrder){

/* (Count, List[ID's])*/
		
		val tmp = data(chr)(pos)
	//	if(tmp._3 >= 5){
		/* Store is Chr -> Pos -> Tuple5[Obs,List[Carriers],ImpropPaired,NumOFBreakpoints(Array[Int]),# SplitReads,BreakDif HashSet[Int]]*/
		if (tmp._4.contains(2)){  // require at least one individual to have both breakpoints
			analysis.write(s"${chr}:${pos}-${pos + 1500}\t${tmp._1/populationSize.toFloat}\t${tmp._3}\t${tmp._5}\t${tmp._1}")
			analysis.write(s"\t" + pedEventType(tmp._2) + s"\t${tmp._7(0)}\t${tmp._7(1)}\t${tmp._7(2)}\t${tmp._7(3)}\t${tmp._6}\t${tmp._9}\t")
			tmp._8.foreach(s => analysis.write(s"${s._1},${s._2._1},${s._2._2}:"))
			tmp._2.foreach(da => analysis.write("\t" + da))
			analysis.write("\n")
		}

	//	}
	}

}
analysis.close
