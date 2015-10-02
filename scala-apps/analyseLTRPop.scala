import scala.collection.mutable.HashMap
import org.apache.commons.io.FileUtils._
import java.io._

/* to filter with Awk | awk '{n=split($0,col,"\t"); y=split(col[8],p,":"); if ((p[1]+p[2]+p[3]+p[4]) >= 5 && col[5]/col[7] >= 2 && col[7] >= 2 && col[5] >= 10 && col[4] < 0.2) {print $0};}'*/

val fwd = org.apache.commons.io.FileUtils.listFiles(new File("."),Array("notPP.fwd.bedgraph"),true).iterator
//val sFwd = org.apache.commons.io.FileUtils.listFiles(new File("."),Array("split.fwd.bedgraph"),true).iterator

/* Chr -> Pos -> Tuple5[Obs,Carriers,Reads,Array[Int],SplitReads]*/

val data = new HashMap[String,HashMap[Int,Tuple5[Int,List[String],Int,Array[Int],Int]]]

/* HashMap[Proband -> (Sire, Dam, List[Children])]*/


val families = new HashMap[String,Tuple3[String,String,List[String]]]

val popIn = new BufferedReader(new FileReader(new File("/scratch/aeonsim/vcfs/ref-peds/Damona-Full-8Jan15.ped")))
val pop = new HashMap[String,Array[String]]

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

def pedEventType (posCarriers: List[String]): String = {
	var fullped, partped, denovo, junk = 0
	for (individual <- posCarriers){
	/* if core family member */
		if (families.contains(individual)){
			val cfam = families(individual)
			if (posCarriers.contains(cfam._1) || posCarriers.contains(cfam._1)){
				var kids = 0
				for (kid <- cfam._3) if (posCarriers.contains(kid)) kids += 1
				if ( kids >= 1 ){
					fullped += 1
				} else {
					partped += 1
				}
			} else {
			/* Not in Parents */
				var kids = 0
				for (kid <- cfam._3) if (posCarriers.contains(kid)) kids += 1
				if ( kids >= 1 ){
					denovo += 1
				} else {
					junk += 1
				}
				
			}

		}


	}
	s"${fullped}:${partped}:${denovo}:${junk}"

}

/* Loop through all files */

while (fwd.hasNext){

	val F = fwd.next
try {
	//val S = sFwd.next
	val fID = F.toString.split("/").last.split("_").apply(0)

	val cF = new BufferedReader(new FileReader(F))
	val cR = new BufferedReader(new FileReader(new File(F.toString.replace("fwd","rev"))))
	val sF = new BufferedReader(new FileReader(new File(F.toString.replace("notPP","split"))))
	val sR = new BufferedReader(new FileReader(new File(F.toString.replace("notPP.fwd","split.rev"))))
	populationSize += 1

/* Loop through all lines */

	while (cF.ready && cR.ready){
		val cfL = cF.readLine.split("\t")
		val crL = cR.readLine.split("\t")
		val csfL = sF.readLine.split("\t")
		val csrL = sR.readLine.split("\t")

/* If item size is 4 analyse otherwise skip */

		if (cfL.size >= 5 && crL.size >= 5){
			val fKey = s"${cfL(0)}:${cfL(1)}-${cfL(2)}" 
			val rKey = s"${crL(0)}:${crL(1)}-${crL(2)}" 

			if (rKey == fKey){

				if (cfL(3).toInt >= 1 || crL(3).toInt >= 1 || csfL(3).toInt >= 1 || csrL(3).toInt >= 1 ) {
					val RepMatePat = cfL(4).split(":").map(_.toInt)

	/* Update Read mate patterns types */
	
					def updatePat (d: Array[Int]) : Array[Int] = {
						var curData = d
						for (i <- 0 to 3){
						curData(i) = curData(i) + data(cfL(0))(cfL(2).toInt)._4.apply(i).toInt
						}
						curData
					}
					
	/* Calc number of Split Reads */
	
					val splitReads = csfL(3).toInt + csrL(3).toInt
					val notPP = cfL(3).toInt + crL(3).toInt

	/*	IF Chromosome is present check for Pos else add Chrom */

					if (data.contains(cfL(0))){

	/*	If Pos is present add data else add pos then data */

						if (data(cfL(0)).contains(cfL(2).toInt)) {
						/* Chr -> Pos -> Tuple5[Obs,Carriers,Reads,Array[Int],SplitReads]*/

							data(cfL(0))(cfL(2).toInt) = (data(cfL(0))(cfL(2).toInt)._1 + 1, fID :: data(cfL(0))(cfL(2).toInt)._2,notPP + data(cfL(0))(cfL(2).toInt)._3,updatePat(RepMatePat),splitReads + data(cfL(0))(cfL(2).toInt)._5)
							
						} else {
							data(cfL(0)) += cfL(2).toInt -> (1 , fID :: Nil,notPP,RepMatePat,splitReads)
						}

					} else {

						data += cfL(0) -> HashMap( cfL(1).toInt -> (1, fID :: Nil,notPP,RepMatePat,splitReads))

					}
				}
				

			} else {
				System.err.println("Unsynchronised Files " + cfL + "\t" + crL + "\t" + fwd.toString)
			}

		}

	}
	cF.close
	cR.close
	sF.close
	sR.close
} catch {
case e: Exception => System.err.println(e + " " + F.toString)

}

} // end of main while


/* Data should be full loaded so need to Analyse */

val analysis = new BufferedWriter(new FileWriter(new File("Analysis.tab")))
val chromOrder = data.keys.filter(s => ! List("chrX","chrM").contains(s)).toList.sortWith(_.slice(3,7).toInt < _.slice(3,7).toInt) ::: List("chrX")

analysis.write(s"CHROM\tSTART\tEND\tPOP%\tnPP_READS\tsplit_READS\tNUM-CARRIERS\tLTR_R_++,+-,-+,--\tFullPed:PartPED:Denovo:Junk\tCARRIERS\n")

for (chr <- chromOrder){
	val tmpDataOrder = data(chr).keys.toArray.sorted

	for (pos <- tmpDataOrder){

/* (Count, List[ID's])*/
		
		val tmp = data(chr)(pos)
	//	if(tmp._3 >= 5){
			analysis.write(s"${chr}\t${pos - 1000}\t${pos}\t${tmp._1/populationSize.toFloat}\t${tmp._3}\t${tmp._5}\t${tmp._1}")
			analysis.write(s"\t${tmp._4(0)}:${tmp._4(1)}:${tmp._4(2)}:${tmp._4(3)}\t" + pedEventType(tmp._2))
			tmp._2.foreach(da => analysis.write("\t" + da))
			analysis.write("\n")
	//	}
	}

}
analysis.close


