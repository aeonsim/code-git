import scala.collection.mutable.HashMap
import org.apache.commons.io.FileUtils._
import java.io._

val fwd = org.apache.commons.io.FileUtils.listFiles(new File("."),Array("notPP.metrics.fwd.bedgraph"),true).iterator
//val rev = org.apache.commons.io.FileUtils.listFiles(new File("."),Array("notPP.metrics.rev.bedgraph"),true).iterator

val data = new HashMap[String,HashMap[Int,Tuple3[Int,List[String],Int]]]

var populationSize = 0

/* Loop through all files */

while (fwd.hasNext){
val F = fwd.next
val fID = F.toString.split("/").last.split("_").apply(0)

//val R = rev.next

val cF = new BufferedReader(new FileReader(F))
val cR = new BufferedReader(new FileReader(new File(F.toString.replace("fwd","rev"))))


//val rID = R.toString.split("/").last.split("_").apply(0)

// Add check to make sure opening the same individual
//if (fID == rID ){

populationSize += 1

/* Loop through all lines */

while (cF.ready && cR.ready){
val cfL = cF.readLine.split("\t")
val crL = cR.readLine.split("\t")

/* If item size is 4 analyse otherwise skip */

if (cfL.size == 4 && crL.size == 4){
val fKey = s"${cfL(0)}:${cfL(1)}-${cfL(2)}" 
val rKey = s"${crL(0)}:${crL(1)}-${crL(2)}" 

if (rKey == fKey){

if (cfL(3).toInt >= 1 || crL(3).toInt >= 1 ) {

if (data.contains(cfL(0))){

if (data(cfL(0)).contains(cfL(2).toInt)) {
data(cfL(0))(cfL(2).toInt) = (data(cfL(0))(cfL(2).toInt)._1 + 1, fID :: data(cfL(0))(cfL(2).toInt)._2,cfL(3).toInt + crL(3).toInt + data(cfL(0))(cfL(2).toInt)._3)
} else {
data(cfL(0)) += cfL(2).toInt -> (1 , fID :: Nil,cfL(3).toInt + crL(3).toInt)
}

} else {

data += cfL(0) -> HashMap( cfL(1).toInt -> (1, fID :: Nil,cfL(3).toInt + crL(3).toInt))

}

}

} else {
System.err.println("Unsynchronised Files " + cfL + "\t" + crL + "\t" + fwd.toString)
}

}

}
cF.close
cR.close

/*} else {
System.err.println("Incorrect paired Files " + fID + "\t" + rID)
}
*/
}

/* Data should be full loaded so need to Analyse */

val analysis = new BufferedWriter(new FileWriter(new File("Analysis.tab")))
val chromOrder = data.keys.filter(s => ! List("chrX","chrM").contains(s)).toList.sortWith(_.slice(3,7).toInt < _.slice(3,7).toInt) ::: List("chrX")

analysis.write(s"CHROM\tSTART\tEND\tPOP%\tREADS\tNUM-CARRIERS\tCARRIERS\n")

for (chr <- chromOrder){
val tmpDataOrder = data(chr).keys.toArray.sorted

for (pos <- tmpDataOrder){

/* (Count, List[ID's])*/

val tmp = data(chr)(pos)
if(tmp._3 >= 5){
analysis.write(s"${chr}\t${pos - 1000}\t${pos}\t${tmp._1/populationSize.toFloat}\t${tmp._3}\t${tmp._1}")
tmp._2.foreach(da => analysis.write("\t" + da))
analysis.write("\n")
}
}

}
analysis.close


