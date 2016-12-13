import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet
import org.apache.commons.io.FileUtils._
import java.io._

//run against folder of model, provide type of input you want (1K-blood 1K-sperm etc) and the number of oocytes to use...

/* to filter with Awk | awk '{n=split($0,col,"\t"); y=split(col[8],p,":"); if ((p[1]+p[2]+p[3]+p[4]) >= 5 && col[5]/col[7] >= 2 && col[7] >= 2 && col[5] >= 10 && col[4] < 0.2) {print $0};}'*/
println("simSampler /path/to/output/ file_suffix sampleNumber")

val fwd = org.apache.commons.io.FileUtils.listFiles(new File(args(0)),Array("tab"),true).iterator
val rand = scala.util.Random

val sims = new HashMap[String,List[String]]

while (fwd.hasNext){

	val F = fwd.next
	val fID = F.toString.split("/").last.split("-").apply(0)
if (F.toString.contains(args(1))){	
	val cEvents = new BufferedReader(new FileReader(F))

	while(cEvents.ready){
		val spermCell = cEvents.readLine
		if (sims.contains(fID)){
				sims(fID) = spermCell :: sims(fID)
			} else {
				sims += fID -> List(spermCell)
			}
	}

	val data = sims(fID).toArray
	val simsSperm = new HashMap[String,Double]

	for (i <- 1 to args(2).toInt){
		//print( i + " ")
		val cSperm = data(rand.nextInt(data.size)).split("\t")

		var n = 2
		while (n < cSperm.size){
			if (simsSperm.contains(cSperm(n))){

				} else {
					simsSperm += cSperm(n) -> cSperm(n + 1).toDouble
				}
			n += 2
		}


	}

	//print("\n")

	//simsSperm.toSeq.sortBy(-_._2).foreach(spm => print(spm._1 + "\t"))
	//print("\n")
	simsSperm.toSeq.sortBy(-_._2).foreach(spm => print(spm._2 + "\t"))
	print("\n")

}

}