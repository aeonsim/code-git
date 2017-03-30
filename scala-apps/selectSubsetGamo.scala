import java.io._
import scala.collection.mutable.HashSet

val ins = new BufferedReader(new FileReader(new File(args(0))))
val max = args(1).toDouble
val num = args(2).toInt

while (ins.ready){
	val tmp = ins.readLine
	val data = new BufferedReader(new FileReader(new File(tmp)))
	val muts = new HashSet[Tuple2[String,Double]]
	var lines : List[Int] = Nil
	for (i <- 1 to num) lines = scala.util.Random.nextInt(1000) :: lines

	var c = 0
	 while (data.ready && c <= lines.sorted.last){
	 	val cline = data.readLine.split("\t")
	 	if (lines.contains(c)){
	 		for (i <- 1 to ((cline.size -1)/2) ){
	 			val id = cline((2*i) - 1)
	 			val ar = cline(2*i).toDouble
	 			if (ar <= max) muts += Tuple2(id,ar)
	 		}

	 	}
	 	c += 1
	 }

	 muts.toArray.sortBy(s => - s._2).foreach(k => print("\t" + k._2))
	 print("\n")


}