import java.io._
import scala.collection.mutable.HashMap
import scala.math.BigDecimal

//val args = Array("Mut_1x_cellCycle_151_bias_T_PGCs_40_NL2865_Pat.tab","NL286_shared_blue.txt")

val in = new BufferedReader(new FileReader(new File(args(0))))
val in_obs = new BufferedReader(new FileReader(new File(args(1))))
var data = new HashMap[Int,Long]
var dataProb = new HashMap[Int,Double]
var dataObs: List[Tuple2[Int,Int]] = Nil
var total = 0L

while (in.ready){
	val line = in.readLine.split("\t")
	if (line(0).contains("Muts_")){
		val mutNum = line(0).split("_")(1).toInt
		if (data.contains(mutNum)){
			for (k <- 1 to (line.size -1)) {
				data(mutNum) += line(k).toLong
				total += line(k).toInt
			}
		} else {
			for (k <- 1 to (line.size -1)) {
				data += mutNum -> line(k).toLong
				total += line(k).toInt
			}
		}
	}
}

for (i <- data){
	dataProb += i._1 -> i._2/total.toDouble
}


while(in_obs.ready){
	val line = in_obs.readLine.split("\t")
	dataObs = (line(0).toInt,line(1).toInt) :: dataObs
}
in_obs.close



val out = new BufferedWriter (new FileWriter(new File(args(0) + "_output.graph.txt")))

out.write("Numer of Spermcells " + dataObs.size + "\n")

for (k <- 0 to data.size -1){
	out.write(k + "\t" + "%.2f".format(data(k)/total.toDouble * 100) + "\t")
	val perc = ((data(k)/total.toDouble * 100)/0.5).toInt
	for (j <- 0 to perc) out.write("#")
	out.write("\n")

}

var likelihood = BigDecimal(1.0)

for (l <- dataObs){
	//println(l + "\t" + m + "\t" + (if (dataProb(l._1) == 0.0) likelihood * BigDecimal(1/numCells) else likelihood * BigDecimal(dataProb(l._1))))
		likelihood = if (dataProb(l._1) == 0.0) likelihood * BigDecimal(1/10000000.0) else likelihood * BigDecimal(dataProb(l._1))
}

out.write("Likelihood is " + likelihood + "\n")
out.close
println(args(0) + " | " + args(1) + " Likelihood is " + likelihood)