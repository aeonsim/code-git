import java.io._
import scala.collection.mutable.HashMap
import scala.math.BigDecimal
import org.apache.commons.io.FileUtils._

/* scala app <input/file/folder> <suffix> <observedARs> <binSize x.x> <Upperlimit x.x>*/

//11_off_0.85-0.85-0.4_PGCS_40_bias_FALSE/1-1K-sperm-5y.tab

//scala likelihoods.scala 11_off_5.1-0.5-0.4_PGCS_4_bias_FALSE/55-1K-sperm-5y.tab - NL288_Blue.txt 0.05
//val args = Array("11_off_5.1-0.5-0.4_PGCS_4_bias_FALSE/55-1K-sperm-5y.tab","1K-sperm-5y.tab","NL288_Blue.txt","0.05","0.3")

val fwd = org.apache.commons.io.FileUtils.iterateFiles(new File(args(0)),null,false)

val in_obs = new BufferedReader(new FileReader(new File(args(2))))


var data = new HashMap[Int,Int]
var dataObs = new HashMap[Int,Int]
var dataProb = new HashMap[Int,Double]

val binsize = args(3).toDouble
val upperLim = args(4).toDouble
var numCells, numMuts, numObsMut = 0

for (i <- 0 to (1/binsize).toInt){
	data += i -> 0
	dataObs += i -> 0
	dataProb += i -> 0.0
}

while(in_obs.ready){
	val line = in_obs.readLine
	dataObs((line.toDouble/binsize).toInt) += 1
	numObsMut += 1
}
in_obs.close

while (fwd.hasNext){

	val F = fwd.next
	if (F.toString.contains(args(1))) {
		val in_data = new BufferedReader(new FileReader(F))


		while (in_data.ready){
			var lc = 2
			val line = in_data.readLine.split("\t")
			numCells += 1
			while (lc <= line.size){
				if (line(lc).toDouble <= upperLim) data((line(lc).toDouble/binsize).toInt) += 1
				lc += 2
				numMuts += 1
			}
		}
		in_data.close
	}
	
}

val out = new BufferedWriter (new FileWriter(new File(args(0) + "output.graph.txt")))

out.write("Numer of Spermcells " + numCells + "\n")
out.write("Average number of de novos " + numMuts/numCells.toDouble + " vs " + numObsMut + " observed de novos\n")



for (k <- 0 to data.size -1){
	out.write("%.2f".format(k * binsize) + "\t" + "%.2f".format(data(k)/numMuts.toDouble * 100) + "\t")
	val perc = ((data(k)/numMuts.toDouble * 100)/0.5).toInt
	for (j <- 0 to perc) out.write("#")
	out.write("\n")

	out.write("%.2f".format(k * binsize) + "\t" + "%.2f".format(dataObs(k)/numObsMut.toDouble * 100) + "\t")
	val percObs = ((dataObs(k)/numObsMut.toDouble * 100)/0.5).toInt
	for (j <- 0 to percObs) out.write("&")
	out.write("\n")

	dataProb(k) = (data(k)/numMuts.toDouble)

}

var likelihood = BigDecimal(1.0)

for (l <- dataObs){
	for (m <- 1 to l._2) {
	//println(l + "\t" + m + "\t" + (if (dataProb(l._1) == 0.0) likelihood * BigDecimal(1/numCells) else likelihood * BigDecimal(dataProb(l._1))))
		likelihood = if (dataProb(l._1) == 0.0) likelihood * BigDecimal(1/(numCells.toDouble+1)) else likelihood * BigDecimal(dataProb(l._1))

	}
}

out.write("Likelihood is " + likelihood + "\n")
out.close
println(args(0) + "*" + args(1) + "\t" + args(2) + "\t" + numCells + " cells " + numMuts/numCells.toDouble + " simulated de novos vs " + numObsMut + " observed de novos, Likelihood is " + likelihood)