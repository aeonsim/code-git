import java.io._
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet
import scala.math.BigDecimal
import org.apache.commons.io.FileUtils._
import org.apache.commons.math3.distribution.BinomialDistribution

/* scala app <input/file/folder> <suffix> <observedARs> <binSize x.x> <Upperlimit x.x> <SampleSize>*/

//11_off_0.85-0.85-0.4_PGCS_40_bias_FALSE/1-1K-sperm-5y.tab

//scala likelihoods.scala 11_off_5.1-0.5-0.4_PGCS_4_bias_FALSE/55-1K-sperm-5y.tab - NL288_Blue.txt 0.05
//val args = Array("11_off_5.1-0.5-0.4_PGCS_4_bias_FALSE/55-1K-sperm-5y.tab","1K-sperm-5y.tab","NL288_Blue.txt","0.05","0.3")

val fwd = org.apache.commons.io.FileUtils.iterateFiles(new File(args(0)),null,false)

val in_obs = new BufferedReader(new FileReader(new File(args(2))))


var data = new HashMap[Int,Int]
var dataObs = new HashMap[Int,Int]
var dataProb = new HashMap[Int,Double]
var dataUniqProb = new HashMap[Int,Double]
var allData = new HashMap[String,Tuple2[Double,Int]]
var allPools = new HashMap[Double,Int]
var dataSampledGam = new HashMap[Int,Int]

val binsize = args(3).toFloat
val upperLim = args(4).toDouble
var numCells, numMuts, numObsMut, uniqueNum = 0
val sampleSize = args(5).toInt
val detectThresh = 0.13
val avgSeqDepth = 24

var biNoms = new HashMap[Double,BinomialDistribution]
biNoms += 0.5 -> new BinomialDistribution(avgSeqDepth,0.5)

def keep(ar: Double, obs: Int) : Tuple2[Boolean,Double] = {
	val lar = if (ar >= 1.0) ar/1000 else ar
	val pWild = 1 - scala.math.pow((0.5 / ((0.5 - lar) + 0.5).toFloat),sampleSize - obs)

	if (biNoms.contains(lar)){
			val reads = biNoms(lar).sample
			val obsAR = reads/avgSeqDepth.toDouble
			return Tuple2(obsAR <= 0.33 && scala.util.Random.nextFloat <= pWild,obsAR)
		} else {
			biNoms += lar -> new BinomialDistribution(avgSeqDepth,ar)
			val reads = biNoms(lar).sample
			val obsAR = reads/avgSeqDepth.toDouble
			//Tuple2(0.05 >= biNoms(0.5).cumulativeProbability(reads) && scala.util.Random.nextFloat <= pWild,obsAR)
			return Tuple2(obsAR <= 0.33 && scala.util.Random.nextFloat <= pWild,obsAR)

		}

}


for (i <- 0 to (1/binsize).toInt){
	data += i -> 0
	dataObs += i -> 0
	dataProb += i -> 0.0
	dataUniqProb += i -> 0.0
	dataSampledGam += i -> 0
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
		var localData : List[Array[String]] = Nil

		while (in_data.ready){
			var lc = 2
			val line = in_data.readLine.split("\t")
			localData = line :: localData
			numCells += 1
			while (lc <= line.size){
				if (allData.contains(s"${line(lc - 1)}_${line(lc)}")) allData(s"${line(lc - 1)}_${line(lc)}") = Tuple2(allData(s"${line(lc - 1)}_${line(lc)}")._1, allData(s"${line(lc - 1)}_${line(lc)}")._2 + 1) else allData += s"${line(lc - 1)}_${line(lc)}" -> Tuple2(line(lc).toDouble,1)
				if (line(lc).toDouble <= upperLim) data((line(lc).toDouble/binsize).toInt) += 1
				lc += 2
				numMuts += 1
			}
		}
		in_data.close

		/* Need to Sample 100 sets of X gametes */
		var x = 0
		while (x < 999){
			x += 1
			var uniqueGam: List[Double] =Nil
			//for (i <- 1 to sampleSize){
				//val gam = localData(scala.util.Random.nextInt(localData.size))
				val gam = localData(x)
				var lc = 2
				while(lc <= gam.size){
					uniqueGam = gam(lc).toDouble :: uniqueGam
					lc += 2
				}
			
			var tokeep : List[Double] = Nil

			for ( j <- uniqueGam){
				val lkeep = keep(j,1)
				if (lkeep._1) {
					dataSampledGam((lkeep._2.toFloat / binsize).toInt) += 1
					uniqueNum += 1
					tokeep = lkeep._2 :: tokeep
				}
			}
			if (x % 5 == 0) {
				tokeep.sortBy(-_).foreach(s => print(s + "\t"))
				print("\n")
			}
		}
	}
	
}

val out = new BufferedWriter (new FileWriter(new File(args(0) + "OM.SM.output.graph.txt")))

out.write("Numer of Spermcells " + numCells + "\n")
out.write("Average number of de novos " + numMuts/numCells.toDouble + " vs " + numObsMut + " observed de novos\n")



for (k <- 0 to data.size -1){
	out.write("Base\t" + "%.2f".format(k * binsize) + "\t" + "%.2f".format(data(k)/numMuts.toDouble * 100) + "\t" + data(k) + "\t")
	val perc = ((data(k)/numMuts.toDouble * 100)/0.5).toInt
	for (j <- 0 to perc) out.write("#")
	out.write("\n")

	out.write( "Biased\t" + "%.2f".format(k * binsize) + "\t" + "%.2f".format(dataSampledGam(k)/uniqueNum.toDouble * 100) + "\t" + dataSampledGam(k) + "\t")
	val percU = ((dataSampledGam(k)/uniqueNum.toDouble * 100)/0.5).toInt
	for (j <- 0 to percU) out.write("=")
	out.write("\n")

	out.write( "Observed\t" + "%.2f".format(k * binsize) + "\t" + "%.2f".format(dataObs(k)/numObsMut.toDouble * 100) + "\t" + dataObs(k) + "\t")
	val percObs = ((dataObs(k)/numObsMut.toDouble * 100)/0.5).toInt
	for (j <- 0 to percObs) out.write("&")
	out.write("\n")

	dataProb(k) = (data(k)/numMuts.toDouble)
	dataUniqProb(k) = (dataSampledGam(k)/uniqueNum.toDouble)

}

var likelihood, uniqLike = BigDecimal(1.0)

for (l <- dataObs){
	for (m <- 1 to l._2) {
	//println(l + "\t" + m + "\t" + (if (dataProb(l._1) == 0.0) likelihood * BigDecimal(1/numCells) else likelihood * BigDecimal(dataProb(l._1))))
		likelihood = if (dataProb(l._1) == 0.0) likelihood * BigDecimal(1/(numCells.toDouble+1)) else likelihood * BigDecimal(dataProb(l._1))
		uniqLike = if (dataUniqProb(l._1) == 0.0) uniqLike * BigDecimal(1/(uniqueNum.toDouble+1)) else uniqLike * BigDecimal(dataUniqProb(l._1))
	}
}

out.write("Likelihood is " + likelihood + "\n")
out.write("Unique Likelihood is " + uniqLike + "\n")
out.close
System.err.println(args(0) + "*" + args(1) + "\t" + args(2) + "\t" + numCells + " cells " + numMuts/numCells.toDouble + " simulated de novos vs " + numObsMut + " observed de novos, Likelihood is " + likelihood + " UniqueLike " + uniqLike)