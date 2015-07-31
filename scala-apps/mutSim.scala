package mutSim

object mutSim {

  def main(args: Array[String]): Unit = {
    import org.apache.commons.math3.distribution.ChiSquaredDistribution
	import java.io._
	
    if(args.size == 0 || args.contains("-h")) {
      println("Error need arguments\nmutSim: mutrate recombRate popSize gens environmental loops maxMut cost? indv?")
      System.exit(1)
    }
    
    val outPutGens = Array(1,5,10,15,20,25,50,100,200,400,500,1000,1500,2000,2500,3000,4000,4500,4600,4900,5000,8000,8500,9000,9500)
    val mutRate = args(0).toDouble
    val recombRate = args(1).toDouble
    val popsize = args(2).toInt
    val generations = args(3).toInt
    val enviro = args(4).toInt
    val loops = args(5).toInt
    val maxMut = if ( args.size >= 7) args(6).toDouble else 1.5
    val cost = if (args.size == 8 && args(7) == "true") true else false
    val outIndv = if (args.size == 9 && args(8) == "true") true else false
    val out = new BufferedWriter(new FileWriter("Mut"+ mutRate+ "Rec" + recombRate + "Pop" + popsize + "Gens" +  generations + "Env" + enviro + ".loop" + loops + ".tab"))
    var data = new Array[Array[Double]](generations + 1)
    
    for (gens <- 0 to generations){
      data(gens) = new Array[Double](loops + 1)
    }
    
    for (lp <- 0 to loops){
    
val rand = scala.util.Random
var curGenNumb = 0

/* Storing Data as LocA1, LocA2, MutA1, MutA2  thus haplotype without Recomb is LocA1-MutA1 or locA2-MutA2 */


val chi = new org.apache.commons.math3.distribution.ChiSquaredDistribution(0.25)
val chiM = new org.apache.commons.math3.distribution.ChiSquaredDistribution(4)


def getMut(loci: Tuple2[Double,Double], mod: Double): Tuple2[Double,Double] = {
	if (rand.nextFloat <= mod){
		if (rand.nextFloat <= mod) (chi.sample, chi.sample) else (chi.sample, loci._2)
	} else {
		if (rand.nextFloat <= mod) (loci._1, chi.sample) else loci
	}
}

def getMutMut(loci: Tuple2[Double,Double], mod: Double): Tuple2[Double,Double] = {
	var mu1 = chiM.sample
	var mu2 = chiM.sample
	if (mu1 > maxMut) mu1 = (1 + 1/mu1)
	if (mu2 > maxMut) mu2 = (1 + 1/mu2)
	if (rand.nextFloat <= (mutRate * mod)){
		if (rand.nextFloat <= (mutRate * mod)) (mu1, mu2) else (mu1, 1)
	} else {
		if (rand.nextFloat <= (mutRate * mod)) (1, mu2) else (1,1)
	}
}

def selection(breedingPop: Array[Double]) : Int = {
var index = 0
val mate = rand.nextFloat()
while (mate > breedingPop(index)) index += 1
index
}

var curGen = new Array[Tuple4[Double,Double,Double,Double]](popsize)
var nextGen = new Array[Tuple4[Double,Double,Double,Double]](popsize)

var gcounter = 0 
while (gcounter < popsize){
//for (indv <- 0 to (popsize -1)){
curGen(gcounter) = (1.0,1.0,1.0,1.0)
gcounter += 1
}

var envCount = 0

while (curGenNumb <= generations){
val outputIndv  = if (outPutGens.contains(curGenNumb) && outIndv) true else false
val outIndivid = if (outputIndv) new BufferedWriter(new FileWriter("Mut"+ mutRate+ "Rec" + recombRate + "Pop" + popsize + "Gens" +  generations + "Env" + enviro + ".loop" + loops + "-Indivduals-" + curGenNumb + ".tab")) else new BufferedWriter(new FileWriter("Mut"+ mutRate+ "Rec" + recombRate + "Pop" + popsize + "Gens" +  generations + "Env" + enviro + ".loop" + loops + "-Indivduals-" + -1 + ".tab"))
var selectionSpace = 0.0
var mutrates = 0.0

if (envCount >= enviro){
		envCount = 0
		gcounter = 0
		  while (gcounter < popsize){
		    val indv = curGen(gcounter)
		    nextGen(gcounter) = (chi.sample,chi.sample,indv._3, indv._4)
		    gcounter += 1
		  }
		  curGen = nextGen
}

gcounter = 0
while (gcounter < popsize){
//for (indv <- curGen){
val indv = curGen(gcounter) 
selectionSpace += indv._1
selectionSpace += indv._2
mutrates += indv._3
mutrates += indv._4
gcounter += 1
}

data(curGenNumb)(lp) = ((mutrates/(2*popsize))*mutRate)
//out.write("Average Mutation rate for generation \t" + curGenNumb + "\t is \t" + (mutrates/(2*popsize))*mutRate + "\n")

var selectionMapping = new Array[Double](popsize)

var counter = 0
selectionMapping(0) = (curGen(counter)._1 + curGen(counter)._2)/selectionSpace

/* Gives array of Selection strengths*/

while (counter <= (popsize -2)){
counter += 1
selectionMapping(counter) = ((curGen(counter)._1 + curGen(counter)._2)/selectionSpace) + selectionMapping(counter - 1)
}

gcounter = 0
//for (indv <- 0 to (popsize -1)){
while (gcounter < popsize){
/* Selection mutator */
	val mate1 = selection(selectionMapping)
	var mate2 = selection(selectionMapping)
	/* NO self breeding */
	while (mate1 == mate2) mate2 = selection(selectionMapping)
	
	val m1MutMod = (curGen(mate1)._3 + curGen(mate1)._4)/2
	val m2MutMod = (curGen(mate2)._3 + curGen(mate2)._4)/2
	
	var m1Sel = getMut((curGen(mate1)._1,curGen(mate1)._2), m1MutMod)
	var m2Sel = getMut((curGen(mate2)._1,curGen(mate2)._2), m2MutMod)
	
	var m1Mut = getMutMut((curGen(mate1)._3,curGen(mate1)._4), m1MutMod)
	var m2Mut = getMutMut((curGen(mate2)._3,curGen(mate2)._4), m2MutMod)
	
	if (rand.nextFloat <= recombRate) m1Sel = (m1Sel._2,m1Sel._1)
	if (rand.nextFloat <= recombRate) m2Sel = (m2Sel._2,m2Sel._1)
	
	var m1Hap = if (rand.nextFloat >= 0.5) (m1Sel._1, curGen(mate1)._3 * m1Mut._1) else (m1Sel._2, curGen(mate1)._4 * m1Mut._2)
	var m2Hap = if (rand.nextFloat >= 0.5) (m2Sel._1, curGen(mate2)._3 * m2Mut._1) else (m2Sel._2, curGen(mate2)._4 * m2Mut._2)
	
	if (cost) {
	nextGen(gcounter) = (m1Hap._1 * (m1Hap._2 + m2Hap._2) * mutRate * mutRate, m2Hap._1 * (m1Hap._2 + m2Hap._2) * mutRate * mutRate, m1Hap._2, m2Hap._2) 
	} else {
		nextGen(gcounter) = (m1Hap._1, m2Hap._1, m1Hap._2, m2Hap._2)
	}
	
	if (outputIndv) outIndivid.write(curGenNumb + "\t" + curGen(gcounter)._1 + "\t" + curGen(gcounter)._2 + "\t" + (curGen(gcounter)._1 + curGen(gcounter)._2) +"\t" + (curGen(gcounter)._3 / 2 * mutRate) + "\t" + (curGen(gcounter)._4 / 2 * mutRate) + "\t" + ((curGen(gcounter)._3 + curGen(gcounter)._4)/2*mutRate) + "\n")	
	gcounter += 1
}
outIndivid.close
curGen = nextGen
nextGen =  new Array[Tuple4[Double,Double,Double,Double]](popsize)
envCount += 1

curGenNumb += 1

}
    }
    for (gens <- 0 to generations){
      out.write(gens + "\t\t\t")
      data(gens).foreach(loop => out.write("\t" + (if (loop > maxMut ) "Dead" else loop )))
      out.write("\n")
    }
    out.close
  }
  
}