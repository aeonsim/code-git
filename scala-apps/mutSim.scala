package mutSim

object mutSim {

  def main(args: Array[String]): Unit = {
    import org.apache.commons.math3.distribution.ChiSquaredDistribution
	import java.io._
	
    if(args.size == 0 || args.contains("-h")) {
      println("Error need arguments\nmutSim: mutrate recombRate popSize gens environmental loops")
      System.exit(1)
    }
    
    val mutRate = args(0).toDouble
    val recombRate = args(1).toDouble
    val popsize = args(2).toInt
    val generations = args(3).toInt
    val enviro = args(4).toInt
    val loops = args(5).toInt
    val out = new BufferedWriter(new FileWriter("Mut"+ mutRate+ "Rec" + recombRate + "Pop" + popsize + "Gens" +  generations + ".loop" + loops + ".tab"))

    var data = new Array[Array[Double]](generations + 1)
    
    for (gens <- 0 to generations){
      data(gens) = new Array[Double](loops + 1)
    }
    
    for (lp <- 0 to loops){
    
val rand = scala.util.Random
var curGenNumb = 0

/* Storing Data as LocA1, LocA2, MutA1, MutA2  thus haplotype without Recomb is LocA1-MutA1 or locA2-MutA2 */


val chi = new org.apache.commons.math3.distribution.ChiSquaredDistribution(0.25)


def getMut(loci: Tuple2[Double,Double], mod: Double): Tuple2[Double,Double] = {
	if (rand.nextFloat <= (mutRate * mod)){
		if (rand.nextFloat <= (mutRate * mod)) (chi.sample, chi.sample) else (chi.sample, loci._2)
	} else {
		if (rand.nextFloat <= (mutRate * mod)) (loci._1, chi.sample) else loci
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

var selectionSpace = 0.0
var mutrates = 0.0

if (envCount == enviro){
		envCount = 0
		gcounter = 0
		  while (gcounter < popsize){
		    val indv = curGen(gcounter)
		    curGen(gcounter) = (chi.sample,chi.sample,indv._3, indv._4)
		    gcounter += 1
		  }
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
	
	var m1Mut = getMut((curGen(mate1)._3,curGen(mate1)._4), m1MutMod)
	var m2Mut = getMut((curGen(mate2)._3,curGen(mate2)._4), m2MutMod)
	
	if (rand.nextFloat <= recombRate) m1Sel = (m1Sel._2,m1Sel._1)
	if (rand.nextFloat <= recombRate) m2Sel = (m2Sel._2,m2Sel._1)
	
	var m1Hap = if (rand.nextFloat >= 0.5) (m1Sel._1, m1Mut._1) else (m1Sel._2, m1Mut._2)
	var m2Hap = if (rand.nextFloat >= 0.5) (m2Sel._1, m2Mut._1) else (m2Sel._2, m2Mut._2)
	
	nextGen(gcounter) = (m1Hap._1, m2Hap._1, m1Hap._2, m2Hap._2)
		
	//if (rand.nextFloat <= recombRate) nextGen(indv) = (nextSel._2,nextSel._1,nextMut._1 * curGen(indv)._3,nextMut._2 * curGen(indv)._4) else nextGen(indv) = (nextSel._1,nextSel._2,nextMut._1 * curGen(indv)._3,nextMut._2 * curGen(indv)._4)
	gcounter += 1
}
curGen = nextGen
nextGen =  new Array[Tuple4[Double,Double,Double,Double]](popsize)
envCount += 1

curGenNumb += 1

}
    }
    for (gens <- 0 to generations){
      out.write(gens + "\t\t\t")
      data(gens).foreach(loop => out.write("\t" + loop))
      out.write("\n")
    }
    out.close
  }
  
}