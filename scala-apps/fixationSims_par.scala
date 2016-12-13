import scala.collection.mutable.HashMap

System.err.println("app.scala <baseRate x.x> <early *X> <popsize> <male2female 1.0> <generations> <mosaic true/false>")
if(args.size == 0) System.exit(0)
val genome = 2000000000
val rand = new scala.util.Random

//val baseRate = 0.4
val multiplier = args(1).toInt
val popsize = args(2).toInt
val replications = 14
val male2female = args(3).toDouble
val generations = args(4).toInt
val mosaic = args(5).toUpperCase.contains("T")
var toRemove : List[Int] = Nil

def mutGen (num: Int): List[Int] = {
	var count = 0
	var result : List[Int] = Nil
		while (count < num) {
			result = rand.nextInt(genome) :: result
			count += 1
		}
	result
}

val baseRate = new org.apache.commons.math3.distribution.PoissonDistribution(args(0).toDouble)
val earlyRate = new org.apache.commons.math3.distribution.PoissonDistribution(args(0).toDouble * multiplier)
val spermRate = new org.apache.commons.math3.distribution.PoissonDistribution(2)
//Arrays of the Mutations carried by each individual
var curPop = Array.fill(popsize)(List(-1))
var nextGen = Array.fill(popsize)(List(-1))
var gametes = new Array[List[List[Int]]] (popsize)
var curGen = 0
val mutAge = new HashMap[Int,Int]
var fixtime : List[Int] = Nil
var mosDegree : collection.mutable.Map[Int,Double] = collection.mutable.Map()
var mosMutFix : List[Int] = Nil
var nmMutFix : List[Int] = Nil
var avgNumMuts: List[Int] = Nil

while (curGen < generations){
	var indv = 0
	//println("Gen Gametes")
	while(indv < popsize){
		var cells: List[List[Int]] = List(mutGen(0))
		var reps = 0
		while (reps < replications) {
			
			if (mosaic){
				//println("mos")
				if (reps < 4){
					cells = cells.map(e => mutGen(earlyRate.sample) ::: e) ::: cells.map(e => mutGen(earlyRate.sample) ::: e)	
				} else {
					if (reps == 13 && indv < (popsize / 2)) cells = cells.map(e => mutGen(spermRate.sample) ::: e) ::: cells.map(e => mutGen(spermRate.sample) ::: e) else cells = cells.map(e => mutGen(baseRate.sample) ::: e) ::: cells.map(e => mutGen(baseRate.sample) ::: e)
				}
				reps += 1
			} else {
				if (reps == (replications - 1)){
					if (indv < (popsize / 2)) cells = cells.map(e => mutGen(spermRate.sample) ::: e) ::: cells.map(e => mutGen(spermRate.sample) ::: e) else cells = cells.map(e => mutGen(baseRate.sample) ::: e) ::: cells.map(e => mutGen(baseRate.sample) ::: e)
				} else {
					cells = cells.map(e => mutGen(baseRate.sample) ::: e) ::: cells.map(e => mutGen(baseRate.sample) ::: e)
				}
				reps += 1

			}


		}

		val mosCount = cells.flatten.groupBy(identity).map(x => (x._1, x._2.size)).map(x => (x._1,x._2/cells.size.toDouble))
		println("avg muts " + cells.map(_.size).sum/cells.size)
		//val tmp : scala.collection.immutable.HashMap[Int,Double] = mosCount.map(x => (x._1,x._2/cells.size.toDouble))
		//val tmpMD = collection.mutable.HashMap() ++ mosCount ++ mosDegree
		val tmpMD = collection.mutable.HashMap() ++ mosCount ++ mosDegree
		mosDegree = tmpMD
		gametes(indv) = cells
		indv += 1
	}

	indv = 0
	//If zero Gen muts = gametes, randomly take a gamete from two individuals and merge
	if (curGen == 0){
		//println( curGen + " 1st Gen Gametes")
			for (indv2 <- (0 to popsize -1).par){
			//while(indv < popsize){
				val father = rand.nextInt(popsize / 2)
				val mother = rand.nextInt(popsize / 2) + (popsize / 2)
				val maleGam = gametes(father)(rand.nextInt(gametes(father).size))
				val femaleGam = gametes(mother)(rand.nextInt(gametes(mother).size))
				nextGen(indv2) = maleGam ::: femaleGam
				//indv += 1
			}
		} else {
			println( curGen + " Gen Gametes")
			for (indv2 <- (0 to popsize -1).par){
			//while(indv < popsize){
				val father = rand.nextInt(popsize / 2)
				val mother = rand.nextInt(popsize / 2) + (popsize / 2)
				val maleGam = (gametes(father)(rand.nextInt(gametes(father).size)) ::: curPop(father).filter(y => rand.nextInt(2) == 1)).filter(x => ! toRemove.contains(x))
				val femaleGam = gametes(mother)(rand.nextInt(gametes(mother).size)) ::: curPop(mother).filter(y => rand.nextInt(2) == 1).filter(x => ! toRemove.contains(x))
				nextGen(indv2) = maleGam ::: femaleGam
				//indv += 1
			}
		}
		toRemove = Nil
		//build list of all muts and number of times seen
		val muts = new HashMap[Int,Int]
		//indv = 0
		//for (indv2 <- (0 to popsize -1).par){
		while(indv < popsize){
			//println("Counting Muts")
			var mutnum = 0
			while (mutnum < nextGen(indv).size){
				if (muts.contains(nextGen(indv)(mutnum))){
						muts(nextGen(indv)(mutnum)) += 1
					} else {
						muts +=  nextGen(indv)(mutnum) -> 1
					}
				mutnum += 1
			}
			indv += 1
		}
		//Check to see if mut is new if not add one gen onto age if is add and set to 1
		for (m <- muts){
			if (mutAge.contains(m._1)){
				mutAge(m._1) += 1
			} else {
				mutAge += m._1 -> 1
			}

			if (m._2 == popsize) {
				println(s"Mutation ${m._1} fixed ${m._2} in ${mutAge(m._1)} was at AR ${mosDegree(m._1)}")
				//println(s"Mutation ${m._1} fixed ${m._2} in ${mutAge(m._1)}")
				if (mosDegree(m._1) >= 0.01) mosMutFix = mutAge(m._1) :: mosMutFix else nmMutFix = mutAge(m._1) :: nmMutFix
				fixtime = mutAge(m._1) :: fixtime
				toRemove = m._1 :: toRemove
				//mosDegree.remove(m._1)
			} else {
				//if ((m._2 / popsize.toDouble) >= 0.9) System.err.println(s"Mutation ${m._1} at ${m._2/popsize.toDouble} in ${mutAge(m._1)}")
			}

		}

		//check to see if mut has reach fixation if so report age and add to ignore

		curPop = nextGen

	curGen += 1
}

System.err.println("Average fix time = " + fixtime.sum/fixtime.size.toDouble + " Num of fixed = " + fixtime.size)
System.err.println("M >= 1% avg fix  = " + mosMutFix.sum/mosMutFix.size.toDouble + " Num of fixed = " + mosMutFix.size.toDouble)
System.err.println("NM <= 1% avg fix = " + nmMutFix.sum/nmMutFix.size.toDouble + " Num of fixed = " + nmMutFix.size.toDouble) 
System.err.println("avg number of muts = " + avgNumMuts.sum/avgNumMuts.size.toDouble )

println(curPop.map(_.size).sum/curPop.size)
